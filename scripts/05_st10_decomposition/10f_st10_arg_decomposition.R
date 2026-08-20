#!/usr/bin/env Rscript
# Reviewer response: decompose ST10 ARG decreasing trend
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(writexl); library(data.table); library(broom)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST10", "reviewer_arg_decomposition")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
    grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
    TRUE ~ paste0("Escherichia_coli_", x))
}

# ---- 1. Load ST10 CARD data ----
cat("Loading ST10 CARD...\n")
card_raw <- fread(config$st_card_burden("ST10"), sep = "\t", header = TRUE,
  data.table = FALSE, colClasses = "character", na.strings = "")
card_raw$genome_id <- make_genome_id(str_remove(card_raw[[1]], "^card/ST10/"))
card_raw$genome_id <- str_remove(card_raw$genome_id, "_card\\.tsv$")

# Count ARGs per genome
arg_genes <- setdiff(colnames(card_raw), c("#FILE", "NUM_FOUND", "genome_id"))
to_bin <- function(x) as.integer(sapply(strsplit(as.character(x), ";"),
  function(v) any(suppressWarnings(as.numeric(v) >= 90), na.rm = TRUE)))
card_bin <- card_raw %>% select(genome_id, all_of(arg_genes)) %>%
  mutate(across(all_of(arg_genes), to_bin))
card_bin$total_arg <- rowSums(card_bin[, arg_genes, drop = FALSE], na.rm = TRUE)

cat("ST10 CARD genomes:", nrow(card_bin), "\n")

# ---- 2. Load ST10 metadata ----
cat("Loading ST10 metadata...\n")
meta_file <- config$st_metadata("ST10")
meta <- read_excel(meta_file, col_types = "text")

# Find key columns
name_col <- {
  candidates <- c("Name", "genome", "strain", "isolate", "assembly", "sample")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[2]
}
year_col <- {
  candidates <- c("Collection Year", "Collection_Year", "year", "Year", "date")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[7]
}
source_col <- {
  candidates <- c("Source Niche", "Source_Niche", "source_niche", "source", "niche", "Host", "host")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[4]
}
country_col <- {
  candidates <- c("Country", "country")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[12]
}
continent_col <- {
  candidates <- c("Continent", "continent")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else NA_character_
}

meta_clean <- meta %>%
  transmute(
    name_raw = .data[[name_col]],
    year = suppressWarnings(as.integer(.data[[year_col]])),
    source_niche = as.character(.data[[source_col]]),
    country = as.character(.data[[country_col]]),
    continent = if (!is.na(continent_col)) as.character(.data[[continent_col]]) else NA_character_
  ) %>%
  mutate(genome_id = make_genome_id(name_raw))

# ---- 3. Merge ----
merged <- card_bin %>% inner_join(meta_clean, by = "genome_id") %>%
  filter(!is.na(year), year >= 2016, year <= 2025)

cat("Merged genomes:", nrow(merged), "\n")

# Classify niche
merged <- merged %>%
  mutate(
    clinical = case_when(
      grepl("human|patient|clinical|blood|urine|wound|hospital|csf|respiratory|sputum",
        source_niche, ignore.case = TRUE) ~ "Clinical",
      is.na(source_niche) | source_niche == "" ~ "Unknown",
      TRUE ~ "Non-clinical"
    ),
    region = case_when(
      grepl("Europe|UK|Germany|France|Spain|Italy|Netherlands|Denmark|Sweden|Norway|Switzerland|Austria|Belgium|Ireland|Portugal|Poland|Czech|Hungary|Romania|Greece|Finland",
        country, ignore.case = TRUE) ~ "Europe",
      grepl("USA|United States|Canada|Mexico",
        country, ignore.case = TRUE) ~ "North America",
      grepl("China|Japan|India|South Korea|Taiwan|Thailand|Vietnam|Indonesia|Malaysia|Philippines|Singapore|Hong Kong|Pakistan|Bangladesh|Sri Lanka",
        country, ignore.case = TRUE) ~ "Asia",
      grepl("Australia|New Zealand",
        country, ignore.case = TRUE) ~ "Oceania",
      grepl("Brazil|Argentina|Chile|Colombia|Peru|Uruguay|Paraguay|Ecuador|Venezuela",
        country, ignore.case = TRUE) ~ "South America",
      grepl("South Africa|Nigeria|Kenya|Ghana|Ethiopia|Tanzania|Uganda|Zambia|Malawi|Mozambique|Botswana|Senegal|Burkina Faso|Benin|Zimbabwe",
        country, ignore.case = TRUE) ~ "Africa",
      TRUE ~ "Other"
    )
  )

# ---- 4. Annual ARG burden trends ----
annual <- merged %>%
  group_by(year) %>%
  summarise(
    n = n(),
    mean_arg = mean(total_arg),
    se_arg = sd(total_arg) / sqrt(n()),
    pct_clinical = mean(clinical == "Clinical") * 100,
    n_regions = n_distinct(region),
    .groups = "drop"
  )
cat("\nAnnual ARG burden:\n")
print(annual, n = Inf)

# ---- 5. Annual by clinical status ----
annual_clinical <- merged %>%
  filter(clinical != "Unknown") %>%
  group_by(year, clinical) %>%
  summarise(
    n = n(),
    mean_arg = mean(total_arg),
    .groups = "drop"
  )
cat("\nAnnual ARG by clinical status:\n")
print(annual_clinical %>% arrange(year, clinical), n = Inf)

# ---- 6. Annual by region ----
annual_region <- merged %>%
  filter(region != "Other") %>%
  group_by(year, region) %>%
  summarise(
    n = n(),
    mean_arg = mean(total_arg),
    .groups = "drop"
  )
cat("\nAnnual ARG by region:\n")
print(annual_region %>% arrange(year, region), n = Inf)

# ---- 7. Test: does clinical % decrease over time? ----
clinical_trend <- merged %>%
  filter(clinical != "Unknown") %>%
  group_by(year) %>%
  summarise(pct_clinical = mean(clinical == "Clinical") * 100, .groups = "drop")
m_clin <- lm(pct_clinical ~ year, data = clinical_trend)
cat("\nClinical % trend:\n")
cat(sprintf("  slope = %.4f, p = %.4e\n", coef(m_clin)[2], summary(m_clin)$coefficients[2,4]))

# ---- 8. Test: does regional composition shift over time? ----
regional_trend <- merged %>%
  filter(region != "Other") %>%
  count(year, region) %>%
  group_by(year) %>%
  mutate(pct = n / sum(n) * 100)

# ---- 9. Adjusted models ----
# Model 1: raw trend
m1 <- lm(total_arg ~ year, data = merged)
# Model 2: adjusted for clinical
m2 <- lm(total_arg ~ year + clinical, data = merged %>% filter(clinical != "Unknown"))
# Model 3: adjusted for region
m3 <- lm(total_arg ~ year + region, data = merged %>% filter(region != "Other"))
# Model 4: adjusted for both
m4 <- lm(total_arg ~ year + clinical + region, data = merged %>%
    filter(clinical != "Unknown", region != "Other"))

model_comp <- bind_rows(
  tidy(m1) %>% filter(term == "year") %>% mutate(model = "raw", dataset = "all"),
  tidy(m2) %>% filter(term == "year") %>% mutate(model = "+ clinical", dataset = "known clinical"),
  tidy(m3) %>% filter(term == "year") %>% mutate(model = "+ region", dataset = "known region"),
  tidy(m4) %>% filter(term == "year") %>% mutate(model = "+ clinical + region", dataset = "both")
) %>% select(model, estimate, std.error, p.value)

cat("\nModel comparison (year coefficient):\n")
print(model_comp)

# ---- 10. Clinical-stratified trends ----
niche_stratified <- merged %>%
  filter(clinical != "Unknown") %>%
  group_by(clinical) %>%
  summarise(
    n = n(),
    slope = tryCatch(coef(lm(total_arg ~ year))[2], error = function(e) NA_real_),
    p = tryCatch(summary(lm(total_arg ~ year))$coefficients[2,4], error = function(e) NA_real_),
    mean_arg = mean(total_arg),
    .groups = "drop"
  )
cat("\nClinical-stratified trends:\n")
print(niche_stratified)

# ---- 11. Region-stratified trends ----
region_stratified <- merged %>%
  filter(region != "Other") %>%
  group_by(region) %>%
  summarise(
    n = n(),
    slope = tryCatch(coef(lm(total_arg ~ year))[2], error = function(e) NA_real_),
    p = tryCatch(summary(lm(total_arg ~ year))$coefficients[2,4], error = function(e) NA_real_),
    mean_arg = mean(total_arg),
    .groups = "drop"
  )
cat("\nRegion-stratified trends:\n")
print(region_stratified)

# ---- 12. Write outputs ----
write_xlsx(list(
  annual_arg = annual,
  annual_by_clinical = annual_clinical %>% arrange(year, clinical),
  annual_by_region = annual_region %>% arrange(year, region),
  clinical_trend = clinical_trend,
  model_comparison = model_comp,
  niche_stratified = niche_stratified,
  region_stratified = region_stratified,
  regional_composition = regional_trend %>% arrange(year, region)
), path = file.path(OUT, "ST10_arg_decomposition.xlsx"))

# ---- 13. Summary ----
sink(file.path(OUT, "ST10_arg_decomposition_summary.txt"))
cat("ST10 ARG Burden Decomposition\n")
cat("==============================\n\n")
cat(sprintf("Genomes: %d (2016-2025)\n", nrow(merged)))
cat(sprintf("Raw ARG trend slope: %.4f (p = %.4e)\n", coef(m1)[2], summary(m1)$coefficients[2,4]))
cat(sprintf("\nAfter adjusting for clinical status: slope = %.4f (p = %.4e)\n",
  filter(model_comp, model == "+ clinical")$estimate,
  filter(model_comp, model == "+ clinical")$p.value))
cat(sprintf("After adjusting for region: slope = %.4f (p = %.4e)\n",
  filter(model_comp, model == "+ region")$estimate,
  filter(model_comp, model == "+ region")$p.value))
cat(sprintf("After adjusting for both: slope = %.4f (p = %.4e)\n",
  filter(model_comp, model == "+ clinical + region")$estimate,
  filter(model_comp, model == "+ clinical + region")$p.value))
cat("\nClinical % trend:\n")
cat(sprintf("  slope = %.4f, p = %.4e\n", coef(m_clin)[2], summary(m_clin)$coefficients[2,4]))
cat("\nInterpretation:\n")
if (summary(m_clin)$coefficients[2,4] < 0.05) {
  cat("  Clinical proportion is changing significantly over time, which may contribute\n")
  cat("  to the ARG burden trend if clinical isolates carry more ARGs.\n")
} else {
  cat("  Clinical proportion is stable. The ARG decrease is not driven by\n")
  cat("  shifting clinical representation.\n")
}
sink()

# ---- 14. Plots ----
p1 <- ggplot(annual, aes(x = year, y = mean_arg)) +
  geom_line() + geom_point(aes(size = n)) +
  geom_ribbon(aes(ymin = mean_arg - se_arg, ymax = mean_arg + se_arg), alpha = 0.2) +
  theme_minimal() + labs(title = "ST10 Mean ARG Burden (2016-2025)", y = "Mean ARG count", x = "")
ggsave(file.path(OUT, "arg_trend_overall.pdf"), p1, width = 6, height = 4)

p2 <- ggplot(annual_clinical %>% filter(clinical != "Unknown"), aes(x = year, y = mean_arg, color = clinical)) +
  geom_line() + geom_point(aes(shape = clinical)) + theme_minimal() +
  labs(title = "ARG Trend by Clinical Status", y = "Mean ARG count", x = "", color = "", shape = "")
ggsave(file.path(OUT, "arg_trend_by_clinical.pdf"), p2, width = 7, height = 4)

p3 <- ggplot(annual_region %>% filter(region != "Other", n >= 10), aes(x = year, y = mean_arg, color = region)) +
  geom_line() + geom_point() + theme_minimal() +
  labs(title = "ARG Trend by Region", y = "Mean ARG count", x = "", color = "")
ggsave(file.path(OUT, "arg_trend_by_region.pdf"), p3, width = 8, height = 5)

p4 <- ggplot(clinical_trend, aes(x = year, y = pct_clinical)) +
  geom_line() + geom_point() + theme_minimal() +
  labs(title = "Clinical Proportion Over Time", y = "% Clinical", x = "")
ggsave(file.path(OUT, "clinical_proportion.pdf"), p4, width = 6, height = 4)

cat("Plots saved.\n")
cat("Done. Output in:", OUT, "\n")
