#!/usr/bin/env Rscript
# ==============================================================================
# Analysis / Figure7.R
# ==============================================================================
# PROVENANCE: verbatim relocation of the original scripts/.../09a_sensitivity_analysis.R.
#
# Reproduces submitted Figure 7 ("Sensitivity analysis: temporal burden increase"; ST69 VFDB Cluster_3 / ST10 VF Cluster_3 / ST10 ResFinder Cluster_6). Exact plot-title text match confirmed against the submitted Figure_07.png.
#
# ANALYSIS IS UNCHANGED -- this is the original script's body, unmodified
# (aside from this provenance header and, where noted above, a filename
# correction to match the manuscript's actual published figure numbering).
# ==============================================================================

#!/usr/bin/env Rscript
# Sensitivity analysis: ST69 VFDB Cluster_3 vs ST10 VF Cluster_3
# Output: Excel + figure showing raw trend, clinical composition, model sensitivity
suppressPackageStartupMessages({
  library(tidyverse); library(broom); library(patchwork); library(writexl)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "sensitivity_analysis")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

# ======================================================================
# Load ST69 VFDB data
# ======================================================================
cat("Loading ST69 VFDB...\n")
master69 <- read_csv(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv"), show_col_types = FALSE,
  guess_max = 10000) %>%
  mutate(genome_id = as.character(genome_id), year = as.integer(year),
    clinical_binary = as.integer(clinical_binary),
    source_niche = as.character(source_niche), country = as.character(country),
    total = total_vf) %>%
  filter(!is.na(clinical_binary), !is.na(year), !is.na(shell_cluster))

to_bin_vfdb <- function(x) as.integer(sapply(strsplit(as.character(x), ";"),
  function(v) any(suppressWarnings(as.numeric(v) >= 90), na.rm = TRUE)))

vfdb_raw <- read_tsv(config$st_vfdb_summary("ST69"), show_col_types = FALSE,
  progress = FALSE, col_types = cols(.default = "c")) %>%
  filter(str_detect(.data[["#FILE"]], "^ST69/")) %>%
  mutate(genome_id = make_genome_id(str_remove(.data[["#FILE"]], "^ST69/")),
         genome_id = str_remove(genome_id, "_vfdb\\.tsv$"))
vf_genes69 <- setdiff(colnames(vfdb_raw), c("#FILE", "NUM_FOUND", "genome_id"))
vfdb_bin69 <- vfdb_raw %>% select(genome_id, all_of(vf_genes69)) %>%
  mutate(across(all_of(vf_genes69), to_bin_vfdb))

df69 <- master69 %>% select(genome_id, shell_cluster, year, clinical_binary,
                             source_niche, country, total) %>%
  inner_join(vfdb_bin69, by = "genome_id")
cd69 <- df69 %>% filter(shell_cluster == "Cluster_3") %>%
  filter(!is.na(source_niche), !is.na(country))
cat("ST69 Cluster_3:", nrow(cd69), "genomes\n")

# ======================================================================
# Load ST10 VF data
# ======================================================================
cat("Loading ST10 VF...\n")
master10 <- read_csv(file.path(config$OUTPUT_DIR, "ST10", "virulencefinder_validation",
  "04_master_shell_cluster_metadata_VF_table.csv"), show_col_types = FALSE) %>%
  mutate(genome_id = as.character(genome_id), year = as.integer(year),
    clinical_binary = as.integer(clinical_binary),
    source_niche = as.character(source_niche), country = as.character(country)) %>%
  filter(!is.na(clinical_binary), !is.na(year), !is.na(shell_cluster))

vfb10 <- read.delim(config$VF_BINARY, header = TRUE, sep = "\t",
  check.names = FALSE, na.strings = "", stringsAsFactors = FALSE) %>%
  filter(st == "ST10") %>%
  mutate(genome_id = make_genome_id(genome)) %>% select(-st, -genome)

genes_vf10 <- setdiff(colnames(vfb10), "genome_id")
for (g in genes_vf10)
  vfb10[[g]] <- as.integer(!is.na(vfb10[[g]]) & vfb10[[g]] != "." & vfb10[[g]] != "0" & vfb10[[g]] != "")
vfb10 <- vfb10 %>% mutate(total = rowSums(across(all_of(genes_vf10))))

df10 <- master10 %>% inner_join(vfb10, by = "genome_id")
cd10 <- df10 %>% filter(shell_cluster == "Cluster_3") %>%
  filter(!is.na(source_niche), !is.na(country))
cat("ST10 Cluster_3:", nrow(cd10), "genomes\n")

# ======================================================================
# Load ST10 ResFinder Cluster_6 data
# ======================================================================
cat("Loading ST10 ResFinder...\n")
rf_raw <- read_tsv("finder_result/resfinder_summary/resfinder_binary_matrix.tsv",
  show_col_types = FALSE) %>%
  filter(st == "ST10") %>%
  mutate(genome_id = as.character(genome))
st10_vfdb <- read_csv(file.path(config$OUTPUT_DIR, "ST10", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv"),
  show_col_types = FALSE) %>%
  select(genome_id, shell_cluster, year, clinical_binary, source_niche, country)
rf_df <- rf_raw %>% inner_join(st10_vfdb, by = "genome_id")
gene_cols <- setdiff(colnames(rf_raw), c("genome_id", "st", "genome", "X.FILE", "X"))
for (g in gene_cols) rf_df[[g]] <- as.integer(!is.na(rf_df[[g]]) & rf_df[[g]] != "0" & rf_df[[g]] != "")
rf_df <- rf_df %>% mutate(total = rowSums(across(all_of(gene_cols))))
cd6 <- rf_df %>% filter(shell_cluster == "Cluster_6") %>%
  filter(!is.na(source_niche), !is.na(country))
cat("ST10 ResFinder Cluster_6:", nrow(cd6), "genomes\n")

# ======================================================================
# Model comparison function
# ======================================================================
run_models <- function(data, label) {
  m1 <- lm(total ~ year, data = data)
  m2 <- lm(total ~ year + country, data = data)
  m3 <- lm(total ~ year + source_niche, data = data)
  m4 <- lm(total ~ year + source_niche + country, data = data)

  extract <- function(m, model_name) {
    s <- coef(summary(m))
    ci <- confint(m)
    r <- rownames(s)[rownames(s) == "year"]
    if (length(r) == 0) return(NULL)
    tibble(dataset = label, model = model_name,
           slope = s[r, "Estimate"],
           se = s[r, "Std. Error"],
           ci_lower = ci[r, "2.5 %"],
           ci_upper = ci[r, "97.5 %"],
           p_value = s[r, "Pr(>|t|)"],
           n = nobs(m))
  }

  bind_rows(
    extract(m1, "year-only"),
    extract(m2, "+ country"),
    extract(m3, "+ niche"),
    extract(m4, "+ niche + country")
  )
}

# ======================================================================
# Sensitivity models (stepwise + stratified)
# ======================================================================
sensitivity <- bind_rows(
  run_models(cd69, "ST69 VFDB Cluster_3"),
  run_models(cd10, "ST10 VF Cluster_3"),
  run_models(cd6, "ST10 ResFinder Cluster_6")
)

# Niche-stratified
niche_trends <- bind_rows(
  cd69 %>% group_by(source_niche) %>%
    group_modify(~ {
      if (nrow(.x) < 10 || n_distinct(.x$year) < 3)
        return(tibble(dataset = "ST69 VFDB", slope = NA, p = NA, n = nrow(.x), mean_vf = mean(.x$total)))
      m <- lm(total ~ year, data = .x)
      s <- coef(summary(m))
      tibble(dataset = "ST69 VFDB", slope = s["year","Estimate"], p = s["year","Pr(>|t|)"],
             n = nrow(.x), mean_vf = mean(.x$total))
    }),
  cd10 %>% group_by(source_niche) %>%
    group_modify(~ {
      if (nrow(.x) < 10 || n_distinct(.x$year) < 3)
        return(tibble(dataset = "ST10 VF", slope = NA, p = NA, n = nrow(.x), mean_vf = mean(.x$total)))
      m <- lm(total ~ year, data = .x)
      s <- coef(summary(m))
      tibble(dataset = "ST10 VF", slope = s["year","Estimate"], p = s["year","Pr(>|t|)"],
             n = nrow(.x), mean_vf = mean(.x$total))
    }),
  cd6 %>% group_by(source_niche) %>%
    group_modify(~ {
      if (nrow(.x) < 10 || n_distinct(.x$year) < 3)
        return(tibble(dataset = "ST10 ResF C6", slope = NA, p = NA, n = nrow(.x), mean_vf = mean(.x$total)))
      m <- lm(total ~ year, data = .x)
      s <- coef(summary(m))
      tibble(dataset = "ST10 ResF C6", slope = s["year","Estimate"], p = s["year","Pr(>|t|)"],
             n = nrow(.x), mean_vf = mean(.x$total))
    })
)

# Country-stratified (top 8 by n)
top_ctry69 <- cd69 %>% count(country, sort = TRUE) %>% slice_max(n, n = 8) %>% pull(country)
top_ctry10 <- cd10 %>% count(country, sort = TRUE) %>% slice_max(n, n = 8) %>% pull(country)
top_ctry6 <- cd6 %>% count(country, sort = TRUE) %>% slice_max(n, n = 8) %>% pull(country)

country_trends <- bind_rows(
  cd69 %>% filter(country %in% top_ctry69) %>%
    group_by(country) %>%
    group_modify(~ {
      if (nrow(.x) < 10 || n_distinct(.x$year) < 3)
        return(tibble(dataset = "ST69 VFDB", slope = NA, p = NA, n = nrow(.x)))
      m <- lm(total ~ year, data = .x)
      s <- coef(summary(m))
      tibble(dataset = "ST69 VFDB", slope = s["year","Estimate"], p = s["year","Pr(>|t|)"],
             n = nrow(.x))
    }),
  cd10 %>% filter(country %in% top_ctry10) %>%
    group_by(country) %>%
    group_modify(~ {
      if (nrow(.x) < 10 || n_distinct(.x$year) < 3)
        return(tibble(dataset = "ST10 VF", slope = NA, p = NA, n = nrow(.x)))
      m <- lm(total ~ year, data = .x)
      s <- coef(summary(m))
      tibble(dataset = "ST10 VF", slope = s["year","Estimate"], p = s["year","Pr(>|t|)"],
             n = nrow(.x))
    }),
  cd6 %>% filter(country %in% top_ctry6) %>%
    group_by(country) %>%
    group_modify(~ {
      if (nrow(.x) < 10 || n_distinct(.x$year) < 3)
        return(tibble(dataset = "ST10 ResF C6", slope = NA, p = NA, n = nrow(.x)))
      m <- lm(total ~ year, data = .x)
      s <- coef(summary(m))
      tibble(dataset = "ST10 ResF C6", slope = s["year","Estimate"], p = s["year","Pr(>|t|)"],
             n = nrow(.x))
    })
)

# Yearly composition
yearly_comp <- bind_rows(
  cd69 %>% group_by(year) %>% summarise(
    dataset = "ST69 VFDB Cluster_3", n = n(),
    mean_vf = mean(total), sd_vf = sd(total),
    pct_clinical = mean(clinical_binary) * 100,
    pct_Human = mean(source_niche == "Human") * 100,
    n_countries = n_distinct(country), .groups = "drop"),
  cd10 %>% group_by(year) %>% summarise(
    dataset = "ST10 VF Cluster_3", n = n(),
    mean_vf = mean(total), sd_vf = sd(total),
    pct_clinical = mean(clinical_binary) * 100,
    pct_Human = mean(source_niche == "Human") * 100,
    n_countries = n_distinct(country), .groups = "drop"),
  cd6 %>% group_by(year) %>% summarise(
    dataset = "ST10 ResF Cluster_6", n = n(),
    mean_vf = mean(total), sd_vf = sd(total),
    pct_clinical = mean(clinical_binary) * 100,
    pct_Human = mean(source_niche == "Human") * 100,
    n_countries = n_distinct(country), .groups = "drop")
)

# Interaction models (clinical ~ year * is_target + country)
interaction_models <- bind_rows(
  {
    df69_model <- df69 %>% filter(!is.na(source_niche), !is.na(country)) %>%
      mutate(is_target = as.integer(shell_cluster == "Cluster_3"))
    m <- glm(clinical_binary ~ year * is_target + country,
             data = df69_model, family = binomial())
    tidy(m) %>% filter(str_detect(term, ":")) %>%
      mutate(dataset = "ST69 VFDB", model = "clinical ~ year * target + country")
  },
  {
    df10_model <- df10 %>% filter(!is.na(source_niche), !is.na(country)) %>%
      mutate(is_target = as.integer(shell_cluster == "Cluster_3"))
    m <- glm(clinical_binary ~ year * is_target + country,
             data = df10_model, family = binomial())
    tidy(m) %>% filter(str_detect(term, ":")) %>%
      mutate(dataset = "ST10 VF", model = "clinical ~ year * target + country")
  },
  {
    df6_model <- rf_df %>% filter(!is.na(source_niche), !is.na(country)) %>%
      mutate(is_target = as.integer(shell_cluster == "Cluster_6"))
    m <- glm(clinical_binary ~ year * is_target + country,
             data = df6_model, family = binomial())
    tidy(m) %>% filter(str_detect(term, ":")) %>%
      mutate(dataset = "ST10 ResF", model = "clinical ~ year * target + country")
  }
)

# Adjusted model full output
adjusted_full <- bind_rows(
  {
    m <- lm(total ~ year + source_niche + country, data = cd69)
    tidy(m) %>% mutate(dataset = "ST69 VFDB Cluster_3", model = "total ~ year + niche + country")
  },
  {
    m <- lm(total ~ year + source_niche + country, data = cd10)
    tidy(m) %>% mutate(dataset = "ST10 VF Cluster_3", model = "total ~ year + niche + country")
  },
  {
    m <- lm(total ~ year + source_niche + country, data = cd6)
    tidy(m) %>% mutate(dataset = "ST10 ResF Cluster_6", model = "total ~ year + niche + country")
  }
)

# ======================================================================
# Write Excel
# ======================================================================
cat("Writing Excel...\n")
list(
  Stepwise_model_comparison = sensitivity,
  Niche_stratified_trends = niche_trends,
  Country_stratified_trends = country_trends,
  Yearly_composition = yearly_comp,
  Interaction_models = interaction_models,
  Adjusted_model_full = adjusted_full
) %>% write_xlsx(file.path(OUT, "sensitivity_analysis.xlsx"))
cat("  ->", file.path(OUT, "sensitivity_analysis.xlsx"), "\n")

# ======================================================================
# Figure: 3 columns (ST69 | ST10 VF | ST10 ResF) x 3 rows
# ======================================================================
cat("Making figure...\n")

# Row 1: Yearly mean burden
p1a <- cd69 %>% group_by(year) %>%
  summarise(mean_vf = mean(total), se = sd(total)/sqrt(n()), .groups = "drop") %>%
  ggplot(aes(x = year, y = mean_vf)) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "grey80", alpha = 0.3) +
  geom_point(size = 2.5) + geom_line(linetype = "dashed", alpha = 0.5) +
  labs(x = NULL, y = "Mean burden", title = "ST69 VFDB Cluster_3") +
  theme_minimal(base_size = 10)

p1b <- cd10 %>% group_by(year) %>%
  summarise(mean_vf = mean(total), se = sd(total)/sqrt(n()), .groups = "drop") %>%
  ggplot(aes(x = year, y = mean_vf)) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "grey80", alpha = 0.3) +
  geom_point(size = 2.5) + geom_line(linetype = "dashed", alpha = 0.5) +
  labs(x = NULL, y = "Mean burden", title = "ST10 VF Cluster_3") +
  theme_minimal(base_size = 10)

p1c <- cd6 %>% group_by(year) %>%
  summarise(mean_vf = mean(total), se = sd(total)/sqrt(n()), .groups = "drop") %>%
  ggplot(aes(x = year, y = mean_vf)) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "grey80", alpha = 0.3) +
  geom_point(size = 2.5) + geom_line(linetype = "dashed", alpha = 0.5) +
  labs(x = NULL, y = "Mean burden", title = "ST10 ResF Cluster_6") +
  theme_minimal(base_size = 10)

# Row 2: Yearly clinical %
p2a <- cd69 %>% group_by(year) %>%
  summarise(pct_clinical = mean(clinical_binary) * 100, .groups = "drop") %>%
  ggplot(aes(x = year, y = pct_clinical)) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", fill = "grey80", alpha = 0.3) +
  geom_point(size = 2.5) + geom_line(linetype = "dashed", alpha = 0.5) +
  labs(x = NULL, y = "Clinical (%)") +
  theme_minimal(base_size = 10)

p2b <- cd10 %>% group_by(year) %>%
  summarise(pct_clinical = mean(clinical_binary) * 100, .groups = "drop") %>%
  ggplot(aes(x = year, y = pct_clinical)) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", fill = "grey80", alpha = 0.3) +
  geom_point(size = 2.5) + geom_line(linetype = "dashed", alpha = 0.5) +
  labs(x = NULL, y = "Clinical (%)") +
  theme_minimal(base_size = 10)

p2c <- cd6 %>% group_by(year) %>%
  summarise(pct_clinical = mean(clinical_binary) * 100, .groups = "drop") %>%
  ggplot(aes(x = year, y = pct_clinical)) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", fill = "grey80", alpha = 0.3) +
  geom_point(size = 2.5) + geom_line(linetype = "dashed", alpha = 0.5) +
  labs(x = NULL, y = "Clinical (%)") +
  theme_minimal(base_size = 10)

# Row 3: Forest plot of year coefficient across models
p3_data <- sensitivity %>%
  mutate(dataset = factor(dataset,
    levels = c("ST69 VFDB Cluster_3", "ST10 VF Cluster_3", "ST10 ResFinder Cluster_6")),
    model = factor(model, levels = c("year-only", "+ country", "+ niche", "+ niche + country")))

p3 <- p3_data %>%
  ggplot(aes(x = slope, y = model, xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(size = 0.8, fatten = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey40"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "ns"),
    name = NULL) +
  labs(x = "Year slope (95% CI)", y = NULL) +
  facet_wrap(~ dataset, ncol = 3, scales = "free_x") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 9))

# Assemble
fig <- (p1a | p1b | p1c) / (p2a | p2b | p2c) / p3 +
  plot_annotation(title = "Sensitivity analysis: temporal burden increase")

ggsave(file.path(OUT, "sensitivity_figure.png"), fig, width = 8, height = 10, dpi = 300)
ggsave(file.path(OUT, "sensitivity_figure.pdf"), fig, width = 8, height = 10)
cat("  ->", file.path(OUT, "sensitivity_figure.png"), "\n")

# ======================================================================
cat("Done.\n")
