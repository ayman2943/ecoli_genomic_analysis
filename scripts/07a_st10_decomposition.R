#!/usr/bin/env Rscript
# Decompose ST10 Cluster_3 VF increase: composition vs within-stratum
suppressPackageStartupMessages(library(tidyverse))
source("config.R"); config$TARGET_ST <- "ST10"; source("config.R", local=TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

master <- read_csv(file.path(config$OUTPUT_DIR, "ST10", "virulencefinder_validation",
  "04_master_shell_cluster_metadata_VF_table.csv"), show_col_types=FALSE) %>%
  mutate(genome_id=as.character(genome_id), year=as.integer(year),
    clinical_binary=as.integer(clinical_binary),
    source_niche=as.character(source_niche), country=as.character(country)) %>%
  filter(!is.na(clinical_binary),!is.na(year),!is.na(shell_cluster))

vfb <- read.delim(config$VF_BINARY, header=TRUE, sep="\t", check.names=FALSE,
  na.strings="", stringsAsFactors=FALSE) %>%
  filter(st=="ST10") %>% mutate(genome_id=make_genome_id(genome)) %>% select(-st,-genome)
for (g in setdiff(colnames(vfb),"genome_id"))
  vfb[[g]] <- as.integer(!is.na(vfb[[g]])&vfb[[g]]!="."&vfb[[g]]!="0"&vfb[[g]]!="")
vfb <- vfb %>% mutate(total=rowSums(across(all_of(setdiff(colnames(vfb),"genome_id")))))

cd <- master %>% inner_join(vfb, by="genome_id") %>%
  filter(shell_cluster=="Cluster_3",!is.na(source_niche),!is.na(country))

# Create stratum = country | niche (collapse small groups)
cd <- cd %>% mutate(stratum = paste0(country, " | ", source_niche))

# 1. Yearly composition (proportion of each stratum)
comp <- cd %>% count(year, stratum) %>% group_by(year) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup()

# 2. Each stratum's overall mean VF (fixed reference)
stratum_means <- cd %>% group_by(stratum) %>%
  summarise(mean_vf = mean(total), total_n = n(), .groups = "drop")

# Join composition with stratum means
comp <- comp %>% left_join(stratum_means, by = "stratum")

# 3. Standardization: expected VF if only composition changed, but each stratum's VF stayed fixed
# For each stratum, compute overall mean VF (pooled across all years) — the "fixed" within value
stratum_fixed <- cd %>% group_by(stratum) %>%
  summarise(fixed_vf = mean(total), total_n = n(), .groups = "drop")

# Expected yearly VF = Σ(stratum_proportion_year × fixed_vf) for ALL strata in that year
expected <- comp %>% left_join(stratum_fixed, by = "stratum") %>%
  group_by(year) %>%
  summarise(expected = weighted.mean(fixed_vf, n), .groups = "drop")

# Actual yearly trend
actual <- cd %>% group_by(year) %>%
  summarise(actual = mean(total), n = n(), .groups = "drop")

compare <- actual %>% left_join(expected, by = "year")

# ===== Figure =====
# Panel A: stacked bar of composition
top_strata <- stratum_means %>% slice_max(total_n, n = 8) %>% pull(stratum)
cd_plot <- cd %>% mutate(stratum_plot = ifelse(stratum %in% top_strata, stratum, "Other"),
  stratum_plot = factor(stratum_plot, levels = c(sort(top_strata), "Other")))

p1 <- cd_plot %>% count(year, stratum_plot) %>% group_by(year) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup() %>%
  ggplot(aes(x = year, y = pct, fill = stratum_plot)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = NULL, y = "Proportion (%)", title = "Composition of ST10 Cluster_3 by year") +
  theme_minimal(base_size = 10) +
  theme(legend.text = element_text(size = 7))

# Panel B: within-stratum yearly VF (strata with >= 10 genomes)
big_strata <- stratum_means %>% filter(total_n >= 10) %>% pull(stratum)
pd <- cd %>% filter(stratum %in% big_strata) %>%
  group_by(year, stratum) %>%
  summarise(mean_vf = mean(total), se = sd(total)/sqrt(n()), n = n(), .groups = "drop")
pd_clean <- pd %>% filter(n >= 3)  # drop years with <3 genomes for plotting

p2 <- pd_clean %>%
  ggplot(aes(x = year, y = mean_vf, color = stratum)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, alpha = 0.5) +
  scale_color_brewer(palette = "Set2", name = NULL) +
  labs(x = NULL, y = "Mean VF burden",
       title = "Within-stratum VF trends (strata with >= 10 genomes)") +
  theme_minimal(base_size = 10) +
  theme(legend.text = element_text(size = 7))

# Panel C: actual vs expected (composition-only)
p3 <- compare %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = actual, color = "Actual"), linewidth = 1) +
  geom_point(aes(y = actual, color = "Actual"), size = 2.5) +
  geom_line(aes(y = expected, color = "Expected (composition only)"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = expected, color = "Expected (composition only)"),
             size = 2.5, shape = 1) +
  scale_color_manual(values = c("Actual" = "red",
    "Expected (composition only)" = "steelblue"), name = NULL) +
  labs(x = "Year", y = "Mean VF burden",
       title = "Actual vs composition-only expected VF") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")

# Assemble
library(patchwork)
fig <- (p1 / p2 / p3) + plot_annotation(
       title = "ST10 VF Cluster_3: compositional analysis",
       subtitle = "If actual and expected lines overlap, the increase is entirely compositional")

outdir <- file.path(config$OUTPUT_DIR, "sensitivity_analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(outdir, "decomposition_compositional.png"), fig, width = 9, height = 12, dpi = 300)
cat("Saved:", file.path(outdir, "decomposition_compositional.png"), "\n")

# Print summary
cat("\n=== Actual vs Standardized (composition-adjusted) ===\n")
print(compare, n = 20)

# LM on actual
m_actual <- lm(actual ~ year, data = compare)
cat(sprintf("\nActual trend: slope=%.4f  p=%.4f\n",
  coef(m_actual)["year"], coef(summary(m_actual))["year","Pr(>|t|)"]))

# LM on expected (should match actual if composition fully explains)
m_exp <- lm(expected ~ year, data = compare)
cat(sprintf("Expected (composition-only) trend: slope=%.4f  p=%.4f\n",
  coef(m_exp)["year"], coef(summary(m_exp))["year","Pr(>|t|)"]))
cat(sprintf("Actual - Expected difference in slope: %.4f\n",
  coef(m_actual)["year"] - coef(m_exp)["year"]))
