#!/usr/bin/env Rscript
#
# VirulenceFinder.R  —  Shell-gene cluster virulence analysis (VirulenceFinder track)
#
# Split out of the original virulence_analysis.R (which ran the VFDB and
# VirulenceFinder tracks back-to-back). This file keeps only the VF half;
# see VFDB.R for the VFDB half. Both files duplicate the shared
# `analyze_vf_track()` function and setup code rather than sourcing each
# other, matching how every script in this repo only ever sources config.R.
#
# Maps VirulenceFinder to shell-gene clusters, detects temporal trends,
# selects the cluster with strongest VF signal, and performs clinical
# enrichment + figures.
#
# Usage:
#   TARGET_ST=ST10 Rscript scripts/Virulence/VirulenceFinder.R
#
suppressPackageStartupMessages({
  library(tidyverse); library(broom); library(writexl); library(readr)
  library(ggplot2); library(scales); library(patchwork); library(forcats)
})
ts <- function(m) cat("[", Sys.time(), "] ", m, "\n", sep = "")

source("config.R")
st <- config$TARGET_ST
OUT_DIR <- file.path(config$OUTPUT_DIR, st, "virulence_analysis")
FIG_DIR <- file.path(OUT_DIR, "figures"); TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
EARLY_YEARS <- c(2016, 2017, 2018); LATE_YEARS <- c(2022, 2023, 2024, 2025)

# -------------------------------------------------------------------
# 1. Load master table (cluster assignments, metadata)
# -------------------------------------------------------------------
master_file <- file.path(config$OUTPUT_DIR, st, "vfdb_analysis",
  grep("04_master.*\\.csv$", list.files(file.path(config$OUTPUT_DIR, st, "vfdb_analysis"), pattern = "04_master"), value = TRUE)[1])
if (is.na(master_file)) stop("Master table not found")
master <- read_csv(master_file, show_col_types = FALSE) %>%
  mutate(
    genome_id = as.character(genome_id), year = as.integer(year),
    clinical_binary = case_when(niche == "clinical" ~ 1L, niche == "non-clinical" ~ 0L, TRUE ~ NA_integer_),
    source_niche = as.character(source_niche), country = as.character(country)
  ) %>%
  filter(!is.na(clinical_binary), !is.na(year), !is.na(shell_cluster))

common <- function(df) {
  df %>% left_join(master %>% select(genome_id, shell_cluster, year, clinical_binary, source_niche, country), by = "genome_id") %>%
    filter(!is.na(clinical_binary), !is.na(year), !is.na(shell_cluster))
}

# -------------------------------------------------------------------
# 2. Analysis function (VirulenceFinder track)
# -------------------------------------------------------------------
analyze_vf_track <- function(label, df) {
  ts(paste("\n===========", label, "==========="))
  genes <- setdiff(colnames(df), c("genome_id", "shell_cluster", "year", "clinical_binary", "source_niche", "country"))
  df <- df %>% mutate(total = rowSums(across(all_of(genes))))
  track_dir <- file.path(OUT_DIR, gsub(" ", "_", tolower(label)))
  dir.create(track_dir, showWarnings = FALSE, recursive = TRUE)

  # Per-cluster temporal trend
  ct <- df %>% group_by(shell_cluster) %>%
    summarise(n = n(), slope = tryCatch(coef(lm(total ~ year))[[2]], error = function(e) NA_real_),
              p = tryCatch(summary(lm(total ~ year))$coefficients[2,4], error = function(e) 1),
              mean_vf = mean(total), .groups = "drop") %>% arrange(p, desc(abs(slope)))
  cat("Temporal VF trends:\n"); print(ct %>% select(shell_cluster, n, slope, p, mean_vf), n = Inf)
  write_csv(ct, file.path(track_dir, "01_cluster_temporal_trends.csv"))

  target <- as.character(ct$shell_cluster[1])
  slope_val <- ct$slope[1]; p_val <- ct$p[1]
  cat(sprintf("Target: %s (slope=%.4f, p=%.4e)\n", target, slope_val, p_val))

  target_df <- df %>% filter(shell_cluster == target)

  # Clinical enrichment (Fisher)
  enrich <- bind_rows(lapply(genes, function(g) {
    d <- df %>% filter(!is.na(.data[[g]]))
    if (nrow(d) < 20 || length(unique(d[[g]])) < 2) return(tibble(gene = g, or = NA_real_, p = NA_real_))
    ft <- tryCatch(fisher.test(table(d$clinical_binary, d[[g]])), error = function(e) NULL)
    if (is.null(ft)) return(tibble(gene = g, or = NA_real_, p = NA_real_))
    tibble(gene = g, or = unname(ft$estimate), p = ft$p.value)
  })) %>% mutate(p_adj = p.adjust(pmax(p, 0, na.rm = TRUE), method = "BH"))
  write_csv(enrich, file.path(track_dir, "02_clinical_enrichment.csv"))

  # Gene-level temporal trends in target cluster
  trends <- bind_rows(lapply(genes, function(g) {
    d <- target_df %>% filter(!is.na(.data[[g]]))
    if (nrow(d) < 20 || n_distinct(d$year) < 3) return(tibble(gene = g, slope = NA_real_, p = NA_real_))
    f <- tryCatch(lm(as.formula(paste0(g, " ~ year")), data = d), error = function(e) NULL)
    if (is.null(f)) return(tibble(gene = g, slope = NA_real_, p = NA_real_))
    s <- summary(f)$coefficients
    tibble(gene = g, slope = s["year", "Estimate"], p = s["year", "Pr(>|t|)"])
  })) %>% mutate(p_adj = p.adjust(pmax(p, 0, na.rm = TRUE), method = "BH"))
  write_csv(trends, file.path(track_dir, "03_target_temporal_trends.csv"))

  # Early vs late prevalence
  n_early <- sum(target_df$year %in% EARLY_YEARS); n_late <- sum(target_df$year %in% LATE_YEARS)
  if (n_early > 0 && n_late > 0 && length(genes) > 0) {
    pe <- target_df %>% filter(year %in% EARLY_YEARS) %>% select(all_of(genes)) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% pivot_longer(everything(), names_to = "gene", values_to = "early")
    pl <- target_df %>% filter(year %in% LATE_YEARS) %>% select(all_of(genes)) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% pivot_longer(everything(), names_to = "gene", values_to = "late")
    prevalence <- pe %>% inner_join(pl, by = "gene") %>% mutate(delta = late - early, n_early = n_early, n_late = n_late)
    write_csv(prevalence, file.path(track_dir, "04_early_vs_late_prevalence.csv"))
  } else { prevalence <- tibble() }

  # Adjusted model: VF ~ year + niche + country (in target cluster)
  model_adj <- tryCatch({
    m <- lm(total ~ year + source_niche + country, data = target_df)
    tidy(m) %>% filter(term == "year") %>% mutate(model = "VF ~ year + niche + country") %>% select(model, estimate, std.error, p.value)
  }, error = function(e) tibble(model = "failed", estimate = NA_real_, std.error = NA_real_, p.value = NA_real_))
  write_csv(model_adj, file.path(track_dir, "05_adjusted_vf_model.csv"))

  # Interaction model: clinical ~ year * is_target
  master2 <- df %>% mutate(is_target = as.integer(shell_cluster == target))
  interact <- tryCatch({
    m <- glm(clinical_binary ~ year * is_target + country, data = master2, family = binomial())
    tidy(m) %>% filter(stringr::str_detect(term, ":")) %>% mutate(model = "clinical ~ year * target + country")
  }, error = function(e) tibble(model = "failed", term = NA, estimate = NA))
  write_csv(interact, file.path(track_dir, "06_interaction_model.csv"))

  # Summary
  clin_genes <- enrich %>% filter(p_adj < 0.05, or > 1, !is.na(p_adj)) %>% pull(gene)
  inc_genes <- trends %>% filter(p_adj < 0.05, slope > 0, !is.na(p_adj)) %>% pull(gene)
  dec_genes <- trends %>% filter(p_adj < 0.05, slope < 0, !is.na(p_adj)) %>% pull(gene)
  cat(sprintf("Clinically enriched: %d | Increasing in %s: %d | Decreasing: %d\n",
      length(clin_genes), target, length(inc_genes), length(dec_genes)))

  # ---- Figures ----
  # F1: Burden by cluster over time
  p1 <- df %>% group_by(shell_cluster, year) %>%
    summarise(m = mean(total), .groups = "drop") %>%
    mutate(is_target = shell_cluster == target) %>%
    ggplot(aes(x = year, y = m, color = shell_cluster, linewidth = is_target)) +
    geom_line(alpha = 0.7) + scale_linewidth_manual(values = c("TRUE" = 1.5, "FALSE" = 0.5), guide = "none") +
    labs(title = paste(label, "- VF burden by cluster"), x = "Year", y = "Mean VF genes") + theme_minimal()
  ggsave(file.path(track_dir, "F1_burden_by_cluster.png"), p1, width = 8, height = 5)

  # F2: Target burden trend
  ty <- target_df %>% group_by(year) %>%
    summarise(m = mean(total), se = sd(total)/sqrt(n()), .groups = "drop")
  p2 <- ty %>% ggplot(aes(x = year, y = m)) +
    geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.2) +
    geom_line() + geom_point() + geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
    labs(title = paste(target, "-", label, "VF burden"), x = "Year", y = "Mean VF genes") + theme_minimal()
  ggsave(file.path(track_dir, "F2_target_temporal.png"), p2, width = 8, height = 5)

  # F3: Clinical fraction of target over time
  cf <- target_df %>% group_by(year) %>%
    summarise(pct_clinical = mean(clinical_binary) * 100, .groups = "drop")
  p3 <- cf %>% ggplot(aes(x = year, y = pct_clinical)) +
    geom_line() + geom_point() + geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
    labs(title = paste(target, "- clinical fraction over time"), x = "Year", y = "% Clinical") + theme_minimal() + ylim(0, 100)
  ggsave(file.path(track_dir, "F3_clinical_fraction.png"), p3, width = 8, height = 5)

  # F4: Top clinically enriched genes
  tenr <- enrich %>% filter(!is.na(or), !is.na(p_adj), or > 0, is.finite(or)) %>% arrange(p_adj) %>% slice_head(n = 20) %>%
    mutate(gene = fct_reorder(gene, or))
  if (nrow(tenr) >= 3) {
    p4 <- tenr %>% ggplot(aes(x = or, y = gene, color = p_adj < 0.05)) +
      geom_point(size = 3) + geom_vline(xintercept = 1, linetype = "dashed") +
      scale_x_log10() + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
      labs(title = paste("Clinically enriched VF genes -", label), x = "Odds ratio", y = "") +
      theme_minimal() + theme(legend.position = "none")
    ggsave(file.path(track_dir, "F4_clinical_enrichment.png"), p4, width = 7, height = 6)
  }

  # F5: Top increasing VF genes in target
  inc <- trends %>% filter(p_adj < 0.05, slope > 0, !is.na(p_adj)) %>%
    arrange(p_adj) %>% slice_head(n = 15) %>% mutate(gene = fct_reorder(gene, slope))
  if (nrow(inc) >= 3) {
    p5 <- inc %>% ggplot(aes(x = slope, y = gene)) +
      geom_col(fill = "coral") + geom_vline(xintercept = 0, linetype = "dashed") +
      labs(title = paste("Increasing VF genes in", target), x = "Year slope (gene ~ year)", y = "") +
      theme_minimal()
    ggsave(file.path(track_dir, "F5_increasing_genes.png"), p5, width = 7, height = 5)
  }

  # F6: Early vs late delta (top 15)
  if (nrow(prevalence) > 0) {
    top_prev <- prevalence %>% arrange(desc(abs(delta))) %>% slice_head(n = 15) %>%
      mutate(gene = fct_reorder(gene, delta))
    p6 <- top_prev %>% ggplot(aes(x = delta, y = gene)) +
      geom_col(aes(fill = delta > 0)) + geom_vline(xintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "darkred"), guide = "none") +
      labs(title = paste(target, "- early vs late prevalence"),
           subtitle = paste0("Early (2016-2018, n=", n_early, ") vs Late (2022-2025, n=", n_late, ")"),
           x = "Delta prevalence", y = "") + theme_minimal()
    ggsave(file.path(track_dir, "F6_early_vs_late.png"), p6, width = 7, height = 6)
  }

  # Combined figure
  if (exists("p2") && exists("p3")) {
    combined <- p2 / p3 + plot_annotation(title = paste(label, "-", target))
    ggsave(file.path(track_dir, "F_combined_temporal.png"), combined, width = 8, height = 8)
  }

  # Write track summary
  summary <- tibble(Metric = c("Target cluster", "Temporal slope", "Temporal p", "Clinically enriched", "Increasing", "Decreasing"),
                    Value = c(target, sprintf("%.4f", slope_val), sprintf("%.4e", p_val),
                              as.character(length(clin_genes)), as.character(length(inc_genes)), as.character(length(dec_genes))))
  write_csv(summary, file.path(track_dir, "00_summary.csv"))

  list(target = target, enrich = enrich, trends = trends, model_adj = model_adj, interact = interact)
}

# -------------------------------------------------------------------
# 3. VirulenceFinder track
# -------------------------------------------------------------------
ts("Reading VF binary matrix...")
vf_file <- config$VF_BINARY
if (file.exists(vf_file)) {
  raw <- read.delim(vf_file, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "", na.strings = "", stringsAsFactors = FALSE)
  raw <- raw %>% filter(st == st)
  genes <- setdiff(colnames(raw), c("st", "genome"))
  df <- raw %>% mutate(genome_id = as.character(genome)) %>% select(genome_id, all_of(genes))
  df <- df %>% mutate(across(all_of(genes), ~ as.integer(!is.na(.x) & .x != "." & .x != "0" & .x != "")))
  df <- common(df)
  if (nrow(df) > 50) {
    res_vf <- analyze_vf_track("VirulenceFinder", df)
  } else { ts("VF: too few genomes"); res_vf <- NULL }
} else { ts("VF file not found") }

# -------------------------------------------------------------------
# 4. Summary
# -------------------------------------------------------------------
cat("\n========================================\n")
cat("  VIRULENCEFINDER ANALYSIS COMPLETE -", st, "\n")
cat("========================================\n")
if (exists("res_vf") && !is.null(res_vf)) {
  cat("VF: target =", res_vf$target, "\n")
}
cat("Output:", OUT_DIR, "\n")
ts("DONE")
