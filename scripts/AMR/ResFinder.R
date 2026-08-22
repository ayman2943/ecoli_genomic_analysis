#!/usr/bin/env Rscript
#
# ResFinder.R  —  Shell-gene cluster AMR analysis (ResFinder track)
#
# Split out of the original resistance_analysis.R (which ran the CARD and
# ResFinder tracks back-to-back) — see CARD.R for the CARD half. This file
# also folds in the two follow-on ResFinder-only reviewer-response scripts
# that depended on this track's output (originally
# resfinder_clusters.R and resfinder_decreasing.R), each preserved verbatim
# in its own `local({ ... })` block below so their variables can't collide
# with each other or with the main track above — the same isolation
# run_pipeline.R already gives every script via `source(path, local=TRUE)`.
#
# Usage:
#   TARGET_ST=ST10 Rscript scripts/AMR/ResFinder.R
#
suppressPackageStartupMessages({
  library(tidyverse); library(broom); library(writexl); library(readr)
  library(ggplot2); library(scales); library(patchwork); library(forcats)
})
ts <- function(m) cat("[", Sys.time(), "] ", m, "\n", sep = "")

source("config.R")
st <- config$TARGET_ST
OUT_DIR <- file.path(config$OUTPUT_DIR, st, "resistance_analysis")
FIG_DIR <- file.path(OUT_DIR, "figures"); TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
EARLY_YEARS <- c(2016, 2017, 2018); LATE_YEARS <- c(2022, 2023, 2024, 2025)

# -------------------------------------------------------------------
# 1. Load master table
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
# 2. Analysis function (ResFinder track)
# -------------------------------------------------------------------
analyze_amr_track <- function(label, df) {
  ts(paste("\n===========", label, "==========="))
  genes <- setdiff(colnames(df), c("genome_id", "shell_cluster", "year", "clinical_binary", "source_niche", "country"))
  df <- df %>% mutate(total = rowSums(across(all_of(genes))))
  track_dir <- file.path(OUT_DIR, gsub(" ", "_", tolower(label)))
  dir.create(track_dir, showWarnings = FALSE, recursive = TRUE)

  # Per-cluster temporal trend
  ct <- df %>% group_by(shell_cluster) %>%
    summarise(n = n(), slope = tryCatch(coef(lm(total ~ year))[[2]], error = function(e) NA_real_),
              p = tryCatch(summary(lm(total ~ year))$coefficients[2,4], error = function(e) 1),
              mean_amr = mean(total), .groups = "drop") %>% arrange(p, desc(abs(slope)))
  cat("Temporal AMR trends:\n"); print(ct %>% select(shell_cluster, n, slope, p, mean_amr), n = Inf)
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

  # Adjusted model: AMR ~ year + niche + country
  model_adj <- tryCatch({
    m <- lm(total ~ year + source_niche + country, data = target_df)
    tidy(m) %>% filter(term == "year") %>% mutate(model = "AMR ~ year + niche + country") %>% select(model, estimate, std.error, p.value)
  }, error = function(e) tibble(model = "failed", estimate = NA, std.error = NA, p.value = NA))
  write_csv(model_adj, file.path(track_dir, "05_adjusted_amr_model.csv"))

  # Interaction: clinical ~ year * is_target
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
  p1 <- df %>% group_by(shell_cluster, year) %>%
    summarise(m = mean(total), .groups = "drop") %>%
    mutate(is_target = shell_cluster == target) %>%
    ggplot(aes(x = year, y = m, color = shell_cluster, linewidth = is_target)) +
    geom_line(alpha = 0.7) + scale_linewidth_manual(values = c("TRUE" = 1.5, "FALSE" = 0.5), guide = "none") +
    labs(title = paste(label, "- AMR burden by cluster"), x = "Year", y = "Mean AMR genes") + theme_minimal()
  ggsave(file.path(track_dir, "F1_burden_by_cluster.png"), p1, width = 8, height = 5)

  ty <- target_df %>% group_by(year) %>%
    summarise(m = mean(total), se = sd(total)/sqrt(n()), .groups = "drop")
  p2 <- ty %>% ggplot(aes(x = year, y = m)) +
    geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.2) +
    geom_line() + geom_point() + geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
    labs(title = paste(target, "-", label, "AMR burden"), x = "Year", y = "Mean AMR genes") + theme_minimal()
  ggsave(file.path(track_dir, "F2_target_temporal.png"), p2, width = 8, height = 5)

  tenr <- enrich %>% filter(!is.na(or), !is.na(p_adj), or > 0, is.finite(or)) %>% arrange(p_adj) %>% slice_head(n = 20) %>%
    mutate(gene = fct_reorder(gene, or))
  if (nrow(tenr) >= 3) {
    p4 <- tenr %>% ggplot(aes(x = or, y = gene, color = p_adj < 0.05)) +
      geom_point(size = 3) + geom_vline(xintercept = 1, linetype = "dashed") +
      scale_x_log10() + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
      labs(title = paste("Clinically enriched AMR genes -", label), x = "Odds ratio", y = "") +
      theme_minimal() + theme(legend.position = "none")
    ggsave(file.path(track_dir, "F3_clinical_enrichment.png"), p4, width = 7, height = 6)
  }

  inc <- trends %>% filter(p_adj < 0.05, slope > 0, !is.na(p_adj)) %>%
    arrange(p_adj) %>% slice_head(n = 15) %>% mutate(gene = fct_reorder(gene, slope))
  if (nrow(inc) >= 3) {
    p5 <- inc %>% ggplot(aes(x = slope, y = gene)) +
      geom_col(fill = "coral") + geom_vline(xintercept = 0, linetype = "dashed") +
      labs(title = paste("Increasing AMR genes in", target), x = "Year slope", y = "") + theme_minimal()
    ggsave(file.path(track_dir, "F4_increasing_genes.png"), p5, width = 7, height = 5)
  }

  dec <- trends %>% filter(p_adj < 0.05, slope < 0, !is.na(p_adj)) %>%
    arrange(p_adj) %>% slice_head(n = 15) %>% mutate(gene = fct_reorder(gene, slope))
  if (nrow(dec) >= 3) {
    p6 <- dec %>% ggplot(aes(x = slope, y = gene)) +
      geom_col(fill = "steelblue") + geom_vline(xintercept = 0, linetype = "dashed") +
      labs(title = paste("Decreasing AMR genes in", target), x = "Year slope", y = "") + theme_minimal()
    ggsave(file.path(track_dir, "F5_decreasing_genes.png"), p6, width = 7, height = 5)
  }

  if (nrow(prevalence) > 0) {
    top_prev <- prevalence %>% arrange(desc(abs(delta))) %>% slice_head(n = 15) %>%
      mutate(gene = fct_reorder(gene, delta))
    p7 <- top_prev %>% ggplot(aes(x = delta, y = gene)) +
      geom_col(aes(fill = delta > 0)) + geom_vline(xintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "darkred"), guide = "none") +
      labs(title = paste(target, "- early vs late AMR prevalence"),
           subtitle = paste0("Early (2016-2018, n=", n_early, ") vs Late (2022-2025, n=", n_late, ")"),
           x = "Delta prevalence", y = "") + theme_minimal()
    ggsave(file.path(track_dir, "F6_early_vs_late.png"), p7, width = 7, height = 6)
  }

  summary <- tibble(Metric = c("Target cluster", "Temporal slope", "Temporal p", "Clinically enriched", "Increasing", "Decreasing"),
                    Value = c(target, sprintf("%.4f", slope_val), sprintf("%.4e", p_val),
                              as.character(length(clin_genes)), as.character(length(inc_genes)), as.character(length(dec_genes))))
  write_csv(summary, file.path(track_dir, "00_summary.csv"))

  list(target = target, enrich = enrich, trends = trends)
}

# -------------------------------------------------------------------
# 3. ResFinder track
# -------------------------------------------------------------------
ts("Reading ResFinder binary matrix...")
resf_file <- file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv")
if (file.exists(resf_file)) {
  raw <- read.delim(resf_file, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "", na.strings = "", stringsAsFactors = FALSE)
  raw <- raw %>% filter(st == st)
  genes <- setdiff(colnames(raw), c("st", "genome"))
  df <- raw %>% mutate(genome_id = as.character(genome)) %>% select(genome_id, all_of(genes))
  df <- df %>% mutate(across(all_of(genes), ~ as.integer(!is.na(.x) & .x != "." & .x != "0" & .x != "")))
  df <- common(df)
  if (nrow(df) > 50) res_resf <- analyze_amr_track("ResFinder", df) else { ts("ResFinder: too few"); res_resf <- NULL }
} else { ts("ResFinder file not found") }

cat("\n========================================\n")
cat("  RESFINDER RESISTANCE ANALYSIS COMPLETE -", st, "\n")
cat("========================================\n")
if (exists("res_resf") && !is.null(res_resf)) cat("ResFinder: target =", res_resf$target, "\n")
cat("Output:", OUT_DIR, "\n")
ts("DONE")

# =====================================================================
# 4. Follow-on: quick per-cluster decreasing-gene scan (ST10, 3 fixed
#    clusters) — originally resfinder_clusters.R. Preserved verbatim in
#    its own scope; runs after the main track above.
# =====================================================================
local({
  suppressPackageStartupMessages({
    library(tidyverse); library(readxl)
  })

  EARLY <- c(2016, 2017, 2018)
  LATE  <- c(2022, 2023, 2024, 2025)

  make_id <- function(x) {
    x <- as.character(x); x <- trimws(x)
    ifelse(grepl("^Escherichia_coli_", x), x,
           ifelse(grepl("^E\\.coli_", x), sub("^E\\.coli_", "Escherichia_coli_", x),
                  paste0("Escherichia_coli_", x)))
  }

  cat("\n=== [ResFinder follow-on 1/2] Per-cluster decreasing-gene scan ===\n")
  meta <- read_xlsx(config$st_metadata("ST10")) %>%
    rename(id = Name) %>%
    mutate(id = make_id(id),
           clinical = .data[["Source Niche"]] == "Human",
           year = as.integer(.data[["Collection Year"]])) %>%
    filter(!is.na(year), !is.na(clinical))

  rf <- read.delim(file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
                    header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE) %>%
    filter(st == "ST10") %>%
    mutate(id = make_id(genome)) %>%
    select(-st, -genome)

  genes <- setdiff(colnames(rf), "id")

  clust <- read.csv("output/ST10/vfdb_analysis/03_shell_gene_cluster_assignments_k9.csv",
                     stringsAsFactors = FALSE) %>% rename(id = 1, shell_cluster = 2)

  rf_long <- rf %>% pivot_longer(-id, names_to = "gene", values_to = "present") %>%
    mutate(present = as.integer(!is.na(present) & present != "" & present != "0")) %>%
    inner_join(meta %>% select(id, clinical, year), by = "id") %>%
    inner_join(clust %>% select(id, shell_cluster), by = "id")

  analyze_cluster <- function(d, clust_name) {
    d <- d %>% filter(shell_cluster == clust_name) %>%
      mutate(period = case_when(year %in% EARLY ~ "early", year %in% LATE ~ "late")) %>%
      filter(!is.na(period))
    n_early <- n_distinct(d$id[d$period == "early"])
    n_late  <- n_distinct(d$id[d$period == "late"])
    cat("  n_early:", n_early, "n_late:", n_late, "\n")

    per_gene <- d %>% group_by(gene) %>%
      summarise(pe = sum(present[period == "early"]), pl = sum(present[period == "late"]),
                .groups = "drop") %>%
      rowwise() %>%
      mutate(fr = list(tryCatch(fisher.test(matrix(c(pl, n_late-pl, pe, n_early-pe), 2, byrow = TRUE)),
                                 error = function(e) NULL))) %>%
      ungroup()

    per_gene$or <- sapply(per_gene$fr, function(x) if (is.null(x)) NA else unname(x$estimate))
    per_gene$p <- sapply(per_gene$fr, function(x) if (is.null(x)) NA else x$p.value)
    per_gene$fr <- NULL

    per_gene <- per_gene %>% mutate(
      p_adj = p.adjust(p, "BH"),
      prev_early = pe / n_early,
      prev_late  = pl / n_late,
      delta_pp = (prev_late - prev_early) * 100,
      direction = case_when(delta_pp < 0 & p_adj < 0.05 ~ "Decreasing",
                            delta_pp > 0 & p_adj < 0.05 ~ "Increasing",
                            TRUE ~ "No change")
    ) %>% arrange(delta_pp)

    dec <- per_gene %>% filter(direction == "Decreasing")
    inc <- per_gene %>% filter(direction == "Increasing")
    cat("  Decreasing:", nrow(dec), "Increasing:", nrow(inc), "\n")
    if (nrow(dec) > 0) {
      cat("  Top decreasing:\n")
      print(dec %>% select(gene, prev_early, prev_late, delta_pp, or, p_adj), row.names = FALSE)
    }
    if (nrow(inc) > 0) {
      cat("  Top increasing:\n")
      print(inc %>% select(gene, prev_early, prev_late, delta_pp, or, p_adj), row.names = FALSE)
    }
    per_gene
  }

  dir.create("output/resfinder_decreasing_analysis", showWarnings = FALSE, recursive = TRUE)
  for (cl in c("Cluster_1", "Cluster_5", "Cluster_8")) {
    cat("\n===== ", cl, " =====\n")
    n <- n_distinct(rf_long$id[rf_long$shell_cluster == cl])
    cat("  Total genomes:", n, "\n")
    r <- analyze_cluster(rf_long, cl)
    write.csv(r, paste0("output/resfinder_decreasing_analysis/", cl, "_gene_early_vs_late.csv"), row.names = FALSE)
  }
  cat("[ResFinder follow-on 1/2] Done.\n")
})

# =====================================================================
# 5. Follow-on: full decreasing-cluster + clinical-overlap analysis
#    (ST10) — originally resfinder_decreasing.R. Depends on the
#    "01_cluster_temporal_trends.csv" written by the main ResFinder
#    track above (section 3). Preserved verbatim in its own scope,
#    including its original hardcoded ST10 output paths.
# =====================================================================
local({
  suppressPackageStartupMessages({
    library(tidyverse); library(readxl); library(ggplot2); library(patchwork); library(scales)
  })

  OUT <- file.path(config$OUTPUT_DIR, "resfinder_decreasing_analysis")
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

  EARLY_YEARS <- c(2016, 2017, 2018)
  LATE_YEARS  <- c(2022, 2023, 2024, 2025)

  make_genome_id <- function(x) {
    x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
    case_when(grepl("^Escherichia_coli_", x) ~ x,
              grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
              TRUE ~ paste0("Escherichia_coli_", x))
  }

  cat("\n=== [ResFinder follow-on 2/2] Cluster temporal trends ===\n")
  clust_trends <- read.csv("output/ST10/resistance_analysis/resfinder/01_cluster_temporal_trends.csv",
                            stringsAsFactors = FALSE) %>%
    arrange(p)
  print(clust_trends, row.names = FALSE)

  down_clusters <- clust_trends %>% filter(slope < 0, p < 0.05) %>% arrange(p)
  TARGET_CLUSTER <- down_clusters$shell_cluster[1]
  cat("\nSelected target cluster:", TARGET_CLUSTER, "\n")
  cat("  slope:", down_clusters$slope[1], "p:", down_clusters$p[1], "n:", down_clusters$n[1], "mean_amr:", down_clusters$mean_amr[1], "\n")

  cat("\n=== Loading data ===\n")
  meta10 <- read_xlsx(config$st_metadata("ST10")) %>%
    rename(genome_id = Name) %>%
    mutate(genome_id = make_genome_id(genome_id),
           source_niche = .data[["Source Niche"]],
           clinical = source_niche == "Human",
           year = as.integer(.data[["Collection Year"]])) %>%
    filter(!is.na(year), !is.na(source_niche))

  resf_binary <- read.delim(file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
                             header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  resf10 <- resf_binary %>% filter(st == "ST10") %>%
    mutate(genome_id = make_genome_id(genome)) %>% select(-st, -genome)
  resf10_genes <- setdiff(colnames(resf10), "genome_id")

  st10_clusters <- read.csv("output/ST10/vfdb_analysis/03_shell_gene_cluster_assignments_k9.csv",
                             stringsAsFactors = FALSE) %>% rename(genome_id = 1, shell_cluster = 2, raw_cluster = 3)

  resf10_long <- resf10 %>% pivot_longer(-genome_id, names_to = "gene", values_to = "present") %>%
    mutate(present = as.integer(!is.na(present) & present != "" & present != "0")) %>%
    left_join(meta10 %>% select(genome_id, clinical, year), by = "genome_id") %>%
    left_join(st10_clusters %>% select(genome_id, shell_cluster), by = "genome_id") %>%
    filter(!is.na(shell_cluster), !is.na(year))

  cat("Loaded", n_distinct(resf10_long$genome_id), "genomes,",
      n_distinct(resf10_long$gene), "ARG genes\n")

  cat("\n=== Per-gene early vs late in", TARGET_CLUSTER, "===\n")

  gene_cluster_analysis <- function(long_df, target_cluster) {
    clust_df <- long_df %>% filter(shell_cluster == target_cluster) %>%
      mutate(period = case_when(year %in% EARLY_YEARS ~ "early",
                                 year %in% LATE_YEARS ~ "late",
                                 TRUE ~ NA_character_)) %>%
      filter(!is.na(period))

    cat("  Genomes in cluster:", n_distinct(clust_df$genome_id), "\n")
    cat("  Early genomes:", n_distinct(clust_df$genome_id[clust_df$period == "early"]), "\n")
    cat("  Late genomes:", n_distinct(clust_df$genome_id[clust_df$period == "late"]), "\n")

    per_gene <- clust_df %>%
      group_by(gene) %>%
      summarise(
        n_early = sum(period == "early"), n_late = sum(period == "late"),
        present_early = sum(present[period == "early"]),
        present_late  = sum(present[period == "late"]),
        absent_early = n_early - present_early,
        absent_late  = n_late - present_late,
        prev_early = present_early / n_early,
        prev_late  = present_late / n_late,
        delta_pp = (prev_late - prev_early) * 100,
        .groups = "drop"
      ) %>% rowwise() %>%
      mutate(fr = {
        mat <- matrix(c(present_late, absent_late, present_early, absent_early), nrow = 2, byrow = TRUE)
        out <- tryCatch(fisher.test(mat), error = function(e) NULL)
        if (is.null(out)) list(tibble(or = NA_real_, p = NA_real_))
        else list(tibble(or = unname(out$estimate), p = out$p.value))
      }) %>% unnest(fr) %>% ungroup() %>%
      mutate(
        p_adj = p.adjust(p, method = "BH"),
        direction = case_when(
          delta_pp < 0 & p_adj < 0.05 ~ "Decreasing",
          delta_pp > 0 & p_adj < 0.05 ~ "Increasing",
          TRUE ~ "No change"
        )
      ) %>% arrange(delta_pp)

    decreasing_genes <- per_gene %>% filter(direction == "Decreasing")
    cat("  Decreasing genes (BH adj p<0.05):", nrow(decreasing_genes), "\n")
    if (nrow(decreasing_genes) > 0) {
      cat("  Top decreasing:\n")
      print(decreasing_genes %>% select(gene, prev_early, prev_late, delta_pp, or, p_adj), row.names = FALSE)
    }
    per_gene
  }

  resf_cluster_genes <- gene_cluster_analysis(resf10_long, TARGET_CLUSTER)

  cat("\n=== Clinical enrichment across ALL ST10 ResFinder genes ===\n")

  clinical_analysis <- function(long_df) {
    clin_df <- long_df %>% filter(!is.na(clinical))
    per_gene <- clin_df %>%
      group_by(gene) %>%
      summarise(
        n_clin = sum(clinical), n_non = sum(!clinical),
        present_clin = sum(present[clinical]),
        present_non = sum(present[!clinical]),
        absent_clin = n_clin - present_clin,
        absent_non = n_non - present_non,
        prev_clin = present_clin / n_clin,
        prev_non = present_non / n_non,
        delta_clin_pp = (prev_clin - prev_non) * 100,
        .groups = "drop"
      ) %>% filter(present_clin >= 3, present_non >= 3) %>%
      rowwise() %>%
      mutate(fr = {
        mat <- matrix(c(present_clin, absent_clin, present_non, absent_non), nrow = 2, byrow = TRUE)
        out <- tryCatch(fisher.test(mat), error = function(e) NULL)
        if (is.null(out)) list(tibble(or = NA_real_, p = NA_real_))
        else list(tibble(or = unname(out$estimate), p = out$p.value))
      }) %>% unnest(fr) %>% ungroup() %>%
      mutate(
        p_adj = p.adjust(p, method = "BH"),
        clin_direction = case_when(
          or > 1 & p_adj < 0.05 ~ "Clinically enriched",
          or < 1 & p_adj < 0.05 ~ "Clinically depleted",
          TRUE ~ "No association"
        )
      ) %>% arrange(or)

    cat("  Clinically enriched:", sum(per_gene$clin_direction == "Clinically enriched"), "\n")
    cat("  Clinically depleted:", sum(per_gene$clin_direction == "Clinically depleted"), "\n")
    per_gene
  }

  resf_clinical <- clinical_analysis(resf10_long)

  cat("\n=== Overlap analysis ===\n")

  overlap <- resf_cluster_genes %>% filter(direction == "Decreasing") %>%
    inner_join(
      resf_clinical %>% filter(clin_direction != "No association") %>%
        select(gene, clin_direction, clin_or = or, clin_p = p, clin_p_adj = p_adj,
               prev_clin, prev_non, delta_clin_pp),
      by = "gene"
    ) %>% arrange(delta_pp)

  clin_depleted <- overlap %>% filter(clin_direction == "Clinically depleted")
  clin_enriched <- overlap %>% filter(clin_direction == "Clinically enriched")

  cat("Genes decreasing in", TARGET_CLUSTER, "AND clinically depleted:", nrow(clin_depleted), "\n")
  if (nrow(clin_depleted) > 0) {
    print(clin_depleted %>% select(gene, prev_early, prev_late, delta_pp, or, p_adj,
                                   prev_clin, prev_non, delta_clin_pp, clin_or), row.names = FALSE)
  }
  cat("\nGenes decreasing in", TARGET_CLUSTER, "AND clinically enriched:", nrow(clin_enriched), "\n")
  if (nrow(clin_enriched) > 0) {
    print(clin_enriched %>% select(gene, prev_early, prev_late, delta_pp, or, p_adj,
                                    prev_clin, prev_non, delta_clin_pp, clin_or), row.names = FALSE)
  }
  cat("\nGenes decreasing in", TARGET_CLUSTER, "with NO clinical association:",
      sum(resf_cluster_genes$direction == "Decreasing") - nrow(overlap), "\n")

  write.csv(resf_cluster_genes, file.path(OUT, paste0(TARGET_CLUSTER, "_gene_level_early_vs_late.csv")), row.names = FALSE)
  write.csv(resf_clinical, file.path(OUT, "ST10_ResFinder_clinical_enrichment.csv"), row.names = FALSE)
  write.csv(overlap, file.path(OUT, paste0(TARGET_CLUSTER, "_decreasing_x_clinical_overlap.csv")), row.names = FALSE)

  cat("\n=== Figures ===\n")

  p_trend <- clust_trends %>%
    mutate(sig = p < 0.05) %>%
    ggplot(aes(x = reorder(shell_cluster, slope), y = slope, fill = sig)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(aes(label = sprintf("p=%.4f", p)), vjust = ifelse(clust_trends$slope >= 0, -0.5, 1.5),
              size = 2.8, color = "grey40") +
    scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "grey70"), name = "p < 0.05") +
    labs(title = "ResFinder cluster temporal trends", x = NULL, y = "Slope (AMR/year)") +
    theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  top_dec <- resf_cluster_genes %>%
    filter(direction == "Decreasing") %>%
    slice_min(delta_pp, n = min(15, sum(resf_cluster_genes$direction == "Decreasing")))

  if (nrow(top_dec) > 0) {
    p_dec <- top_dec %>%
      ggplot(aes(x = delta_pp, y = reorder(gene, delta_pp), color = p_adj < 0.01)) +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(aes(x = 0, xend = delta_pp), linewidth = 0.7) +
      geom_point(size = 2.5) +
      scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "grey50"), guide = "none") +
      labs(title = paste0("Decreasing genes in ", TARGET_CLUSTER),
           subtitle = "Early (2016-18) vs late (2022-25), Fisher BH adj p<0.05",
           x = "Late minus early prevalence (pp)", y = NULL) +
      theme_classic(base_size = 10)
  } else {
    p_dec <- ggplot() + annotate("text", x=0.5, y=0.5, label="No significantly decreasing genes") +
      labs(title = paste0("Decreasing genes in ", TARGET_CLUSTER)) + theme_void()
  }

  if (nrow(overlap) > 0) {
    p_overlap <- overlap %>%
      ggplot(aes(x = delta_pp, y = reorder(gene, delta_pp), color = clin_direction)) +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(aes(x = 0, xend = delta_pp), linewidth = 0.7) +
      geom_point(size = 3) +
      scale_color_manual(values = c("Clinically enriched" = "#2166AC", "Clinically depleted" = "#D73027")) +
      labs(title = "Overlap: cluster-decreasing & clinical association",
           x = "Cluster delta (pp)", y = NULL, color = "Clinical") +
      theme_classic(base_size = 10)
  } else {
    p_overlap <- ggplot() + annotate("text", x=0.5, y=0.5, label="No overlapping genes") +
      labs(title = "Overlap") + theme_void()
  }

  fig <- (p_trend | p_dec) / p_overlap +
    plot_annotation(title = paste0("ResFinder decreasing analysis in ST10 ", TARGET_CLUSTER),
                    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))

  ggsave(file.path(OUT, "ResFinder_decreasing_analysis.png"), fig, width = 14, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(OUT, "ResFinder_decreasing_analysis.pdf"), fig, width = 14, height = 10, bg = "white")

  cat("\n[ResFinder follow-on 2/2] Done. Results in:", OUT, "\n")
})
