#!/usr/bin/env Rscript
#
# ResFinder cluster decreasing analysis
#
# Identifies downward-trending shell clusters and analyzes:
#   1. Per-gene decreasing prevalence within that cluster
#   2. Genes that are clinically depleted (less common in clinical)
#   3. Overlap between cluster-decreasing and clinically depleted genes
#
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(ggplot2); library(patchwork); library(scales)
})
source("config.R")

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

# ===== 1. Identify the most significant downward-trend cluster =====
cat("=== Step 1: Cluster temporal trends ===\n")
clust_trends <- read.csv("output/ST10/resistance_analysis/resfinder/01_cluster_temporal_trends.csv",
                          stringsAsFactors = FALSE) %>%
  arrange(p)
print(clust_trends, row.names = FALSE)

# Pick the one with most significant negative slope (lowest p-value among negative-slope clusters)
down_clusters <- clust_trends %>% filter(slope < 0, p < 0.05) %>% arrange(p)
TARGET_CLUSTER <- down_clusters$shell_cluster[1]
cat("\nSelected target cluster:", TARGET_CLUSTER, "\n")
cat("  slope:", down_clusters$slope[1], "p:", down_clusters$p[1], "n:", down_clusters$n[1], "mean_amr:", down_clusters$mean_amr[1], "\n")

# ===== 2. Load data =====
cat("\n=== Step 2: Loading data ===\n")
meta10 <- read_xlsx(config$st_metadata("ST10")) %>%
  rename(genome_id = Name) %>%
  mutate(genome_id = make_genome_id(genome_id),
         source_niche = .data[["Source Niche"]],
         clinical = source_niche == "Human",
         year = as.integer(.data[["Collection Year"]])) %>%
  filter(!is.na(year), !is.na(source_niche))

# ResFinder binary matrix
resf_binary <- read.delim(file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
                           header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
resf10 <- resf_binary %>% filter(st == "ST10") %>%
  mutate(genome_id = make_genome_id(genome)) %>% select(-st, -genome)
resf10_genes <- setdiff(colnames(resf10), "genome_id")

# Shell cluster assignments (use VFDB clusters for ResFinder)
st10_clusters <- read.csv("output/ST10/vfdb_analysis/03_shell_gene_cluster_assignments_k9.csv",
                           stringsAsFactors = FALSE) %>% rename(genome_id = 1, shell_cluster = 2, raw_cluster = 3)

# Join
resf10_long <- resf10 %>% pivot_longer(-genome_id, names_to = "gene", values_to = "present") %>%
  mutate(present = as.integer(!is.na(present) & present != "" & present != "0")) %>%
  left_join(meta10 %>% select(genome_id, clinical, year), by = "genome_id") %>%
  left_join(st10_clusters %>% select(genome_id, shell_cluster), by = "genome_id") %>%
  filter(!is.na(shell_cluster), !is.na(year))

cat("Loaded", n_distinct(resf10_long$genome_id), "genomes,",
    n_distinct(resf10_long$gene), "ARG genes\n")

# ===== 3. Per-gene decreasing analysis in target cluster =====
cat("\n=== Step 3: Per-gene early vs late in", TARGET_CLUSTER, "===\n")

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

# ===== 4. Clinical enrichment for ResFinder genes =====
cat("\n=== Step 4: Clinical enrichment across ALL ST10 ResFinder genes ===\n")

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

# ===== 5. Overlap: cluster-decreasing AND clinically depleted =====
cat("\n=== Step 5: Overlap analysis ===\n")

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

# ===== 6. Write results =====
write.csv(resf_cluster_genes, file.path(OUT, paste0(TARGET_CLUSTER, "_gene_level_early_vs_late.csv")), row.names = FALSE)
write.csv(resf_clinical, file.path(OUT, "ST10_ResFinder_clinical_enrichment.csv"), row.names = FALSE)
write.csv(overlap, file.path(OUT, paste0(TARGET_CLUSTER, "_decreasing_x_clinical_overlap.csv")), row.names = FALSE)

# ===== 7. Summary figure =====
cat("\n=== Step 7: Figures ===\n")

# Cluster trend figure
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

# Top decreasing genes in target cluster
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

# Overlap figure
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

cat("\n=== DONE ===\n")
cat("Results in:", OUT, "\n")
