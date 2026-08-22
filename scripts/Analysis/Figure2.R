#!/usr/bin/env Rscript
# ==============================================================================
# Analysis / Figure2.R
# ==============================================================================
# Reproduces submitted manuscript Figure 2: "Shell-gene cluster silhouette
# and gene content."
#
# PROVENANCE: extracted from `12_fig02-05_ST69_analysis.R`. That single
# original script internally computed nine draft "Fig2"-"Fig9" panels; only
# two of them are the source of the submitted figures, confirmed by direct
# visual/text match against the submitted Figure_02.png:
#   - its own "FIGURE 2: Silhouette" section         -> top rows of Figure 2
#   - its own "FIGURE 4: Gene content by cluster" section -> bottom row
#   - its own "fig23 <- (fig2 / fig4)" combine step, which produces
#     "Fig2_3_combined.png" with the exact title "Shell-gene cluster
#     silhouette and gene content" -- an EXACT text match to the submitted
#     figure -- IS Figure 2.
#
# Everything else that lived in the original `12_fig02-05_ST69_analysis.R`
# (its own internal "Figure 3" PCA panel, "Figure 6"/"Figure 9" file-copy
# stubs, "Figure 7" preliminary tree mapping, "Figure 8" clinical
# trajectory) is NOT part of submitted Figure 2 (or any other main figure)
# and is preserved separately, in full, as an extra block inside
# `Analysis/Supplementary.R` -- so nothing from the original script is lost,
# it is just not duplicated into this figure-specific extraction.
#
# ANALYSIS IS UNCHANGED -- this is a verbatim extraction of the relevant
# sections, not a rewrite of the underlying logic.
# ==============================================================================

#!/usr/bin/env Rscript
#
# All publication figures (2-9) as 3-panel comparison
# Panel A: ST69 VFDB | Panel B: ST10 VirulenceFinder | Panel C: ST10 ResFinder
#
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(ggplot2); library(patchwork)
  library(scales)
  if (requireNamespace("ape", quietly=TRUE)) library(ape)
  if (requireNamespace("ggtree", quietly=TRUE)) library(ggtree)
})
source("config.R")

OUT <- file.path(config$OUTPUT_DIR, "figures_all_3panel")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_classic(base_size = 10))

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

cluster_colors <- function(n) {
  if (n <= 9) RColorBrewer::brewer.pal(max(n, 3), "Set1") else rainbow(n)
}


cat("Loading metadata...\n")
meta69 <- read_xlsx(config$st_metadata("ST69")) %>%
  rename(genome_id = Name) %>%
  mutate(genome_id = make_genome_id(genome_id),
         source_niche = .data[["Source Niche"]],
         clinical = source_niche == "Human",
         year = as.integer(.data[["Collection Year"]])) %>%
  filter(!is.na(year), !is.na(source_niche))
meta10 <- read_xlsx(config$st_metadata("ST10")) %>%
  rename(genome_id = Name) %>%
  mutate(genome_id = make_genome_id(genome_id),
         source_niche = .data[["Source Niche"]],
         clinical = source_niche == "Human",
         year = as.integer(.data[["Collection Year"]])) %>%
  filter(!is.na(year), !is.na(source_niche))

# ===== Load data per case =====
# ===== Load data per case =====
# Case A: ST69 VFDB
cat("Loading ST69 VFDB master...\n")
st69_master <- read.csv("output/ST69/vfdb_analysis/04_master_shell_cluster_metadata_VFDB_table.csv",
                         stringsAsFactors = FALSE) %>%
  mutate(genome_id = make_genome_id(genome_id))
# Master already has clinical_binary, year, source_niche columns
st69_vf_genes <- setdiff(colnames(st69_master),
                          c("genome_id","shell_cluster","raw_cluster","X.FILE","total_vf",
                            "clinical_binary","year","source_niche","continent","country","niche","genome"))
# Convert VFDB gene columns to binary presence
for (g in st69_vf_genes) {
  st69_master[[g]] <- as.integer(!is.na(st69_master[[g]]) & st69_master[[g]] != "." & st69_master[[g]] != "")
}
st69_master <- st69_master %>% distinct(genome_id, .keep_all = TRUE)
st69_master$year <- as.integer(st69_master$year)
st69_sil <- read.csv("output/ST69/vfdb_analysis/02_silhouette_k4_to_k10.csv", stringsAsFactors = FALSE)
st69_ksel <- read.csv("output/ST69/vfdb_analysis/02B_k_selection_summary.csv", stringsAsFactors = FALSE)

# Case B: ST10 VF
cat("Loading ST10 VF master...\n")
st10_vf_master <- read.csv("output/ST10/virulencefinder_validation/04_master_shell_cluster_metadata_VF_table.csv",
                            stringsAsFactors = FALSE) %>%
  mutate(genome_id = make_genome_id(genome_id))
st10_vf_master$year <- as.integer(st10_vf_master$year)
# clinical_binary already exists; convert to integer
st10_vf_master <- st10_vf_master %>%
  mutate(clinical_binary = as.integer(clinical_binary == 1 | clinical_binary == "Human" | clinical_binary == TRUE))
st10_vf_sil <- read.csv("output/ST10/virulencefinder_validation/02_silhouette_k4_to_k10.csv", stringsAsFactors = FALSE)
st10_vf_ksel <- read.csv("output/ST10/virulencefinder_validation/02B_k_selection_summary.csv", stringsAsFactors = FALSE)

# Case C: ST10 ResFinder
cat("Loading ST10 ResFinder data...\n")
st10_resf_master <- read.csv("output/ST10/vfdb_analysis/04_master_shell_cluster_metadata_VFDB_table.csv",
                              stringsAsFactors = FALSE) %>%
  mutate(genome_id = make_genome_id(genome_id)) %>%
  left_join(meta10 %>% select(genome_id, clinical, year), by = "genome_id") %>%
  mutate(clinical_binary = as.integer(clinical))
# Actually this is VFDB data, not ResFinder. For ResFinder we need the resfinder binary matrix
# and then join with VFDB shell cluster assignments (same pangenome)
st10_shell_clusters <- read.csv("output/ST10/vfdb_analysis/03_shell_gene_cluster_assignments_k9.csv",
                                 stringsAsFactors = FALSE) %>% rename(genome_id = 1, shell_cluster = 2, raw_cluster = 3)
st10_resf_sil <- read.csv("output/ST10/resistance_analysis/resfinder/02_silhouette_k4_to_k10.csv", stringsAsFactors = FALSE)
st10_resf_ksel <- read.csv("output/ST10/resistance_analysis/resfinder/02B_k_selection_summary.csv", stringsAsFactors = FALSE)

resf_binary <- read.delim(file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
                           header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
resf10 <- resf_binary %>% filter(st == "ST10") %>%
  mutate(genome_id = make_genome_id(genome)) %>%
  select(-st, -genome)
resf10_genes <- setdiff(colnames(resf10), "genome_id")
resf10_long <- resf10 %>% pivot_longer(-genome_id, names_to = "gene", values_to = "present") %>%
  mutate(present = as.integer(!is.na(present) & present != "" & present != "0"))
# Compute total ARG per genome
resf10_total <- resf10_long %>% group_by(genome_id) %>%
  summarise(total_arg = sum(present), .groups = "drop")

# Join with metadata and cluster assignments
st10_resf_master <- resf10_total %>% left_join(meta10 %>% select(genome_id, clinical, year), by = "genome_id") %>%
  left_join(st10_shell_clusters %>% select(genome_id, shell_cluster), by = "genome_id") %>%
  filter(!is.na(shell_cluster))

cat("Data loaded.\n")
cat("  ST69:", nrow(st69_master), "genomes,", length(st69_vf_genes), "VFDB genes\n")
cat("  ST10 VF:", nrow(st10_vf_master), "genomes\n")
cat("  ST10 ResF:", nrow(st10_resf_master), "genomes\n")


# ================================================================
# FIGURE 2: Silhouette
# ================================================================
cat("\n--- Figure 2: Silhouette ---\n")

plot_silhouette <- function(sil, ksel, title) {
  # Find best k from ksel (different column names across outputs)
  ksel_cols <- colnames(ksel)
  best_k <- NA
  if ("best_k_primary" %in% ksel_cols) best_k <- ksel$best_k_primary[1]
  else if ("best_k" %in% ksel_cols) best_k <- ksel$best_k[1]
  else if ("selected_k" %in% ksel_cols) best_k <- ksel$selected_k[1]
  if (is.na(best_k)) return(ggplot() + annotate("text", x=0.5, y=0.5, label="Could not determine best k") +
                              labs(title=title) + theme_void())

  # Silhouette data by k
  sil_cols <- colnames(sil)
  k_col <- if ("k" %in% sil_cols) "k" else "K"
  sil_col <- if ("average_silhouette_width" %in% sil_cols) "average_silhouette_width" else
    if ("avg_sil" %in% sil_cols) "avg_sil" else sil_cols[3]
  method_col <- if ("method" %in% sil_cols) "method" else NULL

  # Filter to ward.D2 method if available; else first method
  sil <- sil %>% rename(kk = !!k_col, sil_val = !!sil_col)
  if (!is.null(method_col)) {
    sil <- sil %>% rename(m = !!method_col)
    if ("ward.D2" %in% sil$m) sil <- sil %>% filter(m == "ward.D2")
  }

  if (n_distinct(sil$kk) > 1) {
    p1 <- sil %>% mutate(best = kk == best_k) %>%
      ggplot(aes(x = kk, y = sil_val)) +
      geom_line(color = "grey60") + geom_point(aes(color = best), size = 2.5) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
      labs(x = "k", y = "Mean silhouette") + ylim(0, NA)
  } else {
    # Only one k value available
    p1 <- sil %>% mutate(best = TRUE) %>%
      ggplot(aes(x = factor(kk), y = sil_val)) +
      geom_col(aes(fill = best), show.legend = FALSE) +
      scale_fill_manual(values = c("TRUE" = "red")) +
      labs(x = "Selected k", y = "Mean silhouette") + ylim(0, NA)
  }

  # Per-k silhouette not available (only average per k), so show bar of best_k vs others
  p2 <- sil %>% mutate(label = ifelse(kk == best_k, paste0("k=", kk, " (selected)"), paste0("k=", kk))) %>%
    ggplot(aes(x = sil_val, y = reorder(label, sil_val), fill = kk == best_k)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey60"), guide = "none") +
    labs(x = "Mean silhouette", y = NULL) + xlim(0, NA)

  wrap_elements(p1 | p2) + labs(title = title) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

p2A <- plot_silhouette(st69_sil, st69_ksel, "ST69 VFDB")
p2C <- plot_silhouette(st10_vf_sil, st10_vf_ksel, "ST10 shell clusters")
fig2 <- (p2A | p2C) + plot_layout(ncol = 2) +
  plot_annotation(title = "Shell-gene cluster silhouette analysis",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))
ggsave(file.path(OUT, "Fig2_silhouette_3panel.png"), fig2, width = 12, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig2_silhouette_3panel.pdf"), fig2, width = 12, height = 5, bg = "white")
cat("  Saved Fig2\n")

# ================================================================
# FIGURE 4: Gene content by cluster
# ================================================================
cat("\n--- Figure 4: Content by cluster ---\n")

content_plot <- function(d, val_col, ylab, title) {
  summ <- d %>% group_by(shell_cluster) %>%
    summarise(n = n(), mean_val = mean(.data[[val_col]], na.rm = TRUE),
              total = sum(.data[[val_col]], na.rm = TRUE), .groups = "drop") %>% drop_na()
  ggplot(summ, aes(x = reorder(shell_cluster, desc(n)), y = total, fill = shell_cluster)) +
    geom_col() +
    geom_text(aes(label = n, y = 0), vjust = 1.5, size = 3, color = "white") +
    scale_fill_manual(values = cluster_colors(nrow(summ)), guide = "none") +
    labs(title = title, x = NULL, y = ylab) +
    theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p4A <- content_plot(st69_master, "total_vf", "Total VFDB genes", "ST69 VFDB")
p4B <- content_plot(st10_vf_master, "total_vf", "Total VF genes", "ST10 VirulenceFinder")
p4C <- content_plot(st10_resf_master, "total_arg", "Total ARG count", "ST10 ResFinder")

fig4 <- (p4A | p4B | p4C) + plot_layout(ncol = 3) +
  plot_annotation(title = "Gene content by shell cluster",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))
ggsave(file.path(OUT, "Fig4_content_by_cluster_3panel.png"), fig4, width = 14, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig4_content_by_cluster_3panel.pdf"), fig4, width = 14, height = 5, bg = "white")
cat("  Saved Fig4\n")

# ---- Combined Fig 2 + Fig 4 (silhouette + content) -- THIS IS SUBMITTED FIGURE 2 ----
# ---- Combined Fig 2 + Fig 3 (silhouette + content) ----
fig23 <- (fig2 / fig4) + plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Shell-gene cluster silhouette and gene content",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave(file.path(OUT, "Fig2_3_combined.png"), fig23, width = 12, height = 11, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig2_3_combined.pdf"), fig23, width = 12, height = 11, bg = "white")
cat("  Saved Fig2_3_combined\n")

cat("\n=== Figure 2 complete ===\n")
cat("Output directory:", OUT, "\n")
