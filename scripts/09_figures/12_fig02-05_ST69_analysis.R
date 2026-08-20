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
# FIGURE 3: PCA of gene content
# ================================================================
cat("\n--- Figure 3: PCA ---\n")

pca_plot <- function(gene_mat, clust_df, title) {
  mat <- as.matrix(gene_mat)
  keep <- complete.cases(mat) & apply(mat, 1, var) > 0
  if (sum(keep) < 10) return(NULL)
  mat <- mat[keep, , drop = FALSE]
  vars <- apply(mat, 2, var, na.rm = TRUE)
  mat <- mat[, vars > 0.01, drop = FALSE]
  if (ncol(mat) < 3) return(NULL)
  pca <- prcomp(mat, scale. = TRUE, center = TRUE)
  ve <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    genome_id = rownames(pca$x)) %>%
    left_join(clust_df, by = "genome_id")
  ggplot(df, aes(x = PC1, y = PC2, color = shell_cluster)) +
    geom_point(alpha = 0.6, size = 0.8) +
    stat_ellipse(aes(fill = shell_cluster), geom = "polygon", alpha = 0.06, level = 0.8, show.legend = FALSE) +
    scale_color_manual(values = cluster_colors(n_distinct(df$shell_cluster))) +
    labs(title = title, x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"), color = "Cluster") +
    theme_classic(base_size = 10)
}

# ST69: use VFDB genes as the matrix (full pangenome not needed)
st69_pca_mat <- st69_master %>% select(genome_id, all_of(st69_vf_genes)) %>%
  distinct(genome_id, .keep_all = TRUE) %>%
  column_to_rownames("genome_id")
st69_clust <- st69_master %>% select(genome_id, shell_cluster) %>% distinct(genome_id, .keep_all = TRUE)
p3A <- pca_plot(st69_pca_mat, st69_clust, "ST69 VFDB")
if (is.null(p3A)) p3A <- ggplot() + annotate("text", x=0.5, y=0.5, label="PCA failed") + labs(title="ST69 VFDB") + theme_void()

# ST10 VF: use binary VF matrix
vf_binary <- read.delim(config$VF_BINARY, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
vf10_mat <- vf_binary %>% filter(st == "ST10") %>%
  mutate(genome_id = make_genome_id(genome)) %>% select(-st, -genome) %>%
  distinct(genome_id, .keep_all = TRUE) %>% column_to_rownames("genome_id")
vf10_mat[] <- lapply(vf10_mat, function(x) as.integer(!is.na(x) & x != "" & x != "0"))
st10_vf_clust <- st10_vf_master %>% select(genome_id, shell_cluster)
p3B <- pca_plot(vf10_mat, st10_vf_clust, "ST10 VirulenceFinder")
if (is.null(p3B)) p3B <- ggplot() + annotate("text", x=0.5, y=0.5, label="PCA failed") + labs(title="ST10 VirulenceFinder") + theme_void()

# ST10 ResFinder: use binary ARG matrix
resf10_mat <- resf10 %>% select(genome_id, all_of(resf10_genes)) %>%
  distinct(genome_id, .keep_all = TRUE) %>% column_to_rownames("genome_id")
resf10_mat[] <- lapply(resf10_mat, function(x) as.integer(!is.na(x) & x != "" & x != "0"))
st10_resf_clust <- st10_resf_master %>% select(genome_id, shell_cluster)
p3C <- pca_plot(resf10_mat, st10_resf_clust, "ST10 ResFinder")
if (is.null(p3C)) p3C <- ggplot() + annotate("text", x=0.5, y=0.5, label="PCA failed") + labs(title="ST10 ResFinder") + theme_void()

fig3 <- (p3A | p3B | p3C) + plot_layout(ncol = 3) +
  plot_annotation(title = "Gene content PCA by shell cluster",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))
ggsave(file.path(OUT, "Fig3_pca_3panel.png"), fig3, width = 18, height = 6, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig3_pca_3panel.pdf"), fig3, width = 18, height = 6, bg = "white")
cat("  Saved Fig3\n")

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

# ---- Combined Fig 2 + Fig 3 (silhouette + content) ----
fig23 <- (fig2 / fig4) + plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Shell-gene cluster silhouette and gene content",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave(file.path(OUT, "Fig2_3_combined.png"), fig23, width = 12, height = 11, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig2_3_combined.pdf"), fig23, width = 12, height = 11, bg = "white")
cat("  Saved Fig2_3_combined\n")

# ================================================================
# FIGURE 5: Oaxaca decomposition
# ================================================================
cat("\n--- Figure 5: Oaxaca decomposition ---\n")

oaxaca_data <- function(df, val_col) {
  df <- df %>% filter(!is.na(year), !is.na(.data[[val_col]]), !is.na(shell_cluster))
  yr_cl <- df %>% group_by(year, shell_cluster) %>%
    summarise(n = n(), mean_val = mean(.data[[val_col]], na.rm = TRUE), .groups = "drop") %>%
    group_by(year) %>%
    mutate(total_n = sum(n), prop = n / total_n) %>% ungroup()
  n_clusters <- n_distinct(yr_cl$shell_cluster)
  yr_valid <- yr_cl %>% group_by(year) %>%
    summarise(all_present = n_distinct(shell_cluster) == n_clusters, .groups = "drop") %>%
    filter(all_present)
  if (nrow(yr_valid) < 2) return(NULL)
  ref_year <- min(yr_valid$year)
  ref <- yr_cl %>% filter(year == ref_year) %>%
    select(shell_cluster, ref_prop = prop, ref_mean = mean_val)
  decomp <- yr_cl %>% left_join(ref, by = "shell_cluster") %>%
    mutate(delta_prop = prop - ref_prop, delta_mean = mean_val - ref_mean,
           composition = delta_prop * ref_mean,
           within = ref_prop * delta_mean,
           interaction = delta_prop * delta_mean)
  total_decomp <- decomp %>% group_by(year) %>%
    summarise(composition = sum(composition, na.rm = TRUE),
              within = sum(within, na.rm = TRUE),
              interaction = sum(interaction, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(-year, names_to = "component", values_to = "value") %>%
    mutate(component = str_to_title(component),
           component = case_when(
             component == "Composition" ~ "Composition",
             component == "Within" ~ "Within-cluster",
             component == "Interaction" ~ "Interaction"),
           component = factor(component, levels = c("Composition", "Within-cluster", "Interaction")))
  obs <- df %>% group_by(year) %>%
    summarise(observed = mean(.data[[val_col]], na.rm = TRUE), .groups = "drop")
  ref_obs <- obs$observed[obs$year == ref_year]
  obs <- obs %>% mutate(delta_from_ref = observed - ref_obs)
  list(components = total_decomp, observed = obs, ref_year = ref_year,
       clusters = decomp)
}

plot_obs <- function(d, title) {
  if (is.null(d)) return(NULL)
  d$observed %>%
    ggplot(aes(x = year, y = observed)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.8, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.6, color = "grey40") +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(title = title, x = NULL, y = "Mean burden", subtitle = paste("Ref:", d$ref_year)) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"))
}
plot_comp <- function(d, title) {
  if (is.null(d)) return(NULL)
  cc <- c("Composition" = "#4DAF4A", "Within-cluster" = "#377EB8", "Interaction" = "#E41A1C")
  d$components %>%
    ggplot(aes(x = year, y = value, color = component, group = component)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5, alpha = 0.85) +
    scale_color_manual(values = cc) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(title = title, x = NULL, y = "Change in weighted mean", color = NULL) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 9), legend.key.size = unit(0.8, "lines"),
          legend.position = "right")
}
plot_drv <- function(d, title, show_smooth = TRUE) {
  if (is.null(d)) return(NULL)
  n_cl <- n_distinct(d$clusters$shell_cluster)
  pal <- if (n_cl <= 9) RColorBrewer::brewer.pal(max(n_cl, 3), "Set1") else rainbow(n_cl)
  p <- d$clusters %>%
    ggplot(aes(x = year, y = within, color = shell_cluster, group = shell_cluster)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.3, alpha = 0.8) +
    scale_color_manual(values = pal) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(title = title, x = "Year", y = "Within-cluster contribution", color = "Cluster") +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 9), legend.key.size = unit(0.8, "lines")) +
    guides(color = guide_legend(ncol = 3))
  if (show_smooth) {
    p <- p + geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5)
  }
  p
}

d69 <- oaxaca_data(st69_master, "total_vf")
d10vf <- oaxaca_data(st10_vf_master, "total_vf")
d10res <- oaxaca_data(st10_resf_master, "total_arg")

lbls <- c("ST69 VFDB", "ST10 VirulenceFinder", "ST10 ResFinder")
ds <- list(d69, d10vf, d10res)
plts <- list()
for (i in seq_along(ds)) {
  plts[[paste0("A", i)]] <- plot_obs(ds[[i]], lbls[i])
  plts[[paste0("B", i)]] <- plot_comp(ds[[i]], lbls[i])
  plts[[paste0("C", i)]] <- plot_drv(ds[[i]], lbls[i])
}
rowA <- wrap_plots(plts[c("A1","A2","A3")], ncol = 3) +
  plot_annotation(title = "A  Observed lineage-level mean burden",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))
rowB <- wrap_plots(plts[c("B1","B2","B3")], ncol = 3) +
  plot_annotation(title = "B  Decomposition components",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))
rowC <- wrap_plots(plts[c("C1","C2","C3")], ncol = 3) & theme(legend.position = "bottom")
rowC <- rowC + plot_annotation(title = "C  Per-cluster within-cluster driver effects",
                                theme = theme(plot.title = element_text(size = 11, face = "bold")))
fig5 <- (rowA / rowB / rowC) + plot_layout(heights = c(1, 1, 1.3)) +
  plot_annotation(
    title = "Oaxaca-style decomposition of burden change",
    subtitle = "Composition = cluster shift x ref.mean  |  Within-cluster = ref.prop x delta mean  |  Interaction = joint",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"))
  )
ggsave(file.path(OUT, "Fig5_decomposition_3panel.png"), fig5, width = 12, height = 12, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig5_decomposition_3panel.pdf"), fig5, width = 12, height = 12, bg = "white")
cat("  Saved Fig5\n")

# Version 2: driver panel without smooth lines
plts_ns <- list()
for (i in seq_along(ds)) {
  plts_ns[[paste0("A", i)]] <- plot_obs(ds[[i]], lbls[i])
  plts_ns[[paste0("B", i)]] <- plot_comp(ds[[i]], lbls[i])
  plts_ns[[paste0("C", i)]] <- plot_drv(ds[[i]], lbls[i], show_smooth = FALSE)
}
rowC_ns <- wrap_plots(plts_ns[c("C1","C2","C3")], ncol = 3) & theme(legend.position = "bottom")
rowC_ns <- rowC_ns + plot_annotation(title = "C  Per-cluster within-cluster driver effects (no trend lines)",
                                     theme = theme(plot.title = element_text(size = 11, face = "bold")))
fig5_ns <- (rowA / rowB / rowC_ns) + plot_layout(heights = c(1, 1, 1.3)) +
  plot_annotation(
    title = "Oaxaca-style decomposition of burden change",
    subtitle = "Composition = cluster shift x ref.mean  |  Within-cluster = ref.prop x delta mean  |  Interaction = joint",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"))
  )
ggsave(file.path(OUT, "Fig5_decomposition_3panel_nosmooth.png"), fig5_ns, width = 12, height = 12, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig5_decomposition_3panel_nosmooth.pdf"), fig5_ns, width = 12, height = 12, bg = "white")
cat("  Saved Fig5 (nosmooth)\n")

# ================================================================
# FIGURE 6: Gene prevalence (copy from earlier script)
# ================================================================
cat("\n--- Figure 6 (copy): Gene prevalence ---\n")
f6 <- file.path(config$OUTPUT_DIR, "figures_gene_prevalence_3panel", "Fig_gene_prevalence_3panel.png")
if (file.exists(f6)) file.copy(f6, file.path(OUT, "Fig6_gene_prevalence_3panel.png"), overwrite = TRUE)
f6p <- sub(".png$", ".pdf", f6)
if (file.exists(f6p)) file.copy(f6p, file.path(OUT, "Fig6_gene_prevalence_3panel.pdf"), overwrite = TRUE)
cat("  Copied Fig6\n")

# ================================================================
# FIGURE 7: Tree mapping
# ================================================================
cat("\n--- Figure 7: Tree mapping ---\n")

tree_plot <- function(tree_file, mapping_file, title, layout = "rectangular") {
  if (!requireNamespace("ape", quietly=TRUE) || !requireNamespace("ggtree", quietly=TRUE)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="ape/ggtree not installed") +
             labs(title = title) + theme_void())
  }
  if (!file.exists(tree_file)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Tree not found") +
             labs(title = title) + theme_void())
  }
  tree <- tryCatch(ape::read.tree(tree_file), error = function(e) NULL)
  if (is.null(tree)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Could not read tree") +
             labs(title = title) + theme_void())
  }
  if (file.exists(mapping_file)) {
    mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)
  } else {
    mapping <- NULL
  }
  tips <- data.frame(genome_id = tree$tip.label, stringsAsFactors = FALSE) %>%
    mutate(genome_id = make_genome_id(genome_id))
  if (!is.null(mapping)) {
    if ("shell_cluster" %in% colnames(mapping)) {
      tips <- tips %>% left_join(mapping %>% select(genome_id, shell_cluster), by = "genome_id")
    } else if (ncol(mapping) >= 2) {
      colnames(mapping)[1:2] <- c("genome_id", "shell_cluster")
      tips <- tips %>% left_join(mapping %>% select(genome_id, shell_cluster), by = "genome_id")
    }
  }
  tips <- tips %>% filter(!is.na(shell_cluster))
  if (nrow(tips) < 10) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Insufficient tip data") +
             labs(title = title) + theme_void())
  }
  # Prune tree to only matched tips
  tree <- tryCatch(ape::keep.tip(tree, tips$genome_id), error = function(e) tree)
  tips <- tips %>% slice(match(tree$tip.label, genome_id))
  # Normalize branch lengths so all panels have comparable proportions
  tree$edge.length <- rep(1, nrow(tree$edge))

  all_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                 "#FFFF33", "#A65628", "#F781BF", "#999999")
  names(all_cols) <- paste0("Cluster_", 1:9)
  present <- sort(unique(tips$shell_cluster))
  clust_cols <- all_cols[present]

  p <- ggtree(tree, layout = layout, size = 0.3, color = "grey50") %<+%
    tips + geom_tippoint(aes(color = shell_cluster), size = 1.2, alpha = 0.9) +
    scale_color_manual(values = clust_cols) +
    labs(title = title, color = "Cluster") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 11, face = "bold"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())
  if (layout == "circular") {
    p <- p + theme(plot.margin = margin(10, 10, 10, 10))
  }
  p
}

# ST69: tree from pangenome, mapping from master table
st69_tree <- "pangenome_output/msa_output/phylogeny/ST69_bootstrap.treefile"
st69_map <- "output/ST69/vfdb_analysis/04_master_shell_cluster_metadata_VFDB_table.csv"
p7A <- tree_plot(st69_tree, st69_map, "ST69")
p7A_circ <- tree_plot(st69_tree, st69_map, "ST69", layout = "circular")

# ST10: shell genome clusters (same pangenome structure for VF and ResFinder)
st10_tree <- "output/ST10/vfdb_analysis/pruned_tree.nwk"
st10_map <- "output/ST10/vfdb_analysis/tree_cluster_mapping.csv"
p7B <- tree_plot(st10_tree, st10_map, "ST10")
p7B_circ <- tree_plot(st10_tree, st10_map, "ST10", layout = "circular")

fig7 <- (p7A | p7B) + plot_layout(ncol = 2) +
  plot_annotation(theme = theme(plot.title = element_blank()))
ggsave(file.path(OUT, "Fig7_tree_mapping_2panel.png"), fig7, width = 16, height = 8, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig7_tree_mapping_2panel.pdf"), fig7, width = 16, height = 8, bg = "white")
cat("  Saved Fig7 (linear)\n")

fig7_circ <- (p7A_circ | p7B_circ) + plot_layout(ncol = 2) +
  plot_annotation(theme = theme(plot.title = element_blank()))
ggsave(file.path(OUT, "Fig7_tree_mapping_2panel_circular.png"), fig7_circ, width = 14, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig7_tree_mapping_2panel_circular.pdf"), fig7_circ, width = 14, height = 7, bg = "white")
cat("  Saved Fig7 (circular)\n")

# ================================================================
# FIGURE 8: Clinical trajectory
# ================================================================
cat("\n--- Figure 8: Clinical trajectory ---\n")

clin_traj <- function(d, title) {
  if (!"year" %in% colnames(d) || !"shell_cluster" %in% colnames(d)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Missing columns") +
             labs(title = title) + theme_void())
  }
  if ("clinical_binary" %in% colnames(d)) {
    d <- d %>% mutate(clin_flag = as.integer(clinical_binary == 1 | clinical_binary == "Human" | clinical_binary == TRUE))
  } else if ("clinical" %in% colnames(d)) {
    d <- d %>% mutate(clin_flag = as.integer(clinical))
  } else {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No clinical data") +
             labs(title = title) + theme_void())
  }
  d <- d %>% filter(!is.na(year)) %>%
    group_by(year, shell_cluster) %>%
    summarise(pct_clin = mean(clin_flag, na.rm = TRUE) * 100, n = n(), .groups = "drop")
  # Keep top clusters by total n
  top <- d %>% group_by(shell_cluster) %>% summarise(total = sum(n), .groups = "drop") %>%
    slice_max(total, n = 5)
  d2 <- d %>% filter(shell_cluster %in% top$shell_cluster)
  if (nrow(d2) < 5) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Insufficient data") +
             labs(title = title) + theme_void())
  }
  ggplot(d2, aes(x = year, y = pct_clin, color = shell_cluster, group = shell_cluster)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    scale_color_manual(values = cluster_colors(n_distinct(d2$shell_cluster))) +
    scale_x_continuous(breaks = pretty(d2$year, n = 5)) +
    ylim(0, 100) +
    labs(title = title, x = "Year", y = "Clinical isolates (%)", color = "Cluster") +
    theme_classic(base_size = 10)
}

p8A <- clin_traj(st69_master, "ST69 VFDB")
p8B <- clin_traj(st10_vf_master, "ST10 VirulenceFinder")
p8C <- clin_traj(st10_resf_master, "ST10 ResFinder")

fig8 <- (p8A | p8B | p8C) + plot_layout(ncol = 3) +
  plot_annotation(title = "Clinical trajectory by shell cluster",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))
ggsave(file.path(OUT, "Fig8_clinical_trajectory_3panel.png"), fig8, width = 16, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig8_clinical_trajectory_3panel.pdf"), fig8, width = 16, height = 5, bg = "white")
cat("  Saved Fig8\n")

# ================================================================
# FIGURE 9: Clinical enrichment (copy from earlier script)
# ================================================================
cat("\n--- Figure 9 (copy): Clinical enrichment ---\n")
f9 <- file.path(config$OUTPUT_DIR, "figures_clinical_enrichment_3panel", "Fig_clinical_enrichment_3panel.png")
if (file.exists(f9)) file.copy(f9, file.path(OUT, "Fig9_clinical_enrichment_3panel.png"), overwrite = TRUE)
f9p <- sub(".png$", ".pdf", f9)
if (file.exists(f9p)) file.copy(f9p, file.path(OUT, "Fig9_clinical_enrichment_3panel.pdf"), overwrite = TRUE)
cat("  Copied Fig9\n")

cat("\n=== ALL FIGURES 2-9 COMPLETE ===\n")
cat("Output directory:", OUT, "\n")
