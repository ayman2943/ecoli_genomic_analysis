#!/usr/bin/env Rscript
# ==============================================================================
# Analysis / Figure4.R
# ==============================================================================
# Reproduces submitted manuscript Figure 4: "Oaxaca-style decomposition of
# burden change" (ST69 VFDB / ST10 VirulenceFinder / ST10 ResFinder), with
# subtitle "Composition = cluster shift x ref.mean | Within-cluster =
# ref.prop x delta mean | Interaction = joint".
#
# PROVENANCE: extracted from `12_fig02-05_ST69_analysis.R`'s own
# "FIGURE 5: Oaxaca decomposition" section -- confirmed via an EXACT text
# match of the plot title/subtitle against the submitted Figure_04.png.
#
# The original section produces TWO renderings from the same underlying
# decomposition data: "Fig5_decomposition_3panel.png" (with linear trend
# lines on the driver panel) and "Fig5_decomposition_3panel_nosmooth.png"
# (without). The submitted figure matches the "nosmooth" rendering, but
# BOTH are kept here, verbatim, since they share the same computation and
# dropping either would be a change to the original analysis rather than a
# pure reorganisation.
#
# Everything else in the original `12_fig02-05_ST69_analysis.R` (Figures 2/3
# and the trailing sections) is out of scope for this file -- see
# Analysis/Figure2.R and the full-script backup block in
# Analysis/Supplementary.R.
#
# ANALYSIS IS UNCHANGED -- this is a verbatim extraction of the relevant
# section, not a rewrite of the underlying logic.
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

cat("\n=== Figure 4 complete ===\n")
cat("Output directory:", OUT, "\n")
