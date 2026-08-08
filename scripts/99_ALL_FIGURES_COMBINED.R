#!/usr/bin/env Rscript
#
# 99_ALL_FIGURES_COMBINED.R
#
# Master script: generates ALL publication figures (1-9) + Supplementary S1-S2
# into a single multi-page PDF.
#
# Output: All_Figures.pdf (in config$OUTPUT_DIR)
#

#

suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(ggplot2); library(patchwork)
  library(scales); library(broom); library(writexl)
  if (requireNamespace("ape", quietly=TRUE)) library(ape)
  if (requireNamespace("phangorn", quietly=TRUE)) library(phangorn)
})

# Config is loaded from the repo-root config.R (env-driven paths).
if (file.exists("config.R")) {
  source("config.R")
} else {
  source("../config.R")
}
OUT <- file.path(config$OUTPUT_DIR, "combined_figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

cluster_colors <- function(n) {
  if (n <= 9) RColorBrewer::brewer.pal(max(n, 3), "Set1") else rainbow(n)
}

theme_pub <- theme_classic(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom", plot.title = element_text(face = "bold", size = 11))
theme_set(theme_pub)

pdf_file <- file.path(OUT, "All_Figures.pdf")
pdf(pdf_file, width = 12, height = 8, bg = "white")

cat("============================================\n")
cat(" Generating ALL figures into:\n")
cat(" ", pdf_file, "\n")
cat("============================================\n\n")

# ============================================================
# FIGURE 1: Temporal trends across 5 STs
# ============================================================
cat("--- Figure 1: Temporal trends + Mann-Kendall ---\n")

mk_test <- function(x) {
  n <- length(x); if (n < 3) return(list(tau = NA, p = 1))
  s <- 0; for (i in 1:(n-1)) for (j in (i+1):n) s <- s + sign(x[j] - x[i])
  denom <- n * (n - 1) / 2; tau <- s / denom
  var_s <- n * (n - 1) * (2 * n + 5) / 18
  z <- if (var_s > 0) (s - sign(s)) / sqrt(var_s) else 0
  list(tau = tau, p = 2 * pnorm(-abs(z)))
}

st_list <- c("ST10", "ST69", "ST73", "ST95", "ST131")
all_data <- list()

for (st in st_list) {
  cat("  ", st, "...\n", sep = "")
  meta_file <- config$st_metadata(st)
  if (!file.exists(meta_file)) next
  meta <- read_xlsx(meta_file)
  name_col <- grep("^(Name|genome|strain|isolate|assembly|sample)$", colnames(meta), value = TRUE)[1]
  if (is.na(name_col)) name_col <- colnames(meta)[2]
  year_col <- grep("year", colnames(meta), ignore.case = TRUE, value = TRUE)[1]
  meta <- meta %>% rename(genome_id = all_of(name_col)) %>%
    mutate(genome_id = make_genome_id(genome_id), year = as.integer(.data[[year_col]])) %>%
    filter(!is.na(year)) %>% select(genome_id, year)

  vfdb_file <- config$st_vfdb_summary(st)
  if (file.exists(vfdb_file)) {
    v <- read_tsv(vfdb_file, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
      mutate(genome_id = make_genome_id(str_remove(basename(.data[["#FILE"]]), "_vfdb\\.tsv$")),
             burden = as.numeric(NUM_FOUND)) %>%
      select(genome_id, burden) %>% inner_join(meta, by = "genome_id") %>%
      group_by(year) %>% summarise(m = mean(burden, na.rm = TRUE),
        se = sd(burden, na.rm = TRUE) / sqrt(n()), n = n(), .groups = "drop") %>%
      filter(n >= 5) %>% mutate(st = st, db = "VFDB")
    all_data[[length(all_data) + 1]] <- v
  }

  card_file <- config$st_card_burden(st)
  if (file.exists(card_file)) {
    c <- read_tsv(card_file, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
      mutate(genome_id = make_genome_id(str_remove(basename(.data[["#FILE"]]), "_card\\.tsv$")),
             burden = as.numeric(NUM_FOUND)) %>%
      select(genome_id, burden) %>% inner_join(meta, by = "genome_id") %>%
      group_by(year) %>% summarise(m = mean(burden, na.rm = TRUE),
        se = sd(burden, na.rm = TRUE) / sqrt(n()), n = n(), .groups = "drop") %>%
      filter(n >= 5) %>% mutate(st = st, db = "CARD")
    all_data[[length(all_data) + 1]] <- c
  }

  vf_file <- config$st_vf_summary(st)
  if (file.exists(vf_file)) {
    vf <- read_tsv(vf_file, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
      mutate(genome_id = make_genome_id(str_remove(basename(.data[["#FILE"]]), "_vf\\.tsv$")),
             burden = as.numeric(NUM_FOUND)) %>%
      select(genome_id, burden) %>% inner_join(meta, by = "genome_id") %>%
      group_by(year) %>% summarise(m = mean(burden, na.rm = TRUE),
        se = sd(burden, na.rm = TRUE) / sqrt(n()), n = n(), .groups = "drop") %>%
      filter(n >= 5) %>% mutate(st = st, db = "VF")
    all_data[[length(all_data) + 1]] <- vf
  }

  resf_file <- config$st_resfinder_summary(st)
  if (file.exists(resf_file)) {
    r <- read_tsv(resf_file, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
      mutate(genome_id = make_genome_id(str_remove(basename(.data[["#FILE"]]), "_resfinder\\.tsv$")),
             burden = as.numeric(NUM_FOUND)) %>%
      select(genome_id, burden) %>% inner_join(meta, by = "genome_id") %>%
      group_by(year) %>% summarise(m = mean(burden, na.rm = TRUE),
        se = sd(burden, na.rm = TRUE) / sqrt(n()), n = n(), .groups = "drop") %>%
      filter(n >= 5) %>% mutate(st = st, db = "ResFinder")
    all_data[[length(all_data) + 1]] <- r
  }
}

df <- bind_rows(all_data)
st_colors <- c(ST10 = "#FF6B6B", ST131 = "#4ECDC4", ST69 = "#95E1D3", ST73 = "#F38181", ST95 = "#FFA07A")
st_levels <- c("ST10", "ST69", "ST73", "ST95", "ST131")
db_levels <- c("CARD", "ResFinder", "VFDB", "VF")
db_names <- c(CARD = "CARD (ARG)", ResFinder = "ResFinder (ARG)", VFDB = "VFDB (VF)", VF = "VirulenceFinder (VF)")

df <- df %>% mutate(st = factor(st, levels = st_levels), db = factor(db, levels = db_levels))

fig1_top <- ggplot(df, aes(year, m, color = st, fill = st)) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.1, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
  facet_wrap(~ db, nrow = 1, ncol = 4, labeller = labeller(db = db_names), scales = "free_y") +
  scale_color_manual(values = st_colors, name = "Sequence Type") +
  scale_fill_manual(values = st_colors, guide = "none") +
  labs(x = "Year", y = "Mean gene burden", title = "Temporal trends across ExPEC lineages (2016-2025)") +
  theme(legend.position = "bottom")

mk_results <- tibble()
for (combo in unique(paste(df$st, df$db))) {
  parts <- str_split(combo, " ")[[1]]
  sub_df <- df %>% filter(st == parts[1], db == parts[2]) %>% arrange(year)
  if (nrow(sub_df) >= 3) {
    mk <- mk_test(sub_df$m)
    mk_results <- bind_rows(mk_results, tibble(st = parts[1], db = parts[2],
      tau = mk$tau, p = mk$p))
  }
}
mk_results <- mk_results %>% mutate(p_adj = p.adjust(p, method = "BH"),
  sig = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", p_adj < 0.05 ~ "*", TRUE ~ "ns"))

fig1_bot <- ggplot(mk_results, aes(x = st, y = db, fill = tau)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0("tau=", round(tau, 3), "\n", sig)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "MK tau") +
  labs(x = "Sequence Type", y = "Database", title = "Mann-Kendall trend summary (BH-adjusted)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig1 <- fig1_top / fig1_bot + plot_layout(heights = c(3, 1))
print(fig1)

# ============================================================
# FIGURES 2-5: ST69 3-panel analysis
# ============================================================
cat("\n--- Figures 2-5: ST69 shell-gene clustering + decomposition ---\n")

st69_master <- read.csv(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv"), stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(genome_id = make_genome_id(as.character(genome_id)),
         year = as.integer(year), total = total_vf)
st69_sil <- read.csv(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "02_silhouette_k4_to_k10.csv"), stringsAsFactors = FALSE)
st69_ksel <- read.csv(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "02B_k_selection_summary.csv"), stringsAsFactors = FALSE)

# Figure 2: Silhouette
fig2 <- ggplot(st69_sil, aes(x = factor(k), y = avg_silhouette_width)) +
  geom_col(fill = "#4ECDC4", alpha = 0.8, width = 0.6) +
  geom_point(size = 2) + geom_line(aes(group = 1), linewidth = 0.5) +
  labs(x = "Number of clusters (k)", y = "Average silhouette width",
       title = "Figure 2. ST69 silhouette analysis: k = 4 is optimal") +
  theme_pub
print(fig2)

# Figure 3: VFDB burden by cluster
cluster_summary <- st69_master %>% group_by(shell_cluster) %>%
  summarise(mean_burden = mean(total, na.rm = TRUE),
    clinical_pct = mean(clinical_binary, na.rm = TRUE) * 100, n = n(), .groups = "drop")
fig3 <- ggplot(st69_master, aes(x = reorder(shell_cluster, total, FUN = median), y = total,
  fill = factor(shell_cluster))) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = cluster_colors(length(unique(st69_master$shell_cluster)))) +
  labs(x = "Shell-gene cluster", y = "VFDB burden (total genes)",
       title = "Figure 3. ST69 VFDB burden by shell-gene cluster") +
  theme(legend.position = "none")
print(fig3)

# Figure 4: Virulence system profiles
vf_genes_key <- c("kpsM", "kpsC", "kpsE", "kpsF", "papA", "papB", "papC", "papF", "papG", "papH",
                   "papI", "papJ", "papK", "papX", "sat", "hlyD", "hlyA", "cnf1", "iutA", "iha")
gene_prev <- st69_master %>%
  select(genome_id, shell_cluster, all_of(intersect(vf_genes_key, colnames(st69_master)))) %>%
  pivot_longer(-c(genome_id, shell_cluster), names_to = "gene", values_to = "present") %>%
  mutate(present = as.integer(!is.na(present) & present != "" & present != "0")) %>%
  group_by(shell_cluster, gene) %>% summarise(prevalence = mean(present) * 100, .groups = "drop")
fig4 <- ggplot(gene_prev, aes(x = gene, y = prevalence, fill = factor(shell_cluster))) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = cluster_colors(length(unique(st69_master$shell_cluster)))) +
  labs(x = "Virulence gene", y = "Prevalence (%)", fill = "Cluster",
       title = "Figure 4. Virulence-system profiles across clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
print(fig4)

# Figure 5: Decomposition
decomp_data <- st69_master %>% filter(!is.na(shell_cluster)) %>%
  group_by(year, shell_cluster) %>% summarise(m = mean(total, na.rm = TRUE), .groups = "drop") %>%
  complete(year, shell_cluster, fill = list(m = NA))
decomp_wide <- decomp_data %>% pivot_wider(names_from = shell_cluster, values_from = m)
ref_year <- min(decomp_wide$year, na.rm = TRUE)
ref_vals <- decomp_wide %>% filter(year == ref_year) %>% select(-year) %>% as.numeric()
names(ref_vals) <- colnames(decomp_wide)[-1]
overall_mean <- st69_master %>% group_by(year) %>% summarise(overall = mean(total, na.rm = TRUE), .groups = "drop")
decomp_effects <- decomp_wide %>% left_join(overall_mean, by = "year") %>%
  mutate(across(-c(year, overall), ~ .x - ref_vals[cur_column()],
                .names = "within_{.col}")) %>%
  mutate(composition = overall - rowSums(select(., starts_with("within_")), na.rm = TRUE) -
           sum(ref_vals, na.rm = TRUE) / length(ref_vals))

fig5A <- ggplot(overall_mean, aes(year, overall)) +
  geom_line(color = "#E74C3C", linewidth = 1) + geom_point(size = 2) +
  labs(x = "Year", y = "Mean VFDB burden", title = "A. Annual mean VFDB burden") + theme_pub
within_df <- decomp_effects %>% select(year, starts_with("within_")) %>%
  pivot_longer(-year, names_to = "cluster", values_to = "effect") %>%
  mutate(cluster = str_remove(cluster, "within_"))
fig5B <- ggplot(within_df, aes(year, effect, color = cluster)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "Year", y = "Within-cluster effect", title = "B. Within-cluster enrichment effects") +
  scale_color_manual(values = cluster_colors(length(unique(within_df$cluster)))) + theme_pub
cluster_contrib <- within_df %>% filter(year == max(year)) %>%
  mutate(abs_effect = abs(effect)) %>% arrange(desc(abs_effect))
fig5C <- ggplot(cluster_contrib, aes(x = reorder(cluster, effect), y = effect,
  fill = effect > 0)) +
  geom_col(alpha = 0.8) + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"), guide = "none") +
  labs(x = "Cluster", y = "VFDB unit contribution (2025)",
       title = "C. Cluster-level contributions (2025)") + theme_pub
fig5 <- fig5A / (fig5B | fig5C) + plot_layout(heights = c(1, 1))
print(fig5)

# ============================================================
# FIGURE 6: Gene-level trajectories (Cluster_3)
# ============================================================
cat("\n--- Figure 6: Gene-level VFDB trajectories ---\n")

if (requireNamespace("scales", quietly = TRUE)) {
  cluster3 <- st69_master %>% filter(shell_cluster == "Cluster_3")
  vf_cols <- intersect(vf_genes_key, colnames(cluster3))
  if (length(vf_cols) > 0 && nrow(cluster3) > 0) {
    gene_traj <- cluster3 %>%
      select(genome_id, year, all_of(vf_cols)) %>%
      pivot_longer(-c(genome_id, year), names_to = "gene", values_to = "present") %>%
      mutate(present = as.integer(!is.na(present) & present != "" & present != "0")) %>%
      group_by(year, gene) %>% summarise(prev = mean(present) * 100, .groups = "drop")
    fig6 <- ggplot(gene_traj, aes(year, prev, color = gene)) +
      geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
      labs(x = "Year", y = "Prevalence (%)", color = "Gene",
           title = "Figure 6. Gene-level VFDB trajectories (Cluster_3)") +
      theme_pub + theme(legend.position = "right")
    print(fig6)
  }
}

# ============================================================
# FIGURE 7: Tree mapping + parsimony
# ============================================================
cat("\n--- Figure 7: Phylogeny + parsimony ---\n")

tree_file <- file.path(config$OUTPUT_DIR, "ST69", "ST69_bootstrap.treefile")
if (!file.exists(tree_file)) tree_file <- file.path(config$BASE_DIR, "ST69_bootstrap.treefile")
if (file.exists(tree_file) && requireNamespace("ape", quietly = TRUE)) {
  tree <- read.tree(tree_file)
  shell_assignments <- st69_master %>% select(genome_id, shell_cluster) %>%
    filter(genome_id %in% tree$tip.label) %>% distinct(genome_id, .keep_all = TRUE)
  if (nrow(shell_assignments) > 0) {
    trait <- setNames(shell_assignments$shell_cluster, shell_assignments$genome_id)
    missing_tips <- setdiff(tree$tip.label, names(trait))
    if (length(missing_tips) > 0) trait <- c(trait, setNames(rep("Unknown", length(missing_tips)), missing_tips))
    trait <- trait[tree$tip.label]
    par(mar = c(2, 2, 3, 1))
    plot(tree, type = "fan", tip.color = RColorBrewer::brewer.pal(4, "Set1")[as.integer(factor(trait))],
         cex = 0.3, main = "Figure 7. Shell-gene groups on core-genome phylogeny")
    legend("topright", legend = levels(factor(shell_assignments$shell_cluster)),
           col = RColorBrewer::brewer.pal(4, "Set1"), pch = 19, cex = 0.7, bty = "n")
  }
} else {
  cat("  Tree file not found or ape not installed — skipping Figure 7\n")
}

# ============================================================
# FIGURE 8: Clinical enrichment + sensitivity
# ============================================================
cat("\n--- Figure 8: Clinical trajectory + sensitivity ---\n")

cluster_by_year <- st69_master %>%
  group_by(year, shell_cluster) %>%
  summarise(n = n(), clinical_pct = mean(clinical_binary, na.rm = TRUE) * 100, .groups = "drop")
fig8A <- ggplot(cluster_by_year, aes(year, clinical_pct, color = shell_cluster)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  labs(x = "Year", y = "Clinical representation (%)",
       title = "A. Clinical representation by cluster over time") +
  scale_color_manual(values = cluster_colors(length(unique(st69_master$shell_cluster)))) + theme_pub

c3 <- st69_master %>% filter(shell_cluster == "Cluster_3") %>%
  group_by(year) %>% summarise(n = n(), n_clinical = sum(clinical_binary), .groups = "drop") %>%
  mutate(p_clinical = n_clinical / n)
if (nrow(c3) >= 3) {
  mod <- glm(cbind(n_clinical, n - n_clinical) ~ year, data = c3, family = binomial)
  or <- exp(coef(mod)["year"])
  ci <- exp(confint(mod)["year", ])
  sens_text <- paste0("OR/year = ", round(or, 3), "\n95% CI: ", round(ci[1], 3), "-", round(ci[2], 3))
} else { sens_text <- "" }
fig8B <- ggplot(c3, aes(year, p_clinical)) +
  geom_point(size = 2) + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  annotate("text", x = min(c3$year) + 1, y = max(c3$p_clinical) - 0.05,
           label = sens_text, hjust = 0, size = 3) +
  labs(x = "Year", y = "Clinical proportion", title = "B. Cluster_3 clinical trajectory") + theme_pub
fig8 <- fig8A | fig8B
print(fig8)

# ============================================================
# FIGURE 9: Clinical enrichment top genes + MGE context
# ============================================================
cat("\n--- Figure 9: Clinical enrichment + MGE context ---\n")

vf_cols_present <- intersect(vf_genes_key, colnames(st69_master))
if (length(vf_cols_present) > 0) {
  enrichment <- st69_master %>%
    select(clinical_binary, all_of(vf_cols_present)) %>%
    pivot_longer(-clinical_binary, names_to = "gene", values_to = "present") %>%
    mutate(present = as.integer(!is.na(present) & present != "" & present != "0")) %>%
    group_by(gene) %>%
    summarise(
      prev_clinical = mean(present[clinical_binary == 1], na.rm = TRUE),
      prev_nonclinical = mean(present[clinical_binary == 0], na.rm = TRUE),
      delta = (prev_clinical - prev_nonclinical) * 100,
      .groups = "drop"
    ) %>% arrange(desc(abs(delta))) %>% head(15)
  fig9A <- ggplot(enrichment, aes(x = reorder(gene, delta), y = delta, fill = delta > 0)) +
    geom_col(alpha = 0.8) + coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"), guide = "none") +
    labs(x = "Gene", y = "Prevalence difference (pp)", title = "A. Clinical enrichment (clinical - non-clinical)") +
    theme_pub
  print(fig9A)
} else {
  cat("  No VF columns found — skipping Figure 9\n")
}

# ============================================================
# SUPPLEMENTARY S1: Parsimony permutation test
# ============================================================
cat("\n--- Supplementary S1: Parsimony permutation test ---\n")

if (file.exists(tree_file) && requireNamespace("ape", quietly = TRUE) &&
    requireNamespace("phangorn", quietly = TRUE)) {
  tree <- read.tree(tree_file)
  shell_assignments <- st69_master %>% select(genome_id, shell_cluster) %>%
    filter(genome_id %in% tree$tip.label) %>% distinct(genome_id, .keep_all = TRUE)
  if (nrow(shell_assignments) > 0) {
    trait <- setNames(shell_assignments$shell_cluster, shell_assignments$genome_id)
    missing_tips <- setdiff(tree$tip.label, names(trait))
    if (length(missing_tips) > 0) trait <- c(trait, setNames(rep("Unknown", length(missing_tips)), missing_tips))
    trait <- trait[tree$tip.label]

    parsimony_score <- function(tree, tv) {
      phangorn::parsimony(tree, phangorn::phyDat(tv[!is.na(tv) & tv != "Unknown"], type = "USER"))
    }
    obs <- parsimony_score(tree, trait)
    N_PERM <- 999
    set.seed(42)
    perm_scores <- replicate(N_PERM, {
      perm_trait <- sample(trait)
      names(perm_trait) <- names(trait)
      parsimony_score(tree, perm_trait)
    })
    p_val <- (sum(perm_scores <= obs) + 1) / (N_PERM + 1)

    perm_df <- tibble(score = perm_scores)
    figS1 <- ggplot(perm_df, aes(x = score)) +
      geom_histogram(binwidth = 10, fill = "grey70", color = "grey40", alpha = 0.8) +
      geom_vline(xintercept = obs, color = "#E74C3C", linewidth = 1.2, linetype = "dashed") +
      annotate("text", x = obs, y = Inf, label = paste0("Observed = ", obs, "\np = ", round(p_val, 4)),
               vjust = 2, hjust = -0.1, color = "#E74C3C", fontface = "bold") +
      labs(x = "Parsimony score", y = "Frequency",
           title = paste0("Supplementary S1. Parsimony permutation test (p = ", round(p_val, 4), ")")) +
      theme_pub
    print(figS1)
  }
}

# ============================================================
# SUPPLEMENTARY S2: Allelic conversion scatter
# ============================================================
cat("\n--- Supplementary S2: Allelic conversion ---\n")

s2_file <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_rigorous_analysis",
  "VFDB_vs_PP_rigorous.csv")
if (file.exists(s2_file)) {
  s2_data <- read_csv(s2_file, show_col_types = FALSE)
  if (all(c("VFDB", "PP") %in% colnames(s2_data))) {
    figS2 <- ggplot(s2_data, aes(VFDB, PP)) +
      geom_point(alpha = 0.3, size = 0.8, color = "#4ECDC4") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "lm", se = TRUE, color = "#2C3E50") +
      labs(x = "VFDB identity (%)", y = "PPanGGOLiN identity (%)",
           title = "Supplementary S2. Allelic conversion: VFDB vs PPanGGOLiN") +
      coord_equal() + theme_pub
    print(figS2)
  }
}

# ============================================================
# Close PDF
# ============================================================
dev.off()
cat("\n============================================\n")
cat(" DONE\n")
cat(" Output:", pdf_file, "\n")
cat("============================================\n")
