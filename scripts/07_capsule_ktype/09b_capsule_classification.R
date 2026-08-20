#!/usr/bin/env Rscript
# Reviewer response: Capsule classification (G2 vs G3) on clinical enrichment + temporal increase
# Uses VFDB summary + master join matching the publication
suppressPackageStartupMessages({
  library(tidyverse); library(writexl); library(data.table); library(broom)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_capsule_classification")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
    grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
    TRUE ~ paste0("Escherichia_coli_", x))
}
to_bin_vfdb <- function(x) as.integer(sapply(strsplit(as.character(x), ";"),
  function(v) any(suppressWarnings(as.numeric(v) >= 90), na.rm = TRUE)))

# ---- 1. Load VFDB summary (same source as publication Figure_7) ----
cat("Loading VFDB summary...\n")
vfdb_raw <- read_tsv(config$st_vfdb_summary("ST69"), show_col_types = FALSE,
  progress = FALSE, col_types = cols(.default = "c")) %>%
  filter(str_detect(.data[["#FILE"]], "^ST69/")) %>%
  mutate(genome_id = make_genome_id(str_remove(.data[["#FILE"]], "^ST69/")),
    genome_id = str_remove(genome_id, "_vfdb\\.tsv$"))
genes <- setdiff(colnames(vfdb_raw), c("#FILE", "NUM_FOUND", "genome_id"))
bin <- vfdb_raw %>% select(genome_id, all_of(genes)) %>%
  mutate(across(all_of(genes), to_bin_vfdb))
cat("VFDB summary:", nrow(bin), "genomes\n")

# ---- 2. Load master table for metadata ----
cat("Loading master table...\n")
master <- fread(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv"),
  sep = ",", header = TRUE, data.table = FALSE, na.strings = "")
meta <- master %>% select(genome_id, shell_cluster, year, clinical_binary)
cat("Master table:", nrow(meta), "genomes\n")

# ---- 3. Join (same as publication figure_combined_summary.R) ----
df <- meta %>% inner_join(bin, by = "genome_id") %>%
  mutate(year = as.integer(year), clin = as.integer(clinical_binary)) %>%
  filter(!is.na(clin), !is.na(year), year >= 2016, year <= 2025)
cat("Joined:", nrow(df), "genomes\n")

# ---- 4. Capsule classification ----
kps_cols <- intersect(grep("^(kps|neu)", genes, value = TRUE), colnames(df))
G2_markers <- intersect(c("kpsC", "kpsS"), kps_cols)
G3_markers <- intersect(c("kpsE", "kpsM", "kpsT"), kps_cols)

df$capsule_type <- apply(df[, kps_cols, drop = FALSE], 1, function(r) {
  has_g2 <- any(r[names(r) %in% G2_markers] >= 1)
  has_g3 <- any(r[names(r) %in% G3_markers] >= 1)
  any_kps <- any(r >= 1)
  if (has_g2) return("G2")
  if (has_g3) return("G3")
  if (any_kps) return("Unclassified")
  return("No capsule")
})

caps_dist <- df %>% count(capsule_type) %>% mutate(pct = n / sum(n) * 100)
cat("\n=== Capsule Distribution ===\n")
print(caps_dist)

# ---- 5. Key genes for clinical enrichment + temporal increase ----
key_genes <- intersect(c("kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT","kpsU",
  "papX","papF","papB","papG","sat","hlyA"), genes)

# Clinical enrichment (delta prevalence + OR, matching Figure 7B)
cat("\n=== Clinical Enrichment (all genomes) ===\n")
clin_enrich <- map_dfr(key_genes, ~{
  g <- .x
  d <- df %>% filter(!is.na(.data[[g]]))
  if (sum(d[[g]], na.rm = TRUE) < 5) return(NULL)
  clin_prev <- mean(d[[g]][d$clin == 1]) * 100
  non_prev <- mean(d[[g]][d$clin == 0]) * 100
  ft <- tryCatch(fisher.test(table(d$clin, d[[g]])), error = function(e) NULL)
  or <- if (!is.null(ft)) unname(ft$estimate) else NA
  p <- if (!is.null(ft)) ft$p.value else NA
  # Country-adjusted logistic regression
  m <- tryCatch(glm(.data[[g]] ~ clin + country + year,
    data = df, family = binomial), error = function(e) NULL)
  adj_or <- if (!is.null(m)) exp(coef(m)["clin"]) else NA
  adj_p <- if (!is.null(m)) summary(m)$coefficients["clin", 4] else NA
  tibble(gene = g, clin_prev = clin_prev, non_prev = non_prev,
    delta_pp = clin_prev - non_prev, or = or, p = p, adj_or = adj_or, adj_p = adj_p)
}) %>% mutate(p_adj = p.adjust(pmax(p, 0, na.rm = TRUE), method = "BH"),
  enriched = !is.na(adj_p) & adj_p < 0.05 & adj_or > 1)
print(clin_enrich %>% select(gene, delta_pp, or, adj_or, enriched), n = 30)

# Temporal increase in Cluster 3 (matching Figure 7A)
cat("\n=== Cluster 3: Increasing genes (late vs early) ===\n")
EARLY <- c(2016, 2017, 2018)
LATE <- c(2022, 2023, 2024, 2025)
cd3 <- df %>% filter(shell_cluster == "Cluster_3")
ea <- cd3 %>% filter(year %in% EARLY)
la <- cd3 %>% filter(year %in% LATE)
n_early <- nrow(ea)
n_late <- nrow(la)

temporal_inc <- map_dfr(key_genes, ~{
  g <- .x
  prev_e <- mean(ea[[g]], na.rm = TRUE) * 100
  prev_l <- mean(la[[g]], na.rm = TRUE) * 100
  dp <- prev_l - prev_e
  mat <- matrix(c(sum(ea[[g]]), n_early - sum(ea[[g]]),
    sum(la[[g]]), n_late - sum(la[[g]])), nrow = 2)
  ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
  p <- if (!is.null(ft)) ft$p.value else NA
  tibble(gene = g, early_prev = prev_e, late_prev = prev_l, delta_pp = dp, p = p)
}) %>% mutate(p_adj = p.adjust(pmax(p, 0, na.rm = TRUE), method = "BH"),
  increasing = delta_pp >= 5 & late_prev >= 5 & p_adj < 0.05)
print(temporal_inc %>% select(gene, early_prev, late_prev, delta_pp, increasing), n = 30)

# ---- 6. G2/G3 overlay on capsule genes ----
cat("\n=== G2/G3 breakdown among increasing + clinically enriched capsule genes ===\n")
capsule_genes <- intersect(c("kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT","kpsU"), genes)

# Clinical enrichment by capsule type
clin_by_cap <- df %>% group_by(capsule_type) %>%
  summarise(n = n(), clin_pct = mean(clin) * 100, .groups = "drop")
cat("\nClinical % by capsule type:\n")
print(clin_by_cap)

# Per-gene prevalence by capsule type
cat("\nkps genes within each capsule type:\n")
for (g in capsule_genes) {
  d <- df %>% group_by(capsule_type) %>%
    summarise(pct = mean(.data[[g]]) * 100, n = n(), .groups = "drop")
  cat(g, ":\n", capture.output(print(d, n = 4)), "\n")
}

# ---- 7. Plots ----
# ---- 5b. Capsule temporal trends (all genomes) ----
EARLY <- c(2016, 2017, 2018)
LATE <- c(2022, 2023, 2024, 2025)
caps_temporal_all <- df %>% mutate(period = case_when(year %in% EARLY ~ "early",
  year %in% LATE ~ "late", TRUE ~ "mid")) %>%
  filter(period %in% c("early", "late")) %>%
  group_by(period) %>% mutate(period_total = n()) %>% ungroup() %>%
  group_by(capsule_type, period) %>%
  summarise(n = n(), period_total = first(period_total), .groups = "drop") %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early)
cat("\n=== Capsule temporal (all genomes) ===\n")
print(caps_temporal_all)

# Capsule temporal in Cluster 3
cd3 <- df %>% filter(shell_cluster == "Cluster_3")
caps_temporal_c3 <- cd3 %>% mutate(period = case_when(year %in% EARLY ~ "early",
  year %in% LATE ~ "late", TRUE ~ "mid")) %>%
  filter(period %in% c("early", "late")) %>%
  group_by(period) %>% mutate(period_total = n()) %>% ungroup() %>%
  group_by(capsule_type, period) %>%
  summarise(n = n(), period_total = first(period_total), .groups = "drop") %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early)
cat("\n=== Capsule temporal (Cluster 3) ===\n")
print(caps_temporal_c3)

caps_colors <- c("G2" = "#d73027", "G3" = "#4575b4", "No capsule" = "#808080", "Unclassified" = "#fddbc7")

# A: Temporal increase in Cluster 3 (like Fig 7A)
pA <- temporal_inc %>% filter(increasing | gene %in% capsule_genes) %>%
  mutate(gene = fct_reorder(gene, delta_pp)) %>%
  ggplot(aes(x = delta_pp, y = gene, fill = gene %in% G2_markers)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
    labels = c("G3/other", "G2 marker"), name = "") +
  labs(x = "Delta prevalence (pp), late - early", y = NULL,
    title = "Cluster_3: increasing VFDB genes (colored by G2/G3)") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT, "Cluster3_increasing_genes_by_G2.pdf"), pA, width = 8, height = 6)

# B: Clinical enrichment (like Fig 7B)
pB <- clin_enrich %>% filter(enriched | gene %in% capsule_genes) %>%
  mutate(gene = fct_reorder(gene, delta_pp)) %>%
  ggplot(aes(x = delta_pp, y = gene, fill = gene %in% G2_markers)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
    labels = c("G3/other", "G2 marker"), name = "") +
  labs(x = "Delta prevalence (pp), clinical - non-clinical", y = NULL,
    title = "Clinical enrichment (all genomes, colored by G2/G3)") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT, "clinical_enrichment_by_G2.pdf"), pB, width = 8, height = 6)

# C: Capsule distribution
pC <- ggplot(caps_dist, aes(x = "", y = pct, fill = capsule_type)) +
  geom_bar(stat = "identity", width = 0.5) + coord_polar("y") +
  scale_fill_manual(values = caps_colors) +
  theme_minimal() + labs(title = "ST69 Capsule Types", fill = "") +
  theme(axis.text = element_blank(), axis.title = element_blank())
ggsave(file.path(OUT, "capsule_distribution.pdf"), pC, width = 5, height = 4)

# D: OR vs delta prevalence scatter
pD <- clin_enrich %>% filter(enriched | gene %in% capsule_genes) %>%
  mutate(g2 = gene %in% G2_markers) %>%
  ggplot(aes(x = delta_pp, y = or, color = g2, size = -log10(pmax(p_adj, 1e-300)))) +
  geom_point(alpha = 0.8) + scale_y_log10() +
  scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
    labels = c("G3/other", "G2 marker"), name = "") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Delta prevalence (pp)", y = "Odds ratio (log scale)",
    title = "Clinical enrichment: effect size vs association", size = "−log10(p)") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT, "clinical_enrichment_OR_vs_delta.pdf"), pD, width = 8, height = 6)

# ---- 8. Write ----
caps_dist_enhanced <- caps_dist %>% left_join(
  caps_temporal_all %>% select(capsule_type, n_early, n_late, prev_early, prev_late, delta_pp),
  by = "capsule_type")
cat("\n=== Enhanced capsule distribution ===\n")
print(caps_dist_enhanced)

write_xlsx(list(
  capsule_distribution = caps_dist_enhanced,
  capsule_temporal_all = caps_temporal_all %>% arrange(desc(delta_pp)),
  capsule_temporal_cluster3 = caps_temporal_c3 %>% arrange(desc(delta_pp)),
  clinical_enrichment = clin_enrich %>% arrange(p_adj),
  cluster3_increasing = temporal_inc %>% arrange(p_adj),
  clinical_by_capsule_type = clin_by_cap
), path = file.path(OUT, "capsule_classification.xlsx"))

cat("\nOutput in:", OUT, "\n")
