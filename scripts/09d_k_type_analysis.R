#!/usr/bin/env Rscript
# Reviewer: K-type inference from VirulenceFinder kpsM alleles in ST69
# Addresses question: is temporal increase K52/K54/K96-specific?
suppressPackageStartupMessages({library(tidyverse); library(writexl); library(readxl); library(data.table)})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, config$TARGET_ST, "reviewer_capsule_classification")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load VF binary matrix (sample-level kpsM alleles) ----
cat("Loading VirulenceFinder binary matrix...\n")
vf <- fread(file.path(config$INPUT_DIR,
  "virulencefinder_summary", "virulencefinder_binary_matrix.tsv"),
  sep = "\t", header = TRUE, data.table = FALSE, na.strings = "")
cat("Rows:", nrow(vf), "Cols:", ncol(vf), "\n")
# genome ID column
genome_col <- if ("Genome" %in% colnames(vf)) "Genome" else "genome"
if ("st" %in% colnames(vf)) {
  vf <- vf %>% filter(.data[["st"]] == config$TARGET_ST)
  cat("Filtered to", config$TARGET_ST, "-> Rows:", nrow(vf), "\n")
}

# Also load VF frequency to know total
vf_freq <- fread(file.path(config$INPUT_DIR,
  "virulencefinder_summary", "virulencefinder_gene_frequency.tsv"),
  sep = "\t", header = TRUE, data.table = FALSE)
if ("st" %in% colnames(vf_freq)) {
  vf_freq <- vf_freq %>% filter(.data[["st"]] == config$TARGET_ST)
}

# ---- 2. Extract kpsM allele columns ----
k_cols <- grep("kpsM[I]*[_]", colnames(vf), value = TRUE)
k_cols <- setdiff(k_cols, c("kpsMII", "kpsMIII"))
cat("kpsM K-type alleles:", length(k_cols), "\n")
print(k_cols)

# K-type mapping
k_type_map <- c(
  kpsMII_K1 = "K1", kpsMII_K5 = "K5", kpsMII_K52 = "K52",
  kpsMII_K4 = "K4", kpsMII_K23 = "K23", kpsM_K15 = "K15",
  kpsM_K11 = "K11", kpsM_K19 = "K19", kpsM_K19K23 = "K19/K23",
  kpsMIII_K96 = "K96", kpsMIII_K98 = "K98", kpsMIII_K10 = "K10"
)
present_k <- intersect(k_cols, names(k_type_map))
k_type_map <- k_type_map[present_k]

# Binary columns are character ("1" or "." or "0")
is_present <- function(x) {
  x <- as.character(x); x[is.na(x)] <- "."; x == "1"
}
# Assign K-type: first column name -> type
ktype_df <- vf %>% select(genome, all_of(present_k))
colnames(ktype_df)[1] <- "Genome"

# Per-genome K-type assignment (mutually exclusive, first match)
ktype_df$k_type <- apply(ktype_df[, present_k, drop = FALSE], 1, function(r) {
  hits <- names(r)[is_present(r)]
  if (length(hits) == 0) return(NA_character_)
  unname(k_type_map[hits[1]])
})
# For genomes with no specific allele but generic kpsMII/MIII, assign broadly
has_g2 <- is_present(vf[["kpsMII"]])
has_g3 <- is_present(vf[["kpsMIII"]])
ktype_df$k_type <- ifelse(is.na(ktype_df$k_type) & has_g2, "G2-unknown", ktype_df$k_type)
ktype_df$k_type <- ifelse(is.na(ktype_df$k_type) & has_g3, "G3-unknown", ktype_df$k_type)
ktype_df$k_type <- ifelse(is.na(ktype_df$k_type), "No kpsM allele", ktype_df$k_type)

cat(sprintf("\n=== K-type distribution in %s ===\n", config$TARGET_ST))
k_tab <- table(ktype_df$k_type, useNA = "ifany")
print(data.frame(K_type = names(k_tab), n = as.integer(k_tab), pct = round(as.integer(k_tab)/sum(k_tab)*100, 1)))

# ---- 3. Merge with metadata ----
# VirulenceFinder has "Genome" column - need to match with master
# Check format
cat("\nGenome format (first 3):\n")
print(head(ktype_df$Genome, 3))
# Probably Escherichia_coli_XXXXX format
meta <- fread(file.path(config$OUTPUT_DIR, config$TARGET_ST, "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv"),
  sep = ",", header = TRUE, data.table = FALSE, na.strings = "")
cat("\nMaster genome_id format (first 3):\n")
print(head(meta$genome_id, 3))

df <- ktype_df %>% inner_join(meta, by = c("Genome" = "genome_id")) %>%
  mutate(year = as.integer(year), clin = as.integer(clinical_binary)) %>%
  filter(!is.na(clin), !is.na(year))
cat("Merged:", nrow(df), "genomes\n")

# ---- 4. Temporal trends by K-type (all genomes) ----
EARLY <- c(2016, 2017, 2018)
LATE <- c(2022, 2023, 2024, 2025)
k_temporal <- df %>% filter(year %in% c(EARLY, LATE)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late")) %>%
  group_by(period) %>% mutate(period_total = n()) %>% ungroup() %>%
  group_by(k_type, period) %>%
  summarise(n = n(), period_total = first(period_total), .groups = "drop") %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early) %>%
  arrange(desc(delta_pp))
cat("\n=== K-type temporal trends (all genomes) ===\n")
print(as.data.frame(k_temporal))

# ---- 5. Temporal trends by K-type (Cluster 3 only) ----
c3 <- df %>% filter(shell_cluster == "Cluster_3")
k_temporal_c3 <- c3 %>% filter(year %in% c(EARLY, LATE)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late")) %>%
  group_by(period) %>% mutate(period_total = n()) %>% ungroup() %>%
  group_by(k_type, period) %>%
  summarise(n = n(), period_total = first(period_total), .groups = "drop") %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early) %>%
  arrange(desc(delta_pp))
cat("\n=== K-type temporal trends (Cluster 3) ===\n")
print(as.data.frame(k_temporal_c3))

# ---- 6. Clinical enrichment by K-type ----
k_clin <- df %>%
  group_by(k_type) %>%
  summarise(n = n(), clin_pct = mean(clin) * 100, .groups = "drop")
cat("\n=== K-type clinical enrichment ===\n")
print(as.data.frame(k_clin))

# ---- 7. G2/G3 overlay with K-type ----
# Map K-type to G2/G3
g_map <- c("K1" = "G2", "K5" = "G2", "K52" = "G2", "K4" = "G2", "K15" = "G2",
  "K11" = "G2", "K19" = "G2", "K19/K23" = "G2", "G2-unknown" = "G2",
  "K96" = "G3", "K98" = "G3", "K10" = "G3", "G3-unknown" = "G3",
  "No kpsM allele" = "No kpsM")
df$k_group <- g_map[df$k_type]

# Temporal in C3 by k_group
k_grp_c3 <- c3 %>% filter(year %in% c(EARLY, LATE)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late"),
    k_group = g_map[k_type]) %>%
  group_by(period) %>% mutate(period_total = n()) %>% ungroup() %>%
  group_by(k_group, period) %>%
  summarise(n = n(), period_total = first(period_total), .groups = "drop") %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early)
cat("\n=== K-type group temporal (Cluster 3) ===\n")
print(k_grp_c3)

# Clinical by k_group
k_clin_grp <- df %>% group_by(k_group) %>%
  summarise(n = n(), clin_pct = mean(clin) * 100, .groups = "drop")
cat("\n=== K-type group clinical ===\n")
print(k_clin_grp)

# ---- 8. K-type within G2 genomes in C3 ----
c3_g2 <- c3 %>% filter(k_type %in% c("K52", "K5", "K1", "K15", "K11", "G2-unknown"))
cat("\n=== G2 K-type prevalence in C3 (early vs late) ===\n")
g2_c3 <- c3_g2 %>% filter(year %in% c(EARLY, LATE)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late")) %>%
  group_by(k_type, period) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(period) %>% mutate(period_total = sum(n)) %>% ungroup() %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early)
print(g2_c3)

# ---- 9. Output ----
xlsx_file <- file.path(OUT, "capsule_classification.xlsx")
existing_sheets <- list()
for (s in c("capsule_distribution", "capsule_temporal_all", "capsule_temporal_cluster3",
  "clinical_enrichment", "cluster3_increasing", "clinical_by_capsule_type")) {
  tryCatch({
    existing_sheets[[s]] <- read_xlsx(xlsx_file, sheet = s)
  }, error = function(e) NULL)
}
existing[["k_type_distribution"]] <- k_dist
existing[["k_type_temporal"]] <- k_temporal
existing[["k_type_temporal_cluster3"]] <- k_temporal_c3
existing[["k_type_clinical"]] <- k_clin
existing[["k_type_group_clinical"]] <- k_clin_grp
existing[["k_type_group_c3_temporal"]] <- k_grp_c3
existing[["g2_k_type_c3_temporal"]] <- g2_c3
write_xlsx(existing_sheets, xlsx_file)
cat("\nAppended to:", xlsx_file, "\n")

# Summary for reviewer
cat("\n\n=== SUMMARY FOR REVIEWER ===\n")
cat("K54: NOT detected in ST69 (not in VirulenceFinder DB)\n")
cat("K96 (G3):", round(k_dist$pct[k_dist$k_type == "K96"], 1), "% of ST69\n")
cat("K52 (G2):", round(k_dist$pct[k_dist$k_type == "K52"], 1), "% of ST69\n")
cat("G2 types (K1/K5/K52/K15/K11): combined", 
  round(sum(k_dist$pct[k_dist$k_type %in% c("K1","K5","K52","K15","K11")]), 1), "%\n")
cat("G3 types (K96/K98/K10): combined",
  round(sum(k_dist$pct[k_dist$k_type %in% c("K96","K98","K10")]), 1), "%\n")
