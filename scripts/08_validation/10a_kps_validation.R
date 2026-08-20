#!/usr/bin/env Rscript
# Reviewer response: kps locus-level validation
# Compare VFDB-based kps detection vs PPanGGOLiN gene family-based detection
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(writexl); library(data.table)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_kps_validation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

# ---- 1. Load cluster assignments ----
assignments <- read_csv(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "03_shell_gene_cluster_assignments_k4.csv"), show_col_types = FALSE) %>%
  mutate(genome_id = as.character(genome_id))

# ---- 2. Load VFDB kps genes ----
cat("Loading VFDB summary...\n")
vfdb_raw <- fread(config$st_vfdb_summary("ST69"), sep = "\t", header = TRUE,
  data.table = FALSE, colClasses = "character", na.strings = "")
vfdb_raw$genome_id <- make_genome_id(str_remove(vfdb_raw[[1]], "^ST69/"))
vfdb_raw$genome_id <- str_remove(vfdb_raw$genome_id, "_vfdb\\.tsv$")

# Identify kps genes in VFDB
kps_vfdb_pattern <- "^(kps|neu)"
vf_genes <- setdiff(colnames(vfdb_raw), c("#FILE", "NUM_FOUND", "genome_id"))
kps_genes_vfdb <- grep(kps_vfdb_pattern, vf_genes, value = TRUE)
cat("VFDB kps genes found:", length(kps_genes_vfdb), "\n")
cat("  ", paste(kps_genes_vfdb, collapse = ", "), "\n")

# Binary matrix for VFDB kps
to_bin <- function(x) as.integer(sapply(strsplit(as.character(x), ";"),
  function(v) any(suppressWarnings(as.numeric(v) >= 90), na.rm = TRUE)))
vfdb_kps <- vfdb_raw %>% select(genome_id, all_of(kps_genes_vfdb)) %>%
  mutate(across(all_of(kps_genes_vfdb), to_bin))

cat("VFDB genomes:", nrow(vfdb_kps), "\n")

# ---- 3. Load PPanGGOLiN gene families ----
cat("Loading PPanGGOLiN gene_presence_absence.Rtab...\n")
rtab <- fread(file.path(config$PANGENOME_DIR, "gene_presence_absence.Rtab"),
  sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

# Load VF family annotations to map gene names to family IDs
vf_families <- read_csv(file.path(config$PANGENOME_DIR, "vf_module_rgp_results",
  "tables", "gene_family_annotation_VF_only.csv"), show_col_types = FALSE)
cat("VF-annotated families:", nrow(vf_families), "\n")

# Identify kps families in PPanGGOLiN
kps_families <- vf_families %>% filter(vf_system == "kps") %>%
  distinct(family_id, gene_name)
cat("PPanGGOLiN kps families:", nrow(kps_families), "\n")
print(table(kps_families$gene_name))

# Also check what gene families map to the same VFDB kps gene names
kps_family_map <- kps_families %>%
  mutate(gene_base = tolower(str_extract(gene_name, "^[a-z]+[0-9]*")))

# Get the Rtab columns (genomes)
rtab_genomes <- colnames(rtab)[-1]
common <- intersect(rtab_genomes, vfdb_kps$genome_id)
cat("Common genomes:", length(common), "\n")

# PPanGGOLiN-based kps: which families are present in each genome
rtab_family_ids <- rtab[[1]]
kps_family_ids <- intersect(kps_families$family_id, rtab_family_ids)
cat("KPS families in Rtab:", length(kps_family_ids), "\n")

# Build PPanGGOLiN kps matrix
kps_idx <- which(rtab_family_ids %in% kps_family_ids)
if (length(kps_idx) > 0) {
  pp_kps_mat <- rtab[kps_idx, common, drop = FALSE]
  pp_kps_presence <- apply(pp_kps_mat[, -1], 2, function(x) as.integer(any(as.numeric(x) > 0, na.rm = TRUE)))
} else {
  pp_kps_presence <- setNames(rep(0L, length(common)), common)
}

# ---- 4. Compare per-genome kps prevalence ----
comparison <- tibble(genome_id = common,
  vfdb_kps = rowSums(vfdb_kps[match(common, vfdb_kps$genome_id), kps_genes_vfdb, drop = FALSE], na.rm = TRUE),
  pp_kps = pp_kps_presence[match(common, names(pp_kps_presence))],
  vfdb_kps_binary = as.integer(vfdb_kps > 0))

# Count per-gene prevalence in VFDB
kps_gene_prev_vfdb <- vfdb_kps %>%
  filter(genome_id %in% common) %>%
  summarise(across(all_of(kps_genes_vfdb), ~ mean(.x, na.rm = TRUE) * 100)) %>%
  pivot_longer(everything(), names_to = "gene", values_to = "prevalence_vfdb_pct")

cat("\nVFDB kps gene prevalence:\n")
print(kps_gene_prev_vfdb %>% arrange(desc(prevalence_vfdb_pct)), n = Inf)

# Per-gene family prevalence in PPanGGOLiN
if (length(kps_family_ids) > 0) {
  kps_family_prev <- tibble(
    family_id = kps_family_ids,
    prevalence_pp_pct = rowMeans(rtab[match(kps_family_ids, rtab_family_ids), common, drop = FALSE] > 0, na.rm = TRUE) * 100
  ) %>% left_join(kps_families %>% distinct(family_id, gene_name), by = "family_id")
  cat("\nPPanGGOLiN kps family prevalence:\n")
  print(kps_family_prev %>% arrange(desc(prevalence_pp_pct)), n = Inf)
}

# ---- 5. Kps module analysis ----
kps_modules <- read_csv(file.path(config$PANGENOME_DIR, "vf_module_rgp_results",
  "tables", "VF_families_in_functional_modules.csv"), show_col_types = FALSE) %>%
  filter(vf_system == "kps")

cat("\nkps module summary:\n")
print(table(kps_modules$module_id))

# ---- 6. Group 2 vs Group 3 capsule assessment ----
# kpsD, kpsE, kpsF, kpsM, kpsT are shared across group 2 capsules
# neu genes are specific to K1 (group 2)
# Different kps gene families may indicate different capsule types
kps_family_details <- kps_families %>%
  group_by(gene_name) %>%
  summarise(n_families = n_distinct(family_id), .groups = "drop")

cat("\nkps gene -> family count (allelic diversity):\n")
print(kps_family_details)

# Count distinct PPanGGOLiN families per kps gene
# Multiple families per gene name = allelic variation
multi_family <- kps_family_details %>% filter(n_families > 1)
cat("\nkps genes with multiple PPanGGOLiN families (allelic variation):\n")
if (nrow(multi_family) > 0) print(multi_family)

# ---- 7. Cluster-specific kps prevalence ----
comparison <- comparison %>% left_join(assignments %>% select(genome_id, shell_cluster), by = "genome_id")

kps_by_cluster <- comparison %>%
  filter(!is.na(shell_cluster)) %>%
  group_by(shell_cluster) %>%
  summarise(
    n = n(),
    vfdb_kps_mean = mean(vfdb_kps, na.rm = TRUE),
    vfdb_kps_pct = mean(vfdb_kps_binary, na.rm = TRUE) * 100,
    pp_kps_pct = mean(pp_kps, na.rm = TRUE) * 100,
    .groups = "drop"
  )
cat("\nkps prevalence by cluster:\n")
print(kps_by_cluster)

# ---- 8. VFDB-only kps vs PPanGGOLiN kps comparison ----
# The key question: is VFDB under-calling kps?
discrepancy <- comparison %>%
  mutate(
    vfdb_positive = vfdb_kps_binary == 1,
    pp_positive = pp_kps == 1,
    vfdb_only = vfdb_positive & !pp_positive,
    pp_only = pp_positive & !vfdb_positive,
    both = vfdb_positive & pp_positive
  )

cat("\nDiscrepancy analysis:\n")
cat("  VFDB+ PPanGGOLiN+ (both):", sum(discrepancy$both), "\n")
cat("  VFDB+ PPanGGOLiN- (VFDB only):", sum(discrepancy$vfdb_only), "\n")
cat("  VFDB- PPanGGOLiN+ (PP only):", sum(discrepancy$pp_only), "\n")
cat("  VFDB- PPanGGOLiN- (neither):", sum(!discrepancy$vfdb_positive & !discrepancy$pp_positive), "\n")

# ---- 9. Literature comparison ----
# Gladstone et al. 2026: <3% of ST69 lack kps
pct_lacking_kps <- mean(comparison$pp_kps == 0, na.rm = TRUE) * 100
cat(sprintf("\n%% lacking kps (PPanGGOLiN): %.1f%%\n", pct_lacking_kps))
pct_lacking_kps_vfdb <- mean(comparison$vfdb_kps_binary == 0, na.rm = TRUE) * 100
cat(sprintf("%% lacking kps (VFDB): %.1f%%\n", pct_lacking_kps_vfdb))

# ---- 10. Write outputs ----
write_xlsx(list(
  kps_gene_prevalence_VFDB = kps_gene_prev_vfdb,
  kps_family_prevalence_PPanGGOLiN = if (exists("kps_family_prev")) kps_family_prev else tibble(),
  kps_by_cluster = kps_by_cluster,
  discrepancy = discrepancy %>% select(genome_id, vfdb_kps, pp_kps, vfdb_kps_binary, vfdb_positive, pp_positive, shell_cluster),
  kps_module_summary = as.data.frame(table(kps_modules$module_id)),
  kps_gene_family_details = kps_family_details
), path = file.path(OUT, "kps_validation.xlsx"))

# ---- 11. Summary text ----
sink(file.path(OUT, "kps_validation_summary.txt"))
cat("kps Locus-Level Validation Summary\n")
cat("==================================\n\n")
cat(sprintf("VFDB kps genes detected: %d\n", length(kps_genes_vfdb)))
cat(sprintf("PPanGGOLiN kps families: %d\n", nrow(kps_families)))
cat(sprintf("PPanGGOLiN kps families in Rtab: %d\n", length(kps_family_ids)))
cat(sprintf("Shared genomes: %d\n\n", length(common)))
cat("VFDB kps gene prevalence:\n")
print(kps_gene_prev_vfdb %>% arrange(desc(prevalence_vfdb_pct)), n = Inf)
cat("\n")
cat(sprintf("kps prevalence by method:\n"))
cat(sprintf("  VFDB: %.1f%% genomes carry at least one kps gene\n", mean(comparison$vfdb_kps_binary)*100))
cat(sprintf("  PPanGGOLiN: %.1f%% genomes carry at least one kps family\n", mean(comparison$pp_kps)*100))
cat(sprintf("\nLiterature (Gladstone et al. 2026): <3%% of ST69 lack kps\n"))
cat(sprintf("Our PPanGGOLiN estimate: %.1f%% lack kps\n", pct_lacking_kps))
cat(sprintf("Our VFDB estimate: %.1f%% lack kps\n\n", pct_lacking_kps_vfdb))
cat("Interpretation:\n")
if (pct_lacking_kps < 5) {
  cat("  PPanGGOLiN-based kps detection approaches the literature estimate.\n")
} else if (pct_lacking_kps > 10) {
  cat("  PPanGGOLiN-based detection still shows more genomes lacking kps than expected.\n")
  cat("  This may reflect divergent kps alleles not captured in current gene families.\n")
}
if (pct_lacking_kps_vfdb > 10) {
  cat("  VFDB-based detection misses more kps-positive genomes than PPanGGOLiN.\n")
  cat("  This is consistent with divergent capsular loci (group 2/3) not represented in VFDB.\n")
}
sink()

cat("\nDone. Output in:", OUT, "\n")
