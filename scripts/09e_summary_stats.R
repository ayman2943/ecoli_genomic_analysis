#!/usr/bin/env Rscript
# Reviewer response: summary statistics at each filtering step
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(writexl); library(data.table)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_summary_stats")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
    TRUE ~ paste0("Escherichia_coli_", x))
}

# ---- 1. PPanGGOLiN genomes (gene_presence_absence.Rtab) ----
cat("Loading PPanGGOLiN gene presence/absence...\n")
rtab <- fread(config$GENE_PA_AB, sep = "\t", header = TRUE, data.table = FALSE)
n_panx_genomes <- nrow(rtab)
n_panx_genes <- ncol(rtab) - 1
cat("PPanGGOLiN genomes:", n_panx_genomes, "\n")
cat("PPanGGOLiN gene families:", n_panx_genes, "\n")

# ---- 2. Metadata ----
meta_file <- config$st_metadata("ST69")
meta <- read_excel(meta_file, col_types = "text")
n_meta <- nrow(meta)
cat("Metadata entries:", n_meta, "\n")

# ---- 3. Metadata Name column ----
name_col <- {
  candidates <- c("Name", "genome", "strain", "isolate", "assembly")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[2]
}
cat("Metadata name column:", name_col, "\n")
cat("First few values:", paste(head(as.character(meta[[name_col]]), 3), collapse = ", "), "\n")
# Note: PPanGGOLiN IDs are contig-based (Escherichia_coli_XXXX_00001) while metadata IDs
# are sample-based. These identifiers differ by design and cannot be directly matched.

# ---- 4. Metadata with years ----
year_col <- {
  candidates <- c("Collection Year", "Collection_Year", "year", "Year")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[7]
}
meta_with_year <- meta %>% filter(!is.na(.data[[year_col]]))
n_with_year <- nrow(meta_with_year)
cat("Metadata with year:", n_with_year, "\n")

# ---- 5. VFDB analysis ----
vfdb_file <- config$VFDB_SUMMARY
if (file.exists(vfdb_file)) {
  vfdb <- fread(vfdb_file, sep = "\t", header = TRUE, data.table = FALSE)
  n_vfdb_genomes <- nrow(vfdb)
  n_vf_genes <- sum(grepl("VF", colnames(vfdb))) + sum(grepl("virulence", colnames(vfdb)))
  cat("VFDB summary genomes:", n_vfdb_genomes, "\n")
} else {
  vfdb <- NULL
  n_vfdb_genomes <- NA
  n_vf_genes <- NA
}

# ---- 6. VF binary matrix ----
if (file.exists(config$VF_BINARY)) {
  vf_bin <- fread(config$VF_BINARY, sep = "\t", header = TRUE, data.table = FALSE)
  if ("#FILE" %in% colnames(vf_bin)) n_vf_bin <- nrow(vf_bin) else n_vf_bin <- ncol(vf_bin)
  cat("VF binary matrix entries:", n_vf_bin, "\n")
}

# ---- 7. CARD burden ----
card_file <- config$st_card_burden("ST69")
if (file.exists(card_file)) {
  card <- fread(card_file, sep = "\t", header = TRUE, data.table = FALSE)
  n_card_genomes <- nrow(card)
  cat("CARD ARG summary genomes:", n_card_genomes, "\n")
}

# ---- 8. VFDB shell cluster master table ----
vfdb_master <- file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv")
if (file.exists(vfdb_master)) {
  master <- fread(vfdb_master, sep = ",", header = TRUE, data.table = FALSE)
  n_master_genomes <- nrow(master)
  n_clusters <- master %>% filter(!is.na(shell_cluster)) %>% pull(shell_cluster) %>% n_distinct()
  cat("Master table genomes:", n_master_genomes, "\n")
  cat("Shell clusters:", n_clusters, "\n")
}

# ---- 9. RGP analysis ----
rgp_file <- config$RGP_REGIONS
if (file.exists(rgp_file)) {
  rgp <- fread(rgp_file, sep = "\t", header = TRUE, data.table = FALSE, nThread = 1, nrows = 10)
  cat("RGP regions file exists:", rgp_file, "\n")
}

# ---- 10. Summary table ----
summary <- tribble(
  ~step, ~count,
  "PPanGGOLiN gene families", n_panx_genes,
  "PPanGGOLiN genomes (contigs)", n_panx_genomes,
  "Metadata entries (isolates)", n_meta,
  "Metadata with year", n_with_year,
  "VFDB summary genomes", n_vfdb_genomes,
  "Master table genomes (VFDB + cluster)", if (exists("n_master_genomes")) n_master_genomes else NA,
  "Shell clusters", if (exists("n_clusters")) n_clusters else NA,
  "CARD summary genomes", if (exists("n_card_genomes")) n_card_genomes else NA
)

cat("\n=== Summary Statistics ===\n")
print(summary)

write_xlsx(list(summary = summary), file.path(OUT, "summary_statistics.xlsx"))
cat("\nDone. Output in:", OUT, "\n")
