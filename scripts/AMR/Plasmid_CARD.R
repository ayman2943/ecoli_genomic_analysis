#!/usr/bin/env Rscript
# ==============================================================================
# AMR / Plasmid_CARD.R
# ==============================================================================
# Plasmid-associated AMR gene (CARD) annotation layer.
#
# PROVENANCE: this script did not exist in the original repository -- no
# CARD/AMR equivalent of the plasmid-VFDB analysis (10d_plasmid_context.R /
# Virulence/Plasmid_VFDB.R) had ever been run. It is written here, now that
# scripts/Plasmid_assembly/run_mobsuite.sh and run_abricate_plasmids.sh
# produce the plasmid CARD summary those scripts feed. It mirrors
# Virulence/Plasmid_VFDB.R's structure and logic exactly (same clean_gid()
# helper, same to_bin() presence rule, same chromosomal/plasmid/absent
# classification), applied to CARD instead of VFDB.
#
# ONE DELIBERATE DIFFERENCE from Plasmid_VFDB.R: that script restricts
# classification to a curated "interest" gene list (iuc/pap/sat/kps) that
# came directly from the original 10d_plasmid_context.R -- reproducing an
# existing, manuscript-used analysis choice. No equivalent curated AMR gene
# list exists anywhere in this repository, and inventing one here would not
# be "reproducing" anything -- it would be new analysis. So this script
# instead classifies every CARD gene that appears in BOTH the whole-genome
# master table and the plasmid CARD summary (their column intersection),
# which is the general-purpose, non-curated equivalent.
#
# Requires scripts/Plasmid_assembly/run_mobsuite.sh and
# run_abricate_plasmids.sh to have been run first (produces
# ../Plasmid/plasmid_card_summary/{ST}_plasmid_summary.tsv).
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(data.table)
})
source("config.R")

OUT <- file.path(config$OUTPUT_DIR, config$TARGET_ST, "reviewer_rgp_context")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

################################################################################
# 1. Helper: clean genome IDs (identical to Virulence/Plasmid_VFDB.R)
################################################################################
clean_gid <- function(raw) {
  gid <- basename(raw)
  gid <- sub("^Escherichia_coli_","",gid)
  gid <- sub("_combined_plasmids_vfdb\\.tsv$","",gid)
  gid <- sub("_combined_plasmids_card\\.tsv$","",gid)
  gid <- sub("_vfdb\\.tsv$","",gid)
  gid <- sub("_card\\.tsv$","",gid)
  gid <- sub("\\.tsv$","",gid)
  gid
}

################################################################################
# 2. Load whole-genome master CARD table
################################################################################
cat("Loading data...\n")

card_file <- config$st_card_burden()
if (!file.exists(card_file)) {
  stop("Whole-genome CARD summary not found at ", card_file,
       " -- run scripts/Annotation/run_abricate.sh first.")
}
master <- fread(card_file, sep="\t", header=TRUE, data.table=FALSE,
    na.strings=c("",".","-","NA"), colClasses="character", check.names=FALSE) %>%
  mutate(gid_clean = clean_gid(.data[["#FILE"]]))

# Every column that isn't abricate's own bookkeeping is a CARD gene.
meta_cols <- c("#FILE", "NUM_FOUND", "gid_clean")
card_genes <- setdiff(colnames(master), meta_cols)
cat("Whole-genome CARD genes:", length(card_genes), "across", nrow(master), "genomes\n")

# Parse whole-genome CARD binary from master table
to_bin <- function(x) as.integer(!is.na(x) & x != ".")
master_bin <- master %>%
  mutate(across(all_of(card_genes), to_bin, .names="{.col}_wg"))

################################################################################
# 3. Load plasmid CARD summary
################################################################################
cat("Loading plasmid data...\n")
pl_file <- file.path(config$BASE_DIR,"..","Plasmid","plasmid_card_summary",paste0(config$TARGET_ST, "_plasmid_summary.tsv"))
if (!file.exists(pl_file)) {
  stop("Plasmid CARD summary not found at ", pl_file,
       " -- run scripts/Plasmid_assembly/run_mobsuite.sh then ",
       "run_abricate_plasmids.sh first (requires the external Plasmid dataset).")
}
pl_card <- fread(pl_file,
  sep="\t", header=TRUE, data.table=FALSE, na.strings=c("",".","-","NA"),
  colClasses="character", check.names=FALSE) %>%
  mutate(gid_clean = clean_gid(.data[["#FILE"]]))

# Classify every gene common to both the whole-genome and plasmid tables
# (no curated "interest" subset -- see header).
pl_genes <- intersect(card_genes, setdiff(colnames(pl_card), c("#FILE","NUM_FOUND","gid_clean")))
cat("CARD genes observed on plasmids:", length(pl_genes), "of", length(card_genes), "\n")
pl_bin <- pl_card %>%
  select(gid_clean, all_of(pl_genes)) %>%
  mutate(across(all_of(pl_genes), to_bin, .names="{.col}_pl"))

################################################################################
# 4. Per-gene chromosomal vs. plasmid vs. absent classification (all genomes)
################################################################################
cat("\n=== Per-gene genomic location (all", config$TARGET_ST, ") ===\n")

merged <- master_bin %>%
  select(gid_clean, ends_with("_wg")) %>%
  left_join(pl_bin, by="gid_clean", suffix=c("", ".pl_only"))

gene_loc <- map_dfr(pl_genes, function(g) {
  wg <- paste0(g, "_wg")
  pl <- paste0(g, "_pl")
  if (!pl %in% colnames(merged)) merged[[pl]] <<- 0L

  merged %>%
    summarise(
      n_total = n(),
      n_wg = sum(.data[[wg]] == 1, na.rm=TRUE),
      n_pl = sum(.data[[pl]] == 1, na.rm=TRUE)
    ) %>%
    mutate(
      gene = g,
      pct_wg = round(n_wg/n_total*100, 1),
      pct_pl = round(n_pl/n_total*100, 1),
      pct_chr = round((n_wg - n_pl)/n_total*100, 1),
      pct_abs = round(100 - pct_wg, 1)
    )
})

cat("Gene-level localization (all", config$TARGET_ST, "):\n")
print(gene_loc %>% select(gene, pct_wg, pct_chr, pct_pl, pct_abs), row.names=FALSE)

write_csv(gene_loc, file.path(OUT, "plasmid_card_gene_location.csv"))
cat("\nSaved:", file.path(OUT, "plasmid_card_gene_location.csv"), "\n")
cat("Done.\n")
