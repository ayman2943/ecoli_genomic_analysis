#!/usr/bin/env Rscript
# ==============================================================================
# Virulence / Plasmid_VFDB.R
# ==============================================================================
# Plasmid-associated virulence gene (VFDB) annotation layer.
#
# PROVENANCE: extracted from the original `10d_plasmid_context.R`
# (sections 1-3 only: genome-ID cleanup, master-table loading, and plasmid
# VFDB summary loading/binary classification). This is the general-purpose
# data layer that classifies each "interest" virulence gene as present on a
# plasmid vs. chromosome-only vs. absent, for every ST69 genome.
#
# The remainder of the original `10d_plasmid_context.R` (finding the top-10
# increasing VF genes in Cluster_3, the RGP co-occurrence analysis, and the
# 3-panel combined figure) is FIGURE-SPECIFIC — it reproduces submitted
# Figure 8 exactly, and has been kept intact and unmodified as
# `Analysis/Figure8.R`. That script re-loads the same plasmid summary file
# itself, so it does not depend on this script having been run first; the
# two are independent entry points into the same source data, matching how
# the original single script was structured.
#
# ANALYSIS IS UNCHANGED — this is a verbatim relocation of the original
# code (sections 1-3), not a rewrite of the underlying logic.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(data.table)
})
source("config.R")

OUT <- file.path(config$OUTPUT_DIR, config$TARGET_ST, "reviewer_rgp_context")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

################################################################################
# 1. Helper: clean genome IDs
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
# 2. Load data
################################################################################
cat("Loading data...\n")

# Master table (has cluster + VF binary)
master <- fread(
  file.path(config$st_out_vfdb(config$TARGET_ST), "04_master_shell_cluster_metadata_VFDB_table.csv"),
  sep=",", header=TRUE, data.table=FALSE) %>%
  mutate(gid_clean = clean_gid(.data[["#FILE"]]), year=as.integer(year))

# Interest genes
interest <- c("iucA","iucB","iucC","iucD","iutA","sat",
  "papA","papB","papC","papD","papE","papF","papG","papH","papI","papJ","papK","papX",
  "kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT","kpsU")
gene_sys <- c(
  iucA="iuc",iucB="iuc",iucC="iuc",iucD="iuc",iutA="iuc",
  sat="sat",
  papA="pap",papB="pap",papC="pap",papD="pap",papE="pap",papF="pap",
  papG="pap",papH="pap",papI="pap",papJ="pap",papK="pap",papX="pap",
  kpsC="kps",kpsD="kps",kpsE="kps",kpsF="kps",kpsM="kps",kpsS="kps",kpsT="kps",kpsU="kps")

# Parse whole-genome VF binary from master table
to_bin <- function(x) as.integer(!is.na(x) & x != ".")
master_bin <- master %>%
  mutate(across(all_of(intersect(interest, colnames(master))), to_bin, .names="{.col}_wg"))

################################################################################
# 3. Load plasmid VFDB summary
################################################################################
cat("Loading plasmid data...\n")
pl_file <- file.path(config$BASE_DIR,"..","Plasmid","plasmid_vfdb_summary",paste0(config$TARGET_ST, "_plasmid_summary.tsv"))
if (!file.exists(pl_file)) {
  stop("Plasmid VFDB summary not found at ", pl_file,
       " (requires the external Plasmid dataset). Skipping plasmid-context analysis.")
}
pl_vf <- fread(pl_file,
  sep="\t", header=TRUE, data.table=FALSE, na.strings=c("",".","-","NA"),
  colClasses="character", check.names=FALSE) %>%
  mutate(gid_clean = clean_gid(.data[["#FILE"]]))

# Parse plasmid binary
pl_genes <- intersect(interest, setdiff(colnames(pl_vf), c("#FILE","NUM_FOUND","gid_clean")))
pl_bin <- pl_vf %>%
  select(gid_clean, all_of(pl_genes)) %>%
  mutate(across(all_of(pl_genes), to_bin, .names="{.col}_pl"))

cat("\nDone. Loaded master table (", nrow(master), " genomes) and plasmid VFDB summary (",
    nrow(pl_vf), " genomes) for ", length(interest), " interest genes.\n", sep="")
cat("Outputs of this script are consumed directly by Analysis/Figure8.R, which\n")
cat("re-derives the same objects (master, pl_vf, pl_bin) and carries the\n")
cat("figure-specific top-10/RGP/plotting logic forward.\n")
