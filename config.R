#!/usr/bin/env Rscript
# config.R - central configuration for the E. coli genomic analysis pipeline.
#
# All paths are derived from environment variables (with sensible defaults) so
# the pipeline can run on any machine. Set these in config/pipeline_config.sh or
# export them before launching the analysis:
#
#   ECOLI_BASE_DIR       working directory (default: current working directory)
#   ECOLI_PANGENOME_DIR  external pangenome directory used for non-ST69 lineages
#                        (default: empty - falls back to <BASE>/pangenome_output)
#   TARGET_ST            lineage to analyse (default: ST69)
#
# NOTE: this file replaces the machine-specific 00_config.R used during
# development; field names are preserved for compatibility with every analysis
# script in scripts/.

config <- list()

# ---- Base directories (must be first) ----
config$BASE_DIR <- Sys.getenv("ECOLI_BASE_DIR", unset = getwd())
config$PANGENOME_DIR <- file.path(config$BASE_DIR, "pangenome_output")
config$PPANGGOLIN_DIR <- file.path(config$PANGENOME_DIR, "ppanggolin_output")
config$INPUT_DIR <- file.path(config$BASE_DIR, "finder_result")
config$CARD_VFDB_DIR <- file.path(config$BASE_DIR, "card_vfdb_result")
config$OUTPUT_DIR <- file.path(config$BASE_DIR, "output")

# ---- ST parameter ----
config$TARGET_ST <- toupper(Sys.getenv("TARGET_ST", unset = "ST69"))
config$AVAILABLE_STS <- c("ST10", "ST69", "ST73", "ST95", "ST131")
config$MIN_YEAR_N <- 20
config$MAX_YEAR <- 2025
config$CLUSTER_MAX_K <- 10
config$N_CLUSTERS <- 4
config$EARLY_YEARS <- c(2016, 2017, 2018)
config$LATE_YEARS <- c(2022, 2023, 2024, 2025)

# ---- ST-specific path generators ----
config$st_tree <- function(st = config$TARGET_ST) {
  file.path(config$BASE_DIR, paste0(st, "_bootstrap.treefile"))
}
config$st_metadata <- function(st = config$TARGET_ST) {
  f <- file.path(config$BASE_DIR, "metadata_matched", paste0("matched_", st, ".xlsx"))
  if (file.exists(f)) return(f)
  file.path(config$BASE_DIR, paste0(st, "_metadata.xlsx"))
}
config$st_pangenome <- function(st = config$TARGET_ST) {
  ext <- Sys.getenv("ECOLI_PANGENOME_DIR", unset = "")
  if (nzchar(ext) && st != "ST69") {
    return(file.path(ext, "pangenome_output"))
  }
  file.path(config$PANGENOME_DIR, "ppanggolin_output")
}
config$st_vfdb_summary <- function(st = config$TARGET_ST) {
  file.path(config$CARD_VFDB_DIR, "vfdb_summary", paste0(st, "_vfdb_summary.tsv"))
}
config$st_card_burden <- function(st = config$TARGET_ST) {
  file.path(config$CARD_VFDB_DIR, "card_summary", paste0(st, ".tsv"))
}
config$st_vf_summary <- function(st = config$TARGET_ST) {
  file.path(config$INPUT_DIR, st, "virulencefinder_summary", "virulencefinder_burden.tsv")
}
config$st_resfinder_summary <- function(st = config$TARGET_ST) {
  file.path(config$INPUT_DIR, st, "resfinder_summary", "resfinder_burden.tsv")
}
config$st_out_vf <- function(st = config$TARGET_ST) {
  d <- file.path(config$OUTPUT_DIR, st, "virulencefinder_validation")
  dir.create(d, showWarnings = FALSE, recursive = TRUE); d
}
config$st_out_vfdb <- function(st = config$TARGET_ST) {
  d <- file.path(config$OUTPUT_DIR, st, "vfdb_analysis")
  dir.create(d, showWarnings = FALSE, recursive = TRUE); d
}

# ---- ST-specific input files ----
config$GENE_PA_AB <- file.path(config$st_pangenome(), "gene_presence_absence.Rtab")
config$SHELL_FILE <- file.path(config$st_pangenome(), "partitions", "shell.txt")
config$VF_BINARY <- file.path(config$INPUT_DIR, "virulencefinder_summary", "virulencefinder_binary_matrix.tsv")
config$VF_BURDEN <- file.path(config$INPUT_DIR, "virulencefinder_summary", "virulencefinder_burden.tsv")
config$RESF_BURDEN <- file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_burden.tsv")
config$MODULES_RGP <- file.path(config$PPANGGOLIN_DIR, "modules_RGP_lists.tsv")
config$RGP_REGIONS <- file.path(config$PANGENOME_DIR, "regions_of_genomic_plasticity.tsv")

# ---- Resolve ST-specific file paths and check existence ----
config$TREE_FILE <- config$st_tree()
config$METADATA_FILE <- config$st_metadata()
config$VFDB_SUMMARY <- config$st_vfdb_summary()

# ---- Output directories for this ST ----
config$VF_DIR <- config$st_out_vf()
config$VFDB_DIR <- config$st_out_vfdb()

cat(sprintf("=== CONFIG: TARGET_ST = %s ===\n", config$TARGET_ST))
cat(sprintf("BASE_DIR  = %s\n", config$BASE_DIR))
cat(sprintf("VF_DIR    = %s\n", config$VF_DIR))
cat(sprintf("VFDB_DIR  = %s\n", config$VFDB_DIR))
