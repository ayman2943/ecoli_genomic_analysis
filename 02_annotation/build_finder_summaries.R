#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(stringr)
  library(janitor)
})

options(stringsAsFactors = FALSE)

# =============================================================================
# CONFIG
# =============================================================================
# BASE_DIR: directory containing the {ST}/ analysis_results folders
#           (i.e. the folder that holds ST10/, ST69/, ... with virulence/ and
#           resfinder/ subfolders). For this pipeline that is the
#           analysis_results/ folder produced by run_vf_resfinder.sh.
# OUT_DIR : where the summary matrices are written. Defaults to the
#           finder_result/ folder expected by config.R. Override with:
#           Rscript build_finder_summaries.R <base_dir> <out_dir>
BASE_DIR <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[1] else "."
VF_DIR_NAME <- "virulence"
RF_DIR_NAME <- "resfinder"

OUT_DIR <- if (length(commandArgs(trailingOnly = TRUE)) >= 2) commandArgs(trailingOnly = TRUE)[2] else "finder_result"

# =============================================================================
# FIND ST DIRECTORIES
# =============================================================================
st_dirs <- list.dirs(BASE_DIR, full.names = TRUE, recursive = FALSE) |>
  keep(~ dir.exists(file.path(.x, VF_DIR_NAME)) &&
         dir.exists(file.path(.x, RF_DIR_NAME))) |>
  sort()

if (length(st_dirs) == 0) {
  stop("No ST directories found with virulence/ and resfinder/ folders.")
}

# =============================================================================
# OUTPUT DIRS
# =============================================================================
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(OUT_DIR, "virulencefinder_summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "resfinder_summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "logs"), recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================
clean_genome_id <- function(x) {
  x |>
    basename() |>
    str_remove("\\.(txt|tsv|csv)$") |>
    str_remove("_(VF|RF)$") |>
    str_replace_all("[[:space:]]+", "_")
}

sanitize_gene <- function(x) {
  x |>
    as.character() |>
    str_trim() |>
    str_replace_all("\\s+", "_") |>
    str_replace_all("[^A-Za-z0-9_.()\\-]", "_")
}

safe_fread <- function(path) {
  tryCatch(
    fread(
      file = path,
      sep = "\t",
      header = TRUE,
      fill = TRUE,
      quote = "",
      data.table = FALSE,
      showProgress = FALSE
    ),
    error = function(e) NULL
  )
}

choose_col <- function(df, candidates) {
  nm <- names(df)
  idx <- which(tolower(nm) %in% tolower(candidates))
  if (length(idx) == 0) return(NA_character_)
  nm[idx[1]]
}

make_binary_matrix <- function(df, id_cols = c("st", "genome"), gene_col = "gene") {
  if (nrow(df) == 0) return(tibble(st = character(), genome = character()))

  df |>
    distinct(across(all_of(id_cols)), .data[[gene_col]]) |>
    mutate(value = 1L) |>
    pivot_wider(
      id_cols = all_of(id_cols),
      names_from = .data[[gene_col]],
      values_from = value,
      values_fill = 0
    ) |>
    arrange(across(all_of(id_cols)))
}

# =============================================================================
# EMPTY STRUCTURES
# =============================================================================
empty_vf <- function() {
  tibble(
    st = character(),
    genome = character(),
    gene = character(),
    database = character(),
    identity = numeric(),
    contig = character(),
    contig_position = character(),
    protein_function = character(),
    accession = character(),
    source_file = character()
  )
}

empty_rf <- function() {
  tibble(
    st = character(),
    genome = character(),
    gene = character(),
    identity = numeric(),
    alignment_length = numeric(),
    gene_length = numeric(),
    coverage = numeric(),
    position_in_reference = character(),
    contig = character(),
    position_in_contig = character(),
    phenotype = character(),
    accession = character(),
    source_file = character()
  )
}

# =============================================================================
# VIRULENCEFINDER PARSER
# =============================================================================
parse_vf_file <- function(path, st) {

  genome <- clean_genome_id(path)
  df <- safe_fread(path)

  if (is.null(df) || nrow(df) == 0) return(empty_vf())

  names(df) <- janitor::make_clean_names(names(df))

  gene_col <- choose_col(df, c("virulence_factor", "gene", "gene_name"))
  db_col   <- choose_col(df, c("database"))
  id_col   <- choose_col(df, c("identity"))
  contig_col <- choose_col(df, c("contig"))
  pos_col    <- choose_col(df, c("position_in_contig"))
  func_col   <- choose_col(df, c("protein_function"))
  acc_col    <- choose_col(df, c("accession_number", "accession"))

  if (is.na(gene_col)) return(empty_vf())

  tibble(
    st = st,
    genome = genome,
    gene = sanitize_gene(df[[gene_col]]),
    database = if (!is.na(db_col)) as.character(df[[db_col]]) else NA_character_,
    identity = if (!is.na(id_col)) suppressWarnings(as.numeric(df[[id_col]])) else NA_real_,
    contig = if (!is.na(contig_col)) as.character(df[[contig_col]]) else NA_character_,
    contig_position = if (!is.na(pos_col)) as.character(df[[pos_col]]) else NA_character_,
    protein_function = if (!is.na(func_col)) as.character(df[[func_col]]) else NA_character_,
    accession = if (!is.na(acc_col)) as.character(df[[acc_col]]) else NA_character_,
    source_file = path
  ) |>
    filter(!is.na(gene), gene != "") |>
    distinct()
}

# =============================================================================
# RESFINDER PARSER (FIXED - NO readLines HEADER LOGIC)
# =============================================================================
parse_rf_file <- function(path, st) {

  genome <- clean_genome_id(path)

  df <- safe_fread(path)

  if (is.null(df) || nrow(df) == 0) return(empty_rf())

  names(df) <- janitor::make_clean_names(names(df))

  gene_col  <- choose_col(df, c("resistance_gene", "gene"))
  id_col    <- choose_col(df, c("identity"))
  cov_col   <- choose_col(df, c("coverage"))
  contig_col <- choose_col(df, c("contig"))
  pcontig_col <- choose_col(df, c("position_in_contig"))
  pref_col  <- choose_col(df, c("position_in_reference"))
  pheno_col <- choose_col(df, c("phenotype"))
  acc_col   <- choose_col(df, c("accession_no", "accession"))
  alg_col   <- choose_col(df, c("alignment_length_gene_length"))

  if (is.na(gene_col)) return(empty_rf())

  alg_raw <- if (!is.na(alg_col)) as.character(df[[alg_col]]) else NA_character_

  alignment_length <- suppressWarnings(str_extract(alg_raw, "^\\s*\\d+"))
  alignment_length <- as.numeric(alignment_length)

  gene_length <- suppressWarnings(str_extract(alg_raw, "(?<=/)\\s*\\d+"))
  gene_length <- as.numeric(gene_length)

  tibble(
    st = st,
    genome = genome,
    gene = sanitize_gene(df[[gene_col]]),

    identity = if (!is.na(id_col)) suppressWarnings(as.numeric(df[[id_col]])) else NA_real_,
    coverage = if (!is.na(cov_col)) suppressWarnings(as.numeric(df[[cov_col]])) else NA_real_,

    alignment_length = alignment_length,
    gene_length = gene_length,

    position_in_reference = if (!is.na(pref_col)) as.character(df[[pref_col]]) else NA_character_,
    contig = if (!is.na(contig_col)) as.character(df[[contig_col]]) else NA_character_,
    position_in_contig = if (!is.na(pcontig_col)) as.character(df[[pcontig_col]]) else NA_character_,

    phenotype = if (!is.na(pheno_col)) as.character(df[[pheno_col]]) else NA_character_,
    accession = if (!is.na(acc_col)) as.character(df[[acc_col]]) else NA_character_,

    source_file = path
  ) |>
    filter(!is.na(gene), gene != "") |>
    distinct()
}

# =============================================================================
# MAIN LOOP
# =============================================================================
vf_all <- list()
rf_all <- list()

cat("Found STs:\n")
print(basename(st_dirs))
cat("\n")

for (st_dir in st_dirs) {

  st <- basename(st_dir)

  vf_dir <- file.path(st_dir, VF_DIR_NAME)
  rf_dir <- file.path(st_dir, RF_DIR_NAME)

  vf_files <- list.files(vf_dir, pattern = "\\.(tsv|txt)$", full.names = TRUE, recursive = TRUE)
  rf_files <- list.files(rf_dir, pattern = "\\.(tsv|txt)$", full.names = TRUE, recursive = TRUE)

  cat("Processing:", st, "\n")
  cat(" VF:", length(vf_files), " RF:", length(rf_files), "\n")

  vf_long <- purrr::map_dfr(vf_files, parse_vf_file, st = st)
  rf_long <- purrr::map_dfr(rf_files, parse_rf_file, st = st)

  vf_all[[st]] <- vf_long
  rf_all[[st]] <- rf_long

  cat(" Done:", st, "\n\n")
}

vf_all <- bind_rows(vf_all)
rf_all <- bind_rows(rf_all)

# =============================================================================
# OUTPUTS
# =============================================================================
write_tsv(vf_all, file.path(OUT_DIR, "virulencefinder_summary", "virulencefinder_long.tsv"))
write_tsv(rf_all, file.path(OUT_DIR, "resfinder_summary", "resfinder_long.tsv"))

write_tsv(
  make_binary_matrix(vf_all),
  file.path(OUT_DIR, "virulencefinder_summary", "virulencefinder_binary_matrix.tsv")
)

write_tsv(
  make_binary_matrix(rf_all),
  file.path(OUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv")
)

vf_burden <- vf_all |>
  distinct(st, genome, gene) |>
  count(st, genome, name = "virulencefinder_burden")

rf_burden <- rf_all |>
  distinct(st, genome, gene) |>
  count(st, genome, name = "resfinder_burden")

write_tsv(vf_burden, file.path(OUT_DIR, "virulencefinder_summary", "virulencefinder_burden.tsv"))
write_tsv(rf_burden, file.path(OUT_DIR, "resfinder_summary", "resfinder_burden.tsv"))

vf_gene_freq <- vf_all |>
  distinct(st, genome, gene) |>
  count(st, gene, sort = TRUE, name = "n_genomes")

rf_gene_freq <- rf_all |>
  distinct(st, genome, gene) |>
  count(st, gene, sort = TRUE, name = "n_genomes")

write_tsv(vf_gene_freq, file.path(OUT_DIR, "virulencefinder_summary", "virulencefinder_gene_frequency.tsv"))
write_tsv(rf_gene_freq, file.path(OUT_DIR, "resfinder_summary", "resfinder_gene_frequency.tsv"))

# =============================================================================
# QC
# =============================================================================
qc <- tibble(
  dataset = c("VirulenceFinder", "ResFinder"),
  n_rows = c(nrow(vf_all), nrow(rf_all)),
  n_genomes = c(n_distinct(vf_all$genome), n_distinct(rf_all$genome)),
  n_genes = c(n_distinct(vf_all$gene), n_distinct(rf_all$gene))
)

write_tsv(qc, file.path(OUT_DIR, "logs", "summary_qc.tsv"))

cat("\nDone.\n")
print(qc)