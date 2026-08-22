#!/usr/bin/env Rscript
#
# 00_metadata/filter_metadata.R
#
# Filter a raw EnteroBase metadata export into per-ST input files used by the
# rest of the pipeline:
#   metadata/{ST}_filtered.xlsx          <- rows for the target ST (genomes)
#   metadata_matched/matched_{ST}.xlsx   <- metadata matched to downloaded genomes
#
# Usage:
#   Rscript 00_metadata/filter_metadata.R RAW_ENTERO_EXPORT.xlsx
#   TARGET_ST=ST69 Rscript 00_metadata/filter_metadata.R export.xlsx
#
# The raw export must contain the standard EnteroBase columns, at minimum
# 'Uberstrain' (or 'Name') and a column identifying the ST.

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript filter_metadata.R <raw_entero_export.xlsx>")
}
raw_file <- args[1]
TARGET_ST <- toupper(Sys.getenv("TARGET_ST", unset = "ST69"))

# ---- Read raw export ----
if (!file.exists(raw_file)) stop("Raw export not found: ", raw_file)
raw <- read_excel(raw_file)

# ---- Identify name + ST columns ----
name_col <- intersect(c("Uberstrain", "Name"), colnames(raw))[1]
st_col <- intersect(c("ST", "Sequence Type", "MLST"), colnames(raw))[1]
if (is.na(name_col)) stop("No name column (Uberstrain/Name) found")
if (is.na(st_col)) stop("No ST column found")

# ---- Filter to the target ST ----
df <- raw %>%
  rename(genome_id = all_of(name_col)) %>%
  filter(!is.na(genome_id))

st_match <- str_remove(TARGET_ST, "^ST")
df <- df %>%
  filter(!is.na(.data[[st_col]]) &
         as.character(.data[[st_col]]) == st_match)

cat(sprintf("TARGET_ST=%s : %d genomes retained from %d rows\n",
            TARGET_ST, nrow(df), nrow(raw)))

dir.create("metadata", showWarnings = FALSE)
dir.create("metadata_matched", showWarnings = FALSE)

# ---- Write filtered input ----
filtered_file <- file.path("metadata", paste0(TARGET_ST, "_filtered.xlsx"))
writexl::write_xlsx(df, filtered_file)
cat("Wrote:", filtered_file, "\n")

# ---- Write matched metadata (used by config$st_metadata) ----
matched_file <- file.path("metadata_matched", paste0("matched_", TARGET_ST, ".xlsx"))
writexl::write_xlsx(df, matched_file)
cat("Wrote:", matched_file, "\n")
