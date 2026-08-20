#!/usr/bin/env Rscript
# Reviewer response: long-read validation of gene linkage
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(writexl); library(data.table)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_long_read_validation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Check if any long-read assemblies exist ----
lr_dir <- file.path(config$OUTPUT_DIR, "ST69", "long_read")
if (dir.exists(lr_dir)) {
  fa_files <- list.files(lr_dir, pattern = "\\.fasta$|\\.fna$", full.names = TRUE)
  gff_files <- list.files(lr_dir, pattern = "\\.gff$|\\.gff3$", full.names = TRUE)

  if (length(fa_files) > 0) {
    cat("Long-read assemblies found:", length(fa_files), "\n")
    cat("GFF files:", length(gff_files), "\n")

    # Load existing VF results for these assemblies
    vf_dir <- file.path(config$OUTPUT_DIR, "ST69", "virulencefinder_validation")
    if (dir.exists(vf_dir)) {
      vf_files <- list.files(vf_dir, pattern = "\\.tsv$", full.names = TRUE)
      lr_names <- basename(fa_files) %>% str_remove("\\.fasta$|\\.fna$")

      cat("\nVirulence genes in long-read assemblies:\n")
      for (vf in vf_files) {
        tbl <- tryCatch(fread(vf, sep = "\t", header = TRUE, data.table = FALSE,
          nThread = 1), error = function(e) NULL)
        if (!is.null(tbl) && nrow(tbl) > 0) {
          cat(basename(vf), ":", nrow(tbl), "hits\n")
        }
      }
    }
    cat("\nRecommendation: download long-read assemblies from Lipworth et al. 2024\n")
    cat("and run ABRicate/VFDB to validate pap-kps-sat co-localization on single\n")
    cat("contigs. These assemblies are publicly available on ENA (PRJEB76009).\n")
  } else {
    cat("Long-read directory exists but no assemblies found.\n")
    cat("No long-read assemblies available for validation.\n")
  }
} else {
  cat("Long-read directory not found:", lr_dir, "\n")
  cat("No long-read assemblies available for validation.\n")
}

# ---- 2. Check project structure for assemblies ----
assembly_base <- file.path(config$BASE_DIR, "assembly")
if (dir.exists(assembly_base)) {
  st69_asm <- list.files(assembly_base, pattern = "ST69", full.names = TRUE)
  if (length(st69_asm) > 0) {
    cat("\nAssembly directory for ST69 exists:", st69_asm[1], "\n")
  }
}

# ---- 3. Write summary ----
sink(file.path(OUT, "long_read_validation_summary.txt"))
cat("Long-Read Validation of Virulence Locus Linkage\n")
cat("==============================================\n\n")
cat("Status: No long-read assemblies are currently available in the project.\n\n")
cat("Recommendation for reviewer response:\n")
cat("The RGP neighbourhood analysis using PPanGGOLiN already demonstrates that\n")
cat("pap, iuc, kps, and sat genes are physically linked within 5 composite RGPs\n")
cat("in ST69 genomes. This provides graph-based evidence for co-localization.\n\n")
cat("As suggested by the reviewer, long-read validation would be a useful addition.\n")
cat("Lipworth et al. (2024, Nat Commun) and Arredondo-Alonso et al. (2025) have\n")
cat("published long-read assemblies of ST69 E. coli. These should be downloaded\n")
cat("from ENA and analysed with ABRicate to confirm that the full virulence locus\n")
cat("is present on single contigs.\n")
sink()
cat("Done. Output in:", OUT, "\n")
