#!/usr/bin/env Rscript
################################################################################
# PLASMID ANALYSIS PIPELINE - Main Script
# Plasmid-associated ARG and VF Analysis
################################################################################

# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(pheatmap)
  library(scales)
  library(writexl)
})

# Configuration
theme_set(theme_minimal(base_size = 13))

STs <- c("ST10", "ST131", "ST69", "ST73", "ST95")
OUT_DIR <- "outputs/plasmid_outputs"

# Create output directories
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Source utility functions
source("scripts/utils/common_utils.R")

cat("\n================================================================================\n")
cat("PLASMID-ASSOCIATED GENE ANALYSIS\n")
cat("================================================================================\n\n")

################################################################################
# SECTION 1: DATA LOADING
################################################################################

cat("SECTION 1: Loading plasmid datasets\n")
cat("================================================================================\n\n")

source("scripts/modules/plasmid/01_plasmid_loading.R")
plasmid_data <- load_plasmid_data(STs)

card <- plasmid_data$card
vfdb <- plasmid_data$vfdb

################################################################################
# SECTION 2: ARG BURDEN ANALYSIS
################################################################################

cat("\nSECTION 2: Plasmid ARG burden analysis\n")
cat("================================================================================\n\n")

source("scripts/modules/plasmid/02_arg_burden.R")
arg_results <- analyze_plasmid_arg_burden(card, STs, OUT_DIR)

################################################################################
# SECTION 3: HIGH-RISK ARG DETECTION
################################################################################

cat("\nSECTION 3: High-risk ARG detection (ESBL/Carbapenemase/MCR)\n")
cat("================================================================================\n\n")

source("scripts/modules/plasmid/03_highrisk_detection.R")
highrisk_results <- analyze_highrisk_args(card, STs, OUT_DIR)

################################################################################
# SECTION 4: VF PLASMID ANALYSIS
################################################################################

cat("\nSECTION 4: Plasmid-borne virulence factors\n")
cat("================================================================================\n\n")

source("scripts/modules/plasmid/04_vf_plasmid.R")
vf_results <- analyze_plasmid_vf(vfdb, STs, OUT_DIR)

################################################################################
# SECTION 5: ARG-VF CORRELATION
################################################################################

cat("\nSECTION 5: ARG-VF co-occurrence on plasmids\n")
cat("================================================================================\n\n")

source("scripts/modules/plasmid/05_correlation_analysis.R")
correlation_results <- analyze_arg_vf_correlation(card, vfdb, OUT_DIR)

################################################################################
# SECTION 6: EXPORT SUMMARY
################################################################################

cat("\nSECTION 6: Exporting summary tables\n")
cat("================================================================================\n\n")

# Combine all results
summary_data <- list(
  "ARG_Burden" = arg_results$summary,
  "MDR_Plasmids" = arg_results$mdr_summary,
  "HighRisk_ARGs" = highrisk_results$summary,
  "VF_Burden" = vf_results$summary,
  "ARG_VF_Correlation" = correlation_results$correlation_test,
  "Temporal_ARG_Trends" = arg_results$temporal_trends,
  "Temporal_VF_Trends" = vf_results$temporal_trends
)

write_xlsx(summary_data, file.path(OUT_DIR, "plasmid_summary.xlsx"))
cat(sprintf("  ✓ plasmid_summary.xlsx\n"))

################################################################################
# COMPLETION SUMMARY
################################################################################

cat("\n================================================================================\n")
cat("PLASMID ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat(sprintf("Output directory: %s\n\n", OUT_DIR))
cat(sprintf("Total CARD genomes: %d\n", nrow(card)))
cat(sprintf("Total VFDB genomes: %d\n", nrow(vfdb)))
cat(sprintf("Sequence types: %s\n", paste(STs, collapse = ", ")))

cat("\nKey Findings:\n")
if (exists("correlation_results") && !is.null(correlation_results$correlation_test)) {
  cat(sprintf("  ARG-VF Spearman rho: %.3f (p = %.2e)\n", 
              correlation_results$correlation_test$rho,
              correlation_results$correlation_test$p_value))
}

cat("\nVisualizations generated:\n")
cat("  P1_plasmid_ARG_trend.png\n")
cat("  P2_plasmid_MDR_rate.png\n")
cat("  P3_highrisk_plasmid.png\n")
cat("  P4_plasmid_VF_trend.png\n")
cat("  P5_ARG_VF_correlation.png\n")

cat("\n================================================================================\n")
