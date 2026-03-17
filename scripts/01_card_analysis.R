#!/usr/bin/env Rscript
################################################################################
# AMR ANALYSIS PIPELINE - Main Script
# Comprehensive Antibiotic Resistance Database (CARD) Analysis
################################################################################

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(scales)
  library(RColorBrewer)
  library(viridis)
  library(ggpubr)
  library(gridExtra)
  library(patchwork)
  library(ggsci)
  library(data.table)
  library(parallel)
  library(stringi)
  library(ggridges)
  library(ggtext)
  library(broom)
  library(ggrepel)
  library(vegan)
  library(MASS)
  library(rstatix)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(corrplot)
  library(flextable)
  library(officer)
  library(epitools)
  library(tidytext)
})

# Configuration
n_cores <- max(1, detectCores() - 1)
cat(sprintf("Using %d cores for parallel processing\n\n", n_cores))

STs      <- c("ST10", "ST131", "ST69", "ST73", "ST95")
CARD_DIR <- "data/card_summary"
META_DIR <- "data/metadata"
OUT_DIR  <- "outputs/R_analysis_outputs"

# Create output directories
dir.create(file.path(OUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

# Source utility functions
source("scripts/utils/common_utils.R")
source("scripts/utils/metadata_matching.R")

# Set default theme
theme_set(theme_modern())

################################################################################
# SECTION 1: METADATA MATCHING
################################################################################

cat("================================================================================\n")
cat("SECTION 1: METADATA MATCHING\n")
cat("================================================================================\n\n")

# Source and run metadata matching module
source("scripts/modules/amr/01_metadata_matching.R")
matching_results <- run_metadata_matching(STs, CARD_DIR, META_DIR)

data <- matching_results$data
matching_summary <- matching_results$summary
overall_rate <- with(matching_summary, sum(Matched) / sum(Total) * 100)

################################################################################
# SECTION 2: DATA LOADING & VALIDATION
################################################################################

cat("\n================================================================================\n")
cat("SECTION 2: DATA LOADING & GENE MATRIX CONSTRUCTION\n")
cat("================================================================================\n\n")

source("scripts/modules/amr/02_data_loading.R")
gene_data <- load_and_process_gene_data(data, STs)

# Update data with gene counts
data <- gene_data$data
pres_gene_bin <- gene_data$pres_gene_bin
gene_cols <- gene_data$gene_cols

################################################################################
# SECTION 3: BASE VISUALIZATIONS
################################################################################

cat("\n================================================================================\n")
cat("SECTION 3: GENERATING BASE AMR VISUALIZATIONS\n")
cat("================================================================================\n\n")

source("scripts/modules/amr/03_base_visualizations.R")
generate_base_visualizations(data, OUT_DIR)

################################################################################
# SECTION 4: ENRICHMENT ANALYSIS
################################################################################

cat("\n================================================================================\n")
cat("SECTION 4: ST-SPECIFIC GENE ENRICHMENT ANALYSIS\n")
cat("================================================================================\n\n")

source("scripts/modules/amr/04_enrichment_analysis.R")
fc_results <- run_enrichment_analysis(data, pres_gene_bin, gene_cols, STs, OUT_DIR)

################################################################################
# SECTION 5: TEMPORAL ANALYSIS
################################################################################

cat("\n================================================================================\n")
cat("SECTION 5: TEMPORAL TREND ANALYSIS\n")
cat("================================================================================\n\n")

source("scripts/modules/amr/05_temporal_analysis.R")
temporal_results <- run_temporal_analysis(data, STs, OUT_DIR)

################################################################################
# SECTION 6: CLINICAL VS NON-CLINICAL COMPARISON
################################################################################

cat("\n================================================================================\n")
cat("SECTION 6: CLINICAL VS NON-CLINICAL COMPARISON\n")
cat("================================================================================\n\n")

source("scripts/modules/amr/06_clinical_comparison.R")
clinical_results <- run_clinical_comparison(
  data, pres_gene_bin, gene_cols, STs, OUT_DIR
)

################################################################################
# SECTION 7: DIVERSITY & ORDINATION
################################################################################

cat("\n================================================================================\n")
cat("SECTION 7: DIVERSITY METRICS & ORDINATION\n")
cat("================================================================================\n\n")

source("scripts/modules/amr/07_diversity_analysis.R")
diversity_results <- run_diversity_analysis(data, pres_gene_bin, gene_cols, OUT_DIR)

################################################################################
# SECTION 8: CO-OCCURRENCE NETWORK
################################################################################

cat("\n================================================================================\n")
cat("SECTION 8: GENE CO-OCCURRENCE NETWORK ANALYSIS\n")
cat("================================================================================\n\n")

source("scripts/modules/amr/08_cooccurrence_network.R")
network_results <- run_network_analysis(pres_gene_bin, gene_cols, OUT_DIR)

################################################################################
# SECTION 9: STATISTICAL MODELS
################################################################################

cat("\n================================================================================\n")
cat("SECTION 9: ADVANCED STATISTICAL MODELS\n")
cat("================================================================================\n\n")

source("scripts/modules/amr/09_statistical_models.R")
model_results <- run_statistical_models(data, OUT_DIR)

################################################################################
# SECTION 10: EXPORT TABLES
################################################################################

cat("\n================================================================================\n")
cat("SECTION 10: EXPORTING SUMMARY TABLES\n")
cat("================================================================================\n\n")

source("scripts/modules/amr/10_export_tables.R")
export_summary_tables(
  data, pres_gene_bin, gene_cols, STs,
  matching_summary, fc_results, clinical_results,
  temporal_results, diversity_results, network_results,
  model_results, OUT_DIR
)

################################################################################
# COMPLETION SUMMARY
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("AMR ANALYSIS COMPLETE!\n")
cat(rep("=", 80), "\n\n", sep = "")

cat(sprintf("  Total genomes analyzed  : %d\n", nrow(data)))
cat(sprintf("  Sequence types          : %s\n", paste(unique(data$ST), collapse = ", ")))
cat(sprintf("  Year range              : %d - %d\n",
            min(data$Collection_Year, na.rm = TRUE), 
            max(data$Collection_Year, na.rm = TRUE)))
cat(sprintf("  Overall match rate      : %.1f%%\n", overall_rate))
cat(sprintf("  Gene columns analyzed   : %d\n", length(gene_cols)))

cat("\nOutputs saved to", OUT_DIR, "\n")
cat("  - plots/   : Publication-quality visualizations (PNG, 300 DPI)\n")
cat("  - tables/  : Statistical summaries (Excel, CSV, DOCX)\n")

cat("\n", rep("=", 80), "\n", sep = "")
