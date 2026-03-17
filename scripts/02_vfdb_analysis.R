#!/usr/bin/env Rscript
################################################################################
# VFDB ANALYSIS PIPELINE - Main Script
# Comprehensive Virulence Factor Database Analysis
################################################################################

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(data.table)
  library(scales)
  library(viridis)
  library(patchwork)
  library(RColorBrewer)
  library(ggridges)
  library(vegan)
  library(MASS)
  library(broom)
  library(corrplot)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(rstatix)
  library(ggpubr)
  library(flextable)
  library(officer)
  library(epitools)
  library(ggrepel)
  library(ggsci)
  library(tidytext)
  library(colorspace)
})

# Configuration
n_cores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("Using %d cores for parallel processing\n\n", n_cores))

STs      <- c("ST10", "ST131", "ST69", "ST73", "ST95")
VFDB_DIR <- "data/vfdb_summary"
META_DIR <- "data/metadata"
OUT_DIR  <- "outputs/VFDB_analysis_outputs"

# Create output directories
dir.create(file.path(OUT_DIR, "plots"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "stats"),  recursive = TRUE, showWarnings = FALSE)

# Source utility functions
source("scripts/utils/common_utils.R")
source("scripts/utils/metadata_matching.R")

# Set default theme
theme_set(theme_pub())

################################################################################
# SECTION 1: METADATA MATCHING
################################################################################

cat("================================================================\n")
cat("SECTION 1: METADATA MATCHING\n")
cat("================================================================\n\n")

source("scripts/modules/vfdb/01_vfdb_loading.R")
vfdb_results <- load_vfdb_data(STs, VFDB_DIR, META_DIR)

data <- vfdb_results$data
matching_stats <- vfdb_results$matching_stats

################################################################################
# SECTION 2: VF CLASS ANALYSIS
################################################################################

cat("\n================================================================\n")
cat("SECTION 2: VIRULENCE FACTOR CLASS ANALYSIS\n")
cat("================================================================\n\n")

source("scripts/modules/vfdb/02_vf_class_analysis.R")
class_results <- analyze_vf_classes(data, STs, OUT_DIR)

################################################################################
# SECTION 3: PREVALENCE ANALYSIS
################################################################################

cat("\n================================================================\n")
cat("SECTION 3: GENE PREVALENCE ANALYSIS\n")
cat("================================================================\n\n")

source("scripts/modules/vfdb/03_prevalence_analysis.R")
prevalence_results <- analyze_vf_prevalence(data, STs, OUT_DIR)

################################################################################
# SECTION 4: TEMPORAL TRENDS
################################################################################

cat("\n================================================================\n")
cat("SECTION 4: TEMPORAL TREND ANALYSIS\n")
cat("================================================================\n\n")

source("scripts/modules/vfdb/04_temporal_vf_trends.R")
temporal_results <- analyze_vf_temporal_trends(data, STs, OUT_DIR)

################################################################################
# SECTION 5: CLINICAL COMPARISON
################################################################################

cat("\n================================================================\n")
cat("SECTION 5: CLINICAL VS NON-CLINICAL COMPARISON\n")
cat("================================================================\n\n")

source("scripts/modules/vfdb/05_clinical_vf_comparison.R")
clinical_results <- analyze_vf_clinical_comparison(data, STs, OUT_DIR)

################################################################################
# SECTION 6: DIVERSITY METRICS
################################################################################

cat("\n================================================================\n")
cat("SECTION 6: VF PROFILE DIVERSITY ANALYSIS\n")
cat("================================================================\n\n")

source("scripts/modules/vfdb/06_diversity_metrics.R")
diversity_results <- analyze_vf_diversity(data, STs, OUT_DIR)

################################################################################
# SECTION 7: STATISTICAL TESTS
################################################################################

cat("\n================================================================\n")
cat("SECTION 7: COMPREHENSIVE STATISTICAL TESTING\n")
cat("================================================================\n\n")

source("scripts/modules/vfdb/07_statistical_tests.R")
stats_results <- run_vfdb_statistical_tests(data, STs, OUT_DIR)

################################################################################
# SECTION 8: EXPORT TABLES
################################################################################

cat("\n================================================================\n")
cat("SECTION 8: EXPORTING SUMMARY TABLES\n")
cat("================================================================\n\n")

source("scripts/modules/vfdb/08_export_tables.R")
export_vfdb_tables(
  data, matching_stats, class_results, prevalence_results,
  temporal_results, clinical_results, diversity_results,
  stats_results, OUT_DIR
)

################################################################################
# COMPLETION SUMMARY
################################################################################

cat("\n================================================================\n")
cat("VFDB ANALYSIS COMPLETE\n")
cat("================================================================\n\n")

cat(sprintf("Output directory: %s\n\n", OUT_DIR))
cat(sprintf("Total genomes analyzed: %d\n", nrow(data)))
cat(sprintf("Sequence types: %s\n", paste(STs, collapse = ", ")))
cat(sprintf("VF genes analyzed: %d\n", 
            length(setdiff(colnames(data), 
                          c("Name", "ST", "Collection_Year", "Country", 
                            "Continent", "Source_Niche", "Total_VF_Count", 
                            "VF_Shannon")))))

cat("\nSTATISTICAL TESTS PERFORMED:\n")
cat("  ✓ Kruskal-Wallis + Bonferroni-Dunn post-hoc\n")
cat("  ✓ Jaccard PERMANOVA / adonis2\n")
cat("  ✓ NMDS ordination with stress metric\n")
cat("  ✓ Shannon diversity index\n")
cat("  ✓ Mann-Kendall tau trend tests\n")
cat("  ✓ Linear mixed models (Year × ST interaction)\n")
cat("  ✓ Fisher exact test + OR (clinical vs non-clinical)\n")
cat("  ✓ Log2 fold-change enrichment (per ST)\n")
cat("  ✓ Negative binomial GLM (IRR analysis)\n")
cat("  ✓ Logistic regression (high VF burden)\n")
cat("  ✓ Chi-squared heterogeneity tests\n")
cat("  ✓ Spearman correlation (VF-AMR if available)\n")

cat("\n================================================================\n")
