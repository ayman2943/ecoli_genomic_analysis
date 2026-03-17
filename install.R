#!/usr/bin/env Rscript
################################################################################
# E. coli Genomic Analysis Pipeline - Installation Script
# Automatically installs all required R packages
################################################################################

cat("================================================================================\n")
cat("E. coli Genomic Analysis Pipeline - Package Installation\n")
cat("================================================================================\n\n")

# Check R version
r_version <- getRversion()
min_version <- "4.0.0"

cat(sprintf("Current R version: %s\n", r_version))
cat(sprintf("Minimum required: %s\n\n", min_version))

if (r_version < min_version) {
  stop(sprintf("R version %s or higher is required. Please upgrade R.", min_version))
}

# List of required packages
packages <- c(
  # Core tidyverse
  "tidyverse",
  "readxl",
  "writexl",
  "data.table",
  
  # Visualization
  "ggplot2",
  "scales",
  "RColorBrewer",
  "viridis",
  "ggpubr",
  "gridExtra",
  "patchwork",
  "ggsci",
  "ggridges",
  "ggtext",
  "ggrepel",
  "colorspace",
  "pheatmap",
  
  # Statistical analysis
  "vegan",
  "MASS",
  "rstatix",
  "broom",
  "epitools",
  
  # Network analysis
  "igraph",
  "ggraph",
  "tidygraph",
  "corrplot",
  
  # Text and data
  "stringi",
  "tidytext",
  
  # Reporting
  "flextable",
  "officer"
)

cat("Checking for required packages...\n\n")

# Check which packages are already installed
installed <- installed.packages()[, "Package"]
to_install <- packages[!packages %in% installed]

if (length(to_install) == 0) {
  cat("✓ All required packages are already installed!\n\n")
} else {
  cat(sprintf("Installing %d missing packages:\n", length(to_install)))
  cat(paste("  -", to_install, collapse = "\n"), "\n\n")
  
  # Install missing packages
  install.packages(
    to_install,
    dependencies = TRUE,
    repos = "https://cloud.r-project.org/"
  )
  
  cat("\n✓ Package installation complete!\n\n")
}

# Verify installation
cat("Verifying package installation...\n\n")

failed <- character()
for (pkg in packages) {
  result <- tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    TRUE
  }, error = function(e) {
    FALSE
  })
  
  if (result) {
    cat(sprintf("  ✓ %s\n", pkg))
  } else {
    cat(sprintf("  ✗ %s (FAILED)\n", pkg))
    failed <- c(failed, pkg)
  }
}

cat("\n")

if (length(failed) > 0) {
  cat("================================================================================\n")
  cat("WARNING: The following packages failed to install:\n")
  cat(paste("  -", failed, collapse = "\n"), "\n\n")
  cat("Please try installing them manually:\n")
  cat(sprintf("  install.packages(c(%s))\n", 
              paste0('"', failed, '"', collapse = ", ")))
  cat("================================================================================\n")
} else {
  cat("================================================================================\n")
  cat("SUCCESS: All packages installed and verified!\n")
  cat("================================================================================\n\n")
  
  cat("Next steps:\n")
  cat("  1. Prepare your data (see docs/USAGE.md)\n")
  cat("  2. Run analyses:\n")
  cat("       Rscript scripts/01_card_analysis.R\n")
  cat("       Rscript scripts/02_vfdb_analysis.R\n")
  cat("       Rscript scripts/03_plasmid_analysis.R\n\n")
}

# Print session info for reference
cat("\n================================================================================\n")
cat("Session Information\n")
cat("================================================================================\n\n")
print(sessionInfo())
