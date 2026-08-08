# R package installation for the E. coli ExPEC pipeline.
#
# Run after creating the 'wgs' conda environment:
#   conda activate wgs
#   Rscript install.R
#
# ggtree and treeio come from Bioconductor; all others from CRAN.
# ggtree is only needed for the two optional tree figure scripts (15/16).

options(timeout = 600)
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran <- c(
  "tidyverse", "data.table", "readxl", "writexl",
  "ape", "phangorn", "phytools", "vegan", "cluster",
  "broom", "scales", "ggrepel", "patchwork", "gridExtra",
  "viridis", "janitor", "png"
)

bioc <- c("ggtree")

# ---- CRAN ----
to_install <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) {
  cat("Installing CRAN packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install)
} else {
  cat("All CRAN packages already installed.\n")
}

# ---- Bioconductor (optional: tree figures only) ----
to_install_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install_bioc)) {
  cat("Installing Bioconductor packages:", paste(to_install_bioc, collapse = ", "), "\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(to_install_bioc, update = FALSE, ask = FALSE)
} else {
  cat("All Bioconductor packages already installed.\n")
}

cat("Done.\n")
