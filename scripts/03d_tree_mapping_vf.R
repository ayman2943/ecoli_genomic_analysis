#!/usr/bin/env Rscript
#
# 02b_tree_mapping_vf.R  —  Map VF shell-gene clusters onto core-tree
#
source("config.R")

if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape", repos = "https://cloud.r-project.org")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse", repos = "https://cloud.r-project.org")
library(ape); library(tidyverse)

ST <- config$TARGET_ST
CLUSTER_FILE <- list.files(config$VF_DIR, pattern = "03_shell_gene_cluster_assignments_k.*\\.csv", full.names = TRUE)
if (length(CLUSTER_FILE) == 0) stop("No cluster assignment file found in ", config$VF_DIR)
CLUSTER_FILE <- CLUSTER_FILE[1]
TREE_FILE <- config$TREE_FILE
OUT_DIR <- config$VF_DIR

for (f in c(CLUSTER_FILE, TREE_FILE)) if (!file.exists(f)) stop("Missing: ", f)

clusters <- read_csv(CLUSTER_FILE, show_col_types = FALSE)
tree <- read.tree(TREE_FILE)

clusters$genome_id <- gsub(" ", "_", clusters$genome_id)
tree$tip.label <- gsub(" ", "_", tree$tip.label)

common <- intersect(tree$tip.label, clusters$genome_id)
cat(sprintf("Tree tips: %d  Cluster genomes: %d  Overlap: %d\n",
            length(tree$tip.label), nrow(clusters), length(common)))
if (length(common) < 50) stop("Too few overlapping genomes")

tree_pruned <- keep.tip(tree, common)
cluster_vec <- setNames(clusters$shell_cluster[match(common, clusters$genome_id)], common)

write.tree(tree_pruned, file.path(OUT_DIR, "pruned_tree.nwk"))
write_csv(tibble(genome_id = names(cluster_vec), shell_cluster = cluster_vec),
          file.path(OUT_DIR, "tree_cluster_mapping.csv"))
cat("Tree mapping done.\n")
