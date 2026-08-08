#!/usr/bin/env Rscript

# ============================================================
# Map PPanGGOLiN shell-gene clusters onto core-genome tree
#
# Required files:
#   ST69_bootstrap.treefile
#   ST69_shell_gene_sublineage_validation_CORRECTED/03_shell_gene_cluster_assignments_k6.csv
#
# Optional:
#   ST69_shell_gene_sublineage_validation_CORRECTED/04_master_shell_cluster_metadata_VF_table.csv
#
# Main questions:
#   1. Do shell-gene clusters form monophyletic clades on the core tree?
#   2. Are they compact or scattered across the tree?
#   3. Is shell-cluster structure significantly associated with core phylogeny?
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

packages <- c(
  "ape", "phytools", "tidyverse", "ggplot2",
  "ggtree", "phangorn", "writexl"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("ggtree")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(p, ask = FALSE, update = FALSE)
    } else {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
}

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(tidyverse)
  library(ggplot2)
  library(ggtree)
  library(phangorn)
  library(writexl)
})

cat("Packages loaded\n")

# ============================================================
# Config (env-driven, from repo-root config.R)
# ============================================================
if (file.exists("config.R")) {
  source("config.R")
} else {
  source("../config.R")
}


# ============================================================
# 1. Input files
# ============================================================

TREE_FILE <- config$st_tree()

# Shell-cluster assignment / master tables written by scripts 03c/03d into
# output/{ST}/vfdb_analysis/ (and legacy dev locations kept for back-compat).
cluster_dirs <- c(config$st_out_vfdb(),
                  file.path(config$OUTPUT_DIR, "ST69_shell_gene_sublineage_validation_CORRECTED"),
                  file.path(config$OUTPUT_DIR, "ST69_shell_gene_sublineage_validation_FIXED"))
possible_cluster_files <- c(
  file.path(cluster_dirs[1], "03_shell_gene_cluster_assignments_k4.csv"),
  file.path(cluster_dirs[1], "03_shell_gene_cluster_assignments_k6.csv"),
  file.path(cluster_dirs[2], "03_shell_gene_cluster_assignments_k4.csv"),
  file.path(cluster_dirs[2], "03_shell_gene_cluster_assignments_k6.csv"),
  file.path(cluster_dirs[3], "03_shell_gene_cluster_assignments_k6.csv"),
  "03_shell_gene_cluster_assignments_k4.csv",
  "03_shell_gene_cluster_assignments_k6.csv"
)

CLUSTER_FILE <- possible_cluster_files[file.exists(possible_cluster_files)][1]

possible_master_files <- c(
  file.path(cluster_dirs[1], "04_master_shell_cluster_metadata_VF_table.csv"),
  file.path(cluster_dirs[2], "04_master_shell_cluster_metadata_VF_table.csv"),
  file.path(cluster_dirs[3], "04_master_shell_cluster_metadata_VF_table.csv"),
  "04_master_shell_cluster_metadata_VF_table.csv"
)

MASTER_FILE <- possible_master_files[file.exists(possible_master_files)][1]

OUT_DIR <- paste0(config$TARGET_ST, "_shell_clusters_mapped_to_core_tree")
dir.create(OUT_DIR, showWarnings = FALSE)

if (!file.exists(TREE_FILE)) {
  stop("Tree file not found: ", TREE_FILE)
}

if (is.na(CLUSTER_FILE)) {
  stop("Shell-cluster assignment file not found.")
}

cat("Using tree:", TREE_FILE, "\n")
cat("Using cluster file:", CLUSTER_FILE, "\n")

if (!is.na(MASTER_FILE)) {
  cat("Using master file:", MASTER_FILE, "\n")
} else {
  cat("No master table found. Tree mapping will use cluster assignments only.\n")
}


# ============================================================
# 2. Helper functions
# ============================================================

get_desc_tips <- function(tree, node) {
  ape::extract.clade(tree, node)$tip.label
}

is_pure_node <- function(tree, node, target_tips) {
  tips <- get_desc_tips(tree, node)
  all(tips %in% target_tips)
}

parent_of_node <- function(tree, node) {
  parent <- tree$edge[tree$edge[, 2] == node, 1]
  if (length(parent) == 0) return(NA_integer_)
  parent[1]
}

get_maximal_pure_subclades <- function(tree, target_tips) {
  n_tip <- length(tree$tip.label)
  internal_nodes <- (n_tip + 1):(n_tip + tree$Nnode)

  pure_nodes <- internal_nodes[
    sapply(internal_nodes, function(nd) is_pure_node(tree, nd, target_tips))
  ]

  if (length(pure_nodes) == 0) {
    return(tibble())
  }

  maximal_nodes <- pure_nodes[
    sapply(pure_nodes, function(nd) {
      parent <- parent_of_node(tree, nd)

      if (is.na(parent)) {
        return(TRUE)
      }

      !(parent %in% pure_nodes)
    })
  ]

  map_dfr(maximal_nodes, function(nd) {
    tips <- get_desc_tips(tree, nd)

    tibble(
      node = nd,
      n_tips = length(tips),
      tip_list = paste(tips, collapse = ";")
    )
  }) %>%
    arrange(desc(n_tips))
}

safe_is_monophyletic <- function(tree, tips) {
  if (length(tips) < 2) return(NA)
  tryCatch(
    ape::is.monophyletic(tree, tips),
    error = function(e) NA
  )
}

safe_get_mrca <- function(tree, tips) {
  if (length(tips) < 2) return(NA_integer_)
  tryCatch(
    ape::getMRCA(tree, tips),
    error = function(e) NA_integer_
  )
}

parsimony_score_multistate <- function(tree, trait_vector) {
  trait_vector <- trait_vector[tree$tip.label]
  phy_dat <- phangorn::phyDat(
    as.matrix(trait_vector),
    type = "USER",
    levels = sort(unique(as.character(trait_vector)))
  )
  phangorn::parsimony(tree, phy_dat)
}


# ============================================================
# 3. Load tree
# ============================================================

tree_raw <- ape::read.tree(TREE_FILE)

cat("Raw tree tips:", length(tree_raw$tip.label), "\n")
cat("Raw internal nodes:", tree_raw$Nnode, "\n")

# Rooting is useful for plotting/clade interpretation.
tree <- tryCatch(
  phytools::midpoint.root(tree_raw),
  error = function(e) tree_raw
)

cat("Tree loaded and midpoint-rooted if possible.\n")


# ============================================================
# 4. Load shell cluster assignments
# ============================================================

cluster_assign <- readr::read_csv(CLUSTER_FILE, show_col_types = FALSE)

if (!all(c("genome_id", "shell_cluster") %in% colnames(cluster_assign))) {
  stop("Cluster file must contain columns: genome_id, shell_cluster")
}

cluster_assign <- cluster_assign %>%
  mutate(
    genome_id = as.character(genome_id),
    shell_cluster = as.factor(shell_cluster)
  ) %>%
  distinct(genome_id, .keep_all = TRUE)

cat("Cluster assignment genomes:", nrow(cluster_assign), "\n")
cat("Shell clusters:", paste(levels(cluster_assign$shell_cluster), collapse = ", "), "\n")


# ============================================================
# 5. Optional master table with VF/metadata
# ============================================================

if (!is.na(MASTER_FILE)) {
  master <- readr::read_csv(MASTER_FILE, show_col_types = FALSE) %>%
    mutate(genome_id = as.character(genome_id)) %>%
    distinct(genome_id, .keep_all = TRUE)

  tip_meta <- cluster_assign %>%
    left_join(master, by = c("genome_id", "shell_cluster"))

} else {
  tip_meta <- cluster_assign
}

common_tips <- intersect(tree$tip.label, tip_meta$genome_id)

cat("Shared tree + shell-cluster genomes:", length(common_tips), "\n")

if (length(common_tips) < 50) {
  stop("Too few shared genomes between tree and shell clusters. Check genome IDs.")
}

tree <- ape::keep.tip(tree, common_tips)

tip_meta <- tip_meta %>%
  filter(genome_id %in% tree$tip.label)

# Reorder metadata to tree tip order.
tip_meta <- tip_meta %>%
  slice(match(tree$tip.label, genome_id))

write_csv(
  tip_meta,
  file.path(OUT_DIR, "01_tree_tip_shell_cluster_metadata.csv")
)

cat("Final mapped tree tips:", length(tree$tip.label), "\n")


# ============================================================
# 6. Monophyly and MRCA-purity tests
# ============================================================

clusters <- sort(unique(as.character(tip_meta$shell_cluster)))

cluster_monophyly <- map_dfr(clusters, function(cl) {
  cl_tips <- tip_meta %>%
    filter(shell_cluster == cl) %>%
    pull(genome_id)

  mrca_node <- safe_get_mrca(tree, cl_tips)

  if (is.na(mrca_node)) {
    return(tibble(
      shell_cluster = cl,
      n_cluster_tips = length(cl_tips),
      is_monophyletic = NA,
      mrca_node = NA_integer_,
      mrca_total_tips = NA_integer_,
      mrca_cluster_tips = NA_integer_,
      mrca_noncluster_tips = NA_integer_,
      mrca_purity = NA_real_
    ))
  }

  mrca_tips <- get_desc_tips(tree, mrca_node)

  tibble(
    shell_cluster = cl,
    n_cluster_tips = length(cl_tips),
    is_monophyletic = safe_is_monophyletic(tree, cl_tips),
    mrca_node = mrca_node,
    mrca_total_tips = length(mrca_tips),
    mrca_cluster_tips = sum(mrca_tips %in% cl_tips),
    mrca_noncluster_tips = sum(!mrca_tips %in% cl_tips),
    mrca_purity = mrca_cluster_tips / mrca_total_tips
  )
}) %>%
  arrange(desc(mrca_purity), desc(n_cluster_tips))

write_csv(
  cluster_monophyly,
  file.path(OUT_DIR, "02_cluster_monophyly_and_mrca_purity.csv")
)

cat("\nCluster monophyly / MRCA purity:\n")
print(cluster_monophyly)


# ============================================================
# 7. Fragmentation: maximal pure subclades per shell cluster
# ============================================================

pure_subclades <- map_dfr(clusters, function(cl) {
  cl_tips <- tip_meta %>%
    filter(shell_cluster == cl) %>%
    pull(genome_id)

  sub <- get_maximal_pure_subclades(tree, cl_tips)

  if (nrow(sub) == 0) {
    return(tibble(
      shell_cluster = cl,
      pure_subclade_node = NA_integer_,
      pure_subclade_n_tips = 0,
      tip_list = ""
    ))
  }

  sub %>%
    mutate(
      shell_cluster = cl,
      pure_subclade_node = node,
      pure_subclade_n_tips = n_tips
    ) %>%
    select(shell_cluster, pure_subclade_node, pure_subclade_n_tips, tip_list)
})

pure_subclade_summary <- pure_subclades %>%
  filter(!is.na(pure_subclade_node)) %>%
  group_by(shell_cluster) %>%
  summarise(
    n_pure_subclades = n(),
    largest_pure_subclade = max(pure_subclade_n_tips, na.rm = TRUE),
    total_tips_in_pure_subclades = sum(pure_subclade_n_tips, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(
    cluster_monophyly %>% select(shell_cluster, n_cluster_tips),
    by = "shell_cluster"
  ) %>%
  mutate(
    n_pure_subclades = replace_na(n_pure_subclades, 0L),
    largest_pure_subclade = replace_na(largest_pure_subclade, 0),
    total_tips_in_pure_subclades = replace_na(total_tips_in_pure_subclades, 0),
    largest_pure_subclade_fraction = largest_pure_subclade / n_cluster_tips
  ) %>%
  arrange(desc(largest_pure_subclade_fraction))

write_csv(
  pure_subclades,
  file.path(OUT_DIR, "03_maximal_pure_subclades_by_cluster.csv")
)

write_csv(
  pure_subclade_summary,
  file.path(OUT_DIR, "04_pure_subclade_fragmentation_summary.csv")
)

cat("\nPure-subclade fragmentation summary:\n")
print(pure_subclade_summary)


# ============================================================
# 8. Phylogenetic signal by parsimony permutation
# ============================================================

set.seed(123)

trait <- tip_meta$shell_cluster
names(trait) <- tip_meta$genome_id
trait <- trait[tree$tip.label]

observed_parsimony <- parsimony_score_multistate(tree, trait)

N_PERM <- 999

perm_scores <- replicate(N_PERM, {
  perm_trait <- sample(trait)
  names(perm_trait) <- names(trait)
  parsimony_score_multistate(tree, perm_trait)
})

p_lower <- (sum(perm_scores <= observed_parsimony) + 1) / (N_PERM + 1)

parsimony_result <- tibble(
  test = "multistate_shell_cluster_parsimony",
  observed_parsimony = observed_parsimony,
  permutation_mean = mean(perm_scores),
  permutation_sd = sd(perm_scores),
  empirical_p_lower = p_lower,
  interpretation = case_when(
    p_lower < 0.001 ~ "Very strong phylogenetic structure",
    p_lower < 0.01 ~ "Strong phylogenetic structure",
    p_lower < 0.05 ~ "Moderate phylogenetic structure",
    TRUE ~ "No strong phylogenetic structure"
  )
)

write_csv(
  parsimony_result,
  file.path(OUT_DIR, "05_parsimony_phylogenetic_signal.csv")
)

write_csv(
  tibble(perm_parsimony_score = perm_scores),
  file.path(OUT_DIR, "06_parsimony_permutation_scores.csv")
)

cat("\nParsimony phylogenetic-signal test:\n")
print(parsimony_result)


# ============================================================
# 9. Binary parsimony signal per cluster
# ============================================================

binary_parsimony <- map_dfr(clusters, function(cl) {
  binary_trait <- ifelse(trait == cl, "in_cluster", "other")
  names(binary_trait) <- names(trait)

  obs <- parsimony_score_multistate(tree, binary_trait)

  perms <- replicate(N_PERM, {
    perm_trait <- sample(binary_trait)
    names(perm_trait) <- names(binary_trait)
    parsimony_score_multistate(tree, perm_trait)
  })

  p <- (sum(perms <= obs) + 1) / (N_PERM + 1)

  tibble(
    shell_cluster = cl,
    observed_binary_parsimony = obs,
    permutation_mean = mean(perms),
    permutation_sd = sd(perms),
    empirical_p_lower = p,
    interpretation = case_when(
      p < 0.001 ~ "Very strong tree clustering",
      p < 0.01 ~ "Strong tree clustering",
      p < 0.05 ~ "Moderate tree clustering",
      TRUE ~ "No strong tree clustering"
    )
  )
}) %>%
  arrange(empirical_p_lower)

write_csv(
  binary_parsimony,
  file.path(OUT_DIR, "07_binary_cluster_phylogenetic_signal.csv")
)

cat("\nBinary cluster phylogenetic-signal tests:\n")
print(binary_parsimony)


# ============================================================
# 10. Tree plots
# ============================================================

p_tree <- ggtree(tree, size = 0.12, color = "grey60") %<+% tip_meta +
  geom_tippoint(aes(color = shell_cluster), size = 0.9, alpha = 0.9) +
  theme_tree2() +
  labs(
    title = "PPanGGOLiN shell-gene clusters mapped onto ST69 core-genome tree",
    subtitle = paste0(
      "Shared genomes mapped: ", length(tree$tip.label),
      " | Test whether shell clusters are evolutionary clades"
    ),
    color = "Shell cluster"
  )

ggsave(
  file.path(OUT_DIR, "Figure_1_shell_clusters_on_core_tree.pdf"),
  p_tree,
  width = 14,
  height = 18
)

ggsave(
  file.path(OUT_DIR, "Figure_1_shell_clusters_on_core_tree.png"),
  p_tree,
  width = 14,
  height = 18,
  dpi = 250
)

p_tree_circular <- ggtree(tree, layout = "circular", size = 0.10, color = "grey65") %<+% tip_meta +
  geom_tippoint(aes(color = shell_cluster), size = 0.7, alpha = 0.9) +
  labs(
    title = "Circular tree: shell clusters on ST69 core phylogeny",
    color = "Shell cluster"
  ) +
  theme(legend.position = "right")

ggsave(
  file.path(OUT_DIR, "Figure_2_shell_clusters_on_core_tree_circular.pdf"),
  p_tree_circular,
  width = 12,
  height = 12
)

ggsave(
  file.path(OUT_DIR, "Figure_2_shell_clusters_on_core_tree_circular.png"),
  p_tree_circular,
  width = 12,
  height = 12,
  dpi = 250
)

# Optional VF-colored tree if total_vf is present.
if ("total_vf" %in% colnames(tip_meta)) {
  p_vf <- ggtree(tree, size = 0.12, color = "grey60") %<+% tip_meta +
    geom_tippoint(aes(color = total_vf), size = 0.9, alpha = 0.9) +
    scale_color_gradient(low = "grey80", high = "red3", na.value = "grey90") +
    theme_tree2() +
    labs(
      title = "Total VF burden mapped onto ST69 core-genome tree",
      subtitle = "Compare VF burden with core phylogenetic structure",
      color = "Total VF"
    )

  ggsave(
    file.path(OUT_DIR, "Figure_3_total_VF_on_core_tree.pdf"),
    p_vf,
    width = 14,
    height = 18
  )

  ggsave(
    file.path(OUT_DIR, "Figure_3_total_VF_on_core_tree.png"),
    p_vf,
    width = 14,
    height = 18,
    dpi = 250
  )
}


# ============================================================
# 11. Final interpretation table
# ============================================================

interpretation_table <- cluster_monophyly %>%
  left_join(pure_subclade_summary, by = c("shell_cluster", "n_cluster_tips")) %>%
  left_join(binary_parsimony, by = "shell_cluster") %>%
  mutate(
    final_interpretation = case_when(
      is_monophyletic == TRUE ~
        "Can be described as a core-tree monophyletic evolutionary sublineage",

      is_monophyletic == FALSE &
        mrca_purity >= 0.80 &
        largest_pure_subclade_fraction >= 0.70 &
        empirical_p_lower < 0.05 ~
        "Mostly tree-localized; can be cautiously described as a phylogenetically structured sublineage",

      is_monophyletic == FALSE &
        empirical_p_lower < 0.05 ~
        "Phylogenetically structured but not monophyletic; describe as shell-gene-defined subgroup/profile",

      TRUE ~
        "Not strongly localized on the core tree; describe as accessory-gene profile rather than evolutionary sublineage"
    )
  )

write_csv(
  interpretation_table,
  file.path(OUT_DIR, "08_FINAL_cluster_tree_mapping_interpretation.csv")
)

cat("\nFinal cluster/tree interpretation:\n")
print(interpretation_table)


# ============================================================
# 12. Excel workbook and report
# ============================================================

write_xlsx(
  list(
    tip_metadata = tip_meta,
    monophyly_purity = cluster_monophyly,
    pure_subclades = pure_subclades,
    fragmentation = pure_subclade_summary,
    parsimony_signal = parsimony_result,
    binary_signal = binary_parsimony,
    final_interpretation = interpretation_table
  ),
  path = file.path(OUT_DIR, "ST69_shell_clusters_core_tree_mapping_results.xlsx")
)

report_file <- file.path(OUT_DIR, "TREE_MAPPING_REPORT.txt")

sink(report_file)

cat("ST69 shell clusters mapped onto core-genome tree\n")
cat("================================================\n\n")

cat("Tree file: ", TREE_FILE, "\n", sep = "")
cat("Cluster file: ", CLUSTER_FILE, "\n", sep = "")
cat("Tree tips mapped: ", length(tree$tip.label), "\n\n", sep = "")

cat("Monophyly and MRCA purity:\n")
print(cluster_monophyly)

cat("\nPure-subclade fragmentation:\n")
print(pure_subclade_summary)

cat("\nMultistate parsimony phylogenetic signal:\n")
print(parsimony_result)

cat("\nBinary cluster phylogenetic signal:\n")
print(binary_parsimony)

cat("\nFinal interpretation:\n")
print(interpretation_table)

cat("\nHow to use this:\n")
cat("If clusters are monophyletic or highly tree-localized, you can describe them as evolutionary sublineages.\n")
cat("If clusters are scattered but still phylogenetically structured, describe them as shell-gene-defined accessory-genome subgroups.\n")
cat("If clusters are scattered with weak signal, avoid evolutionary sublineage language.\n")

sink()


# ============================================================
# 13. Final message
# ============================================================

cat("\nDONE.\n")
cat("Outputs saved in:", OUT_DIR, "\n\n")
cat("Most important files:\n")
cat("  02_cluster_monophyly_and_mrca_purity.csv\n")
cat("  04_pure_subclade_fragmentation_summary.csv\n")
cat("  05_parsimony_phylogenetic_signal.csv\n")
cat("  07_binary_cluster_phylogenetic_signal.csv\n")
cat("  08_FINAL_cluster_tree_mapping_interpretation.csv\n")
cat("  TREE_MAPPING_REPORT.txt\n")
cat("  Figure_1_shell_clusters_on_core_tree.pdf/png\n")
cat("  ST69_shell_clusters_core_tree_mapping_results.xlsx\n")