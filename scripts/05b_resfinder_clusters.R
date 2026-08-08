#!/usr/bin/env Rscript
#
# Quick analysis: ResFinder per-gene decreasing for ALL downward-trend clusters
#
suppressPackageStartupMessages({
  library(tidyverse); library(readxl)
})
source("config.R")

EARLY <- c(2016, 2017, 2018)
LATE  <- c(2022, 2023, 2024, 2025)

make_id <- function(x) {
  x <- as.character(x); x <- trimws(x)
  ifelse(grepl("^Escherichia_coli_", x), x,
         ifelse(grepl("^E\\.coli_", x), sub("^E\\.coli_", "Escherichia_coli_", x),
                paste0("Escherichia_coli_", x)))
}

cat("Loading...\n")
meta <- read_xlsx(config$st_metadata("ST10")) %>%
  rename(id = Name) %>%
  mutate(id = make_id(id),
         clinical = .data[["Source Niche"]] == "Human",
         year = as.integer(.data[["Collection Year"]])) %>%
  filter(!is.na(year), !is.na(clinical))

rf <- read.delim(file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
                  header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(st == "ST10") %>%
  mutate(id = make_id(genome)) %>%
  select(-st, -genome)

genes <- setdiff(colnames(rf), "id")

clust <- read.csv("output/ST10/vfdb_analysis/03_shell_gene_cluster_assignments_k9.csv",
                   stringsAsFactors = FALSE) %>% rename(id = 1, shell_cluster = 2)

rf_long <- rf %>% pivot_longer(-id, names_to = "gene", values_to = "present") %>%
  mutate(present = as.integer(!is.na(present) & present != "" & present != "0")) %>%
  inner_join(meta %>% select(id, clinical, year), by = "id") %>%
  inner_join(clust %>% select(id, shell_cluster), by = "id")

analyze_cluster <- function(d, clust_name) {
  d <- d %>% filter(shell_cluster == clust_name) %>%
    mutate(period = case_when(year %in% EARLY ~ "early", year %in% LATE ~ "late")) %>%
    filter(!is.na(period))
  n_early <- n_distinct(d$id[d$period == "early"])
  n_late  <- n_distinct(d$id[d$period == "late"])
  cat("  n_early:", n_early, "n_late:", n_late, "\n")

  per_gene <- d %>% group_by(gene) %>%
    summarise(pe = sum(present[period == "early"]), pl = sum(present[period == "late"]),
              .groups = "drop") %>%
    rowwise() %>%
    mutate(fr = list(tryCatch(fisher.test(matrix(c(pl, n_late-pl, pe, n_early-pe), 2, byrow = TRUE)),
                               error = function(e) NULL))) %>%
    ungroup()

  per_gene$or <- sapply(per_gene$fr, function(x) if (is.null(x)) NA else unname(x$estimate))
  per_gene$p <- sapply(per_gene$fr, function(x) if (is.null(x)) NA else x$p.value)
  per_gene$fr <- NULL

  per_gene <- per_gene %>% mutate(
    p_adj = p.adjust(p, "BH"),
    prev_early = pe / n_early,
    prev_late  = pl / n_late,
    delta_pp = (prev_late - prev_early) * 100,
    direction = case_when(delta_pp < 0 & p_adj < 0.05 ~ "Decreasing",
                          delta_pp > 0 & p_adj < 0.05 ~ "Increasing",
                          TRUE ~ "No change")
  ) %>% arrange(delta_pp)

  dec <- per_gene %>% filter(direction == "Decreasing")
  inc <- per_gene %>% filter(direction == "Increasing")
  cat("  Decreasing:", nrow(dec), "Increasing:", nrow(inc), "\n")
  if (nrow(dec) > 0) {
    cat("  Top decreasing:\n")
    print(dec %>% select(gene, prev_early, prev_late, delta_pp, or, p_adj), row.names = FALSE)
  }
  if (nrow(inc) > 0) {
    cat("  Top increasing:\n")
    print(inc %>% select(gene, prev_early, prev_late, delta_pp, or, p_adj), row.names = FALSE)
  }
  per_gene
}

for (cl in c("Cluster_1", "Cluster_5", "Cluster_8")) {
  cat("\n===== ", cl, " =====\n")
  n <- n_distinct(rf_long$id[rf_long$shell_cluster == cl])
  cat("  Total genomes:", n, "\n")
  r <- analyze_cluster(rf_long, cl)
  write.csv(r, paste0("output/resfinder_decreasing_analysis/", cl, "_gene_early_vs_late.csv"), row.names = FALSE)
}
cat("\nDone.\n")
