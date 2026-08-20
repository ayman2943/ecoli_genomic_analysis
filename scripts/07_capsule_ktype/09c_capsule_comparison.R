#!/usr/bin/env Rscript
# Compare kps prevalence: old (master table, numeric parsing) vs new (VFDB summary, ABRicate parsing)
suppressPackageStartupMessages({library(tidyverse); library(writexl); library(readxl)})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_capsule_classification")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

kps <- c("kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT","kpsU")

# Old approach: master table, numeric >= 90
master <- read_csv(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv"), show_col_types=FALSE)
old <- master %>% select(genome_id, any_of(kps)) %>%
  mutate(across(any_of(kps), ~{
    v <- suppressWarnings(as.numeric(.))
    as.integer(!is.na(v) & v >= 90)
  }))

# New approach: VFDB summary, ABRicate dot parsing
vfdb <- read_tsv(config$st_vfdb_summary("ST69"), show_col_types=FALSE,
  progress=FALSE, col_types=cols(.default="c")) %>%
  filter(str_detect(.data[["#FILE"]], "^ST69/")) %>%
  mutate(genome_id = str_remove(str_remove(.data[["#FILE"]], "^ST69/"), "_vfdb.tsv$"))
new_g <- intersect(kps, colnames(vfdb))
new <- vfdb %>% select(genome_id, all_of(new_g)) %>%
  mutate(across(all_of(new_g), ~as.integer(!is.na(.x) & .x != ".")))

# Join and compute comparison
cmp <- inner_join(old, new, by="genome_id", suffix=c("_old", "_new"))

comparison <- map_dfr(new_g, ~{
  o <- cmp[[paste0(.x, "_old")]]
  n <- cmp[[paste0(.x, "_new")]]
  tibble(gene=.x,
    old_prev = mean(o, na.rm=TRUE) * 100,
    new_prev = mean(n, na.rm=TRUE) * 100,
    diff_pp = old_prev - new_prev,
    n_agree = sum(o == n, na.rm=TRUE),
    n_total = nrow(cmp),
    pct_agree = n_agree / n_total * 100)
})

cat("\n=== kps prevalence comparison: old vs new ===\n")
print(comparison, n=10)

# Also show clinical enrichment overlap
cat("\n=== Clinical enrichment: both methods agree ===\n")
old_clin <- old %>% inner_join(select(master, genome_id, clinical_binary), by="genome_id") %>%
  filter(!is.na(clinical_binary)) %>%
  mutate(clin = as.integer(clinical_binary))
new_clin <- new %>% inner_join(select(master, genome_id, clinical_binary), by="genome_id") %>%
  filter(!is.na(clinical_binary)) %>%
  mutate(clin = as.integer(clinical_binary))

enrich_cmp <- map_dfr(new_g, ~{
  old_g <- old_clin[[.x]]; new_g <- new_clin[[.x]]
  c_old <- old_clin$clin; c_new <- new_clin$clin
  tibble(gene=.x,
    old_delta = (mean(old_g[c_old==1]) - mean(old_g[c_old==0])) * 100,
    new_delta = (mean(new_g[c_new==1]) - mean(new_g[c_new==0])) * 100)
})
print(enrich_cmp, n=10)

# Write to xlsx (combine all sheets from scratch)
xlsx_file <- file.path(OUT, "capsule_classification.xlsx")
all_sheets <- list(
  capsule_distribution = read_xlsx(xlsx_file, sheet="capsule_distribution"),
  clinical_enrichment = read_xlsx(xlsx_file, sheet="clinical_enrichment"),
  cluster3_increasing = read_xlsx(xlsx_file, sheet="cluster3_increasing"),
  clinical_by_capsule_type = read_xlsx(xlsx_file, sheet="clinical_by_capsule_type"),
  kps_method_comparison = comparison,
  kps_clinical_enrichment_comparison = enrich_cmp
)
write_xlsx(all_sheets, xlsx_file)
cat("\nAppended to:", xlsx_file, "\n")
