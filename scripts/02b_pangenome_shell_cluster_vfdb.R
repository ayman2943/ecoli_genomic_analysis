#!/usr/bin/env Rscript
#
# 01a_shell_cluster_vfdb.R — ST69 shell-gene clustering with VFDB validation.
# Pipeline step: 1a (VFDB track). Environment: R 4.x, packages below.
# Inputs: gene_presence_absence.Rtab, partitions/shell.txt, ST69_vfdb_summary.tsv, metadata
# Output: ST69_shell_gene_sublineage_validation_CORRECTED/
#
source("config.R")

packages <- c("tidyverse", "readxl", "vegan", "cluster", "broom", "writexl")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(vegan); library(cluster)
  library(broom); library(writexl)
})
cat("Packages loaded\n")

RTAB_FILE <- config$GENE_PA_AB
SHELL_FILE <- config$SHELL_FILE
VFDB_SUMMARY_FILE <- config$VFDB_SUMMARY
METADATA_FILE <- config$METADATA_FILE
OUT_DIR <- config$VFDB_DIR

dir.create(OUT_DIR, showWarnings = FALSE)
K_RANGE <- 4:10
MIN_YEAR_N <- 20
MAX_YEAR <- 2025
EARLY_YEARS <- c(2016, 2017, 2018)
LATE_YEARS <- c(2022, 2023, 2024, 2025)
PRIMARY_HCLUST_METHOD <- "ward.D2"
SENSITIVITY_HCLUST_METHOD <- "average"
TARGET_ST <- config$TARGET_ST

for (f in c(RTAB_FILE, SHELL_FILE, VFDB_SUMMARY_FILE, METADATA_FILE)) {
  if (!file.exists(f)) stop("Missing: ", f)
}
cat("Using Rtab:", RTAB_FILE, "\nShell:", SHELL_FILE, "\nVFDB:", VFDB_SUMMARY_FILE, "\nMetadata:", METADATA_FILE, "\n")

safe_spearman <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2)
    return(list(rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value)
}
safe_wilcox <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(NULL)
  tryCatch(wilcox.test(x, y), error = function(e) NULL)
}
safe_kruskal <- function(value, group) {
  d <- tibble(value = value, group = group) %>% filter(!is.na(value), !is.na(group))
  if (n_distinct(d$group) < 2) return(NULL)
  tryCatch(kruskal.test(value ~ group, data = d), error = function(e) NULL)
}
find_col <- function(df, patterns, fallback = NA_integer_) {
  cn <- colnames(df)
  for (p in patterns) {
    h <- grep(p, cn, ignore.case = TRUE, value = TRUE)
    if (length(h)) return(h[1])
  }
  if (!is.na(fallback) && fallback <= ncol(df)) return(cn[fallback])
  NA_character_
}
make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}
score_gene_group <- function(df, genes) {
  present <- intersect(genes, colnames(df))
  if (length(present) == 0) return(rep(NA_real_, nrow(df)))
  rowSums(df[, present, drop = FALSE], na.rm = TRUE)
}
top_values <- function(x, n = 5) {
  x <- as.character(x); x <- x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) return("")
  tb <- sort(table(x), decreasing = TRUE)
  paste(names(tb)[seq_len(min(n, length(tb)))], tb[seq_len(min(n, length(tb)))], sep = "=", collapse = "; ")
}

# Read pangenome binary matrix
rtab <- read_tsv(RTAB_FILE, show_col_types = FALSE, progress = FALSE, guess_max = 10000)
if (ncol(rtab) < 3) stop("Rtab too few columns")
colnames(rtab)[1] <- "family_id"
rtab <- rtab %>% mutate(family_id = as.character(family_id)) %>% distinct(family_id, .keep_all = TRUE)
genome_cols_all <- setdiff(colnames(rtab), "family_id")
rtab[genome_cols_all] <- lapply(rtab[genome_cols_all], function(x) {
  x <- suppressWarnings(as.integer(as.character(x)))
  x[is.na(x)] <- 0L
  as.integer(x > 0)
})
cat("Rtab families:", nrow(rtab), "\nRtab genomes:", length(genome_cols_all), "\n")

# Read shell gene partition
read_shell_file <- function(path, rtab_families) {
  lines <- readLines(path, warn = FALSE); lines <- trimws(lines); lines <- lines[lines != ""]
  tokens <- unlist(strsplit(lines, "[,;\\t ]+"))
  tokens <- trimws(tokens); tokens <- tokens[tokens != ""]
  tokens <- tokens[!grepl("^(shell|family|families|gene|genes|partition)$", tokens, ignore.case = TRUE)]
  shell_fams <- intersect(tokens, rtab_families)
  tibble(family_id = shell_fams, partition = "shell") %>% distinct(family_id, .keep_all = TRUE)
}
partition_tbl <- read_shell_file(SHELL_FILE, rtab$family_id)
cat("PPanGGOLiN shell families matched to Rtab:", nrow(partition_tbl), "\n")
if (nrow(partition_tbl) < 100) stop("Too few shell families")
shell_families <- partition_tbl$family_id
SHELL_SOURCE <- "PPanGGOLiN partitions/shell.txt"
write_csv(partition_tbl, file.path(OUT_DIR, "00_shell_families_used.csv"))

# Read VFDB summary
vfdb <- read_tsv(VFDB_SUMMARY_FILE, show_col_types = FALSE, progress = FALSE) %>%
  filter(str_detect(.data[["#FILE"]], paste0("^", TARGET_ST, "/"))) %>%
  mutate(genome = str_remove(.data[["#FILE"]], paste0("^", TARGET_ST, "/")),
         genome = str_remove(genome, "_vfdb\\.tsv$")) %>%
  mutate(genome_id = make_genome_id(genome))
cat("VFDB", TARGET_ST, "genomes:", nrow(vfdb), "\n")
vf_genes <- setdiff(colnames(vfdb), c("#FILE", "NUM_FOUND", "genome", "genome_id"))

cat("VFDB genes:", length(vf_genes), "\n")
vfdb_binary <- vfdb %>% rename(total_vf = NUM_FOUND)
# Metadata
meta_raw <- read_excel(METADATA_FILE, col_names = TRUE)
name_col <- find_col(meta_raw, c("^Name$", "genome", "strain", "isolate", "assembly", "sample"), fallback = 2)
source_col <- find_col(meta_raw, c("Source.*Niche", "source_niche", "^source$", "niche", "host"), fallback = 4)
year_col <- find_col(meta_raw, c("Collection.*Year", "^year$", "collection_year", "date"), fallback = 7)
continent_col <- find_col(meta_raw, c("^Continent$", "continent"), fallback = 11)
country_col <- find_col(meta_raw, c("^Country$", "country"), fallback = 12)
meta <- meta_raw %>%
  transmute(name_raw = .data[[name_col]], source_niche = .data[[source_col]],
            year = .data[[year_col]], continent = .data[[continent_col]],
            country = .data[[country_col]]) %>%
  filter(!is.na(name_raw)) %>%
  mutate(genome_id = make_genome_id(name_raw), source_niche = as.character(source_niche),
         continent = as.character(continent), country = as.character(country),
         year = suppressWarnings(as.integer(year)),
         niche = case_when(source_niche == "Human" ~ "clinical",
                           grepl("human|patient|clinical|blood|urine|wound|hospital|cerebrospinal|respiratory",
                                 source_niche, ignore.case = TRUE) ~ "clinical",
                           is.na(source_niche) | source_niche == "" ~ "unknown",
                           TRUE ~ "non-clinical")) %>%
  select(genome_id, niche, source_niche, year, continent, country) %>%
  distinct(genome_id, .keep_all = TRUE)
cat("Metadata genomes:", nrow(meta), "\n")

# Intersect
rtab_genomes <- setdiff(colnames(rtab), "family_id")
common_genomes <- Reduce(intersect, list(rtab_genomes, vfdb_binary$genome_id, meta$genome_id))
cat("Shared genomes before year filter:", length(common_genomes), "\n")
if (length(common_genomes) < 100) stop("Too few shared genomes")
meta_common <- meta %>% filter(genome_id %in% common_genomes, !is.na(year), year <= MAX_YEAR)
valid_years <- meta_common %>% count(year, name = "n") %>% filter(n >= MIN_YEAR_N) %>% pull(year)
common_genomes <- meta_common %>% filter(year %in% valid_years) %>% pull(genome_id)
cat("Shared genomes after filter:", length(common_genomes), "\nYears:", paste(sort(valid_years), collapse = ", "), "\n")
if (length(common_genomes) < 100) stop("Too few after filtering")

# Shell-gene Jaccard distance + clustering
rtab_shell_common <- rtab %>% filter(family_id %in% shell_families) %>%
  select(family_id, all_of(common_genomes))
shell_matrix_family_by_genome <- rtab_shell_common %>% select(-family_id)
shell_matrix_family_by_genome[] <- lapply(shell_matrix_family_by_genome, function(x) {
  x <- suppressWarnings(as.integer(as.character(x)))
  x[is.na(x)] <- 0L
  as.integer(x > 0)
})
mat_fam_genome <- as.matrix(shell_matrix_family_by_genome)
rownames(mat_fam_genome) <- rtab_shell_common$family_id
present_counts <- rowSums(mat_fam_genome)
keep_families <- present_counts > 0 & present_counts < ncol(mat_fam_genome)
mat_fam_genome <- mat_fam_genome[keep_families, , drop = FALSE]
mat_genome_family <- t(mat_fam_genome)
cat("Variable shell families:", ncol(mat_genome_family), "\n")
write_csv(tibble(statistic = c("shared_genomes", "variable_shell_families", "years_used", "shell_source"),
                 value = c(length(common_genomes), ncol(mat_genome_family),
                           paste(sort(valid_years), collapse = ";"), SHELL_SOURCE)),
          file.path(OUT_DIR, "01_input_summary.csv"))

cat("Computing Jaccard distance...\n")
jaccard_dist <- vegdist(mat_genome_family, method = "jaccard", binary = TRUE)
hc_primary <- hclust(jaccard_dist, method = PRIMARY_HCLUST_METHOD)
hc_sensitivity <- hclust(jaccard_dist, method = SENSITIVITY_HCLUST_METHOD)

silhouette_primary <- map_dfr(K_RANGE, function(k) {
  cl <- cutree(hc_primary, k = k)
  sil <- silhouette(cl, jaccard_dist)
  tibble(method = PRIMARY_HCLUST_METHOD, k = k, average_silhouette_width = mean(sil[, "sil_width"]))
})
silhouette_sensitivity <- map_dfr(K_RANGE, function(k) {
  cl <- cutree(hc_sensitivity, k = k)
  sil <- silhouette(cl, jaccard_dist)
  tibble(method = SENSITIVITY_HCLUST_METHOD, k = k, average_silhouette_width = mean(sil[, "sil_width"]))
})
silhouette_results <- bind_rows(silhouette_primary, silhouette_sensitivity)
write_csv(silhouette_results, file.path(OUT_DIR, "02_silhouette_k4_to_k10.csv"))

best_k_primary <- silhouette_primary$k[which.max(silhouette_primary$average_silhouette_width)]
best_silhouette <- max(silhouette_primary$average_silhouette_width, na.rm = TRUE)
k4_silhouette <- silhouette_primary %>% filter(k == 4) %>% pull(average_silhouette_width)
if (length(k4_silhouette) == 0) k4_silhouette <- NA_real_
k4_delta_from_best <- best_silhouette - k4_silhouette
K4_CLOSE_THRESHOLD <- 0.02
k4_supported <- !is.na(k4_silhouette) && (best_k_primary == 4 || k4_delta_from_best <= K4_CLOSE_THRESHOLD)
k_selection_summary <- tibble(
  best_k_primary = best_k_primary, best_silhouette = best_silhouette,
  k4_silhouette = k4_silhouette, k4_delta_from_best = k4_delta_from_best,
  k4_close_threshold = K4_CLOSE_THRESHOLD, k4_supported = k4_supported,
  selected_k = ifelse(k4_supported, 4, best_k_primary),
  selection_rule = case_when(best_k_primary == 4 ~ "k=4 highest silhouette",
                             k4_supported ~ "k=4 close to best",
                             TRUE ~ "best k selected"))
write_csv(k_selection_summary, file.path(OUT_DIR, "02B_k_selection_summary.csv"))
print(k_selection_summary)
SELECTED_K <- k_selection_summary$selected_k[1]
cat("Best k:", best_k_primary, " k=4 sil:", round(k4_silhouette, 4), " Selected:", SELECTED_K, "\n")

forced_clusters_raw <- cutree(hc_primary, k = SELECTED_K)
if (is.null(names(forced_clusters_raw)))
  names(forced_clusters_raw) <- rownames(mat_genome_family)
cluster_assignments <- tibble(genome_id = names(forced_clusters_raw),
                              raw_cluster = as.integer(forced_clusters_raw))
cluster_key <- cluster_assignments %>% count(raw_cluster, name = "n") %>%
  arrange(desc(n), raw_cluster) %>% mutate(shell_cluster = paste0("Cluster_", row_number()))
cluster_assignments <- cluster_assignments %>%
  left_join(cluster_key %>% select(raw_cluster, shell_cluster), by = "raw_cluster") %>%
  select(genome_id, shell_cluster, raw_cluster)
write_csv(cluster_assignments, file.path(OUT_DIR, paste0("03_shell_gene_cluster_assignments_k", SELECTED_K, ".csv")))
cat("Cluster sizes (k=", SELECTED_K, "):\n", sep = "")
print(table(cluster_assignments$shell_cluster))

# Master table
master <- cluster_assignments %>%
  left_join(vfdb_binary, by = "genome_id") %>%
  left_join(meta, by = "genome_id") %>%
  mutate(clinical_binary = case_when(niche == "clinical" ~ 1, niche == "non-clinical" ~ 0, TRUE ~ NA_real_),
         shell_cluster = factor(shell_cluster))
write_csv(master, file.path(OUT_DIR, "04_master_shell_cluster_metadata_VFDB_table.csv"))
cat("Master table genomes:", nrow(master), "\n")

# Cluster summary
metrics <- c("total_vf", vf_genes)
cluster_summary <- master %>% group_by(shell_cluster) %>%
  summarise(n = n(), mean_total_vf = mean(total_vf, na.rm = TRUE),
            median_total_vf = median(total_vf, na.rm = TRUE),
            sd_total_vf = sd(total_vf, na.rm = TRUE),
            pct_clinical = mean(clinical_binary, na.rm = TRUE) * 100,
            n_clinical = sum(niche == "clinical", na.rm = TRUE),
            n_nonclinical = sum(niche == "non-clinical", na.rm = TRUE),
            mean_year = mean(year, na.rm = TRUE),
            min_year = min(year, na.rm = TRUE),
            max_year = max(year, na.rm = TRUE),
            top_source_niches = top_values(source_niche),
            top_countries = top_values(country),
            .groups = "drop") %>%
  arrange(desc(mean_total_vf))
write_csv(cluster_summary, file.path(OUT_DIR, "05_shell_cluster_VFDB_metadata_summary.csv"))
cat("Done. Output in:", OUT_DIR, "\n")
