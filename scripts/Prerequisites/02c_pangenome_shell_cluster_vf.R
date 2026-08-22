#!/usr/bin/env Rscript
# ==============================================================================
# Prerequisites / 02c_pangenome_shell_cluster_vf.R
# ==============================================================================
# PROVENANCE: verbatim relocation of the original scripts/.../02c_pangenome_shell_cluster_vf.R.
# Shell-gene clustering (silhouette + k-selection) on the VirulenceFinder/ResFinder pangenome (ST10).
#
# ANALYSIS IS UNCHANGED -- this is the original script's body, unmodified.
# ==============================================================================

#!/usr/bin/env Rscript
#
# 01b_shell_cluster_vf.R  —  ST{ST} shell-gene clustering + VirulenceFinder
#
# Uses the common pangenome Rtab + shell partition, VirulenceFinder binary matrix,
# and ST-matched metadata.  Produces cluster assignments, master table, and
# all downstream validation outputs (temporal, clinical, decomposition, verdict).
#
# ST is set via config$TARGET_ST (default ST69, override with TARGET_ST=ST10).
#
source("config.R")

packages <- c("tidyverse", "readxl", "vegan", "cluster", "broom", "writexl", "ggplot2", "scales")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(vegan); library(cluster)
  library(broom); library(writexl); library(ggplot2); library(scales)
})

OUT_DIR <- config$VF_DIR
ST <- config$TARGET_ST
K_RANGE <- 4:config$CLUSTER_MAX_K
MIN_YEAR_N <- config$MIN_YEAR_N
MAX_YEAR <- config$MAX_YEAR
EARLY_YEARS <- config$EARLY_YEARS
LATE_YEARS <- config$LATE_YEARS
PRIMARY_HCLUST_METHOD <- "ward.D2"
SENSITIVITY_HCLUST_METHOD <- "average"

# ---- Verify inputs ----
for (f in c(config$GENE_PA_AB, config$SHELL_FILE, config$VF_BINARY,
            config$VF_BURDEN, config$METADATA_FILE))
  if (!file.exists(f)) stop("Missing: ", f)

# ---- Helper functions ----
safe_spearman <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
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
  for (p in patterns) { h <- grep(p, cn, ignore.case = TRUE, value = TRUE); if (length(h)) return(h[1]) }
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
read_shell_file <- function(path, rtab_families) {
  lines <- readLines(path, warn = FALSE); lines <- trimws(lines); lines <- lines[lines != ""]
  tokens <- unlist(strsplit(lines, "[,;\\t ]+")); tokens <- trimws(tokens); tokens <- tokens[tokens != ""]
  tokens <- tokens[!grepl("^(shell|family|families|gene|genes|partition)$", tokens, ignore.case = TRUE)]
  shell_fams <- intersect(tokens, rtab_families)
  tibble(family_id = shell_fams, partition = "shell") %>% distinct(family_id, .keep_all = TRUE)
}

# 1. Read pangenome matrix —---
rtab <- read_tsv(config$GENE_PA_AB, show_col_types = FALSE, progress = FALSE, guess_max = 10000)
if (ncol(rtab) < 3) stop("Rtab too few columns")
colnames(rtab)[1] <- "family_id"
rtab <- rtab %>% mutate(family_id = as.character(family_id)) %>% distinct(family_id, .keep_all = TRUE)
genome_cols_all <- setdiff(colnames(rtab), "family_id")
rtab[genome_cols_all] <- lapply(rtab[genome_cols_all], function(x) {
  x <- suppressWarnings(as.integer(as.character(x))); x[is.na(x)] <- 0L; as.integer(x > 0)
})
cat(sprintf("Rtab: %d families, %d genomes\n", nrow(rtab), length(genome_cols_all)))

# 2. Shell partition —---
partition_tbl <- read_shell_file(config$SHELL_FILE, rtab$family_id)
cat("Shell families matched:", nrow(partition_tbl), "\n")
if (nrow(partition_tbl) < 100) stop("Too few shell families")
write_csv(partition_tbl, file.path(OUT_DIR, "00_shell_families_used.csv"))

# 3. VirulenceFinder binary —---
vf_all <- read.delim(config$VF_BINARY, header = TRUE, sep = "\t", check.names = FALSE,
                     comment.char = "", na.strings = "", stringsAsFactors = FALSE)
vf_st <- vf_all %>% filter(st == ST)
if (nrow(vf_st) == 0) stop("No ", ST, " in VF binary matrix")
vf_st <- vf_st %>% mutate(genome_id = as.character(genome), .after = "st") %>% select(-genome)
vf_genes <- colnames(vf_st)[3:ncol(vf_st)]
vf_binary <- vf_st
vf_binary[, vf_genes] <- lapply(vf_binary[, vf_genes], function(x)
  as.integer(!is.na(x) & x != "." & x != "0" & x != ""))
vf_binary$total_vf <- rowSums(vf_binary[, vf_genes, drop = FALSE], na.rm = TRUE)

vf_burden <- read.delim(config$VF_BURDEN, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(st == ST) %>% mutate(genome_id = as.character(genome))
cat(sprintf("VF: %d genomes, %d genes\n", nrow(vf_st), length(vf_genes)))

# 4. VF gene groups —---
vf_groups <- list(
  Yersiniabactin = c("fyuA", "irp2"),
  P_fimbriae = c("papC", grep("^papA", vf_genes, value = TRUE)),
  Capsule_kps = c("kpsE", grep("^kpsM", vf_genes, value = TRUE), "neuC"),
  iro_cluster = c("iroN"),
  Type1_fimbriae = c("fimH", "fimF41"),
  Siderophores = c("iucC", "iutA", "ireA", "sitA", "iha"),
  Toxins = c("sat", "senB", "hlyA", "hlyE", "hlyF", "cnf1", "cnf2",
             grep("^cdt", vf_genes, value = TRUE), "astA", "espP", "pic", "pet", "subA", "eatA"),
  Adhesins = c(grep("^afa", vf_genes, value = TRUE), grep("^sfa", vf_genes, value = TRUE),
               grep("^foc", vf_genes, value = TRUE), grep("^lng", vf_genes, value = TRUE),
               grep("^fed", vf_genes, value = TRUE), grep("^fan", vf_genes, value = TRUE),
               grep("^fae", vf_genes, value = TRUE), grep("^csg", vf_genes, value = TRUE),
               "lpfA", "tsh", "tia", "ecpD", "csnB"),
  Invasins = c("ibeA", "iss", "tcpC", "ompT"),
  Colicins = c("cib", "cma", "cvaC", "cba", "cia",
               grep("^colE", vf_genes, value = TRUE), grep("^mch", vf_genes, value = TRUE),
               "mcmA", "mcbA"),
  Protectins = c("traT", "terC", "nlpI", "yghJ"),
  T6SS = intersect(vf_genes, c("hcp", "vgrG", "tssA", "tssB", "tssC", "tssD",
                                "tssE", "tssF", "tssG", "tssH", "tssJ", "tssK", "tssL", "tssM")),
  Eae_nle = intersect(vf_genes, c("eae", "espA", "espB", "espF", "espJ", "nleA", "nleB", "nleC",
                                   "tir", "tccP", "cif"))
)
vf_groups <- lapply(vf_groups, function(g) intersect(g, vf_genes))
vf_groups <- vf_groups[sapply(vf_groups, length) > 0]

vf_scores <- vf_binary %>% select(genome_id, total_vf)
for (g in names(vf_groups)) vf_scores[[g]] <- score_gene_group(vf_binary, vf_groups[[g]])

# kps binary flag (>97% expected for ST69)
kps_genes <- intersect(
  c("kpsE", grep("^kpsM", vf_genes, value = TRUE), "neuC"), colnames(vf_binary))
vf_scores$has_kps <- as.integer(rowSums(vf_binary[, kps_genes, drop = FALSE]) > 0)
cat(sprintf("kps prevalence: %.1f%%\n", mean(vf_scores$has_kps) * 100))

# 5. Metadata —---
meta_raw <- read_excel(config$METADATA_FILE, col_names = TRUE)
name_col <- find_col(meta_raw, c("^Name$", "genome", "strain", "isolate", "assembly", "sample"), 2)
source_col <- find_col(meta_raw, c("Source.*Niche", "source_niche", "^source$", "niche", "host"), 4)
year_col <- find_col(meta_raw, c("Collection.*Year", "^year$", "collection_year", "date"), 7)
continent_col <- find_col(meta_raw, c("^Continent$", "continent"), 11)
country_col <- find_col(meta_raw, c("^Country$", "country"), 12)
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

# 6. Intersect —---
rtab_genomes <- setdiff(colnames(rtab), "family_id")
common <- Reduce(intersect, list(rtab_genomes, vf_scores$genome_id, meta$genome_id))
cat("Shared genomes before filter:", length(common), "\n")
if (length(common) < 100) stop("Too few shared genomes")
meta_common <- meta %>% filter(genome_id %in% common, !is.na(year), year <= MAX_YEAR)
valid_years <- meta_common %>% count(year, name = "n") %>% filter(n >= MIN_YEAR_N) %>% pull(year)
common <- meta_common %>% filter(year %in% valid_years) %>% pull(genome_id)
cat(sprintf("Shared after filter: %d  years: %s\n", length(common),
            paste(sort(valid_years), collapse = ",")))
if (length(common) < 100) stop("Too few after filtering")

# 7. Shell matrix —---
rtab_shell_common <- rtab %>% filter(family_id %in% partition_tbl$family_id) %>%
  select(family_id, all_of(common))
shell_matrix <- rtab_shell_common %>% select(-family_id)
shell_matrix[] <- lapply(shell_matrix, function(x) {
  x <- suppressWarnings(as.integer(as.character(x))); x[is.na(x)] <- 0L; as.integer(x > 0)
})
mat_fam <- as.matrix(shell_matrix); rownames(mat_fam) <- rtab_shell_common$family_id
keep <- rowSums(mat_fam) > 0 & rowSums(mat_fam) < ncol(mat_fam)
mat_fam <- mat_fam[keep, , drop = FALSE]
mat_gen <- t(mat_fam)
cat("Variable shell families:", ncol(mat_gen), "\n")

write_csv(tibble(stat = c("shared_genomes", "variable_shell_families", "years", "shell_source"),
                 val = c(length(common), ncol(mat_gen),
                         paste(sort(valid_years), collapse = ";"), "PPanGGOLiN shell.txt")),
          file.path(OUT_DIR, "01_input_summary.csv"))

# 8. Jaccard + clustering —---
cat("Jaccard distance...\n")
jaccard_dist <- vegdist(mat_gen, method = "jaccard", binary = TRUE)
hc_primary <- hclust(jaccard_dist, method = PRIMARY_HCLUST_METHOD)
hc_sens <- hclust(jaccard_dist, method = SENSITIVITY_HCLUST_METHOD)

sil_prim <- map_dfr(K_RANGE, function(k) {
  cl <- cutree(hc_primary, k = k)
  sil <- silhouette(cl, jaccard_dist)
  tibble(method = PRIMARY_HCLUST_METHOD, k = k, avg_sil = mean(sil[, "sil_width"]))
})
sil_sens <- map_dfr(K_RANGE, function(k) {
  cl <- cutree(hc_sens, k = k)
  sil <- silhouette(cl, jaccard_dist)
  tibble(method = SENSITIVITY_HCLUST_METHOD, k = k, avg_sil = mean(sil[, "sil_width"]))
})
sil_results <- bind_rows(sil_prim, sil_sens)
write_csv(sil_results, file.path(OUT_DIR, "02_silhouette_k4_to_k10.csv"))

best_k <- sil_prim$k[which.max(sil_prim$avg_sil)]
best_sil <- max(sil_prim$avg_sil, na.rm = TRUE)
k4_sil <- sil_prim %>% filter(k == 4) %>% pull(avg_sil)
if (length(k4_sil) == 0) k4_sil <- NA_real_
k4_ok <- !is.na(k4_sil) && (best_k == 4 || (best_sil - k4_sil) <= 0.02)
SELECTED_K <- if (k4_ok) 4 else best_k
cat(sprintf("Best k=%d (sil=%.4f)  k=4 sil=%.4f  selected k=%d\n", best_k, best_sil, k4_sil, SELECTED_K))

write_csv(tibble(best_k = best_k, best_sil = best_sil, k4_sil = k4_sil,
                 k4_supported = k4_ok, selected_k = SELECTED_K),
          file.path(OUT_DIR, "02B_k_selection_summary.csv"))

raw_cl <- cutree(hc_primary, k = SELECTED_K)
if (is.null(names(raw_cl))) names(raw_cl) <- rownames(mat_gen)
assignments <- tibble(genome_id = names(raw_cl), raw_cluster = as.integer(raw_cl))
cl_key <- assignments %>% count(raw_cluster) %>%
  arrange(desc(n), raw_cluster) %>% mutate(shell_cluster = paste0("Cluster_", row_number()))
assignments <- assignments %>% left_join(cl_key, by = "raw_cluster") %>%
  select(genome_id, shell_cluster, raw_cluster)
write_csv(assignments, file.path(OUT_DIR, sprintf("03_shell_gene_cluster_assignments_k%d.csv", SELECTED_K)))
cat("Cluster sizes:\n"); print(table(assignments$shell_cluster))

# 9. Master table —---
master <- assignments %>%
  left_join(vf_scores, by = "genome_id") %>%
  left_join(meta, by = "genome_id") %>%
  mutate(clinical_binary = case_when(niche == "clinical" ~ 1, niche == "non-clinical" ~ 0, TRUE ~ NA_real_),
         shell_cluster = factor(shell_cluster),
         has_kps = replace_na(has_kps, 0L))
write_csv(master, file.path(OUT_DIR, "04_master_shell_cluster_metadata_VF_table.csv"))
cat("Master:", nrow(master), "genomes\n")

# 10. Summary tables —---
metrics <- c("total_vf", names(vf_groups))

cl_summary <- master %>% group_by(shell_cluster) %>%
  summarise(n = n(), mean_vf = mean(total_vf, na.rm = TRUE),
            median_vf = median(total_vf, na.rm = TRUE), sd_vf = sd(total_vf, na.rm = TRUE),
            pct_clinical = mean(clinical_binary, na.rm = TRUE) * 100,
            mean_year = mean(year, na.rm = TRUE), .groups = "drop",
            pct_kps = mean(has_kps, na.rm = TRUE) * 100) %>%
  arrange(desc(mean_vf))
for (g in names(vf_groups)) {
  cl_summary[[paste0("mean_", g)]] <- master %>% group_by(shell_cluster) %>%
    summarise(m = mean(.data[[g]], na.rm = TRUE), .groups = "drop") %>% pull(m)
}
write_csv(cl_summary, file.path(OUT_DIR, "05_shell_cluster_VF_metadata_summary.csv"))

# 11. Global tests —---
global_tests <- map_dfr(metrics, function(m) {
  kt <- safe_kruskal(master[[m]], master$shell_cluster)
  tibble(metric = m, stat = ifelse(is.null(kt), NA_real_, kt$statistic),
         p = ifelse(is.null(kt), NA_real_, kt$p.value))
}) %>% mutate(p_adj = p.adjust(p, "BH")) %>% arrange(p_adj)
write_csv(global_tests, file.path(OUT_DIR, "06_global_VF_tests_across_shell_clusters.csv"))

# 12. Temporal —---
sp <- safe_spearman(master$year, master$total_vf)
lm1 <- lm(total_vf ~ year, data = master)
write_csv(tibble(test = c("spearman", "lm"), est = c(sp$rho, coef(lm1)[["year"]]),
                  p = c(sp$p, summary(lm1)$coef["year", "Pr(>|t|)"])),
          file.path(OUT_DIR, "07_global_temporal_VF_escalation.csv"))

yr_cl <- master %>% count(year, shell_cluster) %>%
  group_by(year) %>% mutate(total_n = sum(n), prop = n / total_n) %>% ungroup()
write_csv(yr_cl, file.path(OUT_DIR, "08_yearly_shell_cluster_proportions.csv"))

cl_prop_trends <- yr_cl %>% group_by(shell_cluster) %>%
  summarise(n_years = n_distinct(year), first_prop = prop[which.min(year)][1],
            last_prop = prop[which.max(year)][1], delta = last_prop - first_prop,
            rho = safe_spearman(year, prop)$rho, p = safe_spearman(year, prop)$p, .groups = "drop") %>%
  mutate(p_adj = p.adjust(p, "BH"))
write_csv(cl_prop_trends, file.path(OUT_DIR, "09_shell_cluster_proportion_temporal_trends.csv"))

within_trends <- map_dfr(levels(master$shell_cluster), function(cl) {
  d <- master %>% filter(shell_cluster == cl)
  map_dfr(metrics, function(m) {
    dd <- d %>% filter(!is.na(year), !is.na(.data[[m]]))
    if (nrow(dd) < 10 || n_distinct(dd$year) < 3)
      return(tibble(cluster = cl, metric = m, n = nrow(dd), rho = NA_real_, p = NA_real_))
    sp2 <- safe_spearman(dd$year, dd[[m]])
    tibble(cluster = cl, metric = m, n = nrow(dd), rho = sp2$rho, p = sp2$p)
  })
}) %>% group_by(metric) %>% mutate(p_adj = p.adjust(p, "BH")) %>% ungroup()
write_csv(within_trends, file.path(OUT_DIR, "10_within_shell_cluster_VF_temporal_trends.csv"))

# 13. Models —---
country_counts <- table(master$country)
major_cty <- names(country_counts[country_counts >= 10])
df_m <- master %>% mutate(cty = ifelse(as.character(country) %in% major_cty, as.character(country), "Other")) %>%
  filter(!is.na(total_vf), !is.na(year), !is.na(shell_cluster), !is.na(niche))
m1 <- lm(total_vf ~ year, data = df_m)
m2 <- lm(total_vf ~ year + shell_cluster, data = df_m)
m3 <- lm(total_vf ~ year + shell_cluster + niche + cty, data = df_m)
write_csv(bind_rows(broom::tidy(m1) %>% mutate(model = "m1"),
                     broom::tidy(m2) %>% mutate(model = "m2"),
                     broom::tidy(m3) %>% mutate(model = "m3")),
          file.path(OUT_DIR, "11_model_coefficients.csv"))
write_csv(tibble(model = c("m1","m2","m3"), AIC = c(AIC(m1),AIC(m2),AIC(m3)),
                  BIC = c(BIC(m1),BIC(m2),BIC(m3)),
                  adjR2 = c(summary(m1)$adj.r.squared, summary(m2)$adj.r.squared,
                            summary(m3)$adj.r.squared)),
          file.path(OUT_DIR, "12_model_fit_stats.csv"))

# 14. Early-late decomposition —---
pf <- master %>% filter(year %in% c(EARLY_YEARS, LATE_YEARS)) %>%
  mutate(period = case_when(year %in% EARLY_YEARS ~ "early", year %in% LATE_YEARS ~ "late"))
p_sum <- pf %>% group_by(period, shell_cluster) %>%
  summarise(n = n(), mv = mean(total_vf, na.rm = TRUE), .groups = "drop") %>%
  group_by(period) %>% mutate(prop = n / sum(n)) %>% ungroup()
all_pc <- expand_grid(period = c("early", "late"), shell_cluster = levels(master$shell_cluster))
overall <- master %>% group_by(shell_cluster) %>% summarise(gmv = mean(total_vf, na.rm = TRUE), .groups = "drop")
p_full <- all_pc %>% left_join(p_sum, by = c("period", "shell_cluster")) %>% left_join(overall, by = "shell_cluster") %>%
  mutate(n = replace_na(n, 0L), prop = ifelse(is.na(prop), 0, prop), mv = ifelse(is.na(mv), gmv, mv))
early <- p_full %>% filter(period == "early") %>% select(shell_cluster, pe = prop, me = mv)
late <- p_full %>% filter(period == "late") %>% select(shell_cluster, pl = prop, ml = mv)
decomp <- early %>% inner_join(late, by = "shell_cluster") %>%
  mutate(dp = pl - pe, dm = ml - me, pb = (pe + pl)/2, mb = (me + ml)/2,
         demog = dp * mb, within = pb * dm, interact = dp * dm / 2)
decomp_summary <- with(decomp, tibble(
  early_vf = sum(pe * me), late_vf = sum(pl * ml),
  delta = sum(pl * ml) - sum(pe * me),
  demog = sum(demog), within = sum(within), interact = sum(interact),
  demog_pct = sum(demog) / (sum(pl * ml) - sum(pe * me)) * 100,
  within_pct = sum(within) / (sum(pl * ml) - sum(pe * me)) * 100,
  interact_pct = sum(interact) / (sum(pl * ml) - sum(pe * me)) * 100))
write_csv(p_full, file.path(OUT_DIR, "15_period_shell_cluster_summary.csv"))
write_csv(decomp, file.path(OUT_DIR, "16_early_late_decomposition_by_shell_cluster.csv"))
write_csv(decomp_summary, file.path(OUT_DIR, "17_early_late_decomposition_summary.csv"))

# 15. Weighted contribution —---
yr_cl_vf <- master %>% group_by(year, shell_cluster) %>%
  summarise(n = n(), mv = mean(total_vf, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>% mutate(prop = n / sum(n), wc = prop * mv) %>% ungroup()
wc_trends <- yr_cl_vf %>% group_by(shell_cluster) %>%
  summarise(first = wc[which.min(year)][1], last = wc[which.max(year)][1],
            delta = last - first, rho = safe_spearman(year, wc)$rho, p = safe_spearman(year, wc)$p,
            .groups = "drop") %>% mutate(p_adj = p.adjust(p, "BH"))
write_csv(yr_cl_vf, file.path(OUT_DIR, "18_yearly_weighted_VF_contributions.csv"))
write_csv(wc_trends, file.path(OUT_DIR, "19_weighted_contribution_trends.csv"))

# 16. Clinical tests —---
niche_tests <- map_dfr(metrics, function(m) {
  d <- master %>% filter(niche %in% c("clinical", "non-clinical"), !is.na(.data[[m]]))
  wt <- safe_wilcox(d %>% filter(niche == "clinical") %>% pull(m),
                    d %>% filter(niche == "non-clinical") %>% pull(m))
  tibble(metric = m, clin_mean = mean(d[[m]][d$niche == "clinical"], na.rm = TRUE),
         nonclin_mean = mean(d[[m]][d$niche == "non-clinical"], na.rm = TRUE),
         delta = clin_mean - nonclin_mean, p = ifelse(is.null(wt), NA_real_, wt$p.value))
}) %>% mutate(p_adj = p.adjust(p, "BH")) %>% arrange(p_adj)
write_csv(niche_tests, file.path(OUT_DIR, "20_clinical_vs_nonclinical_VF_tests.csv"))

# 17. Verdict —---
verdict <- tibble(
  criterion = c("Shell clustering silhouette", "Global VF ↑ over time",
                "VF differs across clusters", "High-VF cluster expands",
                "Within-cluster VF ↑", "Weighted contribution ↑",
                "Demographic component +", "Within-cluster component +",
                "FINAL: sublineage-driven escalation"),
  status = c(if(k4_ok) "Supported" else "Not supported",
             if(sp$p < .05 && sp$rho > 0) "Supported" else "Not supported",
             if(!is.null(global_tests$p[1]) && global_tests$p[1] < .05) "Supported" else "Not supported",
             if(any(cl_prop_trends$p_adj < .05 & cl_prop_trends$delta > 0 &
                     cl_summary$mean_vf[match(cl_prop_trends$shell_cluster, cl_summary$shell_cluster)] >
                     mean(master$total_vf, na.rm = TRUE))) "Supported" else "Not strong",
             if(any(within_trends$p_adj < .05 & within_trends$rho > 0)) "Supported" else "Not strong",
             if(any(wc_trends$p_adj < .05 & wc_trends$delta > 0)) "Supported" else "Not strong",
             if(decomp_summary$demog > 0) "Supported" else "Not supported",
             if(decomp_summary$within > 0) "Supported" else "Not supported",
             if(k4_ok && sp$p < .05 && sp$rho > 0) "Supported with caveats" else "Not fully supported"))
write_csv(verdict, file.path(OUT_DIR, "21_VALIDATION_VERDICT.csv"))
print(verdict)

# 18. Figures —---
p_sil <- sil_results %>%
  ggplot(aes(x = k, y = avg_sil, color = method)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  geom_vline(xintercept = SELECTED_K, linetype = "dashed") +
  theme_classic() + labs(title = paste(ST, "shell clustering silhouette"), x = "k", y = "Avg silhouette")
ggsave(file.path(OUT_DIR, "Figure_1_silhouette.pdf"), p_sil, width = 7, height = 5)
ggsave(file.path(OUT_DIR, "Figure_1_silhouette.png"), p_sil, width = 7, height = 5, dpi = 300)

p_glob <- master %>% group_by(year) %>%
  summarise(n = n(), m = mean(total_vf, na.rm = TRUE), se = sd(total_vf, na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
  ggplot(aes(x = year, y = m)) + geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = .25) +
  geom_line() + geom_point(aes(size = n)) + theme_classic() +
  labs(title = paste(ST, "VF burden over time"), x = "Year", y = "Mean VF")
ggsave(file.path(OUT_DIR, "Figure_2_global_temporal_trend.pdf"), p_glob, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "Figure_2_global_temporal_trend.png"), p_glob, width = 8, height = 6, dpi = 300)

p_prop <- yr_cl %>% ggplot(aes(x = year, y = prop, color = shell_cluster)) +
  geom_line() + geom_point() + theme_classic() + scale_y_continuous(labels = percent_format()) +
  labs(title = paste(ST, "cluster proportions"), x = "Year", y = "Proportion")
ggsave(file.path(OUT_DIR, "Figure_3_cluster_proportions.pdf"), p_prop, width = 9, height = 6)
ggsave(file.path(OUT_DIR, "Figure_3_cluster_proportions.png"), p_prop, width = 9, height = 6, dpi = 300)

decomp_long <- decomp %>% select(shell_cluster, demog, within, interact) %>%
  pivot_longer(-shell_cluster, names_to = "comp", values_to = "val")
p_dec <- decomp_long %>% ggplot(aes(x = reorder(shell_cluster, val), y = val, fill = comp)) +
  geom_col(position = "stack") + coord_flip() + theme_classic() +
  labs(title = paste(ST, "early-late VF decomposition"), x = "Cluster", y = "Contribution")
ggsave(file.path(OUT_DIR, "Figure_4_decomposition.pdf"), p_dec, width = 9, height = 6)
ggsave(file.path(OUT_DIR, "Figure_4_decomposition.png"), p_dec, width = 9, height = 6, dpi = 300)

# 19. Excel —---
write_xlsx(list(master = master, cluster_summary = cl_summary, global_tests = global_tests,
                 temporal = decomp_summary, verdict = verdict),
           path = file.path(OUT_DIR, paste0(ST, "_virulencefinder_results.xlsx")))

cat("\nDONE. Output in:", OUT_DIR, "\n")
