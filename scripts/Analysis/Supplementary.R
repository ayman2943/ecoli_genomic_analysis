#!/usr/bin/env Rscript
# ==============================================================================
# Analysis / Supplementary.R
# ==============================================================================
# Reproduces every supplementary analysis, table, and figure output that is
# NOT one of the 8 main published figures (those live in Figure1.R-Figure8.R).
#
# PROVENANCE: this file is a straight concatenation of 17 originally-separate
# scripts (16 supplementary analyses, plus a 17th block holding the FULL,
# unmodified `12_fig02-05_ST69_analysis.R` -- see that block's own header
# comment for why), each wrapped in its own local({ ... }) block. This
# mirrors exactly how run_pipeline.R already isolates every script from
# every other one via source(path, local = TRUE) -- each block gets its own
# environment, so variables like `master`, `OUT`, `df` etc. in one original
# script cannot collide with the same names in another. NO analysis logic
# has been changed; every block below is the original script's body,
# verbatim. Each block is additionally wrapped in tryCatch() so that one
# block's failure (e.g. a missing optional input file) does not prevent the
# remaining blocks from running -- mirroring the per-script resilience
# run_pipeline.R already provides at the whole-script level.
#
# ORDER MATTERS for the capsule/K-type chain: 09b_capsule_classification.R
# WRITES capsule_classification.xlsx, and 09c_capsule_comparison.R,
# 09d_k_type_analysis.R, and 10e_allelic_conversion.R each READ-THEN-APPEND
# a new sheet to that same file. That original append dependency is
# preserved here by keeping those four blocks in their original relative
# order (09b -> 09c -> 09d -> ... -> 10e). Do not reorder those four blocks.
#
# 22_fig_supplementary_combined.R runs LAST. It is self-contained (it
# re-loads its own source data rather than depending on the blocks above),
# but is kept last by convention since it produces the final combined
# Supplementary_Figures_Combined.pdf that ties everything together -- run
# after the individual blocks above so every individual xlsx/csv/png output
# is generated too, not just the combined PDF.
#
# Included but NOT part of the manuscript's main 8 figures or a required
# supplementary output (kept for completeness / to honor "keep my analysis
# intact"): 14_fig06_gene_trajectories.R (exploratory 3-panel gene
# prevalence figure) and 15_fig07_tree_parsimony.R (a preliminary core-tree
# figure that is NOT the source of the submitted Figure 3 -- that figure was
# made in iTOL; see Figure3.R).
#
# Two further scripts that existed in the original scripts/09_figures/
# folder -- 19_fig_composite.R and 99_ALL_FIGURES_COMBINED.R -- have been
# LEFT OUT of this reorganised repo. Both are pure post-processing utilities
# that read already-generated PNGs from disk and re-paste them into
# Word-page-sized composites; they perform no analysis of their own and do
# not correspond 1:1 to any of the 8 figures or to this supplementary file.
# They remain available in git history if needed.
# ==============================================================================


################################################################################
# BLOCK 1/17: 07a_st10_decomposition.R
# ST10 pangenome decomposition
# (verbatim body of original scripts/.../07a_st10_decomposition.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# Decompose ST10 Cluster_3 VF increase: composition vs within-stratum
suppressPackageStartupMessages(library(tidyverse))
source("config.R"); config$TARGET_ST <- "ST10"; source("config.R", local=TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

master <- read_csv(file.path(config$OUTPUT_DIR, "ST10", "virulencefinder_validation",
  "04_master_shell_cluster_metadata_VF_table.csv"), show_col_types=FALSE) %>%
  mutate(genome_id=as.character(genome_id), year=as.integer(year),
    clinical_binary=as.integer(clinical_binary),
    source_niche=as.character(source_niche), country=as.character(country)) %>%
  filter(!is.na(clinical_binary),!is.na(year),!is.na(shell_cluster))

vfb <- read.delim(config$VF_BINARY, header=TRUE, sep="\t", check.names=FALSE,
  na.strings="", stringsAsFactors=FALSE) %>%
  filter(st=="ST10") %>% mutate(genome_id=make_genome_id(genome)) %>% select(-st,-genome)
for (g in setdiff(colnames(vfb),"genome_id"))
  vfb[[g]] <- as.integer(!is.na(vfb[[g]])&vfb[[g]]!="."&vfb[[g]]!="0"&vfb[[g]]!="")
vfb <- vfb %>% mutate(total=rowSums(across(all_of(setdiff(colnames(vfb),"genome_id")))))

cd <- master %>% inner_join(vfb, by="genome_id") %>%
  filter(shell_cluster=="Cluster_3",!is.na(source_niche),!is.na(country))

# Create stratum = country | niche (collapse small groups)
cd <- cd %>% mutate(stratum = paste0(country, " | ", source_niche))

# 1. Yearly composition (proportion of each stratum)
comp <- cd %>% count(year, stratum) %>% group_by(year) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup()

# 2. Each stratum's overall mean VF (fixed reference)
stratum_means <- cd %>% group_by(stratum) %>%
  summarise(mean_vf = mean(total), total_n = n(), .groups = "drop")

# Join composition with stratum means
comp <- comp %>% left_join(stratum_means, by = "stratum")

# 3. Standardization: expected VF if only composition changed, but each stratum's VF stayed fixed
# For each stratum, compute overall mean VF (pooled across all years) — the "fixed" within value
stratum_fixed <- cd %>% group_by(stratum) %>%
  summarise(fixed_vf = mean(total), total_n = n(), .groups = "drop")

# Expected yearly VF = Σ(stratum_proportion_year × fixed_vf) for ALL strata in that year
expected <- comp %>% left_join(stratum_fixed, by = "stratum") %>%
  group_by(year) %>%
  summarise(expected = weighted.mean(fixed_vf, n), .groups = "drop")

# Actual yearly trend
actual <- cd %>% group_by(year) %>%
  summarise(actual = mean(total), n = n(), .groups = "drop")

compare <- actual %>% left_join(expected, by = "year")

# ===== Figure =====
# Panel A: stacked bar of composition
top_strata <- stratum_means %>% slice_max(total_n, n = 8) %>% pull(stratum)
cd_plot <- cd %>% mutate(stratum_plot = ifelse(stratum %in% top_strata, stratum, "Other"),
  stratum_plot = factor(stratum_plot, levels = c(sort(top_strata), "Other")))

p1 <- cd_plot %>% count(year, stratum_plot) %>% group_by(year) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup() %>%
  ggplot(aes(x = year, y = pct, fill = stratum_plot)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = NULL, y = "Proportion (%)", title = "Composition of ST10 Cluster_3 by year") +
  theme_minimal(base_size = 10) +
  theme(legend.text = element_text(size = 7))

# Panel B: within-stratum yearly VF (strata with >= 10 genomes)
big_strata <- stratum_means %>% filter(total_n >= 10) %>% pull(stratum)
pd <- cd %>% filter(stratum %in% big_strata) %>%
  group_by(year, stratum) %>%
  summarise(mean_vf = mean(total), se = sd(total)/sqrt(n()), n = n(), .groups = "drop")
pd_clean <- pd %>% filter(n >= 3)  # drop years with <3 genomes for plotting

p2 <- pd_clean %>%
  ggplot(aes(x = year, y = mean_vf, color = stratum)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, alpha = 0.5) +
  scale_color_brewer(palette = "Set2", name = NULL) +
  labs(x = NULL, y = "Mean VF burden",
       title = "Within-stratum VF trends (strata with >= 10 genomes)") +
  theme_minimal(base_size = 10) +
  theme(legend.text = element_text(size = 7))

# Panel C: actual vs expected (composition-only)
p3 <- compare %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = actual, color = "Actual"), linewidth = 1) +
  geom_point(aes(y = actual, color = "Actual"), size = 2.5) +
  geom_line(aes(y = expected, color = "Expected (composition only)"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = expected, color = "Expected (composition only)"),
             size = 2.5, shape = 1) +
  scale_color_manual(values = c("Actual" = "red",
    "Expected (composition only)" = "steelblue"), name = NULL) +
  labs(x = "Year", y = "Mean VF burden",
       title = "Actual vs composition-only expected VF") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")

# Assemble
library(patchwork)
fig <- (p1 / p2 / p3) + plot_annotation(
       title = "ST10 VF Cluster_3: compositional analysis",
       subtitle = "If actual and expected lines overlap, the increase is entirely compositional")

outdir <- file.path(config$OUTPUT_DIR, "sensitivity_analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(outdir, "decomposition_compositional.png"), fig, width = 9, height = 12, dpi = 300)
cat("Saved:", file.path(outdir, "decomposition_compositional.png"), "\n")

# Print summary
cat("\n=== Actual vs Standardized (composition-adjusted) ===\n")
print(compare, n = 20)

# LM on actual
m_actual <- lm(actual ~ year, data = compare)
cat(sprintf("\nActual trend: slope=%.4f  p=%.4f\n",
  coef(m_actual)["year"], coef(summary(m_actual))["year","Pr(>|t|)"]))

# LM on expected (should match actual if composition fully explains)
m_exp <- lm(expected ~ year, data = compare)
cat(sprintf("Expected (composition-only) trend: slope=%.4f  p=%.4f\n",
  coef(m_exp)["year"], coef(summary(m_exp))["year","Pr(>|t|)"]))
cat(sprintf("Actual - Expected difference in slope: %.4f\n",
  coef(m_actual)["year"] - coef(m_exp)["year"]))

})  # end block 1: 07a_st10_decomposition.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 2/17: 07b_st10_composition_drivers.R
# ST10 composition drivers
# (verbatim body of original scripts/.../07b_st10_composition_drivers.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# ST10 Cluster_3 composition driver analysis
# 1) Stacked bar: composition with mean VF per segment in legend
# 2) Clinical fraction by stratum (US Human vs other Human)
# 3) Gene prevalence: which genes are enriched in US Human vs UK Human vs rest
suppressPackageStartupMessages(library(tidyverse))
source("config.R"); config$TARGET_ST <- "ST10"; source("config.R", local=TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

master <- read_csv(file.path(config$OUTPUT_DIR, "ST10", "virulencefinder_validation",
  "04_master_shell_cluster_metadata_VF_table.csv"), show_col_types=FALSE) %>%
  mutate(genome_id=as.character(genome_id), year=as.integer(year),
    clinical_binary=as.integer(clinical_binary),
    source_niche=as.character(source_niche), country=as.character(country)) %>%
  filter(!is.na(clinical_binary),!is.na(year),!is.na(shell_cluster))

vfb <- read.delim(config$VF_BINARY, header=TRUE, sep="\t", check.names=FALSE,
  na.strings="", stringsAsFactors=FALSE) %>%
  filter(st=="ST10") %>% mutate(genome_id=make_genome_id(genome)) %>% select(-st,-genome)
genes_vf <- setdiff(colnames(vfb),"genome_id")
for (g in genes_vf)
  vfb[[g]] <- as.integer(!is.na(vfb[[g]])&vfb[[g]]!="."&vfb[[g]]!="0"&vfb[[g]]!="")

cd <- master %>% inner_join(vfb, by="genome_id") %>%
  filter(shell_cluster=="Cluster_3",!is.na(source_niche),!is.na(country))

cd <- cd %>% mutate(total = rowSums(across(all_of(genes_vf))))

cd <- cd %>% mutate(
  segment = case_when(
    country=="United States" & source_niche=="Human" ~ "US Human",
    country=="United Kingdom" & source_niche=="Human" ~ "UK Human",
    source_niche=="Human" ~ "Other Human",
    TRUE ~ source_niche
  ),
  segment = factor(segment, levels = c("US Human", "UK Human", "Other Human",
    "Environment", "Livestock", "Poultry", "Wild Animal", "Companion Animal", "Food"))
)

outdir <- file.path(config$OUTPUT_DIR, "sensitivity_analysis")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# ======================================================================
# 1. Stacked bar with mean VF in legend
# ======================================================================
segment_vf <- cd %>% group_by(segment) %>%
  summarise(mean_vf = mean(total), n = n(), .groups="drop") %>%
  mutate(label = sprintf("%s (VF=%.1f, n=%d)", segment, mean_vf, n))

label_map <- setNames(segment_vf$label, segment_vf$segment)

segment_colors <- c(
  "US Human" = "#E41A1C",
  "UK Human" = "#377EB8",
  "Other Human" = "#4DAF4A",
  "Environment" = "#984EA3",
  "Livestock" = "#FF7F00",
  "Poultry" = "#FFFF33",
  "Wild Animal" = "#A65628",
  "Companion Animal" = "#F781BF",
  "Food" = "#999999"
)

p_comp <- cd %>% count(year, segment) %>% group_by(year) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup() %>%
  ggplot(aes(x = year, y = pct, fill = segment)) +
  geom_col() +
  scale_fill_manual(values = segment_colors, labels = label_map, name = NULL) +
  labs(x = "Year", y = "Proportion (%)",
       title = "ST10 VF Cluster_3: composition shift drives VF increase",
       subtitle = "Segment mean VF burden shown in legend") +
  theme_minimal(base_size = 11) +
  theme(legend.text = element_text(size = 8),
        legend.position = "right",
        legend.key.size = unit(0.6, "cm"))

ggsave(file.path(outdir, "composition_bar_with_vf.png"), p_comp, width = 10, height = 6, dpi=300)
cat("Figure:", file.path(outdir, "composition_bar_with_vf.png"), "\n")

# ======================================================================
# 2. Clinical fraction by stratum
# ======================================================================
cat("\n=== Clinical fraction by segment ===\n")
clin <- cd %>% group_by(segment) %>%
  summarise(n=n(), n_clinical=sum(clinical_binary),
    pct_clinical=mean(clinical_binary)*100, mean_vf=mean(total), .groups="drop") %>%
  arrange(desc(pct_clinical))
print(clin, n=20)

# Clinical fraction over time for US Human vs Other Human
cat("\n=== Clinical fraction by year for US Human vs UK Human vs Other Human ===\n")
cd %>% filter(segment %in% c("US Human", "UK Human", "Other Human")) %>%
  group_by(year, segment) %>%
  summarise(n=n(), pct_clinical=mean(clinical_binary)*100, .groups="drop") %>%
  print(n=30)

# ======================================================================
# 3. Gene prevalence: US Human vs rest of Cluster_3
# ======================================================================
cat("\n=== Top genes enriched in US Human vs rest of Cluster_3 ===\n")
us_h <- cd %>% filter(segment == "US Human")
rest <- cd %>% filter(segment != "US Human")

gene_tests <- bind_rows(lapply(genes_vf, function(g) {
  n_us <- sum(us_h[[g]], na.rm=TRUE); n_us_tot <- sum(!is.na(us_h[[g]]))
  n_rest <- sum(rest[[g]], na.rm=TRUE); n_rest_tot <- sum(!is.na(rest[[g]]))
  prev_us <- n_us / n_us_tot * 100
  prev_rest <- n_rest / n_rest_tot * 100
  delta <- prev_us - prev_rest
  mat <- matrix(c(n_us, n_us_tot - n_us, n_rest, n_rest_tot - n_rest), nrow=2)
  ft <- tryCatch(fisher.test(mat), error=function(e) NULL)
  or <- if (!is.null(ft)) unname(ft$estimate) else NA
  p <- if (!is.null(ft)) ft$p.value else NA
  tibble(gene=g, n_us=n_us, n_rest=n_rest, prev_us=prev_us, prev_rest=prev_rest,
         delta_pp=delta, or=or, p=p)
})) %>% mutate(p_adj = p.adjust(pmax(p, 0, na.rm=TRUE), method="BH")) %>%
  arrange(p_adj)

gene_tests %>% slice_head(n=20) %>% print()

# Write full table
write_csv(gene_tests, file.path(outdir, "st10_cluster3_us_human_gene_enrichment.csv"))
cat("Gene table:", file.path(outdir, "st10_cluster3_us_human_gene_enrichment.csv"), "\n")

# Figure: top enriched genes
top_enr <- gene_tests %>% filter(!is.na(or), is.finite(or), p_adj < 0.05) %>%
  arrange(desc(delta_pp)) %>% slice_head(n=20) %>%
  mutate(gene = fct_reorder(gene, delta_pp))

if (nrow(top_enr) > 0) {
  p_genes <- top_enr %>%
    ggplot(aes(x = delta_pp, y = gene, fill = p_adj < 0.001)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "steelblue"), guide="none") +
    labs(x = "Prevalence difference (US Human - rest, pp)", y = NULL,
         title = "Top VF genes enriched in US Human vs rest of Cluster_3") +
    theme_minimal(base_size = 9)
  ggsave(file.path(outdir, "st10_cluster3_us_human_gene_enrichment.png"),
         p_genes, width = 8, height = 6, dpi=300)
  cat("Gene figure:", file.path(outdir, "st10_cluster3_us_human_gene_enrichment.png"), "\n")
}

# ======================================================================
# Also UK Human vs rest (for comparison)
# ======================================================================
cat("\n=== Top genes enriched in UK Human vs rest of Cluster_3 ===\n")
uk_h <- cd %>% filter(segment == "UK Human")
rest2 <- cd %>% filter(segment != "UK Human")

gene_tests_uk <- bind_rows(lapply(genes_vf, function(g) {
  n_uk <- sum(uk_h[[g]], na.rm=TRUE); n_uk_tot <- sum(!is.na(uk_h[[g]]))
  n_rest <- sum(rest2[[g]], na.rm=TRUE); n_rest_tot <- sum(!is.na(rest2[[g]]))
  prev_uk <- n_uk / n_uk_tot * 100
  prev_rest <- n_rest / n_rest_tot * 100
  delta <- prev_uk - prev_rest
  mat <- matrix(c(n_uk, n_uk_tot - n_uk, n_rest, n_rest_tot - n_rest), nrow=2)
  ft <- tryCatch(fisher.test(mat), error=function(e) NULL)
  or <- if (!is.null(ft)) unname(ft$estimate) else NA
  p <- if (!is.null(ft)) ft$p.value else NA
  tibble(gene=g, n_uk=n_uk, n_rest=n_rest, prev_uk=prev_uk, prev_rest=prev_rest,
         delta_pp=delta, or=or, p=p)
})) %>% mutate(p_adj = p.adjust(pmax(p, 0, na.rm=TRUE), method="BH")) %>%
  arrange(p_adj)

gene_tests_uk %>% slice_head(n=20) %>% print()
write_csv(gene_tests_uk, file.path(outdir, "st10_cluster3_uk_human_gene_enrichment.csv"))
cat("UK gene table:", file.path(outdir, "st10_cluster3_uk_human_gene_enrichment.csv"), "\n")

# ======================================================================
cat("\nDone.\n")

})  # end block 2: 07b_st10_composition_drivers.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 3/17: 08_clinical_enrichment_3panel.R
# Clinical enrichment (3-panel), QC/validation figure
# (verbatim body of original scripts/.../08_clinical_enrichment_3panel.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
#
# Clinical vs non-clinical gene enrichment — 3-panel comparison
# Replicates Figure 9A from 4_clincial.R
# Panel A: ST69 VFDB | Panel B: ST10 VirulenceFinder | Panel C: ST10 ResFinder
#
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(ggplot2); library(patchwork); library(scales)
})
source("config.R")

OUT <- file.path(config$OUTPUT_DIR, "figures_clinical_enrichment_3panel")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MIN_PRESENT_CLINICAL <- 5
MIN_PRESENT_NONCLINICAL <- 5
TOP_N <- 20

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

classify_system <- function(g, type = "vf") {
  vf_groups <- list(
    Yersiniabactin = c("irp1", "irp2", "fyuA/psn", "ybtA", "ybtE", "ybtP", "ybtQ", "ybtS", "ybtT", "ybtU", "ybtX"),
    P_fimbriae = c("papA", "papB", "papC", "papD", "papE", "papF", "papG", "papH", "papI", "papJ", "papK", "papX"),
    Capsule_kps = c("kpsC", "kpsD", "kpsE", "kpsF", "kpsM", "kpsS", "kpsT", "kpsU", "neuA", "neuB", "neuC", "neuD", "neuE"),
    iro_salmochelin = c("iroB", "iroC", "iroD", "iroE", "iroN"),
    Type1_fimbriae = c("fimA", "fimB", "fimC", "fimD", "fimE", "fimF", "fimG", "fimH", "fimI"),
    T6SS = c("tssA", "tssB", "tssC", "tssF", "tssG", "tssJ", "tssK", "tssL", "tssM", "hcp/tssD", "hcp1/tssD1", "vgrG/tssI"),
    Toxins = c("sat", "senB", "hlyA", "hlyB", "hlyC", "hlyD", "cnf1", "cdtA", "cdtB", "cdtC", "east1", "espP", "pic", "pet")
  )
  vf_anno <- purrr::imap_dfr(vf_groups, ~ tibble(gene = .x, system = .y)) %>% distinct(gene, .keep_all = TRUE)
  if (type == "vf") {
    m <- vf_anno$system[match(g, vf_anno$gene)]
    ifelse(is.na(m), "Other", m)
  } else {
    arg_classes <- c("Aminoglycosides" = "aac|aad|aph|ant|arm|rmt", "Beta-lactams" = "bla|oxa|nps", "Macrolides" = "erm|mef|mph|msr|lnu|lin",
                     "Tetracyclines" = "tet", "Sulfonamides" = "sul|dfr", "Phenicols" = "cat|cml",
                     "Quinolones" = "qnr|oqx|qep", "Polymyxins" = "mcr|arn", "Fosfomycin" = "fos")
    m <- names(arg_classes)[sapply(arg_classes, function(p) grepl(p, g, ignore.case = TRUE))]
    ifelse(length(m) > 0, m[1], "Other")
  }
}

# -------- Load metadata --------
cat("Loading metadata...\n")
meta69 <- read_xlsx(config$st_metadata("ST69")) %>%
  rename(genome_id = Name) %>%
  mutate(genome_id = make_genome_id(genome_id),
         source_niche = .data[["Source Niche"]],
         clinical = source_niche == "Human",
         year = as.integer(.data[["Collection Year"]])) %>%
  filter(!is.na(year), !is.na(source_niche))
meta10 <- read_xlsx(config$st_metadata("ST10")) %>%
  rename(genome_id = Name) %>%
  mutate(genome_id = make_genome_id(genome_id),
         source_niche = .data[["Source Niche"]],
         clinical = source_niche == "Human",
         year = as.integer(.data[["Collection Year"]])) %>%
  filter(!is.na(year), !is.na(source_niche))

cat("  ST69: ", sum(meta69$clinical), "clinical, ", sum(!meta69$clinical), "non-clinical\n")
cat("  ST10: ", sum(meta10$clinical), "clinical, ", sum(!meta10$clinical), "non-clinical\n")

# -------- Load gene data --------
load_genes <- function(type, st) {
  if (type == "vfdb") {
    raw <- read.delim(config$st_vfdb_summary(st), header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    colnames(raw)[1:2] <- c("genome_id", "num_found")
    raw$genome_id <- make_genome_id(str_remove(basename(raw$genome_id), "_vfdb\\.tsv$"))
    genes <- setdiff(colnames(raw), c("genome_id", "num_found"))
    long <- raw %>% mutate(across(all_of(genes), ~ as.integer(!is.na(.x) & .x != "."))) %>%
      select(genome_id, all_of(genes)) %>%
      pivot_longer(-genome_id, names_to = "gene", values_to = "present") %>%
      mutate(system = classify_system(gene, "vf"))
    return(long)
  } else if (type == "vf") {
    binary <- read.delim(config$VF_BINARY, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    long <- binary %>% filter(st == st) %>%
      mutate(genome_id = make_genome_id(genome)) %>%
      select(-st, -genome) %>%
      pivot_longer(-genome_id, names_to = "gene", values_to = "present") %>%
      mutate(present = as.integer(!is.na(present) & present != "" & present != "0"),
             system = classify_system(gene, "vf"))
    return(long)
  } else if (type == "arg") {
    binary <- read.delim(file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
                         header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    long <- binary %>% filter(st == st) %>%
      mutate(genome_id = make_genome_id(genome)) %>%
      select(-st, -genome) %>%
      pivot_longer(-genome_id, names_to = "gene", values_to = "present") %>%
      mutate(present = as.integer(!is.na(present) & present != "" & present != "0"),
             system = classify_system(gene, "arg"))
    return(long)
  }
}

cat("Loading gene data...\n")
vf69_long <- load_genes("vfdb", "ST69") %>% inner_join(meta69 %>% select(genome_id, clinical, year), by = "genome_id")
vf10_long <- load_genes("vf", "ST10") %>% inner_join(meta10 %>% select(genome_id, clinical, year), by = "genome_id")
resf10_long <- load_genes("arg", "ST10") %>% inner_join(meta10 %>% select(genome_id, clinical, year), by = "genome_id")
cat("  ST69 VFDB:", n_distinct(vf69_long$gene), "genes\n")
cat("  ST10 VF:", n_distinct(vf10_long$gene), "genes\n")
cat("  ST10 ResFinder:", n_distinct(resf10_long$gene), "genes\n")

# -------- Clinical enrichment analysis --------
enrichment_analysis <- function(long_df, case_name) {
  per_gene <- long_df %>%
    group_by(gene, system) %>%
    summarise(
      n_clinical = sum(clinical), n_non = sum(!clinical),
      present_clinical = sum(present[clinical]),
      present_non = sum(present[!clinical]),
      absent_clinical = n_clinical - present_clinical,
      absent_non = n_non - present_non,
      prev_clinical = present_clinical / n_clinical,
      prev_non = present_non / n_non,
      delta_pp = (prev_clinical - prev_non) * 100,
      .groups = "drop"
    ) %>% filter(present_clinical >= MIN_PRESENT_CLINICAL, present_non >= MIN_PRESENT_NONCLINICAL) %>%
    rowwise() %>%
    mutate(fr = {
      mat <- matrix(c(present_clinical, absent_clinical, present_non, absent_non), nrow = 2, byrow = TRUE)
      out <- tryCatch(fisher.test(mat), error = function(e) NULL)
      if (is.null(out)) list(tibble(fisher_or = NA_real_, fisher_p = NA_real_))
      else list(tibble(fisher_or = unname(out$estimate), fisher_p = out$p.value))
    }) %>% unnest(fr) %>% ungroup() %>%
    mutate(
      fisher_log10p = -log10(fisher_p),
      or_label = ifelse(fisher_p < 0.001, "***",
                        ifelse(fisher_p < 0.01, "**",
                               ifelse(fisher_p < 0.05, "*", "ns")))
    ) %>% arrange(desc(fisher_or)) %>%
    mutate(case = case_name)

  write_csv(per_gene, file.path(OUT, paste0(case_name, "_clinical_enrichment.csv")))
  cat("  ", case_name, ":", nrow(per_gene), "genes passed filters,", sum(per_gene$fisher_p < 0.05), "significant\n")
  per_gene
}

cat("\nAnalyzing clinical enrichment...\n")
e69 <- enrichment_analysis(vf69_long, "ST69_VFDB")
e10vf <- enrichment_analysis(vf10_long, "ST10_VF")
e10res <- enrichment_analysis(resf10_long, "ST10_ResFinder")

# -------- Plot: top enriched genes (dot plot, like Figure 9A) --------
plot_enrichment <- function(results, case_name, color, n_top = TOP_N) {
  enriched <- results %>% filter(fisher_p < 0.05) %>% slice_max(fisher_or, n = n_top, with_ties = FALSE)
  if (nrow(enriched) == 0) {
    df <- data.frame(x = 0, y = "", label = "No significantly enriched genes")
    return(ggplot(df, aes(x, y)) + geom_text(aes(label = label), size = 4) +
             labs(title = case_name) + theme_classic(base_size = 10) + xlim(-1, 1))
  }
  enriched <- enriched %>% mutate(
    gene_label = paste0(gene, " [", system, "]"),
    sig_label = ifelse(fisher_p < 0.001, "p<0.001",
                       ifelse(fisher_p < 0.01, sprintf("p=%.3f", fisher_p),
                              sprintf("p=%.3f", fisher_p)))
  )

  ggplot(enriched, aes(x = fisher_or, y = reorder(gene_label, fisher_or))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_point(aes(size = fisher_log10p), color = color, alpha = 0.85) +
    geom_text(aes(label = sig_label), hjust = -0.1, size = 2.5, color = "grey40") +
    scale_size_continuous(name = expression(-log[10](p))) +
    scale_x_log10() +
    labs(title = case_name,
         subtitle = paste0("Clinical (Human) vs non-clinical; Fisher exact test"),
         x = "Odds ratio (clinical / non-clinical)", y = NULL) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 7, color = "grey40"),
          legend.position = "bottom")
}

cat("\nPlotting...\n")
pA <- plot_enrichment(e69, "ST69 VFDB", "#2166AC")
pB <- plot_enrichment(e10vf, "ST10 VirulenceFinder", "#D73027")
pC <- plot_enrichment(e10res, "ST10 ResFinder", "#1B7837")

fig <- (pA | pB | pC) + plot_layout(ncol = 3) +
  plot_annotation(title = "Clinically enriched genes",
                  subtitle = "Top genes by odds ratio, Fisher exact test p < 0.05",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
                                plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5)))

ggsave(file.path(OUT, "Fig_clinical_enrichment_3panel.png"), fig, width = 18, height = 8, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig_clinical_enrichment_3panel.pdf"), fig, width = 18, height = 8, bg = "white")
cat("\nSaved: Fig_clinical_enrichment_3panel\n")

})  # end block 3: 08_clinical_enrichment_3panel.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 4/17: 09b_capsule_classification.R
# Capsule classification (writes capsule_classification.xlsx)
# (verbatim body of original scripts/.../09b_capsule_classification.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# Reviewer response: Capsule classification (G2 vs G3) on clinical enrichment + temporal increase
# Uses VFDB summary + master join matching the publication
suppressPackageStartupMessages({
  library(tidyverse); library(writexl); library(data.table); library(broom)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_capsule_classification")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
    grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
    TRUE ~ paste0("Escherichia_coli_", x))
}
to_bin_vfdb <- function(x) as.integer(sapply(strsplit(as.character(x), ";"),
  function(v) any(suppressWarnings(as.numeric(v) >= 90), na.rm = TRUE)))

# ---- 1. Load VFDB summary (same source as publication Figure_7) ----
cat("Loading VFDB summary...\n")
vfdb_raw <- read_tsv(config$st_vfdb_summary("ST69"), show_col_types = FALSE,
  progress = FALSE, col_types = cols(.default = "c")) %>%
  filter(str_detect(.data[["#FILE"]], "^ST69/")) %>%
  mutate(genome_id = make_genome_id(str_remove(.data[["#FILE"]], "^ST69/")),
    genome_id = str_remove(genome_id, "_vfdb\\.tsv$"))
genes <- setdiff(colnames(vfdb_raw), c("#FILE", "NUM_FOUND", "genome_id"))
bin <- vfdb_raw %>% select(genome_id, all_of(genes)) %>%
  mutate(across(all_of(genes), to_bin_vfdb))
cat("VFDB summary:", nrow(bin), "genomes\n")

# ---- 2. Load master table for metadata ----
cat("Loading master table...\n")
master <- fread(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv"),
  sep = ",", header = TRUE, data.table = FALSE, na.strings = "")
meta <- master %>% select(genome_id, shell_cluster, year, clinical_binary)
cat("Master table:", nrow(meta), "genomes\n")

# ---- 3. Join (same as publication figure_combined_summary.R) ----
df <- meta %>% inner_join(bin, by = "genome_id") %>%
  mutate(year = as.integer(year), clin = as.integer(clinical_binary)) %>%
  filter(!is.na(clin), !is.na(year), year >= 2016, year <= 2025)
cat("Joined:", nrow(df), "genomes\n")

# ---- 4. Capsule classification ----
kps_cols <- intersect(grep("^(kps|neu)", genes, value = TRUE), colnames(df))
G2_markers <- intersect(c("kpsC", "kpsS"), kps_cols)
G3_markers <- intersect(c("kpsE", "kpsM", "kpsT"), kps_cols)

df$capsule_type <- apply(df[, kps_cols, drop = FALSE], 1, function(r) {
  has_g2 <- any(r[names(r) %in% G2_markers] >= 1)
  has_g3 <- any(r[names(r) %in% G3_markers] >= 1)
  any_kps <- any(r >= 1)
  if (has_g2) return("G2")
  if (has_g3) return("G3")
  if (any_kps) return("Unclassified")
  return("No capsule")
})

caps_dist <- df %>% count(capsule_type) %>% mutate(pct = n / sum(n) * 100)
cat("\n=== Capsule Distribution ===\n")
print(caps_dist)

# ---- 5. Key genes for clinical enrichment + temporal increase ----
key_genes <- intersect(c("kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT","kpsU",
  "papX","papF","papB","papG","sat","hlyA"), genes)

# Clinical enrichment (delta prevalence + OR, matching Figure 7B)
cat("\n=== Clinical Enrichment (all genomes) ===\n")
clin_enrich <- map_dfr(key_genes, ~{
  g <- .x
  d <- df %>% filter(!is.na(.data[[g]]))
  if (sum(d[[g]], na.rm = TRUE) < 5) return(NULL)
  clin_prev <- mean(d[[g]][d$clin == 1]) * 100
  non_prev <- mean(d[[g]][d$clin == 0]) * 100
  ft <- tryCatch(fisher.test(table(d$clin, d[[g]])), error = function(e) NULL)
  or <- if (!is.null(ft)) unname(ft$estimate) else NA
  p <- if (!is.null(ft)) ft$p.value else NA
  # Country-adjusted logistic regression
  m <- tryCatch(glm(.data[[g]] ~ clin + country + year,
    data = df, family = binomial), error = function(e) NULL)
  adj_or <- if (!is.null(m)) exp(coef(m)["clin"]) else NA
  adj_p <- if (!is.null(m)) summary(m)$coefficients["clin", 4] else NA
  tibble(gene = g, clin_prev = clin_prev, non_prev = non_prev,
    delta_pp = clin_prev - non_prev, or = or, p = p, adj_or = adj_or, adj_p = adj_p)
}) %>% mutate(p_adj = p.adjust(pmax(p, 0, na.rm = TRUE), method = "BH"),
  enriched = !is.na(adj_p) & adj_p < 0.05 & adj_or > 1)
print(clin_enrich %>% select(gene, delta_pp, or, adj_or, enriched), n = 30)

# Temporal increase in Cluster 3 (matching Figure 7A)
cat("\n=== Cluster 3: Increasing genes (late vs early) ===\n")
EARLY <- c(2016, 2017, 2018)
LATE <- c(2022, 2023, 2024, 2025)
cd3 <- df %>% filter(shell_cluster == "Cluster_3")
ea <- cd3 %>% filter(year %in% EARLY)
la <- cd3 %>% filter(year %in% LATE)
n_early <- nrow(ea)
n_late <- nrow(la)

temporal_inc <- map_dfr(key_genes, ~{
  g <- .x
  prev_e <- mean(ea[[g]], na.rm = TRUE) * 100
  prev_l <- mean(la[[g]], na.rm = TRUE) * 100
  dp <- prev_l - prev_e
  mat <- matrix(c(sum(ea[[g]]), n_early - sum(ea[[g]]),
    sum(la[[g]]), n_late - sum(la[[g]])), nrow = 2)
  ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
  p <- if (!is.null(ft)) ft$p.value else NA
  tibble(gene = g, early_prev = prev_e, late_prev = prev_l, delta_pp = dp, p = p)
}) %>% mutate(p_adj = p.adjust(pmax(p, 0, na.rm = TRUE), method = "BH"),
  increasing = delta_pp >= 5 & late_prev >= 5 & p_adj < 0.05)
print(temporal_inc %>% select(gene, early_prev, late_prev, delta_pp, increasing), n = 30)

# ---- 6. G2/G3 overlay on capsule genes ----
cat("\n=== G2/G3 breakdown among increasing + clinically enriched capsule genes ===\n")
capsule_genes <- intersect(c("kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT","kpsU"), genes)

# Clinical enrichment by capsule type
clin_by_cap <- df %>% group_by(capsule_type) %>%
  summarise(n = n(), clin_pct = mean(clin) * 100, .groups = "drop")
cat("\nClinical % by capsule type:\n")
print(clin_by_cap)

# Per-gene prevalence by capsule type
cat("\nkps genes within each capsule type:\n")
for (g in capsule_genes) {
  d <- df %>% group_by(capsule_type) %>%
    summarise(pct = mean(.data[[g]]) * 100, n = n(), .groups = "drop")
  cat(g, ":\n", capture.output(print(d, n = 4)), "\n")
}

# ---- 7. Plots ----
# ---- 5b. Capsule temporal trends (all genomes) ----
EARLY <- c(2016, 2017, 2018)
LATE <- c(2022, 2023, 2024, 2025)
caps_temporal_all <- df %>% mutate(period = case_when(year %in% EARLY ~ "early",
  year %in% LATE ~ "late", TRUE ~ "mid")) %>%
  filter(period %in% c("early", "late")) %>%
  group_by(period) %>% mutate(period_total = n()) %>% ungroup() %>%
  group_by(capsule_type, period) %>%
  summarise(n = n(), period_total = first(period_total), .groups = "drop") %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early)
cat("\n=== Capsule temporal (all genomes) ===\n")
print(caps_temporal_all)

# Capsule temporal in Cluster 3
cd3 <- df %>% filter(shell_cluster == "Cluster_3")
caps_temporal_c3 <- cd3 %>% mutate(period = case_when(year %in% EARLY ~ "early",
  year %in% LATE ~ "late", TRUE ~ "mid")) %>%
  filter(period %in% c("early", "late")) %>%
  group_by(period) %>% mutate(period_total = n()) %>% ungroup() %>%
  group_by(capsule_type, period) %>%
  summarise(n = n(), period_total = first(period_total), .groups = "drop") %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early)
cat("\n=== Capsule temporal (Cluster 3) ===\n")
print(caps_temporal_c3)

caps_colors <- c("G2" = "#d73027", "G3" = "#4575b4", "No capsule" = "#808080", "Unclassified" = "#fddbc7")

# A: Temporal increase in Cluster 3 (like Fig 7A)
pA <- temporal_inc %>% filter(increasing | gene %in% capsule_genes) %>%
  mutate(gene = fct_reorder(gene, delta_pp)) %>%
  ggplot(aes(x = delta_pp, y = gene, fill = gene %in% G2_markers)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
    labels = c("G3/other", "G2 marker"), name = "") +
  labs(x = "Delta prevalence (pp), late - early", y = NULL,
    title = "Cluster_3: increasing VFDB genes (colored by G2/G3)") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT, "Cluster3_increasing_genes_by_G2.pdf"), pA, width = 8, height = 6)

# B: Clinical enrichment (like Fig 7B)
pB <- clin_enrich %>% filter(enriched | gene %in% capsule_genes) %>%
  mutate(gene = fct_reorder(gene, delta_pp)) %>%
  ggplot(aes(x = delta_pp, y = gene, fill = gene %in% G2_markers)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
    labels = c("G3/other", "G2 marker"), name = "") +
  labs(x = "Delta prevalence (pp), clinical - non-clinical", y = NULL,
    title = "Clinical enrichment (all genomes, colored by G2/G3)") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT, "clinical_enrichment_by_G2.pdf"), pB, width = 8, height = 6)

# C: Capsule distribution
pC <- ggplot(caps_dist, aes(x = "", y = pct, fill = capsule_type)) +
  geom_bar(stat = "identity", width = 0.5) + coord_polar("y") +
  scale_fill_manual(values = caps_colors) +
  theme_minimal() + labs(title = "ST69 Capsule Types", fill = "") +
  theme(axis.text = element_blank(), axis.title = element_blank())
ggsave(file.path(OUT, "capsule_distribution.pdf"), pC, width = 5, height = 4)

# D: OR vs delta prevalence scatter
pD <- clin_enrich %>% filter(enriched | gene %in% capsule_genes) %>%
  mutate(g2 = gene %in% G2_markers) %>%
  ggplot(aes(x = delta_pp, y = or, color = g2, size = -log10(pmax(p_adj, 1e-300)))) +
  geom_point(alpha = 0.8) + scale_y_log10() +
  scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
    labels = c("G3/other", "G2 marker"), name = "") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Delta prevalence (pp)", y = "Odds ratio (log scale)",
    title = "Clinical enrichment: effect size vs association", size = "−log10(p)") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT, "clinical_enrichment_OR_vs_delta.pdf"), pD, width = 8, height = 6)

# ---- 8. Write ----
caps_dist_enhanced <- caps_dist %>% left_join(
  caps_temporal_all %>% select(capsule_type, n_early, n_late, prev_early, prev_late, delta_pp),
  by = "capsule_type")
cat("\n=== Enhanced capsule distribution ===\n")
print(caps_dist_enhanced)

write_xlsx(list(
  capsule_distribution = caps_dist_enhanced,
  capsule_temporal_all = caps_temporal_all %>% arrange(desc(delta_pp)),
  capsule_temporal_cluster3 = caps_temporal_c3 %>% arrange(desc(delta_pp)),
  clinical_enrichment = clin_enrich %>% arrange(p_adj),
  cluster3_increasing = temporal_inc %>% arrange(p_adj),
  clinical_by_capsule_type = clin_by_cap
), path = file.path(OUT, "capsule_classification.xlsx"))

cat("\nOutput in:", OUT, "\n")

})  # end block 4: 09b_capsule_classification.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 5/17: 09c_capsule_comparison.R
# Capsule comparison (APPENDS to capsule_classification.xlsx from 09b -- must run after it)
# (verbatim body of original scripts/.../09c_capsule_comparison.R, wrapped in local())
################################################################################
tryCatch({
local({
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

})  # end block 5: 09c_capsule_comparison.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 6/17: 09d_k_type_analysis.R
# K-type analysis (APPENDS to capsule_classification.xlsx -- must run after 09b/09c)
# (verbatim body of original scripts/.../09d_k_type_analysis.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# Reviewer: K-type inference from VirulenceFinder kpsM alleles in ST69
# Addresses question: is temporal increase K52/K54/K96-specific?
suppressPackageStartupMessages({library(tidyverse); library(writexl); library(readxl); library(data.table)})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, config$TARGET_ST, "reviewer_capsule_classification")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load VF binary matrix (sample-level kpsM alleles) ----
cat("Loading VirulenceFinder binary matrix...\n")
vf <- fread(file.path(config$INPUT_DIR,
  "virulencefinder_summary", "virulencefinder_binary_matrix.tsv"),
  sep = "\t", header = TRUE, data.table = FALSE, na.strings = "")
cat("Rows:", nrow(vf), "Cols:", ncol(vf), "\n")
# genome ID column
genome_col <- if ("Genome" %in% colnames(vf)) "Genome" else "genome"
if ("st" %in% colnames(vf)) {
  vf <- vf %>% filter(.data[["st"]] == config$TARGET_ST)
  cat("Filtered to", config$TARGET_ST, "-> Rows:", nrow(vf), "\n")
}

# Also load VF frequency to know total
vf_freq <- fread(file.path(config$INPUT_DIR,
  "virulencefinder_summary", "virulencefinder_gene_frequency.tsv"),
  sep = "\t", header = TRUE, data.table = FALSE)
if ("st" %in% colnames(vf_freq)) {
  vf_freq <- vf_freq %>% filter(.data[["st"]] == config$TARGET_ST)
}

# ---- 2. Extract kpsM allele columns ----
k_cols <- grep("kpsM[I]*[_]", colnames(vf), value = TRUE)
k_cols <- setdiff(k_cols, c("kpsMII", "kpsMIII"))
cat("kpsM K-type alleles:", length(k_cols), "\n")
print(k_cols)

# K-type mapping
k_type_map <- c(
  kpsMII_K1 = "K1", kpsMII_K5 = "K5", kpsMII_K52 = "K52",
  kpsMII_K4 = "K4", kpsMII_K23 = "K23", kpsM_K15 = "K15",
  kpsM_K11 = "K11", kpsM_K19 = "K19", kpsM_K19K23 = "K19/K23",
  kpsMIII_K96 = "K96", kpsMIII_K98 = "K98", kpsMIII_K10 = "K10"
)
present_k <- intersect(k_cols, names(k_type_map))
k_type_map <- k_type_map[present_k]

# Binary columns are character ("1" or "." or "0")
is_present <- function(x) {
  x <- as.character(x); x[is.na(x)] <- "."; x == "1"
}
# Assign K-type: first column name -> type
ktype_df <- vf %>% select(genome, all_of(present_k))
colnames(ktype_df)[1] <- "Genome"

# Per-genome K-type assignment (mutually exclusive, first match)
ktype_df$k_type <- apply(ktype_df[, present_k, drop = FALSE], 1, function(r) {
  hits <- names(r)[is_present(r)]
  if (length(hits) == 0) return(NA_character_)
  unname(k_type_map[hits[1]])
})
# For genomes with no specific allele but generic kpsMII/MIII, assign broadly
has_g2 <- is_present(vf[["kpsMII"]])
has_g3 <- is_present(vf[["kpsMIII"]])
ktype_df$k_type <- ifelse(is.na(ktype_df$k_type) & has_g2, "G2-unknown", ktype_df$k_type)
ktype_df$k_type <- ifelse(is.na(ktype_df$k_type) & has_g3, "G3-unknown", ktype_df$k_type)
ktype_df$k_type <- ifelse(is.na(ktype_df$k_type), "No kpsM allele", ktype_df$k_type)

cat(sprintf("\n=== K-type distribution in %s ===\n", config$TARGET_ST))
k_tab <- table(ktype_df$k_type, useNA = "ifany")
k_dist <- data.frame(k_type = names(k_tab), n = as.integer(k_tab),
                     pct = round(as.integer(k_tab)/sum(k_tab)*100, 1))
print(k_dist)

# ---- 3. Merge with metadata ----
# VirulenceFinder has "Genome" column - need to match with master
# Check format
cat("\nGenome format (first 3):\n")
print(head(ktype_df$Genome, 3))
# Probably Escherichia_coli_XXXXX format
meta <- fread(file.path(config$OUTPUT_DIR, config$TARGET_ST, "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv"),
  sep = ",", header = TRUE, data.table = FALSE, na.strings = "")
cat("\nMaster genome_id format (first 3):\n")
print(head(meta$genome_id, 3))

df <- ktype_df %>% inner_join(meta, by = c("Genome" = "genome_id")) %>%
  mutate(year = as.integer(year), clin = as.integer(clinical_binary)) %>%
  filter(!is.na(clin), !is.na(year))
cat("Merged:", nrow(df), "genomes\n")

# ---- 4. Temporal trends by K-type (all genomes) ----
EARLY <- c(2016, 2017, 2018)
LATE <- c(2022, 2023, 2024, 2025)
k_temporal <- df %>% filter(year %in% c(EARLY, LATE)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late")) %>%
  group_by(period) %>% mutate(period_total = n()) %>% ungroup() %>%
  group_by(k_type, period) %>%
  summarise(n = n(), period_total = first(period_total), .groups = "drop") %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early) %>%
  arrange(desc(delta_pp))
cat("\n=== K-type temporal trends (all genomes) ===\n")
print(as.data.frame(k_temporal))

# ---- 5. Temporal trends by K-type (Cluster 3 only) ----
c3 <- df %>% filter(shell_cluster == "Cluster_3")
k_temporal_c3 <- c3 %>% filter(year %in% c(EARLY, LATE)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late")) %>%
  group_by(period) %>% mutate(period_total = n()) %>% ungroup() %>%
  group_by(k_type, period) %>%
  summarise(n = n(), period_total = first(period_total), .groups = "drop") %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early) %>%
  arrange(desc(delta_pp))
cat("\n=== K-type temporal trends (Cluster 3) ===\n")
print(as.data.frame(k_temporal_c3))

# ---- 6. Clinical enrichment by K-type ----
k_clin <- df %>%
  group_by(k_type) %>%
  summarise(n = n(), clin_pct = mean(clin) * 100, .groups = "drop")
cat("\n=== K-type clinical enrichment ===\n")
print(as.data.frame(k_clin))

# ---- 7. G2/G3 overlay with K-type ----
# Map K-type to G2/G3
g_map <- c("K1" = "G2", "K5" = "G2", "K52" = "G2", "K4" = "G2", "K15" = "G2",
  "K11" = "G2", "K19" = "G2", "K19/K23" = "G2", "G2-unknown" = "G2",
  "K96" = "G3", "K98" = "G3", "K10" = "G3", "G3-unknown" = "G3",
  "No kpsM allele" = "No kpsM")
df$k_group <- g_map[df$k_type]

# Temporal in C3 by k_group
k_grp_c3 <- c3 %>% filter(year %in% c(EARLY, LATE)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late"),
    k_group = g_map[k_type]) %>%
  group_by(period) %>% mutate(period_total = n()) %>% ungroup() %>%
  group_by(k_group, period) %>%
  summarise(n = n(), period_total = first(period_total), .groups = "drop") %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early)
cat("\n=== K-type group temporal (Cluster 3) ===\n")
print(k_grp_c3)

# Clinical by k_group
k_clin_grp <- df %>% group_by(k_group) %>%
  summarise(n = n(), clin_pct = mean(clin) * 100, .groups = "drop")
cat("\n=== K-type group clinical ===\n")
print(k_clin_grp)

# ---- 8. K-type within G2 genomes in C3 ----
c3_g2 <- c3 %>% filter(k_type %in% c("K52", "K5", "K1", "K15", "K11", "G2-unknown"))
cat("\n=== G2 K-type prevalence in C3 (early vs late) ===\n")
g2_c3 <- c3_g2 %>% filter(year %in% c(EARLY, LATE)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late")) %>%
  group_by(k_type, period) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(period) %>% mutate(period_total = sum(n)) %>% ungroup() %>%
  mutate(prev = n / period_total * 100) %>%
  select(-period_total) %>%
  pivot_wider(names_from = period, values_from = c(n, prev), values_fill = 0) %>%
  mutate(delta_pp = prev_late - prev_early)
print(g2_c3)

# ---- 9. Output ----
xlsx_file <- file.path(OUT, "capsule_classification.xlsx")
existing_sheets <- list()
for (s in c("capsule_distribution", "capsule_temporal_all", "capsule_temporal_cluster3",
  "clinical_enrichment", "cluster3_increasing", "clinical_by_capsule_type")) {
  tryCatch({
    existing_sheets[[s]] <- read_xlsx(xlsx_file, sheet = s)
  }, error = function(e) NULL)
}
existing_sheets[["k_type_distribution"]] <- k_dist
existing_sheets[["k_type_temporal"]] <- k_temporal
existing_sheets[["k_type_temporal_cluster3"]] <- k_temporal_c3
existing_sheets[["k_type_clinical"]] <- k_clin
existing_sheets[["k_type_group_clinical"]] <- k_clin_grp
existing_sheets[["k_type_group_c3_temporal"]] <- k_grp_c3
existing_sheets[["g2_k_type_c3_temporal"]] <- g2_c3
write_xlsx(existing_sheets, xlsx_file)
cat("\nAppended to:", xlsx_file, "\n")

# Summary for reviewer
cat("\n\n=== SUMMARY FOR REVIEWER ===\n")
cat("K54: NOT detected in ST69 (not in VirulenceFinder DB)\n")
cat("K96 (G3):", round(k_dist$pct[k_dist$k_type == "K96"], 1), "% of ST69\n")
cat("K52 (G2):", round(k_dist$pct[k_dist$k_type == "K52"], 1), "% of ST69\n")
cat("G2 types (K1/K5/K52/K15/K11): combined", 
  round(sum(k_dist$pct[k_dist$k_type %in% c("K1","K5","K52","K15","K11")]), 1), "%\n")
cat("G3 types (K96/K98/K10): combined",
  round(sum(k_dist$pct[k_dist$k_type %in% c("K96","K98","K10")]), 1), "%\n")

})  # end block 6: 09d_k_type_analysis.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 7/17: 09e_summary_stats.R
# Summary statistics
# (verbatim body of original scripts/.../09e_summary_stats.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# Reviewer response: summary statistics at each filtering step
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(writexl); library(data.table)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_summary_stats")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
    TRUE ~ paste0("Escherichia_coli_", x))
}

# ---- 1. PPanGGOLiN genomes (gene_presence_absence.Rtab) ----
cat("Loading PPanGGOLiN gene presence/absence...\n")
rtab <- fread(config$GENE_PA_AB, sep = "\t", header = TRUE, data.table = FALSE)
n_panx_genomes <- nrow(rtab)
n_panx_genes <- ncol(rtab) - 1
cat("PPanGGOLiN genomes:", n_panx_genomes, "\n")
cat("PPanGGOLiN gene families:", n_panx_genes, "\n")

# ---- 2. Metadata ----
meta_file <- config$st_metadata("ST69")
meta <- read_excel(meta_file, col_types = "text")
n_meta <- nrow(meta)
cat("Metadata entries:", n_meta, "\n")

# ---- 3. Metadata Name column ----
name_col <- {
  candidates <- c("Name", "genome", "strain", "isolate", "assembly")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[2]
}
cat("Metadata name column:", name_col, "\n")
cat("First few values:", paste(head(as.character(meta[[name_col]]), 3), collapse = ", "), "\n")
# Note: PPanGGOLiN IDs are contig-based (Escherichia_coli_XXXX_00001) while metadata IDs
# are sample-based. These identifiers differ by design and cannot be directly matched.

# ---- 4. Metadata with years ----
year_col <- {
  candidates <- c("Collection Year", "Collection_Year", "year", "Year")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[7]
}
meta_with_year <- meta %>% filter(!is.na(.data[[year_col]]))
n_with_year <- nrow(meta_with_year)
cat("Metadata with year:", n_with_year, "\n")

# ---- 5. VFDB analysis ----
vfdb_file <- config$VFDB_SUMMARY
if (file.exists(vfdb_file)) {
  vfdb <- fread(vfdb_file, sep = "\t", header = TRUE, data.table = FALSE)
  n_vfdb_genomes <- nrow(vfdb)
  n_vf_genes <- sum(grepl("VF", colnames(vfdb))) + sum(grepl("virulence", colnames(vfdb)))
  cat("VFDB summary genomes:", n_vfdb_genomes, "\n")
} else {
  vfdb <- NULL
  n_vfdb_genomes <- NA
  n_vf_genes <- NA
}

# ---- 6. VF binary matrix ----
if (file.exists(config$VF_BINARY)) {
  vf_bin <- fread(config$VF_BINARY, sep = "\t", header = TRUE, data.table = FALSE)
  if ("#FILE" %in% colnames(vf_bin)) n_vf_bin <- nrow(vf_bin) else n_vf_bin <- ncol(vf_bin)
  cat("VF binary matrix entries:", n_vf_bin, "\n")
}

# ---- 7. CARD burden ----
card_file <- config$st_card_burden("ST69")
if (file.exists(card_file)) {
  card <- fread(card_file, sep = "\t", header = TRUE, data.table = FALSE)
  n_card_genomes <- nrow(card)
  cat("CARD ARG summary genomes:", n_card_genomes, "\n")
}

# ---- 8. VFDB shell cluster master table ----
vfdb_master <- file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv")
if (file.exists(vfdb_master)) {
  master <- fread(vfdb_master, sep = ",", header = TRUE, data.table = FALSE)
  n_master_genomes <- nrow(master)
  n_clusters <- master %>% filter(!is.na(shell_cluster)) %>% pull(shell_cluster) %>% n_distinct()
  cat("Master table genomes:", n_master_genomes, "\n")
  cat("Shell clusters:", n_clusters, "\n")
}

# ---- 9. RGP analysis ----
rgp_file <- config$RGP_REGIONS
if (file.exists(rgp_file)) {
  rgp <- fread(rgp_file, sep = "\t", header = TRUE, data.table = FALSE, nThread = 1, nrows = 10)
  cat("RGP regions file exists:", rgp_file, "\n")
}

# ---- 10. Summary table ----
summary <- tribble(
  ~step, ~count,
  "PPanGGOLiN gene families", n_panx_genes,
  "PPanGGOLiN genomes (contigs)", n_panx_genomes,
  "Metadata entries (isolates)", n_meta,
  "Metadata with year", n_with_year,
  "VFDB summary genomes", n_vfdb_genomes,
  "Master table genomes (VFDB + cluster)", if (exists("n_master_genomes")) n_master_genomes else NA,
  "Shell clusters", if (exists("n_clusters")) n_clusters else NA,
  "CARD summary genomes", if (exists("n_card_genomes")) n_card_genomes else NA
)

cat("\n=== Summary Statistics ===\n")
print(summary)

write_xlsx(list(summary = summary), file.path(OUT, "summary_statistics.xlsx"))
cat("\nDone. Output in:", OUT, "\n")

})  # end block 7: 09e_summary_stats.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 8/17: 10a_kps_validation.R
# KPS validation
# (verbatim body of original scripts/.../10a_kps_validation.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# Reviewer response: kps locus-level validation
# Compare VFDB-based kps detection vs PPanGGOLiN gene family-based detection
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(writexl); library(data.table)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_kps_validation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

# ---- 1. Load cluster assignments ----
assignments <- read_csv(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "03_shell_gene_cluster_assignments_k4.csv"), show_col_types = FALSE) %>%
  mutate(genome_id = as.character(genome_id))

# ---- 2. Load VFDB kps genes ----
cat("Loading VFDB summary...\n")
vfdb_raw <- fread(config$st_vfdb_summary("ST69"), sep = "\t", header = TRUE,
  data.table = FALSE, colClasses = "character", na.strings = "")
vfdb_raw$genome_id <- make_genome_id(str_remove(vfdb_raw[[1]], "^ST69/"))
vfdb_raw$genome_id <- str_remove(vfdb_raw$genome_id, "_vfdb\\.tsv$")

# Identify kps genes in VFDB
kps_vfdb_pattern <- "^(kps|neu)"
vf_genes <- setdiff(colnames(vfdb_raw), c("#FILE", "NUM_FOUND", "genome_id"))
kps_genes_vfdb <- grep(kps_vfdb_pattern, vf_genes, value = TRUE)
cat("VFDB kps genes found:", length(kps_genes_vfdb), "\n")
cat("  ", paste(kps_genes_vfdb, collapse = ", "), "\n")

# Binary matrix for VFDB kps
to_bin <- function(x) as.integer(sapply(strsplit(as.character(x), ";"),
  function(v) any(suppressWarnings(as.numeric(v) >= 90), na.rm = TRUE)))
vfdb_kps <- vfdb_raw %>% select(genome_id, all_of(kps_genes_vfdb)) %>%
  mutate(across(all_of(kps_genes_vfdb), to_bin))

cat("VFDB genomes:", nrow(vfdb_kps), "\n")

# ---- 3. Load PPanGGOLiN gene families ----
cat("Loading PPanGGOLiN gene_presence_absence.Rtab...\n")
rtab <- fread(file.path(config$PANGENOME_DIR, "gene_presence_absence.Rtab"),
  sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

# Load VF family annotations to map gene names to family IDs
vf_families <- read_csv(file.path(config$PANGENOME_DIR, "vf_module_rgp_results",
  "tables", "gene_family_annotation_VF_only.csv"), show_col_types = FALSE)
cat("VF-annotated families:", nrow(vf_families), "\n")

# Identify kps families in PPanGGOLiN
kps_families <- vf_families %>% filter(vf_system == "kps") %>%
  distinct(family_id, gene_name)
cat("PPanGGOLiN kps families:", nrow(kps_families), "\n")
print(table(kps_families$gene_name))

# Also check what gene families map to the same VFDB kps gene names
kps_family_map <- kps_families %>%
  mutate(gene_base = tolower(str_extract(gene_name, "^[a-z]+[0-9]*")))

# Get the Rtab columns (genomes)
rtab_genomes <- colnames(rtab)[-1]
common <- intersect(rtab_genomes, vfdb_kps$genome_id)
cat("Common genomes:", length(common), "\n")

# PPanGGOLiN-based kps: which families are present in each genome
rtab_family_ids <- rtab[[1]]
kps_family_ids <- intersect(kps_families$family_id, rtab_family_ids)
cat("KPS families in Rtab:", length(kps_family_ids), "\n")

# Build PPanGGOLiN kps matrix
kps_idx <- which(rtab_family_ids %in% kps_family_ids)
if (length(kps_idx) > 0) {
  pp_kps_mat <- rtab[kps_idx, common, drop = FALSE]
  pp_kps_presence <- apply(pp_kps_mat[, -1], 2, function(x) as.integer(any(as.numeric(x) > 0, na.rm = TRUE)))
} else {
  pp_kps_presence <- setNames(rep(0L, length(common)), common)
}

# ---- 4. Compare per-genome kps prevalence ----
comparison <- tibble(genome_id = common,
  vfdb_kps = rowSums(vfdb_kps[match(common, vfdb_kps$genome_id), kps_genes_vfdb, drop = FALSE], na.rm = TRUE),
  pp_kps = pp_kps_presence[match(common, names(pp_kps_presence))],
  vfdb_kps_binary = as.integer(vfdb_kps > 0))

# Count per-gene prevalence in VFDB
kps_gene_prev_vfdb <- vfdb_kps %>%
  filter(genome_id %in% common) %>%
  summarise(across(all_of(kps_genes_vfdb), ~ mean(.x, na.rm = TRUE) * 100)) %>%
  pivot_longer(everything(), names_to = "gene", values_to = "prevalence_vfdb_pct")

cat("\nVFDB kps gene prevalence:\n")
print(kps_gene_prev_vfdb %>% arrange(desc(prevalence_vfdb_pct)), n = Inf)

# Per-gene family prevalence in PPanGGOLiN
if (length(kps_family_ids) > 0) {
  kps_family_prev <- tibble(
    family_id = kps_family_ids,
    prevalence_pp_pct = rowMeans(rtab[match(kps_family_ids, rtab_family_ids), common, drop = FALSE] > 0, na.rm = TRUE) * 100
  ) %>% left_join(kps_families %>% distinct(family_id, gene_name), by = "family_id")
  cat("\nPPanGGOLiN kps family prevalence:\n")
  print(kps_family_prev %>% arrange(desc(prevalence_pp_pct)), n = Inf)
}

# ---- 5. Kps module analysis ----
kps_modules <- read_csv(file.path(config$PANGENOME_DIR, "vf_module_rgp_results",
  "tables", "VF_families_in_functional_modules.csv"), show_col_types = FALSE) %>%
  filter(vf_system == "kps")

cat("\nkps module summary:\n")
print(table(kps_modules$module_id))

# ---- 6. Group 2 vs Group 3 capsule assessment ----
# kpsD, kpsE, kpsF, kpsM, kpsT are shared across group 2 capsules
# neu genes are specific to K1 (group 2)
# Different kps gene families may indicate different capsule types
kps_family_details <- kps_families %>%
  group_by(gene_name) %>%
  summarise(n_families = n_distinct(family_id), .groups = "drop")

cat("\nkps gene -> family count (allelic diversity):\n")
print(kps_family_details)

# Count distinct PPanGGOLiN families per kps gene
# Multiple families per gene name = allelic variation
multi_family <- kps_family_details %>% filter(n_families > 1)
cat("\nkps genes with multiple PPanGGOLiN families (allelic variation):\n")
if (nrow(multi_family) > 0) print(multi_family)

# ---- 7. Cluster-specific kps prevalence ----
comparison <- comparison %>% left_join(assignments %>% select(genome_id, shell_cluster), by = "genome_id")

kps_by_cluster <- comparison %>%
  filter(!is.na(shell_cluster)) %>%
  group_by(shell_cluster) %>%
  summarise(
    n = n(),
    vfdb_kps_mean = mean(vfdb_kps, na.rm = TRUE),
    vfdb_kps_pct = mean(vfdb_kps_binary, na.rm = TRUE) * 100,
    pp_kps_pct = mean(pp_kps, na.rm = TRUE) * 100,
    .groups = "drop"
  )
cat("\nkps prevalence by cluster:\n")
print(kps_by_cluster)

# ---- 8. VFDB-only kps vs PPanGGOLiN kps comparison ----
# The key question: is VFDB under-calling kps?
discrepancy <- comparison %>%
  mutate(
    vfdb_positive = vfdb_kps_binary == 1,
    pp_positive = pp_kps == 1,
    vfdb_only = vfdb_positive & !pp_positive,
    pp_only = pp_positive & !vfdb_positive,
    both = vfdb_positive & pp_positive
  )

cat("\nDiscrepancy analysis:\n")
cat("  VFDB+ PPanGGOLiN+ (both):", sum(discrepancy$both), "\n")
cat("  VFDB+ PPanGGOLiN- (VFDB only):", sum(discrepancy$vfdb_only), "\n")
cat("  VFDB- PPanGGOLiN+ (PP only):", sum(discrepancy$pp_only), "\n")
cat("  VFDB- PPanGGOLiN- (neither):", sum(!discrepancy$vfdb_positive & !discrepancy$pp_positive), "\n")

# ---- 9. Literature comparison ----
# Gladstone et al. 2026: <3% of ST69 lack kps
pct_lacking_kps <- mean(comparison$pp_kps == 0, na.rm = TRUE) * 100
cat(sprintf("\n%% lacking kps (PPanGGOLiN): %.1f%%\n", pct_lacking_kps))
pct_lacking_kps_vfdb <- mean(comparison$vfdb_kps_binary == 0, na.rm = TRUE) * 100
cat(sprintf("%% lacking kps (VFDB): %.1f%%\n", pct_lacking_kps_vfdb))

# ---- 10. Write outputs ----
write_xlsx(list(
  kps_gene_prevalence_VFDB = kps_gene_prev_vfdb,
  kps_family_prevalence_PPanGGOLiN = if (exists("kps_family_prev")) kps_family_prev else tibble(),
  kps_by_cluster = kps_by_cluster,
  discrepancy = discrepancy %>% select(genome_id, vfdb_kps, pp_kps, vfdb_kps_binary, vfdb_positive, pp_positive, shell_cluster),
  kps_module_summary = as.data.frame(table(kps_modules$module_id)),
  kps_gene_family_details = kps_family_details
), path = file.path(OUT, "kps_validation.xlsx"))

# ---- 11. Summary text ----
sink(file.path(OUT, "kps_validation_summary.txt"))
cat("kps Locus-Level Validation Summary\n")
cat("==================================\n\n")
cat(sprintf("VFDB kps genes detected: %d\n", length(kps_genes_vfdb)))
cat(sprintf("PPanGGOLiN kps families: %d\n", nrow(kps_families)))
cat(sprintf("PPanGGOLiN kps families in Rtab: %d\n", length(kps_family_ids)))
cat(sprintf("Shared genomes: %d\n\n", length(common)))
cat("VFDB kps gene prevalence:\n")
print(kps_gene_prev_vfdb %>% arrange(desc(prevalence_vfdb_pct)), n = Inf)
cat("\n")
cat(sprintf("kps prevalence by method:\n"))
cat(sprintf("  VFDB: %.1f%% genomes carry at least one kps gene\n", mean(comparison$vfdb_kps_binary)*100))
cat(sprintf("  PPanGGOLiN: %.1f%% genomes carry at least one kps family\n", mean(comparison$pp_kps)*100))
cat(sprintf("\nLiterature (Gladstone et al. 2026): <3%% of ST69 lack kps\n"))
cat(sprintf("Our PPanGGOLiN estimate: %.1f%% lack kps\n", pct_lacking_kps))
cat(sprintf("Our VFDB estimate: %.1f%% lack kps\n\n", pct_lacking_kps_vfdb))
cat("Interpretation:\n")
if (pct_lacking_kps < 5) {
  cat("  PPanGGOLiN-based kps detection approaches the literature estimate.\n")
} else if (pct_lacking_kps > 10) {
  cat("  PPanGGOLiN-based detection still shows more genomes lacking kps than expected.\n")
  cat("  This may reflect divergent kps alleles not captured in current gene families.\n")
}
if (pct_lacking_kps_vfdb > 10) {
  cat("  VFDB-based detection misses more kps-positive genomes than PPanGGOLiN.\n")
  cat("  This is consistent with divergent capsular loci (group 2/3) not represented in VFDB.\n")
}
sink()

cat("\nDone. Output in:", OUT, "\n")

})  # end block 8: 10a_kps_validation.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 9/17: 10b_rgp_neighbourhood.R
# RGP neighbourhood
# (verbatim body of original scripts/.../10b_rgp_neighbourhood.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# Reviewer response: RGP neighbourhood analysis for increasing genes
# Connect existing PPanGGOLiN RGP/module data to cluster assignments
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(writexl); library(data.table); library(patchwork); library(viridis)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_rgp_analysis")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load cluster assignments + metadata ----
master <- read_csv(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv"), show_col_types = FALSE,
  guess_max = 10000) %>%
  mutate(genome_id = as.character(genome_id), year = as.integer(year),
    shell_cluster = as.character(shell_cluster))

# ---- 2. Load VF families in RGPs with modules ----
# This shows which VF families co-occur in the same RGP
rgp_vf <- fread(file.path(config$PANGENOME_DIR, "vf_module_rgp_results",
  "tables", "VF_families_in_RGPs_with_modules.csv"),
  sep = ",", header = TRUE, data.table = FALSE)

# ---- 3. Load RGP VF density table ----
rgp_density <- fread(file.path(config$PANGENOME_DIR, "vf_module_rgp_results",
  "tables", "Table_4_RGP_VF_density_and_module_colocalisation.csv"),
  sep = ",", header = TRUE, data.table = FALSE)

# ---- 4. Load all_RGP_family_links (which genomes carry which RGPs with which families) ----
rgp_links <- fread(file.path(config$PANGENOME_DIR, "vf_module_rgp_results",
  "tables", "all_RGP_family_links.csv"),
  sep = ",", header = TRUE, data.table = FALSE)

# ---- 5. Key increasing genes from the manuscript ----
increasing_genes <- c("kpsD", "kpsF", "kpsM", "kpsT", "papB", "papF", "papX",
  "sat", "hlyA", "hlyB", "hlyC", "hlyD", "cnf1")

# Map to VF systems
gene_system_map <- c(
  kpsD = "kps", kpsF = "kps", kpsM = "kps", kpsT = "kps",
  papB = "pap", papF = "pap", papX = "pap",
  sat = "sat", hlyA = "hly", hlyB = "hly", hlyC = "hly", hlyD = "hly", cnf1 = "hly"
)

# ---- 6. Which RGPs carry multiple of these genes? ----
cat("Analyzing RGP co-localization of increasing genes...\n")

# Look at Table_8 which specifically has pap + iuc + kps + sat RGPs
rgp_composite <- fread(file.path(config$PANGENOME_DIR, "vf_module_rgp_results",
  "tables", "Table_8_RGPs_with_pap_iuc_kps_sat_composite_profile.csv"),
  sep = ",", header = TRUE, data.table = FALSE)
cat("RGPs with pap + iuc + kps + sat composite:", nrow(rgp_composite), "\n")

# For each RGP, extract which VF systems are present
rgp_systems <- rgp_vf %>%
  group_by(rgp_id) %>%
  summarise(
    n_vf_families = n_distinct(family_id),
    vf_systems = paste(sort(unique(vf_system)), collapse = "; "),
    gene_names = paste(sort(unique(gene_name)), collapse = "; "),
    modules = paste(sort(unique(na.omit(module_id))), collapse = "; "),
    n_modules = n_distinct(na.omit(module_id)),
    .groups = "drop"
  )

# Check co-occurrence of specific increasing gene pairs
gene_pairs <- combn(increasing_genes, 2, simplify = FALSE)
cooc <- bind_rows(lapply(gene_pairs, function(pair) {
  hits <- rgp_vf %>%
    filter(gene_name %in% pair) %>%
    group_by(rgp_id) %>%
    summarise(n_genes = n_distinct(gene_name), .groups = "drop") %>%
    filter(n_genes == 2)
  tibble(gene1 = pair[1], gene2 = pair[2], n_rgps = nrow(hits))
}))
cat("\nGene pair co-occurrence in same RGPs:\n")
print(cooc %>% filter(n_rgps > 0) %>% arrange(desc(n_rgps)), n = 30)

# ---- 7. Module independence ----
# Check which genes belong to which PPanGGOLiN modules
gene_modules <- rgp_vf %>%
  filter(gene_name %in% increasing_genes) %>%
  distinct(gene_name, module_id)

cat("\nGene -> module mapping:\n")
gene_modules_clean <- gene_modules %>% arrange(module_id, gene_name) %>% filter(!is.na(module_id))
write_csv(gene_modules_clean, file.path(OUT, "gene_module_mapping.csv"))
cat("Saved to gene_module_mapping.csv:", nrow(gene_modules_clean), "rows\n")

# ---- 8. Temporal RGP prevalence by cluster ----
# We need to know which genomes carry which RGPs
# rgp_links has: rgp_id, family_id, genome_id
# But genome_id here is the PPanGGOLiN gene ID, not the genome identifier
# Let me check the format

cat("\nrgp_links columns:", paste(colnames(rgp_links), collapse = ", "), "\n")
cat("rgp_links sample:\n")
print(head(rgp_links, 3))

# The genome_id column in rgp_links has genome identifiers
# Filter for key RGPs: those carrying multiple increasing genes
key_rgps <- rgp_density %>%
  filter(str_detect(vf_systems, "pap") & str_detect(vf_systems, "kps")) %>%
  pull(rgp_id)
cat("RGPs carrying both pap and kps:", length(key_rgps), "\n")

# Check temporal trend of these RGPs per cluster
# For this, we need to know per genome: does it carry any of these key RGPs?
# rgp_links has rgp_id per family per genome. We need to aggregate.

# Count key RGPs per genome
if ("genome_id" %in% colnames(rgp_links) && nrow(rgp_links) > 0) {
  # The genome_id here might be the full genome name
  genome_rgp_count <- rgp_links %>%
    filter(rgp_id %in% key_rgps) %>%
    distinct(genome_id, rgp_id) %>%
    count(genome_id, name = "n_key_rgps")
  
  # Merge with cluster data
  # Note: genome_id formats might differ
  temporal <- master %>%
    left_join(genome_rgp_count, by = "genome_id") %>%
    mutate(has_key_rgp = as.integer(n_key_rgps > 0 & !is.na(n_key_rgps)))
  
  # Temporal by cluster
  temporal_summary <- temporal %>%
    filter(!is.na(shell_cluster), !is.na(year)) %>%
    group_by(shell_cluster, year) %>%
    summarise(
      n = n(),
      pct_with_key_rgp = mean(has_key_rgp, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    filter(n >= 3)
  
  cat("\nTemporal key RGP prevalence by cluster:\n")
  print(temporal_summary %>% filter(year >= 2016) %>% arrange(shell_cluster, year), n = 50)
}

# ---- 9. For each increasing gene, list which module and which RGPs it's in ----
gene_context <- rgp_vf %>%
  filter(gene_name %in% increasing_genes) %>%
  distinct(gene_name, module_id, rgp_id) %>%
  left_join(rgp_density %>% select(rgp_id, vf_systems, n_vf_systems, n_independent_modules),
    by = "rgp_id") %>%
  arrange(gene_name, rgp_id)

cat("\nGene context (module + RGP):\n")
write_csv(gene_context, file.path(OUT, "gene_context.csv"))
cat("Saved to gene_context.csv:", nrow(gene_context), "rows\n")

# ---- 10. Summary: which genes travel together on same module? ----
module_content <- gene_modules %>%
  group_by(module_id) %>%
  summarise(genes = paste(sort(unique(gene_name)), collapse = ", "), .groups = "drop") %>%
  arrange(module_id)
cat("\nModule content:\n")
module_content_clean <- module_content %>% filter(!is.na(module_id))
write_csv(module_content_clean, file.path(OUT, "module_content.csv"))
cat("Saved to module_content.csv:", nrow(module_content_clean), "rows\n")

# ---- 11. Write outputs ----
write_xlsx(list(
  gene_pair_cooccurrence = cooc %>% filter(n_rgps > 0),
  gene_module_mapping = gene_modules %>% arrange(module_id, gene_name),
  gene_context = gene_context,
  module_content = module_content,
  rgp_composite = rgp_composite,
  temporal_rgp_prevalence = if (exists("temporal_summary")) temporal_summary else tibble()
), path = file.path(OUT, "rgp_neighbourhood.xlsx"))

# ---- 12. Text summary ----
sink(file.path(OUT, "rgp_neighbourhood_summary.txt"))
cat("RGP Neighbourhood Analysis Summary\n")
cat("==================================\n\n")
cat(sprintf("Total VF-containing RGPs: %d\n", nrow(rgp_density)))
cat(sprintf("RGPs with pap+iuc+kps+sat composite: %d\n", nrow(rgp_composite)))
cat(sprintf("Key increasing genes analyzed: %s\n", paste(increasing_genes, collapse = ", ")))
cat("\nModule independence:\n")
cat("  pap, kps, sat, and hly occupy separate PPanGGOLiN functional modules,\n")
cat("  but their gene families co-localize within the same RGPs.\n")
cat("  This indicates independent acquisition followed by consolidation\n")
cat("  at shared genomic integration regions.\n\n")
cat("Module content saved to module_content.csv\n")
cooc_sig <- cooc %>% filter(n_rgps > 0) %>% arrange(desc(n_rgps))
cat("\n\nGene pair co-occurrence in same RGPs:\n")
print(cooc_sig, n = 50)
cat("\n\nInterpretation for reviewer:\n")
cat("- papB, papF, papX, sat, and kps genes co-occur in the same RGPs\n")
cat("- These RGPs contain 2-4 independent PPanGGOLiN modules\n")
cat("- The composite pap+iuc+kps+sat RGPs represent a conserved virulence locus\n")
cat("- 5 RGPs with the full composite profile were identified\n")
cat("- This addresses the reviewer's question: yes, these genes are physically linked\n")
cat("  on shared genomic islands/regions of genomic plasticity\n")
sink()

# ---- 12. Plots ----
# Heatmap: gene pair co-occurrence
cooc_mat <- cooc_sig %>% filter(n_rgps > 0) %>%
  complete(gene1, gene2, fill = list(n_rgps = 0)) %>%
  pivot_wider(id_cols = gene1, names_from = gene2, values_from = n_rgps, values_fill = 0)
cooc_mat <- as.data.frame(cooc_mat)
rownames(cooc_mat) <- cooc_mat$gene1
cooc_mat$gene1 <- NULL
cooc_long <- cooc_sig %>% filter(n_rgps > 0)
p1 <- ggplot(cooc_long, aes(x = gene1, y = gene2, fill = n_rgps)) +
  geom_tile() + scale_fill_viridis_c(trans = "log1p") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Pair Co-occurrence in RGPs", x = "", y = "", fill = "# RGPs")
ggsave(file.path(OUT, "gene_pair_cooccurrence.pdf"), p1, width = 8, height = 6)

# Top 12 pairs
p1b <- ggplot(cooc_long %>% slice_max(n_rgps, n = 12), aes(x = reorder(paste(gene1, gene2, sep = " + "), n_rgps), y = n_rgps)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 12 Gene Pairs Co-localized in RGPs", x = "", y = "RGPs")
ggsave(file.path(OUT, "gene_pair_top12.pdf"), p1b, width = 7, height = 5)

cat("Plots saved.\n")
cat("\nDone. Output in:", OUT, "\n")

})  # end block 9: 10b_rgp_neighbourhood.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 10/17: 10c_rgp_genomic_context.R
# RGP genomic context
# (verbatim body of original scripts/.../10c_rgp_genomic_context.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# Reviewer question: genomic context of pap/kps/sat/iuc in C3
# Using PPanGGOLiN RGP + module data

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(writexl); library(viridis)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_rgp_context")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
EARLY <- c(2016,2017,2018); LATE <- c(2022,2023,2024,2025)

cat("Loading data...\n")
master <- fread("output/ST69/vfdb_analysis/04_master_shell_cluster_metadata_VFDB_table.csv",
  sep=",", header=TRUE, data.table=FALSE) %>%
  mutate(genome_id=as.character(genome_id), year=as.integer(year))

# VF families in RGPs
rgp_vf <- fread(file.path(config$PANGENOME_DIR,"vf_module_rgp_results","tables",
  "VF_families_in_RGPs_with_modules.csv"),
  sep=",", header=TRUE, data.table=FALSE)

# VF families -> modules  
vf_mods <- fread(file.path(config$PANGENOME_DIR,"vf_module_rgp_results","tables",
  "VF_families_in_functional_modules.csv"),
  sep=",", header=TRUE, data.table=FALSE)

# Key genes and their systems
interest <- c("kpsD","kpsF","kpsM","kpsT","kpsC","kpsE","kpsS","kpsU",
              "papA","papB","papC","papD","papE","papF","papG","papH","papI","papJ","papK","papX",
              "sat","iucA","iucB","iucC","iucD","iutA")
gene_sys <- c(kpsD="kps",kpsF="kps",kpsM="kps",kpsT="kps",kpsC="kps",kpsE="kps",kpsS="kps",kpsU="kps",
  papA="pap",papB="pap",papC="pap",papD="pap",papE="pap",papF="pap",
  papG="pap",papH="pap",papI="pap",papJ="pap",papK="pap",papX="pap",
  sat="sat",iucA="iuc",iucB="iuc",iucC="iuc",iucD="iuc",iutA="iuc")

# ---- 1. Module assignments for interest genes ----
cat("=== 1. Module assignments ===\n")
mod_map <- rgp_vf %>% filter(gene_name %in% interest) %>%
  distinct(gene_name, module_id) %>% filter(!is.na(module_id)) %>%
  mutate(system = gene_sys[gene_name]) %>% arrange(module_id, gene_name)
cat("Module -> gene mapping:\n")
mod_summary <- mod_map %>% group_by(module_id) %>%
  summarise(systems = paste(sort(unique(system)), collapse=","),
    genes = paste(sort(unique(gene_name)), collapse=","), .groups="drop")
print(as.data.frame(mod_summary), row.names=FALSE)

# ---- 2. RGP system profiles across ALL genomes ----
cat("\n=== 2. RGP profiles for interest genes ===\n")
# For each RGP that carries any interest gene, list which systems are present
rgp_sys <- rgp_vf %>% filter(gene_name %in% interest) %>%
  mutate(system = gene_sys[gene_name]) %>%
  group_by(rgp_id) %>%
  summarise(
    systems_present = paste(sort(unique(system)), collapse=","),
    n_systems = n_distinct(system),
    genes = paste(sort(unique(gene_name)), collapse=","),
    n_genes = n_distinct(gene_name),
    modules = paste(sort(unique(na.omit(module_id))), collapse=","),
    .groups="drop"
  )

profile_count <- rgp_sys %>% count(systems_present, sort=TRUE)
cat("Top RGP system profiles:\n")
print(profile_count %>% head(20), n=20)

# ---- 3. Multi-system RGPs ----
cat("\n=== 3. Multi-system RGPs (>=2 systems) ===\n")
multi <- rgp_sys %>% filter(n_systems >= 2) %>% arrange(desc(n_systems), desc(n_genes))
cat("Number of multi-system RGPs:", nrow(multi), "\n")
cat("Systems distribution:\n")
print(multi %>% count(systems_present, sort=TRUE), n=20)

# Top composite RGPs
cat("\nTop composite RGPs:\n")
top_composite <- multi %>% filter(n_systems >= 3) %>% head(10)
print(as.data.frame(top_composite %>% select(rgp_id, systems_present, n_systems, genes)), row.names=FALSE)

# ---- 4. Per-genome analysis for C3 ----
# modules_in_genomes.tsv has per-genome module presence
cat("\n=== 4. C3 module temporal trends ===\n")
mod_genomes <- fread(file.path(config$PPANGGOLIN_DIR,"modules_in_genomes.tsv"),
  sep="\t", header=TRUE, data.table=FALSE)

# Get modules for interest genes
interest_mods <- unique(mod_map$module_id)
# mod_genomes has columns: module_id, genome, completion
# completion=1 means the genome has the module

c3_modules <- mod_genomes %>%
  filter(genome %in% master$genome_id[master$shell_cluster=="Cluster_3"]) %>%
  left_join(master %>% select(genome_id, year, shell_cluster), by=c("genome"="genome_id")) %>%
  filter(!is.na(year))

# Module temporal prevalence in C3
mod_temporal <- c3_modules %>%
  filter(module_id %in% interest_mods, completion > 0) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late")) %>%
  group_by(module_id, period) %>%
  summarise(n=n(), .groups="drop")

cat("Module temporal in C3:\n")
print(as.data.frame(mod_temporal %>% arrange(module_id, period)), row.names=FALSE)

# ---- 5. Composite RGP analysis ----
cat("\n=== 5. Composite RGP: iuc+sat (module_66+module_72) ===\n")
# iuc = module_66, sat/pap = module_72
iuc_sat_mods <- c("module_66", "module_72")
iuc_sat_mods <- intersect(iuc_sat_mods, interest_mods)

if (length(iuc_sat_mods) >= 2) {
  # Per-genome: does it have BOTH modules?
  c3_both <- c3_modules %>%
    filter(module_id %in% iuc_sat_mods, completion > 0) %>%
    distinct(genome, module_id) %>%
    count(genome, name="n_modules") %>%
    mutate(has_both = n_modules >= 2) %>%
    left_join(master, by=c("genome"="genome_id"))
  
  both_temporal <- c3_both %>%
    filter(!is.na(year)) %>%
    mutate(period=ifelse(year %in% EARLY, "early", "late")) %>%
    group_by(period) %>% summarise(
      n=n(), n_both=sum(has_both), pct=mean(has_both)*100, .groups="drop")
  
  cat("iuc+sat co-occurrence in C3:\n")
  print(as.data.frame(both_temporal), row.names=FALSE)
}

# ---- 6. Gene-level temporal in C3 (VFDB) ----
cat("\n=== 6. Gene-level temporal in C3 (VFDB) ===\n")
vfdb <- fread(config$st_vfdb_summary("ST69"), sep="\t", header=TRUE, data.table=FALSE) %>%
  filter(str_detect(.data[["#FILE"]],"^ST69/")) %>%
  mutate(genome_id=str_remove(str_remove(.data[["#FILE"]],"^ST69/"),"_vfdb.tsv$"))

to_bin <- function(x) as.integer(!is.na(x) & x != ".")
vf_genes <- intersect(interest, setdiff(colnames(vfdb), c("#FILE","NUM_FOUND","genome_id")))
v_bin <- vfdb %>% select(genome_id, all_of(vf_genes)) %>%
  mutate(across(all_of(vf_genes), to_bin, .names="{.col}_v"))

df <- v_bin %>% inner_join(master, by="genome_id") %>%
  filter(shell_cluster=="Cluster_3", year %in% c(EARLY,LATE)) %>%
  mutate(period=ifelse(year%in%EARLY,"early","late"))

gene_delta <- map_dfr(vf_genes, ~{
  gv <- paste0(.x, "_v")
  d <- df %>% group_by(period) %>% summarise(prev=mean(.data[[gv]],na.rm=T)*100, .groups="drop")
  if (nrow(d)<2) return(NULL)
  tibble(gene=.x, system=gene_sys[.x], early=d$prev[d$period=="early"],
    late=d$prev[d$period=="late"], delta=round(late-early,1))
}) %>% filter(!is.na(delta))

cat("Per-gene temporal change (|delta| >= 3pp):\n")
print(as.data.frame(gene_delta %>% filter(abs(delta)>=3) %>% arrange(system, desc(delta))), row.names=FALSE)

# System aggregate
sys_delta <- gene_delta %>% filter(!is.na(system)) %>%
  group_by(system) %>%
  summarise(n_genes=n(), mean_delta=round(mean(delta),1),
    range=paste0(round(min(delta),1)," to ",round(max(delta),1)), .groups="drop")
cat("\nSystem-level:\n")
print(as.data.frame(sys_delta), row.names=FALSE)

# ---- 7. Plot: RGP composite profile heatmap ----
cat("\n=== 7. Saving plots ===\n")
# Top RGP profiles bar chart
p1 <- profile_count %>% slice_max(n, n=15) %>%
  ggplot(aes(x=reorder(systems_present, n), y=n)) +
  geom_col(fill="steelblue") + coord_flip() +
  labs(title="RGP profiles carrying interest genes", x="systems in RGP", y="RGPs") +
  theme_minimal(base_size=9)
ggsave(file.path(OUT, "rgp_profiles.pdf"), p1, width=8, height=5)

# ---- 8. Summary text ----
sink(file.path(OUT, "rgp_context_summary.txt"))
cat("RGP Genomic Context Summary\n")
cat("===========================\n\n")
cat("Question: Are pap, kps, sat, and iuc genes physically linked?\n\n")

cat("Module organization:\n")
cat("  kps: modules 59 (kpsD,F), 79 (kpsD,M,T), 122,186,376 (kpsM,T)\n")
cat("  pap: modules 72 (papA-H,sat), 511,1003,981 (papB), 441,1071,983 (papF), etc.\n")
cat("  iuc: module 66 (iucA-D,iutA)\n")
cat("  sat: module 72 (shared with pap — co-transcribed)\n")
cat("  => pap+sat share module_72 — they are in the same functional unit\n")
cat("  => kps and iuc are in separate modules\n\n")

cat("RGP co-localization:\n")
cat(sprintf("  RGPs carrying interest genes: %d\n", nrow(rgp_sys)))
cat(sprintf("  RGPs with >=2 interest systems: %d\n", nrow(multi)))
cat(sprintf("  RGPs with pap+iuc: %d\n",
  sum(str_detect(rgp_sys$systems_present,"pap") & str_detect(rgp_sys$systems_present,"iuc"))))
cat(sprintf("  RGPs with pap+kps: %d\n",
  sum(str_detect(rgp_sys$systems_present,"pap") & str_detect(rgp_sys$systems_present,"kps"))))
cat(sprintf("  RGPs with kps+iuc: %d\n",
  sum(str_detect(rgp_sys$systems_present,"kps") & str_detect(rgp_sys$systems_present,"iuc"))))
cat(sprintf("  RGPs with all 4 (pap+kps+iuc+sat): %d\n",
  sum(str_detect(rgp_sys$systems_present,"pap") & str_detect(rgp_sys$systems_present,"kps") &
      str_detect(rgp_sys$systems_present,"iuc") & str_detect(rgp_sys$systems_present,"sat"))))

cat("\nKey finding:\n")
cat("  pap and sat share module_72 -> same functional module (co-transcribed PAI)\n")
cat("  iuc (module_66) and pap+sat (module_72) are INDEPENDENT modules\n")
cat("  but their families CO-LOCALIZE in the same RGPs in many genomes\n")
cat("  => They form a composite genomic island with multiple independently-acquired modules\n\n")

cat("C3 temporal changes:\n")
for (i in 1:nrow(sys_delta)) {
  r <- sys_delta[i,]
  cat(sprintf("  %s: mean %.1fpp (%s)\n", r$system, r$mean_delta, r$range))
}

cat("\nInterpretation for reviewer:\n")
cat("- pap+sat are in the same functional module (module_72) -> physically linked\n")
cat("- iuc (aerobactin) is in a separate module (66) but commonly co-localizes\n")
cat("  in the same RGP as pap+sat -> these form a composite pathogenicity island\n")
cat("- kps is independent but can be found in the same RGP in some genomes\n")
cat("- The temporal increases in C3 affect MULTIPLE independent modules simultaneously\n")
cat("  suggesting a shared genomic island is being gained/lost as a unit\n")
cat("- However, iuc+sat without pap and kps+iuc+sat combinations also exist,\n")
cat("  indicating the island can be partially acquired or variably assembled\n")
cat("- Long-read sequencing would be needed to confirm the exact structure\n")
cat("  and integration site of this composite island\n")
sink()
# ============================================================
# === MGE context: per-genome system profiles, RGP sharing  ===
# ============================================================

cat("\n\n=== MGE CONTEXT ANALYSIS ===\n")
cat("Question: Is a single mobile genetic element driving the co-occurrence?\n\n")

# ---- A. Per-genome system profiles in C3 ----
cat("--- A. Per-genome system profiles in C3 ---\n")

# Define systems from interest genes
system_genes <- list(
  kps = c("kpsD","kpsF","kpsM","kpsT","kpsC","kpsE","kpsS","kpsU"),
  pap = c("papA","papB","papC","papD","papE","papF","papG","papH","papI","papJ","papK","papX"),
  sat = "sat",
  iuc = c("iucA","iucB","iucC","iucD","iutA")
)

# Use VFDB binary data for system presence per genome
# System present if >= 50% of its genes are detected
calc_sys_prev <- function(df, sys_genes) {
  map_dfr(names(sys_genes), ~{
    gs <- intersect(sys_genes[[.x]], colnames(df))
    if (length(gs) == 0) return(NULL)
    vals <- df[[paste0(gs[1], "_v")]]
    if (length(gs) > 1) {
      vals <- rowSums(df[, paste0(gs, "_v"), drop=FALSE], na.rm=TRUE)
    }
    has_sys <- vals >= ceiling(length(gs) * 0.5)
    tibble(genome = df$genome_id, system = .x, present = has_sys)
  })
}

sys_matrix <- calc_sys_prev(v_bin, system_genes)

# Pivot to profile per genome
sys_profiles <- sys_matrix %>%
  filter(!is.na(present)) %>%
  group_by(genome) %>%
  summarise(profile = paste(sort(unique(system[present])), collapse=","),
            n_systems = sum(present), .groups="drop") %>%
  mutate(profile = ifelse(n_systems == 0, "none", profile))

# Add metadata
c3_sys <- sys_profiles %>%
  left_join(master %>% select(genome_id, year, shell_cluster), by=c("genome"="genome_id")) %>%
  filter(shell_cluster == "Cluster_3", !is.na(year)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late"))

# Profile counts in C3 overall
cat("\nC3 system profiles (all years):\n")
c3_profiles <- c3_sys %>% count(profile, sort=TRUE) %>% mutate(pct = round(n/sum(n)*100, 1))
print(as.data.frame(c3_profiles), row.names=FALSE)

# Profile counts by period (early vs late)
cat("\nC3 profile distribution by period:\n")
c3_profile_time <- c3_sys %>%
  group_by(period, profile) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(period) %>%
  mutate(pct = round(n/sum(n)*100, 1)) %>%
  arrange(period, desc(n))
print(as.data.frame(c3_profile_time), row.names=FALSE)

# Which profiles changed most?
cat("\nProfile delta (late - early pp):\n")
profile_delta <- c3_profile_time %>%
  select(period, profile, pct) %>%
  pivot_wider(names_from=period, values_from=pct, values_fill=0) %>%
  mutate(delta = round(late - early, 1)) %>%
  arrange(desc(abs(delta)))
print(as.data.frame(profile_delta), row.names=FALSE)

# ---- B. System co-occurrence in C3 ----
cat("\n--- B. System co-occurrence in C3 ---\n")

# Use sys_matrix directly (cleaner)
c3_sys_mat <- sys_matrix %>%
  filter(genome %in% c3_sys$genome) %>%
  left_join(master %>% select(genome_id, year, shell_cluster), by=c("genome"="genome_id")) %>%
  filter(shell_cluster == "Cluster_3", !is.na(year)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", "late"))

# Pairwise co-occurrence
pair_cooc <- map_dfr(combn(names(system_genes), 2, simplify=FALSE), ~{
  s1 <- .x[1]; s2 <- .x[2]
  d <- c3_sys_mat %>% filter(system %in% c(s1, s2)) %>%
    distinct(genome, system, present) %>%
    pivot_wider(names_from=system, values_from=present)
  n_both <- sum(d[[s1]] & d[[s2]], na.rm=TRUE)
  n_total <- nrow(d)
  tibble(sys1=s1, sys2=s2, n_both=n_both, n=n_total, pct=round(n_both/n_total*100, 1))
})
cat("Pairwise co-occurrence in C3:\n")
print(as.data.frame(pair_cooc), row.names=FALSE)

# ---- C. Module-RGP co-localization (from modules_RGP_lists.tsv) ----
cat("\n--- C. Module-RGP co-localization ---\n")
mod_rgp <- fread(file.path(config$PPANGGOLIN_DIR, "modules_RGP_lists.tsv"),
  sep="\t", header=TRUE, data.table=FALSE)

# Interest modules
interest_mods_all <- unique(mod_map$module_id)
cat("Interest modules:", paste(sort(interest_mods_all), collapse=", "), "\n")

# Filter to interest modules
mod_rgp_interest <- mod_rgp %>%
  filter(str_detect(mod_list, paste(interest_mods_all, collapse="|")))

# For each module, count RGPs
mod_rgp_count <- mod_rgp_interest %>%
  separate_rows(mod_list, sep=",") %>%
  filter(mod_list %in% interest_mods_all) %>%
  group_by(mod_list) %>%
  summarise(n_RGPs = n(), .groups="drop") %>%
  mutate(system = case_when(
    mod_list %in% c("module_66") ~ "iuc",
    mod_list %in% c("module_72") ~ "pap+sat",
    mod_list %in% c("module_59","module_79","module_122","module_186","module_376") ~ "kps",
    mod_list %in% c("module_441","module_511","module_768","module_981","module_983","module_1003","module_1071") ~ "pap",
    TRUE ~ "mixed"
  ))
cat("Module-RGP counts:\n")
print(as.data.frame(mod_rgp_count), row.names=FALSE)

# Do any RGPs carry MULTIPLE interest modules?
multi_mod_rgp <- mod_rgp_interest %>%
  separate_rows(mod_list, sep=",") %>%
  filter(mod_list %in% interest_mods_all) %>%
  group_by(representative_RGP) %>%
  summarise(
    n_mods = n_distinct(mod_list),
    mods = paste(sort(unique(mod_list)), collapse=","),
    n_all_rgps = n_distinct(unlist(strsplit(RGP_list, ","))),
    .groups="drop"
  ) %>%
  filter(n_mods >= 2)

cat("\nRGPs with >=2 interest modules:\n")
cat(sprintf("  Total: %d\n", nrow(multi_mod_rgp)))
cat("  Breakdown by module combo:\n")
multi_mod_summary <- multi_mod_rgp %>%
  count(mods, sort=TRUE) %>% head(10)
print(as.data.frame(multi_mod_summary), row.names=FALSE)

# Key: Do module_66 (iuc) and module_72 (pap+sat) share RGPs?
iuc72_rgp <- multi_mod_rgp %>%
  filter(str_detect(mods, "module_66") & str_detect(mods, "module_72"))
cat(sprintf("\n  RGPs with BOTH module_66 (iuc) AND module_72 (pap+sat): %d\n", nrow(iuc72_rgp)))

# Do ANY kps modules share RGPs with iuc or pap?
kps_mods <- c("module_59","module_79","module_122","module_186","module_376")
for (km in kps_mods) {
  n_iuc <- sum(str_detect(multi_mod_rgp$mods, km) & str_detect(multi_mod_rgp$mods, "module_66"))
  n_pap <- sum(str_detect(multi_mod_rgp$mods, km) & str_detect(multi_mod_rgp$mods, "module_72"))
  if (n_iuc > 0 || n_pap > 0)
    cat(sprintf("  %s shares RGPs with iuc:%d, pap+sat:%d\n", km, n_iuc, n_pap))
}

# ---- D. RGP cluster sharing ----
cat("\n--- D. RGP cluster sharing (shared homologous RGPs across genomes) ---\n")
rgp_cl <- fread(file.path(config$PANGENOME_DIR, "rgp_clusters", "rgp_cluster.tsv"),
  sep="\t", header=TRUE, data.table=FALSE)

# For the RGPs carrying composite profiles, which clusters are they in?
# First, get the list of RGPs for each interest module
mod_rgp_all <- mod_rgp_interest %>%
  separate_rows(RGP_list, sep=",") %>%
  filter(str_detect(mod_list, paste(interest_mods_all, collapse="|"))) %>%
  distinct(representative_RGP, mod_list, RGP = RGP_list)

# Join with rgp clusters
mod_rgp_cluster <- mod_rgp_all %>%
  left_join(rgp_cl, by=c("RGP"="RGPs"), relationship="many-to-many")

# For each module, which clusters
cluster_mod_summary <- mod_rgp_cluster %>%
  filter(!is.na(cluster)) %>%
  group_by(mod_list, cluster) %>%
  summarise(n_genomes = n(), .groups="drop") %>%
  mutate(system = case_when(
    mod_list == "module_66" ~ "iuc",
    mod_list == "module_72" ~ "pap+sat",
    TRUE ~ "other"
  ))

cat("\nTop clusters for module_66 (iuc):\n")
print(as.data.frame(cluster_mod_summary %>% filter(mod_list=="module_66") %>% arrange(desc(n_genomes)) %>% head(10)), row.names=FALSE)
cat("Top clusters for module_72 (pap+sat):\n")
print(as.data.frame(cluster_mod_summary %>% filter(mod_list=="module_72") %>% arrange(desc(n_genomes)) %>% head(10)), row.names=FALSE)

# Do module_66 and module_72 share the SAME clusters?
shared_clusters <- cluster_mod_summary %>%
  filter(mod_list %in% c("module_66","module_72")) %>%
  group_by(cluster) %>%
  summarise(n_mods = n_distinct(mod_list),
            mods = paste(sort(unique(mod_list)), collapse=","),
            n_total = sum(n_genomes),
            .groups="drop") %>%
  filter(n_mods >= 2)
cat(sprintf("\nClusters shared by module_66 AND module_72: %d\n", nrow(shared_clusters)))
if (nrow(shared_clusters) > 0) {
  print(as.data.frame(shared_clusters %>% arrange(desc(n_total))), row.names=FALSE)
}

# ---- E. Kps operon progression (VFDB gene-level) ----
cat("\n--- E. Kps operon progression in C3 (VFDB) ---\n")

# All kps genes
kps_genes_all <- c("kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT","kpsU")
kps_vars <- paste0(kps_genes_all, "_v")

# For each C3 genome, count kps genes present
c3_kps_vf <- df %>%
  filter(shell_cluster == "Cluster_3") %>%
  select(genome_id, year, period, all_of(kps_vars)) %>%
  mutate(
    n_kps = rowSums(across(all_of(kps_vars))),
    kps_complete = n_kps >= 8,
    kps_partial = n_kps >= 1 & n_kps < 8,
    kps_none = n_kps == 0,
    has_mt = .data[["kpsM_v"]] == 1 & .data[["kpsT_v"]] == 1,
    has_core = .data[["kpsC_v"]] == 1 & .data[["kpsD_v"]] == 1 &
                .data[["kpsE_v"]] == 1 & .data[["kpsF_v"]] == 1
  )

c3_kps_prog <- c3_kps_vf %>%
  group_by(period) %>%
  summarise(
    n = n(),
    complete = sum(kps_complete),
    partial = sum(kps_partial),
    none = sum(kps_none),
    complete_pct = round(mean(kps_complete)*100, 1),
    partial_pct = round(mean(kps_partial)*100, 1),
    none_pct = round(mean(kps_none)*100, 1),
    mean_kps_genes = round(mean(n_kps), 1),
    has_mt = round(mean(has_mt)*100, 1),
    has_core = round(mean(has_core)*100, 1),
    .groups="drop"
  )
cat("Kps operon progression in C3:\n")
print(as.data.frame(c3_kps_prog), row.names=FALSE)

# Also: what is the distribution of kps gene counts?
cat("\nKps gene count distribution in C3:\n")
c3_kps_dist <- c3_kps_vf %>%
  group_by(period, n_kps) %>%
  summarise(count=n(), .groups="drop") %>%
  group_by(period) %>% mutate(pct=round(count/sum(count)*100, 1))
print(as.data.frame(c3_kps_dist), row.names=FALSE)

# ---- F. Profile summary for single MGE hypothesis ----
cat("\n\n--- SYNTHESIS: Single MGE hypothesis ---\n")

# Evidence for single MGE:
# 1. Systems co-occurrence rates
cat("1. Pairwise co-occurrence in C3 genomes:\n")
for (i in 1:nrow(pair_cooc)) {
  r <- pair_cooc[i,]
  cat(sprintf("   %s + %s: %.1f%% (%d/%d)\n", r$sys1, r$sys2, r$pct, r$n_both, r$n))
}

# 2. Multi-module RGPs
cat(sprintf("2. RGPs carrying multiple interest modules: %d\n", nrow(multi_mod_rgp)))
cat(sprintf("   module_66 (iuc) + module_72 (pap+sat) in same RGP: %d\n", nrow(iuc72_rgp)))
cat(sprintf("   module_66 (iuc) + module_72 (pap+sat) share same RGP cluster: %d\n", nrow(shared_clusters)))

# 3. C3 temporal profile trajectory
cat("3. C3 profile trajectory:\n")
for (i in 1:nrow(profile_delta)) {
  r <- profile_delta[i,]
  cat(sprintf("   %s: early=%.1f%%  late=%.1f%%  delta=%.1fpp\n",
              r$profile, r$early, r$late, r$delta))
}

# Conclusion text
sink(file.path(OUT, "mge_context_summary.txt"))
cat("MGE Context Summary\n")
cat("===================\n\n")
cat("Question: Is a single mobile genetic element driving the observation\n")
cat("that kps, pap, sat, and iuc all increase simultaneously in C3?\n\n")

cat("PER-GENOME SYSTEM PROFILES IN C3:\n")
for (i in 1:nrow(c3_profiles)) {
  r <- c3_profiles[i,]
  cat(sprintf("  %s: n=%d (%.1f%%)\n", r$profile, r$n, r$pct))
}

cat("\nTEMPORAL PROFILE CHANGE:\n")
for (i in 1:nrow(profile_delta)) {
  r <- profile_delta[i,]
  cat(sprintf("  %s: %.1f%% -> %.1f%% (%.1fpp)\n", r$profile, r$early, r$late, r$delta))
}

cat("\nPAIRWISE CO-OCCURRENCE IN C3:\n")
for (i in 1:nrow(pair_cooc)) {
  r <- pair_cooc[i,]
  cat(sprintf("  %s + %s: %.1f%%\n", r$sys1, r$sys2, r$pct))
}

cat("\nMODULE-RGP CO-LOCALIZATION:\n")
cat(sprintf("  module_66 (iuc) + module_72 (pap+sat) in same RGP: %d RGPs\n", nrow(iuc72_rgp)))
cat(sprintf("  module_66 + module_72 share same RGP cluster: %d clusters\n", nrow(shared_clusters)))

cat("\nKPS OPERON PROGRESSION IN C3 (VFDB):\n")
for (i in 1:nrow(c3_kps_prog)) {
  r <- c3_kps_prog[i,]
  cat(sprintf("  %s: complete=%.1f%% partial=%.1f%% none=%.1f%% (mean %.1f kps genes)\n",
              r$period, r$complete_pct, r$partial_pct, r$none_pct, r$mean_kps_genes))
  cat(sprintf("        has kpsM+T=%.1f%% has core(kpsC,D,E,F)=%.1f%%\n", r$has_mt, r$has_core))
}

cat("\nINTERPRETATION:\n")
cat("- pap+sat share module_72 -> the same functional unit, co-transcribed\n")
cat("- iuc (module_66) is in a SEPARATE module but co-localizes in\n")
cat("  the same RGPs as module_72 in many genomes\n")
cat(sprintf("  => %d RGPs carry BOTH modules (%d share same cluster_775 = shared element)\n",
    nrow(iuc72_rgp), nrow(shared_clusters)))
cat("- kps is in separate modules and rarely co-localizes with iuc/pap\n")
cat("  in the same RGP (only 5 RGPs carry all 4 systems)\n")
cat("- kps operon completion: C3 genomes go from mean ~3.4 to ~5.3 kps\n")
cat("  genes, with the full 8-gene operon increasing from 0% to 10.6%\n\n")

cat("CONCLUSION FOR REVIEWER:\n")
cat("The data are most consistent with TWO independent MGE events:\n\n")
cat("1. A composite genomic island carrying iuc+pap+sat (and sometimes\n")
cat("   additional pap genes) that is increasingly acquired in C3 genomes.\n")
cat("   iuc (aerobactin) and pap/sat (P fimbriae + SAT toxin) are in\n")
cat("   distinct functional modules but frequently co-localize in the\n")
cat("   same RGPs, suggesting they form a composite PAI that can be\n")
cat("   acquired as a unit or partially assembled.\n")
cat("   The shared RGP cluster_775 (406 iuc + 304 pap+sat genomes)\n")
cat("   confirms a homologous element is spreading across genomes.\n\n")
cat("2. Independent kps operon completion within existing G2 capsule\n")
cat("   backgrounds. The kps genes are rarely in the same RGP as\n")
cat("   iuc/pap/sat, and their temporal increase reflects completion\n")
cat("   of the capsule operon (adding kpsC,D,E,F,M,S,T,U) rather than\n")
cat("   acquisition of a new island.\n\n")
cat("The simultaneous increase in C3 is best explained by selection\n")
cat("for a composite virulence plasmid or integrative element carrying\n")
cat("iuc+pap+sat that is sweeping through the population, while kps\n")
cat("operon completion is a separate, concurrent process that happens\n")
cat("to show the same temporal trend.\n\n")
cat("Long-read sequencing (e.g., Lipworth et al. 2024, Arredondo-Alonso\n")
cat("2025) would be needed to resolve the exact structure and\n")
cat("integration site of this composite element.\n")
sink()

# ---- G. Save figures ----
cat("\n--- Saving additional figures ---\n")

# Figure: C3 profile distribution over time
p_profiles <- c3_sys %>%
  mutate(profile = fct_infreq(profile) %>% fct_rev()) %>%
  ggplot(aes(x=profile, fill=period)) +
  geom_bar(position="dodge") +
  coord_flip() +
  scale_fill_viridis_d() +
  labs(title="C3 system profiles over time",
       x="system profile", y="genomes") +
  theme_minimal(base_size=9)
ggsave(file.path(OUT, "c3_profiles_temporal.pdf"), p_profiles, width=8, height=5)

# Figure: Kps module progression
p_kps_prog <- c3_kps_prog %>%
  select(period, complete, partial, none) %>%
  pivot_longer(-period, names_to="kps_status", values_to="count") %>%
  mutate(kps_status = factor(kps_status, levels=c("none","partial","complete"))) %>%
  ggplot(aes(x=period, y=count, fill=kps_status)) +
  geom_col(position="fill") +
  scale_fill_manual(values=c("none"="grey80","partial"="goldenrod","complete"="steelblue")) +
  labs(title="Kps operon status in C3", y="proportion", x="") +
  theme_minimal(base_size=10)
ggsave(file.path(OUT, "c3_kps_progression.pdf"), p_kps_prog, width=5, height=4)

# Figure: Module co-occurrence heatmap
pairs <- pair_cooc %>%
  mutate(label = sprintf("%.0f%%", pct)) %>%
  ggplot(aes(x=sys1, y=sys2, fill=pct)) +
  geom_tile() +
  geom_text(aes(label=label), size=3) +
  scale_fill_viridis_c() +
  labs(title="System co-occurrence in C3 (%)", x="", y="") +
  theme_minimal(base_size=10)
ggsave(file.path(OUT, "c3_system_cooccurrence.pdf"), pairs, width=5, height=4)

cat("Done. All outputs in:", OUT, "\n")

})  # end block 10: 10c_rgp_genomic_context.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 11/17: 10e_allelic_conversion.R
# Allelic conversion (Supplementary Figure S2; APPENDS to capsule_classification.xlsx chain -- run after 09b/09c/09d)
# (verbatim body of original scripts/.../10e_allelic_conversion.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# RIGOROUS allele replacement analysis:
# - Uses BLAST-based mapping ONLY (≥90% identity for high confidence)
# - Separates low-confidence annotations from existing mapping
# - Handles: 1) each PP family maps to ONE best VFDB gene
# - For each VFDB gene: compute temporal dynamics using HIGH-CONF families only
# - Then check if results change when including LOW-CONF families

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(writexl); library(readxl)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_capsule_classification")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
EARLY <- c(2016,2017,2018); LATE <- c(2022,2023,2024,2025)

cat(strrep("=", 60), "\n")
cat("RIGOROUS allele replacement analysis\n")
cat(strrep("=", 60), "\n")

# === 1. Load VFDB ===
cat("1. Loading VFDB...\n")
vfdb <- fread(config$st_vfdb_summary("ST69"), sep="\t", header=TRUE, data.table=FALSE)
vfdb <- vfdb %>% filter(str_detect(.data[["#FILE"]],"^ST69/")) %>%
  mutate(genome_id=str_remove(str_remove(.data[["#FILE"]],"^ST69/"),"_vfdb.tsv$"))
to_bin <- function(x) as.integer(!is.na(x) & x != ".")
all_vf <- setdiff(colnames(vfdb), c("#FILE","NUM_FOUND","genome_id"))
meta <- fread("output/ST69/vfdb_analysis/04_master_shell_cluster_metadata_VFDB_table.csv",
  sep=",", header=TRUE, data.table=FALSE)

v_bin <- vfdb %>% select(genome_id, all_of(all_vf)) %>%
  mutate(across(all_of(all_vf), to_bin, .names="{.col}_v"))
df <- v_bin %>% inner_join(meta, by="genome_id") %>% mutate(year=as.integer(year)) %>%
  filter(!is.na(year))
c3 <- df %>% filter(shell_cluster=="Cluster_3", year %in% c(EARLY,LATE)) %>%
  mutate(period=ifelse(year%in%EARLY,"early","late"))

vfdb_deltas <- map_dfr(all_vf, ~{
  gv <- paste0(.x, "_v")
  if (!gv %in% colnames(c3)) return(NULL)
  d <- c3 %>% group_by(period) %>% summarise(prev=mean(.data[[gv]],na.rm=T)*100, .groups="drop")
  if (nrow(d) < 2) return(NULL)
  tibble(gene=.x, vfdb_early=d$prev[d$period=="early"], vfdb_late=d$prev[d$period=="late"],
    vfdb_delta=round(vfdb_late-vfdb_early,1))
}) %>% filter(!is.na(vfdb_delta))
cat("  VFDB genes in C3:", nrow(vfdb_deltas), "\n")

# === 2. Load PP ===
cat("2. Loading PPanGGOLiN...\n")
rtab <- fread(file.path(config$PANGENOME_DIR,"gene_presence_absence.Rtab"),
  sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
rtab_fams <- rtab[[1]]; rtab_mat <- as.matrix(rtab[,-1])
rownames(rtab_mat) <- rtab_fams; genome_ids <- colnames(rtab)[-1]

# === 3. Build BLAST-based mapping ===
cat("3. Building BLAST-based mapping...\n")
blast <- fread("output/vfdb_vs_pp_blast.txt", sep="\t", header=FALSE,
  col.names=c("vfdb_gene","pp_family","pident","length","mismatch","gapopen",
              "qstart","qend","sstart","send","evalue","bitscore"))

# For each PP family, find its best BLAST hit
# Then classify as HIGH (>=90%) or LOW (<90% but >=80%) confidence
best_hit <- blast %>%
  group_by(pp_family) %>%
  slice_max(pident, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  mutate(confidence = case_when(
    pident >= 90 ~ "high",
    pident >= 80 ~ "low",
    TRUE ~ "poor"
  ))

cat("  PP families with BLAST hits:", nrow(best_hit), "\n")
cat("  High confidence (>=90%):", sum(best_hit$confidence=="high"), "\n")
cat("  Low confidence (80-90%):", sum(best_hit$confidence=="low"), "\n")

# Group families by their best-matching VFDB gene
high_map <- best_hit %>% filter(confidence=="high") %>%
  group_by(vfdb_gene) %>% summarise(families=list(pp_family), .groups="drop")

# === 4. Temporal comparison using ONLY high-confidence mapping ===
cat("4. Computing temporal comparison (high-confidence only)...\n")

results <- map_dfr(vfdb_deltas$gene, ~{
  # Find high-confidence families for this gene
  row <- high_map %>% filter(vfdb_gene == .x)
  fams <- unlist(row$families)
  if (length(fams) == 0) {
    return(tibble(gene=.x, hi_pp_early=NA, hi_pp_late=NA, hi_pp_delta=NA,
                  lo_pp_early=NA, lo_pp_late=NA, lo_pp_delta=NA,
                  n_hi=0, n_lo=0, verdict="No PP data"))
  }
  # High-confidence families
  combined <- colSums(rtab_mat[fams[fams %in% rownames(rtab_mat)],,drop=FALSE] > 0, na.rm=TRUE) > 0
  hi_tib <- tibble(genome_id=genome_ids, present=as.integer(combined)) %>%
    inner_join(meta, by="genome_id") %>% mutate(year=as.integer(year)) %>%
    filter(shell_cluster=="Cluster_3", year %in% c(EARLY,LATE)) %>%
    mutate(period=ifelse(year%in%EARLY,"early","late")) %>%
    group_by(period) %>% summarise(prev=mean(present)*100, .groups="drop")
  
  hi_early <- hi_tib$prev[hi_tib$period=="early"]
  hi_late <- hi_tib$prev[hi_tib$period=="late"]
  hi_delta <- round(hi_late - hi_early, 1)
  
  # Also check if there are any low-confidence families
  low_row <- best_hit %>% filter(vfdb_gene == .x, confidence=="low")
  lo_fams <- setdiff(low_row$pp_family, fams)
  lo_delta <- NA; lo_early <- NA; lo_late <- NA
  
  if (length(lo_fams) > 0) {
    lo_fams <- lo_fams[lo_fams %in% rownames(rtab_mat)]
    if (length(lo_fams) > 0) {
      lo_comb <- colSums(rtab_mat[lo_fams,,drop=FALSE] > 0, na.rm=TRUE) > 0
      lo_tib <- tibble(genome_id=genome_ids, present=as.integer(lo_comb)) %>%
        inner_join(meta, by="genome_id") %>% mutate(year=as.integer(year)) %>%
        filter(shell_cluster=="Cluster_3", year %in% c(EARLY,LATE)) %>%
        mutate(period=ifelse(year%in%EARLY,"early","late")) %>%
        group_by(period) %>% summarise(prev=mean(present)*100, .groups="drop")
      lo_early <- lo_tib$prev[lo_tib$period=="early"]
      lo_late <- lo_tib$prev[lo_tib$period=="late"]
      lo_delta <- round(lo_late - lo_early, 1)
    }
  }
  
  # Verdict based on high-confidence families
  vd <- vfdb_deltas$vfdb_delta[vfdb_deltas$gene==.x]
  if (abs(vd - hi_delta) < 5) {
    verdict <- "REAL"
  } else if (vd > hi_delta + 5 & hi_delta < 5) {
    verdict <- "Allele replacement"
  } else if (vd > hi_delta + 5 & hi_delta >= 5) {
    verdict <- "Partial (both real + replacement)"
  } else if (hi_delta > vd + 5) {
    verdict <- "PP only increase"
  } else {
    verdict <- "Mixed"
  }
  
  tibble(gene=.x,
    vfdb_early=vfdb_deltas$vfdb_early[vfdb_deltas$gene==.x],
    vfdb_late=vfdb_deltas$vfdb_late[vfdb_deltas$gene==.x],
    vfdb_delta=vd,
    n_hi=length(fams), n_lo=length(lo_fams),
    hi_pp_early=hi_early, hi_pp_late=hi_late, hi_pp_delta=hi_delta,
    lo_pp_early=lo_early, lo_pp_late=lo_late, lo_pp_delta=lo_delta,
    verdict=verdict)
})

# === 5. Summary ===
cat("\n=== SUMMARY (high-confidence BLAST mapping only) ===\n")
cat("Total VF genes tested:", nrow(results), "\n")
for (v in c("REAL", "Allele replacement", "Partial (both real + replacement)",
           "PP only increase", "Mixed", "No PP data")) {
  n <- sum(results$verdict == v)
  p <- round(n / nrow(results) * 100, 1)
  if (n > 0) cat(sprintf("  %s: %d (%.1f%%)\n", v, n, p))
}

# Show genes with meaningful change (|VFDB delta| >= 5pp)
cat("\n=== GENES WITH |VFDB delta| >= 5pp ===\n")
meaningful <- results %>% filter(abs(vfdb_delta) >= 5) %>% arrange(desc(vfdb_delta))
for (i in 1:nrow(meaningful)) {
  r <- meaningful[i,]
  hid <- ifelse(is.na(r$hi_pp_delta), "NA", sprintf("%+.1f", r$hi_pp_delta))
  cat(sprintf("  %-12s VFDB=%+.1fpp  PP(high)=%spp  %s\n",
    r$gene, r$vfdb_delta, hid, r$verdict))
}

# Detail: flagged genes
cat("\n=== ALLELE REPLACEMENT FLAGGED ===\n")
flagged <- results %>% filter(verdict=="Allele replacement") %>% arrange(desc(vfdb_delta))
for (i in 1:nrow(flagged)) {
  r <- flagged[i,]
  hid <- ifelse(is.na(r$hi_pp_delta), "NA", sprintf("%+.1f", r$hi_pp_delta))
  cat(sprintf("  %-12s VFDB=%+.1fpp  PP(high)=%spp  high_fams=%d  low_fams=%d\n",
    r$gene, r$vfdb_delta, hid, r$n_hi, r$n_lo))
}

# === 6. Plot ===
cat("\n6. Plotting...\n")
plot_df <- results %>% filter(!is.na(hi_pp_delta)) %>%
  mutate(label = case_when(
    verdict == "REAL" ~ "Real increase",
    verdict == "Allele replacement" ~ "Allele replacement",
    TRUE ~ "Other"),
    label = factor(label, levels=c("Real increase","Allele replacement","Other")))

p <- ggplot(plot_df, aes(x=vfdb_delta, y=hi_pp_delta, color=label)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  geom_hline(yintercept=0, linetype="dotted", color="grey70") +
  geom_vline(xintercept=0, linetype="dotted", color="grey70") +
  geom_point(size=3, alpha=0.8) +
  ggrepel::geom_text_repel(aes(label=gene), size=3, max.overlaps=20) +
  scale_color_manual(values=c("Real increase"="#1b7837",
    "Allele replacement"="#d73027", "Other"="grey60")) +
  labs(x="VFDB delta (pp)", y="PPanGGOLiN delta (pp)",
    title=paste0("Allelic Conversion (>=90% identity, n=", nrow(plot_df), ")"),
    subtitle="Diagonal = agreement; below = VFDB overestimates", color="") +
  coord_fixed(xlim=c(-10,50), ylim=c(-40,50)) +
  theme_minimal(base_size=10)
ggsave(file.path(OUT, "VFDB_vs_PP_rigorous.pdf"), p, width=8, height=7)

# === 7. Detailed breakdown for key genes ===
cat("7. Key gene per-family breakdowns...\n")
key_genes <- c("kpsD","kpsM","kpsT","kpsF","kpsC","kpsE","kpsS","kpsU",
               "papB","papF","fepA")

for (gene in key_genes) {
  # Get BLAST-based families (high + low confidence)
  hi_fams <- best_hit %>% filter(vfdb_gene==gene, confidence=="high") %>% pull(pp_family)
  lo_fams <- best_hit %>% filter(vfdb_gene==gene, confidence=="low") %>% pull(pp_family)
  all_fams <- unique(c(hi_fams, lo_fams))
  all_fams <- all_fams[all_fams %in% rownames(rtab_mat)]
  if (length(all_fams) < 1) next
  
  # Skip genes with no temporal change
  rd <- results %>% filter(gene==!!gene)
  if (nrow(rd) == 0 || (abs(rd$vfdb_delta) < 3 && abs(rd$hi_pp_delta) < 3)) next
  
  cat("  ", gene, "...\n")
  
  # Per-family temporal
  fam_data <- map_dfr(all_fams, function(f) {
    pres <- rtab_mat[f,] > 0
    conf <- ifelse(f %in% hi_fams, "high", "low")
    tibble(genome_id=genome_ids, present=as.integer(pres), family=f, conf=conf)
  }) %>% inner_join(meta, by="genome_id") %>%
    filter(shell_cluster=="Cluster_3") %>% mutate(year=as.integer(year))
  
  # Yearly prevalence per family
  p <- fam_data %>% filter(!is.na(year)) %>%
    group_by(family, year, conf) %>% summarise(prev=mean(present)*100, .groups="drop") %>%
    mutate(family_short = str_extract(family, "[^_]+_[^_]+_[^_]+$"),
           label = paste0(family_short, " (", conf, ")")) %>%
    ggplot(aes(x=year, y=prev, color=label, group=family, linetype=conf)) +
    geom_line(size=1) + geom_point(size=2) +
    scale_x_continuous(breaks=2016:2025) +
    scale_linetype_manual(values=c("high"="solid","low"="dashed")) +
    labs(title=paste0(gene, " - PP family temporal dynamics (rigorous BLAST mapping)"),
      y="Prevalence (%)", x="Year", color="PP family", linetype="Confidence") +
    theme_minimal(base_size=10) +
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1))
  ggsave(file.path(OUT, paste0(gene, "_PP_families_rigorous.pdf")), p, width=8, height=5)
}

# === 8. Save ===
cat("8. Saving...\n")
xf <- file.path(OUT, "capsule_classification.xlsx")
all_sheets <- list()
for (sn in excel_sheets(xf)) { all_sheets[[sn]] <- read_xlsx(xf, sheet=sn) }
all_sheets[["rigorous_mapping"]] <- results
write_xlsx(all_sheets, xf)
write_csv(results, file.path(OUT, "VFDB_vs_PP_rigorous.csv"))
cat("Done.\n")

})  # end block 11: 10e_allelic_conversion.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 12/17: 10f_st10_arg_decomposition.R
# ST10 ARG decomposition
# (verbatim body of original scripts/.../10f_st10_arg_decomposition.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# Reviewer response: decompose ST10 ARG decreasing trend
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(writexl); library(data.table); library(broom)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST10", "reviewer_arg_decomposition")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
    grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
    TRUE ~ paste0("Escherichia_coli_", x))
}

# ---- 1. Load ST10 CARD data ----
cat("Loading ST10 CARD...\n")
card_raw <- fread(config$st_card_burden("ST10"), sep = "\t", header = TRUE,
  data.table = FALSE, colClasses = "character", na.strings = "")
card_raw$genome_id <- make_genome_id(str_remove(card_raw[[1]], "^card/ST10/"))
card_raw$genome_id <- str_remove(card_raw$genome_id, "_card\\.tsv$")

# Count ARGs per genome
arg_genes <- setdiff(colnames(card_raw), c("#FILE", "NUM_FOUND", "genome_id"))
to_bin <- function(x) as.integer(sapply(strsplit(as.character(x), ";"),
  function(v) any(suppressWarnings(as.numeric(v) >= 90), na.rm = TRUE)))
card_bin <- card_raw %>% select(genome_id, all_of(arg_genes)) %>%
  mutate(across(all_of(arg_genes), to_bin))
card_bin$total_arg <- rowSums(card_bin[, arg_genes, drop = FALSE], na.rm = TRUE)

cat("ST10 CARD genomes:", nrow(card_bin), "\n")

# ---- 2. Load ST10 metadata ----
cat("Loading ST10 metadata...\n")
meta_file <- config$st_metadata("ST10")
meta <- read_excel(meta_file, col_types = "text")

# Find key columns
name_col <- {
  candidates <- c("Name", "genome", "strain", "isolate", "assembly", "sample")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[2]
}
year_col <- {
  candidates <- c("Collection Year", "Collection_Year", "year", "Year", "date")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[7]
}
source_col <- {
  candidates <- c("Source Niche", "Source_Niche", "source_niche", "source", "niche", "Host", "host")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[4]
}
country_col <- {
  candidates <- c("Country", "country")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else colnames(meta)[12]
}
continent_col <- {
  candidates <- c("Continent", "continent")
  hits <- intersect(candidates, colnames(meta))
  if (length(hits) > 0) hits[1] else NA_character_
}

meta_clean <- meta %>%
  transmute(
    name_raw = .data[[name_col]],
    year = suppressWarnings(as.integer(.data[[year_col]])),
    source_niche = as.character(.data[[source_col]]),
    country = as.character(.data[[country_col]]),
    continent = if (!is.na(continent_col)) as.character(.data[[continent_col]]) else NA_character_
  ) %>%
  mutate(genome_id = make_genome_id(name_raw))

# ---- 3. Merge ----
merged <- card_bin %>% inner_join(meta_clean, by = "genome_id") %>%
  filter(!is.na(year), year >= 2016, year <= 2025)

cat("Merged genomes:", nrow(merged), "\n")

# Classify niche
merged <- merged %>%
  mutate(
    clinical = case_when(
      grepl("human|patient|clinical|blood|urine|wound|hospital|csf|respiratory|sputum",
        source_niche, ignore.case = TRUE) ~ "Clinical",
      is.na(source_niche) | source_niche == "" ~ "Unknown",
      TRUE ~ "Non-clinical"
    ),
    region = case_when(
      grepl("Europe|UK|Germany|France|Spain|Italy|Netherlands|Denmark|Sweden|Norway|Switzerland|Austria|Belgium|Ireland|Portugal|Poland|Czech|Hungary|Romania|Greece|Finland",
        country, ignore.case = TRUE) ~ "Europe",
      grepl("USA|United States|Canada|Mexico",
        country, ignore.case = TRUE) ~ "North America",
      grepl("China|Japan|India|South Korea|Taiwan|Thailand|Vietnam|Indonesia|Malaysia|Philippines|Singapore|Hong Kong|Pakistan|Bangladesh|Sri Lanka",
        country, ignore.case = TRUE) ~ "Asia",
      grepl("Australia|New Zealand",
        country, ignore.case = TRUE) ~ "Oceania",
      grepl("Brazil|Argentina|Chile|Colombia|Peru|Uruguay|Paraguay|Ecuador|Venezuela",
        country, ignore.case = TRUE) ~ "South America",
      grepl("South Africa|Nigeria|Kenya|Ghana|Ethiopia|Tanzania|Uganda|Zambia|Malawi|Mozambique|Botswana|Senegal|Burkina Faso|Benin|Zimbabwe",
        country, ignore.case = TRUE) ~ "Africa",
      TRUE ~ "Other"
    )
  )

# ---- 4. Annual ARG burden trends ----
annual <- merged %>%
  group_by(year) %>%
  summarise(
    n = n(),
    mean_arg = mean(total_arg),
    se_arg = sd(total_arg) / sqrt(n()),
    pct_clinical = mean(clinical == "Clinical") * 100,
    n_regions = n_distinct(region),
    .groups = "drop"
  )
cat("\nAnnual ARG burden:\n")
print(annual, n = Inf)

# ---- 5. Annual by clinical status ----
annual_clinical <- merged %>%
  filter(clinical != "Unknown") %>%
  group_by(year, clinical) %>%
  summarise(
    n = n(),
    mean_arg = mean(total_arg),
    .groups = "drop"
  )
cat("\nAnnual ARG by clinical status:\n")
print(annual_clinical %>% arrange(year, clinical), n = Inf)

# ---- 6. Annual by region ----
annual_region <- merged %>%
  filter(region != "Other") %>%
  group_by(year, region) %>%
  summarise(
    n = n(),
    mean_arg = mean(total_arg),
    .groups = "drop"
  )
cat("\nAnnual ARG by region:\n")
print(annual_region %>% arrange(year, region), n = Inf)

# ---- 7. Test: does clinical % decrease over time? ----
clinical_trend <- merged %>%
  filter(clinical != "Unknown") %>%
  group_by(year) %>%
  summarise(pct_clinical = mean(clinical == "Clinical") * 100, .groups = "drop")
m_clin <- lm(pct_clinical ~ year, data = clinical_trend)
cat("\nClinical % trend:\n")
cat(sprintf("  slope = %.4f, p = %.4e\n", coef(m_clin)[2], summary(m_clin)$coefficients[2,4]))

# ---- 8. Test: does regional composition shift over time? ----
regional_trend <- merged %>%
  filter(region != "Other") %>%
  count(year, region) %>%
  group_by(year) %>%
  mutate(pct = n / sum(n) * 100)

# ---- 9. Adjusted models ----
# Model 1: raw trend
m1 <- lm(total_arg ~ year, data = merged)
# Model 2: adjusted for clinical
m2 <- lm(total_arg ~ year + clinical, data = merged %>% filter(clinical != "Unknown"))
# Model 3: adjusted for region
m3 <- lm(total_arg ~ year + region, data = merged %>% filter(region != "Other"))
# Model 4: adjusted for both
m4 <- lm(total_arg ~ year + clinical + region, data = merged %>%
    filter(clinical != "Unknown", region != "Other"))

model_comp <- bind_rows(
  tidy(m1) %>% filter(term == "year") %>% mutate(model = "raw", dataset = "all"),
  tidy(m2) %>% filter(term == "year") %>% mutate(model = "+ clinical", dataset = "known clinical"),
  tidy(m3) %>% filter(term == "year") %>% mutate(model = "+ region", dataset = "known region"),
  tidy(m4) %>% filter(term == "year") %>% mutate(model = "+ clinical + region", dataset = "both")
) %>% select(model, estimate, std.error, p.value)

cat("\nModel comparison (year coefficient):\n")
print(model_comp)

# ---- 10. Clinical-stratified trends ----
niche_stratified <- merged %>%
  filter(clinical != "Unknown") %>%
  group_by(clinical) %>%
  summarise(
    n = n(),
    slope = tryCatch(coef(lm(total_arg ~ year))[2], error = function(e) NA_real_),
    p = tryCatch(summary(lm(total_arg ~ year))$coefficients[2,4], error = function(e) NA_real_),
    mean_arg = mean(total_arg),
    .groups = "drop"
  )
cat("\nClinical-stratified trends:\n")
print(niche_stratified)

# ---- 11. Region-stratified trends ----
region_stratified <- merged %>%
  filter(region != "Other") %>%
  group_by(region) %>%
  summarise(
    n = n(),
    slope = tryCatch(coef(lm(total_arg ~ year))[2], error = function(e) NA_real_),
    p = tryCatch(summary(lm(total_arg ~ year))$coefficients[2,4], error = function(e) NA_real_),
    mean_arg = mean(total_arg),
    .groups = "drop"
  )
cat("\nRegion-stratified trends:\n")
print(region_stratified)

# ---- 12. Write outputs ----
write_xlsx(list(
  annual_arg = annual,
  annual_by_clinical = annual_clinical %>% arrange(year, clinical),
  annual_by_region = annual_region %>% arrange(year, region),
  clinical_trend = clinical_trend,
  model_comparison = model_comp,
  niche_stratified = niche_stratified,
  region_stratified = region_stratified,
  regional_composition = regional_trend %>% arrange(year, region)
), path = file.path(OUT, "ST10_arg_decomposition.xlsx"))

# ---- 13. Summary ----
sink(file.path(OUT, "ST10_arg_decomposition_summary.txt"))
cat("ST10 ARG Burden Decomposition\n")
cat("==============================\n\n")
cat(sprintf("Genomes: %d (2016-2025)\n", nrow(merged)))
cat(sprintf("Raw ARG trend slope: %.4f (p = %.4e)\n", coef(m1)[2], summary(m1)$coefficients[2,4]))
cat(sprintf("\nAfter adjusting for clinical status: slope = %.4f (p = %.4e)\n",
  filter(model_comp, model == "+ clinical")$estimate,
  filter(model_comp, model == "+ clinical")$p.value))
cat(sprintf("After adjusting for region: slope = %.4f (p = %.4e)\n",
  filter(model_comp, model == "+ region")$estimate,
  filter(model_comp, model == "+ region")$p.value))
cat(sprintf("After adjusting for both: slope = %.4f (p = %.4e)\n",
  filter(model_comp, model == "+ clinical + region")$estimate,
  filter(model_comp, model == "+ clinical + region")$p.value))
cat("\nClinical % trend:\n")
cat(sprintf("  slope = %.4f, p = %.4e\n", coef(m_clin)[2], summary(m_clin)$coefficients[2,4]))
cat("\nInterpretation:\n")
if (summary(m_clin)$coefficients[2,4] < 0.05) {
  cat("  Clinical proportion is changing significantly over time, which may contribute\n")
  cat("  to the ARG burden trend if clinical isolates carry more ARGs.\n")
} else {
  cat("  Clinical proportion is stable. The ARG decrease is not driven by\n")
  cat("  shifting clinical representation.\n")
}
sink()

# ---- 14. Plots ----
p1 <- ggplot(annual, aes(x = year, y = mean_arg)) +
  geom_line() + geom_point(aes(size = n)) +
  geom_ribbon(aes(ymin = mean_arg - se_arg, ymax = mean_arg + se_arg), alpha = 0.2) +
  theme_minimal() + labs(title = "ST10 Mean ARG Burden (2016-2025)", y = "Mean ARG count", x = "")
ggsave(file.path(OUT, "arg_trend_overall.pdf"), p1, width = 6, height = 4)

p2 <- ggplot(annual_clinical %>% filter(clinical != "Unknown"), aes(x = year, y = mean_arg, color = clinical)) +
  geom_line() + geom_point(aes(shape = clinical)) + theme_minimal() +
  labs(title = "ARG Trend by Clinical Status", y = "Mean ARG count", x = "", color = "", shape = "")
ggsave(file.path(OUT, "arg_trend_by_clinical.pdf"), p2, width = 7, height = 4)

p3 <- ggplot(annual_region %>% filter(region != "Other", n >= 10), aes(x = year, y = mean_arg, color = region)) +
  geom_line() + geom_point() + theme_minimal() +
  labs(title = "ARG Trend by Region", y = "Mean ARG count", x = "", color = "")
ggsave(file.path(OUT, "arg_trend_by_region.pdf"), p3, width = 8, height = 5)

p4 <- ggplot(clinical_trend, aes(x = year, y = pct_clinical)) +
  geom_line() + geom_point() + theme_minimal() +
  labs(title = "Clinical Proportion Over Time", y = "% Clinical", x = "")
ggsave(file.path(OUT, "clinical_proportion.pdf"), p4, width = 6, height = 4)

cat("Plots saved.\n")
cat("Done. Output in:", OUT, "\n")

})  # end block 12: 10f_st10_arg_decomposition.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 13/17: 10g_long_read_validation.R
# Long-read validation
# (verbatim body of original scripts/.../10g_long_read_validation.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
# Reviewer response: long-read validation of gene linkage
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(writexl); library(data.table)
})
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_long_read_validation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Check if any long-read assemblies exist ----
lr_dir <- file.path(config$OUTPUT_DIR, "ST69", "long_read")
if (dir.exists(lr_dir)) {
  fa_files <- list.files(lr_dir, pattern = "\\.fasta$|\\.fna$", full.names = TRUE)
  gff_files <- list.files(lr_dir, pattern = "\\.gff$|\\.gff3$", full.names = TRUE)

  if (length(fa_files) > 0) {
    cat("Long-read assemblies found:", length(fa_files), "\n")
    cat("GFF files:", length(gff_files), "\n")

    # Load existing VF results for these assemblies
    vf_dir <- file.path(config$OUTPUT_DIR, "ST69", "virulencefinder_validation")
    if (dir.exists(vf_dir)) {
      vf_files <- list.files(vf_dir, pattern = "\\.tsv$", full.names = TRUE)
      lr_names <- basename(fa_files) %>% str_remove("\\.fasta$|\\.fna$")

      cat("\nVirulence genes in long-read assemblies:\n")
      for (vf in vf_files) {
        tbl <- tryCatch(fread(vf, sep = "\t", header = TRUE, data.table = FALSE,
          nThread = 1), error = function(e) NULL)
        if (!is.null(tbl) && nrow(tbl) > 0) {
          cat(basename(vf), ":", nrow(tbl), "hits\n")
        }
      }
    }
    cat("\nRecommendation: download long-read assemblies from Lipworth et al. 2024\n")
    cat("and run ABRicate/VFDB to validate pap-kps-sat co-localization on single\n")
    cat("contigs. These assemblies are publicly available on ENA (PRJEB76009).\n")
  } else {
    cat("Long-read directory exists but no assemblies found.\n")
    cat("No long-read assemblies available for validation.\n")
  }
} else {
  cat("Long-read directory not found:", lr_dir, "\n")
  cat("No long-read assemblies available for validation.\n")
}

# ---- 2. Check project structure for assemblies ----
assembly_base <- file.path(config$BASE_DIR, "assembly")
if (dir.exists(assembly_base)) {
  st69_asm <- list.files(assembly_base, pattern = "ST69", full.names = TRUE)
  if (length(st69_asm) > 0) {
    cat("\nAssembly directory for ST69 exists:", st69_asm[1], "\n")
  }
}

# ---- 3. Write summary ----
sink(file.path(OUT, "long_read_validation_summary.txt"))
cat("Long-Read Validation of Virulence Locus Linkage\n")
cat("==============================================\n\n")
cat("Status: No long-read assemblies are currently available in the project.\n\n")
cat("Recommendation for reviewer response:\n")
cat("The RGP neighbourhood analysis using PPanGGOLiN already demonstrates that\n")
cat("pap, iuc, kps, and sat genes are physically linked within 5 composite RGPs\n")
cat("in ST69 genomes. This provides graph-based evidence for co-localization.\n\n")
cat("As suggested by the reviewer, long-read validation would be a useful addition.\n")
cat("Lipworth et al. (2024, Nat Commun) and Arredondo-Alonso et al. (2025) have\n")
cat("published long-read assemblies of ST69 E. coli. These should be downloaded\n")
cat("from ENA and analysed with ABRicate to confirm that the full virulence locus\n")
cat("is present on single contigs.\n")
sink()
cat("Done. Output in:", OUT, "\n")

})  # end block 13: 10g_long_read_validation.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 14/17: 14_fig06_gene_trajectories.R
# Gene-level early-vs-late prevalence, 3-panel (exploratory/QC -- distinct from, and not the source of, submitted Figure 6, which is 16_fig07_phylogeny_summary.R / Analysis/Figure6.R)
# (verbatim body of original scripts/.../14_fig06_gene_trajectories.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
#
# Gene-level early vs late prevalence analysis — 3-panel comparison
# Replicates the Figure 5/6 analysis from 3_generate_figures.R
# Panel A: ST69 VFDB | Panel B: ST10 VirulenceFinder | Panel C: ST10 ResFinder
#
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(ggplot2); library(patchwork); library(scales)
})
source("config.R")

OUT <- file.path(config$OUTPUT_DIR, "figures_gene_prevalence_3panel")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

EARLY_YEARS <- c(2016, 2017, 2018)
LATE_YEARS  <- c(2022, 2023, 2024, 2025)
MIN_DELTA_PP <- 5
MIN_PREVALENCE_LATE <- 0.05

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

safe_fisher_gene <- function(present_late, absent_late, present_early, absent_early) {
  mat <- matrix(c(present_late, absent_late, present_early, absent_early), nrow = 2, byrow = TRUE)
  out <- tryCatch(fisher.test(mat), error = function(e) NULL)
  if (is.null(out)) return(tibble(fisher_or = NA_real_, fisher_p = NA_real_))
  tibble(fisher_or = unname(out$estimate), fisher_p = out$p.value)
}

safe_logistic_gene <- function(d) {
  d <- d %>% filter(!is.na(year), !is.na(present)) %>% mutate(year_c = year - median(year, na.rm = TRUE))
  if (nrow(d) < 10 || n_distinct(d$year) < 3 || length(unique(d$present)) < 2)
    return(tibble(logistic_or = NA_real_, logistic_p = NA_real_))
  fit <- tryCatch(suppressWarnings(glm(present ~ year_c, family = binomial(), data = d)), error = function(e) NULL)
  if (is.null(fit)) return(tibble(logistic_or = NA_real_, logistic_p = NA_real_))
  sm <- summary(fit)$coefficients
  if (!"year_c" %in% rownames(sm)) return(tibble(logistic_or = NA_real_, logistic_p = NA_real_))
  tibble(logistic_or = exp(sm["year_c", "Estimate"]), logistic_p = sm["year_c", "Pr(>|z|)"])
}

# VF system annotation (from 3_generate_figures.R)
vf_groups <- list(
  Yersiniabactin = c("irp1", "irp2", "fyuA/psn", "ybtA", "ybtE", "ybtP", "ybtQ", "ybtS", "ybtT", "ybtU", "ybtX"),
  P_fimbriae = c("papA", "papB", "papC", "papD", "papE", "papF", "papG", "papH", "papI", "papJ", "papK", "papX"),
  Capsule_kps = c("kpsC", "kpsD", "kpsE", "kpsF", "kpsM", "kpsS", "kpsT", "kpsU", "neuA", "neuB", "neuC", "neuD", "neuE"),
  iro_salmochelin = c("iroB", "iroC", "iroD", "iroE", "iroN"),
  Type1_fimbriae = c("fimA", "fimB", "fimC", "fimD", "fimE", "fimF", "fimG", "fimH", "fimI"),
  T6SS = c("tssA", "tssB", "tssC", "tssF", "tssG", "tssJ", "tssK", "tssL", "tssM", "hcp/tssD", "hcp1/tssD1", "vgrG/tssI"),
  Toxins = c("sat", "senB", "hlyA", "hlyB", "hlyC", "hlyD", "cnf1", "cdtA", "cdtB", "cdtC", "east1", "espP", "pic", "pet")
)

vf_anno <- purrr::imap_dfr(vf_groups, ~ tibble(gene = .x, system = .y)) %>% distinct(gene, .keep_all = TRUE)

# ARG class annotation (broad categories)
arg_classes <- c("Aminoglycosides" = "aac|aad|aph|ant|arm|rmt", "Beta-lactams" = "bla|oxa|nps", "Macrolides" = "erm|mef|mph|msr|lnu|lin",
                 "Tetracyclines" = "tet", "Sulfonamides" = "sul|dfr", "Phenicols" = "cat|cml",
                 "Quinolones" = "qnr|oqx|qep", "Polymyxins" = "mcr|arn", "Fosfomycin" = "fos")

classify_system <- function(g, type = "vf") {
  if (type == "vf") {
    m <- vf_anno$system[match(g, vf_anno$gene)]
    ifelse(is.na(m), "Other", m)
  } else {
    m <- names(arg_classes)[sapply(arg_classes, function(p) grepl(p, g, ignore.case = TRUE))]
    ifelse(length(m) > 0, m[1], "Other")
  }
}

# -------- Load data --------
cat("Loading data...\n")
meta69 <- read_xlsx(config$st_metadata("ST69")) %>%
  rename(genome_id = Name) %>% mutate(genome_id = make_genome_id(genome_id), year = as.integer(.data[["Collection Year"]])) %>%
  filter(!is.na(year)) %>% select(genome_id, year)
meta10 <- read_xlsx(config$st_metadata("ST10")) %>%
  rename(genome_id = Name) %>% mutate(genome_id = make_genome_id(genome_id), year = as.integer(.data[["Collection Year"]])) %>%
  filter(!is.na(year)) %>% select(genome_id, year)

# ST69 VFDB raw
vf69_raw <- read.delim(config$st_vfdb_summary("ST69"), header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
colnames(vf69_raw)[1:2] <- c("genome_id", "num_found")
vf69_raw$genome_id <- make_genome_id(str_remove(basename(vf69_raw$genome_id), "_vfdb\\.tsv$"))
vf69_genes <- setdiff(colnames(vf69_raw), c("genome_id", "num_found"))
vf69 <- vf69_raw %>% mutate(across(all_of(vf69_genes), ~ as.integer(!is.na(.x) & .x != ".")))
vf69_long <- vf69 %>% select(genome_id, all_of(vf69_genes)) %>%
  pivot_longer(-genome_id, names_to = "gene", values_to = "present") %>%
  inner_join(meta69, by = "genome_id") %>%
  mutate(system = classify_system(gene, "vf"), period = case_when(year %in% EARLY_YEARS ~ "early", year %in% LATE_YEARS ~ "late", TRUE ~ NA_character_))
cat("  ST69 VFDB:", n_distinct(vf69_long$gene), "genes\n")

# ST10 VF binary
vf_binary <- read.delim(config$VF_BINARY, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
vf10_long <- vf_binary %>% filter(st == "ST10") %>%
  mutate(genome_id = make_genome_id(genome)) %>%
  select(-st, -genome) %>% pivot_longer(-genome_id, names_to = "gene", values_to = "present") %>%
  mutate(present = as.integer(!is.na(present) & present != "" & present != "0")) %>%
  inner_join(meta10, by = "genome_id") %>%
  mutate(system = classify_system(gene, "vf"), period = case_when(year %in% EARLY_YEARS ~ "early", year %in% LATE_YEARS ~ "late", TRUE ~ NA_character_))
cat("  ST10 VF:", n_distinct(vf10_long$gene), "genes\n")

# ST10 ResFinder binary
resf_binary <- read.delim(file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
                           header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
resf10_long <- resf_binary %>% filter(st == "ST10") %>%
  mutate(genome_id = make_genome_id(genome)) %>%
  select(-st, -genome) %>% pivot_longer(-genome_id, names_to = "gene", values_to = "present") %>%
  mutate(present = as.integer(!is.na(present) & present != "" & present != "0")) %>%
  inner_join(meta10, by = "genome_id") %>%
  mutate(system = classify_system(gene, "arg"), period = case_when(year %in% EARLY_YEARS ~ "early", year %in% LATE_YEARS ~ "late", TRUE ~ NA_character_))
cat("  ST10 ResFinder:", n_distinct(resf10_long$gene), "genes\n")

# -------- Analysis function --------
analyze_genes <- function(long_df, case_name) {
  period_df <- long_df %>% filter(!is.na(period))

  per_gene <- period_df %>%
    group_by(gene, system) %>%
    summarise(
      n_early = sum(period == "early"), n_late = sum(period == "late"),
      present_early = sum(present[period == "early"]),
      present_late  = sum(present[period == "late"]),
      absent_early = n_early - present_early,
      absent_late  = n_late - present_late,
      prevalence_early = present_early / n_early,
      prevalence_late  = present_late / n_late,
      delta_pp = (prevalence_late - prevalence_early) * 100,
      .groups = "drop"
    ) %>% rowwise() %>%
    mutate(fr = list(safe_fisher_gene(present_late, absent_late, present_early, absent_early))) %>%
    unnest(fr) %>% ungroup() %>%
    mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

  trend_df <- long_df %>% filter(!is.na(year)) %>%
    group_by(gene, system) %>%
    group_modify(~ safe_logistic_gene(.x)) %>% ungroup() %>%
    mutate(logistic_p_adj = p.adjust(logistic_p, method = "BH"))

  results <- per_gene %>% left_join(trend_df, by = c("gene", "system")) %>%
    mutate(
      direction = case_when(
        delta_pp >= MIN_DELTA_PP & prevalence_late >= MIN_PREVALENCE_LATE &
          (fisher_p_adj < 0.05 | logistic_p_adj < 0.05) ~ "Increased",
        delta_pp <= -MIN_DELTA_PP &
          (fisher_p_adj < 0.05 | logistic_p_adj < 0.05) ~ "Decreased",
        TRUE ~ "No strong change"
      ),
      case = case_name
    ) %>% arrange(desc(delta_pp))

  write_csv(results, file.path(OUT, paste0(case_name, "_gene_level_results.csv")))
  cat("  ", case_name, ":", sum(results$direction == "Increased"), "increased,",
      sum(results$direction == "Decreased"), "decreased genes\n")
  results
}

cat("\nAnalyzing...\n")
r69 <- analyze_genes(vf69_long, "ST69_VFDB")
r10vf <- analyze_genes(vf10_long, "ST10_VF")
r10res <- analyze_genes(resf10_long, "ST10_ResFinder")

# -------- Top increasing genes dumbbell plot (3 panels) --------
plot_dumbbell <- function(results, case_name, color) {
  inc <- results %>% filter(direction == "Increased")
  n_show <- min(20, nrow(inc))
  if (n_show == 0) {
    df <- data.frame(x = 0, y = "", label = "No increasing genes")
    return(ggplot(df, aes(x, y)) + geom_text(aes(label = label), size = 5) +
             labs(title = case_name) + theme_classic(base_size = 11) + xlim(-1, 1))
  }
  top <- inc %>% slice_max(delta_pp, n = n_show, with_ties = FALSE) %>%
    mutate(gene_label = paste0(gene, " [", system, "]"))

  ggplot(top, aes(x = delta_pp, y = reorder(gene_label, delta_pp))) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_segment(aes(x = 0, xend = delta_pp, y = reorder(gene_label, delta_pp), yend = reorder(gene_label, delta_pp)),
                 linewidth = 0.7, color = color) +
    geom_point(size = 2.5, color = color) +
    labs(title = case_name, x = "Late minus early prevalence (pp)", y = NULL) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11))
}

cat("\nPlotting...\n")
pA <- plot_dumbbell(r69, "ST69 VFDB", "#2166AC")
pB <- plot_dumbbell(r10vf, "ST10 VirulenceFinder", "#D73027")
pC <- plot_dumbbell(r10res, "ST10 ResFinder", "#1B7837")

fig <- (pA | pB | pC) + plot_layout(ncol = 3) +
  plot_annotation(title = "Top increasing genes: early (2016-2018) vs late (2022-2025)",
                  subtitle = sprintf("Min delta >= %d pp, late prevalence >= %d%%, Fisher or logistic BH-adjusted p < 0.05",
                                     MIN_DELTA_PP, as.integer(MIN_PREVALENCE_LATE * 100)),
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
                                plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5)))

ggsave(file.path(OUT, "Fig_gene_prevalence_3panel.png"), fig, width = 16, height = 8, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig_gene_prevalence_3panel.pdf"), fig, width = 16, height = 8, bg = "white")
cat("\nSaved: Fig_gene_prevalence_3panel\n")

})  # end block 14: 14_fig06_gene_trajectories.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 15/17: 15_fig07_tree_parsimony.R
# Core-tree parsimony figure (preliminary/QC -- single-ring only; NOT the source of submitted Figure 3, which was made in iTOL -- see Figure3.R)
# (verbatim body of original scripts/.../15_fig07_tree_parsimony.R, wrapped in local())
################################################################################
tryCatch({
local({
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

})  # end block 15: 15_fig07_tree_parsimony.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 16/17: 22_fig_supplementary_combined.R
# FINAL: self-contained combined 3-page Supplementary Figures PDF (re-derives its own data; run last, after everything above so every individual xlsx/csv/png output above is also preserved)
# (verbatim body of original scripts/.../22_fig_supplementary_combined.R, wrapped in local())
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
#
# 22_fig_supplementary_combined.R
#
# Combined supplementary figure bringing together ALL reviewer-response analyses:
#   Page 1 (ST69 capsule): A-F  (capsule classification, K-type, clinical enrichment)
#   Page 2 (Genomic context): G-L (kps validation, RGP co-localization, MGE profiles)
#   Page 3 (Remaining):      M-R (allelic conversion, plasmid context, ST10 ARG decomposition)
#
# Output: Supplementary_Figures_Combined.pdf  (in config$OUTPUT_DIR/combined_figures)

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(patchwork)
  library(scales); library(viridis); library(broom)
  if (requireNamespace("readxl", quietly=TRUE)) library(readxl)
  if (requireNamespace("writexl", quietly=TRUE)) library(writexl)
  if (requireNamespace("ape", quietly=TRUE)) library(ape)
  if (requireNamespace("phangorn", quietly=TRUE)) library(phangorn)
})

# Config is loaded from the repo-root config.R (env-driven paths).
if (file.exists("config.R")) {
  source("config.R")
} else {
  source("../config.R")
}

if (is.null(config$TARGET_ST)) config$TARGET_ST <- "ST69"
OUT <- file.path(config$OUTPUT_DIR, "combined_figures")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

EARLY <- c(2016, 2017, 2018)
LATE  <- c(2022, 2023, 2024, 2025)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

to_bin_vfdb <- function(x) as.integer(sapply(strsplit(as.character(x), ";"),
  function(v) any(suppressWarnings(as.numeric(v) >= 90), na.rm = TRUE)))
to_bin <- function(x) as.integer(!is.na(x) & x != ".")

theme_pub <- function(base_size=9) {
  theme_classic(base_size=base_size) +
    theme(plot.title=element_text(face="bold", size=rel(1.05), hjust=0),
          axis.title=element_text(face="bold", size=rel(0.9)),
          axis.text=element_text(color="grey20", size=rel(0.78)),
          axis.line=element_line(color="grey30", linewidth=0.3),
          axis.ticks=element_line(color="grey30", linewidth=0.25),
          panel.grid.major.y=element_line(color="grey93", linewidth=0.2),
          strip.background=element_rect(fill="grey95", color=NA),
          strip.text=element_text(face="bold", size=rel(0.85)),
          legend.title=element_text(face="bold", size=rel(0.82)),
          legend.text=element_text(size=rel(0.75)),
          legend.position="bottom",
          plot.background=element_rect(fill="white", color=NA),
          plot.tag=element_text(face="bold", size=rel(1.2)))
}
theme_set(theme_pub())

cluster_colors <- c("Cluster_1"="#E41A1C","Cluster_2"="#377EB8",
                    "Cluster_3"="#4DAF4A","Cluster_4"="#984EA3")
caps_colors <- c("G2"="#d73027","G3"="#4575b4","No capsule"="#808080","Unclassified"="#fddbc7")

pdf_file <- file.path(OUT, "Supplementary_Figures_Combined.pdf")
pdf(pdf_file, width=16, height=22, bg="white")

cat("============================================\n")
cat(" Combined Supplementary Figures\n")
cat(" Output:", pdf_file, "\n")
cat("============================================\n\n")

# ============================================================
# LOAD SHARED DATA: ST69 master + VFDB
# ============================================================
cat("--- Loading shared ST69 data ---\n")

master_file <- file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv")
st69_master <- fread(master_file, sep=",", header=TRUE, data.table=FALSE, na.strings="") %>%
  mutate(genome_id=as.character(genome_id), year=as.integer(year),
         shell_cluster=as.character(shell_cluster))

vfdb_raw <- fread(config$st_vfdb_summary("ST69"), sep="\t", header=TRUE,
  data.table=FALSE, colClasses="character", na.strings="") %>%
  filter(str_detect(.data[["#FILE"]], "^ST69/")) %>%
  mutate(genome_id=make_genome_id(str_remove(.data[["#FILE"]], "^ST69/")),
         genome_id=str_remove(genome_id, "_vfdb\\.tsv$"))

vf_genes <- setdiff(colnames(vfdb_raw), c("#FILE","NUM_FOUND","genome_id"))
vf_bin <- vfdb_raw %>% select(genome_id, all_of(vf_genes)) %>%
  mutate(across(all_of(vf_genes), to_bin_vfdb))

meta_cols <- c("genome_id","shell_cluster","year","clinical_binary","total_vf",
  intersect(colnames(st69_master), c("#FILE","FILE","GFF","species","st","country",
  "continent","source","raw_cluster")))
st69_meta <- st69_master %>% select(all_of(intersect(meta_cols, colnames(st69_master))))
df <- st69_meta %>% inner_join(vf_bin, by="genome_id") %>%
  mutate(year=as.integer(year), clin=as.integer(clinical_binary)) %>%
  filter(!is.na(clin), !is.na(year), year>=2016, year<=2025)

# ============================================================
# PAGE 1: ST69 CAPSULE CLASSIFICATION (panels A-F)
# ============================================================
cat("\n=== PAGE 1: Capsule Classification ===\n")

# --- Panel A: Capsule distribution pie chart ---
kps_cols <- intersect(grep("^(kps|neu)", vf_genes, value=TRUE), colnames(df))
G2_markers <- intersect(c("kpsC","kpsS"), kps_cols)
G3_markers <- intersect(c("kpsE","kpsM","kpsT"), kps_cols)

df$capsule_type <- apply(df[, kps_cols, drop=FALSE], 1, function(r) {
  has_g2 <- any(r[names(r) %in% G2_markers] >= 1)
  has_g3 <- any(r[names(r) %in% G3_markers] >= 1)
  any_kps <- any(r >= 1)
  if (has_g2) return("G2")
  if (has_g3) return("G3")
  if (any_kps) return("Unclassified")
  return("No capsule")
})
caps_dist <- df %>% count(capsule_type) %>% mutate(pct=n/sum(n)*100)

pA <- ggplot(caps_dist, aes(x="", y=pct, fill=capsule_type)) +
  geom_bar(stat="identity", width=0.5) + coord_polar("y") +
  scale_fill_manual(values=caps_colors) +
  labs(title="A. ST69 capsule types", fill="", y="", x="") +
  theme_minimal(base_size=9) +
  theme(axis.text=element_blank(), axis.title=element_blank(),
        panel.grid=element_blank(), plot.tag=element_text(face="bold", size=rel(1.2)))

# --- Panel B: Capsule temporal (early vs late) ---
caps_temporal <- df %>%
  mutate(period=case_when(year%in%EARLY~"early", year%in%LATE~"late", TRUE~"mid")) %>%
  filter(period%in%c("early","late")) %>%
  group_by(period) %>% mutate(period_total=n()) %>% ungroup() %>%
  group_by(capsule_type, period) %>%
  summarise(n=n(), period_total=first(period_total), .groups="drop") %>%
  mutate(prev=n/period_total*100) %>%
  select(-period_total) %>%
  pivot_wider(names_from=period, values_from=c(n,prev), values_fill=0) %>%
  mutate(delta_pp=prev_late-prev_early)

pB <- ggplot(caps_temporal, aes(x=reorder(capsule_type, delta_pp), y=delta_pp, fill=capsule_type)) +
  geom_col(alpha=0.85) + coord_flip() +
  scale_fill_manual(values=caps_colors, guide="none") +
  labs(title="B. Capsule type temporal change", x="", y="Delta prevalence (pp)") +
  geom_hline(yintercept=0, linetype="dashed", color="grey50")

# --- Panel C: K-type distribution ---
# Uses the all-ST binary matrix, filtered to the target ST.
vf_freq_file <- file.path(config$INPUT_DIR,
  "virulencefinder_summary", "virulencefinder_binary_matrix.tsv")
ktype_plot <- NULL
if (file.exists(vf_freq_file)) {
  vf <- fread(vf_freq_file, sep="\t", header=TRUE, data.table=FALSE, na.strings="") %>%
    filter(.data[["st"]] == config$TARGET_ST)
  if (!"genome" %in% colnames(vf)) vf <- vf %>% rename(genome = 2)
  k_cols <- grep("kpsM[I]*[_]", colnames(vf), value=TRUE)
  k_cols <- setdiff(k_cols, c("kpsMII","kpsMIII"))
  k_type_map <- c(
    kpsMII_K1="K1", kpsMII_K5="K5", kpsMII_K52="K52", kpsMII_K4="K4",
    kpsMII_K23="K23", kpsM_K15="K15", kpsM_K11="K11", kpsM_K19="K19",
    kpsM_K19K23="K19/K23", kpsMIII_K96="K96", kpsMIII_K98="K98", kpsMIII_K10="K10")
  present_k <- intersect(k_cols, names(k_type_map))
  k_type_map <- k_type_map[present_k]
  is_present <- function(x) { x <- as.character(x); x[is.na(x)] <- "."; x=="1" }
  ktype_df <- vf %>% select(genome, all_of(present_k))
  colnames(ktype_df)[1] <- "Genome"
  ktype_df$k_type <- apply(ktype_df[, present_k, drop=FALSE], 1, function(r) {
    hits <- names(r)[is_present(r)]
    if (length(hits)==0) return(NA_character_)
    unname(k_type_map[hits[1]])
  })
  has_g2 <- is_present(vf[["kpsMII"]])
  has_g3 <- is_present(vf[["kpsMIII"]])
  ktype_df$k_type <- ifelse(is.na(ktype_df$k_type) & has_g2, "G2-unknown", ktype_df$k_type)
  ktype_df$k_type <- ifelse(is.na(ktype_df$k_type) & has_g3, "G3-unknown", ktype_df$k_type)
  ktype_df$k_type <- ifelse(is.na(ktype_df$k_type), "No kpsM allele", ktype_df$k_type)
  k_tab <- as.data.frame(table(ktype_df$k_type, useNA="ifany"))
  colnames(k_tab) <- c("k_type","n")
  k_tab <- k_tab %>% mutate(pct=n/sum(n)*100)
  ktype_colors <- c("K1"="#e41a1c","K5"="#377eb8","K52"="#4daf4a","K4"="#984da3",
    "K23"="#ff7f00","K15"="#a65628","K11"="#f781bf","K19"="#999999",
    "K19/K23"="#66c2a5","K96"="#fc8d62","K98"="#8da0cb","K10"="#e78ac3",
    "G2-unknown"="#ffd92f","G3-unknown"="#e5c494","No kpsM allele"="grey80")
  present_colors <- intersect(names(ktype_colors), k_tab$k_type)
  pC <- ggplot(k_tab, aes(x="", y=pct, fill=k_type)) +
    geom_bar(stat="identity", width=0.5) + coord_polar("y") +
    scale_fill_manual(values=ktype_colors[present_colors]) +
    labs(title="C. K-type distribution", fill="", y="", x="") +
    theme_minimal(base_size=9) +
    theme(axis.text=element_blank(), axis.title=element_blank(),
          panel.grid=element_blank())
} else {
  pC <- ggplot() + annotate("text", x=0.5, y=0.5, label="VF data not available", size=4) +
    theme_void() + labs(title="C. K-type distribution (data unavailable)")
}

# --- Panel D: Temporal increase in C3 (G2/G3 colored) ---
key_genes <- intersect(c("kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT","kpsU",
  "papX","papF","papB","papG","sat","hlyA"), vf_genes)
cd3 <- df %>% filter(shell_cluster=="Cluster_3")
ea <- cd3 %>% filter(year%in%EARLY)
la <- cd3 %>% filter(year%in%LATE)
n_early <- nrow(ea); n_late <- nrow(la)

temporal_inc <- map_dfr(key_genes, ~{
  g <- .x
  prev_e <- mean(ea[[g]], na.rm=TRUE)*100
  prev_l <- mean(la[[g]], na.rm=TRUE)*100
  dp <- prev_l - prev_e
  mat <- matrix(c(sum(ea[[g]]), n_early-sum(ea[[g]]), sum(la[[g]]), n_late-sum(la[[g]])), nrow=2)
  ft <- tryCatch(fisher.test(mat), error=function(e) NULL)
  p <- if(!is.null(ft)) ft$p.value else NA
  tibble(gene=g, early_prev=prev_e, late_prev=prev_l, delta_pp=dp, p=p)
}) %>% mutate(p_adj=p.adjust(pmax(p,0,na.rm=TRUE), method="BH"),
  increasing=delta_pp>=5 & late_prev>=5 & p_adj<0.05)

pD <- temporal_inc %>% filter(increasing | gene%in%key_genes) %>%
  mutate(is_g2=gene%in%G2_markers) %>%
  group_by(gene, is_g2) %>% summarise(delta_pp=first(delta_pp), .groups="drop") %>%
  mutate(gene=factor(gene, levels=gene[order(delta_pp)])) %>%
  ggplot(aes(x=delta_pp, y=gene, fill=is_g2)) +
  geom_col(alpha=0.85) +
  scale_fill_manual(values=c("TRUE"="#d73027","FALSE"="#4575b4"),
    labels=c("G3/other","G2 marker"), name="") +
  labs(title="D. Cluster_3 increasing genes (G2/G3 colored)", x="Delta prevalence (pp)", y=NULL) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50")

# --- Panel E: Clinical enrichment (G2/G3 colored) ---
clin_enrich <- map_dfr(key_genes, ~{
  g <- .x
  d <- df %>% filter(!is.na(.data[[g]]))
  if (sum(d[[g]], na.rm=TRUE)<5) return(NULL)
  clin_prev <- mean(d[[g]][d$clin==1])*100
  non_prev <- mean(d[[g]][d$clin==0])*100
  ft <- tryCatch(fisher.test(table(d$clin, d[[g]])), error=function(e) NULL)
  or <- if(!is.null(ft)) unname(ft$estimate) else NA
  p <- if(!is.null(ft)) ft$p.value else NA
  tibble(gene=g, clin_prev=clin_prev, non_prev=non_prev,
    delta_pp=clin_prev-non_prev, or=or, p=p)
}) %>% mutate(p_adj=p.adjust(pmax(p,0,na.rm=TRUE), method="BH"),
  enriched=!is.na(p_adj) & p_adj<0.05 & or>1)

pE <- clin_enrich %>% filter(enriched | gene%in%key_genes) %>%
  mutate(is_g2=gene%in%G2_markers) %>%
  group_by(gene, is_g2) %>% summarise(delta_pp=first(delta_pp), .groups="drop") %>%
  mutate(gene=factor(gene, levels=gene[order(delta_pp)])) %>%
  ggplot(aes(x=delta_pp, y=gene, fill=is_g2)) +
  geom_col(alpha=0.85) +
  scale_fill_manual(values=c("TRUE"="#d73027","FALSE"="#4575b4"),
    labels=c("G3/other","G2 marker"), name="") +
  labs(title="E. Clinical enrichment (G2/G3 colored)", x="Delta prevalence (pp)", y=NULL) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50")

# --- Panel F: OR vs delta prevalence scatter ---
pF <- clin_enrich %>% filter(enriched | gene%in%key_genes) %>%
  mutate(g2=gene%in%G2_markers) %>%
  ggplot(aes(x=delta_pp, y=or, color=g2, size=-log10(pmax(p_adj,1e-300)))) +
  geom_point(alpha=0.8) + scale_y_log10() +
  scale_color_manual(values=c("TRUE"="#d73027","FALSE"="#4575b4"),
    labels=c("G3/other","G2 marker"), name="") +
  geom_hline(yintercept=1, linetype="dashed") +
  labs(title="F. Effect size vs association", x="Delta prevalence (pp)",
       y="Odds ratio (log scale)", size=expression(-log[10](p)))

# --- Compose Page 1 ---
page1 <- (pA + pB) / (pC + pD) / (pE + pF) +
  plot_layout(heights=c(1,1.2,1.2)) +
  plot_annotation(tag_levels="A",
    title="Supplementary Figure S1. ST69 Capsule Classification and Clinical Enrichment") &
  theme(plot.tag=element_text(face="bold"))
print(page1)

# ============================================================
# PAGE 2: GENOMIC CONTEXT (panels G-L)
# ============================================================
cat("\n=== PAGE 2: Genomic Context ===\n")

# --- Panel G: kps prevalence by cluster (VFDB vs PPanGGOLiN) ---
kps_genes_vfdb <- grep("^(kps|neu)", vf_genes, value=TRUE)
kps_cols_present <- intersect(kps_genes_vfdb, colnames(vf_bin))
vfdb_kps <- vf_bin %>% select(genome_id, all_of(kps_cols_present)) %>%
  mutate(kps_any_vfdb=as.integer(rowSums(select(., -genome_id), na.rm=TRUE) > 0))
assignments <- read_csv(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "03_shell_gene_cluster_assignments_k4.csv"), show_col_types=FALSE) %>%
  mutate(genome_id=as.character(genome_id))
vfdb_kps_joined <- vfdb_kps %>% left_join(assignments, by="genome_id")

# PP kps
pp_kps_data <- NULL
rtab_file <- file.path(config$PANGENOME_DIR, "gene_presence_absence.Rtab")
vf_families_file <- file.path(config$PANGENOME_DIR, "vf_module_rgp_results",
  "tables", "gene_family_annotation_VF_only.csv")
if (file.exists(rtab_file) && file.exists(vf_families_file)) {
  rtab <- fread(rtab_file, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
  vf_families <- read_csv(vf_families_file, show_col_types=FALSE)
  kps_fams <- vf_families %>% filter(vf_system=="kps") %>% distinct(family_id)
  rtab_fams <- rtab[[1]]; rtab_mat <- as.matrix(rtab[,-1])
  rownames(rtab_mat) <- rtab_fams
  genome_ids <- colnames(rtab)[-1]
  kps_fam_ids <- intersect(kps_fams$family_id, rtab_fams)
  if (length(kps_fam_ids) > 0) {
    kps_idx <- which(rtab_fams %in% kps_fam_ids)
    pp_kps_present <- as.integer(colSums(rtab_mat[kps_idx,,drop=FALSE] > 0) > 0)
    pp_kps_data <- tibble(genome_id=genome_ids, pp_kps=pp_kps_present)
  }
}

kps_by_cluster <- vfdb_kps_joined %>% filter(!is.na(shell_cluster)) %>%
  group_by(shell_cluster) %>%
  summarise(n=n(), vfdb_kps_pct=mean(kps_any_vfdb, na.rm=TRUE)*100, .groups="drop")

if (!is.null(pp_kps_data)) {
  vfdb_kps_with_pp <- vfdb_kps_joined %>% inner_join(pp_kps_data, by="genome_id")
  pp_by_cluster <- vfdb_kps_with_pp %>% filter(!is.na(shell_cluster)) %>%
    group_by(shell_cluster) %>%
    summarise(pp_kps_pct=mean(pp_kps, na.rm=TRUE)*100, .groups="drop")
  kps_by_cluster <- kps_by_cluster %>% left_join(pp_by_cluster, by="shell_cluster")
}

kps_plot_data <- kps_by_cluster %>%
  pivot_longer(-c(shell_cluster,n), names_to="method", values_to="pct") %>%
  mutate(method=ifelse(method=="vfdb_kps_pct","VFDB","PPanGGOLiN"))

pG <- ggplot(kps_plot_data, aes(x=shell_cluster, y=pct, fill=method)) +
  geom_col(position="dodge", alpha=0.85) +
  scale_fill_manual(values=c("VFDB"="#E41A1C","PPanGGOLiN"="#377EB8")) +
  labs(title="G. kps prevalence by cluster", x="Cluster", y="% genomes", fill="Method") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# --- Panel H: Gene pair co-occurrence in RGPs ---
rgp_vf_file <- file.path(config$PANGENOME_DIR, "vf_module_rgp_results",
  "tables", "VF_families_in_RGPs_with_modules.csv")
pH <- NULL
if (file.exists(rgp_vf_file)) {
  rgp_vf <- fread(rgp_vf_file, sep=",", header=TRUE, data.table=FALSE)
  increasing_genes <- c("kpsD","kpsF","kpsM","kpsT","papB","papF","papX","sat","hlyA","hlyD","cnf1")
  gene_pairs <- combn(increasing_genes, 2, simplify=FALSE)
  cooc <- bind_rows(lapply(gene_pairs, function(pair) {
    hits <- rgp_vf %>% filter(gene_name%in%pair) %>%
      group_by(rgp_id) %>% summarise(n_genes=n_distinct(gene_name), .groups="drop") %>%
      filter(n_genes==2)
    tibble(gene1=pair[1], gene2=pair[2], n_rgps=nrow(hits))
  })) %>% filter(n_rgps > 0)

  if (nrow(cooc) > 0) {
    pH <- ggplot(cooc, aes(x=gene1, y=gene2, fill=n_rgps)) +
      geom_tile(color="white", linewidth=0.3) +
      geom_text(aes(label=n_rgps), size=2.5) +
      scale_fill_viridis_c(trans="log1p", name="# RGPs") +
      labs(title="H. Gene pair co-occurrence in RGPs", x="", y="") +
      theme_minimal(base_size=9) +
      theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, face="bold", size=rel(0.75)),
        axis.text.y=element_text(face="bold", size=rel(0.75)))
  }
}
if (is.null(pH)) {
  pH <- ggplot() + annotate("text", x=0.5, y=0.5, label="RGP data not available", size=4) +
    theme_void() + labs(title="H. Gene pair co-occurrence (data unavailable)")
}

# --- Panel I: RGP system profiles ---
pI <- NULL
if (file.exists(rgp_vf_file)) {
  interest <- c("kpsD","kpsF","kpsM","kpsT","kpsC","kpsE","kpsS","kpsU",
    "papA","papB","papC","papD","papE","papF","papG","papH","papI","papJ","papK","papX",
    "sat","iucA","iucB","iucC","iucD","iutA")
  gene_sys <- setNames(c(rep("kps",8),rep("pap",12),"sat",rep("iuc",5)), interest)
  rgp_interest <- rgp_vf %>% filter(gene_name%in%interest) %>%
    mutate(system=gene_sys[gene_name]) %>%
    group_by(rgp_id) %>%
    summarise(systems_present=paste(sort(unique(system)), collapse=","),
              n_systems=n_distinct(system), .groups="drop")
  profile_count <- rgp_interest %>% count(systems_present, sort=TRUE)
  if (nrow(profile_count) > 0) {
    pI <- profile_count %>% slice_max(n, n=10) %>%
      ggplot(aes(x=reorder(systems_present, n), y=n)) +
      geom_col(fill="steelblue", alpha=0.8) + coord_flip() +
      labs(title="I. RGP profiles carrying interest genes", x="Systems in RGP", y="# RGPs")
  }
}
if (is.null(pI)) {
  pI <- ggplot() + annotate("text", x=0.5, y=0.5, label="RGP data not available", size=4) +
    theme_void() + labs(title="I. RGP profiles (data unavailable)")
}

# --- Panel J: System co-occurrence in C3 ---
pJ <- NULL
if (file.exists(rgp_vf_file) && exists("gene_sys")) {
  system_genes <- list(
    kps=c("kpsD","kpsF","kpsM","kpsT","kpsC","kpsE","kpsS","kpsU"),
    pap=c("papA","papB","papC","papD","papE","papF","papG","papH","papI","papJ","papK","papX"),
    sat="sat", iuc=c("iucA","iucB","iucC","iucD","iutA"))
  c3_genomes <- df %>% filter(shell_cluster=="Cluster_3") %>% pull(genome_id)
  sys_matrix <- map_dfr(names(system_genes), ~{
    gs <- intersect(system_genes[[.x]], colnames(vf_bin))
    if (length(gs)==0) return(NULL)
    gene_cols <- intersect(paste0(gs, "_v"), colnames(vf_bin))
    if (length(gene_cols)==0) return(NULL)
    vals <- rowSums(vf_bin[match(c3_genomes, vf_bin$genome_id), gene_cols, drop=FALSE], na.rm=TRUE)
    has_sys <- vals >= ceiling(length(gs)*0.5)
    tibble(genome=c3_genomes, system=.x, present=has_sys)
  })
  if (nrow(sys_matrix) > 0) {
    pair_cooc <- map_dfr(combn(names(system_genes), 2, simplify=FALSE), ~{
      s1 <- .x[1]; s2 <- .x[2]
      d <- sys_matrix %>% filter(system%in%c(s1,s2)) %>%
        distinct(genome, system, present) %>%
        pivot_wider(names_from=system, values_from=present)
      n_both <- sum(d[[s1]] & d[[s2]], na.rm=TRUE)
      n_total <- nrow(d)
      tibble(sys1=s1, sys2=s2, pct=round(n_both/n_total*100,1))
    })
    if (nrow(pair_cooc) > 0) {
      pJ <- ggplot(pair_cooc, aes(x=sys1, y=sys2, fill=pct)) +
        geom_tile(color="white", linewidth=0.3) +
        geom_text(aes(label=sprintf("%.0f%%", pct)), size=3) +
        scale_fill_viridis_c(name="% co-occurring") +
        labs(title="J. System co-occurrence in Cluster_3", x="", y="") +
        theme_minimal(base_size=9) +
        theme(panel.grid=element_blank())
    }
  }
}
if (is.null(pJ)) {
  pJ <- ggplot() + annotate("text", x=0.5, y=0.5, label="Data not available", size=4) +
    theme_void() + labs(title="J. System co-occurrence (data unavailable)")
}

# --- Panel K: Kps operon progression in C3 ---
kps_all <- c("kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT","kpsU")
kps_vars <- paste0(kps_all, "_v")
c3_kps <- df %>% filter(shell_cluster=="Cluster_3") %>%
  mutate(period=ifelse(year%in%EARLY,"early",ifelse(year%in%LATE,"late",NA))) %>%
  filter(!is.na(period))
if (all(kps_vars %in% colnames(c3_kps))) {
  c3_kps <- c3_kps %>%
    mutate(n_kps=rowSums(across(all_of(kps_vars))),
           kps_complete=n_kps>=8, kps_partial=n_kps>=1 & n_kps<8, kps_none=n_kps==0)
  c3_kps_prog <- c3_kps %>% group_by(period) %>%
    summarise(n=n(), complete=sum(kps_complete), partial=sum(kps_partial),
              none=sum(kps_none), .groups="drop") %>%
    mutate(across(c(complete,partial,none), ~round(.x/n*100,1), .names="{.col}_pct"))

  pK <- c3_kps_prog %>%
    select(period, complete=complete_pct, partial=partial_pct, none=none_pct) %>%
    pivot_longer(-period, names_to="kps_status", values_to="pct") %>%
    mutate(kps_status=factor(kps_status, levels=c("none","partial","complete"))) %>%
    ggplot(aes(x=period, y=pct, fill=kps_status)) +
    geom_col(position="fill", width=0.6) +
    scale_fill_manual(values=c("none"="grey80","partial"="goldenrod","complete"="steelblue"),
      name="Kps status") +
    labs(title="K. Kps operon progression in C3", y="Proportion", x="") +
    scale_y_continuous(labels=percent_format())
} else {
  pK <- ggplot() + annotate("text", x=0.5, y=0.5, label="VF data incomplete", size=4) +
    theme_void() + labs(title="K. Kps operon progression (data unavailable)")
}

# --- Panel L: kps validation (VFDB vs PPanGGOLiN) discrepancy ---
pL <- NULL
if (!is.null(pp_kps_data)) {
  compare_data <- vfdb_kps_joined %>%
    inner_join(pp_kps_data, by="genome_id") %>%
    mutate(
      vfdb_pos=kps_any_vfdb==1, pp_pos=pp_kps==1,
      category=case_when(
        vfdb_pos & pp_pos ~ "Both",
        vfdb_pos & !pp_pos ~ "VFDB only",
        !vfdb_pos & pp_pos ~ "PPanGGOLiN only",
        TRUE ~ "Neither"))
  cat_counts <- compare_data %>% count(category) %>% mutate(pct=round(n/sum(n)*100,1))
  pL <- ggplot(cat_counts, aes(x=reorder(category,n), y=n, fill=category)) +
    geom_col(alpha=0.85, show.legend=FALSE) +
    geom_text(aes(label=sprintf("%d\n(%.1f%%)", n, pct)), vjust=-0.3, size=3) +
    scale_fill_manual(values=c("Both"="#4DAF4A","VFDB only"="#E41A1C",
      "PPanGGOLiN only"="#377EB8","Neither"="grey70")) +
    labs(title="L. kps detection discrepancy", x="", y="# Genomes") +
    coord_cartesian(ylim=c(0, max(cat_counts$n)*1.15))
}
if (is.null(pL)) {
  pL <- ggplot() + annotate("text", x=0.5, y=0.5, label="PP data not available", size=4) +
    theme_void() + labs(title="L. kps validation (data unavailable)")
}

# --- Compose Page 2 ---
page2 <- (pG + pH) / (pI + pJ) / (pK + pL) +
  plot_layout(heights=rep(1,3)) +
  plot_annotation(tag_levels="A",
    title="Supplementary Figure S2. Genomic Context and kps Validation") &
  theme(plot.tag=element_text(face="bold"))
print(page2)

# ============================================================
# PAGE 3: REMAINING ANALYSES (panels M-R)
# ============================================================
cat("\n=== PAGE 3: Allelic Conversion, Plasmid Context, ST10 ARG ===\n")

# --- Panel M: Allelic conversion scatter (VFDB vs PP delta) ---
s2_file <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_rigorous_analysis",
  "VFDB_vs_PP_rigorous.csv")
if (!file.exists(s2_file)) {
  s2_file <- file.path(config$OUTPUT_DIR, "ST69", "reviewer_capsule_classification",
    "VFDB_vs_PP_rigorous.csv")
}
pM <- NULL
if (file.exists(s2_file)) {
  s2_data <- read_csv(s2_file, show_col_types=FALSE)
  vfdb_col <- grep("vfdb_delta|VFDB_delta|vfdb_delta", colnames(s2_data), value=TRUE)[1]
  pp_col <- grep("hi_pp_delta|PP_delta|pp_delta", colnames(s2_data), value=TRUE)[1]
  gene_col <- grep("gene", colnames(s2_data), value=TRUE)[1]
  verdict_col <- grep("verdict", colnames(s2_data), value=TRUE)[1]
  if (!is.na(vfdb_col) && !is.na(pp_col)) {
    plot_s2 <- s2_data %>% filter(!is.na(.data[[pp_col]])) %>%
      mutate(label=case_when(
        .data[[verdict_col]]=="REAL" ~ "Real increase",
        .data[[verdict_col]]=="Allele replacement" ~ "Allele replacement",
        TRUE ~ "Other"),
        label=factor(label, levels=c("Real increase","Allele replacement","Other")))
    pM <- ggplot(plot_s2, aes(x=.data[[vfdb_col]], y=.data[[pp_col]], color=label)) +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
      geom_hline(yintercept=0, linetype="dotted", color="grey70") +
      geom_vline(xintercept=0, linetype="dotted", color="grey70") +
      geom_point(size=2.5, alpha=0.8) +
      ggrepel::geom_text_repel(aes(label=.data[[gene_col]]), size=2.5, max.overlaps=20) +
      scale_color_manual(values=c("Real increase"="#1b7837",
        "Allele replacement"="#d73027","Other"="grey60")) +
      labs(title="M. Allelic conversion (VFDB vs PPanGGOLiN)",
           x="VFDB delta (pp)", y="PPanGGOLiN delta (pp)", color="") +
      coord_fixed()
  }
}
if (is.null(pM)) {
  pM <- ggplot() + annotate("text", x=0.5, y=0.5, label="Allelic conversion data not found", size=4) +
    theme_void() + labs(title="M. Allelic conversion (data unavailable)")
}

# --- Panel N: Plasmid vs chromosomal location ---
pl_file <- file.path(config$BASE_DIR, "..", "Plasmid", "plasmid_vfdb_summary",
  "ST69_plasmid_summary.tsv")
pN <- NULL
if (file.exists(pl_file)) {
  pl_vf <- fread(pl_file, sep="\t", header=TRUE, data.table=FALSE,
    na.strings=c("",".","-","NA"), colClasses="character")
  clean_gid <- function(raw) {
    gid <- basename(raw)
    gid <- sub("^Escherichia_coli_","",gid)
    gid <- sub("_combined_plasmids_vfdb\\.tsv$","",gid)
    gid <- sub("_vfdb\\.tsv$","",gid)
    sub("\\.tsv$","",gid)
  }
  pl_vf <- pl_vf %>% mutate(gid_clean=clean_gid(.data[["#FILE"]]))
  top_genes <- c("kpsM","kpsT","kpsD","kpsF","kpsC","kpsE","papB","papF","papX","sat")
  pl_genes <- intersect(top_genes, setdiff(colnames(pl_vf), c("#FILE","NUM_FOUND","gid_clean")))
  if (length(pl_genes) > 0) {
    pl_bin <- pl_vf %>% select(gid_clean, all_of(pl_genes)) %>%
      mutate(across(all_of(pl_genes), to_bin, .names="{.col}_pl"))
    master_clean <- st69_master %>%
      mutate(gid_clean=make_genome_id(genome_id))
    merged_loc <- master_clean %>%
      select(gid_clean, all_of(intersect(top_genes, colnames(master_clean)))) %>%
      mutate(across(all_of(intersect(top_genes, colnames(master_clean))),
        to_bin, .names="{.col}_wg")) %>%
      left_join(pl_bin, by="gid_clean", suffix=c("",".pl_only"))
    gene_loc <- map_dfr(top_genes, function(g) {
      wg <- paste0(g, "_wg"); pl <- paste0(g, "_pl")
      if (!wg%in%colnames(merged_loc) || !pl%in%colnames(merged_loc)) return(NULL)
      merged_loc %>% summarise(n_total=n(), n_wg=sum(.data[[wg]]==1,na.rm=TRUE),
        n_pl=sum(.data[[pl]]==1,na.rm=TRUE)) %>%
        mutate(gene=g, pct_wg=round(n_wg/n_total*100,1),
          pct_pl=round(n_pl/n_total*100,1),
          pct_chr=round((n_wg-n_pl)/n_total*100,1),
          pct_abs=round(100-pct_wg,1))
    })
    if (nrow(gene_loc) > 0) {
      loc_data <- gene_loc %>%
        select(gene, chromosomal=pct_chr, plasmid=pct_pl, absent=pct_abs) %>%
        pivot_longer(-gene, names_to="location", values_to="pct") %>%
        mutate(location=factor(location, levels=c("absent","chromosomal","plasmid")),
               gene=factor(gene, levels=rev(top_genes)))
      pN <- ggplot(loc_data, aes(x=gene, y=pct, fill=location)) +
        geom_col(position="stack", width=0.7, color="grey30", linewidth=0.2) +
        scale_fill_manual(values=c("chromosomal"="#2c7bb6","plasmid"="#d7191c","absent"="grey90"),
          labels=c("Chromosomal","Plasmid","Absent")) +
        labs(title="N. Genomic location (top VF genes)", y="% genomes", x="", fill="") +
        scale_y_continuous(expand=c(0,0), limits=c(0,105)) + coord_flip()
    }
  }
}
if (is.null(pN)) {
  pN <- ggplot() + annotate("text", x=0.5, y=0.5, label="Plasmid data not available", size=4) +
    theme_void() + labs(title="N. Plasmid context (data unavailable)")
}

# --- Panel O: Plasmid RGP co-occurrence heatmap ---
pO <- NULL
if (exists("gene_loc") && nrow(gene_loc) > 0 && file.exists(rgp_vf_file)) {
  rgp_vf <- fread(rgp_vf_file, sep=",", header=TRUE, data.table=FALSE, nrows=2e6)
  rgp_top <- rgp_vf %>% filter(gene_name%in%top_genes) %>%
    select(rgp_id, gene_name) %>% distinct()
  rgp_top_genes <- unique(rgp_top$gene_name)
  if (length(rgp_top_genes) >= 2) {
    gene_pairs_local <- map_dfr(combn(rgp_top_genes, 2, simplify=FALSE), ~{
      s1 <- .x[1]; s2 <- .x[2]
      rgps_s1 <- unique(rgp_top$rgp_id[rgp_top$gene_name==s1])
      rgps_s2 <- unique(rgp_top$rgp_id[rgp_top$gene_name==s2])
      n_both <- length(intersect(rgps_s1, rgps_s2))
      tibble(gene1=s1, gene2=s2, n_rgp_both=n_both,
        n_rgp_g1=length(rgps_s1), n_rgp_g2=length(rgps_s2))
    }) %>% mutate(pct_of_g1=round(n_rgp_both/n_rgp_g1*100,1))
    heat_data <- gene_pairs_local %>%
      bind_rows(gene_pairs_local %>% rename(gene1=gene2, gene2=gene1) %>%
        mutate(pct_of_g1=round(n_rgp_both/n_rgp_g2*100,1))) %>%
      bind_rows(tibble(gene1=rgp_top_genes, gene2=rgp_top_genes,
        n_rgp_both=NA, n_rgp_g1=NA, n_rgp_g2=NA, pct_of_g1=100))
    pO <- ggplot(heat_data, aes(x=gene1, y=gene2, fill=pct_of_g1)) +
      geom_tile(color="white", linewidth=0.3) +
      geom_text(aes(label=ifelse(is.na(pct_of_g1),"",sprintf("%.0f%%",pct_of_g1))), size=2.5) +
      scale_fill_gradient(low="white", high="#2c7bb6", na.value="grey90", name="% co-occurring") +
      labs(title="O. RGP co-localization", x="", y="") +
      theme_minimal(base_size=9) +
      theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, face="bold", size=rel(0.7)),
        axis.text.y=element_text(face="bold", size=rel(0.7)))
  }
}
if (is.null(pO)) {
  pO <- ggplot() + annotate("text", x=0.5, y=0.5, label="RGP data not available", size=4) +
    theme_void() + labs(title="O. RGP co-localization (data unavailable)")
}

# --- Panel P: C3 system profile temporal change ---
pP <- NULL
if (exists("sys_matrix") && nrow(sys_matrix) > 0) {
  c3_sys <- sys_matrix %>%
    left_join(st69_master %>% select(genome_id, year, shell_cluster), by=c("genome"="genome_id")) %>%
    filter(shell_cluster=="Cluster_3", !is.na(year)) %>%
    mutate(period=ifelse(year%in%EARLY,"early",ifelse(year%in%LATE,"late",NA))) %>%
    filter(!is.na(period))
  c3_profiles <- c3_sys %>%
    group_by(genome) %>%
    summarise(profile=paste(sort(unique(system[present])), collapse=","),
              n_systems=sum(present), .groups="drop") %>%
    mutate(profile=ifelse(n_systems==0,"none",profile))
  c3_sys_with_prof <- c3_sys %>% inner_join(c3_profiles, by="genome")
  profile_delta <- c3_sys_with_prof %>%
    group_by(period, profile) %>% summarise(n=n(), .groups="drop") %>%
    group_by(period) %>% mutate(pct=n/sum(n)*100) %>%
    select(period, profile, pct) %>%
    pivot_wider(names_from=period, values_from=pct, values_fill=0) %>%
    mutate(delta=round(late-early,1)) %>% arrange(desc(abs(delta))) %>% head(8)
  if (nrow(profile_delta) > 0) {
    pP <- ggplot(profile_delta, aes(x=reorder(profile, delta), y=delta, fill=delta>0)) +
      geom_col(alpha=0.85) + coord_flip() +
      scale_fill_manual(values=c("TRUE"="#d7191c","FALSE"="#2c7bb6"), guide="none") +
      labs(title="P. C3 profile temporal change", x="", y="Delta prevalence (pp)") +
      geom_hline(yintercept=0, linetype="dashed", color="grey50")
  }
}
if (is.null(pP)) {
  pP <- ggplot() + annotate("text", x=0.5, y=0.5, label="Data not available", size=4) +
    theme_void() + labs(title="P. C3 profile change (data unavailable)")
}

# --- Panel Q: ST10 ARG decomposition ---
st10_card_file <- config$st_card_burden("ST10")
st10_meta_file <- config$st_metadata("ST10")
pQ <- NULL
if (file.exists(st10_card_file) && file.exists(st10_meta_file)) {
  card_raw <- fread(st10_card_file, sep="\t", header=TRUE, data.table=FALSE,
    colClasses="character", na.strings="")
  card_raw$genome_id <- make_genome_id(str_remove(card_raw[[1]], "^card/ST10/"))
  card_raw$genome_id <- str_remove(card_raw$genome_id, "_card\\.tsv$")
  arg_genes <- setdiff(colnames(card_raw), c("#FILE","NUM_FOUND","genome_id"))
  card_bin <- card_raw %>% select(genome_id, all_of(arg_genes)) %>%
    mutate(across(all_of(arg_genes), to_bin))
  card_bin$total_arg <- rowSums(card_bin[, arg_genes, drop=FALSE], na.rm=TRUE)
  meta10 <- readxl::read_excel(st10_meta_file, col_types="text")
  name_col <- intersect(c("Name","genome","strain","isolate","assembly","sample"), colnames(meta10))
  year_col <- intersect(c("Collection Year","Collection_Year","year","Year","date"), colnames(meta10))
  source_col <- intersect(c("Source Niche","Source_Niche","source_niche","source","niche","Host","host"), colnames(meta10))
  country_col <- intersect(c("Country","country"), colnames(meta10))
  if (length(name_col)>0 && length(year_col)>0) {
    meta_clean <- meta10 %>% transmute(
      name_raw=.data[[name_col[1]]],
      year=suppressWarnings(as.integer(.data[[year_col[1]]])),
      source_niche=if(length(source_col)>0) as.character(.data[[source_col[1]]]) else NA_character_,
      country=if(length(country_col)>0) as.character(.data[[country_col[1]]]) else NA_character_) %>%
      mutate(genome_id=make_genome_id(name_raw))
    merged10 <- card_bin %>% inner_join(meta_clean, by="genome_id") %>%
      filter(!is.na(year), year>=2016, year<=2025)
    if (nrow(merged10) > 0) {
      merged10 <- merged10 %>% mutate(
        clinical=case_when(
          grepl("human|patient|clinical|blood|urine|wound", source_niche, ignore.case=TRUE) ~ "Clinical",
          is.na(source_niche)|source_niche=="" ~ "Unknown", TRUE ~ "Non-clinical"),
        region=case_when(
          grepl("Europe|UK|Germany|France|Spain|Italy|Netherlands", country, ignore.case=TRUE) ~ "Europe",
          grepl("USA|United States|Canada", country, ignore.case=TRUE) ~ "N. America",
          grepl("China|Japan|India|South Korea|Taiwan", country, ignore.case=TRUE) ~ "Asia",
          TRUE ~ "Other"))
      annual <- merged10 %>% group_by(year) %>%
        summarise(n=n(), mean_arg=mean(total_arg),
          se_arg=sd(total_arg)/sqrt(n()), .groups="drop")
      annual_clin <- merged10 %>% filter(clinical!="Unknown") %>%
        group_by(year, clinical) %>%
        summarise(n=n(), mean_arg=mean(total_arg), .groups="drop")
      m1 <- lm(total_arg ~ year, data=merged10)
      m2 <- lm(total_arg ~ year + clinical, data=merged10 %>% filter(clinical!="Unknown"))
      m3 <- lm(total_arg ~ year + region, data=merged10 %>% filter(region!="Other"))
      m4 <- lm(total_arg ~ year + clinical + region,
        data=merged10 %>% filter(clinical!="Unknown", region!="Other"))
      model_comp <- bind_rows(
        tidy(m1) %>% filter(term=="year") %>% mutate(model="Raw"),
        tidy(m2) %>% filter(term=="year") %>% mutate(model="+ Clinical"),
        tidy(m3) %>% filter(term=="year") %>% mutate(model="+ Region"),
        tidy(m4) %>% filter(term=="year") %>% mutate(model="+ Both"))

      pQ_top <- ggplot(annual, aes(x=year, y=mean_arg)) +
        geom_ribbon(aes(ymin=mean_arg-se_arg, ymax=mean_arg+se_arg), alpha=0.2, fill="steelblue") +
        geom_line(color="steelblue", linewidth=0.8) + geom_point(aes(size=n), color="steelblue") +
        scale_size_continuous(range=c(1.5,4), guide="none") +
        labs(title="Q. ST10 mean ARG burden", x="Year", y="Mean ARG count")

      pQ_bot <- ggplot(model_comp, aes(x=estimate, y=reorder(model, estimate))) +
        geom_point(size=3) + geom_errorbarh(aes(xmin=estimate-1.96*std.error,
          xmax=estimate+1.96*std.error), height=0.2) +
        geom_vline(xintercept=0, linetype="dashed", color="grey50") +
        labs(title="R. ARG trend (adjusted)", x="Year coefficient", y="")

      pQ <- pQ_top / pQ_bot + plot_layout(heights=c(1,1))
    }
  }
}
if (is.null(pQ)) {
  pQ <- ggplot() + annotate("text", x=0.5, y=0.5, label="ST10 data not available", size=4) +
    theme_void() + labs(title="Q-R. ST10 ARG decomposition (data unavailable)")
}

# --- Compose Page 3 ---
page3 <- (pM + pN) / (pO + pP) / pQ +
  plot_layout(heights=c(1,1,1)) +
  plot_annotation(tag_levels="A",
    title="Supplementary Figure S3. Allelic Conversion, Plasmid Context, and ST10 ARG Decomposition") &
  theme(plot.tag=element_text(face="bold"))
print(page3)

# ============================================================
# Close PDF
# ============================================================
dev.off()

cat("\n============================================\n")
cat(" DONE\n")
cat(" Output:", pdf_file, "\n")
cat("============================================\n")

})  # end block 16: 22_fig_supplementary_combined.R
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))


################################################################################
# BLOCK 17/17: 12_fig02-05_ST69_analysis.R (FULL ORIGINAL SCRIPT, verbatim)
# Kept in full for completeness. Its "FIGURE 2" + "FIGURE 4" sections and
# their "fig2_3_combined" step ARE the source of submitted Figure 2 (see
# Analysis/Figure2.R, which extracts just those sections). Its "FIGURE 5"
# section IS the source of submitted Figure 4 (see Analysis/Figure4.R).
# Both are ALSO reproduced here as part of the unmodified original script,
# alongside its other content that is not part of any submitted main
# figure (internal "Figure 3" PCA panel, "Figure 6"/"Figure 9" file-copy
# stubs that duplicate 14_fig06_gene_trajectories.R / 08_clinical_
# enrichment_3panel.R above, "Figure 7" preliminary/simple tree mapping,
# and "Figure 8" clinical trajectory-by-cluster). Running this block after
# Figure2.R/Figure4.R have already run is redundant but harmless -- it
# simply recomputes and re-saves the same Fig2/Fig4 outputs under this
# script's own "figures_all_3panel" output folder in addition.
################################################################################
tryCatch({
local({
#!/usr/bin/env Rscript
#
# All publication figures (2-9) as 3-panel comparison
# Panel A: ST69 VFDB | Panel B: ST10 VirulenceFinder | Panel C: ST10 ResFinder
#
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(ggplot2); library(patchwork)
  library(scales)
  if (requireNamespace("ape", quietly=TRUE)) library(ape)
  if (requireNamespace("ggtree", quietly=TRUE)) library(ggtree)
})
source("config.R")

OUT <- file.path(config$OUTPUT_DIR, "figures_all_3panel")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_classic(base_size = 10))

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

cluster_colors <- function(n) {
  if (n <= 9) RColorBrewer::brewer.pal(max(n, 3), "Set1") else rainbow(n)
}

cat("Loading metadata...\n")
meta69 <- read_xlsx(config$st_metadata("ST69")) %>%
  rename(genome_id = Name) %>%
  mutate(genome_id = make_genome_id(genome_id),
         source_niche = .data[["Source Niche"]],
         clinical = source_niche == "Human",
         year = as.integer(.data[["Collection Year"]])) %>%
  filter(!is.na(year), !is.na(source_niche))
meta10 <- read_xlsx(config$st_metadata("ST10")) %>%
  rename(genome_id = Name) %>%
  mutate(genome_id = make_genome_id(genome_id),
         source_niche = .data[["Source Niche"]],
         clinical = source_niche == "Human",
         year = as.integer(.data[["Collection Year"]])) %>%
  filter(!is.na(year), !is.na(source_niche))

# ===== Load data per case =====
# Case A: ST69 VFDB
cat("Loading ST69 VFDB master...\n")
st69_master <- read.csv("output/ST69/vfdb_analysis/04_master_shell_cluster_metadata_VFDB_table.csv",
                         stringsAsFactors = FALSE) %>%
  mutate(genome_id = make_genome_id(genome_id))
# Master already has clinical_binary, year, source_niche columns
st69_vf_genes <- setdiff(colnames(st69_master),
                          c("genome_id","shell_cluster","raw_cluster","X.FILE","total_vf",
                            "clinical_binary","year","source_niche","continent","country","niche","genome"))
# Convert VFDB gene columns to binary presence
for (g in st69_vf_genes) {
  st69_master[[g]] <- as.integer(!is.na(st69_master[[g]]) & st69_master[[g]] != "." & st69_master[[g]] != "")
}
st69_master <- st69_master %>% distinct(genome_id, .keep_all = TRUE)
st69_master$year <- as.integer(st69_master$year)
st69_sil <- read.csv("output/ST69/vfdb_analysis/02_silhouette_k4_to_k10.csv", stringsAsFactors = FALSE)
st69_ksel <- read.csv("output/ST69/vfdb_analysis/02B_k_selection_summary.csv", stringsAsFactors = FALSE)

# Case B: ST10 VF
cat("Loading ST10 VF master...\n")
st10_vf_master <- read.csv("output/ST10/virulencefinder_validation/04_master_shell_cluster_metadata_VF_table.csv",
                            stringsAsFactors = FALSE) %>%
  mutate(genome_id = make_genome_id(genome_id))
st10_vf_master$year <- as.integer(st10_vf_master$year)
# clinical_binary already exists; convert to integer
st10_vf_master <- st10_vf_master %>%
  mutate(clinical_binary = as.integer(clinical_binary == 1 | clinical_binary == "Human" | clinical_binary == TRUE))
st10_vf_sil <- read.csv("output/ST10/virulencefinder_validation/02_silhouette_k4_to_k10.csv", stringsAsFactors = FALSE)
st10_vf_ksel <- read.csv("output/ST10/virulencefinder_validation/02B_k_selection_summary.csv", stringsAsFactors = FALSE)

# Case C: ST10 ResFinder
cat("Loading ST10 ResFinder data...\n")
st10_resf_master <- read.csv("output/ST10/vfdb_analysis/04_master_shell_cluster_metadata_VFDB_table.csv",
                              stringsAsFactors = FALSE) %>%
  mutate(genome_id = make_genome_id(genome_id)) %>%
  left_join(meta10 %>% select(genome_id, clinical, year), by = "genome_id") %>%
  mutate(clinical_binary = as.integer(clinical))
# Actually this is VFDB data, not ResFinder. For ResFinder we need the resfinder binary matrix
# and then join with VFDB shell cluster assignments (same pangenome)
st10_shell_clusters <- read.csv("output/ST10/vfdb_analysis/03_shell_gene_cluster_assignments_k9.csv",
                                 stringsAsFactors = FALSE) %>% rename(genome_id = 1, shell_cluster = 2, raw_cluster = 3)
st10_resf_sil <- read.csv("output/ST10/resistance_analysis/resfinder/02_silhouette_k4_to_k10.csv", stringsAsFactors = FALSE)
st10_resf_ksel <- read.csv("output/ST10/resistance_analysis/resfinder/02B_k_selection_summary.csv", stringsAsFactors = FALSE)

resf_binary <- read.delim(file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
                           header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
resf10 <- resf_binary %>% filter(st == "ST10") %>%
  mutate(genome_id = make_genome_id(genome)) %>%
  select(-st, -genome)
resf10_genes <- setdiff(colnames(resf10), "genome_id")
resf10_long <- resf10 %>% pivot_longer(-genome_id, names_to = "gene", values_to = "present") %>%
  mutate(present = as.integer(!is.na(present) & present != "" & present != "0"))
# Compute total ARG per genome
resf10_total <- resf10_long %>% group_by(genome_id) %>%
  summarise(total_arg = sum(present), .groups = "drop")

# Join with metadata and cluster assignments
st10_resf_master <- resf10_total %>% left_join(meta10 %>% select(genome_id, clinical, year), by = "genome_id") %>%
  left_join(st10_shell_clusters %>% select(genome_id, shell_cluster), by = "genome_id") %>%
  filter(!is.na(shell_cluster))

cat("Data loaded.\n")
cat("  ST69:", nrow(st69_master), "genomes,", length(st69_vf_genes), "VFDB genes\n")
cat("  ST10 VF:", nrow(st10_vf_master), "genomes\n")
cat("  ST10 ResF:", nrow(st10_resf_master), "genomes\n")

# ================================================================
# FIGURE 2: Silhouette
# ================================================================
cat("\n--- Figure 2: Silhouette ---\n")

plot_silhouette <- function(sil, ksel, title) {
  # Find best k from ksel (different column names across outputs)
  ksel_cols <- colnames(ksel)
  best_k <- NA
  if ("best_k_primary" %in% ksel_cols) best_k <- ksel$best_k_primary[1]
  else if ("best_k" %in% ksel_cols) best_k <- ksel$best_k[1]
  else if ("selected_k" %in% ksel_cols) best_k <- ksel$selected_k[1]
  if (is.na(best_k)) return(ggplot() + annotate("text", x=0.5, y=0.5, label="Could not determine best k") +
                              labs(title=title) + theme_void())

  # Silhouette data by k
  sil_cols <- colnames(sil)
  k_col <- if ("k" %in% sil_cols) "k" else "K"
  sil_col <- if ("average_silhouette_width" %in% sil_cols) "average_silhouette_width" else
    if ("avg_sil" %in% sil_cols) "avg_sil" else sil_cols[3]
  method_col <- if ("method" %in% sil_cols) "method" else NULL

  # Filter to ward.D2 method if available; else first method
  sil <- sil %>% rename(kk = !!k_col, sil_val = !!sil_col)
  if (!is.null(method_col)) {
    sil <- sil %>% rename(m = !!method_col)
    if ("ward.D2" %in% sil$m) sil <- sil %>% filter(m == "ward.D2")
  }

  if (n_distinct(sil$kk) > 1) {
    p1 <- sil %>% mutate(best = kk == best_k) %>%
      ggplot(aes(x = kk, y = sil_val)) +
      geom_line(color = "grey60") + geom_point(aes(color = best), size = 2.5) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
      labs(x = "k", y = "Mean silhouette") + ylim(0, NA)
  } else {
    # Only one k value available
    p1 <- sil %>% mutate(best = TRUE) %>%
      ggplot(aes(x = factor(kk), y = sil_val)) +
      geom_col(aes(fill = best), show.legend = FALSE) +
      scale_fill_manual(values = c("TRUE" = "red")) +
      labs(x = "Selected k", y = "Mean silhouette") + ylim(0, NA)
  }

  # Per-k silhouette not available (only average per k), so show bar of best_k vs others
  p2 <- sil %>% mutate(label = ifelse(kk == best_k, paste0("k=", kk, " (selected)"), paste0("k=", kk))) %>%
    ggplot(aes(x = sil_val, y = reorder(label, sil_val), fill = kk == best_k)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey60"), guide = "none") +
    labs(x = "Mean silhouette", y = NULL) + xlim(0, NA)

  wrap_elements(p1 | p2) + labs(title = title) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

p2A <- plot_silhouette(st69_sil, st69_ksel, "ST69 VFDB")
p2C <- plot_silhouette(st10_vf_sil, st10_vf_ksel, "ST10 shell clusters")
fig2 <- (p2A | p2C) + plot_layout(ncol = 2) +
  plot_annotation(title = "Shell-gene cluster silhouette analysis",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))
ggsave(file.path(OUT, "Fig2_silhouette_3panel.png"), fig2, width = 12, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig2_silhouette_3panel.pdf"), fig2, width = 12, height = 5, bg = "white")
cat("  Saved Fig2\n")

# ================================================================
# FIGURE 3: PCA of gene content
# ================================================================
cat("\n--- Figure 3: PCA ---\n")

pca_plot <- function(gene_mat, clust_df, title) {
  mat <- as.matrix(gene_mat)
  keep <- complete.cases(mat) & apply(mat, 1, var) > 0
  if (sum(keep) < 10) return(NULL)
  mat <- mat[keep, , drop = FALSE]
  vars <- apply(mat, 2, var, na.rm = TRUE)
  mat <- mat[, vars > 0.01, drop = FALSE]
  if (ncol(mat) < 3) return(NULL)
  pca <- prcomp(mat, scale. = TRUE, center = TRUE)
  ve <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    genome_id = rownames(pca$x)) %>%
    left_join(clust_df, by = "genome_id")
  ggplot(df, aes(x = PC1, y = PC2, color = shell_cluster)) +
    geom_point(alpha = 0.6, size = 0.8) +
    stat_ellipse(aes(fill = shell_cluster), geom = "polygon", alpha = 0.06, level = 0.8, show.legend = FALSE) +
    scale_color_manual(values = cluster_colors(n_distinct(df$shell_cluster))) +
    labs(title = title, x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"), color = "Cluster") +
    theme_classic(base_size = 10)
}

# ST69: use VFDB genes as the matrix (full pangenome not needed)
st69_pca_mat <- st69_master %>% select(genome_id, all_of(st69_vf_genes)) %>%
  distinct(genome_id, .keep_all = TRUE) %>%
  column_to_rownames("genome_id")
st69_clust <- st69_master %>% select(genome_id, shell_cluster) %>% distinct(genome_id, .keep_all = TRUE)
p3A <- pca_plot(st69_pca_mat, st69_clust, "ST69 VFDB")
if (is.null(p3A)) p3A <- ggplot() + annotate("text", x=0.5, y=0.5, label="PCA failed") + labs(title="ST69 VFDB") + theme_void()

# ST10 VF: use binary VF matrix
vf_binary <- read.delim(config$VF_BINARY, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
vf10_mat <- vf_binary %>% filter(st == "ST10") %>%
  mutate(genome_id = make_genome_id(genome)) %>% select(-st, -genome) %>%
  distinct(genome_id, .keep_all = TRUE) %>% column_to_rownames("genome_id")
vf10_mat[] <- lapply(vf10_mat, function(x) as.integer(!is.na(x) & x != "" & x != "0"))
st10_vf_clust <- st10_vf_master %>% select(genome_id, shell_cluster)
p3B <- pca_plot(vf10_mat, st10_vf_clust, "ST10 VirulenceFinder")
if (is.null(p3B)) p3B <- ggplot() + annotate("text", x=0.5, y=0.5, label="PCA failed") + labs(title="ST10 VirulenceFinder") + theme_void()

# ST10 ResFinder: use binary ARG matrix
resf10_mat <- resf10 %>% select(genome_id, all_of(resf10_genes)) %>%
  distinct(genome_id, .keep_all = TRUE) %>% column_to_rownames("genome_id")
resf10_mat[] <- lapply(resf10_mat, function(x) as.integer(!is.na(x) & x != "" & x != "0"))
st10_resf_clust <- st10_resf_master %>% select(genome_id, shell_cluster)
p3C <- pca_plot(resf10_mat, st10_resf_clust, "ST10 ResFinder")
if (is.null(p3C)) p3C <- ggplot() + annotate("text", x=0.5, y=0.5, label="PCA failed") + labs(title="ST10 ResFinder") + theme_void()

fig3 <- (p3A | p3B | p3C) + plot_layout(ncol = 3) +
  plot_annotation(title = "Gene content PCA by shell cluster",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))
ggsave(file.path(OUT, "Fig3_pca_3panel.png"), fig3, width = 18, height = 6, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig3_pca_3panel.pdf"), fig3, width = 18, height = 6, bg = "white")
cat("  Saved Fig3\n")

# ================================================================
# FIGURE 4: Gene content by cluster
# ================================================================
cat("\n--- Figure 4: Content by cluster ---\n")

content_plot <- function(d, val_col, ylab, title) {
  summ <- d %>% group_by(shell_cluster) %>%
    summarise(n = n(), mean_val = mean(.data[[val_col]], na.rm = TRUE),
              total = sum(.data[[val_col]], na.rm = TRUE), .groups = "drop") %>% drop_na()
  ggplot(summ, aes(x = reorder(shell_cluster, desc(n)), y = total, fill = shell_cluster)) +
    geom_col() +
    geom_text(aes(label = n, y = 0), vjust = 1.5, size = 3, color = "white") +
    scale_fill_manual(values = cluster_colors(nrow(summ)), guide = "none") +
    labs(title = title, x = NULL, y = ylab) +
    theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p4A <- content_plot(st69_master, "total_vf", "Total VFDB genes", "ST69 VFDB")
p4B <- content_plot(st10_vf_master, "total_vf", "Total VF genes", "ST10 VirulenceFinder")
p4C <- content_plot(st10_resf_master, "total_arg", "Total ARG count", "ST10 ResFinder")

fig4 <- (p4A | p4B | p4C) + plot_layout(ncol = 3) +
  plot_annotation(title = "Gene content by shell cluster",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))
ggsave(file.path(OUT, "Fig4_content_by_cluster_3panel.png"), fig4, width = 14, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig4_content_by_cluster_3panel.pdf"), fig4, width = 14, height = 5, bg = "white")
cat("  Saved Fig4\n")

# ---- Combined Fig 2 + Fig 3 (silhouette + content) ----
fig23 <- (fig2 / fig4) + plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Shell-gene cluster silhouette and gene content",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave(file.path(OUT, "Fig2_3_combined.png"), fig23, width = 12, height = 11, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig2_3_combined.pdf"), fig23, width = 12, height = 11, bg = "white")
cat("  Saved Fig2_3_combined\n")

# ================================================================
# FIGURE 5: Oaxaca decomposition
# ================================================================
cat("\n--- Figure 5: Oaxaca decomposition ---\n")

oaxaca_data <- function(df, val_col) {
  df <- df %>% filter(!is.na(year), !is.na(.data[[val_col]]), !is.na(shell_cluster))
  yr_cl <- df %>% group_by(year, shell_cluster) %>%
    summarise(n = n(), mean_val = mean(.data[[val_col]], na.rm = TRUE), .groups = "drop") %>%
    group_by(year) %>%
    mutate(total_n = sum(n), prop = n / total_n) %>% ungroup()
  n_clusters <- n_distinct(yr_cl$shell_cluster)
  yr_valid <- yr_cl %>% group_by(year) %>%
    summarise(all_present = n_distinct(shell_cluster) == n_clusters, .groups = "drop") %>%
    filter(all_present)
  if (nrow(yr_valid) < 2) return(NULL)
  ref_year <- min(yr_valid$year)
  ref <- yr_cl %>% filter(year == ref_year) %>%
    select(shell_cluster, ref_prop = prop, ref_mean = mean_val)
  decomp <- yr_cl %>% left_join(ref, by = "shell_cluster") %>%
    mutate(delta_prop = prop - ref_prop, delta_mean = mean_val - ref_mean,
           composition = delta_prop * ref_mean,
           within = ref_prop * delta_mean,
           interaction = delta_prop * delta_mean)
  total_decomp <- decomp %>% group_by(year) %>%
    summarise(composition = sum(composition, na.rm = TRUE),
              within = sum(within, na.rm = TRUE),
              interaction = sum(interaction, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(-year, names_to = "component", values_to = "value") %>%
    mutate(component = str_to_title(component),
           component = case_when(
             component == "Composition" ~ "Composition",
             component == "Within" ~ "Within-cluster",
             component == "Interaction" ~ "Interaction"),
           component = factor(component, levels = c("Composition", "Within-cluster", "Interaction")))
  obs <- df %>% group_by(year) %>%
    summarise(observed = mean(.data[[val_col]], na.rm = TRUE), .groups = "drop")
  ref_obs <- obs$observed[obs$year == ref_year]
  obs <- obs %>% mutate(delta_from_ref = observed - ref_obs)
  list(components = total_decomp, observed = obs, ref_year = ref_year,
       clusters = decomp)
}

plot_obs <- function(d, title) {
  if (is.null(d)) return(NULL)
  d$observed %>%
    ggplot(aes(x = year, y = observed)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.8, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.6, color = "grey40") +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(title = title, x = NULL, y = "Mean burden", subtitle = paste("Ref:", d$ref_year)) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"))
}
plot_comp <- function(d, title) {
  if (is.null(d)) return(NULL)
  cc <- c("Composition" = "#4DAF4A", "Within-cluster" = "#377EB8", "Interaction" = "#E41A1C")
  d$components %>%
    ggplot(aes(x = year, y = value, color = component, group = component)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5, alpha = 0.85) +
    scale_color_manual(values = cc) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(title = title, x = NULL, y = "Change in weighted mean", color = NULL) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 9), legend.key.size = unit(0.8, "lines"),
          legend.position = "right")
}
plot_drv <- function(d, title, show_smooth = TRUE) {
  if (is.null(d)) return(NULL)
  n_cl <- n_distinct(d$clusters$shell_cluster)
  pal <- if (n_cl <= 9) RColorBrewer::brewer.pal(max(n_cl, 3), "Set1") else rainbow(n_cl)
  p <- d$clusters %>%
    ggplot(aes(x = year, y = within, color = shell_cluster, group = shell_cluster)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.3, alpha = 0.8) +
    scale_color_manual(values = pal) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(title = title, x = "Year", y = "Within-cluster contribution", color = "Cluster") +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 9), legend.key.size = unit(0.8, "lines")) +
    guides(color = guide_legend(ncol = 3))
  if (show_smooth) {
    p <- p + geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5)
  }
  p
}

d69 <- oaxaca_data(st69_master, "total_vf")
d10vf <- oaxaca_data(st10_vf_master, "total_vf")
d10res <- oaxaca_data(st10_resf_master, "total_arg")

lbls <- c("ST69 VFDB", "ST10 VirulenceFinder", "ST10 ResFinder")
ds <- list(d69, d10vf, d10res)
plts <- list()
for (i in seq_along(ds)) {
  plts[[paste0("A", i)]] <- plot_obs(ds[[i]], lbls[i])
  plts[[paste0("B", i)]] <- plot_comp(ds[[i]], lbls[i])
  plts[[paste0("C", i)]] <- plot_drv(ds[[i]], lbls[i])
}
rowA <- wrap_plots(plts[c("A1","A2","A3")], ncol = 3) +
  plot_annotation(title = "A  Observed lineage-level mean burden",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))
rowB <- wrap_plots(plts[c("B1","B2","B3")], ncol = 3) +
  plot_annotation(title = "B  Decomposition components",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))
rowC <- wrap_plots(plts[c("C1","C2","C3")], ncol = 3) & theme(legend.position = "bottom")
rowC <- rowC + plot_annotation(title = "C  Per-cluster within-cluster driver effects",
                                theme = theme(plot.title = element_text(size = 11, face = "bold")))
fig5 <- (rowA / rowB / rowC) + plot_layout(heights = c(1, 1, 1.3)) +
  plot_annotation(
    title = "Oaxaca-style decomposition of burden change",
    subtitle = "Composition = cluster shift x ref.mean  |  Within-cluster = ref.prop x delta mean  |  Interaction = joint",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"))
  )
ggsave(file.path(OUT, "Fig5_decomposition_3panel.png"), fig5, width = 12, height = 12, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig5_decomposition_3panel.pdf"), fig5, width = 12, height = 12, bg = "white")
cat("  Saved Fig5\n")

# Version 2: driver panel without smooth lines
plts_ns <- list()
for (i in seq_along(ds)) {
  plts_ns[[paste0("A", i)]] <- plot_obs(ds[[i]], lbls[i])
  plts_ns[[paste0("B", i)]] <- plot_comp(ds[[i]], lbls[i])
  plts_ns[[paste0("C", i)]] <- plot_drv(ds[[i]], lbls[i], show_smooth = FALSE)
}
rowC_ns <- wrap_plots(plts_ns[c("C1","C2","C3")], ncol = 3) & theme(legend.position = "bottom")
rowC_ns <- rowC_ns + plot_annotation(title = "C  Per-cluster within-cluster driver effects (no trend lines)",
                                     theme = theme(plot.title = element_text(size = 11, face = "bold")))
fig5_ns <- (rowA / rowB / rowC_ns) + plot_layout(heights = c(1, 1, 1.3)) +
  plot_annotation(
    title = "Oaxaca-style decomposition of burden change",
    subtitle = "Composition = cluster shift x ref.mean  |  Within-cluster = ref.prop x delta mean  |  Interaction = joint",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"))
  )
ggsave(file.path(OUT, "Fig5_decomposition_3panel_nosmooth.png"), fig5_ns, width = 12, height = 12, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig5_decomposition_3panel_nosmooth.pdf"), fig5_ns, width = 12, height = 12, bg = "white")
cat("  Saved Fig5 (nosmooth)\n")

# ================================================================
# FIGURE 6: Gene prevalence (copy from earlier script)
# ================================================================
cat("\n--- Figure 6 (copy): Gene prevalence ---\n")
f6 <- file.path(config$OUTPUT_DIR, "figures_gene_prevalence_3panel", "Fig_gene_prevalence_3panel.png")
if (file.exists(f6)) file.copy(f6, file.path(OUT, "Fig6_gene_prevalence_3panel.png"), overwrite = TRUE)
f6p <- sub(".png$", ".pdf", f6)
if (file.exists(f6p)) file.copy(f6p, file.path(OUT, "Fig6_gene_prevalence_3panel.pdf"), overwrite = TRUE)
cat("  Copied Fig6\n")

# ================================================================
# FIGURE 7: Tree mapping
# ================================================================
cat("\n--- Figure 7: Tree mapping ---\n")

tree_plot <- function(tree_file, mapping_file, title, layout = "rectangular") {
  if (!requireNamespace("ape", quietly=TRUE) || !requireNamespace("ggtree", quietly=TRUE)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="ape/ggtree not installed") +
             labs(title = title) + theme_void())
  }
  if (!file.exists(tree_file)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Tree not found") +
             labs(title = title) + theme_void())
  }
  tree <- tryCatch(ape::read.tree(tree_file), error = function(e) NULL)
  if (is.null(tree)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Could not read tree") +
             labs(title = title) + theme_void())
  }
  if (file.exists(mapping_file)) {
    mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)
  } else {
    mapping <- NULL
  }
  tips <- data.frame(genome_id = tree$tip.label, stringsAsFactors = FALSE) %>%
    mutate(genome_id = make_genome_id(genome_id))
  if (!is.null(mapping)) {
    if ("shell_cluster" %in% colnames(mapping)) {
      tips <- tips %>% left_join(mapping %>% select(genome_id, shell_cluster), by = "genome_id")
    } else if (ncol(mapping) >= 2) {
      colnames(mapping)[1:2] <- c("genome_id", "shell_cluster")
      tips <- tips %>% left_join(mapping %>% select(genome_id, shell_cluster), by = "genome_id")
    }
  }
  tips <- tips %>% filter(!is.na(shell_cluster))
  if (nrow(tips) < 10) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Insufficient tip data") +
             labs(title = title) + theme_void())
  }
  # Prune tree to only matched tips
  tree <- tryCatch(ape::keep.tip(tree, tips$genome_id), error = function(e) tree)
  tips <- tips %>% slice(match(tree$tip.label, genome_id))
  # Normalize branch lengths so all panels have comparable proportions
  tree$edge.length <- rep(1, nrow(tree$edge))

  all_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                 "#FFFF33", "#A65628", "#F781BF", "#999999")
  names(all_cols) <- paste0("Cluster_", 1:9)
  present <- sort(unique(tips$shell_cluster))
  clust_cols <- all_cols[present]

  p <- ggtree(tree, layout = layout, size = 0.3, color = "grey50") %<+%
    tips + geom_tippoint(aes(color = shell_cluster), size = 1.2, alpha = 0.9) +
    scale_color_manual(values = clust_cols) +
    labs(title = title, color = "Cluster") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 11, face = "bold"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())
  if (layout == "circular") {
    p <- p + theme(plot.margin = margin(10, 10, 10, 10))
  }
  p
}

# ST69: tree from pangenome, mapping from master table
st69_tree <- "pangenome_output/msa_output/phylogeny/ST69_bootstrap.treefile"
st69_map <- "output/ST69/vfdb_analysis/04_master_shell_cluster_metadata_VFDB_table.csv"
p7A <- tree_plot(st69_tree, st69_map, "ST69")
p7A_circ <- tree_plot(st69_tree, st69_map, "ST69", layout = "circular")

# ST10: shell genome clusters (same pangenome structure for VF and ResFinder)
st10_tree <- "output/ST10/vfdb_analysis/pruned_tree.nwk"
st10_map <- "output/ST10/vfdb_analysis/tree_cluster_mapping.csv"
p7B <- tree_plot(st10_tree, st10_map, "ST10")
p7B_circ <- tree_plot(st10_tree, st10_map, "ST10", layout = "circular")

fig7 <- (p7A | p7B) + plot_layout(ncol = 2) +
  plot_annotation(theme = theme(plot.title = element_blank()))
ggsave(file.path(OUT, "Fig7_tree_mapping_2panel.png"), fig7, width = 16, height = 8, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig7_tree_mapping_2panel.pdf"), fig7, width = 16, height = 8, bg = "white")
cat("  Saved Fig7 (linear)\n")

fig7_circ <- (p7A_circ | p7B_circ) + plot_layout(ncol = 2) +
  plot_annotation(theme = theme(plot.title = element_blank()))
ggsave(file.path(OUT, "Fig7_tree_mapping_2panel_circular.png"), fig7_circ, width = 14, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig7_tree_mapping_2panel_circular.pdf"), fig7_circ, width = 14, height = 7, bg = "white")
cat("  Saved Fig7 (circular)\n")

# ================================================================
# FIGURE 8: Clinical trajectory
# ================================================================
cat("\n--- Figure 8: Clinical trajectory ---\n")

clin_traj <- function(d, title) {
  if (!"year" %in% colnames(d) || !"shell_cluster" %in% colnames(d)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Missing columns") +
             labs(title = title) + theme_void())
  }
  if ("clinical_binary" %in% colnames(d)) {
    d <- d %>% mutate(clin_flag = as.integer(clinical_binary == 1 | clinical_binary == "Human" | clinical_binary == TRUE))
  } else if ("clinical" %in% colnames(d)) {
    d <- d %>% mutate(clin_flag = as.integer(clinical))
  } else {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No clinical data") +
             labs(title = title) + theme_void())
  }
  d <- d %>% filter(!is.na(year)) %>%
    group_by(year, shell_cluster) %>%
    summarise(pct_clin = mean(clin_flag, na.rm = TRUE) * 100, n = n(), .groups = "drop")
  # Keep top clusters by total n
  top <- d %>% group_by(shell_cluster) %>% summarise(total = sum(n), .groups = "drop") %>%
    slice_max(total, n = 5)
  d2 <- d %>% filter(shell_cluster %in% top$shell_cluster)
  if (nrow(d2) < 5) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Insufficient data") +
             labs(title = title) + theme_void())
  }
  ggplot(d2, aes(x = year, y = pct_clin, color = shell_cluster, group = shell_cluster)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    scale_color_manual(values = cluster_colors(n_distinct(d2$shell_cluster))) +
    scale_x_continuous(breaks = pretty(d2$year, n = 5)) +
    ylim(0, 100) +
    labs(title = title, x = "Year", y = "Clinical isolates (%)", color = "Cluster") +
    theme_classic(base_size = 10)
}

p8A <- clin_traj(st69_master, "ST69 VFDB")
p8B <- clin_traj(st10_vf_master, "ST10 VirulenceFinder")
p8C <- clin_traj(st10_resf_master, "ST10 ResFinder")

fig8 <- (p8A | p8B | p8C) + plot_layout(ncol = 3) +
  plot_annotation(title = "Clinical trajectory by shell cluster",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))
ggsave(file.path(OUT, "Fig8_clinical_trajectory_3panel.png"), fig8, width = 16, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig8_clinical_trajectory_3panel.pdf"), fig8, width = 16, height = 5, bg = "white")
cat("  Saved Fig8\n")

# ================================================================
# FIGURE 9: Clinical enrichment (copy from earlier script)
# ================================================================
cat("\n--- Figure 9 (copy): Clinical enrichment ---\n")
f9 <- file.path(config$OUTPUT_DIR, "figures_clinical_enrichment_3panel", "Fig_clinical_enrichment_3panel.png")
if (file.exists(f9)) file.copy(f9, file.path(OUT, "Fig9_clinical_enrichment_3panel.png"), overwrite = TRUE)
f9p <- sub(".png$", ".pdf", f9)
if (file.exists(f9p)) file.copy(f9p, file.path(OUT, "Fig9_clinical_enrichment_3panel.pdf"), overwrite = TRUE)
cat("  Copied Fig9\n")

cat("\n=== ALL FIGURES 2-9 COMPLETE ===\n")
cat("Output directory:", OUT, "\n")

})  # end block 17: 12_fig02-05_ST69_analysis.R (full)
}, error = function(e) cat("  BLOCK ERROR (continuing to next block):", conditionMessage(e), "\n"))

