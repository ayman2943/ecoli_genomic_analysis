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
