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
