#!/usr/bin/env Rscript
#
# cluster_gene_analysis.R  —  Per-gene driver analysis matching
# publication method: early-vs-late prevalence (Δ≥5pp, late≥5%,
# Fisher or logit BH-adj p<0.05) + clinical enrichment on all genomes.
#
suppressPackageStartupMessages({
  library(tidyverse); library(broom); library(patchwork)
})
ts <- function(m) cat("[", Sys.time(), "] ", m, "\n", sep = "")
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "cluster_gene_analysis")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

EARLY <- c(2016, 2018); LATE <- c(2022, 2025)
# Publication 2016,2017,2018 and 2022,2023,2024,2025
EARLY <- c(2016, 2017, 2018)
LATE  <- c(2022, 2023, 2024, 2025)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

# ====================================================================
# Core analysis function
# ====================================================================
analyze_dataset <- function(label, metadata_df, bin_df, clusters) {
  ts(paste("\n=========== ", label, " ===========", sep = ""))

  # Join metadata + binary matrix
  genes <- setdiff(colnames(bin_df), c("genome_id", "st", "genome", "total"))
  df <- metadata_df %>% inner_join(bin_df, by = "genome_id") %>%
    filter(!is.na(clinical_binary), !is.na(year), !is.na(shell_cluster))
  # Compute total burden if not provided
  if (!"total" %in% colnames(df)) {
    df <- df %>% mutate(total = rowSums(across(all_of(genes))))
  }
  cat("  Joined data:", nrow(df), "genomes x", length(genes), "genes\n")

  track_dir <- file.path(OUT, gsub(" ", "_", tolower(label)))
  dir.create(track_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Clinical enrichment across ALL genomes ----
  ts("  Clinical enrichment (all genomes, NOT per cluster)...")
  enrich_all <- bind_rows(lapply(genes, function(g) {
    d <- df %>% filter(!is.na(.data[[g]]))
    if (nrow(d) < 20 || length(unique(d[[g]])) < 2 ||
        sum(d[[g]] == 1, na.rm = TRUE) < 5) {
      return(tibble(gene = g, clin_or = NA_real_, clin_p = NA_real_))
    }
    ft <- tryCatch(fisher.test(table(d$clinical_binary, d[[g]])),
                   error = function(e) NULL)
    if (is.null(ft)) return(tibble(gene = g, clin_or = NA_real_, clin_p = NA_real_))
    tibble(gene = g, clin_or = unname(ft$estimate), clin_p = ft$p.value)
  })) %>% mutate(clin_p_adj = p.adjust(pmax(clin_p, 0, na.rm = TRUE), method = "BH"),
                  clin_enriched = !is.na(clin_p_adj) & clin_p_adj < 0.05 & clin_or > 1)

  write_csv(enrich_all, file.path(track_dir, "clinical_enrichment_all_genomes.csv"))

  for (cl in clusters) {
    ts(paste("  Cluster:", cl))
    cd <- df %>% filter(shell_cluster == cl)
    cl_safe <- gsub("[^A-Za-z0-9]", "_", cl)
    n_early <- sum(cd$year %in% EARLY)
    n_late  <- sum(cd$year %in% LATE)

    # ---- Early vs late prevalence (publication method) ----
    ea <- cd %>% filter(year %in% EARLY)
    la <- cd %>% filter(year %in% LATE)

    tests <- bind_rows(lapply(genes, function(g) {
      e_tab <- table(factor(ea[[g]], levels = 0:1))
      l_tab <- table(factor(la[[g]], levels = 0:1))
      prev_e <- if (n_early > 0) e_tab["1"] / n_early else NA
      prev_l <- if (n_late  > 0) l_tab["1"] / n_late  else NA
      delta_pp <- (prev_l - prev_e) * 100

      # Fisher
      mat <- matrix(c(e_tab["1"], e_tab["0"], l_tab["1"], l_tab["0"]), nrow = 2)
      fisher <- tryCatch(fisher.test(mat), error = function(e) NULL)
      fisher_p <- if (!is.null(fisher)) fisher$p.value else NA
      fisher_or <- if (!is.null(fisher)) unname(fisher$estimate) else NA

      # Logistic gene ~ year (all years in cluster)
      d <- cd %>% filter(!is.na(.data[[g]]))
      if (nrow(d) >= 20 && n_distinct(d$year) >= 3 &&
          sum(d[[g]], na.rm = TRUE) >= 5) {
        logit <- tryCatch(glm(as.formula(paste0(g, " ~ year")),
                             data = d, family = binomial()),
                          error = function(e) NULL)
        if (!is.null(logit)) {
          s <- tryCatch(summary(logit)$coefficients, error = function(e) NULL)
          if (!is.null(s) && "year" %in% rownames(s)) {
            logit_or <- exp(s["year", "Estimate"])
            logit_p  <- s["year", "Pr(>|z|)"]
          } else { logit_or <- NA; logit_p <- NA }
        } else { logit_or <- NA; logit_p <- NA }
      } else { logit_or <- NA; logit_p <- NA }

      tibble(gene = g, prev_early = prev_e, prev_late = prev_l,
             delta_pp = delta_pp, n_early = n_early, n_late = n_late,
             fisher_or = fisher_or, fisher_p = fisher_p,
             logit_or = logit_or, logit_p = logit_p)
    })) %>%
      mutate(fisher_p_adj = p.adjust(pmax(fisher_p, 0, na.rm = TRUE), method = "BH"),
             logit_p_adj  = p.adjust(pmax(logit_p, 0, na.rm = TRUE),  method = "BH")) %>%
      # Publication criteria: Δ≥5pp, late≥5%, Fisher adj<0.05 OR logit adj<0.05
      mutate(candidate_increasing = !is.na(delta_pp) & delta_pp >= 5 &
               !is.na(prev_late) & prev_late >= 0.05 &
               ((!is.na(fisher_p_adj) & fisher_p_adj < 0.05) |
                (!is.na(logit_p_adj)  & logit_p_adj  < 0.05)),
             candidate_decreasing = !is.na(delta_pp) & delta_pp <= -5 &
               !is.na(prev_early) & prev_early >= 0.05 &
               ((!is.na(fisher_p_adj) & fisher_p_adj < 0.05) |
                (!is.na(logit_p_adj)  & logit_p_adj  < 0.05)))

    # ---- Merge with clinical enrichment ----
    merged <- tests %>%
      left_join(enrich_all, by = "gene") %>%
      arrange(desc(candidate_increasing), desc(delta_pp))

    write_csv(merged, file.path(track_dir, paste0(cl_safe, "_gene_analysis.csv")))

    inc <- sum(merged$candidate_increasing, na.rm = TRUE)
    dec <- sum(merged$candidate_decreasing, na.rm = TRUE)
    enr <- sum(merged$candidate_increasing & merged$clin_enriched, na.rm = TRUE)
    ts(sprintf("  %s: %d increasing, %d decreasing, %d inc+clinically enriched",
               cl, inc, dec, enr))

    # ---- Sampling bias models ----
    cd_model <- cd %>% filter(!is.na(source_niche), !is.na(country))
    adj_model <- tryCatch({
      m <- lm(total ~ year + source_niche + country, data = cd_model)
      tidy(m) %>% filter(term == "year") %>%
        transmute(model = "total ~ year + niche + country",
                  estimate = estimate, std.error = std.error, p.value = p.value)
    }, error = function(e) tibble(model = "failed", estimate = NA, std.error = NA, p.value = NA))

    df_model <- df %>% filter(!is.na(source_niche), !is.na(country)) %>%
      mutate(is_target = as.integer(shell_cluster == cl))
    interaction_model <- tryCatch({
      m <- glm(clinical_binary ~ year * is_target + country,
               data = df_model, family = binomial())
      tidy(m) %>% filter(str_detect(term, ":")) %>%
        transmute(model = "clinical ~ year * target + country",
                  term = term, estimate = estimate, std.error = std.error, p.value = p.value)
    }, error = function(e) tibble(model = "failed", term = NA, estimate = NA, std.error = NA, p.value = NA))

    write_csv(adj_model, file.path(track_dir, paste0(cl_safe, "_adjusted_model.csv")))
    write_csv(interaction_model, file.path(track_dir, paste0(cl_safe, "_interaction_model.csv")))

    ts(sprintf("  Adjusted model p=%.4f | Interaction p=%.4f",
               adj_model$p.value[1] %||% NA, interaction_model$p.value[1] %||% NA))

    # ---- Figure: top candidate increasing + clinical enrichment ----
    candidates <- merged %>% filter(candidate_increasing) %>% arrange(desc(delta_pp))
    if (nrow(candidates) > 0) {
      p1 <- candidates %>% slice_head(n = 25) %>%
        mutate(gene = fct_reorder(gene, delta_pp)) %>%
        ggplot(aes(x = delta_pp, y = gene, fill = clin_enriched)) +
        geom_col() +
        geom_vline(xintercept = 5, linetype = "dashed", color = "grey50") +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "steelblue"),
          labels = c("TRUE" = "Clinically enriched", "FALSE" = "Not enriched"),
          name = NULL) +
        labs(x = bquote(Delta~prevalence~(late-early)~"(pp)"), y = NULL,
             title = paste(label, "-", cl, "candidate increasing genes"),
             subtitle = sprintf("%d increasing, %d also clinically enriched", inc, enr)) +
        theme_minimal(base_size = 9) + theme(legend.position = "bottom")
    } else {
      p1 <- ggplot() + ggtitle(paste(label, cl, "- No candidate increasing genes")) +
        theme_void()
    }

    # Top clinically enriched across all genomes
    top_enr <- enrich_all %>% filter(!is.na(clin_or), is.finite(clin_or)) %>%
      arrange(clin_p_adj) %>% slice_head(n = 20) %>%
      mutate(gene = fct_reorder(gene, clin_or))
    if (nrow(top_enr) >= 3) {
      p2 <- top_enr %>%
        ggplot(aes(x = clin_or, y = gene, color = clin_enriched)) +
        geom_point(size = 2.5) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        scale_x_log10() +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey"), guide = "none") +
        labs(x = "OR (clinical vs non-clinical)", y = NULL,
             title = paste(label, "- Clinical enrichment (all genomes)")) +
        theme_minimal(base_size = 9)
    } else {
      p2 <- ggplot() + ggtitle(paste(label, "- Too few genes")) + theme_void()
    }

    combined <- p1 / p2 + plot_annotation(title = paste(label, "-", cl))
    ggsave(file.path(track_dir, paste0(cl_safe, "_gene_analysis.png")),
           combined, width = 10, height = 14)
  }
}

# ====================================================================
# 1. ST69 VFDB — use master for metadata + VFDB summary TSV for genes
# ====================================================================
ts("\n=========== 1. ST69 VFDB ===========")
master69 <- read_csv(
  file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
            "04_master_shell_cluster_metadata_VFDB_table.csv"),
  show_col_types = FALSE, guess_max = 10000) %>%
  mutate(genome_id = as.character(genome_id),
         year = as.integer(year),
         clinical_binary = as.integer(clinical_binary),
         total = total_vf) %>%
  filter(!is.na(clinical_binary), !is.na(year), !is.na(shell_cluster))

# Load raw VFDB summary and binarize using publication method
to_bin_vfdb <- function(x) as.integer(sapply(strsplit(as.character(x), ";"),
  function(v) any(suppressWarnings(as.numeric(v) >= 90), na.rm = TRUE)))

vfdb_raw <- read_tsv(config$st_vfdb_summary("ST69"),
  show_col_types = FALSE, progress = FALSE, col_types = cols(.default = "c")) %>%
  filter(str_detect(.data[["#FILE"]], "^ST69/")) %>%
  mutate(genome_id = make_genome_id(str_remove(.data[["#FILE"]], "^ST69/")),
         genome_id = str_remove(genome_id, "_vfdb\\.tsv$"))
vf_genes <- setdiff(colnames(vfdb_raw), c("#FILE", "NUM_FOUND", "genome_id"))
vfdb_bin_raw <- vfdb_raw %>% select(genome_id, all_of(vf_genes)) %>%
  mutate(across(all_of(vf_genes), to_bin_vfdb))

cat("ST69 VFDB raw:", nrow(vfdb_bin_raw), "genomes x", length(vf_genes), "genes\n")

meta69 <- master69 %>% select(genome_id, shell_cluster, year, clinical_binary,
                              source_niche, country, total)

analyze_dataset("ST69 VFDB", meta69, vfdb_bin_raw, c("Cluster_3"))

# ====================================================================
# 2. ST10 VF — load master (metadata) + VF binary (genes)
# ====================================================================
ts("\n=========== 2. ST10 VF ===========")
master10 <- read_csv(
  file.path(config$OUTPUT_DIR, "ST10", "virulencefinder_validation",
            "04_master_shell_cluster_metadata_VF_table.csv"),
  show_col_types = FALSE) %>%
  mutate(genome_id = as.character(genome_id),
         year = as.integer(year),
         clinical_binary = as.integer(clinical_binary)) %>%
  filter(!is.na(clinical_binary), !is.na(year), !is.na(shell_cluster))

vfbinary <- read.delim(config$VF_BINARY, header = TRUE, sep = "\t",
  check.names = FALSE, na.strings = "", stringsAsFactors = FALSE) %>%
  filter(st == "ST10") %>%
  mutate(genome_id = make_genome_id(genome)) %>%
  select(-st, -genome)

genes_vf <- setdiff(colnames(vfbinary), "genome_id")
# Convert to binary (handling various formats)
for (g in genes_vf) {
  vfbinary[[g]] <- as.integer(!is.na(vfbinary[[g]]) & vfbinary[[g]] != "." &
    vfbinary[[g]] != "0" & vfbinary[[g]] != "")
}

cat("ST10 VF:", nrow(master10), "genomes in master,", ncol(vfbinary) - 1, "genes\n")

analyze_dataset("ST10 VF", master10, vfbinary, c("Cluster_3"))

# ====================================================================
# 3. ST10 ResFinder — load master (metadata) + ResFinder binary (genes)
# ====================================================================
ts("\n=========== 3. ST10 ResFinder ===========")
resf <- read.delim(
  file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
  header = TRUE, sep = "\t", check.names = FALSE,
  na.strings = "", stringsAsFactors = FALSE) %>%
  filter(st == "ST10") %>%
  mutate(genome_id = make_genome_id(genome)) %>%
  select(-st, -genome)

genes_resf <- setdiff(colnames(resf), "genome_id")
for (g in genes_resf) {
  resf[[g]] <- as.integer(!is.na(resf[[g]]) & resf[[g]] != "." &
    resf[[g]] != "0" & resf[[g]] != "")
}

cat("ST10 ResFinder:", nrow(resf), "genomes x", length(genes_resf), "genes\n")

analyze_dataset("ST10 ResFinder", master10, resf, c("Cluster_1", "Cluster_6"))

# ====================================================================
ts("\n=========== DONE ===========")
cat("Output:", normalizePath(OUT), "\n")
