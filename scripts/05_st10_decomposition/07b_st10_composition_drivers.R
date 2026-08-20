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
