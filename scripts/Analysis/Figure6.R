#!/usr/bin/env Rscript
# ==============================================================================
# Analysis / Figure6.R
# ==============================================================================
# PROVENANCE: verbatim relocation of the original scripts/.../16_fig07_phylogeny_summary.R.
#
# Reproduces submitted Figure 6 (8-panel A-H: ST69 Cluster_3 increasing VFDB genes + clinical enrichment; ST10 VF/ResFinder composition shift, within-stratum trends, and increasing/decreasing genes). The script's own header comment ("Row 1 (ST69): Cluster_3 increasing genes + Clinical enrichment; Row 2 (ST10): VF composition shift + Within-stratum VF trends + AMR decreasing") matches the submitted Figure_06.png panel-for-panel. NOTE: the original filename says "fig07" -- that internal label is WRONG relative to the manuscript's final figure numbering; this is definitively Figure 6, not Figure 7. Renamed here to reflect the confirmed, correct figure number. This script reads config$OUTPUT_DIR/cluster_gene_analysis/ data written by Prerequisites/04b_cluster_gene_analysis.R, which must run first.
#
# ANALYSIS IS UNCHANGED -- this is the original script's body, unmodified
# (aside from this provenance header and, where noted above, a filename
# correction to match the manuscript's actual published figure numbering).
# ==============================================================================

#!/usr/bin/env Rscript
# Combined summary figure:
# Row 1 (ST69): Cluster_3 increasing genes + Clinical enrichment
# Row 2 (ST10): VF composition shift + Within-stratum VF trends + AMR decreasing
suppressPackageStartupMessages(library(tidyverse))
source("config.R")
OUT <- file.path(config$OUTPUT_DIR, "sensitivity_analysis")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

# ======================================================================
# ST69 VFDB data
# ======================================================================
to_bin_vfdb <- function(x) as.integer(sapply(strsplit(as.character(x), ";"),
  function(v) any(suppressWarnings(as.numeric(v) >= 90), na.rm=TRUE)))
master69 <- read_csv(file.path(config$OUTPUT_DIR, "ST69", "vfdb_analysis",
  "04_master_shell_cluster_metadata_VFDB_table.csv"), show_col_types=FALSE, guess_max=10000) %>%
  mutate(genome_id=as.character(genome_id), year=as.integer(year),
    clinical_binary=as.integer(clinical_binary),
    source_niche=as.character(source_niche), country=as.character(country)) %>%
  filter(!is.na(clinical_binary),!is.na(year),!is.na(shell_cluster))
vfdb_raw <- read_tsv(config$st_vfdb_summary("ST69"), show_col_types=FALSE,
  progress=FALSE, col_types=cols(.default="c")) %>%
  filter(str_detect(.data[["#FILE"]], "^ST69/")) %>%
  mutate(genome_id=make_genome_id(str_remove(.data[["#FILE"]], "^ST69/")),
         genome_id=str_remove(genome_id, "_vfdb\\.tsv$"))
genes69 <- setdiff(colnames(vfdb_raw), c("#FILE","NUM_FOUND","genome_id"))
bin69 <- vfdb_raw %>% select(genome_id, all_of(genes69)) %>%
  mutate(across(all_of(genes69), to_bin_vfdb))
meta69 <- master69 %>% select(genome_id, shell_cluster, year, clinical_binary, source_niche, country, total=total_vf)
df69 <- meta69 %>% inner_join(bin69, by="genome_id")
cd69 <- df69 %>% filter(shell_cluster=="Cluster_3")

# G2/G3 capsule classification
kps_cols <- intersect(grep("^(kps|neu)", genes69, value=TRUE), colnames(df69))
G2_markers <- intersect(c("kpsC","kpsS"), kps_cols)
G3_markers <- intersect(c("kpsE","kpsM","kpsT"), kps_cols)
classify_capsule <- function(r) {
  hg2 <- any(r[names(r) %in% G2_markers] >= 1)
  hg3 <- any(r[names(r) %in% G3_markers] >= 1)
  anyk <- any(r >= 1)
  if (hg2) return("G2")
  if (hg3) return("G3")
  if (anyk) return("Unclassified")
  return("No capsule")
}
if (length(kps_cols) > 0) {
  df69$cap_type <- apply(df69[, kps_cols, drop=FALSE], 1, classify_capsule)
} else {
  df69$cap_type <- NA_character_
}
is_g2 <- function(g) g %in% G2_markers
is_g3 <- function(g) g %in% G3_markers
# Gene family classification
classify_gene <- function(g) {
  if (g %in% G2_markers) return("G2 capsule genes (kpsC/S)")
  if (g %in% setdiff(intersect(c("kpsD","kpsE","kpsF","kpsM","kpsU"), genes69), c())) return("G2 operon (kpsD/E/F/M/U)")
  if (g %in% G3_markers) return("G3 capsule (kpsT)")
  if (str_detect(g, "^pap")) return("P fimbriae (pap)")
  if (str_detect(g, "^hly|^sat|^vat|^cnf")) return("Toxins (hly/sat/vat)")
  if (str_detect(g, "^iuc|^iut|^fyu|^irp|^ire|^chu")) return("Siderophore")
  if (str_detect(g, "^fim|^ecp|^afa|^sfa|^hra|^iha")) return("Adhesin (other)")
  if (str_detect(g, "^kps|^neu")) return("Capsule (other)")
  return("Other VF")
}

# ======================================================================
# ST10 VF data
# ======================================================================
master10 <- read_csv(file.path(config$OUTPUT_DIR, "ST10", "virulencefinder_validation",
  "04_master_shell_cluster_metadata_VF_table.csv"), show_col_types=FALSE) %>%
  mutate(genome_id=as.character(genome_id), year=as.integer(year),
    clinical_binary=as.integer(clinical_binary),
    source_niche=as.character(source_niche), country=as.character(country)) %>%
  filter(!is.na(clinical_binary),!is.na(year),!is.na(shell_cluster))
vfb10 <- read.delim(config$VF_BINARY, header=TRUE, sep="\t", check.names=FALSE,
  na.strings="", stringsAsFactors=FALSE) %>%
  filter(st=="ST10") %>% mutate(genome_id=make_genome_id(genome)) %>% select(-st,-genome)
genes_vf10 <- setdiff(colnames(vfb10),"genome_id")
for (g in genes_vf10) vfb10[[g]] <- as.integer(!is.na(vfb10[[g]])&vfb10[[g]]!="."&vfb10[[g]]!="0"&vfb10[[g]]!="")
meta10 <- master10 %>% select(genome_id, shell_cluster, year, clinical_binary, source_niche, country)
df10 <- meta10 %>% inner_join(vfb10 %>% mutate(total_vf=rowSums(across(all_of(genes_vf10)))), by="genome_id")
cd10 <- df10 %>% filter(shell_cluster=="Cluster_3",!is.na(source_niche),!is.na(country))

# ======================================================================
# ST10 ResFinder data
# ======================================================================
resf <- read.delim(file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
  header=TRUE, sep="\t", check.names=FALSE, na.strings="", stringsAsFactors=FALSE) %>%
  filter(st=="ST10") %>% mutate(genome_id=make_genome_id(genome)) %>% select(-st,-genome)
genes_resf <- setdiff(colnames(resf),"genome_id")
for (g in genes_resf) resf[[g]] <- as.integer(!is.na(resf[[g]])&resf[[g]]!="."&resf[[g]]!="0"&resf[[g]]!="")
df_resf <- meta10 %>% inner_join(resf %>% mutate(total_amr=rowSums(across(all_of(genes_resf)))), by="genome_id")
cd6 <- df_resf %>% filter(shell_cluster == "Cluster_6") %>%
  mutate(segment = case_when(
    source_niche == "Human" ~ paste(country, "Human"),
    TRUE ~ source_niche),
    seg_simple = fct_infreq(factor(case_when(
      country=="United Kingdom" & source_niche=="Human" ~ "UK Human",
      source_niche=="Human" ~ "Other Human",
      TRUE ~ source_niche
    ))))
resf_seg_fill <- c("UK Human"="#377EB8", "Other Human"="#4DAF4A", "Environment"="#984EA3",
  "Companion Animal"="#F781BF", "Wild Animal"="#A65628")
# Keep all seg_simple categories as-is
cd6 <- cd6 %>% mutate(seg_top10 = as.character(seg_simple))
seg_amr_v <- cd6 %>% group_by(seg_top10) %>%
  summarise(ma = mean(total_amr), n = n(), .groups = "drop") %>%
  mutate(lab = sprintf("%s (AMR=%.1f, n=%d)", seg_top10, ma, n))
lm6 <- setNames(seg_amr_v$lab, seg_amr_v$seg_top10)

# ======================================================================
# Panel A: ST69 Cluster_3 increasing genes
# ======================================================================
EARLY <- c(2016,2017,2018); LATE <- c(2022,2023,2024,2025)
ea69 <- cd69 %>% filter(year %in% EARLY); la69 <- cd69 %>% filter(year %in% LATE)
n_early <- nrow(ea69); n_late <- nrow(la69)

inc_genes <- bind_rows(lapply(genes69, function(g) {
  prev_e <- mean(ea69[[g]], na.rm=TRUE); prev_l <- mean(la69[[g]], na.rm=TRUE)
  dp <- (prev_l - prev_e) * 100
  mat <- matrix(c(sum(ea69[[g]]), n_early - sum(ea69[[g]]),
                  sum(la69[[g]]), n_late - sum(la69[[g]])), nrow=2)
  ft <- tryCatch(fisher.test(mat), error=function(e) NULL)
  fp <- if (!is.null(ft)) ft$p.value else NA
  tibble(gene=g, prev_early=prev_e, prev_late=prev_l, delta_pp=dp, fisher_p=fp)
})) %>% mutate(p_adj = p.adjust(pmax(fisher_p,0,na.rm=TRUE), method="BH"),
  candidate = delta_pp >= 5 & prev_late >= 0.05 & p_adj < 0.05)

inc_genes <- inc_genes %>% filter(candidate) %>% arrange(desc(delta_pp))

inc_genes <- inc_genes %>% mutate(
  cap_group = map_chr(gene, classify_gene)
)
cap_colors <- c(
  "G2 capsule genes (kpsC/S)" = "#d73027",
  "G2 operon (kpsD/E/F/M/U)" = "#fc8d59",
  "G3 capsule (kpsT)" = "#4575b4",
  "P fimbriae (pap)" = "#1b7837",
  "Toxins (hly/sat/vat)" = "#762a83",
  "Siderophore" = "#9970ab",
  "Adhesin (other)" = "#bababa",
  "Capsule (other)" = "#fee0b6",
  "Other VF" = "grey60"
)

pA <- inc_genes %>% slice_head(n=20) %>%
  mutate(gene = fct_reorder(gene, delta_pp)) %>%
  ggplot(aes(x = delta_pp, y = gene, fill = cap_group)) +
  geom_col() +
  scale_fill_manual(values = cap_colors, name = "") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey50") +
  labs(x = "Delta prevalence (pp), late - early", y = NULL,
       title = "ST69 Cluster_3: increasing VFDB genes") +
  theme_minimal(base_size = 9)

# ======================================================================
# Panel B: ST69 clinical enrichment (all genomes)
# ======================================================================
clin69 <- bind_rows(lapply(genes69, function(g) {
  d <- df69 %>% filter(!is.na(.data[[g]]))
  if (sum(d[[g]], na.rm=TRUE) < 5) return(tibble(gene=g, or=NA, p=NA))
  ft <- tryCatch(fisher.test(table(d$clinical_binary, d[[g]])), error=function(e) NULL)
  if (is.null(ft)) return(tibble(gene=g, or=NA, p=NA))
  tibble(gene=g, or=unname(ft$estimate), p=ft$p.value)
})) %>% mutate(p_adj = p.adjust(pmax(p,0,na.rm=TRUE), method="BH"),
  enriched = !is.na(p_adj) & p_adj < 0.05 & or > 1)

clin69 <- clin69 %>% mutate(
  cap_group = map_chr(gene, classify_gene)
)
top_clin <- clin69 %>% filter(!is.na(or), is.finite(or)) %>%
  arrange(p_adj) %>% slice_head(n=20) %>% mutate(gene = fct_reorder(gene, or))

pB <- top_clin %>%
  ggplot(aes(x = or, y = gene, color = cap_group, shape = enriched)) +
  geom_point(size = 2.5) + geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(values = cap_colors, name = "") +
  scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1), guide = "none") +
  labs(x = "OR (clinical vs non-clinical)", y = NULL,
       title = "ST69: clinical enrichment (all genomes)") +
  theme_minimal(base_size = 9)

# ======================================================================
# Panel C: ST10 composition shift
# ======================================================================
cd10 <- cd10 %>% mutate(
  segment = factor(case_when(
    country=="United States" & source_niche=="Human" ~ "US Human",
    country=="United Kingdom" & source_niche=="Human" ~ "UK Human",
    source_niche=="Human" ~ "Other Human",
    TRUE ~ source_niche
  ), levels = c("US Human","UK Human","Other Human","Environment","Livestock","Poultry",
    "Wild Animal","Companion Animal","Food")))

seg_vf <- cd10 %>% group_by(segment) %>%
  summarise(mvf = mean(total_vf), n = n(), .groups = "drop") %>%
  mutate(lab = sprintf("%s (VF=%.0f, n=%d)", segment, mvf, n))
lab_map <- setNames(seg_vf$lab, seg_vf$segment)

seg_colors <- c("US Human"="#E41A1C","UK Human"="#377EB8","Other Human"="#4DAF4A",
  "Environment"="#984EA3","Livestock"="#FF7F00","Poultry"="#FFFF33",
  "Wild Animal"="#A65628","Companion Animal"="#F781BF","Food"="#999999")

pC <- cd10 %>% count(year, segment) %>% group_by(year) %>%
  mutate(pct = n/sum(n)*100) %>% ungroup() %>%
  ggplot(aes(x = year, y = pct, fill = segment)) +
  geom_col() +
  scale_fill_manual(values = seg_colors, labels = lab_map, name = NULL) +
  labs(x = NULL, y = "Proportion (%)",
       title = "ST10 VF Cluster_3: composition drives VF increase") +
  theme_minimal(base_size = 8) +
  theme(legend.text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

# ======================================================================
# Panel D: ST10 within-stratum VF trends
# ======================================================================
big_seg <- seg_vf %>% filter(n >= 10) %>% pull(segment)
pd <- cd10 %>% filter(segment %in% big_seg) %>%
  group_by(year, segment) %>%
  summarise(mvf = mean(total_vf), n = n(), .groups = "drop") %>%
  filter(n >= 3)

pD <- pd %>%
  ggplot(aes(x = year, y = mvf, color = segment)) +
  geom_point(size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, alpha = 0.4) +
  scale_color_manual(values = seg_colors, name = NULL) +
  labs(x = NULL, y = "Mean VF burden",
       title = "ST10: within-stratum VF trends (flat)") +
  theme_minimal(base_size = 8) +
  theme(legend.text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

# ======================================================================
# Panel E: ST10 ResFinder Cluster_6 decreasing genes (from pre-computed)
# ======================================================================
dec6_raw <- read_csv(file.path(config$OUTPUT_DIR, "cluster_gene_analysis",
  "st10_resfinder", "Cluster_6_gene_analysis.csv"), show_col_types = FALSE)

dec6_genes <- dec6_raw %>% filter(candidate_decreasing) %>% arrange(delta_pp)

pE <- dec6_genes %>%
  mutate(gene = fct_reorder(gene, delta_pp)) %>%
  ggplot(aes(x = delta_pp, y = gene)) +
  geom_col(fill = "#377EB8") +
  geom_vline(xintercept = -5, linetype = "dashed", color = "grey50") +
  labs(x = "Delta prevalence (pp)", y = NULL,
       title = "ST10 ResFinder Cluster_6: decreasing AMR genes",
       subtitle = "Early high-AMR strata replaced by lower-AMR Human genomes") +
  theme_minimal(base_size = 9)

# ======================================================================
# Panel F: ST10 ResFinder Cluster_6 within-stratum AMR trends (simplified)
# ======================================================================
seg_amr_sum <- cd6 %>% group_by(seg_top10) %>%
  summarise(ma = mean(total_amr), n = n(), .groups = "drop") %>% filter(n >= 5)
big_seg6 <- seg_amr_sum %>% pull(seg_top10)

pF <- cd6 %>% filter(seg_top10 %in% big_seg6) %>%
  group_by(year, seg_top10) %>%
  summarise(m = mean(total_amr), n = n(), .groups = "drop") %>% filter(n >= 3) %>%
  ggplot(aes(x = year, y = m, color = seg_top10)) +
  geom_point(size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, alpha = 0.4) +
  scale_color_manual(values = resf_seg_fill, name = NULL) +
  labs(x = NULL, y = "Mean AMR burden",
       title = "ST10 ResF Cluster_6: within-stratum AMR trends (flat)") +
  theme_minimal(base_size = 8) +
  theme(legend.text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

# ======================================================================
# Panel G: ST10 ResFinder Cluster_6 AMR composition (top 10 legend)
# ======================================================================
pG <- cd6 %>% count(year, seg_top10) %>% group_by(year) %>%
  mutate(pct = n/sum(n)*100) %>% ungroup() %>%
  ggplot(aes(x = year, y = pct, fill = seg_top10)) +
  geom_col() +
  scale_fill_manual(values = resf_seg_fill, labels = lm6, name = NULL) +
  labs(x = NULL, y = "Proportion (%)",
       title = "ST10 ResF Cluster_6: composition shift") +
  theme_minimal(base_size = 8) +
  theme(legend.text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

# ======================================================================
# Panel H: ST10 VF Cluster_3 increasing genes
# ======================================================================
ea10 <- cd10 %>% filter(year %in% EARLY); la10 <- cd10 %>% filter(year %in% LATE)
ne10 <- nrow(ea10); nl10 <- nrow(la10)

inc10 <- bind_rows(lapply(genes_vf10, function(g) {
  pe <- mean(ea10[[g]], na.rm=TRUE); pl <- mean(la10[[g]], na.rm=TRUE)
  dp <- (pl-pe)*100
  ft <- tryCatch(fisher.test(matrix(c(sum(ea10[[g]]),ne10-sum(ea10[[g]]),sum(la10[[g]]),nl10-sum(la10[[g]])),nrow=2)), error=function(e) NULL)
  tibble(gene=g, prev_e=pe, prev_l=pl, dp=dp, p=if(!is.null(ft))ft$p.value else NA)
})) %>% mutate(p_a=p.adjust(pmax(p,0,na.rm=TRUE),method="BH")) %>%
  filter(dp>=5&prev_l>=0.05&p_a<0.05) %>% arrange(desc(dp))

pH <- inc10 %>% slice_head(n=20) %>% mutate(gene=fct_reorder(gene,dp)) %>%
  ggplot(aes(x=dp,y=gene)) + geom_col(fill="#4DAF4A") +
  geom_vline(xintercept=5,linetype="dashed",color="grey50") +
  labs(x="Delta prevalence (pp)",y=NULL,title="ST10 VF Cluster_3: increasing VF genes") +
  theme_minimal(base_size=9)

# ======================================================================
# Assemble: 4 rows × 2 columns
# Row 1: A (ST69 inc) | B (ST69 clin)
# Row 2: H (ST10 inc) | C (ST10 comp)
# Row 3: D (ST10 stratum) | G (ResF comp)
# Row 4: E (ResF dec) | F (ResF stratum)
# ======================================================================
library(patchwork)

fig <- (pA | pB) / (pC | pG) / (pD | pF) / (pH | pE) +
  plot_layout(heights = c(1, 1, 1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(OUT, "figure_combined_summary.png"), fig, width = 8, height = 10, dpi = 300)
ggsave(file.path(OUT, "figure_combined_summary.pdf"), fig, width = 8, height = 10)

cat("Saved:", file.path(OUT, "figure_combined_summary.png"), "\n")
cat("Saved:", file.path(OUT, "figure_combined_summary.pdf"), "\n")
