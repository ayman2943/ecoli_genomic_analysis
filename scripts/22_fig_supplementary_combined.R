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
