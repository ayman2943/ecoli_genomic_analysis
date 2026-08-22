#!/usr/bin/env Rscript
# ==============================================================================
# Analysis / Figure8.R
# ==============================================================================
# PROVENANCE: verbatim relocation of the original scripts/.../10d_plasmid_context.R.
#
# Reproduces submitted Figure 8 (Panel A "Genomic location (all ST69)" chromosomal-vs-plasmid stacked bars + Panel B "RGP co-localization" heatmap). This script independently produces BOTH panels itself from the plasmid VFDB summary and RGP module data -- confirmed via a full read of its code, matching the submitted Figure_08.png exactly. It also produces a Panel C ("C3 increase early->late") as part of its own 3-panel combined output; that panel does not appear in the submitted Figure_08.png but is kept here since it is unmodified original analysis code, not something added or removed. This script's sections 1-3 (data loading / interest-gene classification) are ALSO reproduced, in isolation, as Virulence/Plasmid_VFDB.R -- both are independent, self-contained entry points into the same source data, matching how this script was always structured (it re-derives everything itself).
#
# ANALYSIS IS UNCHANGED -- this is the original script's body, unmodified
# (aside from this provenance header and, where noted above, a filename
# correction to match the manuscript's actual published figure numbering).
# ==============================================================================

#!/usr/bin/env Rscript
# Combined MGE context: plasmid vs chromosomal + RGP neighbourhood
# Generates Figure C?? for reviewer response

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(patchwork); library(viridis)
})
source("config.R")

OUT <- file.path(config$OUTPUT_DIR, config$TARGET_ST, "reviewer_rgp_context")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
EARLY <- c(2016,2017,2018); LATE <- c(2022,2023,2024,2025)

theme_pub <- function(base_size=10) {
  theme_classic(base_size=base_size) +
    theme(plot.title=element_text(face="bold",size=rel(1.08),hjust=0.5),
      axis.title=element_text(face="bold",size=rel(0.9)),
      axis.text=element_text(color="grey20",size=rel(0.78)),
      axis.line=element_line(color="grey30",linewidth=0.4),
      axis.ticks=element_line(color="grey30",linewidth=0.3),
      panel.grid.major.y=element_line(color="grey92",linewidth=0.3),
      strip.background=element_rect(fill="grey95",color=NA),
      strip.text=element_text(face="bold",size=rel(0.82)),
      legend.title=element_text(face="bold",size=rel(0.82)),
      legend.text=element_text(size=rel(0.75)),
      plot.background=element_rect(fill="white",color=NA))
}

################################################################################
# 1. Helper: clean genome IDs
################################################################################
clean_gid <- function(raw) {
  gid <- basename(raw)
  gid <- sub("^Escherichia_coli_","",gid)
  gid <- sub("_combined_plasmids_vfdb\\.tsv$","",gid)
  gid <- sub("_combined_plasmids_card\\.tsv$","",gid)
  gid <- sub("_vfdb\\.tsv$","",gid)
  gid <- sub("_card\\.tsv$","",gid)
  gid <- sub("\\.tsv$","",gid)
  gid
}

################################################################################
# 2. Load data
################################################################################
cat("Loading data...\n")

# Master table (has cluster + VF binary)
master <- fread(
  file.path(config$st_out_vfdb(config$TARGET_ST), "04_master_shell_cluster_metadata_VFDB_table.csv"),
  sep=",", header=TRUE, data.table=FALSE) %>%
  mutate(gid_clean = clean_gid(.data[["#FILE"]]), year=as.integer(year))

# Interest genes
interest <- c("iucA","iucB","iucC","iucD","iutA","sat",
  "papA","papB","papC","papD","papE","papF","papG","papH","papI","papJ","papK","papX",
  "kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT","kpsU")
gene_sys <- c(
  iucA="iuc",iucB="iuc",iucC="iuc",iucD="iuc",iutA="iuc",
  sat="sat",
  papA="pap",papB="pap",papC="pap",papD="pap",papE="pap",papF="pap",
  papG="pap",papH="pap",papI="pap",papJ="pap",papK="pap",papX="pap",
  kpsC="kps",kpsD="kps",kpsE="kps",kpsF="kps",kpsM="kps",kpsS="kps",kpsT="kps",kpsU="kps")

# Parse whole-genome VF binary from master table
to_bin <- function(x) as.integer(!is.na(x) & x != ".")
master_bin <- master %>%
  mutate(across(all_of(intersect(interest, colnames(master))), to_bin, .names="{.col}_wg"))

################################################################################
# 3. Load plasmid VFDB summary
################################################################################
cat("Loading plasmid data...\n")
pl_file <- file.path(config$BASE_DIR,"..","Plasmid","plasmid_vfdb_summary",paste0(config$TARGET_ST, "_plasmid_summary.tsv"))
if (!file.exists(pl_file)) {
  stop("Plasmid VFDB summary not found at ", pl_file,
       " (requires the external Plasmid dataset). Skipping plasmid-context analysis.")
}
pl_vf <- fread(pl_file,
  sep="\t", header=TRUE, data.table=FALSE, na.strings=c("",".","-","NA"),
  colClasses="character", check.names=FALSE) %>%
  mutate(gid_clean = clean_gid(.data[["#FILE"]]))

# Parse plasmid binary
pl_genes <- intersect(interest, setdiff(colnames(pl_vf), c("#FILE","NUM_FOUND","gid_clean")))
pl_bin <- pl_vf %>%
  select(gid_clean, all_of(pl_genes)) %>%
  mutate(across(all_of(pl_genes), to_bin, .names="{.col}_pl"))

################################################################################
# 4. Find top 10 increasing VF genes in C3
################################################################################
cat("\n=== Finding top 10 increasing VF genes in C3 ===\n")

# Use the original master (pre-binary conversion) for all VF columns
c3_master <- master %>%
  filter(shell_cluster == "Cluster_3", !is.na(year)) %>%
  mutate(period = ifelse(year %in% EARLY, "early", ifelse(year %in% LATE, "late", NA))) %>%
  filter(!is.na(period))

# Find VF columns (anything not metadata)
all_cols <- colnames(c3_master)
meta_cols <- c("#FILE","FILE","GFF","species","st","year","country","continent",
  "source","shell_cluster","gid_clean","gid")
vf_cols <- setdiff(all_cols, meta_cols)
vf_cols <- vf_cols[!vf_cols %in% c("clinical_binary","genome_id","raw_cluster","total_vf","AAA92657","C_RS05810")]
cat("Number of VF columns analyzed:", length(vf_cols), "\n")

early_total <- sum(c3_master$period == "early")
late_total <- sum(c3_master$period == "late")
cat("Early C3 n:", early_total, "Late C3 n:", late_total, "\n")

# Convert to binary and compute deltas
to_bin_vf <- function(x) as.integer(!is.na(x) & x != ".")
top_deltas <- data.frame(gene=character(), early_pct=numeric(), late_pct=numeric(),
  delta=numeric(), stringsAsFactors=FALSE)
for (gene in vf_cols) {
  vals <- to_bin_vf(c3_master[[gene]])
  early_n <- sum(vals[c3_master$period == "early"] == 1, na.rm=TRUE)
  late_n <- sum(vals[c3_master$period == "late"] == 1, na.rm=TRUE)
  top_deltas <- rbind(top_deltas, data.frame(
    gene=gene,
    early_pct=round(early_n/early_total*100, 1),
    late_pct=round(late_n/late_total*100, 1),
    delta=round((late_n/late_total - early_n/early_total)*100, 1),
    stringsAsFactors=FALSE))
}
top_deltas <- top_deltas[order(top_deltas$delta, decreasing=TRUE),]
top10 <- head(top_deltas, 10)
cat("Top 10:\n")
print(top10, row.names=FALSE)

################################################################################
# 5. Per-gene plasmid prevalence for top 10
################################################################################
cat("\n=== Per-gene plasmid prevalence ===\n")

# Create binary WG columns for top 10 genes
merged <- master %>%
  select(gid_clean, shell_cluster, year, all_of(top10$gene)) %>%
  mutate(across(all_of(top10$gene), to_bin_vf, .names="{.col}_wg")) %>%
  select(gid_clean, shell_cluster, year, ends_with("_wg")) %>%
  left_join(pl_bin, by="gid_clean", suffix=c("",".pl_only"))

# For top-10 genes not in plasmid data, add 0 column
gene_loc <- map_dfr(top10$gene, function(g) {
  wg <- paste0(g, "_wg")
  pl <- paste0(g, "_pl")
  if (!pl %in% colnames(merged)) merged[[pl]] <<- 0L
  
  merged %>%
    summarise(
      n_total = n(),
      n_wg = sum(.data[[wg]] == 1, na.rm=TRUE),
      n_pl = sum(.data[[pl]] == 1, na.rm=TRUE)
    ) %>%
    mutate(
      gene = g,
      pct_wg = round(n_wg/n_total*100, 1),
      pct_pl = round(n_pl/n_total*100, 1),
      pct_chr = round((n_wg - n_pl)/n_total*100, 1),
      pct_abs = round(100 - pct_wg, 1)
    )
})

gene_loc$system <- c("kps","kps","kps","kps","kps","kps","pap","sat","pap","pap")
cat("Gene-level localization (all ST69):\n")
print(gene_loc %>% select(gene, pct_wg, pct_chr, pct_pl, pct_abs), row.names=FALSE)

################################################################################
# 6. Temporal deltas per gene (Panel C data)
################################################################################
cat("\n=== Temporal deltas ===\n")

delta_data <- top10 %>% mutate(
  gene = factor(gene, levels=rev(top10$gene)),
  early_pct = round(early_pct, 1),
  late_pct = round(late_pct, 1),
  delta = round(late_pct - early_pct, 1))
print(delta_data, row.names=FALSE)

################################################################################
# 7. RGP co-occurrence for top genes
################################################################################
cat("\n=== RGP co-occurrence ===\n")

rgp_vf <- fread(file.path(config$PANGENOME_DIR,"vf_module_rgp_results","tables",
  "VF_families_in_RGPs_with_modules.csv"),
  sep=",", header=TRUE, data.table=FALSE, nrows=2e6)

rgp_top <- rgp_vf %>% filter(gene_name %in% top10$gene) %>%
  select(rgp_id, gene_name) %>% distinct()

rgp_top_genes <- unique(rgp_top$gene_name)
cat("Top-10 genes found in RGPs:", paste(rgp_top_genes, collapse=", "), "\n")
cat("Not found in RGPs:", paste(setdiff(top10$gene, rgp_top_genes), collapse=", "), "\n")

# Gene-level co-occurrence matrix
if (length(rgp_top_genes) >= 2) {
  gene_pairs <- map_dfr(combn(rgp_top_genes, 2, simplify=FALSE), ~{
    s1 <- .x[1]; s2 <- .x[2]
    rgps_s1 <- unique(rgp_top$rgp_id[rgp_top$gene_name == s1])
    rgps_s2 <- unique(rgp_top$rgp_id[rgp_top$gene_name == s2])
    n_both <- length(intersect(rgps_s1, rgps_s2))
    tibble(gene1=s1, gene2=s2, n_rgp_both=n_both, n_rgp_g1=length(rgps_s1), n_rgp_g2=length(rgps_s2))
  }) %>%
    mutate(pct_of_g1 = round(n_rgp_both/n_rgp_g1*100, 1))
  
  cat("RGP co-occurrence:\n")
  print(gene_pairs, row.names=FALSE)
}

################################################################################
# 8. Combined figure (top 10 genes)
################################################################################
cat("\n=== Generating combined figure ===\n")

# --- Panel A: Stacked bar per gene (chromosomal vs plasmid) ---
pA_data <- gene_loc %>%
  select(gene, chromosomal=pct_chr, plasmid=pct_pl, absent=pct_abs) %>%
  pivot_longer(-gene, names_to="location", values_to="pct") %>%
  mutate(location = factor(location, levels=c("absent","chromosomal","plasmid")),
    gene = factor(gene, levels=rev(top10$gene)))

pA <- ggplot(pA_data, aes(x=gene, y=pct, fill=location)) +
  geom_col(position="stack", width=0.7, color="grey30", linewidth=0.3) +
  scale_fill_manual(values=c("chromosomal"="#2c7bb6","plasmid"="#d7191c","absent"="grey90"),
    labels=c("chromosomal"="Chromosomal only","plasmid"="Plasmid-associated","absent"="Absent")) +
  labs(title="A. Genomic location (all ST69)", y="% genomes", x="", fill="") +
  scale_y_continuous(expand=c(0,0), limits=c(0,105)) +
  coord_flip() +
  theme_pub()

# --- Panel B: RGP co-occurrence heatmap ---
if (length(rgp_top_genes) >= 2) {
  heat_data <- gene_pairs %>%
    bind_rows(gene_pairs %>% rename(gene1=gene2, gene2=gene1) %>%
      mutate(pct_of_g1 = round(n_rgp_both/n_rgp_g2*100, 1))) %>%
    bind_rows(tibble(gene1=rgp_top_genes, gene2=rgp_top_genes,
      n_rgp_both=NA, n_rgp_g1=NA, n_rgp_g2=NA, pct_of_g1=100))
  
  pB <- ggplot(heat_data, aes(x=gene1, y=gene2, fill=pct_of_g1)) +
    geom_tile(color="white", linewidth=0.5) +
    geom_text(aes(label=ifelse(is.na(pct_of_g1), "", sprintf("%.0f%%", pct_of_g1))), size=3) +
    scale_fill_gradient(low="white", high="#2c7bb6", na.value="grey90", name="% RGPs\nco-occurring") +
    labs(title="B. RGP co-localization", x="", y="") +
    theme_minimal(base_size=9) +
    theme(panel.grid=element_blank(),
      axis.text.x=element_text(angle=45, hjust=1, face="bold", size=rel(0.8)),
      axis.text.y=element_text(face="bold", size=rel(0.8)))
} else {
  pB <- ggplot() + annotate("text", x=0.5, y=0.5, label="No RGP co-occurrence data", size=5) +
    theme_void() + labs(title="B. RGP co-localization")
}

# --- Panel C: Delta bar chart ---
pC <- ggplot(delta_data, aes(x=delta, y=gene, fill=delta>0)) +
  geom_col(alpha=0.85, width=0.6) +
  scale_fill_manual(values=c("TRUE"="#d7191c","FALSE"="#2c7bb6"), guide="none") +
  geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=0.3) +
  labs(title="C. C3 increase (early->late)", x="Percentage points", y="") +
  theme_pub()

# --- Combine ---
p_combined <- (pA + pB + pC) +
  plot_layout(ncol=3, widths=c(1.2, 1, 1.2)) &
  theme(legend.position="bottom")

ggsave(file.path(OUT, "combined_mge_context.pdf"), p_combined, width=16, height=7, dpi=300)
ggsave(file.path(OUT, "combined_mge_context.png"), p_combined, width=16, height=7, dpi=300)
cat("Saved combined figure.\n")

################################################################################
# 9. Summary text (top 10 genes)
################################################################################
sink(file.path(OUT, "plasmid_context_summary.txt"))
cat("Top-10 Increasing VF Genes in C3 - Summary\n")
cat("==========================================\n\n")
cat("Top 10 VF genes by increasing prevalence in C3 (early 2016-2018 -> late 2022-2025):\n\n")

cat("GENE-LEVEL LOCALIZATION (all ST69):\n")
for (i in 1:nrow(gene_loc)) {
  r <- gene_loc[i,]
  cat(sprintf("  %s (%s): WG=%.1f%%  Chromosomal=%.1f%%  Plasmid=%.1f%%  Absent=%.1f%%\n",
    r$gene, r$system, r$pct_wg, r$pct_chr, r$pct_pl, r$pct_abs))
}

cat("\nC3 TEMPORAL DELTAS:\n")
for (i in 1:nrow(delta_data)) {
  r <- delta_data[i,]
  cat(sprintf("  %s: %.1f%% -> %.1f%% (+%.1fpp)\n", r$gene, r$early_pct, r$late_pct, r$delta))
}

cat("\nRGP CO-LOCALIZATION:\n")
if (length(rgp_top_genes) >= 2) {
  cat(sprintf("  Genes in RGPs: %s\n", paste(rgp_top_genes, collapse=", ")))
  for (i in 1:nrow(gene_pairs)) {
    r <- gene_pairs[i,]
    cat(sprintf("  %s + %s: %d shared RGPs (%.1f%% of %s RGPs)\n",
      r$gene1, r$gene2, r$n_rgp_both, r$pct_of_g1, r$gene1))
  }
} else {
  cat("  No top-10 genes found in RGP data.\n")
}

cat("\nINTERPRETATION:\n")
cat("The top 10 increasing VF genes in C3 fall into two functional groups:\n")
cat("  1. kps capsule operon (kpsC/F/E/U/D/M): +36 to +42pp, all >99% chromosomal\n")
cat("  2. pap+sat (papX/B/I, sat): +20 to +28pp, mostly chromosomal\n\n")
cat("The kps genes are not found in RGPs (core capsule locus) while pap/sat\n")
cat("are RGP-associated (PAI). This supports a two-event model:\n")
cat("  - kps operon completion within the existing capsule locus\n")
cat("  - Acquisition/expansion of a pap+sat chromosomal PAI\n")
sink()
cat("\nDone.\n")
