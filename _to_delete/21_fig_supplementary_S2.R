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
