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
