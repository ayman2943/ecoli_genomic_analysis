#!/usr/bin/env Rscript
#
# run_pipeline.R — Master analysis pipeline for E. coli ExPEC
#
# Runs the full set of canonical analysis scripts in dependency order:
#   shell clustering → tree mapping → virulence → resistance → temporal/MK
#   → ST10 decomposition → clinical enrichment → figures → validation
#
# Usage:
#   Rscript run_pipeline.R                      # all steps for ST69
#   TARGET_ST=ST10 Rscript run_pipeline.R       # different lineage
#   Rscript run_pipeline.R ST95                  # positional arg also works
#   Rscript run_pipeline.R --quick               # clustering + mapping only
#   Rscript run_pipeline.R --from 09e            # resume at a given script
#   Rscript run_pipeline.R --skip 15 --skip 16   # skip tree figures (no ggtree)
#
# Pre-processing (bash — run separately, see run_all.sh):
#   download assemblies  →  abricate (VFDB/CARD)  →  VF/ResFinder
#   → build finder summaries (R)  →  pangenome/tree (bash)

args <- commandArgs(trailingOnly = TRUE)
quick_mode  <- "--quick" %in% args
resume_from <- sub("^--from=", "", args[grep("^--from=", args)])
skip_set    <- sub("^--skip=", "", args[grep("^--skip=", args)])
if (length(resume_from) == 0 && "--from" %in% args) {
  resume_from <- args[which(args == "--from") + 1]
}
resume_from <- if (length(resume_from)) resume_from[1] else NULL
skip_set    <- if (length(skip_set)) strsplit(paste(skip_set, collapse = ","), ",")[[1]] else character()

pos <- args[!args %in% c("--quick", "--from", skip_set, resume_from, " ")]
if (length(pos) > 0 && grepl("^ST", pos[1])) Sys.setenv(TARGET_ST = pos[1])

if (file.exists("config.R")) {
  source("config.R")
} else if (file.exists("config/config.R")) {
  source("config/config.R")
} else {
  stop("config.R not found in CWD or config/")
}

bar <- function(x) sprintf("\n%s\n%s\n%s\n", paste0(rep("=", 70), collapse = ""), x,
                          paste0(rep("=", 70), collapse = ""))
cat(bar(sprintf("PIPELINE: %s   Mode: %s   Resume: %s",
                config$TARGET_ST,
                if (quick_mode) "QUICK" else "FULL",
                if (is.null(resume_from)) "from start" else resume_from)))

run_script <- function(label, script, optional = FALSE) {
  if (quick_mode && !grepl("cluster|mapping", label, ignore.case = TRUE)) return(invisible(NULL))
  if (!is.null(resume_from) && script < resume_from) return(invisible(NULL))
  if (script %in% skip_set) { cat(sprintf("\n>>> SKIPPED (--skip): %s %s\n", script, label)); return(invisible(NULL)) }
  path <- file.path("scripts", script)
  cat(sprintf("\n>>> [%s] %s <<<\n", script, label))
  if (!file.exists(path)) {
    if (optional) { cat("  SKIPPED (optional, not present):", path, "\n"); return(invisible(NULL)) }
    stop("Script not found: ", path)
  }
  t0 <- Sys.time()
  tryCatch({
    source(path, local = TRUE)
    cat(sprintf("  DONE in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }, error = function(e) {
    cat(sprintf("  ERROR in %.1f min: %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), e$message))
    if (!optional) stop("Pipeline aborted at ", script, " (", label, ")")
  })
}

# === PHASE 1: Shell-gene clustering ===
run_script("Pangenome shell cluster (VFDB)", "02b_pangenome_shell_cluster_vfdb.R")
run_script("Pangenome shell cluster (VF)",  "02c_pangenome_shell_cluster_vf.R")

if (quick_mode) { cat(bar("QUICK MODE — stopping after clustering.")); quit(save = "no", status = 0) }

# === PHASE 2: Tree mapping ===
run_script("Map clusters to core tree (VFDB)", "03c_tree_mapping_vfdb.R")
run_script("Map clusters to core tree (VF)",  "03d_tree_mapping_vf.R")

# === PHASE 3: Virulence + resistance analysis ===
run_script("Virulence analysis (VFDB + VF)", "04a_virulence_analysis.R")
run_script("Cluster gene analysis",          "04b_cluster_gene_analysis.R")
run_script("Resistance analysis (CARD + ResFinder)", "05a_resistance_analysis.R")
run_script("ResFinder cluster analysis",     "05b_resfinder_clusters.R")
run_script("ResFinder decreasing trend",     "05c_resfinder_decreasing.R")

# === PHASE 4: Temporal trend + Mann-Kendall ===
run_script("Temporal trend + Mann-Kendall",  "06a_temporal_mk.R")

# === PHASE 5: ST10 decomposition (multilineage context) ===
run_script("ST10 pangenome decomposition",   "07a_st10_decomposition.R")
run_script("ST10 composition drivers",       "07b_st10_composition_drivers.R")

# === PHASE 6: Clinical enrichment ===
run_script("Clinical enrichment (3-panel)",  "08_clinical_enrichment_3panel.R")

# === PHASE 7: Supplementary analyses ===
run_script("Sensitivity analysis",           "09a_sensitivity_analysis.R")
run_script("Capsule classification",         "09b_capsule_classification.R")
run_script("Capsule comparison",             "09c_capsule_comparison.R")
run_script("K-type analysis",                "09d_k_type_analysis.R")
run_script("Summary statistics",             "09e_summary_stats.R")

# === PHASE 8: Validation analyses ===
run_script("KPS validation",                 "10a_kps_validation.R")
run_script("RGP neighbourhood",              "10b_rgp_neighbourhood.R")
run_script("RGP genomic context",            "10c_rgp_genomic_context.R")
run_script("Plasmid context",                "10d_plasmid_context.R",     optional = TRUE)
run_script("Allelic conversion",             "10e_allelic_conversion.R")
run_script("ST10 ARG decomposition",         "10f_st10_arg_decomposition.R")
run_script("Long-read validation",           "10g_long_read_validation.R")

# === PHASE 9: Figure scripts ===
run_script("Figure 1: temporal trends",               "11_fig01_temporal_trends.R")
run_script("Figures 2–5: ST69 analysis",              "12_fig02-05_ST69_analysis.R")
run_script("Figure 6: cluster temporal",              "13_fig06_cluster_temporal.R")
run_script("Figure 6: gene trajectories",             "14_fig06_gene_trajectories.R")
run_script("Figure 7: tree parsimony",                "15_fig07_tree_parsimony.R",     optional = TRUE)
run_script("Figure 7: phylogeny summary",             "16_fig07_phylogeny_summary.R",  optional = TRUE)
run_script("Figure 8: clinical sensitivity",          "17_fig08_clinical_sensitivity.R")
run_script("Figure 9: clinical enrichment",           "18_fig09_clinical_enrichment.R")
run_script("Composite figures",                     "19_fig_composite.R",     optional = TRUE)
run_script("KPS validation figure",                   "20_fig_kps_validation.R")
run_script("Supplementary S2 figure",                 "21_fig_supplementary_S2.R")
run_script("Supplementary combined figure",           "22_fig_supplementary_combined.R")

# === PHASE 10: Combined figure output (optional) ===
run_script("Combine all figures",                     "99_ALL_FIGURES_COMBINED.R", optional = TRUE)

cat(bar(sprintf("PIPELINE COMPLETE for %s", config$TARGET_ST)))
cat("Outputs in:", config$OUTPUT_DIR, "\n")
