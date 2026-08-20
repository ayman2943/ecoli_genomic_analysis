#!/usr/bin/env Rscript
#
# run_pipeline.R — Master analysis pipeline for E. coli ExPEC
#
# Runs the full set of canonical analysis scripts in dependency order:
#   shell clustering -> tree mapping -> virulence -> resistance -> temporal/MK
#   -> ST10 decomposition -> clinical enrichment -> capsule/K-type -> validation -> figures
#
# The scripts/ folder is organised into one subfolder per pipeline stage
# (01_pangenome_clustering, 02_tree_mapping, ... 09_figures) so the stage
# a script belongs to is visible from its path, not just its filename.
# See scripts/README.md for the full script -> purpose -> figure/table map.
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
#   download assemblies  ->  abricate (VFDB/CARD)  ->  VF/ResFinder
#   -> build finder summaries (R)  ->  pangenome/tree (bash)

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

# run_script(): `script` is still the bare filename (e.g. "09e_summary_stats.R")
# so --from/--skip matching and console labels are unchanged; `folder` is the
# new scripts/ subfolder it now lives in. The two are joined only to build
# the actual path passed to source().
run_script <- function(label, script, folder, optional = FALSE) {
  if (quick_mode && !grepl("cluster|mapping", label, ignore.case = TRUE)) return(invisible(NULL))
  if (!is.null(resume_from) && script < resume_from) return(invisible(NULL))
  if (script %in% skip_set) { cat(sprintf("\n>>> SKIPPED (--skip): %s %s\n", script, label)); return(invisible(NULL)) }
  path <- file.path("scripts", folder, script)
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
run_script("Pangenome shell cluster (VFDB)", "02b_pangenome_shell_cluster_vfdb.R", "01_pangenome_clustering")
run_script("Pangenome shell cluster (VF)",  "02c_pangenome_shell_cluster_vf.R",   "01_pangenome_clustering")

if (quick_mode) { cat(bar("QUICK MODE — stopping after clustering.")); quit(save = "no", status = 0) }

# === PHASE 2: Tree mapping ===
run_script("Map clusters to core tree (VFDB)", "03c_tree_mapping_vfdb.R", "02_tree_mapping")
run_script("Map clusters to core tree (VF)",  "03d_tree_mapping_vf.R",   "02_tree_mapping")

# === PHASE 3: Virulence + resistance analysis ===
run_script("Virulence analysis (VFDB + VF)", "04a_virulence_analysis.R",           "03_virulence_resistance")
run_script("Cluster gene analysis",          "04b_cluster_gene_analysis.R",        "03_virulence_resistance")
run_script("Resistance analysis (CARD + ResFinder)", "05a_resistance_analysis.R",  "03_virulence_resistance")
run_script("ResFinder cluster analysis",     "05b_resfinder_clusters.R",           "03_virulence_resistance")
run_script("ResFinder decreasing trend",     "05c_resfinder_decreasing.R",         "03_virulence_resistance")

# === PHASE 4: Temporal trend + Mann-Kendall ===
run_script("Temporal trend + Mann-Kendall",  "06a_temporal_mk.R", "04_temporal_trends")

# === PHASE 5: ST10 decomposition (multilineage context) ===
run_script("ST10 pangenome decomposition",   "07a_st10_decomposition.R",        "05_st10_decomposition")
run_script("ST10 composition drivers",       "07b_st10_composition_drivers.R",  "05_st10_decomposition")

# === PHASE 6: Clinical enrichment + sensitivity ===
run_script("Clinical enrichment (3-panel)",  "08_clinical_enrichment_3panel.R", "06_clinical_sensitivity")
run_script("Sensitivity analysis",           "09a_sensitivity_analysis.R",      "06_clinical_sensitivity")

# === PHASE 7: Capsule / K-type classification + summary stats ===
run_script("Capsule classification",         "09b_capsule_classification.R", "07_capsule_ktype")
run_script("Capsule comparison",             "09c_capsule_comparison.R",     "07_capsule_ktype")
run_script("K-type analysis",                "09d_k_type_analysis.R",        "07_capsule_ktype")
run_script("Summary statistics",             "09e_summary_stats.R",          "07_capsule_ktype")

# === PHASE 8: Validation analyses ===
run_script("KPS validation",                 "10a_kps_validation.R",         "08_validation")
run_script("RGP neighbourhood",              "10b_rgp_neighbourhood.R",      "08_validation")
run_script("RGP genomic context",            "10c_rgp_genomic_context.R",    "08_validation")
run_script("Plasmid context",                "10d_plasmid_context.R",       "08_validation", optional = TRUE)
run_script("Allelic conversion",             "10e_allelic_conversion.R",     "08_validation")
run_script("ST10 ARG decomposition",         "10f_st10_arg_decomposition.R", "05_st10_decomposition")
run_script("Long-read validation",           "10g_long_read_validation.R",   "08_validation")

# === PHASE 9: Figure scripts ===
# NOTE (deduplicated): Figure 1 (temporal trends), the Figure 7 sensitivity
# figure, the kps-validation figure, and Supplementary Figure S2 (allelic
# conversion) are produced as a side effect of their analysis-stage scripts
# above (06a_temporal_mk.R, 09a_sensitivity_analysis.R, 10a_kps_validation.R,
# 10e_allelic_conversion.R respectively). Those figures were previously
# ALSO generated a second time from byte-identical duplicate scripts
# (11_/17_/20_/21_) kept under different "figNN" names — pure dead weight
# that doubled runtime and could silently drift out of sync with the
# analysis scripts if only one copy were ever edited. The duplicates have
# been removed from this pipeline.
#
# NOTE (unclear provenance — see scripts/README.md "Figures: open questions"):
# 15_fig07_tree_parsimony.R and 16_fig07_phylogeny_summary.R are both
# labelled "fig07" but neither one's output filename matches Figure 7
# (Sensitivity analysis, produced by 09a_sensitivity_analysis.R above). By
# content, 15_ appears to build the core-genome tree figure (Figure 3) and
# 16_ appears to build the 8-panel ST69/ST10 summary figure (Figure 6) —
# but this is inferred from the code, not confirmed against which script
# was actually run to produce the submitted Figure_03.png/Figure_06.png.
# Kept under their original names and numbers until that is confirmed.
run_script("Figures 2-5: ST69 analysis",              "12_fig02-05_ST69_analysis.R",      "09_figures")
run_script("Figure 6: cluster temporal",              "13_fig06_cluster_temporal.R",      "09_figures")
run_script("Figure 6: gene trajectories",             "14_fig06_gene_trajectories.R",     "09_figures")
run_script("Figure 7: tree parsimony",                "15_fig07_tree_parsimony.R",        "09_figures", optional = TRUE)
run_script("Figure 7: phylogeny summary",             "16_fig07_phylogeny_summary.R",     "09_figures", optional = TRUE)
# NOTE: "Figure 9" (clinical enrichment) is not a standalone figure in the
# published manuscript -- its content was folded into Figure 7's middle
# panel (clinical proportion over time). 08_clinical_enrichment_3panel.R
# above still produces its own 3-panel figure for validation/QC purposes.
run_script("Composite figures",                       "19_fig_composite.R",               "09_figures", optional = TRUE)
run_script("Supplementary combined figure",           "22_fig_supplementary_combined.R",  "09_figures")

# === PHASE 10: Combined figure output (optional) ===
run_script("Combine all figures",                     "99_ALL_FIGURES_COMBINED.R", "09_figures", optional = TRUE)

cat(bar(sprintf("PIPELINE COMPLETE for %s", config$TARGET_ST)))
cat("Outputs in:", config$OUTPUT_DIR, "\n")
