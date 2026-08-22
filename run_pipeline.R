#!/usr/bin/env Rscript
#
# run_pipeline.R — Master analysis pipeline for E. coli ExPEC
#
# Runs the full set of canonical analysis scripts in dependency order:
#   prerequisites (finder summaries -> shell clustering -> tree mapping)
#   -> Virulence (VFDB / VirulenceFinder / plasmid VFDB)
#   -> AMR (CARD / ResFinder / plasmid CARD)
#   -> cluster-gene prerequisite for Figure 6
#   -> Analysis (Figure1.R .. Figure8.R, Supplementary.R)
#
# The scripts/ folder is organised one-script-per-database or
# one-script-per-figure (Virulence/, AMR/, Plasmid_assembly/, Analysis/),
# with a Prerequisites/ folder for the shared infrastructure every other
# folder depends on. See scripts/README.md for the full script -> purpose
# -> figure/table map and the confirmed figure->script evidence trail.
#
# Usage:
#   Rscript run_pipeline.R                      # all steps for ST69
#   TARGET_ST=ST10 Rscript run_pipeline.R       # different lineage
#   Rscript run_pipeline.R ST95                  # positional arg also works
#   Rscript run_pipeline.R --quick               # clustering + mapping only
#   Rscript run_pipeline.R --from AMR/CARD.R     # resume at a given script
#   Rscript run_pipeline.R --skip Analysis/Figure3.R   # skip a script
#
# Pre-processing (bash — run separately, see run_all.sh):
#   download assemblies  ->  abricate (VFDB/CARD)  ->  VF/ResFinder
#   -> Prerequisites/00_build_finder_summaries.R  -> pangenome/tree (bash)

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

# run_script(): `path` is the full path relative to scripts/, e.g.
# "Virulence/VFDB.R" -- used directly for --from/--skip matching, console
# labels, and the source() call.
run_script <- function(label, path, optional = FALSE) {
  if (quick_mode && !grepl("cluster|mapping", label, ignore.case = TRUE)) return(invisible(NULL))
  if (!is.null(resume_from) && path < resume_from) return(invisible(NULL))
  if (path %in% skip_set) { cat(sprintf("\n>>> SKIPPED (--skip): %s %s\n", path, label)); return(invisible(NULL)) }
  full_path <- file.path("scripts", path)
  cat(sprintf("\n>>> [%s] %s <<<\n", path, label))
  if (!file.exists(full_path)) {
    if (optional) { cat("  SKIPPED (optional, not present):", full_path, "\n"); return(invisible(NULL)) }
    stop("Script not found: ", full_path)
  }
  t0 <- Sys.time()
  tryCatch({
    source(full_path, local = TRUE)
    cat(sprintf("  DONE in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }, error = function(e) {
    cat(sprintf("  ERROR in %.1f min: %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), e$message))
    if (!optional) stop("Pipeline aborted at ", path, " (", label, ")")
  })
}

# === PHASE 1: Prerequisites (finder summaries, shell clustering, tree mapping) ===
run_script("Build finder summary matrices", "Prerequisites/00_build_finder_summaries.R")
run_script("Pangenome shell cluster (VFDB)", "Prerequisites/02b_pangenome_shell_cluster_vfdb.R")
run_script("Pangenome shell cluster (VF)",  "Prerequisites/02c_pangenome_shell_cluster_vf.R")

if (quick_mode) { cat(bar("QUICK MODE — stopping after clustering.")); quit(save = "no", status = 0) }

run_script("Map clusters to core tree (VFDB)", "Prerequisites/03c_tree_mapping_vfdb.R")
run_script("Map clusters to core tree (VF)",  "Prerequisites/03d_tree_mapping_vf.R")

# === PHASE 2: Virulence (one script per database) ===
run_script("Virulence: VFDB",           "Virulence/VFDB.R")
run_script("Virulence: VirulenceFinder", "Virulence/VirulenceFinder.R")
run_script("Virulence: Plasmid VFDB",   "Virulence/Plasmid_VFDB.R", optional = TRUE)  # needs external Plasmid dataset

# === PHASE 3: AMR (one script per database) ===
run_script("AMR: CARD",       "AMR/CARD.R")
run_script("AMR: ResFinder",  "AMR/ResFinder.R")  # includes 05b/05c cluster + decreasing-trend follow-ons, see file header
run_script("AMR: Plasmid CARD", "AMR/Plasmid_CARD.R", optional = TRUE)  # placeholder -- no such analysis exists yet, see file header

# === PHASE 4: Plasmid assembly ===
# No R script here -- mob-suite is run externally (see Plasmid_assembly/README.md).
# Its output feeds Virulence/Plasmid_VFDB.R and Analysis/Figure8.R above/below.

# === PHASE 5: Cluster-gene prerequisite for Figure 6 ===
# Organisationally lives in Prerequisites/, but must run AFTER Virulence/VFDB.R
# since it reads that script's master shell-cluster table.
run_script("Cluster gene analysis (for Figure 6)", "Prerequisites/04b_cluster_gene_analysis.R")

# === PHASE 6: Analysis — one script per published figure, plus Supplementary.R ===
run_script("Figure 1: temporal trend + Mann-Kendall", "Analysis/Figure1.R")
run_script("Figure 2: shell-cluster silhouette + gene content", "Analysis/Figure2.R")
run_script("Figure 3: iTOL circular phylogenies (placeholder, see file header)", "Analysis/Figure3.R", optional = TRUE)
run_script("Figure 4: Oaxaca-style decomposition", "Analysis/Figure4.R")
run_script("Figure 5: per-cluster temporal burden trends", "Analysis/Figure5.R")
run_script("Figure 6: 8-panel ST69/ST10 summary", "Analysis/Figure6.R")
run_script("Figure 7: sensitivity analysis", "Analysis/Figure7.R")
run_script("Figure 8: plasmid/RGP context", "Analysis/Figure8.R", optional = TRUE)  # needs external Plasmid dataset
run_script("Supplementary: all remaining analyses/tables/figures", "Analysis/Supplementary.R")

cat(bar(sprintf("PIPELINE COMPLETE for %s", config$TARGET_ST)))
cat("Outputs in:", config$OUTPUT_DIR, "\n")
