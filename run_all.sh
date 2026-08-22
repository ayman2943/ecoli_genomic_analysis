#!/bin/bash
# =============================================================================
# run_all.sh — End-to-end E. coli ExPEC pipeline
#
# 01. Filter EnteroBase metadata        -> metadata/{ST}_filtered.xlsx
# 02. Download assemblies               -> {ST}/Escherichia_coli_*.fna.gz
# 03. abricate VFDB + CARD              -> card_vfdb_result/
# 04. VirulenceFinder + ResFinder       -> analysis_results/{ST}/...
# 05. Build finder summary matrices     -> finder_result/
# 06. Plasmid assembly (optional)       -> mob-suite + abricate -> ../Plasmid/
# 07. Pangenome + core tree             -> pangenome_output/ + {ST}_bootstrap.treefile
# 08. R analysis + figures              -> output/
#
# Usage:
#   bash run_all.sh                  # all steps for ST69 (assumes data staged)
#   bash run_all.sh ST10             # different lineage
#   bash run_all.sh ST69 --download  # only metadata+download
#   bash run_all.sh ST69 --annotate  # only abricate + VF/RF + summaries
#   bash run_all.sh ST69 --plasmid   # only mob-suite + plasmid abricate
#   bash run_all.sh ST69 --pangenome # only pangenome + tree
#   bash run_all.sh ST69 --analyze   # only the R analysis/figures
#
# Environment: activate conda env 'wgs' before running; the R scripts need the
# packages listed in install.R. --plasmid additionally needs mob_suite.
# =============================================================================
set -euo pipefail

ST="${1:-ST69}"
MODE="${2:-full}"

# ---- Config (paths, threads, memory) ----
if [ -f config/pipeline_config.sh ]; then
  source config/pipeline_config.sh
fi
export TARGET_ST="$ST"

case "$MODE" in
  --download) echo "[run_all] MODE=download ($ST)";;
  --annotate) echo "[run_all] MODE=annotate ($ST)";;
  --plasmid)  echo "[run_all] MODE=plasmid ($ST)";;
  --pangenome) echo "[run_all] MODE=pangenome ($ST)";;
  --analyze)   echo "[run_all] MODE=analyze ($ST)";;
  full)        echo "[run_all] MODE=full ($ST)";;
  *) echo "Unknown mode: $MODE (full | --download | --annotate | --plasmid | --pangenome | --analyze)"; exit 1;;
esac

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found"; exit 1; }; }

# ---- 01. Metadata ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--download" ]; then
  if [ -z "${RAW_ENTERO_EXPORT:-}" ]; then
    echo "[run_all] 01: RAW_ENTERO_EXPORT not set — skipping metadata filter."
    echo "           Set RAW_ENTERO_EXPORT=/path/to/EnteroBase_export.xlsx and re-run."
  else
    need Rscript
    echo "[run_all] 01: filtering metadata -> metadata/${ST}_filtered.xlsx"
    Rscript scripts/Setup/filter_metadata.R "$RAW_ENTERO_EXPORT"
  fi
fi

# ---- 02. Download ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--download" ]; then
  if [ -d "$ST" ] && [ -n "$(ls "$ST"/*.fna* 2>/dev/null)" ]; then
    echo "[run_all] 02: assemblies already present in $ST — skipping download."
  else
    need python3
    echo "[run_all] 02: downloading assemblies for $ST"
    python3 scripts/Setup/download_assemblies.py \
      "metadata/${ST}_filtered.xlsx" "$ST" "${SPECIES:-Escherichia_coli}"
    # Decompress for abricate/VF/ResFinder (they read plain .fna)
    echo "[run_all] 02: decompressing .fna.gz"
    for f in "$ST"/*.fna.gz; do
      [ -f "$f" ] && gzip -df "$f"
    done
  fi
fi

if [ "$MODE" = "--download" ]; then exit 0; fi

# ---- 03. abricate VFDB + CARD ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--annotate" ]; then
  need abricate; need parallel
  echo "[run_all] 03: abricate VFDB + CARD for $ST"
  bash scripts/Annotation/run_abricate.sh "$ST"
fi

# ---- 04. VirulenceFinder + ResFinder ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--annotate" ]; then
  need python; need run_resfinder.py
  echo "[run_all] 04: VirulenceFinder + ResFinder for $ST"
  bash scripts/Annotation/run_vf_resfinder.sh "$ST"
fi

# ---- 05. Finder summary matrices ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--annotate" ]; then
  echo "[run_all] 05: building finder_result/ summary matrices"
  Rscript scripts/Prerequisites/00_build_finder_summaries.R "analysis_results" "finder_result"
fi

if [ "$MODE" = "--annotate" ]; then exit 0; fi

# ---- 06. Plasmid assembly (optional: mob-suite + abricate on plasmid contigs) ----
if [ "$MODE" = "--plasmid" ] || { [ "$MODE" = "full" ] && command -v mob_recon >/dev/null 2>&1; }; then
  need abricate; need parallel
  if ! command -v mob_recon >/dev/null 2>&1; then
    echo "ERROR: mob_recon not found. Install with: conda install -c bioconda mob_suite"
    exit 1
  fi
  echo "[run_all] 06: mob-suite + plasmid abricate for $ST"
  bash scripts/Plasmid_assembly/run_mobsuite.sh "$ST"
  bash scripts/Plasmid_assembly/run_abricate_plasmids.sh "$ST"
elif [ "$MODE" = "full" ]; then
  echo "[run_all] 06: mob_recon not found — skipping optional plasmid stage."
  echo "           Install mob_suite and re-run with --plasmid to add it later;"
  echo "           Virulence/Plasmid_VFDB.R, AMR/Plasmid_CARD.R, and Figure8.R"
  echo "           will pick it up automatically once it's present."
fi

if [ "$MODE" = "--plasmid" ]; then exit 0; fi

# ---- 07. Pangenome + tree (needed for clustering/mapping steps) ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--pangenome" ]; then
  if [ ! -f "${ST}_bootstrap.treefile" ]; then
    echo "[run_all] 07: running pangenome + core tree for $ST"
    bash scripts/Pangenome/run_pangenome_tree.sh "$ST" "${PANG_N_GENOMES:-0}"
  else
    echo "[run_all] 07: ${ST}_bootstrap.treefile present — skipping pangenome."
  fi
fi

if [ "$MODE" = "--pangenome" ]; then exit 0; fi

# ---- 08. R analysis + figures ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--analyze" ]; then
  need Rscript
  echo "[run_all] 08: running analysis pipeline for $ST"
  Rscript run_pipeline.R "$ST"
fi

echo "[run_all] DONE ($ST)"
