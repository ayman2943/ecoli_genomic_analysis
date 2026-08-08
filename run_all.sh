#!/usr/bin/env bash
# =============================================================================
# run_all.sh — End-to-end E. coli ExPEC pipeline
#
# 00. Set up environment & config
# 01. Filter EnteroBase metadata        -> metadata/{ST}_filtered.xlsx
# 02. Download assemblies               -> {ST}/Escherichia_coli_*.fna.gz
# 03. abricate VFDB + CARD              -> card_vfdb_result/
# 04. VirulenceFinder + ResFinder       -> analysis_results/{ST}/...
# 05. Build finder summary matrices     -> finder_result/
# 06. Pangenome + core tree (optional)  -> pangenome_output/ + {ST}_bootstrap.treefile
# 07. R analysis + figures              -> output/
#
# Usage:
#   bash run_all.sh                  # all steps for ST69 (assumes data staged)
#   bash run_all.sh ST10             # different lineage
#   bash run_all.sh ST69 --download  # only metadata+download
#   bash run_all.sh ST69 --annotate  # only abricate + VF/RF + summaries
#   bash run_all.sh ST69 --pangenome # only pangenome + tree
#   bash run_all.sh ST69 --analyze   # only the R analysis/figures
#
# Environment: activate conda env 'wgs' before running; the R scripts need the
# packages listed in install.R.
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
  --pangenome) echo "[run_all] MODE=pangenome ($ST)";;
  --analyze)   echo "[run_all] MODE=analyze ($ST)";;
  full)        echo "[run_all] MODE=full ($ST)";;
  *) echo "Unknown mode: $MODE (full | --download | --annotate | --pangenome | --analyze)"; exit 1;;
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
    Rscript 00_metadata/filter_metadata.R "$RAW_ENTERO_EXPORT"
  fi
fi

# ---- 02. Download ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--download" ]; then
  if [ -d "$ST" ] && [ -n "$(ls "$ST"/*.fna* 2>/dev/null)" ]; then
    echo "[run_all] 02: assemblies already present in $ST — skipping download."
  else
    need python3
    echo "[run_all] 02: downloading assemblies for $ST"
    python3 01_download/download_assemblies.py \
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
  bash 02_annotation/run_abricate.sh "$ST"
fi

# ---- 04. VirulenceFinder + ResFinder ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--annotate" ]; then
  need python; need run_resfinder.py
  echo "[run_all] 04: VirulenceFinder + ResFinder for $ST"
  bash 02_annotation/run_vf_resfinder.sh "$ST"
fi

# ---- 05. Finder summary matrices ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--annotate" ]; then
  echo "[run_all] 05: building finder_result/ summary matrices"
  Rscript 02_annotation/build_finder_summaries.R "analysis_results" "finder_result"
fi

if [ "$MODE" = "--annotate" ]; then exit 0; fi

# ---- 06. Pangenome + tree (needed for clustering/mapping steps) ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--pangenome" ]; then
  if [ ! -f "${ST}_bootstrap.treefile" ]; then
    echo "[run_all] 06: running pangenome + core tree for $ST"
    bash 03_pangenome/run_pangenome_tree.sh "$ST" "${PANG_N_GENOMES:-0}"
  else
    echo "[run_all] 06: ${ST}_bootstrap.treefile present — skipping pangenome."
  fi
fi

if [ "$MODE" = "--pangenome" ]; then exit 0; fi

# ---- 07. R analysis + figures ----
if [ "$MODE" = "full" ] || [ "$MODE" = "--analyze" ]; then
  need Rscript
  echo "[run_all] 07: running analysis pipeline for $ST"
  Rscript run_pipeline.R "$ST"
fi

echo "[run_all] DONE ($ST)"
