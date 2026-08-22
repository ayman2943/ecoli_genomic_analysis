#!/bin/bash
# =============================================================================
# 02_annotation/run_abricate.sh
#
# Run abricate (VFDB and CARD) over all downloaded assemblies for one or all STs
# at 80% identity / 80% coverage, then produce per-ST summary tables.
#
# Usage:
#   bash 02_annotation/run_abricate.sh ST69            # single ST
#   bash 02_annotation/run_abricate.sh all             # all 5 STs (default)
#   JOBS=6 bash 02_annotation/run_abricate.sh ST69     # parallel jobs
#
# Requires: abricate (conda env 'wgs')
# Outputs:
#   card_vfdb_result/vfdb_summary/{ST}_vfdb_summary.tsv
#   card_vfdb_result/card_summary/{ST}.tsv
# =============================================================================
set -euo pipefail

JOBS="${JOBS:-6}"
ST_ARG="${1:-all}"

if [ "$ST_ARG" = "all" ]; then ST_LIST=("ST10" "ST69" "ST73" "ST95" "ST131"); else ST_LIST=("$ST_ARG"); fi

check_abricate() {
  if ! command -v abricate >/dev/null 2>&1; then
    echo "ERROR: abricate not found. Activate conda env 'wgs' first."
    exit 1
  fi
}
check_abricate

# CARD full-table output (all columns). abricate's default card output already
# includes every CARD column; we keep the raw per-genome table for downstream.
card_raw() {
  local db="$1" infile="$2" outfile="$3"
  abricate --db "$db" --minid 80 --mincov 80 "$infile" > "$outfile" 2>/dev/null
}

process_st() {
  local st="$1"
  local input_dir="$st"
  if [ ! -d "$input_dir" ] || [ -z "$(ls "$input_dir"/*.fna 2>/dev/null)" ]; then
    echo "[SKIP] $st : no .fna files in $input_dir"
    return
  fi

  echo "=== $st : abricate VFDB + CARD ==="
  local vfdb_dir="card_vfdb_result/vfdb/${st}"
  local card_dir="card_vfdb_result/card/${st}"
  mkdir -p "$vfdb_dir" "$card_dir"

  local files=("$input_dir"/*.fna)
  local total=${#files[@]}

  process_file() {
    local file="$1"
    local base
    base=$(basename "$file" .fna)
    if [[ ! -s "$vfdb_dir/${base}.tsv" ]]; then
      abricate --db vfdb --minid 80 --mincov 80 "$file" > "$vfdb_dir/${base}.tsv" 2>/dev/null
    fi
    if [[ ! -s "$card_dir/${base}.tsv" ]]; then
      card_raw card "$file" "$card_dir/${base}.tsv"
    fi
  }
  export -f process_file
  export vfdb_dir card_dir

  printf '%s\n' "${files[@]}" | \
    parallel --bar --jobs "$JOBS" --halt soon,fail=1 process_file {}

  # Per-ST summary tables (rows = genes, columns = genomes)
  mkdir -p card_vfdb_result/vfdb_summary card_vfdb_result/card_summary
  abricate --summary "$vfdb_dir"/*.tsv > "card_vfdb_result/vfdb_summary/${st}_vfdb_summary.tsv"

  # CARD burden-style summary (transposed, per-genome ARG count is derived downstream)
  abricate --summary "$card_dir"/*.tsv > "card_vfdb_result/card_summary/${st}.tsv"

  echo "  DONE $st : $(ls "$vfdb_dir"/*.tsv | wc -l | tr -d ' ') VFDB, $(ls "$card_dir"/*.tsv | wc -l | tr -d ' ') CARD genomes"
}

for st in "${ST_LIST[@]}"; do
  process_st "$st"
done

echo "=== abricate finished ==="
