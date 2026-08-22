#!/bin/bash
# =============================================================================
# scripts/Plasmid_assembly/run_abricate_plasmids.sh
#
# Run abricate (VFDB and CARD, 80% identity / 80% coverage -- same thresholds
# as scripts/Annotation/run_abricate.sh) over the combined per-genome plasmid
# FASTAs produced by run_mobsuite.sh, then build the per-ST plasmid summary
# tables that scripts/Virulence/Plasmid_VFDB.R, scripts/AMR/Plasmid_CARD.R,
# and scripts/Analysis/Figure8.R read directly:
#
#   ../Plasmid/plasmid_vfdb_summary/{ST}_plasmid_summary.tsv
#   ../Plasmid/plasmid_card_summary/{ST}_plasmid_summary.tsv
#
# The per-genome abricate output filenames
# ("{base}_combined_plasmids_vfdb.tsv" / "..._card.tsv") are exactly what
# those three R scripts' clean_gid() helper already expects to strip -- this
# is not a new convention, it is the one the original analysis code was
# already written against.
#
# Usage:
#   bash scripts/Plasmid_assembly/run_abricate_plasmids.sh ST69        # single ST
#   bash scripts/Plasmid_assembly/run_abricate_plasmids.sh all         # all 5 STs (default)
#   JOBS=6 bash scripts/Plasmid_assembly/run_abricate_plasmids.sh ST69 # parallel jobs
#
# Requires: abricate (conda env 'wgs') -- run AFTER run_mobsuite.sh.
# Input:  ../Plasmid/combined_fasta/{ST}/{base}_combined_plasmids.fasta
# Output: ../Plasmid/plasmid_vfdb/{ST}/{base}_combined_plasmids_vfdb.tsv
#         ../Plasmid/plasmid_card/{ST}/{base}_combined_plasmids_card.tsv
#         ../Plasmid/plasmid_vfdb_summary/{ST}_plasmid_summary.tsv
#         ../Plasmid/plasmid_card_summary/{ST}_plasmid_summary.tsv
# =============================================================================
set -euo pipefail

JOBS="${JOBS:-6}"
ST_ARG="${1:-all}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PLASMID_ROOT="$(dirname "$REPO_ROOT")/Plasmid"

if [ "$ST_ARG" = "all" ]; then ST_LIST=("ST10" "ST69" "ST73" "ST95" "ST131"); else ST_LIST=("$ST_ARG"); fi

check_abricate() {
  if ! command -v abricate >/dev/null 2>&1; then
    echo "ERROR: abricate not found. Activate conda env 'wgs' first."
    exit 1
  fi
}
check_abricate

process_file() {
  local file="$1" vfdb_dir="$2" card_dir="$3"
  local base
  base=$(basename "$file" .fasta)
  if [[ ! -s "$vfdb_dir/${base}_vfdb.tsv" ]]; then
    abricate --db vfdb --minid 80 --mincov 80 "$file" > "$vfdb_dir/${base}_vfdb.tsv" 2>/dev/null
  fi
  if [[ ! -s "$card_dir/${base}_card.tsv" ]]; then
    abricate --db card --minid 80 --mincov 80 "$file" > "$card_dir/${base}_card.tsv" 2>/dev/null
  fi
}
export -f process_file

process_st() {
  local st="$1"
  local input_dir="$PLASMID_ROOT/combined_fasta/$st"
  if [ ! -d "$input_dir" ] || [ -z "$(ls "$input_dir"/*_combined_plasmids.fasta 2>/dev/null)" ]; then
    echo "[SKIP] $st : no combined plasmid FASTAs in $input_dir (run run_mobsuite.sh first)"
    return
  fi

  echo "=== $st : abricate VFDB + CARD on plasmid contigs ==="
  local vfdb_dir="$PLASMID_ROOT/plasmid_vfdb/$st"
  local card_dir="$PLASMID_ROOT/plasmid_card/$st"
  mkdir -p "$vfdb_dir" "$card_dir"

  local files=("$input_dir"/*_combined_plasmids.fasta)
  printf '%s\n' "${files[@]}" | \
    parallel --bar --jobs "$JOBS" --halt soon,fail=1 process_file {} "$vfdb_dir" "$card_dir"

  # Per-ST summary tables -- same abricate --summary format the whole-genome
  # pipeline already uses (wide table: #FILE + one column per gene, "." for
  # absent), which is exactly what clean_gid()/to_bin() in the downstream R
  # scripts expect.
  mkdir -p "$PLASMID_ROOT/plasmid_vfdb_summary" "$PLASMID_ROOT/plasmid_card_summary"
  abricate --summary "$vfdb_dir"/*_vfdb.tsv > "$PLASMID_ROOT/plasmid_vfdb_summary/${st}_plasmid_summary.tsv"
  abricate --summary "$card_dir"/*_card.tsv > "$PLASMID_ROOT/plasmid_card_summary/${st}_plasmid_summary.tsv"

  echo "  DONE $st : $(ls "$vfdb_dir"/*_vfdb.tsv | wc -l | tr -d ' ') genomes with plasmid annotations"
}

for st in "${ST_LIST[@]}"; do
  process_st "$st"
done

echo "=== plasmid abricate finished. Summaries in $PLASMID_ROOT/plasmid_{vfdb,card}_summary/ ==="
echo "Next: Rscript run_pipeline.R  (Virulence/Plasmid_VFDB.R, AMR/Plasmid_CARD.R, Analysis/Figure8.R will pick these up automatically)"
