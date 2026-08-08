#!/bin/bash
# =============================================================================
# 02_annotation/run_vf_resfinder.sh
#
# Run VirulenceFinder and ResFinder over all downloaded assemblies for one or
# all STs at 80% identity / 80% coverage. Skips genomes whose output already
# exists, so it is safe to resume.
#
# Usage:
#   bash 02_annotation/run_vf_resfinder.sh ST69        # single ST
#   bash 02_annotation/run_vf_resfinder.sh all         # all 5 STs (default)
#   JOBS=4 bash 02_annotation/run_vf_resfinder.sh ST69
#
# Requires: conda env 'wgs' (python -m virulencefinder, run_resfinder.py),
#           databases at $HOME/databases/{virulencefinder_db,resfinder_db}
# Outputs:
#   analysis_results/{ST}/virulence/{genome}_VF.tsv
#   analysis_results/{ST}/resfinder/{genome}_RF.txt
# =============================================================================
set -euo pipefail

JOBS="${JOBS:-4}"
ST_ARG="${1:-all}"

if [ "$ST_ARG" = "all" ]; then ST_LIST=("ST10" "ST69" "ST73" "ST95" "ST131"); else ST_LIST=("$ST_ARG"); fi

BASE_OUT="analysis_results"

VF_DB="${VF_DB:-$HOME/databases/virulencefinder_db}"
RF_DB="${RF_DB:-$HOME/databases/resfinder_db}"

ID=0.80
COV=0.80

mkdir -p "$BASE_OUT"

process_vf() {
    local file=$1
    local st=$2
    local base
    base=$(basename "$file" .fna)
    local tmp="$BASE_OUT/$st/vf_tmp/${base}"
    local final="$BASE_OUT/$st/virulence/${base}_VF.tsv"

    if [[ -s "$final" ]]; then
        echo "[VF SKIP] $st / $base"
        return
    fi

    mkdir -p "$tmp"

    echo "[VF RUN ] $st / $base"

    python -m virulencefinder \
        -ifa "$file" -o "$tmp" -p "$VF_DB" -x -l "$COV" -t "$ID" \
        > /dev/null 2>&1

    if [[ -f "$tmp/results_tab.tsv" ]]; then
        mv "$tmp/results_tab.tsv" "$final"
    fi

    rm -rf "$tmp"
    echo "[VF DONE] $st / $base"
}

process_rf() {
    local file=$1
    local st=$2
    local base
    base=$(basename "$file" .fna)
    local tmp="$BASE_OUT/$st/rf_tmp/${base}"
    local final="$BASE_OUT/$st/resfinder/${base}_RF.txt"

    if [[ -s "$final" ]]; then
        echo "[RF SKIP] $st / $base"
        return
    fi

    mkdir -p "$tmp"

    echo "[RF RUN ] $st / $base"

    run_resfinder.py \
        -ifa "$file" -o "$tmp" -db_res "$RF_DB" -acq -l "$COV" -t "$ID" \
        > /dev/null 2>&1

    if [[ -f "$tmp/ResFinder_results_tab.txt" ]]; then
        mv "$tmp/ResFinder_results_tab.txt" "$final"
    fi

    rm -rf "$tmp"
    echo "[RF DONE] $st / $base"
}

export -f process_vf process_rf
export VF_DB RF_DB ID COV BASE_OUT JOBS

for ST in "${ST_LIST[@]}"; do
    INPUT_DIR="$ST"
    if [ ! -d "$INPUT_DIR" ] || [ -z "$(ls "$INPUT_DIR"/*.fna 2>/dev/null)" ]; then
        echo "[WARN] No genomes found in $ST"
        continue
    fi

    echo "=== $ST : VirulenceFinder + ResFinder ==="
    files=("$INPUT_DIR"/*.fna)
    total=${#files[@]}
    mkdir -p "$BASE_OUT/$ST/virulence" "$BASE_OUT/$ST/resfinder"

    echo "Genomes: $total"
    printf '%s\n' "${files[@]}" | parallel --bar -j "$JOBS" process_vf {} "$ST"
    printf '%s\n' "${files[@]}" | parallel --bar -j "$JOBS" process_rf {} "$ST"
done

echo "======================================"
echo "DONE ALL ST TYPES"
echo "======================================"
for ST in "${ST_LIST[@]}"; do
    VF_COUNT=$(ls "$BASE_OUT/$ST/virulence"/*_VF.tsv 2>/dev/null | wc -l | tr -d ' ')
    RF_COUNT=$(ls "$BASE_OUT/$ST/resfinder"/*_RF.txt 2>/dev/null | wc -l | tr -d ' ')
    echo "$ST: VF=$VF_COUNT RF=$RF_COUNT"
done
