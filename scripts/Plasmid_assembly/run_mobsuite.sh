#!/bin/bash
# =============================================================================
# scripts/Plasmid_assembly/run_mobsuite.sh
#
# Run mob-suite (mob_recon) over all downloaded assemblies for one or all STs
# to reconstruct/type plasmid contigs, then concatenate each genome's plasmid
# contigs into one combined per-genome FASTA for downstream plasmid abricate
# annotation (see run_abricate_plasmids.sh).
#
# Naming matches the existing downstream analysis code exactly: each genome's
# combined plasmid FASTA is named "{base}_combined_plasmids.fasta", where
# {base} is the same genome identifier used everywhere else in this pipeline
# (the .fna filename without extension, e.g. "Escherichia_coli_GCF_XXXXXXXXX").
# scripts/Analysis/Figure8.R and scripts/Virulence/Plasmid_VFDB.R's
# clean_gid() helper strips the "_combined_plasmids_vfdb.tsv" /
# "_combined_plasmids_card.tsv" suffix abricate will add to this filename in
# the next step -- so this naming convention is load-bearing, not cosmetic.
#
# Usage:
#   bash scripts/Plasmid_assembly/run_mobsuite.sh ST69          # single ST
#   bash scripts/Plasmid_assembly/run_mobsuite.sh all           # all 5 STs (default)
#   JOBS=6 bash scripts/Plasmid_assembly/run_mobsuite.sh ST69   # parallel jobs
#
# Requires: mob_suite (mob_recon) -- see environment.yml / install docs at
#   https://github.com/phac-nml/mob-suite
# Input:  {ST}/*.fna                              (same assemblies abricate uses)
# Output: ../Plasmid/mob_suite_raw/{ST}/{base}/    (raw mob_recon output per genome)
#         ../Plasmid/combined_fasta/{ST}/{base}_combined_plasmids.fasta
#           (only written for genomes where mob_recon found >=1 plasmid contig;
#           genomes with no predicted plasmids are skipped and logged)
#
# NOTE ON PATHS: like the existing "Plasmid/" dataset referenced from
# scripts/Virulence/Plasmid_VFDB.R and scripts/Analysis/Figure8.R
# (config$BASE_DIR/../Plasmid/...), this script's output is written to a
# sibling "Plasmid/" directory NEXT TO the repo root, not inside it -- it is
# a large generated dataset, not code, so it is not committed to git.
# =============================================================================
set -euo pipefail

JOBS="${JOBS:-4}"
ST_ARG="${1:-all}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PLASMID_ROOT="$(dirname "$REPO_ROOT")/Plasmid"

if [ "$ST_ARG" = "all" ]; then ST_LIST=("ST10" "ST69" "ST73" "ST95" "ST131"); else ST_LIST=("$ST_ARG"); fi

check_mobsuite() {
  if ! command -v mob_recon >/dev/null 2>&1; then
    echo "ERROR: mob_recon not found. Install mob_suite first:"
    echo "  conda install -c bioconda mob_suite"
    echo "  (or: pip install mob_suite -- see https://github.com/phac-nml/mob-suite)"
    exit 1
  fi
}
check_mobsuite

# Threads PER mob_recon call. We already parallelize ACROSS genomes with GNU
# parallel below (JOBS concurrent genomes) -- keep this at 1 so JOBS genomes
# in flight doesn't multiply out to JOBS x MOB_THREADS cores requested at
# once. Override explicitly (e.g. MOB_THREADS=2) if you lower JOBS instead.
MOB_THREADS="${MOB_THREADS:-1}"

process_genome() {
  local st="$1" file="$2" raw_dir="$3" combined_dir="$4"
  local base
  base=$(basename "$file" .fna)
  local out_dir="$raw_dir/$base"

  if [ -s "$combined_dir/${base}_combined_plasmids.fasta" ]; then
    return  # already done
  fi

  local log_file="$raw_dir/${base}.mob_recon.log"
  # mob_recon reconstructs and MOB-types plasmid contigs from a single
  # assembly (typing runs automatically -- there is no separate --run_typer
  # flag). -f/--force lets it overwrite out_dir on a re-run; mob_recon
  # creates out_dir itself.
  if ! mob_recon --infile "$file" --outdir "$out_dir" \
        --num_threads "$MOB_THREADS" --force \
        >"$log_file" 2>&1; then
    echo "  [WARN] mob_recon failed for $base (see $log_file)"
    return
  fi

  # mob_recon writes one FASTA per predicted plasmid: plasmid_AA000001.fasta,
  # plasmid_AA000002.fasta, ... (chromosome.fasta is also written; excluded).
  shopt -s nullglob
  local plasmid_fastas=("$out_dir"/plasmid_*.fasta)
  shopt -u nullglob
  if [ "${#plasmid_fastas[@]}" -eq 0 ]; then
    echo "  [SKIP] $base : no plasmids predicted"
    return
  fi

  cat "${plasmid_fastas[@]}" > "$combined_dir/${base}_combined_plasmids.fasta"
  echo "  [OK] $base : ${#plasmid_fastas[@]} plasmid contig(s)"
}
export -f process_genome

process_st() {
  local st="$1"
  local input_dir="$REPO_ROOT/$st"
  if [ ! -d "$input_dir" ] || [ -z "$(ls "$input_dir"/*.fna 2>/dev/null)" ]; then
    echo "[SKIP] $st : no .fna files in $input_dir"
    return
  fi

  echo "=== $st : mob-suite plasmid reconstruction ==="
  local raw_dir="$PLASMID_ROOT/mob_suite_raw/$st"
  local combined_dir="$PLASMID_ROOT/combined_fasta/$st"
  mkdir -p "$raw_dir" "$combined_dir"

  local files=("$input_dir"/*.fna)
  printf '%s\n' "${files[@]}" | \
    parallel --bar --jobs "$JOBS" --halt soon,fail=1 \
      process_genome "$st" {} "$raw_dir" "$combined_dir"

  local n_with_plasmids
  n_with_plasmids=$(ls "$combined_dir"/*_combined_plasmids.fasta 2>/dev/null | wc -l | tr -d ' ')
  echo "  DONE $st : ${n_with_plasmids}/${#files[@]} genomes had >=1 predicted plasmid"
}

for st in "${ST_LIST[@]}"; do
  process_st "$st"
done

echo "=== mob-suite finished. Combined plasmid FASTAs in $PLASMID_ROOT/combined_fasta/{ST}/ ==="
echo "Next: bash scripts/Plasmid_assembly/run_abricate_plasmids.sh $ST_ARG"
