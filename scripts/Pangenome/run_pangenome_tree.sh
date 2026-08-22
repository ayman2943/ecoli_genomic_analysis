#!/usr/bin/env bash
# =============================================================================
# 03_pangenome/run_pangenome_tree.sh
#
# Build a pangenome + core-genome tree for a given ST from downloaded
# assemblies: Prokka → PPanGGOLiN → core MSA → trimAL → IQ-TREE.
#
# Usage:
#   bash 03_pangenome/run_pangenome_tree.sh ST69            # all downloaded genomes
#   bash 03_pangenome/run_pangenome_tree.sh ST10 500        # stratified subsample N
#   bash 03_pangenome/run_pangenome_tree.sh ST10 500 --prokka-only   # stop early
#
# Requires: conda env 'wgs' (prokka, ppanggolin, trimal, iqtree)
# RAM note: PPanGGOLiN on >1000 genomes needs ≥32 GB RAM.
# Outputs:
#   pangenome_output/ppanggolin_output/   (gene_presence_absence.Rtab, partitions/)
#   {ST}_bootstrap.treefile               (copied to repo root for the analysis)
# =============================================================================
set -euo pipefail

ST="${1:-ST69}"
N_GENOMES="${2:-0}"          # 0 = use ALL genomes for the ST
MODE="${3:-full}"            # full | prokka-only

BASE_DIR="$PWD"
RAW_DIR="${BASE_DIR}/${ST}"
ANNOT_DIR="${BASE_DIR}/annotation/annotations_${ST}"
PANG_DIR="${BASE_DIR}/pangenome_output"
PP_OUT="${PANG_DIR}/ppanggolin_output"
CPU="${CPU:-8}"

if [ ! -d "$RAW_DIR" ] || [ -z "$(ls "$RAW_DIR"/*.fna 2>/dev/null)" ]; then
  echo "ERROR: no .fna files in $RAW_DIR. Run 01_download first."
  exit 1
fi
mkdir -p "$ANNOT_DIR" "$PANG_DIR" "$PP_OUT"

echo "============================================"
echo " ST:        $ST"
echo " Genomes:   $([ "$N_GENOMES" = "0" ] && echo ALL || echo "$N_GENOMES subsample")"
echo " Annotation: $ANNOT_DIR"
echo " Pangenome: $PANG_DIR"
echo " Mode:      $MODE"
echo "============================================"

# ---- 1. Select genomes ----
if [ "$N_GENOMES" != "0" ]; then
  SUBSAMPLE_LIST="${BASE_DIR}/selected_${ST}_${N_GENOMES}.txt"
  SUBSAMPLE_DIR="${BASE_DIR}/${ST}_subsampled_${N_GENOMES}"
  mkdir -p "$SUBSAMPLE_DIR"

  if [ ! -f "$SUBSAMPLE_LIST" ]; then
    echo "[1/6] Stratified subsampling ${N_GENOMES} genomes for ${ST}..."
    Rscript -e '
      suppressPackageStartupMessages({library(readxl); library(tidyverse)})
      args <- commandArgs(TRUE)
      st <- args[1]; n_target <- as.integer(args[2])
      meta <- read_excel(file.path("metadata_matched", paste0("matched_", st, ".xlsx")))
      raw_ids <- sub("\\.fna$", "", list.files(args[3], pattern = "\\.fna$"))
      df <- meta %>%
        filter(.data[[intersect(c("Uberstrain","Name"), colnames(.))[1]]] %in% raw_ids) %>%
        rename(genome_id = intersect(c("Uberstrain","Name"), colnames(.))[1])
      if (nrow(df) < n_target) { selected <- df
      } else {
        year_c <- intersect(c("Collection Year","Year"), colnames(df))[1]
        niche_c <- intersect(c("Source Niche","Source Type"), colnames(df))[1]
        df$year_bin <- if (is.na(year_c)) "unknown" else cut(as.numeric(df[[year_c]]),
          breaks=c(-Inf,2010,2015,2020,Inf), labels=c("pre2010","2010-14","2015-19","2020+"))
        df$niche <- if (is.na(niche_c)) "unknown" else ifelse(grepl("Human|clinical|blood|urine",
          df[[niche_c]], ignore.case=TRUE), "clinical", "other")
        df$stratum <- paste(df$year_bin, df$niche, sep="|")
        cnt <- count(df, stratum) %>% mutate(target = pmax(round(n/sum(n)*n_target), 1))
        set.seed(42)
        selected <- df %>% group_by(stratum) %>%
          slice_sample(n = cnt$target[match(stratum, cnt$stratum)]) %>% ungroup()
      }
      write_lines(paste0(selected$genome_id, ".fna"), args[4])
      cat("Selected:", nrow(selected), "\n")
    ' "$ST" "$N_GENOMES" "$RAW_DIR" "$SUBSAMPLE_LIST"
  fi

  if [ "$(ls "$SUBSAMPLE_DIR"/*.fna 2>/dev/null | wc -l | tr -d ' ')" = "0" ]; then
    while IFS= read -r fname; do cp "${RAW_DIR}/${fname}" "${SUBSAMPLE_DIR}/" 2>/dev/null || true; done < "$SUBSAMPLE_LIST"
  fi
  WORK_DIR="$SUBSAMPLE_DIR"
else
  WORK_DIR="$RAW_DIR"
fi

TOTAL=$(ls "$WORK_DIR"/*.fna 2>/dev/null | wc -l | tr -d ' ')
echo "  Genomes ready: $TOTAL"

# ---- 2. Prokka ----
echo "[2/6] Prokka annotation..."
COUNT=0
for fasta_file in "$WORK_DIR"/*.fna; do
  COUNT=$((COUNT + 1))
  base=$(basename "${fasta_file}" .fna)
  OUTDIR="${ANNOT_DIR}/${base}"
  if [ -f "${OUTDIR}/${base}.gff" ]; then
    echo "  [${COUNT}/${TOTAL}] SKIP ${base}"
    continue
  fi
  echo "  [${COUNT}/${TOTAL}] Annotating ${base}..."
  rm -rf "${OUTDIR}"
  prokka --outdir "${OUTDIR}" --prefix "${base}" --locustag "${base}" \
    --genus Escherichia --species coli --strain "${base}" --kingdom Bacteria \
    --gcode 11 --usegenus --cpus "$CPU" --force "${fasta_file}"
done
echo "  Done."

if [ "$MODE" = "prokka-only" ]; then echo "Stopping after annotation (--prokka-only)."; exit 0; fi

# ---- 3. GFF list ----
GFF_LIST="${ANNOT_DIR}/gff_list_${ST}.tsv"
echo "[3/6] Creating GFF list..."
> "$GFF_LIST"
for d in "${ANNOT_DIR}"/*/; do
  [ -f "${d}$(basename "${d%/}").gff" ] && echo -e "${d}$(basename "${d%/}").gff\t$(basename "${d%/}")" >> "$GFF_LIST"
done
echo "  GFF entries: $(wc -l < "$GFF_LIST")"

# ---- 4. PPanGGOLiN ----
echo "[4/6] Running PPanGGOLiN..."
ppanggolin workflow --anno "$GFF_LIST" --output "$PANG_DIR/ppanggolin_raw" --cpu "$CPU" --ram "${PPANG_RAM:-16}"
echo "  Done."

# ---- 5. Core MSA + trimAL + IQ-TREE ----
echo "[5/6] Core MSA + trimming + tree..."
ppanggolin write_pangenome -p "$PANG_DIR/ppanggolin_raw/pangenome.h5" --fasta "${ST}_core_genome"
ALN=$(ls "${PANG_DIR}/${ST}_core_genome".* 2>/dev/null | head -1)
if [ -n "$ALN" ]; then
  trimal -automated1 -in "$ALN" -out "${PANG_DIR}/${ST}_trimmed.aln"
  iqtree -s "${PANG_DIR}/${ST}_trimmed.aln" -m MFP -bb 1000 -nt AUTO -mem "${IQTREE_MEM:-14G}" -pre "${PANG_DIR}/${ST}_bootstrap"
  cp "${PANG_DIR}/${ST}_bootstrap.treefile" "${BASE_DIR}/" 2>/dev/null || true
  echo "  Tree: ${ST}_bootstrap.treefile"
fi

# ---- 6. Export pangenome tables ----
echo "[6/6] Exporting pangenome tables..."
ppanggolin write -p "$PANG_DIR/ppanggolin_raw/pangenome.h5" --output "$PP_OUT" \
  --Rtab --csv --regions --spots --families_tsv

echo ""
echo "============================================"
echo " PRE-PROCESSING DONE for ${ST}"
echo "============================================"
echo "Next step: TARGET_ST=${ST} Rscript run_pipeline.R"
