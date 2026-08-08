# Usage guide

Detailed instructions for each pipeline stage. Run everything from the repo
root inside the `wgs` conda environment.

## Stage 00 — Metadata filtering

```bash
export RAW_ENTERO_EXPORT=/path/to/EnteroBase_export.xlsx
Rscript 00_metadata/filter_metadata.R "$RAW_ENTERO_EXPORT"
```

Produces `metadata/{TARGET_ST}_filtered.xlsx` and
`metadata_matched/matched_{TARGET_ST}.xlsx` containing only rows whose ST
matches the target lineage. Column requirements: `Uberstrain` or `Name`, plus
an ST column (`ST` / `Sequence Type` / `MLST`).

## Stage 01 — Download assemblies

```bash
python3 01_download/download_assemblies.py metadata/ST69_filtered.xlsx ST69
gzip -d ST69/*.fna.gz
```

- Resolves each row to an assembly via GCF_/GCA_ accession, BioSample, or SRA
  (in that order).
- Downloads by FTP (fast) with `datasets` fallback.
- Requires the NCBI `datasets` CLI and Entrez tools (`esearch`/`elink`/
  `efetch`/`xtract`) — both in `environment.yml`.
- Set `DOWNLOAD_ASSUME_YES=1` to skip the interactive confirmation.
- Without `NCBI_API_KEY` the downloader runs at 3 req/sec (slower).

## Stage 02 — Annotation

### abricate (VFDB + CARD)

```bash
bash 02_annotation/run_abricate.sh ST69        # or "all" for the 5 STs
```

Writes per-genome tables to `card_vfdb_result/{vfdb,card}/{ST}/` and
per-ST summaries to `card_vfdb_result/{vfdb,card}_summary/{ST}*.tsv`.
Thresholds: 80% identity, 80% coverage. Resumable (skips existing outputs).

### VirulenceFinder + ResFinder

```bash
bash 02_annotation/run_vf_resfinder.sh ST69
```

Writes `analysis_results/{ST}/virulence/{genome}_VF.tsv` and
`analysis_results/{ST}/resfinder/{genome}_RF.txt`. Databases default to
`$HOME/databases/{virulencefinder_db,resfinder_db}` (override via `VF_DB`,
`RF_DB`).

### Build summary matrices

```bash
Rscript 02_annotation/build_finder_summaries.R analysis_results finder_result
```

Parses every VF/ResFinder output (under `analysis_results/{ST}/`) into:
- `virulencefinder_summary/{binary_matrix,burden,gene_frequency,long}.tsv`
- `resfinder_summary/{binary_matrix,burden,gene_frequency,long}.tsv`
- `logs/summary_qc.tsv`

These are the matrices consumed by every downstream analysis script.

## Stage 03 — Pangenome + core tree

```bash
bash 03_pangenome/run_pangenome_tree.sh ST69            # all genomes
bash 03_pangenome/run_pangenome_tree.sh ST10 500        # stratified subsample
```

1. Prokka annotation (`annotation/annotations_{ST}/`)
2. PPanGGOLiN workflow (`pangenome_output/ppanggolin_raw/`)
3. Core-genome MSA export + trimAL
4. IQ-TREE (`-m MFP -bb 1000`) → `{ST}_bootstrap.treefile` (copied to root)
5. Pangenome table export → `pangenome_output/ppanggolin_output/`

RAM: ≥32 GB recommended for >1,000 genomes.

## Stage 04 — R analysis

```bash
Rscript run_pipeline.R ST69            # full
Rscript run_pipeline.R ST69 --quick    # clustering + tree mapping only
Rscript run_pipeline.R ST69 --from 09e # resume at a given script
Rscript run_pipeline.R ST69 --skip 15 --skip 16   # skip tree figures
```

The runner sources `config.R` (paths from env, see README) and executes the
numbered scripts in `scripts/` in dependency order. Optional scripts (tree
figures requiring `ggtree`) are skipped gracefully if the package is absent.

### Outputs

- `output/{ST}/virulencefinder_validation/` — VF-based results
- `output/{ST}/vfdb_analysis/` — VFDB-based results
- `output/figures_*/` — manuscript figure panels
- `output/combined_figures/All_Figures.pdf` — combined multi-page PDF (script `99`)
- `output/{ST}/...` — cluster assignments, K-type tables, enrichment results

## Resuming / partial runs

Everything is resumable:

- Download: existing `.fna.gz` are skipped.
- abricate / VF / ResFinder: existing per-genome outputs are skipped.
- R analysis: use `--from <script>` to resume, or `--skip` to omit specific
  scripts (e.g. tree figures).
