# E. coli ExPEC Genomic Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

End-to-end pipeline for the population genomics of extraintestinal pathogenic
*Escherichia coli* (ExPEC) across five sequence types (ST10, ST69, ST73, ST95,
ST131). Reproduces every analysis and figure in the accompanying manuscript.

```
EnteroBase metadata ──▶ download assemblies ──▶ abricate VFDB/CARD
    + VirulenceFinder + ResFinder ──▶ summary matrices ──▶ pangenome + core tree
    ──▶ shell-gene clustering ──▶ tree mapping ──▶ virulence/resistance analysis
    ──▶ temporal trends (Mann-Kendall) ──▶ clinical enrichment ──▶ figures
```

## Table of contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Pipeline stages](#pipeline-stages)
- [Directory layout](#directory-layout)
- [Reproducing the manuscript](#reproducing-the-manuscript)
- [Configuration](#configuration)
- [Outputs](#outputs)
- [FAQ / troubleshooting](#faq--troubleshooting)
- [License](#license)

## Overview

This pipeline analyses a genome collection of ExPEC lineages to dissect
temporal trends, cluster structure, clinical enrichment, capsule (K-type)
diversity, and antimicrobial-resistance/virulence content. All thresholds are
80% identity / 80% coverage unless stated otherwise.

Key analyses:

- **Shell-gene clustering** — PPanGGOLiN shell genes projected onto VFDB and
  VirulenceFinder content; cluster assignment per genome.
- **Tree mapping** — cluster labels mapped onto the core-genome
  phylogeny (`{ST}_bootstrap.treefile` from IQ-TREE).
- **Virulence & resistance** — per-cluster enrichment, temporal dynamics,
  Mann-Kendall trend tests for VFDB, VirulenceFinder, CARD and ResFinder.
- **ST10 decomposition** — compositional drivers of expansion within ST10.
- **Clinical enrichment** — clinical vs non-clinical comparisons with
  sensitivity analysis.
- **Capsule (K-type) classification** — group 2/group 3 capsule assignment
  based on `kpsM` allele classes from the VirulenceFinder binary matrix.
- **Validation** — long-read assemblies, KPS operon structure, RGP context,
  plasmid context.
- **Figures** — all manuscript figures under `output/` (see [Outputs](#outputs)).

## Installation

### 1. Conda environment

```bash
conda env create -f environment.yml
conda activate wgs
```

This provides abricate, prokka, ppanggolin, iqtree, trimal, seqtk,
VirulenceFinder, ResFinder and the NCBI `datasets`/Entrez tools.

### 2. Databases

Download the reference databases (all FASTA-based, no config needed):

```bash
mkdir -p ~/databases
cd ~/databases
# VirulenceFinder
git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db.git
# ResFinder
git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
# abricate databases are installed with the package
abricate --setupdb
```

Set `VF_DB` / `RF_DB` in `config/pipeline_config.sh` if you placed them
elsewhere.

### 3. R packages

```bash
Rscript install.R
```

Installs all CRAN packages plus `ggtree` (Bioconductor, used only by the two
optional tree-figure scripts).

## Quick start

1. Place your EnteroBase export (`RAW_ENTERO_EXPORT`) where the pipeline can
   read it and set the variable in `config/pipeline_config.sh`.
2. Optionally export `NCBI_API_KEY` to speed up downloads.

```bash
bash run_all.sh ST69                  # full run for ST69
bash run_all.sh ST69 --annotate       # skip download, just run tools + summaries
bash run_all.sh ST69 --analyze        # R analysis + figures only
```

Or run stages manually:

```bash
# 01 filter metadata
Rscript 00_metadata/filter_metadata.R "$RAW_ENTERO_EXPORT"
# 02 download (requires NCBI datasets + Enterz tools)
python3 01_download/download_assemblies.py metadata/ST69_filtered.xlsx ST69
gzip -d ST69/*.fna.gz
# 03-04 annotation
bash 02_annotation/run_abricate.sh ST69
bash 02_annotation/run_vf_resfinder.sh ST69
# 05 summary matrices
Rscript 02_annotation/build_finder_summaries.R analysis_results finder_result
# 06 pangenome + tree
bash 03_pangenome/run_pangenome_tree.sh ST69
# 07 analysis + figures
Rscript run_pipeline.R ST69
```

## Pipeline stages

| # | Script | Purpose |
|---|--------|---------|
| 00 | `00_metadata/filter_metadata.R` | Filter EnteroBase export to target ST |
| 01 | `01_download/download_assemblies.py` | Parallel NCBI download of assemblies |
| 02 | `02_annotation/run_abricate.sh` | abricate VFDB + CARD (80/80) |
| 03 | `02_annotation/run_vf_resfinder.sh` | VirulenceFinder + ResFinder (80/80) |
| 04 | `02_annotation/build_finder_summaries.R` | binary/burden/frequency matrices |
| 05 | `03_pangenome/run_pangenome_tree.sh` | Prokka → PPanGGOLiN → MSA → trimAL → IQ-TREE |
| 06 | `run_pipeline.R` | All analyses + figures (see below) |

`run_pipeline.R` runs the numbered scripts in `scripts/` in dependency order:

- `02b`–`02c` shell clustering (VFDB / VF)
- `03c`–`03d` tree mapping (VFDB / VF)
- `04a`–`04b` virulence analysis + cluster genes
- `05a`–`05c` resistance analysis (CARD + ResFinder)
- `06a` temporal trends + Mann-Kendall
- `07a`–`07b` ST10 decomposition
- `08` clinical enrichment
- `09a`–`09e` sensitivity, capsule classification/comparison, K-types, summary stats
- `10a`–`10g` validation (KPS, RGP, plasmid, allelic conversion, long reads)
- `11`–`22`, `99` figures (all manuscript figures + supplementary)

## Directory layout

```
Ecoli_genomic_analysis/
├── 00_metadata/        metadata filtering
├── 01_download/        assembly downloader (sanitized, API key via env)
├── 02_annotation/      abricate, VirulenceFinder/ResFinder, summary builder
├── 03_pangenome/       Prokka / PPanGGOLiN / IQ-TREE pipeline
├── config/             pipeline_config.sh (site settings)
├── scripts/            canonical numbered R analysis scripts (manuscript code)
├── config.R            R configuration (env-driven, see below)
├── run_pipeline.R      master R runner
├── run_all.sh          end-to-end orchestrator
├── environment.yml     conda environment
├── requirements.txt    python deps
├── install.R           R package installer
└── docs/               detailed structure & usage docs
```

## Reproducing the manuscript

All results in the manuscript (tables + figures) come from the R scripts in
`scripts/`. Running `Rscript run_pipeline.R ST69` from the repo root
regenerates everything under `output/`. The three supplementary workbook
sheets are produced by scripts `09d` (K-type assignment) and `09e`
(summary statistics); the accession sheet is built from
`metadata_matched/matched_{ST}.xlsx`.

## Configuration

All settings are environment variables with defaults in `config.R`:

| Variable | Default | Purpose |
|----------|---------|---------|
| `TARGET_ST` | `ST69` | Lineage to analyse |
| `ECOLI_BASE_DIR` | `getwd()` | Working directory containing inputs |
| `ECOLI_PANGENOME_DIR` | — | External pangenome dir for non-ST69 lineages |
| `NCBI_API_KEY` | — | NCBI key to speed downloads (never committed) |
| `JOBS` | 6 | parallel abricate jobs |
| `CPU` | 8 | prokka / ppanggolin threads |
| `PPANG_RAM` | 16 | PPanGGOLiN RAM (GB) |
| `IQTREE_MEM` | 14G | IQ-TREE memory |

## Outputs

- `finder_result/` — VirulenceFinder/ResFinder binary, burden, gene-frequency,
  long-format matrices + QC (flat, all-ST, `st` column).
- `card_vfdb_result/` — abricate VFDB/CARD per-genome tables and summaries.
- `output/{ST}/` — per-analysis tables (virulence, resistance, temporal,
  enrichment, K-type, clusters).
- `output/figures_*/` — manuscript figure panels; `output/combined_figures/`
  holds the combined multi-page `All_Figures.pdf` (script `99`).
- `{ST}_bootstrap.treefile` — core-genome phylogeny from IQ-TREE.

## FAQ / troubleshooting

- **`ggtree` not installed** — tree figure scripts (`15_`, `16_`) are marked
  optional and skipped automatically; run `Rscript install.R` to add ggtree.
- **Downloader asks for confirmation** — set `DOWNLOAD_ASSUME_YES=1` for
  unattended runs.
- **PPanGGOLiN out of memory** — reduce `N_GENOMES` (subsample) or lower
  `PPANG_RAM`; ≥32 GB recommended for >1,000 genomes.
- **No ST directories found** — the summary builder (`build_finder_summaries.R`)
  expects `analysis_results/{ST}/virulence|resfinder/`; run
  `run_vf_resfinder.sh` first.

## License

MIT — see `LICENSE`.

## Citation

If you use this pipeline in your work, please cite the associated manuscript
(after publication) and link to this repository.
