# Quick start

End-to-end reproduction of the manuscript analyses for a target lineage
(default `ST69`).

## 1. Setup (one time)

```bash
# Conda environment with all bioinformatics tools
conda env create -f environment.yml
conda activate wgs

# R packages
Rscript install.R
```

## 2. Configuration

Edit `config/pipeline_config.sh`:

- `RAW_ENTERO_EXPORT` → path to your EnteroBase export `.xlsx`
- `VF_DB` / `RF_DB` → paths to the VirulenceFinder / ResFinder databases
- `JOBS`, `CPU`, `PPANG_RAM`, `IQTREE_MEM` → compute settings

Optional: `export NCBI_API_KEY=...` (put it in `.env` or your shell; never
commit it).

## 3. Full run

```bash
bash run_all.sh ST69          # everything: download → tools → analysis → figures
```

## 4. Stage-by-stage

```bash
# Metadata + download
bash run_all.sh ST69 --download

# Annotation + summary matrices
bash run_all.sh ST69 --annotate

# Pangenome + core tree
bash run_all.sh ST69 --pangenome

# R analysis + figures
bash run_all.sh ST69 --analyze
```

## 5. Verify outputs

```bash
# Check summary matrices were built
head finder_result/virulencefinder_summary/virulencefinder_binary_matrix.tsv
# Check figures (3-panel panels + combined multi-page PDF)
ls output/figures_all_3panel/
ls output/combined_figures/
```

## Troubleshooting

- `abricate: command not found` → `conda activate wgs`
- PPanGGOLiN OOM → lower `PPANG_RAM`, or subsample with
  `bash 03_pangenome/run_pangenome_tree.sh ST69 500`
- `ggtree` missing → tree figures skipped (optional); `Rscript install.R`
- Downloader interactive prompt → `DOWNLOAD_ASSUME_YES=1`
