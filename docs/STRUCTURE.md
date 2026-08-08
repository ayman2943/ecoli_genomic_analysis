# Repository structure

Layout of `Ecoli_genomic_analysis`.

```
Ecoli_genomic_analysis/
├── 00_metadata/
│   └── filter_metadata.R          Filter EnteroBase export → metadata/{ST}_filtered.xlsx
│                                  and metadata_matched/matched_{ST}.xlsx
├── 01_download/
│   └── download_assemblies.py     Parallel NCBI assembly download (FTP + datasets
│                                  fallback). API key via NCBI_API_KEY env/.env only.
├── 02_annotation/
│   ├── run_abricate.sh            abricate VFDB + CARD (80% id / 80% cov) + summaries
│   ├── run_vf_resfinder.sh        VirulenceFinder + ResFinder (80/80), resumable
│   └── build_finder_summaries.R   → finder_result/{virulencefinder,resfinder}_summary/*.tsv
├── 03_pangenome/
│   └── run_pangenome_tree.sh      Prokka → PPanGGOLiN → core MSA → trimAL → IQ-TREE
├── config/
│   ├── pipeline_config.sh         Site-specific settings (paths, threads, memory)
│   └── config.R                   R config (env-driven), kept in sync with repo-root copy
├── scripts/                       Canonical numbered R analysis scripts (manuscript code)
│   ├── 00_config.R                Original dev config (kept for reference)
│   ├── 02b/02c                    Shell-gene clustering (VFDB / VF)
│   ├── 03c/03d                    Tree mapping (VFDB / VF)
│   ├── 04a/04b                    Virulence analysis, cluster gene analysis
│   ├── 05a/05b/05c                Resistance analysis (CARD + ResFinder)
│   ├── 06a                        Temporal trends + Mann-Kendall
│   ├── 07a/07b                    ST10 decomposition, composition drivers
│   ├── 08                         Clinical enrichment (3-panel)
│   ├── 09a–09e                    Sensitivity, capsule classification/comparison,
│   │                              K-type analysis, summary stats
│   ├── 10a–10g                    Validation (KPS, RGP, plasmid, allelic, ARG, long-read)
│   ├── 11–22                      All manuscript figure scripts
│   └── 99                         ALL_FIGURES_COMBINED
├── config.R                       R configuration (used by all scripts via source())
├── run_pipeline.R                 Master R runner (ordered script execution)
├── run_all.sh                     End-to-end bash orchestrator
├── environment.yml                Conda environment (bioinformatics + R)
├── requirements.txt               Python deps for downloader
├── install.R                      R package installer (CRAN + ggtree)
├── README.md                      Overview + usage
├── QUICKSTART.md                  5-step quick start
├── LICENSE                        MIT
└── docs/
    ├── STRUCTURE.md               This file
    └── USAGE.md                   Detailed usage of each stage
```

## Runtime data (gitignored, regenerable)

| Path | Produced by | Consumed by |
|------|-------------|-------------|
| `metadata_matched/matched_{ST}.xlsx` | `00_metadata/` | `config.R` (`METADATA_FILE`) |
| `{ST}/Escherichia_coli_*.fna` | `01_download/` | abricate, VF, ResFinder |
| `card_vfdb_result/` | `02_annotation/run_abricate.sh` | `config.R` (`CARD_VFDB_DIR`) |
| `analysis_results/{ST}/` | `02_annotation/run_vf_resfinder.sh` | `build_finder_summaries.R` |
| `finder_result/` | `build_finder_summaries.R` | `config.R` (`INPUT_DIR`) |
| `pangenome_output/` | `03_pangenome/` | `config.R` (`PANGENOME_DIR`) |
| `{ST}_bootstrap.treefile` | `03_pangenome/` | `config.R` (`TREE_FILE`) |
| `output/{ST}/` | `run_pipeline.R` | manuscript tables (virulence, resistance, K-type, clusters) |
| `output/figures_*/` | scripts `11`–`22` | manuscript figure panels |
| `output/combined_figures/` | script `99` | combined multi-page `All_Figures.pdf` |
