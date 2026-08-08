# Changelog

## [2.0.0] - 2026

### Changed
- Rebuilt repository as an end-to-end reproducible pipeline.
- Added full input stages: metadata filtering, assembly download, abricate
  (VFDB/CARD), VirulenceFinder/ResFinder, finder-summary matrices, pangenome +
  IQ-TREE core tree.
- Replaced the generic placeholder analysis modules with the canonical numbered
  scripts used to generate every manuscript result and figure (`scripts/`).
- Configuration is now environment-driven (`config.R`, `config/pipeline_config.sh`).
- `run_all.sh` orchestrates the whole workflow; `run_pipeline.R` drives the
  analysis/figures stage.
- Removed hardcoded NCBI API key from the downloader (env var only).
- Added `environment.yml`, `install.R`, expanded docs and `.gitignore`.

## [1.0.0] - 2025-03-17

### Added
- Initial release of the modular analysis pipeline.
- AMR (CARD), VFDB and plasmid analysis modules.
- Metadata matching system, visualisations, statistical framework.
