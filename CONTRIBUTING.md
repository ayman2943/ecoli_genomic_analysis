# Contributing to E. coli Genomic Analysis Pipeline

Thank you for considering contributing.

## How to contribute

1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/...`).
3. Make changes, following the existing conventions (numbered scripts,
   env-driven config in `config.R`, no hardcoded paths or secrets).
4. Verify your R scripts parse: `Rscript -e 'parse("scripts/<file>.R")'` and
   shell scripts with `bash -n`.
5. Open a pull request describing the change and the reproduction step.

## Guidelines

- **No secrets** — never commit NCBI API keys or database paths. Use the
  `NCBI_API_KEY` environment variable / `.env`.
- **Keep it reproducible** — inputs should be derivable from EnteroBase
  metadata; outputs regenerable with `run_all.sh`.
- **Document** — update README/QUICKSTART/docs if you change the workflow.
- **Match thresholds** — the pipeline uses 80% identity / 80% coverage unless
  explicitly stated otherwise.

## Reporting bugs

Include the exact command, the failing script name, and the error output
(prefer `Rscript run_pipeline.R ST69 --from <script>` traces).
