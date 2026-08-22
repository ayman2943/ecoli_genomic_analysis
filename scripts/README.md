# scripts/ — the complete, reproducible pipeline

Every script needed to go from a raw EnteroBase metadata export to every
figure and table in the manuscript lives here, in one place, organised by
pipeline stage. Nothing outside `scripts/` (plus `config.R` and
`run_pipeline.R`/`run_all.sh` at the repo root) is required to reproduce the
publication.

`run_all.sh` (repo root) runs the bash/python stages below in order, then
calls `run_pipeline.R`, which runs the R stages in order. See the repo-root
`README.md` for the one-command quick start.

## Folder layout, in pipeline order

| Folder | Stage | Runs via |
|---|---|---|
| `Setup/` | Filter EnteroBase metadata to the target ST; download assemblies from NCBI | `run_all.sh` |
| `Annotation/` | abricate (VFDB + CARD) and VirulenceFinder/ResFinder on whole-genome assemblies | `run_all.sh` |
| `Plasmid_assembly/` | mob-suite plasmid reconstruction, then abricate (VFDB + CARD) on the reconstructed plasmid contigs | `run_all.sh` (optional stage) |
| `Pangenome/` | Prokka -> PPanGGOLiN -> core-genome MSA -> trimAL -> IQ-TREE | `run_all.sh` |
| `Prerequisites/` | Finder-summary matrices, shell-gene clustering, tree-cluster mapping, per-cluster gene stats | `run_pipeline.R` |
| `Virulence/` | Raw per-database virulence gene burden/annotation: VFDB, VirulenceFinder, plasmid VFDB | `run_pipeline.R` |
| `AMR/` | Raw per-database AMR gene burden/annotation: CARD, ResFinder, plasmid CARD | `run_pipeline.R` |
| `Analysis/` | One script per published figure (`Figure1.R`...`Figure8.R`), plus `Supplementary.R` | `run_pipeline.R` |

## Setup/

- `filter_metadata.R` — filters an EnteroBase export down to the target ST,
  writing `metadata/{ST}_filtered.xlsx` and `metadata_matched/matched_{ST}.xlsx`.
- `download_assemblies.py` — downloads the matched assemblies from NCBI
  (FTP, with `datasets` CLI fallback).

## Annotation/

- `run_abricate.sh` — abricate VFDB + CARD (80% identity / 80% coverage) on
  every whole-genome assembly, plus the per-ST summary tables
  (`card_vfdb_result/{vfdb,card}_summary/`) that `Prerequisites/` and `AMR/CARD.R`
  read.
- `run_vf_resfinder.sh` — VirulenceFinder + ResFinder (80/80) on every
  assembly.

## Plasmid_assembly/

- `run_mobsuite.sh` — runs `mob_recon` (mob-suite) on every whole-genome
  assembly to reconstruct and type plasmid contigs, then concatenates each
  genome's plasmid contigs into one combined FASTA
  (`{genome}_combined_plasmids.fasta`) for the next step.
- `run_abricate_plasmids.sh` — runs abricate (VFDB + CARD, same 80/80
  thresholds) on those combined plasmid FASTAs, then builds the per-ST
  summary tables that `Virulence/Plasmid_VFDB.R`, `AMR/Plasmid_CARD.R`, and
  `Analysis/Figure8.R` read directly:
  `../Plasmid/plasmid_vfdb_summary/{ST}_plasmid_summary.tsv` and
  `../Plasmid/plasmid_card_summary/{ST}_plasmid_summary.tsv`.

Run in that order: `run_mobsuite.sh` before `run_abricate_plasmids.sh`. Both
require external tools (`mob_suite`, `abricate`) and write their (large,
regenerable) output to a `Plasmid/` directory that sits next to the repo
root, not inside it — same convention the rest of the pipeline already uses
for `card_vfdb_result/`, `finder_result/`, and `output/`. This stage is
optional: `Virulence/Plasmid_VFDB.R`, `AMR/Plasmid_CARD.R`, and
`Analysis/Figure8.R` are all marked `optional = TRUE` in `run_pipeline.R`
and are skipped (not failed) if this stage hasn't been run yet.

## Pangenome/

- `run_pangenome_tree.sh` — Prokka annotation, PPanGGOLiN pangenome
  clustering, core-gene MSA, trimAL, and IQ-TREE, producing the
  `{ST}_bootstrap.treefile` and pangenome outputs that `Prerequisites/`
  depends on.

## Prerequisites/

Shared infrastructure every folder below depends on: builds the
finder-tool binary matrices (`00_build_finder_summaries.R`), shell-gene
clustering (`02b`/`02c`), tree-cluster mapping (`03c`/`03d`), and per-cluster
gene statistics (`04b`, required by `Analysis/Figure6.R`).

## Virulence/ and AMR/

Each script here is a raw per-database data-loading + statistical-analysis
script (temporal trend, clinical enrichment, early/late deltas) for exactly
one database:

- `Virulence/VFDB.R` and `Virulence/VirulenceFinder.R` are a clean split of
  the original combined virulence-analysis script; same statistical logic
  in both, only the database-specific data-loading section differs.
- `Virulence/Plasmid_VFDB.R` is the plasmid-vs-chromosomal classification
  layer for a curated set of virulence genes (the same `iuc`/`pap`/`sat`/`kps`
  gene list used in the manuscript's plasmid-context analysis), reading the
  `Plasmid/plasmid_vfdb_summary/` table `Plasmid_assembly/` produces.
- `AMR/CARD.R` and `AMR/ResFinder.R` are the equivalent split for AMR.
  `AMR/ResFinder.R` additionally folds in two ResFinder-only follow-on
  analyses, each in its own isolated `local({...})` block.
- `AMR/Plasmid_CARD.R` is the AMR-side equivalent of `Plasmid_VFDB.R`,
  reading the `Plasmid/plasmid_card_summary/` table. It classifies every
  CARD gene common to the whole-genome and plasmid tables rather than a
  curated subset, since (unlike the VFDB list) no curated AMR gene list for
  this analysis exists in the manuscript — see the file's own header.

## Analysis/ — figure -> script mapping (confirmed against the submitted PNGs)

| Figure | Script | Evidence |
|---|---|---|
| 1 | `Figure1.R` | Filename/title match to Figure_01.png, corroborated twice. |
| 2 | `Figure2.R` | Its combine step's title, "Shell-gene cluster silhouette and gene content," is an exact text match to Figure_02.png. |
| 3 | `Figure3.R` (**placeholder — no script**) | Built in iTOL, not any R script (confirmed). See the file for the upstream data files needed to rebuild it manually. |
| 4 | `Figure4.R` | "Oaxaca-style decomposition of burden change" title/subtitle is an exact text match to Figure_04.png. |
| 5 | `Figure5.R` | Content match (per-cluster Kendall-tau trends) to Figure_05.png. |
| 6 | `Figure6.R` | Script's own header matches Figure_06.png panel-for-panel. Depends on `Prerequisites/04b_cluster_gene_analysis.R`. |
| 7 | `Figure7.R` | Exact plot-title match ("Sensitivity analysis: temporal burden increase") to Figure_07.png. |
| 8 | `Figure8.R` | Independently produces both Figure_08.png panels (genomic location + RGP co-localization) itself. Needs the `Plasmid_assembly/` output. |

`Analysis/Supplementary.R` is a concatenation of every other original
analysis/figure script that is not one of the 8 above, each in its own
`local({...})` block, in dependency order.

## Notes

- Every R script only ever `source()`s `config.R` (never another R script,
  except where explicitly folded into a `local({...})` block as documented
  above), and all read/write paths go through `config$OUTPUT_DIR` — run
  `Rscript run_pipeline.R` from the repo root.
- `--from` and `--skip` in `run_pipeline.R` take the folder-qualified path
  (e.g. `--from Virulence/VFDB.R`), not the bare filename.
- No analysis logic was changed anywhere in this reorganisation. Every R
  script under `Prerequisites/`, `Virulence/`, `AMR/`, and `Analysis/` is
  either an unmodified relocation, a mechanical split at an existing
  internal section boundary, or a `local({...})`-wrapped concatenation.
  `Virulence/Plasmid_VFDB.R`, `AMR/Plasmid_CARD.R`, and everything in
  `Plasmid_assembly/`, `Setup/`, `Annotation/`, and `Pangenome/` were
  written fresh (or relocated from outside `scripts/`) to make the
  reorganised repo runnable end-to-end; they read/write data, they do not
  change any statistical method used elsewhere in the pipeline.
