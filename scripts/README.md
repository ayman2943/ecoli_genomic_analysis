# scripts/ — one script per virulence/AMR database, one per figure

This folder was reorganised a second time, at your request, to go from
"one subfolder per pipeline stage" (the previous layout) to "one script per
figure / one script per database" — so a script's name tells you exactly
what manuscript output it reproduces, with no ambiguity.

`run_pipeline.R` (repo root) calls every script below in dependency order.
`00_config.R` stays at the top level because every stage sources it.

## Folder layout

| Folder | Contents | Purpose |
|---|---|---|
| `Prerequisites/` | `00_build_finder_summaries.R`, `02b_pangenome_shell_cluster_vfdb.R`, `02c_pangenome_shell_cluster_vf.R`, `03c_tree_mapping_vfdb.R`, `03d_tree_mapping_vf.R`, `04b_cluster_gene_analysis.R` | Shared infrastructure every other folder depends on: builds the finder-tool binary matrices, shell-gene clustering, tree-cluster mapping, and per-cluster gene statistics. Not named by the four folders you asked for, because these scripts don't correspond to one virulence/AMR database or one figure — they're upstream of all of them. Dropping this folder would break Virulence/, AMR/, and Analysis/Figure6.R. |
| `Virulence/` | `VFDB.R`, `VirulenceFinder.R`, `Plasmid_VFDB.R` | Raw per-database virulence gene burden/annotation, one script per database. |
| `AMR/` | `CARD.R`, `ResFinder.R`, `Plasmid_CARD.R` | Raw per-database AMR gene burden/annotation, one script per database. |
| `Plasmid_assembly/` | `README.md` (placeholder — see below) | mob-suite plasmid assembly/typing. |
| `Analysis/` | `Figure1.R` … `Figure8.R`, `Supplementary.R` | One script per published figure, plus one script for everything supplementary. |

## Virulence/ and AMR/

Each script here is a raw per-database data-loading + statistical-analysis
script (temporal trend, clinical enrichment, early/late deltas — the full
per-database analysis, not just annotation) for exactly one database:

- `Virulence/VFDB.R` and `Virulence/VirulenceFinder.R` are a clean split of
  the original `04a_virulence_analysis.R`, which analysed both databases
  back-to-back in one file. Same statistical logic in both, unchanged —
  only the database-specific data-loading section differs.
- `Virulence/Plasmid_VFDB.R` is the general-purpose plasmid-vs-chromosomal
  classification layer for virulence genes, extracted from
  `10d_plasmid_context.R`'s data-loading sections (that script's
  figure-specific "top 10 genes in Cluster_3" logic stays in
  `Analysis/Figure8.R`, which is the full, unmodified original script).
- `AMR/CARD.R` and `AMR/ResFinder.R` are the equivalent split of
  `05a_resistance_analysis.R`. `AMR/ResFinder.R` additionally folds in the
  two follow-on scripts that only ever ran on ResFinder data
  (`05b_resfinder_clusters.R`, `05c_resfinder_decreasing.R`), each in its
  own `local({...})` block so their variables can't collide with the main
  ResFinder analysis above them in the same file.
- `AMR/Plasmid_CARD.R` is a **placeholder**. After reading every script in
  the repo, no CARD/AMR equivalent of `10d_plasmid_context.R` exists —
  there is no plasmid CARD summary file and no script that classifies AMR
  genes as chromosomal vs. plasmid-associated. Rather than invent an
  analysis that was never run for the manuscript, this file documents the
  gap and stops with an explanatory message if executed. See its header
  for what would be needed to fill it in.

## Plasmid_assembly/

No mob-suite script or output exists anywhere in this repository — per
your note, that step is run on your own machine/hard drive, outside this
reproducibility repo. `Plasmid_assembly/README.md` documents this and
lists the downstream files (`Plasmid/plasmid_vfdb_summary/...`) that its
output feeds into. If you share the actual mob-suite command you run, it
can be dropped in as a real script.

## Analysis/ — figure -> script mapping (confirmed against the submitted PNGs)

The previous round of this repo flagged that several scripts in the old
`09_figures/` folder had self-contradictory internal "figNN" labels that
didn't match the manuscript's actual figure numbers (see the previous
"Figures: open questions" section, preserved below for the record). That
was resolved this round by reading the actual submitted
`Figures/Figure_0N.png` files directly and matching their visual content
and exact title text against the code, rather than trusting each script's
own comments or filename:

| Figure | Script | Evidence |
|---|---|---|
| 1 | `Figure1.R` (= `06a_temporal_mk.R`, unchanged) | Filename/title match, corroborated twice. |
| 2 | `Figure2.R` (extracted from `12_fig02-05_ST69_analysis.R`) | Its `fig2_3_combined` step's title, "Shell-gene cluster silhouette and gene content," is an exact text match to Figure_02.png. |
| 3 | `Figure3.R` (**placeholder — no script**) | **Confirmed by you: made in iTOL**, not any R script. See the file for the upstream data (tree + cluster/clinicality/region/K-type/year annotations) needed to rebuild it manually. |
| 4 | `Figure4.R` (extracted from `12_fig02-05_ST69_analysis.R`) | Its "Oaxaca-style decomposition of burden change" title/subtitle is an exact text match to Figure_04.png. |
| 5 | `Figure5.R` (= `13_fig06_cluster_temporal.R`, renamed) | Content match (per-cluster Kendall-tau trends) to Figure_05.png. Original filename said "fig06" — that was wrong. |
| 6 | `Figure6.R` (= `16_fig07_phylogeny_summary.R`, renamed) | Script's own header ("Row 1: ST69 Cluster_3 genes + clinical enrichment; Row 2: ST10 VF composition + trends + AMR decline") matches Figure_06.png panel-for-panel. Original filename said "fig07" — that was wrong. Depends on `Prerequisites/04b_cluster_gene_analysis.R`. |
| 7 | `Figure7.R` (= `09a_sensitivity_analysis.R`, unchanged) | Exact plot-title match ("Sensitivity analysis: temporal burden increase") to Figure_07.png. |
| 8 | `Figure8.R` (= `10d_plasmid_context.R`, unchanged) | Full read confirms it independently produces both Figure_08.png panels (genomic location + RGP co-localization) itself. Needs the external Plasmid dataset. |

`Analysis/Supplementary.R` is a concatenation of every other original
analysis/figure script that is not one of the 8 above, each in its own
`local({...})` block (17 blocks total, including a full unmodified copy of
`12_fig02-05_ST69_analysis.R` for completeness — see the file's own header
for exactly what that covers and why). The append-chain order
(`09b` -> `09c` -> `09d` -> `10e`, which all write to the same
`capsule_classification.xlsx`) is preserved. Two pure post-processing
utilities from the old `09_figures/` folder — `19_fig_composite.R` and
`99_ALL_FIGURES_COMBINED.R` — were left out entirely: both just re-paste
already-generated PNGs into combined layouts and perform no analysis of
their own; they remain in git history if needed.

## Notes

- Every script only ever `source()`s `config.R` (never another script,
  except where explicitly folded into a `local({...})` block as documented
  above), and all read/write paths go through `config$OUTPUT_DIR` — so
  moving/renaming files here doesn't change any script's behaviour. Run
  `Rscript run_pipeline.R` from the repo root exactly as before.
- `--from` and `--skip` in `run_pipeline.R` now take the folder-qualified
  path (e.g. `--from Virulence/VFDB.R`, `--skip Analysis/Figure3.R`), not
  the bare filename.
- No analysis logic was changed anywhere in this reorganisation — every
  script is either an unmodified relocation, a mechanical split at an
  existing internal section boundary, or a `local({...})`-wrapped
  concatenation. Where a script's own filename previously disagreed with
  the manuscript's actual figure numbering, only the name was corrected;
  the code inside is untouched.

## Previous round's open questions (resolved above; kept for the record)

The scripts in the old `09_figures/` folder were not all written against
the final, 8-figure manuscript structure, and their filenames didn't all
agree with what they actually produced. That ambiguity has now been
resolved by direct comparison against the submitted `Figures/Figure_0N.png`
files (see the table above) and confirmation from you on Figure 3's iTOL
origin. The specific mismatches previously flagged were:

- `12_fig02-05_ST69_analysis.R` covered multiple figures in one script —
  resolved by extracting exactly its Figure-2-producing and
  Figure-4-producing sections into `Figure2.R`/`Figure4.R`, and keeping
  the full original as a Supplementary.R backup block.
- `13_fig06_cluster_temporal.R` was labelled "fig06" but is actually
  **Figure 5** — confirmed and renamed.
- `14_fig06_gene_trajectories.R` predates the current numbering and never
  matched any submitted figure — kept in `Supplementary.R` as exploratory
  content, not treated as a main figure.
- `15_fig07_tree_parsimony.R` was labelled "fig07" but its content
  (single-ring tree, not the real multi-ring Figure 3) doesn't match any
  submitted figure — kept in `Supplementary.R`, not used for Figure 3
  (which was made in iTOL — see `Figure3.R`).
- `16_fig07_phylogeny_summary.R` was labelled "fig07" but is actually
  **Figure 6** — confirmed and renamed.
- `19_fig_composite.R` and `99_ALL_FIGURES_COMBINED.R` were utility/export
  scripts rather than sources of new analysis — left out of this repo
  entirely (see above).
