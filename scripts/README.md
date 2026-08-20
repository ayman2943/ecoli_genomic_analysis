# scripts/ — analysis pipeline, organised by stage

33 R scripts, grouped into one subfolder per pipeline stage so you can find
things by what they do rather than by a bare number. `run_pipeline.R` (repo
root) calls every script below in this same order; `00_config.R` stays at
the top level because every stage sources it.

| Folder | Stage | Scripts | Produces |
|---|---|---|---|
| `00_config.R` | shared config (paths, target ST, thresholds) | — | — |
| `01_pangenome_clustering/` | PPanGGOLiN shell-gene clustering | `02b_pangenome_shell_cluster_vfdb.R`, `02c_pangenome_shell_cluster_vf.R` | Cluster assignments (VFDB- and VirulenceFinder-based) |
| `02_tree_mapping/` | map clusters onto the core-genome tree | `03c_tree_mapping_vfdb.R`, `03d_tree_mapping_vf.R` | Tip-annotated tree data |
| `03_virulence_resistance/` | per-genome VF/AMR burden + per-cluster gene content | `04a_virulence_analysis.R`, `04b_cluster_gene_analysis.R`, `05a_resistance_analysis.R`, `05b_resfinder_clusters.R`, `05c_resfinder_decreasing.R` | Table 1 inputs, Results §1-2 |
| `04_temporal_trends/` | Mann-Kendall trend test, 5 STs x 4 databases | `06a_temporal_mk.R` | **Figure 1**, Results §1 |
| `05_st10_decomposition/` | ST10 Oaxaca-style decomposition + ARG | `07a_st10_decomposition.R`, `07b_st10_composition_drivers.R`, `10f_st10_arg_decomposition.R` | **Figure 4** (ST10 panels), Results §7-8 |
| `06_clinical_sensitivity/` | clinical enrichment + adjusted sensitivity models | `08_clinical_enrichment_3panel.R`, `09a_sensitivity_analysis.R` | **Figure 7**, **Table 2**, Results §4, 5 |
| `07_capsule_ktype/` | capsule group / K-type classification (VirulenceFinder `kpsM`-based, corrected per Reviewer 2 #1) | `09b_capsule_classification.R`, `09c_capsule_comparison.R`, `09d_k_type_analysis.R`, `09e_summary_stats.R` | Supplementary Table S2, Results §9 |
| `08_validation/` | KPS/RGP/plasmid/allelic-conversion + long-read checks | `10a_kps_validation.R`, `10b_rgp_neighbourhood.R`, `10c_rgp_genomic_context.R`, `10d_plasmid_context.R`, `10e_allelic_conversion.R`, `10g_long_read_validation.R` | **Figure 8**, Supplementary Table S4, Results §6 |
| `09_figures/` | figure assembly scripts | see below | Figures 2-8 + Supplementary |

## Figures: open questions

The scripts in `09_figures/` were not all written against the final,
8-figure manuscript structure, and their filenames don't all agree with
what they actually produce. Before treating this folder as the definitive
source of the submitted figures, it's worth confirming against your own
records which script you last ran for each PNG in `Figures/`:

- `12_fig02-05_ST69_analysis.R` — one script covering multiple figures
  (silhouette/clustering, PCA, gene content by cluster, decomposition,
  tree mapping). Its own `ggsave()` calls use internal names
  (`Fig2_silhouette_3panel`, `Fig5_decomposition_3panel`, ...) that don't
  literally match `Figure_0N.png`, so treat it as a source of panels to
  assemble rather than a 1:1 figure generator.
- `13_fig06_cluster_temporal.R` — by content (per-cluster Kendall τ trends)
  this looks like **Figure 5**, not Figure 6.
- `14_fig06_gene_trajectories.R` — its own header says it "replicates the
  Figure 5/6 analysis from `3_generate_figures.R`" — a script that no
  longer exists in this repo, a sign this file predates the current
  numbering and was never fully reconciled with it.
- `15_fig07_tree_parsimony.R` — labelled "fig07" but its content (mapping
  shell-gene clusters onto the core-genome tree) and its own output names
  (`Figure_1_shell_clusters_on_core_tree...`, `Figure_2_...`,
  `Figure_3_...`) point to **Figure 3**, not Figure 7. It also looks
  exploratory — it saves three tree-layout variants rather than one final
  panel.
- `16_fig07_phylogeny_summary.R` — labelled "fig07" and "phylogeny", but
  its own header says "Row 1: ST69 Cluster 3 genes + clinical enrichment;
  Row 2: ST10 VF composition + trends + AMR decline" and it writes to an
  `OUT` folder named `sensitivity_analysis`. That description matches
  **Figure 6** (panels A-H), not Figure 7 or a phylogeny.
- `19_fig_composite.R` and `99_ALL_FIGURES_COMBINED.R` are utility/export
  scripts (Word-ready single-page exports, and an all-in-one combined PDF)
  rather than sources of new analysis.

None of this changes what's in the manuscript — the published figures are
already final. It only means the *script-to-figure* labelling in this
folder is inherited from an earlier round of exploration and doesn't fully
match the current manuscript numbering. If you remember which of `13_`/`14_`
and `15_`/`16_` you actually used last, renaming the "winner" to match its
real figure number (and moving the other to an `archive/` subfolder) would
close this out — happy to do that once you confirm.

## Notes

- Every script only ever `source()`s `config.R` (never another script), and
  all read/write paths go through `config$OUTPUT_DIR` — so moving files
  into subfolders here doesn't change any script's behaviour. Run
  `Rscript run_pipeline.R` from the repo root exactly as before.
- `--from` and `--skip` in `run_pipeline.R` still take the bare script
  filename (e.g. `--from 09e`, `--skip 15`), not the folder path.
