# Plasmid_assembly

This folder is reserved for the **mob-suite** plasmid assembly/typing step.

## Status: not present in this repository

After reviewing every script in the codebase, no mob-suite command, wrapper
script, or mob-suite output file exists anywhere in this repo. Per your
answer ("see the hdd"), the mob-suite step is run locally on your own
machine/hard drive, outside of this reproducibility repository -- its
output (the per-genome plasmid sequences/typing calls) is what feeds the
`Plasmid/` input directory that `Virulence/Plasmid_VFDB.R`,
`AMR/Plasmid_CARD.R`, and `Analysis/Figure8.R` all read from
(`config$BASE_DIR/../Plasmid/plasmid_vfdb_summary/{ST}_plasmid_summary.tsv`).

## What's needed to complete this folder

To make the pipeline fully reproducible end-to-end (assembly -> plasmid
typing -> plasmid VF/AMR annotation), this folder should eventually contain
the actual mob-suite invocation you use, for example something in the shape
of:

```bash
#!/usr/bin/env bash
# Plasmid_assembly/run_mobsuite.sh
# Run mob_recon (mob-suite) on each assembled genome to identify and
# reconstruct plasmid contigs, producing per-genome plasmid FASTA files
# that are subsequently annotated (VFDB/CARD) to build the
# Plasmid/plasmid_vfdb_summary and Plasmid/plasmid_card_summary tables.

for genome in assemblies/*.fasta; do
  sample=$(basename "$genome" .fasta)
  mob_recon --infile "$genome" --outdir "mob_suite_out/${sample}" --run_typer
done
```

Since no such command was available to extract from the existing repo, this
placeholder documents the gap transparently rather than inventing an
untested command. If you'd like this filled in with your actual mob-suite
parameters (database version, thresholds, etc.), just share the command you
run and it can be dropped in here verbatim.

## Downstream dependents

- `Virulence/Plasmid_VFDB.R`
- `AMR/Plasmid_CARD.R` (currently also a placeholder -- see that file)
- `Analysis/Figure8.R` (submitted Figure 8, Panel A "Genomic location" and
  Panel B "RGP co-localization")
