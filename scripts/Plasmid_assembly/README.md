# Plasmid_assembly/

Reconstructs plasmid contigs from each whole-genome assembly (mob-suite),
then annotates them (abricate VFDB + CARD) — the input `Virulence/Plasmid_VFDB.R`,
`AMR/Plasmid_CARD.R`, and `Analysis/Figure8.R` need.

## Scripts, run in order

1. **`run_mobsuite.sh`** — runs `mob_recon` on every assembly in `{ST}/*.fna`
   to reconstruct and type plasmid contigs. For each genome with at least
   one predicted plasmid, its plasmid contigs are concatenated into a
   single `{genome}_combined_plasmids.fasta`.

   ```bash
   bash scripts/Plasmid_assembly/run_mobsuite.sh ST69
   bash scripts/Plasmid_assembly/run_mobsuite.sh all      # all 5 STs
   ```

   Requires `mob_recon` (mob-suite):
   ```bash
   conda install -c bioconda mob_suite
   ```

2. **`run_abricate_plasmids.sh`** — runs abricate (VFDB + CARD, 80% identity
   / 80% coverage — same thresholds as the whole-genome annotation) on each
   combined plasmid FASTA from step 1, then builds the per-ST summary
   tables:

   ```bash
   bash scripts/Plasmid_assembly/run_abricate_plasmids.sh ST69
   ```

   Requires `abricate` (already needed for `Annotation/run_abricate.sh`).

## Output

Both scripts write to a `Plasmid/` directory that sits **next to** the repo
root (`../Plasmid/` relative to this repo), not inside it — like
`card_vfdb_result/`, `finder_result/`, and `output/` elsewhere in this
pipeline, it's large, fully regenerable data, not code, so it isn't
committed to git:

```
../Plasmid/
├── mob_suite_raw/{ST}/{genome}/          raw mob_recon output per genome
├── combined_fasta/{ST}/{genome}_combined_plasmids.fasta
├── plasmid_vfdb/{ST}/{genome}_combined_plasmids_vfdb.tsv
├── plasmid_card/{ST}/{genome}_combined_plasmids_card.tsv
├── plasmid_vfdb_summary/{ST}_plasmid_summary.tsv    <- read by Virulence/Plasmid_VFDB.R, Analysis/Figure8.R
└── plasmid_card_summary/{ST}_plasmid_summary.tsv    <- read by AMR/Plasmid_CARD.R
```

## Why this is optional

Not every genome has a predicted plasmid, and this stage adds a real
compute cost (mob-suite + a second abricate pass) on top of the
whole-genome annotation `Annotation/` already does. `run_pipeline.R` marks
`Virulence/Plasmid_VFDB.R`, `AMR/Plasmid_CARD.R`, and `Analysis/Figure8.R`
as `optional = TRUE` — they run automatically once this stage's output
exists, and are cleanly skipped (not failed) otherwise.
