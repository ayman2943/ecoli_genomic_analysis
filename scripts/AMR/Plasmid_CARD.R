#!/usr/bin/env Rscript
# ==============================================================================
# AMR / Plasmid_CARD.R  --  PLACEHOLDER (no corresponding analysis exists yet)
# ==============================================================================
#
# This folder's naming convention mirrors Virulence/ (VFDB.R, VirulenceFinder.R,
# Plasmid_VFDB.R) with a parallel AMR-side script for each. However, after a
# full read of every script in this repository -- including
# `10d_plasmid_context.R` (the sole plasmid-context script, which is
# VFDB/virulence-gene-specific: it loads `Plasmid/plasmid_vfdb_summary/` and
# classifies only the virulence "interest" gene list defined there) -- no
# script, output file, or intermediate data anywhere in this repo performs
# the equivalent plasmid-localization analysis for CARD/AMR resistance genes.
#
# In short: this repo currently contains a Virulence/Plasmid_VFDB.R analysis
# but NOT an AMR/Plasmid_CARD.R analysis. No `plasmid_card_summary` (or
# similarly named) file exists under `Plasmid/`, and no script reads or
# writes one.
#
# Per your instruction to "keep my analysis intact ... just rewrite those
# scripts to reproduce same images and files" -- this script is intentionally
# left as a placeholder rather than fabricating a new CARD-plasmid analysis
# that was never actually run for the manuscript. Nothing in the submitted
# figures or supplementary material depends on this file.
#
# TO COMPLETE THIS FOLDER: if a plasmid-associated CARD/AMR analysis exists
# (e.g. run separately, or on your local machine/HDD alongside the mob-suite
# step -- see Plasmid_assembly/README.md), point Claude (or fill in below) to:
#   1. The plasmid CARD/AMR summary TSV (analogous to
#      Plasmid/plasmid_vfdb_summary/{ST}_plasmid_summary.tsv, but for CARD),
#   2. The specific AMR "interest" gene list to classify (analogous to the
#      `interest` vector of VF genes in Plasmid_VFDB.R).
# Once supplied, this script would be written to mirror Plasmid_VFDB.R's
# structure exactly (clean_gid() helper, master AMR table load, plasmid CARD
# summary load, binary chromosomal/plasmid/absent classification).
#
# ==============================================================================

stop(
  "AMR/Plasmid_CARD.R is a placeholder. No plasmid-associated CARD/AMR ",
  "analysis exists in the original repository (confirmed by a full read of ",
  "10d_plasmid_context.R, which is VFDB/virulence-only). See the header ",
  "comment in this file for details and next steps."
)
