#!/usr/bin/env bash
# =============================================================================
# config/pipeline_config.sh
#
# Site-specific settings for the E. coli pipeline. Edit the values below for
# your machine, or override any of them on the command line / environment.
# Sourced by run_all.sh.
# =============================================================================

# ---- Target lineage (also set via TARGET_ST env or first arg to run_all.sh) ----
# export TARGET_ST="ST69"

# ---- Paths ----
# export RAW_ENTERO_EXPORT="/path/to/EnteroBase_export.xlsx"
# export SPECIES="Escherichia_coli"

# ---- NCBI API key (OPTIONAL - speed up downloads; NEVER commit real keys) ----
# Put the key in a local .env file or export it in your shell:
#   export NCBI_API_KEY="your_key_here"

# ---- Databases ----
# export VF_DB="$HOME/databases/virulencefinder_db"
# export RF_DB="$HOME/databases/resfinder_db"

# ---- Compute ----
export JOBS="${JOBS:-6}"
export CPU="${CPU:-8}"
export PPANG_RAM="${PPANG_RAM:-16}"
export IQTREE_MEM="${IQTREE_MEM:-14G}"

# ---- External pangenome dir (non-ST69 lineages, if stored elsewhere) ----
# export ECOLI_PANGENOME_DIR="/Volumes/ayman_ssd/Ecoli/All_genome/annotation/pangenome_output"

echo "[pipeline_config] JOBS=$JOBS CPU=$CPU (override via env)"
