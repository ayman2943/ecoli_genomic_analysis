# Prerequisites/

Scripts in this folder are not virulence, AMR, plasmid, or figure scripts
themselves -- they are infrastructure that produces the master shell-gene
cluster tables, tree-cluster mappings, and cluster-gene statistics that the
Virulence/, AMR/, and Analysis/ scripts all depend on. This folder was
added beyond the four you named (Virulence, AMR, Plasmid_assembly,
Analysis) because these scripts are genuine, unavoidable prerequisites --
dropping them would break the pipeline (e.g. Analysis/Figure6.R directly
reads output written by 04b_cluster_gene_analysis.R below).

Run in this order (also reflected in run_pipeline.R):

  1. 00_build_finder_summaries.R   (raw finder-tool output -> binary matrices)
  2. 02b_pangenome_shell_cluster_vfdb.R   (ST69 shell clustering)
  3. 02c_pangenome_shell_cluster_vf.R    (ST10 shell clustering)
  4. 03c_tree_mapping_vfdb.R    (ST69 cluster -> tree mapping)
  5. 03d_tree_mapping_vf.R     (ST10 cluster -> tree mapping)
  6. 04b_cluster_gene_analysis.R  (per-cluster gene statistics; required by Figure6.R)

Then Virulence/, AMR/, and Analysis/ scripts can run (in any order among
themselves, except where noted in their own headers).
