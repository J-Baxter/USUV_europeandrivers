#! /bin/bash
# CLI to generate json for lineage assignment
# First generate tree using nextstrain/augur
nextstrain shell.
augur tree --alignment USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsampled.fasta --output plain_tree.nwk 
augur refine --tree plain_tree.nwk --root mid_point --output-tree plain_tree.nwk
augur ancestral --tree plain_tree.nwk --alignment USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsampled.fasta --output-node-data nt_muts.json --inference joint
augur export v2 --tree plain_tree.nwk --node-data plain_tree.node_data.json nt_muts.json --auspice-config auspice_config.json --output 
exit

# Use CLI autolin to annotate tree with lineages
python annotate_json.py --input USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsampled.json --output USUV_2025Jun24_lineages.json --size 2 --distinction 10 --cutoff 0.99 --levels 3 --floor 1 --labels USUV_2025Jun24_labels.tsv --report USUV_2025Jun24_report.tsv



