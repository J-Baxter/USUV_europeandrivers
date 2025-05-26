#! /bin/bash
# CLI to generate json for lineage assignment
augur import beast --mcc USUV_2025May22_noFLI_NFLG_subsampled_SRD06_RelaxLn_constant_mcc.tree --most-recent-tip-date 2024.779 --output-tree beast_mcc.nwk --output-node-data beast_mcc.json
augur ancestral --tree beast_mcc.nwk --alignment USUV_2025May22_alldata_aligned_formatted_noFLI_NFLG_subsampled.fasta --output-node-data nt_muts.json --inference joint
augur export v2 --tree beast_mcc.nwk --node-data beast_mcc.json nt_muts.json --auspice-config auspice_config.json --output usuv.json

