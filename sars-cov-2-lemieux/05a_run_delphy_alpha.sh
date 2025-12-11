#!/bin/bash
mkdir -p delphy_outputs_alpha_a

time ../delphy \
   --v0-in-fasta delphy_inputs/ma_sars_cov_2.fasta \
   --v0-steps 2000000000 \
   --v0-out-log-file delphy_outputs_alpha_a/ma_sars_cov_2_delphy_alpha.log \
   --v0-log-every 100000 \
   --v0-out-trees-file delphy_outputs_alpha_a/ma_sars_cov_2_delphy_alpha.trees \
   --v0-tree-every 1000000 \
   --v0-out-delphy-file delphy_outputs_alpha_a/ma_sars_cov_2_delphy_alpha.dphy \
   --v0-site-rate-heterogeneity \
   --v0-delphy-snapshot-every 2000000

../delphy_mcc delphy_outputs_alpha_a/ma_sars_cov_2_delphy_alpha.trees delphy_outputs_alpha_a/ma_sars_cov_2_delphy_alpha.mcc
../loganalyser277 -burnin 30 delphy_outputs_alpha_a/ma_sars_cov_2_delphy_alpha.log >delphy_outputs_alpha_a/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 delphy_outputs_alpha_a/ma_sars_cov_2_delphy_alpha.trees > delphy_outputs_alpha_a/tree_ess.json
