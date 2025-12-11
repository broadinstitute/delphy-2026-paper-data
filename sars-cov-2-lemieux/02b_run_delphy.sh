#!/bin/bash
mkdir -p delphy_outputs_b

time ../delphy \
   --v0-in-fasta delphy_inputs/ma_sars_cov_2.fasta \
   --v0-steps 2000000000 \
   --v0-out-log-file delphy_outputs_b/ma_sars_cov_2_delphy.log \
   --v0-log-every 100000 \
   --v0-out-trees-file delphy_outputs_b/ma_sars_cov_2_delphy.trees \
   --v0-tree-every 1000000 \
   --v0-out-delphy-file delphy_outputs_b/ma_sars_cov_2_delphy.dphy \
   --v0-delphy-snapshot-every 2000000

../delphy_mcc delphy_outputs_b/ma_sars_cov_2_delphy.trees delphy_outputs_b/ma_sars_cov_2_delphy.mcc
../loganalyser277 -burnin 30 delphy_outputs_b/ma_sars_cov_2_delphy.log >delphy_outputs_b/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 delphy_outputs_b/ma_sars_cov_2_delphy.trees > delphy_outputs_b/tree_ess.json

