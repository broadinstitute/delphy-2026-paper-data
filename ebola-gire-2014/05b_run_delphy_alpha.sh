#!/bin/bash
mkdir -p delphy_outputs_alpha_b

time ../delphy \
   --v0-in-fasta delphy_inputs/ebola.fasta \
   --v0-steps 1000000000 \
   --v0-out-log-file delphy_outputs_alpha_b/ebola_delphy_alpha.log \
   --v0-log-every 50000 \
   --v0-out-trees-file delphy_outputs_alpha_b/ebola_delphy_alpha.trees \
   --v0-tree-every 500000 \
   --v0-out-delphy-file delphy_outputs_alpha_b/ebola_delphy_alpha.dphy \
   --v0-delphy-snapshot-every 1000000 \
   --v0-site-rate-heterogeneity

../delphy_mcc delphy_outputs_alpha_b/ebola_delphy_alpha.trees delphy_outputs_alpha_b/ebola_delphy_alpha.mcc
../loganalyser277 -burnin 30 delphy_outputs_alpha_b/ebola_delphy_alpha.log >delphy_outputs_alpha_b/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 delphy_outputs_alpha_a/ebola_delphy_alpha.trees > delphy_outputs_alpha_b/tree_ess.json
