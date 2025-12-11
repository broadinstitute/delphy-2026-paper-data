#!/bin/bash
mkdir -p delphy_outputs_b

time ../delphy \
   --v0-in-fasta delphy_inputs/mpox-parker-2025.fasta \
   --v0-steps 1000000000 \
   --v0-out-log-file delphy_outputs_b/mpox-parker-2025.log \
   --v0-log-every 50000 \
   --v0-out-trees-file delphy_outputs_b/mpox-parker-2025.trees \
   --v0-tree-every 500000 \
   --v0-out-delphy-file delphy_outputs_b/mpox-parker-2025.dphy \
   --v0-delphy-snapshot-every 500000 \
   --v0-mpox-hack

../delphy_mcc delphy_outputs_b/mpox-parker-2025.trees delphy_outputs_b/mpox-parker-2025.mcc
../loganalyser277 -burnin 30 delphy_outputs_b/mpox-parker-2025.log >delphy_outputs_b/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 delphy_outputs_b/mpox-parker-2025.trees > delphy_outputs_b/tree_ess.json
