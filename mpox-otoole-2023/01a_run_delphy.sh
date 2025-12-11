#!/bin/bash
mkdir -p delphy_outputs_a

time ../delphy \
   --v0-in-fasta delphy_inputs/mpox-otoole-2023.fasta \
   --v0-steps 200000000 \
   --v0-out-log-file delphy_outputs_a/mpox-otoole-2023.log \
   --v0-log-every 10000 \
   --v0-out-trees-file delphy_outputs_a/mpox-otoole-2023.trees \
   --v0-tree-every 100000 \
   --v0-out-delphy-file delphy_outputs_a/mpox-otoole-2023.dphy \
   --v0-delphy-snapshot-every 100000 \
   --v0-mpox-hack

../delphy_mcc delphy_outputs_a/mpox-otoole-2023.trees delphy_outputs_a/mpox-otoole-2023.mcc
../loganalyser277 -burnin 30 delphy_outputs_a/mpox-otoole-2023.log >delphy_outputs_a/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 delphy_outputs_a/mpox-otoole-2023.trees > delphy_outputs_a/tree_ess.json
