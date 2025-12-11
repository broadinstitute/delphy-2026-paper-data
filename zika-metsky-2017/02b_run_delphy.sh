#!/bin/bash
mkdir -p delphy_outputs_b

time ../delphy \
   --v0-in-fasta delphy_inputs/zika.fasta \
   --v0-steps 2000000000 \
   --v0-out-log-file delphy_outputs_b/zika_delphy.log \
   --v0-log-every 200000 \
   --v0-out-trees-file delphy_outputs_b/zika_delphy.trees \
   --v0-tree-every 1000000 \
   --v0-out-delphy-file delphy_outputs_b/zika_delphy.dphy \
   --v0-delphy-snapshot-every 1000000

../delphy_mcc delphy_outputs_b/zika_delphy.trees delphy_outputs_b/zika_delphy.mcc
../loganalyser277 -burnin 30 delphy_outputs_b/zika_delphy.log >delphy_outputs_b/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 delphy_outputs_b/zika_delphy.trees > delphy_outputs_b/tree_ess.json
