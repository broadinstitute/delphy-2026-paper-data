#!/bin/bash
time ../beastX1050 \
     -threads 8 \
     -beagle \
     -working \
     beastX_run/mpox-parker-2025-beastX.xml

../delphy_mcc beastX_run/Mpox_2poch_combined.trees beastX_run/mpox-parker-2025-beastX.mcc
../loganalyser277 -burnin 30 beastX_run/Mpox_2poch_combined.log >beastX_run/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beastX_run/Mpox_2poch_combined.trees > beastX_run/tree_ess.json
