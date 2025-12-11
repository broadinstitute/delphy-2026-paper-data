#!/bin/bash
time ../beastX1050 \
     -beagle \
     -working \
     beastX_run_b/h3n2.xml

../delphy_mcc beastX_run_b/output.trees beastX_run_b/h3n2.mcc
../loganalyser277 -burnin 30 beastX_run_b/output.log >beastX_run_b/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beastX_run_b/output.trees > beastX_run_b/tree_ess.json
