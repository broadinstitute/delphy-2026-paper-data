#!/bin/bash
time ../beastX1050 \
     -beagle \
     -working \
     beastX_run_a/h3n2.xml

../delphy_mcc beastX_run_a/output.trees beastX_run_a/h3n2.mcc
../loganalyser277 -burnin 30 beastX_run_a/output.log >beastX_run_a/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beastX_run_a/output.trees > beastX_run_a/tree_ess.json
