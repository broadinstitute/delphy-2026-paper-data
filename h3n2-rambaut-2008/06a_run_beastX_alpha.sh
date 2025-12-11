#!/bin/bash
time ../beastX1050 \
     -beagle \
     -working \
     beastX_run_alpha_a/h3n2_alpha.xml

../delphy_mcc beastX_run_alpha_a/output.trees beastX_run_alpha_a/h3n2_alpha.mcc
../loganalyser277 -burnin 30 beastX_run_alpha_a/output.log >beastX_run_alpha_a/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beastX_run_alpha_a/output.trees > beastX_run_alpha_a/tree_ess.json
