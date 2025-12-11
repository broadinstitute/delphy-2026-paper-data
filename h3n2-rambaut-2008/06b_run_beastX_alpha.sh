#!/bin/bash
time ../beastX1050 \
     -beagle \
     -working \
     beastX_run_alpha_b/h3n2_alpha.xml

../delphy_mcc beastX_run_alpha_b/output.trees beastX_run_alpha_b/h3n2_alpha.mcc
../loganalyser277 -burnin 30 beastX_run_alpha_b/output.log >beastX_run_alpha_b/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beastX_run_alpha_b/output.trees > beastX_run_alpha_b/tree_ess.json
