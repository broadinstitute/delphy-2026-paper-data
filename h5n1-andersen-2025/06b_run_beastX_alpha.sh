#!/bin/bash
time ../beastX1050 \
     -beagle \
     -working \
     beastX_run_alpha_b/h5n1-andersen-e756a15-ALL_full_dates_only.xml

../delphy_mcc beastX_run_alpha_b/output.trees beastX_run_alpha_b/h5n1-andersen-e756a15-ALL_full_dates_only_alpha_beastX.mcc
../loganalyser277 -burnin 30 beastX_run_alpha_b/output.log >beastX_run_alpha_b/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beastX_run_alpha_b/output.trees > beastX_run_alpha_b/tree_ess.json
