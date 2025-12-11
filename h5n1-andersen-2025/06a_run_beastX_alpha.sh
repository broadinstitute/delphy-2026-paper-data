#!/bin/bash
time ../beastX1050 \
     -beagle \
     -working \
     beastX_run_alpha_a/h5n1-andersen-e756a15-ALL_full_dates_only.xml

../delphy_mcc beastX_run_alpha_a/output.trees beastX_run_alpha_a/h5n1-andersen-e756a15-ALL_full_dates_only_alpha_beastX.mcc
../loganalyser277 -burnin 30 beastX_run_alpha_a/output.log >beastX_run_alpha_a/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beastX_run_alpha_a/output.trees > beastX_run_alpha_a/tree_ess.json
