#!/bin/bash
time ../beastX1050 \
     -beagle \
     -working \
     beastX_run_alpha/ma_sars_cov_2_alpha.xml

../delphy_mcc beastX_run_alpha/output.trees beastX_run_alpha/ma_sars_cov_2_alpha_beastX.mcc
../loganalyser277 -burnin 30 beastX_run_alpha/output.log >beastX_run_alpha/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beastX_run_alpha/output.trees > beastX_run_alpha/tree_ess.json
