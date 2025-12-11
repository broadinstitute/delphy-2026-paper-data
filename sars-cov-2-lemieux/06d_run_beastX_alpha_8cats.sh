#!/bin/bash
time ../beastX1050 \
     -beagle \
     -working \
     beastX_run_alpha_8cats/ma_sars_cov_2_alpha_8cats.xml

../delphy_mcc beastX_run_alpha_8cats/output.trees beastX_run_alpha_8cats/ma_sars_cov_2_alpha_8cats_beastX.mcc
../loganalyser277 -burnin 30 beastX_run_alpha_8cats/output.log >beastX_run_alpha_8cats/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beastX_run_alpha_8cats/output.trees > beastX_run_alpha_8cats/tree_ess.json
