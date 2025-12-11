#!/bin/bash
time ../beast277 \
     -threads -1 \
     -beagle \
     -working \
     beast2_run/ma_sars_cov_2.xml

../delphy_mcc beast2_run/output.trees beast2_run/ma_sars_cov_2_beast2.mcc
../loganalyser277 -burnin 30 beast2_run/output.log >beast2_run/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beast2_run/output.trees > beast2_run/tree_ess.json
