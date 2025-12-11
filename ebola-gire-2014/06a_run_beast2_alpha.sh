#!/bin/bash
time ../beast277 \
     -threads -1 \
     -beagle \
     -working \
     beast2_run_alpha/ebola_alpha.xml

../delphy_mcc beast2_run_alpha/output.trees beast2_run_alpha/ebola_beast2.mcc
../loganalyser277 -burnin 30 beast2_run_alpha/output.log >beast2_run_alpha/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beast2_run_alpha/output.trees > beast2_run_alpha/tree_ess.json
