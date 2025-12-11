#!/bin/bash
time ../beastX1050 \
     -beagle \
     -working \
     beastX_run/zika.xml

../delphy_mcc beastX_run/output.trees beastX_run/zika_beastX.mcc
../loganalyser277 -burnin 30 beastX_run/output.log >beastX_run/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 beastX_run/output.trees > beastX_run/tree_ess.json
