#!/bin/bash
mkdir -p delphy_outputs_a

# 1610 samples => ~10,000,000,000 steps
#
# Skygrid parameters are chosen to allow direct comparison with BEAST X run.
#
# Ordinarily, we'd use a log-linear Skygrid with much lower prior double-half time (the
# default 1 month instead of 6 months = 0.5 years).  However, in early test runs with a
# staircase and a 1-month double-half time, we observed the usual metastability of
# consecutive population intervals with wildly disparate populations leading to very many
# coalescences pinned to the end of the earlier low-population interval.
#
time ../delphy \
   --v0-in-fasta delphy_inputs/ebola_dudas.fasta \
   --v0-steps 10000000000 \
   --v0-out-log-file delphy_outputs_a/ebola_dudas_delphy.log \
   --v0-log-every 2000000 \
   --v0-out-trees-file delphy_outputs_a/ebola_dudas_delphy.trees \
   --v0-tree-every 10000000 \
   --v0-out-delphy-file delphy_outputs_a/ebola_dudas_delphy.dphy \
   --v0-delphy-snapshot-every 10000000 \
   --v0-pop-model skygrid \
   --v0-skygrid-type staircase \
   --v0-skygrid-prior-double-half-time 0.5 \
   --v0-skygrid-disable-low-pop-barrier \
   --v0-skygrid-num-parameters 25 \
   --v0-skygrid-cutoff 2

../delphy_mcc delphy_outputs_a/ebola_dudas_delphy.trees delphy_outputs_a/ebola_dudas_delphy.mcc
../loganalyser277 -burnin 30 delphy_outputs_a/ebola_dudas_delphy.log >delphy_outputs_a/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 delphy_outputs_a/ebola_dudas_delphy.trees > delphy_outputs_a/tree_ess.json
