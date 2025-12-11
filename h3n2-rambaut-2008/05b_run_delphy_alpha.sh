#!/bin/bash
mkdir -p delphy_outputs_alpha_b

# 165 samples => ~1,000,000,000 steps
#
# Skygrid parameters are chosen to allow direct comparison with BEAST X run, which retraces
# original workshop tutorial run.
# Latest tip is dated 23 Dec 2003.  Cutoff date is 23 Dec 1998 (5 years before), so intervals are
# a little over 1 month wide (49 intervals = 50 parameters).
#
# Ordinarily, we'd use a log-linear Skygrid with much lower prior double-half time (the
# default 1 month instead of 3 months = 0.25 years).  However, in early test runs with a
# staircase and a 1-month double-half time, we observed the usual metastability of
# consecutive population intervals with wildly disparate populations leading to very many
# coalescences pinned to the end of the earlier low-population interval.
#

time ../delphy \
   --v0-in-fasta delphy_inputs/h3n2.fasta \
   --v0-steps 1000000000 \
   --v0-out-log-file delphy_outputs_alpha_b/h3n2_alpha.log \
   --v0-log-every 100000 \
   --v0-out-trees-file delphy_outputs_alpha_b/h3n2_alpha.trees \
   --v0-tree-every 1000000 \
   --v0-out-delphy-file delphy_outputs_alpha_b/h3n2_alpha.dphy \
   --v0-delphy-snapshot-every 1000000 \
   --v0-site-rate-heterogeneity \
   --v0-pop-model skygrid \
   --v0-skygrid-prior-double-half-time 0.25 \
   --v0-skygrid-type staircase \
   --v0-skygrid-disable-low-pop-barrier \
   --v0-skygrid-num-parameters 50 \
   --v0-skygrid-cutoff 5

../delphy_mcc delphy_outputs_alpha_b/h3n2_alpha.trees delphy_outputs_alpha_b/h3n2_alpha.mcc
../loganalyser277 -burnin 30 delphy_outputs_alpha_b/h3n2_alpha.log >delphy_outputs_alpha_b/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 delphy_outputs_alpha_b/h3n2_alpha.trees > delphy_outputs_alpha_b/tree_ess.json
