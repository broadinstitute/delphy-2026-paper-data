#!/bin/bash
mkdir -p delphy_outputs_alpha_b

# 3339 samples => ~20,000,000,000 steps
#
# Skygrid parameters are chosen to allow direct comparison with BEAST X run.
# Latest tip is dated 3 Aug 2025.
# Cutoff date is in early Oct 2023, intervals are 1 month wide => 22 intervals => 23 parameters.
#
# Ordinarily, we'd use a log-linear Skygrid with much lower prior double-half time (the
# default 1 month instead of 3 months = 0.25 years).  However, in early test runs with a
# staircase and a 1-month double-half time, we observed the usual metastability of
# consecutive population intervals with wildly disparate populations leading to very many
# coalescences pinned to the end of the earlier low-population interval.
#
# An early run here showed that errors arising from the discretization of the coalescent prior
# with default settings (~400 cells) were leading to small but observable differences between
# the results of BEAST X and Delphy.  The final runs were thus made with a finer discretization
# (~1000 cells)
#
time ../delphy \
   --v0-in-fasta delphy_inputs/h5n1-andersen-e756a15-ALL_full_dates_only.fasta \
   --v0-steps   20000000000 \
   --v0-out-log-file delphy_outputs_alpha_b/h5n1-andersen-e756a15-ALL_full_dates_only_alpha_b.log \
   --v0-log-every   2000000 \
   --v0-out-trees-file delphy_outputs_alpha_b/h5n1-andersen-e756a15-ALL_full_dates_only_alpha_b.trees \
   --v0-tree-every 20000000 \
   --v0-out-delphy-file delphy_outputs_alpha_b/h5n1-andersen-e756a15-ALL_full_dates_only_alpha_b.dphy \
   --v0-delphy-snapshot-every 20000000 \
   --v0-site-rate-heterogeneity \
   --v0-target-coal-prior-cells 1000 \
   --v0-pop-model skygrid \
   --v0-skygrid-prior-double-half-time 0.25 \
   --v0-skygrid-type staircase \
   --v0-skygrid-disable-low-pop-barrier \
   --v0-skygrid-num-parameters 23 \
   --v0-skygrid-cutoff 1.833

../delphy_mcc delphy_outputs_alpha_b/h5n1-andersen-e756a15-ALL_full_dates_only_alpha_b.trees delphy_outputs_alpha_b/h5n1-andersen-e756a15-ALL_full_dates_only_alpha_b.mcc
../loganalyser277 -burnin 30 delphy_outputs_alpha_b/h5n1-andersen-e756a15-ALL_full_dates_only_alpha_b.log >delphy_outputs_alpha_b/log_analysis.txt
../calc-tree-ess --compact --burnin-pct 30 delphy_outputs_alpha_b/h5n1-andersen-e756a15-ALL_full_dates_only_alpha_b.trees > delphy_outputs_alpha_b/tree_ess.json
