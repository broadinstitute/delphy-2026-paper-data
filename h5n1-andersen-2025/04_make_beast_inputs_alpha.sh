#!/bin/bash

# These generate equivalent input files for BEAST 2 and BEAST X that are statistically equivalent to
# a Delphy run with the parameters below.  Note that "--v0-steps" is the number of Delphy steps, and
# as a heuristic, the "equivalent" BEAST runs use 1/10th as many steps (similarly, the log and
# tree frequencies are reduced by a 1/10th)

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

mkdir -p beastX_run_alpha_a
../delphy \
   --dry-run \
   --v0-in-fasta delphy_inputs/h5n1-andersen-e756a15-ALL_full_dates_only.fasta \
   --v0-steps 20000000000 \
   --v0-log-every 2000000 \
   --v0-tree-every 20000000 \
   --v0-site-rate-heterogeneity \
   --v0-pop-model skygrid \
   --v0-skygrid-prior-double-half-time 0.25 \
   --v0-skygrid-type staircase \
   --v0-skygrid-disable-low-pop-barrier \
   --v0-skygrid-num-parameters 23 \
   --v0-skygrid-cutoff 1.833 \
   --v0-out-beast-xml beastX_run_alpha_a/h5n1-andersen-e756a15-ALL_full_dates_only.xml \
   --v0-out-beast-version X-10.5.0

mkdir -p beastX_run_alpha_b
cp beastX_run_alpha_a/h5n1-andersen-e756a15-ALL_full_dates_only.xml beastX_run_alpha_b/
