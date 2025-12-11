#!/bin/bash

# These generate equivalent input files for BEAST 2 and BEAST X that are statistically equivalent to
# a Delphy run with the parameters below.  Note that "--v0-steps" is the number of Delphy steps, and
# as a heuristic, the "equivalent" BEAST runs use 1/10th as many steps (similarly, the log and
# tree frequencies are reduced by a 1/10th)

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

mkdir -p beastX_run_a
../delphy \
   --dry-run \
   --v0-in-fasta delphy_inputs/ebola_dudas.fasta \
   --v0-steps 10000000000 \
   --v0-log-every 2000000 \
   --v0-tree-every 10000000 \
   --v0-pop-model skygrid \
   --v0-skygrid-type staircase \
   --v0-skygrid-prior-double-half-time 0.5 \
   --v0-skygrid-disable-low-pop-barrier \
   --v0-skygrid-num-parameters 25 \
   --v0-skygrid-cutoff 2 \
   --v0-out-beast-xml beastX_run_a/ebola_dudas.xml \
   --v0-out-beast-version X-10.5.0

mkdir -p beastX_run_b
cp beastX_run_a/ebola_dudas.xml beastX_run_b/
