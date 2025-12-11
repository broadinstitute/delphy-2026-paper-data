#!/bin/bash

# These generate equivalent input files for BEAST 2 and BEAST X that are statistically equivalent to
# a Delphy run with the parameters below.  Note that "--v0-steps" is the number of Delphy steps, and
# as a heuristic, the "equivalent" BEAST runs use 1/10th as many steps (similarly, the log and
# tree frequencies are reduced by a 1/10th)

mkdir -p beast2_run_alpha
../delphy \
   --dry-run \
   --v0-in-fasta delphy_inputs/zika.fasta \
   --v0-steps 2000000000 \
   --v0-log-every 200000 \
   --v0-tree-every 1000000 \
   --v0-site-rate-heterogeneity \
   --v0-out-beast-xml beast2_run_alpha/zika_alpha.xml \
   --v0-out-beast-version 2.7.7

mkdir -p beastX_run_alpha
../delphy \
   --dry-run \
   --v0-in-fasta delphy_inputs/zika.fasta \
   --v0-steps 2000000000 \
   --v0-log-every 200000 \
   --v0-tree-every 1000000 \
   --v0-site-rate-heterogeneity \
   --v0-out-beast-xml beastX_run_alpha/zika_alpha.xml \
   --v0-out-beast-version X-10.5.0
