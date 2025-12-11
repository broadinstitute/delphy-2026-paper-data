#!/bin/bash

# These generate equivalent input files for BEAST 2 and BEAST X that are statistically equivalent to
# a Delphy run with the parameters below.  Note that "--v0-steps" is the number of Delphy steps, and
# as a heuristic, the "equivalent" BEAST runs use 1/10th as many steps (similarly, the log and
# tree frequencies are reduced by a 1/10th)

mkdir -p beast2_run
../delphy \
   --dry-run \
   --v0-in-fasta delphy_inputs/ma_sars_cov_2.fasta \
   --v0-steps 2000000000 \
   --v0-log-every 100000 \
   --v0-tree-every 1000000 \
   --v0-out-beast-xml beast2_run/ma_sars_cov_2.xml \
   --v0-out-beast-version 2.7.7

mkdir -p beastX_run
../delphy \
   --dry-run \
   --v0-in-fasta delphy_inputs/ma_sars_cov_2.fasta \
   --v0-steps 2000000000 \
   --v0-log-every 100000 \
   --v0-tree-every 1000000 \
   --v0-out-beast-xml beastX_run/ma_sars_cov_2.xml \
   --v0-out-beast-version X-10.5.0
