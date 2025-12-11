#!/bin/bash

time ../iqtree2 -s delphy_inputs/h5n1-andersen-e756a15-ALL_full_dates_only.fasta -m HKY+FO --prefix ml/h5n1_andersen
time treetime --tree ml/h5n1_andersen.treefile --dates ml/h5n1_andersen_dates.csv --aln delphy_inputs/h5n1-andersen-e756a15-ALL_full_dates_only.fasta --outdir ml/tt --coalescent skyline --n-skyline 23 --stochastic-resolve  --reconstruct-tip-states

# Tweak TreeTime output to not confuse baltic
sed -i -e 's/ Tree /tree /g' ml/tt/timetree.nexus
