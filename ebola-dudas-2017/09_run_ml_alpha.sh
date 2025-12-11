#!/bin/bash

time ../iqtree2 -s delphy_inputs/ebola_dudas.fasta -m HKY+FO+G4 --prefix ml/ebola_dudas_alpha
time treetime --tree ml/ebola_dudas_alpha.treefile --dates ml/ebola_dudas_dates.csv --aln delphy_inputs/ebola_dudas.fasta --outdir ml/tt_alpha --coalescent skyline --n-skyline 25 --stochastic-resolve  --reconstruct-tip-states

# Tweak TreeTime output to not confuse baltic
sed -i -e 's/ Tree /tree /g' ml/tt_alpha/timetree.nexus
