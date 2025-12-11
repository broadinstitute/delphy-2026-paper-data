#!/bin/bash

time ../iqtree2 -s delphy_inputs/ebola_dudas.fasta -m HKY+FO --prefix ml/ebola_dudas
time treetime --tree ml/ebola_dudas.treefile --dates ml/ebola_dudas_dates.csv --aln delphy_inputs/ebola_dudas.fasta --outdir ml/tt --coalescent skyline --n-skyline 25 --stochastic-resolve  --reconstruct-tip-states

# Tweak TreeTime output to not confuse baltic
sed -i -e 's/ Tree /tree /g' ml/tt/timetree.nexus
