#!/bin/bash

time ../iqtree2 -s delphy_inputs/ebola.fasta -m HKY+FO --prefix ml/ebola
time treetime --tree ml/ebola.treefile --dates ml/ebola_dates.csv --aln delphy_inputs/ebola.fasta --outdir ml/tt --coalescent skyline --n-skyline 2 --stochastic-resolve

# Tweak TreeTime output to not confuse baltic
sed -i -e 's/ Tree /tree /g' ml/tt/timetree.nexus
