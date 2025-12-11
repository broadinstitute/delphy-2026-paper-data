#!/bin/bash

time ../iqtree2 -s delphy_inputs/ebola.fasta -m HKY+FO+G4 --prefix ml/ebola_alpha
time treetime --tree ml/ebola_alpha.treefile --dates ml/ebola_dates.csv --aln delphy_inputs/ebola.fasta --outdir ml/tt_alpha --coalescent skyline --n-skyline 2 --stochastic-resolve

# Tweak TreeTime output to not confuse baltic
sed -i -e 's/ Tree /tree /g' ml/tt_alpha/timetree.nexus
