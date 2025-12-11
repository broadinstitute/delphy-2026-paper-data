#!/bin/bash

time ../iqtree2 -s delphy_inputs/zika.fasta -m HKY+FO+G4 --prefix ml/zika_alpha
time treetime --tree ml/zika_alpha.treefile --dates ml/zika_dates.csv --aln delphy_inputs/zika.fasta --outdir ml/tt_alpha --coalescent skyline --n-skyline 2 --stochastic-resolve

# Tweak TreeTime output to not confuse baltic
sed -i -e 's/ Tree /tree /g' ml/tt_alpha/timetree.nexus
