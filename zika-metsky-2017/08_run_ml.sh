#!/bin/bash

time ../iqtree2 -s delphy_inputs/zika.fasta -m HKY+FO --prefix ml/zika
time treetime --tree ml/zika.treefile --dates ml/zika_dates.csv --aln delphy_inputs/zika.fasta --outdir ml/tt --coalescent skyline --n-skyline 2 --stochastic-resolve

# Tweak TreeTime output to not confuse baltic
sed -i -e 's/ Tree /tree /g' ml/tt/timetree.nexus
