#!/bin/bash

time ../iqtree2 -s delphy_inputs/ma_sars_cov_2.fasta -m HKY+FO --prefix ml/ma_sars_cov_2
time treetime --tree ml/ma_sars_cov_2.treefile --dates ml/ma_sars_cov_2_dates.csv --aln delphy_inputs/ma_sars_cov_2.fasta --outdir ml/tt --coalescent skyline --n-skyline 2 --stochastic-resolve

# Tweak TreeTime output to not confuse baltic
sed -i -e 's/ Tree /tree /g' ml/tt/timetree.nexus
