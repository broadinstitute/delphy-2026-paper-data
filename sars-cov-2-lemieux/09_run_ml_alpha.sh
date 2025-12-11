#!/bin/bash

time ../iqtree2 -s delphy_inputs/ma_sars_cov_2.fasta -m HKY+FO+G4 --prefix ml/ma_sars_cov_2_alpha
time treetime --tree ml/ma_sars_cov_2_alpha.treefile --dates ml/ma_sars_cov_2_dates.csv --aln delphy_inputs/ma_sars_cov_2.fasta --outdir ml/tt_alpha --coalescent skyline --n-skyline 2 --stochastic-resolve

# Tweak TreeTime output to not confuse baltic
sed -i -e 's/ Tree /tree /g' ml/tt_alpha/timetree.nexus
