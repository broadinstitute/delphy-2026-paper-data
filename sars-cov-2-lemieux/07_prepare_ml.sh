#!/bin/bash

mkdir -p ml
../extract_fasta_dates.py <delphy_inputs/ma_sars_cov_2.fasta >ml/ma_sars_cov_2_dates.csv
