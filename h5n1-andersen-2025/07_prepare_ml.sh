#!/bin/bash

mkdir -p ml
mkdir -p ml_alpha
../extract_fasta_dates.py <delphy_inputs/h5n1-andersen-e756a15-ALL_full_dates_only.fasta >ml/h5n1_andersen_dates.csv
