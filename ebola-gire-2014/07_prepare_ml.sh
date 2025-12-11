#!/bin/bash

mkdir -p ml
../extract_fasta_dates.py <delphy_inputs/ebola.fasta >ml/ebola_dates.csv
