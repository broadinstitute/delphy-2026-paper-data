#!/bin/bash

mkdir -p ml
../extract_fasta_dates.py <delphy_inputs/ebola_dudas.fasta >ml/ebola_dudas_dates.csv
