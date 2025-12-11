#!/bin/bash

mkdir -p ml
../extract_fasta_dates.py <delphy_inputs/zika.fasta >ml/zika_dates.csv
