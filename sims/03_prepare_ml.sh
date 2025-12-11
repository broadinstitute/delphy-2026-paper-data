#!/bin/bash

for n in 100 1000  # Only for runs where we output a FASTA file
do
    for sim in exp_${n} const_${n}
    do
        ../extract_fasta_dates.py <${sim}/inputs/${sim}.fasta >${sim}/ml_outputs/${sim}_dates.csv
    done
done

