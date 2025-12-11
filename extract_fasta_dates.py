#!/usr/bin/env python3
#
# Usage: ./extract_fasta_dates.py <input.fasta >out_dates.csv
#
# Extract dates from a Delphy-ready FASTA file for use with TreeTime

import sys

print("name,date")
for line in sys.stdin:
    if line.startswith('>'):
        fasta_id, *rest = line[1:].split()
        fields = fasta_id.split('|')
        if len(fields) < 2:
            sys.stderr.write(f"ERROR: FASTA Id '{fasta_id}' doesn't have date field")
        else:
            date_field = fields[-1]
            print(f"{fasta_id},{date_field}")
