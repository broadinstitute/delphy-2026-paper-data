#!/usr/bin/env python
#
# This vignette retraces Exercise 2 in the BEAST workshop tutorial "Revealing the
# evolutionary dynamics of influenza"
# (https://beast.community/workshop_influenza_phylodynamics), using the included H3N2
# dataset, itself a subset of the data analyzed in Rambaut et al, "The genomic and
# epidemiological dynamics of human influenza A virus", Nature (453), pp. 615--619 (doi:
# 10.1038/nature06945).

import sys
import requests
from pathlib import Path
import subprocess
import shutil
import xml.etree.ElementTree as ET
import re
import csv
import math
import datetime

if len(sys.argv) != 1:
    sys.stderr.write("Usage: ./00_prepare_runs.py\n")
    sys.exit(1)

# Prepare directory structure
# ===========================
raw_path = Path('raw')
raw_path.mkdir(parents=True, exist_ok=True)
scratch_path = Path('scratch')
scratch_path.mkdir(parents=True, exist_ok=True)
delphy_inputs_path = Path('delphy_inputs')
delphy_inputs_path.mkdir(parents=True, exist_ok=True)

# Download H3N2 dataset in FASTA format
# =====================================
print("\nDownloading H3N2 dataset from BEAST workshop tutorials website")
raw_nexus_path = raw_path / 'NewYork.HA.2000-2003.nex'
do_download = True
if raw_nexus_path.exists():
    do_download = False
    print(f'[SKIPPING] {raw_nexus_path} already exists')

if do_download:
    subprocess.run([
        'curl',
        '-O',
        'https://beast.community/tutorials/workshop_influenza_phylodynamics/files/NewYork.HA.2000-2003.nex'],
                   cwd=raw_path.as_posix())


# Extract sequences from NEXUS file
# =================================
# The below is the world's most fragile NEXUS file reader

print("\nReading sequence data from NEXUS file...")

with raw_nexus_path.open('r') as f:
    nexus_lines = [line.strip() for line in f]
assert nexus_lines[:6] == [
    '#NEXUS',
    '',
    'Begin DATA;',
    'Dimensions ntax=165 nchar=1698;',
    'Format datatype=NUCLEOTIDE gap=-;',
    'Matrix',
], repr(nexus_lines[:6])
# Skip nexus_lines[6]
raw_seqs = []
for line in nexus_lines[7:]:
    if not line:
        continue
    if line == ';':
        break

    if m := re.match(r'^\'([^\']+)\'\s+(.*)$', line):
        seq_id = m[1]
        seq = m[2]
        assert all(c in 'ACGTURYSWKMBDHVN.-' for c in seq), seq  # See https://www.bioinformatics.org/sms/iupac.html

        raw_seqs.append((seq_id, seq))

first_seq_id, first_seq = raw_seqs[0]
assert all(len(seq) == len(first_seq) for (seq_id, seq) in raw_seqs)
        
print(f'Found {len(raw_seqs)} aligned sequences of length {len(first_seq)}')

# Extract dates from sequence IDs
# ===============================
# Sequence IDs have the form 'NewYork_NNN_YYYY.YY' (with varying precision on NNN and YYYY.YY)
print("\nExtracting sequence dates")
id_2_date = {}
for seq_id, _ in raw_seqs:
    
    state_part, num_part, date_part = seq_id.split('_')
    date_as_linear_year = float(date_part)
    year = math.floor(date_as_linear_year)
    year_fraction = date_as_linear_year - year
    days_in_year = (datetime.date(year+1, 1, 1) - datetime.date(year, 1, 1)).days
    seq_date = datetime.date(year, 1, 1) + datetime.timedelta(days=int(round(year_fraction * days_in_year)))

    id_2_date[seq_id] = f'{seq_date.year:04d}-{seq_date.month:02d}-{seq_date.day:02d}'

# Prepare fasta with all the sequences
# ====================================
#
print("\nPreparing input FASTA...")
input_fasta_path = delphy_inputs_path / 'h3n2.fasta'
do_reassembly = True
if input_fasta_path.exists():
    do_reassembly = False
    print(f'[SKIPPING] {input_fasta_path} already exist')

if do_reassembly:
    with open(input_fasta_path, 'w') as f:
        for (seq_id, seq) in raw_seqs:
            f.write('>')
            f.write(f'{seq_id}|{id_2_date[seq_id]}')
            f.write('\n')
            f.write(seq)
            f.write('\n')
            
    print(f'Output FASTA {input_fasta_path}')
