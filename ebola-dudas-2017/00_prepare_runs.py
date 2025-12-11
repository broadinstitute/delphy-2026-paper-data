#!/usr/bin/env python

import sys
import requests
from pathlib import Path
import subprocess
import shutil
import xml.etree.ElementTree as ET
import re
import csv

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

# Clone Github repo
# =================
commit_hash = '9db59a4'
print(f"\nCloning Github repo ebov/space-time at commit {commit_hash}")
data_repo_path = raw_path / 'space-time'
do_clone = True
if data_repo_path.exists():
    do_clone = False
    print(f'[SKIPPING] {data_repo_path} already exists')

if do_clone:
    subprocess.run([
        'git',
        'clone',
        'https://github.com/ebov/space-time.git'],
                   cwd=raw_path.as_posix())
    subprocess.run([
        'git',
        'checkout',
        commit_hash],
                   cwd=data_repo_path.as_posix())

# Extract sequences from FASTA file
# =================================

def read_fasta(path):
    seqs = []
    current_seq_id = ''
    current_seq = ''
    with open(path, 'r') as raw_f:
        for line in raw_f:
            if line.startswith('>'):
                if current_seq:
                    seqs.append((current_seq_id, current_seq))
                current_seq_id = line[1:].strip()
                current_seq = ''
            else:
                current_seq += line.strip()
        if current_seq:
            seqs.append((current_seq_id, current_seq))
    return seqs

in_fasta_file_path = (
    data_repo_path /
    'Data' /
    'Makona_1610_genomes_2016-06-23.fasta'
)

print("\nReading sequence data from FASTA file...")
in_f = read_fasta(in_fasta_file_path)

print(f'Found {len(in_f)} aligned sequences')

# Prepare fasta with all the sequences
# ====================================
# We simply sanitize the IDs so that the Delphy UI doesn't label all of them "EBOV"
#
print("\nReconstituting input FASTA with sanitized IDs...")

def id_2_short_id(seq_id):        # The initial pair of pieces of an ID is unique, so we can correlate to the metadata
    pieces = seq_id.split('|')
    return '|'.join(pieces[0:2])

id_2_sane_id = {}
id_2_date = {}  # During curation, Dudas et al imputed exact dates for all sequences
short_id_2_id = {}
for (seq_id, seq) in in_f:
    pieces = seq_id.split('|')
    id_2_sane_id[seq_id] = '-'.join(pieces[:-1])
    id_2_date[seq_id] = pieces[-1]
    
    short_id = id_2_short_id(seq_id)
    assert short_id not in short_id_2_id
    short_id_2_id[short_id] = seq_id
    
input_fasta_path = delphy_inputs_path / 'ebola_dudas.fasta'
do_reassembly = True
if input_fasta_path.exists():
    do_reassembly = False
    print(f'[SKIPPING] {input_fasta_path} already exist')

if do_reassembly:
    with open(input_fasta_path, 'w') as f:
        for (seq_id, seq) in in_f:
            f.write('>')
            f.write(f'{id_2_sane_id[seq_id]}|{id_2_date[seq_id]}')
            f.write('\n')
            f.write(seq)
            f.write('\n')
            
    print(f'Output FASTA with sanitized IDs {input_fasta_path}')

# Read in metadata
# ================
metadata_file_path = (
    data_repo_path
    / 'Data'
    / 'Makona_1610_metadata_2016-06-23.csv'
)
print(f"\nParsing metadata in {metadata_file_path}")

id_2_metadata = {}
with metadata_file_path.open('r') as f:
    ff = csv.reader(f, delimiter=',')
    header = None
    for row in ff:
        if not header:
            header = row
            continue

        (num, final_label, virus, sample_id, sequence_id, sample_lab, sequencing_lab, accession, country, location, platform, prop_ambig, date_field, collection_date, imputed_date, reference, GP82) = row

        short_id = id_2_short_id(final_label)  # IDs in metadata aren't 1-1 identical to IDs in FASTA (!!)
        seq_id = short_id_2_id[short_id]
        assert seq_id in id_2_sane_id, seq_id
        assert seq_id not in id_2_metadata
        
        id_2_metadata[seq_id] = {
            'country': country,
            'location': location,
        }

# Check that we have metadata for everything
# (as of this writing, we do)
assert len(id_2_sane_id) == len(id_2_metadata)

# Prepare metadata file
# =====================
print("\nPreparing metadata...")
input_metadata_path = delphy_inputs_path / 'ebola_dudas_metadata.csv'
do_metadata = True
if input_metadata_path.exists():
    do_metadata = False
    print(f'[SKIPPING] {input_metadata_path} already exist')

country_2_sane_country = {
    'SLE': 'Sierra Leone',
    'GIN': 'Guinea',
    'LBR': 'Liberia',
    '?': '-',
    '': '-',
}

def sane_location(country, location):
    if location in ('?', ''):
        return country_2_sane_country[country]
    else:
        return f'{location} - {country_2_sane_country[country]}'

if do_metadata:
    with open(input_metadata_path, 'w') as f:
        f.write('id,Country,Location\n')
        for (seq_id, _) in in_f:
            sane_id = id_2_sane_id[seq_id]
            country = id_2_metadata[seq_id]["country"]
            location = id_2_metadata[seq_id]["location"]
            f.write(f'{sane_id},{country_2_sane_country[country]},{sane_location(country, location)}\n')

    print(f'Metadata written to {input_metadata_path}')
