#!/usr/bin/env python

import sys
import requests
from pathlib import Path
import subprocess
import shutil
import xml.etree.ElementTree as ET
import re

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

# Download original supplementary data
# ====================================
print("\nDownloading Metsky et al (2017) supplementary data")
si_data_path = raw_path / '41586_2017_BFnature22402_MOESM6_ESM.zip'
do_download = True
if si_data_path.exists():
    do_download = False
    print(f'[SKIPPING] {si_data_path} already exists')

if do_download:
    # Some HTTP trickery makes it nontrivial to download this file simply using the Python `requests` library;
    # instead of digging into why that is, just punt the job to `curl`, which does it correctly
    subprocess.run([
        'curl',
        '-O',
        'https://static-content.springer.com/esm/art%3A10.1038%2Fnature22402/MediaObjects/41586_2017_BFnature22402_MOESM6_ESM.zip'],
                   cwd=raw_path.as_posix())

# Unpack SI zip file and get BEAST file
# =====================================

print("\nUnpacking supplementary data zip file...")
beast_file_path = (
    scratch_path /
    'SupplementaryData' /
    'BEAST input and output' /
    'Phylogenetic analyses and model selection' /
    'SRD06-strict-exponential.xml')
do_unzip = True
if beast_file_path.exists():
    do_unzip = False
    print(f'[SKIPPING] {beast_file_path} already exists')

if do_unzip:
    shutil.unpack_archive(si_data_path, scratch_path)
    if not beast_file_path.exists():
        raise ValueError(f'{beast_file_path} not in supplementary info zip file?')

# Extract sequences from BEAST file
# =================================
# We'll extract the essentials of the run in `SRD06-strict-exponential.xml`.  The original data
# (in `SupplementaryData/Sequences and alignments/alignment-used-for-beast-analyses.no-outgroup.fasta`) is not quite
# aligned to the sequence that defines the coordinates in the paper
# ([KX197192.1](https://www.ncbi.nlm.nih.gov/nuccore/KX197192.1)).  However, the BEAST analyses only use the coding
# regions of this sequence (sites 108-10379, inclusive).  Mercifully, that portion is aligned to KX197192.1 and is what appears in the BEAST files, unscrambled.

print("\nReading sequence data from BEAST XML file...")
beastXml = ET.parse(beast_file_path.as_posix())
xmlTaxa = [taxon for taxon in beastXml.find("taxa")]
xmlAlignments = beastXml.findall("alignment")
print(f'Found {len(xmlAlignments)} alignments and {len(xmlTaxa)} taxa')

# Check that taxa are listed in the same order in <taxa> section and in alignment
taxaIds = [taxon.get('id') for taxon in xmlTaxa]
alignmentIdRefs = [s.find('taxon').get('idref') for s in xmlAlignments[0].findall('sequence')]
if taxaIds != alignmentIdRefs:
    raise ValueError('Taxa are listed in different order in <taxa> section and in alignment')

# Prepare fasta with all the sequences
# ====================================
print("\nReconstituting multiple sequence alignment from BEAST XML file...")

input_fasta_path = delphy_inputs_path / 'zika.fasta'
do_reassembly = True
if input_fasta_path.exists():
    do_reassembly = False
    print(f'[SKIPPING] {input_fasta_path} already exist')

if do_reassembly:
    with open(input_fasta_path, 'w') as f:
        for (taxon, s) in zip(xmlTaxa,
                              xmlAlignments[0].findall('sequence')):
            taxonId = taxon.get('id')
            gbId, geo, date = taxonId.split('|')
            
            # "Align" CDS to reference KX197192.1 (UTRs at both ends are just masked)
            finalSeq = ('N'*107) + re.sub(r"\s+", "", ''.join(s.itertext())) + ('N'*428)
            
            f.write('>')
            f.write(f'{gbId}|{date}')
            f.write('\n')
            f.write(finalSeq)
            f.write('\n')
            
    print(f'Reconstituted multiple sequence alignment written to {input_fasta_path}')

# Prepare metadata file
# =====================
print("\nPreparing metadata...")
input_metadata_path = delphy_inputs_path / 'zika_metadata.csv'
do_metadata = True
if input_metadata_path.exists():
    do_metadata = False
    print(f'[SKIPPING] {input_metadata_path} already exist')

if do_metadata:
    with open(input_metadata_path, 'w') as f:
        f.write('id,Geo\n')
        for taxon in xmlTaxa:
            taxonId = taxon.get('id')
            gbId, geo, date = taxonId.split('|')
        
            f.write(f'{gbId},{geo}\n')

    print(f'Metadata written to {input_metadata_path}')
