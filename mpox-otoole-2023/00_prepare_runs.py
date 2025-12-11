#!/usr/bin/env python3

from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

# Prepare directory structure
# ===========================
raw_path = Path('raw')
raw_path.mkdir(parents=True, exist_ok=True)
delphy_inputs_path = Path('delphy_inputs')
delphy_inputs_path.mkdir(parents=True, exist_ok=True)

# Clone O'Toole et al data repo
# =============================
print("\nCloning O'Toole et al 2023 data repo")
data_repo_path = raw_path / 'apobec3'
do_clone = True
if data_repo_path.exists():
    do_clone = False
    print(f'[SKIPPING] {data_repo_path} already exists')

if do_clone:
    subprocess.run([
        'git',
        'clone',
        'https://github.com/hmpxv/apobec3.git'],
                   cwd=raw_path.as_posix())

beast_file_path = data_repo_path / 'data' / 'apobec3_2partition.epoch.xml'
if not beast_file_path.exists():
    raise ValueError(f"BEAST file not found in O'Toole et al 2023 data repo: {beast_file_path}")
    
# Extract sequences from BEAST file
# =================================
print("\nReading sequence data from BEAST XML file...")
beastXml = ET.parse(beast_file_path.as_posix())
alignments = beastXml.findall("alignment")
apobec3_alignment = [a for a in alignments if a.get('id') == 'apobec3_alignment'][0]
non_apobec3_alignment = [a for a in alignments if a.get('id') == 'non_apobec3_alignment'][0]

included_sequences = []

# Check that sequences are listed in the same order in both alignments
for apo, nonapo in zip(apobec3_alignment.findall('sequence'),
                       non_apobec3_alignment.findall('sequence')):
    apo_id = apo.find("taxon").get("idref")
    non_apo_id = nonapo.find("taxon").get("idref")
    
    if apo_id != non_apo_id:
        raise ValueError(f'Sequences in BEAST files alignments listed in different order (first discrepancy: {apo_id} vs {non_apo_id})')

for d in apobec3_alignment.findall('sequence'):
    theId = d.find("taxon").get('idref')
    theDate = theId.split('|')[-1]
    
    # Skip the 1 GISAID sequence
    if theId.startswith('EPI'):
        print(f'Skipping missing GISAID sequence {theId}')
        continue

    # Skip any sequences before 2017 (i.e., outgroups KJ642615 from 1978 and KJ642617 from 1971, before this spillover)
    yy = int(theDate[0:4])
    if yy < 2017:
        print(f'Skipping pre-spillover sequence {theId}')
        continue
    
    included_sequences.append(theId)
    
    print(f'Found sequence for {theId}, dated to {theDate}')

# Prepare fasta with all the sequences
# ====================================
print("\nReconstituting multiple sequence alignment from 2 partitions in BEAST XML file...")

input_fasta_path = delphy_inputs_path / 'mpox-otoole-2023.fasta'
do_reassembly = True
if input_fasta_path.exists():
    do_reassembly = False
    print(f'[SKIPPING] {input_fasta_path} already exist')

if do_reassembly:
    with open(input_fasta_path, 'w') as f:
        for apo, nonApo in zip(apobec3_alignment.findall('sequence'),
                               non_apobec3_alignment.findall('sequence')):
            apoId = apo.find("taxon").get("idref")
            nonApoId = nonApo.find("taxon").get("idref")
            
            if apoId != nonApoId:
                raise ValueError('APO and non-APO sequences with different ids?')
            
            theId = apoId
            theDate = theId.split('|')[-1]
    
            if apoId not in included_sequences:
                continue
    
            apoSeq = ''.join(t for t in apo.itertext()).strip()
            nonApoSeq = ''.join(t for t in nonApo.itertext()).strip()
            states = []
            for apoState, nonApoState in zip(apoSeq, nonApoSeq):
                if apoState == 'N':
                    states.append(nonApoState)
                else:
                    states.append(apoState)
            
            seq = ''.join(states)
                
            f.write('>')
            f.write(theId)
            f.write('\n')
            f.write(seq)
            f.write('\n')


# Prepare metadata file
# =====================
print("\nPreparing metadata...")
input_metadata_path = delphy_inputs_path / 'mpox-otoole-2023_metadata.csv'
do_metadata = True
if input_metadata_path.exists():
    do_metadata = False
    print(f'[SKIPPING] {input_metadata_path} already exist')

if do_metadata:
    with open(input_metadata_path, 'w') as f:
        f.write('id,Geo\n')
        for theId in included_sequences:
            seqId = theId.split('|')[0]
            geo, date = theId.split('|')[-2:]
            
            f.write(f'{seqId},{geo}\n')

    print(f'Metadata written to {input_metadata_path}')
