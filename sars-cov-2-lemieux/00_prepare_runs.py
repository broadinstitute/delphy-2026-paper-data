#!/usr/bin/env python

from Bio import Entrez
import re
import datetime
import sys, os
from pathlib import Path
import subprocess

if len(sys.argv) != 2:
    sys.stderr.write("Usage: ./00_prepare_runs.py <email-address-for-Entrez>\n")
    sys.exit(1)

# Config
# ======
Entrez.email = sys.argv[1]
path_to_mafft = Path("../mafft")

# Prepare directory structure
# ===========================
raw_path = Path('raw')
raw_path.mkdir(parents=True, exist_ok=True)
scratch_path = Path('scratch')
scratch_path.mkdir(parents=True, exist_ok=True)
delphy_inputs_path = Path('delphy_inputs')
delphy_inputs_path.mkdir(parents=True, exist_ok=True)

# Read in sample IDs & clades
# ===========================
# sample_ids.txt = Sample ids for 772 sequences used in Fig 3A of LeMieux et al (2021) (private communication)
print("\nReading in sample IDs...")
sample_id_2_clade = {}
with open('sample_ids.csv', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue   # Skip comments
        if line.startswith('id,'):
            continue   # Skip header
        stripped_line = line.strip()
        if not stripped_line:
            continue   # Skip empty lines
        
        id, clade = stripped_line.split(',')
        sample_id_2_clade[id] = clade

sample_ids = list(sample_id_2_clade.keys())
print(f'Read {len(sample_ids)} samples')
print(f'First 5 IDs: {sample_ids[:5]}')
print(f'Last 5 IDs: {sample_ids[-5:]}')

# Look up accession IDs in GenBank
# ================================
print("\nLooking up accession IDs in GenBank...")
accession_ids_path = raw_path / 'accession_ids.csv'
do_lookup = True
if accession_ids_path.exists():
    do_lookup = False
    print(f'[SKIPPING] {accession_ids_path} already exists')

def lookup_gb_ids(sample_ids):
    with Entrez.esearch(db="nucleotide", term=' OR '.join(sample_ids), retmax=len(sample_ids)) as handle:
        record = Entrez.read(handle)
        print(f'Errors: {record["ErrorList"]}')
        return record['IdList']

if do_lookup:
    acc_ids = lookup_gb_ids(sample_ids)
    with open(accession_ids_path, 'w') as f:
        for acc_id in acc_ids:
            f.write(f'{acc_id}\n')
else:
    acc_ids = []
    with open(accession_ids_path, 'r') as f:
        for line in f:
            if line.strip():
                acc_ids.append(line.strip())

print(f'Found {len(acc_ids)} accession IDs vs {len(sample_ids)} samples')
print(f'First 5 accession IDs: {acc_ids[:5]}')
print(f'Last 5 accession IDs: {acc_ids[-5:]}')


# Fetch sequences from GenBank
# ============================
print("\nFetching sequences from GenBank...")
raw_genomes_path = raw_path / 'raw_genomes.xml'
do_fetch = True
if raw_genomes_path.exists():
    do_fetch = False
    print(f'[SKIPPING] {raw_genomes_path} already exists')

def fetch_sequences(ids, out_f):
    with Entrez.efetch(db="nucleotide", id=','.join(ids), rettype="gb", retmode="xml") as handle:
        while True:
            lines = handle.readlines(1000)
            if not lines:
                break
            out_f.write(b''.join(lines))

if do_fetch:
    with open(raw_genomes_path, 'wb') as f:
        fetch_sequences(acc_ids, f)

        
# Fetch ref sequence from GenBank
# ===============================
print("\nFetching reference sequence from GenBank...")
raw_reference_path = raw_path / 'raw_reference.xml'
do_fetch_ref = True
if raw_reference_path.exists():
    do_fetch_ref = False
    print(f'[SKIPPING] {raw_reference_path} already exists')

ref_acc_id = 'NC_045512.2'
    
if do_fetch_ref:
    with open(raw_reference_path, 'wb') as f:
        fetch_sequences([ref_acc_id], f)


# Convert raw GenBank genomes to FASTA files suitable for Delphy
# ==============================================================
print("\nPreparing unaligned and untrimmed FASTA...")
scratch_raw_genomes_fasta_path = scratch_path / 'raw_genomes.fasta'
scratch_ref_fasta_path = scratch_path / 'ref.fasta'
do_convert = True
if scratch_raw_genomes_fasta_path.exists() and scratch_ref_fasta_path.exists():
    do_convert = False
    print(f'[SKIPPING] {scratch_raw_genomes_fasta_path} and {scratch_ref_fasta_path} already exist')

def parse_genbank_date(gb_date):
    """
    Parses a date from GenBank to a pair of min and max dates, encoded as ((yyyy, mm, dd), (yyyy, mm, dd)) min-max tuples.
    
    Formats recognized from GenBank:
    * 1978
    * 09-Nov-2017
    * Nov-2017
    * 2019-01-16
    * 2018-08
    """
    
    months = {
        'jan': 1,
        'feb': 2,
        'mar': 3,
        'apr': 4,
        'may': 5,
        'jun': 6,
        'jul': 7,
        'aug': 8,
        'sep': 9,
        'oct': 10,
        'nov': 11,
        'dec': 12
    }
    
    # 1978
    if m := re.match(r'^(\d\d\d\d)$', gb_date):
        yyyy = int(m[1])
        return ((yyyy, 1, 1), (yyyy+1, 1, 1))
        
    # Nov-2017
    if m := re.match(r'^(\w\w\w)-(\d\d\d\d)$', gb_date):
        if m[1].lower() in months:
            yyyy = int(m[2])
            mm = months[m[1].lower()]
            min_date = (yyyy, mm, 1)
            if mm < 12:
                max_date = (yyyy, mm+1, 1)
            else:
                max_date = (yyyy+1, 1, 1)
            return (min_date, max_date)
        
    # 09-Nov-2017
    if m := re.match(r'^(\d\d)-(\w\w\w)-(\d\d\d\d)$', gb_date):
        if m[2].lower() in months:
            yyyy = int(m[3])
            mm = months[m[2].lower()]
            dd = int(m[1])
            return ((yyyy, mm, dd), (yyyy, mm, dd))
        
    # 2019-01-16
    if m := re.match(r'^(\d\d\d\d)-(\d\d)-(\d\d)$', gb_date):
        yyyy = int(m[1])
        mm = int(m[2])
        dd = int(m[3])
        return ((yyyy, mm, dd), (yyyy, mm, dd))
        
    # 2018-08
    if m := re.match(r'^(\d\d\d\d)-(\d\d)$', gb_date):
        yyyy = int(m[1])
        mm = int(m[2])
        min_date = (yyyy, mm, 1)
        if mm < 12:
            max_date = (yyyy, mm+1, 1)
        else:
            max_date = (yyyy+1, 1, 1)
        return (min_date, max_date)
        
    # Nothing
    return ((1900, 1, 1), (1900, 1, 1))

def unparse_date_range(date_range):
    (min_date, max_date) = date_range
    (yyyy1, mm1, dd1) = min_date
    (yyyy2, mm2, dd2) = max_date

    if min_date == max_date:
        return f'{yyyy1:04d}-{mm1:02d}-{dd1:02d}'

    if min_date == (yyyy1, 1, 1) and max_date == (yyyy1+1, 1, 1):
        return f'{yyyy1:04d}'

    if ((min_date == (yyyy1, mm1, 1) and max_date == (yyyy1, mm1+1, 1)) or
        (min_date == (yyyy1, 12, 1) and max_date == (yyyy1+1, 1, 1))):
        return f'{yyyy1:04d}-{mm1:02d}'

    return f'{yyyy1:04d}-{mm1:02d}-{dd1:02d}/{yyyy2:04d}-{mm2:02d}-{dd2:02d}'

def gb_to_fasta(in_gb_xml_f, out_fasta_f, default_assembly_name='UNKNOWN'):
    records = Entrez.parse(in_gb_xml_f)
    for record in records:
        
        acc = record['GBSeq_accession-version']
        
        # Extract collection date if it's there
        raw_collection_date = ''
        collection_date = datetime.date(1900, 1, 1)
        for feature in record['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'source':
                for qual in feature['GBFeature_quals']:
                    if qual['GBQualifier_name'] == 'collection_date':
                        raw_collection_date = qual['GBQualifier_value']
                        collection_date_range = parse_genbank_date(raw_collection_date)

        # Extract Assembly Name (what was used in LeMieux et al)
        m = re.search(r'Assembly Name :: (\S+) ;', record['GBSeq_comment'])
        assembly_name = default_assembly_name
        if m:
            assembly_name = m[1]
        
        seq = record['GBSeq_sequence']
        
        out_fasta_f.write(f'>{assembly_name}|{acc}|{unparse_date_range(collection_date_range)}\n')
        out_fasta_f.write(seq)
        out_fasta_f.write('\n')

if do_convert:
    with open(raw_genomes_path, 'rb') as gb_xml:
        with open(scratch_raw_genomes_fasta_path, 'w') as fasta:
            gb_to_fasta(gb_xml, fasta)

    with open(raw_reference_path, 'rb') as gb_xml:
        with open(scratch_ref_fasta_path, 'w') as fasta:
            gb_to_fasta(gb_xml, fasta, ref_acc_id)

# Align with mafft
# ================
print("\nAligning sequences with mafft")
scratch_aligned_fasta_path = scratch_path / 'aligned.fasta'
do_align = True
if scratch_aligned_fasta_path.exists():
    do_align = False
    print(f'[SKIPPING] {scratch_aligned_fasta_path} already exist')

if do_align:
    subprocess.run([
        path_to_mafft.as_posix(),
        "--version"
    ])

    with open(scratch_aligned_fasta_path, 'w') as f:
        subprocess.run([
            path_to_mafft.as_posix(),
            "--thread", "-1",
            "--auto",
            "--addfragments", scratch_raw_genomes_fasta_path.as_posix(),
            "--keeplength",
            scratch_ref_fasta_path.as_posix(),
        ], stdout=f)

# Trim UTRs as in LeMieux et al (2021)
# ====================================
# (we don't bother with masking out problematic sites as would usually be done with
# this kind of data; the point of this dataset is to compare BEAST to Delphy, not to get a perfect tree)
print("\nTrimming UTRs...")
input_fasta_path = delphy_inputs_path / 'ma_sars_cov_2.fasta'
do_trim = True
if input_fasta_path.exists():
    do_trim = False
    print(f'[SKIPPING] {input_fasta_path} already exist')

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

if do_trim:
    aligned_seqs = read_fasta(scratch_aligned_fasta_path)
    print(f'Alignment has {len(aligned_seqs)} sequences')
    print(f'First sequence is {aligned_seqs[0][0]}, with {len(aligned_seqs[0][1])} bases')

    with open(input_fasta_path, 'w') as out_f:
        for seq_id, seq in aligned_seqs:
            if not seq_id.startswith(ref_acc_id):
                seq = '-'*267 + seq[267:-230] + '-'*230  # Only mask out non-ref seqs
            out_f.write(f'>{seq_id}\n')
            out_f.write(f'{seq}\n')

    print(f'Trimmed and aligned sequences written to {input_fasta_path}')

# Prepare metadata file
# =====================
print("\nPreparing metadata...")
input_metadata_path = delphy_inputs_path / 'metadata.csv'
do_metadata = True
if input_metadata_path.exists():
    do_metadata = False
    print(f'[SKIPPING] {input_metadata_path} already exist')

if do_metadata:
    with open(input_metadata_path, 'w') as f:
        f.write('id,clade\n')
        f.write(f'{ref_acc_id},None\n')
        for sample_id in sorted(sample_id_2_clade.keys()):
            f.write(f'{sample_id},{sample_id_2_clade[sample_id]}\n')

    print(f'Metadata written to {input_metadata_path}')
