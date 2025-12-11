#!/usr/bin/env python

from Bio import Entrez
import datetime
import numpy as np
import scipy as sp
import xml.etree.ElementTree as ET
import re
import shutil
import sys
from pathlib import Path

if len(sys.argv) != 3:
    sys.stderr.write("Usage: ./00_prepare_runs.py <email-address-for-Entrez> <path_to_1259657_file_s3.zip>\n")
    sys.stderr.write("  NOTE: you may need to manually download the original data from https://www.science.org/doi/10.1126/science.1259657, supplementary file S3 (1259657_file_s3.zip)\n");
    sys.exit(1)

# Config
# ======
Entrez.email = sys.argv[1]
path_to_S3_zip = sys.argv[2]

# Prepare directory structure
# ===========================
raw_path = Path('raw')
raw_path.mkdir(parents=True, exist_ok=True)
scratch_path = Path('scratch')
scratch_path.mkdir(parents=True, exist_ok=True)
delphy_inputs_path = Path('delphy_inputs')
delphy_inputs_path.mkdir(parents=True, exist_ok=True)

# Unpack S3 zip file and get BEAST file
# =====================================
print("\nUnpacking supplementary data zip file...")
beast_file_path = scratch_path / 'beast' / '2014_GN.SL_SRD.HKY_strict_ctmc.exp.xml'
do_unzip = True
if beast_file_path.exists():
    do_unzip = False
    print(f'[SKIPPING] {beast_file_path} already exists')

if do_unzip:
    shutil.unpack_archive(path_to_S3_zip, scratch_path)
    if not beast_file_path.exists():
        raise ValueError(f'{beast_file_path} not in supplementary info zip file?')

# Extract sequences from BEAST file
# =================================
# We'll extract the essentials of the run in `2014_GN.SL_SRD.HKY_strict_ctmc.exp.xml`.  The original data is split into two partitions (genic and intergenic), but this particular BEAST run treats both partitions together.  However, the relation to site indices of the reference genome used in the paper (KJ660346) is totally scrambled.

print("\nReading sequence data from BEAST XML file...")
beastXml = ET.parse(beast_file_path.as_posix())
xmlTaxa = [taxon for taxon in beastXml.find("taxa")]
xmlAlignments = beastXml.findall("alignment")
print(f'Found {len(xmlAlignments)} alignments and {len(xmlTaxa)} taxa')

# Check that taxa are listed in the same order in <taxa> section and in alignments
taxaIds = [taxon.get('id') for taxon in xmlTaxa]
alignment1IdRefs = [s.find('taxon').get('idref') for s in xmlAlignments[0].findall('sequence')]
alignment2IdRefs = [s.find('taxon').get('idref') for s in xmlAlignments[1].findall('sequence')]
if taxaIds != alignment1IdRefs:
    raise ValueError("Taxa listed in a different order in <taxa> section and in first alignment?")
if taxaIds != alignment2IdRefs:
    raise ValueError("Taxa listed in a different order in <taxa> section and in second alignment?")

# We try to unscramble the genic vs intergenic partition by aligning everything (manually) to the KJ660346 reference
# genome used in Fig 4 of the paper (in the end, that should make it so that the reference to "position 10,218" in the
# paper maps with position 10,218 here).
#
# From the GenBank entry for KJ660346 (https://www.ncbi.nlm.nih.gov/nuccore/KJ660346), that sequence has 18959 bases,
# and its genic regions (`CDS`) are as follows:
#
# * NP: 470..2689
# * VP35: 3129..4151
# * VP40: 4479..5459
# * GP: 6039..8068
# * VP30: 8509..9375
# * VP24: 10345..11100
# * L: 11581..18219
#
# KJ660346 also appears in the original BEAST file (`EBOV|KJ660346|Kissidougou-C15|Kissidougou_Guinea|2014-03-17`).
# The partitioning process extracted all the genic regions above and concatenated in the first partition, and the
# remaining bases are joined into the second partition.  The actual alignment has a few tiny surprises, but nothing
# showstopping.  See `do_reassembly` below for full details.

# Download GenBank file into KJ660346.2
# =====================================
def fetch_sequences(ids, out_f):
    with Entrez.efetch(db="nucleotide", id=','.join(ids), rettype="gb", retmode="xml") as handle:
        while True:
            lines = handle.readlines(1000)
            if not lines:
                break
            out_f.write(b''.join(lines))

print("\nFetching KJ660346.2 from GenBank...")
raw_reference_path = raw_path / 'KJ660346.2.gb'
do_fetch_ref = True
if raw_reference_path.exists():
    do_fetch_ref = False
    print(f'[SKIPPING] {raw_reference_path} already exists')

ref_acc_id = 'KJ660346.2'

if do_fetch_ref:
    with open(raw_reference_path, 'wb') as f:
        fetch_sequences([ref_acc_id], f)

# Convert GenBank file to FASTA
# =============================
scratch_ref_fasta_path = scratch_path / 'ref.fasta'
do_convert = True
if scratch_ref_fasta_path.exists():
    do_convert = False
    print(f'[SKIPPING] {scratch_ref_fasta_path} already exist')

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
        m = re.search('Assembly Name :: (\S+) ;', record['GBSeq_comment'])
        assembly_name = default_assembly_name
        if m:
            assembly_name = m[1]
        
        seq = record['GBSeq_sequence']
        
        out_fasta_f.write(f'>{assembly_name}|{acc}|{unparse_date_range(collection_date_range)}\n')
        out_fasta_f.write(seq)
        out_fasta_f.write('\n')

if do_convert:
    with open(raw_reference_path, 'rb') as gb_xml:
        with open(scratch_ref_fasta_path, 'w') as fasta:
            gb_to_fasta(gb_xml, fasta, ref_acc_id)

# Prepare fasta with all the sequences
# ====================================
print("\nReconstituting multiple sequence alignment from 2 partitions in BEAST XML file...")

input_fasta_path = delphy_inputs_path / 'ebola.fasta'
do_reassembly = True
if input_fasta_path.exists():
    do_reassembly = False
    print(f'[SKIPPING] {input_fasta_path} already exist')

if do_reassembly:
    # The following map was manually rebuilt by carefully realigning the GenBank reference sequence (GB)
    # to each of the two partitions (genic = P1, intergenic = P2):
    #
    # GB: [    1,   470) -> P2: [    1,   470): CGGACACACA...AATTCCGAGT
    # GB: [  470,  2690) -> P1: [    1,  2221): ATGGATTCTC...TCATCAGTGA
    # GB: [ 2690,  3129) -> P2: [  470,   909): ATGAGCATGT...GCCTAACAAG
    # GB: [ 3129,  4152) -> P1: [ 2221,  3244): ATGACAACTA...CAAAATTTGA
    # GB: [ 4152,  4479) -> P2: [  909,  1236): GCCAATCTCT...TGTTAAAAAT
    # GB: [ 4479,  5460) -> P1: [ 3244,  4225): ATGAGGCGGG...TGAGAAGTAA
    # GB: [ 5460,  6039) -> P2: [ 1236,  1815): TTGCAATAAT...CGACAACACA
    # GB: [ 6039,  6924) -> P1: [ 4225,  5110): ATGGGTGTTA...AACTAAAAAA
    # GB: [ 6924,  8068) -> P1: [ 5111,  6255): ACCTCACTAG...TTGTCTTTTA
    # GB: [ 8068,  8069) -> Gap
    # GB: [ 8069,  8509) -> P2: [ 1815,  2255): TCTTTCTTCA...CAACTCTTAA
    # GB: [ 8509,  9376) -> P1: [ 6256,  7123): ATGGAAGCTT...TACCCCTTAA
    # GB: [ 9376, 10345) -> P2: [ 2255,  3224): TAAGGCTGAC...AGAAAAAACC
    # GB: [10345, 11101) -> P1: [ 7123,  7879): ATGGCCAAAG...TGCTATCTAA
    # GB: [11101, 11581) -> P2: [ 3224,  3704): CTAAGATGGA...TTGATATTAA
    # GB: [11581, 18220) -> P1: [ 7879, 14518): ATGGCTACAC...GTTCGATTGA
    # GB: [18220, 18960) -> P2: [ 3704,  4444): ATAACCGTGC...TGTGTGTCCA
    
    with open(input_fasta_path, 'w') as f:
        for (taxon, s1, s2) in zip(xmlTaxa,
                                   xmlAlignments[0].findall('sequence'),
                                   xmlAlignments[1].findall('sequence')):
            taxonId = taxon.get('id')
            _, gbId, internalId, geo, date = taxonId.split('|')
            p1 = re.sub(r"\s+", "", ''.join(s1.itertext()))
            p2 = re.sub(r"\s+", "", ''.join(s2.itertext()))
            
            # Now realign this to KJ660346 according to the above map
            finalSeq = ''.join([
                p2[    1-1:   470-1],  # [    1,   470)
                p1[    1-1:  2221-1],  # [  470,  2690)
                p2[  470-1:   909-1],  # [ 2690,  3129)
                p1[ 2221-1:  3244-1],  # [ 3129,  4152)
                p2[  909-1:  1236-1],  # [ 4152,  4479)
                p1[ 3244-1:  4225-1],  # [ 4479,  5460)
                p2[ 1236-1:  1815-1],  # [ 5460,  6039)
                p1[ 4225-1:  5110-1],  # [ 6039,  6924)
                p1[ 5111-1:  6255-1],  # [ 6924,  8068)
                "N"*1,                 # [ 8068,  8069)
                p2[ 1815-1:  2255-1],  # [ 8069,  8509)
                p1[ 6256-1:  7123-1],  # [ 8509,  9376)
                p2[ 2255-1:  3224-1],  # [ 9376, 10345)
                p1[ 7123-1:  7879-1],  # [10345, 11101)
                p2[ 3224-1:  3704-1],  # [11101, 11581)
                p1[ 7879-1: 14518-1],  # [11581, 18220)
                p2[ 3704-1:  4444-1],  # [18220, 18960)
            ])
            
            f.write('>')
            f.write(f'{gbId}-{internalId}|{date}')
            f.write('\n')
            f.write(finalSeq)
            f.write('\n')
            
    print(f'Reconstituted multiple sequence alignment written to {input_fasta_path}')

# Prepare metadata file
# =====================
print("\nPreparing metadata...")
input_metadata_path = delphy_inputs_path / 'ebola_metadata.csv'
do_metadata = True
if input_metadata_path.exists():
    do_metadata = False
    print(f'[SKIPPING] {input_metadata_path} already exist')

if do_metadata:
    with open(input_metadata_path, 'w') as f:
        f.write('id,Geo\n')
        for taxon in xmlTaxa:
            taxonId = taxon.get('id')
            _, gbId, internalId, geo, date = taxonId.split('|')
            
            f.write(f'{gbId}-{internalId},{geo}\n')

    print(f'Metadata written to {input_metadata_path}')
