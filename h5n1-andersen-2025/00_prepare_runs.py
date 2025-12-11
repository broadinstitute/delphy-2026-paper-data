#!/usr/bin/env python3
#
# This script gathers all the cattle H5N1 sequences in the Andersen lab's avian-influenza repo,
# cross-references them to dates & geo tags in GenBank where available and applies some
# basic QA filters (excluding anything with a genotype different from B3.13).
# Then sequences are aligned segment by segment against the
# reference in `avian-influenza` (`A/cattle/Texas/24-008749-003/2024(H5N1)`).
# All segments are concatenated, longest to shortest, and Delphy fasta and metadata
# files are generated, both for all matching sequences as well as only those with
# exact dates.
#
# The reason for the two data sources is that SRA dates are very uncertain (most are "2024")
# and geos are a mess, while GenBank has almost complete date & clean geo info for
# sequences, **but with an apparent 3-month lag**

from Bio import Entrez
import csv
import os
from collections import defaultdict
import re
from pathlib import Path
import subprocess
import datetime
import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: ./00_prepare_runs.py <email-address-for-Entrez>\n")
    sys.exit(1)

# Config
# ======
Entrez.email = sys.argv[1]
path_to_mafft = Path("../mafft")

# Basic info about 2024+ H5N1-in-cattle outbreak in the USA
# =========================================================
refseq_genbank_name = 'A/cattle/Texas/24-008749-003/2024' # See README
refseq_genbank_name_in_genbank_mapping = 'A/cattle/TX/24-008749-003-original/2024' # Not sure why this changed

refseq_srr = 'SRR28752635'  # Found in metadata/genbank_mapping.tsv
dominant_genotype = 'B3.13'

# Segments ordered longest-to-shortest
# This is how NextStrain orders them, and it also coincides with GenBank segment numbers:
# genbank_seq         1      2      3     4     5     6     7     8
ordered_segments = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
segments_set = set(ordered_segments)

# Civics 101
us_state_abbrev_2_name = {
    'AL': 'Alabama',
    'AK': 'Alaska',
    'AZ': 'Arizona',
    'AR': 'Arkansas',
    'CA': 'California',
    'CO': 'Colorado',
    'CT': 'Connecticut',
    'DE': 'Delaware',
    'FL': 'Florida',
    'GA': 'Georgia',
    'HI': 'Hawaii',
    'ID': 'Idaho',
    'IL': 'Illinois',
    'IN': 'Indiana',
    'IA': 'Iowa',
    'KS': 'Kansas',
    'KY': 'Kentucky',
    'LA': 'Louisiana',
    'ME': 'Maine',
    'MD': 'Maryland',
    'MA': 'Massachusetts',
    'MI': 'Michigan',
    'MN': 'Minnesota',
    'MS': 'Mississippi',
    'MO': 'Missouri',
    'MT': 'Montana',
    'NE': 'Nebraska',
    'NV': 'Nevada',
    'NH': 'New Hampshire',
    'NJ': 'New Jersey',
    'NM': 'New Mexico',
    'NY': 'New York',
    'NC': 'North Carolina',
    'ND': 'North Dakota',
    'OH': 'Ohio',
    'OK': 'Oklahoma',
    'OR': 'Oregon',
    'PA': 'Pennsylvania',
    'RI': 'Rhode Island',
    'SC': 'South Carolina',
    'SD': 'South Dakota',
    'TN': 'Tennessee',
    'TX': 'Texas',
    'UT': 'Utah',
    'VT': 'Vermont',
    'VA': 'Virginia',
    'WA': 'Washington',
    'WV': 'West Virginia',
    'WI': 'Wisconsin',
    'WY': 'Wyoming',
}

# Helper routines
# ===============
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


# Prepare directory structure
# ===========================
raw_path = Path('raw')
raw_path.mkdir(parents=True, exist_ok=True)
scratch_path = Path('scratch')
scratch_path.mkdir(parents=True, exist_ok=True)
delphy_inputs_path = Path('delphy_inputs')
delphy_inputs_path.mkdir(parents=True, exist_ok=True)

# Clone Andersen lab's avian-influenza data repo
# ==============================================
# Also fix to commit e756a15 --- this repo updates often, and we want this benchmark to be reproducible
commit_hash = 'e756a15'
print(f"\nCloning Andersen lab's avian-influenza repo at commit {commit_hash}")
data_repo_path = raw_path / 'avian-influenza'
do_clone = True
if data_repo_path.exists():
    do_clone = False
    print(f'[SKIPPING] {data_repo_path} already exists')

if do_clone:
    subprocess.run([
        'git',
        'clone',
        'https://github.com/andersen-lab/avian-influenza.git'],
                   cwd=raw_path.as_posix())
    subprocess.run([
        'git',
        'checkout',
        commit_hash],
                   cwd=data_repo_path.as_posix())

# Map out SRRs & segments with sequences in Andersen lab repo
# ===========================================================
print(f"\nEnumerating SRRs and segments in Andersen lab repo")
fasta_dir_path = data_repo_path / 'fasta'
srr_2_seg_2_filename = defaultdict(dict)

pattern = re.compile(r'(SRR[0-9]+)_([A-Z0-9]+)_cns.fa')
for entry in os.scandir(fasta_dir_path.as_posix()):
    
    # Check that every file belongs to a specific SRR and segment, as in `pattern`
    # (as of this writing, all files under `fasta` match the above pattern)
    m = pattern.match(entry.name)
    if not m:
        raise ValueError(f'FASTA file does not match expected pattern: {entry.name}')
        
    srr, segment = m.groups()
    if segment in srr_2_seg_2_filename[srr]:
        raise ValueError(f'Two FASTA files for the same SRR and segment?: {entry.name}')

    srr_2_seg_2_filename[srr][segment] = entry.name

print(f'Found {len(srr_2_seg_2_filename)} unique SRRs under {fasta_dir_path}')

# Check that we have exactly the expected segments for each SRR
# (as of this writing, all SRRs have exactly the 8 segments shown above, no more, no less)
for srr, seg_2_filename in srr_2_seg_2_filename.items():
    if set(seg_2_filename.keys()) != segments_set:
        raise ValueError(f'SRR {srr} has these segments: {list(seg_2_filename.keys())}, expected these: {ordered_segments}')

# Read in SRA metadata
# ====================
sra_metadata_file_path = data_repo_path / 'metadata' / 'SraRunTable_automated.csv'
print(f"\nParsing SRA metadata in {sra_metadata_file_path}")

srr_2_metadata = defaultdict(dict)
with sra_metadata_file_path.open('r') as f:
    ff = csv.reader(f, delimiter=',')
    header = None
    for row in ff:
        if not header:
            header = row
            continue
        
        (srr, assayType, avgSpotLen, bases, bioProject, bioSample, bioSampleModel, theBytes,
         centerName, collectionDate, consent, datastoreFiletype, datastoreProvider, datastoreRegion,
         experiment, geo_loc_name_country, geo_loc_name_country_continent, geo_loc_name,
         host, instrument, isolate, libraryName, libraryLayout, librarySelection,
         librarySource, organism, platform, releaseDate, createDate, version, sampleName,
         sraStudy, serotype, isolation_source, biosampleAccession, is_retracted, retraction_detection_date) = row
        
        assert srr not in srr_2_metadata
        assert srr in srr_2_seg_2_filename, srr
        
        srr_2_metadata[srr] = {
            'collectionDate': collectionDate,
            'releaseDate': releaseDate,
            'geo_loc_name_country': geo_loc_name_country,
            'geo_loc_name_country_continent': geo_loc_name_country_continent,
            'geo_loc_name': geo_loc_name,
            'host': host,
            'releaseDate': releaseDate,
            'is_retracted': is_retracted == 'True'
        }

# Check that we have metadata for everything
# (as of this writing, we do)
assert len(srr_2_seg_2_filename) == len(srr_2_metadata)

# Filter out non-cattle samples right off the bat (they're more than half the data!)
# ==================================================================================
non_cattle_hosts = set()
non_cattle_srrs = []
for srr, metadata in srr_2_metadata.items():
    if metadata['host'].lower() != 'cattle':
        non_cattle_hosts.add(metadata['host'].lower())
        non_cattle_srrs.append(srr)

print(f'Dropping {len(non_cattle_srrs)} SRRs from {len(non_cattle_hosts)} non-cattle hosts, among which: {list(non_cattle_hosts)[:5]}')

for non_cattle_srr in non_cattle_srrs:
    del srr_2_seg_2_filename[non_cattle_srr]
    del srr_2_metadata[non_cattle_srr]

print(f'A total of {len(srr_2_seg_2_filename)} cattle-host SRRs remain')

# Filter off retracted samples too right off the bat (they're often responsible for apparent collisions later)
# ============================================================================================================
retracted_srrs = []
for srr, metadata in srr_2_metadata.items():
    if metadata['is_retracted']:
        retracted_srrs.append(srr)

print(f'Dropping {len(retracted_srrs)} SRRs because they are marked as retracted')

for retracted_srr in retracted_srrs:
    del srr_2_seg_2_filename[retracted_srr]
    del srr_2_metadata[retracted_srr]

print(f'A total of {len(srr_2_seg_2_filename)} cattle-host non-retracted SRRs remain')

# Read Genotype assignment
# ========================
genoflu_results_file_name = data_repo_path / 'metadata' / 'genoflu_results.tsv'
print(f"\nParsing Genoflu genotype assignments in {genoflu_results_file_name}")

srr_2_genotype = defaultdict(list)
with genoflu_results_file_name.open('r') as f:
    ff = csv.reader(f, delimiter='\t')
    header = None
    for row in ff:
        if not header:
            header = row
            continue

        (sra_run, sample_date, fasta_filename, genotype, *rest) = row

        if sra_run not in srr_2_seg_2_filename:
            continue   # We may have dropped this SRR already owning to non-cattle host

        if sra_run in srr_2_genotype:
            raise ValueError(f'SRR {sra_run} appears 2+ times in {genoflu_results_file_name} ?')
        
        srr_2_genotype[sra_run] = genotype

non_dominant_genotype_srrs = []
for srr, genotype in srr_2_genotype.items():
    if genotype != dominant_genotype:
        non_dominant_genotype_srrs.append(srr)

print(f'Dropping {len(non_dominant_genotype_srrs)} SRRs because they are not confidently of genotype "{dominant_genotype}"')

for non_dominant_genotype_srr in non_dominant_genotype_srrs:
    del srr_2_seg_2_filename[non_dominant_genotype_srr]
    del srr_2_metadata[non_dominant_genotype_srr]

print(f'A total of {len(srr_2_seg_2_filename)} cattle-host non-retracted SRRs with genotype {dominant_genotype} remain')


# Read SRR -> GenBank mappings
# ============================
# Read GenBank mappings (not available for every SRR)
genbank_mapping_file_name = data_repo_path / 'metadata' / 'genbank_mapping.tsv'
print(f"\nParsing GenBank mappings in {genbank_mapping_file_name}")

srr_2_genbank_mappings = defaultdict(list)
with genbank_mapping_file_name.open('r') as f:
    ff = csv.reader(f, delimiter='\t')
    header = None
    for row in ff:
        if not header:
            header = row
            continue
        
        (seg_file, seg_seq_name, sra_run, seg, genbank_acc, genbank_seq, genbank_name) = row

        if sra_run not in srr_2_seg_2_filename:
            continue   # We may have dropped this SRR already owning to non-cattle host
        
        # Don't bother gathering GenBank accession IDs for individual segments
        # Amazingly, different versions of the same SRR appear several times in GenBank with
        # different accession IDs and names!
        #
        # For example, the HA segment of SRR28752471 is
        #  - A/cattle/Michigan/24-009027-002/2024     -  PP824612.1
        #  - A/Cattle/Michigan/24-009027-002-v/2024   -  PQ011563.1
        #
        # First we gather all the names as sets, to see how bad the problem is

        # We gather the GenBank accession of the smallest segment to download
        # as little as possible from GenBank below
        if seg == ordered_segments[-1]:
            srr_2_genbank_mappings[sra_run].append((genbank_acc, genbank_name))

print(f'Number of SRRs with matching GenBank accessions: {len(srr_2_genbank_mappings)}')

# When there are multiple GenBank accessions for the same SRR & segment,
# disambiguate using the one with the lexicographically largest accession ID,
# in the hopes that this represents the latest (and most correct) accession
srr_2_genbank_name = {}   # Only apply mapping if it's unique
srr_2_ns_genbank_acc = {}
num_ambiguous_srrs = 0
warnings_on = True
for srr, genbank_mappings in srr_2_genbank_mappings.items():
    chosen_mapping = max(genbank_mappings)  # mappings are tuples (acc, name)

    # OVERRIDE for the reference sequence
    if srr == refseq_srr:
        chosen_mapping = [(acc, name) for (acc, name) in genbank_mappings if name == refseq_genbank_name_in_genbank_mapping][0]
    
    if len(genbank_mappings) > 1:
        if warnings_on:
            print(f'WARNING: {srr} NS segment has multiple GenBank entries, picking lexicographically latest one: {chosen_mapping}')
            print(f'          (available mappings: {genbank_mappings})')
        num_ambiguous_srrs += 1
        if warnings_on and num_ambiguous_srrs > 5:
            warnings_on = False
            print(f'...similar warnings suppressed henceforth')

    srr_2_ns_genbank_acc[srr], srr_2_genbank_name[srr] = chosen_mapping

print(f'Out of {len(srr_2_genbank_mappings)} SRRs that map to GenBank, {num_ambiguous_srrs} have an ambiguous mapping')

# Compile reverse mappings, which are useful below
# In passing, we check that mappings are 1-1, i.e., invertible
ns_genbank_acc_2_srr = {}
for srr, ns_genbank_acc in srr_2_ns_genbank_acc.items():
    assert ns_genbank_acc not in ns_genbank_acc_2_srr, ns_genbank_acc
    ns_genbank_acc_2_srr[ns_genbank_acc] = srr

genbank_name_2_srr = {}
for srr, genbank_name in srr_2_genbank_name.items():
    assert genbank_name not in genbank_name_2_srr, genbank_name
    genbank_name_2_srr[genbank_name] = srr


# Fetch sequences (and metadata!) from GenBank
# ============================================
print("\nFetching sequences and metadata from GenBank...")
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
        fetch_sequences(list(srr_2_ns_genbank_acc.values()), f)


# Extract dates and geos from GenBank where available
# ===================================================
print("\nParsing dates and geos from GenBank (where available)...")
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

def extract_gb_collection_dates_and_geo(in_gb_xml_f):
    result_dates = {}
    result_geos = {}
    
    records = Entrez.parse(in_gb_xml_f)
    for record in records:
        
        acc = record['GBSeq_accession-version']

        srr = ns_genbank_acc_2_srr[acc]
        
        # Extract collection date if it's there
        collection_date_range = None
        geo = None
        for feature in record['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'source':
                for qual in feature['GBFeature_quals']:
                    if qual['GBQualifier_name'] == 'collection_date':
                        raw_collection_date = qual['GBQualifier_value']
                        collection_date_range = parse_genbank_date(raw_collection_date)
                    if qual['GBQualifier_name'] == 'geo_loc_name':
                        geo = qual['GBQualifier_value']

        if collection_date_range is not None:
            result_dates[srr] = unparse_date_range(collection_date_range)
        if geo_loc_name is not None:
            result_geos[srr] = geo

    return (result_dates, result_geos)

with open(raw_genomes_path, 'rb') as gb_xml:
    srr_2_genbank_collection_date, srr_2_genbank_geo = extract_gb_collection_dates_and_geo(gb_xml)

num_srrs_with_gb_dates = len(srr_2_genbank_collection_date)
num_unique_dates = len(set(srr_2_genbank_collection_date.values()))
print(f'GenBank has dates for {num_srrs_with_gb_dates} SRRs ({num_unique_dates} distinct values)')

num_srrs_with_gb_full_dates = len(list(t
                                       for t in srr_2_genbank_collection_date.values()
                                       if len(t) == len('2024-01-01')))
num_unique_full_dates = len(set(t
                                for t in srr_2_genbank_collection_date.values()
                                if len(t) == len('2024-01-01')))
print(f'GenBank has full dates for {num_srrs_with_gb_full_dates} SRRs ({num_unique_full_dates} distinct values)')

num_srrs_with_gb_geos = len(srr_2_genbank_geo)
num_srrs_with_gb_geos_not_USA = len(list(g for g in srr_2_genbank_geo.values() if g != 'USA'))
num_unique_geos = len(set(srr_2_genbank_geo.values()))
print(f"GenBank has geos for {num_srrs_with_gb_geos} SRRs - of which {num_srrs_with_gb_geos_not_USA} SRRs's geo is not USA ({num_unique_geos} distinct values)")

print(f'Unique geos: {sorted(list(set(srr_2_genbank_geo.values())))}')


# Read in reference and sample sequences to prepare for segment-by-segment alignment
# ==================================================================================
assert srr_2_genbank_name[refseq_srr] == refseq_genbank_name_in_genbank_mapping
refseq_file_path = data_repo_path / 'reference' / 'reference.fasta'

print(f'\nReading reference genome in {refseq_file_path}')

seg_2_refseq = {}
for fasta_id, seq in read_fasta(refseq_file_path):
    m = re.match(r'([A-Z0-9]+)\|PP[0-9]+\.1\|'+refseq_genbank_name+r'\(H5N1\)', fasta_id)
    assert m, fasta_id
    seg = m.group(1)
    assert seg in segments_set, seg
    seg_2_refseq[seg] = seq

# Quick sanity checks on reference sequence
seg_2_reflen = {seg: len(refseq) for seg, refseq in seg_2_refseq.items()}
check_ordered_segs = sorted(list(seg_2_reflen.keys()), key=lambda seg: seg_2_reflen[seg], reverse=True)
assert ordered_segments == check_ordered_segs

print(f'Segment lengths in reference: ')
for seg in ordered_segments:
    print(f'- {seg}: {seg_2_reflen[seg]} bases')

# Prepare final IDs, dates and geos for Delphy FASTA and metadata files
# =====================================================================
def sanitize_genbank_name(raw_gb_name):
    # Fasta IDs cannot contain spaces!
    return raw_gb_name.replace(' ', '_')

srr_2_final_long_id = {}
srr_2_final_short_id = {}
for srr, metadata in srr_2_metadata.items():
    # Ideal: 'A/Cattle/blah/blah/blah|SRR12345678'  (among included SRRs, GenBank names are unique!)
    # When no GenBank name: 'SRR12345678'
    if srr in srr_2_genbank_name:
        srr_2_final_long_id[srr] = f'{sanitize_genbank_name(srr_2_genbank_name[srr])}|{srr}'
        srr_2_final_short_id[srr] = f'{sanitize_genbank_name(srr_2_genbank_name[srr])}'
    else:
        srr_2_final_long_id[srr] = srr
        srr_2_final_short_id[srr] = srr

# Assign best estimate for dates (may need to think harder here)
srr_2_final_date = {}
for srr, metadata in srr_2_metadata.items():
    if srr in srr_2_genbank_collection_date:   # GenBank has much higher collection date coverage than SRA
        srr_2_final_date[srr] = srr_2_genbank_collection_date[srr]
    else:
        sra_date = metadata['collectionDate']
        assert sra_date != 'missing', srr  # This *could* happen, but currently doesn't
        srr_2_final_date[srr] = sra_date  # Could refine to be earlier than release date

# Assign best estimate for geos (may need to think harder here)
srr_2_final_geo = {}

genbank_prefix_2_geo_guess = {}
for abbrev, name in us_state_abbrev_2_name.items():
    genbank_prefix_2_geo_guess[f'A/cattle/{name}/'] = f'USA: {abbrev}'

for srr, metadata in srr_2_metadata.items():
    geo_guess = None
    if srr in srr_2_genbank_geo:   # GenBank has much cleaner geo coverage than SRA
        geo_guess = srr_2_genbank_geo[srr]
    else:
        sra_geo = metadata['geo_loc_name']
        geo_guess = sra_geo  # Could refine to be earlier than release date

    if geo_guess == 'USA':
        # Try to guess from GenBank name if available
        if srr in srr_2_genbank_name:
            genbank_name = srr_2_genbank_name[srr]
            for prefix, prefix_geo_guess in genbank_prefix_2_geo_guess.items():
                if genbank_name.lower().startswith(prefix.lower()):
                    geo_guess = prefix_geo_guess
                    break

    srr_2_final_geo[srr] = geo_guess
                    
    
        
# Prepare unaligned FASTAs for each segment
# =========================================
for seg in ordered_segments:
    scratch_unaligned_genomes_fasta_path = scratch_path / f'raw_genomes_{seg}.fasta'
    scratch_ref_fasta_path = scratch_path / f'ref_{seg}.fasta'
    print(f'Preparing unaligned combined FASTA file for segment {seg} in {scratch_unaligned_genomes_fasta_path}'
          + f' (ref genome in {scratch_ref_fasta_path})')

    if scratch_unaligned_genomes_fasta_path.exists() and scratch_ref_fasta_path.exists():
        print(f'[SKIPPING] - both files exist already')
        continue
    
    with scratch_ref_fasta_path.open('w') as fasta:
        fasta.write(f'>{srr_2_final_long_id[refseq_srr]}|{srr_2_final_date[refseq_srr]}\n')
        fasta.write(seg_2_refseq[seg])

    with scratch_unaligned_genomes_fasta_path.open('w') as fasta:
        for srr in srr_2_seg_2_filename.keys():
            # Ref seq will come first because of how mafft works, no need to duplicate it here
            if srr == refseq_srr:
                continue
            fasta.write(f'>{srr_2_final_long_id[srr]}|{srr_2_final_date[srr]}\n')
            srr_fasta_file_path = data_repo_path / 'fasta' / f'{srr}_{seg}_cns.fa'
            srr_fasta = read_fasta(srr_fasta_file_path)
            assert len(srr_fasta) == 1
            fasta_id, seq = srr_fasta[0]
            fasta.write(f'{seq}\n')

# Aligned sequences for each segment (time-consuming!)
# ====================================================
for seg in ordered_segments:
    scratch_aligned_fasta_path = scratch_path / f'aligned_{seg}.fasta'
    scratch_raw_genomes_fasta_path = scratch_path / f'raw_genomes_{seg}.fasta'
    scratch_ref_fasta_path = scratch_path / f'ref_{seg}.fasta'
    print(f"\nAligning sequences for {seg} segment with mafft in {scratch_aligned_fasta_path}")
    
    if scratch_aligned_fasta_path.exists():
        print(f'[SKIPPING] {scratch_aligned_fasta_path} already exist')
        continue

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

# Make final concatenated alignments
# ==================================

# Exclude clear outliers from previous runs (visual inspection)
excluded_ids = set([
    # All samples previously on this list were D1.1 samples from Nevada (2025-06-02)

    # These two CA sequences are supposedly dated to late Jan 2024, but they are separated
    # from the nearest sequences by very many SNPs.  Including them single-handedly pushes the
    # tMRCA to ~ Oct 2023 instead of late Dec 2023, and implies that the California outbreak of
    # mid-late 2024 descends from the root via a ~6-month long branch with no other descendants.
    # All of the above suggests these two sequences have a data problem, so we exclude them.
    'A/cattle/CA/24-035207-001-original/2024',
    'A/cattle/CA/24-035207-002-original/2024',
])

aligned_all_fasta_path = delphy_inputs_path / f'h5n1-andersen-{commit_hash}-ALL.fasta'

print(f"\nConcatenating all aligned segments into {aligned_all_fasta_path}")
if aligned_all_fasta_path.exists():
    print(f'[SKIPPING] File exists')
else:
    seg_fastas = [read_fasta(scratch_path / f'aligned_{seg}.fasta') for seg in ordered_segments]
    with open(aligned_all_fasta_path, 'w') as f:
        for (i,per_seg_id_seq_pairs) in enumerate(zip(*seg_fastas)):
            first_id_line = per_seg_id_seq_pairs[0][0]
            assert all(id_line == first_id_line for (id_line,seg_seq) in per_seg_id_seq_pairs)

            # No amount of heuristics here seems to save the runs with a large fraction
            # of tips with very uncertain dates.  So best not to introduce any arbitrary
            # whatsoever
            #
            #dateStr = first_id_line.strip().split('|')[-1]
            #if len(dateStr) == 4 and dateStr == '2025':
            #    # Cap uncertain date at release date if that improves things
            #    srr = first_id_line.strip().split('|')[-2]
            #    releaseDate = srr_2_metadata[srr]['releaseDate']
            #    if dateStr == releaseDate[:4]:
            #        first_id_line = f"{first_id_line[:first_id_line.rindex('|')]}|{dateStr}-01-01/{releaseDate[:10]}"
                
            if first_id_line.split('|')[0] in excluded_ids:
                print(f'Removing gross outlier {first_id_line.strip()} after visual inspection')
            else:
                f.write(f'>{first_id_line.strip()}\n')
                f.write(f"{''.join(seg_seq.strip() for (id_line,seg_seq) in per_seg_id_seq_pairs)}\n")

                
aligned_all_full_dates_only_fasta_path = delphy_inputs_path / f'h5n1-andersen-{commit_hash}-ALL_full_dates_only.fasta'
print(f"\nConcatenating all aligned segments (only seqs with full dates) into {aligned_all_full_dates_only_fasta_path}")
if aligned_all_full_dates_only_fasta_path.exists():
    print(f'[SKIPPING] File exists')
else:
    seg_fastas = [read_fasta(scratch_path / f'aligned_{seg}.fasta') for seg in ordered_segments]
    with open(aligned_all_full_dates_only_fasta_path, 'w') as f:
        for (i,per_seg_id_seq_pairs) in enumerate(zip(*seg_fastas)):
            first_id_line = per_seg_id_seq_pairs[0][0]
            assert all(id_line == first_id_line for (id_line,seg_seq) in per_seg_id_seq_pairs)
    
            dateStr = first_id_line.strip().split('|')[-1]
            if len(dateStr) == len('2024-01-01'):
                if first_id_line.split('|')[0] in excluded_ids:
                    print(f'Removing gross outlier {first_id_line.strip()} after visual inspection')
                else:
                    f.write(f'>{first_id_line.strip()}\n')
                    f.write(f"{''.join(seg_seq.strip() for (id_line,seg_seq) in per_seg_id_seq_pairs)}\n")

# Prepare metadata file
# =====================
print("\nPreparing metadata...")
input_metadata_path = delphy_inputs_path / f'h5n1-andersen-{commit_hash}_metadata.csv'
do_metadata = True
if input_metadata_path.exists():
    do_metadata = False
    print(f'[SKIPPING] {input_metadata_path} already exist')

# The following might need to be expanded if other variations on US states names start showing up in the metadata
geos_2_clean_geos = {}
geos_2_clean_geos['USA'] = '-'
for abbrev, name in us_state_abbrev_2_name.items():
    geos_2_clean_geos[f'USA: {abbrev}'] = abbrev
    geos_2_clean_geos[f'USA:{abbrev}'] = abbrev

unmapped_geos = set()

if do_metadata:
    with open(input_metadata_path, 'w') as f:
        f.write('id,geo\n')
        for srr in srr_2_metadata.keys():
            shortId = srr_2_final_short_id[srr]
            raw_geo = srr_2_final_geo.get(srr, '-')

            if raw_geo not in geos_2_clean_geos:
                unmapped_geos.add(raw_geo)
                geo = '???'
            else:
                geo = geos_2_clean_geos[raw_geo]
            
            f.write(f'{shortId},{geo}\n')

    print(f'Metadata written to {input_metadata_path}')

    assert len(unmapped_geos) == 0, f'Unmapped geos: {sorted(list(unmapped_geos))}'
