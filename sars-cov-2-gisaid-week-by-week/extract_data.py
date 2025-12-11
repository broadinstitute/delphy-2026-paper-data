#!/usr/bin/env python3

import argparse
import lzma
import csv
import re
import gzip
import datetime
from collections import defaultdict
from pathlib import Path
import subprocess
import random

# Config
# ======
MODE_SUBMITTED_BY_DATE = 'submissionDate'
MODE_COLLECTED_BY_DATE = 'collectionDate'

parser = argparse.ArgumentParser(description="Prepare Delphy runs for first N weeks of SARS-CoV-2")
parser.add_argument('--mode', help='Include samples in run if _this_ date is on or before CDC epi week N',
                    required=True, choices=[MODE_SUBMITTED_BY_DATE, MODE_COLLECTED_BY_DATE])

args = parser.parse_args()
    
path_to_mafft = Path("../mafft")

if args.mode == MODE_SUBMITTED_BY_DATE:
    min_submission_date = datetime.datetime(2019, 12,  1)
    max_submission_date = datetime.datetime(2020,  3, 28)  # end of CDC epi week 2020-13
    
    min_collection_date = datetime.datetime(2019, 12,  1)
    max_collection_date = max_submission_date
    
elif args.mode == MODE_COLLECTED_BY_DATE:
    min_submission_date = datetime.datetime(2019, 12,  1)
    max_submission_date = datetime.datetime(2024, 12, 31)  # Allow anything!

    min_collection_date = datetime.datetime(2019, 12,  1)
    max_collection_date = datetime.datetime(2020,  3,  7)  # end of CDC epi week 2020-10


# The tail ends of a genome are hard to sequence.  These crude parameters are from Lemieux et al 2021.
# To mimic the lack of knowledge at the beginning of the outbreak, we do not mask any other sites
num_initial_masked_sites = 268
num_final_masked_sites = 230

# And we filter out fragment submissions if they weren't filtered out before
min_seq_len = 20000




metadata = {}  # Key: Accession ID, Value: {virusName, collectionDate, location, length, host, submissionDate}
virusName2AccessionId = {}

def sanitizeTabs(s):
    return s.replace("\t", " ")  # Avoid CSV pain

with open('metadata_20200331.tsv', 'r') as f:
    ff = csv.reader(f, delimiter='\t')
    header = None
    for row in ff:
        if not header:
            header = row
            continue
        
        virusName, accessionId, collectionDate, location, additionalLocationInfo, seqLength, host, submissionDate = row
        assert accessionId not in metadata;
        metadata[accessionId] = {
            'virusName': sanitizeTabs(virusName),
            'collectionDate': sanitizeTabs(collectionDate),
            'location': sanitizeTabs(location),
            'length': int(seqLength),
            'host': sanitizeTabs(host),
            'submissionDate': sanitizeTabs(submissionDate)
        }
        virusName2AccessionId[virusName] = accessionId

# Filter out anything with an imprecise collectionDate or a date outside [min_collection_date, max_collection_date]
dateRE = re.compile(r'(\d\d\d\d)\-(\d\d)\-(\d\d)')
def isCollectionDateGood(collectionDate):
    if len(collectionDate) != len('YYYY-MM-DD'):
        return False
    m = dateRE.match(collectionDate)
    if not m:
        return False

    dt = datetime.datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)))
    
    return min_collection_date <= dt <= max_collection_date
        
# Filter out anything with an imprecise submissionDate or a date outside [min_submission_date, max_submission_date]
def isSubmissionDateGood(submissionDate):
    if len(submissionDate) != len('YYYY-MM-DD'):
        return False
    m = dateRE.match(submissionDate)
    if not m:
        return False

    dt = datetime.datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)))
    
    return min_submission_date <= dt <= max_submission_date

# Produce a set of accession IDs which pass very basic filters
interestingAccessionIDs = set()
for accessionId, info in metadata.items():
    if (info['host'] == 'Human' and
        info['length'] >= min_seq_len and
        isCollectionDateGood(info['collectionDate']) and
        isSubmissionDateGood(info['submissionDate'])):
        
        interestingAccessionIDs.add(accessionId)

# Filter out a few known-bad genomes (including them completely distorts the tree)
#
# These were determined by iteratively refining short initial Delphy runs.  Inclusion of one of these sequences
# causes the tMRCA to shoot way back to early 2019 or before, with the offending sample being a single child
# branch with O(10)-O(100) mutations to the sequence.  Many of these were noted at the time (references).
#
# The method here is purposely crude to mimic what would probably be done in the early stages of an outbreak,
# where quality controls, protocol shortcomings and known-difficult sites have not yet been mapped out.
# We do not mask any sites except the tail ends for the same reason.
badAccessionIDs = [
        # References:
        # [1] https://virological.org/t/temporal-signal-and-the-evolutionary-rate-of-2019-n-cov-using-47-genomes-collected-by-feb-01-2020/379 from 3 Feb 2020
        # [2] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7142683/ (Bal et al 2020)
        # [3] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7228400/ (Wang et al 2020)
        # [4] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8769011/ (Lv et al 2021)
        # [5] https://raw.githubusercontent.com/roblanf/sarscov2phylo/master/excluded_sequences.tsv
        # [6] https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/exclude.txt
        # ** = Marked as "Under Investigation" in GISAID
        # NS-YYYY-MM-DD = Name and date in which sequence appeared in NextStrain exclusions list

        # Submitted in epi week 5
        'EPI_ISL_406592',  # [1,2,5,6] - hCoV-19/Guangdong/SZTH-001/2020 = Shenzhen/SZTH-001/2020 in NS-2020-02-17

        # Submitted in epi week 7
        'EPI_ISL_408483',  # [2,5,6] - ** - hCoV-19/Shanghai/IVDC-SH-001/2020 = Shanghai/IVDC-SH-001/2020 in NS-2020-02-09
        'EPI_ISL_408487',  # [2,5,6] - ** - hCoV-19/Henan/IVDC-HeN-002/2020 = Henan/IVDC-HeN-002/2020 in NS-2020-02-09

        # Submitted in epi week 9
        'EPI_ISL_412900',  # [2,5,6] - ** - hCoV-19/Wuhan/HBCDC-HB-04/2019 = Wuhan/HBCDC-HB-04/2019 in NS-2020-06-08

        # Submitted in epi week 10
        'EPI_ISL_406595',  # [1,2,5,6] - hCoV-19/Guangdong/SZTH-004/2020 = Shenzhen/SZTH-004/2020 in NS-2020-02-17

        # Submitted in epi week 11
        'EPI_ISL_416030',  # NONE - ** - hCoV-19/Saudi Arabia/CV1249/2020
        'EPI_ISL_414935',  # [5] - ** - hCoV-19/Shandong/LY002-2/2020
        'EPI_ISL_527871',  # [6] - ** - hCoV-19/Yunnan/KMS1/2020 = Yunnan/KMS1/2020 in NS-2021-01-09
        'EPI_ISL_414934',  # [5,6] - ** - hCoV-19/Shandong/LY001-2/2020 = Shandong/LY001-2/2020 in NS-2021-04-03
        'EPI_ISL_415593',  # [5,6] - ** - hCoV-19/USA/WA-UW65/2020 = USA/WA-UW65/2020 in NS-2020-06-08

        # Submitted in epi week 12
        'EPI_ISL_417919',  # [5] - ** - hCoV-19/Malaysia/186197/2020 = Malaysia/186197/2020 in NS-2020-06-08
        'EPI_ISL_417413',  # [5] - ** - hCoV-19/Turkey/6224-Ankara1034/2020 = Turkey/6224-Ankara1034/2020 in NS-2020-06-08

        # Collected by epi week 3
        # Cluster of Lebanese sequences probably incorrectly annotated with collection date 2020-01-10
        'EPI_ISL_8311702',  # NONE - hCoV-19/Lebanon/LAU-R107/2020
        'EPI_ISL_8311703',  # [6] - hCoV-19/Lebanon/LAU-R126/2020 - Lebanon/LAU-R126/2020 in NS-2023-06-09
        'EPI_ISL_8311704',  # [6] - hCoV-19/Lebanon/LAU-R128/2020 - Lebanon/LAU-R128/2020 in NS-2022-09-18
        'EPI_ISL_8311705',  # [6] - hCoV-19/Lebanon/LAU-R130/2020 - Lebanon/LAU-R130/2020 in NS-2022-09-18
        'EPI_ISL_8311706',  # [6] - hCoV-19/Lebanon/LAU-R135/2020 - Lebanon/LAU-R135/2020 in NS-2022-08-23
        'EPI_ISL_8311707',  # [6] - hCoV-19/Lebanon/LAU-R137/2020 - Lebanon/LAU-R137/2020 in NS-2022-09-18
        'EPI_ISL_8311708',  # [6] - hCoV-19/Lebanon/LAU-R150/2020 - Lebanon/LAU-R150/2020 in NS-2022-06-21
        'EPI_ISL_8311709',  # [6] - hCoV-19/Lebanon/LAU-R158/2020 - Lebanon/LAU-R158/2020 in NS-2022-09-18
        'EPI_ISL_8311711',  # [6] - hCoV-19/Lebanon/LAU-R178/2020 - Lebanon/LAU-R178/2020 in NS-2023-06-09
        'EPI_ISL_8311729',  # [6] - hCoV-19/Lebanon/LAU-R289/2020 - Lebanon/LAU-R289/2020 in NS-2022-09-18
        'EPI_ISL_8311730',  # [6] - hCoV-19/Lebanon/LAU-R293/2020 - Lebanon/LAU-R293/2020 in NS-2022-09-18
        'EPI_ISL_8311731',  # [6] - hCoV-19/Lebanon/LAU-R294/2020 - Lebanon/LAU-R294/2020 in NS-2022-10-24
        'EPI_ISL_8311733',  # [6] - hCoV-19/Lebanon/LAU-R328/2020 - Lebanon/LAU-R328/2020 in NS-2022-10-06
        'EPI_ISL_8311734',  # [6] - hCoV-19/Lebanon/LAU-R330/2020 - Lebanon/LAU-R330/2020 in NS-2022-09-18
        'EPI_ISL_8311741',  # [6] - hCoV-19/Lebanon/LAU-R350/2020 - Lebanon/LAU-R350/2020 in NS-2023-06-17
        'EPI_ISL_8311755',  # [6] - hCoV-19/Lebanon/LAU-R84/2020 - Lebanon/LAU-R84/2020 in NS-2022-09-18
        'EPI_ISL_8311756',  # [6] - hCoV-19/Lebanon/LAU-R89/2020 - Lebanon/LAU-R89/2020 in NS-2022-09-18
        'EPI_ISL_8311758',  # [6] - hCoV-19/Lebanon/LAU-R97/2020 - Lebanon/LAU-R97/2020 in NS-2022-09-18
        'EPI_ISL_8311766',  # NONE - ** - hCoV-19/Lebanon/LAU-R101/2020
        'EPI_ISL_8311767',  # NONE - ** - hCoV-19/Lebanon/LAU-R146/2020
        'EPI_ISL_8311771',  # [6] - ** - hCoV-19/Lebanon/LAU-R104/2020 - Lebanon/LAU-R104/2020 in NS-2022-09-18
        'EPI_ISL_8317204',  # [6] - hCoV-19/Lebanon/LAU-R99/2020 - Lebanon/LAU-R99/2020 in NS-2023-06-18
        'EPI_ISL_8317205',  # [6] - ** - hCoV-19/Lebanon/LAU-R95/2020 - Lebanon/LAU-R95/2020 in NS-2022-06-03
        'EPI_ISL_8317209',  # NONE - ** - hCoV-19/Lebanon/LAU-R295/2020
        'EPI_ISL_8317210',  # [6] - ** - hCoV-19/Lebanon/LAU-R351/2020 - Lebanon/LAU-R351/2020 in NS-2023-06-23
        
        'EPI_ISL_2671842',  # [6] - ** - hCoV-19/Japan/20200409-129/2020 - Japan/20200409-129/2020 in NS-2021-06-25
        'EPI_ISL_2716627',  # [6] - ** - hCoV-19/Sierra Leone/SL10/2020 - SierraLeone/SL10/2020 in NS-2021-06-01
        'EPI_ISL_2716636',  # [6] - ** - hCoV-19/Sierra Leone/SL21/2020 - SierraLeone/SL21/2020 in NS-2021-06-01
        
        'EPI_ISL_4405694',  # [6] - hCoV-19/Argentina/PAIS-A1026/2020 - Argentina/PAIS-A1026/2020 in NS-2021-10-01
        'EPI_ISL_7955525',  # [6] - ** - hCoV-19/Argentina/PAIS-C0160/2020 - Argentina/PAIS-C0160/2020 in NS-2023-06-23
        'EPI_ISL_11575530', # [6] - hCoV-19/USA/CO_c8ftxs_0117/2020 - USA/CO_c8ftxs_20200117/2020 in NS-2023-06-18

        'EPI_ISL_14307752', # NONE - hCoV-19/USA/MT-UMGC-02796/2020
        'EPI_ISL_3804266',  # NONE - hCoV-19/Niger/16250A/2020  (collected date is 2020-01-01; maybe meant 2020-*-* ?)

        # Collected by epi week 4
        'EPI_ISL_2835566',  # [6] - hCoV-19/USA/CA-SEARCH-103149/2020 - USA/CA-SEARCH-103149/2020 in NS-2021-06-09
        'EPI_ISL_17116757', # [6] - hCoV-19/USA/WV-WV064576/2020 - USA/WV064576/2020 in NS-2021-09-01
        'EPI_ISL_17116755', # [6] - ** - hCoV-19/USA/WV-WV064569/2020 - USA/WV064569/2020 in NS-2021-09-01
        'EPI_ISL_19045047', # NONE - ** - hCoV-19/USA/CA-SEARCH-139716/2020
        'EPI_ISL_19037038', # NONE - ** - hCoV-19/Ghana/CRI-01/2020
        'EPI_ISL_19037039', # NONE - ** - hCoV-19/Ghana/CRI-02/2020
        'EPI_ISL_19037040', # NONE - ** - hCoV-19/Ghana/CRI-03/2020
        'EPI_ISL_19037041', # NONE - ** - hCoV-19/Ghana/CRI-04/2020
        'EPI_ISL_19037042', # NONE - ** - hCoV-19/Ghana/CRI-05/2020
        'EPI_ISL_19037043', # NONE - ** - hCoV-19/Ghana/CRI-06/2020
        'EPI_ISL_19037044', # NONE - ** - hCoV-19/Ghana/CRI-08/2020
        'EPI_ISL_19037045', # NONE - ** - hCoV-19/Ghana/CRI-09/2020
        'EPI_ISL_19037046', # NONE - ** - hCoV-19/Ghana/CRI-10/2020
        'EPI_ISL_19037047', # NONE - ** - hCoV-19/Ghana/CRI-11/2020
        'EPI_ISL_19037048', # NONE - ** - hCoV-19/Ghana/CRI-13/2020
        'EPI_ISL_19037049', # NONE - ** - hCoV-19/Ghana/CRI-14/2020
        'EPI_ISL_19037050', # NONE - ** - hCoV-19/Ghana/CRI-15/2020
        'EPI_ISL_19037051', # NONE - ** - hCoV-19/Ghana/CRI-16/2020
        
        # Collected by epi week 5
        'EPI_ISL_17121378', # [6] - ** - hCoV-19/USA/OK-PHL-0022663/2020 - USA/OK-PHL-0022663/2020 in NS-2023-06-18
        'EPI_ISL_2426018',  # [6] - hCoV-19/Norway/7651/2020 - Norway/7651/2020 in NS-2021-06-07
        'EPI_ISL_3364539',  # NONE - ** - hCoV-19/USA/NY-GEO-0231/2020

        # Collected by epi week 6
        'EPI_ISL_17846484', # NONE - ** - hCoV-19/USA/un-KDHE-2352575/2020
        'EPI_ISL_1603195',  # [6] - ** - hCoV-19/Spain/CL-COV00781/2020 - Spain/CL-COV00781/2020 in NS-2021-04-19
        'EPI_ISL_4899903',  # [6] - hCoV-19/Morocco/INH-108/2020 - Morocco/INH-108/2020 in NS-2021-10-13
        'EPI_ISL_4899911',  # [6] - hCoV-19/Morocco/INH-109/2020 - Morocco/INH-109/2020 in NS-2021-10-13
        'EPI_ISL_4899898',  # [6] - hCoV-19/Morocco/INH-107/2020 - Morocco/INH-107/2020 in NS-2021-10-13
        'EPI_ISL_4899870',  # [6] - hCoV-19/Morocco/INH-103/2020 - Morocco/INH-103/2020 in NS-2022-09-18
        'EPI_ISL_4899881',  # [6] - hCoV-19/Morocco/INH-104/2020 - Morocco/INH-104/2020 in NS-2021-10-19
        'EPI_ISL_4899863',  # [6] - hCoV-19/Morocco/INH-101/2020 - Morocco/INH-101/2020 in NS-2021-10-13
        'EPI_ISL_4899888',  # [6] - hCoV-19/Morocco/INH-105/2020 - Morocco/INH-105/2020 in NS-2021-10-14
        'EPI_ISL_4899917',  # [6] - hCoV-19/Morocco/INH-MN908947/2020 - Morocco/INH-MN908947/2020 in NS-2022-10-06
        'EPI_ISL_4899892',  # NONE - hCoV-19/Morocco/INH-106/2020
        'EPI_ISL_1263332',  # [6] - hCoV-19/Belgium/UZA-UA-CV0615326772/2020 - Belgium/UZA-UA-CV0615326772/2020 in NS-2021-03-18

        # Collected by epi week 7
        'EPI_ISL_1265909',  # [6] - ** - hCoV-19/USA/OH-ODH-SC1040172/2020 - USA/OH-ODH-SC1040172/2020 in NS-2021-03-19

        # Collected by epi week 8
        'EPI_ISL_1014733',  # [6] - hCoV-19/Spain/MD-IBV-99018532/2020 - Spain/MD-IBV-99018532/2020 in NS-2021-02-18
        'EPI_ISL_1311841',  # [6] - hCoV-19/Netherlands/ZH-EMC-2080/2020 - Netherlands/ZH-EMC-2080/2020 in NS-2021-03-24
        'EPI_ISL_1311840',  # [6] - hCoV-19/Netherlands/ZH-EMC-2079/2020 - Netherlands/ZH-EMC-2079/2020 in NS-2021-03-24
        'EPI_ISL_10980369', # NONE - hCoV-19/Zambia/MH23_107_4328/2020
        'EPI_ISL_1136974',  # [6] - hCoV-19/USA/NV-NSPHL-338695/2020 - USA/NV-NSPHL-338695/2020 in NS-2021-03-04
        'EPI_ISL_9879582',  # [6] - ** - hCoV-19/Mongolia/UB-109429/2020 - Mongolia/UB-109429/2020 in NS-2022-04-30
        'EPI_ISL_1167830',  # [6] - hCoV-19/Chile/AR-265171/2020 - Chile/AR-265171/2020 in NS-2021-03-08

        # The following appear by epi week 5 but they looked borderline.  By epi week 8, it's becoming clearer that
        # there's something wrong with them
        'EPI_ISL_7946128',  # NONE - hCoV-19/USA/CA-SEARCH-58338/2020
        'EPI_ISL_7946211',  # NONE - hCoV-19/USA/CA-SEARCH-58362/2020
        'EPI_ISL_7946167',  # NONE - hCoV-19/USA/CA-SEARCH-58348/2020

        # Collected by epi week 9
        'EPI_ISL_10980370', # NONE - ** - hCoV-19/Zambia/MH21_105_5506/2020
        'EPI_ISL_417446',   # [6] - ** - hCoV-19/Italy/LOM-UniMI02/2020 - Italy/UniMI02/2020 in NS-2020-06-08
        
        # Collected by epi week 10
        'EPI_ISL_2758215',  # [6] - ** - hCoV-19/India/un-IRSHA-CD210871/2020 - India/un-IRSHA-CD210871/2020 in NS-2021-06-05
        'EPI_ISL_2758214',  # [6] - ** - hCoV-19/India/un-IRSHA-CD210927/2020 - India/un-IRSHA-CD210927/2020 in NS-2021-06-05
        'EPI_ISL_2758213',  # [6] - ** - hCoV-19/India/un-IRSHA-CD210929/2020 - India/un-IRSHA-CD210929/2020 in NS-2021-06-05
]

for badAccessionID in badAccessionIDs:
    interestingAccessionIDs.discard(badAccessionID)

print(f'INFO: {len(interestingAccessionIDs)} sequences pass QC and have precise collection and submission dates:')
print(f'INFO:  - Host is Human')
print(f'INFO:  - Sequence length >= {min_seq_len}')
print(f'INFO:  - Collection date in range [{min_collection_date}, {max_collection_date}]')
print(f'INFO:  - Submission date in range [{min_submission_date}, {max_submission_date}]')

# Prepare input data folder (BAIL if it's there instead of overwriting data)
if args.mode == MODE_SUBMITTED_BY_DATE:
    inputs_path = Path('inputs_by_submission_date')
else:
    inputs_path = Path('inputs_by_collection_date')
inputs_path.mkdir(parents=True, exist_ok=False)  # Bail out if this exists
    
# Read in interesting sequences
def read_fasta(f):
    inSequence = False
    curFastaId = ''
    curSequence = ''
    for line in f:
        if type(line) != str:
            line = line.decode('utf-8')
        if not line.strip():
            continue
        if inSequence and line[0] == '>':
            yield (curFastaId, curSequence)
            inSequence = False
        if inSequence:
            curSequence += line.strip()
        else:
            if line[0] != '>':
                raise ValueError(f"Expecting FASTA ID line, found {line}")
            curFastaId = line[1:].strip()
            curSequence = ''
            inSequence = True
    if inSequence and curSequence:
        yield (curFastaId, curSequence)
    
sequences = {}  # Key: AccessionID, Value: Raw Unaligned Sequence
with lzma.open('20200331.fasta.xz', mode='rt', encoding='utf-8') as f:
    numRead = 0
    for fastaId, seq in read_fasta(f):
        virusName, _ = fastaId.split('|', maxsplit=1)
        accessionId = virusName2AccessionId[virusName]
        if accessionId in interestingAccessionIDs:
            
            # Discard very low quality sequences
            
            if len(seq) != metadata[accessionId]['length']:
                print(f"WARNING: Dropping {accessionId} = {metadata[accessionId]['virusName']}, actual length {len(seq)} doesn't match length in metadata {metadata[accessionId]['length']}")
                interestingAccessionIDs.discard(accessionId)
                continue

            seq = seq.lower()
            if (seq.count('n') + seq.count('-')) > (len(seq) // 10):
                print(f'WARNING: Dropping {accessionId} = {metadata[accessionId]["virusName"]}, more than 10% missing')
                interestingAccessionIDs.discard(accessionId)
                continue
            
            sequences[accessionId] = seq
            numRead += 1
            print(f'({numRead} / {len(interestingAccessionIDs)}) Read {accessionId} = {metadata[accessionId]["virusName"]}')



# Split out sequences into CDC epi weeks according to submission date
def epiWeek(dateStr):
    m = dateRE.match(dateStr)
    if not m:
        raise ValueError(f'Not a precise ISO-8601 date: {dateStr}')
    year = int(m.group(1))
    month = int(m.group(2))
    day = int(m.group(3))

    d = datetime.datetime(year, month, day)

    # Find end of epi week
    weekday = d.weekday()        # 0 = Monday, ..., 6 = Sunday
    weekday = (weekday + 1) % 7  # 0 = Sunday, ..., 6 = Saturday

    endOfEpiWeek = d + datetime.timedelta(days=6-weekday)  # CDC Epi Weeks end on a Saturday

    return endOfEpiWeek.year*100 + ((endOfEpiWeek - datetime.datetime(endOfEpiWeek.year, 1, 1)).days // 7 + 1)
    
interestingAccessionIDsByEpiWeek = defaultdict(set)  # Key: Epi week, e.g., 202043; Value: set of accession IDs
for accessionId in interestingAccessionIDs:
    info = metadata[accessionId]
    if args.mode == MODE_SUBMITTED_BY_DATE:
        week = epiWeek(info['submissionDate'])
    else:
        week = epiWeek(info['collectionDate'])
    interestingAccessionIDsByEpiWeek[week].add(accessionId)

for week in sorted(interestingAccessionIDsByEpiWeek.keys()):
    seq_type = 'submissions' if args.mode == MODE_SUBMITTED_BY_DATE else 'collections'
    print(f'INFO: In epi week {week}, there were {len(interestingAccessionIDsByEpiWeek[week])} {seq_type}')

# Produce unaligned FASTA and metadata files for all the weeks
for week in sorted(interestingAccessionIDsByEpiWeek.keys()):
    print(f'INFO: Preparing data for epi week {week}')
    week_path = inputs_path / f'to_epi_week_{week}'
    week_path.mkdir(parents=True, exist_ok=False)  # Bail out if exists

    unaligned_fasta_path = week_path / f'to_epi_week_{week}_unaligned.fasta'
    metadata_path = week_path / f'to_epi_week_{week}.tsv'
    with open(unaligned_fasta_path, 'w') as ff:
        with open(metadata_path, 'w') as fm:
            fm.write('id\t' + '\t'.join(f'locLevel{n+1}' for n in range(6)) + '\n')

            selectedIds = set()
            for accessionId in interestingAccessionIDs:
                if accessionId not in sequences:
                    print(f'WARNING: Skipping {accessionId}, no sequence found (not just dropped owing to QC)')
                else:
                    if ((args.mode == MODE_SUBMITTED_BY_DATE and epiWeek(metadata[accessionId]['submissionDate']) <= week)
                        or
                        (args.mode == MODE_COLLECTED_BY_DATE and epiWeek(metadata[accessionId]['collectionDate']) <= week)):
                        
                        selectedIds.add(accessionId)
                    

            for accessionId in selectedIds:
                ff.write(f">{accessionId}|{metadata[accessionId]['collectionDate']}")
                ff.write('\n')
                ff.write(sequences[accessionId])
                ff.write('\n')
            
                locByLevel = metadata[accessionId]["location"].split(' / ')
                fm.write(f'{accessionId}\t' + '\t'.join(' / '.join(locByLevel[:n+1]) for n in range(6)) + '\n')

    
    # Leverage NC_045512.2 reference genome in sars-cov-2-lemieux
    ref_fasta_path = Path('..') / 'sars-cov-2-lemieux' / 'scratch' / 'ref.fasta'
    aligned_fasta_path = week_path / f'to_epi_week_{week}_aligned_but_not_masked.fasta'
    print(f'INFO: Aligning sequences with mafft into {aligned_fasta_path.as_posix()}')
    
    subprocess.run([
        path_to_mafft.as_posix(),
        "--version"
    ])

    with open(aligned_fasta_path, 'w') as f:
        subprocess.run([
            path_to_mafft.as_posix(),
            "--thread", "-1",
            "--auto",
            "--addfragments", unaligned_fasta_path.as_posix(),
            "--keeplength",
            ref_fasta_path.as_posix(),
        ], stdout=f)
    
    # Mask tail ends naively (as in LeMieux et al 2021)
    aligned_and_masked_fasta_path = week_path / f'to_epi_week_{week}.fasta'
    print(f'INFO: Naive masking of tail ends into {aligned_and_masked_fasta_path.as_posix()}')
    
    with open(aligned_fasta_path, 'r') as r:
        with open(aligned_and_masked_fasta_path, 'w') as w:
            for fastaId, seq in read_fasta(r):
                # NOTE: NC_045512.2 isn't in GISAID (it differs in two sites from EPI_ISL_406798)
                # and I'm not sure of its collection date, so remove it from the dataset
                if fastaId.startswith('NC_045512.2'):
                    #fastaId = 'NC_045512.2|2019-12-26'
                    continue
                w.write(f'>{fastaId}\n')
                w.write('-' * num_initial_masked_sites)
                w.write(seq[num_initial_masked_sites:-num_final_masked_sites])
                w.write('-' * num_final_masked_sites)
                w.write('\n')
    
