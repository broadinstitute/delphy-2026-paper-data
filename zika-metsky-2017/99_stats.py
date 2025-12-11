#!/usr/bin/env python3

import sys

def load_fasta(f):
    result = {}
    current_id = ''
    current_seq = ''
    for line in f:
        if not line.strip():
            continue
        if line.startswith('>'):
            if current_seq:
                result[current_id] = current_seq
            current_id = line[1:].strip()
            current_seq = ''
        else:
            current_seq += line.strip()
    if current_seq:
        result[current_id] = current_seq
    return result

with open('delphy_inputs/zika.fasta', 'r') as f:
    samples = load_fasta(f)

print(f'Number of sequences: {len(samples)}')
for seq in samples.values():
    L = len(seq)
    break
print(f'Sequence length: {L}')

# Date range
dates = {}
for sample_id, seq in samples.items():
    dates[sample_id] = sample_id[-len('2020-01-01'):]
print(f'Date range: {min(dates.values())} to {max(dates.values())}')

# Missing data
num_missings = []
for sample_id, seq in samples.items():
    num_missing = sum((c == '-' or c == 'n') for c in seq.lower())
    if num_missing == 5285:
        print(sample_id)
    num_missings.append(num_missing)

num_missings.sort()

print(f'Missing data: between {min(num_missings)} bases = {100*min(num_missings)/L:.3} % '
      + f'and {max(num_missings)} bases = {100*max(num_missings)/L:.3} % per sequence')

avg_missing = sum(num_missings) / len(num_missings)
print(f'Avg missing: {avg_missing} bases = {100*avg_missing / L:.3} %')
