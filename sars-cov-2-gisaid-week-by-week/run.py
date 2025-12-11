#!/usr/bin/env python3

import argparse
from pathlib import Path
import datetime
import subprocess

MODE_SUBMITTED_BY_DATE = 'submissionDate'
MODE_COLLECTED_BY_DATE = 'collectionDate'

parser = argparse.ArgumentParser(description="Driver for Delphy runs covering first N weeks of SARS-CoV-2")
parser.add_argument('--mode', help='Include samples in run if _this_ date is on or before CDC epi week N',
                    required=True, choices=[MODE_SUBMITTED_BY_DATE, MODE_COLLECTED_BY_DATE])
parser.add_argument('--batch', help='Batch ID (a or b)',
                    required=True, choices=['a', 'b'])

args = parser.parse_args()

path_to_delphy = Path("../delphy")

# Prepare folders (BAIL if it's there instead of overwriting data)
batch = args.batch
if args.mode == MODE_SUBMITTED_BY_DATE:
    inputs_path = Path('inputs_by_submission_date')
    outputs_path = Path(f'outputs_by_submission_date_{batch}')
else:
    inputs_path = Path('inputs_by_collection_date')
    outputs_path = Path(f'outputs_by_collection_date_{batch}')

if not inputs_path.exists():
    raise SystemExit(f'Inputs folder {inputs_path.as_posix()} not found')
outputs_path.mkdir(parents=True, exist_ok=False)  # Bail out if this exists

# Go over each epi week
with open(outputs_path / 'all_timing.tsv', mode='wt', encoding='utf-8') as all_timing_f:
    all_timing_f.write(
        '\t'.join([
            'Epi week',
            'Number of sequences',
            'Start time',
            'End time',
            'Total time (s)',
            'Total steps',
            'Number of threads',
            'Steps / s'
        ]) + '\n')
    
    epiweek_dirs = [d for d in inputs_path.iterdir()]
    epiweek_dirs.sort(key=lambda d : d.name)
    for epiweek_inputs_path in epiweek_dirs:
        epiweek_str = epiweek_inputs_path.name[-len('202001'):]
        epiweek = int(epiweek_str)
        if epiweek < 202001:
            # Epi weeks before 2020 are too sparse
            continue
        
        start_time = datetime.datetime.now()
        print(f'- {start_time.isoformat()}: Running epi week {epiweek}')
    
        epiweek_outputs_path = outputs_path / f'to_epi_week_{epiweek}'
        epiweek_outputs_path.mkdir(parents=True, exist_ok=False)  # Bail out if this exists
        
        # Count sequences
        input_fasta_path = epiweek_inputs_path / f'to_epi_week_{epiweek}.fasta'
        num_seqs = 0
        with open(input_fasta_path, mode='rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('>'):
                    num_seqs += 1
        print(f'  Found {num_seqs} sequences')
    
        # Decide on number of steps, sampling rate and number of threads
        # Rough heuristics:
        #  - 5,000,000 steps per sequence
        #  - 200 samples in .dphy file & trees file
        #  - 10,000 samples in log file
        #  - At least 100 samples per thread
        #  - At least 1 thread, no more than 32 threads
    
        steps_per_seq = 5_000_000
        num_steps = num_seqs * steps_per_seq
        steps_per_sample = num_steps // 200
        steps_per_log = num_steps // 10_000
        num_threads = max(1, min(32, num_seqs // 100))
    
        # Build delphy command line
        output_log_path = epiweek_outputs_path / f'to_epi_week_{epiweek}.log'
        output_trees_path = epiweek_outputs_path / f'to_epi_week_{epiweek}.trees'
        output_dphy_path = epiweek_outputs_path / f'to_epi_week_{epiweek}.dphy'
        delphy_cli = [
            path_to_delphy.as_posix(),
            "--v0-in-fasta", input_fasta_path.as_posix(),
            "--v0-threads", str(num_threads),
            "--v0-steps", str(num_steps),
            "--v0-out-log-file", output_log_path.as_posix(),
            "--v0-log-every", str(steps_per_log),
            "--v0-out-trees-file", output_trees_path.as_posix(),
            "--v0-tree-every", str(steps_per_sample),
            "--v0-out-delphy-file", output_dphy_path.as_posix(),
            "--v0-delphy-snapshot-every", str(steps_per_sample),
        ]
    
        # Record command
        cmd_path = epiweek_outputs_path / 'run.sh'
        with open(cmd_path, mode='wt', encoding='utf-8') as f:
            f.write(" ".join(delphy_cli) + '\n')  # Ignore quoting subtleties
    
        # Run Delphy
        subprocess.run(delphy_cli)
    
        # Finish timing info
        end_time = datetime.datetime.now()
    
        time_for_run = end_time - start_time
        steps_per_second = num_steps / time_for_run.total_seconds()
    
        # Record timing
        timing_path = epiweek_outputs_path / 'timing.tsv'
        timing_str = (
            '\t'.join([
                f'{epiweek}',
                f'{num_seqs}',
                f'{start_time.isoformat()}',
                f'{end_time.isoformat()}',
                f'{time_for_run.total_seconds()}',
                f'{num_steps}',
                f'{num_threads}',
                f'{steps_per_second}'
            ]) + '\n'
        )
        with open(timing_path, mode='wt', encoding='utf-8') as f:
            f.write(timing_str)
        all_timing_f.write(timing_str)
    
        print(f'  {end_time.isoformat()}: Finished {epiweek}, took {time_for_run.total_seconds()} s = {steps_per_second} steps / s')
