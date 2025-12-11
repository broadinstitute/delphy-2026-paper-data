#!/usr/bin/env python3

import argparse
from pathlib import Path
import datetime
import subprocess
import json

parser = argparse.ArgumentParser(description="Driver for Delphy runs for simulated data")
parser.add_argument('--sim', help='Name of simulation (e.g., const_10000)', required=True)
parser.add_argument('--rep', help='Name of replica (e.g., a or b)', required=True)
parser.add_argument('--coal-cells', help='Target number of coal prior cells (default = 400)', required=False, type=int)

args = parser.parse_args()

path_to_delphy = Path("../delphy")

sim_path = Path(args.sim)
ground_truth_path = sim_path / 'ground_truth'
inputs_path = sim_path / 'inputs'
delphy_outputs_path = sim_path / f'delphy_outputs_{args.rep}'
delphy_outputs_path.mkdir(parents=True, exist_ok=True)

if not sim_path.exists():
    raise SystemExit(f'Simulation folder {sim_path.as_posix()} not found')
input_maple_path = inputs_path / f'{args.sim}.maple'

start_time = datetime.datetime.now()
print(f'- {start_time.isoformat()}: Running simulation {args.sim} rep {args.rep}')

with open(ground_truth_path / f'{args.sim}_info.json', 'rt', encoding='utf-8') as f:
    info = json.load(f)

num_seqs = info['sampling_strategy']['num_samples']
print(f'  Simulation has {num_seqs} sequences')
    
# Decide on number of steps, sampling rate and number of threads
# Rough heuristics:
#  - 5,000,000 steps per sequence
#  - 200 samples in .dphy file & trees file
#  - 10,000 samples in log file
#  - At least 20 sequences per thread (~40 nodes per partition)
#  - At least 1 thread, no more than 2*96 threads

steps_per_seq = 5_000_000
num_steps = num_seqs * steps_per_seq
steps_per_sample = num_steps // 200
steps_per_log = num_steps // 10_000
num_threads = max(1, min(2*96, num_seqs // 20))

# Build delphy command line
output_log_path = delphy_outputs_path / f'{args.sim}.log'
output_trees_path = delphy_outputs_path / f'{args.sim}.trees'
output_dphy_path = delphy_outputs_path / f'{args.sim}.dphy'
delphy_cli = [
    path_to_delphy.as_posix(),
    "--v0-in-maple", input_maple_path.as_posix(),
    "--v0-threads", str(num_threads),
    "--v0-steps", str(num_steps),
    "--v0-out-log-file", output_log_path.as_posix(),
    "--v0-log-every", str(steps_per_log),
    "--v0-out-trees-file", output_trees_path.as_posix(),
    "--v0-tree-every", str(steps_per_sample),
    "--v0-out-delphy-file", output_dphy_path.as_posix(),
    "--v0-delphy-snapshot-every", str(steps_per_sample),
]
if args.coal_cells:
    delphy_cli.extend([
        "--v0-target-coal-prior-cells", str(args.coal_cells),
    ])

# Record command
cmd_path = delphy_outputs_path / 'run.sh'
with open(cmd_path, mode='wt', encoding='utf-8') as f:
    f.write(" ".join(delphy_cli) + '\n')  # Ignore quoting subtleties
    
# Run Delphy
subprocess.run(delphy_cli)

# Finish timing info
end_time = datetime.datetime.now()

time_for_run = end_time - start_time
steps_per_second = num_steps / time_for_run.total_seconds()
    
# Record timing
timing_path = delphy_outputs_path / 'timing.tsv'
timing_str = (
    '\t'.join([
        f'{args.sim}',
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
    
print(f'  {end_time.isoformat()}: Finished {args.sim} rep {args.rep}, took {time_for_run.total_seconds()} s = {steps_per_second} steps / s')
