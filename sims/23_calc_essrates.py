#!/usr/bin/env python

import pandas as pd
import subprocess
from io import StringIO
import math
from pathlib import Path
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

def analyse_log(log_filename, burnin=0.30):
    result = subprocess.run(
        [
            '../loganalyser2',
            '-burnin', str(math.floor(burnin*100)),
            log_filename,
        ],
        stdout=subprocess.PIPE
    )
    raw_analysis = result.stdout.decode('utf-8')
    return pd.read_fwf(StringIO(raw_analysis), index_col=0)

def analyse_sim(sim_name, burnin=0.30, log_filename=None, timing_filename=None):
    if not log_filename:
        log_filename = f'./{sim_name}/delphy_outputs_a/{sim_name}.log'
    if not timing_filename:
        timing_filename = f'./{sim_name}/delphy_outputs_a/timing.tsv'

    analysis = analyse_log(log_filename, burnin)
    with open(timing_filename, 'r') as f:
        for line in f:
            _, num_seqs, start_time, end_time, wallclock_seconds, num_steps, num_threads, steps_per_second = line.strip().split('\t')
            num_seqs = int(num_seqs)
            wallclock_seconds = float(wallclock_seconds)
            num_steps = int(num_steps)
            num_threads = int(num_threads)
            break

    ess_posterior = analysis['ESS']['posterior_for_Delphy']

    print(f'For {sim_name}, num_steps={num_steps}, num_threads={num_threads}, ESS={ess_posterior}, wallclock_seconds={wallclock_seconds}, wallclock_hours={wallclock_seconds/60.0/60.0}, ESS / hour = {ess_posterior / (wallclock_seconds/60.0/60.0)}')

analyse_sim('exp_100')
analyse_sim('exp_1000')
analyse_sim('exp_10000')
analyse_sim('exp_100000')

analyse_sim('const_100')
analyse_sim('const_1000')
analyse_sim('const_10000')
analyse_sim('const_100000')
