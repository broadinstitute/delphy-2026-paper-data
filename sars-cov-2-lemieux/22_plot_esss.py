#!/usr/bin/env python

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from utils import *

import subprocess
from io import StringIO
import math
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

def compare_logs(log_and_wallclock_minutes: list[tuple[Log, float]],
                 out_pdf_filename,
                 show_alpha=False):

    ess_tables = []
    for log, wallclock_minutes in log_and_wallclock_minutes:
        ess_tables.append(read_in_ess_table(log.log_analysis_filename, log.tree_ess_filename, log.log_format))

    esss = []
    for (ess_table, (log, wallclock_minutes)) in zip(ess_tables, log_and_wallclock_minutes):
        this_esss = {
            'Posterior': ess_table['joint'],
            'Likelihood': ess_table['likelihood'],
            'Prior': ess_table['prior'],
            'Tree Height': ess_table['TreeHeight'],
            'Mutation Rate': ess_table['clockRate'],
            r'HKY $\pi_A$': ess_table['freqParameter.1'],
            r'HKY $\pi_C$': ess_table['freqParameter.2'],
            r'HKY $\pi_G$': ess_table['freqParameter.3'],
            r'HKY $\pi_T$': ess_table['freqParameter.4'],
            r'HKY $\kappa$': ess_table['kappa'],
            'Coalescent\nExponential': ess_table['CoalescentExponential'],
            'Final Pop Size': ess_table['ePopSize'],
            'Pop Growth Rate': ess_table['growthRate'],
            'Tree Topology': ess_table['TreeESS'],
        }
        if show_alpha:
            this_esss[r'Site Rate $\alpha$'] = ess_table['gammaShape']
        
        esss.append(this_esss)

    sorted_names = sorted(list(esss[0].keys()), key=lambda x: esss[0][x])

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)
    
    df = pd.DataFrame(
        [
            [f'{i+1}. {n}'] +
            [this_esss[n] / wallclock_minutes for this_esss, (log, wallclock_minutes) in zip(esss, log_and_wallclock_minutes)]
            for (i, n) in enumerate(sorted_names)
        ],
        columns=['Observable'] + [log.label for log, _ in log_and_wallclock_minutes])
    
    # view data 
    df.plot(x='Observable', 
            kind='bar', 
            stacked=False, 
            ax=ax,
            width=0.7,
            color=[log.color for log, _ in log_and_wallclock_minutes])

    ax.set_xlabel(None)
    ax.semilogy();
    ax.set_ylabel('Effective Samples per minute');
    #ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter(2));
    plt.xticks(rotation=40, ha='right');
    plt.legend(loc='upper right')
    ax.set_ylim([1e-3, 1e7])
    
    left, bottom, width, height = [0.19, 0.68, 0.35, 0.25]
    ax2 = fig.add_axes([left, bottom, width, height])

    ref_esss, (ref_log, ref_wallclock_minutes) = esss[0], log_and_wallclock_minutes[0]
    df2 = pd.DataFrame(
        [
            [f'{i+1}'] +
            [(ref_esss[n]/ref_wallclock_minutes)/(this_esss[n]/wallclock_minutes)
             for this_esss, (log, wallclock_minutes) in zip(esss, log_and_wallclock_minutes)][1:]
            for (i, n) in enumerate(sorted_names)
        ],
        columns=['Observable'] + ['Speedup ' + log.label for log, _ in log_and_wallclock_minutes[1:]])
    df2.plot(x='Observable', 
             kind='bar', 
             stacked=False, 
             ax=ax2,
             width=0.7,
             legend=False,
             color=[log.color for log, _ in log_and_wallclock_minutes[1:]])
    
    ax2.set_xlabel(None)
    ax2.set_ylabel('Speedup')
    ax2.set_ylim([10, 20000]);
    ax2.semilogy();
    plt.xticks(rotation=90);
    
    plt.savefig(out_pdf_filename, bbox_inches='tight')
    print(f'Plot saved to {out_pdf_filename}')


Path('plots').mkdir(parents=True, exist_ok=True)

# Non-alpha datasets
delphy_a_log = Log(
    filename='delphy_outputs_a/ma_sars_cov_2_delphy.log',
    log_analysis_filename='delphy_outputs_a/log_analysis.txt',
    tree_ess_filename='delphy_outputs_a/tree_ess.json',
    log_format=delphy_log_format,
    label='Delphy',
    color='green'
)
delphy_b_log = Log(
    filename='delphy_outputs_b/ma_sars_cov_2_delphy.log',
    log_analysis_filename='delphy_outputs_b/log_analysis.txt',
    tree_ess_filename='delphy_outputs_b/tree_ess.json',
    log_format=delphy_log_format,
    label='Delphy',
    color='green'
)
beast2_log = Log(
    filename='beast2_run/output.log',
    log_analysis_filename='beast2_run/log_analysis.txt',
    tree_ess_filename='beast2_run/tree_ess.json',
    log_format=beast2_log_format,
    label='BEAST 2',
    color='orange'
)
beastX_log = Log(
    filename='beastX_run/output.log',
    log_analysis_filename='beastX_run/log_analysis.txt',
    tree_ess_filename='beastX_run/tree_ess.json',
    log_format=beastX_log_format,
    label='BEAST X',
    color='blue'
)

compare_logs(
    [(delphy_a_log,   10 + 57.318/60.0),
     (beast2_log,   1323 + 35.510/60.0),
     (beastX_log,    760 + 27.714/60.0)],
    'plots/ESSDelphyAVsBeast2VsBeastX.pdf'
)
compare_logs(
    [(delphy_b_log,   10 + 56.218/60.0),
     (beast2_log,   1323 + 35.510/60.0),
     (beastX_log,    760 + 27.714/60.0)],
    'plots/ESSDelphyBVsBeast2VsBeastX.pdf'
)


# Alpha datasets
delphy_alpha_a_log = Log(
    filename='delphy_outputs_alpha_a/ma_sars_cov_2_delphy_alpha.log',
    log_analysis_filename='delphy_outputs_alpha_a/log_analysis.txt',
    tree_ess_filename='delphy_outputs_alpha_a/tree_ess.json',
    log_format=delphy_log_format,
    label=r'Delphy (K=$\infty$)',
    color='green'
)
delphy_alpha_b_log = Log(
    filename='delphy_outputs_alpha_b/ma_sars_cov_2_delphy_alpha.log',
    log_analysis_filename='delphy_outputs_alpha_b/log_analysis.txt',
    tree_ess_filename='delphy_outputs_alpha_b/tree_ess.json',
    log_format=delphy_log_format,
    label=r'Delphy (K=$\infty$)',
    color='green'
)
beast2_alpha_log = Log(
    filename='beast2_run_alpha/output.log',
    log_analysis_filename='beast2_run_alpha/log_analysis.txt',
    tree_ess_filename='beast2_run_alpha/tree_ess.json',
    log_format=beast2_log_format,
    label='BEAST 2 (K=4)',
    color='orange'
)
beastX_alpha_log = Log(
    filename='beastX_run_alpha/output.log',
    log_analysis_filename='beastX_run_alpha/log_analysis.txt',
    tree_ess_filename='beastX_run_alpha/tree_ess.json',
    log_format=beastX_log_format,
    label='BEAST X (K=4)',
    color='blue'
)
beastX_alpha_2cats_log = Log(
    filename='beastX_run_alpha_2cats/output.log',
    log_analysis_filename='beastX_run_alpha_2cats/log_analysis.txt',
    tree_ess_filename='beastX_run_alpha_2cats/tree_ess.json',
    log_format=beastX_log_format,
    label='BEAST X (K=2)',
    color='blue'
)
beastX_alpha_8cats_log = Log(
    filename='beastX_run_alpha_8cats/output.log',
    log_analysis_filename='beastX_run_alpha_8cats/log_analysis.txt',
    tree_ess_filename='beastX_run_alpha_8cats/tree_ess.json',
    log_format=beastX_log_format,
    label='BEAST X (K=8)',
    color='blue'
)

compare_logs(
    [(delphy_alpha_a_log,         13 + 26.296/60.0),
     (beast2_alpha_log,         3140 + 21.372/60.0),
     (beastX_alpha_log,         1563 + 57.460/60.0),
     #(beastX_alpha_2cats_log,   1043 +  5.921/60.0),
     #(beastX_alpha_8cats_log,   2764 +  6.714/60.0),
     ],
    'plots/ESSDelphyAVsBeast2VsBeastXAlpha.pdf'
)
compare_logs(
    [(delphy_alpha_b_log,         13 + 31.223/60.0),
     (beast2_alpha_log,         3140 + 21.372/60.0),
     (beastX_alpha_log,         1563 + 57.460/60.0),
     #(beastX_alpha_2cats_log,   1043 +  5.921/60.0),
     #(beastX_alpha_8cats_log,   2764 +  6.714/60.0),
     ],
    'plots/ESSDelphyBVsBeast2VsBeastXAlpha.pdf'
)
