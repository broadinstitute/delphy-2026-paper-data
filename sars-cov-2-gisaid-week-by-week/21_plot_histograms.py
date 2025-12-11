#!/usr/bin/env python

import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import math
from pathlib import Path
import baltic as bt

# Conversions for the humans
# ==========================
epiweek_2_hyphenated = {
    '202001': '2020-01',
    '202002': '2020-02',
    '202003': '2020-03',
    '202004': '2020-04',
    '202005': '2020-05',
    '202006': '2020-06',
    '202007': '2020-07',
    '202008': '2020-08',
    '202009': '2020-09',
    '202010': '2020-10',
    '202011': '2020-11',
    '202012': '2020-12',
    '202013': '2020-13',
}
def extract_log_params(log_filename, burnin):
    # Extremely brittle way of extracting distributions from the log files
    raw_data = pd.read_table(log_filename, comment='#')
    num_pts = len(raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    data = raw_data[burnin_pts:]
    print(f'Read in log file {log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(data)} usable')

    result = {}
    for key in [
            'TreeHeight',
            'clockRate',
            'growthRate',
    ]:
        result[key] = list(data[key])

    return result

def mean_and_95HPD(xs):
    mean = sum(xs)/len(xs)
    
    sorted_xs = sorted(xs)
    hpd_size = int(len(xs)*0.95)
    hpd = (sorted_xs[0], sorted_xs[-1])
    best_hpd_width = sorted_xs[-1] - sorted_xs[0]
    for i in range(len(sorted_xs)-hpd_size):
        hpd_i_width = sorted_xs[i+hpd_size] - sorted_xs[i]
        if hpd_i_width < best_hpd_width:
            best_hpd_width = hpd_i_width
            hpd = (sorted_xs[i], sorted_xs[i+hpd_size])

    return mean, hpd

def do_plots(epiweek_2_log_file, epiweek_2_tree_file, out_pdf_filename):
    fig,ax = plt.subplots(3, 1, figsize=(11,8.5), facecolor='w')
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.1, hspace=0.3)

    epiweek_2_logs = {}
    for epiweek, log_filename in epiweek_2_log_file.items():
        epiweek_2_logs[epiweek] = extract_log_params(log_filename, burnin=0.30)

    epiweek_2_trees = {}
    epiweek_2_maxdate = {}
    for epiweek, tree_filename in epiweek_2_tree_file.items():
        print(f'Reading {tree_filename}')
        epiweek_2_trees[epiweek] = mcc_tree = bt.loadNexus(tree_filename)
        leaves = mcc_tree.traverse_tree()
        max_date = max(leaf.name.split('|')[-1] for leaf in leaves)  # e.g., 2020-03-01
        epiweek_2_maxdate[epiweek] = max_date

    kwargs = dict(showmeans=True, showextrema=False, quantiles=[[0.025,0.975] for log in epiweek_2_logs.values()])
    
    ax = plt.subplot(3,1,1)
    for y in [0.5, 1.0, 1.5, 2.0]:
        ax.axhline(y, color='#DDDDDD', zorder=-1, lw=0.5)
    ax.violinplot([[mu*1e3 for mu in log['clockRate']] for log in epiweek_2_logs.values()], **kwargs)
    ax.set_ylim([0.0, 2.5])
    ax.set_xticks(np.arange(1, len(epiweek_2_logs)+1))
    ax.set_xticklabels(list(epiweek_2_hyphenated[w] for w in epiweek_2_logs.keys()))
    ax.set_ylabel('Clock rate (x10${}^{-3}$/site/year)')
    
    ax = plt.subplot(3,1,2)
    for y in [0.5, 1.0, 1.5]:
        ax.axhline(y, color='#DDDDDD', zorder=-1, lw=0.5)
    ax.violinplot([[12*np.log(2)/g for g in log['growthRate']] for log in epiweek_2_logs.values()], **kwargs)
    ax.set_ylim([0.0, 2.0])
    ax.set_xticks(np.arange(1, len(epiweek_2_logs)+1))
    ax.set_xticklabels(list(epiweek_2_hyphenated[w] for w in epiweek_2_logs.keys()))
    ax.set_ylabel('Doubling time (months)')

    for w, log in epiweek_2_logs.items():
        mean_dbl, hpd_dbl = mean_and_95HPD([365*np.log(2)/g for g in log['growthRate']])
        print(f'Epiweek {epiweek_2_hyphenated[w]}: {mean_dbl} days (95% HPD {hpd_dbl[0]}-{hpd_dbl[1]})')
    
    ax = plt.subplot(3,1,3)
    ax.violinplot([[bt.decimalDate(epiweek_2_maxdate[epiweek]) - h for h in log['TreeHeight']]
                   for epiweek, log in epiweek_2_logs.items()],
                  **kwargs)
    ax.set_ylim([bt.decimalDate('2019-09-01'), bt.decimalDate('2020-01-01')])
    yticks = {
        '1 Sep': '2019-09-01',
        '1 Oct': '2019-10-01',
        '1 Nov': '2019-11-01',
        '1 Dec': '2019-12-01',
        '1 Jan': '2020-01-01',
       }
    ax.set_yticks([bt.decimalDate(datestr) for label, datestr in yticks.items()])
    ax.set_yticklabels(list(yticks.keys()))
    for label, datestr in yticks.items():
        ax.axhline(bt.decimalDate(datestr), color='#DDDDDD', zorder=-1, lw=0.5)
    
    ax.set_xticks(np.arange(1, len(epiweek_2_logs)+1))
    ax.set_xticklabels(list(epiweek_2_hyphenated[w] for w in epiweek_2_logs.keys()))
    ax.set_ylabel('tMRCA (2019)')
    ax.set_xlabel('CDC Epiweek')
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')




Path('plots').mkdir(parents=True, exist_ok=True)

for run in ('a',):#, 'b'):
    epiweek_2_log_file = {
        epiweek: f'./outputs_by_submission_date_{run}/to_epi_week_{epiweek}/to_epi_week_{epiweek}.log'
        for epiweek in [
                '202002',
                '202003',
                '202004',
                '202005',
                '202006',
                '202007',
                '202008',
                '202009',
                '202010',
                '202011',
                '202012',
                '202013',
        ]
       }
    epiweek_2_tree_file = {
        epiweek: f'./outputs_by_submission_date_{run}/to_epi_week_{epiweek}/to_epi_week_{epiweek}.trees'
        for epiweek in epiweek_2_log_file.keys()
       }
    do_plots(epiweek_2_log_file, epiweek_2_tree_file, f'plots/BySubmissionDate{run.upper()}-CompoundDistrs.pdf')
    
    # epiweek_2_log_file = {
    #     epiweek: f'./outputs_by_collection_date_{run}/to_epi_week_{epiweek}/to_epi_week_{epiweek}.log'
    #     for epiweek in [
    #             '202002',
    #             '202003',
    #             '202004',
    #             '202005',
    #             '202006',
    #             '202007',
    #             '202008',
    #             '202009',
    #             '202010',
    #     ]
    #    }
    # epiweek_2_tree_file = {
    #     epiweek: f'./outputs_by_collection_date_{run}/to_epi_week_{epiweek}/to_epi_week_{epiweek}.trees'
    #     for epiweek in epiweek_2_log_file.keys()
    #    }
    # do_plots(epiweek_2_log_file, epiweek_2_tree_file, f'plots/ByCollectionDate{run.upper()}-CompoundDistrs.pdf')
