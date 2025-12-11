#!/usr/bin/env python

import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import math
from pathlib import Path
import datetime

def plot_logs(delphy_log_filename_and_colors_and_labels, out_pdf_filename, ground_truth_nwk, burnin=0.10):
    delphy_raw_datas = []
    delphy_datas = []
    for (delphy_log_filename, _, _) in delphy_log_filename_and_colors_and_labels:
        delphy_raw_datas.append(pd.read_table(delphy_log_filename, comment='#'))
        num_pts = len(delphy_raw_datas[-1])
        burnin_pts = math.floor(burnin * num_pts)
        delphy_datas.append(delphy_raw_datas[-1][burnin_pts:])
        print(f'Read in Delphy log in {delphy_log_filename}')
        print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(delphy_datas[-1])} usable')

    # The ground truth trees end something like this: ...683)NODE_199|2024-01-24:0;
    # The following is an extremely brittle way to extract the date of the root node
    # (really, I should have included this in sapling's info.json...)
    with open(ground_truth_nwk, 'r') as f:
        nwk_str = f.read()
    root_date_start = nwk_str.rfind('|')+1
    root_date_end = nwk_str.rfind(':')
    root_date = nwk_str[root_date_start:root_date_end]
    print(f'Read in root date of {root_date} from {ground_truth_nwk}')

    fig,ax = plt.subplots(3, 4, figsize=(11,8.5), facecolor='w')
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, hspace=0.3)

    def make_subplot(i, colname, title, actual_value=None, adjuster=lambda x: x):
        plt.subplot(3,4,i)
        plt.title(title)
        plt.xlabel(None)

        for (delphy_data, (_, data_color, _)) in zip(delphy_datas, delphy_log_filename_and_colors_and_labels):
            this_ax = sb.kdeplot(adjuster(delphy_data[colname]), fill=True, color=data_color);
        this_ax.set(xlabel=None)
        this_ax.get_yaxis().set_visible(False)
        if actual_value is not None:
            this_ax.axvline(actual_value, color='green')
            
        return this_ax

    make_subplot( 1, 'clockRate', r'Mutation rate (x$10^{-3}$/site/year)', 1.0, lambda mu: mu*1000)
    make_subplot( 2, 'TreeHeight', r'Tree Height (years)',
                  (datetime.date(2024, 7, 31) - datetime.date.fromisoformat(root_date)).days / 365.0)
    this_ax = make_subplot( 3, 'kappa', r'HKY $\kappa$ parameter', 5.0)
    #fig.legend(this_ax.get_legend_handles_labels()[0], ['BEAST2', 'Delphy'], loc='lower right')
    fig.delaxes(ax[0][3])
    
    make_subplot( 5, 'freqParameter.1', r'HKY $\pi_A$', 0.30)
    make_subplot( 6, 'freqParameter.2', r'HKY $\pi_C$', 0.18)
    make_subplot( 7, 'freqParameter.3', r'HKY $\pi_G$', 0.20)
    make_subplot( 8, 'freqParameter.4', r'HKY $\pi_T$', 0.32)

    this_ax = make_subplot( 9, 'CoalescentExponential', r'Coalescent Prior')
    this_ax.set_xticks([-170000, -140000, -110000])
    make_subplot(10, 'ePopSize', r'Final Pop Size (years)',
                 6.0 if delphy_log_filename.startswith('exp') else 2.0)
    make_subplot(11, 'growthRate', r'Pop Growth Rate (1/year)',
                 10.0 if delphy_log_filename.startswith('exp') else 0.0)
    fig.delaxes(ax[2][3])

    fig.legend([label for (_,_,label) in delphy_log_filename_and_colors_and_labels] + ['Ground Truth'],
               loc='lower right', 
               bbox_to_anchor=(0.95, 0.05))
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')

delphy_log_filename_and_colors_and_labels = []
colors = {
    625: '#FF0000',
    1250: '#BB4400',
    2500: '#888800',
    5000: '#44BB00',
    10000: '#00FF00',
}
for coal_bins in [625, 1250, 2500, 5000, 10000]:
    sim = 'exp_100000'
    rep = f'c_{coal_bins}bins' if coal_bins != 10000 else 'a'
    delphy_log_filename_and_colors_and_labels.append(
        (f'{sim}/delphy_outputs_{rep}/{sim}.log', colors[coal_bins], f'{coal_bins} bins'))
    
plot_logs(
    delphy_log_filename_and_colors_and_labels,
    f'{sim}/{sim}_bybin_distrs.pdf',
    f'{sim}/ground_truth/{sim}_tree.nwk',
    burnin=0.30
)
