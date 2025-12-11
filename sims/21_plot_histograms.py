#!/usr/bin/env python

import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import math
from pathlib import Path
import datetime

def extract_ml_results(iqtree_log_filename, tt_folder_name, look_for_alpha=False):
    """Extremely brittle IQ-Tree log file parser, enough for these plots though"""

    results = {}

    if not Path(iqtree_log_filename).exists():
        return {}  # Nothing!
    
    with open(iqtree_log_filename, 'r') as f:
        lines = f.readlines()

    # Extract HKY kappa
    # -----------------
    # Looking for a section like this:
    #
    # Rate parameter R:
    # 
    #   A-C: 1.0000
    #   A-G: 9.0197
    #   A-T: 1.0000
    #   C-G: 1.0000
    #   C-T: 9.0197
    #   G-T: 1.0000
    #
    rate_line_idx = [i
                     for i in range(len(lines))
                     if lines[i] == 'Rate parameter R:\n'][0]
    results['kappa'] = float(lines[rate_line_idx+3].strip().split(':')[1])

    # Extract HKY stationary state frequencies
    # ----------------------------------------
    # Looking for a section like this:
    #
    # State frequencies: (estimated with maximum likelihood)
    # 
    #   pi(A) = 0.3192
    #   pi(C) = 0.214
    #   pi(G) = 0.198
    #   pi(T) = 0.2688
    #
    state_freqs_line_idx = [i
                            for i in range(len(lines))
                            if lines[i].startswith('State frequencies')][0]
    results['freqParameter.1'] = float(lines[state_freqs_line_idx+2].split('=')[1])
    results['freqParameter.2'] = float(lines[state_freqs_line_idx+3].split('=')[1])
    results['freqParameter.3'] = float(lines[state_freqs_line_idx+4].split('=')[1])
    results['freqParameter.4'] = float(lines[state_freqs_line_idx+5].split('=')[1])

    # Extract site-rate heterogeneity parameter alpha
    # -----------------------------------------------
    # Looking for a section like this:
    #
    # Model of rate heterogeneity: Gamma with 4 categories
    # Gamma shape alpha: 998.4
    #
    if look_for_alpha:
        alpha_line_idx = [i
                          for i in range(len(lines))
                          if lines[i].startswith('Gamma shape alpha')][0]
        results['gammaShape'] = float(lines[alpha_line_idx].split(':')[1])

    # Extract mutation rate
    # ---------------------
    # Looking in tt/molecular_clock.txt for a line like this:
    #
    # --rate:	1.785e-03
    with (Path(tt_folder_name) / 'molecular_clock.txt').open('r') as f:
        for line in f:
            if 'rate' in line:
                results['clockRate'] = float(line.split(':')[1])

    # Extract tree height
    # -------------------
    # Read all node dates from tt/dates.tsv and extract rate
    with (Path(tt_folder_name) / 'dates.tsv').open('r') as f:
        dates = []
        for line in f:
            if line.startswith('#'):
                continue
            dates.append(float(line.split()[-1]))

        results['TreeHeight'] = max(dates) - min(dates)

    # Extract coalescent prior parameters
    # -----------------------------------
    # Assume tt/skyline.tsv has two useful lines, like this:
    #
    # #Skyline assuming 50.0 gen/year and approximate confidence bounds (+/- 2.000000 standard deviations of the LH)
    # #date 	N_e 	lower 	upper
    # 2014.175	6.230e+00	2.804e+00	1.384e+01
    # 2014.462	1.254e+01	9.726e+00	1.617e+01
    with (Path(tt_folder_name) / 'skyline.tsv').open('r') as f:
        lines = f.readlines()
        assert len(lines) >= 4
        min_date, min_Ne, *_ = map(float, lines[2].split())
        max_date, max_Ne, *_ = map(float, lines[3].split())

        rho = 1./50.  # Default TreeTime assumption = 50 generations per year; rho = generation time in years
        min_Ne_rho = min_Ne*rho
        max_Ne_rho = max_Ne*rho

        results['ePopSize'] = max_Ne_rho
        results['growthRate'] = math.log(max_Ne/min_Ne) / (max_date - min_date)
                
    return results

def plot_logs(delphy_log_filename, out_pdf_filename, ground_truth_nwk, burnin=0.10, ml_results={}):
    delphy_raw_data = pd.read_table(delphy_log_filename, comment='#')
    num_pts = len(delphy_raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    delphy_data = delphy_raw_data[burnin_pts:]
    print(f'Read in Delphy log in {delphy_log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(delphy_data)} usable')

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
        
        this_ax = sb.kdeplot(adjuster(delphy_data[colname]), fill=True, color='green');
        this_ax.set(xlabel=None)
        this_ax.get_yaxis().set_visible(False)
        if actual_value is not None:
            this_ax.axvline(actual_value, color='green')

        if colname in ml_results:
            (minx, maxx) = this_ax.get_xaxis().get_data_interval()
            (miny, maxy) = this_ax.get_yaxis().get_data_interval()
            value = adjuster(ml_results[colname])
            if 0.5*minx <= value <= 2*maxx:
                this_ax.axvline(value, color='red')
            else:
                ypos = 0.5*(miny+maxy)
                if value < minx:
                    this_ax.annotate("",
                                     xy=(minx, ypos),
                                     xytext=(minx+0.15*(maxx-minx), ypos),
                                     arrowprops=dict(arrowstyle="->", color='red'))
                else:
                    this_ax.annotate("",
                                     xy=(maxx, ypos),
                                     xytext=(maxx-0.15*(maxx-minx), ypos),
                                     arrowprops=dict(arrowstyle="->", color='red'))
            
        return this_ax

    make_subplot( 1, 'clockRate', r'Mutation rate (x$10^{-3}$/site/year)', 1.0, lambda mu: mu*1000)
    make_subplot( 2, 'TreeHeight', r'Tree Height (years)',
                  (datetime.date(2024, 7, 31) - datetime.date.fromisoformat(root_date)).days / 365.0)
    this_ax = make_subplot( 3, 'kappa', r'HKY $\kappa$ parameter', 5.0)
    fig.legend(this_ax.get_legend_handles_labels()[0], ['BEAST2', 'Delphy'], loc='lower right')
    fig.delaxes(ax[0][3])
    
    make_subplot( 5, 'freqParameter.1', r'HKY $\pi_A$', 0.30)
    make_subplot( 6, 'freqParameter.2', r'HKY $\pi_C$', 0.18)
    make_subplot( 7, 'freqParameter.3', r'HKY $\pi_G$', 0.20)
    make_subplot( 8, 'freqParameter.4', r'HKY $\pi_T$', 0.32)

    make_subplot( 9, 'CoalescentExponential', r'Coalescent Prior')
    make_subplot(10, 'ePopSize', r'Final Pop Size (years)',
                 6.0 if delphy_log_filename.startswith('exp') else 2.0)
    make_subplot(11, 'growthRate', r'Pop Growth Rate (1/year)',
                 10.0 if delphy_log_filename.startswith('exp') else 0.0)
    fig.delaxes(ax[2][3])
    
    fig.legend(['Delphy', 'Ground truth'] + (['TreeTime'] if ml_results else []), loc='lower right', 
           bbox_to_anchor=(0.95, 0.05))
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')

for sim in ['exp_100',
            'exp_1000',
            'exp_10000',
            'exp_100000',
            'const_100',
            'const_1000',
            'const_10000',
            'const_100000',
            ]:
    for rep in ('a', 'b'):
        plot_logs(
            f'{sim}/delphy_outputs_{rep}/{sim}.log',
            f'{sim}/delphy_outputs_{rep}/{sim}_distrs.pdf',
            f'{sim}/ground_truth/{sim}_tree.nwk',
            burnin=0.30,
            ml_results=extract_ml_results(f'{sim}/ml_outputs/{sim}.iqtree', f'{sim}/ml_outputs', False)
        )
