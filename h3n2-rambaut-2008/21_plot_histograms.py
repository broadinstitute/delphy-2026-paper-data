#!/usr/bin/env python

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from utils import *

import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import math
import baltic as bt
from pathlib import Path
from matplotlib.ticker import FuncFormatter

def compare_logs(
        logs: list[Log],
        out_pdf_filename,
        burnin=0.30,
        show_alpha=False
):

    datasets = []
    ess_tables = []
    for log in logs:
        raw_data = pd.read_table(log.filename, comment='#')
        num_pts = len(raw_data)
        burnin_pts = math.floor(burnin * num_pts)
        data = raw_data[burnin_pts:]
        print(f'Read in log in {log.filename}')
        print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(data)} usable')
        datasets.append(data)
        ess_tables.append(read_in_ess_table(log.log_analysis_filename, log.tree_ess_filename, log.log_format))

    fig = plt.figure(figsize=(11,8.5), facecolor='w')
    plt.tight_layout()
    gs = fig.add_gridspec(3, 4)
    plt.subplots_adjust(top=0.9, hspace=0.3)

    def make_subplot(i, j, colname, title, adjuster=lambda x: x):
        def translate_colname(colname, log_filename, log_data, log_format):
            #print(f'Looking for {colname} in {log_filename} [-> {log_format[colname]}]')
            candidates = [x for x in log_data.columns if x.startswith(log_format[colname])]
            if len(candidates) == 0:
                raise ValueError(f'Column {colname} not found in {log_filename}')
            if len(candidates) > 1:
                raise ValueError(f'Column {colname} ambiguous in {log_filename} (matches: {candidates})')
            return candidates[0]
            
        this_ax = fig.add_subplot(gs[i, j])
        plt.title(title)
        plt.xlabel(None)

        for log, data, ess_table in zip(logs, datasets, ess_tables):

            # sb.kdeplot will estimate the Gaussian kernel bandwidth using Scott's rule assuming that all samples
            # are independent.  But they are not!  In Scott's rule, the bandwidth scales as n^{-1/5},
            # and we really want it to scale as ESS^{-1/5}, so we have to add a bandwidth-adjustment
            # factor of (n/ESS)^{1/5}
            
            log_colname = translate_colname(colname, log.filename, data, log.log_format)
            n = len(data[log_colname])
            ess = ess_table[colname]
            bw_adjust = (n / ess) ** (1./5.)
            this_ax = sb.kdeplot(adjuster(data[log_colname]), fill=True, color=log.color, bw_adjust=bw_adjust);
            this_ax.set(xlabel=None)
            this_ax.get_yaxis().set_visible(False)

        #if colname == 'gammaShape':
        #    this_ax.set(xlim=[-0.05, 0.50])  # Manual adjustment
            
        return this_ax

    def make_legend():
        fig.legend([log.label for log in logs], loc='lower right', 
                   bbox_to_anchor=(0.95, 0.05))
    
    if show_alpha:
        # Override this_ax to set up the legend properly right below
        this_ax = make_subplot(0, 3, 'gammaShape', r'Heterogeneity Spread $\alpha$')
        make_legend()  # Otherwise, we make it below
    
    make_subplot(0, 0, 'clockRate', r'Mutation rate (x$10^{-3}$/site/year)', lambda mu: mu*1000)
    make_subplot(0, 1, 'TreeHeight', r'Tree Height (years)')
    make_subplot(0, 2, 'kappa', r'HKY $\kappa$ parameter')
        
    make_subplot(1, 0, 'freqParameter.1', r'HKY $\pi_A$')
    make_subplot(1, 1, 'freqParameter.2', r'HKY $\pi_C$')
    make_subplot(1, 2, 'freqParameter.3', r'HKY $\pi_G$')
    make_subplot(1, 3, 'freqParameter.4', r'HKY $\pi_T$')

    #make_subplot(2, 0, 'CoalescentExponential', r'Coalescent Prior')
    #make_subplot(2, 1, 'ePopSize', r'Final Pop Size (years)')
    #make_subplot(2, 2, 'growthRate', r'Pop Growth Rate (1/year)')
    #fig.delaxes(ax[2][3])

    # Skygrid plots
    t0 = bt.decimalDate("2025-08-05")  # Hard-coded, but ok
    def extract_skygrid(data, t0):
        gamma_cols = [x for x in data.columns if x.startswith('skygrid.logPopSize')]
        num_params = len(gamma_cols)
        assert num_params >= 2
        M = num_params - 1

        for cutoff in data['skygrid.cutOff']:
            K = cutoff
            break

        dt = K / M

        ts, lo, median, hi = [], [], [], []
        for i, col in enumerate(gamma_cols):
            ts.append(t0 - i*dt)
            gammas = data[col]
            this_median, hpd = median_and_95HPD(gammas)
            lo.append(hpd[0])
            median.append(this_median)
            hi.append(hpd[1])

        return (ts, np.exp(lo), np.exp(median), np.exp(hi))

    skygrid_results = []
    for log, data in zip(logs, datasets):
        [ts, lo, median, hi] = extract_skygrid(data, t0)
        skygrid_results.append([ts, lo, median, hi])

    ax1 = fig.add_subplot(gs[2, 0:3])
    plt.title('Effective Population Size (years)')
    for log, data, (ts, lo, median, hi) in zip(logs, datasets, skygrid_results):
        ax1.fill_between(ts, lo, hi, color=log.color, alpha=0.2)
        ax1.plot(ts, median, color=log.color, linestyle='-')
    ax1.set_yscale('log')
    ax1.margins(x=0, y=0)
    def log_fmt(y, pos):
        # Format with 1 significant figure, no scientific notation
        rounded = float(f"{y:.1g}")
        if math.floor(rounded) == rounded:
            return '%s' % int(rounded)
        else:
            return '%s' % rounded
    ax1.yaxis.set_major_formatter(FuncFormatter(log_fmt))
    
    if not show_alpha:
        make_legend()    # Otherwise, we make it above
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')

Path('plots').mkdir(parents=True, exist_ok=True)

# Non-alpha datasets
delphy_a_log = Log(
    filename='delphy_outputs_a/h3n2.log',
    log_analysis_filename='delphy_outputs_a/log_analysis.txt',
    tree_ess_filename='delphy_outputs_a/tree_ess.json',
    log_format=delphy_log_format,
    label='Delphy',
    color='green'
)
delphy_b_log = Log(
    filename='delphy_outputs_b/h3n2.log',
    log_analysis_filename='delphy_outputs_b/log_analysis.txt',
    tree_ess_filename='delphy_outputs_b/tree_ess.json',
    log_format=delphy_log_format,
    label='Delphy',
    color='green'
)
beastX_a_log = Log(
    filename='beastX_run_a/output.log',
    log_analysis_filename='beastX_run_a/log_analysis.txt',
    tree_ess_filename='beastX_run_a/tree_ess.json',
    log_format=beastX_log_format,
    label='BEAST X',
    color='blue'
)
beastX_b_log = Log(
    filename='beastX_run_b/output.log',
    log_analysis_filename='beastX_run_b/log_analysis.txt',
    tree_ess_filename='beastX_run_b/tree_ess.json',
    log_format=beastX_log_format,
    label='BEAST X',
    color='blue'
)

compare_logs(
    [delphy_a_log, beastX_a_log],
    'plots/DelphyAVsBeastXA.pdf',
    burnin=0.30
)
compare_logs(
    [delphy_b_log, beastX_b_log],
    'plots/DelphyBVsBeastXB.pdf',
    burnin=0.30
)

# Alpha datasets
delphy_alpha_a_log = Log(
    filename='delphy_outputs_alpha_a/h3n2_alpha.log',
    log_analysis_filename='delphy_outputs_alpha_a/log_analysis.txt',
    tree_ess_filename='delphy_outputs_alpha_a/tree_ess.json',
    log_format=delphy_log_format,
    label=r'Delphy (K=$\infty$)',
    color='green'
)
delphy_alpha_b_log = Log(
    filename='delphy_outputs_alpha_b/h3n2_alpha.log',
    log_analysis_filename='delphy_outputs_alpha_b/log_analysis.txt',
    tree_ess_filename='delphy_outputs_alpha_b/tree_ess.json',
    log_format=delphy_log_format,
    label=r'Delphy (K=$\infty$)',
    color='green'
)
beastX_alpha_a_log = Log(
    filename='beastX_run_alpha_a/output.log',
    log_analysis_filename='beastX_run_alpha_a/log_analysis.txt',
    tree_ess_filename='beastX_run_alpha_a/tree_ess.json',
    log_format=beastX_log_format,
    label='BEAST X (K=4)',
    color='blue'
)
beastX_alpha_b_log = Log(
    filename='beastX_run_alpha_b/output.log',
    log_analysis_filename='beastX_run_alpha_b/log_analysis.txt',
    tree_ess_filename='beastX_run_alpha_b/tree_ess.json',
    log_format=beastX_log_format,
    label='BEAST X (K=4)',
    color='blue'
)

compare_logs(
    [delphy_alpha_a_log, beastX_alpha_a_log],
    'plots/DelphyAVsBeastXAAlpha.pdf',
    burnin=0.30,
    show_alpha=True
)
compare_logs(
    [delphy_alpha_b_log, beastX_alpha_b_log],
    'plots/DelphyBVsBeastXBAlpha.pdf',
    burnin=0.30,
    show_alpha=True
)
