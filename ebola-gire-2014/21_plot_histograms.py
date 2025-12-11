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
        show_alpha=False,
        ml_results={}
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

    fig = plt.figure(figsize=(9,6.5), facecolor='w')
    plt.subplots_adjust(top=0.9, hspace=0.6, wspace=0.1)
    plt.tight_layout()
    gs = fig.add_gridspec(3, 4, hspace=0.6, top=0.9)

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

    def make_legend():
        fig.legend([log.label for log in logs] + (['ML' + (' (K=4)' if show_alpha else '')] if ml_results else []), loc='lower right', 
                   bbox_to_anchor=(0.90, 0.10), fontsize=10)
    
    if show_alpha:
        # Override this_ax to set up the legend properly right below
        this_ax = make_subplot(0, 3, 'gammaShape', 'Heterogeneity\nSpread $\\alpha$')
        make_legend()  # Otherwise, we make it below
    
    make_subplot(0, 0, 'clockRate', 'Mutation rate\n(x$10^{-3}$/site/year)', lambda mu: mu*1000)
    make_subplot(0, 1, 'TreeHeight', 'Tree Height\n(years)')
    make_subplot(0, 2, 'kappa', r'HKY $\kappa$')
        
    make_subplot(1, 0, 'freqParameter.1', 'HKY $\\pi_A$')
    make_subplot(1, 1, 'freqParameter.2', 'HKY $\\pi_C$')
    make_subplot(1, 2, 'freqParameter.3', 'HKY $\\pi_G$')
    make_subplot(1, 3, 'freqParameter.4', 'HKY $\\pi_T$')

    make_subplot(2, 0, 'CoalescentExponential', r'Coalescent Prior')
    make_subplot(2, 1, 'ePopSize', 'Final Pop Size\n(years)')
    make_subplot(2, 2, 'growthRate', 'Pop Growth Rate\n(1/year)')

    if not show_alpha:
        make_legend()    # Otherwise, we make it above
        
    plt.savefig(out_pdf_filename, bbox_inches='tight')
    print(f'Plot saved to {out_pdf_filename}')

Path('plots').mkdir(parents=True, exist_ok=True)

# Non-alpha datasets
ml_results = extract_ml_results('ml/ebola.iqtree', 'ml/tt', False)
delphy_a_log = Log(
    filename='delphy_outputs_a/ebola_delphy.log',
    log_analysis_filename='delphy_outputs_a/log_analysis.txt',
    tree_ess_filename='delphy_outputs_a/tree_ess.json',
    log_format=delphy_log_format,
    label='Delphy',
    color='green'
)
delphy_b_log = Log(
    filename='delphy_outputs_b/ebola_delphy.log',
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
    [delphy_a_log, beast2_log, beastX_log],
    'plots/DelphyAVsBeast2VsBeastX.pdf',
    burnin=0.30,
    ml_results=ml_results
)
compare_logs(
    [delphy_b_log, beast2_log, beastX_log],
    'plots/DelphyBVsBeast2VsBeastX.pdf',
    burnin=0.30,
    ml_results=ml_results
)

# Alpha datasets
ml_results_alpha = extract_ml_results('ml/ebola_alpha.iqtree', 'ml/tt_alpha', True)
delphy_alpha_a_log = Log(
    filename='delphy_outputs_alpha_a/ebola_delphy_alpha.log',
    log_analysis_filename='delphy_outputs_alpha_a/log_analysis.txt',
    tree_ess_filename='delphy_outputs_alpha_a/tree_ess.json',
    log_format=delphy_log_format,
    label=r'Delphy (K=$\infty$)',
    color='green'
)
delphy_alpha_b_log = Log(
    filename='delphy_outputs_alpha_b/ebola_delphy_alpha.log',
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
compare_logs(
    [delphy_alpha_a_log, beast2_alpha_log, beastX_alpha_log],
    'plots/DelphyAVsBeast2VsBeastXAlpha.pdf',
    burnin=0.30,
    show_alpha=True,
    ml_results=ml_results_alpha
)
compare_logs(
    [delphy_alpha_b_log, beast2_alpha_log, beastX_alpha_log],
    'plots/DelphyBVsBeast2VsBeastXAlpha.pdf',
    burnin=0.30,
    show_alpha=True,
    ml_results=ml_results_alpha
)

