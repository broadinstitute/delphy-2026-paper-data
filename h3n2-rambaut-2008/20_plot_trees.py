#!/usr/bin/env python

import baltic as bt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import csv
from scipy import stats
import pandas as pd
import math
import numpy as np

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from utils import *

def extract_tMRCA_distr(log_filename, burnin, t0):
    # Extremely brittle way of extracting a tMRCA distribution from the log files
    raw_data = pd.read_table(log_filename, comment='#')
    num_pts = len(raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    data = raw_data[burnin_pts:]
    print(f'Read in log file {log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(data)} usable')

    return [t0 - h for h in data['rootHeight']]


def plot_tree(mcc_tree, out_pdf_filename, tMRCAs, tMRCA_ess):
    leaves = mcc_tree.traverse_tree()
    max_date = max(leaf.name.split('|')[-1] for leaf in leaves)  # as of this writing, "2015-10-24"
    mcc_tree.setAbsoluteTime(bt.decimalDate(max_date))

    fig,ax = plt.subplots(figsize=(6,6),facecolor='w')

    x_attr=lambda k: k.absoluteTime
    c_func=lambda k: 'black'
    s_func=lambda k: 15 if (k.branchType == 'leaf') else 1
    
    mcc_tree.sortBranches(descending=True)
    #for k in mcc_tree.getInternal():
    #    k.children.reverse()
    mcc_tree.drawTree()
    
    mcc_tree.plotTree(ax,
                      x_attr=x_attr,
                      colour='black',
                      width=1)
    mcc_tree.plotPoints(
        ax,
        x_attr=x_attr,
        size=s_func,
        colour=c_func,
        zorder=100,
        outline=False)
    
    # stats.gaussian_kde will estimate the Gaussian kernel bandwidth using Scott's rule assuming that all samples
    # are independent.  But they are not!  In Scott's rule, the bandwidth scales as n^{-1/5},
    # and we really want it to scale as ESS^{-1/5}
            
    mean_tMRCA = np.mean(tMRCAs)  # Slightly different from root height because tMRCA was sampled much more often than trees
    kde = stats.gaussian_kde(tMRCAs, bw_method=tMRCA_ess ** (-1./5.))
    tt = np.linspace(bt.decimalDate('1995-01-01'), bt.decimalDate('2000-01-01'), 200)
    ax.plot(tt, -20+60*kde(tt), color='#888888', lw=0.5)
    ax.fill_between(tt, -20+60*kde(tt), -20, color='#000000', alpha=0.2)
    
    ax.set_xlim(1995, 2004.1);
    ax.set_ylim(-20, mcc_tree.ySpan+20);
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);
    ax.get_xaxis().set_ticks([
        bt.decimalDate("1995-01-01"),
        bt.decimalDate("1996-01-01"),
        bt.decimalDate("1997-01-01"),
        bt.decimalDate("1998-01-01"),
        bt.decimalDate("1999-01-01"),
        bt.decimalDate("2000-01-01"),
        bt.decimalDate("2001-01-01"),
        bt.decimalDate("2002-01-01"),
        bt.decimalDate("2003-01-01"),
        bt.decimalDate("2004-01-01"),
       ], minor=True)    # Dirty trick to place tick marks at start of year but label mid-year
    ax.get_xaxis().set_ticks([
        bt.decimalDate("1995-07-01"),
        bt.decimalDate("1996-07-01"),
        bt.decimalDate("1997-07-01"),
        bt.decimalDate("1998-07-01"),
        bt.decimalDate("1999-07-01"),
        bt.decimalDate("2000-07-01"),
        bt.decimalDate("2001-07-01"),
        bt.decimalDate("2002-07-01"),
        bt.decimalDate("2003-07-01"),
       ])    # Dirty trick to place tick marks at start of year but label mid-year
    ax.tick_params(axis='x', which='major', length=0)
    ax.tick_params(axis='x', which='minor', length=5)
    ax.set_xticklabels([
        "1995",
        "1996",
        "1997",
        "1998",
        "1999",
        "2000",
        "2001",
        "2002",
        "2003",
       ]);
    ax.set_xlabel('Time (years)')
    
    plt.axvspan(bt.decimalDate("1996-01-01"), bt.decimalDate("1997-01-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("1998-01-01"), bt.decimalDate("1999-01-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2000-01-01"), bt.decimalDate("2001-01-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2002-01-01"), bt.decimalDate("2003-01-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2004-01-01"), bt.decimalDate("2005-01-01"), color='black', alpha=0.05)
    
    # plt.axvspan(bt.decimalDate("2014-01-01"), bt.decimalDate("2014-02-01"), color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2014-03-01"), bt.decimalDate("2014-04-01"), color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2014-05-01"), bt.decimalDate("2014-06-01"), color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2014-07-01"), bt.decimalDate("2014-08-01"), color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2014-09-01"), bt.decimalDate("2014-10-01"), color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2014-11-01"), bt.decimalDate("2014-12-01"), color='black', alpha=0.05)

    # plt.axvspan(bt.decimalDate("2015-01-01"), bt.decimalDate("2015-02-01"), color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2015-03-01"), bt.decimalDate("2015-04-01"), color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2015-05-01"), bt.decimalDate("2015-06-01"), color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2015-07-01"), bt.decimalDate("2015-08-01"), color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2015-09-01"), bt.decimalDate("2015-10-01"), color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2015-11-01"), bt.decimalDate("2015-12-01"), color='black', alpha=0.05)

    def annotate_node(node, text, diamond=False, **kwargs):
        x = kwargs.get('x', x_attr(node))
        y = node.y

        # Allow kwargs to override defaults
        allkwargs = dict(horizontalalignment='right', verticalalignment='bottom') | kwargs
        if 'x' in allkwargs:
            del allkwargs['x']
        plt.text(x-0.03, y, text, **allkwargs)
        if diamond:
            ax.plot(x, y, marker='D', color='black')

    annotate_node(mcc_tree.root, 'tMRCA', diamond=True, x=mean_tMRCA, rotation='vertical', verticalalignment='top')
    
    ax.plot([mean_tMRCA, mean_tMRCA], [mcc_tree.root.y, -20], color='black', linestyle='--', lw=0.75)
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')

Path('plots').mkdir(parents=True, exist_ok=True)

data = [
    # data_dir                  mcc_filename      log_filename      plot_filename
    ('delphy_outputs_a',       'h3n2.mcc',       'h3n2.log',       'DelphyAMcc.pdf'),
    ('delphy_outputs_b',       'h3n2.mcc',       'h3n2.log',       'DelphyBMcc.pdf'),
    ('beastX_run_a',           'h3n2.mcc',       'output.log',     'BeastXAMcc.pdf'),
    ('beastX_run_b',           'h3n2.mcc',       'output.log',     'BeastXBMcc.pdf'),
    ('delphy_outputs_alpha_a', 'h3n2_alpha.mcc', 'h3n2_alpha.log', 'DelphyAAlphaMcc.pdf'),
    ('delphy_outputs_alpha_b', 'h3n2_alpha.mcc', 'h3n2_alpha.log', 'DelphyBAlphaMcc.pdf'),
    ('beastX_run_alpha_a',     'h3n2_alpha.mcc', 'output.log',     'BeastXAAlphaMcc.pdf'),
    ('beastX_run_alpha_b',     'h3n2_alpha.mcc', 'output.log',     'BeastXBAlphaMcc.pdf'),
]

for (data_dir, mcc_filename, log_filename, plot_filename) in data:
    mcc_tree = bt.loadNexus(f'{data_dir}/{mcc_filename}')
    leaves = mcc_tree.traverse_tree()
    max_date = max(leaf.name.split('|')[-1] for leaf in leaves)  # e.g., 2015-10-24
    
    tMRCAs = extract_tMRCA_distr(f'{data_dir}/{log_filename}', burnin=0.30, t0=bt.decimalDate(max_date))
    ess_table = read_in_ess_table(f'{data_dir}/log_analysis.txt', None,
                                  delphy_log_format if data_dir.startswith('delphy') else beastX_log_format)
    plot_tree(mcc_tree, f'plots/{plot_filename}', tMRCAs, ess_table['TreeHeight'])
