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

def mrca_node_of(tree, leaf_names):
    leaf_nodes = [node for node in tree.Objects if node.is_leaf() and node.name.split('|')[0] in leaf_names]
    if len(leaf_nodes) == 1:
        return leaf_nodes[0]
    else:
        return tree.commonAncestor(leaf_nodes)

def extract_tMRCA_distr(log_filename, burnin):
    # Extremely brittle way of extracting a tMRCA distribution from the log files
    raw_data = pd.read_table(log_filename, comment='#')
    num_pts = len(raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    data = raw_data[burnin_pts:]
    print(f'Read in log file {log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(data)} usable')

    return list(data['age(root)'])


def plot_tree(mcc_filename, out_pdf_filename, tMRCAs):
    mcc_tree = bt.loadNexus(mcc_filename)
    mcc_tree.traverse_tree()
    mcc_tree.setAbsoluteTime(bt.decimalDate("2022-08-15"))  # Hard-coded, but ok

    fig,ax = plt.subplots(figsize=(5,4),facecolor='w')

    x_attr=lambda k: k.absoluteTime
    c_func='black'
    s_func=lambda k: 20 if (k.branchType == 'leaf') else 1
    
    mcc_tree.sortBranches(descending=True)
    #for k in mcc_tree.getInternal():  # Match ordering in O'Toole et al 2023
    #    k.children.reverse()
    mcc_tree.drawTree()

    mcc_tree.plotTree(ax,
                      x_attr=x_attr,
                      colour='black',
                      width=0.5)
    mcc_tree.plotPoints(
        ax,
        x_attr=x_attr,
        size=s_func,
        colour=c_func,
        zorder=100,
        outline=False)

    mean_tMRCA = np.mean(tMRCAs)  # Slightly different from root height because tMRCA was sampled much more often than trees
    kde = stats.gaussian_kde(tMRCAs)
    tt = np.linspace(2013.5, 2018.0, 200)
    ax.plot(tt, -5+10*kde(tt), color='#888888', lw=0.5)
    ax.fill_between(tt, -5+10*kde(tt), -5, color='#000000', alpha=0.2)

    ax.set_xlim(2014.5, 2022.7)
    ax.set_ylim(-5, mcc_tree.ySpan+1);
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([
        -5+10*0.00,
        -5+10*0.50,
        -5+10*1.00,
        -5+10*1.50,
    ]);
    ax.get_yaxis().set_ticklabels([
        '0',
        '50',
        '100',
        '150',
    ]);
    ax.set_ylabel('Density', loc='bottom')
    ax.yaxis.set_label_coords(-0.12, 0.12)  # Trial and error
    
    ax.get_xaxis().set_ticks([
        bt.decimalDate("2015-01-01"),
        bt.decimalDate("2016-01-01"),
        bt.decimalDate("2017-01-01"),
        bt.decimalDate("2018-01-01"),
        bt.decimalDate("2019-01-01"),
        bt.decimalDate("2020-01-01"),
        bt.decimalDate("2021-01-01"),
        bt.decimalDate("2022-01-01"),
        #bt.decimalDate("2023-01-01"),
       ])
    ax.set_xticklabels([
        "2015",
        "2016",
        "2017",
        "2018",
        "2019",
        "2020",
        "2021",
        "2022",
        #"2023",
       ]);
    ax.set_xlabel('Time (years)')

    plt.axvspan(2014.0, 2015.0, color='black', alpha=0.05)
    plt.axvspan(2016.0, 2017.0, color='black', alpha=0.05)
    plt.axvspan(2018.0, 2019.0, color='black', alpha=0.05)
    plt.axvspan(2020.0, 2021.0, color='black', alpha=0.05)
    plt.axvspan(2022.0, 2023.0, color='black', alpha=0.05)

    def annotate_node(node, text, diamond=False, **kwargs):
        x = kwargs.get('x', x_attr(node))
        y = node.y

        # Allow kwargs to override defaults
        allkwargs = dict(horizontalalignment='right', verticalalignment='bottom') | kwargs
        if 'x' in allkwargs:
            del allkwargs['x']
        plt.text(x-0.1, y, text, **allkwargs)
        if diamond:
            ax.plot(x, y, marker='D', color='black')

    annotate_node(mrca_node_of(mcc_tree, ['ON563414']), 'B.1')
    annotate_node(mrca_node_of(mcc_tree, ['ON676708', 'ON563414']), 'A.1.1')
    annotate_node(mrca_node_of(mcc_tree, ['ON676708', 'OP612681']), 'A.1')
    annotate_node(mrca_node_of(mcc_tree, ['ON674051', 'ON676707']), 'A.2')

    annotate_node(mcc_tree.root, 'tMRCA', diamond=True, x=mean_tMRCA, rotation='vertical', verticalalignment='top')

    ax.plot([mean_tMRCA, mean_tMRCA], [mcc_tree.root.y, -5], color='black', linestyle='--', lw=0.75)
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')


Path('plots').mkdir(parents=True, exist_ok=True)

for run in ('a', 'b'):
    tMRCAs = extract_tMRCA_distr(f'./delphy_outputs_{run}/mpox-otoole-2023.log', burnin=0.30)
    plot_tree(f'./delphy_outputs_{run}/mpox-otoole-2023.mcc', f'plots/{run.upper()}DelphyMcc.pdf', tMRCAs)
