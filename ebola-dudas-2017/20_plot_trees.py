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

# Read in metadata
# ===========================
print("\nReading in metadata...")
sample_id_2_country = {}
with open('delphy_inputs/ebola_dudas_metadata.csv', 'r') as f:
    r = csv.DictReader(f)
    for record in r:
        sample_id_2_country[record['id']] = record['Country']

# Country color map
# $ cut -d',' -f2 delphy_inputs/ebola_dudas_metadata.csv | sort | uniq -c | sort -nr
#    1031 Sierra Leone
#     368 Guinea
#     209 Liberia
#       2 -
#       1 Country
        
# Colors to match Fig 1 of Holmes et al (2017) - doi: 10.1038/nature19790
country_2_color = {
    'Sierra Leone': '#5275b4',
    'Guinea':       '#5da56f',
    'Liberia':      '#da877a',
    '-':            'lightgrey',
}

def country_of(leaf_name):
    sample_id, *rest = leaf_name.split('|')
    return sample_id_2_country.get(sample_id, "-")

def color_of(leaf_name):
    return country_2_color[country_of(leaf_name)]

def simple_parsimony(tree, tip_state_assigner):
    # Post-order traversal implementing Fitch algorithm
    index_2_fs = {}
    def go(node):
        if not node.is_leaf():
            for child in node.children:
                go(child)
        
        if node.is_leaf():
            fs = tip_state_assigner(node)
            if isinstance(fs, str):
                fs = frozenset([fs])  # Singleton set
            fs = frozenset(fs)  # Make lists, sets and frozensets into just frozensets
        else:
            child_fss = [index_2_fs[child.index] for child in node.children]

            # Fitch rule: intersection if nonempty, otherwise union
            fs = frozenset.intersection(*child_fss)
            if len(fs) == 0:
                fs = frozenset.union(*child_fss)

        index_2_fs[node.index] = fs
    go(tree.root)

    # Pre-order traversal to assign states (pick "lowest" state on ties according to natural sort order of keys)
    index_2_state = {}
    def gogo(node, parent_state):
        node_fs = index_2_fs[node.index]
        
        if parent_state in node_fs:
            node_state = parent_state
        else:
            node_state = sorted(list(node_fs))[0]
            if len(node_fs) > 1:
                print(f'WARNING: Arbitrarily breaking tie between {node_fs} in favor of {node_state}')
        index_2_state[node.index] = node_state

        if not node.is_leaf():
            for child in node.children:
                gogo(child, node_state)
    gogo(tree.root, None)

    return index_2_state


def extract_tMRCA_distr(log_filename, burnin, t0):
    # Extremely brittle way of extracting a tMRCA distribution from the log files
    raw_data = pd.read_table(log_filename, comment='#')
    num_pts = len(raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    data = raw_data[burnin_pts:]
    print(f'Read in log file {log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(data)} usable')

    return [t0 - h for h in data['rootHeight']]


def plot_tree(mcc_tree, out_pdf_filename, tMRCAs, tMRCA_ess, legend=False):
    leaves = mcc_tree.traverse_tree()
    max_date = max(leaf.name.split('|')[-1] for leaf in leaves)  # as of this writing, "2015-10-24"
    mcc_tree.setAbsoluteTime(bt.decimalDate(max_date))

    index_2_country = simple_parsimony(mcc_tree, lambda node: country_of(node.name))
    
    fig,ax = plt.subplots(figsize=(6,6),facecolor='w', constrained_layout=True)

    x_attr=lambda k: k.absoluteTime
    c_func=lambda k: 'black' if k.branchType != 'leaf' else color_of(k.name)
    s_func=lambda k: 15 if (k.branchType == 'leaf') else 1
    
    mcc_tree.sortBranches(descending=True)
    #for k in mcc_tree.getInternal():
    #    k.children.reverse()
    mcc_tree.drawTree()
    
    mcc_tree.plotTree(ax,
                      x_attr=x_attr,
                      colour=lambda node: country_2_color[index_2_country[node.index]],
                      width=1)
    mcc_tree.plotPoints(
        ax,
        x_attr=x_attr,
        size=s_func,
        colour=c_func,
        zorder=100,
        outline=True)
    
    target_func=lambda k: k.branchType=='node' and 'posterior' in k.traits and k.traits['posterior'] > 0.95 ## only target high-posterior-support nodes

    mcc_tree.plotPoints(
        ax,
        x_attr=x_attr,
        target=target_func,
        size=5,
        colour="white",
        outline=True,
        outline_size=10,
        outline_colour='black',
        zorder=120
    )

    if tMRCAs:
        # stats.gaussian_kde will estimate the Gaussian kernel bandwidth using Scott's rule assuming that all samples
        # are independent.  But they are not!  In Scott's rule, the bandwidth scales as n^{-1/5},
        # and we really want it to scale as ESS^{-1/5}
        mean_tMRCA = np.mean(tMRCAs)  # Slightly different from root height because tMRCA was sampled much more often than trees
        kde = stats.gaussian_kde(tMRCAs, bw_method=tMRCA_ess ** (-1./5.))
        tt = np.linspace(bt.decimalDate('2013-09-01'), bt.decimalDate('2014-03-01'), 200)
        ax.plot(tt, -20+60*kde(tt), color='#888888', lw=0.5)
        ax.fill_between(tt, -20+60*kde(tt), -20, color='#000000', alpha=0.2)
    
    ax.set_xlim(bt.decimalDate("2013-08-01"), bt.decimalDate("2015-11-01"));
    ax.set_ylim(-20, mcc_tree.ySpan+150);
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);
    ax.get_xaxis().set_ticks([
        bt.decimalDate("2013-10-01"),
        bt.decimalDate("2014-01-01"),
        bt.decimalDate("2014-04-01"),
        bt.decimalDate("2014-07-01"),
        bt.decimalDate("2014-10-01"),
        bt.decimalDate("2015-01-01"),
        bt.decimalDate("2015-04-01"),
        bt.decimalDate("2015-07-01"),
        bt.decimalDate("2015-10-01"),
       ])
    ax.set_xticklabels([
        "1 Oct\n",
        "1 Jan\n2014",
        "1 Apr\n",
        "1 Jul\n",
        "1 Oct\n",
        "1 Jan\n2015",
        "1 Apr\n",
        "1 Jul\n",
        "1 Oct\n",
       ]);
    #ax.set_xlabel('Time (years)')
    
    plt.axvspan(bt.decimalDate("2013-10-01"), bt.decimalDate("2013-11-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2013-12-01"), bt.decimalDate("2014-01-01"), color='black', alpha=0.05)
    
    plt.axvspan(bt.decimalDate("2014-02-01"), bt.decimalDate("2014-03-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2014-04-01"), bt.decimalDate("2014-05-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2014-06-01"), bt.decimalDate("2014-07-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2014-08-01"), bt.decimalDate("2014-09-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2014-10-01"), bt.decimalDate("2014-11-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2014-12-01"), bt.decimalDate("2015-01-01"), color='black', alpha=0.05)

    plt.axvspan(bt.decimalDate("2015-02-01"), bt.decimalDate("2015-03-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2015-04-01"), bt.decimalDate("2015-05-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2015-06-01"), bt.decimalDate("2015-07-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2015-08-01"), bt.decimalDate("2015-09-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2015-10-01"), bt.decimalDate("2015-11-01"), color='black', alpha=0.05)

    # Legend code adapted from https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html
    if legend:
        custom_lines = []
        custom_labels = []
        for country, color in country_2_color.items():
            if country == '-':
                continue    # Suppress extremely rare but confusing category in legend
            custom_lines.append(Line2D([0], [0], color=color, lw=3))
            custom_labels.append(country)
        ax.legend(custom_lines, custom_labels, bbox_to_anchor=(0.69, 0.79), framealpha=0.5)
    
    def annotate_node(node, text, diamond=False, **kwargs):
        x = kwargs.get('x', x_attr(node))
        y = node.y

        # Allow kwargs to override defaults
        allkwargs = dict(horizontalalignment='right', verticalalignment='bottom') | kwargs
        if 'x' in allkwargs:
            del allkwargs['x']
        plt.text(x-0.03, y-300.05, text, **allkwargs)
        if diamond:
            ax.plot(x, y, marker='D', color='black')

    if tMRCAs:
        annotate_node(mcc_tree.root,
                      '',   #'tMRCA',
                      diamond=True, x=mean_tMRCA, rotation='vertical', verticalalignment='top')
        
        ax.plot([mean_tMRCA, mean_tMRCA], [mcc_tree.root.y, -20], color='black', linestyle='--', lw=0.75)
    
    plt.savefig(out_pdf_filename, bbox_inches='tight', pad_inches=0)
    print(f'Plot saved to {out_pdf_filename}')

Path('plots').mkdir(parents=True, exist_ok=True)

data = [
    # data_dir                  mcc_filename                    log_filename                    plot_filename
    ('delphy_outputs_a',       'ebola_dudas_delphy.mcc',       'ebola_dudas_delphy.log',       'DelphyAMcc.pdf'),
    ('delphy_outputs_b',       'ebola_dudas_delphy.mcc',       'ebola_dudas_delphy.log',       'DelphyBMcc.pdf'),
    ('beastX_run_a',           'ebola_dudas_beastX.mcc',       'output.log',                   'BeastXAMcc.pdf'),
    ('beastX_run_b',           'ebola_dudas_beastX.mcc',       'output.log',                   'BeastXBMcc.pdf'),
    ('delphy_outputs_alpha_a', 'ebola_dudas_alpha_delphy.mcc', 'ebola_dudas_alpha_delphy.log', 'DelphyAAlphaMcc.pdf'),
    ('delphy_outputs_alpha_b', 'ebola_dudas_alpha_delphy.mcc', 'ebola_dudas_alpha_delphy.log', 'DelphyBAlphaMcc.pdf'),
    ('beastX_run_alpha_a',     'ebola_dudas_alpha_beastX.mcc', 'output.log',                   'BeastXAAlphaMcc.pdf'),
    ('beastX_run_alpha_b',     'ebola_dudas_alpha_beastX.mcc', 'output.log',                   'BeastXBAlphaMcc.pdf'),
]

for (data_dir, mcc_filename, log_filename, plot_filename) in data:
    mcc_tree = bt.loadNexus(f'{data_dir}/{mcc_filename}')
    leaves = mcc_tree.traverse_tree()
    max_date = max(leaf.name.split('|')[-1] for leaf in leaves)  # e.g., 2015-10-24
    
    tMRCAs = extract_tMRCA_distr(f'{data_dir}/{log_filename}', burnin=0.30, t0=bt.decimalDate(max_date))
    ess_table = read_in_ess_table(f'{data_dir}/log_analysis.txt', None,
                                  delphy_log_format if data_dir.startswith('delphy') else beastX_log_format)
    plot_tree(mcc_tree, f'plots/{plot_filename}', tMRCAs, ess_table['TreeHeight'], legend=data_dir.startswith('delphy'))

plot_tree(bt.loadNexus('./ml/tt/timetree.nexus'), 'plots/ML.pdf', None, None, legend=False)
plot_tree(bt.loadNexus('./ml/tt_alpha/timetree.nexus'), 'plots/MLAlpha.pdf', None, None, legend=False)
