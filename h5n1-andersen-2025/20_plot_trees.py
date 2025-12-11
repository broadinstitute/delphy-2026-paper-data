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

# Oops!  Somehow 4 sequences without exact dates slipped through; ignore them
BAD_LEAVES = [
    'A/cattle/ID/25-014338-001-original/2025',
    'A/cattle/CA/24-033172-001-original/2024',
    'A/cattle/CA/24-033172-003-original/2024',
    'A/cattle/CA/24-033172-002-original/2024',
]

# Marked in outliers.tsv by TimeTree
OUTLIERS = [
    'A/cattle/CA/25-010022-002-original/2025|SRR33029766|2025-03-18',
    'A/cattle/CA/25-015059-001-original/2025|SRR33764531|2025-05-06',
]

# Read in metadata
# ===========================
print("\nReading in metadata...")
sample_id_2_state = {}
with open('delphy_inputs/h5n1-andersen-e756a15_metadata.csv', 'r') as f:
    r = csv.DictReader(f)
    for record in r:
        sample_id_2_state[record['id']] = record['geo']

# State color map
# $ cut -d',' -f2 delphy_inputs/h5n1-andersen-ebfdf65_metadata.csv | sort | uniq -c | sort -nr
#     962 -          # Almost all (all?) of these are SRRs that have no GenBank counterpart
#     808 CA
#     271 CO
#     240 TX
#     142 ID
#     108 MI
#      94 MN
#      56 IA
#      45 WY
#      42 SD
#      40 NM
#      20 OH
#      11 KS
#       4 OK
#       4 NC
#       3 UT
#       1 geo

state_2_color = {
    'CA': 'blue',
    'CO': 'orange',
    'TX': 'red',
    'ID': 'green',    
    'MI': 'yellow',
    'MN': 'cyan',
    'IA': 'darkgrey',
    'WY': 'grey',
    'SD': 'magenta',
    'NM': 'teal',
    'OH': 'purple',
    'KS': 'lightgreen',
    'OK': 'pink',
    'NC': 'brown',
    'UT': 'black',
    '-':  'lightgrey',
}

def state_of(leaf_name):
    sample_id, *rest = leaf_name.split('|')
    return sample_id_2_state.get(sample_id, "-")

def color_of(leaf_name):
    return state_2_color[state_of(leaf_name)]

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

    # Sadly, one of the "outliers" detected by TimeTree, whose "apparent date" is far later than
    # its stated date, _is_ the latest tip in the tree, so we need to be more careful here
    for leaf in leaves:
        if leaf.name not in OUTLIERS:
            leaf_date_str = leaf.name.split('|')[-1]
            leaf_date = bt.decimalDate(leaf_date_str)
            abs_time_delta = leaf_date - leaf.height
            break

    for k in mcc_tree.Objects:
        k.absoluteTime = abs_time_delta + k.height  # moral equivalent of mcc_tree.setAbsoluteTime(max_date)
    
    #max_date = max(leaf.name.split('|')[-1] for leaf in leaves if leaf.name not in BAD_LEAVES)  # e.g., 2020-03-01
    #mcc_tree.setAbsoluteTime(bt.decimalDate(max_date))

    index_2_cat = simple_parsimony(mcc_tree, lambda node: state_of(node.name))
    
    fig,ax = plt.subplots(figsize=(6,6),facecolor='w', constrained_layout=True)

    x_attr=lambda k: k.absoluteTime
    c_func=lambda k: 'black' if k.branchType != 'leaf' else color_of(k.name)
    s_func=lambda k: 15 if (k.branchType == 'leaf') else 1
    
    mcc_tree.sortBranches(descending=True)
    mcc_tree.drawTree()
    
    mcc_tree.plotTree(ax,
                      x_attr=x_attr,
                      colour=lambda node: state_2_color[index_2_cat[node.index]],
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

    target_func=lambda k: k.branchType=='leaf' and k.name in OUTLIERS

    mcc_tree.plotPoints(
        ax,
        x_attr=x_attr,
        target=target_func,
        size=35,
        colour='white',
        outline=True,
        outline_size=50,
        outline_colour='red',
        zorder=90
    )

    if tMRCAs:
        # stats.gaussian_kde will estimate the Gaussian kernel bandwidth using Scott's rule assuming that all samples
        # are independent.  But they are not!  In Scott's rule, the bandwidth scales as n^{-1/5},
        # and we really want it to scale as ESS^{-1/5}
        mean_tMRCA = np.mean(tMRCAs)  # Slightly different from root height because tMRCA was sampled much more often than trees
        kde = stats.gaussian_kde(tMRCAs, bw_method=tMRCA_ess ** (-1./5.))
        tt = np.linspace(bt.decimalDate('2023-10-01'), bt.decimalDate('2024-03-01'), 200)
        ax.plot(tt, -20+60*kde(tt), color='#888888', lw=0.5)
        ax.fill_between(tt, -20+60*kde(tt), -20, color='#000000', alpha=0.2)
    
    ax.set_xlim(bt.decimalDate("2023-10-01"), bt.decimalDate("2025-11-01"))
    ax.set_ylim(-20, mcc_tree.ySpan+300);
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);
    # ax.get_yaxis().set_ticks([
    #     -5+10*0.00,
    #     -5+10*0.50,
    #     -5+10*1.00,
    #     -5+10*1.50,
    # ]);
    # ax.get_yaxis().set_ticklabels([
    #     '0',
    #     '50',
    #     '100',
    #     '150',
    # ]);
    # ax.set_ylabel('Density', loc='bottom')
    # ax.yaxis.set_label_coords(-0.12, 0.12)  # Trial and error

    ax.get_xaxis().set_ticks([
        bt.decimalDate("2023-10-01"),
        bt.decimalDate("2024-01-01"),
        bt.decimalDate("2024-04-01"),
        bt.decimalDate("2024-07-01"),
        bt.decimalDate("2024-10-01"),
        bt.decimalDate("2025-01-01"),
        bt.decimalDate("2025-04-01"),
        bt.decimalDate("2025-07-01"),
        bt.decimalDate("2025-10-01"),
    ])
    ax.set_xticklabels([
        "1 Oct\n",
        "1 Jan\n2024",
        "1 Apr\n",
        "1 Jul\n",
        "1 Oct\n",
        "1 Jan\n2025",
        "1 Apr\n",
        "1 Jul\n",
        "1 Oct\n",
    ]);
    ax.set_xlabel('Time (years)')

    plt.axvspan(bt.decimalDate("2023-10-01"), bt.decimalDate("2023-11-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2023-12-01"), bt.decimalDate("2024-01-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-02-01"), bt.decimalDate("2024-03-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-04-01"), bt.decimalDate("2024-05-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-06-01"), bt.decimalDate("2024-07-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-08-01"), bt.decimalDate("2024-09-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-10-01"), bt.decimalDate("2024-11-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-12-01"), bt.decimalDate("2025-01-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2025-02-01"), bt.decimalDate("2025-03-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2025-04-01"), bt.decimalDate("2025-05-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2025-06-01"), bt.decimalDate("2025-07-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2025-08-01"), bt.decimalDate("2025-09-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2025-10-01"), bt.decimalDate("2025-11-01"), color='black', alpha=0.05)

    # Legend code adapted from https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html
    if legend:
        custom_lines = []
        custom_labels = []
        for cat, color in state_2_color.items():
            if cat == '-':
                continue    # Suppress extremely rare but confusing category in legend
            custom_lines.append(Line2D([0], [0], color=color, lw=3))
            custom_labels.append(cat)
        ax.legend(custom_lines, custom_labels, bbox_to_anchor=(0.85,0.98), framealpha=0.5)
    
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

    #annotate_node(mrca_node_of(mcc_tree, ['ON563414']), 'B.1')
    #annotate_node(mrca_node_of(mcc_tree, ['ON676708', 'ON563414']), 'A.1.1')
    #annotate_node(mrca_node_of(mcc_tree, ['ON676708', 'OP612681']), 'A.1')
    #annotate_node(mrca_node_of(mcc_tree, ['ON674051', 'ON676707']), 'A.2')

    if tMRCAs:
        annotate_node(mcc_tree.root,
                      '',   #'tMRCA',
                      diamond=True, x=mean_tMRCA, rotation='vertical', verticalalignment='top')
        
        ax.plot([mean_tMRCA, mean_tMRCA], [mcc_tree.root.y, -20], color='black', linestyle='--', lw=0.75)
    
    plt.savefig(out_pdf_filename, bbox_inches='tight', pad_inches=0)
    print(f'Plot saved to {out_pdf_filename}')

Path('plots').mkdir(parents=True, exist_ok=True)

prefix = 'h5n1-andersen-e756a15-ALL_full_dates_only'
data = [
    # data_dir                   mcc_filename             log_filename            plot_filename
    ('delphy_outputs_a',       f'{prefix}_a.mcc',       f'{prefix}_a.log',       'DelphyAMcc.pdf'),
    ('delphy_outputs_b',       f'{prefix}_b.mcc',       f'{prefix}_b.log',       'DelphyBMcc.pdf'),
    ('beastX_run_a',           f'{prefix}_beastX.mcc',   'output.log',           'BeastXAMcc.pdf'),
    ('beastX_run_b',           f'{prefix}_beastX.mcc',   'output.log',           'BeastXBMcc.pdf'),
    ('delphy_outputs_alpha_a', f'{prefix}_alpha_a.mcc', f'{prefix}_alpha_a.log', 'DelphyAAlphaMcc.pdf'),
    ('delphy_outputs_alpha_b', f'{prefix}_alpha_b.mcc', f'{prefix}_alpha_b.log', 'DelphyBAlphaMcc.pdf'),
    ('beastX_run_alpha_a',     f'{prefix}_alpha_beastX.mcc',   'output.log',           'BeastXAAlphaMcc.pdf'),
    ('beastX_run_alpha_b',     f'{prefix}_alpha_beastX.mcc',   'output.log',           'BeastXBAlphaMcc.pdf'),
]

for (data_dir, mcc_filename, log_filename, plot_filename) in data:
    mcc_tree = bt.loadNexus(f'{data_dir}/{mcc_filename}')
    leaves = mcc_tree.traverse_tree()

    # Both Delphy nor BEAST absolutely respect the stated date of every leaf, so there's no need
    # to account for "outliers" as in the ML tree
    max_date = max(leaf.name.split('|')[-1] for leaf in leaves if leaf.name not in BAD_LEAVES)  # e.g., 2015-10-24
    
    tMRCAs = extract_tMRCA_distr(f'{data_dir}/{log_filename}', burnin=0.30, t0=bt.decimalDate(max_date))
    ess_table = read_in_ess_table(f'{data_dir}/log_analysis.txt', None,
                                  delphy_log_format if data_dir.startswith('delphy') else beastX_log_format)
    plot_tree(mcc_tree, f'plots/{plot_filename}', tMRCAs, ess_table['TreeHeight'], legend=data_dir.startswith('delphy'))

plot_tree(bt.loadNexus('./ml/tt/timetree.nexus'), 'plots/ML.pdf', None, None, legend=False)
plot_tree(bt.loadNexus('./ml/tt_alpha/timetree.nexus'), 'plots/MLAlpha.pdf', None, None, legend=False)
