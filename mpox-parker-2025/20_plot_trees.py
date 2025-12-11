#!/usr/bin/env python
#
# These plots are adapted from (from https://github.com/andersen-lab/Mpox_West_Africa/blob/main/Scripts/Figure_3/Tree_DTA.ipynb)

import baltic as bt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import csv
from scipy import stats
import pandas as pd
import math
import numpy as np
import random

# Read in metadata
# ===========================
print("\nReading in metadata...")
sample_id_2_region = {}
sample_id_2_state = {}
state_2_region = {}
with open('delphy_inputs/mpox-parker-2025_metadata.csv', 'r') as f:
    r = csv.DictReader(f)
    for record in r:
        sample_id_2_state[record['id']] = record['State']
        sample_id_2_region[record['id']] = record['Region']

        if record['State'] in state_2_region:
            assert state_2_region[record['State']] == record['Region']
        else:
            state_2_region[record['State']] = record['Region']

print(sorted(list(state_2_region.items())));
            
# Mapping of state to region and pseudo-region
def region_of(leaf_name):
    sample_id, *rest = leaf_name.split('|')
    region = sample_id_2_region.get(sample_id, "Rest")
    state = sample_id_2_state.get(sample_id, "Rest")
    return region
        
def state_of(leaf_name):
    sample_id, *rest = leaf_name.split('|')
    region = sample_id_2_region.get(sample_id, "Rest")
    state = sample_id_2_state.get(sample_id, "Rest")
    return state
        
def merge_states_to_pseudo_regions(state):
    if state == 'Rivers':
        return 'Rivers'
    else:
        return state_2_region[state]
        
# Pseudo-region color map
region_2_color = {
    'SS':      '#16166B',
    'SW':      '#138808',
    'SE':      '#0070BB',
    'NW':      '#FFD700',
    'NE':      '#FE5A1D',
    'NC':      '#960018',
    'Rest':    '#C0C0C0',
    'Nigeria': 'black',
}
pseudo_region_2_color = {
    'Rivers': '#273BE2',
} | region_2_color
#del pseudo_region_2_color['Rest']
#del pseudo_region_2_color['Nigeria']

def color_of(leaf_name, color_palette=region_2_color, mapper=region_of):
    return color_palette[mapper(leaf_name)]

def mrca_node_of(tree, leaf_names):
    leaf_nodes = [node for node in tree.Objects if node.is_leaf() and node.name.split('|')[0] in leaf_names]
    if len(leaf_nodes) == 1:
        return leaf_nodes[0]
    else:
        return tree.commonAncestor(leaf_nodes)

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
            if len(node_fs) == 1:
                node_state = list(node_fs)[0]
            else:
                node_state = random.choice(list(node_fs))
                print(f'WARNING: Arbitrarily breaking tie between {node_fs} in favor of {node_state}')
        index_2_state[node.index] = node_state

        if not node.is_leaf():
            for child in node.children:
                gogo(child, node_state)
    gogo(tree.root, None)

    return index_2_state


def extract_tMRCA_distr(log_filename, burnin):
    # Extremely brittle way of extracting a tMRCA distribution from the log files
    raw_data = pd.read_table(log_filename, comment='#')
    num_pts = len(raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    data = raw_data[burnin_pts:]
    print(f'Read in log file {log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(data)} usable')

    return list(data['age(root)'])


def plot_tree(mcc_filename, out_pdf_filename, tMRCAs, tMRCAs_BEAST, color_palette=region_2_color, mapper=region_of, cat_merger=lambda x: x):
    mcc_tree = bt.loadNexus(mcc_filename)
    leaves = mcc_tree.traverse_tree()
    max_date = max(leaf.name.split('|')[-1] for leaf in leaves)  # e.g., 2020-03-01
    mcc_tree.setAbsoluteTime(bt.decimalDate(max_date))

    index_2_cat = simple_parsimony(mcc_tree, lambda node: mapper(node.name))
    for index, cat in list(index_2_cat.items()):
        index_2_cat[index] = cat_merger(cat)
    
    fig,ax = plt.subplots(figsize=(6,6),facecolor='w')

    x_attr=lambda k: k.absoluteTime
    c_func=lambda node: color_palette[index_2_cat[node.index]]
    s_func=lambda k: 15 if (k.branchType == 'leaf') else 1
    
    mcc_tree.sortBranches(descending=True)
    #for k in mcc_tree.getInternal():  # Match ordering in Parker et al 2025
    #    k.children.reverse()
    mcc_tree.drawTree()

    mcc_tree.plotTree(ax,
                      x_attr=x_attr,
                      colour=c_func,
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
    
    mean_tMRCA = np.mean(tMRCAs)  # Slightly different from root height because tMRCA was sampled much more often than trees
    kde = stats.gaussian_kde(tMRCAs)
    tt = np.linspace(2014, 2017.0, 200)
    ax.plot(tt, -5+60*kde(tt), color='#888888', lw=0.5)
    ax.fill_between(tt, -5+60*kde(tt), -5, color='green', alpha=0.2)

    mean_tMRCA_BEAST = np.mean(tMRCAs_BEAST)  # Slightly different from root height because tMRCA was sampled much more often than trees
    kde = stats.gaussian_kde(tMRCAs_BEAST)
    tt = np.linspace(2014, 2017.0, 200)
    ax.plot(tt, -5+60*kde(tt), color='#888888', lw=0.5, ls='dashed')
    ax.fill_between(tt, -5+60*kde(tt), -5, color='blue', alpha=0.1)
    
    ax.set_xlim(2014.7, 2024.3)
    
    ax.set_ylim(-5, mcc_tree.ySpan+1);
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);

    ax.get_xaxis().set_ticks([
        bt.decimalDate("2015-01-01"),
        bt.decimalDate("2016-01-01"),
        bt.decimalDate("2017-01-01"),
        bt.decimalDate("2018-01-01"),
        bt.decimalDate("2019-01-01"),
        bt.decimalDate("2020-01-01"),
        bt.decimalDate("2021-01-01"),
        bt.decimalDate("2022-01-01"),
        bt.decimalDate("2023-01-01"),
        bt.decimalDate("2024-01-01"),
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
        "2023",
        "2024",
    ]);
    #ax.set_xlabel('Time (years)')

    plt.axvspan(2015.0, 2016.0, color='black', alpha=0.05)
    plt.axvspan(2017.0, 2018.0, color='black', alpha=0.05)
    plt.axvspan(2019.0, 2020.0, color='black', alpha=0.05)
    plt.axvspan(2021.0, 2022.0, color='black', alpha=0.05)
    plt.axvspan(2023.0, 2024.0, color='black', alpha=0.05)

    # Legend code adapted from https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html
    custom_lines = []
    custom_labels = []
    for cat, color in color_palette.items():
        custom_lines.append(Line2D([0], [0], color=color, lw=3))
        custom_labels.append(cat)
    ax.legend(custom_lines, custom_labels, bbox_to_anchor=(0.25,0.40), framealpha=0.5, prop={'size': 9})
    
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

    annotate_node(mcc_tree.root, 'tMRCA hMPXV-1', diamond=True, x=mean_tMRCA, rotation='vertical', verticalalignment='top')
    
    ax.plot([mean_tMRCA, mean_tMRCA], [mcc_tree.root.y, -5], color='black', linestyle='--', lw=0.75)
    
    plt.savefig(out_pdf_filename, bbox_inches='tight', pad_inches=0)
    print(f'Plot saved to {out_pdf_filename}')


Path('plots').mkdir(parents=True, exist_ok=True)

tMRCAs_BEAST = extract_tMRCA_distr(f'./beastX_run/Mpox_2poch_combined.log', burnin=0.30)
for run in ('a', 'b'):
    tMRCAs = extract_tMRCA_distr(f'./delphy_outputs_{run}/mpox-parker-2025.log', burnin=0.30)
    plot_tree(f'./delphy_outputs_{run}/mpox-parker-2025.mcc', f'plots/{run.upper()}DelphyMcc-ByRegion.pdf', tMRCAs, tMRCAs_BEAST,
              mapper=region_of, color_palette=region_2_color)
    plot_tree(f'./delphy_outputs_{run}/mpox-parker-2025.mcc', f'plots/{run.upper()}DelphyMcc-ByPseudoRegion.pdf', tMRCAs, tMRCAs_BEAST,
              mapper=state_of, cat_merger=merge_states_to_pseudo_regions, color_palette=pseudo_region_2_color)

def shorten_id(long_id):
    if long_id.split('|')[0] == 'unpub':   # e.g., unpub|TRM076|Nigeria|Rivers|2022-09-27
        return long_id.split('|')[1]
    elif long_id.startswith('PP'):         # e.g., PP852976|MPXV|TRM081|Nigeria|Bayelsa|2022-10-21
        return long_id.split('|')[2]
    else:
        return long_id.split('|')[0]  # e.g., OP535319|MPXV|Nigeria|Cross-River|2017-12

plot_tree(f'./beastX_run/mpox-parker-2025-beastX.mcc', f'plots/BeastMcc-ByRegion.pdf', tMRCAs, tMRCAs_BEAST,
          mapper=lambda n : region_of(shorten_id(n)), color_palette=region_2_color)
plot_tree(f'./beastX_run/mpox-parker-2025-beastX.mcc', f'plots/BeastMcc-ByPseudoRegion.pdf', tMRCAs, tMRCAs_BEAST,
          mapper=lambda n : state_of(shorten_id(n)), cat_merger=merge_states_to_pseudo_regions, color_palette=pseudo_region_2_color)
