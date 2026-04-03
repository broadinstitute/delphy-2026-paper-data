#!/usr/bin/env python

import baltic as bt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

# Read in sample IDs & clades
# ===========================
# sample_ids.txt = Sample ids for 772 sequences used in Fig 3A of LeMieux et al (2021) (private communication)
print("\nReading in sample IDs...")
sample_id_2_internal_clade = {}
with open('sample_ids.csv', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue   # Skip comments
        if line.startswith('id,'):
            continue   # Skip header
        stripped_line = line.strip()
        if not stripped_line:
            continue   # Skip empty lines
        
        id, clade = stripped_line.split(',')
        sample_id_2_internal_clade[id] = clade

sample_ids = list(sample_id_2_internal_clade.keys())
print(f'Read {len(sample_ids)} samples')
print(f'First 5 IDs: {sample_ids[:5]}')
print(f'Last 5 IDs: {sample_ids[-5:]}')

# Colors from https://github.com/JacobLemieux/sarscov2pub/blob/ec90f30ca5c115f8dfb2abfab3604723ff8a521e/scripts/main_figures.R#L204
cbPalette = [
    "#E69F00",   # 0 == orange
    "#56B4E9",   # 1 == light blue
    "#009E73",   # 2 == dark green
    "#F0E442",   # 3 == yellow
    "#0072B2",   # 4 == dark blue
    "#D55E00",   # 5 == dark orange
    "#CC79A7",   # 6 == violet
    "#999999",   # 7 == grey
]

# Simplified geo labels
# We use shorter labels than in the original paper here

# internal_2_external_clade = {
#     'Clade_1': "C20099T (BHCHP)",
#     'Clade_2': "G3892T (SNF)",
#     'Clade_3': "C2416T (Conference, BHCHP)",
#     'Clade_4': "G105T (BHCHP)",
#     'Clade_5': "G28899T",
# }

# external_clade_2_color = {
#     'C2416T (Conference, BHCHP)': cbPalette[2],  # dark green
#     'G105T (BHCHP)':              cbPalette[4],  # dark blue
#     'G3892T (SNF)':               cbPalette[3],  # yellow
#     'C20099T (BHCHP)':            cbPalette[0],  # orange
#     'G28899T':                    cbPalette[1],  # light blue
#     'Other':                      'black',
# }

internal_2_external_clade = {
    'Clade_1': "C20099T",
    'Clade_2': "G3892T",
    'Clade_3': "C2416T",
    'Clade_4': "G105T",
    'Clade_5': "G28899T",
}

external_clade_2_color = {
    'C2416T':  cbPalette[2],  # dark green
    'G105T':   cbPalette[4],  # dark blue
    'G3892T':  cbPalette[3],  # yellow
    'C20099T': cbPalette[0],  # orange
    'G28899T': cbPalette[1],  # light blue
    'Other':   'black',
}

def external_clade_of(leaf_name):
    sample_id, *rest = leaf_name.split('|')
    geo = sample_id_2_internal_clade.get(sample_id, "?")
    return internal_2_external_clade.get(geo, 'Other')

def color_of(leaf_name):
    return external_clade_2_color[external_clade_of(leaf_name)]

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


def plot_tree(mcc_filename, out_pdf_filename, legend=False):
    mcc_tree = bt.loadNexus(mcc_filename)
    mcc_tree.traverse_tree()
    mcc_tree.setAbsoluteTime(bt.decimalDate("2020-05-09"))  # Hard-coded, but ok

    index_2_external_clade = simple_parsimony(mcc_tree, lambda node: external_clade_of(node.name))
    
    fig,ax = plt.subplots(figsize=(4,5),facecolor='w')

    x_attr=lambda k: k.absoluteTime
    c_func=lambda k: 'black' if k.branchType != 'leaf' else color_of(k.name)
    s_func=lambda k: 20 if (k.branchType == 'leaf' and external_clade_of(k.name) != "Other") else 1
    
    mcc_tree.sortBranches(descending=True)
    for k in mcc_tree.getInternal():
        k.children.reverse()
    mcc_tree.drawTree()
    
    mcc_tree.plotTree(ax,
                      x_attr=x_attr,
                      colour=lambda node: external_clade_2_color[index_2_external_clade[node.index]],
                      width=0.3)
    mcc_tree.plotPoints(
        ax,
        x_attr=x_attr,
        size=s_func,
        colour=c_func,
        zorder=100,
        outline=False)

    target_func=lambda k: k.branchType=='node' and 'posterior' in k.traits and k.traits['posterior'] >= 0.95 ## only target high-posterior-support nodes

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

    
    ax.set_xlim(2019.9, 2020.37);
    ax.set_ylim(-5, mcc_tree.ySpan+75);
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);
    ax.get_xaxis().set_ticks([
        bt.decimalDate("2020-01-01"),
        bt.decimalDate("2020-02-06"),
        bt.decimalDate("2020-03-14"),
        bt.decimalDate("2020-04-19"),
       ])
    ax.set_xticklabels([
        "1 Jan\n2020",
        "6 Feb\n2020",
        "14 March\n2020",
        "19 April\n2020",
       ]);
    
    plt.axvspan(bt.decimalDate("2020-01-01"), bt.decimalDate("2020-02-01"),
                color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2020-03-01"), bt.decimalDate("2020-04-01"),
                color='black', alpha=0.05)

    if legend:
        # Legend code adapted from https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html
        custom_lines = []
        custom_labels = []
        for external_clade, color in external_clade_2_color.items():
            custom_lines.append(Line2D([0], [0], color=color, lw=3))
            custom_labels.append(external_clade)
        ax.legend(custom_lines, custom_labels, bbox_to_anchor=(0.35, 0.85), framealpha=0.5)
    
    plt.savefig(out_pdf_filename, bbox_inches='tight', pad_inches=0)
    print(f'Plot saved to {out_pdf_filename}')

Path('plots').mkdir(parents=True, exist_ok=True)

plot_tree('./delphy_outputs_a/ma_sars_cov_2_delphy.mcc', 'plots/DelphyAMcc.pdf', legend=True)
plot_tree('./delphy_outputs_b/ma_sars_cov_2_delphy.mcc', 'plots/DelphyBMcc.pdf')
plot_tree('./delphy_outputs_alpha_a/ma_sars_cov_2_delphy_alpha.mcc', 'plots/DelphyAMccAlpha.pdf', legend=True)
plot_tree('./delphy_outputs_alpha_b/ma_sars_cov_2_delphy_alpha.mcc', 'plots/DelphyBMccAlpha.pdf')
plot_tree('./beastX_run/ma_sars_cov_2_beastX.mcc', 'plots/BeastXMcc.pdf')
plot_tree('./beastX_run_alpha/ma_sars_cov_2_alpha_beastX.mcc', 'plots/BeastXMccAlpha.pdf')
plot_tree('./beast2_run/ma_sars_cov_2_beast2.mcc', 'plots/Beast2Mcc.pdf')
plot_tree('./beast2_run_alpha/ma_sars_cov_2_alpha_beast2.mcc', 'plots/Beast2MccAlpha.pdf')
plot_tree('./ml/tt/timetree.nexus', 'plots/ML.pdf')
plot_tree('./ml/tt_alpha/timetree.nexus', 'plots/MLAlpha.pdf')
