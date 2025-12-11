#!/usr/bin/env python

import math
import baltic as bt
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess

def plots_for(sim, rep):  # e.g., sim == "exp_100", rep == "a"
    
    # Baltic doesn't like inner node names (!)
    orig_tree_filename = f'{sim}/ground_truth/{sim}_tree.nwk'
    sanitized_tree_filename = f'{sim}/ground_truth/{sim}_tree.sanitized_nwk'
    mcc_tree_filename = f'{sim}/delphy_outputs_{rep}/{sim}.mcc'
    ml_tree_filename = f'{sim}/ml_outputs/timetree.nexus'
    with open(orig_tree_filename, mode='rt', encoding='utf-8') as fin, open(sanitized_tree_filename, mode='wt', encoding='utf-8') as fout:
        subprocess.run([
            'sed',
            's/NODE_[^:]*://g'
        ], stdin=fin, stdout=fout)

    real_tree = date_tree(bt.loadNewick(sanitized_tree_filename))
    mcc_tree = date_tree(bt.loadNexus(mcc_tree_filename))
    min_t = min([real_tree.root.absoluteTime,
                 mcc_tree.root.absoluteTime])

    if Path(ml_tree_filename).exists():
        ml_tree = date_tree(bt.loadNexus(ml_tree_filename))
        min_t = min([min_t,
                     ml_tree.root.absoluteTime])
    else:
        ml_tree = None
        
    tip_colors = plot_tree(real_tree, min_t, {}, f'{sim}/ground_truth/{sim}_tree.pdf')
    plot_tree(mcc_tree, min_t, tip_colors, f'{sim}/delphy_outputs_{rep}/{sim}_{rep}_mcc.pdf')
    if ml_tree is not None:
        plot_tree(ml_tree, min_t, tip_colors, f'{sim}/ml_outputs/{sim}_ml.pdf')

def date_tree(tree):
    tree.traverse_tree();
    tree.setAbsoluteTime(bt.decimalDate('2024-07-31'))  # This is a slight approximation
    return tree
    
def plot_tree(tree, min_t, tip_colors, out_pdf_filename):
    # If tip_colors is falsy, we'll create a dictionary of tip names to tip colors
    # that ranges from green to blue in the vertical order of the tips (when the tree is ladderized)
    # and return it

    fig,ax = plt.subplots(figsize=(5,7),facecolor='w')

    if not tip_colors:
        total_tips = len(tree.getExternal())
        i = 0
        tip_colors = {}
        for tip in tree.getExternal():
            c = i / total_tips
            tip_colors[tip.name] = [0, c, (1-c)]
            i = i+1
    
    x_attr=lambda k: k.absoluteTime
    c_func=lambda k: 'black' if k.branchType != 'leaf' else tip_colors[k.name]
    s_func=lambda k: 20 if (k.branchType == 'leaf') else 1
    
    tree.sortBranches(descending=True)
    for k in tree.getInternal():
        k.children.reverse()
    tree.drawTree()
        
    tree.plotTree(ax,x_attr=x_attr,colour='#5B5B5C',width=0.5)
    tree.plotPoints(
        ax,
        x_attr=x_attr,
        size=s_func,
        colour=c_func,
        zorder=100,
        outline=False)
    
    ax.set_ylim(-5, tree.ySpan+5);
    ax.set_xlim(min(min_t, bt.decimalDate('2023-12-01')), bt.decimalDate('2024-08-01'));
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);
    if min_t >= bt.decimalDate('2023-09-01'):
        ax.get_xaxis().set_ticks([
            bt.decimalDate("2023-12-01"),
            bt.decimalDate("2024-02-01"),
            bt.decimalDate("2024-04-01"),
            bt.decimalDate("2024-06-01"),
            bt.decimalDate('2024-08-01')
        ])
        ax.set_xticklabels([
            "1 Dec\n2023",
            "1 Feb\n2024",
            "1 Apr\n2024",
            "1 Jun\n2024",
            "1 Aug\n2024",
        ])
    else:
        ticks = []
        labels = []
        year = int(math.floor(min_t))
        while year <= 2024:
            ticks.append(year)
            labels.append(f"1 Jan\n{year}")
            year += 1
        
        ax.get_xaxis().set_ticks(ticks)
        ax.set_xticklabels(labels)
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')

    return tip_colors

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
        plots_for(sim, rep)

for coal_bins in [625, 1250, 2500, 5000]:
   plots_for('exp_100000', f'c_{coal_bins}bins')
