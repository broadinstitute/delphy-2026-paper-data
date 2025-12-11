#!/usr/bin/env python

import baltic as bt
import matplotlib.pyplot as plt
from pathlib import Path

def plot_tree(mcc_tree, out_pdf_filename, legend=False):
    leaves = mcc_tree.traverse_tree()
    max_date = max(leaf.name.split('|')[-1] for leaf in leaves)
    mcc_tree.setAbsoluteTime(bt.decimalDate(max_date))

    fig,ax = plt.subplots(figsize=(4,5),facecolor='w')

    x_attr=lambda k: k.absoluteTime
    c_func=lambda k: 'black'
    s_func=lambda k: 20 if (k.branchType == 'leaf') else 1
    
    mcc_tree.sortBranches(descending=True)
    for k in mcc_tree.getInternal():
        k.children.reverse()
    mcc_tree.drawTree()
    
    mcc_tree.plotTree(ax,
                      x_attr=x_attr,
                      colour='#5B5B5C',
                      width=0.5)
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

    ax.set_xlim(bt.decimalDate("2014-03-01"), bt.decimalDate("2014-07-05"));
    ax.set_ylim(-5, mcc_tree.ySpan+5);
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);
    ax.get_xaxis().set_ticks([
        bt.decimalDate("2014-03-01"),
        bt.decimalDate("2014-04-01"),
        bt.decimalDate("2014-05-01"),
        bt.decimalDate("2014-06-01"),
        bt.decimalDate("2014-07-01"),
       ])
    ax.set_xticklabels([
        "1 Mar\n2014",
        "1 Apr\n2014",
        "1 May\n2014",
        "1 Jun\n2014",
        "1 Jul\n2014",
       ]);
    
    plt.axvspan(bt.decimalDate("2014-03-01"), bt.decimalDate("2014-04-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2014-05-01"), bt.decimalDate("2014-06-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2014-07-01"), bt.decimalDate("2014-08-01"), color='black', alpha=0.05)

    plt.savefig(out_pdf_filename, bbox_inches='tight', pad_inches=0)
    print(f'Plot saved to {out_pdf_filename}')

Path('plots').mkdir(parents=True, exist_ok=True)

plot_tree(bt.loadNexus('./delphy_outputs_a/ebola_delphy.mcc'), 'plots/DelphyAMcc.pdf', legend=True)
plot_tree(bt.loadNexus('./delphy_outputs_alpha_a/ebola_delphy_alpha.mcc'), 'plots/DelphyAMccAlpha.pdf', legend=True)
plot_tree(bt.loadNexus('./delphy_outputs_b/ebola_delphy.mcc'), 'plots/DelphyBMcc.pdf', legend=True)
plot_tree(bt.loadNexus('./delphy_outputs_alpha_b/ebola_delphy_alpha.mcc'), 'plots/DelphyBMccAlpha.pdf', legend=True)
plot_tree(bt.loadNexus('./beast2_run/ebola_beast2.mcc'), 'plots/Beast2Mcc.pdf', legend=False)
plot_tree(bt.loadNexus('./beast2_run_alpha/ebola_alpha_beast2.mcc'), 'plots/Beast2MccAlpha.pdf', legend=False)
plot_tree(bt.loadNexus('./beastX_run/ebola_beastX.mcc'), 'plots/BeastXMcc.pdf', legend=False)
plot_tree(bt.loadNexus('./beastX_run_alpha/ebola_alpha_beastX.mcc'), 'plots/BeastXMccAlpha.pdf', legend=False)
plot_tree(bt.loadNexus('./ml/tt/timetree.nexus'), 'plots/ML.pdf', legend=False)
plot_tree(bt.loadNexus('./ml/tt_alpha/timetree.nexus'), 'plots/MLAlpha.pdf', legend=False)
