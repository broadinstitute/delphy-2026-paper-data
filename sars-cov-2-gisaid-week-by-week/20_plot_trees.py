#!/usr/bin/env python

import csv
import baltic as bt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import subprocess

# Conversions for the humans
# ==========================
epiweek_2_hyphenated = {
    '202001': '2020-01',
    '202002': '2020-02',
    '202003': '2020-03',
    '202004': '2020-04',
    '202005': '2020-05',
    '202006': '2020-06',
    '202007': '2020-07',
    '202008': '2020-08',
    '202009': '2020-09',
    '202010': '2020-10',
    '202011': '2020-11',
    '202012': '2020-12',
    '202013': '2020-13',
}
epiweek_2_end_date = {  # See https://www.cmmcp.org/sites/g/files/vyhlif2966/f/uploads/epiweekcalendar2020.pdf
    '202001': '4 Jan 2020',
    '202002': '11 Jan 2020',
    '202003': '18 Jan 2020',
    '202004': '25 Jan 2020',
    '202005': '1 Feb 2020',
    '202006': '8 Feb 2020',
    '202007': '15 Feb 2020',
    '202008': '22 Feb 2020',
    '202009': '29 Feb 2020',
    '202010': '7 Mar 2020',
    '202011': '14 Mar 2020',
    '202012': '21 Mar 2020',
    '202013': '28 Mar 2020',
}

# Read in metadata
# ===========================
print("\nReading in metadata...")
def sanitizeTabs(s):
    return s.replace("\t", " ")  # Avoid CSV pain

sample_id_2_geo0 = {}
with open('metadata_20200331.tsv', 'r') as f:
    ff = csv.reader(f, delimiter='\t')
    header = None
    for row in ff:
        if not header:
            header = row
            continue
        
        virusName, accessionId, collectionDate, location, additionalLocationInfo, seqLength, host, submissionDate = row

        sample_id_2_geo0[f'{accessionId}'] = \
            sanitizeTabs(location).split(' / ')[0]

geo0_to_simple_geo0 = {
    'Asia':          'Asia',
    'Europe':        'Europe',
    'North America': 'North America',
    'South America': 'Other',
    'Oceania':       'Other',
    'Africa':        'Other',
    'Other':         'Other',
}
        
simple_geo0_2_color = {
    'Asia':          '#cccccc',
    'Europe':        'red',
    'North America': 'black',
    'Other':         'orange',
    #'Asia':          '#E4E1DC',#'#cccccc',
    #'Europe':        '#00CC66',#'red',
    #'North America': '#DC0073',#'orange',
    #'Other':         '#192A51',#'black',
}

def simple_geo0_of(leaf_name):
    sample_id, *rest = leaf_name.split('|')
    geo0 = sample_id_2_geo0.get(sample_id, "Other")
    simple_geo_0 = geo0_to_simple_geo0.get(geo0, "Other")
    return simple_geo_0

def color_of(leaf_name):
    return simple_geo0_2_color[simple_geo0_of(leaf_name)]

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


def plot_tree(mcc_filename, out_pdf_filename, title):
    mcc_tree = bt.loadNexus(mcc_filename)
    leaves = mcc_tree.traverse_tree()
    max_date = max(leaf.name.split('|')[-1] for leaf in leaves)  # e.g., 2020-03-01
    mcc_tree.setAbsoluteTime(bt.decimalDate(max_date))

    index_2_simple_geo0 = simple_parsimony(mcc_tree, lambda node: simple_geo0_of(node.name))
    
    fig,ax = plt.subplots(figsize=(5,7),facecolor='w')

    x_attr=lambda k: k.absoluteTime
    c_func=lambda k: 'black' if k.branchType != 'leaf' else color_of(k.name)
    s_func=lambda k: 1 if k.branchType == 'leaf' else 0.5
    w_func=lambda k: 0.1 if index_2_simple_geo0[k.index] == 'Asia' else 0.5
    
    mcc_tree.sortBranches(descending=True)
    for k in mcc_tree.getInternal():
        k.children.reverse()
    mcc_tree.drawTree()

    ymin, ymax = -0.05*mcc_tree.ySpan, 1.01*mcc_tree.ySpan
    for (weeknum,d) in [
            (47, "2019-11-17"),  # start of CDC epi week 2019-47
            (48, "2019-11-24"),  # start of CDC epi week 2019-48
            (49, "2019-12-01"),  # start of CDC epi week 2019-49
            (50, "2019-12-08"),  # start of CDC epi week 2019-50
            (51, "2019-12-15"),  # start of CDC epi week 2019-51
            (52, "2019-12-22"),  # start of CDC epi week 2019-52
            ( 1, "2019-12-29"),  # start of CDC epi week 2020-01
            ( 2, "2020-01-05"),  # start of CDC epi week 2020-02
            ( 3, "2020-01-12"),  # start of CDC epi week 2020-03
            ( 4, "2020-01-19"),  # start of CDC epi week 2020-04
            ( 5, "2020-01-26"),  # start of CDC epi week 2020-05
            ( 6, "2020-02-02"),  # start of CDC epi week 2020-06
            ( 7, "2020-02-09"),  # start of CDC epi week 2020-07
            ( 8, "2020-02-16"),  # start of CDC epi week 2020-08
            ( 9, "2020-02-23"),  # start of CDC epi week 2020-09
            (10, "2020-03-01"),  # start of CDC epi week 2020-10
            (11, "2020-03-08"),  # start of CDC epi week 2020-11
            (12, "2020-03-15"),  # start of CDC epi week 2020-12
            (13, "2020-03-22"),  # start of CDC epi week 2020-13
            (14, "2020-03-29"),  # start of CDC epi week 2020-14
    ]:
        plt.vlines(bt.decimalDate(d), ymin, ymax, linestyle='--', color='#dddddd', lw=0.25)
        if weeknum != 14:
            plt.text(bt.decimalDate(d) + 3.5/365, ymin, str(weeknum), horizontalalignment='center', verticalalignment='bottom',
                     color='lightgrey', fontsize=7)
        
    mcc_tree.plotTree(ax,
                      x_attr=x_attr,
                      colour=lambda node: simple_geo0_2_color[index_2_simple_geo0[node.index]],
                      width=w_func)
    mcc_tree.plotPoints(
        ax,
        x_attr=x_attr,
        size=s_func,
        colour=c_func,
        zorder=100,
        outline=False)
    
    ax.set_xlim(bt.decimalDate("2019-11-17"),   # Start of CDC epi week 2019-47
                bt.decimalDate("2020-03-29"));  # Start of CDC epi week 2020-14
    ax.set_ylim(ymin, ymax)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);
    ax.get_xaxis().set_ticks([
        bt.decimalDate("2019-12-01"),  # start of CDC epi week 2019-49
        #bt.decimalDate("2019-12-08"),  # start of CDC epi week 2019-50
        #bt.decimalDate("2019-12-15"),  # start of CDC epi week 2019-51
        #bt.decimalDate("2019-12-22"),  # start of CDC epi week 2019-51
        #bt.decimalDate("2019-12-29"),  # start of CDC epi week 2020-01
        bt.decimalDate("2020-01-05"),  # start of CDC epi week 2020-02
        #bt.decimalDate("2020-01-12"),  # start of CDC epi week 2020-03
        #bt.decimalDate("2020-01-19"),  # start of CDC epi week 2020-04
        #bt.decimalDate("2020-01-26"),  # start of CDC epi week 2020-05
        bt.decimalDate("2020-02-02"),  # start of CDC epi week 2020-06
        #bt.decimalDate("2020-02-09"),  # start of CDC epi week 2020-07
        #bt.decimalDate("2020-02-16"),  # start of CDC epi week 2020-08
        #bt.decimalDate("2020-02-23"),  # start of CDC epi week 2020-09
        bt.decimalDate("2020-03-01"),  # start of CDC epi week 2020-10
        #bt.decimalDate("2020-03-08"),  # start of CDC epi week 2020-11
        #bt.decimalDate("2020-03-15"),  # start of CDC epi week 2020-12
        #bt.decimalDate("2020-03-22"),  # start of CDC epi week 2020-13
        #bt.decimalDate("2020-03-29"),  # start of CDC epi week 2020-14
       ])
    ax.set_xticklabels([
        "1 Dec\n2019",
        "5 Jan\n2020",
        #"12 Jan\n2020",
        #"19 Jan\n2020",
        #"26 Jan\n2020",
        "2 Feb\n2020",
        #"9 Feb\n2020",
        #"16 Feb\n2020",
        #"23 Feb\n2020",
        "1 Mar\n2020",
        #"8 Mar\n2020",
        #"15 Mar\n2020",
        #"22 Mar\n2020",
        #"29 Mar\n2020",
       ]);

    # plt.axvspan(bt.decimalDate("2019-11-24"), bt.decimalDate("2019-12-01"),
    #             color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2019-12-08"), bt.decimalDate("2019-12-15"),
    #             color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2019-12-22"), bt.decimalDate("2019-12-29"),
    #             color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2020-01-05"), bt.decimalDate("2020-01-12"),
    #             color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2020-01-19"), bt.decimalDate("2020-01-26"),
    #             color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2020-02-02"), bt.decimalDate("2020-02-09"),
    #             color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2020-02-16"), bt.decimalDate("2020-02-23"),
    #             color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2020-03-01"), bt.decimalDate("2020-03-08"),
    #             color='black', alpha=0.05)
    # plt.axvspan(bt.decimalDate("2020-03-15"), bt.decimalDate("2020-03-22"),
    #             color='black', alpha=0.05)

    # Legend code adapted from https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html
    custom_lines = []
    custom_labels = []
    for simple_geo0, color in simple_geo0_2_color.items():
        custom_lines.append(Line2D([0], [0], color=color, lw=3))
        custom_labels.append(simple_geo0)
    ax.legend(custom_lines, custom_labels, bbox_to_anchor=(0.3, 0.98), framealpha=0.5)

    plt.title(title.replace('NNN', str(len(leaves))))
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')


def compound_plot(epiweek_2_tree_file, out_pdf_filename):
    
    fig,ax = plt.subplots(3, 4, figsize=(11,8.5), facecolor='w')
    plt.tight_layout()
    plt.subplots_adjust(top=0.95, bottom=0.05, hspace=0.2, wspace=0.3)

    for i in range(12):
        if i >= len(epiweek_2_tree_file):
            this_ax = plt.subplot(3,4,i+1)
            fig.delaxes(this_ax)
    
    for i, (epiweek, mcc_filename) in enumerate(epiweek_2_tree_file.items()):
        print(f'Reading {mcc_filename}')
        mcc_tree = bt.loadNexus(mcc_filename)
        leaves = mcc_tree.traverse_tree()
        max_date = max(leaf.name.split('|')[-1] for leaf in leaves)  # e.g., 2020-03-01
        mcc_tree.setAbsoluteTime(bt.decimalDate(max_date))
    
        index_2_simple_geo0 = simple_parsimony(mcc_tree, lambda node: simple_geo0_of(node.name))

        ax = plt.subplot(3,4,i+1)
        plt.title(f'Epiweek {epiweek_2_hyphenated[epiweek]} (N={len(leaves)})\n(ending on {epiweek_2_end_date[epiweek]})')
        plt.xlabel(None)
        
        x_attr=lambda k: k.absoluteTime
        c_func=lambda k: 'black' if k.branchType != 'leaf' else color_of(k.name)
        s_func=lambda k: 1 if k.branchType == 'leaf' else 0.5
        w_func=lambda k: 0.1 if index_2_simple_geo0[k.index] == 'Asia' else 0.5
        
        mcc_tree.sortBranches(descending=True)
        for k in mcc_tree.getInternal():
            k.children.reverse()
            mcc_tree.drawTree()
            
        ymin, ymax = -0.05*mcc_tree.ySpan, 1.01*mcc_tree.ySpan
        for (weeknum,d) in [
                (47, "2019-11-17"),  # start of CDC epi week 2019-47
                (48, "2019-11-24"),  # start of CDC epi week 2019-48
                (49, "2019-12-01"),  # start of CDC epi week 2019-49
                (50, "2019-12-08"),  # start of CDC epi week 2019-50
                (51, "2019-12-15"),  # start of CDC epi week 2019-51
                (52, "2019-12-22"),  # start of CDC epi week 2019-52
                ( 1, "2019-12-29"),  # start of CDC epi week 2020-01
                ( 2, "2020-01-05"),  # start of CDC epi week 2020-02
                ( 3, "2020-01-12"),  # start of CDC epi week 2020-03
                ( 4, "2020-01-19"),  # start of CDC epi week 2020-04
                ( 5, "2020-01-26"),  # start of CDC epi week 2020-05
                ( 6, "2020-02-02"),  # start of CDC epi week 2020-06
                ( 7, "2020-02-09"),  # start of CDC epi week 2020-07
                ( 8, "2020-02-16"),  # start of CDC epi week 2020-08
                ( 9, "2020-02-23"),  # start of CDC epi week 2020-09
                (10, "2020-03-01"),  # start of CDC epi week 2020-10
                (11, "2020-03-08"),  # start of CDC epi week 2020-11
                (12, "2020-03-15"),  # start of CDC epi week 2020-12
                (13, "2020-03-22"),  # start of CDC epi week 2020-13
                (14, "2020-03-29"),  # start of CDC epi week 2020-14
        ]:
            plt.vlines(bt.decimalDate(d), ymin, ymax, linestyle='--', color='#dddddd', lw=0.25)
            if weeknum != 14:
                plt.text(bt.decimalDate(d) + 3.5/365, ymin, str(weeknum), horizontalalignment='center', verticalalignment='bottom',
                         color='lightgrey', fontsize=6)
        
        mcc_tree.plotTree(ax,
                          x_attr=x_attr,
                          colour=lambda node: simple_geo0_2_color[index_2_simple_geo0[node.index]],
                          width=w_func)
        mcc_tree.plotPoints(
            ax,
            x_attr=x_attr,
            size=s_func,
            colour=c_func,
            zorder=100,
            outline=False)
        
        ax.set_xlim(bt.decimalDate("2019-11-17"),   # Start of CDC epi week 2019-47
                    bt.decimalDate("2020-03-29"));  # Start of CDC epi week 2020-14
        ax.set_ylim(ymin, ymax)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.get_yaxis().set_ticks([]);
        if i < 8:
            ax.get_xaxis().set_ticks([])
        else:
            ax.get_xaxis().set_ticks([
                bt.decimalDate("2019-12-01"),  # start of CDC epi week 2019-49
                #bt.decimalDate("2019-12-08"),  # start of CDC epi week 2019-50
                #bt.decimalDate("2019-12-15"),  # start of CDC epi week 2019-51
                #bt.decimalDate("2019-12-22"),  # start of CDC epi week 2019-51
                #bt.decimalDate("2019-12-29"),  # start of CDC epi week 2020-01
                bt.decimalDate("2020-01-05"),  # start of CDC epi week 2020-02
                #bt.decimalDate("2020-01-12"),  # start of CDC epi week 2020-03
                #bt.decimalDate("2020-01-19"),  # start of CDC epi week 2020-04
                #bt.decimalDate("2020-01-26"),  # start of CDC epi week 2020-05
                bt.decimalDate("2020-02-02"),  # start of CDC epi week 2020-06
                #bt.decimalDate("2020-02-09"),  # start of CDC epi week 2020-07
                #bt.decimalDate("2020-02-16"),  # start of CDC epi week 2020-08
                #bt.decimalDate("2020-02-23"),  # start of CDC epi week 2020-09
                bt.decimalDate("2020-03-01"),  # start of CDC epi week 2020-10
                #bt.decimalDate("2020-03-08"),  # start of CDC epi week 2020-11
                #bt.decimalDate("2020-03-15"),  # start of CDC epi week 2020-12
                #bt.decimalDate("2020-03-22"),  # start of CDC epi week 2020-13
                #bt.decimalDate("2020-03-29"),  # start of CDC epi week 2020-14
            ])
            ax.set_xticklabels([
                "1 Dec\n2019",
                "5 Jan\n2020",
                #"12 Jan\n2020",
                #"19 Jan\n2020",
                #"26 Jan\n2020",
                "2 Feb\n2020",
                #"9 Feb\n2020",
                #"16 Feb\n2020",
                #"23 Feb\n2020",
                "1 Mar\n2020",
                #"8 Mar\n2020",
                #"15 Mar\n2020",
                #"22 Mar\n2020",
                #"29 Mar\n2020",
            ])

    # Legend code adapted from https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html
    custom_lines = []
    custom_labels = []
    for simple_geo0, color in simple_geo0_2_color.items():
        custom_lines.append(Line2D([0], [0], color=color, lw=3))
        custom_labels.append(simple_geo0)
    fig.legend(custom_lines, custom_labels, bbox_to_anchor=(0.235, 0.88), framealpha=0.5, fontsize=8)

    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')


Path('plots').mkdir(parents=True, exist_ok=True)

# for run in ('a', 'b'):
#     for epiweek in ['202002',
#                     '202003',
#                     '202004',
#                     '202005',
#                     '202006',
#                     '202007',
#                     '202008',
#                     '202009',
#                     '202010',
#                     '202011',
#                     '202012',
#                     '202013',
#                     ]:
        
#         plot_tree(f'./outputs_by_submission_date_{run}/to_epi_week_{epiweek}/to_epi_week_{epiweek}.mcc',
#                   f'plots/BySubmissionDate{run.upper()}-ToEpiWeek{epiweek}.pdf',
#                   f'Submitted up to CDC epi week {epiweek} (N = NNN)')
        
#         if epiweek <= '202010':
#             plot_tree(f'./outputs_by_collection_date_{run}/to_epi_week_{epiweek}/to_epi_week_{epiweek}.mcc',
#                       f'plots/ByCollectionDate{run.upper()}-ToEpiWeek{epiweek}.pdf',
#                       f'Collected up to CDC epi week {epiweek} (N = NNN)')

# # Assume that pdftk is installed!
# for run in ('a', 'b'):
#     submission_date_epiweeks = ['202002',
#                                 '202003',
#                                 '202004',
#                                 '202005',
#                                 '202006',
#                                 '202007',
#                                 '202008',
#                                 '202009',
#                                 '202010',
#                                 '202011',
#                                 '202012',
#                                 '202013',
#                                 ]
#     collection_date_epiweeks = [w for w in submission_date_epiweeks if w <= '202010']

#     # Assume pdftk is installed
#     outfile = f'plots/BySubmissionDate{run.upper()}-All.pdf'
#     print(f'Collating plots into {outfile}')
#     pdftk_cmdline = (
#         ['pdftk'] +
#         [f'plots/BySubmissionDate{run.upper()}-ToEpiWeek{epiweek}.pdf' for epiweek in submission_date_epiweeks] +
#         ['cat',
#          'output',
#          outfile])
#     subprocess.run(pdftk_cmdline)
    
#     outfile = f'plots/ByCollectionDate{run.upper()}-All.pdf'
#     print(f'Collating plots into {outfile}')
#     pdftk_cmdline = (
#         ['pdftk'] +
#         [f'plots/ByCollectionDate{run.upper()}-ToEpiWeek{epiweek}.pdf' for epiweek in collection_date_epiweeks] +
#         ['cat',
#          'output',
#          outfile])
#     subprocess.run(pdftk_cmdline)


# Compound plot for paper
for run in ('a','b'):
    epiweek_2_tree_file = {
        epiweek: f'./outputs_by_submission_date_{run}/to_epi_week_{epiweek}/to_epi_week_{epiweek}.mcc'
        for epiweek in [
                '202002',
                '202003',
                '202004',
                '202005',
                '202006',
                '202007',
                '202008',
                '202009',
                '202010',
                '202011',
                '202012',
                '202013',
        ]
    }
    
    compound_plot(epiweek_2_tree_file, f'plots/BySubmissionDate{run.upper()}-Compound.pdf')

    epiweek_2_tree_file = {
        epiweek: f'./outputs_by_collection_date_{run}/to_epi_week_{epiweek}/to_epi_week_{epiweek}.mcc'
        for epiweek in [
                '202002',
                '202003',
                '202004',
                '202005',
                '202006',
                '202007',
                '202008',
                '202009',
                '202010',
        ]
    }
    
    compound_plot(epiweek_2_tree_file, f'plots/ByCollectionDate{run.upper()}-Compound.pdf')
