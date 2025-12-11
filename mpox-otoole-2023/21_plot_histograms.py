#!/usr/bin/env python

import baltic as bt
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import math
from pathlib import Path
from scipy import stats

t0 = bt.decimalDate("2022-08-15")  # Hard-coded, but ok

def extract_log_params(log_filename, burnin):
    # Extremely brittle way of extracting distributions from the log files
    raw_data = pd.read_table(log_filename, comment='#')
    num_pts = len(raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    data = raw_data[burnin_pts:]
    print(f'Read in log file {log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(data)} usable')

    result = {}
    for key in [
            'exponential.popSize',
            'exponential.growthRate',
            'apobec3.clock.rate',
            'non_apobec3.clock.rate',
    ]:
        result[key] = list(data[key])

    return result

def plot_growth_curves(logs, out_pdf_filename):
    years = np.linspace(2014.0, 2023.0, 1000)
    pop_min = []
    pop_2p5 = []
    pop_50 = []
    pop_97p5 = []
    pop_max = []
    for year in years:
        pops = sorted([N0 * np.exp(g*(year-t0))
                       for (N0, g) in zip(logs['exponential.popSize'], logs['exponential.growthRate'])])
        pop_min.append(pops[0])
        pop_2p5.append(pops[int(0.025*len(pops))])
        pop_50.append(pops[int(0.50*len(pops))])
        pop_97p5.append(pops[int(0.975*len(pops))])
        pop_max.append(pops[-1])

    fig,ax = plt.subplots(figsize=(5,4),facecolor='w')

    #ax.fill_between(years, pop_min, pop_max, color='#ff0000', alpha=0.2)
    ax.fill_between(years, pop_2p5, pop_97p5, color='#000000', alpha=0.2)
    ax.plot(years, pop_50, color='black', linestyle='--')
    
    ax.set_xlim(2014.5, 2022.7)
    ax.set_ylim(5e-5, 3.2e3);
    ax.set_yscale('log')
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.get_yaxis().set_ticks([
        1e-4,
        1e-3,
        1e-2,
        1e-1,
        1e0,
        1e1,
        1e2,
        1e3,
    ]);
    ax.get_yaxis().set_ticklabels([
        '10$^{-4}$',
        '10$^{-3}$',
        '10$^{-2}$',
        '10$^{-1}$',
        '10$^{0}$',
        '10$^{1}$',
        '10$^{2}$',
        '10$^{3}$',
    ]);
    ax.set_ylabel('Effective population size (years)')
    
    ax.get_xaxis().set_ticks([
        bt.decimalDate("2015-01-01"),
        bt.decimalDate("2016-01-01"),
        bt.decimalDate("2017-01-01"),
        bt.decimalDate("2018-01-01"),
        bt.decimalDate("2019-01-01"),
        bt.decimalDate("2020-01-01"),
        bt.decimalDate("2021-01-01"),
        bt.decimalDate("2022-01-01"),
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
       ]);
    ax.set_xlabel('Time (years)')

    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')

def plot_clock_rates(logs, out_pdf_filename):
    fig,ax = plt.subplots(figsize=(4,4),facecolor='w')

    mean_apobec = np.sum(logs['apobec3.clock.rate']) / len(logs['apobec3.clock.rate'])
    mean_non_apobec = np.sum(logs['non_apobec3.clock.rate']) / len(logs['non_apobec3.clock.rate'])
    
    bplot = ax.boxplot([logs['apobec3.clock.rate'],
                        logs['non_apobec3.clock.rate']],
                       widths=0.7,
                       patch_artist=True,
                       tick_labels=['APOBEC3', 'non-APOBEC3'])

    for patch, color in zip(bplot['boxes'], ['#926568', '#f4e7c8']):
        patch.set_facecolor(color)
    for patch, median, color in zip(bplot['boxes'], bplot['medians'], ['#ded2d2', '#dab965']):
        patch.set_edgecolor(color)
        median.set_color(color)
    for whisker, color in zip(bplot['whiskers'], ['#ded2d2', '#ded2d2', '#dab965', '#dab965']):
        whisker.set_color(color)
    for cap, color in zip(bplot['caps'], ['#926568', '#926568', '#f4e7c8', '#f4e7c8']):
        cap.set_color(color)
    for flier, color in zip(bplot['fliers'], ['#926568', '#f4e7c8']):
        flier.set_color(color)
        flier.set_marker('*')
        flier.set_markersize(1.0)
    
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    
    ax.set_ylim(0, 17.6e-5);
    ax.get_yaxis().set_ticks([
           0e-5,
         2.5e-5,
         5.0e-5,
         7.5e-5,
        10.0e-5,
        12.5e-5,
        15.0e-5,
        17.5e-5,
    ]);
    ax.get_yaxis().set_ticklabels([
        '0.000000',
        '0.000025',
        '0.000050',
        '0.000075',
        '0.000100',
        '0.000125',
        '0.000150',
        '0.000175',
    ])
    ax.set_ylabel('Clock rate')
    ax.set_xlabel('Clock')

    plt.subplots_adjust(left=0.3, right=0.9, top=0.9, bottom=0.1)

    plt.text(1, mean_apobec+0.000040, f'{mean_apobec*1e4:.1f}' + ' x 10$^{{-4}}$', horizontalalignment='center')
    plt.text(2, mean_non_apobec+0.000008, f'{mean_non_apobec*1e6:.1f}' + ' x 10$^{{-6}}$', horizontalalignment='center')
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')
    
def plot_doubling_time_dist(logs, out_pdf_filename):
    fig,ax = plt.subplots(figsize=(5,2.5),facecolor='w')

    # e^g = 2^(1/tdbl)  => tdbl = ln(2) / g
    doubling_times = np.log(2) / np.array(logs['exponential.growthRate'])
    times = np.linspace(0, 5, 200)
    kde = stats.gaussian_kde(doubling_times)
    ax.fill_between(times, kde(times), 0, color='#e5d7d8')
    ax.plot(times, kde(times), color='#9e666a', lw=1.5)

    mean_t_dbl = np.mean(doubling_times)
    plt.vlines(mean_t_dbl, 0, 1.0, linestyle='--', color='black', lw=1)
    plt.text(mean_t_dbl, 1.1, f'Mean: {mean_t_dbl:.2f} years', horizontalalignment='center')
    
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    
    ax.set_ylabel('Density')
    ax.set_ylim(0, 1.6);
    ax.get_yaxis().set_ticks([
        0.0,
        0.5,
        1.0,
        1.5
    ]);
    ax.get_yaxis().set_ticklabels([
        '0.0',
        '0.5',
        '1.0',
        '1.5',
    ])
    
    ax.set_xlabel('Time (years)')
    ax.set_xlim(0, 5)
    ax.get_xaxis().set_ticks([
        0,
        1,
        2,
        3,
        4,
        5,
    ]);
    ax.get_xaxis().set_ticklabels([
        '0',
        '1',
        '2',
        '3',
        '4',
        '5',
    ])

    plt.subplots_adjust(left=0.12, right=0.93, top=0.9, bottom=0.2)
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')
    
Path('plots').mkdir(parents=True, exist_ok=True)

for run in ('a', 'b'):
    logs = extract_log_params(f'delphy_outputs_{run}/mpox-otoole-2023.log', burnin=0.3)
    plot_growth_curves(logs, f'plots/{run.upper()}GrowthCurve.pdf')
    plot_clock_rates(logs, f'plots/{run.upper()}ClockRates.pdf')
    plot_doubling_time_dist(logs, f'plots/{run.upper()}DoublingTimeDistr.pdf')
