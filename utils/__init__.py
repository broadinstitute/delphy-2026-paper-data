## Horrendously quick-n-dirty module to share many common tidbits in the analysis and plotting scripts

from typing import NamedTuple
import pandas as pd
import json
from pathlib import Path
import math

# Log formats = internal observable name -> name in log
delphy_log_format = {
    'joint': 'posterior_for_Delphy',
    'prior': 'prior_for_Delphy',
    'likelihood': 'likelihood_really_logG',
    'age(root)': 'age(root)',
    'treeLength': 'treeLength',
    'clockRate': 'clock.rate',
    'TreeHeight': 'rootHeight',
    'kappa': 'kappa',
    'gammaShape': 'alpha',
    'freqParameter.1': 'frequencies1',
    'freqParameter.2': 'frequencies2',
    'freqParameter.3': 'frequencies3',
    'freqParameter.4': 'frequencies4',
    'CoalescentExponential': 'coalescent',
    'ePopSize': 'exponential.popSize',
    'growthRate': 'exponential.growthRate',
    'skygrid.precision': 'skygrid.precision',
}
beast2_log_format = {
    'joint': 'posterior',
    'prior': 'prior',
    'likelihood': 'likelihood',
    # No 'age(root)
    'treeLength': 'Tree.treeLength',
    'clockRate': 'clockRate',
    'TreeHeight': 'Tree.height',
    'kappa': 'kappa',
    'gammaShape': 'gammaShape',
    'freqParameter.1': 'freqParameter.1',
    'freqParameter.2': 'freqParameter.2',
    'freqParameter.3': 'freqParameter.3',
    'freqParameter.4': 'freqParameter.4',
    'CoalescentExponential': 'Coalescent',
    'ePopSize': 'ePopSize',
    'growthRate': 'growthRate',
}
beastX_log_format = {
    'joint': 'joint',
    'prior': 'prior',
    'likelihood': 'likelihood',
    'age(root)': 'age(root)',
    'treeLength': 'treeLength',
    'clockRate': 'clock.rate',
    'TreeHeight': 'rootHeight',
    'kappa': 'kappa',
    'gammaShape': 'alpha',
    'freqParameter.1': 'frequencies1',
    'freqParameter.2': 'frequencies2',
    'freqParameter.3': 'frequencies3',
    'freqParameter.4': 'frequencies4',
    'CoalescentExponential': 'coalescent',
    'ePopSize': 'exponential.popSize',
    'growthRate': 'exponential.growthRate',
    'skygrid.precision': 'skygrid.precision',
}

for i in range(51):
    delphy_log_format[f'skygrid.logPopSize{i}'] = f'skygrid.logPopSize{i}'
    beastX_log_format[f'skygrid.logPopSize{i}'] = f'skygrid.logPopSize{i}'

class Log(NamedTuple):
    filename: str
    log_analysis_filename: str
    tree_ess_filename: str
    log_format: dict[str, str]
    label: str
    color: str

def read_in_ess_table(log_analysis_filename, tree_ess_filename, log_format):
    print(f'Reading in {log_analysis_filename} and {tree_ess_filename}')
        
    log_analysis = pd.read_fwf(log_analysis_filename, index_col=0)
    result = {}
    for k, log_k in log_format.items():
        if log_k in log_analysis['ESS']:
            result[k] = log_analysis['ESS'][log_k]

    
    if tree_ess_filename:
        with open(tree_ess_filename, 'r') as f:
            tree_ess = json.load(f)
        result['TreeESS'] = tree_ess['chain_results'][0]['effective_sample_size']

    return result

def extract_ml_results(iqtree_log_filename, tt_folder_name, look_for_alpha=False):
    """Extremely brittle IQ-Tree log file parser, enough for these plots though"""

    results = {}
    
    with open(iqtree_log_filename, 'r') as f:
        lines = f.readlines()

    # Extract HKY kappa
    # -----------------
    # Looking for a section like this:
    #
    # Rate parameter R:
    # 
    #   A-C: 1.0000
    #   A-G: 9.0197
    #   A-T: 1.0000
    #   C-G: 1.0000
    #   C-T: 9.0197
    #   G-T: 1.0000
    #
    rate_line_idx = [i
                     for i in range(len(lines))
                     if lines[i] == 'Rate parameter R:\n'][0]
    results['kappa'] = float(lines[rate_line_idx+3].strip().split(':')[1])

    # Extract HKY stationary state frequencies
    # ----------------------------------------
    # Looking for a section like this:
    #
    # State frequencies: (estimated with maximum likelihood)
    # 
    #   pi(A) = 0.3192
    #   pi(C) = 0.214
    #   pi(G) = 0.198
    #   pi(T) = 0.2688
    #
    state_freqs_line_idx = [i
                            for i in range(len(lines))
                            if lines[i].startswith('State frequencies')][0]
    results['freqParameter.1'] = float(lines[state_freqs_line_idx+2].split('=')[1])
    results['freqParameter.2'] = float(lines[state_freqs_line_idx+3].split('=')[1])
    results['freqParameter.3'] = float(lines[state_freqs_line_idx+4].split('=')[1])
    results['freqParameter.4'] = float(lines[state_freqs_line_idx+5].split('=')[1])

    # Extract site-rate heterogeneity parameter alpha
    # -----------------------------------------------
    # Looking for a section like this:
    #
    # Model of rate heterogeneity: Gamma with 4 categories
    # Gamma shape alpha: 998.4
    #
    if look_for_alpha:
        alpha_line_idx = [i
                          for i in range(len(lines))
                          if lines[i].startswith('Gamma shape alpha')][0]
        results['gammaShape'] = float(lines[alpha_line_idx].split(':')[1])

    # Extract mutation rate
    # ---------------------
    # Looking in tt/molecular_clock.txt for a line like this:
    #
    # --rate:	1.785e-03
    #
    # Depending on TreeTime options, there may be an associated range, e.g.,
    #
    # --rate:	5.678e-04 +/- 1.03e-04 (one std-dev)
    #
    with (Path(tt_folder_name) / 'molecular_clock.txt').open('r') as f:
        for line in f:
            if 'rate' in line:
                results['clockRate'] = float(line.split(':')[1].split()[0])

    # Extract tree height
    # -------------------
    # Read all node dates from tt/dates.tsv and extract rate
    with (Path(tt_folder_name) / 'dates.tsv').open('r') as f:
        dates = []
        for line in f:
            if line.startswith('#'):
                continue
            date_str = line.split()[-1]
            if '--' in date_str:
                continue   # Missing dates ?
            dates.append(float(date_str))

        results['TreeHeight'] = max(dates) - min(dates)

    # Extract coalescent prior parameters
    # -----------------------------------
    # We interpret tt/skyline.tsv in two ways.  One is to read
    # all the raw data (rescaled to units of years) into results['skyline']
    # (a dict of `linear year` -> (N_e, lower, upper)).  The other is to
    # take the first and last points of this skyline as the endpoints of
    # an exponential growth curve and derive the final pop size & growth rate;
    # this is to compare directly against an exponential growth coalescent
    # when not using Skygrid in Delphy.
    #
    # tt/skyline.tsv looks like this:
    #
    # #Skyline assuming 50.0 gen/year and approximate confidence bounds (+/- 2.000000 standard deviations of the LH)
    # #date 	N_e 	lower 	upper
    # 2014.175	6.230e+00	2.804e+00	1.384e+01
    # 2014.462	1.254e+01	9.726e+00	1.617e+01
    skyline = {}
    with (Path(tt_folder_name) / 'skyline.tsv').open('r') as f:
        lines = f.readlines()
        assert len(lines) >= 4
        rho = 1./50.  # Default TreeTime assumption = 50 generations per year; rho = generation time in years

        for line in lines[2:]:
            if not line.strip():
                continue
            linear_year, Ne, lower, upper = map(float, line.split())
            Ne *= rho
            lower *= rho
            upper *= rho
            skyline[linear_year] = (Ne, lower, upper)

        results['skyline'] = skyline
        
        min_date = min(skyline.keys())
        max_date = max(skyline.keys())
        
        min_Ne_rho, *_ = skyline[min_date]
        max_Ne_rho, *_ = skyline[max_date]

        results['ePopSize'] = max_Ne_rho
        results['growthRate'] = math.log(max_Ne_rho/min_Ne_rho) / (max_date - min_date)
                
    return results


def median_and_95HPD(xs):
    sorted_xs = sorted(xs)

    L = len(sorted_xs)
    if L % 2 == 0:
        median = 0.5*(sorted_xs[L // 2 - 1] + sorted_xs[L // 2])
    else:
        median = sorted_xs[L // 2]
        
    hpd_size = int(len(xs)*0.95)
    hpd = (sorted_xs[0], sorted_xs[-1])
    best_hpd_width = sorted_xs[-1] - sorted_xs[0]
    for i in range(len(sorted_xs)-hpd_size):
        hpd_i_width = sorted_xs[i+hpd_size] - sorted_xs[i]
        if hpd_i_width < best_hpd_width:
            best_hpd_width = hpd_i_width
            hpd = (sorted_xs[i], sorted_xs[i+hpd_size])


    return median, hpd

def mrca_node_of(tree, leaf_names):
    leaf_nodes = [node for node in tree.Objects if node.is_leaf() and node.name.split('|')[0] in leaf_names]
    if len(leaf_nodes) == 1:
        return leaf_nodes[0]
    else:
        return tree.commonAncestor(leaf_nodes)

