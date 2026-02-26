#!/usr/bin/env python3
"""Analyze WCSS results: run loganalyser, check ESS, compute coverage and ranks.

All results are written to TSV files under analyses/.  This script produces
no plots; see 04_plot.py for visualization.

Uses BEAST 2's loganalyser for posterior summary statistics (mean, HPD, ESS).
Raw posterior samples are read only for Rank Uniformity Validation.

Usage:
    ./03_analyze.py --n 200
    ./03_analyze.py --n 200 --burnin 30
"""

import argparse
import io
import json
import math
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import binom


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

LOGANALYSER = "../../loganalyser277"

# (display_name, log_column_name / loganalyser column prefix)
PARAMS = [
    ("kappa",      "kappa"),
    ("pi_A",       "frequencies1"),
    ("pi_C",       "frequencies2"),
    ("pi_G",       "frequencies3"),
    ("pi_T",       "frequencies4"),
    ("rootHeight", "rootHeight"),
]

# Observables whose ESS should be ignored in the ESS check
# (fixed parameters will have near-zero ESS by construction).
ESS_IGNORE = {"meanRate"}  # mutation rate is fixed in this study

ESS_THRESHOLD = 200


# ---------------------------------------------------------------------------
# Step 1: Run loganalyser
# ---------------------------------------------------------------------------

def run_loganalyser(script_dir, analyses_dir, n, burnin_pct):
    """Run BEAST 2's loganalyser, save raw TSV, return results as DataFrame."""
    log_files = [
        os.path.join("sims", f"sim_{i:03d}", "delphy.log")
        for i in range(n)
    ]
    cmd = [LOGANALYSER, "-oneline", "-burnin", str(burnin_pct)] + log_files
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=script_dir)
    if result.returncode != 0:
        raise RuntimeError(f"loganalyser failed:\n{result.stderr}")

    # Parse TSV output — find the header line starting with "sample\t"
    lines = result.stdout.strip().split("\n")
    header_idx = None
    for idx, line in enumerate(lines):
        if line.startswith("sample\t"):
            header_idx = idx
            break
    if header_idx is None:
        raise RuntimeError(
            f"Could not find loganalyser header in output:\n{result.stdout[:500]}"
        )

    tsv_text = "\n".join(lines[header_idx:])

    # Save raw TSV for later inspection
    tsv_path = os.path.join(analyses_dir, "loganalyser_output.tsv")
    with open(tsv_path, "w") as f:
        f.write(tsv_text + "\n")
    print(f"  Saved {tsv_path}")

    df = pd.read_csv(io.StringIO(tsv_text), sep="\t")
    return df


def check_ess(la_df, analyses_dir, ignore_low_ess=False):
    """Check ESS values, save TSV, and fail if any non-ignored observable is low.

    Returns True if all non-ignored observables pass, False otherwise.
    """
    ess_cols = sorted(c for c in la_df.columns if c.endswith(".ESS"))
    n = len(la_df)

    # Save TSV with all observables
    tsv_path = os.path.join(analyses_dir, "ess_check.tsv")
    with open(tsv_path, "w") as f:
        f.write("Observable\tLow_ESS_runs\tN\tFraction\t"
                "Mean_ESS\tStd_ESS\tIgnored\n")
        for col in ess_cols:
            obs_name = col[:-4]  # strip ".ESS"
            low_count = int((la_df[col] < ESS_THRESHOLD).sum())
            mean_ess = la_df[col].mean()
            std_ess = la_df[col].std()
            ignored = "yes" if obs_name in ESS_IGNORE else "no"
            f.write(f"{obs_name}\t{low_count}\t{n}\t"
                    f"{low_count/n:.3f}\t{mean_ess:.1f}\t{std_ess:.1f}\t"
                    f"{ignored}\n")
    print(f"  Saved {tsv_path}")

    # Print human-readable summary (non-ignored with at least one low-ESS run)
    print(f"\nESS check (threshold={ESS_THRESHOLD}):")
    any_low = False
    max_frac = 0.0
    max_frac_obs = ""
    for col in ess_cols:
        obs_name = col[:-4]
        if obs_name in ESS_IGNORE:
            continue
        low_count = int((la_df[col] < ESS_THRESHOLD).sum())
        frac = low_count / n
        if frac > max_frac:
            max_frac = frac
            max_frac_obs = obs_name
        if low_count > 0:
            if not any_low:
                print(f"{'Observable':<30} {'Low ESS runs':>14} {'Fraction':>10}")
                print(f"{'-'*30} {'-'*14} {'-'*10}")
                any_low = True
            print(f"{obs_name:<30} {low_count:>14} {frac:>10.3f}")

    if not any_low:
        print("  All non-ignored observables have ESS >= "
              f"{ESS_THRESHOLD} across all {n} replicates.")
        return True

    if max_frac < 0.05:
        print(f"\n  WARNING: Low ESS detected in a few runs"
              f" (observable with highest low-ESS fraction: {max_frac_obs}"
              f" at {max_frac:.1%}), but below 5% threshold.  Continuing.")
        return True

    if ignore_low_ess:
        print("\n  WARNING: Low ESS detected but continuing (--ignore-low-ess).")
        return True

    print(f"\n  ERROR: Some non-ignored observables have ESS < {ESS_THRESHOLD}"
          f" in more than 5% of replicates.")
    print("  To fix, increase the number of MCMC steps and rerun:")
    print("    ./01_generate.py --n <N> --steps <more_steps>")
    print("    make -C sims clean && ./02_run.py")
    print("    ./03_analyze.py --n <N>")
    print("  Or pass --ignore-low-ess to continue anyway.")
    return False


# ---------------------------------------------------------------------------
# Step 2: Read true parameters
# ---------------------------------------------------------------------------

def read_true_params(sim_dir):
    """Read true parameter values from Sapling's info JSON."""
    info_path = os.path.join(sim_dir, "sim_info.json")
    with open(info_path) as f:
        info = json.load(f)

    return {
        "kappa": info["subst_model"]["kappa"],
        "pi_A": info["subst_model"]["pi"][0],
        "pi_C": info["subst_model"]["pi"][1],
        "pi_G": info["subst_model"]["pi"][2],
        "pi_T": info["subst_model"]["pi"][3],
        "rootHeight": info["tree_stats"]["tree_height"],
    }


# ---------------------------------------------------------------------------
# Step 3: Coverage analysis
# ---------------------------------------------------------------------------

def compute_coverage(true_vals, la_df, params):
    """Compute 95% HPD coverage using loganalyser output."""
    n = len(true_vals)
    results = {}

    for name, col in params:
        in_hpd = 0
        for i in range(n):
            lo = la_df[f"{col}.95%HPDlo"].iloc[i]
            hi = la_df[f"{col}.95%HPDup"].iloc[i]
            if lo <= true_vals[i][name] <= hi:
                in_hpd += 1
        results[name] = in_hpd / n

    return results


# ---------------------------------------------------------------------------
# Step 4: Rank Uniformity Validation (RUV)
# ---------------------------------------------------------------------------

def compute_normalized_ranks(true_vals, script_dir, n, burnin_frac, params):
    """Compute normalized rank of true value among posterior samples."""
    param_names = [name for name, _ in params]
    ranks = {name: [] for name in param_names}

    for i in range(n):
        log_path = os.path.join(script_dir, "sims", f"sim_{i:03d}", "delphy.log")
        raw = pd.read_table(log_path, comment="#")
        burnin_rows = math.floor(burnin_frac * len(raw))
        data = raw.iloc[burnin_rows:]

        for name, col in params:
            true_val = true_vals[i][name]
            samples = data[col].values
            rank = np.sum(samples < true_val)
            L = len(samples)
            ranks[name].append(rank / L)

    return ranks


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Analyze WCSS results")
    parser.add_argument("--n", type=int, default=200,
                        help="Number of simulation replicates (default: 200)")
    parser.add_argument("--burnin", type=int, default=30,
                        help="Burnin percentage (default: 30)")
    parser.add_argument("--ignore-low-ess", action="store_true",
                        help="Continue even if some observables have low ESS")
    args = parser.parse_args()

    n = args.n
    burnin_pct = args.burnin
    burnin_frac = burnin_pct / 100.0
    script_dir = os.path.dirname(os.path.abspath(__file__))
    analyses_dir = os.path.join(script_dir, "analyses")
    os.makedirs(analyses_dir, exist_ok=True)

    param_names = [name for name, _ in PARAMS]

    # Step 1: Run loganalyser and ESS check
    print(f"Running loganalyser on {n} replicates (burnin={burnin_pct}%)...")
    la_df = run_loganalyser(script_dir, analyses_dir, n, burnin_pct)
    print(f"  loganalyser returned {len(la_df)} rows")
    if not check_ess(la_df, analyses_dir, args.ignore_low_ess):
        sys.exit(1)

    # Step 2: Read true parameters and save to TSV
    print("\nReading true parameters...")
    true_vals = []
    for i in range(n):
        sim_dir = os.path.join(script_dir, "sims", f"sim_{i:03d}")
        true_vals.append(read_true_params(sim_dir))

    true_path = os.path.join(analyses_dir, "true_params.tsv")
    with open(true_path, "w") as f:
        f.write("replicate\tkappa\tpi_A\tpi_C\tpi_G\tpi_T\trootHeight\n")
        for i, tv in enumerate(true_vals):
            f.write(f"{i}\t{tv['kappa']}\t{tv['pi_A']}\t{tv['pi_C']}\t"
                    f"{tv['pi_G']}\t{tv['pi_T']}\t{tv['rootHeight']}\n")
    print(f"  Saved {true_path}")

    # Step 3: Coverage analysis
    print()
    coverage = compute_coverage(true_vals, la_df, PARAMS)
    lo_binom = binom.ppf(0.025, n, 0.95) / n
    hi_binom = binom.ppf(0.975, n, 0.95) / n

    # Print human-readable table
    print(f"{'Parameter':<14} {'Coverage':>10}")
    print(f"{'-'*14} {'-'*10}")
    for name in param_names:
        print(f"{name:<14} {coverage[name]:>10.2f}")
    print(f"\nExpected coverage: 0.95.  "
          f"Binomial 95% interval for N={n}: [{lo_binom:.3f}, {hi_binom:.3f}].")

    # Save TSV
    cov_path = os.path.join(analyses_dir, "coverage_summary.txt")
    with open(cov_path, "w") as f:
        f.write("Parameter\tCoverage\tN\tExpected\tBinom_2.5%\tBinom_97.5%\n")
        for name in param_names:
            f.write(f"{name}\t{coverage[name]:.2f}\t{n}\t"
                    f"0.950\t{lo_binom:.3f}\t{hi_binom:.3f}\n")
    print(f"  Saved {cov_path}")

    # Step 4: RUV — compute and save normalized ranks
    print(f"\nComputing normalized ranks...")
    ranks = compute_normalized_ranks(true_vals, script_dir, n, burnin_frac,
                                     PARAMS)

    ranks_path = os.path.join(analyses_dir, "ranks.tsv")
    with open(ranks_path, "w") as f:
        f.write("replicate\t" + "\t".join(param_names) + "\n")
        for i in range(n):
            vals = "\t".join(str(ranks[name][i]) for name in param_names)
            f.write(f"{i}\t{vals}\n")
    print(f"  Saved {ranks_path}")

    print("\nDone.  Run 04_plot.py to generate plots.")


if __name__ == "__main__":
    main()
