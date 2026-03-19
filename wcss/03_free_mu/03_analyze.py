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
import concurrent.futures
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
CALC_CLADE_COVERAGE = "../../tree_ess/target/release/calc-clade-coverage"

# (display_name, log_column_name / loganalyser column prefix)
PARAMS = [
    ("mu",         "meanRate"),
    ("kappa",      "kappa"),
    ("pi_A",       "frequencies1"),
    ("pi_C",       "frequencies2"),
    ("pi_G",       "frequencies3"),
    ("pi_T",       "frequencies4"),
    ("rootHeight", "rootHeight"),
]

# Observables whose ESS should be ignored in the ESS check
ESS_IGNORE = set()  # no fixed parameters to ignore in this study

ESS_THRESHOLD_LOW = 200
ESS_THRESHOLD_VERY_LOW = 150


# ---------------------------------------------------------------------------
# Step 0: Digest log files
# ---------------------------------------------------------------------------

def digest_log_file(src_path, dst_path):
    """Produce delphy-digested.log from delphy.log (symlink, no-op)."""
    if os.path.islink(dst_path) or os.path.exists(dst_path):
        os.remove(dst_path)
    os.symlink(os.path.basename(src_path), dst_path)


# ---------------------------------------------------------------------------
# Step 1: Run loganalyser
# ---------------------------------------------------------------------------

def run_loganalyser(script_dir, analyses_dir, n, burnin_pct):
    """Run BEAST 2's loganalyser, save raw TSV, return results as DataFrame."""
    log_files = [
        os.path.join("sims", f"sim_{i:03d}", "delphy-digested.log")
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
    """Check ESS values, save TSV, and return (ok, excluded_indices).

    Two-tier check:
    - "low ESS": ESS < 200.  Error if > 10% of runs.
    - "very low ESS": ESS < 150.  Error if > 5% of runs.
      Runs with very low ESS on any non-ignored observable are excluded
      from downstream analysis.

    Returns (True, excluded_set) on success, (False, set()) on error.
    """
    ess_cols = sorted(c for c in la_df.columns if c.endswith(".ESS"))
    n = len(la_df)

    # Save TSV with all observables
    tsv_path = os.path.join(analyses_dir, "ess_check.tsv")
    with open(tsv_path, "w") as f:
        f.write("Observable\tLow_ESS_runs\tVery_low_ESS_runs\tN\t"
                "Low_frac\tVery_low_frac\tMean_ESS\tStd_ESS\tIgnored\n")
        for col in ess_cols:
            obs_name = col[:-4]  # strip ".ESS"
            low_count = int((la_df[col] < ESS_THRESHOLD_LOW).sum())
            very_low_count = int((la_df[col] < ESS_THRESHOLD_VERY_LOW).sum())
            mean_ess = la_df[col].mean()
            std_ess = la_df[col].std()
            ignored = "yes" if obs_name in ESS_IGNORE else "no"
            f.write(f"{obs_name}\t{low_count}\t{very_low_count}\t{n}\t"
                    f"{low_count/n:.3f}\t{very_low_count/n:.3f}\t"
                    f"{mean_ess:.1f}\t{std_ess:.1f}\t{ignored}\n")
    print(f"  Saved {tsv_path}")

    # Print human-readable summary (non-ignored with at least one low-ESS run)
    print(f"\nESS check (low={ESS_THRESHOLD_LOW}, "
          f"very low={ESS_THRESHOLD_VERY_LOW}):")
    any_low = False
    max_low_frac = 0.0
    max_very_low_frac = 0.0
    for col in ess_cols:
        obs_name = col[:-4]
        if obs_name in ESS_IGNORE:
            continue
        low_count = int((la_df[col] < ESS_THRESHOLD_LOW).sum())
        very_low_count = int((la_df[col] < ESS_THRESHOLD_VERY_LOW).sum())
        low_frac = low_count / n
        very_low_frac = very_low_count / n
        max_low_frac = max(max_low_frac, low_frac)
        max_very_low_frac = max(max_very_low_frac, very_low_frac)
        if low_count > 0:
            if not any_low:
                print(f"{'Observable':<30} {'<200':>8} {'<150':>8} "
                      f"{'%<200':>8} {'%<150':>8}")
                print(f"{'-'*30} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
                any_low = True
            print(f"{obs_name:<30} {low_count:>8} {very_low_count:>8} "
                  f"{low_frac:>8.1%} {very_low_frac:>8.1%}")

    if not any_low:
        print(f"  All non-ignored observables have ESS >= {ESS_THRESHOLD_LOW} "
              f"across all {n} replicates.")
        return (True, set())

    # Check error conditions
    error = False
    if max_low_frac > 0.10:
        error = True
        print(f"\n  ERROR: Some non-ignored observables have ESS < "
              f"{ESS_THRESHOLD_LOW} in > 10% of replicates.")
    if max_very_low_frac > 0.05:
        error = True
        print(f"\n  ERROR: Some non-ignored observables have ESS < "
              f"{ESS_THRESHOLD_VERY_LOW} in > 5% of replicates.")

    if error and not ignore_low_ess:
        print("  To fix, increase the number of MCMC steps and rerun:")
        print("    ./01_generate.py --n <N> --steps <more_steps>")
        print("    make -C sims clean && ./02_run.py")
        print("    ./03_analyze.py --n <N>")
        print("  Or pass --ignore-low-ess to continue anyway.")
        return (False, set())

    # Compute excluded replicates: any row where a non-ignored observable
    # has ESS < ESS_THRESHOLD_VERY_LOW
    excluded = set()
    for col in ess_cols:
        obs_name = col[:-4]
        if obs_name in ESS_IGNORE:
            continue
        for idx in la_df.index[la_df[col] < ESS_THRESHOLD_VERY_LOW]:
            excluded.add(idx)

    if error:
        print(f"\n  WARNING: Low ESS detected but continuing "
              f"(--ignore-low-ess).")
    if excluded:
        print(f"  Excluding {len(excluded)} replicate(s) with very low ESS "
              f"from downstream analysis.")
    else:
        print(f"\n  WARNING: Some runs have low ESS but none below "
              f"{ESS_THRESHOLD_VERY_LOW}.  Continuing with all replicates.")

    return (True, excluded)


# ---------------------------------------------------------------------------
# Step 2: Read true parameters
# ---------------------------------------------------------------------------

def read_true_params(sim_dir):
    """Read true parameter values from Sapling's info JSON."""
    info_path = os.path.join(sim_dir, "sim_info.json")
    with open(info_path) as f:
        info = json.load(f)

    return {
        "mu": info["subst_model"]["mu"],
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

def compute_normalized_ranks(true_vals, script_dir, included, burnin_frac,
                             params):
    """Compute normalized rank of true value among posterior samples."""
    param_names = [name for name, _ in params]
    ranks = {name: [] for name in param_names}

    for j, i in enumerate(included):
        log_path = os.path.join(script_dir, "sims", f"sim_{i:03d}", "delphy-digested.log")
        raw = pd.read_table(log_path, comment="#")
        burnin_rows = math.floor(burnin_frac * len(raw))
        data = raw.iloc[burnin_rows:]

        for name, col in params:
            true_val = true_vals[j][name]
            samples = data[col].values
            rank = np.sum(samples < true_val)
            L = len(samples)
            ranks[name].append(rank / L)

    return ranks


# ---------------------------------------------------------------------------
# Step 5: Clade coverage
# ---------------------------------------------------------------------------

def _run_one_clade_coverage(i, script_dir, burnin_pct, num_bins, header):
    """Run calc-clade-coverage for a single replicate.  Returns (i, stdout)."""
    sim_tag = f"sim_{i:03d}"
    true_tree = os.path.join("sims", sim_tag, "sim.nwk")
    posterior_trees = os.path.join("sims", sim_tag, "delphy.trees")

    cmd = [
        CALC_CLADE_COVERAGE,
        "--burnin-pct", str(burnin_pct),
        "--num-bins", str(num_bins),
        "--replicate", sim_tag,
        true_tree,
        posterior_trees,
    ]
    if header:
        cmd.append("--header")

    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=script_dir)
    if result.returncode != 0:
        raise RuntimeError(
            f"calc-clade-coverage failed for {sim_tag}:\n{result.stderr}")

    return (i, result.stdout.strip())


def compute_clade_coverage(script_dir, analyses_dir, included, burnin_pct,
                           num_bins):
    """Run calc-clade-coverage per replicate, aggregate, and save results."""
    n = len(included)
    raw_lines = [None] * n
    num_workers = os.cpu_count()
    done = 0

    with concurrent.futures.ProcessPoolExecutor(
            max_workers=num_workers) as executor:
        futures = {
            executor.submit(
                _run_one_clade_coverage, i, script_dir, burnin_pct,
                num_bins, header=(j == 0)): j
            for j, i in enumerate(included)
        }
        for future in concurrent.futures.as_completed(futures):
            j = futures[future]
            _, line = future.result()
            raw_lines[j] = line
            done += 1
            if done % 50 == 0 or done == n:
                print(f"  {done}/{n} replicates processed")

    # Save raw per-replicate TSV
    raw_path = os.path.join(analyses_dir, "clade_coverage_raw.tsv")
    with open(raw_path, "w") as f:
        f.write("\n".join(raw_lines) + "\n")
    print(f"  Saved {raw_path}")

    # Parse and aggregate
    raw_df = pd.read_csv(io.StringIO("\n".join(raw_lines)), sep="\t")

    # Extract bin boundaries from column names (totals_0_10, totals_10_20, ...).
    # Preserve the column order from calc-clade-coverage output, which is
    # already in ascending bin order.
    total_cols = [c for c in raw_df.columns if c.startswith("totals_")]
    bins = []
    for col in total_cols:
        parts = col.replace("totals_", "").split("_")
        lo, hi = int(parts[0]), int(parts[1])
        bins.append((lo, hi))

    agg_rows = []
    for lo, hi in bins:
        total = int(raw_df[f"totals_{lo}_{hi}"].sum())
        true_hits = int(raw_df[f"true_hits_{lo}_{hi}"].sum())
        frac = true_hits / total if total > 0 else float("nan")
        agg_rows.append((lo, hi, total, true_hits, frac))

    # Save aggregated TSV
    agg_path = os.path.join(analyses_dir, "clade_coverage.tsv")
    with open(agg_path, "w") as f:
        f.write("bin_lo\tbin_hi\ttotal\ttrue_hits\tfraction\n")
        for lo, hi, total, true_hits, frac in agg_rows:
            f.write(f"{lo}\t{hi}\t{total}\t{true_hits}\t{frac:.6f}\n")
    print(f"  Saved {agg_path}")

    # Print human-readable summary
    print(f"\nClade coverage (aggregated across {n} replicates):")
    print(f"  {'Bin':>10} {'Total':>8} {'True hits':>10} {'Fraction':>10}")
    print(f"  {'-'*10} {'-'*8} {'-'*10} {'-'*10}")
    for lo, hi, total, true_hits, frac in agg_rows:
        print(f"  {lo:>3}-{hi:<3}%   {total:>8} {true_hits:>10} {frac:>10.4f}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Analyze WCSS results")
    parser.add_argument("--n", type=int, default=200,
                        help="Number of simulation replicates (default: 200)")
    parser.add_argument("--burnin", type=int, default=30,
                        help="Burnin percentage (default: 30)")
    parser.add_argument("--num-coverage-bins", type=int, default=20,
                        help="Number of bins for clade coverage (default: 20)")
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

    # Step 0: Digest log files
    print("Digesting log files...")
    for i in range(n):
        sim_dir = os.path.join(script_dir, "sims", f"sim_{i:03d}")
        digest_log_file(
            os.path.join(sim_dir, "delphy.log"),
            os.path.join(sim_dir, "delphy-digested.log"))
    print(f"  Digested {n} log files")

    # Step 1: Run loganalyser and ESS check
    print(f"\nRunning loganalyser on {n} replicates (burnin={burnin_pct}%)...")
    la_df = run_loganalyser(script_dir, analyses_dir, n, burnin_pct)
    print(f"  loganalyser returned {len(la_df)} rows")
    ok, excluded = check_ess(la_df, analyses_dir, args.ignore_low_ess)
    if not ok:
        sys.exit(1)

    # Build included replicate list and filter la_df
    included = sorted(set(range(n)) - excluded)
    n_eff = len(included)
    if excluded:
        print(f"\n  {len(excluded)} replicate(s) excluded, "
              f"{n_eff} remaining for analysis.")
        # Save excluded replicates
        excl_path = os.path.join(analyses_dir, "excluded_replicates.tsv")
        ess_cols = sorted(c for c in la_df.columns if c.endswith(".ESS"))
        with open(excl_path, "w") as f:
            f.write("replicate\tlow_ess_observables\n")
            for i in sorted(excluded):
                low_obs = []
                for col in ess_cols:
                    obs_name = col[:-4]
                    if obs_name in ESS_IGNORE:
                        continue
                    if la_df[col].iloc[i] < ESS_THRESHOLD_VERY_LOW:
                        low_obs.append(obs_name)
                f.write(f"sim_{i:03d}\t{','.join(low_obs)}\n")
        print(f"  Saved {excl_path}")

    la_df = la_df.iloc[included].reset_index(drop=True)

    # Save filtered loganalyser output so that all downstream TSV files
    # in analyses/ are row-aligned.
    filtered_path = os.path.join(analyses_dir,
                                 "loganalyser_output_filtered.tsv")
    la_df.to_csv(filtered_path, sep="\t", index=False)
    print(f"  Saved {filtered_path} ({n_eff} rows)")

    # Step 2: Read true parameters and save to TSV
    print("\nReading true parameters...")
    true_vals = []
    for i in included:
        sim_dir = os.path.join(script_dir, "sims", f"sim_{i:03d}")
        true_vals.append(read_true_params(sim_dir))

    true_path = os.path.join(analyses_dir, "true_params.tsv")
    with open(true_path, "w") as f:
        f.write("replicate\t" + "\t".join(param_names) + "\n")
        for i, tv in zip(included, true_vals):
            vals = "\t".join(str(tv[name]) for name in param_names)
            f.write(f"sim_{i:03d}\t{vals}\n")
    print(f"  Saved {true_path}")

    # Step 3: Coverage analysis
    print()
    coverage = compute_coverage(true_vals, la_df, PARAMS)
    lo_binom = binom.ppf(0.025, n_eff, 0.95) / n_eff
    hi_binom = binom.ppf(0.975, n_eff, 0.95) / n_eff

    # Print human-readable table
    print(f"{'Parameter':<14} {'Coverage':>10}")
    print(f"{'-'*14} {'-'*10}")
    for name in param_names:
        print(f"{name:<14} {coverage[name]:>10.2f}")
    print(f"\nExpected coverage: 0.95.  "
          f"Binomial 95% interval for N={n_eff}: "
          f"[{lo_binom:.3f}, {hi_binom:.3f}].")

    # Save TSV
    cov_path = os.path.join(analyses_dir, "coverage_summary.txt")
    with open(cov_path, "w") as f:
        f.write("Parameter\tCoverage\tN\tExpected\tBinom_2.5%\tBinom_97.5%\n")
        for name in param_names:
            f.write(f"{name}\t{coverage[name]:.2f}\t{n_eff}\t"
                    f"0.950\t{lo_binom:.3f}\t{hi_binom:.3f}\n")
    print(f"  Saved {cov_path}")

    # Step 4: RUV — compute and save normalized ranks
    print(f"\nComputing normalized ranks...")
    ranks = compute_normalized_ranks(true_vals, script_dir, included,
                                     burnin_frac, PARAMS)

    ranks_path = os.path.join(analyses_dir, "ranks.tsv")
    with open(ranks_path, "w") as f:
        f.write("replicate\t" + "\t".join(param_names) + "\n")
        for j, i in enumerate(included):
            vals = "\t".join(str(ranks[name][j]) for name in param_names)
            f.write(f"sim_{i:03d}\t{vals}\n")
    print(f"  Saved {ranks_path}")

    # Step 5: Clade coverage
    print(f"\nComputing clade coverage...")
    compute_clade_coverage(script_dir, analyses_dir, included, burnin_pct,
                           args.num_coverage_bins)

    print("\nDone.  Run 04_plot.py to generate plots.")


if __name__ == "__main__":
    main()
