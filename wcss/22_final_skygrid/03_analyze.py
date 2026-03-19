#!/usr/bin/env python3
"""Analyze WCSS results: run loganalyser, check ESS, compute coverage and ranks.

All results are written to TSV files under analyses/.  This script produces
no plots; see 04_plot.py for visualization.

Uses BEAST 2's loganalyser for posterior summary statistics (mean, HPD, ESS).
Raw posterior samples are read only for Rank Uniformity Validation and
tip-date posterior analysis.

Usage:
    ./03_analyze.py
    ./03_analyze.py --n 10 --ignore-low-ess --force-include-all-replicates
"""

import argparse
import concurrent.futures
import io
import json
import math
import os
import subprocess
import sys
from datetime import date

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
    ("alpha",      "alpha"),
    ("tau",        "skygrid.precision"),
    ("gamma_0",    "skygrid.logPopSize5"),   # oldest knot (reversed order)
    ("gamma_1",    "skygrid.logPopSize4"),
    ("gamma_2",    "skygrid.logPopSize3"),
    ("gamma_3",    "skygrid.logPopSize2"),
    ("gamma_4",    "skygrid.logPopSize1"),   # most recent knot
    ("N_bar",      "N_bar"),                 # augmented column
    ("kappa",      "kappa"),
    ("pi_A",       "frequencies1"),
    ("pi_C",       "frequencies2"),
    ("pi_G",       "frequencies3"),
    ("pi_T",       "frequencies4"),
    ("rootHeight", "rootHeight"),
]

# Observables whose ESS should be ignored in the ESS check
ESS_IGNORE = {"skygrid.isloglinear", "skygrid.cutOff"}

ESS_THRESHOLD_LOW = 200
ESS_THRESHOLD_VERY_LOW = 150

# Columns used to compute N_bar in the augmented log file
GAMMA_COLS = [
    "skygrid.logPopSize1", "skygrid.logPopSize2", "skygrid.logPopSize3",
    "skygrid.logPopSize4", "skygrid.logPopSize5",
]


# ---------------------------------------------------------------------------
# Step 0: Digest log files
# ---------------------------------------------------------------------------

def digest_log_file(src_path, dst_path):
    """Produce delphy-digested.log: strip age(...) columns and add N_bar."""
    with open(src_path) as fin, open(dst_path, "w") as fout:
        keep_cols = None
        gamma_indices = None
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if keep_cols is None:
                # Header line: identify columns to keep (strip age(...))
                keep_cols = [i for i, f in enumerate(fields)
                             if not f.startswith("age(")]
                kept_headers = [fields[i] for i in keep_cols]
                gamma_indices = [kept_headers.index(c) for c in GAMMA_COLS]
                fout.write("\t".join(kept_headers) + "\tN_bar\n")
                continue
            if not line.strip():
                continue
            kept_fields = [fields[i] for i in keep_cols]
            gammas = [float(kept_fields[i]) for i in gamma_indices]
            nbar = math.exp(sum(gammas) / len(gammas))
            fout.write("\t".join(kept_fields) + f"\t{nbar}\n")


# ---------------------------------------------------------------------------
# Step 0b: Tip-date log production
# ---------------------------------------------------------------------------

def produce_tips_log(src_path, dst_path):
    """Produce delphy-tips.log with derived tip calendar years.

    For each age(TIP_XXX|...) column in delphy.log, compute:
        tipYear = rootHeight + age(root) - age(TIP)

    With the fixed Delphy (v1.3.1+), all three log columns use
    calendar-aware conversions, so:
        tipYear = to_linear_year(tip.t)
    """
    with open(src_path) as fin, open(dst_path, "w") as fout:
        header_fields = None
        age_col_indices = None
        root_height_idx = None
        age_root_idx = None
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if header_fields is None:
                header_fields = fields
                age_col_indices = [i for i, f in enumerate(fields)
                                   if f.startswith("age(") and "|" in f]
                root_height_idx = fields.index("rootHeight")
                age_root_idx = fields.index("age(root)")
                # Write header: sample + tipYear columns
                out_cols = [fields[0]]  # sample
                for idx in age_col_indices:
                    # age(TIP_XXX|YYYY-MM) -> tipYear(TIP_XXX|YYYY-MM)
                    out_cols.append("tipYear(" + fields[idx][4:])
                fout.write("\t".join(out_cols) + "\n")
                continue
            if not fields[0].strip():
                continue
            root_height = float(fields[root_height_idx])
            age_root = float(fields[age_root_idx])
            out_fields = [fields[0]]  # sample
            for idx in age_col_indices:
                age_tip = float(fields[idx])
                tip_year = root_height + age_root - age_tip
                out_fields.append(f"{tip_year}")
            fout.write("\t".join(out_fields) + "\n")


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


def check_ess(la_df, analyses_dir, ignore_low_ess=False,
              force_include_all=False):
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
    if force_include_all:
        if error:
            print(f"\n  WARNING: Low ESS detected but keeping all replicates "
                  f"(--force-include-all-replicates).")
        return (True, excluded)
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
    """Read true parameter values from Sapling's info JSON and wcss_true_params."""
    info_path = os.path.join(sim_dir, "sim_info.json")
    with open(info_path) as f:
        info = json.load(f)

    gamma_k = info["pop_model"]["gamma_k"]

    # Read supplementary params (tau) not stored by Sapling
    wcss_path = os.path.join(sim_dir, "wcss_true_params.json")
    with open(wcss_path) as f:
        wcss = json.load(f)

    return {
        "mu": info["subst_model"]["mu"],
        "alpha": info["subst_model"]["site_rate_heterogeneity_alpha"],
        "tau": wcss["tau"],
        "gamma_0": gamma_k[0],    # log-years, no conversion needed
        "gamma_1": gamma_k[1],
        "gamma_2": gamma_k[2],
        "gamma_3": gamma_k[3],
        "gamma_4": gamma_k[4],
        "N_bar": math.exp(sum(gamma_k) / len(gamma_k)),  # years
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
        log_path = os.path.join(script_dir, "sims", f"sim_{i:03d}",
                                "delphy-digested.log")
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
# Step 6: Tip-date posterior analysis
# ---------------------------------------------------------------------------

def date_to_linear_year(d):
    """Convert a date to a fractional year, calendar-aware."""
    y_start = date(d.year, 1, 1)
    y_end = date(d.year + 1, 1, 1)
    days_in_year = (y_end - y_start).days
    return d.year + (d - y_start).days / days_in_year


def parse_complete_maple_dates(maple_path):
    """Parse tip names and dates from a COMPLETE MAPLE file.

    Returns dict mapping tip_name_prefix (e.g. 'TIP_042') to date object.
    """
    tip_dates = {}
    with open(maple_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            name_full = line[1:].strip()
            pipe_idx = name_full.rfind("|")
            if pipe_idx < 0:
                continue
            tip_prefix = name_full[:pipe_idx]
            date_str = name_full[pipe_idx + 1:]
            if len(date_str) == 10:  # YYYY-MM-DD
                tip_dates[tip_prefix] = date.fromisoformat(date_str)
    return tip_dates


def _classify_tip_col(col_name):
    """Classify a tipYear(...) column as 'month' or 'year' uncertainty.

    Returns (tip_prefix, uncertainty_type) or None if unrecognized.
    """
    inner = col_name[8:-1]  # strip "tipYear(" and ")"
    pipe_idx = inner.rfind("|")
    if pipe_idx < 0:
        return None
    tip_prefix = inner[:pipe_idx]
    date_part = inner[pipe_idx + 1:]
    if len(date_part) == 7:  # YYYY-MM
        return (tip_prefix, "month")
    elif len(date_part) == 4:  # YYYY
        return (tip_prefix, "year")
    return None


def _run_one_tip_loganalyser(i, script_dir, burnin_pct):
    """Run loganalyser on a single replicate's delphy-tips.log."""
    tips_log = os.path.join("sims", f"sim_{i:03d}", "delphy-tips.log")
    cmd = [LOGANALYSER, "-oneline", "-burnin", str(burnin_pct), tips_log]
    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=script_dir)
    if result.returncode != 0:
        raise RuntimeError(
            f"loganalyser failed for sim_{i:03d} tips:\n{result.stderr}")
    return (i, result.stdout)


def run_tip_loganalyser(script_dir, n, burnin_pct):
    """Run loganalyser per-replicate on delphy-tips.log files.

    Returns:
        tip_la_results: dict mapping replicate index to list of
            (tip_prefix, uncertainty_type, mean, hpd_lo, hpd_hi, ess)
        min_tip_ess_month: list of min ESS across month-uncertain tips
            per replicate (NaN if no month-uncertain tips)
        min_tip_ess_year: list of min ESS across year-uncertain tips
            per replicate (NaN if no year-uncertain tips)
    """
    tip_la_results = {}
    min_tip_ess_month = [float("nan")] * n
    min_tip_ess_year = [float("nan")] * n
    num_workers = os.cpu_count()
    done = 0

    with concurrent.futures.ProcessPoolExecutor(
            max_workers=num_workers) as executor:
        futures = {
            executor.submit(
                _run_one_tip_loganalyser, i, script_dir, burnin_pct): i
            for i in range(n)
        }
        for future in concurrent.futures.as_completed(futures):
            i = futures[future]
            _, stdout = future.result()
            done += 1
            if done % 50 == 0 or done == n:
                print(f"  {done}/{n} replicates processed")

            # Parse loganalyser output
            lines = stdout.strip().split("\n")
            header_idx = None
            for idx, line in enumerate(lines):
                if line.startswith("sample\t"):
                    header_idx = idx
                    break
            if header_idx is None:
                tip_la_results[i] = []
                continue

            tsv_text = "\n".join(lines[header_idx:])
            la = pd.read_csv(io.StringIO(tsv_text), sep="\t")
            if len(la) == 0:
                tip_la_results[i] = []
                continue

            tips = []
            month_ess_vals = []
            year_ess_vals = []

            # Find all tipYear columns from the loganalyser output
            tip_cols = set()
            for col in la.columns:
                if col.endswith(".mean"):
                    base = col[:-5]
                    if base.startswith("tipYear("):
                        tip_cols.add(base)

            for base in tip_cols:
                classified = _classify_tip_col(base)
                if classified is None:
                    continue
                tip_prefix, utype = classified
                mean = la[f"{base}.mean"].iloc[0]
                hpd_lo = la[f"{base}.95%HPDlo"].iloc[0]
                hpd_hi = la[f"{base}.95%HPDup"].iloc[0]
                ess = la[f"{base}.ESS"].iloc[0]
                tips.append((tip_prefix, utype, mean, hpd_lo, hpd_hi, ess))
                if utype == "month":
                    month_ess_vals.append(ess)
                else:
                    year_ess_vals.append(ess)

            tip_la_results[i] = tips
            if month_ess_vals:
                min_tip_ess_month[i] = min(month_ess_vals)
            if year_ess_vals:
                min_tip_ess_year[i] = min(year_ess_vals)

    return tip_la_results, min_tip_ess_month, min_tip_ess_year


def analyze_tip_dates(script_dir, analyses_dir, included, burnin_frac,
                      tip_la_results):
    """Analyze tip-date posteriors for uncertain tips across all replicates."""
    month_results = []
    year_results = []

    for j, i in enumerate(included):
        sim_dir = os.path.join(script_dir, "sims", f"sim_{i:03d}")
        tips_log_path = os.path.join(sim_dir, "delphy-tips.log")
        complete_maple = os.path.join(sim_dir, "sim-COMPLETE.maple")

        tips = tip_la_results.get(i, [])
        if not tips:
            continue

        # Read true dates from COMPLETE MAPLE
        true_dates = parse_complete_maple_dates(complete_maple)

        # Read raw samples from delphy-tips.log for rank computation
        raw = pd.read_table(tips_log_path, comment="#")
        burnin_rows = math.floor(burnin_frac * len(raw))
        data = raw.iloc[burnin_rows:]

        for tip_prefix, utype, la_mean, la_hpd_lo, la_hpd_hi, _ in tips:
            if tip_prefix not in true_dates:
                print(f"  WARNING: No true date for {tip_prefix} in "
                      f"sim_{i:03d}, skipping")
                continue
            true_date = true_dates[tip_prefix]
            true_calendar_year = date_to_linear_year(true_date)

            # Find the matching tipYear column for rank computation
            matching_cols = [c for c in data.columns
                            if c.startswith(f"tipYear({tip_prefix}|")]
            if not matching_cols:
                continue
            samples = data[matching_cols[0]].values
            rank = np.sum(samples < true_calendar_year)
            normalized_rank = rank / len(samples)

            result = {
                "replicate": f"sim_{i:03d}",
                "tip_name": tip_prefix,
                "true_date": true_calendar_year,
                "posterior_mean": la_mean,
                "hpd_lo": la_hpd_lo,
                "hpd_hi": la_hpd_hi,
                "normalized_rank": normalized_rank,
            }

            if utype == "month":
                month_results.append(result)
            else:
                year_results.append(result)

        if (j + 1) % 50 == 0 or (j + 1) == len(included):
            print(f"  {j+1}/{len(included)} replicates processed")

    # Save results
    for label, results in [("month", month_results), ("year", year_results)]:
        path = os.path.join(analyses_dir, f"tip_date_ranks_{label}.tsv")
        with open(path, "w") as f:
            f.write("replicate\ttip_name\ttrue_date\tposterior_mean\t"
                    "hpd_lo\thpd_hi\tnormalized_rank\n")
            for r in results:
                f.write(f"{r['replicate']}\t{r['tip_name']}\t"
                        f"{r['true_date']:.6f}\t{r['posterior_mean']:.6f}\t"
                        f"{r['hpd_lo']:.6f}\t{r['hpd_hi']:.6f}\t"
                        f"{r['normalized_rank']:.6f}\n")
        print(f"  Saved {path} ({len(results)} tips)")

    # Coverage summary
    cov_path = os.path.join(analyses_dir, "tip_date_coverage_summary.tsv")
    with open(cov_path, "w") as f:
        f.write("uncertainty_type\tn_tips\tn_covered\tcoverage\t"
                "expected\tbinom_2.5%\tbinom_97.5%\n")
        for label, results in [("month", month_results),
                               ("year", year_results)]:
            n_tips = len(results)
            if n_tips == 0:
                f.write(f"{label}\t0\t0\tnan\t0.950\tnan\tnan\n")
                continue
            n_covered = sum(1 for r in results
                           if r["hpd_lo"] <= r["true_date"] <= r["hpd_hi"])
            coverage = n_covered / n_tips
            lo_binom = binom.ppf(0.025, n_tips, 0.95) / n_tips
            hi_binom = binom.ppf(0.975, n_tips, 0.95) / n_tips
            f.write(f"{label}\t{n_tips}\t{n_covered}\t{coverage:.4f}\t"
                    f"0.950\t{lo_binom:.4f}\t{hi_binom:.4f}\n")
    print(f"  Saved {cov_path}")

    # Print summary
    print(f"\nTip-date coverage:")
    print(f"  {'Type':<8} {'N':>6} {'Covered':>8} {'Coverage':>10}")
    print(f"  {'-'*8} {'-'*6} {'-'*8} {'-'*10}")
    for label, results in [("month", month_results),
                           ("year", year_results)]:
        n_tips = len(results)
        if n_tips == 0:
            print(f"  {label:<8} {0:>6} {0:>8} {'N/A':>10}")
            continue
        n_covered = sum(1 for r in results
                       if r["hpd_lo"] <= r["true_date"] <= r["hpd_hi"])
        coverage = n_covered / n_tips
        print(f"  {label:<8} {n_tips:>6} {n_covered:>8} {coverage:>10.4f}")


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
    parser.add_argument("--force-include-all-replicates", action="store_true",
                        help="Do not exclude any replicates regardless of ESS")
    args = parser.parse_args()

    n = args.n
    burnin_pct = args.burnin
    burnin_frac = burnin_pct / 100.0
    script_dir = os.path.dirname(os.path.abspath(__file__))
    analyses_dir = os.path.join(script_dir, "analyses")
    os.makedirs(analyses_dir, exist_ok=True)

    param_names = [name for name, _ in PARAMS]

    # Step 0: Digest log files (strip age columns, add N_bar)
    print("Digesting log files...")
    for i in range(n):
        sim_dir = os.path.join(script_dir, "sims", f"sim_{i:03d}")
        digest_log_file(
            os.path.join(sim_dir, "delphy.log"),
            os.path.join(sim_dir, "delphy-digested.log"))
    print(f"  Digested {n} log files")

    # Step 0b: Produce tip-date log files and run per-replicate loganalyser
    print("\nProducing tip-date log files...")
    for i in range(n):
        sim_dir = os.path.join(script_dir, "sims", f"sim_{i:03d}")
        produce_tips_log(
            os.path.join(sim_dir, "delphy.log"),
            os.path.join(sim_dir, "delphy-tips.log"))
    print(f"  Produced {n} tip-date log files")

    print("Running per-replicate loganalyser on tip-date logs...")
    tip_la_results, min_tip_ess_month, min_tip_ess_year = \
        run_tip_loganalyser(script_dir, n, burnin_pct)

    # Step 1: Run loganalyser and ESS check
    print(f"\nRunning loganalyser on {n} replicates (burnin={burnin_pct}%)...")
    la_df = run_loganalyser(script_dir, analyses_dir, n, burnin_pct)
    print(f"  loganalyser returned {len(la_df)} rows")

    # Append tip-date ESS columns for check_ess
    la_df["minTipESS_month.ESS"] = min_tip_ess_month
    la_df["minTipESS_year.ESS"] = min_tip_ess_year

    ok, excluded = check_ess(la_df, analyses_dir, args.ignore_low_ess,
                             args.force_include_all_replicates)
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

    # Step 6: Tip-date analysis
    print("\nAnalyzing tip-date posteriors...")
    analyze_tip_dates(script_dir, analyses_dir, included, burnin_frac,
                      tip_la_results)

    print("\nDone.  Run 04_plot.py to generate plots.")


if __name__ == "__main__":
    main()
