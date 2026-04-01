#!/usr/bin/env python3
"""Generate all simulation inputs and a Makefile for the final Skygrid WCSS.

Variant of study 22b with no year-long tip-date uncertainty.

Usage:
    ./01_generate.py
    ./01_generate.py --n 10 --steps 50000000
"""

import argparse
import json
import math
import os
import subprocess
from datetime import date, timedelta

import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

NUM_TIPS = 500
NUM_SITES = 30000
T0 = "2026-01-01"
TIP_DATE_START = date(2025, 1, 1)
TIP_DATE_END = date(2025, 12, 31)

# Inverse-Gamma prior on N_bar (in years): mean=0.25, stddev=0.05
# InvGamma: mean = beta/(alpha-1), var = beta^2/((alpha-1)^2*(alpha-2))
# => alpha = 2 + mean^2/var = 27, beta = mean*(alpha-1) = 6.5
NBAR_PRIOR_MEAN = 0.25
NBAR_PRIOR_STDDEV = 0.05
NBAR_PRIOR_ALPHA = 2.0 + NBAR_PRIOR_MEAN**2 / NBAR_PRIOR_STDDEV**2  # 27.0
NBAR_PRIOR_BETA = NBAR_PRIOR_MEAN * (NBAR_PRIOR_ALPHA - 1.0)         # 6.5

# Skygrid configuration
NUM_KNOTS = 5
SKYGRID_FIRST_KNOT_DATE = "2025-01-01"
SKYGRID_LAST_KNOT_DATE = "2026-01-01"

# Gamma prior on GMRF precision tau: mean=1.5, stddev=0.25
TAU_PRIOR_ALPHA = 36.0
TAU_PRIOR_BETA = 24.0

# Gamma prior on mu (in subst/site/year), specified via mean and stddev
MU_PRIOR_MEAN = 1e-3
MU_PRIOR_STDDEV = 1e-4

# Delphy's hardcoded kappa prior: log(kappa) ~ Normal(mean=1.0, std=1.25)
KAPPA_MEAN_LOG = 1.0
KAPPA_SIGMA_LOG = 1.25

# Exponential prior on site-rate heterogeneity alpha
ALPHA_PRIOR_MEAN = 1.0

# Missing data parameters
MISSING_DATA_MEAN_NUM_GAPS = 3.0
MISSING_DATA_MEAN_GAP_LENGTH = 500.0
MISSING_DATA_MEAN_NUM_MISSING_SITES = 3.0

# Tip-date uncertainty parameters
P_TIP_DATE_UNCERTAIN_UPTO_MONTH = 0.15
P_TIP_DATE_UNCERTAIN_UPTO_YEAR = 0.0   # <-- Changed from 0.05 in study 22b

# Delphy run parameters (defaults; --steps overrides STEPS)
DEFAULT_STEPS = 3_000_000_000

# Default master seed for replicate parameter sampling
DEFAULT_MASTER_SEED = 2025

# Path to binaries (relative to the 22c_skygrid_no_year_uncertainty directory)
SAPLING_REL = "../../sapling"
DELPHY_REL = "../../delphy"


# ---------------------------------------------------------------------------
# Step 1: Generate tips and run Sapling for each replicate
# ---------------------------------------------------------------------------

def generate_tips(rng, out_path):
    """Write tips.txt with NUM_TIPS tips, dates uniform over 2025."""
    num_days = (TIP_DATE_END - TIP_DATE_START).days  # 364
    day_offsets = rng.integers(0, num_days + 1, size=NUM_TIPS)
    with open(out_path, "w") as f:
        for i in range(NUM_TIPS):
            tip_date = TIP_DATE_START + timedelta(days=int(day_offsets[i]))
            f.write(f"TIP_{i:03d}|{tip_date.isoformat()}\n")


def sample_prior_params(rng):
    """Sample tau, gamma, mu, kappa, pi, and alpha from their priors."""
    # tau ~ Gamma(alpha, beta)
    tau = rng.gamma(TAU_PRIOR_ALPHA, 1.0 / TAU_PRIOR_BETA)

    # N_bar ~ InvGamma(mean=1, stddev=0.2)
    inv_nbar = rng.gamma(NBAR_PRIOR_ALPHA, 1.0 / NBAR_PRIOR_BETA)
    N_bar = 1.0 / inv_nbar  # years

    # Sample gamma_k offsets from GMRF random walk with precision tau
    gamma = [0.0] * NUM_KNOTS
    for k in range(1, NUM_KNOTS):
        gamma[k] = gamma[k-1] + rng.normal(0, 1.0 / math.sqrt(tau))

    # Shift gamma_k's so that exp(mean(gamma_k)) = N_bar (in years)
    mean_gamma = sum(gamma) / NUM_KNOTS
    target_mean = math.log(N_bar)
    gamma = [g + (target_mean - mean_gamma) for g in gamma]

    # Warn if any gamma_k is very low (approaching default barrier at ~-5.9)
    min_gamma = min(gamma)
    if min_gamma < -5.0:
        print(f"    WARNING: min gamma_k = {min_gamma:.2f} "
              f"(N = {math.exp(min_gamma):.4f} yr), "
              f"approaching low-gamma barrier territory")

    # mu ~ Gamma(mean=MU_PRIOR_MEAN, stddev=MU_PRIOR_STDDEV)
    mu_alpha = (MU_PRIOR_MEAN / MU_PRIOR_STDDEV) ** 2
    mu_beta = MU_PRIOR_MEAN / MU_PRIOR_STDDEV ** 2
    mu = rng.gamma(mu_alpha, 1.0 / mu_beta)

    pi = rng.dirichlet([1.0, 1.0, 1.0, 1.0])
    log_kappa = rng.normal(KAPPA_MEAN_LOG, KAPPA_SIGMA_LOG)
    kappa = np.exp(log_kappa)

    # alpha ~ Exponential(mean=ALPHA_PRIOR_MEAN)
    alpha = rng.exponential(ALPHA_PRIOR_MEAN)

    return tau, gamma, N_bar, mu, kappa, pi, alpha


def run_sapling(i, tau, gamma, N_bar, mu, kappa, pi, alpha, replicate_seed,
                tips_file, script_dir):
    """Run Sapling to simulate data for replicate i."""
    sim_dir = os.path.join(script_dir, "sims", f"sim_{i:03d}")

    gamma_str = ",".join(str(g) for g in gamma)

    cmd = [
        os.path.join(script_dir, SAPLING_REL),
        "--tip-file", tips_file,
        "--t0", T0,
        "--skygrid-first-knot-date", SKYGRID_FIRST_KNOT_DATE,
        "--skygrid-last-knot-date", SKYGRID_LAST_KNOT_DATE,
        "--skygrid-gamma", gamma_str,
        "--skygrid-type", "log-linear",
        "--mu", str(mu),
        "--hky-kappa", str(kappa),
        "--hky-pi-A", str(pi[0]),
        "--hky-pi-C", str(pi[1]),
        "--hky-pi-G", str(pi[2]),
        "--hky-pi-T", str(pi[3]),
        "--site-rate-heterogeneity-alpha", str(alpha),
        "--missing-data-mean-num-gaps", str(MISSING_DATA_MEAN_NUM_GAPS),
        "--missing-data-mean-gap-length", str(MISSING_DATA_MEAN_GAP_LENGTH),
        "--missing-data-mean-num-missing-sites", str(MISSING_DATA_MEAN_NUM_MISSING_SITES),
        "--p-tip-date-uncertain-upto-month", str(P_TIP_DATE_UNCERTAIN_UPTO_MONTH),
        "--p-tip-date-uncertain-upto-year", str(P_TIP_DATE_UNCERTAIN_UPTO_YEAR),
        "--num-sites", str(NUM_SITES),
        "--seed", str(replicate_seed),
        "--out-maple", os.path.join(sim_dir, "sim.maple"),
        "--out-info", os.path.join(sim_dir, "sim_info.json"),
        "--out-newick", os.path.join(sim_dir, "sim.nwk"),
        "--out-nexus", os.path.join(sim_dir, "sim.nexus"),
    ]

    gamma_display = ", ".join(f"{g:.3f}" for g in gamma)
    print(f"  Replicate {i:03d}: tau={tau:.4f}, "
          f"gamma=({gamma_display}), N_bar={N_bar:.4f}, mu={mu:.6e}, "
          f"kappa={kappa:.4f}, alpha={alpha:.4f}, "
          f"pi=({pi[0]:.4f}, {pi[1]:.4f}, {pi[2]:.4f}, {pi[3]:.4f})")
    subprocess.run(cmd, check=True)

    # Write supplementary true params (tau is not stored by Sapling)
    wcss_path = os.path.join(sim_dir, "wcss_true_params.json")
    with open(wcss_path, "w") as f:
        json.dump({"tau": tau}, f)


# ---------------------------------------------------------------------------
# Step 2: Generate Makefile
# ---------------------------------------------------------------------------

def generate_makefile(n, steps, script_dir):
    """Generate a Makefile in sims/ to drive all Delphy runs."""
    makefile_path = os.path.join(script_dir, "sims", "Makefile")
    log_every = steps // 1000
    tree_every = steps // 1000

    # Delphy path is relative to sims/, one level deeper than script_dir
    delphy_from_sims = os.path.join("..", DELPHY_REL)

    lines = []
    lines.append("# Auto-generated by 01_generate.py — do not edit")
    lines.append(f"DELPHY = {delphy_from_sims}")
    lines.append("")
    lines.append("SIM_DIRS := $(wildcard sim_[0-9]*)")
    lines.append("")
    lines.append("all: $(addsuffix /.done,$(SIM_DIRS))")
    lines.append("")
    lines.append("sim_%/.done: sim_%/sim.maple")
    lines.append(f"\t$(DELPHY) \\")
    lines.append(f"\t  --v0-in-maple $< \\")
    lines.append(f"\t  --v0-steps {steps} \\")
    lines.append(f"\t  --v0-out-log-file sim_$*/delphy.log \\")
    lines.append(f"\t  --v0-log-every {log_every} \\")
    lines.append(f"\t  --v0-out-trees-file sim_$*/delphy.trees \\")
    lines.append(f"\t  --v0-tree-every {tree_every} \\")
    lines.append(f"\t  --v0-out-delphy-file sim_$*/delphy.dphy \\")
    lines.append(f"\t  --v0-delphy-snapshot-every {tree_every} \\")
    lines.append(f"\t  --v0-mu-prior-mean {MU_PRIOR_MEAN:g} \\")
    lines.append(f"\t  --v0-mu-prior-stddev {MU_PRIOR_STDDEV:g} \\")
    lines.append(f"\t  --v0-pop-model skygrid \\")
    lines.append(f"\t  --v0-skygrid-num-parameters {NUM_KNOTS} \\")
    lines.append(f"\t  --v0-skygrid-first-knot-date {SKYGRID_FIRST_KNOT_DATE} \\")
    lines.append(f"\t  --v0-skygrid-last-knot-date {SKYGRID_LAST_KNOT_DATE} \\")
    lines.append(f"\t  --v0-skygrid-type log-linear \\")
    lines.append(f"\t  --v0-skygrid-nbar-prior-mean {NBAR_PRIOR_MEAN:g} \\")
    lines.append(f"\t  --v0-skygrid-nbar-prior-stddev {NBAR_PRIOR_STDDEV:g} \\")
    lines.append(f"\t  --v0-skygrid-infer-prior-smoothness \\")
    lines.append(f"\t  --v0-skygrid-tau-prior-alpha {TAU_PRIOR_ALPHA:g} \\")
    lines.append(f"\t  --v0-skygrid-tau-prior-beta {TAU_PRIOR_BETA:g} \\")
    lines.append(f"\t  --v0-skygrid-disable-low-pop-barrier \\")
    lines.append(f"\t  --v0-site-rate-heterogeneity \\")
    lines.append(f"\t&& touch $@")
    lines.append("")
    lines.append("clean:")
    lines.append("\trm -f sim_*/delphy.* sim_*/delphy-digested.log sim_*/delphy-tips.log sim_*/.done")
    lines.append("")
    lines.append(".PHONY: all clean")
    lines.append("")

    with open(makefile_path, "w") as f:
        f.write("\n".join(lines))

    print(f"Wrote {makefile_path}")


# ---------------------------------------------------------------------------
# Step 3: Mutation count statistics and plots
# ---------------------------------------------------------------------------

def collect_mutation_counts(script_dir, n):
    """Read num_mutations from each replicate's sim_info.json."""
    records = []
    for i in range(n):
        info_path = os.path.join(script_dir, "sims", f"sim_{i:03d}",
                                 "sim_info.json")
        with open(info_path) as f:
            info = json.load(f)
        records.append({
            "replicate": f"sim_{i:03d}",
            "num_mutations": info["tree_stats"]["num_mutations"],
        })
    return records


def save_mutation_counts(records, script_dir):
    """Save mutation counts to TSV."""
    tsv_path = os.path.join(script_dir, "sims", "mutation_counts.tsv")
    with open(tsv_path, "w") as f:
        f.write("replicate\tnum_mutations\n")
        for r in records:
            f.write(f"{r['replicate']}\t{r['num_mutations']}\n")
    print(f"  Saved {tsv_path}")


def print_mutation_count_stats(records):
    """Print summary statistics for mutation counts."""
    counts = np.array([r["num_mutations"] for r in records])
    print(f"\nMutation count statistics ({len(counts)} replicates):")
    print(f"  Min:    {np.min(counts)}")
    print(f"  p5:     {np.percentile(counts, 5):.0f}")
    print(f"  Median: {np.median(counts):.0f}")
    print(f"  Mean:   {np.mean(counts):.1f}")
    print(f"  p95:    {np.percentile(counts, 95):.0f}")
    print(f"  Max:    {np.max(counts)}")
    print(f"  Target (2 * {NUM_TIPS} tips): {2 * NUM_TIPS}")


def plot_mutation_counts(records, script_dir):
    """Generate histogram and eCDF of mutation counts."""
    plots_dir = os.path.join(script_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    counts = np.array([r["num_mutations"] for r in records])
    target = 2 * NUM_TIPS

    # Histogram
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(counts, bins=30, edgecolor="black", alpha=0.7)
    ax.axvline(target, color="red", linestyle="--",
               label=f"2 x {NUM_TIPS} tips = {target}")
    ax.set_xlabel("Number of mutations")
    ax.set_ylabel("Count")
    ax.set_title("Mutation counts across replicates")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "mutation_counts_histogram.pdf"))
    plt.close(fig)
    print(f"  Saved {os.path.join(plots_dir, 'mutation_counts_histogram.pdf')}")

    # eCDF
    fig, ax = plt.subplots(figsize=(6, 4))
    sorted_counts = np.sort(counts)
    ecdf = np.arange(1, len(sorted_counts) + 1) / len(sorted_counts)
    ax.step(sorted_counts, ecdf, where="post")
    ax.axvline(target, color="red", linestyle="--",
               label=f"2 x {NUM_TIPS} tips = {target}")
    ax.set_xlabel("Number of mutations")
    ax.set_ylabel("Cumulative probability")
    ax.set_title("Mutation counts eCDF")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "mutation_counts_ecdf.pdf"))
    plt.close(fig)
    print(f"  Saved {os.path.join(plots_dir, 'mutation_counts_ecdf.pdf')}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate WCSS simulation inputs and Makefile")
    parser.add_argument("--n", type=int, default=200,
                        help="Number of simulation replicates (default: 200)")
    parser.add_argument("--steps", type=int, default=DEFAULT_STEPS,
                        help=f"MCMC steps per replicate (default: {DEFAULT_STEPS})")
    parser.add_argument("--master-seed", type=int, default=DEFAULT_MASTER_SEED,
                        help=f"Master seed for parameter sampling (default: {DEFAULT_MASTER_SEED})")
    args = parser.parse_args()

    n = args.n
    steps = args.steps
    master_seed = args.master_seed

    # All paths relative to this script's directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    print(f"Generating WCSS with {n} replicates (master seed: {master_seed})")
    print()

    # Ensure sims/ directory exists
    sims_dir = os.path.join(script_dir, "sims")
    os.makedirs(sims_dir, exist_ok=True)

    # Step 1: generate tips, sample parameters, and run Sapling per replicate
    print("Running Sapling simulations:")
    for i in range(n):
        rng = np.random.default_rng(master_seed + i)
        sim_dir = os.path.join(sims_dir, f"sim_{i:03d}")
        os.makedirs(sim_dir, exist_ok=True)
        tips_file = os.path.join(sim_dir, "tips.txt")
        generate_tips(rng, tips_file)
        tau, gamma, N_bar, mu, kappa, pi, alpha = sample_prior_params(rng)
        replicate_seed = int(rng.integers(0, 2**31))
        run_sapling(i, tau, gamma, N_bar, mu, kappa, pi, alpha,
                    replicate_seed, tips_file, script_dir)
    print()

    # Step 2: Makefile
    generate_makefile(n, steps, script_dir)
    print()

    # Step 3: Mutation count statistics and plots
    print("Collecting mutation count statistics...")
    records = collect_mutation_counts(script_dir, n)
    save_mutation_counts(records, script_dir)
    print_mutation_count_stats(records)
    print()
    print("Generating mutation count plots...")
    plot_mutation_counts(records, script_dir)
    print()

    print("Done!  Next step:")
    print(f"  ./02_run.py   # run Delphy on all {n} replicates")


if __name__ == "__main__":
    main()
