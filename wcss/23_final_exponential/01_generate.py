#!/usr/bin/env python3
"""Generate all simulation inputs and a Makefile for the final exponential WCSS.

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
from scipy.stats import laplace


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

NUM_TIPS = 500
NUM_SITES = 30000
T0 = "2026-01-01"
TIP_DATE_START = date(2025, 1, 1)
TIP_DATE_END = date(2025, 12, 31)

# Inverse-Gamma prior on n0 (in years): mean=2.5, stddev=0.5
# InvGamma: mean = beta/(alpha-1), var = beta^2/((alpha-1)^2*(alpha-2))
# => alpha = 2 + mean^2/var = 27, beta = mean*(alpha-1) = 65.0
N0_PRIOR_MEAN = 2.5
N0_PRIOR_STDDEV = 0.5
N0_PRIOR_ALPHA = 2.0 + N0_PRIOR_MEAN**2 / N0_PRIOR_STDDEV**2  # 27.0
N0_PRIOR_BETA = N0_PRIOR_MEAN * (N0_PRIOR_ALPHA - 1.0)         # 65.0

# Truncated Laplace prior on g (growth rate, year^{-1})
G_PRIOR_MU = 2.0
G_PRIOR_SCALE = 0.2
G_MIN = 0.5

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
P_TIP_DATE_UNCERTAIN_UPTO_YEAR = 0.05

# Delphy run parameters (defaults; --steps overrides STEPS)
DEFAULT_STEPS = 3_000_000_000

# Default master seed for replicate parameter sampling
DEFAULT_MASTER_SEED = 2025

# Path to binaries (relative to the 23_final_exponential directory)
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
    """Sample n0, g, mu, kappa, pi, and alpha from their priors."""
    # n0 ~ InvGamma(alpha, beta): sample 1/n0 ~ Gamma(alpha, 1/beta)
    inv_n0 = rng.gamma(N0_PRIOR_ALPHA, 1.0 / N0_PRIOR_BETA)
    n0 = 1.0 / inv_n0

    # g ~ Truncated Laplace(G_PRIOR_MU, G_PRIOR_SCALE) on [G_MIN, inf)
    u_lo = laplace.cdf(G_MIN, loc=G_PRIOR_MU, scale=G_PRIOR_SCALE)
    u = rng.uniform(u_lo, 1.0)
    g = laplace.ppf(u, loc=G_PRIOR_MU, scale=G_PRIOR_SCALE)

    # mu ~ Gamma(mean=MU_PRIOR_MEAN, stddev=MU_PRIOR_STDDEV)
    mu_alpha = (MU_PRIOR_MEAN / MU_PRIOR_STDDEV) ** 2
    mu_beta = MU_PRIOR_MEAN / MU_PRIOR_STDDEV ** 2
    mu = rng.gamma(mu_alpha, 1.0 / mu_beta)

    pi = rng.dirichlet([1.0, 1.0, 1.0, 1.0])
    log_kappa = rng.normal(KAPPA_MEAN_LOG, KAPPA_SIGMA_LOG)
    kappa = np.exp(log_kappa)

    # alpha ~ Exponential(mean=ALPHA_PRIOR_MEAN)
    alpha = rng.exponential(ALPHA_PRIOR_MEAN)

    return n0, g, mu, kappa, pi, alpha


def run_sapling(i, n0, g, mu, kappa, pi, alpha, replicate_seed, tips_file,
                script_dir):
    """Run Sapling to simulate data for replicate i."""
    sim_dir = os.path.join(script_dir, "sims", f"sim_{i:03d}")

    cmd = [
        os.path.join(script_dir, SAPLING_REL),
        "--tip-file", tips_file,
        "--t0", T0,
        "--exp-pop-n0", str(n0),
        "--exp-pop-g", str(g),
        "--mu", str(mu),
        "--hky-kappa", str(kappa),
        "--hky-pi-A", str(pi[0]),
        "--hky-pi-C", str(pi[1]),
        "--hky-pi-G", str(pi[2]),
        "--hky-pi-T", str(pi[3]),
        "--site-rate-heterogeneity-alpha", str(alpha),
        "--missing-data-mean-num-gaps", str(MISSING_DATA_MEAN_NUM_GAPS),
        "--missing-data-mean-gap-length", str(MISSING_DATA_MEAN_GAP_LENGTH),
        "--missing-data-mean-num-missing-sites",
        str(MISSING_DATA_MEAN_NUM_MISSING_SITES),
        "--p-tip-date-uncertain-upto-month",
        str(P_TIP_DATE_UNCERTAIN_UPTO_MONTH),
        "--p-tip-date-uncertain-upto-year",
        str(P_TIP_DATE_UNCERTAIN_UPTO_YEAR),
        "--num-sites", str(NUM_SITES),
        "--seed", str(replicate_seed),
        "--out-maple", os.path.join(sim_dir, "sim.maple"),
        "--out-info", os.path.join(sim_dir, "sim_info.json"),
        "--out-newick", os.path.join(sim_dir, "sim.nwk"),
        "--out-nexus", os.path.join(sim_dir, "sim.nexus"),
    ]

    print(f"  Replicate {i:03d}: n0={n0:.4f}, g={g:.4f}, mu={mu:.6e}, "
          f"kappa={kappa:.4f}, alpha={alpha:.4f}, "
          f"pi=({pi[0]:.4f}, {pi[1]:.4f}, {pi[2]:.4f}, {pi[3]:.4f})")
    subprocess.run(cmd, check=True)


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
    lines.append(f"\t  --v0-pop-n0-prior-mean {N0_PRIOR_MEAN:g} \\")
    lines.append(f"\t  --v0-pop-n0-prior-stddev {N0_PRIOR_STDDEV:g} \\")
    lines.append(f"\t  --v0-pop-g-prior-mu {G_PRIOR_MU:g} \\")
    lines.append(f"\t  --v0-pop-g-prior-scale {G_PRIOR_SCALE:g} \\")
    lines.append(f"\t  --v0-pop-growth-rate-min {G_MIN:g} \\")
    lines.append(f"\t  --v0-site-rate-heterogeneity \\")
    lines.append(f"\t&& touch $@")
    lines.append("")
    lines.append("clean:")
    lines.append("\trm -f sim_*/delphy.* sim_*/delphy-digested.log "
                 "sim_*/delphy-tips.log sim_*/.done")
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
                        help=f"MCMC steps per replicate "
                             f"(default: {DEFAULT_STEPS})")
    parser.add_argument("--master-seed", type=int, default=DEFAULT_MASTER_SEED,
                        help=f"Master seed for parameter sampling "
                             f"(default: {DEFAULT_MASTER_SEED})")
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
        n0, g, mu, kappa, pi, alpha = sample_prior_params(rng)
        replicate_seed = int(rng.integers(0, 2**31))
        run_sapling(i, n0, g, mu, kappa, pi, alpha, replicate_seed,
                    tips_file, script_dir)
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
