#!/usr/bin/env python3
"""Generate all simulation inputs and a Makefile for the site-rate-heterogeneity WCSS.

Usage:
    ./01_generate.py
    ./01_generate.py --n 10 --steps 20000000
"""

import argparse
import json
import os
import subprocess
from datetime import date, timedelta

import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

NUM_TIPS = 200
NUM_SITES = 30000
T0 = "2026-01-01"
TIP_DATE_START = date(2025, 1, 1)
TIP_DATE_END = date(2025, 12, 31)

# Inverse-Gamma prior on n0 (in years): mean=3, stddev=1
# InvGamma: mean = beta/(alpha-1), var = beta^2/((alpha-1)^2*(alpha-2))
# => alpha = 2 + mean^2/var = 11, beta = mean*(alpha-1) = 30
N0_PRIOR_MEAN = 3.0
N0_PRIOR_STDDEV = 1.0
N0_PRIOR_ALPHA = 2.0 + N0_PRIOR_MEAN**2 / N0_PRIOR_STDDEV**2  # 11.0
N0_PRIOR_BETA = N0_PRIOR_MEAN * (N0_PRIOR_ALPHA - 1.0)         # 30.0

# Exponential prior on g (in e-foldings/year), constrained to g >= 0
G_PRIOR_MEAN = 1.0

# Gamma(alpha, beta) prior on mu (in subst/site/year)
# Exponential with mean 1e-3
MU_PRIOR_ALPHA = 1.0
MU_PRIOR_BETA = 1000.0  # rate parameter in years

# Delphy's hardcoded kappa prior: log(kappa) ~ Normal(mean=1.0, std=1.25)
KAPPA_MEAN_LOG = 1.0
KAPPA_SIGMA_LOG = 1.25

# Exponential prior on site-rate heterogeneity alpha
ALPHA_PRIOR_MEAN = 1.0

# Delphy run parameters (defaults; --steps overrides STEPS)
DEFAULT_STEPS = 1_000_000_000

# Fixed seed for tip date generation (reproducibility)
TIP_SEED = 42

# Default master seed for replicate parameter sampling
DEFAULT_MASTER_SEED = 2025

# Path to binaries (relative to the 05_site_rate_heterogeneity directory)
SAPLING_REL = "../../sapling"
DELPHY_REL = "../../delphy"


# ---------------------------------------------------------------------------
# Step 1: Generate tips file
# ---------------------------------------------------------------------------

def generate_tips_file(out_path):
    """Write tips.txt with 200 tips, dates uniform over 2025."""
    rng = np.random.default_rng(TIP_SEED)
    num_days = (TIP_DATE_END - TIP_DATE_START).days  # 364
    day_offsets = rng.integers(0, num_days + 1, size=NUM_TIPS)

    with open(out_path, "w") as f:
        for i in range(NUM_TIPS):
            tip_date = TIP_DATE_START + timedelta(days=int(day_offsets[i]))
            f.write(f"TIP_{i:03d}|{tip_date.isoformat()}\n")

    print(f"Wrote {out_path} ({NUM_TIPS} tips)")


# ---------------------------------------------------------------------------
# Step 2: Run Sapling for each replicate
# ---------------------------------------------------------------------------

def sample_prior_params(rng):
    """Sample n0, g, mu, kappa, pi, and alpha from their priors."""
    # n0 ~ InvGamma(alpha, beta): sample 1/n0 ~ Gamma(alpha, 1/beta)
    inv_n0 = rng.gamma(N0_PRIOR_ALPHA, 1.0 / N0_PRIOR_BETA)
    n0 = 1.0 / inv_n0

    # g ~ Exponential(mean=G_PRIOR_MEAN), constrained to g >= 0
    g = rng.exponential(G_PRIOR_MEAN)

    mu = rng.gamma(MU_PRIOR_ALPHA, 1.0 / MU_PRIOR_BETA)
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
    os.makedirs(sim_dir, exist_ok=True)

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
# Step 3: Generate Makefile
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
    lines.append(f"\t  --v0-mu-prior-alpha {MU_PRIOR_ALPHA:g} \\")
    lines.append(f"\t  --v0-mu-prior-beta {MU_PRIOR_BETA:g} \\")
    lines.append(f"\t  --v0-pop-n0-prior-mean {N0_PRIOR_MEAN:g} \\")
    lines.append(f"\t  --v0-pop-n0-prior-stddev {N0_PRIOR_STDDEV:g} \\")
    lines.append(f"\t  --v0-pop-g-prior-exponential-with-mean {G_PRIOR_MEAN:g} \\")
    lines.append(f"\t  --v0-site-rate-heterogeneity \\")
    lines.append(f"\t&& touch $@")
    lines.append("")
    lines.append("clean:")
    lines.append("\trm -f sim_*/delphy.* sim_*/delphy-digested.log sim_*/.done")
    lines.append("")
    lines.append(".PHONY: all clean")
    lines.append("")

    with open(makefile_path, "w") as f:
        f.write("\n".join(lines))

    print(f"Wrote {makefile_path}")


# ---------------------------------------------------------------------------
# Step 4: Mutation count statistics and plots
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

    # Step 1: tips file
    tips_file = os.path.join(sims_dir, "tips.txt")
    generate_tips_file(tips_file)
    print()

    # Step 2: simulate each replicate
    print("Running Sapling simulations:")
    for i in range(n):
        rng = np.random.default_rng(master_seed + i)
        n0, g, mu, kappa, pi, alpha = sample_prior_params(rng)
        replicate_seed = int(rng.integers(0, 2**31))
        run_sapling(i, n0, g, mu, kappa, pi, alpha, replicate_seed,
                    tips_file, script_dir)
    print()

    # Step 3: Makefile
    generate_makefile(n, steps, script_dir)
    print()

    # Step 4: Mutation count statistics and plots
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
