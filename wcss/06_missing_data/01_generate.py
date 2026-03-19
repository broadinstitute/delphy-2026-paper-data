#!/usr/bin/env python3
"""Generate all simulation inputs and a Makefile for the missing-data WCSS.

Usage:
    ./01_generate.py
    ./01_generate.py --n 10 --steps 20000000
"""

import argparse
import os
import subprocess
from datetime import date, timedelta

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

# Missing data parameters
MISSING_DATA_MEAN_NUM_GAPS = 3.0
MISSING_DATA_MEAN_GAP_LENGTH = 500.0
MISSING_DATA_MEAN_NUM_MISSING_SITES = 3.0

# Delphy run parameters (defaults; --steps overrides STEPS)
DEFAULT_STEPS = 1_000_000_000

# Fixed seed for tip date generation (reproducibility)
TIP_SEED = 42

# Default master seed for replicate parameter sampling
DEFAULT_MASTER_SEED = 2025

# Path to binaries (relative to the 06_missing_data directory)
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
        "--missing-data-mean-num-gaps", str(MISSING_DATA_MEAN_NUM_GAPS),
        "--missing-data-mean-gap-length", str(MISSING_DATA_MEAN_GAP_LENGTH),
        "--missing-data-mean-num-missing-sites", str(MISSING_DATA_MEAN_NUM_MISSING_SITES),
        "--num-sites", str(NUM_SITES),
        "--seed", str(replicate_seed),
        "--out-maple", os.path.join(sim_dir, "sim.maple"),
        "--out-info", os.path.join(sim_dir, "sim_info.json"),
        "--out-newick", os.path.join(sim_dir, "sim.nwk"),
        "--out-nexus", os.path.join(sim_dir, "sim.nexus"),
        "--out-fasta", os.path.join(sim_dir, "sim.fasta"),
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

    print("Done!  Next step:")
    print(f"  ./02_run.py   # run Delphy on all {n} replicates")


if __name__ == "__main__":
    main()
