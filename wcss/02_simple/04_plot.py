#!/usr/bin/env python3
"""Plot WCSS results from pre-computed TSV files.

Reads all data from analyses/ and produces plots in plots/.
Performs no computation — see 03_analyze.py for that.

Usage:
    ./04_plot.py
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import ConnectionPatch


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# (display_name, loganalyser column prefix)
PARAMS = [
    ("kappa",      "kappa"),
    ("pi_A",       "frequencies1"),
    ("pi_C",       "frequencies2"),
    ("pi_G",       "frequencies3"),
    ("pi_T",       "frequencies4"),
    ("rootHeight", "rootHeight"),
]


# ---------------------------------------------------------------------------
# Individual plots
# ---------------------------------------------------------------------------

def plot_rank_histogram(ranks, name, n, plots_dir):
    """Plot rank histogram for one parameter."""
    fig, ax = plt.subplots(figsize=(5, 3.5))
    num_bins = max(5, n // 10)
    ax.hist(ranks, bins=num_bins, range=(0, 1), edgecolor="black", alpha=0.7)
    ax.axhline(n / num_bins, color="red", linestyle="--",
               label=f"Expected ({n/num_bins:.1f})")
    ax.set_xlabel("Normalized rank")
    ax.set_ylabel("Count")
    ax.set_title(f"Rank histogram: {name}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, f"rank_histogram_{name}.pdf"))
    plt.close(fig)


def plot_ecdf(ranks, name, plots_dir):
    """Plot ECDF of normalized ranks for one parameter."""
    fig, ax = plt.subplots(figsize=(5, 4))
    sorted_ranks = np.sort(ranks)
    N = len(sorted_ranks)
    ecdf = np.arange(1, N + 1) / N
    ax.step(sorted_ranks, ecdf, where="post", label="ECDF")
    ax.plot([0, 1], [0, 1], "r--", label="Uniform")
    # 95% confidence band
    x_grid = np.linspace(0, 1, 200)
    band = 1.96 * np.sqrt(x_grid * (1 - x_grid) / N)
    ax.fill_between(x_grid, x_grid - band, x_grid + band,
                    color="red", alpha=0.15)
    ax.set_xlabel("Normalized rank")
    ax.set_ylabel("Cumulative probability")
    ax.set_title(f"ECDF: {name}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, f"ecdf_{name}.pdf"))
    plt.close(fig)


def plot_clade_coverage(cc_df, plots_dir):
    """Plot standalone clade coverage bar chart."""
    fig, ax = plt.subplots(figsize=(6, 4.5))
    midpoints = (cc_df["bin_lo"].values + cc_df["bin_hi"].values) / 2.0 / 100.0
    width = (cc_df["bin_hi"].values[0] - cc_df["bin_lo"].values[0]) / 100.0
    ax.bar(midpoints, cc_df["fraction"].values, width=width * 0.85,
           edgecolor="black", alpha=0.7)
    ax.plot([0, 1], [0, 1], "r--", alpha=0.7, label="x = y")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Posterior clade support")
    ax.set_ylabel("Fraction present in true tree")
    ax.set_title("Clade coverage")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "clade_coverage.pdf"))
    plt.close(fig)


def plot_scatter(true_df, la_df, name, col, n, plots_dir):
    """Plot true vs posterior mean with HPD error bars."""
    fig, ax = plt.subplots(figsize=(5, 4.5))

    true_arr = true_df[name].values[:n]
    mean_arr = la_df[f"{col}.mean"].values[:n]
    lo_arr = la_df[f"{col}.95%HPDlo"].values[:n]
    hi_arr = la_df[f"{col}.95%HPDup"].values[:n]
    covers = (lo_arr <= true_arr) & (true_arr <= hi_arr)

    for i in range(n):
        color = "lightblue" if covers[i] else "red"
        ax.plot([true_arr[i], true_arr[i]],
                [lo_arr[i], hi_arr[i]],
                color=color, linewidth=1.5, alpha=0.6)
    ax.scatter(true_arr, mean_arr, s=10, color="black", zorder=5)

    # Identity line
    all_vals = np.concatenate([true_arr, mean_arr, lo_arr, hi_arr])
    lo_bound = np.min(all_vals)
    hi_bound = np.max(all_vals)
    margin = 0.05 * (hi_bound - lo_bound)
    ax.plot([lo_bound - margin, hi_bound + margin],
            [lo_bound - margin, hi_bound + margin],
            "r--", alpha=0.5, label="y = x")

    ax.set_xlabel(f"True {name}")
    ax.set_ylabel(f"Posterior mean {name}")
    ax.set_title(f"True vs inferred: {name}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, f"scatter_{name}.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Summary plot
# ---------------------------------------------------------------------------

def plot_summary(true_df, la_df, ranks_df, coverage_df, cc_df, params, n,
                 plots_dir):
    """Generate combined multi-panel summary figure."""
    num_params = len(params)
    num_rows = num_params + 1  # extra row for clade coverage
    fig, axes = plt.subplots(num_rows, 3, figsize=(14, 3.2 * num_rows))

    # Build coverage lookup
    coverage = dict(zip(coverage_df["Parameter"], coverage_df["Coverage"]))

    for row, (name, col) in enumerate(params):
        # Column 1: Rank histogram
        ax = axes[row, 0]
        num_bins = max(5, n // 10)
        ax.hist(ranks_df[name].values, bins=num_bins, range=(0, 1),
                edgecolor="black", alpha=0.7)
        ax.axhline(n / num_bins, color="red", linestyle="--")
        ax.set_title(f"{name}")
        ax.set_ylabel("Count")

        # Column 2: ECDF
        ax = axes[row, 1]
        sorted_ranks = np.sort(ranks_df[name].values)
        ecdf = np.arange(1, len(sorted_ranks) + 1) / len(sorted_ranks)
        ax.step(sorted_ranks, ecdf, where="post")
        ax.plot([0, 1], [0, 1], "r--")
        # 95% confidence band
        x_grid = np.linspace(0, 1, 200)
        band = 1.96 * np.sqrt(x_grid * (1 - x_grid) / n)
        ax.fill_between(x_grid, x_grid - band, x_grid + band,
                        color="red", alpha=0.15)
        cov = coverage[name]
        ax.set_title(f"Coverage: {cov:.0%}")
        ax.set_ylabel("CDF")

        # Column 3: Scatter with HPD bars (from loganalyser)
        ax = axes[row, 2]
        true_arr = true_df[name].values[:n]
        mean_arr = la_df[f"{col}.mean"].values[:n]
        hpd_lo_arr = la_df[f"{col}.95%HPDlo"].values[:n]
        hpd_hi_arr = la_df[f"{col}.95%HPDup"].values[:n]
        covers = (hpd_lo_arr <= true_arr) & (true_arr <= hpd_hi_arr)

        for i in range(n):
            color = "lightblue" if covers[i] else "red"
            ax.plot([true_arr[i], true_arr[i]],
                    [hpd_lo_arr[i], hpd_hi_arr[i]],
                    color=color, linewidth=1.5, alpha=0.6)
        ax.scatter(true_arr, mean_arr, s=10, color="black", zorder=5)

        all_vals = np.concatenate([true_arr, hpd_lo_arr, hpd_hi_arr])
        lo_bound = np.min(all_vals)
        hi_bound = np.max(all_vals)
        margin = 0.05 * (hi_bound - lo_bound) if hi_bound > lo_bound else 0.1
        ax.plot([lo_bound - margin, hi_bound + margin],
                [lo_bound - margin, hi_bound + margin],
                "r--", alpha=0.5)
        ax.set_title(f"{name}")
        ax.set_ylabel("Posterior mean")

        # Zoomed inset for kappa
        if name == "kappa":
            zoom_max = 5.0
            inset_ax = ax.inset_axes([0.05, 0.55, 0.4, 0.4])
            for i in range(n):
                if true_arr[i] <= zoom_max:
                    color = "lightblue" if covers[i] else "red"
                    inset_ax.plot([true_arr[i], true_arr[i]],
                                 [hpd_lo_arr[i], hpd_hi_arr[i]],
                                 color=color, linewidth=1.0, alpha=0.6)
            mask = true_arr <= zoom_max
            if np.any(mask):
                inset_ax.scatter(true_arr[mask], mean_arr[mask],
                                 s=6, color="black", zorder=5)
            zoomed_hi = hpd_hi_arr[true_arr <= zoom_max]
            y_max = (np.max(zoomed_hi) * 1.05
                     if len(zoomed_hi) > 0 else zoom_max)
            inset_ax.plot([0, zoom_max], [0, zoom_max], "r--", alpha=0.5)
            inset_ax.set_xlim(0, zoom_max)
            inset_ax.set_ylim(0, y_max)
            inset_ax.tick_params(labelsize=6)
            rect, connectors = ax.indicate_inset_zoom(
                inset_ax, edgecolor="gray")
            for c in connectors:
                c.set_visible(False)
            # Manual guide lines
            rect_x0, rect_x1 = inset_ax.get_xlim()
            rect_y0, rect_y1 = inset_ax.get_ylim()
            for (main_pt, inset_pt) in [
                ((rect_x0, rect_y1), (0, 0)),
                ((rect_x1, rect_y1), (1, 0)),
            ]:
                con = ConnectionPatch(
                    xyA=main_pt, coordsA=ax.transData,
                    xyB=inset_pt, coordsB=inset_ax.transAxes,
                    color="gray", linewidth=0.7)
                ax.add_artist(con)

    # Clade coverage row (bottom, spanning all 3 columns)
    cc_row = num_params
    # Hide the two right-hand axes
    axes[cc_row, 1].set_visible(False)
    axes[cc_row, 2].set_visible(False)
    # Make the left axis span all 3 columns
    pos0 = axes[cc_row, 0].get_position()
    pos2 = axes[cc_row, 2].get_position()
    axes[cc_row, 0].set_position([pos0.x0, pos0.y0,
                                  pos2.x1 - pos0.x0, pos0.height])
    ax = axes[cc_row, 0]
    midpoints = (cc_df["bin_lo"].values + cc_df["bin_hi"].values) / 2.0 / 100.0
    width = (cc_df["bin_hi"].values[0] - cc_df["bin_lo"].values[0]) / 100.0
    ax.bar(midpoints, cc_df["fraction"].values, width=width * 0.85,
           edgecolor="black", alpha=0.7)
    ax.plot([0, 1], [0, 1], "r--", alpha=0.7)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Posterior clade support")
    ax.set_ylabel("Fraction present in true tree")
    ax.set_title("Clade coverage")

    # Column headers
    axes[0, 0].text(0.5, 1.15, "Rank histogram",
                    transform=axes[0, 0].transAxes,
                    ha="center", fontweight="bold")
    axes[0, 1].text(0.5, 1.15, "ECDF",
                    transform=axes[0, 1].transAxes,
                    ha="center", fontweight="bold")
    axes[0, 2].text(0.5, 1.15, "True vs inferred",
                    transform=axes[0, 2].transAxes,
                    ha="center", fontweight="bold")

    fig.tight_layout()

    # Re-apply the clade coverage spanning after tight_layout
    pos0 = axes[cc_row, 0].get_position()
    # Get the right edge from a parameter row's column 2
    pos2 = axes[0, 2].get_position()
    axes[cc_row, 0].set_position([pos0.x0, pos0.y0,
                                  pos2.x1 - pos0.x0, pos0.height])

    fig.savefig(os.path.join(plots_dir, "wcss_summary.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    analyses_dir = os.path.join(script_dir, "analyses")
    plots_dir = os.path.join(script_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # Read all inputs from analyses/
    print("Reading analysis results...")
    la_df = pd.read_csv(os.path.join(analyses_dir,
                                     "loganalyser_output_filtered.tsv"),
                        sep="\t")
    true_df = pd.read_csv(os.path.join(analyses_dir, "true_params.tsv"),
                          sep="\t")
    ranks_df = pd.read_csv(os.path.join(analyses_dir, "ranks.tsv"),
                           sep="\t")
    coverage_df = pd.read_csv(os.path.join(analyses_dir,
                                           "coverage_summary.txt"),
                              sep="\t")
    cc_df = pd.read_csv(os.path.join(analyses_dir, "clade_coverage.tsv"),
                         sep="\t")

    n = len(true_df)

    # Individual per-parameter plots
    print("Generating individual plots...")
    for name, col in PARAMS:
        plot_rank_histogram(ranks_df[name].values, name, n, plots_dir)
        plot_ecdf(ranks_df[name].values, name, plots_dir)
        plot_scatter(true_df, la_df, name, col, n, plots_dir)

    # Standalone clade coverage plot
    print("Generating clade coverage plot...")
    plot_clade_coverage(cc_df, plots_dir)

    # Summary multi-panel figure
    print("Generating summary plot...")
    plot_summary(true_df, la_df, ranks_df, coverage_df, cc_df, PARAMS, n,
                 plots_dir)

    print(f"Plots saved to {plots_dir}/")


if __name__ == "__main__":
    main()
