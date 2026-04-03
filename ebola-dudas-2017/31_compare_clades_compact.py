#!/usr/bin/env python
"""Compact clade comparison plots (support + tMRCAs stacked) for each comparison."""

import gzip
import json
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path


def load_comparison(json_gz_path):
    with gzip.open(json_gz_path, "rt") as f:
        return json.load(f)


def plot_clade_panel(data, label_a, label_b, out_path, target_height=4.0):
    """Create a compact 2-row figure: support scatter (top), dates scatter (bottom)."""
    clades = data["clades"]

    support_a = np.array([c["file_a"]["support"] for c in clades])
    support_b = np.array([c["file_b"]["support"] for c in clades])
    support_a_se = np.array([c["file_a"]["support_se"] for c in clades])
    support_b_se = np.array([c["file_b"]["support_se"] for c in clades])
    clade_sizes = np.array([len(c["tips_in_clade"]) for c in clades])
    num_clades = len(clades)

    max_area = 80
    sizes = max_area * clade_sizes / clade_sizes.max()

    both_mask = (support_a > 0) & (support_b > 0)
    date_a = np.array([c["file_a"]["mean_date"] for c in clades])
    date_b = np.array([c["file_b"]["mean_date"] for c in clades])
    date_a_se = np.array([c["file_a"]["date_se"] for c in clades])
    date_b_se = np.array([c["file_b"]["date_se"] for c in clades])

    # Compute shared date ticks and limits upfront
    lo = min(date_a[both_mask].min(), date_b[both_mask].min())
    hi = max(date_a[both_mask].max(), date_b[both_mask].max())
    margin = (hi - lo) * 0.05
    date_lo, date_hi = lo - margin, hi + margin
    date_ticks = ticker.MaxNLocator(nbins=5).tick_values(date_lo, date_hi)
    date_ticks = [t for t in date_ticks if date_lo <= t <= date_hi]

    # Auto-detect date format: use %.1f if range < 1.5 years, else %.0f
    date_fmt = "%.1f" if (hi - lo) < 1.5 else "%.0f"

    fig_width = target_height * 0.55
    fig, axes = plt.subplots(2, 1, figsize=(fig_width, target_height),
                             gridspec_kw={"hspace": 0.0})

    # -- Support scatter (top) --
    ax = axes[0]
    ax.errorbar(
        support_a, support_b,
        xerr=support_a_se, yerr=support_b_se,
        fmt="none", elinewidth=0.4, capsize=0, ecolor="steelblue", alpha=0.25,
    )
    ax.scatter(
        support_a, support_b, s=sizes,
        alpha=0.4, color="steelblue", edgecolors="none",
    )
    ax.plot([0, 1], [0, 1], "k--", linewidth=0.6, alpha=0.4)
    ax.set_xlabel("")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", direction="in")
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels([])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"])
    ax.set_ylabel("")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=9)
    ax.text(0.05, 0.90, "Support", transform=ax.transAxes,
            fontsize=11, va="top", ha="left")

    # -- Date scatter (bottom) --
    ax = axes[1]
    ax.errorbar(
        date_a[both_mask], date_b[both_mask],
        xerr=date_a_se[both_mask], yerr=date_b_se[both_mask],
        fmt="none", elinewidth=0.4, capsize=0, ecolor="coral", alpha=0.25,
    )
    ax.scatter(
        date_a[both_mask], date_b[both_mask], s=sizes[both_mask],
        alpha=0.4, color="coral", edgecolors="none",
    )
    ax.plot([date_lo, date_hi], [date_lo, date_hi],
            "k--", linewidth=0.6, alpha=0.4)
    ax.set_xlabel(label_a, fontsize=11)
    ax.tick_params(axis="x", direction="in")
    ax.set_ylabel("")
    ax.set_xticks(date_ticks)
    ax.set_xticklabels([])
    ax.set_yticks(date_ticks)
    ax.set_xlim(date_lo, date_hi)
    ax.set_ylim(date_lo, date_hi)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=9)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(date_fmt))
    ax.text(0.05, 0.90, "tMRCAs", transform=ax.transAxes,
            fontsize=11, va="top", ha="left")
    ax.text(0.97, 0.08, f"{num_clades} clades", transform=ax.transAxes,
            fontsize=9, va="bottom", ha="right")

    # Single y-axis label on the right, vertically centered
    fig.text(0.97, 0.5, label_b, fontsize=11, va="center", ha="center",
             rotation=-90)

    fig.subplots_adjust(left=0.12, right=0.90, top=0.93, bottom=0.08)
    fig.savefig(out_path, bbox_inches="tight")
    print(f"Saved {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    here = Path(__file__).parent
    plots_dir = here / "plots"
    plots_dir.mkdir(exist_ok=True)
    analysis_dir = here / "analysis"

    comparisons = [
        ("BEAST X A", "BEAST X B", "clade_comparison_beastX_a_vs_b"),
        ("BEAST X Alpha A", "BEAST X Alpha B", "clade_comparison_beastX_alpha_a_vs_b"),
        ("Delphy A", "Delphy B", "clade_comparison_delphy_a_vs_b"),
        ("Delphy Alpha A", "Delphy Alpha B", "clade_comparison_delphy_alpha_a_vs_b"),
        ("BEAST X A", "Delphy A", "clade_comparison_beastX_a_vs_delphy_a"),
        ("BEAST X Alpha A", "Delphy Alpha A", "clade_comparison_beastX_alpha_a_vs_delphy_alpha_a"),
        ("BEAST X B", "Delphy B", "clade_comparison_beastX_b_vs_delphy_b"),
    ]

    for label_a, label_b, name in comparisons:
        json_gz_path = analysis_dir / f"{name}.json.gz"
        if not json_gz_path.exists():
            print(f"Skipping {name} (no cached data)")
            continue

        print(f"Plotting {name}...")
        data = load_comparison(json_gz_path)
        plot_clade_panel(data, label_a, label_b, plots_dir / f"{name}_compact.pdf")
