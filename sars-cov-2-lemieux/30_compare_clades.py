#!/usr/bin/env python
"""Compare clade support and dates between two .trees files using compare-clades."""

import json
import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def run_compare_clades(file_a, file_b, burnin_pct=30, min_support=0.01):
    """Run compare-clades and return parsed JSON."""
    cmd = [
        str(Path(__file__).parent.parent / "tree_ess" / "target" / "release" / "compare-clades"),
        "--burnin-pct", str(burnin_pct),
        "--min-support", str(min_support),
        file_a, file_b,
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"compare-clades failed with exit code {result.returncode}")
    print(result.stderr, end="", file=sys.stderr)
    return json.loads(result.stdout)


def plot_comparison(data, label_a, label_b, out_prefix):
    """Make support and date scatter plots as separate PDFs."""
    clades = data["clades"]
    ess_a = data["file_a"]["tree_ess"]
    ess_b = data["file_b"]["tree_ess"]

    support_a = np.array([c["file_a"]["support"] for c in clades])
    support_b = np.array([c["file_b"]["support"] for c in clades])
    support_a_se = np.array([c["file_a"]["support_se"] for c in clades])
    support_b_se = np.array([c["file_b"]["support_se"] for c in clades])
    clade_sizes = np.array([len(c["tips_in_clade"]) for c in clades])

    # Scale dot area proportional to clade size
    min_area, max_area = 5, 150
    sizes = min_area + (max_area - min_area) * (clade_sizes - 1) / max(clade_sizes.max() - 1, 1)

    # -- Support scatter --
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.errorbar(
        support_a, support_b,
        xerr=support_a_se, yerr=support_b_se,
        fmt="none", elinewidth=0.5, capsize=0, ecolor="steelblue", alpha=0.3,
    )
    ax.scatter(
        support_a, support_b, s=sizes,
        alpha=0.4, color="steelblue", edgecolors="none",
    )
    ax.plot([0, 1], [0, 1], "k--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel(f"{label_a} support")
    ax.set_ylabel(f"{label_b} support")
    ax.set_title(f"Clade posterior support\n(ESS: {label_a}={ess_a:.0f}, {label_b}={ess_b:.0f})")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    fig.tight_layout()
    support_path = f"{out_prefix}_support.pdf"
    fig.savefig(support_path)
    print(f"Saved {support_path}")
    plt.close(fig)

    # -- Date scatter --
    # Filter to clades present in both files (support > 0)
    both_mask = (support_a > 0) & (support_b > 0)
    date_a = np.array([c["file_a"]["mean_date"] for c in clades])
    date_b = np.array([c["file_b"]["mean_date"] for c in clades])
    date_a_se = np.array([c["file_a"]["date_se"] for c in clades])
    date_b_se = np.array([c["file_b"]["date_se"] for c in clades])

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.errorbar(
        date_a[both_mask], date_b[both_mask],
        xerr=date_a_se[both_mask], yerr=date_b_se[both_mask],
        fmt="none", elinewidth=0.5, capsize=0, ecolor="coral", alpha=0.3,
    )
    ax.scatter(
        date_a[both_mask], date_b[both_mask], s=sizes[both_mask],
        alpha=0.4, color="coral", edgecolors="none",
    )
    lo = min(date_a[both_mask].min(), date_b[both_mask].min())
    hi = max(date_a[both_mask].max(), date_b[both_mask].max())
    margin = (hi - lo) * 0.05
    ax.plot([lo - margin, hi + margin], [lo - margin, hi + margin], "k--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel(f"{label_a} mean MRCA date")
    ax.set_ylabel(f"{label_b} mean MRCA date")
    ax.set_title(f"Clade MRCA dates\n({sum(both_mask)} clades present in both)")
    ax.set_aspect("equal")
    fig.tight_layout()
    dates_path = f"{out_prefix}_dates.pdf"
    fig.savefig(dates_path)
    print(f"Saved {dates_path}")
    plt.close(fig)


if __name__ == "__main__":
    here = Path(__file__).parent
    plots_dir = here / "plots"
    plots_dir.mkdir(exist_ok=True)
    analysis_dir = here / "analysis"
    analysis_dir.mkdir(exist_ok=True)

    comparisons = [
        (
            str(here / "beast2_run" / "output.trees"),
            str(here / "beastX_run" / "output.trees"),
            "BEAST 2", "BEAST X",
            "clade_comparison_beast2_vs_beastX",
        ),
        (
            str(here / "beast2_run_alpha" / "output.trees"),
            str(here / "beastX_run_alpha" / "output.trees"),
            "BEAST 2 Alpha", "BEAST X Alpha",
            "clade_comparison_beast2_alpha_vs_beastX_alpha",
        ),
        (
            str(here / "delphy_outputs_a" / "ma_sars_cov_2_delphy.trees"),
            str(here / "delphy_outputs_b" / "ma_sars_cov_2_delphy.trees"),
            "Delphy A", "Delphy B",
            "clade_comparison_delphy_a_vs_b",
        ),
        (
            str(here / "delphy_outputs_alpha_a" / "ma_sars_cov_2_delphy_alpha.trees"),
            str(here / "delphy_outputs_alpha_b" / "ma_sars_cov_2_delphy_alpha.trees"),
            "Delphy Alpha A", "Delphy Alpha B",
            "clade_comparison_delphy_alpha_a_vs_b",
        ),
        (
            str(here / "beastX_run" / "output.trees"),
            str(here / "delphy_outputs_a" / "ma_sars_cov_2_delphy.trees"),
            "BEAST X", "Delphy A",
            "clade_comparison_beastX_vs_delphy_a",
        ),
        (
            str(here / "beastX_run_alpha" / "output.trees"),
            str(here / "delphy_outputs_alpha_a" / "ma_sars_cov_2_delphy_alpha.trees"),
            "BEAST X Alpha", "Delphy Alpha A",
            "clade_comparison_beastX_alpha_vs_delphy_alpha_a",
        ),
    ]

    for file_a, file_b, label_a, label_b, name in comparisons:
        print(f"\n{'='*60}")
        print(f"Comparing {label_a} vs {label_b}")
        print(f"{'='*60}")

        json_path = analysis_dir / f"{name}.json"
        out_prefix = plots_dir / name

        # Run compare-clades (or load cached JSON)
        if json_path.exists():
            print(f"Loading cached {json_path}")
            with open(json_path) as f:
                data = json.load(f)
        else:
            data = run_compare_clades(file_a, file_b)
            with open(json_path, "w") as f:
                json.dump(data, f)
            print(f"Cached to {json_path}")

        plot_comparison(data, label_a, label_b, out_prefix)
