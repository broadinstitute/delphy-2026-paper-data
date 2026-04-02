#!/usr/bin/env python3
"""Compact single-page WCSS summary plot.

Generates a single PDF with three panels:

  A. Coverage bar chart — 95% HPD coverage for each model parameter and
     tip-date uncertainty type, with per-bar binomial confidence intervals.

  B. Clade coverage calibration — posterior clade support vs fraction of
     clades present in the true tree.

  C. Per-parameter scatter + diagnostics — for each parameter, a scatter
     plot of posterior mean vs true value with 95% HPD bars (blue = covers
     true value, red = does not), plus a true-rank histogram and eCDF to
     check calibration.  A legend cell in the first row explains the
     scatter plot elements.

Handles both skygrid and exponential population models, optional site-rate
heterogeneity (alpha), and tip-date uncertainty (month and/or year).
Gracefully skips tip-date types with no uncertain tips.

Output: plots/wcss_compact_summary.pdf

Usage:
    ./plot_compact_summary.py <study_dir>
    ./plot_compact_summary.py 22b_final_skygrid_lower_n0
"""

import argparse
import math
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Parameter display name -> loganalyser column prefix
# ---------------------------------------------------------------------------

PARAM_COL_MAP = {
    "mu":         "meanRate",
    "alpha":      "alpha",
    "n0":         "exponential.popSize",
    "g":          "exponential.growthRate",
    "tau":        "skygrid.precision",
    "gamma_0":    "skygrid.logPopSize5",
    "gamma_1":    "skygrid.logPopSize4",
    "gamma_2":    "skygrid.logPopSize3",
    "gamma_3":    "skygrid.logPopSize2",
    "gamma_4":    "skygrid.logPopSize1",
    "N_bar":      "N_bar",
    "kappa":      "kappa",
    "pi_A":       "frequencies1",
    "pi_C":       "frequencies2",
    "pi_G":       "frequencies3",
    "pi_T":       "frequencies4",
    "rootHeight": "rootHeight",
}

DISPLAY_LABELS = {
    "mu":         r"$\mu$",
    "alpha":      r"$\alpha$",
    "n0":         r"$n_0$",
    "g":          r"$g$",
    "tau":        r"$\tau$",
    "gamma_0":    r"$\gamma_0$",
    "gamma_1":    r"$\gamma_1$",
    "gamma_2":    r"$\gamma_2$",
    "gamma_3":    r"$\gamma_3$",
    "gamma_4":    r"$\gamma_4$",
    "N_bar":      r"$\bar{N}$",
    "kappa":      r"$\kappa$",
    "pi_A":       r"$\pi_A$",
    "pi_C":       r"$\pi_C$",
    "pi_G":       r"$\pi_G$",
    "pi_T":       r"$\pi_T$",
    "rootHeight": r"$\mathrm{height}$",
    "td_month":   r"$t_{\mathrm{tip}}\ \mathrm{(month)}$",
    "td_year":    r"$t_{\mathrm{tip}}\ \mathrm{(year)}$",
}

# Row templates: each sub-list is forced onto one row (filtered by available params)
SCATTER_ROW_TEMPLATES = [
    ["mu", "rootHeight", "alpha"],
    ["kappa", "pi_A", "pi_C", "pi_G", "pi_T"],
    ["gamma_0", "gamma_1", "gamma_2", "gamma_3", "gamma_4"],
    ["tau", "N_bar", "n0", "g", "td_month", "td_year"],
]

# Scale factors for display (multiply data by this before plotting)
DISPLAY_SCALE = {
    "mu": 1e3,
}

NCOLS = 5
NBINS_HIST = 10
SPINE_LW = 0.4


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def draw_coverage_bars(ax, coverage_df, td_cov_df, param_names, n):
    """Horizontal bar chart of 95% HPD coverage with per-bar CI markers."""
    names = []
    coverages = []
    sample_sizes = []
    cov_lookup = dict(zip(coverage_df["Parameter"], coverage_df["Coverage"]))
    for name in param_names:   # Ensure names are added in the order in `param_names`, not `coverage_df`
        if name in cov_lookup:
            names.append(name)
            coverages.append(cov_lookup[name])
            sample_sizes.append(n)

    # Add tip-date coverages (skip types with no tips)
    if td_cov_df is not None:
        for _, row in td_cov_df.iterrows():
            n_tips = int(row.get("n_tips", 0))
            if n_tips == 0 or np.isnan(row["coverage"]):
                continue
            names.append(f"td_{row['uncertainty_type']}")
            coverages.append(row["coverage"])
            sample_sizes.append(n_tips)

    coverages = np.array(coverages)
    sample_sizes = np.array(sample_sizes)

    # Space bars evenly across a fixed range so the chart fills the axes
    # regardless of number of parameters
    n_bars = len(names)
    max_bars = 17  # reference: study 22b
    spacing = max_bars / n_bars
    y_pos = np.arange(n_bars) * spacing
    bar_height = 0.7 * spacing

    display_names = [DISPLAY_LABELS.get(name, name) for name in names]

    ax.barh(y_pos, coverages, height=bar_height, color="#1f77b4", alpha=0.7,
            edgecolor="black", linewidth=0.3)

    # Per-bar CI markers
    min_ci_lo = 1.0
    for i, ns in enumerate(sample_sizes):
        ci_lo, ci_hi = stats.binom.interval(0.95, ns, 0.95)
        ci_lo /= ns
        ci_hi /= ns
        min_ci_lo = min(min_ci_lo, ci_lo)
        ax.plot([ci_lo, ci_hi], [y_pos[i], y_pos[i]],
                color="red", linewidth=0.8, alpha=0.5, solid_capstyle="round")
        ax.plot(0.95, y_pos[i], "|", color="red", markersize=3, alpha=0.7)

    # Auto x-axis lower bound (must show all bars and all CI markers)
    min_cov = min(coverages.min(), min_ci_lo)
    x_lo = max(0.0, math.floor((min_cov - 0.02) * 50) / 50)
    x_hi = 1.0
    ax.set_xlim(x_lo, x_hi)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(display_names, fontsize=5.5)
    ax.set_xlabel("95% HPD coverage", fontsize=7)
    ax.tick_params(axis="x", labelsize=6)
    ax.invert_yaxis()


def draw_clade_coverage(ax, cc_df):
    """Clade coverage calibration bar chart."""
    midpoints = (cc_df["bin_lo"].values + cc_df["bin_hi"].values) / 2.0 / 100.0
    width = (cc_df["bin_hi"].values[0] - cc_df["bin_lo"].values[0]) / 100.0
    ax.bar(midpoints, cc_df["fraction"].values, width=width * 0.85,
           edgecolor="black", linewidth=0.3, alpha=0.7)
    ax.plot([0, 1], [0, 1], "r--", alpha=0.7, linewidth=0.8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Posterior clade support", fontsize=7)
    ax.set_ylabel("Fraction in true tree", fontsize=7)
    ax.tick_params(labelsize=6)


def draw_rank_histogram(ax, ranks, ylim_max=2.0):
    """Draw a normalized rank histogram (counts/expected) on ax."""
    counts, edges = np.histogram(ranks, bins=NBINS_HIST, range=(0, 1))
    expected = len(ranks) / NBINS_HIST
    normalized = counts / expected
    midpoints = (edges[:-1] + edges[1:]) / 2.0
    width = 1.0 / NBINS_HIST

    ax.bar(midpoints, normalized, width=width * 0.85, color="#1f77b4",
           alpha=0.6, edgecolor="none")
    ax.axhline(1.0, color="red", linewidth=0.5, alpha=0.7,
               linestyle="--")
    n = len(ranks)
    ci_lo, ci_hi = stats.binom.interval(0.95, n, 1.0 / NBINS_HIST)
    ax.axhspan(ci_lo / expected, ci_hi / expected,
               color="red", alpha=0.08)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, ylim_max)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_LW)
        spine.set_color("grey")


def draw_ecdf(ax, ranks, n):
    """Draw an ECDF on ax with red confidence band, no tick labels."""
    x_grid = np.linspace(0, 1, 200)
    band = 1.96 * np.sqrt(x_grid * (1 - x_grid) / n)
    ax.fill_between(x_grid, x_grid - band, x_grid + band,
                    color="red", alpha=0.12)
    ax.plot([0, 1], [0, 1], "k--", linewidth=0.4, alpha=0.4)

    sorted_ranks = np.sort(ranks)
    ecdf = np.arange(1, len(sorted_ranks) + 1) / len(sorted_ranks)
    ax.step(sorted_ranks, ecdf, where="post", color="#1f77b4",
            linewidth=0.5)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_LW)
        spine.set_color("grey")


def draw_scatter(ax, true_arr, mean_arr, lo_arr, hi_arr, display_name,
                 point_size=2, lw=0.8):
    """Scatter plot with HPD bars and internal label (upper left)."""
    covers = (lo_arr <= true_arr) & (true_arr <= hi_arr)

    for i in range(len(true_arr)):
        color = "#a6cee3" if covers[i] else "#e31a1c"
        ax.plot([true_arr[i], true_arr[i]],
                [lo_arr[i], hi_arr[i]],
                color=color, linewidth=lw, alpha=0.4)
    ax.scatter(true_arr, mean_arr, s=point_size, color="black", zorder=5,
               linewidths=0)

    all_vals = np.concatenate([true_arr, lo_arr, hi_arr])
    lo_bound = np.nanmin(all_vals)
    hi_bound = np.nanmax(all_vals)
    margin = 0.05 * (hi_bound - lo_bound) if hi_bound > lo_bound else 0.1
    ax.plot([lo_bound - margin, hi_bound + margin],
            [lo_bound - margin, hi_bound + margin],
            "r--", alpha=0.4, linewidth=0.5)

    # Parameter name inside the plot (upper left, pushed down)
    ax.text(0.03, 0.90, display_name, transform=ax.transAxes,
            fontsize=6, ha="left", va="top")

    ax.tick_params(labelsize=5, pad=1)
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_LW)
        spine.set_color("grey")


def draw_legend_cell(ax_hist, ax_ecdf, ax_scatter):
    """Draw a key/legend explaining the scatter plot and diagnostic elements."""
    
    # --- Mini-plot labels (histogram and ECDF) ---
    ax_hist.set_xlim(0, 1)
    ax_hist.set_ylim(0, 1)
    ax_hist.set_xticks([])
    ax_hist.set_yticks([])
    for spine in ax_hist.spines.values():
        spine.set_linewidth(SPINE_LW)
        spine.set_color("grey")
    ax_hist.text(0.5, 0.5, "True\nRank\nHistogram", transform=ax_hist.transAxes,
                 fontsize=4.5, ha="center", va="center", color="grey")

    ax_ecdf.set_xlim(0, 1)
    ax_ecdf.set_ylim(0, 1)
    ax_ecdf.set_xticks([])
    ax_ecdf.set_yticks([])
    for spine in ax_ecdf.spines.values():
        spine.set_linewidth(SPINE_LW)
        spine.set_color("grey")
    ax_ecdf.text(0.5, 0.5, "True\nRank\neCDF", transform=ax_ecdf.transAxes,
                 fontsize=4.5, ha="center", va="center", color="grey")

    # --- Scatter key ---
    ax_scatter.set_xlim(0, 1)
    ax_scatter.set_ylim(0, 1)
    ax_scatter.set_xticks([])
    ax_scatter.set_yticks([])
    ax_scatter.set_xlabel("True value", fontsize=5.5, labelpad=2)
    ax_scatter.set_ylabel("Posterior\ndistribution", fontsize=5.5,
                          labelpad=2)
    for spine in ax_scatter.spines.values():
        spine.set_linewidth(SPINE_LW)
        spine.set_color("grey")

    # Diagonal line (no label)
    ax_scatter.plot([0.05, 0.95], [0.05, 0.95], "r--", alpha=0.4,
                    linewidth=0.5)

    # Example 1: blue (HPD covers true value) — shifted right for label room
    x1, y1 = 0.40, 0.74
    y1_lo, y1_hi = 0.30, 0.90
    ax_scatter.plot([x1, x1], [y1_lo, y1_hi],
                    color="#a6cee3", linewidth=1.5, alpha=0.7)
    ax_scatter.scatter([x1], [y1], s=8, color="black", zorder=5,
                       linewidths=0)

    # Example 2: red (HPD does not cover true value) — shifted right
    x2, y2 = 0.80, 0.50
    y2_lo, y2_hi = 0.30, 0.60
    ax_scatter.plot([x2, x2], [y2_lo, y2_hi],
                    color="#e31a1c", linewidth=1.5, alpha=0.7)
    ax_scatter.scatter([x2], [y2], s=8, color="black", zorder=5,
                       linewidths=0)

    # Labels for examples
    ax_scatter.text(x1, y1_lo - 0.04, "True in\n95% HPD",
                    fontsize=4, ha="center", va="top", color="#1f77b4")
    ax_scatter.text(x2, y2_lo - 0.04, "True not in\n95% HPD",
                    fontsize=4, ha="center", va="top", color="#e31a1c")

    # 95% HPD double-ended arrow to the left of the blue bar
    arrow_x = x1 - 0.16
    ax_scatter.arrow(arrow_x, y1_lo, 0.0, y1_hi-y1_lo, color="grey", lw=0.7,
                     head_length = 0.03, head_width = 0.02,
                     length_includes_head = True)
    ax_scatter.arrow(arrow_x, y1_hi, 0.0, y1_lo-y1_hi, color="grey", lw=0.7,
                     head_length = 0.03, head_width = 0.02,
                     length_includes_head = True)
    ax_scatter.text(arrow_x - 0.07, (y1_lo + y1_hi) / 2, "95%\nHPD",
                    fontsize=4, ha="right", va="center", color="grey")

    # Mean label — to the right of the dot
    ax_scatter.annotate("Mean", xy=(x1, y1),
                        xytext=(x1 + 0.05, y1),
                        fontsize=4.5, ha="left", va="center",
                        color="grey")


# ---------------------------------------------------------------------------
# Other helpers
# ---------------------------------------------------------------------------

def read_optional(analyses_dir, filename):
    path = os.path.join(analyses_dir, filename)
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t")
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compact single-page WCSS summary.")
    parser.add_argument("study_dir", help="Path to study directory")
    args = parser.parse_args()

    study_dir = args.study_dir.rstrip("/")
    analyses_dir = os.path.join(study_dir, "analyses")
    plots_dir = os.path.join(study_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # --- Read data ---
    la_df = pd.read_csv(
        os.path.join(analyses_dir, "loganalyser_output_filtered.tsv"),
        sep="\t")
    true_df = pd.read_csv(
        os.path.join(analyses_dir, "true_params.tsv"), sep="\t")
    ranks_df = pd.read_csv(
        os.path.join(analyses_dir, "ranks.tsv"), sep="\t")
    coverage_df = pd.read_csv(
        os.path.join(analyses_dir, "coverage_summary.txt"), sep="\t")
    cc_df = pd.read_csv(
        os.path.join(analyses_dir, "clade_coverage.tsv"), sep="\t")

    td_month_df = read_optional(analyses_dir, "tip_date_ranks_month.tsv")
    td_year_df = read_optional(analyses_dir, "tip_date_ranks_year.tsv")
    td_cov_df = read_optional(analyses_dir, "tip_date_coverage_summary.tsv")

    # Map out observables in the study
    available = {c for c in ranks_df.columns
                 if c != "replicate" and c in PARAM_COL_MAP}
    unknown = {c for c in ranks_df.columns
               if c != "replicate" and c not in PARAM_COL_MAP}
    for name in sorted(unknown):
        print(f"WARNING: ranks column '{name}' not in PARAM_COL_MAP, skipping",
              file=sys.stderr)

    # Tip-date pseudo-params (included in templates)
    if td_month_df is not None and len(td_month_df) > 0:
        available.add("td_month")
    if td_year_df is not None and len(td_year_df) > 0:
        available.add("td_year")

    # Build rows from templates
    scatter_rows = []
    for template in SCATTER_ROW_TEMPLATES:
        row = [p for p in template if p in available]
        if row:
            scatter_rows.append(row)

    param_names = [p for row in scatter_rows for p in row
                   if not p.startswith("td_")]

    n = len(true_df)  # number of replicates in this study
    nrows_scatter = len(scatter_rows)

    # --- Pre-compute histogram ylim ---
    # Normalize all histograms by expected count; find global max ratio
    all_max_ratio = 0.0
    for row_params in scatter_rows:
        for name in row_params:
            if name == "td_month":
                ranks = td_month_df["normalized_rank"].values
            elif name == "td_year":
                ranks = td_year_df["normalized_rank"].values
            else:
                ranks = ranks_df[name].values
            counts, _ = np.histogram(ranks, bins=NBINS_HIST, range=(0, 1))
            expected = len(ranks) / NBINS_HIST  # `len(ranks)` is much larger than `n` for tip-date ranks
            all_max_ratio = max(all_max_ratio, counts.max() / expected)
    hist_ylim = max(2.0, all_max_ratio * 1.1)

    # --- Figure layout ---
    # Fixed top height so clade coverage plot is same size across studies
    top_height = 2.04
    spacer_height = 0.55  # fixed absolute gap between top and scatter panels
    bot_height = 1.4 * nrows_scatter
    fig_height = top_height + spacer_height + bot_height
    fig = plt.figure(figsize=(7.5, fig_height))

    total = top_height + spacer_height + bot_height
    top_frac = top_height / total
    bot_frac = bot_height / total

    # --- Top section: equal horizontal split ---
    gs_top = fig.add_gridspec(
        1, 2, wspace=0.30,
        top=1.0, bottom=1.0 - top_frac,
        left=0.10, right=0.95)
    ax_cov = fig.add_subplot(gs_top[0])
    draw_coverage_bars(ax_cov, coverage_df, td_cov_df, param_names, n)

    ax_cc = fig.add_subplot(gs_top[1])
    draw_clade_coverage(ax_cc, cc_df)

    # --- Bottom section: scatter groups with tight inter-group spacing ---
    gs_bot = fig.add_gridspec(
        nrows_scatter, 1,
        hspace=0.25,
        top=bot_frac, bottom=0.0,
        left=0.08, right=0.95)

    for row_idx, row_params in enumerate(scatter_rows):
        gs_group = gs_bot[row_idx].subgridspec(
            2, NCOLS,
            height_ratios=[0.45, 0.75],
            hspace=0.05, wspace=0.35)

        # Assign column indices sequentially
        col_map = {name: i for i, name in enumerate(row_params)}
        occupied_cols = set(col_map.values())

        for name in row_params:
            col_idx = col_map[name]

            if name == "td_month":
                td_df = td_month_df
                ranks = td_df["normalized_rank"].values
                true_arr = td_df["true_date"].values - 2025.0
                mean_arr = td_df["posterior_mean"].values - 2025.0
                lo_arr = td_df["hpd_lo"].values - 2025.0
                hi_arr = td_df["hpd_hi"].values - 2025.0
                n_ecdf = len(ranks)
                ps, plw = 0.125, 0.2
            elif name == "td_year":
                td_df = td_year_df
                ranks = td_df["normalized_rank"].values
                true_arr = td_df["true_date"].values - 2025.0
                mean_arr = td_df["posterior_mean"].values - 2025.0
                lo_arr = td_df["hpd_lo"].values - 2025.0
                hi_arr = td_df["hpd_hi"].values - 2025.0
                n_ecdf = len(ranks)
                ps, plw = 0.125, 0.2
            else:
                la_col = PARAM_COL_MAP[name]
                scale = DISPLAY_SCALE.get(name, 1.0)
                true_arr = true_df[name].values[:n] * scale
                mean_arr = la_df[f"{la_col}.mean"].values[:n] * scale
                lo_arr = la_df[f"{la_col}.95%HPDlo"].values[:n] * scale
                hi_arr = la_df[f"{la_col}.95%HPDup"].values[:n] * scale
                ranks = ranks_df[name].values
                n_ecdf = n
                ps, plw = 1, 0.6

            display = DISPLAY_LABELS.get(name, name)

            # Diagnostic row
            gs_diag = gs_group[0, col_idx].subgridspec(
                1, 2, wspace=0.06)
            ax_hist = fig.add_subplot(gs_diag[0, 0])
            ax_ecdf = fig.add_subplot(gs_diag[0, 1])

            draw_rank_histogram(ax_hist, ranks, ylim_max=hist_ylim)
            draw_ecdf(ax_ecdf, ranks, n_ecdf)

            # Save first histogram axes for panel label placement
            if row_idx == 0 and col_idx == 0:
                first_hist_ax = ax_hist

            # Scatter row
            ax_scatter = fig.add_subplot(gs_group[1, col_idx])
            draw_scatter(ax_scatter, true_arr, mean_arr, lo_arr, hi_arr,
                         display, point_size=ps, lw=plw)

            # Scale annotation for scaled parameters
            if name in DISPLAY_SCALE:
                exp = int(round(math.log10(DISPLAY_SCALE[name])))
                ax_scatter.text(0.97, 0.03,
                                r"$\times 10^{%d}$" % (-exp),
                                transform=ax_scatter.transAxes,
                                fontsize=5, ha="right", va="bottom",
                                color="grey")

            # Offset annotation for tip-date parameters
            if name in ("td_month", "td_year"):
                ax_scatter.text(0.03, 0.66, r"$+ \, 2025$",
                                transform=ax_scatter.transAxes,
                                fontsize=5, ha="left", va="top",
                                color="grey")

        # Hide unused cells; place legend in the rightmost slot of the first row
        is_first_row = (row_idx == 0)
        legend_col = NCOLS - 1 if is_first_row and NCOLS - 1 not in occupied_cols else None

        for col_idx in range(NCOLS):
            if col_idx in occupied_cols:
                continue
            gs_diag = gs_group[0, col_idx].subgridspec(1, 2, wspace=0.06)
            if col_idx == legend_col:
                ax_lh = fig.add_subplot(gs_diag[0, 0])
                ax_le = fig.add_subplot(gs_diag[0, 1])
                ax_ls = fig.add_subplot(gs_group[1, col_idx])
                draw_legend_cell(ax_lh, ax_le, ax_ls)
            else:
                fig.add_subplot(gs_diag[0, 0]).set_visible(False)
                fig.add_subplot(gs_diag[0, 1]).set_visible(False)
                fig.add_subplot(gs_group[1, col_idx]).set_visible(False)

    # --- Panel labels (placed in figure coordinates for alignment) ---
    label_props = dict(fontsize=10, fontweight="bold", ha="left", va="top")
    label_pad = 0.065  # gap between label and axes left edge (figure coords)

    pos_cov = ax_cov.get_position()
    pos_cc = ax_cc.get_position()
    pos_hist0 = first_hist_ax.get_position()

    label_x_ac = pos_cov.x0 - label_pad   # shared x for A and C
    label_y_ab = pos_cov.y1               # shared y for A and B

    fig.text(label_x_ac, label_y_ab, "A", **label_props)
    fig.text(pos_cc.x0 - label_pad, label_y_ab, "B", **label_props)
    fig.text(label_x_ac, pos_hist0.y1, "C", **label_props)

    # --- Save ---
    out_path = os.path.join(plots_dir, "wcss_compact_summary.pdf")
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
