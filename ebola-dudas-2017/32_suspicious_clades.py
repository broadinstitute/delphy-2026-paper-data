#!/usr/bin/env python
"""Identify suspicious clades in Delphy A vs B with large support differences."""

import argparse
import gzip
import json
from pathlib import Path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find clades with support differing by more than a threshold.")
    parser.add_argument("--threshold", type=float, default=0.2,
                        help="Min |support_a - support_b| to flag (default: 0.2)")
    args = parser.parse_args()
    threshold = args.threshold
    threshold_str = f"{threshold:.2f}".replace(".", "p")

    here = Path(__file__).parent
    analysis_dir = here / "analysis"

    # Load clade comparison data
    with gzip.open(analysis_dir / "clade_comparison_delphy_a_vs_b.json.gz", "rt") as f:
        data = json.load(f)

    # Find suspicious clades: |support_a - support_b| > threshold
    suspicious_clades = []
    for clade in sorted(data["clades"],
                        key=lambda c: len(c["tips_in_clade"]), reverse=True):
        support_a = clade["file_a"]["support"]
        support_b = clade["file_b"]["support"]
        if abs(support_a - support_b) > threshold:
            suspicious_clades.append(clade)

    # Write compressed JSON
    out_json = analysis_dir / f"suspicious_clades_delphy_a_vs_b_{threshold_str}.json.gz"
    with gzip.open(out_json, "wt") as f:
        json.dump(suspicious_clades, f)
    print(f"Wrote {len(suspicious_clades)} suspicious clades to {out_json}")

    # Filter to clades with >20 tips and produce report
    big_suspicious = [c for c in suspicious_clades if len(c["tips_in_clade"]) > 20]

    # Group by clusters of similar size (gaps > 3 tips separate clusters),
    # then refine by splitting into connected components based on tip overlap
    # (two clades are connected if they share > 50% of the smaller's tips).
    big_suspicious.sort(key=lambda c: len(c["tips_in_clade"]), reverse=True)

    # Step 1: initial size-based clusters
    size_clusters = []  # list of lists of indices into big_suspicious
    for i, c in enumerate(big_suspicious):
        size = len(c["tips_in_clade"])
        if not size_clusters or len(big_suspicious[size_clusters[-1][-1]]["tips_in_clade"]) - size > 3:
            size_clusters.append([i])
        else:
            size_clusters[-1].append(i)

    # Step 2: refine each size cluster by tip overlap
    def overlap_components(indices):
        """Split indices into connected components by tip overlap."""
        tip_sets = {i: set(big_suspicious[i]["tips_in_clade"]) for i in indices}
        # Build adjacency: connected if overlap > 50% of smaller set
        adj = {i: set() for i in indices}
        for i in indices:
            for j in indices:
                if i < j:
                    overlap = len(tip_sets[i] & tip_sets[j])
                    min_size = min(len(tip_sets[i]), len(tip_sets[j]))
                    if overlap > 0.5 * min_size:
                        adj[i].add(j)
                        adj[j].add(i)
        # BFS to find components
        visited = set()
        components = []
        for i in indices:
            if i not in visited:
                comp = []
                queue = [i]
                while queue:
                    node = queue.pop()
                    if node in visited:
                        continue
                    visited.add(node)
                    comp.append(node)
                    queue.extend(adj[node] - visited)
                components.append(comp)
        return components

    refined_clusters = []
    for sc in size_clusters:
        if len(sc) == 1:
            refined_clusters.append(sc)
        else:
            refined_clusters.extend(overlap_components(sc))

    # Build groups dict from refined clusters
    groups = {}
    cluster_keys = []
    for comp in refined_clusters:
        sizes = [len(big_suspicious[i]["tips_in_clade"]) for i in comp]
        lo, hi = min(sizes), max(sizes)
        key = f"{lo}-{hi}" if lo != hi else str(lo)
        # Deduplicate keys if needed (e.g., two separate groups both "44")
        orig_key = key
        suffix = 2
        while key in groups:
            key = f"{orig_key}({suffix})"
            suffix += 1
        cluster_keys.append(key)
        groups[key] = [big_suspicious[i] for i in sorted(comp)]

    # Compute core tips (intersection) and extra tips per group
    group_cores = {}
    group_extras = {}
    for key, members in groups.items():
        tip_sets = [set(c["tips_in_clade"]) for c in members]
        core = tip_sets[0].intersection(*tip_sets[1:]) if len(tip_sets) > 1 else tip_sets[0]
        group_cores[key] = core
        group_extras[key] = {c["clade_fp"]: sorted(set(c["tips_in_clade"]) - core)
                             for c in members}

    # Labels for highlighted clades (top 3 groups at threshold 0.1)
    highlight_labels = {
        # Group 1 (X): sizes 1463-1465
        3056168053120098155:  "X_1",
        615860167271339789:   "X_2",
        4208303076419715594:  "X_3",
        1766726551689611884:  "X_4",
        # Group 2 (Y): sizes 1098-1100
        1133435413859384573:  "Y_1",
        3210674569445890789:  "Y_2",
        11774686227717232978: "Y_3",
        # Group 3 (Z): sizes 752-755
        14832638268288838707: "Z_1",
        15512943217666836143: "Z_2",
        9406187989066844029:  "Z_3",
        8906950833791099648:  "Z_4",
        3340453893401912018:  "Z_5",
    }

    # Write TSV report and print to screen
    report_path = analysis_dir / f"suspicious_clades_delphy_a_vs_b_{threshold_str}_report.tsv"
    header = "group_num\tgroup\tlabel\tfingerprint\tsize\tsupport_a\tsupport_b\tdate_a\tdate_b\tcore_size\textra_count\textra_tips"
    print(header)
    bucket_keys = cluster_keys
    with open(report_path, "w") as f:
        f.write(header + "\n")
        for group_num, bucket in enumerate(bucket_keys, 1):
            core_size = len(group_cores[bucket])
            for c in groups.get(bucket, []):
                fp = c["clade_fp"]
                size = len(c["tips_in_clade"])
                sa = c["file_a"]["support"]
                sb = c["file_b"]["support"]
                da = c["file_a"]["mean_date"]
                db = c["file_b"]["mean_date"]
                da_s = f"{da:.4f}" if da is not None else ""
                db_s = f"{db:.4f}" if db is not None else ""
                lbl = highlight_labels.get(fp, "")
                extras = group_extras[bucket][fp]
                extra_str = ",".join(extras) if extras else ""
                line = f"{group_num}\t{bucket}\t{lbl}\t{fp}\t{size}\t{sa:.4f}\t{sb:.4f}\t{da_s}\t{db_s}\t{core_size}\t{len(extras)}\t{extra_str}"
                print(line)
                f.write(line + "\n")

    print(f"\nWrote {len(big_suspicious)} clades (>20 tips) to {report_path}")
