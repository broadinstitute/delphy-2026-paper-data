# Remove the `size` field from `InnerNodeClade`

## Context

The `size` field in `InnerNodeClade` is only used in one place:
`extract_nontrivial_clades` in `calc_clade_coverage.rs` filters out
the root clade via `size < num_tips`.  But the root clade can be
identified more directly: its fingerprint is the XOR of all tip
fingerprints, which is readily available from `tip_fps`.

## Changes to src/clades.rs

Remove `size` from `InnerNodeClade`:

```rust
InnerNodeClade {
    subclade1: CladeFp,
    subclade2: CladeFp,
}
```

Remove the `size()` method from `CladeDefinition`.

Keep `is_complete()` as-is — it doesn't depend on `size`.

Simplify `analyze_tree_clades`: remove all `size` tracking
(`child_clade_size`, `size: 0`, `size: child_clade_size`,
`size: size + child_clade_size`).

## Changes to src/bin/calc_clade_coverage.rs

Replace `extract_nontrivial_clades` with a simpler `extract_clades`
that returns *all* clade fingerprints from the tree (tips, inner
nodes, and root):

```rust
fn extract_clades(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
) -> Vec<CladeFp> {
    let mut clade_map = CladeMap::new();
    analyze_tree_clades(tree, tip_fps, &mut clade_map);
    clade_map.into_keys().collect()
}
```

Move the filtering to the caller.  Compute `root_fp` once in `main()`:

```rust
let root_fp = tip_fps.values().fold(CladeFp::empty(), |acc, fp| acc.union(fp));
```

Then at each call site, filter out tip and root fingerprints:

- Posterior loop (line 159): filter each tree's clades to exclude
  tips and root before counting.
- True tree clades (line 172): same filter.

The filter is: keep `fp` if `!tip_fps.values().any(|t| *t == fp) && fp != root_fp`.
But since `tip_fps` is a `HashMap<String, CladeFp>`, a more efficient
approach is to collect tip fps into a `HashSet<CladeFp>` once, then
filter with `!tip_fp_set.contains(&fp) && fp != root_fp`.

Remove the `CladeDefinition` import — `extract_clades` no longer
matches on variants.

Update tests accordingly — they currently call
`extract_nontrivial_clades` directly and expect pre-filtered results.

## Changes to src/bin/compare_clades.rs

Update the test's expected `InnerNodeClade` literals to remove
the `size` field.

## Verification

`cargo test` — all tests should pass.
