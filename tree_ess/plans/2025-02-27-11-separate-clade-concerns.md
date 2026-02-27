# Separate clade cataloging from clade time computation

## Context

`analyze_tree_clades` interleaves two independent concerns in a single
traversal: (1) building `CladeDefinition`s and populating the shared
`CladeMap`, and (2) computing `time_from_root` for each clade.
`calc_clade_coverage` only needs concern (1) and ignores the returned
per-tree info.  Separating them into two functions makes each simpler
and avoids wasted work when times aren't needed.

## Changes to src/clades.rs

### Replace `analyze_tree_clades` with `catalog_tree_clades`

Strip out all time tracking.  The function populates the shared
`CladeMap` and returns a sorted vector of clade fingerprints found
in this tree.  An `include_trivial` flag controls whether tip
singletons and the root clade are included in the returned list:

```rust
pub fn catalog_tree_clades(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
    clade_map: &mut CladeMap,
    include_trivial: bool,
) -> Vec<CladeFp>
```

When `include_trivial` is false, the root clade (XOR of all tip fps)
and individual tip fps are excluded from the returned vector.  The
`CladeMap` is always fully populated regardless of the flag.

The traversal logic is the same as today minus the `time_from_root`
/ `time_from_root_stack` / `per_tree_info` variables.

### Add `calc_clade_times_from_root`

A new function that computes `time_from_root` for every clade in a
tree.  It does *not* build `CladeDefinition`s — it only needs to
track each node's accumulated XOR fingerprint and branch-length
distance from the root:

```rust
pub fn calc_clade_times_from_root(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
) -> HashMap<CladeFp, f64>
```

The values are `time_from_root` directly as `f64`.  Remove
`CladeInTreeInfo` since it is no longer needed.

The traversal is simpler: maintain a stack of `(CladeFp, f64)` pairs
(node fingerprint, time from root).  On Enter, push a new entry
(tip fp for tips, empty for inner nodes).  On Exit, record
`(fp, time_from_root)` into the result, then XOR the node's fp into
the parent's.

### Remove `is_complete`

Remove the `is_complete()` method and the
`assert!(node_clade_defn.is_complete())` call in
`catalog_tree_clades`.  It only validated that the traversal logic
was correct.

## Changes to src/bin/compare_clades.rs

In `process_tree`, replace:

```rust
let per_tree_info = analyze_tree_clades(tree, tip_fps, clade_map);
```

with two calls:

```rust
let sorted_clade_fps = catalog_tree_clades(tree, tip_fps, clade_map, true);
let clade_times_from_root = calc_clade_times_from_root(tree, tip_fps);
```

`sorted_clade_fps` replaces the locally-built sorted vector.

Convert `clade_times_from_root` (a `HashMap<CladeFp, f64>`) into
`clade_fps_2_dates` (a `HashMap<CladeFp, f64>`) by adding `root_date`
to each time-from-root value.

Remove `CladeInTreeInfo` entirely.  Change `TreeResult` to hold
`clade_fps_2_dates: HashMap<CladeFp, f64>` instead of
`clade_fps_2_info: HashMap<CladeFp, CladeInTreeInfo>`.  Update the
accumulation loop to use `date` directly from the map instead of
`info.date`.

Update imports accordingly.

## Changes to src/bin/calc_clade_coverage.rs

Replace `extract_clades` with a direct call to
`catalog_tree_clades(tree, tip_fps, &mut clade_map, false)`, which
returns a sorted vector of nontrivial clade fps directly.  This
eliminates the need for the local `extract_clades` wrapper, the
`tip_fp_set`, and the `root_fp` computation.

Update imports accordingly.

## What is NOT changing

- `CladeMap`, `CladeDefinition`, `tips_in_clade` are unchanged.
- No behavioral changes.

## Verification

`cargo test` — all tests should pass.
