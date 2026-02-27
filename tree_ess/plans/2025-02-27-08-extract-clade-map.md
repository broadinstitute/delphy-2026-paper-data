# Extract CladeDefinition, CladeMap, and tree clade processing into the library

## Motivation

`compare_clades.rs` defines `CladeDefinition`, `CladeMap`, and
`tips_in_clade`, which are general-purpose clade structures.  Its
`process_tree` function interleaves two concerns: building the clade
structure and computing per-clade dates.  Meanwhile,
`calc_clade_coverage.rs` has `extract_nontrivial_clades`, which
reimplements the same tree traversal just to collect clade fingerprints.

By extracting a library function that performs the tree traversal,
both binaries can share the traversal logic.

## Rename clade_fp.rs to clades.rs

Rename the `clade_fp` module to `clades` to reflect its broader scope
as the home for all clade-related types and algorithms.  Update
`lib.rs` and all import paths accordingly.

## Design

The `CladeMap` is tree-independent: it records only the nesting
structure and sizes of clades, which are the same across all trees
containing a given clade.  A single `CladeMap` is shared across many
trees and only grows (new clades are added, existing ones are kept).

Per-tree information (specifically `time_from_root`) is returned
separately as a `HashMap<CladeFp, CladeInTreeInfo>`.

The library function signature:

```rust
pub struct CladeInTreeInfo {
    pub time_from_root: f64,
}

pub fn analyze_tree_clades(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
    clade_map: &mut CladeMap,
) -> HashMap<CladeFp, CladeInTreeInfo> {
    ...
}
```

This performs the same traversal as the current `process_tree`:
- Builds `CladeDefinition`s and inserts into the shared `clade_map`
  (only if a clade has never been seen before).
- Tracks `time_from_root` (accumulated branch length from root) for
  each clade in this particular tree.
- Returns the per-tree info map covering all clades in the tree
  (tips, inner nodes, and root).

## Changes to src/clades.rs (renamed from clade_fp.rs)

### Move CladeDefinition (unchanged)

Move `CladeDefinition` from `compare_clades.rs` as-is:

```rust
#[derive(Debug, Serialize, Clone, Eq, PartialEq)]
pub enum CladeDefinition {
    TipClade {
        name: String,
        fp: CladeFp,
    },
    InnerNodeClade {
        subclade1: CladeFp,
        subclade2: CladeFp,
        size: usize,
    },
}
```

With its existing methods: `fp()`, `size()`, `is_complete()`.

### Move CladeMap and tips_in_clade

```rust
pub type CladeMap = HashMap<CladeFp, CladeDefinition>;

pub fn tips_in_clade(fp: &CladeFp, clade_map: &CladeMap) -> Vec<String> {
    ...
}
```

### Add CladeInTreeInfo and analyze_tree_clades

```rust
pub struct CladeInTreeInfo {
    pub time_from_root: f64,
}

pub fn analyze_tree_clades(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
    clade_map: &mut CladeMap,
) -> HashMap<CladeFp, CladeInTreeInfo> {
    ...
}
```

## Changes to src/dates.rs

Move `find_root_date` from `compare_clades.rs` to the `dates` module.
Change the return type from `f64` to `Option<f64>`, returning `None`
if no tip has an exact date (instead of panicking):

```rust
pub fn find_root_date(
    tree: &NewickTree,
    exact_tip_dates: &HashMap<String, f64>,
) -> Option<f64> {
    ...
}
```

## Changes to compare_clades.rs

### Remove local definitions

Remove the local `CladeDefinition`, `CladeMap`, `tips_in_clade`,
`find_root_date`, and their section comments.  Import them from the
library.

### Simplify process_tree

Replace the inline traversal with:

1. Call `analyze_tree_clades(tree, tip_fps, clade_map)` to get the
   per-tree info and populate the shared clade_map.
2. Call `find_root_date(tree, exact_tip_dates)` to get the root date
   (unwrap with a panic message).
3. Build the local `clade_fps_2_info` HashMap by iterating the
   per-tree info and computing `date = root_date + info.time_from_root`
   for each clade.
4. Build `sorted_clade_fps` from the per-tree info keys.

The date computation stays entirely in the caller — the library's
`CladeInTreeInfo` provides `time_from_root`, and `process_tree`
converts to absolute dates using the root date.

### Update the test

The `test_process_tree_basic` test continues to check the global
`clade_map` entries.  Since `CladeDefinition` is unchanged, the
expected values are unchanged.

## Changes to calc_clade_coverage.rs

Replace the local `extract_nontrivial_clades` with a thin wrapper
around `analyze_tree_clades`:

```rust
fn extract_nontrivial_clades(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
) -> Vec<CladeFp> {
    let num_tips = tip_fps.len();
    let mut clade_map = CladeMap::new();
    analyze_tree_clades(tree, tip_fps, &mut clade_map);
    clade_map.into_iter()
        .filter_map(|(fp, defn)| match defn {
            CladeDefinition::InnerNodeClade { size, .. }
                if size < num_tips => Some(fp),
            _ => None,
        })
        .collect()
}
```

Filters out `TipClade`s and the root (the `InnerNodeClade` whose
`size == num_tips`).

Remove unused `TraversalAction` import (the traversal now happens
inside `analyze_tree_clades`).

## What is NOT changing

- The `eprintln!` progress messages stay in the binaries.
- `TreeResult`, `CladeAccumulator`, and the rest of the per-file
  processing logic stay local in `compare_clades.rs`.
- `CladeDefinition` is unchanged (no `time_from_root`; that belongs
  in the per-tree `CladeInTreeInfo`).
- No behavioral changes.
