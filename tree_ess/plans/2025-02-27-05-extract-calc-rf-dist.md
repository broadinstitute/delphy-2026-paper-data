# Extract calc_rf_dist into the clade_fp library module

## Motivation

`calc_tree_ess.rs` and `compare_clades.rs` have identical copies of
`calc_rf_dist`, which computes the Robinson-Foulds distance between two
trees represented as sorted slices of `CladeFp`.  Since it operates on
`CladeFp` values, it belongs in `src/clade_fp.rs` alongside `CladeFp`
itself.

## Current state

Both binaries have:

```rust
fn calc_rf_dist(sorted_split_fps_a: &[CladeFp], sorted_split_fps_b: &[CladeFp]) -> u64 {
    use EitherOrBoth::{Both, Left, Right};
    let mut distance = 0;
    for pair in itertools::merge_join_by(
        sorted_split_fps_a,
        sorted_split_fps_b,
        |a, b| a.cmp(b),
    ) {
        match pair {
            Both(_, _) => (),
            Left(_) | Right(_) => distance += 1,
        }
    }
    distance
}
```

## Changes

### 1. Add calc_rf_dist to src/clade_fp.rs

Add the function as a public free function, with `use itertools` and
`use itertools::EitherOrBoth`.

### 2. Update calc_tree_ess.rs

1. Remove the local `calc_rf_dist` function.
2. Add `calc_rf_dist` to the `use tree_ess::clade_fp::{...}` import.
3. Remove `use itertools;` and `use itertools::{EitherOrBoth, ...}` if
   no longer needed locally.  `Itertools` is still used (for
   `collect_vec()`), but `EitherOrBoth` and bare `itertools::` are not.

### 3. Update compare_clades.rs

1. Remove the local `calc_rf_dist` function and its comment
   (`// -- RF distance (from calc_tree_ess.rs) --`).
2. Add `calc_rf_dist` to the `use tree_ess::clade_fp::{...}` import.
3. Remove `use itertools::EitherOrBoth;` if no longer needed locally.
   `Itertools` is still used (for `sorted()` and `collect_vec()`).

### 4. Add a test

Add a test to the `clade_fp::tests` module exercising `calc_rf_dist`:
identical inputs should give distance 0, completely disjoint inputs
should count every element, and partially overlapping inputs should
count only the differences.

## What is NOT changing

- `calc_clade_coverage.rs` does not use RF distance; no changes there.
- No behavioral changes.
