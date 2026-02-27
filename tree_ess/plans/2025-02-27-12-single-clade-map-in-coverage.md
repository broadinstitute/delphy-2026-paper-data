# Use a single growing CladeMap in calc_clade_coverage

## Context

`calc_clade_coverage.rs` currently creates a new `CladeMap` for each
posterior tree and another for the true tree.  Since `CladeMap` entries
are never removed (only inserted via `or_insert`), a single growing
`CladeMap` works just as well and avoids redundant re-insertion of
shared clades.

## Changes to src/bin/calc_clade_coverage.rs

- Create one `CladeMap` before the posterior tree loop.
- Pass `&mut clade_map` to every `catalog_tree_clades` call (posterior
  and true tree).
- Remove the two per-scope `let mut clade_map = CladeMap::new()` calls.

## What is NOT changing

- No changes to `src/clades.rs` or any other file.
- No behavioral changes.

## Verification

`cargo test` — all tests should pass.
