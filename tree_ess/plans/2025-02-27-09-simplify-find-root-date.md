# Simplify find_root_date to derive dates directly from tip names

## Context

`find_root_date` currently takes an `exact_tip_dates: &HashMap<String, f64>`
parameter, but this map is built by calling `parse_tip_date` on each tip
name — information that's already available in the tree's tip nodes.
By calling `parse_tip_date` internally, `find_root_date` can take just
a `&NewickTree` and the entire `exact_tip_dates` map can be eliminated
from `compare_clades.rs`.

## Changes to src/dates.rs

Simplify the signature:

```rust
pub fn find_root_date(tree: &NewickTree) -> Option<f64>
```

Instead of looking up `exact_tip_dates.get(node.name.as_str())`, call
`parse_tip_date(&node.name)` directly.  If it returns `Some(date)`,
compute `date - root_to_node_dist` and return `Some(root_date)`.

## Changes to src/bin/compare_clades.rs

- Remove the `exact_tip_dates` map and its construction loop (lines 313-321).
- Remove the emptiness check and `eprintln!` for tip date count.
- Remove `exact_tip_dates` from `process_tree` and `process_trees` signatures.
- Call `find_root_date(tree)` instead of `find_root_date(tree, exact_tip_dates)`.
- Update the test: remove the `exact_tip_dates` construction and the
  parameter from `process_tree(...)`.

## What is NOT changing

- `parse_tip_date` remains a public function in `dates.rs`.
- `find_root_date` still returns `Option<f64>`.
- The caller still `.expect()`s the result.

## Verification

`cargo test` — all tests should pass.
