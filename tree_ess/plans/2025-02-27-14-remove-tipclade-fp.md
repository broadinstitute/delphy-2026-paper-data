# Remove redundant `fp` field from `TipClade`

## Context

`CladeDefinition::TipClade { name, fp }` stores `fp`, which is always
equal to the `CladeMap` key under which the entry lives.  This is
redundant — the fingerprint is already known from the map key at every
use site.  Removing it tightens the data model.

The `CladeDefinition::fp()` method currently returns `*fp` for
`TipClade` and `subclade1.union(subclade2)` for `InnerNodeClade`.
After removing `fp` from `TipClade`, the method can no longer derive
the fingerprint for tips.  Since the only callsite (`catalog_tree_clades`)
already knows the fingerprint from the traversal state, `fp()` can be
removed entirely.

Additionally, `catalog_tree_clades` currently tracks a
`node_clade_defn_stack: Vec<CladeDefinition>` to push/pop the parent's
clade definition on enter/exit.  Converting this to a
`work_stack: Vec<(CladeFp, CladeDefinition)>` (mirroring
`calc_clade_times_from_root`'s `work_stack: Vec<(CladeFp, f64)>`)
keeps the current node's fingerprint in a local `fp` variable,
avoiding the need for `CladeDefinition::fp()` and the redundant
`tip_fps` re-lookup.

## Changes to src/clades.rs

### Remove `fp` from `TipClade` (line 65)

```rust
TipClade {
    name: String,
}
```

### Remove `CladeDefinition::fp()` method (lines 75-84)

Delete the entire `impl CladeDefinition` block containing `fp()`.

### Restructure `catalog_tree_clades` (lines 122-194)

Replace `node_clade_defn` + `node_clade_defn_stack` with
`fp` + `clade_def` + `work_stack: Vec<(CladeFp, CladeDefinition)>`.
Also rename `node_clade_defn` to `clade_def` throughout.

```rust
pub fn catalog_tree_clades(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
    clade_map: &mut CladeMap,
    include_trivial: bool,
) -> Vec<CladeFp> {
    let mut fp = CladeFp::empty();
    let mut clade_def = CladeDefinition::InnerNodeClade {
        subclade1: CladeFp::empty(),
        subclade2: CladeFp::empty(),
    };
    let mut work_stack: Vec<(CladeFp, CladeDefinition)> = vec![];

    let mut clade_fps: Vec<CladeFp> = Vec::new();

    for (action, node_ref) in tree.traversal_iter() {
        let node = node_ref.borrow();
        match action {
            TraversalAction::Enter => {
                // Shift focus from parent to node
                work_stack.push((fp, clade_def));

                if node.is_tip() {
                    fp = *tip_fps.get(node.name.as_str()).expect("Unknown tip name");
                    clade_def = CladeDefinition::TipClade {
                        name: node.name.clone(),
                    };
                } else {
                    fp = CladeFp::empty();
                    clade_def = CladeDefinition::InnerNodeClade {
                        subclade1: CladeFp::empty(),
                        subclade2: CladeFp::empty(),
                    };
                }
            }
            TraversalAction::Exit => {
                clade_map.entry(fp).or_insert(clade_def);

                let is_trivial = node.is_tip() || node_ref == tree.root;
                if include_trivial || !is_trivial {
                    clade_fps.push(fp);
                }

                // Shift focus from node to parent, merging this node's fp
                let (parent_fp, parent_clade_def) = work_stack.pop().unwrap();
                clade_def = match parent_clade_def {
                    CladeDefinition::TipClade {} => unreachable!(),
                    CladeDefinition::InnerNodeClade {
                        subclade1,
                        subclade2,
                    } => {
                        if subclade1.is_empty() {
                            CladeDefinition::InnerNodeClade {
                                subclade1: fp,
                                subclade2: CladeFp::empty(),
                            }
                        } else if subclade2.is_empty() {
                            CladeDefinition::InnerNodeClade {
                                subclade1,
                                subclade2: fp,
                            }
                        } else {
                            panic!("Non-bifurcating tree?")
                        }
                    }
                };
                fp = parent_fp.union(&fp);
            }
        }
    }

    clade_fps.sort();
    clade_fps
}
```

Key differences from the old version:
- `node_clade_defn` renamed to `clade_def`.
- `fp` is a local `CladeFp` variable, set directly from `tip_fps` on
  enter (tip) or to `CladeFp::empty()` (inner).  No `fp()` call needed.
- `work_stack` holds `(CladeFp, CladeDefinition)`, paralleling
  `calc_clade_times_from_root`'s `(CladeFp, f64)` work stack.
- On exit, the node fp is just `fp` — no method call or match.
- Parent fp merge: `fp = parent_fp.union(&fp)`.

### Update `tips_in_clade` (line 93)

Remove the now-unnecessary `..` from the `TipClade` match arm:

```rust
Some(CladeDefinition::TipClade { name }) => result.push(name.clone()),
```

### Update `unreachable!()` match (line 168)

Remove the `..`:

```rust
CladeDefinition::TipClade {} => unreachable!(),
```

## Changes to src/bin/compare_clades.rs

### Update test `TipClade` literals (lines 445, 452, 459)

Remove `fp: fp_x` from each `TipClade { name: ..., fp: ... }` literal
in the `test_process_tree` test:

```rust
TipClade {
    name: name_a.to_string(),
}
```

## What is NOT changing

- `CladeFp`, `CladeMap`, `assign_tip_fps`, `calc_rf_dist` — unchanged.
- `catalog_tree_clades` signature — unchanged.
- `calc_clade_coverage.rs` — no `TipClade` references outside of
  `catalog_tree_clades`.

## Verification

`cargo test` — all tests should pass.
