# Extract assign_tip_fps into the clade_fp library module

## Motivation

All three binaries have identical copies of `assign_tip_fps`:

```rust
fn assign_tip_fps(tree: &NewickTree, rng: &mut dyn Rng) -> HashMap<String, CladeFp> {
    let mut result = HashMap::new();
    for node_ref in tree.any_order_iter() {
        let node = node_ref.borrow();
        if node.is_tip() {
            assert!(!node.name.is_empty(), "Tip with no name!");
            result.insert(node.name.clone(), CladeFp::random(rng));
        }
    }
    result
}
```

This function is a natural companion to `CladeFp` and belongs alongside
it in `src/clade_fp.rs`.

## Changes

### 1. Add assign_tip_fps to src/clade_fp.rs

Add the function as a public free function, with the necessary imports
(`HashMap`, `NewickTree`, `NodeLike`, `TreeLike`):

```rust
pub fn assign_tip_fps(tree: &NewickTree, rng: &mut dyn Rng) -> HashMap<String, CladeFp> {
    ...
}
```

### 2. Update each binary

In all three binaries:

1. Remove the local `assign_tip_fps` function and its associated comment.
2. Add `use tree_ess::clade_fp::assign_tip_fps;`.
3. Change `use rand::{Rng, SeedableRng};` to `use rand::SeedableRng;`.
   After this change, no binary uses the `Rng` trait directly — each only
   constructs a `Pcg64Mcg` (needs `SeedableRng`) and passes `&mut rng`
   to `assign_tip_fps` (coercion to `&mut dyn Rng` happens automatically).

## What is NOT changing

- No behavioral changes.
- No other functions are being moved.
