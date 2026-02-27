# Refactor calc_tree_ess to use CladeFp newtype

## Motivation

`calc_tree_ess.rs` represents clade fingerprints as bare `u64` values
throughout: in `assign_tip_fingerprints`, `calc_sorted_splits`,
`calc_rf_dist`, the `Sample` struct, and tests.  The other two binaries
(`compare_clades.rs` and `calc_clade_coverage.rs`) already wrap these in a
`CladeFp(u64)` newtype that provides `empty()`, `random(rng)`, and
`union()` methods.

This refactoring brings `calc_tree_ess.rs` in line with the other binaries
by introducing a local `CladeFp` type.  For now, this duplicates the
`CladeFp` definition, laying the groundwork for a later extraction into a
shared library module.

## Current state in calc_tree_ess.rs

- `assign_tip_fingerprints(tree) -> HashMap<String, u64>` uses
  `rand::random()` (thread-local RNG, non-deterministic)
- `calc_sorted_splits(tree, tip_fps) -> Vec<u64>` uses `0u64` for empty
  and `^=` for XOR union
- `calc_rf_dist(a: &Vec<u64>, b: &Vec<u64>) -> u64` compares bare u64s
- `Sample.sorted_split_fingerprints: Vec<u64>`
- Tests use bare `u64` literals like `0b0001`

## Target state in calc_tree_ess.rs (matching compare_clades/calc_clade_coverage)

```rust
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone, Copy)]
struct CladeFp(u64);
impl CladeFp {
    fn empty() -> CladeFp { CladeFp(0) }
    fn random(rng: &mut dyn Rng) -> CladeFp { CladeFp(rng.next_u64()) }
    fn union(&self, other: &CladeFp) -> CladeFp { CladeFp(self.0 ^ other.0) }
}
```

(No `Serialize` or `is_empty()` needed in this binary.)

## Changes

### 1. Add CladeFp struct

Add the `CladeFp` struct shown above, along with `use rand::Rng;`.

### 2. Rename and update assign_tip_fingerprints

Rename `assign_tip_fingerprints` to `assign_tip_fps` (matching the other
binaries) and change its signature to take an explicit `rng: &mut dyn Rng`
parameter:

```rust
fn assign_tip_fps(tree: &NewickTree, rng: &mut dyn Rng) -> HashMap<String, CladeFp> {
    ...
    result.insert(node.name.clone(), CladeFp::random(rng));
    ...
}
```

At the call site in `main`, use a deterministic `Pcg64Mcg` seeded RNG,
matching the pattern in `compare_clades` and `calc_clade_coverage`:

```rust
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;
...
let mut rng = Pcg64Mcg::seed_from_u64(1);
...
tip_name_2_fingerprint_opt.get_or_insert_with(|| assign_tip_fps(&trees[0].1, &mut rng));
```

### 3. Update calc_sorted_splits

Change signature and body to use `CladeFp`:

```rust
fn calc_sorted_splits(
    tree: &NewickTree,
    tip_name_2_fingerprint: &HashMap<String, CladeFp>,
) -> Vec<CladeFp> {
    let mut partial_fingerprints: Vec<CladeFp> = Vec::new();
    let mut result: Vec<CladeFp> = Vec::new();
    ...
        TraversalAction::Enter => {
            if node_ref.borrow().is_tip() {
                partial_fingerprints.push(
                    *tip_name_2_fingerprint.get(...).unwrap(),
                );
            } else {
                partial_fingerprints.push(CladeFp::empty());
            }
        }
        TraversalAction::Exit => {
            let fingerprint = partial_fingerprints.pop().expect("Improper nesting?");
            ...
            if let Some(parent_fingerprint) = partial_fingerprints.last_mut() {
                *parent_fingerprint = parent_fingerprint.union(&fingerprint);
            }
        }
    ...
}
```

### 4. Update calc_rf_dist

Change parameter types from `&Vec<u64>` to `&[CladeFp]`:

```rust
fn calc_rf_dist(
    sorted_split_fingerprints_a: &[CladeFp],
    sorted_split_fingerprints_b: &[CladeFp],
) -> u64 { ... }
```

The body is unchanged because `CladeFp` derives `Ord`.

### 5. Update Sample struct

```rust
struct Sample {
    global_sample_num: usize,
    state: u64,
    sorted_split_fingerprints: Vec<CladeFp>,
}
```

No further changes needed: the field is populated from `calc_sorted_splits`
and consumed by `calc_rf_dist`, both of which are updated above.

### 6. Update tests

Change bare `u64` literals to `CladeFp(...)`:

```rust
let tip_name_2_fingerprint = HashMap::from([
    (String::from("A"), CladeFp(0b0001)),
    (String::from("B"), CladeFp(0b0010)),
    ...
]);
```

And expected values:

```rust
assert_eq!(
    splits,
    vec![CladeFp(0b0011), CladeFp(0b1100), CladeFp(0b1111)]
);
```

Similarly for `calc_rs_dist_test`:

```rust
let sorted_split_fingerprints_a = vec![CladeFp(1), CladeFp(3), CladeFp(4), CladeFp(7)];
let sorted_split_fingerprints_b = vec![CladeFp(2), CladeFp(3), CladeFp(5), CladeFp(7)];
```

### 7. Update imports

- Remove `rand::random` (no longer used directly).
- Add `use rand::{Rng, SeedableRng};` and `use rand_pcg::Pcg64Mcg;`.

## What is NOT changing

- No new library module is created; `CladeFp` is defined locally in
  `calc_tree_ess.rs` (duplicating the other binaries for now).
- No changes to other files.
