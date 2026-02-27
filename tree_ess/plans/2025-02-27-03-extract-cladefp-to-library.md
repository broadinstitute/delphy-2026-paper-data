# Extract CladeFp into a shared library module

## Motivation

All three binaries define their own local `CladeFp(u64)` newtype with
identical core methods (`empty()`, `random(rng)`, `union()`).  The only
differences are:

- `compare_clades.rs` derives `Serialize` and has an `is_empty()` method.
- `calc_tree_ess.rs` and `calc_clade_coverage.rs` have neither.

This plan extracts `CladeFp` into a shared `src/clade_fp.rs` library
module and updates all three binaries to use it.

## New file: src/clade_fp.rs

Create a new module with the union of all features used across binaries:

```rust
use rand::Rng;
use serde::Serialize;

/// A clade fingerprint represents a set of tips.
/// Singletons are represented directly as a random fingerprint.
/// Larger sets are represented as the XOR of all of their tips' fingerprints.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone, Copy, Serialize)]
pub struct CladeFp(u64);

impl CladeFp {
    pub fn new(value: u64) -> CladeFp {
        CladeFp(value)
    }

    pub fn empty() -> CladeFp {
        CladeFp(0)
    }

    pub fn random(rng: &mut dyn Rng) -> CladeFp {
        CladeFp(rng.next_u64())
    }

    pub fn union(&self, other: &CladeFp) -> CladeFp {
        CladeFp(self.0 ^ other.0)
    }

    pub fn is_empty(&self) -> bool {
        self.0 == 0
    }
}
```

The inner `u64` field is private to enforce encapsulation.  A `new()`
constructor is provided for tests that need to create values with specific
bit patterns (e.g., `CladeFp::new(0b0011)`).  No `value()` accessor is
needed: nothing outside the `CladeFp` impl reads the raw `u64`.

## Changes to src/lib.rs

Add `pub mod clade_fp;`.

## Changes to each binary

In all three binaries:

1. Remove the local `CladeFp` struct and `impl` block.
2. Remove the associated comment block (e.g., `// -- Clade fingerprint ...`).
3. Add `use tree_ess::clade_fp::CladeFp;`.
4. Remove `use rand::Rng;` if it is no longer needed locally (it will
   still be needed in `compare_clades.rs` for `assign_tip_fps`).

Additionally, all test code that constructs `CladeFp` values directly
(e.g., `CladeFp(0b0011)`) must change to use `CladeFp::new(0b0011)`,
since the inner field is now private.

No other code changes are needed: all non-test call sites already use
`CladeFp::empty()`, `CladeFp::random(rng)`, and `.union()`, which work
identically with the library version.

## What is NOT changing

- `assign_tip_fps` stays local in each binary for now (that's a separate
  refactoring step).
- No behavioral changes.
