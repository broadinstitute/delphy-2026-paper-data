# Plan: Extract `BurninSpec` into library

## Background

All three binaries (`calc_tree_ess`, `compare_clades`, `calc_clade_coverage`)
contain identical copies of a `BurninSpec` enum that specifies burn-in as
either a fraction or a tree count, along with a `first_sample_idx()` method
that converts it to a starting index.  They also each contain identical
logic to parse the `--burnin-pct` / `--burnin-trees` CLI options into a
`BurninSpec`.

## Goal

Eliminate this duplication by moving `BurninSpec` and the option-parsing
logic into a shared library module.

## New file: `src/burnin.rs`

```rust
use std::cmp;

#[derive(Debug)]
pub enum BurninSpec {
    Fract(f64),
    Trees(usize),
}

impl BurninSpec {
    /// Construct from CLI options.  `pct` and `trees` are mutually
    /// exclusive (enforced by clap's `group`).  If neither is given,
    /// defaults to 10% burn-in.
    pub fn from_options(pct: Option<f64>, trees: Option<usize>) -> Self {
        if let Some(pct) = pct {
            assert!((0.0..=100.0).contains(&pct));
            BurninSpec::Fract(pct / 100.0)
        } else if let Some(trees) = trees {
            BurninSpec::Trees(trees)
        } else {
            BurninSpec::Fract(0.10)
        }
    }

    pub fn first_sample_idx(&self, num_trees: usize) -> usize {
        match *self {
            BurninSpec::Fract(pct) => {
                assert!((0.0..=1.0).contains(&pct));
                (num_trees as f64 * pct).floor() as usize
            }
            BurninSpec::Trees(burnin_trees) => cmp::min(num_trees, burnin_trees),
        }
    }
}
```

## Update `src/lib.rs`

Add `pub mod burnin;`.

## Update each binary

In all three binaries (`calc_tree_ess.rs`, `compare_clades.rs`,
`calc_clade_coverage.rs`):

1. Delete the local `BurninSpec` enum and its `impl` block.
2. Delete the local burn-in parsing block (`let burnin_spec = if let
   Some(pct) = ...`).
3. Add `use tree_ess::burnin::BurninSpec;`.
4. Replace the deleted parsing block with:
   ```rust
   let burnin_spec = BurninSpec::from_options(args.burnin_pct, args.burnin_trees);
   ```

The clap `Args` struct in each binary keeps its `burnin_pct` and
`burnin_trees` fields unchanged — those are per-binary CLI definitions
and can't be shared without pulling in clap as a library dependency.

## Verification

`cargo test` — all existing tests must pass unchanged.
