# Extract Frechet correlation ESS computation into a library module

## Motivation

Both `calc_tree_ess.rs` and `compare_clades.rs` implement the Frechet
correlation ESS from Magee et al 2024.  The core algorithm is identical:
given an all-pairs distance matrix, compute Var[tau], E[Delta_s^2] at
each lag, autocorrelation rho_s, autocorrelation time (ACT), and the
Frechet correlation ESS.

The two binaries differ in how they obtain the distance matrix:

- `calc_tree_ess` builds a global cross-chain RF distance matrix, then
  slices out per-chain submatrices.  It also reports intermediate values
  (rho_s, ACT) in its JSON output.
- `compare_clades` has a standalone `compute_ess()` function that builds
  its own RF distance matrix from sorted clade fingerprints and returns
  only the ESS value.

The shared core is: **given a distance matrix, compute rho_s, ACT, and
the Frechet correlation ESS**.

## New file: src/ess.rs

Create a new module `ess` (named broadly in anticipation of adding
other ESS methods) with a function that takes a distance matrix view
and returns a result struct containing the intermediate values:

```rust
use ndarray::{Array1, ArrayView2};

pub struct FrechetEssResult {
    pub rho_s: Vec<f64>,
    pub auto_correlation_time: f64,
    pub frechet_ess: f64,
}

/// Calculate the Frechet correlation ESS from an all-pairs distance matrix.
///
/// Implements the method from Magee et al, "How Trustworthy Is Your Tree?
/// Bayesian Phylogenetic Effective Sample Size Through the Lens of Monte
/// Carlo Error", Bayesian Analysis 19(2), 565--593 (2024).
///
/// We simplify their expressions because Var(tau_t) and Var(tau_{t+s}) are
/// statistically identical (why should the variance at different points in
/// the chain be different when averaged over all possible chains?), and we
/// can thus use all pairwise distances to get a slightly better estimate
/// for them (Eq (15) instead of the equations following Eq (18)).
/// With that simplification, we get,
///
///   rho_s = 1 - (1/2) E[Delta_s^2] / Var[tau],
///
/// where
///
///   Var[tau] = 1/(n (n-1)) sum_{i=1}^n sum_{j=i+1}^n d^2(x_i, x_j)
///   E[Delta_s^2] = 1/(n-s) sum_{i=1}^{n-s} d^2(x_i, x_{i+s})
pub fn calc_frechet_ess(d_ij: &ArrayView2<f64>) -> FrechetEssResult {
    ...
}
```

The Var[tau] computation iterates only over pairs `0 <= i < j < n`
and divides by `n*(n-1)`.

The function handles the edge cases:
- n <= 1: return ESS = n, ACT = 1, empty rho_s
- var_tau == 0 (all distances zero): return ESS = n

## Changes to src/lib.rs

Add `pub mod ess;`.

## Changes to compare_clades.rs

Replace the local `compute_ess()` with a local `calc_tree_ess()` that:
1. Builds the RF distance matrix (this part stays local).
2. Calls `calc_frechet_ess()` from the library.
3. Returns the Frechet ESS value.

A "tree ESS" is a Frechet correlation ESS computed against a
Robinson-Foulds distance matrix, so the name `calc_tree_ess` is
appropriate for the local wrapper.  This also avoids having two
functions named `calc_frechet_ess` (one local, one in the library).

```rust
fn calc_tree_ess(sorted_clade_fps: &[Vec<CladeFp>]) -> f64 {
    ...build rf_dist matrix...
    calc_frechet_ess(&rf_dist.view()).frechet_ess
}
```

The `tree_ess` field in `FileInfo` and the call site remain unchanged
(the field name is part of the JSON output contract).

## Changes to calc_tree_ess.rs

Replace the inline Frechet ESS computation block (the var_tau /
mean_delta_s_2 / rho_s / ACT / ESS section) with a call to the library
function:

```rust
let frechet_ess_result = calc_frechet_ess(&prod_d_ij);
```

Then populate `ChainResults` from the returned struct:

```rust
chain.results = Some(ChainResults {
    ...
    rho_s: frechet_ess_result.rho_s,
    auto_correlation_time: frechet_ess_result.auto_correlation_time,
    effective_sample_size: frechet_ess_result.frechet_ess,
});
```

The explanatory comment about Frechet covariance (referencing Eq (18)
of Magee et al) moves to the doc comment on `calc_frechet_ess`.

## What is NOT changing

- The RF distance matrix construction in both binaries stays local
  (they build it differently).
- The `eprintln!` progress messages stay in the binaries.
- The `tree_ess` field name in `FileInfo` (JSON output contract).
- No behavioral changes.
