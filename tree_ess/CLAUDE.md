# tree_ess — Rust tools for tree topology analysis

## Overview

A Rust crate providing two binary tools for analyzing BEAST/Delphy `.trees`
files (NEXUS format with Newick tree blocks), plus shared library code for
parsing and tree manipulation.

## Binaries

### `calc-tree-ess`

Computes tree topology ESS using the Frechet correlation ESS method
(Magee et al 2024).

```
calc-tree-ess [OPTIONS] <files>...

OPTIONS:
  --burnin-pct <f64>      Burn-in as percentage of trees
  --burnin-trees <usize>  Burn-in as number of trees
  -c, --compact           Omit full distance matrix from output
```

- Input: one or more `.trees` files
- Output: JSON to stdout with per-chain ESS, ACT, autocorrelation values,
  and optionally the full Robinson-Foulds distance matrix

### `compare-clades`

Compares clade support and MRCA dates across two `.trees` files.

```
compare-clades [OPTIONS] <file_a> <file_b>

OPTIONS:
  --burnin-pct <f64>      Burn-in as percentage of trees
  --burnin-trees <usize>  Burn-in as number of trees
  --min-support <f64>     Only output clades with support >= threshold
```

- Input: two `.trees` files
- Output: JSON to stdout with per-clade support frequencies, mean MRCA
  dates, and standard errors (adjusted for tree ESS)

## Library modules

### `nexus_reader` (src/nexus_reader.rs)

Parses BEAST/Delphy NEXUS `.trees` files: taxa block, translate table,
and tree statements.  Returns `Vec<(state_number, NewickTree)>`.
Supports filtering by minimum state and state interval.

### `newick` (src/newick/)

Full Newick parser: lexer (src/newick/lexer.rs) and recursive-descent
parser (src/newick/parser.rs).  Handles quoted labels, branch lengths,
and BEAST-style `[&key=value;...]` attributes on both nodes and branches.
Produces `NewickTree` / `NewickNode` structures.

### `trees` (src/trees.rs)

Generic tree traversal framework via `TreeLike` and `NodeLike` traits.
Provides preorder, postorder, and enter/exit traversal iterators.
Includes `SimpleTree<T>` (N-ary) and `BinaryTree<T>` (optimized binary)
concrete implementations.

### `refs` (src/refs/)

Abstraction over node references.  Default: `Rc<RefCell<T>>` (safe).
With feature `fast_unsafe_light_refs`: `NonNull<T>` (zero overhead,
unsafe).  Both implement the `Pool<T>` trait for alloc/free.

## Key algorithms

### Clade fingerprinting

Each tip is assigned a random 64-bit fingerprint.  Inner clades get the
XOR of their descendant tips' fingerprints.  This gives O(1) clade
identity comparison and O(n) computation per tree.

### Robinson-Foulds distance

Computed by sorted-merge of sorted split fingerprint arrays.  Counts
splits present in one tree but not the other.

### Frechet correlation ESS

From the all-to-all RF distance matrix:
- Var(tau) = mean of d^2 over all pairs
- E[Delta_s^2] = mean of d^2(x_i, x_{i+s}) at lag s
- rho_s = 1 - 0.5 * E[Delta_s^2] / Var(tau)
- ACT = 1 + 2 * sum of rho_s until consecutive sum goes negative
- ESS = n / ACT

### Standard errors (in compare-clades)

- Support SE: sqrt(p(1-p) / ESS), where ESS is tree topology ESS
- Date SE: sd / sqrt(min(ESS, clade_occurrence_count))

## Tip date convention

Tip names are expected in `name|YYYY-MM-DD` format (date after last `|`).
Dates are converted to decimal years accounting for leap years.

## Test data

`testdata/twotrees.trees`: NEXUS file with 2 trees over 4 tips (A, B, C, D),
useful for testing parsing and RF distance (RF distance between them = 2).

## Build

```
cargo build --release
# Binaries in target/release/calc-tree-ess and target/release/compare-clades
```
