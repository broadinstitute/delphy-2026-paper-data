# Plan: `calc-clade-coverage` tool

## Goal

Add a new binary to the `tree-ess` crate that computes clade coverage
statistics for one replicate of a well-calibrated simulation study.
Given a true tree and a set of posterior trees from inference, it bins
posterior clades by their posterior support and reports how many of those
clades are actually present in the true tree.

The tool is designed to be called once per replicate, appending one TSV
line per call.  A downstream script aggregates these lines across all
replicates and computes the final bar chart.

## Background

See `wcss/plans/01_intro.md`, section "Clade Coverage Validation".  In
brief: the true tree is one sample from the posterior, so asking "is this
clade present in the true tree?" is a Bernoulli trial with success
probability equal to the clade's posterior frequency.  Binning clades by
posterior support and checking the observed success rate against the
expected rate validates that the posterior clade probabilities are
well-calibrated.

References:
- <https://github.com/rbouckaert/DeveloperManual> (section "Clade coverage")
- <https://github.com/christiaanjs/beast-validation/blob/master/src/beastvalidation/experimenter/CladeCoverageCalculator.java>

## CLI interface

```
calc-clade-coverage [OPTIONS] --replicate <string> <true_tree> <posterior_trees>

ARGUMENTS:
  <true_tree>        NEXUS file containing the true (generating) tree
  <posterior_trees>   BEAST/Delphy .trees file with posterior tree samples

OPTIONS:
  --replicate <string>     Replicate identifier (required; included as
                           first column of output for traceability)
  --burnin-pct <f64>       Percentage of posterior trees to discard as
                           burn-in (default: 10)
  --burnin-trees <usize>   Number of posterior trees to discard as burn-in
  --num-bins <usize>       Number of equal-width bins spanning [0, 1]
                           (default: 10, i.e., 10% bins: 0-10%, ..., 90-100%)
  --header                 Prepend a TSV header line before the data line
  -h, --help
  -V, --version
```

The two burn-in options are mutually exclusive (same `group = "burn-in"`
pattern as `calc-tree-ess` and `compare-clades`).

## Input files

### True tree (`<true_tree>`)

A NEXUS file containing exactly one tree.  In the WCSS workflow, this is
the Sapling-generated `.nexus` file (e.g., `sims/sim_000/sim.nexus`).
Parse it with `NexusReader` and assert that it contains exactly one tree.

### Posterior trees (`<posterior_trees>`)

A BEAST/Delphy `.trees` file containing many trees (one per MCMC sample).
Parse with `NexusReader`, apply burn-in, then use the post-burn-in trees
for analysis.

## Processing steps

### Step 1: Parse inputs

1. Parse the true tree file with `NexusReader`.  Assert exactly one tree.
2. Parse the posterior trees file with `NexusReader`.
3. Extract tip names from the first posterior tree.  Assign deterministic
   clade fingerprints using `Pcg64Mcg` seeded from a fixed seed (same
   pattern as `compare-clades`).
4. Verify that the true tree has the same set of tips.

### Step 2: Compute posterior clade support

For each post-burn-in posterior tree:
- Traverse the tree computing clade fingerprints (XOR of tip fingerprints,
  same algorithm as `calc_sorted_splits` in `calc-tree-ess` and
  `process_tree` in `compare-clades`).
- For each clade encountered, increment its count in a
  `HashMap<CladeFp, usize>`.

After processing all posterior trees, each clade's posterior support is
`count / num_post_burnin_trees`.

Note: tip clades (singletons) and the root clade (all tips) always have
support 1.0.  These are not interesting for validation.  Exclude clades
of size 1 (tips) and size N (root) from the analysis.

### Step 3: Merge true tree clades into posterior support map

Traverse the true tree computing clade fingerprints (excluding tips and
root).  For each true clade, look it up in the posterior support HashMap
from Step 2.  If the true clade is not present in the HashMap (i.e., it
was never seen in any posterior tree), insert it with a count of 0.

This ensures the HashMap contains every clade from either the posterior
or the true tree, with accurate counts for all of them.

### Step 4: Bin and count

The number of bins K is given by `--num-bins` (default 10).  The bin width
is `1.0 / K`.  With the default K = 10, the bins are [0%, 10%),
[10%, 20%), ..., [90%, 100%].  The last bin is closed on the right:
[90%, 100%].

Accumulate `totals` and `true_hits` separately:

1. **Totals**: iterate over all clades in the posterior support HashMap.
   For each clade, compute its posterior support and bin index, then
   increment `totals[b]`.

2. **True hits**: iterate over the true tree clades (from Step 3).  For
   each true clade, look up its posterior support in the HashMap (all
   lookups succeed since Step 3 ensured every true clade is present),
   compute its bin index, and increment `true_hits[b]`.

This is efficient: the totals pass touches all M posterior clades (no
lookups), and the true_hits pass does only N lookups (N << M), all of
which succeed.

### Step 5: Output

If `--header` is specified, print the header line first.  Then print one
TSV data line.

The first column is `replicate` (the value of `--replicate`).  Then, for each bin b = 0..K-1, three columns:
- `totals_{LO}_{HI}`: total number of clades with support in [LO%, HI%)
- `true_hits_{LO}_{HI}`: number of those clades present in the true tree
- `frac_{LO}_{HI}`: `true_hits / totals` (or `NaN` if totals is 0)

Where LO and HI are the bin boundaries as integers (e.g., `totals_0_10`,
`true_hits_0_10`, `frac_0_10` for the first bin).

Example header (default 10% bins):
```
replicate	totals_0_10	true_hits_0_10	frac_0_10	totals_10_20	true_hits_10_20	frac_10_20	...	totals_90_100	true_hits_90_100	frac_90_100
```

Example data line:
```
sim_042	342	3	0.008772	58	7	0.120690	...	87	82	0.942529
```

All fractions are printed to 6 decimal places.

## Code organization

### New file: `src/bin/calc_clade_coverage.rs`

Following the pattern of the other binaries:
- `extern crate tree_ess;`
- `clap::Parser` for CLI args
- `BurninSpec` (same as the other binaries — duplicated for now, as the
  other binaries do)
- Reuse library code: `NexusReader`, `AllocPool`, `TreeLike`,
  `TraversalAction`, `NodeLike`
- Clade fingerprinting via tree traversal (same XOR algorithm)
- No JSON output, no serde dependency — just plain TSV to stdout
- Progress messages to stderr

### Cargo.toml addition

```toml
[[bin]]
name = "calc-clade-coverage"
path = "src/bin/calc_clade_coverage.rs"
```

## Tests

### Unit tests (in `calc_clade_coverage.rs`)

1. **Identical trees**: true tree == posterior tree (single posterior
   sample, no burn-in).  Every non-trivial clade should be a true hit.
   All clades land in the 90-100% bin (support = 1.0).  Fraction = 1.0.

2. **Completely different trees**: construct a true tree and a posterior
   tree with no shared non-trivial clades.  All true_hits should be 0.

3. **Bin boundary**: a clade with support exactly 0.1 should land in the
   10-20% bin, not the 0-10% bin.  A clade with support exactly 1.0
   should land in the 90-100% bin (capped).

### Integration test

Add a small `.nexus` true tree and a `.trees` file with a few posterior
samples to `testdata/`, then run the binary and check the output against
expected values.

## Usage in WCSS

The WCSS `03_analyze.py` script will call this tool once per replicate:

```bash
# First replicate: include header
calc-clade-coverage --burnin-pct 30 --header --replicate sim_000 \
    sims/sim_000/sim.nexus sims/sim_000/delphy.trees \
    > analyses/clade_coverage.tsv

# Subsequent replicates: append
for i in $(seq 1 199); do
    calc-clade-coverage --burnin-pct 30 \
        --replicate sim_$(printf '%03d' $i) \
        sims/sim_$(printf '%03d' $i)/sim.nexus \
        sims/sim_$(printf '%03d' $i)/delphy.trees \
        >> analyses/clade_coverage.tsv
done
```

Then a plotting script reads `clade_coverage.tsv`, sums the `totals_*`
and `true_hits_*` columns across all rows, computes the aggregate
fractions, and plots the bar chart with the x=y diagonal.
