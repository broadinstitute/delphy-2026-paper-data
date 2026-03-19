# Plan: WCSS with Missing Data

## Goal

Extend `05_site_rate_heterogeneity` to simulate realistic missing data
in the input sequences.  Real sequencing data contains gaps (e.g., from
amplicon dropout) and isolated missing sites (e.g., from low coverage).
This study validates that Delphy's inference remains well-calibrated
when the input data has missing sites masked as `N`.

The underlying model parameters (tree, mu, alpha, n0, g, kappa, pi) are
drawn from the same priors as in `05_site_rate_heterogeneity`.  The only
difference is that Sapling now also simulates missing data before
producing the MAPLE file that Delphy reads.

## Missing data parameters

Using Sapling's new `--missing-data-*` options (from the MissingData
plan):

- `--missing-data-mean-num-gaps 3` — average of 3 contiguous gaps per
  tip (Poisson-distributed)
- `--missing-data-mean-gap-length 500` — average gap length of 500
  sites (exponentially distributed)
- `--missing-data-mean-num-missing-sites 3` — average of 3 isolated
  missing sites per tip (Poisson-distributed)

With L = 30,000 sites, the expected fraction of missing data per tip is
roughly `(3 * 500 + 3) / 30000 ≈ 5%`, which is a realistic level.

## How it works

When missing data is active, Sapling produces:

- `sim.maple` — masked MAPLE file (missing sites replaced with `N`
  ranges).  **This is what Delphy reads.**
- `sim-COMPLETE.maple` — complete MAPLE file (no missing data).
- `sim.fasta` / `sim-COMPLETE.fasta` — same idea for FASTA.
- `sim_info.json` — includes a `missing_data` section with per-tip gap
  and missing site details.

The true tree and parameter values are unchanged by missing data (Sapling
draws missing data from the RNG after mutation simulation, so the
`-COMPLETE` files at a given seed are identical to a run without missing
data at the same seed).

## Differences from `05_site_rate_heterogeneity`

1. **Sapling is called with missing data options** (3 gaps of mean
   length 500, plus 3 isolated missing sites per tip).
2. **Delphy reads the masked `sim.maple`** (already the default — no
   change to the Makefile rule).

That's it.  All priors, parameters, analysis, and plotting are
identical.

## Configuration

- **Directory:** `wcss/06_missing_data/`
- **Tips:** Same as before (200 tips, dates uniform over 2025)
- **N:** 200 (final), 10 (debug)
- **Steps:** 1,000,000,000 (same as `05_site_rate_heterogeneity`)
- **Missing data:**
  - Mean number of gaps per tip: 3
  - Mean gap length: 500 sites
  - Mean number of isolated missing sites per tip: 3
- **Sampled parameters:** Same as `05_site_rate_heterogeneity`
  (alpha, n0, g, mu, kappa, pi — all drawn from their priors)
- **Fixed parameters:**
  - Genome size L = 30,000 sites

## Deliverables

Files in `wcss/06_missing_data/`:

0. **`00_plan.md`** — Symlink to `../plans/06_missing_data.md`
1. **`01_generate.py`** — Generate simulation inputs, Makefile, and
   mutation count diagnostics
2. **`02_run.py`** — Run Delphy via `make -jN`
3. **`03_analyze.py`** — Run loganalyser, check ESS, compute
   coverage/ranks
4. **`04_plot.py`** — Produce all plots from TSV files

---

## Part 1: `01_generate.py` — Generate

Copy from `05_site_rate_heterogeneity/01_generate.py` with these
changes:

### Configuration changes

Add:
```python
MISSING_DATA_MEAN_NUM_GAPS = 3.0
MISSING_DATA_MEAN_GAP_LENGTH = 500.0
MISSING_DATA_MEAN_NUM_MISSING_SITES = 3.0
```

### Step 2 changes: Run Sapling

`run_sapling()` adds three arguments to the sapling command:
```
"--missing-data-mean-num-gaps", str(MISSING_DATA_MEAN_NUM_GAPS),
"--missing-data-mean-gap-length", str(MISSING_DATA_MEAN_GAP_LENGTH),
"--missing-data-mean-num-missing-sites", str(MISSING_DATA_MEAN_NUM_MISSING_SITES),
```

No other changes to `run_sapling()` or `sample_prior_params()`.

### All other steps

No changes.  The Makefile rule and all parameters are identical to
`05_site_rate_heterogeneity`.

---

## Part 2: `02_run.py` — Run

Identical to `05_site_rate_heterogeneity/02_run.py`.  Copy verbatim.

---

## Part 3: `03_analyze.py` — Analyze

Identical to `05_site_rate_heterogeneity/03_analyze.py`.  Copy verbatim.

---

## Part 4: `04_plot.py` — Plot

Identical to `05_site_rate_heterogeneity/04_plot.py`.  Copy verbatim.

---

## Directory structure at runtime

Same as `05_site_rate_heterogeneity`, except:

- Each `sim_NNN/` also contains `sim-COMPLETE.maple` and
  `sim-COMPLETE.fasta` (the complete versions without masking).

---

## Execution workflow

```bash
cd wcss/06_missing_data

# --- Final run ---
./01_generate.py
./02_run.py
./03_analyze.py
./04_plot.py

# --- Debug run ---
./01_generate.py --n 10 --steps 20000000
./02_run.py
./03_analyze.py --n 10 --ignore-low-ess --force-include-all-replicates
./04_plot.py
```

---

## Potential concerns

- **ESS impact of missing data:** Missing data reduces the effective
  amount of information in the sequences, which could widen posteriors
  and potentially affect mixing.  At ~5% missing data, this should be
  a modest effect.  If ESS is noticeably worse than
  `05_site_rate_heterogeneity`, increase steps.

- **Coverage should still be ~95%:** Missing data doesn't introduce
  model misspecification — Delphy correctly handles `N` sites by
  marginalizing over all possible bases.  Wider posteriors due to less
  data are expected, but coverage should remain well-calibrated.

- **Seed stability:** Sapling draws missing data after mutation
  simulation, so the underlying tree and mutations are identical to
  what would be produced without missing data at the same seed.  The
  `-COMPLETE` files can be used to verify this.
