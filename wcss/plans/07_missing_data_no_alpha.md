# Plan: WCSS with Missing Data (No Site-Rate Heterogeneity)

## Goal

Validate Delphy's inference under missing data without site-rate
heterogeneity.  This is the same as `04_free_exp_pop` (free mu, n0, g,
kappa, pi) but with simulated missing data added to the input
sequences.

Preliminary results from `05_site_rate_heterogeneity` suggest that for
this configuration (200 tips, L=30000, exponential population model),
the data is nearly uninformative about the site-rate heterogeneity
parameter `alpha`, making it difficult to validate.  This study
therefore focuses on validating the core model parameters under missing
data, with site-rate heterogeneity disabled.

## Differences from `04_free_exp_pop`

1. **Sapling is called with missing data options:**
   - `--missing-data-mean-num-gaps 3`
   - `--missing-data-mean-gap-length 500`
   - `--missing-data-mean-num-missing-sites 3`

   This gives ~5% missing data per tip sequence.

2. **Delphy reads the masked `sim.maple`** (with `N` ranges).  No
   change to the Delphy invocation — site-rate heterogeneity remains
   disabled (no `--v0-site-rate-heterogeneity` flag).

That's it.  All priors, parameters, analysis, and plotting are
identical to `04_free_exp_pop`.

## Configuration

- **Directory:** `wcss/07_missing_data_no_alpha/`
- **Tips:** Same as before (200 tips, dates uniform over 2025)
- **N:** 200 (final), 10 (debug)
- **Steps:** 1,000,000,000
- **Missing data:**
  - Mean number of gaps per tip: 3
  - Mean gap length: 500 sites
  - Mean number of isolated missing sites per tip: 3
- **Sampled parameters (drawn fresh per replicate):**
  - Effective population size n0 ~ InvGamma(mean=3 years, stddev=1 year)
  - Exponential growth rate g ~ Exponential(mean=1 e-folding/year),
    constrained to g >= 0
  - Mutation rate mu ~ Gamma(alpha=1, beta=1000)
    (i.e., Exponential with mean 1e-3 subst/site/year)
  - Stationary frequencies (pi_A, pi_C, pi_G, pi_T) ~ Dirichlet(1,1,1,1)
  - HKY kappa ~ LogNormal(meanlog=1.0, sdlog=1.25)
- **Fixed parameters:**
  - Genome size L = 30,000 sites
  - No site-rate heterogeneity (all nu_l = 1)

## Deliverables

Files in `wcss/07_missing_data_no_alpha/`:

0. **`00_plan.md`** — Symlink to `../plans/07_missing_data_no_alpha.md`
1. **`01_generate.py`** — Generate simulation inputs, Makefile, and
   mutation count diagnostics
2. **`02_run.py`** — Run Delphy via `make -jN`
3. **`03_analyze.py`** — Run loganalyser, check ESS, compute
   coverage/ranks
4. **`04_plot.py`** — Produce all plots from TSV files

---

## Part 1: `01_generate.py` — Generate

Copy from `04_free_exp_pop/01_generate.py` with these changes:

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

No `--site-rate-heterogeneity-alpha` argument (site-rate heterogeneity
is disabled).

### All other steps

No changes.  The Makefile rule is identical to `04_free_exp_pop` (no
`--v0-site-rate-heterogeneity` flag).

---

## Part 2: `02_run.py` — Run

Identical to `04_free_exp_pop/02_run.py`.  Copy verbatim.

---

## Part 3: `03_analyze.py` — Analyze

Identical to `04_free_exp_pop/03_analyze.py`.  Copy verbatim.

(PARAMS list: mu, n0, g, kappa, pi_A–pi_T, rootHeight.  No alpha.)

---

## Part 4: `04_plot.py` — Plot

Identical to `04_free_exp_pop/04_plot.py`.  Copy verbatim.

---

## Execution workflow

```bash
cd wcss/07_missing_data_no_alpha

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

- **ESS impact of missing data:** With ~5% missing data and no
  site-rate heterogeneity, the model is simpler than
  `06_missing_data`.  ESS should be comparable to `04_free_exp_pop`.

- **Coverage should still be ~95%:** Delphy correctly handles `N`
  sites by marginalizing over all possible bases.  Wider posteriors
  from less data are expected but coverage should remain calibrated.
