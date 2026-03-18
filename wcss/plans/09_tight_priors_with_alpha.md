# Plan: WCSS with Tighter Priors and Site-Rate Heterogeneity

## Goal

Combine the tighter priors on mu and g from `08_tight_priors` with
the site-rate heterogeneity from `05_site_rate_heterogeneity`.  This
validates Delphy's inference of alpha alongside all other parameters,
using priors that keep mutation counts in a reasonable range.

Missing data is included (same as `07_missing_data_no_alpha` and
`08_tight_priors`).  Site-rate heterogeneity is enabled.

## Differences from `08_tight_priors`

1. **`alpha` is drawn from the prior** per replicate:
   `alpha ~ Exponential(mean=1.0)`.
2. **Sapling is called with `--site-rate-heterogeneity-alpha`** set to
   the sampled value of `alpha`.
3. **Delphy is called with `--v0-site-rate-heterogeneity`** in the
   Makefile to enable inference of `alpha` and `nu_l`.
4. **`alpha` is added to PARAMS** in `03_analyze.py` and `04_plot.py`.
5. **`read_true_params()`** reads `alpha` from
   `info["subst_model"]["site_rate_heterogeneity_alpha"]`.

Everything else (tips, genome size, n0 prior, kappa prior, missing
data, configurable mu/g priors, mutation count diagnostics) is
inherited from `08_tight_priors`.

## Configuration

- **Directory:** `wcss/09_tight_priors_with_alpha/`
- **Tips:** Same as before (200 tips, dates uniform over 2025)
- **N:** 200 (final), 10 (debug)
- **Steps:** 1,000,000,000
- **Missing data:** Same as `08_tight_priors` (3 gaps of mean
  length 500, 3 isolated missing sites)
- **Sampled parameters (drawn fresh per replicate):**
  - Mutation rate mu ~ Gamma(mean, stddev), configurable via CLI.
    Default: mean=5e-4, stddev=1e-4.
  - Growth rate g ~ Truncated Laplace(mu, scale, g_min, g_max),
    configurable via CLI.
    Default: Laplace(mu=2, scale=0.2, g_min=0.5).
  - Site-rate heterogeneity alpha ~ Exponential(mean=1.0)
  - n0 ~ InvGamma(mean=3 years, stddev=1 year)
  - pi ~ Dirichlet(1,1,1,1)
  - kappa ~ LogNormal(meanlog=1.0, sdlog=1.25)
- **Fixed parameters:**
  - Genome size L = 30,000 sites

## Deliverables

Files in `wcss/09_tight_priors_with_alpha/`:

0. **`00_plan.md`** — Symlink to `../plans/09_tight_priors_with_alpha.md`
1. **`01_generate.py`** — Generate simulation inputs, Makefile, and
   mutation count statistics/plots
2. **`02_run.py`** — Run Delphy via `make -jN`
3. **`03_analyze.py`** — Run loganalyser, check ESS, compute
   coverage/ranks
4. **`04_plot.py`** — Produce all plots from TSV files

---

## Part 1: `01_generate.py` — Generate

Copy from `08_tight_priors/01_generate.py` with these changes:

### Configuration changes

Add:
```python
ALPHA_PRIOR_MEAN = 1.0
```

### Default mu/g priors

Update the defaults to match the chosen values from `08_tight_priors`:
```python
DEFAULT_MU_PRIOR_MEAN = 5e-4
DEFAULT_MU_PRIOR_STDDEV = 1e-4

DEFAULT_G_PRIOR_MU = 2.0
DEFAULT_G_PRIOR_SCALE = 0.2
DEFAULT_G_MIN = 0.5
DEFAULT_G_MAX = float("inf")
```

### Sampling changes

`sample_prior_params()` additionally draws:
```python
alpha = rng.exponential(ALPHA_PRIOR_MEAN)
```

Returns `n0, g, mu, kappa, pi, alpha`.

### Step 2 changes: Run Sapling

`run_sapling()` gains an `alpha` parameter and passes
`--site-rate-heterogeneity-alpha {alpha}` to sapling.
The print statement includes `alpha` in its output.

### Makefile generation

`generate_makefile()` adds `--v0-site-rate-heterogeneity` to the
Delphy invocation (boolean flag, no value).

---

## Part 2: `02_run.py` — Run

Identical to `08_tight_priors/02_run.py`.  Copy verbatim.

---

## Part 3: `03_analyze.py` — Analyze

Copy from `08_tight_priors/03_analyze.py` with these changes:

### Configuration changes

Add `("alpha", "alpha")` to the `PARAMS` list:
```python
PARAMS = [
    ("mu",         "meanRate"),
    ("alpha",      "alpha"),
    ("n0",         "exponential.popSize"),
    ("g",          "exponential.growthRate"),
    ("kappa",      "kappa"),
    ("pi_A",       "frequencies1"),
    ("pi_C",       "frequencies2"),
    ("pi_G",       "frequencies3"),
    ("pi_T",       "frequencies4"),
    ("rootHeight", "rootHeight"),
]
```

### Step 2 changes: Read true parameters

`read_true_params()` additionally reads `alpha`:
```python
"alpha": info["subst_model"]["site_rate_heterogeneity_alpha"],
```

---

## Part 4: `04_plot.py` — Plot

Copy from `08_tight_priors/04_plot.py` with these changes:

### Configuration changes

Add `("alpha", "alpha")` to the `PARAMS` list (same position as
in `03_analyze.py`).

---

## Execution workflow

```bash
cd wcss/09_tight_priors_with_alpha

# --- Full run with default priors (from 08_tight_priors tuning) ---
./01_generate.py
./02_run.py
./03_analyze.py
./04_plot.py

# --- Debug run ---
./01_generate.py --n 10 --steps 20000000
./02_run.py
./03_analyze.py --n 10 --ignore-low-ess --force-include-all-replicates
./04_plot.py

# --- Custom priors (if needed) ---
./01_generate.py \
    --mu-prior-mean 5e-4 --mu-prior-stddev 1e-4 \
    --g-prior-mu 2.0 --g-prior-scale 0.2 --g-min 0.5
```

---

## Potential concerns

- **ESS for alpha:** The alpha move uses 10 scaling sub-moves per MCMC
  step, but alpha controls 30,000 latent `nu_l` variables.  If ESS for
  `alpha` is consistently low, increase steps to 2B.

- **Interaction between alpha and mu:** With site-rate heterogeneity,
  the effective per-site mutation rate varies.  The tighter mu prior
  (stddev=1e-4 around mean=5e-4) should still produce reasonable
  mutation counts since the Gamma(alpha, alpha) rate modifiers have
  mean 1.

- **Prior tightening vs. WCSS validity:** Same as `08_tight_priors` —
  the WCSS remains valid as long as the simulation priors match the
  Delphy priors, which they do by construction.
