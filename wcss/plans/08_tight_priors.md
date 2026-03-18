# Plan: WCSS with Tighter Priors on Mutation Rate and Growth Rate

## Goal

Delphy works best when the number of mutations is comparable to the
number of branches (2 * num_tips = 400 for our 200-tip trees).  With
the current wide priors on mu and g, many replicates produce mutation
counts that are far too high (thousands) or nearly zero, pushing
Delphy out of its sweet spot.

This study tightens the priors on mu and g to keep the mutation counts
in a more reasonable range.  It also collects and plots mutation count
statistics from Sapling to help tune the priors.

Missing data is included (same as `07_missing_data_no_alpha`).
Site-rate heterogeneity is disabled.

## Differences from `07_missing_data_no_alpha`

1. **Configurable prior on mu:** Specified via `--mu-prior-mean` and
   `--mu-prior-stddev` (matching Delphy's `--v0-mu-prior-mean` /
   `--v0-mu-prior-stddev` parametrization).  Internally converted to
   Gamma shape/rate for sampling:
   `alpha = (mean/stddev)^2`, `beta = mean/stddev^2`.

2. **Configurable prior on g:** The full bounded Laplace prior is
   exposed via `--g-prior-mu`, `--g-prior-scale`, `--g-min`, and
   `--g-max`.  A convenience shorthand `--g-prior-exponential-with-mean`
   is also available (sets mu=0, scale=|mean|, and the appropriate
   bound), matching Delphy's `--v0-pop-g-prior-exponential-with-mean`.

3. **Mutation count statistics:** After all Sapling simulations,
   `01_generate.py` reads `num_mutations` from each replicate's
   `sim_info.json`, prints summary statistics, saves them to
   `sims/mutation_counts.tsv`, and generates a histogram and eCDF
   plot in `plots/`.

Everything else (tips, genome size, n0 prior, kappa prior, missing
data, Delphy invocation, analysis, and plotting) is unchanged.

## Configuration

- **Directory:** `wcss/08_tight_priors/`
- **Tips:** Same as before (200 tips, dates uniform over 2025)
- **N:** 200 (final), 10 (debug)
- **Steps:** 1,000,000,000
- **Missing data:** Same as `07_missing_data_no_alpha` (3 gaps of mean
  length 500, 3 isolated missing sites)
- **Sampled parameters (drawn fresh per replicate):**
  - Mutation rate mu ~ Gamma(alpha, beta), where alpha and beta are
    derived from `--mu-prior-mean` and `--mu-prior-stddev`.
    Default: mean=1e-3, stddev=1e-3 (i.e., Exponential with mean 1e-3,
    same as previous studies).
  - Growth rate g ~ Truncated Laplace(mu, scale, g_min, g_max).
    Default: Exponential with mean 1.0 (i.e., mu=0, scale=1, g_min=0,
    g_max=+inf, same as previous studies).
  - n0 ~ InvGamma(mean=3 years, stddev=1 year)
  - pi ~ Dirichlet(1,1,1,1)
  - kappa ~ LogNormal(meanlog=1.0, sdlog=1.25)
- **Fixed parameters:**
  - Genome size L = 30,000 sites
  - No site-rate heterogeneity

## Deliverables

Files in `wcss/08_tight_priors/`:

0. **`00_plan.md`** — Symlink to `../plans/08_tight_priors.md`
1. **`01_generate.py`** — Generate simulation inputs, Makefile, and
   mutation count statistics/plots
2. **`02_run.py`** — Run Delphy via `make -jN`
3. **`03_analyze.py`** — Run loganalyser, check ESS, compute
   coverage/ranks
4. **`04_plot.py`** — Produce all plots from TSV files

---

## Part 1: `01_generate.py` — Generate

Copy from `07_missing_data_no_alpha/01_generate.py` with these changes:

### Configuration changes

Replace the fixed prior constants with defaults:

```python
# Mu prior: Gamma distribution specified via mean and stddev
# Default: Exponential with mean 1e-3 (alpha=1, beta=1000)
DEFAULT_MU_PRIOR_MEAN = 1e-3
DEFAULT_MU_PRIOR_STDDEV = 1e-3

# G prior: Truncated Laplace(mu, scale) on [g_min, g_max]
# Default: Exponential with mean 1.0 (mu=0, scale=1, g_min=0)
DEFAULT_G_PRIOR_MU = 0.0
DEFAULT_G_PRIOR_SCALE = 1.0
DEFAULT_G_MIN = 0.0
DEFAULT_G_MAX = float("inf")
```

### CLI changes

Add arguments for the mu prior:

```python
parser.add_argument("--mu-prior-mean", type=float,
                    default=DEFAULT_MU_PRIOR_MEAN,
                    help="Mean of the Gamma prior on mu (subst/site/year)")
parser.add_argument("--mu-prior-stddev", type=float,
                    default=DEFAULT_MU_PRIOR_STDDEV,
                    help="Std dev of the Gamma prior on mu (subst/site/year)")
```

Add arguments for the g prior (two modes, mutually exclusive):

```python
# Direct Laplace parametrization
parser.add_argument("--g-prior-mu", type=float, default=None,
                    help="Location of the Laplace prior on g (e-foldings/year)")
parser.add_argument("--g-prior-scale", type=float, default=None,
                    help="Scale of the Laplace prior on g (e-foldings/year)")
parser.add_argument("--g-min", type=float, default=None,
                    help="Lower bound on g (e-foldings/year)")
parser.add_argument("--g-max", type=float, default=None,
                    help="Upper bound on g (e-foldings/year)")

# Shorthand for exponential prior
parser.add_argument("--g-prior-exponential-with-mean", type=float,
                    default=None,
                    help="Shorthand: Exponential prior on g with this mean")
```

If `--g-prior-exponential-with-mean` is given, it sets mu=0,
scale=|mean|, and g_min=0 (if mean > 0) or g_max=0 (if mean < 0).
If neither mode is specified, the defaults apply (Exponential with
mean 1.0).

### Sampling changes

`sample_prior_params()` takes the prior parameters as arguments:

- **mu:** Convert mean/stddev to Gamma alpha/beta:
  `alpha = (mean/stddev)^2`, `beta = mean/stddev^2`.
  Sample `mu = rng.gamma(alpha, 1/beta)`.

- **g:** Sample from the truncated Laplace using inverse transform
  sampling via `scipy.stats.laplace`.  Compute
  `u_lo = cdf(g_min, loc=mu, scale=scale)` and
  `u_hi = cdf(g_max, loc=mu, scale=scale)`, draw
  `u ~ Uniform(u_lo, u_hi)`, then `g = ppf(u, loc=mu, scale=scale)`.
  This is exact and efficient for any truncation bounds.  (Rejection
  sampling was considered but discarded since it can become
  arbitrarily slow for tight bounds.)

### Makefile generation

`generate_makefile()` receives the prior parameters and generates
the appropriate Delphy arguments:

- For mu: `--v0-mu-prior-mean` and `--v0-mu-prior-stddev`.
- For g: always use the direct Laplace parametrization
  (`--v0-pop-g-prior-mu` and `--v0-pop-g-prior-scale`), regardless
  of how the user specified the prior on the CLI.  Only emit
  `--v0-pop-growth-rate-min` if `g_min` is finite, and
  `--v0-pop-growth-rate-max` if `g_max` is finite.  (Delphy's
  defaults for these are BEAUti values that don't match ours, so
  explicit values are always needed.)

Note: `01_generate.py` adds `scipy` as a dependency (for
`scipy.stats.laplace` in the truncated Laplace sampling).

### Step 4 (new): Mutation count statistics

After all Sapling simulations, read `num_mutations` from each
replicate's `sim_info.json` and:

1. Print summary statistics (min, max, mean, median, p5, p95).
2. Save to `sims/mutation_counts.tsv` (columns: replicate,
   num_mutations, mu, g).
3. Generate two plots in `plots/`:
   - `plots/mutation_counts_histogram.pdf` — histogram of mutation
     counts with a vertical line at 2 * NUM_TIPS = 400.
   - `plots/mutation_counts_ecdf.pdf` — eCDF of mutation counts with
     a vertical line at 2 * NUM_TIPS = 400.

These plots are generated directly by `01_generate.py` (not
`04_plot.py`) since they characterize the simulation inputs, not
the inference results.

---

## Part 2: `02_run.py` — Run

Identical to `07_missing_data_no_alpha/02_run.py`.  Copy verbatim.

---

## Part 3: `03_analyze.py` — Analyze

Identical to `07_missing_data_no_alpha/03_analyze.py`.  Copy verbatim.

---

## Part 4: `04_plot.py` — Plot

Identical to `07_missing_data_no_alpha/04_plot.py`.  Copy verbatim.

---

## Execution workflow

```bash
cd wcss/08_tight_priors

# --- Explore priors (generate only, inspect mutation count plots) ---
./01_generate.py --n 200 --steps 20000000
# Look at plots/mutation_counts_histogram.pdf
# Adjust priors as needed, e.g.:
./01_generate.py --n 200 --steps 20000000 \
    --mu-prior-mean 5e-4 --mu-prior-stddev 3e-4 \
    --g-prior-exponential-with-mean 0.5
# Re-run until satisfied with mutation count distribution

# --- Full run with chosen priors ---
./01_generate.py \
    --mu-prior-mean <M> --mu-prior-stddev <S> \
    --g-prior-exponential-with-mean <G>
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

## Prior tuning results

We iterated over prior choices using 10 debug replicates (`--n 10`)
to find values that keep mutation counts near the target of
2 × 200 tips = 400.

| Run | mu mean | mu stddev | g prior | Mutation counts (min–p5–median–mean–p95–max) | Verdict |
|-----|---------|-----------|---------|----------------------------------------------|---------|
| 1 | 1e-3 | 1e-3 | Exp(mean=1) (defaults) | 5–156–921–1051–2046–2062 | Far too wide; range spans 5 to 2000+, some replicates nearly mutation-free |
| 2 | 1e-3 | 1e-4 | Laplace(mu=1, scale=0.2, g_min=0.5) | 661–706–978–963–1215–1251 | Tight g helps, but mu still too high — all replicates above target |
| 3 | 1e-3 | 1e-4 | Laplace(mu=2, scale=0.2, g_min=0.5) | 493–496–715–692–889–915 | Higher g shrinks trees, reducing mutations, but still ~2× target |
| 4 | 5e-4 | 1e-4 | Laplace(mu=2, scale=0.2, g_min=0.5) | 172–188–374–371–540–552 | Centered on target (median 374 vs target 400), reasonable spread |

**Chosen priors for final run:**
- mu ~ Gamma(mean=5e-4, stddev=1e-4)
- g ~ Laplace(mu=2, scale=0.2, g_min=0.5)

---

## Potential concerns

- **Prior tightening vs. WCSS validity:** The WCSS is still valid as
  long as the Delphy priors match the simulation priors.  Since we
  pass the same prior parameters to both sapling (for sampling) and
  delphy (in the Makefile), the test remains self-consistent.

- **Choosing good priors:** The mutation count plots from
  `01_generate.py` make it easy to iterate: run with different prior
  values, check the histogram, and pick values that concentrate most
  replicates in a reasonable range (say, 100–1000 mutations).

