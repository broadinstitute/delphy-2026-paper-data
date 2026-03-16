# Plan: WCSS with Tip-Date Uncertainty

## Goal

Validate Delphy's inference when some tip dates are uncertain (specified
at month or year precision instead of exact day).  Sapling simulates
data with full dates, then masks 15% of tip dates to month precision
(YYYY-MM) and 5% to year precision (YYYY) in the output MAPLE file.
Delphy receives the masked MAPLE file and must infer the exact tip
dates as part of its MCMC.

This builds on `09_tight_priors_with_alpha`, so it includes tight
priors on mu and g, site-rate heterogeneity, and missing data.

This uses a new feature in Sapling:

- **Sapling:** Tip-date uncertainty (`--p-tip-date-uncertain-upto-month`,
  `--p-tip-date-uncertain-upto-year`).

## Differences from `09_tight_priors_with_alpha`

1. **Sapling is invoked with tip-date uncertainty options:**
   `--p-tip-date-uncertain-upto-month 0.15` and
   `--p-tip-date-uncertain-upto-year 0.05`.  This means 15% of tips
   get month-level dates, 5% get year-level dates, and 80% keep exact
   dates.
2. **Sapling produces `-COMPLETE` files** (with full dates and full
   sequences) alongside the normal files (with masked dates and missing
   data).  Delphy receives the normal `sim.maple` (with uncertain dates).
3. **Delphy invocation is unchanged.**  Delphy automatically detects
   uncertain tip dates in the MAPLE file and samples exact dates during
   MCMC.
4. **`03_analyze.py` adds tip-date posterior analysis.**  For each
   replicate, it identifies which tips have uncertain dates, reads
   their true dates from the `-COMPLETE` MAPLE file, reads their
   posterior date samples from the Delphy log file, computes ranks
   and coverage.  Month-uncertain and year-uncertain tips are analyzed
   separately.
5. **`04_plot.py` adds tip-date rows** to the summary figure: one row
   for month-uncertain tips and one for year-uncertain tips.

Everything else (tips, genome size, n0 prior, kappa prior, mu/g priors,
alpha prior, missing data, mutation count diagnostics, analyzed
parameters) is inherited from `09_tight_priors_with_alpha`.

## Configuration

- **Directory:** `wcss/11_tip_date_uncertainty/`
- **Tips:** Same as `09_tight_priors_with_alpha` (200 tips, dates
  uniform over 2025)
- **N:** 200 (final), 10 (debug)
- **Steps:** 1,000,000,000 (5,000,000 per tip x 200 tips);
  debug: 20,000,000 (100,000 per tip x 200 tips)
- **Tip-date uncertainty:**
  - 15% of tips masked to month precision (YYYY-MM)
  - 5% of tips masked to year precision (YYYY)
  - 80% of tips keep exact dates (YYYY-MM-DD)
- **Missing data:** Same as `09_tight_priors_with_alpha` (3 gaps of
  mean length 500, 3 isolated missing sites)
- **Sampled parameters:** Same as `09_tight_priors_with_alpha`:
  - mu ~ Gamma(mean=5e-4, stddev=1e-4)
  - g ~ Truncated Laplace(mu=2, scale=0.2, g_min=0.5)
  - n0 ~ InvGamma(mean=3 yr, stddev=1 yr)
  - alpha ~ Exponential(mean=1)
  - pi ~ Dirichlet(1,1,1,1)
  - kappa ~ LogNormal(meanlog=1.0, sdlog=1.25)
- **Fixed parameters:**
  - Genome size L = 30,000 sites

## Analyzed parameters

Same as `09_tight_priors_with_alpha`:

| Display name | Log column             | True value source                        |
|--------------|------------------------|------------------------------------------|
| mu           | meanRate               | subst_model.mu                           |
| alpha        | alpha                  | subst_model.site_rate_heterogeneity_alpha |
| n0           | exponential.popSize    | pop_model.n0                             |
| g            | exponential.growthRate  | pop_model.g                              |
| kappa        | kappa                  | subst_model.kappa                        |
| pi_A         | frequencies1           | subst_model.pi[0]                        |
| pi_C         | frequencies2           | subst_model.pi[1]                        |
| pi_G         | frequencies3           | subst_model.pi[2]                        |
| pi_T         | frequencies4           | subst_model.pi[3]                        |
| rootHeight   | rootHeight             | tree_stats.tree_height                   |

## Tip-date posterior analysis

### Overview

In addition to the standard WCSS parameter analysis, this study validates
Delphy's inference of uncertain tip dates.  For each replicate, some tips
have dates masked to month or year precision.  Delphy samples exact dates
during MCMC, and we validate these posteriors against the true dates.

### Identifying uncertain tips

Sapling's `sim_info.json` records only the masking *probabilities*
(not a per-tip breakdown).  Instead, identify uncertain tips by
scanning the Delphy log file header for columns matching `age(...)`.
Delphy emits an `age(TIP_XXX|DATE)` column for each tip whose date
range spans more than one day (`t_min != t_max`).

The uncertainty type is determined from the date format embedded in
the column name:

- `age(TIP_XXX|YYYY-MM)` (7-char date after `|`) → **month-uncertain**
- `age(TIP_XXX|YYYY)` (4-char date after `|`) → **year-uncertain**

The two categories are analyzed separately.

### Reading true tip dates

True (exact) tip dates come from the `-COMPLETE` MAPLE file
(`sim-COMPLETE.maple`) produced by Sapling alongside the masked file.
Parse the tip names and dates from the COMPLETE file's header.

### Reading posterior tip dates from Delphy log files

Delphy's log file includes columns `age(TIP_XXX|...)` for each tip
with an uncertain date.  The value is `(beast_t0 - tip_t) / 365.0`,
where `beast_t0 = max(tip_t)` across all tips (varies per MCMC step
if the latest tip itself has an uncertain date).

To convert to an absolute calendar date, use the identity:

```
tip_calendar_year = rootHeight + age(root) - age(TIP)
```

where `rootHeight` and `age(root)` are also log columns.  This
cancels out `beast_t0` and gives the tip date in calendar years
(e.g., 2025.456).

The true tip date is similarly converted to a calendar year for
comparison.

### Custom analysis (not using loganalyser)

Because different replicates have different sets of uncertain tips
(and therefore different `age(TIP_XXX|...)` columns), loganalyser
cannot be used for tip-date analysis.  Instead, `03_analyze.py`
reads the raw posterior samples directly:

1. **For each replicate**, identify uncertain tips and their log
   columns.
2. **Read the raw log file**, apply burnin, extract posterior samples
   for each uncertain tip's `age(TIP_XXX|...)`, `rootHeight`, and
   `age(root)`.
3. **Convert to calendar years** using the formula above, for both
   the posterior samples and the true value.
4. **Compute normalized rank** of the true value within the posterior
   samples.
5. **Compute 95% HPD** from the posterior samples and check whether
   the true value is covered.

### Aggregation across replicates

All uncertain tip ranks are pooled across replicates, separated by
uncertainty type.  With N=200 replicates, ~200 * 0.15 * 200 = 6000
month-uncertain tip-date samples and ~200 * 0.05 * 200 = 2000
year-uncertain tip-date samples are available (exact counts depend
on Sapling's random masking choices).

### Output files

- `analyses/tip_date_ranks_month.tsv` — one row per month-uncertain
  tip across all replicates: replicate, tip_name, true_date,
  posterior_mean, hpd_lo, hpd_hi, normalized_rank
- `analyses/tip_date_ranks_year.tsv` — same for year-uncertain tips
- `analyses/tip_date_coverage_summary.tsv` — columns:
  `uncertainty_type`, `n_tips`, `n_covered`, `coverage`, `expected`,
  `binom_2.5%`, `binom_97.5%`.  One row for "month", one for "year".

### Plots (in `04_plot.py`)

For each uncertainty type (month, year):

1. **Rank histogram** — histogram of normalized ranks; should be
   uniform.
2. **ECDF** — empirical CDF of normalized ranks with 95% confidence
   band around the diagonal.
3. **Scatter plot** — true tip date vs posterior mean, with 95% HPD
   error bars (blue if covered, red if not).

These appear as two extra rows at the bottom of the summary figure
(before the clade coverage row): one row for month-uncertain tips,
one for year-uncertain tips.

## Deliverables

Files in `wcss/11_tip_date_uncertainty/`:

0. **`00_plan.md`** -- Symlink to `../plans/11_tip_date_uncertainty.md`
1. **`01_generate.py`** -- Generate simulation inputs, Makefile, and
   mutation count statistics/plots
2. **`02_run.py`** -- Run Delphy via `make -jN`
3. **`03_analyze.py`** -- Run loganalyser, check ESS, compute
   coverage/ranks for standard parameters, plus tip-date posterior
   analysis (ranks, HPD, coverage by uncertainty type)
4. **`04_plot.py`** -- Produce all plots from TSV files, including
   tip-date rank histograms, ECDFs, scatter plots, and summary figure
   with tip-date rows

---

## Part 1: `01_generate.py` -- Generate

Copy from `09_tight_priors_with_alpha/01_generate.py` with these
changes:

### Configuration changes

Add:
```python
P_TIP_DATE_UNCERTAIN_UPTO_MONTH = 0.15
P_TIP_DATE_UNCERTAIN_UPTO_YEAR = 0.05
```

### Sapling invocation changes

`run_sapling()` adds two flags to the Sapling command:
```
--p-tip-date-uncertain-upto-month 0.15
--p-tip-date-uncertain-upto-year 0.05
```

No other changes to `run_sapling()`.  The normal `sim.maple` output
(with masked dates) is what Delphy will read.

### Makefile generation

The Delphy invocation is identical to `09_tight_priors_with_alpha`.

The `clean` target adds `sim_*/delphy-digested.log` to the `rm`
command (see Part 3 for why these files exist).

---

## Part 2: `02_run.py` -- Run

Identical to `09_tight_priors_with_alpha/02_run.py`.  Copy verbatim.

---

## Part 3: `03_analyze.py` -- Analyze

Copy from `09_tight_priors_with_alpha/03_analyze.py` with these
additions:

### Loganalyser preprocessing: `delphy-digested.log`

Delphy's raw log files include `age(TIP_XXX|...)` columns for each
tip with an uncertain date.  Because different replicates mask
different tips, these columns vary across replicates.  loganalyser
requires uniform columns across all input files, so it cannot process
the raw log files directly.

Solution: before running loganalyser, strip all `age(...)` columns
from each replicate's `delphy.log` into a `delphy-digested.log` file.
loganalyser then runs on the digested files (which have uniform
columns).  The digested files are kept for reproducibility and are
cleaned up by `make clean`.

```python
def _strip_age_columns(src_path, dst_path):
    """Copy src_path to dst_path, dropping any age(...) columns."""
```

### New: `analyze_tip_dates()`

After the standard parameter analysis (Steps 1–5), add a new step
that analyzes tip-date posteriors:

```python
def analyze_tip_dates(script_dir, analyses_dir, included, burnin_frac):
```

For each replicate `i` in `included`:

1. Read the Delphy log file header to find `age(...)` columns.
   Classify each as month-uncertain or year-uncertain based on the
   date format in the column name (see "Identifying uncertain tips").
   Extract the tip name (everything between `age(` and `|`).
2. Read true tip dates from `sim-COMPLETE.maple` (parse tip names
   and dates from the MAPLE header).
3. Apply burnin to the log file, extract posterior samples for each
   uncertain tip's `age(TIP_XXX|...)`, `rootHeight`, and `age(root)`.
4. Convert to calendar years using the formula above, for both
   the posterior samples and the true value.
5. Compute normalized rank of the true value within the posterior
   samples.
6. Compute 95% HPD: sort posterior samples, find the shortest
   interval containing 95% of samples.  Check whether the true
   value falls within.

Collect results into two lists (month, year) and save:

- `tip_date_ranks_month.tsv`
- `tip_date_ranks_year.tsv`
- `tip_date_coverage_summary.tsv`

### Reading true tip dates from COMPLETE MAPLE

The `-COMPLETE` MAPLE file has the same format as regular MAPLE but
with exact dates.  Parse the `>TIP_XXX|YYYY-MM-DD` lines.

**Important:** The true date must be converted to a calendar year
using the same coordinate system as Delphy's posterior samples.
Delphy uses `to_linear_year(t) = 2020.0 + t / 365.0`, where `t` is
days since 2020-01-01.  The `rootHeight + age(root) - age(TIP)`
formula produces values in this same system.  So the true date
conversion is:

```python
from datetime import date
days = (tip_date - date(2020, 1, 1)).days
true_calendar_year = 2020.0 + days / 365.0
```

Do NOT use `year + (day_of_year - 1) / days_in_year`, which accounts
for leap years differently and causes a systematic error of up to
~5 days for dates 5 years from the epoch.

### Integration into main()

After the existing Step 5 (clade coverage), add:

```python
# Step 6: Tip-date analysis
print("\nAnalyzing tip-date posteriors...")
analyze_tip_dates(script_dir, analyses_dir, included, burnin_frac)
```

---

## Part 4: `04_plot.py` -- Plot

Copy from `09_tight_priors_with_alpha/04_plot.py` with these
additions:

### New tip-date plots

Add standalone plots for each uncertainty type:

- `rank_histogram_tip_dates_month.pdf`
- `rank_histogram_tip_dates_year.pdf`
- `ecdf_tip_dates_month.pdf`
- `ecdf_tip_dates_year.pdf`
- `scatter_tip_dates_month.pdf`
- `scatter_tip_dates_year.pdf`

### Summary figure changes

The summary figure gains two extra rows between the parameter rows
and the clade coverage row:

- Row N+0: **Month-uncertain tip dates** — rank histogram, ECDF,
  scatter plot (using data from `tip_date_ranks_month.tsv`)
- Row N+1: **Year-uncertain tip dates** — rank histogram, ECDF,
  scatter plot (using data from `tip_date_ranks_year.tsv`)

These rows use the same 3-column layout as parameter rows.  The
scatter plot shows true tip date (calendar year) on x-axis vs
posterior mean on y-axis, with 95% HPD error bars.

The tip-date rows have many more data points than standard parameter
rows (thousands vs hundreds), so the rank histogram bins and scatter
point sizes should be adjusted accordingly.  Use a smaller point size
(s=3) and more bins (50) for the rank histograms.

---

## Execution workflow

```bash
cd wcss/11_tip_date_uncertainty

# --- Final run (200 replicates, 1B steps — the defaults) ---
./01_generate.py
./02_run.py
./03_analyze.py
./04_plot.py

# --- Debug run (10 short replicates, for pipeline testing) ---
./01_generate.py --n 10 --steps 20000000
./02_run.py
./03_analyze.py --n 10 --ignore-low-ess --force-include-all-replicates
./04_plot.py
```

---

## Potential concerns

- **ESS for tip dates:**  Delphy samples exact tip dates via uniform
  proposals within the uncertain range (month or year).  With 200 tips
  and 20% uncertain, ~40 tip dates are sampled.  ESS should be fine
  since each tip date is a low-dimensional parameter.

- **Coverage for rootHeight:**  With uncertain tip dates, the root
  height depends on the sampled tip dates.  The true rootHeight from
  `sim_info.json` uses the exact dates, so coverage should still be
  valid — Delphy's posterior for rootHeight should cover the true value
  if inference is correct.

- **Non-independence of tip-date ranks:**  Ranks from tips within
  the same replicate share the same posterior chain (same `rootHeight`
  and `age(root)` samples), so they are not fully independent.  The
  effective sample size for the ECDF confidence band lies between
  N_replicates (200) and N_total_tips (~6000 month / ~2000 year).
  The plotted confidence band (based on N_total_tips) will be slightly
  too narrow, but this is a known limitation and the visual diagnostic
  is still useful.

- **Seed stability:**  Per the Sapling tip-date uncertainty plan,
  tip-date masking draws from the RNG *after* missing data simulation.
  This means the coalescent tree, mutations, and missing data pattern
  are identical to what `09_tight_priors_with_alpha` would produce
  with the same seed.  Only the date masking differs.
