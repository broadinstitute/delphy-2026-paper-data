# Plan: Final WCSS Studies — Overview

## Goal

The preparatory WCSS studies (02–11) incrementally validated Delphy's
inference by adding one feature at a time.  The three final studies
(21–23) are the publishable validation runs.  They use larger trees,
per-replicate tip dates, and longer chains, and they cover the two
population models that Delphy supports (exponential and Skygrid) plus
a high-mutation stress test.

## Studies

| Study | Pop model  | Alpha | Missing data | Tip-date uncertainty | Mutation regime | Based on |
|-------|-----------|-------|-------------|---------------------|----------------|----------|
| 21    | Constant  | No    | No          | No                  | High (wide mu prior) | 03 |
| 22    | Skygrid   | Yes   | Yes         | Yes                 | Densely sampled (~1 mut/branch) | 10 + 11 |
| 23    | Exponential | Yes | Yes         | Yes                 | Densely sampled (~1 mut/branch) | 10 + 11 |

### Study 21 — High-mutation stress test

Like `03_free_mu`: constant population (fixed n0 = 1.0 year), free
mu with a wide Exponential(mean=1e-3) prior, HKY with free kappa and
pi.  No site-rate heterogeneity, no missing data, no tip-date
uncertainty.

With 500 tips, L = 30,000, and mu drawn from a wide prior, many
replicates will have far more mutations than the ~2 × 500 = 1000
expected for a typical densely sampled genomic epi dataset.  This
study validates that Delphy's inference remains well-calibrated even
in high-mutation regimes where branches carry many mutations.

### Studies 22 & 23 — Full-featured validation

These two studies are identical except for the population model.  Both
include all key features that Delphy supports:

- Site-rate heterogeneity (free alpha)
- Missing data (~5% per tip)
- Tip-date uncertainty (15% month-uncertain, 5% year-uncertain)
- Free substitution model (HKY with free mu, kappa, pi)

The priors on mu and the population model parameters are tuned so that
the distribution of total mutation counts is tightly concentrated
around 2 × 500 = 1000, giving approximately 1 mutation per branch on
average.  This is the sweet spot where we expect Delphy's
reformulation of Bayesian phylogenetics to be most effective.

**Study 22** uses a **Skygrid** population model (5 knots, log-linear
interpolation), building on study 10's implementation but adding
tip-date uncertainty.

**Study 23** uses an **exponential** population model (free n0 and g),
building on study 11's implementation.

## Common changes from preparatory studies

### 1. More tips: 500 instead of 200

```python
NUM_TIPS = 500
```

Larger trees stress-test mixing and increase statistical power for
the WCSS calibration check.

### 2. Per-replicate tip dates

In studies 02–11, a single `tips.txt` was generated once and shared
across all replicates.  In the final studies, each replicate gets its
own tip dates, generated from the same distribution (uniform over
2025) but with independent randomness.

This eliminates a subtle conditioning artifact: sharing tip dates
means all replicates have the same tip-date configuration, so the
WCSS only validates calibration conditional on that specific set of
dates.  Per-replicate tip dates validate calibration marginally over
the tip-date distribution.

Implementation: `generate_tips()` moves inside the per-replicate loop
and writes to `sim_NNN/tips.txt` instead of `sims/tips.txt`.  The
`run_sapling()` call uses the per-replicate tips file.

### 3. Longer chains

With 500 tips, the state space is larger and mixing is slower.  All
final studies use at least 2× the chain length of their preparatory
counterparts:

| Study | Steps | Rationale |
|-------|-------|-----------|
| 21    | 2,000,000,000 | 4× study 03's 500M (2.5× more tips, wider mu) |
| 22    | 3,000,000,000 | 3× study 10's 1B (2.5× more tips, added tip-date uncertainty) |
| 23    | 3,000,000,000 | 3× study 11's 1B (2.5× more tips) |

These are starting points.  If ESS is consistently low, increase
further.

### 4. Same number of replicates

```python
DEFAULT_N = 200
```

Same as the preparatory studies.  200 replicates give good resolution
on the WCSS calibration.

## Tuning mutation counts for studies 22 & 23

The mutation rate prior will be a tight Gamma centered on
mu = 1e-3 subst/site/year with a small standard deviation (~1e-4),
reflecting a real virus like SARS-CoV-2.  The population model priors
will be tuned so that, in combination with this mu prior, the
distribution of total mutation counts is concentrated around
2 × 500 = 1000.  The exact prior parameters will be determined by
short pilot runs and documented in the individual study plans.

## Directory structure

```
wcss/
  plans/
    20_final_wcss_overview.md    (this file)
    21_final_high_mutation.md
    22_final_skygrid.md
    23_final_exponential.md
  21_final_high_mutation/
    00_plan.md -> ../plans/21_final_high_mutation.md
    01_generate.py
    02_run.py
    03_analyze.py
    04_plot.py
  22_final_skygrid/
    00_plan.md -> ../plans/22_final_skygrid.md
    01_generate.py
    02_run.py
    03_analyze.py
    04_plot.py
  23_final_exponential/
    00_plan.md -> ../plans/23_final_exponential.md
    01_generate.py
    02_run.py
    03_analyze.py
    04_plot.py
```

## Self-contained plans

Each individual study plan (21, 22, 23) is entirely self-contained.
It may reference the preparatory study it evolved from for context,
but a reader should be able to understand exactly what the study does
— including all prior distributions, model parameters, and Delphy CLI
flags — from the individual plan alone.

## Execution order

1. **Study 21** first — simplest model, quickest to validate.
2. **Studies 22 & 23** in parallel — they are independent.

For each study:
```bash
cd wcss/2N_final_*/

# Debug run (10 replicates, short chains)
./01_generate.py --n 10 --steps <short>
./02_run.py
./03_analyze.py --n 10 --ignore-low-ess --force-include-all-replicates
./04_plot.py

# Final run (200 replicates, full chains)
./01_generate.py
./02_run.py
./03_analyze.py
./04_plot.py
```
