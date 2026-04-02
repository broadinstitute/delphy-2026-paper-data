# Well-Calibrated Simulation Studies (WCSS) for Delphy

This directory contains a series of well-calibrated simulation studies
that validate Delphy's MCMC inference, following the methodology of
Mendes, Bouckaert, Carvalho & Drummond (2025).  Each study simulates
many datasets by drawing parameters from Delphy's priors and
generating sequence data with Sapling, then runs Delphy on each
dataset and checks that the true generating values are recovered at
the expected rates (~95% coverage for 95% HPD intervals, uniform rank
distributions, well-calibrated clade probabilities).  See
[plans/01_intro.md](plans/01_intro.md) for background on the
methodology.

The studies are organized as a progression, where each one builds on
the previous by relaxing a restriction or adding a model feature.

Studies 02-07 use wide default priors, which produce a broad range of
mutation counts per replicate — sometimes hundreds of mutations on
single branches, well outside Delphy's sweet spot.  Starting from 08,
tighter priors on mu and g concentrate runs around densely sampled
trees (0-2 mutations per branch), which is the regime Delphy targets.
The earlier wide-prior studies are retained to demonstrate correctness
in that regime.  Study 10 (Skygrid) branches off from 06 since it
uses a different population model.

## Common setup

All studies share the same basic design:

- **200 replicates** per study
- **200 tips** (studies 02-11) or **500 tips** (studies 21-23),
  dates uniform over 2025, genome size L = 30,000 sites
- **HKY substitution model** with kappa ~ LogNormal(1.0, 1.25) and
  pi ~ Dirichlet(1,1,1,1)
- **4 scripts per study:**
  1. `01_generate.py` — draw parameters from priors, simulate data
     with Sapling, generate a Makefile for Delphy runs
  2. `02_run.py` — run Delphy via `make -jN`
  3. `03_analyze.py` — run loganalyser, check ESS, compute coverage,
     ranks, and clade coverage
  4. `04_plot.py` — produce all figures from analysis TSVs

Detailed plans for each study are in `plans/`.

## Cloud infrastructure

Each study runs 200 replicates of Delphy, which can take many CPU-hours.
Three scripts (`cloud_setup.py`, `cloud_run.py`, `cloud_teardown.py`)
automate offloading runs to Hetzner Cloud instances.  `cloud_setup.py`
provisions an instance and uploads sim inputs; `cloud_run.py` launches
`make -jN` via a detached process and monitors progress with periodic
rsync of results; `cloud_teardown.py` deletes the instance.  Work can
be split across multiple instances by partitioning the sim range.  See
[plans/cloud-runs.md](plans/cloud-runs.md) for details.

## Studies

### 02_simple — Baseline

Fix mu, population size, and growth rate.  Only infer kappa, pi, and
tree topology.  Establishes the WCSS infrastructure.

### 03_free_mu — Free mutation rate

Adds: **mu ~ Gamma(1, 1000)** (Exponential with mean 1e-3).  Delphy
now also infers the mutation rate.

### 04_free_exp_pop — Free exponential population model

Adds: **n0 ~ InvGamma(mean=3, stddev=1)** and
**g ~ Exponential(mean=1)**.  Population model changes from constant
to exponential with both parameters inferred.  Steps doubled to 1B
due to slower mixing.

### 05_site_rate_heterogeneity — Gamma site rates

Adds: **alpha ~ Exponential(mean=1)**.  Sapling generates
per-site rate modifiers nu_l ~ Gamma(alpha, alpha); Delphy infers
alpha and the latent nu_l.

### 06_missing_data — Missing data (full model)

Adds: **simulated missing data** (~5% per tip: 3 gaps of mean length
500, plus 3 isolated missing sites).  All model features from 05 are
retained.  Validates that Delphy handles N-masked sites correctly.

### 07_missing_data_no_alpha — Missing data (no site-rate heterogeneity)

A variant of 06 with **site-rate heterogeneity disabled**.  Same
missing data as 06, but the model matches 04_free_exp_pop.  Motivated
by the observation that with 200 tips and L=30000, the data may be
nearly uninformative about alpha, making it hard to validate.

### 08_tight_priors — Tighter priors on mu and g

Replaces the wide default priors on mu and g with **tighter,
configurable priors** to keep mutation counts near the target of
2 x 200 = 400.  Chosen values after tuning:
mu ~ Gamma(mean=5e-4, stddev=1e-4),
g ~ Laplace(mu=2, scale=0.2, g_min=0.5).
Missing data included, site-rate heterogeneity disabled.
Adds mutation count diagnostics (histogram and eCDF) to
`01_generate.py`.

### 09_tight_priors_with_alpha — Tight priors + site-rate heterogeneity

Combines 08's tighter mu/g priors with **alpha ~ Exponential(mean=1)**
and site-rate heterogeneity enabled.  Validates alpha inference
alongside the tighter priors.

### 10_skygrid — Skygrid population model

Replaces the exponential population model with a **Skygrid** (5
knots, log-linear interpolation).  Adds tau (GMRF precision) and
gamma_0..gamma_4 (log-population at each knot) as inferred
parameters, plus N_bar as a derived parameter.  Branches off from 06
(missing data + site-rate heterogeneity + wide priors).

### 11_tip_date_uncertainty — Uncertain tip dates

Adds: **tip-date uncertainty** (15% of tips masked to month precision,
5% to year precision).  Delphy samples exact tip dates during MCMC
using a bounded exponential proposal matched to the branch likelihood.
Builds on 09 (tight priors + site-rate heterogeneity + missing data).
Validates tip-date posterior coverage and rank uniformity separately
for month-uncertain and year-uncertain tips.

### Final studies (21-23): 500 tips

Studies 21-23 scale up to 500 tips with per-replicate tip dates and
longer chains.  They represent the final validation configurations
intended for the paper.

### 21_final_high_mutation — High-mutation stress test

Uses **500 tips** with a **wide prior on mu** (Exponential, mean
1e-3), so many replicates produce far more mutations than the ~1000
expected for a densely sampled tree.  The model is deliberately
simple — constant population, no site-rate heterogeneity, no missing
data, no tip-date uncertainty — to isolate the effect of the
high-mutation regime on Delphy's tree rearrangement moves.

### 22_final_skygrid — Skygrid with all key features

Uses **500 tips** with a **Skygrid population model** (5 knots,
log-linear interpolation, GMRF precision tau, N_bar ~ InvGamma with
mean 0.75 year) combined with tight mu/g priors, site-rate
heterogeneity, missing data, and tip-date uncertainty.

### 22b_final_skygrid_lower_n0 — Skygrid with lower N_bar

A variant of study 22 with **N_bar prior mean 0.25 year** (vs 0.75).
Reduces deep-branch mutations to improve MCMC convergence.  Study 22
excluded 29.5% of replicates due to low ESS; this variant achieves
better mixing within the same step budget.

### 22c_skygrid_no_year_uncertainty — Skygrid without year-uncertain tips

A variant of study 22b with **year-level tip-date uncertainty disabled**
(P_TIP_DATE_UNCERTAIN_UPTO_YEAR = 0.0).  Tests whether year-uncertain
tips cause the low rootHeight coverage seen in 22b.  Month-level
uncertainty is kept at 0.15.

### 23_final_exponential — Exponential with all key features

Uses **500 tips** with an **exponential population model** (n0 ~
InvGamma, g ~ truncated Laplace) combined with tight mu/g priors,
site-rate heterogeneity, missing data, and tip-date uncertainty.
Companion to study 22, differing only in the population model.

### 23b_exponential_no_year_uncertainty — Exponential without year-uncertain tips

A variant of study 23 with **year-level tip-date uncertainty disabled**,
matching the same change made in study 22c for the skygrid model.
Month-level uncertainty is kept at 0.15.

## Compact summary plots

`plot_compact_summary.py` generates a single-page PDF per study with
three panels: (A) 95% HPD coverage bar chart, (B) clade coverage
calibration, and (C) per-parameter scatter plots with rank histograms
and eCDFs.  Run `99_plot_all_compact_summaries.sh` to regenerate all
plots.

## Files not in this repository

The following file types are excluded from git to keep the repository
within GitHub's size limits:

- **`.log`** and **`.trees`** — MCMC trace logs and sampled tree
  files.  These will be released via Zenodo alongside the paper.

`.dphy` snapshot files and `.fasta` sequence alignments are no longer
produced (removed to save ~0.5 TiB of storage).

See `.gitignore` for the exact exclusion patterns.
