# Plan: Final WCSS — Skygrid with Lower N_bar (Variant of Study 22)

## Goal

This is a variant of study 22 (`22_final_skygrid`) with a
substantially lower N_bar prior mean (0.25 year vs 0.75 year).  The
purpose is to improve MCMC convergence by reducing the number of
mutations on deep branches near the root, which is the tree property
that most slows down Delphy's MCMC mixing.

## Rationale

Study 22 with NBAR_PRIOR_MEAN=0.75 produced trees with a mutation
count distribution centered around 1000 (mean=1061, median=1006).
While this is close to the target of 2 x 500 = 1000, many replicates
had severe convergence problems:

- **59 of 200 replicates** (29.5%) were excluded due to very low ESS
  (< 150) on at least one parameter.
- The worst-affected observable was `numMuts`, with 58 replicates
  (29%) having ESS < 150.
- Other commonly affected observables included `alpha` (23 excluded),
  `rootHeight` (28 excluded), `treeLength` (21 excluded), and
  `likelihood_really_logG` (23 excluded).
- After exclusion, only 141 replicates remained for analysis.

The root cause is that larger N_bar values produce trees with longer
total branch length, and thus more mutations on branches near the
root.  These deep, heavily-mutated branches are where Delphy's SPR
moves are least effective, leading to poor mixing within the 3B step
budget.

## Changes from study 22

The only change is to the N_bar prior:

| Parameter | Study 22 | This variant |
|-----------|----------|-------------|
| NBAR_PRIOR_MEAN | 0.75 year | 0.25 year |
| NBAR_PRIOR_STDDEV | 0.15 year | 0.05 year |
| NBAR_PRIOR_ALPHA | 27.0 | 27.0 |
| NBAR_PRIOR_BETA | 19.5 year | 6.5 year |

The CV remains 0.2 (same relative precision).  All other parameters,
priors, and settings are identical to study 22.

## Results comparison

### Mutation counts

| Variant | Mutation count mean | Mutation count median |
|---------|--------------------:|----------------------:|
| Study 22 (N_bar=0.75) | 1,061 | 1,006 |
| This variant (N_bar=0.25) | 561 | 549 |

### Convergence

| Metric | Study 22 | This variant |
|--------|----------|-------------|
| Excluded replicates | 59 (29.5%) | 39 (19.5%) |
| Replicates used for analysis | 141 | 161 |
| numMuts ESS < 150 | 58 (29%) | 25 (12.5%) |
| alpha ESS < 150 | 23 (11.5%) | 3 (1.5%) |
| rootHeight ESS < 150 | 28 (14%) | 8 (4%) |
| treeLength ESS < 150 | 21 (10.5%) | 3 (1.5%) |
| Mean numMuts ESS | 387 | 479 |

The lower N_bar substantially reduces convergence failures across
all observables, though some replicates still have low ESS on
`numMuts` and `minTipESS_year`.

### Low-population barrier

With NBAR_PRIOR_MEAN = 0.25, the GMRF walk centers around
ln(0.25) ~ -1.39, leaving 4.5 log-units of margin to the barrier
at gamma = ln(1/365) ~ -5.90.  Monte Carlo estimation (10^7 samples)
of the full GMRF random walk with tau drawn from Gamma(36, 24) and
N_bar drawn from InvGamma(27, 6.5) finds 4 barrier hits
(P ~ 4×10⁻⁶), i.e. ~0.00 expected hits per 200 replicates.

## Configuration summary

| Setting | Value |
|---------|-------|
| NUM_TIPS | 500 |
| NUM_SITES | 30,000 |
| DEFAULT_N | 200 |
| DEFAULT_STEPS | 3,000,000,000 |
| NUM_KNOTS | 5 |
| SKYGRID_FIRST_KNOT_DATE | 2025-01-01 |
| SKYGRID_LAST_KNOT_DATE | 2026-01-01 |
| NBAR_PRIOR_MEAN | 0.25 year |
| NBAR_PRIOR_STDDEV | 0.05 year |
| NBAR_PRIOR_ALPHA | 27.0 |
| NBAR_PRIOR_BETA | 6.5 year |
| TAU_PRIOR_ALPHA | 36.0 |
| TAU_PRIOR_BETA | 24.0 |
| MU_PRIOR_MEAN | 1e-3 subst/site/year |
| MU_PRIOR_STDDEV | 1e-4 subst/site/year |
| ALPHA_PRIOR_MEAN | 1.0 |
| KAPPA_MEAN_LOG | 1.0 (hardcoded in Delphy) |
| KAPPA_SIGMA_LOG | 1.25 (hardcoded in Delphy) |
| MISSING_DATA_MEAN_NUM_GAPS | 3.0 |
| MISSING_DATA_MEAN_GAP_LENGTH | 500.0 |
| MISSING_DATA_MEAN_NUM_MISSING_SITES | 3.0 |
| P_TIP_DATE_UNCERTAIN_UPTO_MONTH | 0.15 |
| P_TIP_DATE_UNCERTAIN_UPTO_YEAR | 0.05 |
| TIP_DATE_START | 2025-01-01 |
| TIP_DATE_END | 2025-12-31 |
| DEFAULT_MASTER_SEED | 2025 |

## Delphy invocation

```
delphy \
  --v0-in-maple sim_NNN/sim.maple \
  --v0-steps 3000000000 \
  --v0-out-log-file sim_NNN/delphy.log \
  --v0-log-every 3000000 \
  --v0-out-trees-file sim_NNN/delphy.trees \
  --v0-tree-every 3000000 \
  --v0-out-delphy-file sim_NNN/delphy.dphy \
  --v0-delphy-snapshot-every 3000000 \
  --v0-mu-prior-mean 0.001 \
  --v0-mu-prior-stddev 0.0001 \
  --v0-pop-model skygrid \
  --v0-skygrid-num-parameters 5 \
  --v0-skygrid-first-knot-date 2025-01-01 \
  --v0-skygrid-last-knot-date 2026-01-01 \
  --v0-skygrid-type log-linear \
  --v0-skygrid-nbar-prior-mean 0.25 \
  --v0-skygrid-nbar-prior-stddev 0.05 \
  --v0-skygrid-infer-prior-smoothness \
  --v0-skygrid-tau-prior-alpha 36 \
  --v0-skygrid-tau-prior-beta 24 \
  --v0-skygrid-disable-low-pop-barrier \
  --v0-site-rate-heterogeneity
```

## Implementation

The scripts (`01_generate.py`, `02_run.py`, `03_analyze.py`,
`04_plot.py`) are copied from study 22 with only the N_bar prior
constants changed.  The analysis and plotting pipelines are
identical.
