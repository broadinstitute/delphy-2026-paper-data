# Plan: Final WCSS — Skygrid without Year-Uncertain Tips

## Goal

This is a variant of study 22b (`22b_final_skygrid_lower_n0`) with
year-long tip-date uncertainty disabled (`P_TIP_DATE_UNCERTAIN_UPTO_YEAR = 0.0`).
The purpose is to test whether the low rootHeight coverage (0.86) in
study 22b is caused by year-uncertain tips being misinferred early,
which pushes the root height too high.

## Changes from study 22b

| Parameter | Study 22b | This variant |
|-----------|-----------|-------------|
| P_TIP_DATE_UNCERTAIN_UPTO_YEAR | 0.05 | 0.0 |

All other parameters, priors, and settings are identical to study 22b.
Month-level tip-date uncertainty is kept at 0.15.

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
| TAU_PRIOR_ALPHA | 36.0 |
| TAU_PRIOR_BETA | 24.0 |
| MU_PRIOR_MEAN | 1e-3 subst/site/year |
| MU_PRIOR_STDDEV | 1e-4 subst/site/year |
| ALPHA_PRIOR_MEAN | 1.0 |
| P_TIP_DATE_UNCERTAIN_UPTO_MONTH | 0.15 |
| P_TIP_DATE_UNCERTAIN_UPTO_YEAR | 0.0 |
| DEFAULT_MASTER_SEED | 2025 |

## Implementation

The scripts (`01_generate.py`, `02_run.py`, `03_analyze.py`,
`04_plot.py`) are copied from study 22b with only
`P_TIP_DATE_UNCERTAIN_UPTO_YEAR` changed from 0.05 to 0.0.
Because the Sapling seed sequence is the same but the tip-date
uncertainty parameter differs, the simulated trees and sequences
will NOT be identical to 22b (different tips will be marked as
uncertain, and the resulting coalescent simulations will differ).
