# Plan: Final WCSS — Exponential without Year-Uncertain Tips

## Goal

This is a variant of study 23 (`23_final_exponential`) with
year-long tip-date uncertainty disabled (`P_TIP_DATE_UNCERTAIN_UPTO_YEAR = 0.0`).
The purpose is to test whether removing year-uncertain tips improves
rootHeight coverage in the exponential population model, paralleling
the 22c study for the skygrid model.

## Changes from study 23

| Parameter | Study 23 | This variant |
|-----------|----------|-------------|
| P_TIP_DATE_UNCERTAIN_UPTO_YEAR | 0.05 | 0.0 |
| DEFAULT_MASTER_SEED | 2025 | 9999 |

All other parameters, priors, and settings are identical to study 23.
Month-level tip-date uncertainty is kept at 0.15.

## Configuration summary

| Setting | Value |
|---------|-------|
| NUM_TIPS | 500 |
| NUM_SITES | 30,000 |
| DEFAULT_N | 200 |
| DEFAULT_STEPS | 3,000,000,000 |
| N0_PRIOR_MEAN | 2.5 year |
| N0_PRIOR_STDDEV | 0.5 year |
| G_PRIOR_MU | 2.0 year^{-1} |
| G_PRIOR_SCALE | 0.2 |
| G_MIN | 0.5 |
| MU_PRIOR_MEAN | 1e-3 subst/site/year |
| MU_PRIOR_STDDEV | 1e-4 subst/site/year |
| ALPHA_PRIOR_MEAN | 1.0 |
| P_TIP_DATE_UNCERTAIN_UPTO_MONTH | 0.15 |
| P_TIP_DATE_UNCERTAIN_UPTO_YEAR | 0.0 |
| DEFAULT_MASTER_SEED | 9999 |

## Implementation

The scripts (`01_generate.py`, `02_run.py`, `03_analyze.py`,
`04_plot.py`) are copied from study 23 with
`P_TIP_DATE_UNCERTAIN_UPTO_YEAR` changed from 0.05 to 0.0 and
`DEFAULT_MASTER_SEED` changed from 2025 to 9999.
