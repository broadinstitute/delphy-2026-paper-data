# Plan: WCSS with Free Exponential Population Model

## Goal

Extend `03_free_mu` to also validate inference of the exponential
population model parameters `n0` (effective population size at the time
of the last tip) and `g` (exponential growth rate).  In `03_free_mu`,
the population model was constant with fixed size N_e = 1.0 year; here,
we use an exponential population curve `N(t) = n0 * exp(g * (t - t0))`
with both `n0` and `g` sampled from their priors.

This uses two new Delphy features:

- **Flexible Inverse-Gamma prior on n0** (`--v0-pop-n0-prior-mean` /
  `--v0-pop-n0-prior-stddev`), from plan `2026-03-05-01-flexible-n0-prior.md`.
- **Flexible Laplace prior on g with bounds** (`--v0-pop-g-prior-exponential-with-mean`),
  from plan `2026-03-05-03-flexible-g-prior.md`.

## Differences from `03_free_mu`

1. **Population model changes from constant to exponential** with both
   parameters sampled.
2. **`n0` is drawn from the prior** per replicate:
   `n0 ~ InvGamma(alpha, beta)` with mean 3 years, stddev 1 year.
   (This implies alpha = 11, beta = 30.)
3. **`g` is drawn from the prior** per replicate:
   `g ~ Exponential(mean=1)` (in e-foldings/year), constrained to g >= 0.
   In Delphy, this is a truncated Laplace with mu=0, scale=1, g_min=0,
   configured via `--v0-pop-g-prior-exponential-with-mean 1`.
4. **`n0` and `g` are added to the analysis and plotting** alongside mu,
   kappa, pi, and rootHeight.
5. **`exponential.popSize` and `exponential.growthRate` are added to ESS checking** (they are
   no longer fixed).

Everything else (tips, genome size, mu prior, site rate heterogeneity,
clade coverage, etc.) is identical to `03_free_mu`.

## Configuration

- **Directory:** `wcss/04_free_exp_pop/`
- **Tips:** Same as before (200 tips, dates uniform over 2025)
- **N (number of simulation replicates):** 10 (debug), then 200 (final)
- **Steps:** 100,000 per tip = 20,000,000 total (debug), then
  5,000,000 per tip = 1,000,000,000 total (final; doubled from the
  initial 500M because the extra free pop parameters slow mixing)
- **Sampled parameters (drawn fresh per replicate):**
  - Effective population size n0 ~ InvGamma(mean=3 years, stddev=1 year)
  - Exponential growth rate g ~ Exponential(mean=1 e-folding/year),
    constrained to g >= 0
  - Mutation rate mu ~ Gamma(alpha=1, beta=1000)
    (i.e., Exponential with mean 1e-3 subst/site/year)
  - Stationary frequencies (pi_A, pi_C, pi_G, pi_T) ~ Dirichlet(1,1,1,1)
  - HKY kappa ~ LogNormal(meanlog=1.0, sdlog=1.25)
- **Fixed parameters:**
  - No site rate heterogeneity
  - Genome size L = 30,000 sites

## Deliverables

Files in `wcss/04_free_exp_pop/`:

0. **`00_plan.md`** -- Symlink to `../plans/04_free_exp_pop.md` (this plan)
1. **`01_generate.py`** -- Generate simulation inputs and Makefile
2. **`02_run.py`** -- Run Delphy via `make -jN`
3. **`03_analyze.py`** -- Run loganalyser, check ESS, compute coverage/ranks
4. **`04_plot.py`** -- Produce all plots from TSV files

---

## Part 1: `01_generate.py` -- Generate

Copy from `03_free_mu/01_generate.py` with these changes:

### Configuration changes

- **Remove:** `CONST_POP_N0 = 1.0`
- **Add:**
  - `N0_PRIOR_MEAN = 3.0` (years)
  - `N0_PRIOR_STDDEV = 1.0` (years)
  - `G_PRIOR_MEAN = 1.0` (e-foldings/year; Exponential prior on g >= 0)
- `DEFAULT_STEPS = 1_000_000_000` (5,000,000 per tip x 200 tips)

### Step 2 changes: Sample parameters

- `sample_prior_params(rng)` now also draws `n0` and `g`:
  - `n0 ~ InvGamma(mean=3, stddev=1)`:
    Compute alpha = 2 + mean^2/var = 11, beta = mean*(alpha-1) = 30.
    Sample: `inv_n0 = rng.gamma(alpha, 1.0/beta)`, then `n0 = 1.0/inv_n0`.
    In numpy: `n0 = 1.0 / rng.gamma(11.0, 1.0/30.0)`.
  - `g ~ Exponential(mean=1)`:
    In numpy: `g = rng.exponential(1.0)`.
- `run_sapling()` uses `--exp-pop-n0 {n0} --exp-pop-g {g}` instead of
  `--const-pop-n0 1.0`.
- The print statement includes `n0` and `g` in its output.

### Step 3 changes: Makefile

The Delphy invocation in the Makefile pattern rule changes:
- **Remove:**
  - `--v0-init-final-pop-size 1.0`
  - `--v0-fix-final-pop-size`
  - `--v0-init-pop-growth-rate 0.0`
  - `--v0-fix-pop-growth-rate`
- **Add:**
  - `--v0-pop-n0-prior-mean 3 --v0-pop-n0-prior-stddev 1`
  - `--v0-pop-g-prior-exponential-with-mean 1`

The `log_every` and `tree_every` intervals are computed as
`steps // 1000` (giving 1000 samples per run).

Resulting Makefile rule (with default 1B steps):
```makefile
sim_%/.done: sim_%/sim.maple
	$(DELPHY) \
	  --v0-in-maple $< \
	  --v0-steps 1000000000 \
	  --v0-out-log-file sim_$*/delphy.log \
	  --v0-log-every 1000000 \
	  --v0-out-trees-file sim_$*/delphy.trees \
	  --v0-tree-every 1000000 \
	  --v0-out-delphy-file sim_$*/delphy.dphy \
	  --v0-delphy-snapshot-every 1000000 \
	  --v0-mu-prior-alpha 1 \
	  --v0-mu-prior-beta 1000 \
	  --v0-pop-n0-prior-mean 3 \
	  --v0-pop-n0-prior-stddev 1 \
	  --v0-pop-g-prior-exponential-with-mean 1 \
	&& touch $@
```

---

## Part 2: `02_run.py` -- Run

Identical to `03_free_mu/02_run.py`.  Copy verbatim.

---

## Part 3: `03_analyze.py` -- Analyze

Copy from `03_free_mu/03_analyze.py` with these changes:

### Configuration changes

- Add `("n0", "exponential.popSize")` and `("g", "exponential.growthRate")` to the `PARAMS`
  list, after `mu` and before `kappa`:
  ```python
  PARAMS = [
      ("mu",         "meanRate"),
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
- `ESS_IGNORE` remains `set()` (all parameters are inferred).

### Step 2 changes: Read true parameters

- `read_true_params()` additionally reads `n0` and `g` from
  `info['pop_model']`:
  ```python
  "n0": info["pop_model"]["n0"],
  "g": info["pop_model"]["g"],
  ```
  (Sapling writes `{"type": "exp", "n0": ..., "g": ...}` for
  exponential pop models.)

### New option: `--force-include-all-replicates`

`03_analyze.py` adds a `--force-include-all-replicates` flag that
skips all replicate exclusion based on ESS.  This is needed for debug
runs where ESS is expected to be low across the board.  Without it,
short debug runs can exclude all replicates, causing downstream
crashes.

### All other steps

No structural changes needed.  Coverage, ranks, and clade coverage all
operate on the `PARAMS` list, so adding the two new entries
automatically includes them in all downstream analyses.

---

## Part 4: `04_plot.py` -- Plot

Copy from `03_free_mu/04_plot.py` with these changes:

### Configuration changes

- Add `("n0", "exponential.popSize")` and `("g", "exponential.growthRate")` to the `PARAMS`
  list (same positions as in `03_analyze.py`).

### All other code

No structural changes.  The summary figure gains two extra rows for
`n0` and `g`.

---

## Directory structure at runtime

```
wcss/04_free_exp_pop/
  00_plan.md -> ../plans/04_free_exp_pop.md
  01_generate.py
  02_run.py
  03_analyze.py
  04_plot.py
  analyses/
    clade_coverage.tsv
    clade_coverage_raw.tsv
    coverage_summary.txt
    ess_check.tsv
    excluded_replicates.tsv
    loganalyser_output.tsv
    loganalyser_output_filtered.tsv
    ranks.tsv
    true_params.tsv
  plots/
    clade_coverage.pdf
    ecdf_mu.pdf
    ecdf_n0.pdf
    ecdf_g.pdf
    ecdf_kappa.pdf
    ...
    rank_histogram_mu.pdf
    rank_histogram_n0.pdf
    rank_histogram_g.pdf
    rank_histogram_kappa.pdf
    ...
    scatter_mu.pdf
    scatter_n0.pdf
    scatter_g.pdf
    scatter_kappa.pdf
    ...
    wcss_summary.pdf
  sims/
    Makefile
    tips.txt
    sim_000/
      sim.maple
      sim_info.json
      sim.nwk
      sim.nexus
      sim.fasta
      delphy.log
      delphy.trees
      delphy.dphy
      .done
    sim_001/
      ...
    ...
```

## Execution workflow

```bash
cd wcss/04_free_exp_pop

# --- Final run (200 replicates, 500M steps — the defaults) ---
./01_generate.py
./02_run.py
./03_analyze.py
./04_plot.py

# If step 3 fails due to low ESS, increase the number of steps and rerun:
./01_generate.py --steps 2000000000
make -C sims clean && ./02_run.py
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

- **ESS for pop parameters:** Mixing for `n0` and `g` is slower than
  for substitution model parameters.  An initial run at 500M steps
  showed 10-15% of replicates with low ESS on a few observables, so
  the default was doubled to 1B.  If ESS is still low, increase
  further to 2B.

- **Correlation between n0 and g:** With an exponential population
  model, `n0` and `g` can be correlated in the posterior.  This is
  expected and doesn't affect WCSS calibration — the WCSS validates
  marginal coverage for each parameter independently, and a correctly
  implemented sampler will be well-calibrated regardless of posterior
  correlations.  However, strong correlations may slow mixing, which
  would show up as low ESS.

- **Fixed tip dates:** We provide fixed tip times via `--tip-file`
  (uniform over 2025), the same as in previous studies.  This is
  correct because Delphy's posterior treats tip dates as fixed data,
  not as samples from a prior.