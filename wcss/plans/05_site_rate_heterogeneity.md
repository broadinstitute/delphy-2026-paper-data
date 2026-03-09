# Plan: WCSS with Site-Rate Heterogeneity

## Goal

Extend `04_free_exp_pop` to also validate inference of the site-rate
heterogeneity parameter `alpha` (the Gamma shape parameter).  In
`04_free_exp_pop`, site-rate heterogeneity was disabled (all sites had
rate 1.0); here, each replicate draws `alpha` from its prior and uses
it to generate Gamma-distributed per-site rate modifiers in Sapling.
Delphy then infers the posterior distribution of `alpha` alongside all
other parameters.

This uses:

- **Sapling's `--site-rate-heterogeneity-alpha` option** (commits
  `8c7a0bf` and `50e83d6`), which draws `nu_l ~ Gamma(alpha, alpha)`
  for each site and uses thinning + alias sampling for mutations.
- **Delphy's `--v0-site-rate-heterogeneity` flag**, which enables the
  alpha move (scaling proposals) and Gibbs sampling of per-site rates
  `nu_l`.

## Delphy's prior on alpha

Verified in `core/run.cpp` (lines 485-487 and 1177-1184):

```
alpha ~ Exponential(mean = 1.0)
```

i.e., `pi(alpha) = exp(-alpha)` for `alpha > 0`.  The `nu_l` are then
drawn as `nu_l ~ Gamma(shape=alpha, rate=alpha)` with mean 1 and
variance `1/alpha`.

## Differences from `04_free_exp_pop`

1. **`alpha` is drawn from the prior** per replicate:
   `alpha ~ Exponential(mean=1.0)`.
2. **Sapling is called with `--site-rate-heterogeneity-alpha`** set to
   the sampled value of `alpha`.  Sapling internally draws
   `nu_l ~ Gamma(alpha, alpha)` for each site.
3. **Delphy is called with `--v0-site-rate-heterogeneity`** to enable
   inference of `alpha` and `nu_l`.
4. **`alpha` is added to analysis and plotting** alongside mu, n0, g,
   kappa, pi, and rootHeight.
5. **`alpha` is added to ESS checking** (it is no longer fixed).

Everything else (tips, genome size, population model priors, mu prior,
clade coverage, etc.) is identical to `04_free_exp_pop`.

## Configuration

- **Directory:** `wcss/05_site_rate_heterogeneity/`
- **Tips:** Same as before (200 tips, dates uniform over 2025)
- **N (number of simulation replicates):** 10 (debug), then 200 (final)
- **Steps:** 1,000,000,000 (same as `04_free_exp_pop`; may need to
  increase if alpha mixing is slow)
- **Sampled parameters (drawn fresh per replicate):**
  - Site-rate heterogeneity alpha ~ Exponential(mean=1.0)
  - Effective population size n0 ~ InvGamma(mean=3 years, stddev=1 year)
  - Exponential growth rate g ~ Exponential(mean=1 e-folding/year),
    constrained to g >= 0
  - Mutation rate mu ~ Gamma(alpha=1, beta=1000)
    (i.e., Exponential with mean 1e-3 subst/site/year)
  - Stationary frequencies (pi_A, pi_C, pi_G, pi_T) ~ Dirichlet(1,1,1,1)
  - HKY kappa ~ LogNormal(meanlog=1.0, sdlog=1.25)
- **Fixed parameters:**
  - Genome size L = 30,000 sites

## Deliverables

Files in `wcss/05_site_rate_heterogeneity/`:

0. **`00_plan.md`** -- Symlink to `../plans/05_site_rate_heterogeneity.md`
1. **`01_generate.py`** -- Generate simulation inputs and Makefile
2. **`02_run.py`** -- Run Delphy via `make -jN`
3. **`03_analyze.py`** -- Run loganalyser, check ESS, compute coverage/ranks
4. **`04_plot.py`** -- Produce all plots from TSV files

---

## Part 1: `01_generate.py` -- Generate

Copy from `04_free_exp_pop/01_generate.py` with these changes:

### Configuration changes

- **Add:**
  - `ALPHA_PRIOR_MEAN = 1.0` (Exponential prior on alpha)

### Step 2 changes: Sample parameters

- `sample_prior_params(rng)` now also draws `alpha`:
  - `alpha ~ Exponential(mean=1.0)`:
    In numpy: `alpha = rng.exponential(1.0)`.
- Returns `n0, g, mu, kappa, pi, alpha`.

### Step 2 changes: Run Sapling

- `run_sapling()` gains an `alpha` parameter and passes
  `--site-rate-heterogeneity-alpha {alpha}` to sapling.
- The print statement includes `alpha` in its output.

### Step 3 changes: Makefile

The Delphy invocation in the Makefile pattern rule adds:
- `--v0-site-rate-heterogeneity`

(This is a boolean flag; no value needed.)

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
	  --v0-site-rate-heterogeneity \
	&& touch $@
```

---

## Part 2: `02_run.py` -- Run

Identical to `04_free_exp_pop/02_run.py`.  Copy verbatim.

---

## Part 3: `03_analyze.py` -- Analyze

Copy from `04_free_exp_pop/03_analyze.py` with these changes:

### Configuration changes

- Add `("alpha", "alpha")` to the `PARAMS` list, after `mu` and
  before `n0`:
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
- `ESS_IGNORE` remains `set()` (all parameters are inferred).

### Step 2 changes: Read true parameters

- `read_true_params()` additionally reads `alpha` from
  `info['subst_model']`:
  ```python
  "alpha": info["subst_model"]["site_rate_heterogeneity_alpha"],
  ```
  (Sapling writes `{"site_rate_heterogeneity_alpha": ...}` in the
  `subst_model` object when alpha > 0.)

### All other steps

No structural changes needed.  Coverage, ranks, and clade coverage all
operate on the `PARAMS` list, so adding the new entry automatically
includes alpha in all downstream analyses.

---

## Part 4: `04_plot.py` -- Plot

Copy from `04_free_exp_pop/04_plot.py` with these changes:

### Configuration changes

- Add `("alpha", "alpha")` to the `PARAMS` list (same position as
  in `03_analyze.py`).

### All other code

No structural changes.  The summary figure gains one extra row for
`alpha`.

---

## Directory structure at runtime

```
wcss/05_site_rate_heterogeneity/
  00_plan.md -> ../plans/05_site_rate_heterogeneity.md
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
    ecdf_alpha.pdf
    ecdf_n0.pdf
    ecdf_g.pdf
    ecdf_kappa.pdf
    ...
    rank_histogram_mu.pdf
    rank_histogram_alpha.pdf
    rank_histogram_n0.pdf
    rank_histogram_g.pdf
    rank_histogram_kappa.pdf
    ...
    scatter_mu.pdf
    scatter_alpha.pdf
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
cd wcss/05_site_rate_heterogeneity

# --- Final run (200 replicates, 1B steps -- the defaults) ---
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

- **ESS for alpha:** The alpha move uses 10 scaling sub-moves per MCMC
  step, which should mix reasonably well.  However, alpha controls
  30,000 latent `nu_l` variables, and the Gibbs resampling of all `nu_l`
  after each alpha move could slow things down.  If ESS for `alpha`
  is consistently low, increase steps to 2B.

- **Very small alpha values:** Since `alpha ~ Exponential(mean=1)`,
  small values like `alpha < 0.1` are plausible (~10% probability).
  Very small alpha means extreme site-rate variation (variance =
  1/alpha), which could make inference harder and produce very long
  tails in the posterior.  This is expected behavior and the WCSS should
  still be well-calibrated.

- **Alpha naming:** Delphy's log output uses `alpha` as the column name
  (while BEAST XML output uses `gammaShape`).  Our display name matches
  the log column, so no mapping is needed.
