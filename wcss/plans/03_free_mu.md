# Plan: WCSS with Free Mutation Rate (Gamma Prior)

## Goal

Extend the simple WCSS (`02_simple`) to also validate inference of the
mutation rate `mu`.  In `02_simple`, `mu` was fixed at its true value;
here, we let Delphy infer it under a Gamma(alpha=1, beta=1000) prior,
i.e., an exponential distribution with mean 1e-3 subst/site/year.

This uses the new `--v0-mu-prior-alpha` / `--v0-mu-prior-beta` CLI
options added in Delphy commit `7ed312b` ("Add flexible Gamma prior on
mutation rate mu").

## Differences from `02_simple`

1. **Mutation rate is now a sampled parameter**, not fixed.
2. **`mu` is drawn from the prior** per replicate: `mu ~ Gamma(1, 1000)`
   (shape=1, rate=1000 in units of years, so mean=1e-3 subst/site/year).
3. **Delphy is not told `--v0-fix-mutation-rate`**; instead, it receives
   `--v0-mu-prior-alpha 1 --v0-mu-prior-beta 1000` to match the
   generating prior.
4. **`mu` is added to the analysis and plotting** alongside kappa, pi,
   and rootHeight.
5. **`meanRate` is removed from `ESS_IGNORE`** since it is now inferred.

Everything else (tips, genome size, population model, site rate
heterogeneity, clade coverage, etc.) is identical to `02_simple`.

## Configuration

- **Directory:** `wcss/03_free_mu/`
- **Tips:** Same as `02_simple` (200 tips, dates uniform over 2025)
- **N (number of simulation replicates):** 200
- **Steps:** 2,500,000 per tip = 500,000,000 total
- **Fixed parameters (same across all replicates):**
  - Effective population size N_e = 1.0 year (constant population model)
  - No site rate heterogeneity
  - Genome size L = 30,000 sites
- **Sampled parameters (drawn fresh per replicate):**
  - Mutation rate mu ~ Gamma(alpha=1, beta=1000)
    (i.e., Exponential with mean 1e-3 subst/site/year)
  - Stationary frequencies (pi_A, pi_C, pi_G, pi_T) ~ Dirichlet(1,1,1,1)
  - HKY kappa ~ LogNormal(meanlog=1.0, sdlog=1.25)

## Deliverables

Files in `wcss/03_free_mu/`:

0. **`00_plan.md`** — Symlink to `../plans/03_free_mu.md` (this plan)
1. **`01_generate.py`** — Generate simulation inputs and Makefile
2. **`02_run.py`** — Run Delphy via `make -jN`
3. **`03_analyze.py`** — Run loganalyser, check ESS, compute coverage/ranks
4. **`04_plot.py`** — Produce all plots from TSV files

---

## Part 1: `01_generate.py` — Generate

Copy from `02_simple/01_generate.py` with these changes:

### Configuration changes

- `MU = 1e-3` — **remove** the fixed constant.  Replace with:
  - `MU_PRIOR_ALPHA = 1.0`
  - `MU_PRIOR_BETA = 1000.0`  (rate parameter in years)
  This gives an Exponential(rate=1000) prior, i.e., mean = 1e-3.
- `DEFAULT_STEPS = 500_000_000` (2,500,000 per tip x 200 tips).

### Step 2 changes: Sample parameters

- `sample_prior_params(rng)` now also draws
  `mu ~ Gamma(shape=MU_PRIOR_ALPHA, scale=1/MU_PRIOR_BETA)`.
  In numpy: `rng.gamma(MU_PRIOR_ALPHA, 1.0 / MU_PRIOR_BETA)`.
- `run_sapling()` passes the sampled `mu` value via `--mu {mu}` instead
  of the fixed constant.
- The print statement includes `mu` in its output.

### Step 3 changes: Makefile

The Delphy invocation in the Makefile pattern rule changes:
- **Remove:** `--v0-fix-mutation-rate`
- **Remove:** `--v0-init-mutation-rate 1e-3` (Delphy now defaults to the
  prior mean when a proper Gamma prior is specified; commit `bf51b39`)
- **Add:** `--v0-mu-prior-alpha 1 --v0-mu-prior-beta 1000`
- Everything else (fix pop size, fix growth rate, log/tree/dphy output)
  stays the same.

Resulting Makefile rule:
```makefile
sim_%/.done: sim_%/sim.maple
	$(DELPHY) \
	  --v0-in-maple $< \
	  --v0-steps 500000000 \
	  --v0-out-log-file sim_$*/delphy.log \
	  --v0-log-every 500000 \
	  --v0-out-trees-file sim_$*/delphy.trees \
	  --v0-tree-every 500000 \
	  --v0-out-delphy-file sim_$*/delphy.dphy \
	  --v0-delphy-snapshot-every 500000 \
	  --v0-mu-prior-alpha 1 \
	  --v0-mu-prior-beta 1000 \
	  --v0-init-final-pop-size 1.0 \
	  --v0-fix-final-pop-size \
	  --v0-init-pop-growth-rate 0.0 \
	  --v0-fix-pop-growth-rate \
	&& touch $@
```

---

## Part 2: `02_run.py` — Run

Identical to `02_simple/02_run.py`.  Copy verbatim.

---

## Part 3: `03_analyze.py` — Analyze

Copy from `02_simple/03_analyze.py` with these changes:

### Configuration changes

- Add `("mu", "meanRate")` to the `PARAMS` list (first position, before
  kappa, so that mu appears first in all output tables and plots).
- Change `ESS_IGNORE` from `{"meanRate"}` to `set()` (empty set — no
  observables are expected to have degenerate ESS now).

### Step 2 changes: Read true parameters

- `read_true_params()` additionally reads `"mu"` from
  `info['subst_model']['mu']`.
- `true_params.tsv` gains a `mu` column:
  `replicate`, `mu`, `kappa`, `pi_A`, `pi_C`, `pi_G`, `pi_T`, `rootHeight`.

### All other steps

No structural changes needed.  Coverage, ranks, and clade coverage all
operate on the `PARAMS` list, so adding `("mu", "meanRate")` to `PARAMS`
automatically includes `mu` in all downstream analyses.

---

## Part 4: `04_plot.py` — Plot

Copy from `02_simple/04_plot.py` with these changes:

### Configuration changes

- Add `("mu", "meanRate")` to the `PARAMS` list (first position).

### All other code

No structural changes.  The summary figure gains one extra row for `mu`.
The zoomed inset for kappa continues to apply (the `if name == "kappa"`
check still works).

---

## Directory structure at runtime

```
wcss/03_free_mu/
  00_plan.md -> ../plans/03_free_mu.md
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
    ecdf_kappa.pdf
    ...
    rank_histogram_mu.pdf
    rank_histogram_kappa.pdf
    ...
    scatter_mu.pdf
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
cd wcss/03_free_mu

# 1. Generate all simulation inputs + Makefile
./01_generate.py --n 200

# 2. Run all Delphy jobs
./02_run.py                # uses half the CPUs by default
./02_run.py --jobs 20      # or specify parallelism explicitly

# 3. Analyze results
./03_analyze.py --n 200

# 4. Plot results
./04_plot.py

# If step 3 fails due to low ESS, increase the number of steps and rerun:
./01_generate.py --n 200 --steps 1000000000
make -C sims clean && ./02_run.py
./03_analyze.py --n 200
./04_plot.py
```

---

## Notes from initial run

More replicates than expected had low ESS (compared to `02_simple`).
On investigation, the affected replicates tend to be those with high
sampled mutation rates and long branches near the root carrying many
mutations.  This is a known area where Delphy converges slowly — the
tree rearrangement moves struggle when branches accumulate many
mutations.  This is not a bug in the WCSS setup; the excluded
replicates are consistent with expected MCMC mixing difficulties.