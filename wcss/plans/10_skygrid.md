# Plan: WCSS with Skygrid Population Model

## Goal

Validate Delphy's inference of a Skygrid population model with 5 knots
(log-linear interpolation), including inference of the GMRF precision
parameter tau.  The simulation generates a random Skygrid population
curve by first sampling tau from its prior, then sampling gamma_k offsets
from the GMRF with that tau, then shifting to achieve a sampled N_bar.
Delphy infers all Skygrid parameters (gamma_0..gamma_4 and tau) from the
simulated data.

This builds on `06_missing_data`, so it also includes site-rate
heterogeneity (alpha) and missing data in the simulated sequences.

This uses new features in both Delphy and Sapling:

- **Delphy:** Skygrid N_bar prior (`--v0-skygrid-nbar-prior-mean` /
  `--v0-skygrid-nbar-prior-stddev`), explicit knot dates
  (`--v0-skygrid-first-knot-date` / `--v0-skygrid-last-knot-date`).
- **Sapling:** Skygrid population model (`--skygrid-first-knot-date`,
  `--skygrid-last-knot-date`, `--skygrid-gamma`, `--skygrid-type`).

## Differences from `06_missing_data`

1. **Population model changes from exponential to Skygrid** with 5 knots,
   log-linear interpolation.
2. **`n0` and `g` are replaced by gamma_0..gamma_4** (the log-population
   values at the 5 knots) **and tau** (the GMRF precision).
3. **tau is sampled from a tight Gamma prior** per replicate, then used to
   generate the gamma_k random walk.
4. **Population curve is sampled from the GMRF prior**: first sample tau
   from Gamma(36, 24), then sample N_bar from InvGamma(mean=1 yr,
   stddev=0.2 yr), then sample gamma_k offsets using the sampled tau, and
   shift to achieve the desired N_bar.
5. **Sapling is invoked with `--skygrid-*` options** instead of
   `--exp-pop-*`.
6. **Delphy is invoked with Skygrid options**: 5 knots, explicit dates,
   log-linear, N_bar prior, tau inference enabled, low-gamma barrier
   disabled.
7. **N_bar is included as a derived parameter** via augmented log files
   (each replicate's log file gets an extra `N_bar` column before analysis).
8. **`skygrid.isloglinear` and `skygrid.cutOff` are added to `ESS_IGNORE`**
   since they are constant.

## Configuration

- **Directory:** `wcss/10_skygrid/`
- **Tips:** Same as `06_missing_data` (200 tips, dates uniform over 2025)
- **N (number of simulation replicates):** 10 (debug), then 200 (final)
- **Steps:** 1,000,000,000 (5,000,000 per tip x 200 tips);
  debug: 20,000,000 (100,000 per tip x 200 tips)
- **Population model:** Skygrid, 5 knots, log-linear
  - Knot dates: 2025-01-01 through 2026-01-01, evenly spaced
    (4 intervals of ~91.25 days each)
  - In practice, we pass `--v0-skygrid-first-knot-date 2025-01-01` and
    `--v0-skygrid-last-knot-date 2026-01-01` and Delphy computes the
    internal knot positions.  Likewise for Sapling.
  - Low-gamma barrier: **disabled** (`--v0-skygrid-disable-low-pop-barrier`)
- **Sampled parameters (drawn fresh per replicate):**
  - GMRF precision tau ~ Gamma(alpha=36, beta=24)
    (mean=1.5, stddev=0.25)
  - N_bar ~ InvGamma(mean=1 yr, stddev=0.2 yr)
    => alpha = 2 + 1^2/0.04 = 27, beta = 1 * 26 = 26
  - gamma_k offsets ~ GMRF random walk with the sampled tau
  - Mutation rate mu ~ Gamma(alpha=1, beta=1000)
    (i.e., Exponential with mean 1e-3 subst/site/year)
  - Stationary frequencies (pi_A, pi_C, pi_G, pi_T) ~ Dirichlet(1,1,1,1)
  - HKY kappa ~ LogNormal(meanlog=1.0, sdlog=1.25)
  - Site-rate heterogeneity alpha ~ Exponential(mean=1)
- **Fixed parameters:**
  - Genome size L = 30,000 sites
- **Missing data:** Same as `06_missing_data` (mean 3 gaps of length 500,
  mean 3 missing sites)

## Log file column ordering for Skygrid

Delphy's log file outputs `skygrid.logPopSize` values in **BEAST
convention**: reversed knot order and in **log-years** (not log-days).
Specifically (`beasty_output.cpp:40-41`):

```cpp
// gamma = ln(N(t)), where N(t) is measured in days in Delphy,
// in years in BEAST (hence -ln(365.0) and (M-k))
os << pop_model.gamma(pop_model.M() - k) - std::log(365.0) << "\t";
```

So for 5 knots (M=4):

| Log column          | Internal knot | Meaning                    |
|---------------------|---------------|----------------------------|
| skygrid.logPopSize1 | gamma(4)      | Most recent knot (gamma_4) |
| skygrid.logPopSize2 | gamma(3)      | gamma_3                    |
| skygrid.logPopSize3 | gamma(2)      | gamma_2                    |
| skygrid.logPopSize4 | gamma(1)      | gamma_1                    |
| skygrid.logPopSize5 | gamma(0)      | Oldest knot (gamma_0)      |

All values in the log file are in **log-years**.  Sapling also outputs
gamma_k in log-years.  **No unit conversion is needed** when comparing
true values to posterior samples.

## Augmented log files for N_bar

N_bar = exp(mean(logPopSize1..5)) is a derived quantity not present in
Delphy's raw log file.  To keep all downstream analysis uniform (using
loganalyser for HPD, ESS, etc.), we **augment each replicate's log file**
with an extra `N_bar` column before running loganalyser.

In `03_analyze.py`, before calling loganalyser:
1. For each replicate, read `delphy.log`.
2. Compute `N_bar = exp(mean(skygrid.logPopSize1..5))` for each sample
   row.  (Values are in log-years, so N_bar is in years.)
3. Write `delphy-digested.log` in the same directory, with the `N_bar`
   column appended.
4. Run loganalyser on the digested log files instead of the originals.
5. The rank computation (`compute_normalized_ranks`) also reads from the
   digested log files, so it can access the `N_bar` column.

This way, N_bar is treated identically to all other parameters throughout
the entire pipeline (loganalyser, coverage, ranks, and plots).

## Analyzed parameters

All parameters below (including N_bar) are available as columns in the
augmented log files and are processed uniformly by loganalyser:

| Display name | Log column            | True value source                     |
|--------------|-----------------------|---------------------------------------|
| mu           | meanRate              | subst_model.mu                        |
| alpha        | alpha                 | subst_model.site_rate_heterogeneity_alpha |
| tau          | skygrid.precision     | wcss_true_params.json: tau            |
| gamma_0      | skygrid.logPopSize5   | pop_model.gamma_k[0]                  |
| gamma_1      | skygrid.logPopSize4   | pop_model.gamma_k[1]                  |
| gamma_2      | skygrid.logPopSize3   | pop_model.gamma_k[2]                  |
| gamma_3      | skygrid.logPopSize2   | pop_model.gamma_k[3]                  |
| gamma_4      | skygrid.logPopSize1   | pop_model.gamma_k[4]                  |
| N_bar        | N_bar (augmented)     | exp(mean(gamma_k))                    |
| kappa        | kappa                 | subst_model.kappa                     |
| pi_A         | frequencies1          | subst_model.pi[0]                     |
| pi_C         | frequencies2          | subst_model.pi[1]                     |
| pi_G         | frequencies3          | subst_model.pi[2]                     |
| pi_T         | frequencies4          | subst_model.pi[3]                     |
| rootHeight   | rootHeight            | tree_stats.tree_height                |

Note the reversed mapping: `gamma_0` (oldest knot) maps to
`skygrid.logPopSize5`, and `gamma_4` (most recent) maps to
`skygrid.logPopSize1`.

## Storing true tau

Sapling does not know about tau (it's a GMRF prior parameter, not a
simulation parameter).  The true tau for each replicate is stored in a
supplementary file `sim_XXX/wcss_true_params.json`:

```json
{"tau": 1.234}
```

Written by `01_generate.py` after running Sapling.  Read by
`03_analyze.py` in `read_true_params()`.

## Deliverables

Files in `wcss/10_skygrid/`:

0. **`00_plan.md`** -- Symlink to `../plans/10_skygrid.md` (this plan)
1. **`01_generate.py`** -- Generate simulation inputs, Makefile, and
   mutation count diagnostics
2. **`02_run.py`** -- Run Delphy via `make -jN`
3. **`03_analyze.py`** -- Run loganalyser, check ESS, compute coverage/ranks
4. **`04_plot.py`** -- Produce all plots from TSV files

---

## Part 1: `01_generate.py` -- Generate

Copy from `06_missing_data/01_generate.py` with these changes:

### Configuration changes

- **Remove:** `N0_PRIOR_MEAN`, `N0_PRIOR_STDDEV`, `N0_PRIOR_ALPHA`,
  `N0_PRIOR_BETA`, `G_PRIOR_MEAN`
- **Add:**
  - `NBAR_PRIOR_MEAN = 1.0` (years)
  - `NBAR_PRIOR_STDDEV = 0.2` (years)
  - `NBAR_PRIOR_ALPHA = 2.0 + NBAR_PRIOR_MEAN**2 / NBAR_PRIOR_STDDEV**2`
    (= 27.0)
  - `NBAR_PRIOR_BETA = NBAR_PRIOR_MEAN * (NBAR_PRIOR_ALPHA - 1.0)`
    (= 26.0)
  - `NUM_KNOTS = 5`
  - `SKYGRID_FIRST_KNOT_DATE = "2025-01-01"`
  - `SKYGRID_LAST_KNOT_DATE = "2026-01-01"`
  - `TAU_PRIOR_ALPHA = 36.0`
  - `TAU_PRIOR_BETA = 24.0`
    (Gamma prior on tau: mean=1.5, stddev=0.25)

### Step 2 changes: Sample parameters

- `sample_prior_params(rng)` replaces n0/g sampling with Skygrid sampling:
  ```python
  import math

  # tau ~ Gamma(alpha, beta)
  tau = rng.gamma(TAU_PRIOR_ALPHA, 1.0 / TAU_PRIOR_BETA)

  # N_bar ~ InvGamma(mean=1, stddev=0.2)
  inv_nbar = rng.gamma(NBAR_PRIOR_ALPHA, 1.0 / NBAR_PRIOR_BETA)
  N_bar = 1.0 / inv_nbar  # years

  # Sample gamma_k offsets from GMRF random walk with precision tau
  gamma = [0.0] * NUM_KNOTS
  for k in range(1, NUM_KNOTS):
      gamma[k] = gamma[k-1] + rng.normal(0, 1.0 / math.sqrt(tau))

  # Shift gamma_k's so that exp(mean(gamma_k)) = N_bar (in years)
  mean_gamma = sum(gamma) / NUM_KNOTS
  target_mean = math.log(N_bar)
  gamma = [g + (target_mean - mean_gamma) for g in gamma]

  # Warn if any gamma_k is very low (approaching default barrier at ~-5.9)
  min_gamma = min(gamma)
  if min_gamma < -5.0:
      print(f"    WARNING: min gamma_k = {min_gamma:.2f} "
            f"(N = {math.exp(min_gamma):.4f} yr), "
            f"approaching low-gamma barrier territory")
  ```
- Returns `(tau, gamma, N_bar, mu, kappa, pi, alpha)` instead of
  `(n0, g, mu, kappa, pi, alpha)`.

### Sapling invocation changes

- Replace `--exp-pop-n0` / `--exp-pop-g` with:
  ```
  --skygrid-first-knot-date 2025-01-01
  --skygrid-last-knot-date 2026-01-01
  --skygrid-gamma "g0,g1,g2,g3,g4"   (comma-separated, in log-years)
  --skygrid-type log-linear
  ```
- Keep all other Sapling flags from `06_missing_data` (site-rate
  heterogeneity, missing data parameters, etc.).
- After running Sapling, write `sim_XXX/wcss_true_params.json` with
  `{"tau": ...}`.
- The print statement shows `tau=..., gamma=(...), N_bar=...`.

### Step 3 changes: Makefile

The Delphy invocation in the Makefile pattern rule changes:
- **Remove:**
  - `--v0-pop-n0-prior-mean ...`
  - `--v0-pop-n0-prior-stddev ...`
  - `--v0-pop-g-prior-exponential-with-mean ...`
- **Add:**
  - `--v0-pop-model skygrid`
  - `--v0-skygrid-num-parameters 5`
  - `--v0-skygrid-first-knot-date 2025-01-01`
  - `--v0-skygrid-last-knot-date 2026-01-01`
  - `--v0-skygrid-type log-linear`
  - `--v0-skygrid-nbar-prior-mean 1`
  - `--v0-skygrid-nbar-prior-stddev 0.2`
  - `--v0-skygrid-infer-prior-smoothness`
  - `--v0-skygrid-tau-prior-alpha 36`
  - `--v0-skygrid-tau-prior-beta 24`
  - `--v0-skygrid-disable-low-pop-barrier`
- **Keep:** `--v0-site-rate-heterogeneity` (from `06_missing_data`)

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
	  --v0-pop-model skygrid \
	  --v0-skygrid-num-parameters 5 \
	  --v0-skygrid-first-knot-date 2025-01-01 \
	  --v0-skygrid-last-knot-date 2026-01-01 \
	  --v0-skygrid-type log-linear \
	  --v0-skygrid-nbar-prior-mean 1 \
	  --v0-skygrid-nbar-prior-stddev 0.2 \
	  --v0-skygrid-infer-prior-smoothness \
	  --v0-skygrid-tau-prior-alpha 36 \
	  --v0-skygrid-tau-prior-beta 24 \
	  --v0-skygrid-disable-low-pop-barrier \
	  --v0-site-rate-heterogeneity \
	&& touch $@
```

---

## Part 2: `02_run.py` -- Run

Identical to `06_missing_data/02_run.py`.  Copy verbatim.

---

## Part 3: `03_analyze.py` -- Analyze

Copy from `06_missing_data/03_analyze.py` with these changes:

### Configuration changes

- Replace PARAMS:
  ```python
  PARAMS = [
      ("mu",         "meanRate"),
      ("alpha",      "alpha"),
      ("tau",        "skygrid.precision"),
      ("gamma_0",    "skygrid.logPopSize5"),   # oldest knot (reversed order)
      ("gamma_1",    "skygrid.logPopSize4"),
      ("gamma_2",    "skygrid.logPopSize3"),
      ("gamma_3",    "skygrid.logPopSize2"),
      ("gamma_4",    "skygrid.logPopSize1"),   # most recent knot
      ("N_bar",      "N_bar"),                 # augmented column
      ("kappa",      "kappa"),
      ("pi_A",       "frequencies1"),
      ("pi_C",       "frequencies2"),
      ("pi_G",       "frequencies3"),
      ("pi_T",       "frequencies4"),
      ("rootHeight", "rootHeight"),
  ]
  ```
- `ESS_IGNORE = {"skygrid.isloglinear", "skygrid.cutOff"}`
  (both are constant; tau is inferred, so its ESS matters)

### New step: Augment log files

Before running loganalyser, augment each replicate's log file with an
`N_bar` column:

```python
GAMMA_COLS = [
    "skygrid.logPopSize1", "skygrid.logPopSize2", "skygrid.logPopSize3",
    "skygrid.logPopSize4", "skygrid.logPopSize5",
]

def digest_log_file(src_path, dst_path):
    """Produce delphy-digested.log: add computed N_bar column."""
    with open(src_path) as fin, open(dst_path, "w") as fout:
        gamma_indices = None
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            if gamma_indices is None:
                # First non-comment line is the header
                headers = line.rstrip("\n").split("\t")
                gamma_indices = [headers.index(c) for c in GAMMA_COLS]
                fout.write(line.rstrip("\n") + "\tN_bar\n")
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            gammas = [float(fields[i]) for i in gamma_indices]
            nbar = math.exp(sum(gammas) / len(gammas))
            fout.write(line.rstrip("\n") + f"\t{nbar}\n")
```

Both loganalyser and `compute_normalized_ranks` read from
`delphy-digested.log` (not `delphy.log`), so the `N_bar` column is
available to both.

### Read true parameters

- `read_true_params()` reads gamma_k from `info["pop_model"]["gamma_k"]`
  (already in log-years, matching the log file) and tau from
  `wcss_true_params.json`:
  ```python
  gamma_k = info["pop_model"]["gamma_k"]

  # Read supplementary params (tau) not stored by Sapling
  wcss_path = os.path.join(sim_dir, "wcss_true_params.json")
  with open(wcss_path) as f:
      wcss = json.load(f)

  true = {
      "mu": info["subst_model"]["mu"],
      "alpha": info["subst_model"]["site_rate_heterogeneity_alpha"],
      "tau": wcss["tau"],
      "gamma_0": gamma_k[0],    # log-years, no conversion needed
      "gamma_1": gamma_k[1],
      "gamma_2": gamma_k[2],
      "gamma_3": gamma_k[3],
      "gamma_4": gamma_k[4],
      "N_bar": math.exp(sum(gamma_k) / len(gamma_k)),  # years
      "kappa": info["subst_model"]["kappa"],
      "pi_A": info["subst_model"]["pi"][0],
      "pi_C": info["subst_model"]["pi"][1],
      "pi_G": info["subst_model"]["pi"][2],
      "pi_T": info["subst_model"]["pi"][3],
      "rootHeight": info["tree_stats"]["tree_height"],
  }
  ```

### All other steps

No structural changes needed.  Since N_bar is a column in the augmented
log files, loganalyser handles its HPD, ESS, and mean.  The coverage,
rank, and plotting code all treat it identically to any other PARAMS
entry.

---

## Part 4: `04_plot.py` -- Plot

Copy from `06_missing_data/04_plot.py` with these changes:

### Configuration changes

- Replace PARAMS to match `03_analyze.py` (tau, gamma_0..4, N_bar instead
  of n0/g).

### All other code

No structural changes.  The summary figure has 15 parameter rows
(mu, alpha, tau, gamma_0..4, N_bar, kappa, pi_A..T, rootHeight) + 1
clade coverage row.

---

## GMRF sampling details

The GMRF prior on the Skygrid log-population curve is:

```
log pi_GMRF = sum_{k=1}^M [(1/2) log(tau) - (tau/2) (gamma_k - gamma_{k-1})^2]
```

So `gamma_k - gamma_{k-1} ~ Normal(0, 1/sqrt(tau))`.

With tau sampled from Gamma(36, 24) (mean=1.5, stddev=0.25):
- At tau=1.5: each step has stddev 0.82 log units, so population changes
  by factor ~2.3 per quarter.
- At tau=1.0: each step has stddev 1.0 log units (factor ~2.7 per
  quarter).
- At tau=2.0: each step has stddev 0.71 log units (factor ~2.0 per
  quarter).

The overall level is then set by shifting all gamma_k's so that
`exp(mean(gamma_k)) = N_bar`, where N_bar is drawn from the InvGamma
prior (mean=1 yr, stddev=0.2 yr).

### Low-population barrier analysis

The low-gamma barrier (disabled in this WCSS) sits at N = 1 day =
1/365 yr, i.e., gamma = ln(1/365) ≈ -5.90 in log-years.  Starting from
gamma_0 = 0 (N = 1 yr), the probability that the 4-step random walk
hits this barrier is negligible:

- **At fixed tau = 1.5 (the prior mean):** P(min gamma_k < -5.90) ≈
  0.016%, or ~0.03 expected hits per 200 replicates.
- **With tau drawn from Gamma(36, 24):** P(min gamma_k < -5.90) ≈
  0.032%, or ~0.06 expected hits per 200 replicates (estimated by Monte
  Carlo with 10^6 samples).

Both estimates are before the shift to achieve N_bar, which further
reduces the probability by centering the walk.  In practice, the
low-gamma barrier is essentially never reached.

---

## Directory structure at runtime

```
wcss/10_skygrid/
  00_plan.md -> ../plans/10_skygrid.md
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
    mutation_counts_histogram.pdf
    mutation_counts_ecdf.pdf
    clade_coverage.pdf
    ecdf_mu.pdf
    ecdf_alpha.pdf
    ecdf_tau.pdf
    ecdf_gamma_0.pdf
    ...
    ecdf_N_bar.pdf
    rank_histogram_mu.pdf
    rank_histogram_alpha.pdf
    rank_histogram_tau.pdf
    rank_histogram_gamma_0.pdf
    ...
    rank_histogram_N_bar.pdf
    scatter_mu.pdf
    scatter_alpha.pdf
    scatter_tau.pdf
    scatter_gamma_0.pdf
    ...
    scatter_N_bar.pdf
    wcss_summary.pdf
  sims/
    Makefile
    tips.txt
    mutation_counts.tsv
    sim_000/
      sim.maple
      sim_info.json
      wcss_true_params.json
      sim.nwk
      sim.nexus
      sim.fasta
      delphy.log
      delphy-digested.log
      delphy.trees
      delphy.dphy
      .done
    sim_001/
      ...
    ...
```

## Execution workflow

```bash
cd wcss/10_skygrid

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

- **ESS for Skygrid gamma_k and tau:** Mixing for the gamma_k values
  relies on both HMC moves and the zero-mode Gibbs move.  Tau is updated
  via a Gibbs move.  ESS should be reasonable at 1B steps, but may need
  to be increased if mixing is slow.

- **Correlation between gamma_k values:**  Adjacent gamma_k's are
  correlated by the GMRF prior.  This is expected and doesn't affect WCSS
  calibration (the WCSS validates marginal coverage independently).

- **Disabled low-gamma barrier:**  The barrier is a soft constraint that
  prevents log-population from going too low.  Disabling it means the
  posterior can explore very small population sizes, which is fine for the
  WCSS since the true curves are generated from the same prior.  This also
  avoids the barrier introducing a mismatch between the simulation prior
  and the inference prior.  A warning is printed during generation if any
  sampled gamma_k is very low (below -5.0 in log-years, i.e., N < 0.007
  years ≈ 2.5 days).  As shown in the barrier analysis above, the
  probability of reaching the barrier is negligible (~0.03% per replicate).

- **Tau prior tightness:**  The Gamma(36, 24) prior on tau (mean=1.5,
  stddev=0.25, 95% interval ≈ [1.04, 2.03]) is fairly tight.  This
  prevents extreme tau values that would produce either nearly flat curves
  (high tau) or wildly oscillating ones (low tau).
