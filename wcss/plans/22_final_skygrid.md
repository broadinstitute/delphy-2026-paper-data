# Plan: Final WCSS — Skygrid with All Key Features

## Goal

Validate that Delphy's inference is well-calibrated when using a
Skygrid population model together with all key features: site-rate
heterogeneity, missing data, and tip-date uncertainty.  This study
uses 500 tips with tight priors tuned so that the total number of
mutations is concentrated around 2 x 500 = 1000, giving approximately
1 mutation per branch on average.  This is the regime where Delphy's
reformulation of Bayesian phylogenetics is most effective.

This study evolved from preparatory studies `10_skygrid` (Skygrid
population model) and `11_tip_date_uncertainty` (tip-date uncertainty
with exponential population), scaled up to 500 tips with per-replicate
tip dates and longer chains.  Study 23 is identical to this study
except for the population model (exponential instead of Skygrid).

## Model

### Tree model

- **Tips:** 500, with dates drawn uniformly over 2025-01-01 to
  2025-12-31 (inclusive), generated independently per replicate.
- **Population model:** Skygrid with 5 knots spanning 2025-01-01 to
  2026-01-01, log-linear interpolation.  The population size at each
  knot k is exp(gamma_k) years, where the gamma_k follow a Gaussian
  Markov Random Field (GMRF) with precision tau.  The mean population
  size N_bar = exp(mean(gamma_k)) is drawn from an Inverse-Gamma
  prior.
- **Reference date (t0):** 2026-01-01.

### Substitution model

- **Genome:** L = 30,000 sites.
- **Model:** HKY (Hasegawa-Kishino-Yano).
- **Mutation rate:** mu ~ Gamma(mean=1e-3 subst/site/year,
  stddev=1e-4 subst/site/year).  Equivalently,
  Gamma(shape=100, rate=100,000 year^{-1}).  This is a tight prior:
  the coefficient of variation is 0.1, concentrating mu near 1e-3.
- **Transition/transversion ratio:** kappa ~ LogNormal(meanlog=1.0,
  sdlog=1.25).  This is Delphy's hardcoded prior.
- **Stationary frequencies:** (pi_A, pi_C, pi_G, pi_T) ~
  Dirichlet(1, 1, 1, 1).  This is Delphy's hardcoded prior.

### Site-rate heterogeneity

- **Model:** Continuous Gamma.  Each site l has a rate multiplier
  nu_l ~ Gamma(shape=alpha, rate=alpha), so that E[nu_l] = 1.
- **Shape parameter:** alpha ~ Exponential(mean=1.0).

### Missing data

Sapling generates missing data per tip with these parameters:

| Parameter | Value | Sapling flag |
|-----------|-------|-------------|
| Mean number of gaps | 3.0 | `--missing-data-mean-num-gaps` |
| Mean gap length | 500.0 sites | `--missing-data-mean-gap-length` |
| Mean number of isolated missing sites | 3.0 | `--missing-data-mean-num-missing-sites` |

This produces approximately 5% missing data per tip on average.

### Tip-date uncertainty

A fraction of tips have uncertain dates:

| Uncertainty level | Fraction | Sapling flag |
|-------------------|----------|-------------|
| Month-uncertain (date truncated to YYYY-MM) | 15% | `--p-tip-date-uncertain-upto-month 0.15` |
| Year-uncertain (date truncated to YYYY) | 5% | `--p-tip-date-uncertain-upto-year 0.05` |

The remaining 80% of tips have exact dates.

## Simulation parameters

Per replicate, the following are drawn from their priors and passed
to Sapling:

| Parameter | Prior | Sapling flag |
|-----------|-------|-------------|
| N_bar | InvGamma(mean=0.75 year, stddev=0.15 year) [alpha=27, beta=19.5] | (derived; see gamma_k) |
| tau | Gamma(shape=36, rate=24) [mean=1.5, stddev=0.25] | (derived; see gamma_k) |
| gamma_k (k=0..4) | GMRF random walk with precision tau, shifted so exp(mean(gamma_k)) = N_bar | `--skygrid-gamma` |
| mu | Gamma(mean=1e-3, stddev=1e-4) [shape=100, rate=100,000 year^{-1}] subst/site/year | `--mu` |
| kappa | LogNormal(meanlog=1.0, sdlog=1.25) | `--hky-kappa` |
| (pi_A, pi_C, pi_G, pi_T) | Dirichlet(1,1,1,1) | `--hky-pi-{A,C,G,T}` |
| alpha | Exponential(mean=1.0) | `--site-rate-heterogeneity-alpha` |

### GMRF prior implications

The GMRF prior on the Skygrid log-population curve is:

    log pi_GMRF = sum_{k=1}^{M} [(1/2) log(tau) - (tau/2) (gamma_k - gamma_{k-1})^2]

So gamma_k - gamma_{k-1} ~ Normal(0, 1/sqrt(tau)).  With 5 knots
spanning one year (4 intervals of ~91 days each), each step
represents approximately one quarter.

With tau sampled from Gamma(36, 24) (mean=1.5, stddev=0.25):
- At tau=1.5: each step has stddev 0.82 log units, so population
  changes by a factor of ~2.3 per quarter.
- At tau=1.0: each step has stddev 1.0 log units (factor ~2.7 per
  quarter).
- At tau=2.0: each step has stddev 0.71 log units (factor ~2.0 per
  quarter).

The overall level is set by shifting all gamma_k so that
exp(mean(gamma_k)) = N_bar, where N_bar is drawn from the InvGamma
prior (mean=0.75 year, stddev=0.15 year).

### Low-population barrier analysis

The low-gamma barrier (disabled in this WCSS) sits at N = 1 day =
1/365 year, i.e., gamma = ln(1/365) ~ -5.90 in log-years.  With
N_bar = 0.75 year (the prior mean), the walk centers around
ln(0.75) ~ -0.29, leaving 5.6 log-units of margin to the barrier.
Monte Carlo estimation (10^7 samples; see Appendix A) with tau
drawn from Gamma(36, 24) and N_bar drawn from InvGamma(27, 19.5)
finds zero barrier hits, confirming that the barrier is essentially
never reached.  The existing tau prior does not need tightening.

Fixed Sapling parameters:

| Parameter | Value | Sapling flag |
|-----------|-------|-------------|
| L | 30,000 sites | `--num-sites 30000` |
| t0 | 2026-01-01 | `--t0 2026-01-01` |
| Skygrid first knot | 2025-01-01 | `--skygrid-first-knot-date 2025-01-01` |
| Skygrid last knot | 2026-01-01 | `--skygrid-last-knot-date 2026-01-01` |
| Skygrid type | log-linear | `--skygrid-type log-linear` |

### Gamma_k sampling procedure

1. Sample tau ~ Gamma(shape=36, rate=24).
2. Sample N_bar ~ InvGamma(alpha=27, beta=19.5) by drawing
   1/N_bar ~ Gamma(alpha=27, rate=19.5).
3. Build gamma_k via GMRF random walk: gamma_0 = 0, then
   gamma_k = gamma_{k-1} + Normal(0, 1/sqrt(tau)) for k=1..4.
4. Shift all gamma_k so that exp(mean(gamma_k)) = N_bar.
5. Pass the gamma_k values as a comma-separated string to
   `--skygrid-gamma`.

### Supplementary true parameters

Sapling does not store tau in `sim_info.json`.  The script writes
`wcss_true_params.json` with `{"tau": ...}` alongside each
replicate's simulation outputs.

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
  --v0-skygrid-nbar-prior-mean 0.75 \
  --v0-skygrid-nbar-prior-stddev 0.15 \
  --v0-skygrid-infer-prior-smoothness \
  --v0-skygrid-tau-prior-alpha 36 \
  --v0-skygrid-tau-prior-beta 24 \
  --v0-skygrid-disable-low-pop-barrier \
  --v0-site-rate-heterogeneity
```

The `log_every` and `tree_every` intervals are `steps // 1000`,
giving 1000 posterior samples per run.

Note: tip-date uncertainty requires no extra Delphy flags — Delphy
detects uncertain dates from the MAPLE input automatically.

## Inferred parameters

The WCSS validates Delphy's posterior for these parameters:

| Display name | Delphy log column | True value source |
|-------------|-------------------|-------------------|
| mu | meanRate | `sim_info.json -> subst_model.mu` |
| alpha | alpha | `sim_info.json -> subst_model.site_rate_heterogeneity_alpha` |
| tau | skygrid.precision | `wcss_true_params.json -> tau` |
| gamma_0 | skygrid.logPopSize5 | `sim_info.json -> pop_model.gamma_k[0]` |
| gamma_1 | skygrid.logPopSize4 | `sim_info.json -> pop_model.gamma_k[1]` |
| gamma_2 | skygrid.logPopSize3 | `sim_info.json -> pop_model.gamma_k[2]` |
| gamma_3 | skygrid.logPopSize2 | `sim_info.json -> pop_model.gamma_k[3]` |
| gamma_4 | skygrid.logPopSize1 | `sim_info.json -> pop_model.gamma_k[4]` |
| N_bar | N_bar | computed: exp(mean(gamma_k)) |
| kappa | kappa | `sim_info.json -> subst_model.kappa` |
| pi_A | frequencies1 | `sim_info.json -> subst_model.pi[0]` |
| pi_C | frequencies2 | `sim_info.json -> subst_model.pi[1]` |
| pi_G | frequencies3 | `sim_info.json -> subst_model.pi[2]` |
| pi_T | frequencies4 | `sim_info.json -> subst_model.pi[3]` |
| rootHeight | rootHeight | `sim_info.json -> tree_stats.tree_height` |

Note: Delphy's Skygrid log columns are numbered in reverse order:
`skygrid.logPopSize1` is the most recent knot (gamma_4) and
`skygrid.logPopSize5` is the oldest (gamma_0).

N_bar is a derived quantity not directly in the log file.  It is
computed by the log digestion step as exp(mean(gamma_k)).

ESS is checked for all inferred parameters.  The following
observables are excluded from the ESS check because they are
constant: `skygrid.isloglinear`, `skygrid.cutOff`.

### Tip-date posteriors

In addition to the continuous parameters above, the WCSS also
validates tip-date inference.  For each tip with an uncertain date,
the script computes coverage and normalized ranks of the true tip
date within its posterior.  Results are reported separately for
month-uncertain and year-uncertain tips.

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
| NBAR_PRIOR_MEAN | 0.75 year |
| NBAR_PRIOR_STDDEV | 0.15 year |
| NBAR_PRIOR_ALPHA | 27.0 |
| NBAR_PRIOR_BETA | 19.5 year |
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

## Deliverables

Files in `wcss/22_final_skygrid/`:

0. **`00_plan.md`** -- Symlink to `../plans/22_final_skygrid.md`
1. **`01_generate.py`** -- Generate simulation inputs, Makefile, and
   mutation count diagnostics
2. **`02_run.py`** -- Run Delphy via `make -jN`
3. **`03_analyze.py`** -- Run loganalyser, check ESS, compute
   coverage/ranks, analyze tip-date posteriors
4. **`04_plot.py`** -- Produce all plots from TSV files

---

## Part 1: `01_generate.py` -- Generate

This script performs three steps:
1. For each replicate: generate tip dates, draw model parameters from
   their priors, and run Sapling to simulate data.
2. Generate a Makefile to drive all Delphy runs.
3. Collect mutation count statistics and produce diagnostic plots.

### Tip date generation

Each replicate gets its own tip dates.  There is no shared
`tips.txt` -- each replicate writes to `sim_NNN/tips.txt`.

```python
def generate_tips(rng, out_path):
    """Write tips.txt with NUM_TIPS tips, dates uniform over 2025."""
    num_days = (TIP_DATE_END - TIP_DATE_START).days  # 364
    day_offsets = rng.integers(0, num_days + 1, size=NUM_TIPS)
    with open(out_path, "w") as f:
        for i in range(NUM_TIPS):
            tip_date = TIP_DATE_START + timedelta(days=int(day_offsets[i]))
            f.write(f"TIP_{i:03d}|{tip_date.isoformat()}\n")
```

### Prior sampling

```python
def sample_prior_params(rng):
    """Sample tau, gamma, mu, kappa, pi, and alpha from their priors."""
    # tau ~ Gamma(shape=36, rate=24)
    tau = rng.gamma(TAU_PRIOR_ALPHA, 1.0 / TAU_PRIOR_BETA)

    # N_bar ~ InvGamma(alpha=27, beta=19.5)
    inv_nbar = rng.gamma(NBAR_PRIOR_ALPHA, 1.0 / NBAR_PRIOR_BETA)
    N_bar = 1.0 / inv_nbar

    # gamma_k via GMRF random walk
    gamma = [0.0] * NUM_KNOTS
    for k in range(1, NUM_KNOTS):
        gamma[k] = gamma[k-1] + rng.normal(0, 1.0 / math.sqrt(tau))
    mean_gamma = sum(gamma) / NUM_KNOTS
    target_mean = math.log(N_bar)
    gamma = [g + (target_mean - mean_gamma) for g in gamma]

    # mu ~ Gamma(mean=MU_PRIOR_MEAN, stddev=MU_PRIOR_STDDEV)
    mu_alpha = (MU_PRIOR_MEAN / MU_PRIOR_STDDEV) ** 2
    mu_beta = MU_PRIOR_MEAN / MU_PRIOR_STDDEV ** 2
    mu = rng.gamma(mu_alpha, 1.0 / mu_beta)

    pi = rng.dirichlet([1.0, 1.0, 1.0, 1.0])
    log_kappa = rng.normal(KAPPA_MEAN_LOG, KAPPA_SIGMA_LOG)
    kappa = np.exp(log_kappa)

    # alpha ~ Exponential(mean=ALPHA_PRIOR_MEAN)
    alpha = rng.exponential(ALPHA_PRIOR_MEAN)

    return tau, gamma, N_bar, mu, kappa, pi, alpha
```

### Sapling invocation

```python
cmd = [
    sapling_path,
    "--tip-file", tips_file,
    "--t0", T0,
    "--skygrid-first-knot-date", SKYGRID_FIRST_KNOT_DATE,
    "--skygrid-last-knot-date", SKYGRID_LAST_KNOT_DATE,
    "--skygrid-gamma", gamma_str,    # comma-separated gamma_k values
    "--skygrid-type", "log-linear",
    "--mu", str(mu),
    "--hky-kappa", str(kappa),
    "--hky-pi-A", str(pi[0]),
    "--hky-pi-C", str(pi[1]),
    "--hky-pi-G", str(pi[2]),
    "--hky-pi-T", str(pi[3]),
    "--site-rate-heterogeneity-alpha", str(alpha),
    "--missing-data-mean-num-gaps", str(MISSING_DATA_MEAN_NUM_GAPS),
    "--missing-data-mean-gap-length", str(MISSING_DATA_MEAN_GAP_LENGTH),
    "--missing-data-mean-num-missing-sites", str(MISSING_DATA_MEAN_NUM_MISSING_SITES),
    "--p-tip-date-uncertain-upto-month", str(P_TIP_DATE_UNCERTAIN_UPTO_MONTH),
    "--p-tip-date-uncertain-upto-year", str(P_TIP_DATE_UNCERTAIN_UPTO_YEAR),
    "--num-sites", str(NUM_SITES),
    "--seed", str(replicate_seed),
    "--out-maple", ...,
    "--out-info", ...,
    "--out-newick", ...,
    "--out-nexus", ...,
    "--out-fasta", ...,
]
```

### Supplementary true parameters

After each Sapling run, the script writes `wcss_true_params.json`
containing `{"tau": tau}` to the replicate directory.

### Main loop

```python
for i in range(n):
    rng = np.random.default_rng(master_seed + i)
    sim_dir = os.path.join(sims_dir, f"sim_{i:03d}")
    os.makedirs(sim_dir, exist_ok=True)
    tips_file = os.path.join(sim_dir, "tips.txt")
    generate_tips(rng, tips_file)
    tau, gamma, N_bar, mu, kappa, pi, alpha = sample_prior_params(rng)
    replicate_seed = int(rng.integers(0, 2**31))
    run_sapling(i, tau, gamma, N_bar, mu, kappa, pi, alpha,
                replicate_seed, tips_file, script_dir)
```

### Makefile

```makefile
DELPHY = ../../../delphy

SIM_DIRS := $(wildcard sim_[0-9]*)

all: $(addsuffix /.done,$(SIM_DIRS))

sim_%/.done: sim_%/sim.maple
	$(DELPHY) \
	  --v0-in-maple $< \
	  --v0-steps 3000000000 \
	  --v0-out-log-file sim_$*/delphy.log \
	  --v0-log-every 3000000 \
	  --v0-out-trees-file sim_$*/delphy.trees \
	  --v0-tree-every 3000000 \
	  --v0-out-delphy-file sim_$*/delphy.dphy \
	  --v0-delphy-snapshot-every 3000000 \
	  --v0-mu-prior-mean 0.001 \
	  --v0-mu-prior-stddev 0.0001 \
	  --v0-pop-model skygrid \
	  --v0-skygrid-num-parameters 5 \
	  --v0-skygrid-first-knot-date 2025-01-01 \
	  --v0-skygrid-last-knot-date 2026-01-01 \
	  --v0-skygrid-type log-linear \
	  --v0-skygrid-nbar-prior-mean 0.75 \
	  --v0-skygrid-nbar-prior-stddev 0.15 \
	  --v0-skygrid-infer-prior-smoothness \
	  --v0-skygrid-tau-prior-alpha 36 \
	  --v0-skygrid-tau-prior-beta 24 \
	  --v0-skygrid-disable-low-pop-barrier \
	  --v0-site-rate-heterogeneity \
	&& touch $@

clean:
	rm -f sim_*/delphy.* sim_*/delphy-digested.log sim_*/delphy-tips.log sim_*/.done

.PHONY: all clean
```

### Mutation count diagnostics

After all Sapling runs, the script reads `num_mutations` from each
replicate's `sim_info.json`, saves a TSV summary to
`sims/mutation_counts.tsv`, prints summary statistics, and produces
two diagnostic plots: a histogram and an eCDF of mutation counts
across replicates (saved to `plots/`).

### CLI arguments

```
--n N              Number of replicates (default: 200)
--steps STEPS      MCMC steps per replicate (default: 3,000,000,000)
--master-seed S    Master seed for parameter sampling (default: 2025)
```

---

## Part 2: `02_run.py` -- Run

Runs `make -jN` in the `sims/` directory, where N defaults to half
the available CPUs.  This script is identical across all WCSS studies.

### CLI arguments

```
--jobs N     Number of parallel jobs (default: half of CPU count)
```

---

## Part 3: `03_analyze.py` -- Analyze

This script performs six steps:
1. **Digest log files:** Create `delphy-digested.log` by stripping
   per-replicate `age(...)` columns (needed because each replicate
   has different uncertain tips) and adding a computed `N_bar` column.
2. **Produce tip-date log files:** Create `delphy-tips.log` with
   derived `tipYear(...)` columns for uncertain tips, and run
   per-replicate loganalyser on them.
3. **Run loganalyser:** Extract posterior summaries (mean, ESS, 95%
   HPD) for each replicate, with 30% burn-in.
4. **Check ESS:** Flag replicates with low ESS (< 200) and exclude
   those with very low ESS (< 150) on any parameter.  Minimum
   tip-date ESS (per month/year category) is included in the check.
5. **Compute coverage and ranks:** For each parameter, check whether
   the true value falls within the 95% HPD interval and compute its
   normalized rank within the posterior samples.
6. **Analyze tip-date posteriors:** For each uncertain tip, compute
   coverage and normalized ranks of the true tip date.

### Configuration

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

ESS_IGNORE = {"skygrid.isloglinear", "skygrid.cutOff"}
ESS_THRESHOLD_LOW = 200
ESS_THRESHOLD_VERY_LOW = 150
```

### Log digestion

The `digest_log_file` function combines two transformations:
1. Strip `age(...)` columns (per-replicate tip-date columns that
   differ across replicates, preventing loganalyser from running
   on all files at once).
2. Add a computed `N_bar` column: `N_bar = exp(mean(gamma_k))`,
   where gamma_k are read from the `skygrid.logPopSize*` columns.

### Tip-date log production

The `produce_tips_log` function extracts `age(TIP_XXX|...)` columns
from each replicate's `delphy.log` and computes
`tipYear = rootHeight + age(root) - age(TIP)`, which gives the
calendar year of the tip.  These are written to `delphy-tips.log`.

### True parameter extraction

```python
def read_true_params(sim_dir):
    info_path = os.path.join(sim_dir, "sim_info.json")
    with open(info_path) as f:
        info = json.load(f)

    gamma_k = info["pop_model"]["gamma_k"]

    wcss_path = os.path.join(sim_dir, "wcss_true_params.json")
    with open(wcss_path) as f:
        wcss = json.load(f)

    return {
        "mu": info["subst_model"]["mu"],
        "alpha": info["subst_model"]["site_rate_heterogeneity_alpha"],
        "tau": wcss["tau"],
        "gamma_0": gamma_k[0],
        "gamma_1": gamma_k[1],
        "gamma_2": gamma_k[2],
        "gamma_3": gamma_k[3],
        "gamma_4": gamma_k[4],
        "N_bar": math.exp(sum(gamma_k) / len(gamma_k)),
        "kappa": info["subst_model"]["kappa"],
        "pi_A": info["subst_model"]["pi"][0],
        "pi_C": info["subst_model"]["pi"][1],
        "pi_G": info["subst_model"]["pi"][2],
        "pi_T": info["subst_model"]["pi"][3],
        "rootHeight": info["tree_stats"]["tree_height"],
    }
```

### Clade coverage

In addition to continuous-parameter coverage, the script computes
clade coverage: for each replicate, it reads the true tree and checks
what fraction of true clades appear in the posterior tree samples
(from `delphy.trees`), binned by posterior support threshold.

### Output files

All written to `analyses/`:

- `true_params.tsv` -- true parameter values per replicate
- `loganalyser_output.tsv` -- raw loganalyser output (all replicates)
- `ess_check.tsv` -- ESS values per parameter per replicate
- `excluded_replicates.tsv` -- replicates excluded due to very low ESS
- `loganalyser_output_filtered.tsv` -- loganalyser output after
  excluding bad replicates
- `coverage_summary.txt` -- per-parameter 95% HPD coverage rates
- `ranks.tsv` -- normalized ranks per parameter per replicate
- `clade_coverage.tsv` -- clade coverage by support threshold
- `clade_coverage_raw.tsv` -- per-replicate clade coverage data
- `tip_date_ranks_month.tsv` -- normalized ranks for month-uncertain tips
- `tip_date_ranks_year.tsv` -- normalized ranks for year-uncertain tips
- `tip_date_coverage_summary.tsv` -- tip-date coverage by uncertainty type

### CLI arguments

```
--n N                          Number of replicates (default: 200)
--burnin PCT                   Burn-in percentage (default: 30)
--ignore-low-ess               Don't abort on low ESS warnings
--force-include-all-replicates Skip all ESS-based exclusion
```

---

## Part 4: `04_plot.py` -- Plot

Reads the TSV files produced by `03_analyze.py` and generates all
plots.  Uses the same PARAMS list as `03_analyze.py`.

### Plots produced (in `plots/`)

For each parameter (mu, alpha, tau, gamma_0-gamma_4, N_bar, kappa,
pi_A-pi_T, rootHeight):
- **eCDF plot** (`ecdf_{name}.pdf`): empirical CDF of the normalized
  rank, compared to the Uniform(0,1) reference line.
- **Rank histogram** (`rank_histogram_{name}.pdf`): histogram of
  normalized ranks with a 95% binomial confidence band.
- **Scatter plot** (`scatter_{name}.pdf`): true value vs posterior
  mean, with 95% HPD error bars.

Combined plots:
- **`wcss_summary.pdf`**: multi-panel figure with eCDF, rank
  histogram, and scatter for every parameter.
- **`clade_coverage.pdf`**: clade coverage vs posterior support
  threshold.

Tip-date plots:
- **`tip_date_ecdf_month.pdf`**: eCDF of normalized ranks for
  month-uncertain tips.
- **`tip_date_ecdf_year.pdf`**: eCDF of normalized ranks for
  year-uncertain tips.
- **`tip_date_rank_histogram_month.pdf`**: rank histogram for
  month-uncertain tips.
- **`tip_date_rank_histogram_year.pdf`**: rank histogram for
  year-uncertain tips.

---

## Directory structure at runtime

```
wcss/22_final_skygrid/
  00_plan.md -> ../plans/22_final_skygrid.md
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
    tip_date_ranks_month.tsv
    tip_date_ranks_year.tsv
    tip_date_coverage_summary.tsv
  plots/
    mutation_counts_histogram.pdf
    mutation_counts_ecdf.pdf
    clade_coverage.pdf
    ecdf_mu.pdf
    ecdf_alpha.pdf
    ecdf_tau.pdf
    ecdf_gamma_0.pdf
    ...
    ecdf_rootHeight.pdf
    rank_histogram_mu.pdf
    ...
    scatter_mu.pdf
    ...
    wcss_summary.pdf
    tip_date_ecdf_month.pdf
    tip_date_ecdf_year.pdf
    tip_date_rank_histogram_month.pdf
    tip_date_rank_histogram_year.pdf
  sims/
    Makefile
    mutation_counts.tsv
    sim_000/
      tips.txt
      sim.maple
      sim-COMPLETE.maple
      sim_info.json
      sim.nwk
      sim.nexus
      sim.fasta
      wcss_true_params.json
      delphy.log
      delphy-digested.log
      delphy-tips.log
      delphy.trees
      delphy.dphy
      .done
    sim_001/
      tips.txt
      ...
    ...
```

Note: there is no shared `sims/tips.txt`.  Each `sim_NNN/tips.txt`
contains the per-replicate tip dates.  `sim-COMPLETE.maple` is
produced by Sapling and contains the original exact dates for all
tips (including those with uncertain dates), used for tip-date
posterior validation.

---

## Execution workflow

```bash
cd wcss/22_final_skygrid

# --- Final run ---
./01_generate.py
./02_run.py
./03_analyze.py
./04_plot.py

# --- Debug run ---
./01_generate.py --n 10 --steps 50000000
./02_run.py
./03_analyze.py --n 10 --ignore-low-ess --force-include-all-replicates
./04_plot.py
```

---

## Potential concerns

- **Mutation count tuning:** The N_bar prior mean was chosen by
  pilot runs (200 replicates each) to place the mutation count
  distribution near the target of 2 x 500 = 1000:

  | N_bar prior mean | Mutation count mean | Mutation count median |
  |------------------|--------------------:|----------------------:|
  | 1.0 year         |               1,285 |                 1,211 |
  | 0.5 year         |                 832 |                   793 |
  | 0.75 year (chosen) |             1,061 |                 1,006 |

- **Skygrid low-population barrier:** The
  `--v0-skygrid-disable-low-pop-barrier` flag disables Delphy's
  default barrier that prevents gamma_k from going below ~-5.9.  This
  is necessary because the WCSS draws gamma_k from the prior, and the
  barrier would create a mismatch between the simulation prior and the
  inference prior.

- **ESS with Skygrid:** The Skygrid model has more parameters than
  the constant or exponential population models (5 gamma_k values
  plus tau), which may lead to slower mixing.  The 3B step chain
  length (3x study 10's 1B) should compensate.

- **Per-replicate tip dates:** Each replicate draws its own tip dates,
  so the model parameters for a given replicate index will differ from
  studies 10 and 11.  This is intentional.

- **Log digestion complexity:** Unlike study 10 (which only adds
  N_bar) or study 11 (which only strips age columns), this study's
  digest_log_file must do both: strip age columns and add N_bar.

---

## Appendix A: Low-population barrier Monte Carlo estimation

The following script estimates the probability that the GMRF random
walk hits the low-gamma barrier at gamma = ln(1/365) ~ -5.90, given
the priors on tau and N_bar.  With 10^7 samples, it finds zero hits,
confirming that the barrier is never reached in practice.

```python
import numpy as np
import math

rng = np.random.default_rng(42)
N_MC = 10_000_000
NUM_KNOTS = 5
BARRIER = math.log(1.0 / 365.0)  # ~ -5.90

TAU_ALPHA, TAU_BETA = 36.0, 24.0
NBAR_ALPHA, NBAR_BETA = 27.0, 19.5

hits = 0
for _ in range(N_MC):
    tau = rng.gamma(TAU_ALPHA, 1.0 / TAU_BETA)
    inv_nbar = rng.gamma(NBAR_ALPHA, 1.0 / NBAR_BETA)
    N_bar = 1.0 / inv_nbar

    gamma = [0.0] * NUM_KNOTS
    for k in range(1, NUM_KNOTS):
        gamma[k] = gamma[k-1] + rng.normal(0, 1.0 / math.sqrt(tau))
    mean_gamma = sum(gamma) / NUM_KNOTS
    target_mean = math.log(N_bar)
    gamma = [g + (target_mean - mean_gamma) for g in gamma]

    if min(gamma) < BARRIER:
        hits += 1

prob = hits / N_MC
print(f"P(min gamma_k < {BARRIER:.2f}) = {prob:.6f} ({prob*100:.4f}%)")
print(f"Expected hits per 200 replicates: {prob * 200:.2f}")
```
