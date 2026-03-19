# Plan: Final WCSS — High-Mutation Stress Test

## Goal

Validate that Delphy's inference remains well-calibrated in a
high-mutation regime.  This study uses 500 tips with a wide prior on
the mutation rate mu, so many replicates will have far more mutations
than the ~2 × 500 = 1000 expected for a densely sampled genomic epi
dataset.  Branches routinely carry many mutations, stress-testing
Delphy's tree rearrangement moves.

The model is deliberately simple — constant population, no site-rate
heterogeneity, no missing data, no tip-date uncertainty — so that any
calibration failures can be attributed to the high-mutation regime
rather than model complexity.

This study evolved from preparatory study `03_free_mu`, scaled up to
500 tips with per-replicate tip dates and longer chains.

## Model

### Tree model

- **Tips:** 500, with dates drawn uniformly over 2025-01-01 to
  2025-12-31 (inclusive), generated independently per replicate.
- **Population model:** Constant, with fixed effective population
  size n0 = 1.0 year.
- **Reference date (t0):** 2026-01-01.

### Substitution model

- **Genome:** L = 30,000 sites.
- **Model:** HKY (Hasegawa–Kishino–Yano).
- **Mutation rate:** mu ~ Gamma(shape=1, rate=1000 year⁻¹), i.e.,
  Exponential with mean 1e-3 subst/site/year.  This is a wide prior:
  the coefficient of variation is 1.0, so mu values ranging from
  ~1e-4 to ~3e-3 subst/site/year are common.
- **Transition/transversion ratio:** kappa ~ LogNormal(meanlog=1.0,
  sdlog=1.25).  This is Delphy's hardcoded prior.
- **Stationary frequencies:** (pi_A, pi_C, pi_G, pi_T) ~
  Dirichlet(1, 1, 1, 1).  This is Delphy's hardcoded prior.

### Features not included

- No site-rate heterogeneity (all sites have rate 1.0).
- No missing data.
- No tip-date uncertainty (all tip dates are known exactly).

## Simulation parameters

Per replicate, the following are drawn from their priors and passed
to Sapling:

| Parameter | Prior | Sapling flag |
|-----------|-------|-------------|
| mu | Gamma(shape=1, rate=1000 year⁻¹) [subst/site/year] | `--mu` |
| kappa | LogNormal(meanlog=1.0, sdlog=1.25) | `--hky-kappa` |
| (pi_A, pi_C, pi_G, pi_T) | Dirichlet(1,1,1,1) | `--hky-pi-{A,C,G,T}` |

Fixed Sapling parameters:

| Parameter | Value | Sapling flag |
|-----------|-------|-------------|
| n0 | 1.0 year | `--const-pop-n0 1.0` |
| L | 30,000 sites | `--num-sites 30000` |
| t0 | 2026-01-01 | `--t0 2026-01-01` |

## Delphy invocation

```
delphy \
  --v0-in-maple sim_NNN/sim.maple \
  --v0-steps 2000000000 \
  --v0-out-log-file sim_NNN/delphy.log \
  --v0-log-every 2000000 \
  --v0-out-trees-file sim_NNN/delphy.trees \
  --v0-tree-every 2000000 \
  --v0-out-delphy-file sim_NNN/delphy.dphy \
  --v0-delphy-snapshot-every 2000000 \
  --v0-mu-prior-alpha 1 \
  --v0-mu-prior-beta 1000 \
  --v0-init-final-pop-size 1.0 \
  --v0-fix-final-pop-size \
  --v0-init-pop-growth-rate 0.0 \
  --v0-fix-pop-growth-rate
```

The `log_every` and `tree_every` intervals are `steps // 1000`,
giving 1000 posterior samples per run.

## Inferred parameters

The WCSS validates Delphy's posterior for these parameters:

| Display name | Delphy log column | True value source |
|-------------|-------------------|-------------------|
| mu | meanRate | `sim_info.json → subst_model.mu` |
| kappa | kappa | `sim_info.json → subst_model.kappa` |
| pi_A | frequencies1 | `sim_info.json → subst_model.pi[0]` |
| pi_C | frequencies2 | `sim_info.json → subst_model.pi[1]` |
| pi_G | frequencies3 | `sim_info.json → subst_model.pi[2]` |
| pi_T | frequencies4 | `sim_info.json → subst_model.pi[3]` |
| rootHeight | rootHeight | `sim_info.json → tree_stats.tree_height` |

ESS is checked for all inferred parameters (no parameters are
excluded from the ESS check).

## Configuration summary

| Setting | Value |
|---------|-------|
| NUM_TIPS | 500 |
| NUM_SITES | 30,000 |
| DEFAULT_N | 200 |
| DEFAULT_STEPS | 2,000,000,000 |
| CONST_POP_N0 | 1.0 year |
| MU_PRIOR_ALPHA | 1.0 |
| MU_PRIOR_BETA | 1000.0 year⁻¹ |
| KAPPA_MEAN_LOG | 1.0 (hardcoded in Delphy) |
| KAPPA_SIGMA_LOG | 1.25 (hardcoded in Delphy) |
| TIP_DATE_START | 2025-01-01 |
| TIP_DATE_END | 2025-12-31 |
| DEFAULT_MASTER_SEED | 2025 |

## Deliverables

Files in `wcss/21_final_high_mutation/`:

0. **`00_plan.md`** — Symlink to `../plans/21_final_high_mutation.md`
1. **`01_generate.py`** — Generate simulation inputs, Makefile, and
   mutation count diagnostics
2. **`02_run.py`** — Run Delphy via `make -jN`
3. **`03_analyze.py`** — Run loganalyser, check ESS, compute
   coverage/ranks
4. **`04_plot.py`** — Produce all plots from TSV files

---

## Part 1: `01_generate.py` — Generate

This script performs three steps:
1. For each replicate: generate tip dates, draw model parameters from
   their priors, and run Sapling to simulate data.
2. Generate a Makefile to drive all Delphy runs.
3. Collect mutation count statistics and produce diagnostic plots.

### Tip date generation

Each replicate gets its own tip dates.  There is no shared
`tips.txt` — each replicate writes to `sim_NNN/tips.txt`.

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
    """Sample mu, kappa and pi from their priors."""
    # mu ~ Gamma(shape=1, rate=1000 year^-1) = Exponential(mean=1e-3 subst/site/year)
    mu = rng.gamma(MU_PRIOR_ALPHA, 1.0 / MU_PRIOR_BETA)
    # (pi_A, pi_C, pi_G, pi_T) ~ Dirichlet(1, 1, 1, 1)
    pi = rng.dirichlet([1.0, 1.0, 1.0, 1.0])
    # kappa ~ LogNormal(meanlog=1.0, sdlog=1.25)
    log_kappa = rng.normal(KAPPA_MEAN_LOG, KAPPA_SIGMA_LOG)
    kappa = np.exp(log_kappa)
    return mu, kappa, pi
```

### Sapling invocation

```python
cmd = [
    sapling_path,
    "--tip-file", tips_file,
    "--t0", T0,
    "--const-pop-n0", str(CONST_POP_N0),
    "--mu", str(mu),
    "--hky-kappa", str(kappa),
    "--hky-pi-A", str(pi[0]),
    "--hky-pi-C", str(pi[1]),
    "--hky-pi-G", str(pi[2]),
    "--hky-pi-T", str(pi[3]),
    "--num-sites", str(NUM_SITES),
    "--seed", str(replicate_seed),
    "--out-maple", os.path.join(sim_dir, "sim.maple"),
    "--out-info", os.path.join(sim_dir, "sim_info.json"),
    "--out-newick", os.path.join(sim_dir, "sim.nwk"),
    "--out-nexus", os.path.join(sim_dir, "sim.nexus"),
    "--out-fasta", os.path.join(sim_dir, "sim.fasta"),
]
```

### Main loop

```python
for i in range(n):
    rng = np.random.default_rng(master_seed + i)
    sim_dir = os.path.join(sims_dir, f"sim_{i:03d}")
    os.makedirs(sim_dir, exist_ok=True)
    tips_file = os.path.join(sim_dir, "tips.txt")
    generate_tips(rng, tips_file)
    mu, kappa, pi = sample_prior_params(rng)
    replicate_seed = int(rng.integers(0, 2**31))
    run_sapling(i, mu, kappa, pi, replicate_seed, tips_file, script_dir)
```

### Makefile

```makefile
DELPHY = ../../../delphy

SIM_DIRS := $(wildcard sim_[0-9]*)

all: $(addsuffix /.done,$(SIM_DIRS))

sim_%/.done: sim_%/sim.maple
	$(DELPHY) \
	  --v0-in-maple $< \
	  --v0-steps 2000000000 \
	  --v0-out-log-file sim_$*/delphy.log \
	  --v0-log-every 2000000 \
	  --v0-out-trees-file sim_$*/delphy.trees \
	  --v0-tree-every 2000000 \
	  --v0-out-delphy-file sim_$*/delphy.dphy \
	  --v0-delphy-snapshot-every 2000000 \
	  --v0-mu-prior-alpha 1 \
	  --v0-mu-prior-beta 1000 \
	  --v0-init-final-pop-size 1.0 \
	  --v0-fix-final-pop-size \
	  --v0-init-pop-growth-rate 0.0 \
	  --v0-fix-pop-growth-rate \
	&& touch $@

clean:
	rm -f sim_*/delphy.* sim_*/delphy-digested.log sim_*/.done

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
--steps STEPS      MCMC steps per replicate (default: 2,000,000,000)
--master-seed S    Master seed for parameter sampling (default: 2025)
```

---

## Part 2: `02_run.py` — Run

Runs `make -jN` in the `sims/` directory, where N defaults to half
the available CPUs.  This script is identical across all WCSS studies.

### CLI arguments

```
--jobs N     Number of parallel jobs (default: half of CPU count)
```

---

## Part 3: `03_analyze.py` — Analyze

This script performs three steps:
1. **Digest log files:** Create `delphy-digested.log` symlinks (no
   column transformations needed for this study).
2. **Run loganalyser:** Extract posterior summaries (mean, ESS, 95%
   HPD) for each replicate, with 30% burn-in.
3. **Check ESS:** Flag replicates with low ESS (< 200) and exclude
   those with very low ESS (< 150) on any parameter.
4. **Compute coverage and ranks:** For each parameter, check whether
   the true value falls within the 95% HPD interval and compute its
   normalized rank within the posterior samples.

### Configuration

```python
PARAMS = [
    ("mu",         "meanRate"),
    ("kappa",      "kappa"),
    ("pi_A",       "frequencies1"),
    ("pi_C",       "frequencies2"),
    ("pi_G",       "frequencies3"),
    ("pi_T",       "frequencies4"),
    ("rootHeight", "rootHeight"),
]

ESS_IGNORE = set()
ESS_THRESHOLD_LOW = 200
ESS_THRESHOLD_VERY_LOW = 150
```

### True parameter extraction

```python
def read_true_params(info):
    return {
        "mu":         info["subst_model"]["mu"],
        "kappa":      info["subst_model"]["kappa"],
        "pi_A":       info["subst_model"]["pi"][0],
        "pi_C":       info["subst_model"]["pi"][1],
        "pi_G":       info["subst_model"]["pi"][2],
        "pi_T":       info["subst_model"]["pi"][3],
        "rootHeight": info["tree_stats"]["tree_height"],
    }
```

### Clade coverage

In addition to continuous-parameter coverage, the script also
computes clade coverage: for each replicate, it reads the true tree
and checks what fraction of true clades appear in the posterior tree
samples (from `delphy.trees`), binned by posterior support threshold.

### Output files

All written to `analyses/`:

- `true_params.tsv` — true parameter values per replicate
- `loganalyser_output.tsv` — raw loganalyser output (all replicates)
- `ess_check.tsv` — ESS values per parameter per replicate
- `excluded_replicates.tsv` — replicates excluded due to very low ESS
- `loganalyser_output_filtered.tsv` — loganalyser output after
  excluding bad replicates
- `coverage_summary.txt` — per-parameter 95% HPD coverage rates
- `ranks.tsv` — normalized ranks per parameter per replicate
- `clade_coverage.tsv` — clade coverage by support threshold
- `clade_coverage_raw.tsv` — per-replicate clade coverage data

### CLI arguments

```
--n N                          Number of replicates (default: 200)
--burnin PCT                   Burn-in percentage (default: 30)
--ignore-low-ess               Don't abort on low ESS warnings
--force-include-all-replicates Skip all ESS-based exclusion
```

---

## Part 4: `04_plot.py` — Plot

Reads the TSV files produced by `03_analyze.py` and generates all
plots.  Uses the same PARAMS list as `03_analyze.py`.

### Plots produced (in `plots/`)

For each parameter (mu, kappa, pi_A–pi_T, rootHeight):
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

---

## Directory structure at runtime

```
wcss/21_final_high_mutation/
  00_plan.md -> ../plans/21_final_high_mutation.md
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
    ecdf_kappa.pdf
    ecdf_pi_A.pdf
    ecdf_pi_C.pdf
    ecdf_pi_G.pdf
    ecdf_pi_T.pdf
    ecdf_rootHeight.pdf
    rank_histogram_mu.pdf
    rank_histogram_kappa.pdf
    rank_histogram_pi_A.pdf
    rank_histogram_pi_C.pdf
    rank_histogram_pi_G.pdf
    rank_histogram_pi_T.pdf
    rank_histogram_rootHeight.pdf
    scatter_mu.pdf
    scatter_kappa.pdf
    scatter_pi_A.pdf
    scatter_pi_C.pdf
    scatter_pi_G.pdf
    scatter_pi_T.pdf
    scatter_rootHeight.pdf
    wcss_summary.pdf
  sims/
    Makefile
    mutation_counts.tsv
    sim_000/
      tips.txt
      sim.maple
      sim_info.json
      sim.nwk
      sim.nexus
      sim.fasta
      delphy.log
      delphy-digested.log
      delphy.trees
      delphy.dphy
      .done
    sim_001/
      tips.txt
      ...
    ...
```

Note: there is no shared `sims/tips.txt`.  Each `sim_NNN/tips.txt`
contains the per-replicate tip dates.

---

## Execution workflow

```bash
cd wcss/21_final_high_mutation

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

- **ESS in high-mutation replicates:** With a wide Exponential prior
  on mu, some replicates will have very high mutation rates and trees
  with thousands of mutations.  Delphy's SPR moves struggle when
  branches carry many mutations (the candidate region weights become
  very peaked).  Low ESS in these replicates is expected and was
  observed in preparatory study 03.  The 4× increase in chain length
  (2B vs 500M) should help, but some replicates may still need to be
  excluded.

- **Per-replicate tip dates:** Each replicate draws its own tip
  dates, so the model parameters for a given replicate index will
  differ from study 03.  This is intentional — the replicates are
  not meant to be comparable across studies.
