# Well-Calibrated Simulation Study — Background

## Mendes, Bouckaert, Carvalho & Drummond (2025) — "How to Validate a Bayesian Evolutionary Model"

<https://doi.org/10.1093/sysbio/syae064>

The paper describes two complementary protocols for validating that a Bayesian phylogenetic model is correctly implemented:

### 1. Coverage Validation
- Repeat *n* times (typically n=100): draw parameters from the prior, simulate data under those parameters, then run MCMC inference on the simulated data.
- For each parameter, check whether the true (generating) value falls within the 95% HPD interval of the posterior. If the model is correct, this should happen ~95% of the time.
- **Scenario 1** (correct model): coverage ~= 95%, posterior means scatter around the identity line (true = inferred).
- **Scenario 2** (bug in likelihood): coverage collapses (e.g., 5%), posteriors don't track true values.
- **Scenario 3** (model misspecification): coverage can look fine for some parameters but fail for others.

### 2. Rank Uniformity Validation (RUV)
- For each simulation replicate, compute the rank of the true parameter value among the posterior samples.
- If the model is correct, these ranks should be **uniformly distributed** on (1, L+1).
- Visualized as a histogram (should be flat) and an ECDF plot (should follow the diagonal).
- RUV is more sensitive than coverage alone — it can detect subtle biases even when coverage looks OK.
- They extend RUV to **tree space** using Robinson-Foulds distance to a reference tree, enabling validation of tree topology inference (Figure 8).

### Practical Guidelines
- n=100 simulations is a good starting point; more is better but costlier.
- Use ESS as a ceiling for the number of posterior samples *L* to thin to.
- Data sets should be large enough that posteriors are distinguishable from priors (otherwise bugs can hide).
- Histogram bins should be ~L/20 wide so each bin has >=1 expected count.

---

## Bouckaert — DeveloperManual (GitHub)

<https://github.com/rbouckaert/DeveloperManual>

This is the practical companion, showing how to implement these studies in BEAST 2:

### Validation Hierarchy
1. **Unit tests** — direct correctness of likelihoods, operators
2. **Prior sampling** — verify MCMC sampling from the prior matches a `DirectSimulator`
3. **Parameter recovery (fixed tree)** — simulate data on a known tree, infer parameters back
4. **Full well-calibrated simulation study** — simulate tree + parameters from prior, simulate data, run full MCMC inference, repeat n times
5. **Robustness testing** — what happens under model violations

### Workflow
1. Use `DirectSimulator` to draw (tree, parameters) from the prior and simulate a sequence alignment.
2. Run MCMC on each simulated alignment (same model, same priors).
3. Collect posteriors, compute coverage and rank statistics.
4. Tools: `TraceKSStats` for KS tests, `NumericalScoreFunctionStatistics` for score-function validation, Tracer for visual inspection.

---

## The bottom line for us

A well-calibrated simulation study for Delphy would involve:
1. **Simulating** many (e.g., 100) datasets by drawing trees and parameters from Delphy's prior and then generating sequence data.
2. **Running Delphy** on each simulated dataset.
3. **Checking** that the true generating values fall within the 95% HPD ~95% of the time, and that parameter ranks are uniformly distributed.

This would validate that Delphy's MCMC machinery is correctly sampling from the posterior.

---

## Clade Coverage Validation

- <https://github.com/rbouckaert/DeveloperManual> (section "Clade coverage")
- <https://github.com/christiaanjs/beast-validation/blob/master/src/beastvalidation/experimenter/CladeCoverageCalculator.java>

Coverage and RUV validate scalar parameters, but a phylogenetic MCMC also
infers tree topology.  Clade coverage validation checks whether the
posterior clade probabilities are well-calibrated — that is, whether a
clade assigned posterior probability *p* is truly present in the generating
tree about *p* fraction of the time.

### Key insight

In each simulation replicate, we draw model parameters and a tree from the
prior, then simulate data.  This is equivalent to first drawing data from
P(Data) and then drawing a single model (including tree) from the posterior
P(Model | Data).  The true tree is therefore *one sample from the
posterior*.

This means that for any clade with posterior frequency *p*, asking "is this
clade present in the true tree?" is a Bernoulli trial with success
probability *p*.  If the model is well-calibrated, clades with posterior
frequency ~75% should appear in the true tree ~75% of the time.

To test this, we group all such Bernoulli trials by their success
probability into bins (e.g., 70–80%) and check whether the observed
success rate matches the bin's expected rate.  Since a single replicate
provides few trials per bin, we pool trials across many replicates to
gather enough statistics.

### Procedure

Run *N* simulation replicates.  Each replicate has:
- A **true tree** (from the simulation), which defines a set of true clades.
- A **posterior sample of trees** (from MCMC inference), which defines a
  set of observed clades, each with a posterior frequency (the fraction of
  posterior trees containing that clade).

For each replicate independently:

1. Collect all clades that appear in any posterior tree and compute each
   clade's **posterior frequency** (count / number of posterior trees).
   Bin each clade by its posterior frequency into one of *K* equal-width
   bins (e.g., 0–10%, 10–20%, ..., 90–100%).  Increment a global
   **totals** counter for the appropriate bin.

2. For each clade in the **true tree**, look it up in that replicate's
   posterior clades.  If found, bin it by its posterior frequency and
   increment a global **true hits** counter for that bin.  If not found
   in the posterior at all, count it in the 0% bin.

The totals and true-hits counters accumulate across all *N* replicates.
For each bin *b*, the bar height is:

> (true hits in bin *b*) / (total clades in bin *b*) x 100%

This is the observed success rate of the Bernoulli trials in that bin.

### Interpretation

Plot these fractions as a bar chart with a diagonal *x = y* reference
line:

- **Well-calibrated:** bars track the diagonal.  Clades with 70–80%
  posterior frequency are truly present ~75% of the time.
- **Over-confident:** bars fall below the diagonal.  The model assigns
  higher clade probabilities than warranted.
- **Under-confident:** bars rise above the diagonal.  The model is too
  conservative in its clade support.
