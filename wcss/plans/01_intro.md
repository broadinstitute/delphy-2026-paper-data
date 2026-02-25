# Well-Calibrated Simulation Study — Background

## Mendes, Bouckaert, Carvalho & Drummond (2025) — "How to Validate a Bayesian Evolutionary Model"

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
