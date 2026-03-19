# Plan: Harmonize WCSS study scripts

## Motivation

Studies 02–11 were written sequentially over several weeks.  As the
framework evolved, later studies gained features that weren't backported
to earlier ones.  The result: if you diff two studies' scripts, you see
a mix of genuine conceptual differences (different parameters, priors,
analysis steps) and accidental structural drift (missing CLI flags,
inconsistent file naming, hardcoded headers).

This plan eliminates the accidental drift so that a diff between any
two studies' scripts reveals only the intentional differences.

## Changes

### 1. Uniform log digestion (`delphy-digested.log`)

**Current state:**
- Studies 02–09: loganalyser and `compute_normalized_ranks` read
  `delphy.log` directly.
- Study 10: augments `delphy.log` into `delphy_augmented.log` (adds
  computed `N_bar` column), then loganalyser and
  `compute_normalized_ranks` read `delphy_augmented.log`.
- Study 11: strips `age(...)` columns from `delphy.log` into
  `delphy-digested.log`, then loganalyser reads the digested files.
  `compute_normalized_ranks` still reads `delphy.log` (it needs the
  `age(...)` columns for tip-date analysis later).

**Problem:** `run_loganalyser` and `compute_normalized_ranks` refer to
different filenames across studies, obscuring what is otherwise
identical code.

**Fix:** All studies get a "Step 0: Digest log files" that produces
`delphy-digested.log` in each sim directory.  `run_loganalyser` and
`compute_normalized_ranks` always read `delphy-digested.log`.

- **Studies 02–09:** Digestion is a no-op — create a symbolic link
  `delphy-digested.log -> delphy.log` to avoid duplicating data.
- **Study 10:** Digestion = augmentation (add `N_bar` column).
  Output is a regular file `delphy-digested.log` (renamed from
  `delphy_augmented.log`).
- **Study 11:** Digestion = stripping `age(...)` columns.  Output is
  a regular file `delphy-digested.log` (already named this way).

The digest function signature is the same everywhere:

```python
def digest_log_file(src_path, dst_path):
    """Produce delphy-digested.log from delphy.log."""
```

Studies 02–09:
```python
def digest_log_file(src_path, dst_path):
    """Produce delphy-digested.log from delphy.log (symlink, no-op)."""
    if os.path.islink(dst_path) or os.path.exists(dst_path):
        os.remove(dst_path)
    os.symlink(os.path.basename(src_path), dst_path)
```

Study 10:
```python
def digest_log_file(src_path, dst_path):
    """Produce delphy-digested.log: add computed N_bar column."""
    # (existing augment_log_file body, renamed)
```

Study 11:
```python
def digest_log_file(src_path, dst_path):
    """Produce delphy-digested.log: strip age(...) columns."""
    # (existing _strip_age_columns body, renamed)
```

The call site in `main()` is the same everywhere:

```python
# Step 0: Digest log files
print("Digesting log files...")
for i in range(n):
    sim_dir = os.path.join(script_dir, "sims", f"sim_{i:03d}")
    digest_log_file(
        os.path.join(sim_dir, "delphy.log"),
        os.path.join(sim_dir, "delphy-digested.log"))
print(f"  Digested {n} log files")
```

`run_loganalyser` everywhere changes from reading `delphy.log` (or
`delphy_augmented.log` or `delphy-digested.log`) to always reading
`delphy-digested.log`:

```python
log_files = [
    os.path.join("sims", f"sim_{i:03d}", "delphy-digested.log")
    for i in range(n)
]
```

`compute_normalized_ranks` likewise reads `delphy-digested.log`
everywhere (including study 10, where it currently reads
`delphy_augmented.log`).

Makefile `clean` targets everywhere add `delphy-digested.log`:

```makefile
clean:
	rm -f sim_*/delphy.* sim_*/delphy-digested.log sim_*/.done
```

(Study 11 already has this.  Studies 02–10 need the addition.
For symlinks, `rm -f` removes the link, not the target.)

### 2. Add `--force-include-all-replicates` to studies 02–03

**Current state:** Studies 02 and 03's `check_ess()` has signature:

```python
def check_ess(la_df, analyses_dir, ignore_low_ess=False):
```

Studies 04–11 have:

```python
def check_ess(la_df, analyses_dir, ignore_low_ess=False,
              force_include_all=False):
```

with the corresponding `--force-include-all-replicates` CLI flag and
the `if force_include_all:` block inside `check_ess()`.

**Fix:** Copy the `force_include_all` parameter, CLI flag, and internal
logic from study 04 into studies 02 and 03.  Also update the usage
comment at the top of the file.

### 3. Dynamic `true_params.tsv` header in study 02

**Current state:** Study 02 hardcodes the TSV header:

```python
f.write("replicate\tkappa\tpi_A\tpi_C\tpi_G\tpi_T\trootHeight\n")
```

and writes each row with explicit field references:

```python
f.write(f"sim_{i:03d}\t{tv['kappa']}\t{tv['pi_A']}\t{tv['pi_C']}\t"
        f"{tv['pi_G']}\t{tv['pi_T']}\t{tv['rootHeight']}\n")
```

Studies 03–11 use the dynamic form driven by `param_names`:

```python
f.write("replicate\t" + "\t".join(param_names) + "\n")
# ...
vals = "\t".join(str(tv[name]) for name in param_names)
f.write(f"sim_{i:03d}\t{vals}\n")
```

**Fix:** Replace the hardcoded form in study 02 with the dynamic form.

### 4. Backport mutation count analysis to studies 02–07 and 10

**Current state:** Studies 08, 09, and 11 have a "Step 4: Mutation
count statistics" in `01_generate.py` that:

1. Reads `num_mutations` from each replicate's `sim_info.json`.
2. Prints summary statistics (mean, std, min, max, quartiles).
3. Saves a TSV (`sims/mutation_counts.tsv`).
4. Generates a histogram and eCDF plot (`plots/mutation_counts_*.pdf`).

This is useful for sanity-checking simulated data regardless of which
parameters are free, so it should be present in all studies.

**Fix:** Add the mutation count analysis block to `01_generate.py` in
studies 02–07 and 10.  The code is identical across studies — it only
reads `sim_info.json` fields that all studies produce.

## Summary of changes by study

| Study | Digest (1) | force-include (2) | Dynamic header (3) | Mut counts (4) |
|-------|-----------|-------------------|-------------------|---------------|
| 02    | Add symlink | Add | Fix | Add |
| 03    | Add symlink | Add | Already dynamic | Add |
| 04    | Add symlink | Already present | Already dynamic | Add |
| 05    | Add symlink | Already present | Already dynamic | Add |
| 06    | Add symlink | Already present | Already dynamic | Add |
| 07    | Add symlink | Already present | Already dynamic | Add |
| 08    | Add symlink | Already present | Already dynamic | Already present |
| 09    | Add symlink | Already present | Already dynamic | Already present |
| 10    | Rename augmented -> digested | Already present | Already dynamic | Add |
| 11    | Rename _strip -> digest | Already present | Already dynamic | Already present |

## What NOT to change

- **`02_run.py`** — already identical across all studies.
- **`04_plot.py`** — differences are genuinely study-specific (study 10
  has extra skygrid params, study 11 has tip-date scatter plots).
- **`PARAMS`, `ESS_IGNORE`, `read_true_params()`** — these vary per
  study by design.
- **Sapling invocation, Makefile Delphy flags** — study-specific.

## Execution order

Changes 1–3 are independent and can be done in any order.  Change 4
(mutation counts) requires re-running `01_generate.py` in studies that
didn't previously have it, which regenerates the Makefile and sim
inputs.  Since the data will be regenerated anyway for the final paper
runs, this is fine.

Suggested commit grouping:
1. One commit for the digest harmonization (all studies).
2. One commit for `--force-include-all-replicates` + dynamic header
   (studies 02–03 only).
3. One commit for mutation count backport (studies 02–07, 10).
