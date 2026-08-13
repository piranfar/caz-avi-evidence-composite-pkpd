# hujam-caz-avi — software package

**What R6 delivers.** Not a rewrite — the analysis code in `code/` is unchanged and remains
regression-tested against the frozen v16 outputs. R6 packages it: an installable module, a formal
interface contract separating the reusable decision layer from the drug-pair-specific engine, a
conformance checker, and a second test suite covering the decision layer's own identities.

**License:** MIT (code, this file) / CC BY 4.0 (data summaries and documentation), matching the
parent repository. See `LICENSE`.

---

## What this is

Two things, cleanly separable:

1. **A verified simulation engine for one drug pair** (`code/model2_engine.py`) — continuous-infusion
   ceftazidime/avibactam, reproducing the manuscript's frozen primary-model outputs to
   floating-point precision. This part is **not generic** and should not be genericised: its
   parameters, dose fractions and regimen table belong to this drug pair, and a different pair needs
   its own, regardless of how the code is organised.

2. **A decision-analytic layer that is generic** (`code/model2_hujam.py` and friends) — parameter
   uncertainty, prespecified scenario sets, value of information (EVPI/EVPPI), misselection and
   regret against a fixed-threshold decision rule, and a triage rule for partial monitoring. This
   layer only ever reads two keys from the engine's output — `joint_pta` and `exceedance` — and is
   therefore portable to any drug pair whose engine exposes the same two-function contract.

That contract is written down and enforced in `code/interface.py`, not just implied by convention.

---

## Install

```bash
pip install -e .
```

Installs under the import name `hujam` (declared in `pyproject.toml`), backed by the `code/`
directory — nothing is duplicated or renamed on disk, so every path referenced in
`MODEL1_REPORT.md`, `MODEL2_REPORT.md` and the audit logs still resolves.

```bash
pip install -e ".[figures]"   # adds matplotlib, needed only for model1_finalise.py's plots
```

**This package depends on the primary model from the parent repository** — `code/model2_engine.py`
and `code/test_model1.py` import `reproduce_primary_run.py` directly (read-only) to reuse its
verified pharmacokinetic constants rather than duplicate them. It is not meant to be extracted and
used standalone without that dependency satisfied.

That model exists in **two different locations depending on which repository you cloned**, and this
package finds either automatically, preferring the first:

1. `revision_support/reproduce_primary_run.py` — the local development layout of the full IJAA
   submission package. Its frozen reference table (`revision_support/outputs/primary_pta_results.csv`)
   is matched to floating-point precision (tolerance 0.0) — verified repeatedly throughout this
   project.
2. `src/cazavi/reproduce_primary_run.py` — this package's actual layout in
   `caz-avi-evidence-composite-pkpd` on GitHub, confirmed byte-identical to (1). That repository ships
   two copies of the primary PTA table. `data/reference/primary_pta_results.csv` was found to differ
   from a fresh run of the repo's own bundled code by up to 0.355 percentage points — a pre-existing
   staleness in that one file, not in this package, and not the copy that repo's own `README.md`
   points readers to for its figures and tables. `data/processed/primary_pta_results.csv` — the
   results-facing copy the parent repo's README actually cites — matches a fresh run of the same code
   to floating-point exactness on every regimen and metric. This package's regression test therefore
   reads `data/processed/`, not `data/reference/`, and uses the same exact-zero tolerance in both
   layouts; there is no fudge factor. (The stale `data/reference/` copy is a separate, low-severity
   cleanup item on that repository's own `main` branch, outside the scope of this package.)

Both layouts are tested: `test_model1.py` (119 checks) and `test_model2.py` (44 checks) pass
identically in each, including a real end-to-end check of `pip install -e .` run from inside a clone
of `caz-avi-evidence-composite-pkpd` itself, imported from a completely unrelated directory
afterward.

Tested with Python 3.14, NumPy 2.5, SciPy 1.18 on Windows. Requires Python ≥ 3.10.

## Quickstart

```python
import hujam.model2_engine as E
import hujam.interface as I

# 1. The regression test: with zero uncertainty this must reproduce v16 exactly
ok, worst = E.verify_against_frozen()

# 2. Confirm the engine satisfies the contract the decision layer depends on
I.check_engine(E)

# 3. Run the decision layer
import hujam.model2_hujam as H
result = H.run(n_draws=500, n_per_class=2000,
                rho_scenario="C1_cojutti", target_scenario="T4_uniform")
optimal, summary = H.optimality_and_regret(result)
```

Or run any module as a script directly, exactly as it was developed and tested:

```bash
cd code
python model2_engine.py            # regression test against the frozen v16 outputs
python interface.py                # conformance check, printed
python test_model1.py              # 119 checks: dose conversion, PTA/CFR, Model 1 structure
python test_model2.py              # 44 checks: decision-layer identities, samplers, interface
python model2_hujam.py             # the full decision-analytic pipeline (21 scenario combinations)
```

## Adapting this to a different drug pair

1. Write your own engine module exposing `draw_population(n_per_class, seed)` and
   `evaluate(population, params, regimens)`, returning `{regimen: {"joint_pta": array, "exceedance":
   float, ...}}` for every regimen — see the docstrings in `code/interface.py` for the full contract.
2. Run `interface.check_engine(your_module)`. It checks determinism, required keys, value ranges and
   consistent array lengths across regimens — the things that silently break the decision layer if
   they are wrong.
3. Reuse `model2_hujam.utility`, `optimality_and_regret`, `evppi`, and `convergence_check` unchanged.
   They only ever read `joint_pta` and `exceedance`.
4. Write your own version of `limiting_probability` if your drug pair has more than two components
   or a differently structured target — the current one is specific to a two-component product and
   reads two extra keys (`caz_pta`, `avi_pta`) that only this engine provides.
5. `model2_monitoring.py` and `model2_triage.py` implement *this pair's* specific
   predict-one-component-from-the-other logic and are not claimed as portable; they are a worked
   example of what to build for your own pair's assay-substitution question, not a drop-in.

## What is tested, and what that does and does not establish

| Suite | Checks | What it establishes |
|---|---|---|
| `test_model1.py` | 119 | Dose conversion, infusion rate, free/total concentration, clearance transformation, correlated random effects, renal classes, MIC weighting, PTA, CFR, limiting component, exposure constraint, ELF transformation, Model 1's structural model and variance structure |
| `test_model2.py` | 44 | Engine/decision-toolkit conformance, the zero-uncertainty regression against frozen v16 outputs, all seven target-distribution samplers against their documented support and moments (six analyst-specified, one evidence-derived from Coleman 2014), both correlation scenario samplers, the utility function's monotonicity, the regret identity (zero for the best-in-hindsight choice, non-negative everywhere), and EVPI ≥ EVPPI ≥ 0 |

**What passing tests do not establish:** correctness of the underlying pharmacology (that rests on
independent parameter verification against the source publications, documented in
`RESULT_PROVENANCE_MATRIX.csv`), or that Model 1's clearance-correlation estimate generalises beyond
the CRRT population it was fitted in (that limitation is stated throughout `MODEL1_REPORT.md` and is
not something a test suite can resolve).

## A correctness fix this packaging effort found

Writing `test_model2.py`'s sampler checks caught a real bug: the C1 (published-value) correlation
scenario sampled `Normal(0.94, 0.238*0.94)` on the correlation scale directly and clipped to a valid
range, which clipped **39.5% of draws** to the 0.999 boundary and pulled the sampled mean down to
0.877 — a full 0.06 below the value the scenario is meant to represent. Fixed by sampling on the
Fisher-z scale instead (the same approach the C2 scenario already used), which removes the clipping
entirely and recovers a median of 0.940. Every Model 2 result that used this scenario was rerun; see
`MODEL2_REPORT.md` for the corrected numbers and the magnitude of the change.

## Reproducibility

- Every stochastic step uses an explicit, recorded seed (`MASTER_SEED = 20260811` in
  `model2_hujam.py`; per-analysis seeds documented in each module).
- `code/model2_engine.py` is regression-tested against `revision_support/outputs/primary_pta_results.csv`
  at import time by anyone who runs it — not merely by a test that is easy to skip.
- Outputs are frozen, machine-readable CSVs under `outputs/`, each traceable to the script and log
  that produced it (`audit/log_*.txt`).

## Citation

If you use this package, cite the parent repository (`CITATION.cff` at its root) and this directory
specifically as the source of the joint clearance-correlation model and the decision-analytic layer.
