"""The engine/decision-toolkit contract — what R6 actually generalizes.

WHY THIS FILE EXISTS
    model2_engine.py is a CAZ-AVI-SPECIFIC implementation: its Params fields
    (cl0_caz, exp_avi, ...), its dose fractions, its MIC grid and its regimen table
    all belong to this one drug pair, and none of that should be genericised —
    doing so would touch verified, regression-tested code for no real benefit,
    since a different drug pair needs its own parameters regardless of naming.

    What IS generic, and what this file documents and checks, is the DECISION
    LAYER built on top of it. Reading model2_hujam.py's `run()` and `utility()`
    line by line shows they only ever read two keys from the engine's per-regimen
    output:

        result["joint_pta"]     array of joint attainment, aligned to some grid
        result["exceedance"]    scalar percentage exceeding a safety screen

    Every function that operates ONLY on those two keys is portable to a new
    drug pair unchanged: `utility`, `optimality_and_regret`'s misselection/regret/
    EVPI computation, `evppi`, `convergence_check`. A new drug pair needs to
    write its own engine module exposing `draw_population` and `evaluate` with
    the signatures below, and can then reuse the rest of model2_hujam.py without
    modification.

    Two functions are CAZ-AVI-SPECIFIC EXTENSIONS, not part of the portable core:
        - `limiting_probability` (model2_hujam.py) additionally reads
          result["avi_pta"] and result["caz_pta"] to decide which component
          limits attainment. A new drug pair with more than two components, or a
          single combined target, would need its own version of this function.
        - model2_monitoring.py and model2_triage.py implement the specific
          "predict component B from component A" logic for this drug pair's
          assay-substitution question. That logic is inherently pair-specific
          and is not claimed as portable.

USAGE
    from interface import check_engine
    import model2_engine as E
    check_engine(E)   # raises AssertionError with a specific message on failure
"""
from __future__ import annotations

from typing import Any, Protocol, runtime_checkable

import numpy as np

REQUIRED_RESULT_KEYS = ("joint_pta", "exceedance")
EXTENSION_RESULT_KEYS = ("caz_pta", "avi_pta", "ekfc_class", "daily_g")


@runtime_checkable
class JointAttainmentEngine(Protocol):
    """The two functions a drug-pair-specific engine module must provide.

    This is documentation enforced by `check_engine`, not an import-time
    constraint — model2_engine.py is not required to subclass anything, it only
    has to expose module-level functions matching these signatures.
    """

    def draw_population(self, n_per_class: int, seed: int) -> Any:
        """Return a population object. Opaque to the decision layer: it is
        drawn once per outer-loop run and passed unchanged to `evaluate`, so its
        internal structure is entirely up to the engine. Must be DETERMINISTIC
        given (n_per_class, seed) — the decision layer relies on common random
        numbers across parameter draws to isolate the parameter's effect from
        simulation noise."""
        ...

    def evaluate(self, population: Any, params: Any,
                 regimens: list[str] | None = None) -> dict[str, dict[str, Any]]:
        """Return {regimen_label: result_dict} for one parameter draw.

        Each result_dict MUST contain:
            joint_pta   array-like, joint attainment (%) over some fixed grid
                        (e.g. MIC values), the SAME grid and length for every
                        regimen and every call, so results can be combined with
                        a fixed weights vector.
            exceedance  float, percentage of the population breaching a safety
                        or exposure screen. Use 0.0 if the drug pair has none.

        May contain anything else the drug-pair-specific extensions need (see
        module docstring)."""
        ...


def check_engine(engine_module, n_per_class=200, seed=1, verbose=True):
    """Runtime conformance check for a candidate engine module.

    Not a full test suite — it checks the CONTRACT (required keys, shapes,
    determinism), not correctness of the underlying pharmacology, which only the
    engine's own domain-specific tests can do.
    """
    checks = []

    def check(name, condition):
        checks.append((name, bool(condition)))
        if verbose:
            print(f"  {'PASS' if condition else 'FAIL':4}  {name}")

    check("exposes draw_population", hasattr(engine_module, "draw_population"))
    check("exposes evaluate", hasattr(engine_module, "evaluate"))
    check("exposes a default/BASE parameter set", hasattr(engine_module, "BASE"))

    pop1 = engine_module.draw_population(n_per_class, seed)
    pop2 = engine_module.draw_population(n_per_class, seed)
    check("draw_population is deterministic given the same seed",
          _population_equal(pop1, pop2))

    res = engine_module.evaluate(pop1, engine_module.BASE)
    check("evaluate returns a non-empty mapping of regimen -> result", len(res) > 0)

    regimen0 = next(iter(res))
    d = res[regimen0]
    for key in REQUIRED_RESULT_KEYS:
        check(f"result['{regimen0}'] has the required key '{key}'", key in d)

    if "joint_pta" in d:
        arr = np.asarray(d["joint_pta"])
        check("joint_pta is array-like and 1-dimensional", arr.ndim == 1)
        check("joint_pta values lie in [0, 100]",
              bool(np.all((arr >= 0) & (arr <= 100))))
        lengths = {len(np.asarray(res[r]["joint_pta"])) for r in res}
        check("joint_pta has the same length for every regimen", len(lengths) == 1)

    if "exceedance" in d:
        check("exceedance is a scalar in [0, 100]",
              0.0 <= float(d["exceedance"]) <= 100.0)

    n_fail = sum(1 for _, ok in checks if not ok)
    if n_fail:
        raise AssertionError(f"{n_fail} conformance check(s) failed: "
                             f"{[n for n, ok in checks if not ok]}")
    if verbose:
        print(f"\n  {len(checks)}/{len(checks)} conformance checks passed — "
              f"{engine_module.__name__} satisfies the JointAttainmentEngine contract")
    return True


def _population_equal(a, b):
    """Structural equality for the population objects this project's engines
    return (dict of class -> (array, array)); a new engine returning a
    different structure should extend this rather than have check_engine fail
    on a false negative."""
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            return False
        return all(_population_equal(a[k], b[k]) for k in a)
    if isinstance(a, tuple) and isinstance(b, tuple):
        return len(a) == len(b) and all(_population_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        return bool(np.array_equal(a, b))
    return a == b


if __name__ == "__main__":
    import model2_engine as E
    print("Checking model2_engine against the JointAttainmentEngine contract\n")
    check_engine(E)
