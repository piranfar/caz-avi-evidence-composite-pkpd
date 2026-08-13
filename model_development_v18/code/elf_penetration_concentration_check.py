"""Is the ELF penetration ratio the project applies valid at the concentrations it applies it to?

The submission package applies FIXED ELF/plasma penetration ratios of 0.52 (ceftazidime)
and 0.42 (avibactam), cited to Dimelow 2018 (Drugs R D 18:221-30, doi 10.1007/s40268-018-0241-0).
Those numbers are in the source and are quoted correctly.

But Dimelow's ratios are NOT constants. They are the value of a non-linear plasma-ELF
relationship evaluated at ONE plasma concentration each:

  ceftazidime  52%  at a plasma concentration of 15.3 mg/L
  avibactam    42%  at a plasma concentration of  2.4 mg/L

Ceftazidime's plasma-ELF link is a saturable Michaelis-Menten function; avibactam's is a
power function with exponent < 1. In both cases the penetration ratio FALLS as plasma
concentration rises. The project's continuous-infusion regimens run at steady-state
concentrations far above 15.3 and 2.4 mg/L -- the ceftazidime neurotoxicity screen alone
sits at 104 mg/L -- so the fixed ratios are being read off the curve at the wrong place.

This script rebuilds both published functions from the parameters reported in the paper,
CHECKS them against every numeric checkpoint the paper states, and only then reports the
penetration ratio at the concentrations this project actually simulates.

Nothing outside model_development_v18/ is read or written. This is a diagnostic, not a fix:
the ELF scenario code lives in revision_support/ and is deliberately left untouched.
"""

from __future__ import annotations

import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "data_external", "Dimelow2018_ELF_concentration_dependence")

# --- Ceftazidime: saturable (Michaelis-Menten) plasma-ELF link -----------------------------
# "The maximum possible ELF concentration was predicted to be 45.4 mg/l; half this ELF
#  concentration (22.7 mg/l) would be realised for a plasma concentration of 71.7 mg/l."
CAZ_EMAX = 45.4   # mg/L, maximum attainable ELF concentration
CAZ_EC50 = 71.7   # mg/L, plasma concentration giving half of EMAX in ELF


def caz_elf(cp: float) -> float:
    """ELF ceftazidime concentration (mg/L) for a plasma concentration cp (mg/L)."""
    return CAZ_EMAX * cp / (CAZ_EC50 + cp)


def caz_ratio(cp: float) -> float:
    """ELF/plasma penetration ratio for ceftazidime at plasma concentration cp."""
    if cp <= 0:
        return CAZ_EMAX / CAZ_EC50          # analytic limit as cp -> 0
    return caz_elf(cp) / cp


# --- Avibactam: power plasma-ELF link ------------------------------------------------------
# "At plasma avibactam concentrations relevant for efficacy (~1 mg/l), penetration into ELF
#  was 47%" -> EPR(1 mg/L) = 0.47.  The exponent is not printed, so it is solved from the
# paper's own second checkpoint (42% at 2.4 mg/L) and then validated against its third.
AVI_EPR1 = 0.47


def _solve_avi_pow() -> float:
    """Recover the power-model exponent from EPR(1)=0.47 and ratio(2.4 mg/L)=0.42."""
    # ratio(cp) = EPR1 * cp**(POW-1);  0.42 = 0.47 * 2.4**(POW-1)
    return 1.0 + math.log(0.42 / AVI_EPR1) / math.log(2.4)


AVI_POW = _solve_avi_pow()


def avi_ratio(cp: float) -> float:
    """ELF/plasma penetration ratio for avibactam at plasma concentration cp."""
    if cp <= 0:
        return float("nan")                 # power model diverges at zero
    return AVI_EPR1 * cp ** (AVI_POW - 1.0)


# --- Validation against every number the paper states --------------------------------------
# (label, function, plasma concentration, published value, tolerance)
CHECKPOINTS = [
    ("CAZ ELF conc at plasma 70 mg/L",  lambda: caz_elf(70.0),    22.5,  0.15),
    ("CAZ ratio at plasma 70 mg/L",     lambda: caz_ratio(70.0),  0.321, 0.005),
    ("CAZ ratio at plasma 15.3 mg/L",   lambda: caz_ratio(15.3),  0.52,  0.01),
    ("CAZ ELF conc at plasma 15.3 mg/L", lambda: caz_elf(15.3),    8.0,   0.15),
    ("CAZ ratio as plasma -> 0",        lambda: caz_ratio(0.0),   0.633, 0.005),
    ("CAZ ELF conc at plasma 71.7 mg/L", lambda: caz_elf(71.7),    22.7,  0.10),
    ("AVI ratio at plasma 1 mg/L",      lambda: avi_ratio(1.0),   0.47,  0.005),
    ("AVI ratio at plasma 2.4 mg/L",    lambda: avi_ratio(2.4),   0.42,  0.005),
    ("AVI ELF conc at plasma 12 mg/L",  lambda: avi_ratio(12.0) * 12.0, 4.0, 0.15),
    ("AVI ratio at plasma 12 mg/L",     lambda: avi_ratio(12.0),  0.332, 0.012),
]


def validate() -> bool:
    print("Reproducing Dimelow 2018 from its own published parameters")
    print("-" * 72)
    ok = True
    for label, fn, published, tol in CHECKPOINTS:
        got = fn()
        good = abs(got - published) <= tol
        ok = ok and good
        print(f"  [{'ok' if good else 'FAIL'}] {label:<36} model={got:8.4f}  paper={published:8.4f}")
    print("-" * 72)
    print(f"  avibactam power exponent solved from the paper's checkpoints: {AVI_POW:.4f}")
    return ok


# --- What the project actually applies, and where ------------------------------------------
# The fixed ratios used in revision_support/add_icu_elf_scenario.py and the therapeutic-window
# analysis. The ceiling is the ceftazidime neurotoxicity screen used throughout the package.
APPLIED_CAZ = 0.52
APPLIED_AVI = 0.42
NEUROTOX_SCREEN = 104.0   # mg/L ceftazidime

# Plasma concentrations spanning the range the continuous-infusion regimens actually produce.
CAZ_GRID = [8.0, 15.3, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 104.0, 120.0, 150.0]
AVI_GRID = [1.0, 2.4, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0]


def main() -> int:
    if not validate():
        print("\nParameterisation does not reproduce the source. Stopping rather than extrapolating.")
        return 1

    os.makedirs(OUT_DIR, exist_ok=True)

    print("\nCeftazidime: ratio applied by the project vs Dimelow's own model")
    print("-" * 72)
    caz_rows = []
    for cp in CAZ_GRID:
        r = caz_ratio(cp)
        rel = (APPLIED_CAZ - r) / r * 100.0
        caz_rows.append(dict(drug="ceftazidime", plasma_mg_L=f"{cp:.1f}",
                             dimelow_ratio=f"{r:.4f}", applied_ratio=f"{APPLIED_CAZ:.2f}",
                             overstatement_pct=f"{rel:.1f}"))
        mark = "  <- neurotoxicity screen" if cp == NEUROTOX_SCREEN else ""
        print(f"  plasma {cp:6.1f} mg/L   Dimelow {r:6.1%}   applied {APPLIED_CAZ:5.1%}"
              f"   overstates by {rel:6.1f}%{mark}")

    print("\nAvibactam: ratio applied by the project vs Dimelow's own model")
    print("-" * 72)
    avi_rows = []
    for cp in AVI_GRID:
        r = avi_ratio(cp)
        rel = (APPLIED_AVI - r) / r * 100.0
        avi_rows.append(dict(drug="avibactam", plasma_mg_L=f"{cp:.1f}",
                             dimelow_ratio=f"{r:.4f}", applied_ratio=f"{APPLIED_AVI:.2f}",
                             overstatement_pct=f"{rel:.1f}"))
        print(f"  plasma {cp:6.1f} mg/L   Dimelow {r:6.1%}   applied {APPLIED_AVI:5.1%}"
              f"   overstates by {rel:6.1f}%")

    path = os.path.join(OUT_DIR, "elf_penetration_vs_concentration.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["drug", "plasma_mg_L", "dimelow_ratio",
                                           "applied_ratio", "overstatement_pct"])
        w.writeheader()
        for row in caz_rows + avi_rows:
            w.writerow(row)

    print(f"\nwrote {os.path.relpath(path, HERE)}")

    print("\nBottom line")
    print("-" * 72)
    print(f"  The 0.52 ceftazidime ratio is the value at plasma 15.3 mg/L.")
    print(f"  At the 104 mg/L neurotoxicity screen the same model gives {caz_ratio(104.0):.1%}.")
    print(f"  The package's 'conservative' 0.30 scenario is closer to Dimelow's own model")
    print(f"  at continuous-infusion concentrations than its 'central estimate' 0.52 is.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
