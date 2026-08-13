"""MODEL 2, LAYER 2 — between-study heterogeneity across published population PK models.

WHAT THIS REPLACES
    The manuscript runs four published population PK models one at a time and reports
    the range of the results (population joint CFR 70.2% to 94.5%). A range across
    four models is a description, not an uncertainty statement: it has no coverage
    interpretation, and it gets wider simply by adding another model.

    This layer instead treats the typical clearance as drawn from a distribution
    across studies,

        log CL(study s, class c, analyte a) = mu(c, a) + u_s ,   u_s ~ N(0, tau^2)

    and integrates joint attainment over u. The result is a prediction interval.

ESTIMATION
    A two-way layout with study as a random main effect and class-by-analyte cell as
    a fixed effect, tau^2 recovered by the ANOVA method of moments. The estimand is
    the STUDY-level effect, shared across analytes, because a model that predicts
    high clearance predicts it for both components -- which is what the data show.
    The study-by-cell interaction is retained as the residual term and reported
    separately.

    Restricted to the three renal classes all four models cover (61-90, 91-120,
    121-150), so the layout is balanced.

WHY THIS LAYER IS REPORTED AS A SCENARIO, NOT A RESULT
    Three reasons, and all three must travel with any number it produces.

    1. FOUR STUDIES. Standard meta-analysis guidance is that a between-study
       variance estimated from fewer than about five studies is unreliable. tau
       here has roughly three degrees of freedom.
    2. THE STUDIES ARE NOT EXCHANGEABLE. Critically ill adults on continuous
       infusion (Italy), adults with carbapenem-resistant K. pneumoniae (China),
       pooled phase 1-3 registrational subjects, and adults with cystic fibrosis and
       preserved renal function. Treating these as draws from one distribution is a
       strong assumption; the cystic fibrosis cohort is arguably not from the same
       population at all.
    3. DIFFERENT RENAL DESCRIPTORS. EKFC, Cockcroft-Gault and creatinine clearance
       are compared as if interchangeable at the same class boundary. They are not
       exactly interchangeable.

    Every result is therefore reported WITH and WITHOUT this layer, and with a
    leave-one-study-out check.
"""
from __future__ import annotations

import csv
import io
import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

import model2_engine as E   # noqa: E402

SRC = os.path.join(E.V16, "outputs", "structural_typical_clearances.csv")
COMMON_CLASSES = ("61–90", "91–120", "121–150")


def load_typical():
    raw = open(SRC, "rb").read()
    for enc in ("utf-8", "cp1252"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    rows = list(csv.DictReader(io.StringIO(text)))
    data = {}
    for r in rows:
        if r["ekfc_class"] not in COMMON_CLASSES:
            continue
        data[(r["model"], r["ekfc_class"], "caz")] = float(r["typical_cl_caz_l_h"])
        data[(r["model"], r["ekfc_class"], "avi")] = float(r["typical_cl_avi_l_h"])
    models = sorted({k[0] for k in data})
    return data, models


def estimate_tau(data, models, verbose=True):
    """ANOVA method of moments for the study main effect on the log scale."""
    cells = [(c, a) for c in COMMON_CLASSES for a in ("caz", "avi")]
    Y = np.array([[np.log(data[(m, c, a)]) for (c, a) in cells] for m in models])
    S, J = Y.shape

    grand = Y.mean()
    row = Y.mean(axis=1)          # study means
    col = Y.mean(axis=0)          # cell means

    ms_study = J * float(((row - grand) ** 2).sum()) / (S - 1)
    resid = Y - row[:, None] - col[None, :] + grand
    ms_resid = float((resid ** 2).sum()) / ((S - 1) * (J - 1))
    tau2 = max((ms_study - ms_resid) / J, 0.0)
    tau = float(np.sqrt(tau2))

    if verbose:
        print(f"    studies {S}, cells {J} (3 renal classes x 2 analytes)")
        print(f"    MS study {ms_study:.5f}   MS residual {ms_resid:.5f}")
        print(f"    tau (study, log scale) {tau:.4f}"
              f"   -> between-study CV {100*np.sqrt(np.exp(tau**2)-1):.1f}%")
        print(f"    residual (study x cell) SD {np.sqrt(ms_resid):.4f}"
              f"   -> CV {100*np.sqrt(np.exp(ms_resid)-1):.1f}%")
        for m, r in zip(models, row):
            print(f"      {m:22} mean log CL deviation from grand mean {r-grand:+.4f}"
                  f"   (x{np.exp(r-grand):.2f})")
    return tau, float(np.sqrt(ms_resid)), models


def leave_one_out(data, models):
    out = []
    for drop in models:
        keep = [m for m in models if m != drop]
        tau, res, _ = estimate_tau(data, keep, verbose=False)
        out.append({"dropped_study": drop, "n_studies": len(keep),
                    "tau_log": round(tau, 4),
                    "between_study_cv_pct": round(100 * np.sqrt(np.exp(tau ** 2) - 1), 1)})
    return out


def main():
    print("=" * 76)
    print("MODEL 2, LAYER 2 — between-study heterogeneity")
    print("=" * 76)
    data, models = load_typical()
    print("\n  Estimating tau from the published typical clearances")
    tau, resid_sd, _ = estimate_tau(data, models)

    print("\n  Leave-one-study-out")
    loo = leave_one_out(data, models)
    for r in loo:
        print(f"    without {r['dropped_study']:22} tau {r['tau_log']:.4f}"
              f"   CV {r['between_study_cv_pct']:5.1f}%   (n={r['n_studies']})")
    spread = max(r["tau_log"] for r in loo) - min(r["tau_log"] for r in loo)
    print(f"\n    tau ranges over {spread:.4f} across leave-one-out refits"
          f" — with four studies this is the expected instability, not a defect")

    os.makedirs(E.OUT if hasattr(E, "OUT") else "../outputs", exist_ok=True)
    out_dir = os.path.join(os.path.dirname(E.HERE), "outputs")
    path = os.path.join(out_dir, "model2_heterogeneity_tau.csv")
    rows = [{"quantity": "tau_study_log_scale", "value": round(tau, 4),
             "interpretation": f"between-study CV "
                               f"{100*np.sqrt(np.exp(tau**2)-1):.1f}%",
             "note": "ANOVA method of moments, 4 studies, 3 shared renal classes, "
                     "2 analytes. UNRELIABLE with fewer than about five studies."},
            {"quantity": "residual_study_by_cell_sd", "value": round(resid_sd, 4),
             "interpretation": f"CV {100*np.sqrt(np.exp(resid_sd**2)-1):.1f}%",
             "note": "study-by-cell interaction, retained as the residual term"}]
    rows += [{"quantity": f"tau_without_{r['dropped_study']}", "value": r["tau_log"],
              "interpretation": f"between-study CV {r['between_study_cv_pct']}%",
              "note": f"leave-one-study-out, n={r['n_studies']}"} for r in loo]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"\n  wrote {os.path.relpath(path, os.path.dirname(E.HERE))}")
    return tau


if __name__ == "__main__":
    main()
