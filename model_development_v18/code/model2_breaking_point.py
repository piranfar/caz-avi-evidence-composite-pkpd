"""How wrong would each input have to be before the recommendation changes?

WHY THIS IS DIFFERENT FROM A SENSITIVITY ANALYSIS
    A sensitivity analysis asks how much the OUTPUT moves when an input moves. That
    produces a tornado plot, and a reader still has to judge whether the movement
    matters. This asks the inverse and more useful question:

        for each input, what is the SMALLEST error that would change the decision,
        and is an error that large plausible?

    "Ceftazidime clearance would have to be wrong by 60% — nine times its published
    standard error, and larger than the spread across every published model — before
    the recommended regimen changes" is a far stronger sentence than a tornado plot,
    and it is what a reader actually wants to know.

TWO DECISIONS ARE TESTED
    A  which regimen to use in each renal-function class
    B  whether to measure avibactam or infer it from ceftazidime

    Each input is moved on its own, holding the others at their point estimates, and
    the scan reports the first crossing in each direction. Where an input has a
    published relative standard error the required error is also expressed in
    standard errors, and where the between-study heterogeneity applies it is
    expressed against that too. Those two comparators are what make the answer
    interpretable rather than merely numerical.
"""
from __future__ import annotations

import csv
import os
import sys
from dataclasses import replace

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

import model2_engine as E          # noqa: E402
import model2_hujam as H           # noqa: E402
import model2_monitoring as MON    # noqa: E402
import reproduce_primary_run as P  # noqa: E402

OUT = H.OUT
N_PER_CLASS = 20000
TAU_BETWEEN = 0.4244          # between-study SD on the log clearance scale, Layer 2

# published relative standard errors, Cojutti 2024 Table 2
RSE = dict(E.RSE)
RSE["rho"] = 0.238

# multiplicative parameters, scanned as a factor on the point estimate
MULTIPLICATIVE = ("cl0_caz", "cl0_avi", "omega_caz", "omega_avi",
                  "exp_caz", "exp_avi", "fu_caz", "fu_avi")


def baseline_choice(pop, weights, regimens, pr):
    res = E.evaluate(pop, pr, regimens)
    out = {}
    for cls, idx in group_by_class(regimens).items():
        u = [H.utility(res[regimens[j]], weights) for j in idx]
        out[cls] = regimens[idx[int(np.argmax(u))]]
    return out


def group_by_class(regimens):
    g = {}
    for j, r in enumerate(regimens):
        g.setdefault(P.REGIMENS[r][0], []).append(j)
    return g


def scan_regimen_decision(pop, weights, regimens, param, factors):
    """Smallest factor on `param`, in each direction, that changes the choice."""
    base = baseline_choice(pop, weights, regimens, E.BASE)
    flips = {cls: {"down": None, "up": None} for cls in base}
    for f in factors:
        if param in MULTIPLICATIVE:
            pr = replace(E.BASE, **{param: getattr(E.BASE, param) * f})
        elif param == "rho":
            pr = replace(E.BASE, rho=float(np.clip(E.BASE.rho * f, -0.999, 0.999)))
        elif param == "avi_target":
            pr = replace(E.BASE, avi_target=E.BASE.avi_target * f)
        elif param == "tox_threshold":
            pr = replace(E.BASE, tox_threshold=E.BASE.tox_threshold * f)
        else:
            raise KeyError(param)
        choice = baseline_choice(pop, weights, regimens, pr)
        for cls in base:
            if choice[cls] != base[cls]:
                side = "down" if f < 1 else "up"
                if flips[cls][side] is None:
                    flips[cls][side] = (f, choice[cls])
    return base, flips


def summarise_regimen(param, base, flips):
    rows = []
    rse = RSE.get(param)
    for cls, d in flips.items():
        for side in ("down", "up"):
            hit = d[side]
            if hit is None:
                rows.append({
                    "decision": "A: regimen choice", "parameter": param,
                    "ekfc_class": cls, "direction": side,
                    "baseline_choice": base[cls], "flips_to": "",
                    "required_factor": "", "required_change_pct": "",
                    "in_published_standard_errors": "",
                    "vs_between_study_sd": "",
                    "verdict": "no change across the whole scanned range"})
                continue
            f, to = hit
            pct = 100 * (f - 1)
            n_se = abs(np.log(f)) / rse if rse else None
            n_tau = abs(np.log(f)) / TAU_BETWEEN if param in ("cl0_caz", "cl0_avi") else None
            rows.append({
                "decision": "A: regimen choice", "parameter": param,
                "ekfc_class": cls, "direction": side,
                "baseline_choice": base[cls], "flips_to": to,
                "required_factor": round(f, 3),
                "required_change_pct": round(pct, 1),
                "in_published_standard_errors": round(n_se, 1) if n_se else "",
                "vs_between_study_sd": round(n_tau, 2) if n_tau else "",
                "verdict": ""})
    return rows


def scan_monitoring(pop, rng, assay_cv):
    """Correlation at which inferring becomes as good as measuring."""
    out = []
    for rho in np.arange(0.30, 0.995, 0.02):
        pr = replace(E.BASE, rho=float(rho), avi_target=4.0)
        m = MON.classify(pop, pr, rng, assay_cv, assay_cv)
        out.append((float(rho), m["accuracy_measure"] - m["accuracy_infer"]))
    cross = None
    for i in range(len(out) - 1):
        if (out[i][1] > 0) != (out[i + 1][1] > 0):
            x0, y0 = out[i]
            x1, y1 = out[i + 1]
            cross = x0 + (0 - y0) * (x1 - x0) / (y1 - y0)
            break
    return out, cross


def main():
    pop = E.draw_population(N_PER_CLASS, H.MASTER_SEED)
    dists = E.load_mic_distributions()
    weights = dists[H.PRIMARY_DIST]["w"]
    regimens = list(P.REGIMENS)

    print("=" * 80)
    print("BREAKING-POINT ANALYSIS — how wrong would each input have to be?")
    print("=" * 80)

    base = baseline_choice(pop, weights, regimens, E.BASE)
    print("\n  Baseline choice at the point estimates "
          "(net benefit, lambda = 1, target 4 mg/L)")
    for cls, r in base.items():
        print(f"    {cls:9} {r}")

    factors = np.concatenate([np.arange(0.30, 1.00, 0.01),
                              np.arange(1.01, 4.01, 0.01)])
    factors = np.sort(factors)
    down = factors[factors < 1][::-1]      # walk outward from 1 in each direction
    up = factors[factors > 1]

    rows = []
    print("\n  Decision A — which regimen\n")
    print(f"  {'parameter':12} {'class':9} {'flip at':>9} {'change':>9} "
          f"{'in SEs':>8} {'vs tau':>8}  becomes")
    print("  " + "-" * 74)
    for param in ("cl0_caz", "cl0_avi", "exp_caz", "exp_avi", "omega_caz",
                  "omega_avi", "rho", "fu_caz", "fu_avi", "avi_target",
                  "tox_threshold"):
        b, flips = scan_regimen_decision(pop, weights, regimens, param,
                                         np.concatenate([down, up]))
        r = summarise_regimen(param, b, flips)
        rows += r
        shown = [x for x in r if x["required_factor"] != ""]
        if not shown:
            print(f"  {param:12} {'—':9} {'never':>9}  no error in the scanned "
                  f"range (0.30x to 4x) changes any choice")
            continue
        for x in sorted(shown, key=lambda z: abs(float(z["required_change_pct"])))[:2]:
            print(f"  {param:12} {x['ekfc_class']:9} "
                  f"{x['required_factor']:>9.2f} {x['required_change_pct']:>8.1f}% "
                  f"{str(x['in_published_standard_errors']):>8} "
                  f"{str(x['vs_between_study_sd']):>8}  {x['flips_to']}")

    print("\n  Decision B — measure avibactam, or infer it\n")
    rng = np.random.default_rng(H.MASTER_SEED + 5)
    for cv in (0.0, 0.10, 0.20, 0.30):
        curve, cross = scan_monitoring(pop, rng, cv)
        lo = curve[0][1]
        hi = curve[-1][1]
        if cross is None:
            msg = (f"measuring wins at every correlation from 0.30 to 0.99 "
                   f"(gain {hi:.2f} to {lo:.2f} pp)")
        else:
            msg = f"inferring catches up only at rho = {cross:.3f}"
        print(f"    assay CV {100*cv:3.0f}%  {msg}")
        rows.append({
            "decision": "B: measure or infer", "parameter": "rho",
            "ekfc_class": "all (selected regimens)", "direction": "up",
            "baseline_choice": "measure avibactam",
            "flips_to": "infer" if cross else "",
            "required_factor": round(cross / E.BASE.rho, 3) if cross else "",
            "required_change_pct": round(100 * (cross / E.BASE.rho - 1), 1) if cross else "",
            "in_published_standard_errors": "", "vs_between_study_sd": "",
            "verdict": (f"crossover at rho = {cross:.3f}" if cross
                        else "measuring wins across the entire plausible range"),
            "assay_cv_pct": round(100 * cv)})

    os.makedirs(OUT, exist_ok=True)
    fields = ["decision", "parameter", "ekfc_class", "direction", "baseline_choice",
              "flips_to", "required_factor", "required_change_pct",
              "in_published_standard_errors", "vs_between_study_sd", "verdict",
              "assay_cv_pct"]
    path = os.path.join(OUT, "model2_breaking_points.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\n  wrote {os.path.relpath(path, os.path.dirname(E.HERE))} "
          f"({len(rows)} rows)")
    print("\n  Comparators: 'in SEs' is the required error divided by the parameter's")
    print("  published relative standard error; 'vs tau' divides it by the")
    print("  between-study SD of 0.424 on the log clearance scale, i.e. by the spread")
    print("  across the four published population PK models.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
