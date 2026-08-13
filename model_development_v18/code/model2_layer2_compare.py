"""MODEL 2, LAYER 2 — the effect of including between-study heterogeneity.

Every result is reported WITH and WITHOUT the layer, and with a leave-one-study-out
check on tau, because with four non-exchangeable studies the estimate is fragile.
See the header of `model2_heterogeneity.py` for why this is a scenario, not a result.
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

import model2_engine as E          # noqa: E402
import model2_hujam as H           # noqa: E402
import model2_heterogeneity as HET  # noqa: E402


def compare(n_draws, n_per_class, taus, rho="C1_cojutti", target="T4_uniform"):
    rows, summary = [], []
    for label, tau in taus:
        r = H.run(n_draws, n_per_class, rho, target, tau_between=tau)
        att = H.integrated_attainment(r)
        _, s = H.optimality_and_regret(r)
        widths = [x["prediction_interval_width_pp"] for x in att]
        for x in att:
            x["tau_scenario"] = label
            x["tau_between"] = tau
            rows.append(x)
        summary.append({
            "tau_scenario": label, "tau_between": round(tau, 4),
            "rho_scenario": rho, "target_scenario": target,
            "median_prediction_interval_width_pp": round(float(np.median(widths)), 2),
            "max_prediction_interval_width_pp": round(float(np.max(widths)), 2),
            "mean_misselection_pct": round(
                float(np.mean([x["p_misselection_pct"] for x in s])), 1),
            "mean_evpi_pp": round(float(np.mean([x["evpi_pp"] for x in s])), 3),
            "mean_expected_regret_pp": round(
                float(np.mean([x["expected_regret_pp"] for x in s])), 3),
        })
        print(f"    {label:34} tau {tau:.4f}   PI width "
              f"{np.median(widths):6.2f} pp   misselection "
              f"{np.mean([x['p_misselection_pct'] for x in s]):5.1f}%   EVPI "
              f"{np.mean([x['evpi_pp'] for x in s]):6.3f} pp")
    return rows, summary


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    n_draws = 300 if args.quick else 2000
    n_per_class = 2000 if args.quick else 4000

    print("=" * 80)
    print("MODEL 2, LAYER 2 — with and without between-study heterogeneity")
    print("=" * 80)
    data, models = HET.load_typical()
    tau, resid, _ = HET.estimate_tau(data, models, verbose=False)
    loo = HET.leave_one_out(data, models)
    tau_lo = min(r["tau_log"] for r in loo)
    tau_hi = max(r["tau_log"] for r in loo)
    print(f"\n  tau = {tau:.4f} (between-study CV "
          f"{100*np.sqrt(np.exp(tau**2)-1):.1f}%), leave-one-out range "
          f"{tau_lo:.4f} to {tau_hi:.4f}\n")

    taus = [("off (primary analysis)", 0.0),
            ("estimated tau", tau),
            (f"leave-one-out low ({tau_lo:.3f})", tau_lo),
            (f"leave-one-out high ({tau_hi:.3f})", tau_hi)]
    rows, summary = compare(n_draws, n_per_class, taus)

    out = os.path.join(os.path.dirname(E.HERE), "outputs")
    for data_rows, name in ((rows, "model2_layer2_attainment.csv"),
                            (summary, "model2_layer2_summary.csv")):
        path = os.path.join(out, name)
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data_rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(data_rows)
        print(f"  wrote {os.path.relpath(path, os.path.dirname(E.HERE))} "
              f"({len(data_rows)} rows)")

    base = summary[0]["median_prediction_interval_width_pp"]
    est = summary[1]["median_prediction_interval_width_pp"]
    print(f"\n  Including between-study heterogeneity widens the median prediction")
    print(f"  interval from {base:.1f} to {est:.1f} percentage points "
          f"({est/base:.1f}-fold).")
    print("  Report both. With four non-exchangeable studies, tau is a scenario.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
