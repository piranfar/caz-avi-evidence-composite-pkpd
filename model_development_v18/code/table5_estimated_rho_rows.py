"""Non-RRT continuous-infusion comparator: classification of avibactam attainment from a single
ceftazidime concentration, at the source-model correlation and at the CRRT estimate and its bounds.

The classifier is the one defined in target_checks.py (second_assay_proper):
a population of critically ill adults without renal replacement therapy, simulated from Cojutti et al
(J Antimicrob Chemother 2024;79:2801-8), receiving continuous infusion, with the measured ceftazidime
steady-state concentration carrying 0%, 10% or 20% error. Only the assumed correlation changes
between rows. The rho = 0.94 rows are the comparator used in the CRRT analysis; the rows at the CRRT
estimate and its profile-likelihood bounds (read from outputs/model1_final_parameters.csv) show what
the CRRT value would imply if it applied here, which the CRRT analysis does not assume.

The rho = 0.94, CV = 0 row is a control and must reproduce the published values.
Values are written at full precision (6 decimals) and rounded once, when they are reported.

Writes model_development_v18/outputs/table5_estimated_rho_rows.csv.
"""

from __future__ import annotations

import csv
import os
import sys
from math import erf, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# Local layout: revision_support/; public repository layout: src/cazavi/.
for _p in (os.path.join(HERE, '..', '..', 'revision_support'), os.path.join(HERE, '..', '..', 'src'),
           os.path.join(HERE, '..', '..', 'src', 'cazavi')):
    sys.path.insert(0, _p)

from cazavi_analyses import SELECTED_REGIMENS, _cholesky, draw_population   # noqa: E402
from reproduce_primary_run import (                                          # noqa: E402
    AVI_CT, AVI_FRACTION, CAZ_FRACTION, CL0_AVI, CL0_CAZ, EKFC_REF, EXP_AVI, EXP_CAZ,
    FU_AVI, N_PER_CLASS, OMEGA_AVI, OMEGA_CAZ, PRIMARY_SEED, REGIMENS,
)

OUT = os.path.join(HERE, '..', 'outputs')
DECIMALS = 6


def crrt_rho():
    with open(os.path.join(OUT, 'model1_final_parameters.csv'), encoding='utf-8') as fh:
        r = next(x for x in csv.DictReader(fh) if x['parameter'] == 'corr_CL_caz_avi')
    return float(r['estimate']), float(r['ci_low']), float(r['ci_high'])


def clearances(renal, z, rho):
    eta = z @ _cholesky(OMEGA_CAZ, OMEGA_AVI, rho).T
    return (CL0_CAZ * (renal / EKFC_REF) ** EXP_CAZ * np.exp(eta[:, 0]),
            CL0_AVI * (renal / EKFC_REF) ** EXP_AVI * np.exp(eta[:, 1]))


def css(regimen, cl_caz, cl_avi):
    _cls, dose_g, interval = REGIMENS[regimen]
    return (dose_g * 1000.0 * CAZ_FRACTION / interval / cl_caz,
            dose_g * 1000.0 * AVI_FRACTION / interval / cl_avi)


def operating_characteristics(rho, role, cv, rng):
    """The classifier exactly as defined in target_checks.second_assay_proper."""
    pop = draw_population(N_PER_CLASS, PRIMARY_SEED)
    tp = fp = tn = fn = 0
    for regimen in SELECTED_REGIMENS:
        cls = REGIMENS[regimen][0]
        renal, z = pop[cls]
        cl_caz, cl_avi = clearances(renal, z, rho)
        c_caz, c_avi = css(regimen, cl_caz, cl_avi)
        truth = c_avi * FU_AVI >= AVI_CT

        obs = (c_caz * np.exp(rng.normal(0.0, np.sqrt(np.log(1 + cv ** 2)), c_caz.size))
               if cv > 0 else c_caz)

        typ_caz = CL0_CAZ * (renal / EKFC_REF) ** EXP_CAZ
        _cls, dose_g, interval = REGIMENS[regimen]
        cl_caz_hat = dose_g * 1000.0 * CAZ_FRACTION / interval / obs
        eta_caz_hat = np.log(cl_caz_hat / typ_caz)
        mu = rho * (OMEGA_AVI / OMEGA_CAZ) * eta_caz_hat
        sd = OMEGA_AVI * np.sqrt(max(1.0 - rho ** 2, 1e-12))

        typ_avi = CL0_AVI * (renal / EKFC_REF) ** EXP_AVI
        rate_avi = dose_g * 1000.0 * AVI_FRACTION / interval
        thresh = np.log((rate_avi * FU_AVI / AVI_CT) / typ_avi)
        p_attain = 0.5 * (1.0 + np.vectorize(erf)((thresh - mu) / sd / sqrt(2.0)))
        pred = p_attain >= 0.5

        tp += int(np.sum(pred & truth)); fp += int(np.sum(pred & ~truth))
        tn += int(np.sum(~pred & ~truth)); fn += int(np.sum(~pred & truth))

    n = tp + fp + tn + fn
    r = lambda v: round(v, DECIMALS)
    return {
        'rho': rho if rho == 0.94 else r(rho), 'role': role, 'assay_cv_pct': round(100 * cv),
        'accuracy_pct': r(100.0 * (tp + tn) / n),
        'sensitivity_pct': r(100.0 * tp / max(tp + fn, 1)),
        'specificity_pct': r(100.0 * tn / max(tn + fp, 1)),
        'ppv_pct': r(100.0 * tp / max(tp + fp, 1)),
        'npv_pct': r(100.0 * tn / max(tn + fn, 1)),
        'false_reassurance_pct': r(100.0 * fp / n),
        'prevalence_attaining_pct': r(100.0 * (tp + fn) / n),
    }


def main() -> int:
    est, lo, hi = crrt_rho()
    rhos = ((0.94, 'source-model value'), (hi, 'CRRT estimate, upper profile-likelihood bound'),
            (est, 'CRRT estimate'), (lo, 'CRRT estimate, lower profile-likelihood bound'))
    rows = []
    print('%8s %-46s %4s %9s %6s %6s %6s %6s %7s' % (
        'rho', 'role', 'CV', 'accuracy', 'sens', 'spec', 'PPV', 'NPV', 'falseRe'))
    for rho, role in rhos:
        rng = np.random.default_rng(PRIMARY_SEED + 31)
        for cv in (0.0, 0.10, 0.20):
            r = operating_characteristics(rho, role, cv, rng)
            rows.append(r)
            print('%8.4f %-46s %3d%% %8.2f%% %5.2f%% %5.2f%% %5.2f%% %5.2f%% %6.3f%%' % (
                r['rho'], r['role'], r['assay_cv_pct'], r['accuracy_pct'], r['sensitivity_pct'],
                r['specificity_pct'], r['ppv_pct'], r['npv_pct'], r['false_reassurance_pct']))

    ctrl = next(r for r in rows if r['rho'] == 0.94 and r['assay_cv_pct'] == 0)
    published = dict(accuracy_pct=94.1, sensitivity_pct=97.2, specificity_pct=77.0,
                     ppv_pct=95.8, npv_pct=83.6, false_reassurance_pct=3.6)
    drift = {k: (ctrl[k], v) for k, v in published.items() if round(ctrl[k], 1) != v}
    print()
    print('control row rho=0.94, CV=0 reproduces the published values:',
          'yes' if not drift else 'NO -> %s' % drift)

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, 'table5_estimated_rho_rows.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    print('wrote', os.path.relpath(path, HERE))
    return 0 if not drift else 1


if __name__ == '__main__':
    sys.exit(main())
