"""CRRT virtual-population analysis: target attainment and ceftazidime-only classification.

Simulates a virtual CRRT population from the fitted joint model (Model 1: 4 random effects, the
clearance correlation rho and the volume correlation r_V, read from outputs/model1_joint_popk_parameters.csv
and outputs/model1_final_parameters.csv) and asks, within that population:

  1. What joint target attainment do common CRRT regimens give?
  2. How well does a single measured ceftazidime trough classify avibactam target attainment, at the
     fitted clearance correlation, at its profile-likelihood bounds, and at the non-RRT value of 0.94?

The classifier is the Bayes posterior P(avibactam attains | observed ceftazidime trough, regimen)
computed from the model itself on a grid of observed values; a patient is classified as attaining when
it is >= 0.5. Measurement error on the ceftazidime trough is the model's own proportional residual error
for ceftazidime. A single ceftazidime measurement is used, so the residual correlation between the 2
drugs does not enter the classifier. No assay CV scenario is added.

Population and dosing: intermittent 8- or 12-hourly 2-h infusion, steady state, one compartment per
drug, exactly as fitted. Parameter uncertainty is not propagated here (see crrt_parameter_uncertainty.py);
the correlation is varied explicitly.

The three original regimens use the random stream seeded with SEED. The 2.5 g every 12 h regimen (the
RRT regimen proposed by Hu et al., Drug Des Devel Ther 2026) uses its own stream seeded with SEED + 1.

Values are written at full precision (6 decimals) and rounded once, when they are reported.

Outputs: outputs/crrt_attainment.csv, outputs/crrt_classifier.csv
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "outputs")


def _rows(name):
    with open(os.path.join(OUT, name), encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def load_model1():
    """Model 1 estimates (joint_popk_nlme.py) and the profile-likelihood interval of rho
    (model1_finalise.py). Read at import, so no estimate is copied into this file by hand."""
    est = {r["parameter"]: float(r["estimate"]) for r in _rows("model1_joint_popk_parameters.csv")}
    fin = {r["parameter"]: r for r in _rows("model1_final_parameters.csv")}
    rho = fin["corr_CL_caz_avi"]
    assert abs(float(rho["estimate"]) - est["corr_CL_caz_avi"]) < 1e-6, "Model 1 outputs come from different fits"
    return dict(
        th=dict(cl_caz=est["CL ceftazidime (L/h)"], cl_avi=est["CL avibactam (L/h)"],
                v_caz=est["V ceftazidime (L)"], v_avi=est["V avibactam (L)"]),
        om=dict(cl_caz=est["omega_CL_caz"], cl_avi=est["omega_CL_avi"],
                v_caz=est["omega_V_caz"], v_avi=est["omega_V_avi"]),
        r_v=est["corr_V_caz_avi"], sig_caz=est["sigma_prop_caz"],
        rho=(float(rho["estimate"]), float(rho["ci_low"]), float(rho["ci_high"])))


MODEL1 = load_model1()
TH, OM = MODEL1["th"], MODEL1["om"]
R_V = MODEL1["r_v"]
SIG_CAZ = MODEL1["sig_caz"]
RHO_EST, RHO_LO, RHO_HI = MODEL1["rho"]

FU_CAZ, FU_AVI = 0.85, 0.92          # unbound fractions, as in the primary analysis
T_INF = 2.0                          # h
AVI_CT = 4.0                         # mg/L free, primary avibactam target
AVI_CT_REG = 1.0                     # mg/L free, registrational target (50% fT>CT)
CAZ_K = 4.0                          # ceftazidime fCmin/MIC >= 4
TOX = 104.0                          # mg/L total ceftazidime trough, exposure screen
MICS = (2, 4, 8, 16)
N = 100_000
SEED = 20261002

# (label, caz mg, avi mg, tau h)
REGIMENS = [
    ("2.5 g q8h", 2000.0, 500.0, 8.0),
    ("1.25 g q8h", 1000.0, 250.0, 8.0),
    ("0.94 g q12h", 750.0, 187.5, 12.0),
]
# Simulated in a separate random stream (seed SEED + 1); see the module docstring.
EXTRA_REGIMENS = [
    ("2.5 g q12h", 2000.0, 500.0, 12.0),
]
EXTRA_SEED = SEED + 1
DECIMALS = 6


def css(t, cl, v, dose, tau, t_inf=T_INF):
    """Steady-state one-compartment concentration, repeated IV infusion (as Model 1)."""
    k = cl / v
    rate = dose / t_inf
    c0 = (rate / cl) * (1 - np.exp(-k * t_inf)) * np.exp(-k * (tau - t_inf)) / (1 - np.exp(-k * tau))
    t_in = np.minimum(t, t_inf)
    during = (rate / cl) * (1 - np.exp(-k * t_in)) + c0 * np.exp(-k * t_in)
    c_end = (rate / cl) * (1 - np.exp(-k * t_inf)) + c0 * np.exp(-k * t_inf)
    return np.where(t <= t_inf, during, c_end * np.exp(-k * np.maximum(t - t_inf, 0)))


def population(rng, n, rho, om_scale=1.0, th=None, om=None, r_v=None):
    """Virtual patients. Clearance deviates correlated by rho, volume deviates by r_V.
    th, om and r_v default to the Model 1 estimates; crrt_parameter_uncertainty.py passes draws."""
    th, om = th or TH, om or OM
    r_v = R_V if r_v is None else r_v
    z1 = rng.standard_normal(n)
    z2 = rho * z1 + np.sqrt(1 - rho ** 2) * rng.standard_normal(n)
    z3 = rng.standard_normal(n)
    z4 = r_v * z3 + np.sqrt(1 - r_v ** 2) * rng.standard_normal(n)
    s = om_scale
    return dict(
        cl_caz=th["cl_caz"] * np.exp(s * om["cl_caz"] * z1),
        cl_avi=th["cl_avi"] * np.exp(s * om["cl_avi"] * z2),
        v_caz=th["v_caz"] * np.exp(s * om["v_caz"] * z3),
        v_avi=th["v_avi"] * np.exp(s * om["v_avi"] * z4),
    )


def exposures(p, caz_mg, avi_mg, tau, fu_avi=FU_AVI):
    cmin_caz = css(tau, p["cl_caz"], p["v_caz"], caz_mg, tau)
    cmin_avi = css(tau, p["cl_avi"], p["v_avi"], avi_mg, tau)
    tg = np.linspace(0, tau, 97)[:, None]
    avi_prof = fu_avi * css(tg, p["cl_avi"][None, :], p["v_avi"][None, :], avi_mg, tau)
    ft_avi1 = (avi_prof >= AVI_CT_REG).mean(axis=0)
    return cmin_caz, cmin_avi, ft_avi1


def posterior_curve(rng, rho, caz_mg, avi_mg, tau, om_scale, grid, n_prior=400_000, fu_avi=FU_AVI,
                    sig_caz=None, **pop):
    """P(avibactam fCmin >= target | observed log ceftazidime trough), on a grid."""
    sig = SIG_CAZ if sig_caz is None else sig_caz
    p = population(rng, n_prior, rho, om_scale, **pop)
    cmin_caz, cmin_avi, _ = exposures(p, caz_mg, avi_mg, tau)
    att = (fu_avi * cmin_avi >= AVI_CT).astype(float)
    lc = np.log(cmin_caz)
    out = np.empty(grid.size)
    for i, g in enumerate(grid):
        u = (g - lc) / sig
        w = np.exp(-0.5 * (u * u - np.min(u * u)))      # shifted for numerical stability; cancels in the ratio
        out[i] = (w * att).sum() / w.sum()
    return out


def classify(rng, rho_true, rho_model, caz_mg, avi_mg, tau, om_scale=1.0):
    p = population(rng, N, rho_true, om_scale)
    cmin_caz, cmin_avi, _ = exposures(p, caz_mg, avi_mg, tau)
    truth = FU_AVI * cmin_avi >= AVI_CT
    obs = np.log(cmin_caz) + SIG_CAZ * rng.standard_normal(N)
    grid = np.linspace(obs.min() - 0.05, obs.max() + 0.05, 300)
    curve = posterior_curve(rng, rho_model, caz_mg, avi_mg, tau, om_scale, grid)
    pred = np.interp(obs, grid, curve) >= 0.5
    return summary(truth, pred)


def summary(truth, pred):
    n = truth.size
    tp = np.sum(pred & truth); fp = np.sum(pred & ~truth)
    tn = np.sum(~pred & ~truth); fn = np.sum(~pred & truth)
    f = lambda a, b: 100 * a / b if b else float("nan")
    return dict(prevalence_attain=f(truth.sum(), n), accuracy=f(tp + tn, n),
                sensitivity=f(tp, tp + fn), specificity=f(tn, tn + fp),
                ppv=f(tp, tp + fp), npv=f(tn, tn + fn), false_reassurance=f(fp, n),
                missed_by_population_prior=f((~truth).sum(), n))


def attainment_row(rng, om_scale, lab, cz, av, tau):
    p = population(rng, N, RHO_EST, om_scale)
    cmin_caz, cmin_avi, ft1 = exposures(p, cz, av, tau)
    avi4 = FU_AVI * cmin_avi >= AVI_CT
    row = dict(omega_scale=om_scale, regimen=lab,
               caz_fcmin_median=np.median(FU_CAZ * cmin_caz),
               avi_fcmin_median=np.median(FU_AVI * cmin_avi),
               avi_fcmin_p5=np.percentile(FU_AVI * cmin_avi, 5),
               avi_attain_fcmin4=100 * avi4.mean(),
               avi_attain_50ft1=100 * (ft1 >= 0.5).mean(),
               caz_total_cmin_gt104=100 * (cmin_caz > TOX).mean())
    for m in MICS:
        cz_ok = FU_CAZ * cmin_caz >= CAZ_K * m
        row[f"caz_pta_mic{m}"] = 100 * cz_ok.mean()
        row[f"joint_pta_mic{m}"] = 100 * (cz_ok & avi4).mean()
    return row


MISSPECIFIED = "misspecified: model 0.94, truth estimate"
SCENARIOS = [("estimate", RHO_EST, RHO_EST), ("upper bound", RHO_HI, RHO_HI),
             ("lower bound", RHO_LO, RHO_LO), ("non-RRT value", 0.94, 0.94),
             (MISSPECIFIED, RHO_EST, 0.94)]


def write(path, rows_):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_[0]), lineterminator="\n")
        w.writeheader()
        w.writerows({k: (round(float(v), DECIMALS) if isinstance(v, (float, np.floating)) else v)
                     for k, v in r.items()} for r in rows_)


def main():
    print(f"Model 1: rho {RHO_EST:.6f} ({RHO_LO:.6f} to {RHO_HI:.6f}); r_V {R_V:.6f}; sigma_caz {SIG_CAZ:.6f}")
    rng = np.random.default_rng(SEED)
    att_rows = [attainment_row(rng, om, lab, cz, av, tau)
                for om in (1.0, 1.5, 2.0) for lab, cz, av, tau in REGIMENS]
    cls_rows = []
    for om_scale in (1.0, 1.5, 2.0):
        for lab, cz, av, tau in REGIMENS:
            for name, rt, rm in SCENARIOS:
                r = classify(rng, rt, rm, cz, av, tau, om_scale)
                cls_rows.append(dict(omega_scale=om_scale, regimen=lab, scenario=name,
                                     rho_true=rt, rho_model=rm, **r))

    # Additional regimen in its own stream.
    rng2 = np.random.default_rng(EXTRA_SEED)
    att_rows += [attainment_row(rng2, om, lab, cz, av, tau)
                 for om in (1.0, 1.5, 2.0) for lab, cz, av, tau in EXTRA_REGIMENS]
    for om_scale in (1.0, 1.5, 2.0):
        for lab, cz, av, tau in EXTRA_REGIMENS:
            for name, rt, rm in SCENARIOS:
                r = classify(rng2, rt, rm, cz, av, tau, om_scale)
                cls_rows.append(dict(omega_scale=om_scale, regimen=lab, scenario=name,
                                     rho_true=rt, rho_model=rm, **r))

    write(os.path.join(OUT, "crrt_attainment.csv"), att_rows)
    write(os.path.join(OUT, "crrt_classifier.csv"), cls_rows)

    for r in att_rows:
        print({k: (round(float(v), 2) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    for r in cls_rows:
        print({k: (round(float(v), 2) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
