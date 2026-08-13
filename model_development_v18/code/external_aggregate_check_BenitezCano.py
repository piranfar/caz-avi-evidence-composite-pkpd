"""External aggregate-level predictive check of the primary model against an independent
ICU continuous-infusion cohort (Benitez-Cano 2026, Crit Care 30:305, CC BY-NC-ND).

Observed values are group summary statistics published in that article. No individual
patient data are used, because none are published. This is therefore an EXTERNAL
AGGREGATE-LEVEL PREDICTIVE CHECK, not external clinical validation.
"""
import numpy as np, sys
sys.stdout.reconfigure(encoding="utf-8")

# --- primary model, exactly as in reproduce_primary_run.py ---
CL0_CAZ, EXP_CAZ = 5.0, 0.70
CL0_AVI, EXP_AVI = 5.9, 0.89
EKFC_REF = 70.0
CV_CAZ, CV_AVI = 0.6792, 0.7691
RHO = 0.94
OM_CAZ = np.sqrt(np.log(1 + CV_CAZ**2))
OM_AVI = np.sqrt(np.log(1 + CV_AVI**2))
FU_CAZ, FU_AVI = 0.85, 0.92

# --- observed cohort (Benitez-Cano 2026), CAZ/AVI arm ---
N_OBS = 15
EGFR_MED = 63.0                 # mL/min/1.73 m2, CKD-EPI, median (Table 1)
DAILY_CAZ_MG, DAILY_AVI_MG = 6000.0, 1500.0   # 2 g / 0.5 g q8h infused over 8 h = CI
R_CAZ, R_AVI = DAILY_CAZ_MG / 24.0, DAILY_AVI_MG / 24.0
# observed plasma AUC(0-8,ss), median [IQR width], mg*h/L  ->  Css = AUC/8
OBS = {
    "ceftazidime": dict(auc_med=648.0, auc_iqr=368.8),
    "avibactam":   dict(auc_med=85.6,  auc_iqr=37.6),
}
for v in OBS.values():
    v["css_med"] = v["auc_med"] / 8.0
    v["css_iqr"] = v["auc_iqr"] / 8.0

print("=" * 78)
print("EXTERNAL AGGREGATE-LEVEL PREDICTIVE CHECK")
print("Primary model (Cojutti-derived) vs Benitez-Cano 2026 ICU cohort, n = 15")
print("Regimen 6 g / 1.5 g per day by continuous infusion; cohort median eGFR 63 mL/min/1.73 m2")
print("=" * 78)

# --- 1. deterministic check at the cohort median renal function -----------------
# Under CI at steady state Css = R/CL, and CL is lognormal about CL_typ,
# so the MEDIAN predicted Css is exactly R/CL_typ. No simulation needed.
cl_caz_typ = CL0_CAZ * (EGFR_MED / EKFC_REF) ** EXP_CAZ
cl_avi_typ = CL0_AVI * (EGFR_MED / EKFC_REF) ** EXP_AVI
pred = {"ceftazidime": R_CAZ / cl_caz_typ, "avibactam": R_AVI / cl_avi_typ}

print("\n1. Predicted vs observed MEDIAN total steady-state concentration\n")
print(f"{'analyte':14} {'CL_typ':>8} {'pred Css':>10} {'obs Css':>10} {'obs/pred':>9} {'PE %':>8}")
print("-" * 63)
rows = []
for a, cl in (("ceftazidime", cl_caz_typ), ("avibactam", cl_avi_typ)):
    p, o = pred[a], OBS[a]["css_med"]
    pe = 100.0 * (p - o) / o
    rows.append((a, cl, p, o, o / p, pe))
    print(f"{a:14} {cl:8.3f} {p:10.1f} {o:10.1f} {o/p:9.2f} {pe:+8.1f}")

# --- 2. simulation with interindividual variability ----------------------------
# eGFR spread is not fully identified from a median plus an IQR width, so it is
# drawn lognormal about the reported median with the reported IQR width, and the
# result is reported alongside a fixed-eGFR variant as a bracket.
rng = np.random.default_rng(20260707)
NSIM = 200_000
cov = np.array([[OM_CAZ**2, RHO * OM_CAZ * OM_AVI],
                [RHO * OM_CAZ * OM_AVI, OM_AVI**2]])
eta = rng.multivariate_normal(np.zeros(2), cov, size=NSIM)

# lognormal eGFR with median 63 and interquartile RANGE of 53 mL/min
sigma_egfr = np.log((EGFR_MED + 53 / 2) / (EGFR_MED - 53 / 2)) / (2 * 0.6744898)
egfr = EGFR_MED * np.exp(rng.normal(0, sigma_egfr, NSIM))
egfr = np.clip(egfr, 5, 200)

print(f"\n2. With interindividual variability (N = {NSIM:,}; eGFR lognormal, median 63, IQR 53)\n")
print(f"{'analyte':14} {'pred median':>12} {'pred IQR':>20} {'obs median':>11} {'obs IQR width':>14}")
print("-" * 76)
sim = {}
for a, cl0, ex, om, R in (("ceftazidime", CL0_CAZ, EXP_CAZ, 0, R_CAZ),
                          ("avibactam",   CL0_AVI, EXP_AVI, 1, R_AVI)):
    cl = cl0 * (egfr / EKFC_REF) ** ex * np.exp(eta[:, om])
    css = R / cl
    sim[a] = css
    q1, q2, q3 = np.percentile(css, [25, 50, 75])
    print(f"{a:14} {q2:12.1f} {q1:9.1f} - {q3:8.1f} {OBS[a]['css_med']:11.1f} {OBS[a]['css_iqr']:14.1f}")
    print(f"{'':14} {'':12} {'IQR width ' + format(q3-q1, '.1f'):>20}")

# --- 3. what clearance would reproduce the observed exposure? ------------------
print("\n3. Clearance implied by the observed exposure, and the source's own estimate\n")
print(f"{'analyte':14} {'model CL_typ':>13} {'implied CL':>11} {'reported CL':>12} {'ratio m/r':>10}")
print("-" * 64)
REPORTED = {"ceftazidime": 2.86, "avibactam": 6.08}   # Benitez-Cano Table 2 (Monolix, 1-cpt)
for a, cl in (("ceftazidime", cl_caz_typ), ("avibactam", cl_avi_typ)):
    R = R_CAZ if a == "ceftazidime" else R_AVI
    implied = R / OBS[a]["css_med"]
    print(f"{a:14} {cl:13.2f} {implied:11.2f} {REPORTED[a]:12.2f} {cl/REPORTED[a]:10.2f}")

# --- 4. consequence for target attainment --------------------------------------
print("\n4. Consequence for attainment in this cohort (free concentrations)\n")
for a, fu, tgt, lbl in (("ceftazidime", FU_CAZ, 4 * 8.0, "fCss >= 4 x MIC at MIC 8 mg/L"),
                        ("avibactam",   FU_AVI, 4.0,     "fCss >= 4 mg/L")):
    pred_att = 100.0 * np.mean(sim[a] * fu >= tgt)
    scaled = sim[a] * (OBS[a]["css_med"] / np.median(sim[a]))
    obs_att = 100.0 * np.mean(scaled * fu >= tgt)
    print(f"{a:14} {lbl:32} model {pred_att:5.1f}%   exposure-corrected {obs_att:5.1f}%")

print("\nNOTE: the exposure-corrected column rescales the simulated distribution so its median")
print("matches the observed median. It is an illustration of the direction and size of the")
print("bias, not a re-estimated model.")
