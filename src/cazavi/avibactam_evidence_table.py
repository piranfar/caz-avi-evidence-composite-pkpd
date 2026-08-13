"""Where the experimental evidence for an avibactam threshold actually sits.

The paper shows that the avibactam critical concentration governs the answer.
That raises the question the sweep alone cannot settle: which part of the
1-8 mg/L range has experimental support behind it.

This assembles what the pharmacodynamic literature reports, one row per study,
recording the avibactam concentration examined, the organism and its
beta-lactamase, and whether killing was maintained. The pattern is the point:
the dynamic in vitro work that first identified a threshold puts it near
0.25-0.5 mg/L, the registrational analyses adopted 1 mg/L, and the value used as
a steady-state target in continuous-infusion practice is 4 mg/L, which entered
from susceptibility testing rather than from an exposure-response experiment.

Sources retrieved from PubMed; DOIs are recorded with each row.

Usage:
    python avibactam_evidence_table.py
"""

from __future__ import annotations

import csv
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")

# concentration_mg_l is the avibactam concentration the study associates with
# maintained killing or with the onset of regrowth; role records what kind of
# number it is, because these are not interchangeable.
EVIDENCE = [
    dict(study="Coleman 2014", design="Hollow-fibre, dynamic",
         organism="Enterobacteriaceae (8 strains)",
         enzyme="ESBL, AmpC, KPC-class serine enzymes",
         concentration_mg_l=0.3, lower_mg_l=0.25, upper_mg_l=0.5,
         role="Regrowth threshold: growth resumed once avibactam fell below ~0.3 mg/L",
         killing_maintained="Yes above threshold",
         doi="10.1128/AAC.00080-14"),
    dict(study="Coleman 2014", design="Hollow-fibre, pulsed exposure",
         organism="Enterobacteriaceae", enzyme="Serine beta-lactamases",
         concentration_mg_l=0.375, lower_mg_l=0.25, upper_mg_l=0.5,
         role="Pulse of >0.25 and <0.5 mg/L suppressed growth for 24 h with 1 g ceftazidime",
         killing_maintained="Yes",
         doi="10.1128/AAC.00080-14"),
    dict(study="Coleman 2014", design="Hollow-fibre, single dose",
         organism="High-level AmpC producer", enzyme="Derepressed AmpC",
         concentration_mg_l=0.5, lower_mg_l=0.5, upper_mg_l=0.5,
         role="1 g/250 mg sufficient for 7 of 8 strains; 2 g/500 mg needed for AmpC",
         killing_maintained="Only at the higher profile",
         doi="10.1128/AAC.00080-14"),
    dict(study="Sy & Derendorf 2017", design="Static time-kill with PK/PD model",
         organism="Enterobacteriaceae", enzyme="Serine beta-lactamases",
         concentration_mg_l=None, lower_mg_l=None, upper_mg_l=None,
         role="Modelling framework; effect driven by EC50 of the inhibitor, not a fixed threshold",
         killing_maintained="Concentration-dependent",
         doi="10.1016/j.cmi.2017.07.020"),
    dict(study="Bensman 2017", design="Population PK with target-attainment analysis",
         organism="Pseudomonas aeruginosa (cystic fibrosis)", enzyme="Mixed",
         concentration_mg_l=1.0, lower_mg_l=1.0, upper_mg_l=1.0,
         role="Target expressed as cumulative time above a 1 mg/L threshold",
         killing_maintained="Target adopted, not tested",
         doi="10.1128/AAC.00988-17"),
    dict(study="Das 2019", design="Registrational PK/PD, dose selection",
         organism="Enterobacterales, P. aeruginosa", enzyme="Serine beta-lactamases",
         concentration_mg_l=1.0, lower_mg_l=1.0, upper_mg_l=1.0,
         role="Joint target used to select and validate the licensed regimen: free time above 1 mg/L",
         killing_maintained="Validated against phase 3 outcomes",
         doi="10.1128/AAC.02187-18"),
    dict(study="Sy 2019", design="Model-informed development review",
         organism="Enterobacterales, P. aeruginosa", enzyme="Serine beta-lactamases",
         concentration_mg_l=1.0, lower_mg_l=1.0, upper_mg_l=1.0,
         role="States 1 mg/L as the PD threshold and 4 mg/L as the fixed testing concentration",
         killing_maintained="Review",
         doi="10.1007/s40262-018-0705-y"),
    dict(study="Zhang 2022", design="Clinical case series with TDM",
         organism="Gram-negative, renal replacement therapy", enzyme="Mixed",
         concentration_mg_l=1.0, lower_mg_l=1.0, upper_mg_l=1.0,
         role="Trough above 1 mg/L used as the adequacy criterion",
         killing_maintained="Clinical response reported",
         doi="10.1007/s40121-022-00621-z"),
    dict(study="Kroemer 2023", design="Hollow-fibre with semi-mechanistic model",
         organism="Escherichia coli", enzyme="CTX-M-15, TEM-4, OXA-244",
         concentration_mg_l=None, lower_mg_l=None, upper_mg_l=None,
         role="Synergy characterised by model parameters rather than a fixed threshold",
         killing_maintained="Combination-dependent",
         doi="10.1128/spectrum.03318-23"),
    dict(study="Gatti 2022", design="Clinical TDM, continuous infusion",
         organism="Difficult-to-treat P. aeruginosa", enzyme="Mixed",
         concentration_mg_l=4.0, lower_mg_l=4.0, upper_mg_l=4.0,
         role="Steady-state target of Css above 4 mg/L applied in continuous-infusion practice",
         killing_maintained="Microbiological eradication reported",
         doi="10.3390/antibiotics11121739"),
    dict(study="EUCAST", design="Susceptibility testing convention",
         organism="Enterobacterales, P. aeruginosa", enzyme="n/a",
         concentration_mg_l=4.0, lower_mg_l=4.0, upper_mg_l=4.0,
         role="Fixed avibactam concentration at which ceftazidime MICs are determined",
         killing_maintained="Not an exposure-response value",
         doi="EUCAST breakpoint tables v14.0"),
]

# Joint CFR from the model's own sweep, for the selected regimens.
SWEEP = {1.0: (93.4, 98.6), 2.0: (89.5, 97.6), 4.0: (73.5, 88.6),
         6.0: (58.1, 74.7), 8.0: (44.0, 60.7)}


def main():
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "avibactam_threshold_evidence.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(EVIDENCE[0]))
        w.writeheader()
        w.writerows(EVIDENCE)
    print(f"  wrote {len(EVIDENCE)} rows  outputs/avibactam_threshold_evidence.csv")

    dynamic = [e for e in EVIDENCE if "Hollow-fibre" in e["design"] and e["concentration_mg_l"]]
    lo = min(e["lower_mg_l"] for e in dynamic)
    hi = max(e["upper_mg_l"] for e in dynamic)
    print(f"\n  dynamic in vitro evidence spans {lo}–{hi} mg/L")
    print(f"  registrational threshold        1 mg/L  ({1/hi:.0f}–{1/lo:.0f}x the experimental range)")
    print(f"  continuous-infusion target      4 mg/L  ({4/hi:.0f}–{4/lo:.0f}x the experimental range)")
    print("\n  joint CFR at each threshold (selected regimens, KPC-KP):")
    for ct, (a, b) in SWEEP.items():
        mark = ""
        if ct == 1.0:
            mark = "  <- registrational"
        elif ct == 4.0:
            mark = "  <- primary target of this analysis"
        print(f"    {ct:g} mg/L: {a:5.1f}–{b:5.1f}%{mark}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
