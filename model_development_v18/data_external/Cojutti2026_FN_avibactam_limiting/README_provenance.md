# Independent clinical evidence that avibactam is the limiting component

## Source

Cojutti PG, et al. **Likelihood of aggressive PK/PD target attainment of continuous-infusion
beta-lactams during the first week of treatment of febrile neutropenia: findings from a 1-year
prospective, monocentric study in onco-haematological patients.**
*J Antimicrob Chemother* 2026;81(7):dkag183.
[doi:10.1093/jac/dkag183](https://doi.org/10.1093/jac/dkag183) — PMID 42234486.
OUP, "all rights reserved". Retrieved 12 August 2026 through institutional access during the
systematic JAC sweep. No PDF archived; extracted numbers only.

## Why this one matters more than its size suggests

The cohort is 256 patients, but **only 17 (6.6%) received ceftazidime/avibactam** — 10 targeted,
7 empirical. That is a small subgroup. What makes it valuable is not n; it is that **both analytes
were measured in every one of them, and the paper reports which component failed.**

| | n | CAZ fCss (mg/L) | AVI fCss (mg/L) | joint target attained |
|---|---|---|---|---|
| targeted | 10 | 21.0 (16.6–21.4) | 6.4 (4.4–10.3) | 8/10 (80%) |
| empirical | 7 | 16.2 (9.1–29.8) | 3.3 (2.0–5.6) | **1/7 (14.3%)** |

Medians with min–max range, as printed.

**The paper's own attribution of failure:**

> targeted: "Non-attainment was always due to insufficient concentrations of the BLI in both BL/BLIc,
> namely tazobactam (n = 10) or **avibactam (n = 3)**."
>
> empirical: "Target non-attainment was ... due to insufficient concentrations of the BLI for
> piperacillin/tazobactam and for ceftazidime/avibactam, namely tazobactam (n = 43/92) and
> **avibactam (n = 6/7)**."

## What this supports, and what it does not

**It is direct empirical support for the project's central premise.** Model 2's
`limiting_probability()` asks which component constrains the joint target. This paper answers that
question with measured concentrations in real patients and gets **avibactam, overwhelmingly** — 6 of
7 empirical failures and every targeted failure. A monitoring strategy that measured ceftazidime
alone would have missed all of them.

That is also the sharpest available answer to round 3 of the
[JAC exchange](../JAC_exchange_measure_one_or_both/): Fresan et al. argue avibactam measurement is
not useful because the product is a fixed 4:1 combination. Here, in an independent centre and a
different population, avibactam is the component that fails. The fixed ratio does not make the
avibactam concentration predictable; it makes the failures invisible to single-analyte TDM.

**It is not confirmation of Model 2's numbers.** This is a different population (febrile neutropenia,
onco-haematological), a different target denominator (EUCAST clinical breakpoint for empirical
therapy), and n=17. It corroborates the *premise* — that the inhibitor is often limiting — not any
quantitative output.

## A second, independent citation for C_T = 4 mg/L

The project's avibactam threshold of 4 mg/L has until now rested on the Bologna chain (Gatti 2023).
This paper states the basis independently:

> "The target concentration ratio corresponded to the fixed BLI target concentration defined by the
> EUCAST for testing the in vitro standard susceptibility of each BL/BLIc, namely, **4 mg/L for
> tazobactam or avibactam**."

So C_T = 4 mg/L is an EUCAST susceptibility-testing concentration, not a Bologna-specific choice.
Worth saying in the manuscript — it removes a single-source dependency.

## Unbound fractions — a third published pair

Ceftazidime protein binding 10%, **avibactam 7%** ⇒ fu = 0.90 / 0.93, cited to Sy 2019. The published
pairs now number three: 0.85/0.92 (Dimelow, used by this project), 0.90/0.93 (Gatti 2023 and here),
0.90/0.95 (Falcone 2021). **All three fall inside the project's existing sensitivity ranges**
(fu_CAZ 0.80–0.90, fu_AVI 0.87–0.97), so this is a point in the project's favour, not a gap.

## Discrepancy recorded as printed

For targeted therapy the table gives 8/10 attained (2 failures), but the text attributes avibactam
insufficiency to n = 3. These cannot both be counts of the same thing. Recorded as-printed without
resolution — the same convention used for the Fresan 2023 mortality discrepancy. Do not cite the
n = 3 as a failure count for ceftazidime/avibactam without checking the source table.

## Other findings worth carrying

- **Augmented renal clearance is the dominant risk for non-attainment**: OR 12.29 (3.31–45.57),
  P < 0.001. The project's renal-function strata are the relevant axis, and this says the top of that
  range is where attainment breaks.
- Attainment was **not** associated with improved clinical outcome in either arm (73.4% vs 71.7%
  empirical; 67.5% vs 60.0% targeted). A useful caution: this is another study where the PK/PD target
  is achieved or missed without a detectable outcome signal, which is exactly why the project must not
  translate value-of-information results into lives.
- Dosing was **eGFR-banded continuous infusion** (2.5 g LD then 2.5 g q8h if eGFR > 50), i.e. the same
  regimen family the project simulates.

## Reuse conditions

OUP all-rights-reserved. Cite and link; no PDF redistributed. Quotations above are short and
attributed. The extracted table is in `cazavi_target_attainment.csv`.
