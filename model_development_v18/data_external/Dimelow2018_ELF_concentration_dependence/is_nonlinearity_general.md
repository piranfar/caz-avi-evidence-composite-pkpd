# Is the ELF non-linearity general, or specific to ceftazidime/avibactam?

**Checked 12 August 2026, because `NOVELTY_STRATEGY.md` Part II §5 flagged it as the prerequisite for
any error-class paper. Answer: NOT general. The broad claim does not survive. A narrower one does.**

Literature retrieved from PubMed.

---

## The claim under test

Part II §2.2 proposed: *the fixed ELF penetration ratio is a systematically biased input to lung
target attainment wherever the plasma-ELF relationship is non-linear.* The conditional clause was
always there, but the ranked table scored it 9–10/10 on the implicit assumption that non-linearity is
common. **That assumption is wrong.**

## Evidence AGAINST generality — drugs where a constant ratio fits

| drug | finding | source |
|---|---|---|
| **imipenem** | ELF distribution described by a **time-independent penetration coefficient of 0.44** (RSE 14%), pooled across healthy volunteers, elderly and renal impairment | van Hasselt 2016, *Br J Clin Pharmacol* [10.1111/bcp.12901](https://doi.org/10.1111/bcp.12901) |
| **lefamulin** | ELF "well described using **first-order rate constants** into and out of the ELF compartment" — linear, despite non-linear *plasma* protein binding | Zhang 2019, *JAC* [10.1093/jac/dkz088](https://doi.org/10.1093/jac/dkz088) |
| **ceftaroline** | two-compartment model, penetration ≈23%, no non-linearity reported | Riccobene 2016, *AAC* [10.1128/AAC.02755-15](https://doi.org/10.1128/AAC.02755-15) |
| **sulbactam/durlobactam** | four-compartment model with **linear kinetics**; ELF sub-models per analyte | Cammarata 2024, *AAC* [10.1128/aac.00485-24](https://doi.org/10.1128/aac.00485-24) |

**imipenem is the decisive counter-example.** van Hasselt explicitly fitted and reported a
*time-independent* coefficient across three populations. A fixed ratio is not a shortcut there; it is
the model that fits.

## Evidence FOR failure of the fixed ratio — but for varied and specific reasons

| drug | how the fixed ratio fails | source |
|---|---|---|
| **ceftazidime / avibactam** | plasma-ELF link **saturable** (Michaelis-Menten) for CAZ, **power** for AVI; ratio falls with concentration | Dimelow 2018 [10.1007/s40268-018-0241-0](https://doi.org/10.1007/s40268-018-0241-0) |
| **polymyxin B** | **saturable binding** in ELF; total-drug ELF AUC high, but modelled **unbound** ELF AUC only ~16.7% of plasma total | Jiao 2023, *AAC* [10.1128/aac.00197-23](https://doi.org/10.1128/aac.00197-23) |
| **cefiderocol** | ratio **population-dependent** — 34% in pneumonia patients, **1.4× that of healthy subjects**, with delayed distribution | Kawaguchi 2022, *J Clin Pharmacol* [10.1002/jcph.1986](https://doi.org/10.1002/jcph.1986) |
| **piperacillin / tazobactam** | concentrations "**unpredictable** and **negatively correlated with pulmonary permeability**"; penetration **49.3%** (pip) vs **121.2%** (tazo) | Felton 2014, *Clin Pharmacol Ther* [10.1038/clpt.2014.131](https://doi.org/10.1038/clpt.2014.131) |
| **ceftobiprole** | median 30%, **range 15–45%** — threefold spread across 12 patients | Roger 2025, *JAC* [10.1093/jac/dkaf267](https://doi.org/10.1093/jac/dkaf267) |

Separately, Rouby 2024 (*J Intensive Med* [10.1016/j.jointm.2024.07.006](https://doi.org/10.1016/j.jointm.2024.07.006))
argues that for nebulised aminoglycosides ELF concentrations "grossly overestimate lung interstitial
fluid concentrations" because of bronchial contamination during bronchoscopy, and that lung
microdialysis is the only accurate technique. That is a challenge to the ELF *measurement*, not to
the modelling of it, and is out of scope here — but it should not be forgotten.

## The verdict

**The broad claim is dead.** Non-linearity is not a general property of plasma-ELF relationships;
several β-lactams are well described by a constant coefficient, and one (imipenem) was explicitly
tested and found time-independent. Writing "the field's fixed ratios are systematically biased" would
be false, and the imipenem paper alone would refute it.

**What survives is narrower, better evidenced, and closer to this project's actual subject:**

> For **β-lactam / β-lactamase-inhibitor combinations**, the two components do not penetrate ELF
> alike, and the joint PK/PD target is routinely evaluated by applying one fixed ratio per component —
> ratios usually taken from healthy volunteers, at concentrations unlike those continuous infusion
> produces. That practice is **untested in most published analyses** and **demonstrably wrong in at
> least two combinations**: ceftazidime/avibactam, where the relationship is saturable, and
> piperacillin/tazobactam, where the components differ by 2.5-fold (49.3% vs 121.2%) and penetration
> tracks pulmonary permeability rather than dose.

Two supporting observations make this more than a technicality:

1. **The inhibitor is usually the limiting component** (`Cojutti2026_FN_avibactam_limiting/`), so an
   error in *its* ratio propagates straight into the joint target — the component least likely to have
   been characterised carefully is the one that decides attainment.
2. **Where the link function has actually been tested against alternatives, non-linearity has won
   about as often as it has lost** — Dimelow tested three and chose saturable and power; van Hasselt
   tested and chose constant. The problem is not that the field is wrong; it is that **most analyses
   never test at all**, and a single ratio is reported as though it were a property of the drug rather
   than a fit to one concentration range in one population.

## What this changes for the project

- **The standalone "error-class" paper as originally conceived should not be written.** Its central
  claim would be false as stated.
- A narrower paper on **fixed-ratio ELF scaling in BL/BLI combinations** is defensible and is a better
  fit for this project's data and expertise. Realistic novelty **7–8/10**, not 9–10/10.
- **§5.2 of `MODEL2_REPORT.md` is unaffected.** That section is about *this package's own* use of
  Dimelow's ratios at concentrations where they do not hold, which is established independently of
  whether other drugs behave the same way. It remains correct and remains the more consequential
  finding, because it changes a stated conclusion.
- Nothing here rescues or damages any other route.

## Honest limits of this check

- **This is a targeted search, not a systematic review.** Two PubMed queries, roughly 35 abstracts
  screened, five drugs on each side. A systematic review could shift the balance, and the counts here
  should not be quoted as prevalences.
- **Mostly abstracts, not full texts.** Where an abstract does not mention non-linearity, that is weak
  evidence of a linear fit and no evidence at all that alternatives were tested — the same
  "nobody tests" problem, now applying to this check as well.
- **"Linear fits well" is sample-dependent.** Imipenem's constant coefficient was fitted over the
  concentration range those studies sampled; it says nothing about behaviour far outside it. The
  asymmetry cuts both ways and is not resolved here.
