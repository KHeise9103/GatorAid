# eICU Electrolyte Replacement Analysis

## Project purpose

This analysis reproduces the GatorAId electrolyte-replacement framework in the independent **eICU dataset** for potassium, magnesium, and phosphorus. The reporting structure mirrors the Mayo analysis using the same general 3 x 2 framework:

**Below / Within / Above reference range x Repletion / No repletion.**

The summarized eICU output is an external validation/replication analysis. Medication exposure is based on **medication order/start timing**, not Mayo-style MAR-confirmed administration.

**Reporting convention:** `NR` means the value was not reported in the current summarized eICU output. NR cells are intentionally left unfilled; no values are imputed.

---

# 1. Overall eICU analysis populations

| Measure | Potassium | Magnesium | Phosphorus |
|---|---:|---:|---:|
| Final eligible ICU stays | 711 | 737 | 737 |
| Stays with electrolyte labs | NR | 532 | 347 |
| Electrolyte lab events | 4,336 | 2,244 | 1,281 |
| Replacement-linked events | 556 | 266 | 64 |
| Overall replacement rate | **12.8%** | **11.9%** | **5.0%** |

---

# 2. Potassium

## Cohort summary

| Measure | Result |
|---|---:|
| Final eligible ICU stays | 711 |
| Potassium lab events | 4,336 |
| Replacement-linked lab events | 556 |
| Overall replacement rate | 12.8% |
| Loop-exposed potassium events +/-24 h | 897 / 4,336 (20.7%) |

## Primary 3 x 2 replacement behavior

| Baseline potassium group | Repletion | % within group | Mean Delta K |
|---|---|---:|---:|
| Below <3.5 | Repletion | **35.5%** | +0.468 |
| Below <3.5 | No repletion | 64.5%* | +0.286 |
| Within 3.5-5.0 | Repletion | **8.0%** | +0.179 |
| Within 3.5-5.0 | No repletion | 92.0%* | -0.090 |
| Above >5.0 | Repletion | **0.6%** | -0.500** |
| Above >5.0 | No repletion | 99.4%* | -0.571 |

\* No-repletion percentages are arithmetic complements of the reported repletion percentages.  
\** Only one treated Above-range potassium event had paired follow-up.

## Follow-up outcomes

| Baseline group | Repletion | Normal follow-up | Remained low | Follow-up high |
|---|---|---:|---:|---:|
| Below | Repletion | **58.4%** | 39.6% | 2.0% |
| Below | No repletion | 46.4% | NR | NR |
| Within | Repletion | NR | NR | NR |
| Within | No repletion | NR | NR | NR |
| Above | Repletion | NR | NR | NR |
| Above | No repletion | NR | NR | NR |

## Paired potassium change

| Baseline group | Repletion | Paired n | Mean Delta K | Paired Cohen's d | Effect |
|---|---|---:|---:|---:|---|
| Below | Repletion | NR | +0.468 | NR | NR |
| Below | No repletion | NR | +0.286 | NR | NR |
| Within | Repletion | NR | +0.179 | NR | NR |
| Within | No repletion | NR | -0.090 | NR | NR |
| Above | Repletion | 1 paired follow-up | -0.500 | NR | Not reliable |
| Above | No repletion | NR | -0.571 | NR | NR |

## Direct effect of repletion on Delta K

| Baseline group | Repletion Delta K | No-repletion Delta K | Cohen's d | Effect |
|---|---:|---:|---:|---|
| Below | +0.468 | +0.286 | **0.412** | Small |
| Within | +0.179 | -0.090 | **0.638** | Moderate |
| Above | -0.500 | -0.571 | Not reliable | Insufficient treated n |

## Six-group ANOVA

- Groups compared: **6**
- F statistic: **133.24**
- p-value: **<0.001**
- eta-squared: **0.193**
- Interpretation: **Large overall effect**
- Paired observation count: **NR**

## Loop-diuretic analysis and replacement threshold

- Loop-exposed potassium events +/-24 h: **897 / 4,336 (20.7%)**
- Relevant loop terms were **not statistically significant** in the eICU overall model.
- Detailed odds ratios / 95% CI: **NR**
- Estimated 50% replacement threshold, no loop: **2.88 mmol/L**
- Estimated 50% replacement threshold, loop: **3.11 mmol/L**

Interpretation: the estimated potassium threshold shifted upward with loop exposure, but because the loop terms were not statistically significant, the threshold difference is **exploratory**.

---

# 3. Magnesium

## Cohort and analytic sample

| Measure | Result |
|---|---:|
| Final eligible ICU stays | 737 |
| Stays with magnesium labs | 532 |
| Magnesium lab events | 2,244 |
| Replacement-linked Mg lab events | 266 |
| No-repletion Mg lab events | 1,978 |
| Overall replacement rate | **11.9%** |

## Treatment patterns and Delta Mg

The current summarized eICU output does not report the cell counts or within-category percentages for the six magnesium groups. Reported mean changes are:

| Baseline Mg group | Repletion | Mean Delta Mg |
|---|---|---:|
| Below | Repletion | +0.598 |
| Below | No repletion | +0.672 |
| Within | Repletion | +0.317 |
| Within | No repletion | +0.043 |
| Above | Repletion | NR |
| Above | No repletion | NR |

Follow-up counts, baseline means, follow-up means, and detailed outcome proportions are NR in the summarized output; treated Above-range follow-up was described as limited.

## Primary 3 x 2 ANOVA

- Groups compared: **6**
- F statistic: **83.88**
- p-value: **<0.001**
- eta-squared: **0.277**
- Effect size: **Large**
- Paired observations: **NR**

## Direct repletion effect on Delta Mg

| Baseline Mg group | Repletion Delta Mg | No-repletion Delta Mg | Difference | Cohen's d | Holm p | Effect |
|---|---:|---:|---:|---:|---:|---|
| Below | +0.598 | +0.672 | **-0.074** | NR | 0.565 | Not significant |
| Within | +0.317 | +0.043 | **+0.274** | **0.955** | <0.001 | Large |
| Above | NR | NR | NR | NR | NR | Limited treated observations |

The low-magnesium treatment contrast was not significant after Holm correction, whereas the Within-range treatment-associated difference was large and statistically significant.

## Loop-diuretic model and replacement threshold

- Loop-exposed Mg events +/-24 h: **NR**
- Relevant loop terms were **not statistically significant** in the eICU overall model.
- Detailed odds ratios / 95% CI: **NR**
- Estimated replacement threshold, no loop: **1.285 mg/dL**
- Estimated replacement threshold, loop: **1.519 mg/dL**
- Reportability: **Exploratory**

Interpretation: the estimated threshold shifted upward with loop exposure, but the loop terms were not statistically significant. The Below-range magnesium treatment contrast did not clearly replicate the Mayo finding in this smaller eICU sample.

---

# 4. Phosphorus

## Cohort summary

| Measure | Result |
|---|---:|
| Final eligible ICU stays | 737 |
| Stays with phosphorus labs | 347 |
| Phosphorus lab events | 1,281 |
| Replacement-linked lab events | 64 |
| No-repletion lab events | 1,217 |
| Overall replacement rate | **5.0%** |

## Primary treatment pattern

| Baseline phosphorus group | Repletion | Mean Delta P |
|---|---|---:|
| Below <3.0 | Repletion | +0.903 |
| Below <3.0 | No repletion | +0.289 |
| Within 3.0-4.5 | Repletion | NR |
| Within 3.0-4.5 | No repletion | NR |
| Above >4.5 | Repletion | NR - insufficient treated paired follow-up |
| Above >4.5 | No repletion | NR |

The summarized eICU output does not report the six-group cell counts, within-category percentages, baseline means, or follow-up means for phosphorus.

## Follow-up normalization

| Baseline group | Repletion | Normal follow-up |
|---|---|---:|
| Below | Repletion | **38.9%** |
| Below | No repletion | **22.5%** |
| Within | Repletion | NR |
| Within | No repletion | NR |
| Above | Repletion | NR |
| Above | No repletion | NR |

For Below-range phosphorus, the absolute normalization difference was **+16.4 percentage points** in favor of repletion. A p-value and Cohen's h were not reported for this normalization comparison.

## Direct treatment effect on Delta P

| Baseline group | Repletion Delta P | No-repletion Delta P | Difference | Cohen's d | Holm p |
|---|---:|---:|---:|---:|---:|
| Below | +0.903 | +0.289 | **+0.614** | **0.797** | <0.001 |
| Within | NR | NR | NR | NR | NR |
| Above | NR | NR | NR | NR | Insufficient treated paired follow-up |

The Below-range phosphorus treatment effect was large-to-moderate in magnitude and statistically significant.

## Six-group ANOVA

- Groups compared: **6**
- F statistic: **33.62**
- p-value: **<0.001**
- eta-squared: **0.192**
- Effect size: **Large**
- Paired observations: **NR**

## Loop-diuretic and threshold analysis

- Loop-exposed phosphorus events +/-24 h: **NR**
- Detailed odds ratios / 95% CI: **NR**
- 50% replacement threshold: **No validated threshold supported by the observed eICU distribution**
- No-loop threshold: **No validated threshold**
- Loop threshold: **No validated threshold**

Interpretation: phosphorus showed a clear treatment-associated signal despite the smaller sample, but no eICU phosphorus clinical replacement threshold should be reported from the current model.

---

# 5. Mayo vs eICU replication summary

| Electrolyte | Mayo finding | eICU finding | Replication |
|---|---|---|---|
| Potassium | Replacement improves Delta K / normalization | Same direction; d=0.412 | **Supported** |
| Magnesium | Replacement improves low-Mg outcomes | Low-Mg Delta difference not significant | **Not clearly replicated** |
| Phosphorus | Strong Delta P treatment effect; d=1.037 | Delta P difference +0.614; d=0.797 | **Supported** |

Across all three electrolytes, the overall six-group models were significant with large effect sizes (**eta-squared 0.192-0.277**). Potassium and phosphorus showed directional replication of the Mayo low-electrolyte treatment signal. Magnesium showed a large overall category-by-treatment effect but did not show a significant Below-range treatment contrast in the summarized eICU analysis.

---

# 6. Key findings

1. **Potassium:** replacement occurred in **35.5%** of Below-range events, **8.0%** of Within-range events, and **0.6%** of Above-range events. Low-K repletion was associated with a larger mean increase in potassium and a higher reported normalization rate. The direct low-K effect size was **d=0.412 (small)**.

2. **Magnesium:** the overall six-group model showed a large effect (**eta-squared=0.277**), but the Below-range treatment contrast was not significant (**Holm p=0.565**). The Within-range treatment-associated difference was large (**d=0.955, Holm p<0.001**).

3. **Phosphorus:** Below-range repletion was associated with a substantially larger mean increase (**+0.903 vs +0.289 mg/dL; difference +0.614; d=0.797; Holm p<0.001**) and higher normalization (**38.9% vs 22.5%**).

4. **Loop-diuretic thresholds were exploratory.** Potassium and magnesium threshold estimates shifted upward with loop exposure, but the relevant loop terms were not statistically significant. No validated phosphorus threshold was supported.

5. **Cross-dataset replication:** potassium and phosphorus supported the direction of the Mayo low-electrolyte treatment signal; magnesium did not clearly replicate the low-Mg treatment contrast.

---

# 7. Important limitations

- eICU medication exposure is based on **medication order/start timing**, not MAR-confirmed administration. The eICU analysis should therefore be described as an external validation/replication rather than an exact medication-administration replication of Mayo.
- Many detailed cell counts, baseline means, follow-up means, paired observation counts, and model coefficients are **NR** in the current summarized eICU output. No values should be inferred or imputed for these cells.
- The Above-range potassium treated group had only **one paired follow-up**, so its treatment-effect estimate is not reliable.
- Above-range magnesium had limited treated follow-up, and several detailed comparisons are not reportable from the current summary.
- Phosphorus had a small replacement-linked sample (**64 events**) and insufficient treated paired follow-up in the Above range.
- Loop-diuretic threshold estimates for potassium and magnesium are exploratory because the relevant loop terms were not statistically significant.
- Statistical significance in the six-group models should be interpreted alongside effect sizes and the availability of paired follow-up data.

---

# 8. Repository note

For a public or shared GitHub repository, include **analysis code, documentation, mappings, and aggregate outputs**. Keep patient-level source data and any restricted row-level derivatives outside the repository according to the applicable data-use requirements.

