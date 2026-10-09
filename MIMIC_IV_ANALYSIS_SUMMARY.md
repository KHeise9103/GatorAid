# MIMIC-IV Electrolyte Replacement Analysis

## Project purpose

This analysis reproduces the GatorAId electrolyte-replacement framework in **MIMIC-IV v3.1** for potassium, magnesium, and phosphorus. The core design mirrors the Mayo/eICU reporting structure:

**Below / Within / Above reference range × Repletion / No repletion**, with a 24-hour clinical decision window and a 24-hour follow-up window.

The primary analysis uses **actual eMAR administrations** rather than medication orders whenever possible. Raw MIMIC data are credentialed and are **not included in this repository**.

---

## Software and implementation

- Database engine: DuckDB 1.5.6
- Analysis environment: Jupyter / Python
- Source dataset: MIMIC-IV v3.1
- Primary tables used: `icustays`, `labevents`, `emar`, `diagnoses_icd`, `procedureevents`, and supporting dictionaries
- Primary electrolyte laboratory sources:
  - Potassium: chemistry potassium, itemid **50971**
  - Magnesium: itemid **50960**
  - Phosphate/phosphorus: itemid **50970**
- Primary medication exposure source: eMAR administration records
- Accepted administration statuses:
  - `Administered`
  - `Partial Administered`
  - `Delayed Administered`
  - `Administered in Other Location`

### Key derived files

- `gatoraid_electrolyte_eligibility_final.parquet`
- `gatoraid_analysis_cohort_final.parquet`
- `gatoraid_final_analysis_labs.parquet`
- `gatoraid_final_analysis_repletion.parquet`
- `gatoraid_final_24h_workflow_classified.parquet`
- `gatoraid_final_lab_workflow_dedup.parquet`
- `gatoraid_episode_followup_final_raw.parquet`
- `MIMIC_IV_GatorAId_final_primary_results.csv`

---

# 1. Cohort construction

## Sequential general STROBE flow

| Step | ICU stays at step | Excluded | Remaining |
|---|---:|---:|---:|
| Initial ICU stays | 94,458 | 0 | 94,458 |
| ICU stay ≤48 h | 94,458 | 48,107 | 46,351 |
| Missing ICU LOS | 46,351 | 14 | 46,337 |
| Age <18 | 46,337 | 0 | 46,337 |
| Pregnant | 46,337 | 2 | 46,335 |
| Burn, 2nd/3rd degree | 46,335 | 21 | 46,314 |
| CKD stage 3–5 | 46,314 | 3,670 | 42,644 |
| ESRD | 42,644 | 2,542 | 40,102 |
| AKI with implemented RRT | 40,102 | 1,553 | 38,549 |
| Acute spinal cord injury/transection | 38,549 | 7 | **38,542** |

Of the **38,542** stays remaining after LOS and general exclusions:

- **23,099** had no qualifying electrolyte repletion
- **15,443** had at least one qualifying electrolyte repletion and entered electrolyte-specific branches

## Electrolyte-specific exclusions

### Potassium

Starting from 15,443 stays:

| Exclusion | n excluded | Remaining |
|---|---:|---:|
| Gordon/Liddle syndrome | 0 | 15,443 |
| Gitelman syndrome | 0 | 15,443 |
| Bartter syndrome | 4 | 15,439 |
| 21-hydroxylase deficiency | 1 | 15,438 |
| Fanconi syndrome | 2 | 15,436 |
| Platelets >1,000,000 cells/µL | 26 | 15,410 |
| WBC >120,000 cells/µL | 17 | 15,393 |
| RTA type 4 | 11 | **15,382** |

### Magnesium

| Exclusion | n excluded | Remaining |
|---|---:|---:|
| Bartter syndrome | 4 | 15,439 |
| Gitelman syndrome | 0 | **15,439** |

### Phosphorus

| Exclusion | n excluded | Remaining |
|---|---:|---:|
| Fanconi syndrome | 2 | 15,441 |
| Rickets | 0 | 15,441 |
| Vitamin D deficiency | 398 | 15,043 |
| Bone turnover disease | 3 | **15,040** |

---

# 2. Primary laboratory definitions

## Record-level manuscript filters

- Potassium: **2.5 < K < 7.0**
- Magnesium: **1.0 < Mg < 4.0**
- Phosphorus: **0.5 < P < 4.0**

## Clinical categories used in the 3×2 analysis

### Potassium

- Below: **<3.51**
- Within: **3.51–4.59**
- Source ambiguity: **4.60–5.00**
- Above: **>5.00**

The manuscript source ranges leave a 4.60–5.00 gap. The primary analysis does **not** force these records into either Within or Above; they are retained as `SOURCE_AMBIGUITY`.

### Magnesium

- Below: **<1.4 mg/dL**
- Within: **1.4–2.2 mg/dL**
- Above: **>2.2 mg/dL**

### Phosphorus

- Below: **<3.0 mg/dL**
- Within: **3.0–4.5 mg/dL**
- Above: **>4.5 mg/dL**

Important limitation: the primary manuscript record filter is **P <4.0**, so the primary MIMIC analysis cannot estimate baseline Above >4.5 phosphorus behavior. Above-range phosphorus is therefore reported as **NR / not estimable** in the primary report.

---

# 3. Laboratory filtering and exact-duplicate handling

| Electrolyte | Final analysis stays before lab filter | Primary lab records before filter | Records passing filter | Records removed | Stays with qualifying records |
|---|---:|---:|---:|---:|---:|
| Potassium | 15,382 | 180,565 | 179,905 | 660 | 15,379 |
| Magnesium | 15,439 | 172,220 | 171,674 | 546 | 15,434 |
| Phosphorus | 15,040 | 161,240 | 115,155 | 46,085 | 14,489 |

Exact duplicate clinical lab observations were removed before final decision-level analysis:

- Potassium: **13** duplicate records removed → **179,892** final decision records
- Magnesium: **17** removed → **171,657**
- Phosphorus: **4** removed → **115,151**

---

# 4. Medication/repletion definitions

## Primary medication classes

### Potassium

Primary potassium replacement was based on administered potassium chloride products, including oral/powder and replacement formulations.

### Magnesium

Primary magnesium replacement was based on administered **magnesium sulfate**.

### Phosphorus

Primary phosphorus replacement included phosphorus, sodium phosphate, potassium phosphate, and phosphate replacement products.

## Why eMAR was used

A medication-order/POE audit showed that order identifiers can represent long-lived or reusable protocols, especially for potassium and magnesium. Collapsing by POE order would therefore merge distinct clinical administrations. The primary analysis uses **actual eMAR administrations**, then collapses them into clinical replacement episodes based on the closest preceding lab.

---

# 5. 24-hour event-matching logic

For each administered replacement:

1. Find the closest preceding qualifying electrolyte lab within 24 hours.
2. Treat the unique preceding lab as the clinical decision episode.
3. When multiple labs precede one replacement, the closest lab is selected and earlier labs are classified as `ML→1R` ignored observations.
4. For response analysis, find the first subsequent primary electrolyte lab within 24 hours.
5. Follow-up uses the unrestricted primary lab stream so post-treatment values are not truncated by the baseline STROBE filter.

## Single-lab versus multiple-lab replacement episodes

| Electrolyte | Replacement episodes | 1L→1R | ML→1R |
|---|---:|---:|---:|
| Potassium | 39,760 | 26.85% | **73.15%** |
| Magnesium | 26,530 | 39.34% | **60.66%** |
| Phosphorus | 10,138 | 37.70% | **62.30%** |

Multiple preceding labs were the dominant pattern for all three electrolytes.

---

# 6. Overall MIMIC-IV analysis populations

| Measure | Potassium | Magnesium | Phosphorus |
|---|---:|---:|---:|
| Final eligible ICU stays | 15,382 | 15,439 | 15,040 |
| Stays with qualifying electrolyte labs | 15,379 | 15,434 | 14,489 |
| Final decision-level lab events | 179,892 | 171,657 | 115,151 |
| Replacement-linked events | 39,760 | 26,530 | 10,138 |
| Overall replacement rate | **22.1%** | **15.5%** | **8.8%** |

---

# 7. Primary 3×2 replacement behavior

## Potassium

| Baseline group | Repletion | n | % within group |
|---|---|---:|---:|
| Below <3.51 | Repletion | 12,908 | **53.4%** |
| Below <3.51 | No repletion | 11,266 | 46.6% |
| Within 3.51–4.59 | Repletion | 26,171 | **25.1%** |
| Within 3.51–4.59 | No repletion | 78,112 | 74.9% |
| Above >5.0 | Repletion | 112 | **1.3%** |
| Above >5.0 | No repletion | 8,556 | 98.7% |

Additional source-gap records, 4.60–5.00: **20,309 (11.3% of K decision records)**.

## Magnesium

| Baseline group | Repletion | n | % within group |
|---|---|---:|---:|
| Below <1.4 | Repletion | 644 | **75.4%** |
| Below <1.4 | No repletion | 210 | 24.6% |
| Within 1.4–2.2 | Repletion | 25,017 | **24.5%** |
| Within 1.4–2.2 | No repletion | 77,108 | 75.5% |
| Above >2.2 | Repletion | 869 | **1.8%** |
| Above >2.2 | No repletion | 47,403 | 98.2% |

## Phosphorus

| Baseline group | Repletion | n | % within group |
|---|---|---:|---:|
| Below <3.0 | Repletion | 9,542 | **19.2%** |
| Below <3.0 | No repletion | 40,194 | **80.8%** |
| Within 3.0–4.5* | Repletion | 596 | **1.0%** |
| Within 3.0–4.5* | No repletion | 57,774 | 99.0% |
| Above >4.5 | NR | NR | NR |

\*Because the primary baseline filter is P <4.0, observed Within-range phosphorus values are restricted to the portion below 4.0.

---

# 8. Follow-up outcomes within 24 hours

## Potassium

| Baseline group | Repletion | Follow-up available | Remained low | Normal | K 4.60–5.00 gap | High |
|---|---|---:|---:|---:|---:|---:|
| Below <3.51 | Repletion | 90.4% | 32.6% | **60.9%** | 4.5% | 2.0% |
| Below <3.51 | No repletion | 74.3% | 36.0% | **56.7%** | 5.2% | 2.1% |
| Within 3.51–4.59 | Repletion | 89.3% | 14.0% | 77.2% | 7.0% | 1.8% |
| Within 3.51–4.59 | No repletion | 73.4% | 10.2% | 79.4% | 8.7% | 1.7% |
| Above >5.0 | Repletion | 89.3% | 8.0% | 42.0% | 34.0% | 16.0% |
| Above >5.0 | No repletion | 89.3% | 1.1% | 18.9% | 30.0% | 50.1% |

## Magnesium

| Baseline group | Repletion | Follow-up available | Normal | Remained low | High |
|---|---|---:|---:|---:|---:|
| Below <1.4 | Repletion | 95.3% | **55.2%** | 2.9% | **41.9%** |
| Below <1.4 | No repletion | 85.2% | 52.5% | 5.0% | 42.5% |
| Within 1.4–2.2 | Repletion | 87.7% | 57.9% | 0.2% | 41.8% |
| Within 1.4–2.2 | No repletion | 71.0% | 90.8% | 0.3% | 9.0% |
| Above >2.2 | Repletion | 94.2% | 11.6% | 0.0% | 88.4% |
| Above >2.2 | No repletion | 81.8% | 29.2% | 0.1% | 70.7% |

## Phosphorus

| Baseline group | Repletion | Follow-up available | Normal | Remained low | High |
|---|---|---:|---:|---:|---:|
| Below <3.0 | Repletion | 89.4% | **33.4%** | **63.5%** | 3.1% |
| Below <3.0 | No repletion | 71.9% | 30.0% | 69.0% | 0.9% |
| Within 3.0–4.5 | Repletion | 87.4% | 66.2% | 23.2% | 10.6% |
| Within 3.0–4.5 | No repletion | 73.3% | 72.0% | 23.9% | 4.1% |

---

# 9. Paired change and direct repletion effects

## Potassium

### Paired within-group change

| Baseline group | Repletion | n paired | Mean ΔK | Paired Cohen's d |
|---|---|---:|---:|---:|
| Below <3.51 | Repletion | 11,672 | +0.492 | 0.897 |
| Below <3.51 | No repletion | 8,368 | +0.477 | 0.808 |
| Within 3.51–4.59 | Repletion | 23,374 | +0.169 | 0.377 |
| Within 3.51–4.59 | No repletion | 57,310 | −0.024 | −0.058 |
| Above >5.0 | Repletion | 100 | −0.883 | −1.047 |
| Above >5.0 | No repletion | 7,639 | −0.428 | −0.636 |

### Direct effect of repletion on ΔK

| Baseline group | Repletion ΔK | No-repletion ΔK | Difference | Cohen's d | Effect |
|---|---:|---:|---:|---:|---|
| Below <3.51 | +0.492 | +0.477 | **+0.015** | 0.026 | Negligible |
| Within 3.51–4.59 | +0.169 | −0.024 | **+0.193** | 0.451 | Small |
| Above >5.0 | −0.883 | −0.428 | **−0.455** | −0.674 | Moderate |

Six-group ANOVA on ΔK:

- Paired observations: **108,463**
- F = **5406.44**
- p < **0.001**
- η² = **0.200**
- Overall effect: **Large**

## Magnesium

### Paired within-group change

| Baseline group | Repletion | n paired | Mean ΔMg | Paired Cohen's d |
|---|---|---:|---:|---:|
| Below <1.4 | Repletion | 614 | +0.976 | 0.854 |
| Below <1.4 | No repletion | 179 | +0.907 | 1.842 |
| Within 1.4–2.2 | Repletion | 21,930 | +0.483 | 0.514 |
| Within 1.4–2.2 | No repletion | 54,730 | −0.003 | −0.014 |
| Above >2.2 | Repletion | 819 | +0.221 | 0.322 |
| Above >2.2 | No repletion | 38,760 | −0.132 | −0.540 |

### Direct effect of repletion on ΔMg

| Baseline group | Repletion ΔMg | No-repletion ΔMg | Difference | Cohen's d | Holm p | Effect |
|---|---:|---:|---:|---:|---:|---|
| Below <1.4 | +0.976 | +0.907 | **+0.068** | 0.066 | 0.247 | Negligible |
| Within 1.4–2.2 | +0.483 | −0.003 | **+0.487** | 0.897 | <0.001 | Large |
| Above >2.2 | +0.221 | −0.132 | **+0.354** | 1.351 | <0.001 | Large |

Six-group ANOVA on ΔMg:

- Paired observations: **117,032**
- F = **5595.18**
- p < **0.001**
- η² = **0.193**
- Overall effect: **Large**

## Phosphorus

### Direct effect of repletion on ΔP

| Baseline group | Repletion ΔP | No-repletion ΔP | Difference | Cohen's d | Holm p | Effect |
|---|---:|---:|---:|---:|---:|---|
| Below <3.0 | +0.841 | +0.214 | **+0.627** | 0.728 | <0.001 | Moderate |
| Within 3.0–4.5 | +0.186 | −0.034 | **+0.221** | 0.344 | <0.001 | Small |
| Above >4.5 | NR | NR | NR | NR | NR | NR |

Observed-group ANOVA on ΔP:

- Paired observations: **80,320**
- Observed groups compared: **4**
- F = **3320.09**
- p < **0.001**
- η² = **0.110**
- Overall effect: **Moderate**

### Phosphorus normalization comparison

| Baseline group | Repletion normalization | No-repletion normalization | Absolute difference | p-value | Cohen's h | Effect |
|---|---:|---:|---:|---:|---:|---|
| Below <3.0 | 33.4% | 30.0% | **+3.3 pp** | <0.001 | 0.072 | Negligible |
| Within 3.0–4.5 | 66.2% | 72.0% | **−5.8 pp** | 0.003 | −0.125 | Negligible |

---

# 10. Potassium hemolysis audit

A raw `labevents` comment audit was performed because the manuscript specifically requested hemolysis counts.

Naive matching on the substring `hemoly` was rejected because many comments were generic advisories (for example, statements that hemolysis can falsely elevate potassium) or explicitly stated that the specimen was not hemolyzed.

The final specimen-level definition required explicit language that the **specimen/sample itself was hemolyzed**.

Results:

- Explicitly hemolyzed potassium specimens: **600**
- ICU stays with at least one explicitly hemolyzed K specimen: **475**
- Patients: **468**
- All 600 explicit-positive specimens had numeric K results

Sensitivity excluding the 600 explicit-positive hemolyzed specimens:

| K metric | Primary | Excluding hemolysis |
|---|---:|---:|
| Low K replaced | 53.40% | 53.39% |
| Within-range K replaced | 25.10% | 25.13% |
| Above-range K replaced | 1.29% | 1.29% |

Conclusion: explicit hemolysis did **not materially alter** the potassium behavioral results.

---

# 11. Medication-definition sensitivity analyses

## Magnesium oxide sensitivity

Primary Mg replacement used magnesium sulfate. Adding administered magnesium oxide produced only modest changes:

| Mg metric | Primary sulfate | Sulfate + magnesium oxide |
|---|---:|---:|
| Low Mg replaced | 75.41% | 76.51% |
| Within-range Mg replaced | 24.50% | 25.54% |
| Above-range Mg replaced | 1.80% | 2.18% |

Conclusion: magnesium findings are robust to inclusion of magnesium oxide.

## Potassium phosphate sensitivity

Primary potassium replacement treated potassium phosphate as a phosphorus replacement, not a potassium replacement. Counting potassium phosphate as a K intervention increased replacement rates but did not reverse the interpretation:

| K metric | Primary K | K + potassium phosphate |
|---|---:|---:|
| Low K replaced | 53.39% | 57.70% |
| Within-range K replaced | 25.10% | 26.69% |
| Source-gap replaced | 3.10% | 3.29% |
| Above-range K replaced | 1.29% | 1.35% |

Conclusion: potassium-phosphate assignment is a meaningful but non-interpretation-changing sensitivity.

---

# 12. Source-range sensitivity analyses

## Potassium 4.60–5.00 source gap

The primary analysis leaves **20,309** records unclassified because the manuscript defines Within as 3.51–4.59 and Above as >5.0.

Sensitivity examples:

- Treating the gap as Above: above-range replacement ≈ **2.52%**
- Treating the gap as Within: within-range replacement ≈ **21.8%**

Primary reporting retains the gap as `SOURCE_AMBIGUITY` rather than imposing an arbitrary interpretation.

## Phosphorus expanded-range sensitivity

Because the primary STROBE filter requires `0.5 < P < 4`, an expanded analysis using `0.5 < P < 20` was performed to check whether the truncated primary range distorted behavior:

- Low P replaced: **19.18%**
- Within-range P replaced: **0.88%**
- Above-range P replaced: **0.19%**

The main low-phosphorus finding remained essentially unchanged.

---

# 13. Potassium source validation

Chemistry potassium (itemid 50971) was selected as the primary potassium source. Whole-blood potassium was evaluated as a sensitivity/QA source.

Paired whole-blood vs chemistry K measurements within ±15 minutes:

- Paired observations: **40,237**
- Mean whole-blood − chemistry difference: **−0.244 mmol/L**
- Median difference: **−0.200 mmol/L**
- Mean absolute error: **0.292 mmol/L**
- Correlation: **0.875**

Chemistry potassium was retained as the primary measurement stream.

---

# 14. Timing QA

At the administration-record level, most replacement events had both a preceding lab and a subsequent lab within 24 hours.

Median timing among administrations with matches:

| Electrolyte | Median hours from prior lab to administration | Median hours from administration to follow-up lab |
|---|---:|---:|
| Potassium | 2.60 | 7.55 |
| Magnesium | 2.33 | 8.55 |
| Phosphorus | 4.60 | 7.95 |

The POE/order audit showed substantial order-to-administration lags for potassium and magnesium and confirmed that orders can function as standing/reusable protocols. This supported the decision to use actual administrations rather than order IDs as the primary intervention signal.

---

# 15. Key findings

1. **Potassium:** approximately half of low potassium values were followed by replacement (**53.4%**), while one-quarter of within-range values were also followed by replacement (**25.1%**). The direct low-K ΔK difference between repletion and no repletion was negligible (**+0.015 mmol/L; d=0.026**).

2. **Magnesium:** low magnesium was usually replaced (**75.4%**), but within-range magnesium replacement was also common (**24.5%**). Low-Mg ΔMg did not differ meaningfully by replacement status, whereas within-range and above-range replacement groups showed large treatment-associated differences.

3. **Phosphorus:** low phosphorus was replaced much less often (**19.2%**); **80.8%** of low-P decisions had no replacement. When replacement occurred, the mean improvement was substantially larger (**ΔP difference +0.627 mg/dL; d=0.728**), although the absolute normalization-rate advantage was small (**+3.3 percentage points**).

4. **Multiple testing before replacement was common**, especially for potassium (73.15% of replacement episodes had multiple qualifying labs in the preceding 24 hours).

5. **Explicit potassium hemolysis was uncommon and did not materially change the primary results.**

6. Medication-definition and source-range sensitivity analyses changed numerical estimates modestly but did **not reverse the main clinical interpretation**.

---

# 16. Important limitations

- MIMIC-IV is observational; replacement and no-replacement groups are not randomized.
- Follow-up availability differs by treatment status and baseline category, creating potential selection bias in paired analyses.
- The source manuscript contains an unresolved potassium range gap (4.60–5.00), retained explicitly as source ambiguity.
- The source manuscript's primary phosphorus record filter (`P <4`) conflicts with its nominal/above reference categories (Within through 4.5, Above >4.5). The primary MIMIC analysis preserves the manuscript filter and reports baseline Above phosphorus as NR.
- Medication assignment can be ambiguous for combination products such as potassium phosphate; sensitivity analyses were therefore performed.
- Magnesium oxide was excluded from the primary magnesium definition because it may represent scheduled/chronic oral supplementation; inclusion was assessed in sensitivity analysis.
- Statistical significance is expected to be common with the large MIMIC sample. Effect sizes and absolute differences should be emphasized over p-values alone.

---

# 17. Repository/data-use note

This repository should contain **analysis code, data dictionaries/mappings, aggregate results, and documentation only**. Do not commit raw MIMIC-IV patient-level files or derived row-level data that are restricted by the PhysioNet data-use agreement.

