# GatorAId clean Python analysis bundle

This folder is a clean, repository-ready analysis layer for the three databases currently in scope:

- **MCP / Mayo Clinical Platform** — `mcp_analysis_clean.py`
- **eICU** — `eicu_analysis_clean.py`
- **MIMIC-IV** — `mimic_iv_analysis_clean.py`
- Shared analytic functions — `gatoraid_core.py`

AMD / Amsterdam is intentionally not included.

## Why there is a shared core

The primary analysis should be identical across databases. The shared core prevents small differences in code from changing the definition of:

1. Below / Within / Above electrolyte categories.
2. Replacement linkage to the most recent preceding electrolyte lab within 24 hours.
3. Follow-up linkage to the next lab within 24 hours.
4. Category share and treatment rate.
5. Baseline Repletion vs No Repletion statistics.
6. Paired pre/post statistics.
7. First-dose median and IQR when a database-specific normalized dose is available.

## Current validation status

### MCP

- Potassium: corrected electrolyte-specific exclusions and current primary table logic are incorporated.
- Magnesium: corrected Bartter/Gitelman exclusions and corrected magnesium first-dose conversion are incorporated. Magnesium lab values remain **mg/dL**; normalized first replacement dose is **mEq Mg**.
- Phosphorus: corrected phosphate-specific exclusions are incorporated: Fanconi syndrome (E72.09), rickets, vitamin D deficiency, and bone-turnover disease. The downstream phosphorus table must still be rerun after the corrected cohort. The final phosphate first-dose conversion has **not** yet been QA-finalized, so the clean code deliberately uses only explicit mmol doses rather than inventing conversions.

### eICU

The project sources document the final eICU results and that medication exposure was based on medication order/start timing. The original raw extraction script was not available in the project files used for this bundle. Therefore the clean eICU file is an analysis adapter for a finalized eICU lab/medication export rather than a newly invented raw SQL extraction.

### MIMIC-IV

The project includes a BigQuery extraction starter. Its source-backed lab and input-event candidate queries are preserved in the clean MIMIC file. Medication candidates still require final dictionary QA before being treated as true replacement events.

## Expected standardized columns

The shared engine ultimately needs two tables.

### Labs

- `encounter_id`
- `lab_time`
- `lab_value`
- optional `lab_event_id`

### Replacement medications

- `encounter_id`
- `med_time`
- `medication_name`
- `normalized_first_dose`
- `normalized_dose_unit`

Database wrappers map native columns to this schema.

## Example

```python
from mcp_analysis_clean import run_mcp

results = run_mcp(labs_df, meds_df, "magnesium")

results["main_table"]
results["baseline_stats"]
results["paired_stats"]
```

## Dependencies

```text
pandas
numpy
scipy
```

For MIMIC-IV BigQuery extraction, also install `google-cloud-bigquery` and authenticate to PhysioNet MIMIC-IV in Google Cloud.

## Distribution figures (MCP and MIMIC-IV)

The bundle now includes publication-style histogram figures modeled on the
histogram-plus-normal-reference example discussed during manuscript QA:

- `gatoraid_distribution_figures.py` — shared plotting functions.
- `mcp_make_distribution_figures.py` — MCP figure runner.
- `mimic_make_distribution_figures.py` — MIMIC-IV figure runner.

Each electrolyte is saved as its own figure (PNG and PDF):

1. Baseline distribution with the manuscript clinical thresholds marked.
2. Follow-up minus baseline distribution for No Repletion.
3. Follow-up minus baseline distribution for Repletion.

The overlaid normal curve is descriptive only and is explicitly labeled as not
a formal normality test. Figures are generated from actual event-level values;
the scripts do not simulate histograms from means/SDs.

MCP phosphorus figures should be generated only after the corrected phosphorus
exclusion cohort has been rerun through the downstream event analysis.
