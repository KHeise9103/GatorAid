"""GatorAId — MIMIC-IV clean analysis wrapper.

This file keeps the source-backed MIMIC-IV BigQuery extraction pattern used in
the project and routes standardized lab/replacement events through the common
GatorAId primary analysis engine.

The medication dictionaries must still be reviewed before a final manuscript
run: candidate text matches can include non-repletion products.  The analysis
engine itself does not treat every potassium/magnesium/phosphate text match as
replacement unless it is present in the supplied final medication frame.
"""

from __future__ import annotations

from typing import Dict
import pandas as pd

from gatoraid_core import run_event_analysis


def bigquery_lab_sql() -> str:
    return r"""
WITH lab_ids AS (
  SELECT itemid, label
  FROM `physionet-data.mimiciv_hosp.d_labitems`
  WHERE REGEXP_CONTAINS(LOWER(label), r'(potassium|magnesium|phosphorus|phosphate)')
), icu AS (
  SELECT subject_id, hadm_id, stay_id, intime, outtime
  FROM `physionet-data.mimiciv_icu.icustays`
)
SELECT
  l.subject_id,
  l.hadm_id,
  i.stay_id,
  l.charttime,
  l.itemid,
  d.label,
  l.valuenum,
  l.valueuom,
  l.ref_range_lower,
  l.ref_range_upper,
  l.flag
FROM `physionet-data.mimiciv_hosp.labevents` l
JOIN lab_ids d USING (itemid)
JOIN icu i
  ON l.subject_id = i.subject_id
 AND l.hadm_id = i.hadm_id
 AND l.charttime BETWEEN i.intime AND i.outtime
WHERE l.valuenum IS NOT NULL
ORDER BY l.subject_id, l.hadm_id, i.stay_id, l.charttime
"""


def bigquery_replacement_candidate_sql() -> str:
    return r"""
WITH replacement_ids AS (
  SELECT itemid, label
  FROM `physionet-data.mimiciv_icu.d_items`
  WHERE REGEXP_CONTAINS(
    LOWER(label),
    r'(potassium|magnesium|phosphate|phosphorus|kcl|mgso4)'
  )
)
SELECT
  ie.subject_id,
  ie.hadm_id,
  ie.stay_id,
  ie.starttime,
  ie.endtime,
  ie.itemid,
  d.label,
  ie.amount,
  ie.amountuom,
  ie.rate,
  ie.rateuom,
  ie.ordercategoryname,
  ie.statusdescription
FROM `physionet-data.mimiciv_icu.inputevents` ie
JOIN replacement_ids d USING (itemid)
ORDER BY ie.subject_id, ie.hadm_id, ie.stay_id, ie.starttime
"""


def extract_bigquery(client) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the project BigQuery extraction using an authenticated client."""
    labs = client.query(bigquery_lab_sql()).to_dataframe()
    meds = client.query(bigquery_replacement_candidate_sql()).to_dataframe()
    return labs, meds


def prepare_mimic_frames(
    labs: pd.DataFrame,
    meds: pd.DataFrame,
    *,
    lab_value_col: str = "valuenum",
    med_name_col: str = "label",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Map a finalized MIMIC ICU-stay event export to the shared schema.

    The final medication frame passed here should already contain only true
    replacement events after medication-dictionary QA.
    """
    l = labs.rename(columns={
        "stay_id": "encounter_id",
        "charttime": "lab_time",
        lab_value_col: "lab_value",
    }).copy()
    m = meds.rename(columns={
        "stay_id": "encounter_id",
        "starttime": "med_time",
        med_name_col: "medication_name",
    }).copy()

    if "normalized_first_dose" not in m.columns:
        # Preserve explicit amount only when the caller has already established
        # that amountuom is the desired normalized unit.  Otherwise unresolved.
        m["normalized_first_dose"] = pd.NA
    if "normalized_dose_unit" not in m.columns:
        m["normalized_dose_unit"] = pd.NA
    return l, m


def run_mimic_iv(labs: pd.DataFrame, meds: pd.DataFrame, electrolyte: str) -> Dict[str, pd.DataFrame]:
    labs_std, meds_std = prepare_mimic_frames(labs, meds)
    return run_event_analysis(labs_std, meds_std, electrolyte)


if __name__ == "__main__":
    print(__doc__)
    print("Use extract_bigquery(client) or pass finalized MIMIC lab/medication frames to run_mimic_iv().")
