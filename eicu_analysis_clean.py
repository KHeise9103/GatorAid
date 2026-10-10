"""GatorAId — eICU clean analysis wrapper.

The current project record documents that the eICU analysis used the same
Below / Within / Above x Repletion / No Repletion framework as Mayo, but eICU
replacement exposure was based on medication order/start timing rather than
Mayo MAR-confirmed administration.

The original raw eICU extraction script is not present in the project files
used to build this clean bundle.  Therefore this file does not invent a new raw
SQL query.  Instead it provides a transparent adapter for standardized eICU
lab and medication exports and runs the exact shared primary analysis.

Required standardized inputs
----------------------------
Labs:
    patientunitstayid, lab_time, lab_value
Meds:
    patientunitstayid, med_time, medication_name
Optional normalized dose columns:
    normalized_first_dose, normalized_dose_unit

If your eICU export uses offsets (for example labresultoffset or drugstartoffset),
convert them to a common numeric timeline first; see standardize_eicu_offsets().
"""

from __future__ import annotations

from typing import Dict
import pandas as pd

from gatoraid_core import run_event_analysis


def standardize_eicu_offsets(
    labs: pd.DataFrame,
    meds: pd.DataFrame,
    *,
    stay_col: str = "patientunitstayid",
    lab_offset_col: str = "labresultoffset",
    med_offset_col: str = "drugstartoffset",
    lab_value_col: str = "labresult",
    med_name_col: str = "drugname",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert eICU minute offsets to a shared synthetic datetime timeline.

    A synthetic origin is sufficient because the primary GatorAId linkage uses
    within-stay elapsed time rather than calendar date.
    """
    l = labs.copy()
    m = meds.copy()
    for c in [stay_col, lab_offset_col, lab_value_col]:
        if c not in l.columns:
            raise KeyError(f"Missing eICU lab column: {c}")
    for c in [stay_col, med_offset_col, med_name_col]:
        if c not in m.columns:
            raise KeyError(f"Missing eICU medication column: {c}")

    origin = pd.Timestamp("2000-01-01")
    l = l.rename(columns={stay_col: "encounter_id", lab_value_col: "lab_value"})
    m = m.rename(columns={stay_col: "encounter_id", med_name_col: "medication_name"})
    l["lab_time"] = origin + pd.to_timedelta(pd.to_numeric(l[lab_offset_col], errors="coerce"), unit="m")
    m["med_time"] = origin + pd.to_timedelta(pd.to_numeric(m[med_offset_col], errors="coerce"), unit="m")
    return l, m


def prepare_eicu_frames(
    labs: pd.DataFrame,
    meds: pd.DataFrame,
    *,
    stay_col: str = "patientunitstayid",
    lab_time_col: str = "lab_time",
    lab_value_col: str = "lab_value",
    med_time_col: str = "med_time",
    med_name_col: str = "medication_name",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    l = labs.rename(columns={stay_col: "encounter_id", lab_time_col: "lab_time", lab_value_col: "lab_value"}).copy()
    m = meds.rename(columns={stay_col: "encounter_id", med_time_col: "med_time", med_name_col: "medication_name"}).copy()
    if "normalized_first_dose" not in m.columns:
        m["normalized_first_dose"] = pd.NA
    if "normalized_dose_unit" not in m.columns:
        m["normalized_dose_unit"] = pd.NA
    return l, m


def run_eicu(labs: pd.DataFrame, meds: pd.DataFrame, electrolyte: str) -> Dict[str, pd.DataFrame]:
    labs_std, meds_std = prepare_eicu_frames(labs, meds)
    return run_event_analysis(labs_std, meds_std, electrolyte)


if __name__ == "__main__":
    print(__doc__)
    print("Import run_eicu() after loading the final eligible eICU lab and medication exports.")
