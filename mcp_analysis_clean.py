"""GatorAId — Mayo Clinical Platform (MCP) clean analysis wrapper.

Source basis
------------
This file consolidates the analysis rules used in the three MCP notebooks for
potassium, magnesium, and phosphorus.  It includes the corrected electrolyte-
specific exclusions established during manuscript QA.

Important status note
---------------------
Potassium and magnesium primary manuscript outputs have been QA'd against the
current MCP reruns.  The phosphorus cohort exclusions are corrected here, but
downstream phosphorus manuscript values should be regenerated before being
reported as final.

This script is an analysis layer: it accepts standardized lab and medication
frames and does not embed Mayo credentials or PHI table names.  Secure MCP
extraction remains in the credentialed notebooks/environment.
"""

from __future__ import annotations

from typing import Dict, Tuple
import numpy as np
import pandas as pd

from gatoraid_core import run_event_analysis, export_results

STUDY_START = "2019-01-01"
STUDY_END = "2026-01-01"


def _contains(s: pd.Series, pattern: str, *, case: bool = True) -> pd.Series:
    return s.fillna("").astype(str).str.contains(pattern, regex=True, case=case, na=False)


def mcp_exclusion_flags(problem_df: pd.DataFrame, electrolyte: str) -> pd.DataFrame:
    """Return encounter-level general and electrolyte-specific exclusion flags.

    Expected columns:
      PATIENT_ENCOUNTER_ID, DIAGNOSIS_NAME, ICD9_CODE_LIST, ICD10_CODE_LIST
    """
    required = ["PATIENT_ENCOUNTER_ID", "DIAGNOSIS_NAME", "ICD9_CODE_LIST", "ICD10_CODE_LIST"]
    missing = [c for c in required if c not in problem_df.columns]
    if missing:
        raise KeyError(f"Missing problem-list columns: {missing}")

    x = problem_df.copy()
    enc = x["PATIENT_ENCOUNTER_ID"].astype(str)
    dx = x["DIAGNOSIS_NAME"].fillna("").astype(str)
    icd9 = x["ICD9_CODE_LIST"].fillna("").astype(str)
    icd10 = x["ICD10_CODE_LIST"].fillna("").astype(str)

    general = (
        _contains(icd10, r"\bZ00\.129\b")
        | _contains(icd10, r"\bZ33\.1\b")
        | _contains(icd10, r"\bT20|\bT21|\bT22|\bT23|\bT24|\bT25")
        | _contains(icd10, r"\bN18\.6\b|\bN18\.30\b|\bN18\.4\b|\bN18\.5\b|\bN17\.9\b")
        | _contains(icd10, r"\bZ49\.02\b")
        | _contains(icd9, r"\b99\.2\b|\b39\.95\b")
        | _contains(icd10, r"\bS14\.109A\b|\bS24\.109A\b|\bS34\.109A\b")
    )

    key = electrolyte.lower().strip()
    if key in {"k", "potassium"}:
        specific = (
            _contains(icd10, r"\bI15\.1\b")
            | (_contains(icd10, r"\bI15\.1\b") & _contains(icd10, r"\bE87\.5\b"))
            | _contains(icd10, r"\bN15\.8\b")
            | _contains(icd10, r"\bE26\.81\b")
            | _contains(icd10, r"\bE25\.0\b")
            | _contains(icd10, r"\bE72\.09\b")
            | _contains(icd10, r"\bD75\.83\b")
            | _contains(icd10, r"\bD72\.82\b")
            | _contains(icd10, r"\bN25\.89\b")
        )
    elif key in {"mg", "magnesium"}:
        specific = _contains(icd10, r"\bE26\.81\b|\bN15\.8\b")
    elif key in {"p", "pi", "phos", "phosphate", "phosphorus"}:
        # Fanconi syndrome is E72.09. Do not exclude Fanconi anemia (D61.03).
        fanconi = _contains(icd10, r"\bE72\.09\b")
        rickets = dx.str.contains(r"\bricket", case=False, regex=True, na=False)
        vitamin_d = (
            dx.str.contains(r"vitamin\s*d|vit\s*d|hypovitaminosis\s*d", case=False, regex=True, na=False)
            | _contains(icd10, r"\bE55\.9\b")
        )
        bone_turnover = dx.str.contains(
            r"bone.*turnover|turnover.*bone", case=False, regex=True, na=False
        )
        specific = fanconi | rickets | vitamin_d | bone_turnover
    else:
        raise ValueError(f"Unsupported electrolyte: {electrolyte!r}")

    flags = pd.DataFrame({
        "encounter_id": enc,
        "general_exclusion": general,
        "specific_exclusion": specific,
    })
    out = flags.groupby("encounter_id", as_index=False).max()
    out["any_exclusion"] = out["general_exclusion"] | out["specific_exclusion"]
    return out


def mcp_exclusion_detail(problem_df: pd.DataFrame, electrolyte: str) -> pd.DataFrame:
    """Criterion-level encounter counts for manuscript/STROBE QA."""
    x = problem_df.copy()
    dx = x["DIAGNOSIS_NAME"].fillna("").astype(str)
    icd10 = x["ICD10_CODE_LIST"].fillna("").astype(str)
    enc = x["PATIENT_ENCOUNTER_ID"].astype(str)

    def n(mask):
        return enc[mask].nunique()

    key = electrolyte.lower().strip()
    if key in {"k", "potassium"}:
        criteria = {
            "Gordon syndrome / I15.1": _contains(icd10, r"\bI15\.1\b"),
            "Hyperkalemic hypertension / I15.1 + E87.5": _contains(icd10, r"\bI15\.1\b") & _contains(icd10, r"\bE87\.5\b"),
            "Liddle syndrome / I15.1": _contains(icd10, r"\bI15\.1\b"),
            "Gitelman / N15.8": _contains(icd10, r"\bN15\.8\b"),
            "Bartter / E26.81": _contains(icd10, r"\bE26\.81\b"),
            "21-hydroxylase deficiency / E25.0": _contains(icd10, r"\bE25\.0\b"),
            "Fanconi syndrome / E72.09": _contains(icd10, r"\bE72\.09\b"),
            "Platelet criterion / D75.83": _contains(icd10, r"\bD75\.83\b"),
            "High WBC criterion / D72.82": _contains(icd10, r"\bD72\.82\b"),
            "RTA type 4 / N25.89": _contains(icd10, r"\bN25\.89\b"),
        }
    elif key in {"mg", "magnesium"}:
        criteria = {
            "Bartter / E26.81": _contains(icd10, r"\bE26\.81\b"),
            "Gitelman / N15.8": _contains(icd10, r"\bN15\.8\b"),
        }
    else:
        criteria = {
            "Fanconi syndrome / E72.09": _contains(icd10, r"\bE72\.09\b"),
            "Rickets": dx.str.contains(r"\bricket", case=False, regex=True, na=False),
            "Vitamin D deficiency": (
                dx.str.contains(r"vitamin\s*d|vit\s*d|hypovitaminosis\s*d", case=False, regex=True, na=False)
                | _contains(icd10, r"\bE55\.9\b")
            ),
            "Bone turnover disease": dx.str.contains(r"bone.*turnover|turnover.*bone", case=False, regex=True, na=False),
        }
    return pd.DataFrame({"criterion": list(criteria), "n_encounters": [n(m) for m in criteria.values()]})


def normalize_mcp_potassium_dose(meds: pd.DataFrame) -> pd.DataFrame:
    """Normalize explicit potassium first-dose values to mEq when already documented as mEq."""
    x = meds.copy()
    dose = pd.to_numeric(x.get("dose_numeric", x.get("DOSE_DESCRIPTION")), errors="coerce")
    unit = x.get("DOSE_UNIT", pd.Series(index=x.index, dtype=object)).fillna("").astype(str).str.lower()
    x["normalized_first_dose"] = np.where(unit.str.contains("meq"), dose, np.nan)
    x["normalized_dose_unit"] = "mEq K"
    return x


def normalize_mcp_magnesium_dose(meds: pd.DataFrame) -> pd.DataFrame:
    """Normalize MCP first magnesium dose to mEq elemental magnesium.

    Uses 1 g magnesium sulfate = 8.12 mEq Mg and elemental magnesium
    equivalent weight 12.1525 mg/mEq. Continuous-rate / concentration-only
    records remain unresolved rather than being imputed.
    """
    x = meds.copy()
    name = x.get("MEDICATION_NAME", x.get("medication_name", "")).fillna("").astype(str).str.upper()
    unit = x.get("DOSE_UNIT", pd.Series(index=x.index, dtype=object)).fillna("").astype(str).str.lower().str.strip()
    dose = pd.to_numeric(x.get("dose_numeric", x.get("DOSE_DESCRIPTION")), errors="coerce")
    meq = pd.Series(np.nan, index=x.index, dtype="float64")

    sulfate = name.str.contains("MAGNESIUM SULFATE", na=False)
    meq.loc[sulfate & unit.eq("g") & dose.gt(0)] = dose[sulfate & unit.eq("g") & dose.gt(0)] * 8.12
    meq.loc[sulfate & unit.eq("mg") & dose.gt(0)] = dose[sulfate & unit.eq("mg") & dose.gt(0)] / 1000.0 * 8.12

    # Fixed-dose MgSO4 product names when administered quantity is absent.
    fixed = {
        r"1 GRAM/100 ML": 1.0,
        r"2 GRAM/50 ML": 2.0,
        r"4 GRAM/50 ML": 4.0,
        r"4 GRAM/100 ML": 4.0,
        r"IVPB 1 G": 1.0,
        r"IVPB 2 G": 2.0,
        r"IVPB 3 G": 3.0,
        r"IVPB 4 G": 4.0,
    }
    for pattern, grams in fixed.items():
        mask = sulfate & name.str.contains(pattern, regex=True, na=False) & meq.isna() & (~dose.gt(0))
        meq.loc[mask] = grams * 8.12

    # Elemental magnesium-containing oral products.
    eqwt = 12.1525  # mg elemental Mg per mEq
    # 241.3 mg elemental Mg as magnesium oxide 400 mg.
    mask = name.str.contains(r"241\.3 MG.*MAGNESIUM OXIDE 400 MG", regex=True, na=False)
    meq.loc[mask & meq.isna()] = 241.3 / eqwt

    # Products explicitly named as elemental magnesium strengths.
    for pattern, mg_elemental in [
        (r"MAGNESIUM 64 MG.*CHLORIDE", 64.0),
        (r"MAGNESIUM 71\.5 MG.*CHLORIDE", 71.5),
        (r"54 MG ELEMENTAL MAGNESIUM", 54.0),
        (r"MAGNESIUM 100 MG.*CITRATE", 100.0),
    ]:
        mask = name.str.contains(pattern, regex=True, na=False)
        meq.loc[mask & meq.isna()] = mg_elemental / eqwt

    # Magnesium gluconate 500 mg product documented as 27 mg elemental Mg.
    mask = name.str.contains(r"MAGNESIUM GLUCONATE 500 MG", regex=True, na=False)
    meq.loc[mask & meq.isna()] = 27.0 / eqwt

    # Magnesium oxide products: use 241.3/400 elemental fraction when the
    # documented dose is MgO mass rather than elemental Mg.
    oxide = name.str.contains("MAGNESIUM", na=False) & name.str.contains("OXIDE", na=False)
    oxide_fraction = 241.3 / 400.0
    oxide_mass = oxide & unit.eq("mg") & dose.gt(0) & meq.isna()
    meq.loc[oxide_mass] = (dose[oxide_mass] * oxide_fraction) / eqwt

    # Explicit unit 'mg of magnesium' means elemental Mg.
    elemental_unit = unit.str.contains("mg of magnesium", na=False) & dose.gt(0)
    meq.loc[elemental_unit] = dose[elemental_unit] / eqwt

    x["normalized_first_dose"] = meq
    x["normalized_dose_unit"] = "mEq Mg"
    return x


def normalize_mcp_phosphorus_dose(meds: pd.DataFrame) -> pd.DataFrame:
    """Normalize only explicit mmol phosphorus doses.

    The final MCP phosphorus first-dose mapping has not yet been QA-finalized.
    To avoid silently inventing conversions, only explicit mmol doses are used;
    all other units remain unresolved until the final phosphate dose audit.
    """
    x = meds.copy()
    dose = pd.to_numeric(x.get("dose_numeric", x.get("DOSE_DESCRIPTION")), errors="coerce")
    unit = x.get("DOSE_UNIT", pd.Series(index=x.index, dtype=object)).fillna("").astype(str).str.lower()
    x["normalized_first_dose"] = np.where(unit.str.contains("mmol"), dose, np.nan)
    x["normalized_dose_unit"] = "mmol P"
    return x


def prepare_mcp_frames(labs: pd.DataFrame, meds: pd.DataFrame, electrolyte: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Map MCP notebook column names to the standardized core schema."""
    key = electrolyte.lower().strip()
    value_candidates = {
        "potassium": ["k_value", "lab_value"],
        "magnesium": ["mg_value", "lab_value"],
        "phosphorus": ["phos_value", "lab_value"],
    }
    if key in {"k"}: key = "potassium"
    if key in {"mg"}: key = "magnesium"
    if key in {"p", "pi", "phos", "phosphate"}: key = "phosphorus"

    value_col = next((c for c in value_candidates[key] if c in labs.columns), None)
    if value_col is None:
        raise KeyError(f"Could not find electrolyte value column for {key}")

    l = labs.rename(columns={
        "PATIENT_ENCOUNTER_ID": "encounter_id",
        "lab_time": "lab_time",
        value_col: "lab_value",
    }).copy()

    m = meds.rename(columns={
        "PATIENT_ENCOUNTER_ID": "encounter_id",
        "med_time": "med_time",
        "MEDICATION_NAME": "medication_name",
    }).copy()

    if key == "potassium":
        m = normalize_mcp_potassium_dose(m)
    elif key == "magnesium":
        m = normalize_mcp_magnesium_dose(m)
    else:
        m = normalize_mcp_phosphorus_dose(m)
    return l, m


def run_mcp(labs: pd.DataFrame, meds: pd.DataFrame, electrolyte: str) -> Dict[str, pd.DataFrame]:
    labs_std, meds_std = prepare_mcp_frames(labs, meds, electrolyte)
    return run_event_analysis(labs_std, meds_std, electrolyte)


if __name__ == "__main__":
    print(__doc__)
    print("Import run_mcp(), mcp_exclusion_flags(), and mcp_exclusion_detail() in the secure MCP environment.")
