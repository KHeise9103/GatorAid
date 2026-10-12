"""GatorAId — eICU clean analysis wrapper.

This module exposes the finalized eICU event-level analysis in the same result
shape used by the shared GatorAId code and distribution-figure helpers.

Unlike the MIMIC-IV and MCP wrappers, this file intentionally does *not* rebuild
medication parsing, prior-lab matching, or follow-up matching from raw source
files.  Those steps were extensively QA'd in the eICU notebook and are stored
in finalized DuckDB tables.  The wrapper reads those frozen tables, standardizes
their columns, and calculates manuscript-ready statistics without modifying the
database.

Final event sources
-------------------
Potassium
    eicu_k_decision_followup_corrected
    Corrected categories: <3.5, 3.5-5.0, >5.0 mmol/L.

Magnesium
    eicu_decision_followup, filtered to electrolyte='magnesium'
    Categories: <1.4, 1.4-2.2, >2.2 mg/dL.

Phosphorus
    eicu_p_decision_followup_corrected
    Corrected categories: <3.0, 3.0-4.5, >4.5 mg/dL.

Only true analytic decisions are retained: Repletion and No repletion.  Ignored
episodes are excluded.  The finalized phosphorus stream is read directly so
values above 20 mg/dL that survived the source-specific QA are not silently
removed by a generic plausibility filter.

Example
-------
import duckdb
from eicu_clean_analysis import run_all_eicu

con = duckdb.connect("gatoraid_eicu.duckdb", read_only=True)
results_by_electrolyte = run_all_eicu(con)

# Compatible with the shared distribution-figure helper:
# results_by_electrolyte["potassium"]["events"]

No local database path, credentials, or PHI are embedded in this file.
"""

from __future__ import annotations

from math import asin, sqrt
from typing import Dict, Mapping

import numpy as np
import pandas as pd
from scipy import stats

from gatoraid_core import (
    GROUP_ORDER,
    add_followup_outcomes,
    cohen_d_independent,
    followup_timing_stats,
    manuscript_table,
)


# -----------------------------------------------------------------------------
# Frozen QA targets from the finalized eICU analysis.
# These are event counts, not patient/stay counts.
# -----------------------------------------------------------------------------

EXPECTED_EVENT_COUNTS: Mapping[str, int] = {
    "potassium": 300_003,
    "magnesium": 173_555,
    "phosphorus": 129_649,
}

EXPECTED_FOLLOWUP_COUNTS: Mapping[str, int] = {
    "potassium": 221_650,
    "magnesium": 101_466,
    "phosphorus": 70_173,
}

# Final analytic ICU-stay cohorts used in the manuscript worksheet.
EXPECTED_FINAL_STAY_COUNTS: Mapping[str, int] = {
    "potassium": 34_475,
    "magnesium": 34_534,
    "phosphorus": 34_534,
}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def _key(electrolyte: str) -> str:
    key = electrolyte.lower().strip()
    aliases = {
        "k": "potassium",
        "potassium": "potassium",
        "mg": "magnesium",
        "magnesium": "magnesium",
        "p": "phosphorus",
        "pi": "phosphorus",
        "phos": "phosphorus",
        "phosphate": "phosphorus",
        "phosphorus": "phosphorus",
    }
    if key not in aliases:
        raise ValueError(f"Unsupported electrolyte: {electrolyte!r}")
    return aliases[key]


def _p_display(p: float) -> str:
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def _two_prop_z(x1: int, n1: int, x0: int, n0: int) -> tuple[float, float]:
    """Two-sided pooled two-proportion z test.

    Group 1 is Repletion; group 0 is No Repletion.
    """
    if n1 == 0 or n0 == 0:
        return np.nan, np.nan

    p1 = x1 / n1
    p0 = x0 / n0
    pooled = (x1 + x0) / (n1 + n0)
    se = np.sqrt(pooled * (1.0 - pooled) * ((1.0 / n1) + (1.0 / n0)))

    if not np.isfinite(se) or se == 0:
        return np.nan, np.nan

    z = (p1 - p0) / se
    p = 2.0 * stats.norm.sf(abs(z))
    return float(z), float(p)


def _cohens_h(p1: float, p0: float) -> float:
    """Signed Cohen's h: Repletion minus No Repletion."""
    return 2.0 * asin(sqrt(p1)) - 2.0 * asin(sqrt(p0))


def _h_effect_label(h: float) -> str:
    if pd.isna(h):
        return "NA"
    a = abs(float(h))
    if a < 0.20:
        return "Negligible"
    if a < 0.50:
        return "Small"
    if a < 0.80:
        return "Medium"
    return "Large"


def _status_series(s: pd.Series) -> pd.Series:
    return s.astype(str).map({
        "Repletion": "Repletion",
        "No repletion": "No Repletion",
        "No Repletion": "No Repletion",
    })


# -----------------------------------------------------------------------------
# Finalized event extraction
# -----------------------------------------------------------------------------


def extract_eicu_events(con, electrolyte: str) -> pd.DataFrame:
    """Read finalized analytic episodes for one electrolyte from DuckDB.

    Parameters
    ----------
    con
        Open DuckDB connection containing the finalized eICU analysis tables.
    electrolyte
        potassium / magnesium / phosphorus, with common aliases accepted.

    Returns
    -------
    pandas.DataFrame
        Standardized event-level data compatible with the shared GatorAId
        manuscript and distribution-figure helpers.
    """
    key = _key(electrolyte)

    if key == "potassium":
        sql = """
        SELECT
            index_labid AS lab_event_id,
            CAST(patientunitstayid AS VARCHAR) AS encounter_id,
            lab_category,
            baseline_value AS lab_value,
            decision_status,
            followup_value,
            followup_hours_from_anchor AS hours_to_followup
        FROM eicu_k_decision_followup_corrected
        WHERE decision_status IN ('Repletion', 'No repletion')
        """

    elif key == "magnesium":
        sql = """
        SELECT
            index_labid AS lab_event_id,
            CAST(patientunitstayid AS VARCHAR) AS encounter_id,
            lab_category,
            baseline_value AS lab_value,
            decision_status,
            followup_value,
            followup_hours_from_anchor AS hours_to_followup
        FROM eicu_decision_followup
        WHERE electrolyte = 'magnesium'
          AND decision_status IN ('Repletion', 'No repletion')
        """

    else:  # phosphorus
        sql = """
        SELECT
            index_labid AS lab_event_id,
            CAST(patientunitstayid AS VARCHAR) AS encounter_id,
            lab_category,
            index_value AS lab_value,
            decision_status,
            followup_value,
            minutes_to_followup / 60.0 AS hours_to_followup
        FROM eicu_p_decision_followup_corrected
        WHERE decision_status IN ('Repletion', 'No repletion')
        """

    x = con.execute(sql).fetchdf()

    x["lab_value"] = pd.to_numeric(x["lab_value"], errors="coerce")
    x["followup_value"] = pd.to_numeric(x["followup_value"], errors="coerce")
    x["hours_to_followup"] = pd.to_numeric(x["hours_to_followup"], errors="coerce")

    x["repletion_status"] = _status_series(x["decision_status"])
    if x["repletion_status"].isna().any():
        bad = sorted(x.loc[x["repletion_status"].isna(), "decision_status"].astype(str).unique())
        raise ValueError(f"Unexpected decision status values: {bad}")

    x["repleted"] = x["repletion_status"].eq("Repletion")
    x["electrolyte_group"] = pd.Categorical(
        x["lab_category"], categories=GROUP_ORDER, ordered=True
    )
    x["has_followup_24h"] = x["followup_value"].notna()
    x["delta"] = x["followup_value"] - x["lab_value"]

    # The shared outcome helper applies only clinical Below/Within/Above
    # thresholds.  It does not re-run generic plausibility filtering.
    x = add_followup_outcomes(x, key)
    return x


# -----------------------------------------------------------------------------
# QA / cohort counts
# -----------------------------------------------------------------------------


def validate_eicu_events(events: pd.DataFrame, electrolyte: str) -> None:
    """Raise if event counts differ from the frozen finalized analysis."""
    key = _key(electrolyte)
    n_events = len(events)
    n_follow = int(events["followup_value"].notna().sum())

    expected_events = EXPECTED_EVENT_COUNTS[key]
    expected_follow = EXPECTED_FOLLOWUP_COUNTS[key]

    if n_events != expected_events:
        raise AssertionError(
            f"{key}: event count {n_events:,} != frozen target {expected_events:,}"
        )
    if n_follow != expected_follow:
        raise AssertionError(
            f"{key}: follow-up count {n_follow:,} != frozen target {expected_follow:,}"
        )


def final_analytic_cohort_counts(con) -> pd.DataFrame:
    """Return the finalized eICU stay-level analytic cohort counts."""
    out = con.execute("""
        SELECT
            (SELECT COUNT(*) FROM eicu_cohort_k_final)  AS potassium,
            (SELECT COUNT(*) FROM eicu_cohort_mg_final) AS magnesium,
            (SELECT COUNT(*) FROM eicu_cohort_p_final)  AS phosphorus
    """).fetchdf()

    for key, expected in EXPECTED_FINAL_STAY_COUNTS.items():
        observed = int(out.loc[0, key])
        if observed != expected:
            raise AssertionError(
                f"{key}: final stay count {observed:,} != frozen target {expected:,}"
            )
    return out


# -----------------------------------------------------------------------------
# eICU manuscript statistics
# -----------------------------------------------------------------------------


def baseline_welch_stats(events: pd.DataFrame) -> pd.DataFrame:
    """Repletion vs No Repletion at baseline within each lab category."""
    rows = []

    for group in GROUP_ORDER:
        g = events[events["electrolyte_group"].astype(str).eq(group)]
        rep = g.loc[g["repletion_status"].eq("Repletion"), "lab_value"].dropna().astype(float)
        no = g.loc[g["repletion_status"].eq("No Repletion"), "lab_value"].dropna().astype(float)

        if len(rep) > 1 and len(no) > 1:
            test = stats.ttest_ind(rep, no, equal_var=False)
            d = cohen_d_independent(rep, no)
        else:
            test = type("Test", (), {"statistic": np.nan, "pvalue": np.nan})()
            d = np.nan

        rows.append({
            "group": group,
            "n_repletion": len(rep),
            "n_no_repletion": len(no),
            "mean_repletion": rep.mean(),
            "sd_repletion": rep.std(ddof=1),
            "mean_no_repletion": no.mean(),
            "sd_no_repletion": no.std(ddof=1),
            "mean_difference_repletion_minus_no": rep.mean() - no.mean(),
            "welch_t": test.statistic,
            "p": test.pvalue,
            "p_formatted": _p_display(test.pvalue),
            "cohen_d_signed": d,
        })

    return pd.DataFrame(rows)


def paired_pre_post_stats(events: pd.DataFrame) -> pd.DataFrame:
    """Paired baseline-to-follow-up tests and paired Cohen dz."""
    rows = []

    for group in GROUP_ORDER:
        for status in ("No Repletion", "Repletion"):
            s = events[
                events["electrolyte_group"].astype(str).eq(group)
                & events["repletion_status"].eq(status)
            ].dropna(subset=["lab_value", "followup_value"])

            pre = s["lab_value"].astype(float)
            post = s["followup_value"].astype(float)
            diff = post - pre

            if len(s) > 1:
                test = stats.ttest_rel(post, pre, nan_policy="omit")
                sd_delta = diff.std(ddof=1)
                dz = diff.mean() / sd_delta if sd_delta > 0 else np.nan
            else:
                test = type("Test", (), {"statistic": np.nan, "pvalue": np.nan})()
                dz = np.nan

            rows.append({
                "group": group,
                "repletion_status": status,
                "n_pairs": len(s),
                "baseline_mean": pre.mean(),
                "baseline_sd": pre.std(ddof=1),
                "followup_mean": post.mean(),
                "followup_sd": post.std(ddof=1),
                "delta_mean": diff.mean(),
                "delta_sd": diff.std(ddof=1),
                "paired_t": test.statistic,
                "p": test.pvalue,
                "p_formatted": _p_display(test.pvalue),
                "cohen_dz": dz,
            })

    return pd.DataFrame(rows)


def between_group_followup_change_stats(events: pd.DataFrame) -> pd.DataFrame:
    """Welch tests for follow-up value and change, Repletion vs No Repletion."""
    rows = []

    for group in GROUP_ORDER:
        g = events[
            events["electrolyte_group"].astype(str).eq(group)
            & events["followup_value"].notna()
        ].copy()

        for col, label in (
            ("followup_value", "Follow-up value"),
            ("delta", "Change from baseline"),
        ):
            rep = g.loc[g["repletion_status"].eq("Repletion"), col].dropna().astype(float)
            no = g.loc[g["repletion_status"].eq("No Repletion"), col].dropna().astype(float)

            if len(rep) > 1 and len(no) > 1:
                test = stats.ttest_ind(rep, no, equal_var=False)
                d = cohen_d_independent(rep, no)
            else:
                test = type("Test", (), {"statistic": np.nan, "pvalue": np.nan})()
                d = np.nan

            rows.append({
                "group": group,
                "metric": label,
                "n_repletion": len(rep),
                "n_no_repletion": len(no),
                "mean_repletion": rep.mean(),
                "sd_repletion": rep.std(ddof=1),
                "mean_no_repletion": no.mean(),
                "sd_no_repletion": no.std(ddof=1),
                "mean_difference_repletion_minus_no": rep.mean() - no.mean(),
                "welch_t": test.statistic,
                "p": test.pvalue,
                "p_formatted": _p_display(test.pvalue),
                "cohen_d_signed": d,
            })

    return pd.DataFrame(rows)


def followup_normal_range_stats(events: pd.DataFrame) -> pd.DataFrame:
    """Compare follow-up within-normal proportions and signed Cohen h."""
    rows = []

    for group in GROUP_ORDER:
        g = events[
            events["electrolyte_group"].astype(str).eq(group)
            & events["followup_value"].notna()
        ]

        rep = g[g["repletion_status"].eq("Repletion")]
        no = g[g["repletion_status"].eq("No Repletion")]

        n1, n0 = len(rep), len(no)
        x1 = int(rep["followup_normal"].sum())
        x0 = int(no["followup_normal"].sum())

        p1 = x1 / n1 if n1 else np.nan
        p0 = x0 / n0 if n0 else np.nan
        z, p = _two_prop_z(x1, n1, x0, n0)
        h = _cohens_h(p1, p0) if n1 and n0 else np.nan

        rows.append({
            "group": group,
            "repletion_n": n1,
            "repletion_within_normal_n": x1,
            "repletion_pct": 100.0 * p1 if n1 else np.nan,
            "no_repletion_n": n0,
            "no_repletion_within_normal_n": x0,
            "no_repletion_pct": 100.0 * p0 if n0 else np.nan,
            "absolute_difference_pp": 100.0 * (p1 - p0) if n1 and n0 else np.nan,
            "z": z,
            "p": p,
            "p_formatted": _p_display(p),
            "cohens_h": h,
            "effect": _h_effect_label(h),
        })

    return pd.DataFrame(rows)


def followup_availability_stats(events: pd.DataFrame) -> pd.DataFrame:
    """Compare 24-hour follow-up availability and signed Cohen h."""
    rows = []

    for group in GROUP_ORDER:
        g = events[events["electrolyte_group"].astype(str).eq(group)]
        rep = g[g["repletion_status"].eq("Repletion")]
        no = g[g["repletion_status"].eq("No Repletion")]

        n1, n0 = len(rep), len(no)
        x1 = int(rep["has_followup_24h"].sum())
        x0 = int(no["has_followup_24h"].sum())

        p1 = x1 / n1 if n1 else np.nan
        p0 = x0 / n0 if n0 else np.nan
        z, p = _two_prop_z(x1, n1, x0, n0)
        h = _cohens_h(p1, p0) if n1 and n0 else np.nan

        rows.append({
            "group": group,
            "repletion_n": n1,
            "repletion_followup_n": x1,
            "repletion_followup_pct": 100.0 * p1 if n1 else np.nan,
            "no_repletion_n": n0,
            "no_repletion_followup_n": x0,
            "no_repletion_followup_pct": 100.0 * p0 if n0 else np.nan,
            "absolute_difference_pp": 100.0 * (p1 - p0) if n1 and n0 else np.nan,
            "z": z,
            "p": p,
            "p_formatted": _p_display(p),
            "cohens_h": h,
            "effect": _h_effect_label(h),
        })

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Public runners
# -----------------------------------------------------------------------------


def run_eicu(con, electrolyte: str, *, validate: bool = True) -> Dict[str, pd.DataFrame]:
    """Return finalized eICU results for one electrolyte.

    The returned dictionary deliberately includes an ``events`` frame so it can
    be passed directly to ``save_figure_set`` from the shared distribution-
    figure module.
    """
    key = _key(electrolyte)
    events = extract_eicu_events(con, key)

    if validate:
        validate_eicu_events(events, key)

    return {
        "events": events,
        "main_table": manuscript_table(events),
        "baseline_stats": baseline_welch_stats(events),
        "paired_stats": paired_pre_post_stats(events),
        "between_group_stats": between_group_followup_change_stats(events),
        "normalization_stats": followup_normal_range_stats(events),
        "followup_availability_stats": followup_availability_stats(events),
        "followup_timing": followup_timing_stats(events),
    }


def run_all_eicu(con, *, validate: bool = True) -> Dict[str, Dict[str, pd.DataFrame]]:
    """Run the finalized wrapper for potassium, magnesium, and phosphorus."""
    return {
        electrolyte: run_eicu(con, electrolyte, validate=validate)
        for electrolyte in ("potassium", "magnesium", "phosphorus")
    }


if __name__ == "__main__":
    print(__doc__)
    print(
        "Open the finalized eICU DuckDB externally, then call "
        "run_eicu(con, electrolyte) or run_all_eicu(con)."
    )
