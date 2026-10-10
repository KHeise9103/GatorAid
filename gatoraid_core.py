"""Shared GatorAId electrolyte analysis engine.

This module contains the database-agnostic analysis used by the MCP, eICU,
and MIMIC-IV wrappers in this folder.  The goal is to keep the analytic rules
identical across databases while leaving source-specific extraction and
exclusion logic in the database wrapper.

Primary analysis:
    * classify baseline laboratory values as Below / Within / Above
    * link a replacement event to the most recent preceding lab within 24 h
    * identify the next follow-up lab within 24 h
    * calculate category and treatment frequencies
    * calculate baseline and follow-up descriptive statistics
    * independent large-sample z tests + pooled-SD Cohen d at baseline
    * paired pre/post z tests + paired Cohen d
    * Mann-Whitney comparison of follow-up timing

The module intentionally does not invent medication-dose conversions.  Each
wrapper must normalize dose before calling the shared engine.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm


@dataclass(frozen=True)
class ElectrolyteSpec:
    key: str
    display_name: str
    unit: str
    low: float
    high: float
    plausible_min: float
    plausible_max: float


SPECS: Dict[str, ElectrolyteSpec] = {
    "potassium": ElectrolyteSpec(
        key="potassium", display_name="Potassium", unit="mEq/L",
        low=3.5, high=5.0, plausible_min=1.5, plausible_max=9.0,
    ),
    "magnesium": ElectrolyteSpec(
        key="magnesium", display_name="Magnesium", unit="mg/dL",
        low=1.4, high=2.2, plausible_min=0.4, plausible_max=10.0,
    ),
    "phosphorus": ElectrolyteSpec(
        key="phosphorus", display_name="Phosphorus", unit="mg/dL",
        low=3.0, high=4.5, plausible_min=0.1, plausible_max=20.0,
    ),
}

GROUP_ORDER = ["Below", "Within", "Above"]
STATUS_ORDER = ["No Repletion", "Repletion"]


def get_spec(electrolyte: str) -> ElectrolyteSpec:
    key = electrolyte.lower().strip()
    if key in {"phosphate", "phos", "pi"}:
        key = "phosphorus"
    if key in {"k", "potassium"}:
        key = "potassium"
    if key in {"mg", "magnesium"}:
        key = "magnesium"
    if key not in SPECS:
        raise ValueError(f"Unsupported electrolyte: {electrolyte!r}")
    return SPECS[key]


def format_p(p: float) -> str:
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def interpret_d(d: float) -> str:
    if pd.isna(d):
        return "NA"
    x = abs(float(d))
    if x < 0.2:
        return "Very small"
    if x < 0.5:
        return "Small"
    if x < 0.8:
        return "Moderate"
    return "Large"


def cohen_d_independent(a: Iterable[float], b: Iterable[float]) -> float:
    a = pd.Series(a, dtype="float64").dropna()
    b = pd.Series(b, dtype="float64").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    s1 = a.std(ddof=1)
    s2 = b.std(ddof=1)
    denom = len(a) + len(b) - 2
    if denom <= 0:
        return np.nan
    pooled = np.sqrt(((len(a) - 1) * s1**2 + (len(b) - 1) * s2**2) / denom)
    if pooled == 0 or pd.isna(pooled):
        return np.nan
    return float((a.mean() - b.mean()) / pooled)


def paired_z_test(pre: Iterable[float], post: Iterable[float]):
    pre = pd.Series(pre, dtype="float64")
    post = pd.Series(post, dtype="float64")
    keep = pre.notna() & post.notna()
    diff = (post[keep] - pre[keep]).astype(float)
    n = len(diff)
    if n < 2:
        return np.nan, np.nan, np.nan
    sd = diff.std(ddof=1)
    mean_delta = diff.mean()
    if sd == 0 or pd.isna(sd):
        return np.nan, np.nan, float(mean_delta)
    z = mean_delta / (sd / np.sqrt(n))
    p = 2 * (1 - norm.cdf(abs(z)))
    return float(z), float(p), float(mean_delta)


def classify_values(values: pd.Series, electrolyte: str) -> pd.Categorical:
    spec = get_spec(electrolyte)
    values = pd.to_numeric(values, errors="coerce")
    labels = np.select(
        [values < spec.low, values > spec.high],
        ["Below", "Above"],
        default="Within",
    )
    labels = pd.Series(labels, index=values.index).where(values.notna())
    return pd.Categorical(labels, categories=GROUP_ORDER, ordered=True)


def clean_labs(
    labs: pd.DataFrame,
    electrolyte: str,
    *,
    encounter_col: str = "encounter_id",
    time_col: str = "lab_time",
    value_col: str = "lab_value",
    event_id_col: str = "lab_event_id",
) -> pd.DataFrame:
    spec = get_spec(electrolyte)
    x = labs.copy()
    required = [encounter_col, time_col, value_col]
    missing = [c for c in required if c not in x.columns]
    if missing:
        raise KeyError(f"Missing lab columns: {missing}")

    x[encounter_col] = x[encounter_col].astype(str)
    x[time_col] = pd.to_datetime(x[time_col], errors="coerce", utc=True).dt.tz_localize(None)
    x[value_col] = pd.to_numeric(x[value_col], errors="coerce")
    x = x.dropna(subset=[encounter_col, time_col, value_col]).copy()
    x = x[x[value_col].between(spec.plausible_min, spec.plausible_max, inclusive="both")].copy()
    x = x.sort_values([encounter_col, time_col]).reset_index(drop=True)

    if event_id_col not in x.columns:
        x[event_id_col] = np.arange(len(x), dtype=np.int64)
    else:
        if x[event_id_col].duplicated().any():
            raise ValueError(f"{event_id_col} must uniquely identify lab events")

    x["electrolyte_group"] = classify_values(x[value_col], electrolyte)
    return x


def link_repletion_to_prior_lab(
    labs: pd.DataFrame,
    meds: pd.DataFrame,
    *,
    encounter_col: str = "encounter_id",
    lab_time_col: str = "lab_time",
    med_time_col: str = "med_time",
    event_id_col: str = "lab_event_id",
    med_name_col: str = "medication_name",
    dose_col: str = "normalized_first_dose",
    dose_unit_col: str = "normalized_dose_unit",
    window_hours: float = 24.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Link each replacement administration to its most recent prior lab.

    If multiple administrations occur at the exact first-replacement timestamp,
    normalized doses are summed only when every simultaneous component has a
    resolved normalized dose.
    """
    l = labs.copy()
    m = meds.copy()

    required_l = [encounter_col, lab_time_col, event_id_col]
    required_m = [encounter_col, med_time_col]
    for c in required_l:
        if c not in l.columns:
            raise KeyError(f"Missing lab column: {c}")
    for c in required_m:
        if c not in m.columns:
            raise KeyError(f"Missing medication column: {c}")

    l[encounter_col] = l[encounter_col].astype(str)
    m[encounter_col] = m[encounter_col].astype(str)
    l[lab_time_col] = pd.to_datetime(l[lab_time_col], errors="coerce")
    m[med_time_col] = pd.to_datetime(m[med_time_col], errors="coerce")
    l = l.dropna(subset=[lab_time_col]).copy()
    m = m.dropna(subset=[med_time_col]).copy()

    lsort = l.sort_values([lab_time_col, encounter_col]).reset_index(drop=True)
    msort = m.sort_values([med_time_col, encounter_col]).reset_index(drop=True)

    linked = pd.merge_asof(
        msort,
        lsort[[event_id_col, encounter_col, lab_time_col]],
        left_on=med_time_col,
        right_on=lab_time_col,
        by=encounter_col,
        direction="backward",
        tolerance=pd.Timedelta(hours=window_hours),
    ).dropna(subset=[event_id_col]).copy()

    if linked.empty:
        out = l.copy()
        out["first_repletion_time"] = pd.NaT
        out["repleted"] = False
        out["repletion_status"] = "No Repletion"
        out["first_dose"] = np.nan
        out["first_dose_unit"] = None
        return out, linked

    linked[event_id_col] = linked[event_id_col].astype(l[event_id_col].dtype)
    first_time = linked.groupby(event_id_col)[med_time_col].transform("min")
    first_rows = linked.loc[linked[med_time_col].eq(first_time)].copy()

    if dose_col not in first_rows.columns:
        first_rows[dose_col] = np.nan
    if dose_unit_col not in first_rows.columns:
        first_rows[dose_unit_col] = None
    if med_name_col not in first_rows.columns:
        first_rows[med_name_col] = None

    def collapse_first(g: pd.DataFrame) -> pd.Series:
        doses = pd.to_numeric(g[dose_col], errors="coerce")
        fully_resolved = doses.notna().all()
        units = g[dose_unit_col].dropna().astype(str).unique().tolist()
        return pd.Series({
            "first_repletion_time": g[med_time_col].min(),
            "first_repletion_medication": " + ".join(g[med_name_col].dropna().astype(str).unique()),
            "n_first_admins": len(g),
            "first_dose_fully_resolved": bool(fully_resolved),
            "first_dose": float(doses.sum()) if fully_resolved else np.nan,
            "first_dose_unit": units[0] if len(units) == 1 else ("mixed" if units else None),
        })

    first = first_rows.groupby(event_id_col, sort=False).apply(collapse_first, include_groups=False).reset_index()
    out = l.merge(first, on=event_id_col, how="left")
    out["repleted"] = out["first_repletion_time"].notna()
    out["repletion_status"] = np.where(out["repleted"], "Repletion", "No Repletion")
    return out, linked


def link_next_followup(
    events: pd.DataFrame,
    *,
    encounter_col: str = "encounter_id",
    lab_time_col: str = "lab_time",
    value_col: str = "lab_value",
    event_id_col: str = "lab_event_id",
    window_hours: float = 24.0,
) -> pd.DataFrame:
    base = events.copy()
    base[encounter_col] = base[encounter_col].astype(str)
    base[lab_time_col] = pd.to_datetime(base[lab_time_col], errors="coerce")

    right = base[[encounter_col, lab_time_col, value_col]].rename(
        columns={lab_time_col: "followup_time", value_col: "followup_value"}
    )

    left = base.sort_values([lab_time_col, encounter_col]).reset_index(drop=True)
    right = right.sort_values(["followup_time", encounter_col]).reset_index(drop=True)

    matched = pd.merge_asof(
        left,
        right,
        left_on=lab_time_col,
        right_on="followup_time",
        by=encounter_col,
        direction="forward",
        allow_exact_matches=False,
        tolerance=pd.Timedelta(hours=window_hours),
    )

    keep = matched[[event_id_col, "followup_time", "followup_value"]]
    out = base.merge(keep, on=event_id_col, how="left")
    out["has_followup_24h"] = out["followup_time"].notna()
    out["hours_to_followup"] = (
        out["followup_time"] - out[lab_time_col]
    ).dt.total_seconds() / 3600.0
    out["delta"] = out["followup_value"] - out[value_col]
    return out


def add_followup_outcomes(events: pd.DataFrame, electrolyte: str) -> pd.DataFrame:
    spec = get_spec(electrolyte)
    x = events.copy()
    x["followup_normal"] = x["followup_value"].between(spec.low, spec.high, inclusive="both")
    x["followup_low"] = x["followup_value"] < spec.low
    x["followup_high"] = x["followup_value"] > spec.high
    x["overcorrected"] = (
        x["electrolyte_group"].astype(str).eq("Below")
        & x["repleted"]
        & (x["followup_value"] > spec.high)
    )
    x["failed_to_correct"] = (
        x["electrolyte_group"].astype(str).eq("Below")
        & (x["followup_value"] < spec.low)
    )
    return x


def manuscript_table(events: pd.DataFrame) -> pd.DataFrame:
    rows = []
    total_labs = len(events)

    for group in GROUP_ORDER:
        g = events[events["electrolyte_group"].astype(str).eq(group)].copy()
        group_total = len(g)
        pct_group = 100 * group_total / total_labs if total_labs else np.nan
        pct_treated = 100 * g["repleted"].mean() if group_total else np.nan

        for status in STATUS_ORDER:
            s = g[g["repletion_status"].eq(status)].copy()
            n = len(s)
            n_follow = s["followup_value"].notna().sum()
            first_dose = pd.to_numeric(s.get("first_dose"), errors="coerce") if "first_dose" in s else pd.Series(dtype=float)

            rows.append({
                "group": group,
                "repletion_status": status,
                "n_labs": n,
                "pct_of_all_labs_category": pct_group,
                "pct_treated_labs_category": pct_treated,
                "pct_within_category": 100 * n / group_total if group_total else np.nan,
                "baseline_mean": s["lab_value"].mean(),
                "baseline_sd": s["lab_value"].std(),
                "n_with_followup": int(n_follow),
                "followup_rate": n_follow / n if n else np.nan,
                "median_time_followup_hours": s["hours_to_followup"].median(),
                "followup_mean": s["followup_value"].mean(),
                "followup_sd": s["followup_value"].std(),
                "delta_mean": s["delta"].mean(),
                "delta_sd": s["delta"].std(),
                "first_dose_n_resolved": int(first_dose.notna().sum()) if len(first_dose) else 0,
                "first_dose_median": first_dose.median() if len(first_dose) else np.nan,
                "first_dose_q25": first_dose.quantile(0.25) if len(first_dose) else np.nan,
                "first_dose_q75": first_dose.quantile(0.75) if len(first_dose) else np.nan,
            })

    return pd.DataFrame(rows)


def baseline_repletion_stats(events: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group in GROUP_ORDER:
        rep = events[
            events["electrolyte_group"].astype(str).eq(group)
            & events["repletion_status"].eq("Repletion")
        ]["lab_value"].dropna()
        no = events[
            events["electrolyte_group"].astype(str).eq(group)
            & events["repletion_status"].eq("No Repletion")
        ]["lab_value"].dropna()

        if len(rep) > 1 and len(no) > 1:
            se = np.sqrt(rep.var(ddof=1) / len(rep) + no.var(ddof=1) / len(no))
            z = (rep.mean() - no.mean()) / se if se > 0 else np.nan
            p = 2 * (1 - norm.cdf(abs(z))) if pd.notna(z) else np.nan
        else:
            z = p = np.nan
        d = cohen_d_independent(rep, no)
        rows.append({
            "group": group,
            "n_repletion": len(rep),
            "n_no_repletion": len(no),
            "mean_repletion": rep.mean(),
            "mean_no_repletion": no.mean(),
            "z": z,
            "p": p,
            "p_formatted": format_p(p),
            "cohen_d_signed": d,
            "cohen_d_magnitude": abs(d) if pd.notna(d) else np.nan,
            "effect_size": interpret_d(d),
        })
    return pd.DataFrame(rows)


def paired_pre_post_stats(events: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group in GROUP_ORDER:
        for status in STATUS_ORDER:
            s = events[
                events["electrolyte_group"].astype(str).eq(group)
                & events["repletion_status"].eq(status)
            ].dropna(subset=["lab_value", "followup_value"])
            pre = s["lab_value"]
            post = s["followup_value"]
            diff = post - pre
            z, p, _ = paired_z_test(pre, post)
            if len(s) > 1:
                t_stat, t_p = stats.ttest_rel(post, pre, nan_policy="omit")
                d = diff.mean() / diff.std(ddof=1) if diff.std(ddof=1) > 0 else np.nan
            else:
                t_stat = t_p = d = np.nan
            rows.append({
                "group": group,
                "repletion_status": status,
                "n": len(s),
                "baseline_mean": pre.mean(),
                "baseline_sd": pre.std(),
                "followup_mean": post.mean(),
                "followup_sd": post.std(),
                "delta_mean": diff.mean(),
                "delta_sd": diff.std(),
                "z": z,
                "z_p": p,
                "z_p_formatted": format_p(p),
                "paired_t": t_stat,
                "paired_t_p": t_p,
                "paired_t_p_formatted": format_p(t_p),
                "cohen_d_paired": d,
                "effect_size": interpret_d(d),
            })
    return pd.DataFrame(rows)


def followup_timing_stats(events: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group in GROUP_ORDER:
        rep = events[
            events["electrolyte_group"].astype(str).eq(group)
            & events["repletion_status"].eq("Repletion")
        ]["hours_to_followup"].dropna()
        no = events[
            events["electrolyte_group"].astype(str).eq(group)
            & events["repletion_status"].eq("No Repletion")
        ]["hours_to_followup"].dropna()
        if len(rep) and len(no):
            u, p = stats.mannwhitneyu(rep, no, alternative="two-sided")
            rbc = 2 * u / (len(rep) * len(no)) - 1
        else:
            u = p = rbc = np.nan
        rows.append({
            "group": group,
            "n_repletion": len(rep),
            "n_no_repletion": len(no),
            "repletion_median_hours": rep.median(),
            "no_repletion_median_hours": no.median(),
            "mann_whitney_u": u,
            "p": p,
            "p_formatted": format_p(p),
            "rank_biserial": rbc,
        })
    return pd.DataFrame(rows)


def run_event_analysis(
    labs: pd.DataFrame,
    meds: pd.DataFrame,
    electrolyte: str,
    *,
    encounter_col: str = "encounter_id",
    lab_time_col: str = "lab_time",
    value_col: str = "lab_value",
    med_time_col: str = "med_time",
    med_name_col: str = "medication_name",
    dose_col: str = "normalized_first_dose",
    dose_unit_col: str = "normalized_dose_unit",
) -> dict[str, pd.DataFrame]:
    """Run the primary event-level GatorAId analysis on standardized frames."""
    labs_clean = clean_labs(
        labs, electrolyte,
        encounter_col=encounter_col,
        time_col=lab_time_col,
        value_col=value_col,
    ).rename(columns={encounter_col: "encounter_id", lab_time_col: "lab_time", value_col: "lab_value"})

    meds_std = meds.rename(columns={encounter_col: "encounter_id", med_time_col: "med_time"}).copy()
    if med_name_col != "medication_name" and med_name_col in meds_std.columns:
        meds_std = meds_std.rename(columns={med_name_col: "medication_name"})
    if dose_col != "normalized_first_dose" and dose_col in meds_std.columns:
        meds_std = meds_std.rename(columns={dose_col: "normalized_first_dose"})
    if dose_unit_col != "normalized_dose_unit" and dose_unit_col in meds_std.columns:
        meds_std = meds_std.rename(columns={dose_unit_col: "normalized_dose_unit"})

    events, linked_meds = link_repletion_to_prior_lab(labs_clean, meds_std)
    events = link_next_followup(events)
    events = add_followup_outcomes(events, electrolyte)

    return {
        "events": events,
        "linked_medications": linked_meds,
        "main_table": manuscript_table(events),
        "baseline_stats": baseline_repletion_stats(events),
        "paired_stats": paired_pre_post_stats(events),
        "followup_timing": followup_timing_stats(events),
    }


def export_results(results: dict[str, pd.DataFrame], output_dir: str | Path, prefix: str) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, df in results.items():
        if isinstance(df, pd.DataFrame):
            df.to_csv(out / f"{prefix}_{name}.csv", index=False)
