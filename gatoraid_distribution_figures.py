"""Publication-style distribution figures for GatorAId.

Creates histogram + normal-reference-curve figures from the *actual event-level*
analysis data produced by :mod:`gatoraid_core`.  No values are simulated from
summary statistics.

Primary figure types
--------------------
1. Baseline electrolyte distribution with clinical Below / Within / Above
   thresholds marked.
2. Change (follow-up - baseline) distribution, generated separately for
   Repletion and No Repletion.

The overlaid normal curve is descriptive only.  It is the normal density with
mean and SD estimated from the observed values and is not a normality test.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from gatoraid_core import get_spec


def _finite_numeric(values: Iterable[float]) -> np.ndarray:
    x = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)
    return x[np.isfinite(x)]


def _fd_bins(values: np.ndarray, *, minimum: int = 12, maximum: int = 80) -> int:
    """Freedman-Diaconis bin count, bounded for publication readability."""
    if len(values) < 2:
        return minimum
    q25, q75 = np.quantile(values, [0.25, 0.75])
    iqr = q75 - q25
    if iqr <= 0:
        return min(maximum, max(minimum, int(np.sqrt(len(values)))))
    width = 2 * iqr * (len(values) ** (-1 / 3))
    if width <= 0:
        return minimum
    n = int(np.ceil((values.max() - values.min()) / width))
    return min(maximum, max(minimum, n))


def _save(fig: plt.Figure, stem: Path, *, dpi: int = 350) -> dict[str, Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"png": png, "pdf": pdf}


def _normal_reference(ax: plt.Axes, values: np.ndarray) -> None:
    if len(values) < 3:
        return
    mu = float(np.mean(values))
    sd = float(np.std(values, ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        return
    x = np.linspace(values.min(), values.max(), 600)
    ax.plot(x, norm.pdf(x, loc=mu, scale=sd), linewidth=2.0, label="Normal reference")


def baseline_distribution_figure(
    events: pd.DataFrame,
    electrolyte: str,
    dataset_label: str,
    *,
    value_col: str = "lab_value",
) -> plt.Figure:
    """Histogram of baseline labs with manuscript clinical thresholds."""
    spec = get_spec(electrolyte)
    if value_col not in events.columns:
        raise KeyError(f"Missing {value_col!r} in events")

    values = _finite_numeric(events[value_col])
    if len(values) == 0:
        raise ValueError("No numeric baseline values available for plotting")

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.hist(values, bins=_fd_bins(values), density=True, alpha=0.55)
    _normal_reference(ax, values)

    ax.axvline(spec.low, linestyle="--", linewidth=1.4, label=f"Lower threshold ({spec.low:g})")
    ax.axvline(spec.high, linestyle="--", linewidth=1.4, label=f"Upper threshold ({spec.high:g})")

    ymin, ymax = ax.get_ylim()
    text_y = ymax * 0.94
    xmin, xmax = ax.get_xlim()
    ax.text((xmin + spec.low) / 2, text_y, "Below", ha="center", va="top")
    ax.text((spec.low + spec.high) / 2, text_y, "Within", ha="center", va="top")
    ax.text((spec.high + xmax) / 2, text_y, "Above", ha="center", va="top")

    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
    ax.set_title(f"{dataset_label}: {spec.display_name} baseline distribution")
    ax.set_xlabel(f"Serum {spec.display_name.lower()} ({spec.unit})")
    ax.set_ylabel("Density")
    ax.text(
        0.01,
        -0.17,
        f"N={len(values):,}; mean={mean:.3f}; SD={sd:.3f}. Normal curve is a descriptive reference, not a normality test.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
    )
    ax.grid(axis="y", alpha=0.18)
    ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    return fig


def delta_distribution_figure(
    events: pd.DataFrame,
    electrolyte: str,
    dataset_label: str,
    status: str,
    *,
    delta_col: str = "delta",
    status_col: str = "repletion_status",
) -> plt.Figure:
    """Histogram of follow-up minus baseline values for one treatment status."""
    spec = get_spec(electrolyte)
    if delta_col not in events.columns:
        raise KeyError(f"Missing {delta_col!r} in events")
    if status_col not in events.columns:
        raise KeyError(f"Missing {status_col!r} in events")

    subset = events.loc[events[status_col].astype(str).eq(status), delta_col]
    values = _finite_numeric(subset)
    if len(values) == 0:
        raise ValueError(f"No numeric delta values available for {status}")

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.hist(values, bins=_fd_bins(values), density=True, alpha=0.55)
    _normal_reference(ax, values)
    ax.axvline(0, linestyle="--", linewidth=1.4, label="No change")

    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
    ax.set_title(f"{dataset_label}: {spec.display_name} change — {status}")
    ax.set_xlabel(f"Follow-up - baseline {spec.display_name.lower()} ({spec.unit})")
    ax.set_ylabel("Density")
    ax.text(
        0.01,
        -0.17,
        f"N={len(values):,}; mean change={mean:+.3f}; SD={sd:.3f}. Follow-up window: 24 hours.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
    )
    ax.grid(axis="y", alpha=0.18)
    ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    return fig


def save_baseline_distribution(
    events: pd.DataFrame,
    electrolyte: str,
    dataset_label: str,
    output_dir: str | Path,
    *,
    prefix: str,
) -> dict[str, Path]:
    fig = baseline_distribution_figure(events, electrolyte, dataset_label)
    stem = Path(output_dir) / f"{prefix}_{get_spec(electrolyte).key}_baseline_distribution"
    return _save(fig, stem)


def save_delta_distribution(
    events: pd.DataFrame,
    electrolyte: str,
    dataset_label: str,
    status: str,
    output_dir: str | Path,
    *,
    prefix: str,
) -> dict[str, Path]:
    fig = delta_distribution_figure(events, electrolyte, dataset_label, status)
    status_slug = status.lower().replace(" ", "_")
    stem = Path(output_dir) / f"{prefix}_{get_spec(electrolyte).key}_delta_{status_slug}"
    return _save(fig, stem)


def save_figure_set(
    results_by_electrolyte: dict[str, dict[str, pd.DataFrame]],
    dataset_label: str,
    output_dir: str | Path,
    *,
    prefix: str,
    include_delta: bool = True,
) -> list[dict[str, object]]:
    """Save separate publication-ready figures for K, Mg, and phosphorus.

    ``results_by_electrolyte`` should map electrolyte name to the dictionary
    returned by ``run_event_analysis`` / ``run_mcp`` / ``run_mimic_iv``.
    """
    made: list[dict[str, object]] = []
    for electrolyte in ("potassium", "magnesium", "phosphorus"):
        if electrolyte not in results_by_electrolyte:
            continue
        results = results_by_electrolyte[electrolyte]
        events = results["events"]
        paths = save_baseline_distribution(
            events, electrolyte, dataset_label, output_dir, prefix=prefix
        )
        made.append({"electrolyte": electrolyte, "figure": "baseline", **paths})

        if include_delta:
            for status in ("No Repletion", "Repletion"):
                try:
                    paths = save_delta_distribution(
                        events, electrolyte, dataset_label, status, output_dir, prefix=prefix
                    )
                except ValueError:
                    continue
                made.append({
                    "electrolyte": electrolyte,
                    "figure": f"delta_{status.lower().replace(' ', '_')}",
                    **paths,
                })
    return made
