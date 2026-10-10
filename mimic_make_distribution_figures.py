"""Generate GatorAId distribution figures from finalized MIMIC-IV analyses.

The figures use actual event-level values.  Run after the MIMIC-IV medication
QA is finalized and after creating the three result dictionaries with
``run_mimic_iv``.

Example
-------
results_by_electrolyte = {
    "potassium": potassium_results,
    "magnesium": magnesium_results,
    "phosphorus": phosphorus_results,
}

from mimic_make_distribution_figures import make_mimic_figures
make_mimic_figures(results_by_electrolyte, "figures/mimic_iv")
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

from gatoraid_distribution_figures import save_figure_set


def make_mimic_figures(
    results_by_electrolyte: dict[str, dict[str, pd.DataFrame]],
    output_dir: str | Path = "figures/mimic_iv",
    *,
    include_delta: bool = True,
):
    return save_figure_set(
        results_by_electrolyte,
        dataset_label="MIMIC-IV",
        output_dir=output_dir,
        prefix="mimic_iv",
        include_delta=include_delta,
    )


if __name__ == "__main__":
    print(__doc__)
