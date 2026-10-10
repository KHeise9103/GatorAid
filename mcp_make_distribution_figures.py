"""Generate GatorAId distribution figures from finalized MCP event analyses.

Run this only inside the secure MCP environment after applying the corrected
MCP exclusions and creating the three result dictionaries with ``run_mcp``.
No PHI or Mayo credentials are embedded here.

Example
-------
results_by_electrolyte = {
    "potassium": potassium_results,
    "magnesium": magnesium_results,
    "phosphorus": phosphorus_results,
}

from mcp_make_distribution_figures import make_mcp_figures
make_mcp_figures(results_by_electrolyte, "figures/mcp")
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

from gatoraid_distribution_figures import save_figure_set


def make_mcp_figures(
    results_by_electrolyte: dict[str, dict[str, pd.DataFrame]],
    output_dir: str | Path = "figures/mcp",
    *,
    include_delta: bool = True,
):
    return save_figure_set(
        results_by_electrolyte,
        dataset_label="Mayo Clinical Platform (MCP)",
        output_dir=output_dir,
        prefix="mcp",
        include_delta=include_delta,
    )


if __name__ == "__main__":
    print(__doc__)
