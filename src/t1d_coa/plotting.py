"""Shared figure styling so plots look consistent across the pipeline."""

from __future__ import annotations

# Cell-type → color, used in every per-compartment figure.
COMPARTMENT_COLORS = {
    "beta cell": "#2ca02c",
    "alpha cell": "#9467bd",
    "delta cell": "#8c564b",
    "leukocyte": "#d62728",
    "duct epithelial cell": "#7f7f7f",
    "endothelial cell": "#17becf",
    "mesenchymal cell": "#bcbd22",
    "acinar cell": "#ff7f0e",
    "PP cell": "#e377c2",
    "epsilon cell": "#1f77b4",
}


def style_axes(ax) -> None:
    """Drop top/right spines for the publication-ready look."""
    ax.spines[["top", "right"]].set_visible(False)
