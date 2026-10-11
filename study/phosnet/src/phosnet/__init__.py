"""phosnet: study-effort-aware evaluation of phosphosite data and predictors.

Provides a study-effort index (SEI) for proteins, matched comparisons that
hold SEI constant, and SEI-stratified evaluation of phosphosite predictors.
"""
from .sei import compute_sei_primary, build_sei_table, build_sei_secondary_table, sei_deciles
from .matching import matched_comparison, standardized_mean_difference
from .evaluation import compute_tool_metrics, bootstrap_auroc_ci, bootstrap_gap_ci

__all__ = [
    "compute_sei_primary",
    "build_sei_table",
    "build_sei_secondary_table",
    "sei_deciles",
    "matched_comparison",
    "standardized_mean_difference",
    "compute_tool_metrics",
    "bootstrap_auroc_ci",
    "bootstrap_gap_ci",
]

__version__ = "0.1.0"
