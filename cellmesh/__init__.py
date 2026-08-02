"""
CELL MESH: Metabolite-mediated Event Scoring with Sensor Hierarchies.
一个用于单细胞代谢物通信分析的工具包，完全基于 metabolite availability 算法
"""
from .core import CellMeshResult, run_cell_mesh
from .database import load_cell_mesh_database, load_default_priors
from .io import read_anndata, read_example_data
from .score import compute_metabolite_availability
from .plotting import (
    plot_communication_network,
    plot_event_dotplot,
    plot_metabolite_secretion_violin,
    plot_receptor_expression_violin,
    plot_sample_event_scores,
    plot_significant_event_counts,
)
from .config import MIN_CELL_COUNT, DATA_DIR

__version__ = "0.4.0"
__all__ = [
    "CellMeshResult",
    "run_cell_mesh",
    "load_cell_mesh_database",
    "load_default_priors",
    "compute_metabolite_availability",
    "plot_significant_event_counts",
    "plot_communication_network",
    "plot_event_dotplot",
    "plot_sample_event_scores",
    "plot_metabolite_secretion_violin",
    "plot_receptor_expression_violin",
    "read_anndata",
    "read_example_data",
    "MIN_CELL_COUNT",
    "DATA_DIR"
]
