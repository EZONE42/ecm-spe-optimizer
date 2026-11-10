"""ECM-SPE: Solid Polymer Electrolyte Composition Optimization via Equivalent Circuit Modeling.

Combines 4-resistance equivalent circuit modeling (bulk membrane, cathode porous transmission line,
Bruce-Vincent concentration polarization, and lithium interface ASR) with machine learning
extrapolation analysis to discover optimal SPE compositions and identify true cell-level bottlenecks.
"""

__version__ = "1.0.0"
__author__ = "Jiwon Lee"

from .ecm import (
    CellParameters,
    de_levie_tlm,
    concentration_polarization_resistance,
    cell_resistance_breakdown,
    optimize_composition,
)
from .transference import (
    circle_fit,
    bruce_vincent_tplus,
    interfacial_asr,
)

__all__ = [
    "CellParameters",
    "de_levie_tlm",
    "concentration_polarization_resistance",
    "cell_resistance_breakdown",
    "optimize_composition",
    "circle_fit",
    "bruce_vincent_tplus",
    "interfacial_asr",
]
