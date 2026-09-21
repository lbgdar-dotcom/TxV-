"""txv -- design and QC of IVT mRNA constructs for multi-antigen cancer vaccines.

Typical use::

    from txv.epitopes import Antigen
    from txv.pipeline import design_construct, write_outputs

    report = design_construct("TXV-001", [Antigen("KRAS_G12D", "LVVVGADGVGKSALTIQ")])
    print(report.to_text())
    write_outputs(report, "out/")

See the README for the construct architecture and for what each module does and
does not claim.
"""

__version__ = "0.1.0"

from .constructs import Construct, ConstructBuilder, ConstructSpec, build_construct
from .epitopes import AnchorMotifScorer, Antigen, Cassette
from .codon import CodonOptimizer, OptimizerConfig
from .lnp import formulate
from .parts import PartRegistry, default_registry
from .pipeline import design_construct, load_antigens, write_outputs
from .qc import QCThresholds, run_qc

__all__ = [
    "Antigen", "Cassette", "AnchorMotifScorer",
    "CodonOptimizer", "OptimizerConfig",
    "Construct", "ConstructBuilder", "ConstructSpec", "build_construct",
    "PartRegistry", "default_registry",
    "design_construct", "load_antigens", "write_outputs",
    "run_qc", "QCThresholds", "formulate",
    "__version__",
]
