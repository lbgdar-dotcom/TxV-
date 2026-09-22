"""Matched construct series for the two experiments the panel cannot answer.

The P0-P8 panel settles *where* an antigen goes. It does not answer the two
questions that decide whether a multi-antigen vaccine is actually potent:

**1. How many antigens can one transcript present before presentation decays?**
Published routing work is mostly one or two antigens. The design question for a
real multi-antigen vaccine is whether there is a *presentation budget* -- a
point past which adding antigens dilutes per-epitope pMHC rather than adding
breadth. :func:`multiplex_ladder` builds the series that measures it.

**2. Does uridine depletion trade expression against adjuvanticity?**
Depleting uridine raises expression and lowers innate sensing. For a
prophylactic vaccine that is unambiguously good. For a *cancer* vaccine it may
not be: the same sensing that you suppress is part of what matures the
dendritic cell and licenses it to prime CD8 T cells. So the potency-optimal
uridine content may not be the expression-optimal one.
:func:`uridine_ladder` builds the series that separates them.

Both functions return matched series -- one variable moved, everything else
held -- because that is the only form in which the answer is interpretable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

from .codon import CodonOptimizer, OptimizerConfig, uridine_floor
from .constructs import Construct, ConstructBuilder, ConstructSpec
from .epitopes import Antigen, Cassette, junction_risk
from .parts import PartRegistry, default_registry
from .qc import QCReport, run_qc


@dataclass
class Variant:
    """One arm of a matched series."""

    name: str
    construct: Construct
    qc: QCReport
    #: The variable this arm moves, and its realised value.
    axis: str = ""
    value: float | int = 0
    #: Covariates that moved along with it and must be reported, not hidden.
    covariates: dict[str, float] = field(default_factory=dict)


@dataclass
class VariantSeries:
    axis: str
    variants: list[Variant]
    note: str = ""

    def to_text(self) -> str:
        lines = [f"Series: {self.axis}", f"  {self.note}" if self.note else ""]
        if not self.variants:
            return "\n".join(l for l in lines if l)
        covariate_names = list(self.variants[0].covariates)
        header = f"  {'arm':<16}{self.axis:>12}"
        header += "".join(f"{c:>10}" for c in covariate_names)
        header += f"{'QC':>7}"
        lines.append(header)
        for variant in self.variants:
            row = f"  {variant.name:<16}{variant.value:>12.4g}"
            row += "".join(
                f"{variant.covariates[c]:>10.4g}" for c in covariate_names
            )
            row += f"{'PASS' if variant.qc.passed else 'FAIL':>7}"
            lines.append(row)
        return "\n".join(l for l in lines if l)

    def to_dict(self) -> dict:
        return {
            "axis": self.axis,
            "note": self.note,
            "arms": [
                {
                    "name": v.name, "value": v.value,
                    "covariates": v.covariates,
                    "construct": v.construct.summary(),
                    "qc_passed": v.qc.passed,
                }
                for v in self.variants
            ],
        }


def _routed(route: str, name: str):
    """Default spec/registry for a series: a real route, no placeholders."""
    from .pvax1_ag import routed_spec

    return routed_spec(route, name=name)


def _build(name, cassette, spec, registry, config):
    builder = ConstructBuilder(
        registry=registry, optimizer=CodonOptimizer(config) if config else None
    )
    spec.name = name
    return builder.build(spec, cassette)


# ---------------------------------------------------------------------------
# Axis 1: uridine content
# ---------------------------------------------------------------------------

def uridine_ladder(
    antigens: Sequence[Antigen],
    weights: Sequence[float] = (0.0, 0.25, 0.5, 1.0, 2.0),
    linker: str = "GGSGGGGSGG",
    route: str = "ctla4",
    spec: ConstructSpec | None = None,
    registry: PartRegistry | None = None,
    base_config: OptimizerConfig | None = None,
    name_prefix: str = "U",
) -> VariantSeries:
    """Matched ORFs encoding the identical protein at descending uridine content.

    Each arm sweeps ``w_uridine``; everything else in the optimiser is held.
    The encoded protein is byte-identical across arms, so any difference in a
    cell assay is attributable to the nucleotide sequence, not the antigen.

    **Uridine cannot be varied independently of codon adaptation** -- depleting
    U forces synonymous choices away from the preferred codon. That confound is
    real and unavoidable, so the series reports CAI and GC per arm as
    covariates rather than pretending the axis is clean. Read the result as a
    two-variable response surface, and if an arm's CAI has collapsed, that arm
    is confounded and should be excluded, not explained away.
    """
    if spec is None or registry is None:
        default_spec, default_registry_ = _routed(route, name_prefix)
        spec = spec or default_spec
        registry = registry or default_registry_
    base = base_config or OptimizerConfig()
    cassette = Cassette(list(antigens), linker=linker)

    variants: list[Variant] = []
    for weight in weights:
        config = OptimizerConfig(**{**base.__dict__, "w_uridine": weight})
        name = f"{name_prefix}{weight:g}".replace(".", "p")
        construct = _build(
            name, cassette, spec or ConstructSpec(name=name), registry, config
        )
        transcript = construct.transcript
        uridine = transcript.count("T") / len(transcript)
        variants.append(Variant(
            name=name,
            construct=construct,
            qc=run_qc(construct),
            axis="w_uridine",
            value=weight,
            covariates={
                "U_frac": round(uridine, 4),
                "ORF_U": construct.orf.count("T"),
                "floor_U": uridine_floor(construct.orf_protein)[0],
                "CAI": round(construct.optimization.cai, 4),
                "GC": round(construct.summary()["gc"], 4),
            },
        ))
    return VariantSeries(
        axis="w_uridine",
        variants=variants,
        note="Identical protein, descending uridine. Expression and innate "
             "sensing move in opposite directions across this axis; CAI is a "
             "confound and is reported per arm. floor_U is the minimum "
             "uridines the encoded protein admits: an arm sitting on it cannot "
             "be pushed lower without changing the protein, so more weight "
             "buys nothing.",
    )


# ---------------------------------------------------------------------------
# Axis 2: antigen multiplicity
# ---------------------------------------------------------------------------

def multiplex_ladder(
    antigens: Sequence[Antigen],
    sizes: Sequence[int] = (1, 2, 4, 6, 8),
    anchor: Antigen | None = None,
    linker: str = "GGSGGGGSGG",
    route: str = "ctla4",
    spec: ConstructSpec | None = None,
    registry: PartRegistry | None = None,
    config: OptimizerConfig | None = None,
    name_prefix: str = "N",
) -> VariantSeries:
    """Nested cassettes of increasing antigen count, on one fixed route.

    Two properties make the series interpretable, and both are deliberate:

    * **Nesting.** Arm *k* contains every antigen of arm *k-1*. Without this,
      a drop in per-epitope presentation could just be a different antigen set.
    * **An anchor antigen** present in every arm, at the same cassette position
      (first, right after the signal peptide). It is the internal standard: its
      pMHC signal is what you normalise the others against, so a global drop in
      translation is distinguishable from a genuine per-epitope dilution.

    Pass ``anchor`` explicitly; otherwise the first antigen is used and is
    carried into every arm.
    """
    if spec is None or registry is None:
        default_spec, default_registry_ = _routed(route, name_prefix)
        spec = spec or default_spec
        registry = registry or default_registry_
    pool = list(antigens)
    if not pool:
        raise ValueError("no antigens supplied")
    anchor = anchor or pool[0]
    rest = [a for a in pool if a.name != anchor.name]

    variants: list[Variant] = []
    for size in sorted(sizes):
        if size < 1:
            raise ValueError("cassette size must be at least 1")
        if size > len(rest) + 1:
            raise ValueError(
                f"cassette size {size} exceeds the {len(rest) + 1} antigens available"
            )
        members = [anchor] + rest[: size - 1]
        cassette = Cassette(members, linker=linker)
        name = f"{name_prefix}{size}"
        construct = _build(
            name, cassette, spec or ConstructSpec(name=name), registry, config
        )
        variants.append(Variant(
            name=name,
            construct=construct,
            qc=run_qc(construct),
            axis="n_antigens",
            value=size,
            covariates={
                "ORF_nt": len(construct.orf),
                "tx_nt": construct.transcript_length,
                "junction_risk": round(junction_risk(cassette), 3),
                "CAI": round(construct.optimization.cai, 4),
            },
        ))
    return VariantSeries(
        axis="n_antigens",
        variants=variants,
        note=f"Nested cassettes sharing the anchor antigen {anchor.name!r} in "
             "position 1. Normalise per-epitope readouts to the anchor so a "
             "global translation drop is separable from per-epitope dilution.",
    )


def series_outputs(series: VariantSeries, out_dir) -> dict[str, str]:
    """Write GenBank and mRNA FASTA for every arm, plus the series table."""
    import json
    from pathlib import Path

    from .export import to_fasta, to_genbank

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, str] = {}
    for variant in series.variants:
        construct = variant.construct
        for suffix, body in (
            (".gb", to_genbank(construct)),
            (".mrna.fasta", to_fasta(construct, "mrna")),
        ):
            path = out / f"{construct.name}{suffix}"
            path.write_text(body)
            written[path.name] = str(path)
    table = out / f"{series.axis}_series.json"
    table.write_text(json.dumps(series.to_dict(), indent=2))
    written[table.name] = str(table)
    return written


__all__ = [
    "Variant", "VariantSeries", "uridine_ladder", "multiplex_ladder",
    "series_outputs",
]
