"""Assembly of a complete IVT mRNA DNA template from parts and a payload.

Template layout, 5' -> 3'::

    [T7 promoter] [5' UTR] [Kozak] [==== ORF ====] [3' UTR] [poly(A)] [linearisation site]
                                    |    |      |
                       signal peptide  cassette  trafficking domain + stop

The whole ORF -- signal peptide, every bead, every linker and the trafficking
domain -- is codon-optimised in a **single pass**, not element by element. That
matters: motifs and GC excursions are created at element *boundaries*, and an
element-wise optimiser cannot see them. The price is that you cannot swap one
bead without re-optimising the construct, which is the correct trade for a
drug substance anyway.

Coordinates on :class:`Feature` are 0-based half-open indices into
:attr:`Construct.template`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .codon import CodonOptimizer, OptimizationResult, OptimizerConfig
from .epitopes import Antigen, Cassette
from .parts import Part, PartKind, PartRegistry, Provenance, default_registry
from .seqops import clean, gc_fraction, to_rna, translate


@dataclass(frozen=True)
class Feature:
    name: str
    start: int
    end: int
    kind: str
    note: str = ""

    def __len__(self) -> int:
        return self.end - self.start

    def slice(self, seq: str) -> str:
        return seq[self.start : self.end]


@dataclass
class Construct:
    """An assembled DNA template plus its annotation and design provenance."""

    name: str
    template: str
    features: list[Feature]
    orf_protein: str
    optimization: OptimizationResult
    cassette: Cassette
    transcript_start: int
    transcript_end: int
    placeholders_used: list[str] = field(default_factory=list)
    design_notes: list[str] = field(default_factory=list)

    # -- derived views ------------------------------------------------------
    @property
    def transcript(self) -> str:
        """The DNA-alphabet sequence of the run-off transcript."""
        return self.template[self.transcript_start : self.transcript_end]

    @property
    def mrna(self) -> str:
        """The transcript in RNA alphabet (no cap or tail chemistry implied)."""
        return to_rna(self.transcript)

    @property
    def orf(self) -> str:
        feature = self.feature("ORF")
        return feature.slice(self.template)

    @property
    def transcript_length(self) -> int:
        return self.transcript_end - self.transcript_start

    def feature(self, name: str) -> Feature:
        for feature in self.features:
            if feature.name == name:
                return feature
        raise KeyError(f"no feature {name!r} in {self.name}")

    def features_of_kind(self, kind: str) -> list[Feature]:
        return [f for f in self.features if f.kind == kind]

    def summary(self) -> dict[str, object]:
        return {
            "name": self.name,
            "template_bp": len(self.template),
            "transcript_nt": self.transcript_length,
            "orf_nt": len(self.orf),
            "protein_aa": len(self.orf_protein),
            "n_antigens": len(self.cassette.antigens),
            "gc": round(gc_fraction(self.transcript), 4),
            "uridine_fraction": round(self.transcript.count("T") / self.transcript_length, 4),
            "cai": round(self.optimization.cai, 4),
            "placeholders": list(self.placeholders_used),
        }


@dataclass
class ConstructSpec:
    """Everything that defines a construct except the payload sequences."""

    name: str
    promoter: str = "T7_promoter"
    utr5: str = "UTR5_placeholder"
    kozak: str = "kozak_strong"
    signal_peptide: str | None = "SP_tPA"
    trafficking: str | None = "MITD_placeholder"
    utr3: str = "UTR3_placeholder"
    polya: str = "polyA_100"
    linearization_site: str | None = "BspQI_site"
    stop_codon: str = "TGA"
    #: Append a second stop in the same frame. Cheap insurance against
    #: read-through on a long, highly expressed ORF.
    tandem_stop: bool = True
    #: Extra linker inserted between the signal peptide and the first bead so
    #: signal peptidase has room to cut cleanly.
    sp_spacer: str | None = "linker_GS5"
    #: Linker between the last bead and the trafficking domain.
    traffic_spacer: str | None = "linker_GS5"


class ConstructBuilder:
    """Turns a :class:`ConstructSpec` plus a :class:`Cassette` into a template."""

    def __init__(
        self,
        registry: PartRegistry | None = None,
        optimizer: CodonOptimizer | None = None,
    ) -> None:
        self.registry = registry or default_registry()
        self.optimizer = optimizer or CodonOptimizer()

    def build(self, spec: ConstructSpec, cassette: Cassette) -> Construct:
        if not cassette.antigens:
            raise ValueError("cassette has no antigens")

        placeholders: list[str] = []
        notes: list[str] = []

        def part(name: str | None) -> Part | None:
            if name is None:
                return None
            got = self.registry.get(name)
            if got.is_placeholder:
                placeholders.append(got.name)
            return got

        promoter = part(spec.promoter)
        utr5 = part(spec.utr5)
        kozak = part(spec.kozak)
        sp = part(spec.signal_peptide)
        sp_spacer = part(spec.sp_spacer) if spec.signal_peptide else None
        traffic = part(spec.trafficking)
        traffic_spacer = part(spec.traffic_spacer) if spec.trafficking else None
        utr3 = part(spec.utr3)
        polya = part(spec.polya)
        linearizer = part(spec.linearization_site)

        # ---- protein-level ORF layout ---------------------------------------
        segments: list[tuple[str, str, str]] = []  # (name, protein, kind)
        if sp:
            segments.append((sp.name, sp.protein or "", "signal_peptide"))
            if sp_spacer:
                segments.append((sp_spacer.name, sp_spacer.protein or "", "linker"))
        for index, antigen in enumerate(cassette.antigens):
            if index:
                link = cassette.linker_at(index - 1)
                segments.append((f"linker_{index - 1}", link, "linker"))
            segments.append((antigen.name, antigen.sequence, antigen.kind))
        if traffic:
            if traffic_spacer:
                segments.append((traffic_spacer.name, traffic_spacer.protein or "", "linker"))
            segments.append((traffic.name, traffic.protein or "", "trafficking"))

        orf_protein = "".join(seg for _, seg, _ in segments)
        if not orf_protein.startswith("M"):
            orf_protein = "M" + orf_protein
            segments.insert(0, ("start_Met", "M", "start"))
            notes.append("Prepended an initiator methionine: no element supplied one.")

        # ---- one-pass codon optimisation ------------------------------------
        stop = clean(spec.stop_codon)
        if stop not in ("TAA", "TAG", "TGA"):
            raise ValueError(f"not a stop codon: {spec.stop_codon}")
        result = self.optimizer.optimize(orf_protein, add_stop=stop)
        orf_dna = result.dna
        if spec.tandem_stop:
            orf_dna += "TAA" if stop != "TAA" else "TGA"
            notes.append("Tandem stop codon added downstream of the primary stop.")

        # Map protein segment boundaries onto ORF nucleotide coordinates.
        orf_segment_spans: list[tuple[str, int, int, str]] = []
        aa_cursor = 0
        for seg_name, seg_protein, seg_kind in segments:
            if not seg_protein:
                continue
            start = aa_cursor * 3
            end = (aa_cursor + len(seg_protein)) * 3
            orf_segment_spans.append((seg_name, start, end, seg_kind))
            aa_cursor += len(seg_protein)

        # ---- template assembly ----------------------------------------------
        pieces: list[tuple[str, str, str, str]] = []  # (name, dna, kind, note)
        if promoter:
            pieces.append((promoter.name, promoter.dna or "", "promoter", promoter.note))
        pieces.append((utr5.name, utr5.dna or "", "utr5", utr5.note))
        if kozak:
            pieces.append((kozak.name, kozak.dna or "", "kozak", kozak.note))
        pieces.append(("ORF", orf_dna, "orf", "Codon-optimised in a single pass."))
        pieces.append((utr3.name, utr3.dna or "", "utr3", utr3.note))
        pieces.append((polya.name, polya.dna or "", "polya", polya.note))
        if linearizer:
            pieces.append((linearizer.name, linearizer.dna or "", "linearization",
                           linearizer.note))

        template = ""
        features: list[Feature] = []
        orf_offset = 0
        for piece_name, dna, kind, note in pieces:
            start = len(template)
            template += dna
            features.append(Feature(piece_name, start, len(template), kind, note))
            if piece_name == "ORF":
                orf_offset = start

        for seg_name, start, end, seg_kind in orf_segment_spans:
            features.append(
                Feature(seg_name, orf_offset + start, orf_offset + end, seg_kind)
            )
        stop_start = orf_offset + aa_cursor * 3
        features.append(
            Feature("stop", stop_start, orf_offset + len(orf_dna), "stop",
                    "Tandem stop" if spec.tandem_stop else "")
        )

        # The 3'-terminal G of a class III T7 promoter is transcript position +1.
        if promoter:
            promoter_feature = features[0]
            transcript_start = promoter_feature.end - 1
        else:
            transcript_start = 0
        polya_feature = next(f for f in features if f.kind == "polya")
        transcript_end = polya_feature.end

        features.sort(key=lambda f: (f.start, -len(f)))

        if placeholders:
            notes.append(
                "Placeholder parts in use ("
                + ", ".join(sorted(set(placeholders)))
                + "): replace before any in-vivo or GMP work."
            )

        return Construct(
            name=spec.name,
            template=template,
            features=features,
            orf_protein=translate(orf_dna, stop_at_stop=True).rstrip("*"),
            optimization=result,
            cassette=cassette,
            transcript_start=transcript_start,
            transcript_end=transcript_end,
            placeholders_used=sorted(set(placeholders)),
            design_notes=notes,
        )


def build_construct(
    name: str,
    antigens: list[Antigen],
    linker: str = "GGSGGGGSGG",
    spec: ConstructSpec | None = None,
    registry: PartRegistry | None = None,
    optimizer_config: OptimizerConfig | None = None,
) -> Construct:
    """Convenience wrapper: antigens in, annotated construct out."""
    spec = spec or ConstructSpec(name=name)
    spec.name = name
    builder = ConstructBuilder(
        registry=registry,
        optimizer=CodonOptimizer(optimizer_config) if optimizer_config else None,
    )
    return builder.build(spec, Cassette(antigens, linker=linker))


__all__ = [
    "Construct", "ConstructBuilder", "ConstructSpec", "Feature", "build_construct",
]
