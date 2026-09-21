"""End-to-end design pipeline: antigen list in, registered construct out.

The stages, in order, and why that order:

1. **Deduplicate.** A neoantigen call set routinely contains the same peptide
   twice under different transcript IDs. Registering it twice wastes cassette
   length, which is the scarcest resource in a multi-antigen design.
2. **Choose the linker** against the junction scorer, over the antigen set as
   given. Linker choice dominates junctional risk, so it is settled first.
3. **Order the beads** with the chosen linker. Ordering is a smaller effect
   than linker choice, which is why it comes second.
4. **Assemble and codon-optimise** in a single pass over the whole ORF.
5. **QC.** Nothing is registered before this.
6. **Export and register**, with registration refused on a QC failure.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

from .codon import CodonOptimizer, OptimizerConfig
from .constructs import Construct, ConstructBuilder, ConstructSpec
from .epitopes import (
    CLASS_I_LENGTHS,
    Antigen,
    Cassette,
    JunctionScorer,
    choose_linker,
    dedupe_antigens,
    optimize_order,
)
from .export import to_fasta, to_genbank
from .parts import PartRegistry, default_registry
from .qc import QCReport, QCThresholds, run_qc

DEFAULT_LINKER_CANDIDATES = ("GGSGGGGSGG", "GGGGS", "GPGPG", "AAY", "RAKR")


@dataclass
class DesignReport:
    construct: Construct
    qc: QCReport
    dropped_antigens: list[str] = field(default_factory=list)
    linker_table: dict[str, float] = field(default_factory=dict)
    chosen_linker: str = ""
    order_before: float = 0.0
    order_after: float = 0.0
    reordered: bool = False

    def to_text(self) -> str:
        lines = [f"Design: {self.construct.name}"]
        if self.dropped_antigens:
            lines.append(f"  dropped {len(self.dropped_antigens)} redundant antigen(s):")
            lines += [f"    - {d}" for d in self.dropped_antigens]
        if self.linker_table:
            ranked = sorted(self.linker_table.items(), key=lambda kv: kv[1])
            lines.append("  linker junctional risk (lower is better):")
            lines += [
                f"    {name:<14}{risk:8.2f}" + ("   <- chosen" if name == self.chosen_linker else "")
                for name, risk in ranked
            ]
        if self.reordered:
            lines.append(
                f"  reordered beads: junctional risk {self.order_before:.2f} "
                f"-> {self.order_after:.2f}"
            )
        lines.append(
            "  final order: "
            + " | ".join(a.name for a in self.construct.cassette.antigens)
        )
        summary = self.construct.summary()
        lines.append(
            f"  transcript {summary['transcript_nt']} nt, ORF {summary['orf_nt']} nt, "
            f"{summary['protein_aa']} aa, GC {summary['gc']:.1%}, "
            f"U {summary['uridine_fraction']:.1%}, CAI {summary['cai']:.3f}"
        )
        lines.append("")
        lines.append(self.qc.to_text(show_pass=False))
        return "\n".join(lines)

    def to_dict(self) -> dict[str, object]:
        return {
            "design": {
                "dropped_antigens": self.dropped_antigens,
                "linker_table": self.linker_table,
                "chosen_linker": self.chosen_linker,
                "reordered": self.reordered,
                "junction_risk_before": self.order_before,
                "junction_risk_after": self.order_after,
                "order": [a.name for a in self.construct.cassette.antigens],
            },
            "construct": self.construct.summary(),
            "qc": self.qc.to_dict(),
        }


def design_construct(
    name: str,
    antigens: Sequence[Antigen],
    spec: ConstructSpec | None = None,
    registry: PartRegistry | None = None,
    optimizer_config: OptimizerConfig | None = None,
    thresholds: QCThresholds | None = None,
    scorer: JunctionScorer | None = None,
    lengths: Sequence[int] = CLASS_I_LENGTHS,
    linker: str | None = None,
    linker_candidates: Sequence[str] = DEFAULT_LINKER_CANDIDATES,
    reorder: bool = True,
    fixed_first: bool = False,
    fixed_last: bool = False,
    dedupe: bool = True,
    allow_placeholders: bool = True,
) -> DesignReport:
    """Run the full design pipeline for one construct."""
    if not antigens:
        raise ValueError("no antigens supplied")

    dropped: list[str] = []
    working = list(antigens)
    if dedupe:
        working, dropped = dedupe_antigens(working)
    if not working:
        raise ValueError("all antigens were removed as redundant")

    cassette = Cassette(working, linker=linker or DEFAULT_LINKER_CANDIDATES[0])
    table: dict[str, float] = {}
    if linker is None:
        chosen, table = choose_linker(cassette, linker_candidates, scorer, lengths)
        cassette = Cassette(working, linker=chosen)
    else:
        chosen = linker

    before = after = 0.0
    reordered = False
    if reorder and len(working) >= 3:
        cassette, before, after = optimize_order(
            cassette, scorer, lengths,
            fixed_first=fixed_first, fixed_last=fixed_last,
        )
        reordered = after < before

    spec = spec or ConstructSpec(name=name)
    spec.name = name
    builder = ConstructBuilder(
        registry=registry or default_registry(),
        optimizer=CodonOptimizer(optimizer_config) if optimizer_config else None,
    )
    construct = builder.build(spec, cassette)
    report = run_qc(
        construct, thresholds=thresholds, scorer=scorer,
        allow_placeholders=allow_placeholders,
    )
    return DesignReport(
        construct=construct,
        qc=report,
        dropped_antigens=dropped,
        linker_table=table,
        chosen_linker=chosen,
        order_before=before,
        order_after=after,
        reordered=reordered,
    )


# ---------------------------------------------------------------------------
# Antigen input
# ---------------------------------------------------------------------------

def load_antigens(path: str | Path) -> list[Antigen]:
    """Read antigens from ``.json`` or ``.csv``/``.tsv``.

    CSV/TSV columns: ``name``, ``sequence``, and optionally ``kind``, ``hla``
    (semicolon-separated), ``mutation_offset``, ``note``.
    """
    path = Path(path)
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text())
        rows = payload["antigens"] if isinstance(payload, dict) else payload
        return [_antigen_from_mapping(r) for r in rows]
    delimiter = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","
    with path.open(newline="") as handle:
        return [
            _antigen_from_mapping(row)
            for row in csv.DictReader(handle, delimiter=delimiter)
            if row.get("name") and row.get("sequence")
        ]


def _antigen_from_mapping(row: dict) -> Antigen:
    hla = row.get("hla") or ""
    if isinstance(hla, str):
        hla_tuple = tuple(h.strip() for h in hla.replace(",", ";").split(";") if h.strip())
    else:
        hla_tuple = tuple(hla)
    offset = row.get("mutation_offset")
    return Antigen(
        name=str(row["name"]).strip(),
        sequence=str(row["sequence"]).strip(),
        kind=str(row.get("kind") or "neoepitope").strip(),
        hla=hla_tuple,
        mutation_offset=int(offset) if offset not in (None, "", "NA") else None,
        note=str(row.get("note") or ""),
    )


def write_outputs(report: DesignReport, out_dir: str | Path) -> dict[str, str]:
    """Write GenBank, FASTA (template/mRNA/protein), and the JSON report."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    construct = report.construct
    written: dict[str, str] = {}

    targets = {
        f"{construct.name}.gb": to_genbank(construct),
        f"{construct.name}.template.fasta": to_fasta(construct, "template"),
        f"{construct.name}.mrna.fasta": to_fasta(construct, "mrna"),
        f"{construct.name}.protein.fasta": to_fasta(construct, "protein"),
        f"{construct.name}.report.json": json.dumps(report.to_dict(), indent=2),
        f"{construct.name}.report.txt": report.to_text() + "\n",
    }
    for filename, content in targets.items():
        (out / filename).write_text(content)
        written[filename] = str(out / filename)
    return written


__all__ = [
    "DesignReport", "design_construct", "load_antigens", "write_outputs",
    "DEFAULT_LINKER_CANDIDATES",
]
