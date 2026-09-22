"""Command-line entry point: ``txv <command>``."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .codon import OptimizerConfig
from .constructs import ConstructSpec
from .epitopes import AnchorMotifScorer
from .lnp import COMPOSITIONS, dose_series, formulate
from .parts import default_registry, load_parts_file
from .pipeline import (
    DEFAULT_LINKER_CANDIDATES,
    design_construct,
    load_antigens,
    write_outputs,
)
from .qc import QCThresholds, run_qc


def _registry(args):
    registry = default_registry()
    if getattr(args, "parts", None):
        registry = load_parts_file(args.parts, registry)
    return registry


def cmd_design(args) -> int:
    antigens = load_antigens(args.antigens)
    spec = ConstructSpec(
        name=args.name,
        utr5=args.utr5,
        utr3=args.utr3,
        signal_peptide=None if args.no_signal_peptide else args.signal_peptide,
        trafficking=None if args.no_trafficking else args.trafficking,
        polya=args.polya,
    )
    report = design_construct(
        name=args.name,
        antigens=antigens,
        spec=spec,
        registry=_registry(args),
        optimizer_config=OptimizerConfig(
            target_gc=args.target_gc,
            w_uridine=args.uridine_weight,
            beam_width=args.beam_width,
        ),
        thresholds=QCThresholds(),
        scorer=AnchorMotifScorer(),
        linker=args.linker,
        linker_candidates=tuple(args.linker_candidates or DEFAULT_LINKER_CANDIDATES),
        reorder=not args.no_reorder,
        allow_placeholders=not args.strict,
    )
    print(report.to_text())

    if args.out:
        written = write_outputs(report, args.out)
        print(f"\nwrote {len(written)} file(s) to {args.out}/")

    if args.register or args.dry_run_register:
        from .benchling import BenchlingConfig, make_registrar

        config = BenchlingConfig.from_env(
            folder_id=args.folder_id, schema_id=args.schema_id,
            registry_id=args.registry_id,
        )
        registrar = make_registrar(
            config,
            dry_run=True if args.dry_run_register else None,
            out_dir=args.out or "benchling_dry_run",
        )
        try:
            record = registrar.register(report.construct, report.qc)
        except ValueError as exc:
            print(f"\nregistration refused: {exc}", file=sys.stderr)
            return 2
        print(f"\nbenchling: {json.dumps(record)}")

    return 0 if report.qc.passed else 1


def cmd_formulate(args) -> int:
    if args.fasta:
        text = Path(args.fasta).read_text()
        sequence = "".join(
            line.strip() for line in text.splitlines() if not line.startswith(">")
        )
    else:
        sequence = args.sequence
    stocks = json.loads(args.stocks) if args.stocks else None
    result = formulate(
        sequence=sequence,
        mrna_micrograms=args.mrna_ug,
        np_ratio=args.np,
        composition=args.composition,
        construct_name=args.name,
        lipid_stocks_mg_per_ml=stocks,
        flow_rate_ratio=args.frr,
    )
    print(result.to_text())
    if args.doses:
        weights = [float(w) for w in args.doses.split(",")]
        print("\n  dose plan")
        for row in dose_series(args.mg_per_kg, weights, args.concentration):
            print(
                f"    {row['body_weight_kg']:>8.3f} kg -> "
                f"{row['dose_ug']:>8.1f} ug in {row['volume_ul']:>8.1f} uL"
            )
    if args.json:
        Path(args.json).write_text(json.dumps(result.to_dict(), indent=2))
        print(f"\nwrote {args.json}")
    return 0


def cmd_parts(args) -> int:
    registry = _registry(args)
    for name in registry.names():
        part = registry.get(name)
        body = part.dna if part.dna else part.protein
        shown = body if len(body) <= 48 else body[:45] + "..."
        flag = "  [PLACEHOLDER]" if part.is_placeholder else ""
        print(f"{name:<24}{part.kind.value:<16}{part.provenance.value:<12}{shown}{flag}")
    placeholders = registry.placeholders()
    if placeholders:
        print(
            f"\n{len(placeholders)} placeholder part(s) are not real sequences; "
            "supply your own with --parts before any in-vivo work."
        )
    return 0


def cmd_panel(args) -> int:
    """Build the audited P0-P8 routing panel as IVT mRNA templates."""
    from .export import to_fasta, to_genbank
    from .pvax1_ag import OPEN_DECISIONS, PANEL, build_panel_construct

    names = args.only or [p.name for p in PANEL]
    registry = _registry(args)
    failures = 0
    for name in names:
        spec = ConstructSpec(name=name, utr5=args.utr5, utr3=args.utr3,
                             polya=args.polya)
        construct = build_panel_construct(name, spec=spec, registry=registry)
        report = run_qc(construct, scorer=AnchorMotifScorer(),
                        allow_placeholders=not args.strict)
        panel = next(p for p in PANEL if p.name == name)
        print(
            f"{name}  {len(construct.orf_protein):>4} aa  "
            f"{construct.transcript_length:>5} nt  "
            f"CAI {construct.optimization.cai:.3f}  "
            f"U {construct.transcript.count('T') / construct.transcript_length:.1%}  "
            f"{'PASS' if report.passed else 'FAIL'}  {panel.route}"
        )
        if not report.passed:
            failures += 1
            print(report.to_text(show_pass=False))
        if args.out:
            out = Path(args.out)
            out.mkdir(parents=True, exist_ok=True)
            (out / f"{name}.gb").write_text(to_genbank(construct))
            (out / f"{name}.mrna.fasta").write_text(to_fasta(construct, "mrna"))
            (out / f"{name}.protein.fasta").write_text(to_fasta(construct, "protein"))

    if args.out:
        print(f"\nwrote {len(names) * 3} file(s) to {args.out}/")
    print("\nOpen decisions carried over from the audit:")
    for decision in OPEN_DECISIONS:
        print(f"  - {decision}")
    return 1 if failures else 0


def cmd_order(args) -> int:
    """Produce synthesis-ready fragments for the panel or a designed construct."""
    from .order import (
        FragmentMode, PolyAMode, fragment_genbank, order_fasta, order_table,
        synthesis_fragment,
    )
    from .pvax1_ag import PANEL, build_panel_construct

    names = args.only or [p.name for p in PANEL]
    registry = _registry(args)
    fragments = []
    records: dict[str, str] = {}
    for name in names:
        spec = ConstructSpec(name=name, utr5=args.utr5, utr3=args.utr3,
                             polya=args.polya)
        construct = build_panel_construct(name, spec=spec, registry=registry)
        fragment = synthesis_fragment(construct, FragmentMode(args.mode),
                                      PolyAMode(args.polya_mode))
        fragments.append(fragment)
        records[name] = fragment_genbank(fragment, construct)

    flagged = 0
    for fragment in fragments:
        status = "ok" if fragment.risk.ok else "CHECK"
        print(f"{fragment.name:<6}{fragment.length:>6} bp  "
              f"GC {fragment.risk.gc:>6.1%}  {status}")
        for flag in fragment.risk.flags:
            flagged += 1
            print(f"         ! {flag}")
    if fragments:
        print()
        for note in fragments[0].notes:
            print(f"  note: {note}")

    if args.out:
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        (out / "ordering_table.csv").write_text(order_table(fragments))
        (out / "gblocks.fasta").write_text(order_fasta(fragments))
        for name, record in records.items():
            (out / f"{name}_gblock.gb").write_text(record)
        print(f"\nwrote ordering_table.csv, gblocks.fasta and "
              f"{len(records)} GenBank record(s) to {out}/")
    return 1 if flagged else 0


def cmd_series(args) -> int:
    """Build a matched variant series for one of the two experimental axes."""
    from .variants import multiplex_ladder, series_outputs, uridine_ladder

    antigens = load_antigens(args.antigens)
    if args.axis == "uridine":
        series = uridine_ladder(antigens, route=args.route)
    else:
        sizes = tuple(int(s) for s in args.sizes.split(",")) if args.sizes else (1, 2, 4, 6, 8)
        series = multiplex_ladder(antigens, sizes=sizes, route=args.route)
    print(series.to_text())
    failures = [v.name for v in series.variants if not v.qc.passed]
    if failures:
        print(f"\nQC failed for: {', '.join(failures)}")
    if args.out:
        written = series_outputs(series, args.out)
        print(f"\nwrote {len(written)} file(s) to {args.out}/")
    return 1 if failures else 0


def cmd_compositions(args) -> int:
    for name, comp in COMPOSITIONS.items():
        ratios = " : ".join(f"{k} {v}" for k, v in comp.molar_ratios.items())
        print(f"{name}\n  {ratios}\n  {comp.note}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="txv",
        description="Design IVT mRNA constructs for multi-antigen cancer vaccines.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    design = sub.add_parser("design", help="design a construct from an antigen list")
    design.add_argument("antigens", help="antigen table (.json, .csv or .tsv)")
    design.add_argument("--name", required=True, help="construct name")
    design.add_argument("--out", help="output directory")
    design.add_argument("--parts", help="JSON parts file overriding the built-ins")
    design.add_argument("--utr5", default="UTR5_hAg")
    design.add_argument("--utr3", default="UTR3_AES_mtRNR1")
    design.add_argument("--signal-peptide", default="SP_tPA")
    design.add_argument("--no-signal-peptide", action="store_true")
    design.add_argument("--trafficking", default="MITD_placeholder")
    design.add_argument("--no-trafficking", action="store_true")
    design.add_argument("--polya", default="polyA_120")
    design.add_argument("--linker", help="force a linker instead of selecting one")
    design.add_argument("--linker-candidates", nargs="*")
    design.add_argument("--no-reorder", action="store_true")
    design.add_argument("--target-gc", type=float, default=0.60)
    design.add_argument("--uridine-weight", type=float, default=0.45)
    design.add_argument("--beam-width", type=int, default=40)
    design.add_argument(
        "--strict", action="store_true",
        help="fail QC when placeholder parts are still in the construct",
    )
    design.add_argument("--register", action="store_true", help="register in Benchling")
    design.add_argument(
        "--dry-run-register", action="store_true",
        help="write the Benchling payload without sending it",
    )
    design.add_argument("--folder-id")
    design.add_argument("--schema-id")
    design.add_argument("--registry-id")
    design.set_defaults(func=cmd_design)

    form = sub.add_parser("formulate", help="LNP formulation arithmetic")
    source = form.add_mutually_exclusive_group(required=True)
    source.add_argument("--sequence")
    source.add_argument("--fasta")
    form.add_argument("--name", default="construct")
    form.add_argument("--mrna-ug", type=float, required=True)
    form.add_argument("--np", type=float, default=6.0, help="N/P ratio")
    form.add_argument("--composition", default="sm102_standard",
                      choices=sorted(COMPOSITIONS))
    form.add_argument("--stocks", help='JSON, e.g. \'{"SM-102": 10, "DSPC": 10}\'')
    form.add_argument("--frr", type=float, default=3.0,
                      help="aqueous:organic flow rate ratio")
    form.add_argument("--doses", help="comma-separated body weights in kg")
    form.add_argument("--mg-per-kg", type=float, default=0.5)
    form.add_argument("--concentration", type=float, default=0.1,
                      help="drug product strength, ug/uL")
    form.add_argument("--json", help="write the result as JSON to this path")
    form.set_defaults(func=cmd_formulate)

    parts = sub.add_parser("parts", help="list available construct parts")
    parts.add_argument("--parts", help="JSON parts file overriding the built-ins")
    parts.set_defaults(func=cmd_parts)

    panel = sub.add_parser(
        "panel", help="build the audited P0-P8 routing panel as IVT mRNA"
    )
    panel.add_argument("--only", nargs="*", help="subset, e.g. --only P3 P7")
    panel.add_argument("--out", help="output directory")
    panel.add_argument("--parts", help="JSON parts file overriding the built-ins")
    panel.add_argument("--utr5", default="UTR5_hAg")
    panel.add_argument("--utr3", default="UTR3_AES_mtRNR1")
    panel.add_argument("--polya", default="polyA_120")
    panel.add_argument("--strict", action="store_true")
    panel.set_defaults(func=cmd_panel)

    order = sub.add_parser("order", help="synthesis-ready fragments to order")
    order.add_argument("--only", nargs="*", help="subset, e.g. --only P3 P7")
    order.add_argument("--mode", default="orf", choices=["orf", "cassette"],
                       help="orf: drop into the existing backbone rails. "
                            "cassette: also replace the backbone UTRs")
    order.add_argument("--polya-mode", default="pcr_added",
                       choices=["pcr_added", "encoded"],
                       help="pcr_added omits the A-tract from the order")
    order.add_argument("--out", help="output directory")
    order.add_argument("--parts", help="JSON parts file overriding the built-ins")
    order.add_argument("--utr5", default="UTR5_hAg")
    order.add_argument("--utr3", default="UTR3_AES_mtRNR1")
    order.add_argument("--polya", default="polyA_120")
    order.set_defaults(func=cmd_order)

    series = sub.add_parser(
        "series", help="matched variant series for the experimental axes"
    )
    series.add_argument("antigens", help="antigen table (.json, .csv or .tsv)")
    series.add_argument("--axis", default="multiplex",
                        choices=["multiplex", "uridine"])
    series.add_argument("--route", default="ctla4",
                        choices=["ctla4", "lamp1", "cytosolic"])
    series.add_argument("--sizes", help="multiplex only, e.g. 1,2,4,6,8")
    series.add_argument("--out", help="output directory")
    series.set_defaults(func=cmd_series)

    comps = sub.add_parser("compositions", help="list LNP lipid compositions")
    comps.set_defaults(func=cmd_compositions)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
