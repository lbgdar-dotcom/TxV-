"""End-to-end example: neoantigen + TAA cassette -> QC -> LNP -> Benchling.

Run from the repository root::

    PYTHONPATH=src python examples/design_multiantigen.py
"""

from pathlib import Path

from txv.benchling import BenchlingConfig, make_registrar
from txv.epitopes import AnchorMotifScorer
from txv.lnp import formulate
from txv.pipeline import design_construct, load_antigens, write_outputs

OUT = Path("out")


def main() -> None:
    antigens = load_antigens("examples/antigens_example.csv")
    print(f"{len(antigens)} antigens loaded\n")

    # A real design would pass a NetMHCpan/MHCflurry-backed scorer here; the
    # anchor proxy only ranks options relative to each other.
    report = design_construct(
        "TXV-MULTI-001",
        antigens,
        scorer=AnchorMotifScorer(),
        allow_placeholders=True,   # the built-in UTRs are stand-ins
    )
    print(report.to_text())

    written = write_outputs(report, OUT)
    print(f"\nwrote {len(written)} design file(s) to {OUT}/")

    formulation = formulate(
        sequence=report.construct.transcript,
        mrna_micrograms=100.0,
        np_ratio=6.0,
        composition="sm102_standard",
        construct_name=report.construct.name,
        lipid_stocks_mg_per_ml={
            "SM-102": 10.0, "DSPC": 10.0,
            "cholesterol": 20.0, "DMG-PEG2000": 10.0,
        },
    )
    print("\n" + formulation.to_text())

    # Dry run: writes the exact payload that would be sent, and sends nothing.
    registrar = make_registrar(
        BenchlingConfig.from_env(folder_id="lib_example"),
        dry_run=True,
        out_dir=OUT,
    )
    print("\nbenchling:", registrar.register(report.construct, report.qc))


if __name__ == "__main__":
    main()
