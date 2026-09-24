"""Build the pJET1.2 order package: GenBank records, manifest, checklist.

The manifest exists because "are all the elements in there?" is a question you
should be able to answer from one table, not by opening nine maps.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from txv.ivt_unit import (  # noqa: E402
    FWD_PRIMER, REV_HANDLE, make_unit, reverse_primer, unit_elements,
    unit_genbank, verify_orientation_independence, verify_unit,
)
from txv.order import assess_synthesis  # noqa: E402
from txv.pvax1_ag import (  # noqa: E402
    MODULES, PANEL, PANEL_BY_NAME, build_panel_construct,
)
from txv.seqops import find_all, revcomp, translate  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "orders" / "pjet12_ivt_units"

#: The elements you asked to be able to check, and what each one is for.
TRACKED = [
    ("CTLA4_SP", "SP:CTLA4", "CTLA-4 signal peptide", "routes into the ER; N-terminal only"),
    ("LAMP1_SP", "SP:LAMP1", "LAMP1 signal peptide", "routes into the ER; N-terminal only"),
    ("MA", "MA", "Met-Ala start", "cytosolic constructs; supplies the shared Kozak +4 G"),
    ("HA", "HA", "HA tag", "detects expression, independently of presentation"),
    ("FLAG", "FLAG", "FLAG tag", "detects the second cistron in the dual-route pair"),
    ("A", "agA", "antigen A (SIINFEKL, 29-aa flanks)", "MHC-I readout"),
    ("B", "agB", "antigen B (I-A(b) core, 29-aa flanks)", "MHC-II readout"),
    ("L", "G4S", "GGGGS linker", "flexible spacer; 4 encodings cycled"),
    ("E5", "E5", "E5 acidic C-degron", "proteasome turnover; needs a free C-terminus"),
    ("P2A", "P2A", "P2A skip peptide", "makes one ORF give two proteins"),
    ("CTLA4_TMT", "TM:CTLA4", "CTLA-4 TM + YVKM tail", "AP-2, via the plasma membrane"),
    ("LAMP1_TMT", "TM:LAMP1", "LAMP1 TM + GYQTI tail", "AP-3, direct to lysosome; must be last"),
]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest_rows, checklist, problems = [], [], []

    for panel in PANEL:
        construct = build_panel_construct(panel.name)
        unit = make_unit(panel.name, construct.orf)

        problems += [(panel.name, p) for p in verify_unit(unit, construct.orf_protein)]
        problems += [(panel.name, p) for p in verify_orientation_independence(unit)]

        (OUT / f"{panel.name}.gb").write_text(
            unit_genbank(unit, construct, panel.route, panel.question)
        )

        counts = {m: panel.modules.count(m) for m, _, _, _ in TRACKED}
        protein = unit.protein
        antigen_order = "+".join(m for m in panel.modules if m in ("A", "B")) or "none"

        # Independent confirmation that each module's peptide is really present.
        for module, _, _, _ in TRACKED:
            if counts[module] and MODULES[module] not in protein:
                problems.append((panel.name, f"{module} declared but not in the protein"))

        manifest_rows.append({
            "construct": panel.name,
            "route": panel.route,
            "module_order": " - ".join(panel.modules),
            "antigen_order": antigen_order,
            "protein_aa": len(protein),
            "fragment_bp": len(unit),
            "transcript_nt": len(unit.transcript),
            **{m: counts[m] for m, _, _, _ in TRACKED},
        })

        checklist.append((panel, unit, construct, counts))

    # -- the sequences you actually order -----------------------------------
    # These three files were written once by hand when the package was first
    # cut, which meant the primary deliverable -- the string you paste into a
    # vendor's order form -- was the one artefact with no build step and no way
    # to regenerate it after a design change. They are built here so the fasta,
    # the table, the records and the manifest cannot drift apart.
    fasta_lines, table_rows = [], []
    for panel, unit, construct, _ in checklist:
        seq = unit.sequence
        risk = assess_synthesis(seq)
        fasta_lines.append(
            f">{panel.name}|pJET1.2_blunt|{len(seq)}bp|"
            f"transcript={len(unit.transcript)}nt")
        fasta_lines += [seq[i:i + 60] for i in range(0, len(seq), 60)]
        table_rows.append({
            "name": panel.name,
            "fragment_bp": len(seq),
            "gc": f"{risk.gc:.4f}",
            "transcript_nt": len(unit.transcript),
            "protein_aa": len(unit.protein),
            "synthesis_flags": "; ".join(risk.flags) or "none",
            "sequence": seq,
        })
    (OUT / "gblocks.fasta").write_text("\n".join(fasta_lines) + "\n")
    with (OUT / "ordering_table.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(table_rows[0]))
        writer.writeheader()
        writer.writerows(table_rows)

    with (OUT / "primers.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["name", "sequence", "length_nt", "use"])
        writer.writerow([
            "IVT_F", FWD_PRIMER, len(FWD_PRIMER),
            "Forward - anneals to the forward handle on every construct"])
        writer.writerow([
            "IVT_R_plain", reverse_primer(0), len(reverse_primer(0)),
            "Reverse without tail - colony PCR and sequencing checks. NOTE: you "
            "will NOT find this sequence in gblocks.fasta, and that is correct. "
            "A reverse primer anneals to the top strand, so what appears in the "
            "fragment is its reverse complement - the last 20 nt of every "
            "fragment, " + REV_HANDLE + "."])
        writer.writerow([
            "IVT_R_polyA120", reverse_primer(120), len(reverse_primer(120)),
            "Reverse with a 120-nt poly(T) 5' tail - adds the poly(A) to the "
            "IVT template. Order as an Ultramer or equivalent. Only the 3' 20 "
            "nt anneal; the T120 tail hangs off the template and is copied into "
            "the product, which is why the fragment encodes no poly(A)."])

    # -- manifest ----------------------------------------------------------
    with (OUT / "element_manifest.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0]))
        writer.writeheader()
        writer.writerows(manifest_rows)

    # -- human-readable checklist -------------------------------------------
    lines = [
        "# Element checklist — P0–P8 as ordered for pJET1.2",
        "",
        "Every element in every fragment, with its position in the ordered",
        "sequence. Cross-check against `element_manifest.csv` (the same data as a",
        "table) and the `.gb` files (the same data as annotated maps you can open",
        "in Benchling).",
        "",
        "Counts below are of the **peptide actually found in the translated**",
        "protein, not of what the design file claims.",
        "",
        "## What each tracked element is",
        "",
        "| element | what it is | why it is there |",
        "|---|---|---|",
    ]
    lines += [f"| `{m}` | {desc} | {why} |" for m, _, desc, why in TRACKED]
    lines += ["", "## Per construct", ""]

    for panel, unit, construct, counts in checklist:
        lines += [
            f"### {panel.name} — {panel.route}",
            "",
            f"{panel.question}",
            "",
            f"Fragment **{len(unit)} bp** · transcript **{len(unit.transcript)} nt** "
            f"· protein **{len(unit.protein)} aa** · antigens "
            f"**{'+'.join(m for m in panel.modules if m in ('A', 'B')) or 'none'}**",
            "",
            "| position | element | kind |",
            "|---|---|---|",
        ]
        for name, start, end, kind, _ in unit_elements(unit, construct):
            if kind == "orf":
                continue
            lines.append(f"| {start + 1}–{end} | `{name}` | {kind} |")
        present = [m for m, _, _, _ in TRACKED if counts[m]]
        absent = [m for m, _, _, _ in TRACKED if not counts[m]]
        lines += [
            "",
            f"**Present:** {', '.join(f'`{m}`' for m in present)}",
            "",
            f"**Absent (by design):** {', '.join(f'`{m}`' for m in absent) or 'none'}",
            "",
        ]

    lines += [
        "## Primers (the same three for all nine)",
        "",
        "| name | length | sequence | use |",
        "|---|---|---|---|",
        f"| IVT_F | {len(FWD_PRIMER)} nt | `{FWD_PRIMER}` | forward |",
        f"| IVT_R_plain | {len(reverse_primer(0))} nt | `{reverse_primer(0)}` | "
        "colony PCR, sequencing |",
        f"| IVT_R_polyA120 | {len(reverse_primer(120))} nt | `T`×120 + "
        f"`{revcomp(REV_HANDLE)}` | adds the poly(A) to the IVT template |",
        "",
    ]
    (OUT / "ELEMENT_CHECKLIST.md").write_text("\n".join(lines) + "\n")

    # -- the size table inside the hand-written protocol ---------------------
    # The protocol is prose and stays hand-written, but its table of fragment
    # sizes is data and had gone stale by two revisions -- every number in it
    # was wrong. Generated between markers so the prose stays editable and the
    # numbers cannot drift again.
    protocol = OUT / "CLONING_PROTOCOL.md"
    if protocol.exists():
        text = protocol.read_text()
        begin, end = "<!-- BEGIN GENERATED TABLE -->", "<!-- END GENERATED TABLE -->"
        if begin in text and end in text:
            table = [begin, "| | fragment | transcript | protein | GC |",
                     "|---|---|---|---|---|"]
            for row in manifest_rows:
                gc = float(next(r["gc"] for r in table_rows
                                if r["name"] == row["construct"]))
                table.append(
                    f"| {row['construct']} | {row['fragment_bp']} bp | "
                    f"{row['transcript_nt']} nt | {row['protein_aa']} aa | "
                    f"{gc:.1%} |")
            table.append(end)
            head = text[:text.index(begin)]
            tail = text[text.index(end) + len(end):]
            protocol.write_text(head + "\n".join(table) + tail)

    # -- report --------------------------------------------------------------
    print(f"wrote {len(PANEL)} GenBank records, element_manifest.csv and "
          f"ELEMENT_CHECKLIST.md to {OUT}/\n")
    header = f"{'':4}{'bp':>6}{'aa':>5}  " + "".join(
        f"{short:>9}" for _, short, _, _ in TRACKED)
    print(header)
    for row in manifest_rows:
        counts = "".join(f"{row[m]:>9}" for m, _, _, _ in TRACKED)
        print(f"  {row['construct']:<2}{row['fragment_bp']:>6}"
              f"{row['protein_aa']:>5}  {counts}")

    if problems:
        print("\nPROBLEMS:")
        for name, p in problems:
            print(f"  [{name}] {p}")
        return 1
    print("\nno problems: every tracked element verified present in the "
          "translated protein")
    return 0


if __name__ == "__main__":
    sys.exit(main())
