"""Nine pJET1.2 maps for checking hand-built clones against.

Deliberately minimal. The vector is used exactly as supplied: same sequence,
same numbering, every feature left where it is. The only change is the gene
fragment dropped in at the blunt site, annotated so the design can be read off
the map.

What this does NOT do, on purpose:

* no rotation -- vector positions stay as they are in your file, so a
  coordinate quoted from one map means the same thing in the other;
* no reverse-complementing and no flipped-orientation files;
* no poly(A) -- these maps show the gene fragment as it arrives from the
  vendor, which is what gets pasted, so a hand-built map and this one can be
  compared base for base. Adding the tail is a later step (see
  CLONING_PROTOCOL.md).

The point is to be comparable, not clever.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from txv.genbank_io import GenBankFeature, read_genbank, write_genbank  # noqa: E402
from txv.ivt_unit import make_unit, unit_elements  # noqa: E402
from txv.plasmid import insert_blunt  # noqa: E402
from txv.pvax1_ag import PANEL, build_panel_construct  # noqa: E402
from txv.seqops import find_all, translate  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
VECTOR = ROOT / "backbone" / "pJET1.2.gb"
OUT = ROOT / "orders" / "pjet12_ivt_units" / "benchling_maps"

ECORV, ECORV_OFFSET = "GATATC", 3

_KEY = {"handle": "primer_bind", "promoter": "promoter", "utr5": "5'UTR",
        "utr3": "3'UTR", "kozak": "regulatory", "orf": "CDS",
        "signal_peptide": "sig_peptide", "skip_peptide": "misc_feature"}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    vector = read_genbank(VECTOR)

    sites = find_all(vector.sequence, ECORV)
    if len(sites) != 1:
        print(f"expected one EcoRV site, found {len(sites)}")
        return 1
    site = sites[0] + ECORV_OFFSET

    print(f"{vector.name}: {len(vector)} bp, inserting between "
          f"{site} and {site + 1}\n")

    problems = []
    for panel in PANEL:
        construct = build_panel_construct(panel.name)
        unit = make_unit(panel.name, construct.orf)

        build = insert_blunt(
            vector, unit.sequence, site, f"pJET1.2_{panel.name}",
            definition=f"pJET1.2 with the {panel.name} gene fragment at the "
                       f"Eco32I blunt site ({panel.route})",
        )
        for name, start, end, kind, note in unit_elements(unit, construct):
            build.record.features.append(GenBankFeature(
                _KEY.get(kind, "misc_feature"),
                build.insert_start + start, build.insert_start + end, 1,
                {"label": name, "note": note},
            ))
        build.record.features.sort(key=lambda f: (f.start, f.end))
        (OUT / f"{build.record.name}.gb").write_text(write_genbank(build.record))

        # -- the checks a hand-built map should also pass -------------------
        seq = build.record.sequence
        if len(seq) != len(vector) + len(unit):
            problems.append((panel.name, "length is not vector + fragment"))
        if seq[:site] != vector.sequence[:site]:
            problems.append((panel.name, "sequence before the site changed"))
        if seq[build.insert_end:] != vector.sequence[site:]:
            problems.append((panel.name, "sequence after the site changed"))
        if build.insert != unit.sequence:
            problems.append((panel.name, "the fragment is not intact"))
        if find_all(seq, ECORV):
            problems.append((panel.name, "the Eco32I site survived; the "
                             "fragment did not land inside it"))
        if len(find_all(seq, "AGATCT")) != 2:
            problems.append((panel.name, "expected exactly two BglII sites"))
        if len(find_all(seq, "GCCACCATGGCA")) != 1:
            problems.append((panel.name, "Kozak/start is not unique"))
        if not any("Eco47I" in d for d in build.disrupted):
            problems.append((panel.name, "eco47IR is not disrupted"))
        orf = next(f for f in build.record.features
                   if f.qualifiers.get("label") == "ORF")
        protein = translate(seq[orf.start:orf.end], stop_at_stop=True).rstrip("*")
        if protein != construct.orf_protein:
            problems.append((panel.name, "the ORF does not translate as designed"))

        print(f"  {panel.name}  {len(seq):>5} bp = {len(vector)} + "
              f"{len(unit):<4}  insert {build.insert_start + 1}-"
              f"{build.insert_end}  ORF {len(protein)} aa  "
              f"disrupted: {', '.join(build.disrupted)}")

    print(f"\nwrote {len(PANEL)} maps to "
          f"{OUT.relative_to(ROOT)}/")
    if problems:
        print("\nPROBLEMS:")
        for name, problem in problems:
            print(f"  [{name}] {problem}")
        return 1
    print("\nvector sequence either side of the insert is byte-identical to "
          "the parent; Eco32I destroyed; two BglII sites; ORFs translate as "
          "designed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
