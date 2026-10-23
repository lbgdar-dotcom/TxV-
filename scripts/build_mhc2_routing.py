"""Emit the MHC II routing ORFs as GenBank files for Benchling import.

Start codon to stop codon only -- no UTRs, no promoter, no vector. The
constructs are inserted by hand.

Shared modules keep the canonical encodings from the P0-P8 panel, so a module
that appears in both studies is the same DNA in both and a difference between
them cannot be codon usage.

Constructs whose modules are not all resolved are skipped with a note rather
than emitted from a guessed sequence.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from txv.genbank_io import GenBankFeature, GenBankRecord, write_genbank  # noqa: E402
from txv.mhc2_routing import EXTERNAL, PANEL, module_protein  # noqa: E402
from txv.pvax1_ag import canonical_encodings, linker_variants  # noqa: E402

SPECIES = "mouse"
OUT = ROOT / "orders" / "mhc2_routing"
STOP = "TGA"


def build(construct, canon, linkers):
    """DNA plus one feature per module, in ORF coordinates.

    Linkers cycle through the panel's four encodings for the same GGGGS, which
    is how the P0-P8 panel avoids turning a repeated linker into a long direct
    repeat that recombines during plasmid propagation.
    """
    dna, feats, at, li = [], [], 0, 0
    feats.append(GenBankFeature("CDS", 0, 0, 1, {"label": construct.name,
                                                 "note": construct.route}))
    for mod in construct.modules:
        if mod == "L":
            seg = linkers[li % len(linkers)]
            li += 1
        else:
            seg = canon[mod]
        feats.append(GenBankFeature(
            "misc_feature", at, at + len(seg), 1,
            {"label": mod, "note": module_protein(mod)}))
        dna.append(seg)
        at += len(seg)
    dna.append(STOP)
    feats.append(GenBankFeature("misc_feature", at, at + 3, 1,
                                {"label": "stop"}))
    seq = "".join(dna)
    feats[0] = GenBankFeature("CDS", 0, len(seq), 1,
                              {"label": construct.name,
                               "note": construct.route,
                               "translation": module_protein_all(construct)})
    return GenBankRecord(name=construct.name, sequence=seq, features=feats,
                         is_circular=False,
                         definition=f"{construct.name}: {construct.route}")


def module_protein_all(construct):
    return "".join(module_protein(m) for m in construct.modules)


def main() -> int:
    canon = canonical_encodings(species=SPECIES)
    linkers = linker_variants(species=SPECIES)
    OUT.mkdir(parents=True, exist_ok=True)
    built, skipped = [], []

    for c in PANEL:
        if not c.buildable:
            skipped.append(c)
            continue
        rec = build(c, canon, linkers)
        (OUT / f"{c.name}.gb").write_text(write_genbank(rec))
        prot = module_protein_all(c)
        built.append((c, len(rec.sequence), len(prot)))
        print(f"  {c.name}  {len(rec.sequence):4d} bp  {len(prot):3d} aa  "
              f"{c.route}")

    if skipped:
        print("\nnot built -- sequence required:")
        for c in skipped:
            missing = [m for m in c.external if not EXTERNAL[m].resolved]
            print(f"  {c.name}: {c.route}")
            for m in missing:
                e = EXTERNAL[m]
                print(f"      {m}: {e.protein} ({e.organism})")
                print(f"         region: {e.region}")

    lines = ["# MHC II routing panel -- ORFs for Benchling", "",
             "Start codon to stop codon. No UTRs, no promoter, no vector.",
             "Codon usage: mouse. Shared modules carry the same canonical",
             "encodings as the P0-P8 panel, so the two studies are directly",
             "comparable.", "",
             "| construct | bp | aa | route | tag |", "|---|---|---|---|---|"]
    for c, nt, aa in built:
        lines.append(f"| {c.name} | {nt} | {aa} | {c.route} | {c.tag_position} |")
    if skipped:
        lines += ["", "## Not built", "",
                  "These need a sequence this package does not hold; it was",
                  "left empty rather than recalled, because a gBlock ordered",
                  "from a remembered sequence is an expensive way to find out",
                  "it was wrong.", ""]
        for c in skipped:
            lines.append(f"**{c.name}** -- {c.route}")
            for m in [m for m in c.external if not EXTERNAL[m].resolved]:
                e = EXTERNAL[m]
                lines += [f"- `{m}`: {e.protein} ({e.organism})",
                          f"  - region: {e.region}",
                          f"  - why: {e.why}"]
            if c.notes:
                lines.append(f"- note: {c.notes}")
            lines.append("")
    (OUT / "README.md").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {len(built)} GenBank file(s) and README.md to {OUT}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
