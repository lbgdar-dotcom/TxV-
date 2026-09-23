"""Build the nine finished pJET1.2 plasmid maps and verify them.

You buy pJET1.2 pre-linearised, order the nine fragments, and ligate. This
produces the maps of what you get, so the construct can be checked in SnapGene
or Benchling before anything is sequenced, and works out the three things you
actually need at the bench: whether the selection will work, whether the stock
sequencing primers cover the insert, and how to tell which way round it went in.

Blunt cloning has no orientation control, so both orientations are built. That
is not a hedge -- both are real products, roughly half the colonies each.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from txv.genbank_io import write_genbank  # noqa: E402
from txv.ivt_unit import (  # noqa: E402
    FWD_PRIMER, REV_HANDLE, make_unit, reverse_primer, simulate_pcr,
)
from txv.plasmid import insert_blunt  # noqa: E402
from txv.pvax1_ag import PANEL, build_panel_construct  # noqa: E402
from txv.seqops import find_all, revcomp  # noqa: E402
from txv.snapgene import read_primers, read_snapgene  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
VECTOR_GB = ROOT / "backbone" / "pJET1.2.gb"
OUT = ROOT / "orders" / "pjet12_ivt_units" / "plasmids"

#: Eco32I (EcoRV) recognition site; the kit ships the vector cut here.
ECORV = "GATATC"
#: EcoRV is blunt and cuts GAT^ATC, three bases into its site.
ECORV_OFFSET = 3
#: Excision enzyme. pJET1.2 carries a site either side of the cloning site.
BGLII, BGLII_OFFSET = "AGATCT", 1
#: A Sanger read is reliable for roughly this many bases past the primer.
SANGER_READ = 850


def cut_positions(sequence: str, site: str, offset: int) -> list[int]:
    return [i + offset for i in find_all(sequence, site)]


def digest_fragments(sequence: str, cuts: list[int]) -> list[tuple[int, int, int]]:
    """Circular-digest fragments as ``(start, end, length)``, in cut order.

    The wrap-around fragment is returned with ``end`` past the sequence length
    so a span can be tested against it without special-casing the origin.
    """
    cuts = sorted(cuts)
    if len(cuts) < 2:
        return [(0, len(sequence), len(sequence))]
    out = [(cuts[i], cuts[i + 1], cuts[i + 1] - cuts[i])
           for i in range(len(cuts) - 1)]
    out.append((cuts[-1], cuts[0] + len(sequence),
                len(sequence) - cuts[-1] + cuts[0]))
    return out


def fragment_carrying(fragments, start: int, end: int, length: int):
    """The digest fragment that contains the span ``start..end``.

    Taking the largest fragment instead is wrong and quietly so: in pJET1.2 the
    backbone piece is bigger than the insert-carrying one, so ``max`` reports
    the same size for every construct regardless of insert length -- which is
    what gave the bug away.
    """
    for lo, hi, size in fragments:
        for offset in (0, length):
            if lo <= start + offset and end + offset <= hi:
                return size
    return None


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    vector = read_snapgene(sys.argv[1]) if len(sys.argv) > 1 else None
    if vector is None:
        if not VECTOR_GB.exists():
            print(f"pass the pJET1.2 .dna file, or put a map at {VECTOR_GB}")
            return 1
        from txv.genbank_io import read_genbank
        vector = read_genbank(VECTOR_GB)
    else:
        VECTOR_GB.write_text(write_genbank(vector))
        print(f"wrote {VECTOR_GB.relative_to(ROOT)}")

    stock = read_primers(sys.argv[1]) if len(sys.argv) > 1 else {}
    seq_fwd = next((s for n, s in stock.items() if "forward" in n.lower()), None)
    seq_rev = next((s for n, s in stock.items() if "reverse" in n.lower()), None)

    sites = find_all(vector.sequence, ECORV)
    if len(sites) != 1:
        print(f"expected one EcoRV site in the vector, found {len(sites)}")
        return 1
    site = sites[0] + ECORV_OFFSET

    print(f"\n{vector.name}: {len(vector)} bp, EcoRV blunt cut at {site}")
    print(f"stock primers: {', '.join(stock) or 'none found'}\n")

    rows, problems = [], []
    for panel in PANEL:
        unit = make_unit(panel.name, build_panel_construct(panel.name).orf)

        for orientation, tag in ((1, "fwd"), (-1, "rev")):
            build = insert_blunt(
                vector, unit.sequence, site,
                f"pJET1.2_{panel.name}" + ("" if orientation == 1 else "_rev"),
                orientation=orientation,
                definition=f"pJET1.2 carrying the {panel.name} IVT transcription "
                           f"unit ({panel.route})",
            )
            (OUT / f"{build.record.name}.gb").write_text(
                write_genbank(build.record))

            if orientation != 1:
                continue

            seq = build.record.sequence
            # -- the selection has to be broken for the cloning to work -------
            if "Eco47I/T7" not in " ".join(build.disrupted):
                problems.append((panel.name, "eco47IR is not disrupted; "
                                 "positive selection would not work"))
            # -- and the resistance and origin have to survive ----------------
            for keep in ("AmpR", "ori"):
                feature = next((f for f in build.record.features
                                if f.qualifiers.get("label") == keep), None)
                if feature is None:
                    problems.append((panel.name, f"{keep} missing"))
                elif keep in build.disrupted:
                    problems.append((panel.name, f"{keep} disrupted"))

            # -- excision -----------------------------------------------------
            bgl = cut_positions(seq, BGLII, BGLII_OFFSET)
            fragments = digest_fragments(seq, bgl)
            insert_carrying = fragment_carrying(
                fragments, build.insert_start, build.insert_end, len(seq))
            if len(bgl) >= 2 and insert_carrying is None:
                problems.append((panel.name, "no BglII fragment contains the "
                                 "insert; it cannot be excised in one piece"))
            backbone_piece = max(f[2] for f in fragments) if len(bgl) > 1 else None

            # -- sequencing coverage -----------------------------------------
            covered = None
            if seq_fwd and seq_rev:
                f_hits = find_all(seq, seq_fwd, both_strands=True)
                r_hits = find_all(seq, revcomp(seq_rev), both_strands=True)
                if f_hits and r_hits:
                    f_end = f_hits[0] + len(seq_fwd)
                    reach_f = max(0, f_end + SANGER_READ - build.insert_start)
                    r_start = r_hits[0]
                    reach_r = max(0, build.insert_end - (r_start - SANGER_READ))
                    covered = reach_f + reach_r >= len(unit)
                    if not covered:
                        problems.append((
                            panel.name,
                            f"stock primers reach {reach_f}+{reach_r} nt of a "
                            f"{len(unit)} nt insert; an internal primer is needed"))

            # -- telling the orientation apart --------------------------------
            fwd_products = simulate_pcr(seq, FWD_PRIMER, seq_rev or REV_HANDLE)
            rev_build = insert_blunt(vector, unit.sequence, site, "tmp",
                                     orientation=-1)
            rev_products = simulate_pcr(rev_build.record.sequence,
                                        FWD_PRIMER, seq_rev or REV_HANDLE)

            rows.append({
                "construct": panel.name,
                "plasmid_bp": len(build),
                "insert_bp": len(unit),
                "insert_at": build.insert_start + 1,
                "disrupted": "; ".join(build.disrupted),
                "BglII_sites": len(bgl),
                "BglII_insert_fragment_bp": insert_carrying or "",
                "BglII_backbone_fragment_bp": backbone_piece or "",
                "sequencing_covers_insert": "yes" if covered else "no",
                "orientation_PCR_fwd_bp": ";".join(str(len(p)) for p in fwd_products) or "none",
                "orientation_PCR_rev_bp": ";".join(str(len(p)) for p in rev_products) or "none",
            })
            print(f"  {panel.name}  {len(build):>5} bp  insert {len(unit):>4} bp at "
                  f"{build.insert_start + 1}  BglII: {insert_carrying}+{backbone_piece} bp  "
                  f"seq covers: {'yes' if covered else 'NO'}  "
                  f"orientation PCR fwd/rev: "
                  f"{[len(p) for p in fwd_products]}/{[len(p) for p in rev_products]}")

    # -- the protocol's verification table, generated for the same reason the
    # -- size table is: numbers next to prose go stale silently.
    protocol = OUT.parent / "CLONING_PROTOCOL.md"
    begin = "<!-- BEGIN GENERATED VERIFICATION TABLE -->"
    end = "<!-- END GENERATED VERIFICATION TABLE -->"
    if protocol.exists() and begin in protocol.read_text():
        text = protocol.read_text()
        table = [begin,
                 "| | plasmid | BglII insert band | BglII backbone band | "
                 "orientation PCR (designed) | (flipped) |",
                 "|---|---|---|---|---|---|"]
        for row in rows:
            table.append(
                f"| {row['construct']} | {row['plasmid_bp']} bp | "
                f"{row['BglII_insert_fragment_bp']} bp | "
                f"{row['BglII_backbone_fragment_bp']} bp | "
                f"{row['orientation_PCR_fwd_bp']} bp | "
                f"{row['orientation_PCR_rev_bp']} |")
        table.append(end)
        protocol.write_text(text[:text.index(begin)] + "\n".join(table)
                            + text[text.index(end) + len(end):])

    with (OUT / "plasmid_summary.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nwrote {len(rows) * 2} plasmid maps and plasmid_summary.csv to "
          f"{OUT.relative_to(ROOT)}/")
    if problems:
        print("\nPROBLEMS:")
        for name, problem in problems:
            print(f"  [{name}] {problem}")
        return 1
    print("\nno problems: selection disrupted, AmpR and ori intact, insert "
          "excisable, insert sequenceable with the stock primers")
    return 0


if __name__ == "__main__":
    sys.exit(main())
