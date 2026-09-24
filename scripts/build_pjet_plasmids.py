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

from txv.genbank_io import GenBankFeature, write_genbank  # noqa: E402
from txv.ivt_unit import (  # noqa: E402
    FWD_PRIMER, REV_HANDLE, make_unit, reverse_primer, simulate_pcr,
    unit_elements,
)
from txv.plasmid import (  # noqa: E402
    insert_blunt, revcomp_record, rotate_record,
)
from txv.pvax1_ag import PANEL, build_panel_construct  # noqa: E402
from txv.ivt_unit import T7_CORE  # noqa: E402
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

#: pJET1.2's stock sequencing primers, as read from the SnapGene file the first
#: time and kept here so a rebuild from the GenBank map still knows them.
PJET_STOCK_PRIMERS = {
    "pJET1.2 forward sequencing primer": "CGACTCACTATAGGGAGAGCGGC",
    "pJET1.2 reverse sequencing primer": "AAGAACATCGATTTTCCATGGCAG",
}
#: A Sanger read is reliable for roughly this many bases past the primer.
SANGER_READ = 850

#: Length of the poly(A) carried in the finished plasmid.
#:
#: The tail is **not** synthesised into the gene fragment. A 120-nt A-tract
#: fails vendor screening outright, and so does the BNT162b2-style segmented
#: A30-linker-A70, whose 70-nt run is still past what a gBlock will build.
#: So the fragment ships tail-free and the tail is added by PCR with
#: IVT_R_polyA120 **before cloning**, and that product is what goes into
#: pJET1.2. The plasmid then carries the poly(A) encoded, exactly as
#: pVax1_AG_eGFP does, every miniprep has it, and the Ultramer is spent once
#: rather than on every prep.
POLYA = 120

#: How much sequence may follow the poly(A) in the run-off transcript. The
#: tail should be at or near the 3' end; pVax1_AG manages ~5 nt by cutting
#: BsaI one base past it, and a dozen is the same order.
MAX_NT_AFTER_POLYA = 20

#: GenBank feature kinds by element kind, for transplanting the fragment's own
#: annotation into the finished plasmid.
_KEY = {"handle": "primer_bind", "promoter": "promoter", "utr5": "5'UTR",
        "utr3": "3'UTR", "kozak": "regulatory", "orf": "CDS",
        "signal_peptide": "sig_peptide", "skip_peptide": "misc_feature",
        "degron": "misc_feature", "tag": "misc_feature",
        "antigen": "misc_feature", "linker": "misc_feature",
        "transmembrane": "misc_feature", "start": "misc_feature",
        "stop": "misc_feature"}


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



def read_span(sequence: str, primer: str, read_len: int) -> tuple[set[int], str]:
    """Positions a Sanger read from ``primer`` covers, and which way it reads.

    Circular-safe, which is the whole point. Once a map is rotated so the
    insert starts at position 1, the vector's forward sequencing primer sits
    near the end of the record and reads *across the origin* into the insert.
    Linear arithmetic on that gives a reach larger than the plasmid and
    declares every insert covered, which is how a broken check looks like a
    passing one.

    Returns ``(positions, direction)`` with positions taken modulo the record
    length, or an empty set if the primer does not bind.
    """
    n = len(sequence)
    top = find_all(sequence, primer)
    if top:
        start = top[0] + len(primer)          # reads 3' of the binding site
        return {(start + k) % n for k in range(read_len)}, "forward"
    bottom = find_all(sequence, revcomp(primer))
    if bottom:
        start = bottom[0]                     # reads back from the 5' side
        return {(start - 1 - k) % n for k in range(read_len)}, "reverse"
    return set(), "absent"


def released_fragment(sequence: str, cuts: list[int], start: int, end: int) -> str:
    """The circular-digest fragment containing ``start..end``, as sequence.

    After the map is rotated so the insert begins at position 1, the fragment
    that carries it wraps the origin: its two cuts are near the end and the
    middle of the record, and a plain ``sequence[a:b]`` returns the backbone
    instead. Slicing the doubled sequence handles the wrap.
    """
    cuts = sorted(cuts)
    doubled = sequence + sequence
    n = len(sequence)
    for i, lo in enumerate(cuts):
        hi = cuts[i + 1] if i + 1 < len(cuts) else cuts[0] + n
        for offset in (0, n):
            if lo <= start + offset and end + offset <= hi:
                return doubled[lo:hi]
    return ""


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



def _annotate_insert(build, unit, construct, panel, orientation: int) -> None:
    """Transplant the fragment's own element map into the finished plasmid.

    Without this the plasmid shows one anonymous "insert" block, which is no
    use for the thing these maps exist for: looking at a construct and seeing
    the signal peptide, the antigens and their order, the degron and the P2A
    where they actually sit.

    In the flipped orientation every element is mirrored within the insert, so
    coordinates are reflected and the strand inverted rather than the elements
    being dropped -- a reverse-orientation clone is a real product and its map
    should be readable too.
    """
    lo, hi = build.insert_start, build.insert_end
    span = hi - lo

    for name, start, end, kind, note in unit_elements(unit, construct):
        if orientation == 1:
            f_start, f_end, strand = lo + start, lo + end, 1
        else:
            f_start, f_end, strand = lo + span - end, lo + span - start, -1
        build.record.features.append(GenBankFeature(
            _KEY.get(kind, "misc_feature"), f_start, f_end, strand,
            {"label": name, "note": note},
        ))

    # The poly(A) is not part of the fragment, so unit_elements does not know
    # about it; it comes from the primer during the pre-cloning PCR.
    tail_len = span - len(unit)
    if tail_len > 0:
        if orientation == 1:
            t_start, t_end, strand = hi - tail_len, hi, 1
        else:
            t_start, t_end, strand = lo, lo + tail_len, -1
        build.record.features.append(GenBankFeature(
            "polyA_signal", t_start, t_end, strand,
            {"label": f"poly(A) {tail_len}",
             "note": f"A{tail_len}, added by IVT_R_polyA{tail_len} in the PCR "
                     "before cloning, not synthesised into the gene fragment. "
                     "Encoded here, so every miniprep carries it."},
        ))

    build.record.features.sort(key=lambda f: (f.start, f.end))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    from txv.genbank_io import read_genbank

    source = Path(sys.argv[1]) if len(sys.argv) > 1 else VECTOR_GB
    if not source.exists():
        print(f"pass the pJET1.2 map, or put one at {VECTOR_GB}")
        return 1
    # Accept either the SnapGene original or the GenBank this script writes,
    # so a rebuild does not depend on still having the .dna to hand.
    if source.suffix.lower() == ".dna":
        vector = read_snapgene(source)
        VECTOR_GB.write_text(write_genbank(vector))
        print(f"wrote {VECTOR_GB.relative_to(ROOT)}")
    else:
        vector = read_genbank(source)

    stock = read_primers(source) if source.suffix.lower() == ".dna" else {}
    if not stock:
        # The stock primers are a property of the vector, not of the file
        # format it arrived in.
        stock = dict(PJET_STOCK_PRIMERS)
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
        construct = build_panel_construct(panel.name)
        unit = make_unit(panel.name, construct.orf)

        # What is actually cloned is the PCR product, not the gene fragment:
        # the fragment amplified with IVT_R_polyA120 carries the poly(A).
        cloned = simulate_pcr(unit.sequence, FWD_PRIMER, reverse_primer(POLYA),
                              circular=False)
        if len(cloned) != 1 or not cloned[0].endswith("A" * POLYA):
            problems.append((panel.name, "the poly(A) PCR does not give a "
                             "single tailed product"))
            continue
        insert_seq = cloned[0]

        for orientation, tag in ((1, "fwd"), (-1, "rev")):
            build = insert_blunt(
                vector, insert_seq, site,
                f"pJET1.2_{panel.name}" + ("" if orientation == 1 else "_rev"),
                orientation=orientation,
                definition=f"pJET1.2 carrying the {panel.name} IVT transcription "
                           f"unit + A{POLYA} ({panel.route})",
            )
            _annotate_insert(build, unit, construct, panel, orientation)

            # Present every map 5' first, reading downstream. Rotating puts the
            # insert at position 1, which most viewers draw at twelve o'clock;
            # for the flipped clone the record is first shown from its other
            # strand, so the cassette reads forward there too instead of
            # mirror-imaged. Neither changes the molecule.
            record = build.record
            if orientation == -1:
                record = revcomp_record(record, record.name)
            start = next(f.start for f in record.features
                         if f.qualifiers.get("label") == "insert")
            record = rotate_record(record, start, record.name)

            out_dir = OUT if orientation == 1 else OUT / "reverse_orientation"
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / f"{record.name}.gb").write_text(write_genbank(record))

            if orientation != 1:
                continue

            # Verify the record that was actually written, not the one before
            # orientation. The two differ only in numbering, and every size
            # below is rotation-invariant -- but checking a different object
            # from the one shipped is how a map and its paperwork drift apart,
            # which is the failure this whole package is built to avoid.
            seq = record.sequence
            insert = next(f for f in record.features
                          if f.qualifiers.get("label") == "insert")
            build.insert_start, build.insert_end = insert.start, insert.end
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

            # -- run-off transcription ----------------------------------------
            # pVax1_AG carries a BsaI site 1 nt past its encoded poly(A) and is
            # linearised there for run-off. pJET1.2 gets the same behaviour for
            # free: its two BglII sites bracket the insert, so one digest both
            # releases the insert for the diagnostic gel and produces the IVT
            # template, with the poly(A) near-terminal. No extra enzyme site and
            # no tailed primer once the plasmid exists.
            runoff_nt = runoff_after = None
            runoff_ok = False
            if len(bgl) == 2:
                released = released_fragment(seq, bgl, build.insert_start,
                                             build.insert_end)
                t7_hits = find_all(released, T7_CORE)
                if len(t7_hits) == 1:
                    transcript = released[t7_hits[0] + len(T7_CORE):]
                    runoff_nt = len(transcript)
                    tail = "A" * POLYA
                    if transcript.startswith("AGG") and tail in transcript:
                        runoff_after = len(transcript) - (
                            transcript.rfind(tail) + POLYA)
                        runoff_ok = runoff_after <= MAX_NT_AFTER_POLYA
            if not runoff_ok:
                problems.append((
                    panel.name,
                    "BglII run-off does not give a usable IVT template "
                    f"(transcript {runoff_nt} nt, {runoff_after} nt past the "
                    "poly(A))"))

            # -- sequencing coverage -----------------------------------------
            covered = None
            if seq_fwd and seq_rev:
                span_f, dir_f = read_span(seq, seq_fwd, SANGER_READ)
                span_r, dir_r = read_span(seq, seq_rev, SANGER_READ)
                want = {(build.insert_start + k) % len(seq)
                        for k in range(build.insert_end - build.insert_start)}
                reached = want & (span_f | span_r)
                covered = len(reached) == len(want)
                # The shortest Sanger read that still closes the insert. A
                # "yes" that depends on an 850-nt read is worth less than the
                # number, since read length is a property of whoever runs the
                # sequencing, not of the design.
                lo, hi = 1, SANGER_READ * 2
                while lo < hi:
                    mid = (lo + hi) // 2
                    union = (read_span(seq, seq_fwd, mid)[0]
                             | read_span(seq, seq_rev, mid)[0])
                    if want <= union:
                        hi = mid
                    else:
                        lo = mid + 1
                min_read = lo
                roles = f"fwd reads {dir_f}, rev reads {dir_r}"
                if not covered:
                    problems.append((
                        panel.name,
                        f"stock primers cover {len(reached)} of "
                        f"{len(want)} insert bases ({len(want) - len(reached)} "
                        "unread); an internal sequencing primer is needed"))

            # -- telling the orientation apart --------------------------------
            fwd_products = simulate_pcr(seq, FWD_PRIMER, seq_rev or REV_HANDLE)
            rev_build = insert_blunt(vector, insert_seq, site, "tmp",
                                     orientation=-1)
            rev_products = simulate_pcr(rev_build.record.sequence,
                                        FWD_PRIMER, seq_rev or REV_HANDLE)

            rows.append({
                "construct": panel.name,
                "plasmid_bp": len(build),
                "insert_bp": len(insert_seq),
                "insert_at": build.insert_start + 1,
                "disrupted": "; ".join(build.disrupted),
                "BglII_sites": len(bgl),
                "BglII_insert_fragment_bp": insert_carrying or "",
                "BglII_backbone_fragment_bp": backbone_piece or "",
                "runoff_transcript_nt": runoff_nt or "",
                "nt_after_polyA": runoff_after if runoff_after is not None else "",
                "sequencing_covers_insert": "yes" if covered else "no",
                "min_sanger_read_nt": min_read,
                "orientation_PCR_fwd_bp": ";".join(str(len(p)) for p in fwd_products) or "none",
                "orientation_PCR_rev_bp": ";".join(str(len(p)) for p in rev_products) or "none",
            })
            print(f"  {panel.name}  {len(build):>5} bp  insert {len(insert_seq):>4} bp at "
                  f"{build.insert_start + 1}  BglII: {insert_carrying}+{backbone_piece} bp  "
                  f"run-off {runoff_nt} nt (+{runoff_after} past A{POLYA})  "
                  f"needs {min_read}-nt reads  "
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
                 "run-off transcript | nt after poly(A) | "
                 "min Sanger read | orientation PCR (designed) | (flipped) |",
                 "|---|---|---|---|---|---|---|---|---|"]
        for row in rows:
            table.append(
                f"| {row['construct']} | {row['plasmid_bp']} bp | "
                f"{row['BglII_insert_fragment_bp']} bp | "
                f"{row['BglII_backbone_fragment_bp']} bp | "
                f"{row['runoff_transcript_nt']} nt | "
                f"{row['nt_after_polyA']} | "
                f"{row['min_sanger_read_nt']} nt | "
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
          "excisable, sequenceable, and BglII run-off gives a capped-"
          "compatible IVT template with a near-terminal poly(A)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
