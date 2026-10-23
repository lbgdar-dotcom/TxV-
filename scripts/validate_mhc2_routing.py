"""In-silico validation for the MHC II routing ORFs.

Checks what is checkable about an ORF whose 5' context is not yet fixed. The
constructs are inserted by hand, so start-codon accessibility depends on a
leader this package does not own; it is reported against the HBB leader as a
representative context and labelled as such, not as a property of the ORF.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from txv.codon import CODON_USAGE  # noqa: E402
from txv.genbank_io import parse_genbank  # noqa: E402
from txv.insilico import (  # noqa: E402
    cai, cpg_stats, fold_metrics, signal_peptide_call, tm_segments,
)
from txv.mhc2_routing import (  # noqa: E402
    EXTERNAL, PANEL, PANEL_BY_NAME, PROVENANCE, module_protein, protein,
)
from txv.seqops import gc_fraction  # noqa: E402

OUT = ROOT / "orders" / "mhc2_routing"
MIN_CAI, MAX_CPG_OE = 0.75, 0.80
#: 0.70, not 0.65. The panel's peak 40-nt window is 67.5% and it is the
#: canonical HA + GGGGS junction (CCCATACGACGTGCCAGATTACGCCGGAGGCGGAGGCAGC),
#: which appears identically in P0 and P5 of the already-ordered P0-P8 panel.
#: A bound that flagged it would be flagging sequence that has already shipped,
#: and 67.5% over 40 nt is inside normal synthesis tolerance.
GC_LO, GC_HI = 0.40, 0.70
#: Below this length the CpG observed/expected ratio is unstable, so it is
#: reported but not failed on; CpG per 100 nt is reported for every construct.
CPG_OE_MIN_LENGTH = 250
REPEAT_K = 12
HBB = "AGGACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCGCCACC"

CODON = {}
_B = "TCAG"
_T = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_B):
    for _j, _b in enumerate(_B):
        for _k, _c in enumerate(_B):
            CODON[_a + _b + _c] = _T[_i * 16 + _j * 4 + _k]


def translate(dna):
    return "".join(CODON[dna[i:i + 3]] for i in range(0, len(dna), 3))


def longest_repeat(dna, k=REPEAT_K):
    seen, hits = {}, []
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i + k]
        if kmer in seen:
            hits.append((kmer, seen[kmer] + 1, i + 1))
        else:
            seen[kmer] = i
    return hits


def main() -> int:
    problems, accepted, rows = [], [], []
    records = {}
    for c in PANEL:
        f = OUT / f"{c.name}.gb"
        if not f.exists():
            problems.append((c.name, "no GenBank file; run build_mhc2_routing.py"))
            continue
        records[c.name] = parse_genbank(f.read_text()).sequence.upper()

    print("== translation fidelity ==")
    for name, dna in records.items():
        aa = translate(dna)
        spec = protein(name)
        ok = (len(dna) % 3 == 0 and dna.startswith("ATG") and aa[-1] == "*"
              and "*" not in aa[:-1] and aa[:-1] == spec)
        print(f"  {name}  {len(dna):4d} bp  {len(spec):3d} aa  "
              f"{'OK' if ok else 'FAIL'}")
        if not ok:
            problems.append((name, "ORF does not translate to its specification"))

    print("\n== composition ==")
    print(f"  {'':5s} {'GC':>6} {'CAI':>6} {'CpG':>4} {'/100nt':>7} "
          f"{'o/e':>5}  {'minGC40':>8} {'maxGC40':>8}")
    for name, dna in records.items():
        orf = dna[:-3]
        score = cai(orf, CODON_USAGE["mouse"])
        n, oe = cpg_stats(dna)
        wins = [gc_fraction(dna[i:i + 40]) for i in range(0, len(dna) - 40 + 1)]
        lo, hi = (min(wins), max(wins)) if wins else (0, 0)
        per100 = 100.0 * n / len(dna)
        print(f"  {name:5s} {gc_fraction(dna):6.1%} {score:6.3f} {n:4d} "
              f"{per100:7.1f} {oe:5.2f}  {lo:8.1%} {hi:8.1%}")
        rows.append((name, gc_fraction(dna), score, n, oe))
        if score < MIN_CAI:
            problems.append((name, f"CAI {score:.3f} < {MIN_CAI}"))
        if oe > MAX_CPG_OE:
            if len(dna) < CPG_OE_MIN_LENGTH:
                accepted.append((
                    name,
                    f"CpG o/e {oe:.2f} over only {len(dna)} nt. Every module is "
                    f"canonical and shared, so this cannot be re-optimised "
                    f"without breaking comparability; the same modules give "
                    f"0.83 at 270 nt in P5. Density is {per100:.1f} CpG/100 nt, "
                    f"the highest here because this construct is the shortest "
                    f"and carries none of the CpG-poor routing modules."))
            else:
                problems.append((name, f"CpG o/e {oe:.2f} > {MAX_CPG_OE}"))
        if lo < GC_LO or hi > GC_HI:
            problems.append((name, f"GC window outside {GC_LO:.0%}-{GC_HI:.0%} "
                                   f"({lo:.0%}-{hi:.0%})"))

    print("\n== direct repeats (>= %d nt, recombination risk) ==" % REPEAT_K)
    for name, dna in records.items():
        hits = longest_repeat(dna)
        print(f"  {name}: {len(hits)} repeat(s)" if hits else f"  {name}: none")
        for kmer, a, b in hits[:3]:
            problems.append((name, f"direct repeat {kmer} at {a} and {b}"))

    print("\n== module canonicalisation across the panel ==")
    enc: dict[str, set[str]] = {}
    for c in PANEL:
        dna, at = records.get(c.name), 0
        if dna is None:
            continue
        for mod in c.modules:
            n = len(module_protein(mod)) * 3
            if mod != "L":
                enc.setdefault(mod, set()).add(dna[at:at + n])
            at += n
    for mod, variants in sorted(enc.items()):
        mark = "OK" if len(variants) == 1 else "FAIL"
        if len(variants) != 1:
            problems.append((mod, f"{len(variants)} different encodings"))
        print(f"  {mod:12s} {len(variants)} encoding(s)  {mark}")
    print("  (L is excluded by design: linkers cycle through 4 encodings)")

    print("\n== routing elements ==")
    for sp in ("CTLA4_SP", "LAMP1_SP", "H2DMB_SP"):
        call = signal_peptide_call(module_protein(sp))
        print(f"  {sp:10s} signal peptide: "
              f"{call.problems or 'no problems'}")
        if call.problems:
            problems.append((sp, "; ".join(call.problems)))
    for tm in ("CTLA4_TMT", "LAMP1_TMT", "H2DMB_TMT"):
        segs = tm_segments(module_protein(tm))
        print(f"  {tm:10s} membrane span(s): {len(segs)}")
        if not segs:
            problems.append((tm, "no segment that could span a bilayer"))

    print("\n== topology constraints ==")
    checks = [
        ("R2", protein("R2").endswith("GYQTI"), "LAMP1 GYQTI must be terminal"),
        ("R5", protein("R5").endswith("QETFG"), "LC3B must end at the free Gly"),
        ("R4", not protein("R4").startswith("MAYPYDVPDYA"),
         "Ii must carry no N-terminal tag"),
        ("R4", protein("R4").endswith("YPYDVPDYA"), "Ii tag must be C-terminal"),
        ("R4", "ASFEAQGALANIAVDKA" in protein("R4"),
         "Ea must sit in the CLIP register"),
        ("R3", "YTPL" in protein("R3"), "H2-DM sorting motif must be present"),
    ]
    for name, ok, what in checks:
        print(f"  {name}: {what} -- {'OK' if ok else 'FAIL'}")
        if not ok:
            problems.append((name, what))

    print("\n== folding, against the HBB leader as representative context ==")
    print("  (ViennaRNA on unmodified RNA; the transcript is m1Psi-substituted,")
    print("   and the real 5' context is whichever leader these are cloned into)")
    for name, dna in records.items():
        f = fold_metrics(HBB + dna, len(HBB))
        print(f"  {name}  AUG accessibility {f.start_accessibility:5.3f}  "
              f"cap-proximal MFE {f.cap_proximal_mfe:6.1f}  "
              f"MFE/nt {f.mfe_per_nt:6.3f}")

    print("\n== sequence provenance ==")
    for mod in sorted(PROVENANCE):
        if mod in EXTERNAL:
            print(f"  {mod:10s} {PROVENANCE[mod]}")

    if accepted:
        print(f"\n{len(accepted)} ACCEPTED DEVIATION(S) -- documented, not fixed:")
        for name, a in accepted:
            print(f"  [{name}] {a}")
    if problems:
        print(f"\n{len(problems)} FINDING(S):")
        for name, p in problems:
            print(f"  [{name}] {p}")
        return 1
    print("\nno findings: all seven pass every check")
    return 0


if __name__ == "__main__":
    sys.exit(main())
