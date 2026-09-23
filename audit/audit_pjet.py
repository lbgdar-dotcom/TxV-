"""Independent pre-order audit of the pJET1.2 transcription units.

Like ``audit_panel.py``, this imports nothing from the design library. It reads
the files a vendor would receive, restates the design from scratch, and
compares. A bug in ``txv`` cannot cancel itself out here.

Run:  python audit/audit_pjet.py
Exit: 0 if nothing blocks ordering, 1 otherwise.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _dir(root: Path) -> Path:
    return root / "orders" / "pjet12_ivt_units"


# ---------------------------------------------------------------------------
# The design, restated independently.
# ---------------------------------------------------------------------------

MODULES = {
    "MA": "MA",
    "CTLA4_SP": "MACLGLRRYKAQLQLPSRTWPFVALLTLLFIPVFS",
    "CTLA4_TMT": "FLLWILVAVSLGLFFYSFLVSAVSLSKMLKKRSPLTTGVYVKMPPTEPECEKQFQPYFIPIN",
    "LAMP1_SP": "MAAPGARRPLLLLLLAGLAHGASA",
    "LAMP1_TMT": "MLIPIAVGGALAGLVLIVLIAYLIGRKRSHAGYQTI",
    "A": "EVSGLEQLESIINFEKLTEWTSSNVMEER",
    "B": "EEFAKFASFEAQGALANIAVDKANLDVMK",
    "HA": "YPYDVPDYA",
    "FLAG": "DYKDDDDK",
    "L": "GGGGS",
    "E5": "EEEEE",
    "P2A": "GSGATNFSLLKQAGDVEENPGP",
}

PANEL = {
    "P0": ("CTLA4_SP", "HA", "L", "CTLA4_TMT", "E5"),
    "P1": ("CTLA4_SP", "HA", "L", "A", "L", "CTLA4_TMT", "E5"),
    "P2": ("CTLA4_SP", "HA", "L", "B", "L", "CTLA4_TMT", "E5"),
    "P3": ("CTLA4_SP", "HA", "L", "A", "L", "B", "L", "CTLA4_TMT", "E5"),
    "P4": ("CTLA4_SP", "HA", "L", "B", "L", "A", "L", "CTLA4_TMT", "E5"),
    "P5": ("MA", "HA", "L", "A", "L", "B", "L", "E5"),
    "P6": ("LAMP1_SP", "HA", "L", "A", "L", "B", "L", "LAMP1_TMT"),
    "P7": ("MA", "HA", "L", "A", "L", "E5", "P2A", "LAMP1_SP", "FLAG", "L", "B",
           "L", "LAMP1_TMT"),
    "P8": ("MA", "HA", "L", "B", "L", "E5", "P2A", "LAMP1_SP", "FLAG", "L", "A",
           "L", "LAMP1_TMT"),
}

#: Protein lengths from the audited architecture table, plus the one Kozak
#: alanine added to P5/P7/P8 so the whole panel shares a Kozak context.
AUDITED_AA = {"P0": 116, "P1": 150, "P2": 150, "P3": 184, "P4": 184,
              "P5": 89, "P6": 142, "P7": 184, "P8": 184}

T7_CORE = "TAATACGACTCACTATA"
#: Quoted from the IVT kit manual's "Template requirements" section:
#: "Minimum T7 promotor sequences: 5'-TAATACGACTCACTATAAGG"
#: (VENI all-in-one mRNA Synthesis Kit with Cap1 Analog, Leish Bio). The kit's
#: cap1 analog is the AG-initiating trinucleotide, so a transcript that does
#: not begin AG will not cap.
KIT_MIN_PROMOTER = "TAATACGACTCACTATAAGG"
LEADER = "AGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACC"
#: BNT162b2 3' UTR (AES + mtRNR1), 296 nt -- restated here from the published
#: sequence, not imported, so a corrupted constant in txv cannot pass unseen.
UTR3 = ("CTCGAGCTGGTACTGCATGCACGCAATGCTAGCTGCCCCTTTCCCGTCCTGGGTACCCCGAG"
        "TCTCCCCCGACCTCGGGTCCCAGGTATGCTCCCACCTCCACCTGCCCCACTCACCACCTCTG"
        "CTAGTTCCAGACACCTCCCAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCCACAC"
        "CCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACT"
        "AACCCCAGGGTTGGTCAATTTCGTGCCAGCCACACCCTGGAGCTAGCA")
FWD_HANDLE = "CAGTCACGTAGCATCGACTG"
REV_HANDLE = "GTCATGCAGTCGATCACTGA"

EPITOPES = {"SIINFEKL": "A", "ASFEAQGALANIAVDKA": "B"}

CODON = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _x in enumerate(_B):
    for _j, _y in enumerate(_B):
        for _k, _z in enumerate(_B):
            CODON[_x + _y + _z] = _AA[_i * 16 + _j * 4 + _k]

_COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def tr(dna: str) -> str:
    return "".join(CODON[dna[i:i + 3]] for i in range(0, len(dna) - len(dna) % 3, 3))


def gc(s: str) -> float:
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def occ(hay: str, needle: str, both: bool = False) -> list[int]:
    pat = needle if not both else f"{needle}|{rc(needle)}"
    return [m.start() for m in re.finditer(f"(?=({pat}))", hay)]


BLOCKER, WARN = "BLOCKER", "WARNING"

#: Ceiling on coding-sequence GC. Not a synthesis limit -- vendors will build
#: well above this -- but the band the panel was optimised into, so a drift
#: above it means an encoding changed without the optimiser config changing.
GC_CEILING = 0.60
#: How far apart the constructs' GC may be. Routing is the variable under test;
#: GC should not be a second one.
GC_SPREAD = 0.06


@dataclass
class Report:
    findings: list = field(default_factory=list)
    checks: int = 0

    def check(self, ok, name, msg, sev=BLOCKER):
        self.checks += 1
        if not ok:
            self.findings.append((sev, name, msg))

    @property
    def blockers(self):
        return [f for f in self.findings if f[0] == BLOCKER]

    @property
    def warnings(self):
        return [f for f in self.findings if f[0] == WARN]


def read_fasta(path: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split("|")[0].strip(), []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def gb_sequence(path: Path) -> str:
    body = path.read_text().split("ORIGIN", 1)[1].split("//")[0]
    return "".join(re.findall(r"[acgtnACGTN]+", body)).upper()


def gb_labels(path: Path) -> list[str]:
    return re.findall(r'/label="([^"]+)"', path.read_text())


def main(root: Path | None = None, quiet: bool = False) -> int:
    d = _dir(root or ROOT)
    r = Report()

    fasta = read_fasta(d / "gblocks.fasta")
    with (d / "ordering_table.csv").open() as fh:
        table = {row["name"]: row for row in csv.DictReader(fh)}
    with (d / "element_manifest.csv").open() as fh:
        manifest = {row["construct"]: row for row in csv.DictReader(fh)}

    encodings: dict[str, set[str]] = defaultdict(set)
    linker_encodings: set[str] = set()
    kozak_contexts: dict[str, list[str]] = defaultdict(list)
    orf_gc: dict[str, float] = {}

    for name in sorted(PANEL):
        modules = PANEL[name]
        protein = "".join(MODULES[m] for m in modules)

        frag = fasta.get(name, "")
        r.check(bool(frag), name, "missing from gblocks.fasta")
        r.check(table.get(name, {}).get("sequence") == frag, name,
                "ordering_table.csv and gblocks.fasta disagree")
        if not frag:
            continue

        # -- layout ---------------------------------------------------------
        r.check(frag.startswith(FWD_HANDLE), name, "does not start with IVT_F handle")
        r.check(frag.endswith(REV_HANDLE), name, "does not end with IVT_R handle")
        r.check(len(occ(frag, FWD_HANDLE, True)) == 1, name, "IVT_F handle not unique")
        r.check(len(occ(frag, REV_HANDLE, True)) == 1, name, "IVT_R handle not unique")

        t7 = occ(frag, T7_CORE)
        r.check(len(t7) == 1, name, f"{len(t7)} T7 promoters, expected 1")
        r.check(len(occ(frag, KIT_MIN_PROMOTER)) == 1, name,
                f"the IVT kit's minimum promoter {KIT_MIN_PROMOTER} is not "
                "present exactly once; the transcript would not cap")
        if t7:
            r.check(frag[t7[0] + 17:t7[0] + 20] == "AGG", name,
                    f"transcript starts {frag[t7[0]+17:t7[0]+20]}, not AGG")
            r.check(t7[0] == len(FWD_HANDLE), name, "T7 does not follow the handle")

        r.check(LEADER in frag, name, "5' leader absent or altered")
        r.check(UTR3 in frag, name, "3' UTR absent or altered")
        if LEADER not in frag:
            continue
        r.check("ATG" not in LEADER, name, "upstream AUG in the leader")

        orf_start = frag.index(LEADER) + len(LEADER)
        r.check(frag[orf_start - 6:orf_start + 3] == "GCCACCATG", name,
                "Kozak/ATG junction is not GCCACCATG")
        r.check(frag[orf_start + 3] == "G", name,
                f"Kozak +4 is {frag[orf_start + 3]}, not G")
        kozak_contexts[frag[orf_start - 6:orf_start + 6]].append(name)

        # -- ORF ------------------------------------------------------------
        aa = tr(frag[orf_start:])
        r.check("*" in aa, name, "no in-frame stop codon")
        if "*" not in aa:
            continue
        body = aa[:aa.index("*")]
        orf_end = orf_start + (len(body) + 1) * 3
        orf = frag[orf_start:orf_end]

        r.check(body == protein, name, "translated protein does not match the design")
        r.check(len(body) == AUDITED_AA[name], name,
                f"protein is {len(body)} aa, audited value {AUDITED_AA[name]}")
        r.check(len(orf) % 3 == 0, name, "ORF is not a multiple of 3")
        orf_gc[name] = gc(orf)
        r.check(orf.count("TGA") + orf.count("TAA") + orf.count("TAG") >= 1, name,
                "no stop codon present")
        r.check(frag[orf_end:orf_end + len(UTR3)] == UTR3, name,
                "3' UTR does not follow the stop codon")
        r.check(frag[orf_end + len(UTR3):] == REV_HANDLE, name,
                "reverse handle does not follow the 3' UTR")

        # -- every module, at its expected position -------------------------
        cursor = 0
        for module in modules:
            pep = MODULES[module]
            seg = orf[cursor * 3:(cursor + len(pep)) * 3]
            r.check(tr(seg) == pep, name,
                    f"module {module} not intact at its expected position")
            (linker_encodings if module == "L" else encodings[module]).add(seg)
            cursor += len(pep)

        for epi, owner in EPITOPES.items():
            if owner in modules:
                r.check(epi in body, name, f"epitope {epi} missing from antigen {owner}")

        # -- positional constraints -----------------------------------------
        if "LAMP1_TMT" in modules:
            r.check(body.endswith("GYQTI"), name,
                    "LAMP1 construct does not end in GYQTI; AP-3 sorting fails")
        if "E5" in modules and "P2A" not in modules:
            r.check(body.endswith("EEEEE"), name,
                    "E5 does not terminate the protein, so it has no free C-terminus")
        if "P2A" in modules:
            r.check(modules.index("E5") < modules.index("P2A"), name,
                    "E5 is not upstream of P2A")
            r.check(MODULES["P2A"] in body, name, "P2A skip peptide not intact")
        for i, m in enumerate(modules):
            if m.endswith("_SP"):
                r.check(i == 0 or modules[i - 1] == "P2A", name,
                        f"{m} is not at the N-terminus of a polypeptide")
        r.check("HA" in modules, name, "no HA tag")
        r.check(("FLAG" in modules) == (name in ("P7", "P8")), name,
                "FLAG tag on the wrong construct")

        # -- the manifest must agree with the molecule ----------------------
        row = manifest.get(name, {})
        r.check(bool(row), name, "missing from element_manifest.csv")
        if row:
            r.check(int(row["protein_aa"]) == len(body), name,
                    "manifest protein length disagrees with the sequence")
            r.check(int(row["fragment_bp"]) == len(frag), name,
                    "manifest fragment length disagrees with the sequence")
            r.check(row["module_order"] == " - ".join(modules), name,
                    "manifest module order disagrees with the design")
            for module in set(modules):
                if module in row:
                    r.check(int(row[module]) == modules.count(module), name,
                            f"manifest count for {module} is wrong")

        # -- the annotated record must describe the same molecule ------------
        gb = d / f"{name}.gb"
        r.check(gb.exists(), name, "GenBank record missing")
        if gb.exists():
            r.check(gb_sequence(gb) == frag, name,
                    "the GenBank record is not the sequence being ordered")
            labels = gb_labels(gb)
            for module in set(modules):
                r.check(any(l == module or l.startswith(f"{module}_") for l in labels),
                        name, f"{module} is not annotated in the GenBank record")
            for required in ("T7_promoter", "5'UTR", "Kozak", "ORF", "3'UTR",
                             "IVT_F_handle", "IVT_R_handle"):
                r.check(required in labels, name, f"{required} not annotated")

        # -- manufacturability ----------------------------------------------
        r.check(125 <= len(frag) <= 3000, name, f"{len(frag)} nt outside 125-3000")
        r.check(0.25 <= gc(frag) <= 0.75, name, f"GC {gc(frag):.1%} outside 25-75%")
        windows = [gc(frag[i:i + 100]) for i in range(0, max(1, len(frag) - 99), 10)]
        r.check(all(0.25 <= w <= 0.75 for w in windows), name,
                f"a 100-nt window hits {max(windows):.0%}/{min(windows):.0%} GC", WARN)
        runs = [(m.group(0)[0], len(m.group(0)))
                for m in re.finditer(r"(A+|C+|G+|T+)", frag) if len(m.group(0)) >= 8]
        r.check(not runs, name, f"homopolymer run(s): {runs}", WARN)
        seen, reps = set(), []
        for i in range(len(frag) - 15 + 1):
            k = frag[i:i + 15]
            if k in seen:
                reps.append(k)
            seen.add(k)
        r.check(not reps, name, f"15-nt direct repeat(s): {reps}", WARN)
        r.check("A" * 20 not in frag, name,
                "an encoded poly(A) tract is present; the tail should come from "
                "the reverse primer")

    # -- cross-construct ------------------------------------------------------
    for module, encs in sorted(encodings.items()):
        r.check(len(encs) == 1, "panel",
                f"module {module} has {len(encs)} DNA encodings across the panel; "
                "a difference between arms could then be codon usage, not routing")
    r.check(len(linker_encodings) >= 2, "panel",
            "all GGGGS linkers share one encoding, creating direct repeats", WARN)

    # One Kozak context for the whole panel. The -6..-1 element is supplied by
    # the shared 5' UTR, but +4 is the first base of codon 2, so it belongs to
    # whichever module starts the ORF -- which is exactly how three constructs
    # ended up with a different +4 from the other six. A panel that compares
    # routes cannot afford an initiation difference that tracks the route.
    r.check(len(kozak_contexts) == 1, "panel",
            "constructs do not share one Kozak context: "
            + "; ".join(f"{ctx} ({', '.join(names)})"
                        for ctx, names in sorted(kozak_contexts.items())))

    # GC is held to a tighter band than the 25-75% synthesis limit above,
    # because the point is comparability and vendor success, not legality.
    if orf_gc:
        hi, lo = max(orf_gc.values()), min(orf_gc.values())
        worst = max(orf_gc, key=orf_gc.get)
        r.check(hi <= GC_CEILING, "panel",
                f"{worst} ORF GC is {hi:.1%}, above the {GC_CEILING:.0%} ceiling")
        r.check(hi - lo <= GC_SPREAD, "panel",
                f"ORF GC spans {lo:.1%}-{hi:.1%} ({hi - lo:.1%}), wider than the "
                f"{GC_SPREAD:.0%} the panel allows", WARN)

    swap = {"A": "B", "B": "A"}
    r.check(tuple(swap.get(m, m) for m in PANEL["P7"]) == PANEL["P8"], "panel",
            "P8 is not the antigen swap of P7")
    r.check(len("".join(MODULES[m] for m in PANEL["P3"]))
            == len("".join(MODULES[m] for m in PANEL["P4"])), "panel",
            "P3 and P4 differ in length")
    for antigen in ("A", "B"):
        r.check(MODULES[antigen] not in "".join(MODULES[m] for m in PANEL["P0"]),
                "panel", "the P0 scaffold control carries an antigen")

    if quiet:
        return 1 if r.blockers else 0

    width = 78
    print("=" * width)
    print("PRE-ORDER AUDIT — P0-P8 transcription units for pJET1.2")
    print("=" * width)
    print(f"{r.checks} checks across {len(PANEL)} constructs, "
          "re-derived from the design\n")
    if not r.findings:
        print("  No findings. Nothing blocks ordering.")
    for sev in (BLOCKER, WARN):
        group = [f for f in r.findings if f[0] == sev]
        if not group:
            continue
        print(f"  {sev} ({len(group)})")
        for _, n, m in group:
            print(f"    [{n}] {m}")
        print()
    print("-" * width)
    print(f"{'SAFE TO ORDER' if not r.blockers else 'DO NOT ORDER'}: "
          f"{len(r.blockers)} blocker(s), {len(r.warnings)} warning(s)")
    return 1 if r.blockers else 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]) if len(sys.argv) > 1 else None))
