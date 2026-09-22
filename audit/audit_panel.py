"""Independent pre-order audit of the P0-P8 gBlocks and plasmid maps.

This deliberately does NOT import the design library. It reads the files that
would actually be sent to a vendor and cloned, re-derives what they should
contain from the audited architecture restated below, and compares. A bug in
``txv`` cannot cancel itself out here, because nothing in this file came from
``txv``.

Run:  python audit/audit_panel.py
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


def _paths(root: Path):
    order = root / "orders" / "pvax1_ag_panel"
    return order, order / "plasmids", root / "backbone" / "pvax1_ag_egfp.gb"


ORDER_DIR, PLASMID_DIR, BACKBONE = _paths(ROOT)

# ---------------------------------------------------------------------------
# The design, restated independently from the audited architecture table.
# ---------------------------------------------------------------------------

MODULES = {
    "M": "M",
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
    "P5": ("M", "HA", "L", "A", "L", "B", "L", "E5"),
    "P6": ("LAMP1_SP", "HA", "L", "A", "L", "B", "L", "LAMP1_TMT"),
    "P7": ("M", "HA", "L", "A", "L", "E5", "P2A", "LAMP1_SP", "FLAG", "L", "B",
           "L", "LAMP1_TMT"),
    "P8": ("M", "HA", "L", "B", "L", "E5", "P2A", "LAMP1_SP", "FLAG", "L", "A",
           "L", "LAMP1_TMT"),
}

#: From construct_architecture: protein aa, ORF nt, finished plasmid bp.
AUDITED = {
    "P0": (116, 351, 3486), "P1": (150, 453, 3588), "P2": (150, 453, 3588),
    "P3": (184, 555, 3690), "P4": (184, 555, 3690), "P5": (88, 267, 3402),
    "P6": (142, 429, 3564), "P7": (183, 552, 3687), "P8": (183, 552, 3687),
}

LEFT_ARM = "GAAGAAATATAAGAGCCACC"
RIGHT_ARM = "GCTGCCTTCTGCGGGGCTTG"
INSERT_START_1B = 699          # eGFP CDS start in the parent vector
INSERT_END_1B = 1418           # eGFP CDS end (inclusive)
BACKBONE_BP = 3855

#: Epitopes that must survive intact, or the experiment measures nothing.
EPITOPES = {"SIINFEKL": "A", "ASFEAQGALANIAVDKA": "B"}

#: Enzymes whose sites must not appear inside an insert.
FORBIDDEN_IN_INSERT = {
    "BsaI": "GGTCTC", "BamHI": "GGATCC", "EcoRI": "GAATTC",
    "HindIII": "AAGCTT", "XhoI": "CTCGAG", "NheI": "GCTAGC",
    "BspQI": "GCTCTTC", "NotI": "GCGGCCGC",
}

#: Cryptic splice signals matter here because the same plasmid also drives the
#: ORF from its CMV promoter, i.e. through a nucleus.
SPLICE_DONOR = "GGTAAGT|GGTGAGT"
POLYA_SIGNAL = "AATAAA"

CODON_TABLE = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _x in enumerate(_B):
    for _j, _y in enumerate(_B):
        for _k, _z in enumerate(_B):
            CODON_TABLE[_x + _y + _z] = _AA[_i * 16 + _j * 4 + _k]

_COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def tr(dna: str) -> str:
    return "".join(CODON_TABLE[dna[i:i + 3]] for i in range(0, len(dna) - len(dna) % 3, 3))


def gc(s: str) -> float:
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def occurrences(hay: str, needle: str, both: bool = False) -> list[int]:
    pat = needle if not both else f"{needle}|{rc(needle)}"
    return [m.start() for m in re.finditer(f"(?=({pat}))", hay)]


def expected_protein(name: str) -> str:
    return "".join(MODULES[m] for m in PANEL[name])


# ---------------------------------------------------------------------------

BLOCKER, WARN, INFO = "BLOCKER", "WARNING", "INFO"


@dataclass
class Report:
    findings: list[tuple[str, str, str]] = field(default_factory=list)
    checks: int = 0

    def check(self, ok: bool, construct: str, message: str, severity: str = BLOCKER):
        self.checks += 1
        if not ok:
            self.findings.append((severity, construct, message))
        return ok

    def note(self, construct: str, message: str, severity: str = INFO):
        self.findings.append((severity, construct, message))

    @property
    def blockers(self):
        return [f for f in self.findings if f[0] == BLOCKER]

    @property
    def warnings(self):
        return [f for f in self.findings if f[0] == WARN]


def read_genbank_sequence(path: Path) -> str:
    body = path.read_text().split("ORIGIN", 1)[1].split("//")[0]
    return "".join(re.findall(r"[acgtnACGTN]+", body)).upper()


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


def main(root: Path | None = None, quiet: bool = False) -> int:
    global ORDER_DIR, PLASMID_DIR, BACKBONE
    if root is not None:
        ORDER_DIR, PLASMID_DIR, BACKBONE = _paths(root)
    r = Report()

    backbone = read_genbank_sequence(BACKBONE)
    r.check(len(backbone) == BACKBONE_BP, "backbone",
            f"parent vector is {len(backbone)} bp, expected {BACKBONE_BP}")

    fasta = read_fasta(ORDER_DIR / "gblocks.fasta")
    with (ORDER_DIR / "ordering_table.csv").open() as fh:
        table = {row["name"]: row for row in csv.DictReader(fh)}

    encodings: dict[str, set[str]] = defaultdict(set)
    linker_encodings: set[str] = set()

    for name in sorted(PANEL):
        aa_exp, orf_exp, plasmid_exp = AUDITED[name]
        protein = expected_protein(name)

        # -- the two shipped representations must agree ---------------------
        frag = fasta.get(name, "")
        r.check(bool(frag), name, "missing from gblocks.fasta")
        r.check(table.get(name, {}).get("sequence") == frag, name,
                "ordering_table.csv and gblocks.fasta disagree on the sequence")
        if not frag:
            continue

        # -- fragment layout ------------------------------------------------
        r.check(frag.startswith(LEFT_ARM), name, "fragment does not start with the left homology arm")
        r.check(frag.endswith(RIGHT_ARM), name, "fragment does not end with the right homology arm")
        r.check(len(occurrences(frag, LEFT_ARM)) == 1, name, "left arm is not unique in the fragment")
        r.check(len(occurrences(frag, RIGHT_ARM)) == 1, name, "right arm is not unique in the fragment")
        orf = frag[len(LEFT_ARM):len(frag) - len(RIGHT_ARM)]

        # -- ORF biology ----------------------------------------------------
        r.check(orf.startswith("ATG"), name, f"ORF starts {orf[:3]}, not ATG")
        r.check(len(orf) % 3 == 0, name, f"ORF is {len(orf)} nt, not a multiple of 3")
        r.check(len(orf) == orf_exp, name, f"ORF is {len(orf)} nt, audited value is {orf_exp}")
        prot = tr(orf)
        r.check(prot.count("*") == 1, name, f"{prot.count('*')} stop codons, expected exactly 1")
        r.check(prot.endswith("*"), name, "ORF does not end in a stop codon")
        body = prot.rstrip("*")
        r.check("*" not in body, name, "premature stop codon inside the ORF")
        r.check(body == protein, name, "translated protein does not match the audited design")
        r.check(len(body) == aa_exp, name, f"protein is {len(body)} aa, audited value is {aa_exp}")

        # -- module integrity ------------------------------------------------
        cursor = 0
        for module in PANEL[name]:
            pep = MODULES[module]
            seg = orf[cursor * 3:(cursor + len(pep)) * 3]
            r.check(tr(seg) == pep, name, f"module {module} is not intact at its expected position")
            (linker_encodings if module == "L" else encodings[module]).add(seg)
            cursor += len(pep)

        for epi, owner in EPITOPES.items():
            if owner in PANEL[name]:
                r.check(epi in body, name, f"epitope {epi} missing from antigen {owner}")

        # -- positional constraints -------------------------------------------
        if "LAMP1_TMT" in PANEL[name]:
            r.check(body.endswith("GYQTI"), name,
                    "LAMP1 construct does not end in GYQTI; AP-3 sorting would fail")
        if name in ("P0", "P1", "P2", "P3", "P4", "P5"):
            r.check(body.endswith("EEEEE"), name,
                    "E5 does not terminate the protein, so the C-degron has no free C-terminus")
        if name in ("P7", "P8"):
            r.check(PANEL[name].index("E5") < PANEL[name].index("P2A"), name,
                    "E5 is not upstream of P2A as the design specifies")
        for i, module in enumerate(PANEL[name]):
            if module.endswith("_SP"):
                r.check(i == 0 or PANEL[name][i - 1] == "P2A", name,
                        f"signal peptide {module} is not at the N-terminus of a polypeptide")
        r.check("HA" in PANEL[name], name, "no HA tag: expression could not be distinguished from presentation")
        r.check(("FLAG" in PANEL[name]) == (name in ("P7", "P8")), name,
                "FLAG tag present on a construct that should not carry one (or missing from P7/P8)")
        r.check(("P2A" in PANEL[name]) == (name in ("P7", "P8")), name,
                "P2A present on a construct that should not carry one (or missing from P7/P8)")

        # -- cloning ----------------------------------------------------------
        for enzyme, site in FORBIDDEN_IN_INSERT.items():
            hits = occurrences(orf, site, both=True)
            sev = BLOCKER if enzyme in ("BsaI", "BamHI") else WARN
            r.check(not hits, name, f"{enzyme} site ({site}) inside the insert at {hits}", sev)

        r.check(len(occurrences(backbone, LEFT_ARM)) == 1, name, "left arm is not unique in the backbone")
        r.check(len(occurrences(backbone, RIGHT_ARM)) == 1, name, "right arm is not unique in the backbone")

        # -- nuclear-transcription hazards (the CMV route) ---------------------
        r.check(not occurrences(orf, POLYA_SIGNAL), name,
                f"AATAAA inside the ORF at {occurrences(orf, POLYA_SIGNAL)}: would truncate the CMV transcript")
        donors = [m.start() for m in re.finditer(SPLICE_DONOR, orf)]
        r.check(not donors, name, f"cryptic splice donor in the ORF at {donors}", WARN)

        # -- manufacturability -------------------------------------------------
        runs = [(m.group(0)[0], len(m.group(0)), m.start())
                for m in re.finditer(r"(A+|C+|G+|T+)", frag) if len(m.group(0)) >= 8]
        r.check(not runs, name, f"homopolymer run(s) of 8+ nt in the fragment: {runs}", WARN)

        seen, reps = {}, []
        for i in range(len(orf) - 15 + 1):
            k = orf[i:i + 15]
            if k in seen:
                reps.append((k, seen[k], i))
            seen[k] = i
        r.check(not reps, name, f"direct repeat(s) of 15+ nt in the ORF: {reps}", WARN)

        windows = [gc(frag[i:i + 100]) for i in range(0, max(1, len(frag) - 99), 10)]
        r.check(0.25 <= gc(frag) <= 0.75, name, f"fragment GC {gc(frag):.1%} outside 25-75%")
        r.check(all(0.25 <= w <= 0.75 for w in windows), name,
                f"a 100-nt window is at {max(windows):.0%}/{min(windows):.0%} GC", WARN)
        r.check(125 <= len(frag) <= 3000, name, f"fragment is {len(frag)} nt, outside the 125-3000 nt range")

        # -- the assembled plasmid ---------------------------------------------
        gb = PLASMID_DIR / f"pVax1_AG_{name}.gb"
        r.check(gb.exists(), name, "plasmid map missing")
        if gb.exists():
            plasmid = read_genbank_sequence(gb)
            r.check(len(plasmid) == plasmid_exp, name,
                    f"plasmid is {len(plasmid)} bp, audited value is {plasmid_exp}")

            # Simulate the assembly from the FRAGMENT, independently.
            start = backbone.index(LEFT_ARM) + len(LEFT_ARM)
            end = backbone.index(RIGHT_ARM)
            rebuilt = backbone[:start] + orf + backbone[end:]
            r.check(rebuilt == plasmid, name,
                    "the shipped plasmid map is not what this gBlock would actually assemble into")
            r.check(start == INSERT_START_1B - 1 and end == INSERT_END_1B, name,
                    f"insert window is {start + 1}..{end}, expected "
                    f"{INSERT_START_1B}..{INSERT_END_1B}")

            r.check(plasmid[start - 6:start + 3] == "GCCACCATG", name,
                    f"5' junction reads {plasmid[start - 6:start + 3]}, not GCCACCATG")
            r.check(plasmid[start + len(orf):start + len(orf) + 20] == RIGHT_ARM, name,
                    "3' junction does not run straight into the 3' UTR")

            bsa = occurrences(plasmid, "GGTCTC", both=True)
            r.check(len(bsa) == 2, name, f"{len(bsa)} BsaI sites in the plasmid, expected 2")
            r.check(not any(start <= p < start + len(orf) for p in bsa), name,
                    "a BsaI site lies inside the insert; the IVT template would be cut in half")

            t7 = occurrences(plasmid, "TAATACGACTCACTATA")
            r.check(len(t7) == 1, name, f"{len(t7)} T7 promoters found, expected 1")
            if t7:
                plus1 = t7[0] + 17
                r.check(plasmid[plus1:plus1 + 3] == "AGG", name,
                        f"transcript starts {plasmid[plus1:plus1 + 3]}, not AGG (CleanCap AG needs AG)")
            utr5 = plasmid[start - 44:start]
            r.check("ATG" not in utr5, name, "upstream AUG in the 5' UTR would capture scanning ribosomes")
            r.check("A" * 120 in plasmid, name, "the 120-nt poly(A) tract is not intact")
            for element in ("CMV promoter region", "KanR", "ori"):
                pass  # covered by the exact-rebuild check above

    # -- cross-construct -------------------------------------------------------
    for module, encs in sorted(encodings.items()):
        r.check(len(encs) == 1, "panel",
                f"module {module} has {len(encs)} different DNA encodings across the panel; "
                "a difference between arms could then be codon usage rather than routing")
    r.check(len(linker_encodings) >= 2, "panel",
            "all GGGGS linkers share one encoding, creating exact direct repeats", WARN)
    r.check(all(tr(e) == "GGGGS" for e in linker_encodings), "panel",
            "a linker encoding does not translate to GGGGS")

    r.check(expected_protein("P3") != expected_protein("P4"), "panel", "P3 and P4 are identical")
    r.check(len(expected_protein("P3")) == len(expected_protein("P4")), "panel",
            "P3 and P4 differ in length; they are meant to be an order control")
    r.check(len(expected_protein("P7")) == len(expected_protein("P8")), "panel",
            "P7 and P8 differ in length; they are meant to be a swap control")
    swap = {"A": "B", "B": "A"}
    r.check(tuple(swap.get(m, m) for m in PANEL["P7"]) == PANEL["P8"], "panel",
            "P8 is not the antigen swap of P7")
    r.check(MODULES["A"] not in expected_protein("P0")
            and MODULES["B"] not in expected_protein("P0"), "panel",
            "the P0 scaffold control carries an antigen")

    # -- report ----------------------------------------------------------------
    if quiet:
        return 1 if r.blockers else 0
    width = 78
    print("=" * width)
    print("PRE-ORDER AUDIT — P0-P8 gBlocks and plasmid maps")
    print("=" * width)
    print(f"{r.checks} checks across {len(PANEL)} constructs, "
          f"re-derived from the audited architecture\n")

    if not r.findings:
        print("  No findings. Nothing blocks ordering.")
    for severity in (BLOCKER, WARN, INFO):
        group = [f for f in r.findings if f[0] == severity]
        if not group:
            continue
        print(f"  {severity} ({len(group)})")
        for _, construct, message in group:
            print(f"    [{construct}] {message}")
        print()

    print("-" * width)
    verdict = "SAFE TO ORDER" if not r.blockers else "DO NOT ORDER"
    print(f"{verdict}: {len(r.blockers)} blocker(s), {len(r.warnings)} warning(s)")
    return 1 if r.blockers else 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]) if len(sys.argv) > 1 else None))
