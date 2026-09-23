"""Self-contained IVT transcription units for blunt cloning into pJET1.2.

The orientation problem, and why it disappears
----------------------------------------------
pJET1.2/blunt is a **blunt-end** vector. Blunt ends carry no sequence
information, so T4 ligase joins the insert in either direction with roughly
equal probability. **There is no way to control orientation at the ligation
step.** That is not a flaw in the design; it is what "blunt" means.

The usual answer is directional cloning: cut the insert and the vector with two
*different* enzymes, giving two non-complementary overhangs, so only one
arrangement can anneal. pJET1.2 does not offer that, and adding it would throw
away the thing that makes the kit good -- the lethal ``eco47IR`` gene, which is
disrupted by any insert and gives >99% recombinants with no screening.

So instead of fighting the vector, this module makes orientation **irrelevant**:
each fragment is a complete transcription unit --

    [fwd handle] T7 promoter · 5' leader · Kozak · ORF · 3' UTR [rev handle]

-- and the IVT template is recovered by **PCR with your own primer pair**, not
by cutting the plasmid. A PCR product is double-stranded and its T7 promoter
sits on a defined strand of it, so it transcribes correctly whichever way the
insert happened to land. The plasmid is then nothing but a way to store and
amplify the sequence in *E. coli*, which is all you asked it to be.

Two consequences worth stating plainly:

* **No restriction enzyme is used to insert the fragment.** Blunt ligation uses
  none. The whole question of which enzymes to pick goes away.
* **The vector's own sequence never enters the transcript.** Nothing in the
  design depends on pJET1.2's internals, which is why this works without
  having its map in hand.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .seqops import clean, find_all, gc_fraction, revcomp, translate

#: T7 class III promoter, positions -17..-1. The next base is +1.
T7_CORE = "TAATACGACTCACTATA"

#: Transcript from +1 to just before the ATG: begins AGG (which is what
#: CleanCap AG requires) and ends in the strong Kozak.
LEADER = "AGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACC"

#: From the stop codon to where the encoded poly(A) began in pVax1_AG.
UTR3 = ("GCTGCCTTCTGCGGGGCTTGCCTTCTGGCCATGCCCTTCTTCTCTCCCTTGCACCTGTACC"
        "TCTTGGTCTTTGAATAAAGCCTGAGTAGGAAGGGATCC")

#: Primer-binding handles, identical on all nine fragments so one primer pair
#: amplifies every construct. The forward handle also gives T7 RNA polymerase
#: the upstream duplex it needs -- a promoter sitting flush with the end of a
#: PCR product transcribes poorly.
FWD_HANDLE = "CAGTCACGTAGCATCGACTG"
REV_HANDLE = "GTCATGCAGTCGATCACTGA"

#: The actual oligos. The forward primer IS the forward handle; the reverse
#: primer is the reverse complement of the reverse handle, because a primer
#: anneals to the strand opposite the one it is named after -- a distinction
#: that is easy to lose and gives a PCR that simply does not amplify.
FWD_PRIMER = FWD_HANDLE


def reverse_primer(polya: int = 120) -> str:
    """The IVT reverse primer, carrying the poly(T) that adds the tail.

    T's at the primer's 5' end become A's at the 3' end of the template's top
    strand, which is the strand T7 transcribes. With ``polya=0`` this is just
    the plain amplification primer.
    """
    return "T" * polya + revcomp(REV_HANDLE)


@dataclass
class IVTUnit:
    """One orderable fragment: a complete transcription unit."""

    name: str
    sequence: str
    orf: str
    protein: str
    notes: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def transcript(self) -> str:
        """What T7 makes from this unit, before the PCR-added poly(A)."""
        start = self.sequence.index(T7_CORE) + len(T7_CORE)
        return self.sequence[start:self.sequence.index(REV_HANDLE)]

    @property
    def amplicon(self) -> str:
        """The PCR product: handle to handle, inclusive."""
        return self.sequence


def make_unit(name: str, orf: str,
              fwd_handle: str = FWD_HANDLE,
              rev_handle: str = REV_HANDLE) -> IVTUnit:
    """Build the complete, self-contained fragment for one cassette."""
    orf = clean(orf)
    sequence = fwd_handle + T7_CORE + LEADER + orf + UTR3 + rev_handle
    protein = translate(orf, stop_at_stop=True).rstrip("*")
    notes = [
        "Complete transcription unit: T7 promoter, 5' leader, Kozak, ORF and "
        "3' UTR all travel with the fragment.",
        "Blunt-ligated into pJET1.2 in either orientation; the template is "
        "recovered by PCR with the handle primers, so orientation does not "
        "affect the product.",
        "No poly(A) is encoded. The reverse primer carries the poly(T) tail, "
        "which is what adds it to the template.",
    ]
    return IVTUnit(name, sequence, orf, protein, notes)


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def simulate_pcr(template: str, forward: str, reverse: str,
                 circular: bool = True) -> list[str]:
    """Amplicons produced by one primer pair on a (usually circular) template.

    ``forward`` and ``reverse`` are the **oligos**, not their binding sites: the
    reverse primer's site on the top strand is its reverse complement. Both
    strands are searched, because the point of the exercise is that the insert
    may have landed either way round.
    """
    template = clean(template)
    search = template + template if circular else template
    products = []
    for strand in (search, revcomp(search)):
        for f in find_all(strand, forward):
            r = strand.find(revcomp(reverse), f)
            if r == -1:
                continue
            end = r + len(reverse)
            if end - f <= len(template):
                products.append(strand[f:end])
    # De-duplicate, treating a sequence and its reverse complement as one
    # molecule -- a PCR product is double-stranded.
    unique: list[str] = []
    for p in products:
        if p not in unique and revcomp(p) not in unique:
            unique.append(p)
    return unique


def verify_unit(unit: IVTUnit, expected_protein: str | None = None) -> list[str]:
    """Re-derive everything that must be true of the fragment."""
    problems: list[str] = []
    s = unit.sequence

    for label, handle in (("forward", FWD_HANDLE), ("reverse", REV_HANDLE)):
        hits = find_all(s, handle, both_strands=True)
        if len(hits) != 1:
            problems.append(f"{label} handle occurs {len(hits)} times, expected 1")

    if not s.startswith(FWD_HANDLE):
        problems.append("fragment does not begin with the forward handle")
    if not s.endswith(REV_HANDLE):
        problems.append("fragment does not end with the reverse handle")

    t7 = find_all(s, T7_CORE)
    if len(t7) != 1:
        problems.append(f"{len(t7)} T7 promoters found, expected 1")
    if t7:
        plus1 = t7[0] + len(T7_CORE)
        if s[plus1:plus1 + 3] != "AGG":
            problems.append(
                f"transcript begins {s[plus1:plus1 + 3]}, not AGG "
                "(CleanCap AG needs AG)"
            )

    # A verifier must report a defect, not raise on it: anything located by
    # string search here is exactly what might be broken.
    if LEADER not in s:
        problems.append("5' leader is absent or altered")
        return problems
    leader_start = s.index(LEADER)
    if find_all(LEADER, "ATG"):
        problems.append("upstream AUG in the leader")
    orf_start = leader_start + len(LEADER)
    if s[orf_start:orf_start + 3] != "ATG":
        problems.append("ORF does not start immediately after the Kozak")
    if s[orf_start - 6:orf_start + 3] != "GCCACCATG":
        problems.append("Kozak/ATG junction is not GCCACCATG")

    orf_end = orf_start + len(unit.orf)
    if (orf_end - orf_start) % 3:
        problems.append("ORF length is not a multiple of 3")
    protein = translate(s[orf_start:], stop_at_stop=True)
    if not protein.endswith("*"):
        problems.append("no in-frame stop codon")
    elif orf_start + len(protein) * 3 != orf_end:
        problems.append("the first in-frame stop is not at the end of the ORF")
    if expected_protein is not None and protein.rstrip("*") != expected_protein:
        problems.append("translated protein does not match the design")

    if s[orf_end:orf_end + len(UTR3)] != UTR3:
        problems.append("3' UTR does not follow the stop codon")

    return problems


def verify_orientation_independence(unit: IVTUnit,
                                    mock_vector: str | None = None) -> list[str]:
    """Prove the fragment works whichever way it ligates.

    Inserts the fragment into a mock circular vector forwards and backwards,
    runs the same PCR on both, and requires the same molecule out.
    """
    problems: list[str] = []
    vector = mock_vector or ("TTGCCACGATCGATTACGGATCCAAGTGACTGATCGGTACCAT" * 8)

    forward_plasmid = vector[:100] + unit.sequence + vector[100:]
    reverse_plasmid = vector[:100] + revcomp(unit.sequence) + vector[100:]

    rev = reverse_primer(polya=0)
    a = simulate_pcr(forward_plasmid, FWD_PRIMER, rev)
    b = simulate_pcr(reverse_plasmid, FWD_PRIMER, rev)

    if len(a) != 1:
        problems.append(f"forward orientation gives {len(a)} amplicons, expected 1")
    if len(b) != 1:
        problems.append(f"reverse orientation gives {len(b)} amplicons, expected 1")
    if a and b:
        same = a[0] == b[0] or a[0] == revcomp(b[0])
        if not same:
            problems.append("the two orientations give different amplicons")
        if a[0] != unit.sequence and revcomp(a[0]) != unit.sequence:
            problems.append("the amplicon is not the designed fragment")
    return problems


def screen_handles(forward: str = FWD_HANDLE, reverse: str = REV_HANDLE,
                   context: str = "") -> list[str]:
    """Sanity-check the primer handles."""
    problems = []
    for label, h in (("forward", forward), ("reverse", reverse)):
        if not 0.40 <= gc_fraction(h) <= 0.60:
            problems.append(f"{label} handle GC {gc_fraction(h):.0%} outside 40-60%")
        if re.search(r"(A|C|G|T)\1{3,}", h):
            problems.append(f"{label} handle has a 4+ homopolymer run")
        if h in context or revcomp(h) in context:
            problems.append(f"{label} handle also occurs in the construct body")
    if sum(a != b for a, b in zip(forward, revcomp(reverse))) < 6:
        problems.append("the two handles are near-complementary and will form a primer dimer")
    return problems


__all__ = [
    "T7_CORE", "LEADER", "UTR3", "FWD_HANDLE", "REV_HANDLE",
    "FWD_PRIMER", "reverse_primer",
    "IVTUnit", "make_unit", "simulate_pcr",
    "verify_unit", "verify_orientation_independence", "screen_handles",
]
