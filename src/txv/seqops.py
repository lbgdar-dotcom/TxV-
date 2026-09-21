"""Low-level nucleotide and translation helpers.

Everything in this package works on DNA alphabet (T, not U) because the
deliverable of an IVT programme is a *DNA template* -- a plasmid or a PCR
amplicon -- that a polymerase transcribes. Conversion to RNA is presentation
only (:func:`to_rna`).
"""

from __future__ import annotations

import re
from typing import Iterable

DNA_ALPHABET = frozenset("ACGT")

CODON_TO_AA = {}
_BASES = "TCAG"
_AA_ORDER = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)
for _i, _b1 in enumerate(_BASES):
    for _j, _b2 in enumerate(_BASES):
        for _k, _b3 in enumerate(_BASES):
            CODON_TO_AA[_b1 + _b2 + _b3] = _AA_ORDER[_i * 16 + _j * 4 + _k]

AA_TO_CODONS: dict[str, list[str]] = {}
for _codon, _aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)

STOP_CODONS = ("TAA", "TAG", "TGA")
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def clean(seq: str) -> str:
    """Uppercase and strip whitespace; raise on non-DNA characters."""
    out = re.sub(r"\s+", "", seq).upper().replace("U", "T")
    bad = sorted(set(out) - DNA_ALPHABET - {"N"})
    if bad:
        raise ValueError(f"non-DNA characters in sequence: {bad}")
    return out


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def base_fraction(seq: str, base: str) -> float:
    if not seq:
        return 0.0
    return seq.count(base) / len(seq)


def gc_windows(seq: str, window: int = 100, step: int = 10) -> list[tuple[int, float]]:
    """GC fraction in a sliding window. Returns (start_index, gc) pairs."""
    if len(seq) < window:
        return [(0, gc_fraction(seq))] if seq else []
    return [(i, gc_fraction(seq[i : i + window])) for i in range(0, len(seq) - window + 1, step)]


def longest_homopolymer(seq: str) -> tuple[str, int, int]:
    """Return (base, length, start) of the longest single-base run."""
    best = ("", 0, 0)
    i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if j - i > best[1]:
            best = (seq[i], j - i, i)
        i = j
    return best


def homopolymer_runs(seq: str, min_length: int) -> list[tuple[str, int, int]]:
    """All single-base runs of at least ``min_length``, as (base, length, start)."""
    return [
        (m.group(0)[0], len(m.group(0)), m.start())
        for m in re.finditer(r"(A+|C+|G+|T+)", seq)
        if len(m.group(0)) >= min_length
    ]


def find_all(seq: str, motif: str, both_strands: bool = False) -> list[int]:
    """0-based start positions of ``motif`` (IUPAC-aware) in ``seq``."""
    pattern = iupac_regex(motif)
    hits = [m.start() for m in re.finditer(f"(?={pattern})", seq)]
    if both_strands:
        rc = iupac_regex(revcomp_iupac(motif))
        hits += [m.start() for m in re.finditer(f"(?={rc})", seq)]
    return sorted(set(hits))


_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}
_IUPAC_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def iupac_regex(motif: str) -> str:
    return "".join(_IUPAC[c] for c in motif.upper())


def revcomp_iupac(motif: str) -> str:
    return motif.upper().translate(_IUPAC_COMPLEMENT)[::-1]


def translate(seq: str, stop_at_stop: bool = True) -> str:
    """Translate frame 1. ``*`` marks a stop codon."""
    out = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa = CODON_TO_AA[seq[i : i + 3]]
        if aa == "*" and stop_at_stop:
            out.append("*")
            break
        out.append(aa)
    return "".join(out)


def codons(seq: str) -> list[str]:
    return [seq[i : i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]


def has_internal_stop(cds: str) -> bool:
    """True if a stop codon occurs before the final codon of an in-frame CDS.

    Translation must not halt at the first stop here -- the whole point is to
    look past it.
    """
    return "*" in translate(cds, stop_at_stop=False)[:-1]


def kmers(seq: str, k: int) -> Iterable[tuple[int, str]]:
    for i in range(0, max(0, len(seq) - k + 1)):
        yield i, seq[i : i + k]
