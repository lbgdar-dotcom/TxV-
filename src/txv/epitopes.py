"""Multi-antigen cassette assembly and junctional-epitope control.

A "string-of-beads" cassette concatenates antigen segments into a single open
reading frame. The characteristic failure mode is the **junctional neoepitope**:
the fusion of the C-terminus of one bead to the N-terminus of the next creates
peptides that exist nowhere in the patient's proteome or tumour. Those peptides
compete for presentation, can dominate the response, and are not a therapeutic
target. Cassette design is therefore mostly (a) choosing linkers and (b)
choosing an *order* that minimises junctional binders.

This module gives you:

* :class:`Antigen` / :class:`Cassette` -- the payload data model.
* :func:`scan_junctions` -- enumerate every peptide that spans a junction.
* :class:`JunctionScorer` -- the plug point for a real MHC binding predictor.
* :class:`AnchorMotifScorer` -- a deliberately crude built-in proxy so the
  pipeline runs without a predictor installed.
* :func:`optimize_order` -- 2-opt reordering of beads against a scorer.

.. warning::
   :class:`AnchorMotifScorer` is **not** an MHC binding predictor. It matches
   position-2 and C-terminal anchor residues for a handful of supertypes. Use
   it to smoke-test a pipeline; use NetMHCpan, MHCflurry or your in-house
   predictor (wrapped as a :class:`JunctionScorer`) for anything that informs a
   real design.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Callable, Protocol, Sequence

AA = set("ACDEFGHIKLMNPQRSTVWY")

# Class I peptides are overwhelmingly 8-11mers; class II are longer and have a
# more permissive, open-ended groove.
CLASS_I_LENGTHS = (8, 9, 10, 11)
CLASS_II_LENGTHS = (15,)


@dataclass(frozen=True)
class Antigen:
    """One bead of the cassette.

    ``sequence`` is the amino-acid payload: a minimal epitope, a longer
    "synthetic long peptide" style segment carrying a neoepitope in its native
    flanks, or a whole tumour-associated antigen.
    """

    name: str
    sequence: str
    kind: str = "neoepitope"   # neoepitope | tumor_associated | full_length | helper
    hla: tuple[str, ...] = ()
    #: Index of the mutated residue within ``sequence`` (0-based), if known.
    mutation_offset: int | None = None
    #: Exact DNA for this bead, when its encoding was chosen deliberately and
    #: must not be re-derived by the optimiser.
    pinned_dna: str | None = None
    note: str = ""

    def __post_init__(self) -> None:
        seq = self.sequence.strip().upper()
        object.__setattr__(self, "sequence", seq)
        bad = sorted(set(seq) - AA)
        if bad:
            raise ValueError(f"antigen {self.name!r}: bad residues {bad}")
        if not seq:
            raise ValueError(f"antigen {self.name!r}: empty sequence")
        if self.mutation_offset is not None and not 0 <= self.mutation_offset < len(seq):
            raise ValueError(f"antigen {self.name!r}: mutation_offset out of range")
        if self.pinned_dna is not None:
            from .seqops import clean, translate

            pinned = clean(self.pinned_dna)
            object.__setattr__(self, "pinned_dna", pinned)
            if translate(pinned, stop_at_stop=False) != seq:
                raise ValueError(
                    f"antigen {self.name!r}: pinned_dna encodes "
                    f"{translate(pinned, stop_at_stop=False)!r}, not {seq!r}"
                )

    def __len__(self) -> int:
        return len(self.sequence)


@dataclass
class Junction:
    """A boundary between two beads, with the linker that sits in it."""

    left: str
    right: str
    linker: str
    index: int

    @property
    def context(self) -> str:
        """The stretch of sequence a junctional peptide can be drawn from."""
        return self.left + self.linker + self.right


@dataclass
class Cassette:
    """An ordered polypeptide payload plus the linkers between beads."""

    antigens: list[Antigen]
    linker: str = "GGSGGGGSGG"
    #: Optional per-junction linker overrides, keyed by junction index.
    linker_overrides: dict[int, str] = field(default_factory=dict)

    def linker_at(self, index: int) -> str:
        return self.linker_overrides.get(index, self.linker)

    def sequence(self) -> str:
        parts: list[str] = []
        for i, antigen in enumerate(self.antigens):
            if i:
                parts.append(self.linker_at(i - 1))
            parts.append(antigen.sequence)
        return "".join(parts)

    def layout(self) -> list[tuple[str, int, int, str]]:
        """(name, start, end, kind) spans in the cassette, 0-based half-open."""
        spans: list[tuple[str, int, int, str]] = []
        pos = 0
        for i, antigen in enumerate(self.antigens):
            if i:
                link = self.linker_at(i - 1)
                spans.append((f"linker_{i - 1}", pos, pos + len(link), "linker"))
                pos += len(link)
            spans.append((antigen.name, pos, pos + len(antigen), antigen.kind))
            pos += len(antigen)
        return spans

    def junctions(self) -> list[Junction]:
        return [
            Junction(
                left=self.antigens[i].sequence,
                right=self.antigens[i + 1].sequence,
                linker=self.linker_at(i),
                index=i,
            )
            for i in range(len(self.antigens) - 1)
        ]

    def __len__(self) -> int:
        return len(self.sequence())


class JunctionScorer(Protocol):
    """Scores a peptide's risk of being a presented junctional epitope.

    Return a value in [0, 1]; higher means more likely to be presented. Wrap a
    real predictor like this::

        class NetMHCpanScorer:
            def __init__(self, alleles): self.alleles = alleles
            def __call__(self, peptide: str) -> float:
                # run predictor, convert %rank to risk, e.g. 1.0 if rank < 0.5
                ...
    """

    def __call__(self, peptide: str) -> float: ...


class AnchorMotifScorer:
    """Crude anchor-residue proxy for class I presentation. **Not a predictor.**

    Scores the fraction of configured supertypes whose position-2 and
    C-terminal anchor preferences a peptide satisfies. Its only virtues are
    that it needs no model and is monotone in the thing you care about.
    """

    #: supertype -> (allowed P2 residues, allowed C-terminal residues)
    SUPERTYPES: dict[str, tuple[str, str]] = {
        "A2": ("LMIVQAT", "VLIMAT"),
        "A3": ("LMIVQAST", "KRY"),
        "A24": ("YF", "FLIW"),
        "B7": ("P", "LMIVFA"),
        "B27": ("R", "KRLFY"),
        "B44": ("E", "FWYLIM"),
    }

    def __init__(self, supertypes: Sequence[str] | None = None,
                 lengths: Sequence[int] = CLASS_I_LENGTHS) -> None:
        keys = supertypes or list(self.SUPERTYPES)
        unknown = set(keys) - set(self.SUPERTYPES)
        if unknown:
            raise ValueError(f"unknown supertypes: {sorted(unknown)}")
        self.supertypes = {k: self.SUPERTYPES[k] for k in keys}
        self.lengths = tuple(lengths)

    def __call__(self, peptide: str) -> float:
        if len(peptide) not in self.lengths or len(peptide) < 3:
            return 0.0
        hits = sum(
            1
            for p2, pc in self.supertypes.values()
            if peptide[1] in p2 and peptide[-1] in pc
        )
        return hits / len(self.supertypes)


@dataclass
class JunctionHit:
    junction_index: int
    peptide: str
    score: float
    start_in_context: int


def scan_junctions(
    cassette: Cassette,
    scorer: JunctionScorer | None = None,
    lengths: Sequence[int] = CLASS_I_LENGTHS,
    threshold: float = 0.15,
) -> list[JunctionHit]:
    """Enumerate peptides that genuinely span a junction and score them.

    A peptide counts as junctional if it is not contained in a single bead.
    Peptides lying wholly inside a bead are that bead's own epitopes and must
    not be penalised; peptides lying wholly inside the linker encode no novel
    antigen. Everything else -- bead into linker, linker into bead, or bead
    straight into bead when the linker is short -- is sequence that exists
    nowhere outside this construct, and is what you are screening for.

    The default ``threshold`` of 0.15 means "flag a peptide that matches the
    anchor pattern of at least one supertype" for the six-supertype built-in
    scorer. Set it to suit the dynamic range of whatever scorer you pass.
    """
    scorer = scorer or AnchorMotifScorer(lengths=lengths)
    hits: list[JunctionHit] = []
    for junction in cassette.junctions():
        context = junction.context
        left_end = len(junction.left)
        right_start = left_end + len(junction.linker)
        for k in lengths:
            for start in range(0, max(0, len(context) - k + 1)):
                end = start + k
                if end <= left_end or start >= right_start:
                    continue  # wholly inside one bead: that bead's own epitope
                if start >= left_end and end <= right_start:
                    continue  # wholly inside the linker: not a novel antigen
                peptide = context[start:end]
                score = scorer(peptide)
                if score >= threshold:
                    hits.append(JunctionHit(junction.index, peptide, score, start))
    return sorted(hits, key=lambda h: (-h.score, h.junction_index))


def junction_risk(
    cassette: Cassette,
    scorer: JunctionScorer | None = None,
    lengths: Sequence[int] = CLASS_I_LENGTHS,
    threshold: float = 0.15,
) -> float:
    """Total junctional risk: the summed score of all supra-threshold peptides."""
    return sum(h.score for h in scan_junctions(cassette, scorer, lengths, threshold))


def optimize_order(
    cassette: Cassette,
    scorer: JunctionScorer | None = None,
    lengths: Sequence[int] = CLASS_I_LENGTHS,
    threshold: float = 0.15,
    fixed_first: bool = False,
    fixed_last: bool = False,
    max_passes: int = 40,
) -> tuple[Cassette, float, float]:
    """Reorder beads to minimise junctional risk.

    Exhaustive up to 7 beads; 2-opt with repeated passes beyond that (the
    objective is a path cost over an asymmetric junction matrix, so this is a
    small open TSP and 2-opt gets very close on realistic cassette sizes).

    Returns ``(cassette, risk_before, risk_after)``. ``fixed_first`` /
    ``fixed_last`` pin the terminal beads, which you want when the first bead
    must abut the signal peptide or the last must abut a trafficking domain.
    """
    scorer = scorer or AnchorMotifScorer(lengths=lengths)

    def risk(order: list[Antigen]) -> float:
        trial = Cassette(order, cassette.linker, dict(cassette.linker_overrides))
        return junction_risk(trial, scorer, lengths, threshold)

    original = list(cassette.antigens)
    before = risk(original)
    if len(original) < 3:
        return cassette, before, before

    lo = 1 if fixed_first else 0
    hi = len(original) - 1 if fixed_last else len(original)

    if hi - lo <= 7:
        best_order, best = original, before
        for perm in itertools.permutations(original[lo:hi]):
            candidate = original[:lo] + list(perm) + original[hi:]
            score = risk(candidate)
            if score < best:
                best_order, best = candidate, score
    else:
        best_order, best = list(original), before
        for _ in range(max_passes):
            improved = False
            for i in range(lo, hi - 1):
                for j in range(i + 1, hi):
                    candidate = best_order.copy()
                    candidate[i : j + 1] = reversed(candidate[i : j + 1])
                    score = risk(candidate)
                    if score < best - 1e-12:
                        best_order, best, improved = candidate, score, True
            if not improved:
                break

    return (
        Cassette(best_order, cassette.linker, dict(cassette.linker_overrides)),
        before,
        best,
    )


def choose_linker(
    cassette: Cassette,
    candidates: Sequence[str],
    scorer: JunctionScorer | None = None,
    lengths: Sequence[int] = CLASS_I_LENGTHS,
    threshold: float = 0.15,
) -> tuple[str, dict[str, float]]:
    """Pick the single uniform linker with the lowest junctional risk.

    Returns the winner and the full risk table, so a reviewer can see what the
    alternatives cost rather than just the answer.
    """
    table: dict[str, float] = {}
    for linker in candidates:
        trial = Cassette(list(cassette.antigens), linker, dict(cassette.linker_overrides))
        table[linker] = junction_risk(trial, scorer, lengths, threshold)
    winner = min(table, key=lambda k: (table[k], len(k)))
    return winner, table


def dedupe_antigens(antigens: Sequence[Antigen]) -> tuple[list[Antigen], list[str]]:
    """Drop exact-duplicate and fully contained payload sequences.

    Returns the kept antigens and a list of human-readable drop reasons.
    """
    kept: list[Antigen] = []
    dropped: list[str] = []
    for antigen in sorted(antigens, key=lambda a: -len(a)):
        container = next((k for k in kept if antigen.sequence in k.sequence), None)
        if container is not None:
            reason = "identical to" if len(antigen) == len(container) else "contained in"
            dropped.append(f"{antigen.name}: {reason} {container.name}")
            continue
        kept.append(antigen)
    order = {a.name: i for i, a in enumerate(antigens)}
    kept.sort(key=lambda a: order[a.name])
    return kept, dropped


__all__ = [
    "Antigen", "Cassette", "Junction", "JunctionHit", "JunctionScorer",
    "AnchorMotifScorer", "CLASS_I_LENGTHS", "CLASS_II_LENGTHS",
    "scan_junctions", "junction_risk", "optimize_order", "choose_linker",
    "dedupe_antigens",
]
