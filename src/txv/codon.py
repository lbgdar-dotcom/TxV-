"""Codon optimisation for IVT mRNA, as a constrained sequence search.

Naive "pick the most frequent codon everywhere" optimisation is a poor fit for
mRNA drug substance, because the objectives that matter pull against each other:

* **Translation efficiency** wants common, well-charged codons.
* **Innate immunogenicity and IVT fidelity** want *low uridine*. Uridine content
  drives RIG-I/TLR sensing, and long U runs are where T7 slips and produces
  aborted or +1 products. Modified nucleotides (m1Ψ) reduce but do not remove
  the incentive.
* **Manufacturability** wants a GC band that is high enough for template
  stability but not so high that it forms structure the polymerase stalls on,
  plus no restriction sites used for linearisation, no cryptic splice sites, no
  long homopolymers, and no long direct repeats (which recombine in the
  plasmid and confound sequencing).

So this is a constrained optimisation over synonymous codons, and it is solved
here with a **beam search** over the coding sequence, left to right. The score
of a partial sequence is a weighted sum of the terms above, evaluated
incrementally per codon, plus hard-ish penalties for forbidden motifs that
straddle the codon boundary just written. A final repair pass re-synonymises
locally around any forbidden motif the beam did not manage to avoid.

Weights live in :class:`OptimizerConfig` so the trade-off is an explicit,
reviewable parameter of a design rather than a property of the tool.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

from .seqops import (
    AA_TO_CODONS,
    CODON_TO_AA,
    clean,
    gc_fraction,
    iupac_regex,
    revcomp_iupac,
    translate,
)

# Relative synonymous codon usage in *Homo sapiens*: fraction of that amino
# acid's codons, rounded to two decimals from the standard genome-wide tables.
# Override with your own counts via OptimizerConfig.usage if you have a tissue-
# or expression-system-specific table.
HUMAN_CODON_USAGE: dict[str, float] = {
    "GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11,
    "CGT": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21, "AGA": 0.20, "AGG": 0.20,
    "AAT": 0.46, "AAC": 0.54,
    "GAT": 0.46, "GAC": 0.54,
    "TGT": 0.45, "TGC": 0.55,
    "CAA": 0.25, "CAG": 0.75,
    "GAA": 0.42, "GAG": 0.58,
    "GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
    "CAT": 0.41, "CAC": 0.59,
    "ATT": 0.36, "ATC": 0.48, "ATA": 0.16,
    "TTA": 0.07, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.41,
    "AAA": 0.42, "AAG": 0.58,
    "ATG": 1.00,
    "TTT": 0.45, "TTC": 0.55,
    "CCT": 0.28, "CCC": 0.33, "CCA": 0.27, "CCG": 0.11,
    "TCT": 0.18, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06, "AGT": 0.15, "AGC": 0.24,
    "ACT": 0.24, "ACC": 0.36, "ACA": 0.28, "ACG": 0.12,
    "TGG": 1.00,
    "TAT": 0.43, "TAC": 0.57,
    "GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47,
    "TAA": 0.30, "TAG": 0.24, "TGA": 0.47,
}

# Motifs excluded from coding regions by default. IUPAC codes allowed; each is
# screened on both strands.
DEFAULT_FORBIDDEN: tuple[str, ...] = (
    # Type IIS / common cloning and linearisation sites.
    "GCTCTTC",   # BspQI / SapI
    "GGTCTC",    # BsaI
    "CGTCTC",    # BsmBI / Esp3I
    "GAATTC",    # EcoRI
    "GGATCC",    # BamHI
    "AAGCTT",    # HindIII
    "CTCGAG",    # XhoI
    "GCGGCCGC",  # NotI
    # Cryptic splice donor consensus (exon|GTAAGT). Splicing is not supposed to
    # happen to a cytoplasmically delivered transcript, but these motifs are
    # screened out because the DNA template is propagated and sometimes
    # expressed in mammalian cells during characterisation.
    "GGTAAGT",
    "GGTGAGT",
    # Branch-point-like / 3' splice acceptor core in a pyrimidine context.
    "TTTTTTTTTTNCAG",
    # T7 class II transcription pause / termination-like element.
    "ATCTGTT",
    # Poly(A) signal -- avoid premature polyadenylation of the DNA template.
    "AATAAA",
)


@dataclass
class OptimizerConfig:
    """Objective weights and constraints. All weights are in arbitrary units;
    only their ratios matter."""

    usage: dict[str, float] = field(default_factory=lambda: dict(HUMAN_CODON_USAGE))
    #: Codons below this share of their amino acid are dropped unless the amino
    #: acid has no alternative. Removes the slowest-decoded codons outright.
    min_codon_usage: float = 0.10
    #: Target GC fraction of the coding sequence.
    target_gc: float = 0.60
    gc_window: int = 60
    #: Penalty weight on squared deviation of windowed GC from target.
    w_gc: float = 40.0
    #: Penalty per uridine (thymine) placed. Raise to push uridine depletion
    #: harder at the cost of codon adaptation; 0.0 disables U depletion.
    w_uridine: float = 0.45
    #: Weight on codon adaptiveness, -log(usage / max usage for that residue).
    w_usage: float = 1.0
    #: Flat penalty per forbidden-motif occurrence created.
    w_forbidden: float = 500.0
    #: Runs of a single base at or above this length are penalised.
    max_homopolymer: int = 5
    w_homopolymer: float = 30.0
    #: Direct repeats of this length seen twice are penalised (recombination
    #: and sequencing-assembly risk in the plasmid template).
    repeat_k: int = 14
    w_repeat: float = 25.0
    forbidden: tuple[str, ...] = DEFAULT_FORBIDDEN
    beam_width: int = 40
    seed: int = 0

    def __post_init__(self) -> None:
        if self.beam_width < 1:
            raise ValueError("beam_width must be >= 1")
        if not 0.0 <= self.target_gc <= 1.0:
            raise ValueError("target_gc must be a fraction in [0, 1]")


@dataclass
class OptimizationResult:
    dna: str
    protein: str
    cai: float
    gc: float
    uridine_fraction: float
    residual_forbidden: list[tuple[str, int]]
    n_repair_edits: int

    @property
    def ok(self) -> bool:
        return not self.residual_forbidden


def _relative_adaptiveness(usage: dict[str, float]) -> dict[str, float]:
    """w_i = f_i / max(f) within each amino acid's synonymous family (Sharp & Li)."""
    weights: dict[str, float] = {}
    for aa, syn in AA_TO_CODONS.items():
        best = max(usage.get(c, 0.0) for c in syn) or 1.0
        for codon in syn:
            weights[codon] = max(usage.get(codon, 0.0) / best, 1e-6)
    return weights


def codon_adaptation_index(dna: str, usage: dict[str, float] | None = None) -> float:
    """Geometric mean relative adaptiveness over the sense codons of a CDS."""
    weights = _relative_adaptiveness(usage or HUMAN_CODON_USAGE)
    sense = [
        dna[i : i + 3]
        for i in range(0, len(dna) - len(dna) % 3, 3)
        if CODON_TO_AA[dna[i : i + 3]] != "*"
    ]
    sense = [c for c in sense if CODON_TO_AA[c] not in ("M", "W")]
    if not sense:
        return 0.0
    return math.exp(sum(math.log(weights[c]) for c in sense) / len(sense))


def uridine_floor(protein: str, usage: dict[str, float] | None = None,
                  min_codon_usage: float = 0.0) -> tuple[int, int]:
    """Minimum achievable uridines in a CDS encoding ``protein``.

    Many residues have no U-free codon -- Phe, Tyr, Cys, Trp and Ile cannot
    avoid one, and Leu/Ser/Val can only sometimes. So uridine depletion has a
    hard floor set by the amino-acid sequence, and an optimiser that appears to
    "stop responding" to a higher uridine weight has simply reached it.

    Returns ``(floor_uridines, codons)``; divide by ``codons * 3`` for the
    floor as a fraction of the CDS.
    """
    usage = usage or HUMAN_CODON_USAGE
    total = 0
    for aa in protein.rstrip("*"):
        options = [c for c in AA_TO_CODONS[aa]
                   if usage.get(c, 0.0) >= min_codon_usage] or AA_TO_CODONS[aa]
        total += min(c.count("T") for c in options)
    return total, len(protein.rstrip("*"))


class CodonOptimizer:
    """Beam-search codon optimiser."""

    def __init__(self, config: OptimizerConfig | None = None) -> None:
        self.config = config or OptimizerConfig()
        self._weights = _relative_adaptiveness(self.config.usage)
        self._forbidden_res = [
            (motif, re.compile(f"(?={iupac_regex(motif)})"))
            for motif in self.config.forbidden
        ]
        self._forbidden_rc_res = [
            (motif, re.compile(f"(?={iupac_regex(revcomp_iupac(motif))})"))
            for motif in self.config.forbidden
        ]
        self._max_motif = max((len(m) for m in self.config.forbidden), default=0)

    # -- public API ---------------------------------------------------------
    def optimize(
        self,
        protein: str,
        add_stop: str | None = "TGA",
        pinned: dict[int, str] | None = None,
    ) -> OptimizationResult:
        """Optimise a protein sequence into a coding DNA sequence.

        ``protein`` may end in ``*``; otherwise ``add_stop`` (a stop codon, or
        ``None`` for a translational fusion) is appended.

        ``pinned`` maps a residue index to the exact codon to use there, for
        elements whose encoding was chosen deliberately and must not be
        re-derived. Pinned codons still supply context to the search, so the
        rest of the sequence is optimised *around* them rather than in
        ignorance of them -- which is the whole point of doing this in one pass.
        """
        protein = protein.strip().upper()
        if protein.endswith("*"):
            protein, add_stop = protein[:-1], add_stop or "TGA"
        if not protein:
            raise ValueError("empty protein sequence")
        unknown = set(protein) - set(AA_TO_CODONS)
        if unknown:
            raise ValueError(f"unknown residues: {sorted(unknown)}")

        pinned = dict(pinned or {})
        for index, codon in pinned.items():
            if not 0 <= index < len(protein):
                raise ValueError(f"pinned index {index} is outside the protein")
            if CODON_TO_AA[codon] != protein[index]:
                raise ValueError(
                    f"pinned codon {codon} at {index} encodes "
                    f"{CODON_TO_AA[codon]}, not {protein[index]}"
                )

        dna = self._beam_search(protein, pinned)
        if add_stop:
            dna += clean(add_stop)
        dna, edits = self._repair(dna, protein, pinned)

        return OptimizationResult(
            dna=dna,
            protein=translate(dna, stop_at_stop=False).rstrip("*"),
            cai=codon_adaptation_index(dna, self.config.usage),
            gc=gc_fraction(dna),
            uridine_fraction=dna.count("T") / len(dna),
            residual_forbidden=self.find_forbidden(dna),
            n_repair_edits=edits,
        )

    def find_forbidden(self, dna: str) -> list[tuple[str, int]]:
        """All forbidden-motif occurrences, both strands, as (motif, position)."""
        hits: list[tuple[str, int]] = []
        for motif, rx in self._forbidden_res:
            hits += [(motif, m.start()) for m in rx.finditer(dna)]
        for motif, rx in self._forbidden_rc_res:
            hits += [(f"{motif}(rc)", m.start()) for m in rx.finditer(dna)]
        return sorted(set(hits), key=lambda h: h[1])

    def choices_for(self, aa: str) -> list[str]:
        """Synonymous codons for ``aa`` after the rare-codon cut."""
        syn = AA_TO_CODONS[aa]
        kept = [c for c in syn if self.config.usage.get(c, 0.0) >= self.config.min_codon_usage]
        return sorted(kept or syn)

    # -- search -------------------------------------------------------------
    def _beam_search(self, protein: str, pinned: dict[int, str] | None = None) -> str:
        cfg = self.config
        pinned = pinned or {}
        # Beam entries: (cost, sequence). Kept small; sequences are short enough
        # that carrying the full string is cheaper than reconstructing it.
        beam: list[tuple[float, str]] = [(0.0, "")]
        for index, aa in enumerate(protein):
            options = [pinned[index]] if index in pinned else self.choices_for(aa)
            candidates: list[tuple[float, str]] = []
            for cost, seq in beam:
                for codon in options:
                    candidates.append((cost + self._step_cost(seq, codon), seq + codon))
            candidates.sort(key=lambda c: c[0])
            # Deduplicate on the suffix that can still affect future costs.
            seen: set[str] = set()
            beam = []
            ctx = max(self._max_motif, cfg.repeat_k, cfg.gc_window, cfg.max_homopolymer)
            for cost, seq in candidates:
                key = seq[-ctx:]
                if key in seen:
                    continue
                seen.add(key)
                beam.append((cost, seq))
                if len(beam) >= cfg.beam_width:
                    break
        return min(beam, key=lambda c: c[0])[1]

    def _step_cost(self, prefix: str, codon: str) -> float:
        cfg = self.config
        cost = cfg.w_usage * -math.log(self._weights[codon])
        cost += cfg.w_uridine * codon.count("T")

        seq = prefix + codon
        window = seq[-cfg.gc_window :]
        if len(window) >= cfg.gc_window:
            cost += cfg.w_gc * (gc_fraction(window) - cfg.target_gc) ** 2

        # Only motifs overlapping the newly written codon can be new.
        tail = seq[-(self._max_motif + 2) :]
        for motif, rx in self._forbidden_res:
            if rx.search(tail):
                cost += cfg.w_forbidden
        for motif, rx in self._forbidden_rc_res:
            if rx.search(tail):
                cost += cfg.w_forbidden

        run_tail = seq[-(cfg.max_homopolymer + 3) :]
        if re.search(r"(A|C|G|T)\1{%d,}" % (cfg.max_homopolymer - 1), run_tail):
            cost += cfg.w_homopolymer

        if len(seq) >= cfg.repeat_k:
            kmer = seq[-cfg.repeat_k :]
            if kmer in seq[: -cfg.repeat_k + 1]:
                cost += cfg.w_repeat
        return cost

    # -- repair -------------------------------------------------------------
    def _repair(self, dna: str, protein: str,
                pinned: dict[int, str] | None = None) -> tuple[str, int]:
        """Remove residual forbidden motifs by local re-synonymisation.

        For each remaining hit, try every synonymous substitution at each codon
        the motif overlaps and keep the first that removes the hit without
        creating a new one. Protein sequence is invariant by construction.
        """
        pinned = pinned or {}
        edits = 0
        for _ in range(8):
            hits = self.find_forbidden(dna)
            if not hits:
                break
            motif, pos = hits[0]
            span = len(motif.removesuffix("(rc)"))
            first = max(0, pos // 3)
            last = min(len(dna) // 3 - 1, (pos + span) // 3)
            improved = False
            baseline = len(hits)
            for ci in range(first, last + 1):
                if ci in pinned:
                    continue  # deliberately chosen; not ours to re-derive
                aa = CODON_TO_AA[dna[ci * 3 : ci * 3 + 3]]
                if aa == "*":
                    continue
                for alt in self.choices_for(aa):
                    if alt == dna[ci * 3 : ci * 3 + 3]:
                        continue
                    trial = dna[: ci * 3] + alt + dna[ci * 3 + 3 :]
                    if len(self.find_forbidden(trial)) < baseline:
                        dna, edits, improved = trial, edits + 1, True
                        break
                if improved:
                    break
            if not improved:
                break  # motif is unavoidable for this peptide; QC will flag it
        assert translate(dna, stop_at_stop=False).rstrip("*") == protein.rstrip("*"), (
            "repair changed the encoded protein"
        )
        return dna, edits


__all__ = [
    "CodonOptimizer", "OptimizerConfig", "OptimizationResult",
    "HUMAN_CODON_USAGE", "DEFAULT_FORBIDDEN", "codon_adaptation_index",
    "uridine_floor",
]
