"""In-silico validation of the finished transcription units.

The auditors in ``audit/`` answer "is the molecule the one we designed?" -- they
check composition, frame, order and paperwork. They cannot answer "will the
molecule behave?", because that question is about folding and translation, not
about bookkeeping. This module covers the second question with published
algorithms run on the actual sequences.

What is here, and what each thing is actually evidence for:

* **Secondary structure** (ViennaRNA, Lorenz et al. 2011, *Algorithms Mol Biol*
  6:26). MFE and ensemble free energy over the whole transcript, plus the
  unpaired probability of the start codon and the 5' end. Structure over the
  cap-proximal region and the AUG impedes 43S scanning, so this is the one
  structural number with a direct mechanistic link to how much protein a
  construct makes.
* **Codon adaptation index** (Sharp & Li 1987, *NAR* 15:1281). Geometric mean of
  relative codon adaptiveness. A panel-level comparison, not an absolute score:
  what matters is that the constructs match each other.
* **CpG content.** CpG dinucleotides in a transcript are substrate for ZAP-
  mediated decay and contribute to innate sensing. Counted as observed/expected
  so it is not just a restatement of GC.
* **Hydropathy** (Kyte & Doolittle 1982, *JMB* 157:105). Used two ways: to check
  that each signal peptide has the n/h/c architecture a signal peptide needs,
  and that each transmembrane module contains a segment long and hydrophobic
  enough to span a bilayer. This is what the design's whole routing logic rests
  on, so it is worth checking rather than assuming.

**What is deliberately NOT here.** No MHC binding prediction. The class-I
readout rests on SIINFEKL/H-2K(b) and the class-II readout on an I-A(b) core,
both of which are characterised experimentally in the literature; a neural-net
affinity prediction would add a number without adding information, and the
predictors that cover mouse alleles well are not the ones that are easy to
pin to a version. Epitope *integrity* is checked by the auditors instead.
Cleavage-site prediction (SignalP, TargetP) is likewise absent: those are
licence-encumbered and not reproducible from this repo, so the hydropathy
architecture check below is what stands in for them, and it is weaker.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .codon import HUMAN_CODON_USAGE
from .seqops import CODON_TO_AA

#: Kyte & Doolittle 1982 hydropathy scale.
KD: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

#: Residues accepted at the -1 position of a signal peptidase I site. The
#: "(-3,-1) rule": von Heijne 1983, *Eur J Biochem* 133:17.
SPI_MINUS_ONE = frozenset("AGSCT")
#: -3 tolerates the small set plus short aliphatics.
SPI_MINUS_THREE = frozenset("AGSCTVILM")

#: A transmembrane helix must span ~30 A of bilayer, which takes about 19
#: residues of alpha helix.
TM_MIN_LENGTH = 19
#: Mean KD over the window above which a segment is called hydrophobic. 1.6 is
#: the conventional cut for a 19-residue window.
TM_THRESHOLD = 1.6


def _rna(dna: str) -> str:
    return dna.upper().replace("T", "U")


# --- secondary structure ---------------------------------------------------

@dataclass
class FoldMetrics:
    """ViennaRNA numbers for one transcript."""

    length: int
    mfe: float                 #: kcal/mol
    mfe_per_nt: float
    ensemble_free_energy: float
    ensemble_diversity: float
    start_accessibility: float  #: mean unpaired probability over the AUG
    cap_proximal_mfe: float     #: MFE of the first CAP_WINDOW nt alone
    structure: str = field(repr=False, default="")


#: The 5' window that has to be unwound for the 43S complex to load. Kept
#: separate from the whole-transcript fold because a transcript can be
#: comfortably folded overall and still have a hairpin sitting on the cap.
CAP_WINDOW = 70
#: Window centred on the AUG over which unpaired probability is averaged.
START_WINDOW = 15
#: RNAplfold parameters for start-codon accessibility: window and the maximum
#: span a base pair may cover.
PLFOLD_WINDOW, PLFOLD_MAX_SPAN = 80, 40


def start_accessibility(transcript: str, orf_start: int) -> float:
    """Mean unpaired probability across the start-codon window.

    Computed by **local** folding (RNAplfold; Bernhart et al. 2006, *Algorithms
    Mol Biol* 1:3), which bounds how far apart two paired bases may be.

    Global MFE folding was used here first and was wrong. Over a 700-nt
    transcript it reported five different start-codon accessibilities for
    P0-P4, whose first 45 ORF bases are byte-identical -- the variation came
    entirely from long-range pairing to different downstream sequence. That is
    the least reliable part of an equilibrium folding model, and it is also not
    what a scanning 43S complex encounters on a transcript that is being
    translated as it folds. Bounding the span removes the artefact: the local
    measure returns exactly three values for the panel's three distinct 5'
    ends, which is the answer the sequences say it should return.
    """
    import RNA

    seq = _rna(transcript)
    up = RNA.pfl_fold_up(seq, 1, PLFOLD_WINDOW, PLFOLD_MAX_SPAN)
    flank = (START_WINDOW - 3) // 2
    lo, hi = max(1, orf_start - flank), min(len(seq), orf_start + 3 + flank)
    return sum(up[i][1] for i in range(lo, hi)) / (hi - lo)



def fold_metrics(transcript: str, orf_start: int) -> FoldMetrics:
    """Fold one transcript and report the numbers that bear on translation.

    ``orf_start`` is the 0-based index of the A of the AUG within
    ``transcript``. Requires the ``ViennaRNA`` package.
    """
    import RNA

    seq = _rna(transcript)
    fc = RNA.fold_compound(seq)
    structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    _, efe = fc.pf()
    diversity = fc.mean_bp_distance()

    accessibility = start_accessibility(transcript, orf_start)

    cap = RNA.fold(seq[:CAP_WINDOW])[1]

    return FoldMetrics(
        length=len(seq),
        mfe=mfe,
        mfe_per_nt=mfe / len(seq),
        ensemble_free_energy=efe,
        ensemble_diversity=diversity,
        start_accessibility=accessibility,
        cap_proximal_mfe=cap,
        structure=structure,
    )


# --- codon and dinucleotide statistics -------------------------------------

def cai(orf: str, usage: dict[str, float] | None = None) -> float:
    """Codon adaptation index (Sharp & Li 1987), stop codon excluded.

    Relative adaptiveness of a codon is its usage divided by the highest usage
    among its synonyms; CAI is the geometric mean over the coding sequence.
    Met and Trp have one codon each and contribute 1.0 by definition.
    """
    usage = usage or HUMAN_CODON_USAGE
    best: dict[str, float] = {}
    for codon, freq in usage.items():
        aa = CODON_TO_AA[codon]
        best[aa] = max(best.get(aa, 0.0), freq)

    logs = []
    for i in range(0, len(orf) - 2, 3):
        codon = orf[i:i + 3].upper()
        aa = CODON_TO_AA.get(codon)
        if aa in (None, "*"):
            continue
        w = usage.get(codon, 0.0) / best[aa]
        logs.append(math.log(max(w, 1e-6)))
    return math.exp(sum(logs) / len(logs)) if logs else float("nan")


def cpg_stats(seq: str) -> tuple[int, float]:
    """CpG count and observed/expected ratio.

    Observed/expected normalises out base composition, so a GC-rich sequence
    does not automatically look CpG-rich. Values near 1.0 mean CpG occurs at
    the rate composition predicts; mammalian transcripts typically run well
    below 1.0, and CpG above that is a decay and innate-sensing liability.
    """
    seq = seq.upper()
    observed = sum(1 for i in range(len(seq) - 1) if seq[i:i + 2] == "CG")
    c, g = seq.count("C"), seq.count("G")
    expected = c * g / len(seq) if len(seq) else 0.0
    return observed, (observed / expected if expected else float("nan"))


# --- hydropathy ------------------------------------------------------------

def hydropathy(protein: str, window: int = 19) -> list[float]:
    """Kyte-Doolittle profile, one mean per window, indexed by window start."""
    if len(protein) < window:
        return []
    scores = [KD.get(a, 0.0) for a in protein]
    return [sum(scores[i:i + window]) / window
            for i in range(len(protein) - window + 1)]


def tm_segments(protein: str, window: int = TM_MIN_LENGTH,
                threshold: float = TM_THRESHOLD) -> list[tuple[int, int, float]]:
    """Segments hydrophobic enough and long enough to be a TM helix.

    A window-mean hydropathy screen, which is what TMHMM-class tools do before
    their topology model. It over-calls hydrophobic non-TM stretches, so it is
    used here only to confirm that modules *designed* to be transmembrane
    contain such a segment -- a negative is informative, a positive is not proof.
    """
    profile = hydropathy(protein, window)
    out, start = [], None
    for i, value in enumerate(profile):
        if value >= threshold and start is None:
            start = i
        elif value < threshold and start is not None:
            out.append((start, i + window - 1, max(profile[start:i])))
            start = None
    if start is not None:
        out.append((start, len(protein) - 1, max(profile[start:])))

    # Window starts that are disjoint can still cover overlapping residues,
    # because each window is `window` residues wide. Merge them so the output
    # is spans of protein rather than spans of the profile.
    merged: list[tuple[int, int, float]] = []
    for lo, hi, peak in out:
        if merged and lo <= merged[-1][1]:
            prev = merged[-1]
            merged[-1] = (prev[0], max(prev[1], hi), max(prev[2], peak))
        else:
            merged.append((lo, hi, peak))
    return merged


@dataclass
class SignalPeptideCall:
    """Whether a peptide has the architecture signal peptidase I needs."""

    length: int
    n_region_charge: float
    h_region_start: int
    h_region_length: int
    h_region_max: float
    minus_one: str
    minus_three: str
    problems: list[str]

    @property
    def ok(self) -> bool:
        return not self.problems


def signal_peptide_call(peptide: str) -> SignalPeptideCall:
    """Check a signal peptide against the classical n/h/c architecture.

    A signal peptide is a positively charged n-region, a hydrophobic h-region of
    roughly 7-15 residues, and a c-region ending in small residues at -1 and -3.
    This is the von Heijne description that SignalP's early versions encoded;
    it is a much weaker test than a trained predictor and will not find a
    cleavage site on its own, but it does catch a peptide that has stopped being
    a signal peptide at all -- which is the failure this panel needs to exclude.
    """
    problems: list[str] = []
    h_profile = hydropathy(peptide, window=7)
    if not h_profile:
        return SignalPeptideCall(len(peptide), 0.0, -1, 0, float("nan"), "", "",
                                 ["too short to have an h-region"])

    h_start = max(range(len(h_profile)), key=lambda i: h_profile[i])
    h_max = h_profile[h_start]
    h_len = sum(1 for v in h_profile if v >= 1.6)

    n_region = peptide[:h_start]
    charge = sum(1 for a in n_region if a in "KR") - sum(
        1 for a in n_region if a in "DE")

    minus_one = peptide[-1] if peptide else ""
    minus_three = peptide[-3] if len(peptide) >= 3 else ""

    if h_max < 1.6:
        problems.append(f"no hydrophobic h-region (max window KD {h_max:.2f})")
    if h_len < 7:
        problems.append(f"h-region only {h_len} residues, expected >= 7")
    if charge < 0:
        problems.append(f"n-region net charge {charge:+.0f}, expected >= 0")
    if minus_one not in SPI_MINUS_ONE:
        problems.append(
            f"-1 residue {minus_one!r} is not small; signal peptidase I needs "
            f"one of {''.join(sorted(SPI_MINUS_ONE))}")
    if minus_three not in SPI_MINUS_THREE:
        problems.append(
            f"-3 residue {minus_three!r} is outside "
            f"{''.join(sorted(SPI_MINUS_THREE))}")

    return SignalPeptideCall(len(peptide), charge, h_start, h_len, h_max,
                             minus_one, minus_three, problems)
