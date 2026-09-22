"""Quality control for an assembled construct.

Each check returns a :class:`Check` with a severity. The intent is that ``FAIL``
means "this will not work or will not manufacture", ``WARN`` means "a human
should look at this and may well accept it", and ``PASS`` is recorded too so
the report doubles as evidence the check ran.

Two of these checks are explicitly *proxies* and say so in their output: the
5'-end structure check (stem length and GC, not a folding free energy) and the
junctional-epitope check (whatever scorer you supplied). Nothing here computes
a thermodynamic ensemble; if you need one, run ViennaRNA over
``construct.transcript`` and add it as a custom check.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Sequence

from .codon import CodonOptimizer, OptimizerConfig
from .constructs import Construct
from .epitopes import CLASS_I_LENGTHS, JunctionScorer, scan_junctions
from .seqops import (
    find_all,
    gc_fraction,
    gc_windows,
    homopolymer_runs,
    revcomp,
    translate,
)


class Severity(str, Enum):
    PASS = "pass"
    WARN = "warn"
    FAIL = "fail"


@dataclass
class Check:
    name: str
    severity: Severity
    message: str
    detail: dict[str, object] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.severity is not Severity.FAIL


@dataclass
class QCReport:
    construct: str
    checks: list[Check]

    @property
    def failures(self) -> list[Check]:
        return [c for c in self.checks if c.severity is Severity.FAIL]

    @property
    def warnings(self) -> list[Check]:
        return [c for c in self.checks if c.severity is Severity.WARN]

    @property
    def passed(self) -> bool:
        return not self.failures

    def to_text(self, show_pass: bool = True) -> str:
        glyph = {Severity.PASS: "PASS", Severity.WARN: "WARN", Severity.FAIL: "FAIL"}
        lines = [f"QC report: {self.construct}"]
        for check in self.checks:
            if check.severity is Severity.PASS and not show_pass:
                continue
            lines.append(f"  [{glyph[check.severity]}] {check.name}: {check.message}")
        lines.append(
            f"  -> {len(self.failures)} failure(s), {len(self.warnings)} warning(s)"
        )
        return "\n".join(lines)

    def to_dict(self) -> dict[str, object]:
        return {
            "construct": self.construct,
            "passed": self.passed,
            "checks": [
                {
                    "name": c.name,
                    "severity": c.severity.value,
                    "message": c.message,
                    "detail": c.detail,
                }
                for c in self.checks
            ],
        }


@dataclass
class QCThresholds:
    gc_min: float = 0.45
    gc_max: float = 0.70
    window_gc_min: float = 0.30
    window_gc_max: float = 0.80
    gc_window: int = 100
    uridine_max: float = 0.25
    max_homopolymer: int = 9
    #: Poly(A) is a homopolymer by design, so it is excluded from that check.
    cap_proximal_nt: int = 40
    cap_proximal_gc_max: float = 0.65
    cap_proximal_stem_max: int = 7
    repeat_k: int = 20
    transcript_nt_max: int = 6000
    junction_risk_max: float = 2.0
    min_cai: float = 0.75


def _cap_proximal_stem(seq: str, min_stem: int = 4, max_loop: int = 30) -> tuple[int, int, int]:
    """Longest perfect self-complementary stem within ``seq``.

    Crude structure proxy: finds the longest pair of reverse-complementary
    substrings separated by a loop of at least three nucleotides. Real folding
    is an ensemble calculation -- this only catches the blunt cases.
    Returns (stem_length, left_start, right_start).
    """
    best = (0, -1, -1)
    n = len(seq)
    for i in range(n):
        for j in range(i + min_stem + 3, n):
            length = 0
            while (
                i + length < j - length - 3
                and j + 1 - length > 0
                and seq[i + length] == revcomp(seq[j - length])
            ):
                length += 1
            if length > best[0] and length >= min_stem:
                best = (length, i, j - length + 1)
    return best


def check_orf(construct: Construct) -> list[Check]:
    orf = construct.orf
    checks = []
    checks.append(
        Check("orf.start", Severity.PASS if orf.startswith("ATG") else Severity.FAIL,
              "ORF begins with ATG" if orf.startswith("ATG")
              else f"ORF begins with {orf[:3]!r}, not ATG")
    )
    in_frame = len(orf) % 3 == 0
    checks.append(
        Check("orf.frame", Severity.PASS if in_frame else Severity.FAIL,
              f"ORF length {len(orf)} nt is a multiple of 3" if in_frame
              else f"ORF length {len(orf)} nt is not a multiple of 3")
    )
    protein = translate(orf, stop_at_stop=False)
    internal = [i for i, aa in enumerate(protein[:-1]) if aa == "*"]
    # The tandem stop puts a legitimate second stop at the tail; only stops
    # before the first one are defects.
    first_stop = protein.find("*")
    premature = [i for i in internal if first_stop != -1 and i < first_stop]
    checks.append(
        Check("orf.internal_stop",
              Severity.PASS if not premature else Severity.FAIL,
              "no premature stop codon" if not premature
              else f"premature stop codon(s) at codon {premature}")
    )
    checks.append(
        Check("orf.terminates", Severity.PASS if "*" in protein else Severity.FAIL,
              "ORF terminates in a stop codon" if "*" in protein
              else "ORF has no stop codon")
    )
    return checks


def check_utr5(construct: Construct) -> list[Check]:
    checks = []
    try:
        utr5 = construct.features_of_kind("utr5")[0]
    except IndexError:
        return [Check("utr5.present", Severity.FAIL, "no 5' UTR feature")]
    leader = construct.template[construct.transcript_start : utr5.end]
    uaugs = find_all(leader, "ATG")
    checks.append(
        Check("utr5.upstream_aug",
              Severity.PASS if not uaugs else Severity.FAIL,
              "no upstream AUG in the 5' UTR" if not uaugs
              else f"{len(uaugs)} upstream AUG(s) at {uaugs} will capture scanning "
                   "ribosomes and suppress the main ORF",
              {"positions": uaugs})
    )
    orf = construct.feature("ORF")
    context = construct.template[max(0, orf.start - 6) : orf.start + 4]
    strong = len(context) == 10 and context[3] in "AG" and context[-1] == "G"
    checks.append(
        Check("utr5.kozak", Severity.PASS if strong else Severity.WARN,
              f"Kozak context {context!r} is strong (purine at -3, G at +4)" if strong
              else f"Kozak context {context!r} is not the strong consensus "
                   "(gccRccATGG); initiation may leak to downstream AUGs",
              {"context": context})
    )
    return checks


def check_composition(construct: Construct, t: QCThresholds) -> list[Check]:
    transcript = construct.transcript
    checks = []

    gc = gc_fraction(transcript)
    in_band = t.gc_min <= gc <= t.gc_max
    checks.append(
        Check("composition.gc", Severity.PASS if in_band else Severity.WARN,
              f"transcript GC {gc:.1%} " + ("is within "
              f"{t.gc_min:.0%}-{t.gc_max:.0%}" if in_band else
              f"is outside {t.gc_min:.0%}-{t.gc_max:.0%}"),
              {"gc": round(gc, 4)})
    )

    # Windowed GC, evaluated on the coding region only: the poly(A) tail would
    # otherwise dominate the low-GC tail of the distribution by design.
    orf = construct.feature("ORF")
    coding = orf.slice(construct.template)
    excursions = [
        (start, round(value, 3))
        for start, value in gc_windows(coding, t.gc_window, step=10)
        if value < t.window_gc_min or value > t.window_gc_max
    ]
    checks.append(
        Check("composition.gc_windows",
              Severity.PASS if not excursions else Severity.WARN,
              f"no {t.gc_window}-nt window in the ORF outside "
              f"{t.window_gc_min:.0%}-{t.window_gc_max:.0%}" if not excursions
              else f"{len(excursions)} window(s) outside "
                   f"{t.window_gc_min:.0%}-{t.window_gc_max:.0%}; GC-rich stretches "
                   "are where T7 stalls and produces truncated species",
              {"excursions": excursions[:20]})
    )

    u = transcript.count("T") / len(transcript)
    checks.append(
        Check("composition.uridine", Severity.PASS if u <= t.uridine_max else Severity.WARN,
              f"uridine content {u:.1%} " + ("is at or below "
              f"{t.uridine_max:.0%}" if u <= t.uridine_max
              else f"exceeds {t.uridine_max:.0%}; raises innate sensing and T7 "
                   "slippage risk even with m1-pseudouridine substitution"),
              {"uridine_fraction": round(u, 4)})
    )

    cai = construct.optimization.cai
    checks.append(
        Check("composition.cai", Severity.PASS if cai >= t.min_cai else Severity.WARN,
              f"codon adaptation index {cai:.3f} " +
              ("meets" if cai >= t.min_cai else "is below") + f" the {t.min_cai} floor",
              {"cai": round(cai, 4)})
    )
    return checks


def check_manufacturability(construct: Construct, t: QCThresholds) -> list[Check]:
    checks = []
    template = construct.template

    # Homopolymers, excluding the poly(A) tail, which is one by construction.
    polya = construct.features_of_kind("polya")
    masked = list(template)
    for feature in polya:
        for i in range(feature.start, feature.end):
            masked[i] = "N"
    runs = [r for r in homopolymer_runs("".join(masked), t.max_homopolymer)
            if r[0] != "N"]
    checks.append(
        Check("mfg.homopolymer",
              Severity.PASS if not runs else Severity.WARN,
              f"no single-base run of {t.max_homopolymer}+ nt outside the poly(A) tail"
              if not runs else
              f"{len(runs)} run(s) of {t.max_homopolymer}+ nt: T7 slips on these and "
              "they are hard to sequence through",
              {"runs": runs[:20]})
    )

    # Direct repeats: plasmid recombination and assembly ambiguity.
    seen: dict[str, int] = {}
    repeats: list[tuple[str, int, int]] = []
    for i in range(len(template) - t.repeat_k + 1):
        kmer = template[i : i + t.repeat_k]
        if "A" * t.repeat_k == kmer:
            continue
        if kmer in seen:
            repeats.append((kmer, seen[kmer], i))
        else:
            seen[kmer] = i
    checks.append(
        Check("mfg.direct_repeats",
              Severity.PASS if not repeats else Severity.WARN,
              f"no exact direct repeat of {t.repeat_k}+ nt" if not repeats
              else f"{len(repeats)} exact repeat(s) of {t.repeat_k} nt; these "
                   "recombine during plasmid propagation and confound assembly. "
                   "Common cause: identical linkers between every bead",
              {"repeats": [(a, b) for _, a, b in repeats[:20]]})
    )

    # The linearisation site must cut exactly once in the whole template.
    for feature in construct.features_of_kind("linearization"):
        site = feature.slice(template)
        hits = find_all(template, site, both_strands=True)
        unique = len(hits) == 1
        checks.append(
            Check(f"mfg.linearization.{feature.name}",
                  Severity.PASS if unique else Severity.FAIL,
                  f"{feature.name} site {site} is unique in the template" if unique
                  else f"{feature.name} site {site} occurs {len(hits)} times "
                       f"at {hits}; run-off linearisation would cut internally",
                  {"positions": hits})
        )

    # Cryptic poly(A) signal, attributed to its feature. AAUAAA only directs
    # cleavage and polyadenylation in the nucleus, so in a cytoplasmically
    # delivered transcript it is inert -- and the validated BNT162b2 UTRs each
    # contain one. It is a warning only when the optimiser put it in the ORF.
    transcript = construct.transcript
    orf = construct.feature("ORF")
    paus = find_all(transcript, "AATAAA")
    in_orf = [p for p in paus
              if orf.start <= construct.transcript_start + p < orf.end]
    elsewhere = [(p, _owning_feature(construct, p)) for p in paus if p not in in_orf]
    if in_orf:
        message = (f"AAUAAA at {in_orf} inside the ORF: premature "
                   "polyadenylation risk if the template is ever transcribed in "
                   "a nucleus")
    elif elsewhere:
        message = ("AAUAAA present only in fixed parts ("
                   + ", ".join(f"{owner}@{p}" for p, owner in elsewhere[:6])
                   + "); inert for a cytoplasmic transcript")
    else:
        message = "no AAUAAA in the transcript"
    checks.append(
        Check("mfg.cryptic_polya",
              Severity.WARN if in_orf else Severity.PASS, message,
              {"in_orf": in_orf, "elsewhere": elsewhere})
    )
    return checks


def check_cap_proximal(construct: Construct, t: QCThresholds) -> list[Check]:
    head = construct.transcript[: t.cap_proximal_nt]
    gc = gc_fraction(head)
    stem, left, right = _cap_proximal_stem(head)
    checks = [
        Check("cap.gc", Severity.PASS if gc <= t.cap_proximal_gc_max else Severity.WARN,
              f"cap-proximal GC {gc:.1%} over the first {t.cap_proximal_nt} nt "
              + ("is acceptable" if gc <= t.cap_proximal_gc_max else
                 "is high; structure here impedes 43S loading. Proxy measure, not a "
                 "folding energy -- confirm with a folding tool"),
              {"gc": round(gc, 4)}),
        Check("cap.stem",
              Severity.PASS if stem <= t.cap_proximal_stem_max else Severity.WARN,
              f"longest perfect cap-proximal stem is {stem} bp"
              + ("" if stem <= t.cap_proximal_stem_max else
                 f" (at {left}/{right}); a stable 5' hairpin blocks scanning. "
                 "Proxy measure -- confirm with a folding tool"),
              {"stem_bp": stem}),
        Check("cap.start_base",
              Severity.PASS if construct.transcript[:1] == "G" else Severity.WARN,
              f"transcript begins with {construct.transcript[:1]}"
              + ("" if construct.transcript[:1] == "G" else
                 "; T7 initiates far more efficiently from G at +1"),
              {}),
    ]
    return checks


def _owning_feature(construct: Construct, transcript_pos: int) -> str:
    """Name the smallest annotated feature containing a transcript position."""
    absolute = construct.transcript_start + transcript_pos
    owners = [
        f for f in construct.features
        if f.start <= absolute < f.end and f.kind not in ("orf",)
    ]
    if not owners:
        return "ORF" if construct.feature("ORF").start <= absolute < \
            construct.feature("ORF").end else "?"
    return min(owners, key=len).name


def check_motifs(construct: Construct, optimizer: CodonOptimizer | None = None) -> list[Check]:
    """Screen the transcript for forbidden motifs, attributed to their feature.

    Attribution matters because a validated UTR is a fixed part you chose, not
    something the optimiser can rewrite. The BNT162b2 UTRs, for instance, carry
    XhoI and NheI cloning scars and an AAUAAA inside mtRNR1 -- all in the
    clinical sequence. A hit there is a fact to know before you pick a
    restriction enzyme; a hit in the ORF is a miss by the optimiser. Only the
    latter is a warning.
    """
    optimizer = optimizer or CodonOptimizer()
    # Mask the poly(A) tail: a long A-tract matches pyrimidine-tract motifs on
    # the reverse strand by construction, which is noise, not a finding.
    masked = list(construct.transcript)
    offset = construct.transcript_start
    for feature in construct.features_of_kind("polya"):
        for i in range(feature.start - offset, feature.end - offset):
            if 0 <= i < len(masked):
                masked[i] = "N"
    transcript = "".join(masked)

    orf = construct.feature("ORF")
    in_orf: list[tuple[str, int, str]] = []
    in_fixed: list[tuple[str, int, str]] = []
    for motif, position in optimizer.find_forbidden(transcript):
        owner = _owning_feature(construct, position)
        absolute = offset + position
        record = (motif, position, owner)
        (in_orf if orf.start <= absolute < orf.end else in_fixed).append(record)

    checks = [
        Check("motif.forbidden",
              Severity.PASS if not in_orf else Severity.WARN,
              "no forbidden motif in the ORF" if not in_orf
              else f"{len(in_orf)} forbidden motif(s) survive in the ORF: "
                   + ", ".join(f"{m}@{p}" for m, p, _ in in_orf[:8]),
              {"hits": in_orf[:40]})
    ]
    if in_fixed:
        checks.append(
            Check("motif.fixed_parts", Severity.PASS,
                  f"{len(in_fixed)} motif(s) present in fixed parts you chose "
                  "(not optimiser-editable): "
                  + ", ".join(f"{m}@{o}" for m, _, o in in_fixed[:8])
                  + ". Confirm none clashes with your cloning or linearisation "
                    "enzymes",
                  {"hits": in_fixed[:40]})
        )
    return checks


def check_placeholders(construct: Construct) -> list[Check]:
    used = construct.placeholders_used
    return [
        Check("design.placeholders",
              Severity.PASS if not used else Severity.FAIL,
              "no placeholder parts in the construct" if not used
              else "construct contains placeholder parts that are not real "
                   "sequences: " + ", ".join(used),
              {"parts": used})
    ]


def check_junctions(
    construct: Construct,
    t: QCThresholds,
    scorer: JunctionScorer | None = None,
    lengths: Sequence[int] = CLASS_I_LENGTHS,
) -> list[Check]:
    hits = scan_junctions(construct.cassette, scorer, lengths)
    risk = sum(h.score for h in hits)
    return [
        Check("immuno.junctional_epitopes",
              Severity.PASS if risk <= t.junction_risk_max else Severity.WARN,
              f"junctional risk {risk:.2f} across {len(hits)} flagged peptide(s) "
              + ("is within" if risk <= t.junction_risk_max else "exceeds")
              + f" the {t.junction_risk_max} budget. Score is from the supplied "
                "scorer -- with the built-in anchor proxy, treat it as a relative "
                "ranking only, not a binding prediction",
              {"risk": round(risk, 3),
               "top": [(h.junction_index, h.peptide, round(h.score, 3))
                       for h in hits[:15]]})
    ]


def check_length(construct: Construct, t: QCThresholds) -> list[Check]:
    n = construct.transcript_length
    return [
        Check("mfg.transcript_length",
              Severity.PASS if n <= t.transcript_nt_max else Severity.WARN,
              f"transcript is {n} nt"
              + ("" if n <= t.transcript_nt_max else
                 f", above the {t.transcript_nt_max} nt working limit; yield and "
                 "integrity of long IVT products fall off and LNP encapsulation "
                 "efficiency drops"),
              {"nt": n})
    ]


def run_qc(
    construct: Construct,
    thresholds: QCThresholds | None = None,
    scorer: JunctionScorer | None = None,
    optimizer: CodonOptimizer | None = None,
    allow_placeholders: bool = False,
    extra: Sequence[Callable[[Construct], list[Check]]] = (),
) -> QCReport:
    """Run the full battery. ``extra`` takes your own checks."""
    t = thresholds or QCThresholds()
    checks: list[Check] = []
    checks += check_orf(construct)
    checks += check_utr5(construct)
    checks += check_composition(construct, t)
    checks += check_cap_proximal(construct, t)
    checks += check_manufacturability(construct, t)
    checks += check_motifs(construct, optimizer)
    checks += check_length(construct, t)
    checks += check_junctions(construct, t, scorer)
    for check in check_placeholders(construct):
        if allow_placeholders and check.severity is Severity.FAIL:
            check = Check(check.name, Severity.WARN,
                          check.message + " (downgraded: allow_placeholders=True)",
                          check.detail)
        checks.append(check)
    for fn in extra:
        checks += fn(construct)
    return QCReport(construct.name, checks)


__all__ = [
    "Check", "QCReport", "QCThresholds", "Severity", "run_qc",
    "check_orf", "check_utr5", "check_composition", "check_cap_proximal",
    "check_manufacturability", "check_motifs", "check_junctions",
    "check_placeholders", "check_length",
]
