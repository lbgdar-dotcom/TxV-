"""Build and verify the finished plasmid: backbone + insert, fully annotated.

This is the in-silico counterpart of the bench step. It takes the parent vector
map, swaps the old ORF for a new one at the homology rails, remaps every
downstream feature, carries the insert's per-module annotation across, and then
**re-derives** the things that have to be true of the product rather than
assuming the swap went as intended.

The verification is the point. A plasmid map that was produced by string
surgery and never checked is exactly as trustworthy as the surgery, which is to
say not very -- so :func:`verify_plasmid` recomputes the reading frame, the
transcription unit, the linearisation sites and the regulatory elements from
the assembled sequence itself.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .constructs import Construct
from .genbank_io import GenBankFeature, GenBankRecord
from .seqops import find_all, translate

#: Feature kinds from a :class:`~txv.constructs.Construct` worth carrying into
#: the plasmid map, and the GenBank key to record them under.
_INSERT_KEYS = {
    "orf": "CDS", "signal_peptide": "sig_peptide", "trafficking": "misc_feature",
    "linker": "misc_feature", "neoepitope": "misc_feature", "tag": "misc_feature",
    "degron": "misc_feature", "skip_peptide": "misc_feature",
    "start": "misc_feature", "stop": "misc_feature",
    "tumor_associated": "misc_feature", "helper": "misc_feature",
    "full_length": "misc_feature",
}


@dataclass
class PlasmidBuild:
    record: GenBankRecord
    insert_start: int
    insert_end: int
    removed_span: tuple[int, int]
    removed_sequence: str
    dropped_features: list[str] = field(default_factory=list)

    @property
    def insert(self) -> str:
        return self.record.sequence[self.insert_start : self.insert_end]

    def __len__(self) -> int:
        return len(self.record.sequence)


def build_plasmid(
    backbone: GenBankRecord,
    construct: Construct,
    name: str,
    left_rail: str = "GAAGAAATATAAGAGCCACC",
    right_rail: str = "GCTGCCTTCTGCGGGGCTTG",
) -> PlasmidBuild:
    """Swap ``construct``'s ORF into ``backbone`` between the homology rails.

    Features wholly upstream of the insertion keep their coordinates; features
    wholly downstream are shifted. Features that *overlap* the replaced region
    -- the old ORF and any primer binding inside it -- are dropped and named in
    :attr:`PlasmidBuild.dropped_features`, because silently keeping a primer
    annotation that no longer anneals is how a map starts lying to you.
    """
    sequence = backbone.sequence
    left_hits, right_hits = find_all(sequence, left_rail), find_all(sequence, right_rail)
    for label, hits, rail in (("left", left_hits, left_rail),
                              ("right", right_hits, right_rail)):
        if len(hits) != 1:
            raise ValueError(
                f"{label} rail {rail} occurs {len(hits)} times in {backbone.name} "
                "(expected exactly 1); the assembly would be ambiguous"
            )

    start = left_hits[0] + len(left_rail)
    end = right_hits[0]
    if end <= start:
        raise ValueError(
            "the right rail precedes the left rail in this linear representation; "
            "rotate the backbone so the insertion site does not span the origin"
        )

    insert = construct.orf
    delta = len(insert) - (end - start)
    new_sequence = sequence[:start] + insert + sequence[end:]

    features: list[GenBankFeature] = []
    dropped: list[str] = []
    for feature in backbone.features:
        if feature.end <= start:
            features.append(feature)
        elif feature.start >= end:
            features.append(GenBankFeature(
                feature.key, feature.start + delta, feature.end + delta,
                feature.strand, dict(feature.qualifiers),
            ))
        elif feature.start <= start and feature.end >= end and len(feature) == len(sequence):
            features.append(GenBankFeature(
                feature.key, feature.start, feature.end + delta, feature.strand,
                dict(feature.qualifiers),
            ))
        else:
            dropped.append(feature.label)

    offset = start - construct.feature("ORF").start
    for feature in construct.features:
        key = _INSERT_KEYS.get(feature.kind)
        if key is None:
            continue
        qualifiers = {"label": feature.name}
        if feature.kind == "orf":
            qualifiers["label"] = f"{name}_ORF"
            qualifiers["translation"] = construct.orf_protein
            qualifiers["codon_start"] = "1"
        if feature.note:
            qualifiers["note"] = feature.note[:400]
        features.append(GenBankFeature(
            key, feature.start + offset, feature.end + offset, 1, qualifiers,
        ))

    record = GenBankRecord(
        name=name,
        sequence=new_sequence,
        features=features,
        is_circular=backbone.is_circular,
        definition=(
            f"{name}: {backbone.name} with eGFP replaced by the "
            f"{construct.name} antigen cassette."
        ),
    )
    return PlasmidBuild(
        record=record, insert_start=start, insert_end=start + len(insert),
        removed_span=(start, end), removed_sequence=sequence[start:end],
        dropped_features=dropped,
    )


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

@dataclass
class PlasmidCheck:
    name: str
    ok: bool
    detail: str


def verify_plasmid(
    build: PlasmidBuild,
    construct: Construct,
    backbone: GenBankRecord,
    linearisation_site: str = "GGTCTC",
    expected_linearisation_sites: int = 2,
) -> list[PlasmidCheck]:
    """Re-derive what must be true of the assembled plasmid from its sequence."""
    sequence = build.record.sequence
    checks: list[PlasmidCheck] = []

    def check(name, ok, detail):
        checks.append(PlasmidCheck(name, bool(ok), detail))

    # Reading frame, recomputed from the sequence rather than trusted.
    orf = sequence[build.insert_start : build.insert_end]
    protein = translate(orf, stop_at_stop=True).rstrip("*")
    check("orf.identity", protein == construct.orf_protein,
          f"{len(protein)} aa, matches the designed protein: "
          f"{protein == construct.orf_protein}")
    check("orf.start", orf.startswith("ATG"), f"ORF begins {orf[:3]}")
    check("orf.frame", len(orf) % 3 == 0, f"{len(orf)} nt, multiple of 3")
    check("orf.terminates", "*" in translate(orf, stop_at_stop=False),
          "ORF carries its own stop codon")

    # The Kozak must survive the junction, and there must be no new AUG before it.
    kozak = sequence[build.insert_start - 6 : build.insert_start + 3]
    check("junction.kozak", kozak == "GCCACCATG", f"5' junction reads {kozak}")

    # Transcription unit. The recognised T7 element is the 17 nt -17..-1
    # (TAATACGACTCACTATA); the next base is +1, the first base of the
    # transcript. In this vector that base is A and the transcript begins AGG,
    # which is what CleanCap AG requires -- hence "pVax1_AG".
    t7_core = "TAATACGACTCACTATA"
    t7 = find_all(sequence, t7_core)
    check("t7.present", len(t7) == 1, f"T7 -17..-1 at {[p + 1 for p in t7]}")
    if t7:
        plus_one = t7[0] + len(t7_core)
        start5 = sequence[plus_one : plus_one + 3]
        check("t7.start_AGG", start5.startswith("AG"),
              f"transcript +1 at {plus_one + 1}, begins {start5} "
              "(CleanCap AG needs AG)")

    utr5 = sequence[build.insert_start - 44 : build.insert_start]
    check("utr5.no_uaug", "ATG" not in utr5,
          f"no upstream AUG in the 44-nt 5' UTR ({len(find_all(utr5, 'ATG'))} found)")

    # Linearisation sites must be preserved, and the insert must add none.
    sites = find_all(sequence, linearisation_site, both_strands=True)
    check("linearisation.count", len(sites) == expected_linearisation_sites,
          f"{linearisation_site} at {[p + 1 for p in sites]} "
          f"(expected {expected_linearisation_sites})")
    check("linearisation.not_in_insert",
          not any(build.insert_start <= p < build.insert_end for p in sites),
          "no linearisation site inside the insert")

    # Poly(A) tract intact.
    runs = [len(r) for r in sequence.split("G") if set(r) == {"A"}]
    longest_a = max((len(m) for m in _a_runs(sequence)), default=0)
    check("polya.intact", longest_a >= 120,
          f"longest A-tract {longest_a} nt")

    # Nuclear-transcription hazards, which matter because this vector also
    # drives the same ORF from its CMV promoter.
    paus = [p for p in find_all(sequence, "AATAAA")
            if build.insert_start <= p < build.insert_end]
    check("cmv.no_polya_signal_in_orf", not paus,
          f"AATAAA inside the ORF: {[p + 1 for p in paus] or 'none'}")

    # Backbone elements that must survive untouched.
    for label in ("CMV promoter", "NeoR/KanR", "ori"):
        try:
            original = backbone.find(label).slice(backbone.sequence)
        except KeyError:
            continue
        check(f"backbone.{label.split()[0].lower()}",
              original in sequence, f"{label} present and unchanged")

    check("size.arithmetic",
          len(sequence) == len(backbone.sequence)
          - len(build.removed_sequence) + len(build.insert),
          f"{len(backbone.sequence)} - {len(build.removed_sequence)} + "
          f"{len(build.insert)} = {len(sequence)} bp")
    return checks


def _a_runs(sequence: str):
    import re

    return [m.group(0) for m in re.finditer(r"A+", sequence)]





# ---------------------------------------------------------------------------
# Annotation repair
# ---------------------------------------------------------------------------

def annotate_transcription_unit(
    record: GenBankRecord,
    orf_start: int,
    orf_end: int,
    t7_core: str = "TAATACGACTCACTATA",
    polya_label: str = "poly(A) 120",
    utr3_label: str = "3' UTR",
) -> list[str]:
    """Re-derive the transcription unit's annotation from the sequence.

    Inherited vector maps tend to under-label the untranslated regions, because
    a feature gets drawn once around whatever was being worked on and never
    revisited. That is harmless until someone reads the map to answer a
    question it was never annotated to answer -- "does this transcript contain
    a poly(A) signal?" being exactly such a question.

    This corrects three things in place and returns a list of what it changed:

    * the **3' UTR** feature is extended to its true extent -- from the stop
      codon to the start of the encoded poly(A) tract -- rather than wherever
      the original annotation happened to stop;
    * a **transcript** feature is added spanning T7 +1 to the end of the
      poly(A), so the map shows what the IVT product actually is;
    * any **AAUAAA** inside the transcript is labelled, with its context, since
      whether one is native and harmless or cryptic and damaging depends
      entirely on where it sits.
    """
    changes: list[str] = []
    sequence = record.sequence

    try:
        polya = record.find(polya_label)
    except KeyError:
        return changes

    # -- 3' UTR: stop codon to poly(A) ------------------------------------
    try:
        utr3 = record.find(utr3_label)
    except KeyError:
        utr3 = None
    if utr3 is not None and (utr3.start != orf_end or utr3.end != polya.start):
        old = (utr3.start + 1, utr3.end)
        utr3.start, utr3.end = orf_end, polya.start
        utr3.qualifiers["note"] = (
            "Corrected to the full transcribed 3' UTR: stop codon to poly(A). "
            f"Previously annotated {old[0]}..{old[1]}, which stopped short and "
            "left part of the UTR -- including its poly(A) signal -- unlabelled."
        )
        changes.append(
            f"3' UTR {old[0]}..{old[1]} -> {utr3.start + 1}..{utr3.end} "
            f"({utr3.end - utr3.start} nt)"
        )

    # -- the transcript itself --------------------------------------------
    hits = find_all(sequence, t7_core)
    if len(hits) == 1 and not record.find_all("transcript (T7 run-off)"):
        plus_one = hits[0] + len(t7_core)
        record.features.append(GenBankFeature(
            "misc_RNA", plus_one, polya.end, 1,
            {
                "label": "transcript (T7 run-off)",
                "note": (
                    f"IVT product: {polya.end - plus_one} nt from T7 +1 "
                    f"(begins {sequence[plus_one:plus_one + 3]}) through the "
                    "encoded poly(A). Cap chemistry requiring an AG start is "
                    "compatible with this +1."
                ),
            },
        ))
        changes.append(
            f"added transcript (T7 run-off) {plus_one + 1}..{polya.end} "
            f"({polya.end - plus_one} nt)"
        )

    # -- poly(A) signals, labelled with their context ----------------------
    for position in find_all(sequence, "AATAAA"):
        if not (hits and hits[0] + len(t7_core) <= position < polya.end):
            continue
        if any(f.start == position and f.key == "polyA_signal"
               for f in record.features):
            continue
        in_orf = orf_start <= position < orf_end
        where = "INSIDE THE ORF" if in_orf else (
            "in the 3' UTR" if position >= orf_end else "in the 5' UTR/leader"
        )
        record.features.append(GenBankFeature(
            "polyA_signal", position, position + 6, 1,
            {
                "label": f"AATAAA ({where})",
                "note": (
                    "Cryptic poly(A) signal in the coding sequence -- would "
                    "truncate the mRNA when this ORF is transcribed in the "
                    "nucleus from the CMV promoter. This should not be here."
                    if in_orf else
                    "Native poly(A) signal of the untranslated region. Acts "
                    "only on nuclear transcripts (i.e. the CMV route); inert "
                    "for a cytoplasmically delivered IVT mRNA."
                ),
            },
        ))
        changes.append(f"labelled AAUAAA at {position + 1} ({where})")

    return changes


__all__ = [
    "PlasmidBuild", "PlasmidCheck",
    "build_plasmid", "verify_plasmid", "annotate_transcription_unit",
]
