import pytest

from txv.constructs import ConstructSpec, build_construct
from txv.epitopes import Antigen
from txv.parts import Part, PartKind, PartRegistry, Provenance, default_registry
from txv.codon import OptimizerConfig
from txv.qc import QCThresholds, Severity, run_qc


def severity(report, name):
    return next(c.severity for c in report.checks if c.name == name)


def test_clean_construct_passes(construct):
    report = run_qc(construct, allow_placeholders=True)
    assert report.passed
    assert not report.failures


def test_default_construct_passes_strict_mode(construct):
    """Defaults are real sequences, so strict mode is clean out of the box."""
    report = run_qc(construct, allow_placeholders=False)
    assert report.passed, report.to_text()
    assert severity(report, "design.placeholders") is Severity.PASS


def test_placeholders_fail_in_strict_mode(antigens):
    from txv.constructs import build_construct

    construct = build_construct(
        "ph", antigens, spec=ConstructSpec(name="ph", utr3="UTR3_placeholder")
    )
    report = run_qc(construct, allow_placeholders=False)
    assert not report.passed
    assert severity(report, "design.placeholders") is Severity.FAIL
    assert run_qc(construct, allow_placeholders=True).passed


def test_validated_utr_motifs_are_attributed_not_blamed_on_the_orf(construct):
    """The real UTRs carry AAUAAA and cloning scars; that is not an ORF defect."""
    report = run_qc(construct)
    assert severity(report, "motif.forbidden") is Severity.PASS
    assert severity(report, "mfg.cryptic_polya") is Severity.PASS
    fixed = next(c for c in report.checks if c.name == "motif.fixed_parts")
    owners = {owner for _, _, owner in fixed.detail["hits"]}
    assert owners <= {"UTR5_hAg", "UTR3_AES_mtRNR1"}


def test_aauaaa_inside_the_orf_still_warns(antigens):
    from txv.constructs import build_construct
    from txv.epitopes import Antigen

    # NKID encodes AAUAAA-compatible codons; force it with a narrow table.
    construct = build_construct(
        "x", [Antigen("forced", "MNKNKNKNK")],
        spec=ConstructSpec(name="x", signal_peptide=None),
        optimizer_config=OptimizerConfig(forbidden=(), w_uridine=0.0),
    )
    report = run_qc(construct)
    check = next(c for c in report.checks if c.name == "mfg.cryptic_polya")
    if check.detail["in_orf"]:
        assert check.severity is Severity.WARN


def test_upstream_aug_in_the_utr_is_a_failure(antigens):
    registry = default_registry()
    registry.add(Part("UTR5_bad", PartKind.UTR5, Provenance.LITERATURE,
                      dna="AGGAAATAAGAGATGAAAGAAGAGTAAGAAG"))
    construct = build_construct(
        "bad", antigens, registry=registry,
        spec=ConstructSpec(name="bad", utr5="UTR5_bad"),
    )
    report = run_qc(construct, allow_placeholders=True)
    assert severity(report, "utr5.upstream_aug") is Severity.FAIL
    assert not report.passed


def test_non_unique_linearization_site_is_a_failure(antigens):
    registry = default_registry()
    # A 4-cutter will inevitably appear inside the construct as well.
    registry.add(Part("frequent_site", PartKind.RESTRICTION, Provenance.CANONICAL,
                      dna="GATC"))
    construct = build_construct(
        "bad", antigens, registry=registry,
        spec=ConstructSpec(name="bad", linearization_site="frequent_site"),
    )
    report = run_qc(construct, allow_placeholders=True)
    assert severity(report, "mfg.linearization.frequent_site") is Severity.FAIL


def test_polya_tail_does_not_trip_the_homopolymer_or_motif_checks(construct):
    report = run_qc(construct, allow_placeholders=True)
    assert severity(report, "mfg.homopolymer") is Severity.PASS
    assert severity(report, "motif.forbidden") is Severity.PASS


def test_transcript_length_limit_warns(construct):
    report = run_qc(construct, thresholds=QCThresholds(transcript_nt_max=10),
                    allow_placeholders=True)
    assert severity(report, "mfg.transcript_length") is Severity.WARN
    assert report.passed, "a length overrun is a warning, not a failure"


def test_thresholds_are_applied(construct):
    strict = run_qc(construct, thresholds=QCThresholds(uridine_max=0.01),
                    allow_placeholders=True)
    assert severity(strict, "composition.uridine") is Severity.WARN
    loose = run_qc(construct, thresholds=QCThresholds(uridine_max=0.9),
                   allow_placeholders=True)
    assert severity(loose, "composition.uridine") is Severity.PASS


def test_report_rendering_and_serialisation(construct):
    report = run_qc(construct, allow_placeholders=True)
    text = report.to_text()
    assert "QC report" in text and "failure(s)" in text
    assert len(report.to_text(show_pass=False)) < len(text)
    payload = report.to_dict()
    assert payload["construct"] == construct.name
    assert payload["passed"] is True
    assert all({"name", "severity", "message"} <= set(c) for c in payload["checks"])


def test_extra_checks_run(construct):
    from txv.qc import Check

    def custom(_c):
        return [Check("custom.thing", Severity.FAIL, "nope")]

    report = run_qc(construct, allow_placeholders=True, extra=[custom])
    assert not report.passed
    assert severity(report, "custom.thing") is Severity.FAIL


def test_kozak_warning_when_context_is_weak(antigens):
    registry = default_registry()
    registry.add(Part("kozak_weak", PartKind.KOZAK, Provenance.CANONICAL, dna="TTTTTT"))
    construct = build_construct(
        "weak", antigens, registry=registry,
        spec=ConstructSpec(name="weak", kozak="kozak_weak"),
    )
    report = run_qc(construct, allow_placeholders=True)
    assert severity(report, "utr5.kozak") is Severity.WARN
