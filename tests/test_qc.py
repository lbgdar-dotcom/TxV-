import pytest

from txv.constructs import ConstructSpec, build_construct
from txv.epitopes import Antigen
from txv.parts import Part, PartKind, PartRegistry, Provenance, default_registry
from txv.qc import QCThresholds, Severity, run_qc


def severity(report, name):
    return next(c.severity for c in report.checks if c.name == name)


def test_clean_construct_passes(construct):
    report = run_qc(construct, allow_placeholders=True)
    assert report.passed
    assert not report.failures


def test_placeholders_fail_in_strict_mode(construct):
    report = run_qc(construct, allow_placeholders=False)
    assert not report.passed
    assert severity(report, "design.placeholders") is Severity.FAIL


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
