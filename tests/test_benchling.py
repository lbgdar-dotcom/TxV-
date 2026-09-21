import json

import pytest

from txv.benchling import (
    BenchlingConfig,
    DryRunRegistrar,
    construct_payload,
    guard_qc,
    make_registrar,
)
from txv.benchling.client import ANNOTATION_TYPE
from txv.qc import QCReport, Check, Severity, run_qc


@pytest.fixture
def config():
    return BenchlingConfig(folder_id="lib_x", schema_id="ts_x")


def test_payload_carries_sequence_annotations_and_fields(construct, config):
    qc = run_qc(construct, allow_placeholders=True)
    payload = construct_payload(construct, config, qc)
    assert payload["name"] == construct.name
    assert payload["bases"] == construct.template
    assert payload["isCircular"] is False
    assert payload["folderId"] == "lib_x"
    assert payload["schemaId"] == "ts_x"
    assert len(payload["annotations"]) == len(construct.features)
    assert payload["fields"]["QC status"]["value"] == "PASS"
    assert payload["fields"]["Antigen count"]["value"] == len(
        construct.cassette.antigens
    )


def test_annotation_coordinates_match_the_features(construct, config):
    payload = construct_payload(construct, config)
    for annotation, feature in zip(payload["annotations"], construct.features):
        assert annotation["start"] == feature.start
        assert annotation["end"] == feature.end
        assert annotation["strand"] == 1
        assert annotation["name"] == feature.name
        assert annotation["type"] == ANNOTATION_TYPE.get(feature.kind,
                                                         "misc_feature")


def test_annotations_stay_inside_the_sequence(construct, config):
    payload = construct_payload(construct, config)
    for annotation in payload["annotations"]:
        assert 0 <= annotation["start"] < annotation["end"] <= len(payload["bases"])


def test_registry_fields_only_appear_when_a_registry_is_configured(construct):
    without = construct_payload(construct, BenchlingConfig(folder_id="lib_x"))
    assert "registryId" not in without
    with_registry = construct_payload(
        construct, BenchlingConfig(folder_id="lib_x", registry_id="src_x")
    )
    assert with_registry["registryId"] == "src_x"
    assert with_registry["namingStrategy"] == "NEW_IDS"


def test_extra_and_default_fields_are_merged(construct):
    config = BenchlingConfig(folder_id="lib_x",
                             default_fields={"Programme": "TXV"})
    payload = construct_payload(construct, config, extra_fields={"Batch": "B1"})
    assert payload["fields"]["Programme"]["value"] == "TXV"
    assert payload["fields"]["Batch"]["value"] == "B1"


def test_dry_run_writes_payload_genbank_and_qc(construct, config, tmp_path):
    qc = run_qc(construct, allow_placeholders=True)
    registrar = DryRunRegistrar(config, tmp_path)
    record = registrar.register(construct, qc)
    assert record["mode"] == "dry-run"
    payload_path = tmp_path / f"{construct.name}.payload.json"
    assert json.loads(payload_path.read_text())["bases"] == construct.template
    assert (tmp_path / f"{construct.name}.gb").read_text().startswith("LOCUS")
    assert json.loads((tmp_path / f"{construct.name}.qc.json").read_text())["passed"]
    assert registrar.sent == [record]


def test_dry_run_register_many(construct, config, tmp_path):
    registrar = DryRunRegistrar(config, tmp_path)
    records = registrar.register_many([(construct, None)])
    assert len(records) == 1 and len(registrar.sent) == 1


def test_guard_refuses_a_failing_qc_report(construct):
    failing = QCReport(construct.name,
                       [Check("x.y", Severity.FAIL, "broken")])
    with pytest.raises(ValueError, match="refusing to register"):
        guard_qc(construct, failing)
    # A passing report, and no report at all, are both allowed.
    guard_qc(construct, run_qc(construct, allow_placeholders=True))
    guard_qc(construct, None)


def test_make_registrar_defaults_to_dry_run_without_credentials(tmp_path):
    registrar = make_registrar(BenchlingConfig(folder_id="lib_x"), out_dir=tmp_path)
    assert isinstance(registrar, DryRunRegistrar)


def test_make_registrar_honours_an_explicit_dry_run(tmp_path):
    config = BenchlingConfig(tenant="https://x.benchling.com", api_key="k",
                             folder_id="lib_x")
    assert config.has_credentials
    assert isinstance(make_registrar(config, dry_run=True, out_dir=tmp_path),
                      DryRunRegistrar)


def test_config_from_env(monkeypatch):
    monkeypatch.setenv("BENCHLING_TENANT", "https://x.benchling.com")
    monkeypatch.setenv("BENCHLING_API_KEY", "secret")
    monkeypatch.setenv("BENCHLING_FOLDER_ID", "lib_env")
    config = BenchlingConfig.from_env()
    assert config.folder_id == "lib_env"
    assert config.has_credentials
    assert config.missing() == []
    # Explicit overrides win; None never clobbers the environment.
    assert BenchlingConfig.from_env(folder_id="lib_override").folder_id == "lib_override"
    assert BenchlingConfig.from_env(folder_id=None).folder_id == "lib_env"


def test_missing_lists_what_is_absent():
    gaps = BenchlingConfig().missing()
    assert "BENCHLING_TENANT" in gaps
    assert "BENCHLING_FOLDER_ID" in gaps


def test_live_registrar_refuses_incomplete_configuration():
    from txv.benchling import BenchlingRegistrar

    with pytest.raises(ValueError, match="incomplete"):
        BenchlingRegistrar(BenchlingConfig())
