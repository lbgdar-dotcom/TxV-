import json

import pytest

from txv.epitopes import Antigen
from txv.pipeline import design_construct, load_antigens, write_outputs


def test_design_end_to_end(antigens):
    report = design_construct("TXV-X", antigens)
    assert report.construct.name == "TXV-X"
    assert report.qc.passed
    assert report.chosen_linker
    assert set(report.linker_table) >= {"GGSGGGGSGG", "AAY"}
    assert "TXV-X" in report.to_text()


def test_design_picks_the_lowest_risk_linker(antigens):
    report = design_construct("X", antigens)
    assert report.linker_table[report.chosen_linker] == min(
        report.linker_table.values()
    )


def test_forcing_a_linker_skips_selection(antigens):
    report = design_construct("X", antigens, linker="AAY")
    assert report.chosen_linker == "AAY"
    assert report.linker_table == {}
    assert report.construct.cassette.linker == "AAY"


def test_reordering_can_be_disabled(antigens):
    report = design_construct("X", antigens, reorder=False, linker="AAY")
    assert [a.name for a in report.construct.cassette.antigens] == [
        a.name for a in antigens
    ]


def test_dedupe_is_reported(antigens):
    duplicated = list(antigens) + [Antigen("dupe", antigens[2].sequence)]
    report = design_construct("X", duplicated)
    assert report.dropped_antigens
    assert len(report.construct.cassette.antigens) == len(antigens)


def test_strict_mode_fails_on_placeholders(antigens):
    report = design_construct("X", antigens, allow_placeholders=False)
    assert not report.qc.passed


def test_empty_antigen_list_is_rejected():
    with pytest.raises(ValueError):
        design_construct("X", [])


def test_report_serialisation(antigens):
    payload = design_construct("X", antigens).to_dict()
    assert payload["design"]["order"]
    assert payload["construct"]["name"] == "X"
    assert payload["qc"]["passed"] is True


def test_write_outputs(antigens, tmp_path):
    report = design_construct("X", antigens)
    written = write_outputs(report, tmp_path)
    assert set(written) == {
        "X.gb", "X.template.fasta", "X.mrna.fasta", "X.protein.fasta",
        "X.report.json", "X.report.txt",
    }
    assert (tmp_path / "X.gb").read_text().startswith("LOCUS")
    assert json.loads((tmp_path / "X.report.json").read_text())["construct"]["name"] == "X"


def test_load_antigens_from_csv(tmp_path):
    path = tmp_path / "a.csv"
    path.write_text(
        "name,sequence,kind,hla,mutation_offset,note\n"
        "one,SLLMWITQC,tumor_associated,HLA-A*02:01,,first\n"
        "two,KVAELVHFL,neoepitope,HLA-A*02:01;HLA-A*11:01,3,second\n"
    )
    antigens = load_antigens(path)
    assert [a.name for a in antigens] == ["one", "two"]
    assert antigens[0].mutation_offset is None
    assert antigens[1].mutation_offset == 3
    assert antigens[1].hla == ("HLA-A*02:01", "HLA-A*11:01")


def test_load_antigens_from_json(tmp_path):
    path = tmp_path / "a.json"
    path.write_text(json.dumps({"antigens": [
        {"name": "one", "sequence": "SLLMWITQC"},
        {"name": "two", "sequence": "KVAELVHFL", "kind": "helper"},
    ]}))
    antigens = load_antigens(path)
    assert len(antigens) == 2
    assert antigens[1].kind == "helper"
    assert antigens[0].kind == "neoepitope"


def test_load_antigens_from_a_bare_json_list(tmp_path):
    path = tmp_path / "a.json"
    path.write_text(json.dumps([{"name": "one", "sequence": "SLLMWITQC"}]))
    assert len(load_antigens(path)) == 1


def test_csv_rows_without_a_sequence_are_skipped(tmp_path):
    path = tmp_path / "a.csv"
    path.write_text("name,sequence\none,SLLMWITQC\n,\nthree,\n")
    assert [a.name for a in load_antigens(path)] == ["one"]


def test_example_antigen_table_loads():
    antigens = load_antigens("examples/antigens_example.csv")
    assert len(antigens) == 8
    assert any(a.kind == "helper" for a in antigens)
