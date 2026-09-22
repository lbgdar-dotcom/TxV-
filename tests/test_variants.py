import json

import pytest

from txv.codon import uridine_floor
from txv.seqops import translate
from txv.variants import multiplex_ladder, series_outputs, uridine_ladder


def test_uridine_ladder_holds_the_protein_constant(antigens):
    series = uridine_ladder(antigens[:3], weights=(0.0, 0.5, 2.0))
    proteins = {v.construct.orf_protein for v in series.variants}
    assert len(proteins) == 1, "arms must encode the identical protein"
    dnas = {v.construct.orf for v in series.variants}
    assert len(dnas) > 1, "arms must differ in nucleotide sequence"


def test_uridine_ladder_is_monotone_and_passes_qc(antigens):
    series = uridine_ladder(antigens[:3], weights=(0.0, 0.5, 2.0))
    fractions = [v.covariates["U_frac"] for v in series.variants]
    assert fractions == sorted(fractions, reverse=True)
    assert all(v.qc.passed for v in series.variants)


def test_uridine_ladder_reports_the_floor_it_cannot_cross(antigens):
    series = uridine_ladder(antigens[:3], weights=(0.0, 4.0))
    for variant in series.variants:
        floor = variant.covariates["floor_U"]
        assert floor == uridine_floor(variant.construct.orf_protein)[0]
        assert variant.covariates["ORF_U"] >= floor


def test_uridine_ladder_reports_cai_as_a_covariate(antigens):
    series = uridine_ladder(antigens[:3], weights=(0.0, 2.0))
    assert all("CAI" in v.covariates for v in series.variants)
    assert "confound" in series.note


def test_multiplex_ladder_is_nested(antigens):
    series = multiplex_ladder(antigens, sizes=(1, 2, 3, 5))
    sets = [
        {a.name for a in v.construct.cassette.antigens} for v in series.variants
    ]
    for smaller, larger in zip(sets, sets[1:]):
        assert smaller < larger, "each arm must contain the previous one"


def test_multiplex_ladder_keeps_the_anchor_first_in_every_arm(antigens):
    series = multiplex_ladder(antigens, sizes=(1, 2, 4))
    for variant in series.variants:
        assert variant.construct.cassette.antigens[0].name == antigens[0].name


def test_explicit_anchor_is_honoured(antigens):
    series = multiplex_ladder(antigens, sizes=(1, 3), anchor=antigens[2])
    for variant in series.variants:
        assert variant.construct.cassette.antigens[0].name == antigens[2].name
        names = [a.name for a in variant.construct.cassette.antigens]
        assert len(names) == len(set(names)), "anchor must not be duplicated"


def test_multiplex_arms_grow_and_pass_qc(antigens):
    series = multiplex_ladder(antigens, sizes=(1, 2, 4))
    lengths = [v.covariates["ORF_nt"] for v in series.variants]
    assert lengths == sorted(lengths)
    assert all(v.qc.passed for v in series.variants)
    assert all(v.value == n for v, n in zip(series.variants, (1, 2, 4)))


def test_multiplex_ladder_rejects_impossible_sizes(antigens):
    with pytest.raises(ValueError):
        multiplex_ladder(antigens, sizes=(len(antigens) + 5,))
    with pytest.raises(ValueError):
        multiplex_ladder(antigens, sizes=(0,))
    with pytest.raises(ValueError):
        multiplex_ladder([], sizes=(1,))


def test_series_arms_use_a_real_route_with_no_placeholders(antigens):
    for series in (
        uridine_ladder(antigens[:3], weights=(0.5,)),
        multiplex_ladder(antigens, sizes=(2,)),
    ):
        for variant in series.variants:
            assert variant.construct.placeholders_used == []


def test_route_choice_changes_the_protein(antigens):
    ctla4 = multiplex_ladder(antigens, sizes=(2,), route="ctla4")
    lamp1 = multiplex_ladder(antigens, sizes=(2,), route="lamp1")
    cytosolic = multiplex_ladder(antigens, sizes=(2,), route="cytosolic")
    proteins = [s.variants[0].construct.orf_protein for s in (ctla4, lamp1, cytosolic)]
    assert len(set(proteins)) == 3
    assert proteins[1].endswith("GYQTI"), "LAMP1 tail must end the protein"


def test_series_rendering_and_serialisation(antigens):
    series = multiplex_ladder(antigens, sizes=(1, 2))
    text = series.to_text()
    assert "n_antigens" in text and "N1" in text
    payload = series.to_dict()
    assert payload["axis"] == "n_antigens"
    assert len(payload["arms"]) == 2
    assert json.dumps(payload)


def test_series_outputs_writes_every_arm(antigens, tmp_path):
    series = multiplex_ladder(antigens, sizes=(1, 2))
    written = series_outputs(series, tmp_path)
    assert (tmp_path / "N1.gb").read_text().startswith("LOCUS")
    assert (tmp_path / "N2.mrna.fasta").exists()
    assert (tmp_path / "n_antigens_series.json").exists()
    assert len(written) == 5
