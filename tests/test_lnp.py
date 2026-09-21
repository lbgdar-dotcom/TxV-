import pytest

from txv.lnp import (
    COMPOSITIONS,
    LIPIDS,
    LipidComposition,
    M1PSI_DELTA,
    dose_series,
    formulate,
    phosphate_moles,
    rna_molecular_weight,
)


def test_molecular_weight_is_additive_over_residues():
    single = rna_molecular_weight("A", five_prime="p")
    double = rna_molecular_weight("AA", five_prime="p")
    assert double - single == pytest.approx(329.21, abs=0.01)


def test_triphosphate_and_hydroxyl_ends_differ_by_three_phosphates():
    ppp = rna_molecular_weight("ACGU", five_prime="ppp")
    p = rna_molecular_weight("ACGU", five_prime="p")
    oh = rna_molecular_weight("ACGU", five_prime="oh")
    assert ppp - p == pytest.approx(2 * 79.98, abs=0.01)
    assert p - oh == pytest.approx(79.98, abs=0.01)
    with pytest.raises(ValueError):
        rna_molecular_weight("ACGU", five_prime="cap9")


def test_dna_and_rna_alphabets_agree():
    assert rna_molecular_weight("ACGT") == rna_molecular_weight("ACGU")


def test_modified_uridine_adds_mass_per_uridine():
    plain = rna_molecular_weight("AUUUG")
    modified = rna_molecular_weight("AUUUG", modified_u_delta=M1PSI_DELTA)
    assert modified - plain == pytest.approx(3 * M1PSI_DELTA, abs=0.01)


def test_phosphate_moles_counts_one_per_nucleotide():
    sequence = "ACGU" * 100
    mw = rna_molecular_weight(sequence)
    moles_rna = 100e-6 / mw
    assert phosphate_moles(100.0, sequence) == pytest.approx(moles_rna * 400)
    with pytest.raises(ValueError):
        phosphate_moles(100.0)


def test_np_ratio_drives_ionizable_lipid_stoichiometry():
    sequence = "ACGU" * 250
    low = formulate(sequence, 100.0, np_ratio=3.0)
    high = formulate(sequence, 100.0, np_ratio=6.0)
    low_sm = next(a for a in low.lipids if a.lipid == "SM-102")
    high_sm = next(a for a in high.lipids if a.lipid == "SM-102")
    assert high_sm.micromoles == pytest.approx(2 * low_sm.micromoles)
    # N/P is exactly moles ionisable N over moles phosphate.
    assert high_sm.micromoles / high.phosphate_micromoles == pytest.approx(6.0)


def test_molar_percentages_are_preserved():
    result = formulate("ACGU" * 250, 100.0, composition="sm102_standard")
    total = sum(a.micromoles for a in result.lipids)
    for amount in result.lipids:
        assert amount.micromoles / total * 100 == pytest.approx(amount.molar_percent)


def test_masses_follow_from_moles_and_molecular_weight():
    result = formulate("ACGU" * 250, 100.0)
    for amount in result.lipids:
        assert amount.micrograms == pytest.approx(
            amount.micromoles * LIPIDS[amount.lipid].mw
        )
    assert result.total_lipid_micrograms == pytest.approx(
        sum(a.micrograms for a in result.lipids)
    )
    assert result.lipid_to_mrna_mass_ratio == pytest.approx(
        result.total_lipid_micrograms / 100.0
    )


def test_stock_volumes_and_mixing_volumes():
    stocks = {"SM-102": 10.0, "DSPC": 10.0, "cholesterol": 20.0, "DMG-PEG2000": 10.0}
    result = formulate("ACGU" * 250, 100.0, lipid_stocks_mg_per_ml=stocks,
                       flow_rate_ratio=3.0)
    for amount in result.lipids:
        # ug / (mg/mL) -> uL
        assert amount.stock_volume_ul == pytest.approx(
            amount.micrograms / stocks[amount.lipid]
        )
    assert result.organic_volume_ul == pytest.approx(
        sum(a.stock_volume_ul for a in result.lipids)
    )
    assert result.aqueous_volume_ul == pytest.approx(3 * result.organic_volume_ul)


def test_without_stocks_there_are_no_mixing_volumes():
    result = formulate("ACGU" * 250, 100.0)
    assert result.organic_volume_ul is None
    assert result.aqueous_volume_ul is None


def test_all_shipped_compositions_are_valid():
    for name, comp in COMPOSITIONS.items():
        assert sum(comp.molar_ratios.values()) == pytest.approx(100.0, abs=0.5)
        assert comp.ionizable.ionizable_nitrogens >= 1
        assert {LIPIDS[n].role for n in comp.molar_ratios} >= {"ionizable", "sterol"}


def test_bad_composition_is_rejected():
    with pytest.raises(ValueError):
        LipidComposition("bad", {"SM-102": 50.0, "DSPC": 10.0})
    with pytest.raises(ValueError):
        LipidComposition("unknown", {"NOT-A-LIPID": 100.0})


def test_composition_without_an_ionizable_lipid_is_rejected():
    comp = LipidComposition("nolipid", {"DSPC": 50.0, "cholesterol": 50.0})
    with pytest.raises(ValueError):
        _ = comp.ionizable


def test_dose_series_arithmetic():
    rows = dose_series(0.5, [0.020, 0.025], concentration_ug_per_ul=0.1)
    assert rows[0]["dose_ug"] == pytest.approx(10.0)
    assert rows[0]["volume_ul"] == pytest.approx(100.0)
    assert rows[1]["dose_ug"] == pytest.approx(12.5)


def test_result_serialisation_round_trips_the_numbers():
    result = formulate("ACGU" * 250, 100.0, construct_name="X")
    payload = result.to_dict()
    assert payload["construct"] == "X"
    assert len(payload["lipids"]) == 4
    assert payload["total_lipid_ug"] == pytest.approx(result.total_lipid_micrograms)
    assert "N/P" in result.to_text()
