"""Reading a SnapGene vector, and dropping a blunt fragment into it."""

import pytest

from txv.genbank_io import GenBankFeature, GenBankRecord, read_genbank
from txv.plasmid import insert_blunt
from txv.seqops import find_all, revcomp
from txv.snapgene import read_snapgene

VECTOR = "backbone/pJET1.2.gb"


@pytest.fixture(scope="module")
def vector():
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / VECTOR
    if not path.exists():
        pytest.skip("pJET1.2 map not generated")
    return read_genbank(path)


# --- the vector itself -----------------------------------------------------

def test_vector_is_the_expected_size_and_topology(vector):
    assert len(vector) == 2974
    assert vector.is_circular


def test_vector_has_exactly_one_blunt_cloning_site(vector):
    """pJET1.2's whole premise is a single Eco32I site to drop a blunt end into."""
    assert len(find_all(vector.sequence, "GATATC")) == 1


def test_vector_has_the_two_flanking_excision_sites(vector):
    bgl = find_all(vector.sequence, "AGATCT")
    assert len(bgl) == 2
    ecorv = find_all(vector.sequence, "GATATC")[0]
    assert bgl[0] < ecorv < bgl[1]


def test_the_lethal_gene_spans_the_cloning_site(vector):
    """Positive selection only works if inserting breaks eco47IR."""
    lethal = next(f for f in vector.features
                  if "Eco47I" in f.qualifiers.get("label", ""))
    ecorv = find_all(vector.sequence, "GATATC")[0]
    assert lethal.start < ecorv < lethal.end


# --- blunt insertion -------------------------------------------------------

@pytest.fixture
def toy():
    record = GenBankRecord(
        name="toy", sequence="AAAA" + "GATATC" + "TTTT", is_circular=True,
        features=[
            GenBankFeature("misc_feature", 0, 4, 1, {"label": "before"}),
            GenBankFeature("misc_feature", 10, 14, 1, {"label": "after"}),
            GenBankFeature("CDS", 0, 14, 1, {"label": "spanning"}),
        ],
    )
    return record


def test_insertion_lengthens_by_exactly_the_fragment(toy):
    build = insert_blunt(toy, "CCCCC", 7, "x")
    assert len(build) == len(toy) + 5
    assert build.insert == "CCCCC"


def test_features_before_the_site_do_not_move(toy):
    build = insert_blunt(toy, "CCCCC", 7, "x")
    before = next(f for f in build.record.features
                  if f.qualifiers.get("label") == "before")
    assert (before.start, before.end) == (0, 4)


def test_features_after_the_site_shift_by_the_insert(toy):
    build = insert_blunt(toy, "CCCCC", 7, "x")
    after = next(f for f in build.record.features
                 if f.qualifiers.get("label") == "after")
    assert (after.start, after.end) == (15, 19)


def test_a_feature_spanning_the_site_is_flagged_disrupted(toy):
    build = insert_blunt(toy, "CCCCC", 7, "x")
    assert "spanning" in build.disrupted
    spanning = next(f for f in build.record.features
                    if f.qualifiers.get("label") == "spanning")
    assert "DISRUPTED" in spanning.qualifiers["note"]
    assert spanning.end == 19


def test_reverse_orientation_inserts_the_reverse_complement(toy):
    build = insert_blunt(toy, "CCCCCA", 7, "x", orientation=-1)
    assert build.insert == revcomp("CCCCCA")
    assert build.orientation == -1


def test_an_out_of_range_site_is_rejected(toy):
    with pytest.raises(ValueError):
        insert_blunt(toy, "CCC", 999, "x")


def test_an_invalid_orientation_is_rejected(toy):
    with pytest.raises(ValueError):
        insert_blunt(toy, "CCC", 7, "x", orientation=0)


# --- the real thing --------------------------------------------------------

@pytest.mark.parametrize("name", ["P0", "P5", "P7"])
def test_the_finished_plasmid_carries_the_fragment_intact(vector, name):
    from txv.ivt_unit import make_unit
    from txv.pvax1_ag import build_panel_construct

    unit = make_unit(name, build_panel_construct(name).orf)
    site = find_all(vector.sequence, "GATATC")[0] + 3
    build = insert_blunt(vector, unit.sequence, site, f"pJET1.2_{name}")

    assert build.insert == unit.sequence
    assert len(build) == len(vector) + len(unit)
    assert any("Eco47I" in d for d in build.disrupted)
    for keep in ("AmpR", "ori"):
        assert keep not in build.disrupted


@pytest.mark.parametrize("name", ["P0", "P5", "P7"])
def test_bglii_excises_the_whole_fragment(vector, name):
    """The insert must come out in one piece, or the digest cannot check it."""
    from txv.ivt_unit import make_unit
    from txv.pvax1_ag import build_panel_construct

    unit = make_unit(name, build_panel_construct(name).orf)
    assert find_all(unit.sequence, "AGATCT", both_strands=True) == []

    site = find_all(vector.sequence, "GATATC")[0] + 3
    build = insert_blunt(vector, unit.sequence, site, f"pJET1.2_{name}")
    cuts = sorted(i + 1 for i in find_all(build.record.sequence, "AGATCT"))
    assert len(cuts) == 2
    assert cuts[0] < build.insert_start and build.insert_end < cuts[1]
    assert cuts[1] - cuts[0] == len(unit) + 46


def test_the_orientation_pcr_distinguishes_the_two_products(vector):
    """A band means one orientation, no band the other. That is the test."""
    from txv.ivt_unit import FWD_PRIMER, make_unit, simulate_pcr
    from txv.pvax1_ag import build_panel_construct
    from txv.snapgene import read_primers

    unit = make_unit("P0", build_panel_construct("P0").orf)
    site = find_all(vector.sequence, "GATATC")[0] + 3
    rev = "AAGAACATCGATTTTCCATGGCAG"      # pJET1.2 reverse sequencing primer

    forward = insert_blunt(vector, unit.sequence, site, "f", orientation=1)
    reverse = insert_blunt(vector, unit.sequence, site, "r", orientation=-1)

    assert simulate_pcr(forward.record.sequence, FWD_PRIMER, rev)
    assert not simulate_pcr(reverse.record.sequence, FWD_PRIMER, rev)


# --- the encoded poly(A) and run-off transcription -------------------------

@pytest.fixture(scope="module")
def plasmid_dir():
    from pathlib import Path
    d = Path(__file__).resolve().parents[1] / "orders" / "pjet12_ivt_units" / "plasmids"
    if not (d / "pJET1.2_P0.gb").exists():
        pytest.skip("plasmid maps not generated")
    return d


@pytest.mark.parametrize("name", ["P0", "P5", "P7"])
def test_the_plasmid_carries_an_encoded_polya(plasmid_dir, name):
    record = read_genbank(plasmid_dir / f"pJET1.2_{name}.gb")
    tail = next(f for f in record.features
                if f.qualifiers.get("label", "").startswith("poly(A)"))
    assert set(record.sequence[tail.start:tail.end]) == {"A"}
    assert tail.end - tail.start == 120


@pytest.mark.parametrize("name", ["P0", "P5", "P7"])
def test_the_reverse_orientation_polya_sits_on_the_bottom_strand(plasmid_dir, name):
    """A -1 strand poly(A) reads as poly(T) on the strand the file lists."""
    record = read_genbank(plasmid_dir / f"pJET1.2_{name}_rev.gb")
    tail = next(f for f in record.features
                if f.qualifiers.get("label", "").startswith("poly(A)"))
    assert tail.strand == -1
    assert set(record.sequence[tail.start:tail.end]) == {"T"}


@pytest.mark.parametrize("name", ["P0", "P5", "P7"])
def test_the_gene_fragment_itself_encodes_no_polya(name):
    """The tail is added by PCR. A synthesised A120 would be rejected."""
    from txv.ivt_unit import make_unit
    from txv.pvax1_ag import build_panel_construct
    unit = make_unit(name, build_panel_construct(name).orf)
    assert "A" * 20 not in unit.sequence


@pytest.mark.parametrize("name", [p.name for p in __import__(
    "txv.pvax1_ag", fromlist=["PANEL"]).PANEL])
def test_bglii_runoff_gives_a_usable_ivt_template(plasmid_dir, name):
    """One digest must serve as both the diagnostic and the IVT template."""
    from txv.ivt_unit import T7_CORE
    record = read_genbank(plasmid_dir / f"pJET1.2_{name}.gb")
    seq = record.sequence
    cuts = sorted(i + 1 for i in find_all(seq, "AGATCT"))
    assert len(cuts) == 2

    released = seq[cuts[0]:cuts[1]]
    t7 = find_all(released, T7_CORE)
    assert len(t7) == 1, "the released fragment must carry exactly one promoter"

    transcript = released[t7[0] + len(T7_CORE):]
    assert transcript.startswith("AGG"), "must stay CleanCap AG compatible"
    assert "A" * 120 in transcript
    after = len(transcript) - (transcript.rfind("A" * 120) + 120)
    assert 0 <= after <= 20, f"{after} nt after the poly(A) is too much"


@pytest.mark.parametrize("name", ["P0", "P7"])
def test_the_plasmid_map_annotates_the_modules_not_just_a_block(plasmid_dir, name):
    """The point of these maps is seeing the design, not an opaque insert."""
    record = read_genbank(plasmid_dir / f"pJET1.2_{name}.gb")
    labels = {f.qualifiers.get("label", "") for f in record.features}
    for required in ("T7_promoter", "5'UTR", "Kozak", "ORF", "HA", "3'UTR"):
        assert required in labels, required
    assert any(l.startswith("poly(A)") for l in labels)
    # and the vector's own features must survive
    for required in ("AmpR", "ori", "Eco47I/T7"):
        assert required in labels, required
