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
def test_the_flipped_clone_shows_a_readable_polya_too(plasmid_dir, name):
    """Shown from its other strand, the flipped clone's tail reads as A's.

    Built on the strand it was made on, a flipped clone's poly(A) is a -1
    feature spanning T's -- correct, and unreadable. These maps are drawn from
    the opposite strand precisely so that stops being the case.
    """
    record = read_genbank(
        plasmid_dir / "reverse_orientation" / f"pJET1.2_{name}_rev.gb")
    tail = next(f for f in record.features
                if f.qualifiers.get("label", "").startswith("poly(A)"))
    assert tail.strand == 1
    assert set(record.sequence[tail.start:tail.end]) == {"A"}
    assert tail.end - tail.start == 120


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

    # The insert now starts at position 1, so the fragment carrying it wraps
    # the origin: a plain slice between the two cuts returns the backbone.
    doubled = seq + seq
    released = doubled[cuts[1]:cuts[0] + len(seq)]
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
    # eco47IR is split by the rotation that puts the insert at position 1;
    # it is the feature the insert disrupted, so it straddles the new origin.
    for required in ("AmpR", "ori"):
        assert required in labels, required
    assert any(l.startswith("Eco47I/T7") for l in labels)


# --- map orientation -------------------------------------------------------
# A map is for reading a design off. If the cassette runs backwards it is not
# serving that purpose, whichever strand the record happens to be built on.

CASSETTE = ["IVT_F_handle", "T7_promoter", "5'UTR", "ORF", "3'UTR",
            "IVT_R_handle", "poly(A) 120"]


def _labelled(record):
    return {f.qualifiers.get("label", ""): f for f in record.features}


@pytest.mark.parametrize("name", ["P0", "P5", "P7"])
def test_the_designed_map_starts_at_the_inserts_five_prime_end(plasmid_dir, name):
    record = read_genbank(plasmid_dir / f"pJET1.2_{name}.gb")
    assert _labelled(record)["IVT_F_handle"].start == 0


@pytest.mark.parametrize("name", ["P0", "P5", "P7"])
def test_the_cassette_reads_downstream_in_order(plasmid_dir, name):
    record = read_genbank(plasmid_dir / f"pJET1.2_{name}.gb")
    found = _labelled(record)
    starts = [found[label].start for label in CASSETTE]
    assert starts == sorted(starts), "elements are out of order on the map"
    assert all(found[label].strand == 1 for label in CASSETTE)


@pytest.mark.parametrize("name", ["P0", "P5", "P7"])
def test_the_flipped_clone_map_also_reads_forwards(plasmid_dir, name):
    """Shown from its other strand, a flipped clone is readable too."""
    record = read_genbank(
        plasmid_dir / "reverse_orientation" / f"pJET1.2_{name}_rev.gb")
    found = _labelled(record)
    starts = [found[label].start for label in CASSETTE]
    assert starts == sorted(starts)
    assert all(found[label].strand == 1 for label in CASSETTE)


@pytest.mark.parametrize("name", ["P0", "P5", "P7"])
def test_the_two_clones_differ_in_the_vector_strand_not_the_cassette(
        plasmid_dir, name):
    """The real difference must be the visible one."""
    designed = _labelled(read_genbank(plasmid_dir / f"pJET1.2_{name}.gb"))
    flipped = _labelled(read_genbank(
        plasmid_dir / "reverse_orientation" / f"pJET1.2_{name}_rev.gb"))
    for backbone in ("AmpR", "ori"):
        assert designed[backbone].strand == -flipped[backbone].strand


# --- rotation and strand-flip primitives -----------------------------------

def test_rotation_preserves_the_molecule(vector):
    from txv.plasmid import rotate_record
    rotated = rotate_record(vector, 500)
    assert len(rotated) == len(vector)
    assert rotated.sequence in vector.sequence + vector.sequence
    assert rotated.sequence != vector.sequence


def test_rotation_splits_a_feature_that_crosses_the_new_origin(vector):
    from txv.plasmid import rotate_record
    lethal = next(f for f in vector.features
                  if "Eco47I" in f.qualifiers.get("label", ""))
    rotated = rotate_record(vector, lethal.start + 100)
    halves = [f for f in rotated.features
              if "Eco47I" in f.qualifiers.get("label", "")]
    assert len(halves) == 2
    assert {f.qualifiers["label"] for f in halves} == {
        "Eco47I/T7 (1 of 2)", "Eco47I/T7 (2 of 2)"}


def test_rotation_refuses_a_linear_record(vector):
    from txv.plasmid import rotate_record
    linear = GenBankRecord(name="x", sequence=vector.sequence,
                           is_circular=False)
    with pytest.raises(ValueError):
        rotate_record(linear, 10)


def test_rotating_by_zero_is_a_no_op(vector):
    from txv.plasmid import rotate_record
    assert rotate_record(vector, 0) is vector


def test_strand_flip_is_its_own_inverse(vector):
    from txv.plasmid import revcomp_record
    twice = revcomp_record(revcomp_record(vector))
    assert twice.sequence == vector.sequence
    assert {(f.start, f.end, f.strand) for f in twice.features} == \
           {(f.start, f.end, f.strand) for f in vector.features}
