"""Regression lock on the audited P0-P8 build.

REFERENCE_PROTEINS are the translations from the verified build. If a module in
:mod:`txv.pvax1_ag` is edited, these tests fail -- which is the intent.
"""

import pytest

from txv.pvax1_ag import (
    MODULES,
    OPEN_DECISIONS,
    PANEL,
    PANEL_BY_NAME,
    build_panel_construct,
    module_registry,
    panel_cassette,
    panel_proteins,
)
from txv.qc import run_qc

REFERENCE_PROTEINS = {
    "P0": "MACLGLRRYKAQLQLPSRTWPFVALLTLLFIPVFSYPYDVPDYAGGGGSFLLWILVAVSLGLFFYS"
          "FLVSAVSLSKMLKKRSPLTTGVYVKMPPTEPECEKQFQPYFIPINEEEEE",
    "P1": "MACLGLRRYKAQLQLPSRTWPFVALLTLLFIPVFSYPYDVPDYAGGGGSEVSGLEQLESIINFEKL"
          "TEWTSSNVMEERGGGGSFLLWILVAVSLGLFFYSFLVSAVSLSKMLKKRSPLTTGVYVKMPPTEPE"
          "CEKQFQPYFIPINEEEEE",
    "P2": "MACLGLRRYKAQLQLPSRTWPFVALLTLLFIPVFSYPYDVPDYAGGGGSEEFAKFASFEAQGALAN"
          "IAVDKANLDVMKGGGGSFLLWILVAVSLGLFFYSFLVSAVSLSKMLKKRSPLTTGVYVKMPPTEPE"
          "CEKQFQPYFIPINEEEEE",
    "P3": "MACLGLRRYKAQLQLPSRTWPFVALLTLLFIPVFSYPYDVPDYAGGGGSEVSGLEQLESIINFEKL"
          "TEWTSSNVMEERGGGGSEEFAKFASFEAQGALANIAVDKANLDVMKGGGGSFLLWILVAVSLGLFF"
          "YSFLVSAVSLSKMLKKRSPLTTGVYVKMPPTEPECEKQFQPYFIPINEEEEE",
    "P4": "MACLGLRRYKAQLQLPSRTWPFVALLTLLFIPVFSYPYDVPDYAGGGGSEEFAKFASFEAQGALAN"
          "IAVDKANLDVMKGGGGSEVSGLEQLESIINFEKLTEWTSSNVMEERGGGGSFLLWILVAVSLGLFF"
          "YSFLVSAVSLSKMLKKRSPLTTGVYVKMPPTEPECEKQFQPYFIPINEEEEE",
    "P5": "MYPYDVPDYAGGGGSEVSGLEQLESIINFEKLTEWTSSNVMEERGGGGSEEFAKFASFEAQGALAN"
          "IAVDKANLDVMKGGGGSEEEEE",
    "P6": "MAAPGARRPLLLLLLAGLAHGASAYPYDVPDYAGGGGSEVSGLEQLESIINFEKLTEWTSSNVMEE"
          "RGGGGSEEFAKFASFEAQGALANIAVDKANLDVMKGGGGSMLIPIAVGGALAGLVLIVLIAYLIGR"
          "KRSHAGYQTI",
    "P7": "MYPYDVPDYAGGGGSEVSGLEQLESIINFEKLTEWTSSNVMEERGGGGSEEEEEGSGATNFSLLKQ"
          "AGDVEENPGPMAAPGARRPLLLLLLAGLAHGASADYKDDDDKGGGGSEEFAKFASFEAQGALANIA"
          "VDKANLDVMKGGGGSMLIPIAVGGALAGLVLIVLIAYLIGRKRSHAGYQTI",
    "P8": "MYPYDVPDYAGGGGSEEFAKFASFEAQGALANIAVDKANLDVMKGGGGSEEEEEGSGATNFSLLKQ"
          "AGDVEENPGPMAAPGARRPLLLLLLAGLAHGASADYKDDDDKGGGGSEVSGLEQLESIINFEKLTE"
          "WTSSNVMEERGGGGSMLIPIAVGGALAGLVLIVLIAYLIGRKRSHAGYQTI",
}

#: Protein lengths from the audited construct_architecture table.
REFERENCE_LENGTHS = {
    "P0": 116, "P1": 150, "P2": 150, "P3": 184, "P4": 184,
    "P5": 88, "P6": 142, "P7": 183, "P8": 183,
}

#: ORF and finished-plasmid lengths from the same table. These are the
#: strongest external check available: they close only if the module content,
#: the stop-codon convention and the insertion site are all right.
REFERENCE_ORF_NT = {
    "P0": 351, "P1": 453, "P2": 453, "P3": 555, "P4": 555,
    "P5": 267, "P6": 429, "P7": 552, "P8": 552,
}
REFERENCE_PLASMID_BP = {
    "P0": 3486, "P1": 3588, "P2": 3588, "P3": 3690, "P4": 3690,
    "P5": 3402, "P6": 3564, "P7": 3687, "P8": 3687,
}


@pytest.mark.parametrize("name", sorted(REFERENCE_PROTEINS))
def test_panel_reconstructs_the_audited_translation(name):
    assert panel_proteins()[name] == REFERENCE_PROTEINS[name]


@pytest.mark.parametrize("name", sorted(REFERENCE_LENGTHS))
def test_panel_lengths_match_the_audit(name):
    assert len(PANEL_BY_NAME[name]) == REFERENCE_LENGTHS[name]


def test_panel_covers_p0_to_p8():
    assert [p.name for p in PANEL] == [f"P{i}" for i in range(9)]


def test_signal_peptides_are_n_terminal():
    """SRP reads the first thing out of the ribosome, so an SP works nowhere else."""
    for panel in PANEL:
        for index, module in enumerate(panel.modules):
            if module.endswith("_SP"):
                # An SP is either the first module, or the first module after a
                # P2A -- which starts a new polypeptide.
                assert index == 0 or panel.modules[index - 1] == "P2A", (
                    f"{panel.name}: {module} at index {index} is not at the "
                    "N-terminus of a polypeptide"
                )


def test_lamp1_tail_is_always_last():
    """GYQTI must end the protein or AP-3 (mu3A) binding fails."""
    for panel in PANEL:
        if "LAMP1_TMT" in panel.modules:
            assert panel.modules[-1] == "LAMP1_TMT", panel.name
            assert panel.protein().endswith("GYQTI"), panel.name


def test_ctla4_tail_keeps_its_native_residues_after_yvkm():
    """YVKM is internal; that is what lets the CTLA-4 arm sit upstream of P2A."""
    tail = MODULES["CTLA4_TMT"]
    assert "YVKM" in tail
    assert len(tail) - (tail.index("YVKM") + 4) == 19


def test_dual_route_constructs_put_lamp1_downstream_of_p2a():
    for name in ("P7", "P8"):
        modules = PANEL_BY_NAME[name].modules
        assert "P2A" in modules
        assert modules.index("LAMP1_SP") > modules.index("P2A")
        assert modules.index("LAMP1_TMT") > modules.index("P2A")


def test_p7_and_p8_are_antigen_swaps_of_each_other():
    p7, p8 = PANEL_BY_NAME["P7"], PANEL_BY_NAME["P8"]
    swap = {"A": "B", "B": "A"}
    assert tuple(swap.get(m, m) for m in p7.modules) == p8.modules


def test_p3_and_p4_are_order_controls():
    p3, p4 = PANEL_BY_NAME["P3"], PANEL_BY_NAME["P4"]
    assert p3.antigens == ("A", "B") and p4.antigens == ("B", "A")
    assert len(p3) == len(p4)


def test_scaffold_control_carries_no_antigen():
    assert PANEL_BY_NAME["P0"].antigens == ()
    for antigen in ("A", "B"):
        assert MODULES[antigen] not in PANEL_BY_NAME["P0"].protein()


def test_every_panel_member_carries_a_detection_tag():
    for panel in PANEL:
        assert "HA" in panel.modules, panel.name


@pytest.mark.parametrize("name", ["P0", "P5", "P6", "P7"])
def test_panel_builds_into_a_valid_ivt_construct(name):
    construct = build_panel_construct(name)
    assert construct.orf_protein == REFERENCE_PROTEINS[name]
    assert construct.orf.startswith("ATG")
    assert len(construct.orf) % 3 == 0
    report = run_qc(construct, allow_placeholders=True)
    assert report.passed, report.to_text()


def test_panel_construct_annotates_every_module():
    construct = build_panel_construct("P7")
    labels = {f.name for f in construct.features}
    for module in ("M", "HA", "A", "E5", "P2A", "LAMP1_SP", "FLAG", "B",
                   "LAMP1_TMT"):
        assert module in labels or any(
            label.startswith(f"{module}_") for label in labels
        ), module


def test_panel_cassette_uses_no_inter_module_linker():
    cassette = panel_cassette("P3")
    assert cassette.linker == ""
    assert cassette.sequence() == REFERENCE_PROTEINS["P3"]


def test_module_registry_namespaces_and_keeps_the_generic_parts():
    registry = module_registry()
    assert "pvax1_CTLA4_TMT" in registry
    assert "T7_promoter" in registry
    assert registry.get("pvax1_HA").protein == "YPYDVPDYA"
    assert registry.get("pvax1_backbone_R").dna == "GGTGGCTCTTATATTTCTTCTTACT"


def test_assembly_overlap_is_the_reverse_complement_of_the_backbone_primer():
    """insert_left_overlap must anneal to where backbone_R leaves off."""
    from txv.pvax1_ag import BACKBONE_OLIGOS
    from txv.seqops import revcomp

    left = BACKBONE_OLIGOS["insert_left_overlap"][0]
    backbone_r = BACKBONE_OLIGOS["backbone_R"][0]
    assert left.endswith("GCCACC"), "left overlap must end in the Kozak"
    # backbone_R anneals across the 5' UTR boundary *and* the Kozak, so the
    # insert's left overlap is exactly the 3' end of its reverse complement.
    # If that stops being true, the assembly loses its 20-bp homology.
    assert revcomp(backbone_r).endswith(left)
    assert len(left) == 20, "HiFi/In-Fusion assembly expects a 20-bp overlap"
    assert (
        BACKBONE_OLIGOS["insert_right_overlap"][0]
        == BACKBONE_OLIGOS["backbone_F"][0]
    )


def test_open_decisions_are_recorded():
    assert OPEN_DECISIONS
    joined = " ".join(OPEN_DECISIONS)
    assert "A120" in joined
    assert "AAGCTT" in joined or "HindIII" in joined


# --- the degron ------------------------------------------------------------

def test_the_panel_uses_e5_not_cl1():
    """E5 is the element CVGBM carries and the one that won CureVac's screen."""
    assert MODULES["E5"] == "EEEEE"
    for panel in PANEL:
        assert "CL1" not in panel.modules, panel.name
    with_degron = [p.name for p in PANEL if "E5" in p.modules]
    assert with_degron == ["P0", "P1", "P2", "P3", "P4", "P5", "P7", "P8"]
    # P6 is the one construct with no degron at all.
    assert "E5" not in PANEL_BY_NAME["P6"].modules


def test_cl1_is_retained_only_as_a_documented_alternative():
    from txv.pvax1_ag import MODULE_NOTES

    assert MODULES["CL1"] == "ACKNWFSSLSHFVIHL"
    assert "NOT USED" in MODULE_NOTES["CL1"]
    assert "withdrawn" in MODULE_NOTES["CL1"]


def test_e5_has_a_free_c_terminus_only_where_it_can_work():
    """A C-degron needs a free C-terminus; upstream of P2A it does not have one."""
    for name in ("P0", "P1", "P2", "P3", "P4", "P5"):
        assert PANEL_BY_NAME[name].protein().endswith("EEEEE"), name
    for name in ("P7", "P8"):
        panel = PANEL_BY_NAME[name]
        assert not panel.protein().endswith("EEEEE")
        assert panel.modules.index("E5") < panel.modules.index("P2A")
    assert "weak" in PANEL_BY_NAME["P7"].question


def test_e5_codons_are_pinned_against_their_own_repetitiveness():
    from txv.pvax1_ag import PINNED_DNA

    assert PINNED_DNA["E5"] == "GAAGAGGAAGAGGAG"
    construct = build_panel_construct("P5")
    assert construct.feature("E5").slice(construct.template) == PINNED_DNA["E5"]
    # The pin exists to break periodicity: no 12-nt repeat anywhere in the ORF.
    orf, seen, repeats = construct.orf, {}, []
    for i in range(len(orf) - 12 + 1):
        kmer = orf[i : i + 12]
        if kmer in seen:
            repeats.append(kmer)
        seen[kmer] = i
    assert not repeats, repeats


@pytest.mark.parametrize("name", sorted(REFERENCE_ORF_NT))
def test_orf_length_matches_the_audited_table(name):
    assert len(build_panel_construct(name).orf) == REFERENCE_ORF_NT[name]


def test_single_stop_codon_matches_the_audited_convention():
    construct = build_panel_construct("P3")
    assert len(construct.feature("stop")) == 3
