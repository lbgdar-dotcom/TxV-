"""Negative controls for the pre-order audit.

An audit that passes is only meaningful if it is capable of failing. Each test
copies the real artifacts, introduces one realistic defect, and asserts the
audit blocks on it. If a mutation stops being caught, this suite fails -- which
is the point.
"""

from __future__ import annotations

import csv
import re
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_panel  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(
    not (ROOT / "orders" / "pvax1_ag_panel" / "gblocks.fasta").exists(),
    reason="generated order package not present",
)


@pytest.fixture
def sandbox(tmp_path):
    """A writable copy of everything the audit reads."""
    for rel in ("orders/pvax1_ag_panel", "backbone"):
        shutil.copytree(ROOT / rel, tmp_path / rel)
    return tmp_path


def _fasta(root: Path) -> Path:
    return root / "orders" / "pvax1_ag_panel" / "gblocks.fasta"


def mutate_fragment(root: Path, name: str, fn):
    """Rewrite one construct's sequence in BOTH shipped representations."""
    path = _fasta(root)
    out, current, buf = [], None, []

    def flush():
        if current is None:
            return
        seq = "".join(buf)
        if current.startswith(f">{name}|"):
            seq = fn(seq)
        out.append(current)
        out.extend(seq[i:i + 60] for i in range(0, len(seq), 60))

    for line in path.read_text().splitlines():
        if line.startswith(">"):
            flush()
            current, buf = line, []
        else:
            buf.append(line.strip())
    flush()
    path.write_text("\n".join(out) + "\n")

    table = root / "orders" / "pvax1_ag_panel" / "ordering_table.csv"
    rows = list(csv.DictReader(table.open()))
    for row in rows:
        if row["name"] == name:
            row["sequence"] = fn(row["sequence"])
            row["length_bp"] = str(len(row["sequence"]))
    with table.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)



def epitope_dna_span(fragment: str, peptide: str) -> tuple[int, int]:
    """Locate a peptide in the fragment by translating, not by guessing codons."""
    arm = len(audit_panel.LEFT_ARM)
    orf = fragment[arm: len(fragment) - len(audit_panel.RIGHT_ARM)]
    aa = audit_panel.tr(orf)
    i = aa.find(peptide)
    assert i >= 0, f"{peptide} not found in the translated ORF"
    return arm + i * 3, arm + (i + len(peptide)) * 3


def audit_blocks(root: Path) -> bool:
    return audit_panel.main(root, quiet=True) == 1


# --- the control: unmutated artifacts must pass ---------------------------

def test_unmutated_artifacts_pass(sandbox):
    assert not audit_blocks(sandbox), "the real artifacts should pass"


# --- sequence-level defects ------------------------------------------------

def test_catches_single_base_deletion(sandbox):
    """One lost base frameshifts everything downstream."""
    mutate_fragment(sandbox, "P3", lambda s: s[:120] + s[121:])
    assert audit_blocks(sandbox)


def test_catches_single_base_insertion(sandbox):
    mutate_fragment(sandbox, "P5", lambda s: s[:90] + "A" + s[90:])
    assert audit_blocks(sandbox)


def test_catches_premature_stop_codon(sandbox):
    def stop_it(s):
        i = 20 + 30 * 3  # well inside the ORF, in frame
        return s[:i] + "TGA" + s[i + 3:]
    mutate_fragment(sandbox, "P3", stop_it)
    assert audit_blocks(sandbox)


def test_catches_a_damaged_epitope(sandbox):
    """SIINFEKL is the readout; a point mutation in it makes the arm meaningless."""
    def break_epitope(s):
        start, _ = epitope_dna_span(s, "SIINFEKL")
        # SIINFEKL -> SAINFEKL: one anchor residue, enough to abolish binding
        return s[:start + 3] + "GCG" + s[start + 6:]
    mutate_fragment(sandbox, "P1", break_epitope)
    assert audit_blocks(sandbox)


def test_catches_a_lost_start_codon(sandbox):
    mutate_fragment(sandbox, "P6", lambda s: s[:20] + "GTG" + s[23:])
    assert audit_blocks(sandbox)


# --- positional constraints ------------------------------------------------

def test_catches_residues_appended_after_GYQTI(sandbox):
    """Anything after GYQTI abolishes AP-3 binding, silently."""
    def append_after_tail(s):
        arm = audit_panel.RIGHT_ARM
        orf = s[: -len(arm)]
        return orf[:-3] + "GCCGCC" + orf[-3:] + arm
    mutate_fragment(sandbox, "P6", append_after_tail)
    assert audit_blocks(sandbox)


def test_catches_e5_losing_its_free_c_terminus(sandbox):
    def append_after_degron(s):
        arm = audit_panel.RIGHT_ARM
        orf = s[: -len(arm)]
        return orf[:-3] + "GGCGGC" + orf[-3:] + arm
    mutate_fragment(sandbox, "P5", append_after_degron)
    assert audit_blocks(sandbox)


# --- cloning defects -------------------------------------------------------

def test_catches_a_bsai_site_in_the_insert(sandbox):
    """A BsaI site inside the cassette cuts the IVT template in half."""
    def insert_site(s):
        i = 20 + 60
        return s[:i] + "GGTCTC" + s[i + 6:]
    mutate_fragment(sandbox, "P3", insert_site)
    assert audit_blocks(sandbox)


def test_catches_a_broken_homology_arm(sandbox):
    mutate_fragment(sandbox, "P2", lambda s: "GAAGAAATATAAGAGCCAAA" + s[20:])
    assert audit_blocks(sandbox)


def test_catches_a_plasmid_map_that_does_not_match_its_gblock(sandbox):
    """The map and the fragment must describe the same molecule."""
    gb = sandbox / "orders" / "pvax1_ag_panel" / "plasmids" / "pVax1_AG_P3.gb"
    text = gb.read_text()
    head, body = text.split("ORIGIN", 1)
    # flip one base in the sequence block
    body = body.replace("a", "t", 1)
    gb.write_text(head + "ORIGIN" + body)
    assert audit_blocks(sandbox)


def test_catches_the_two_shipped_files_disagreeing(sandbox):
    table = sandbox / "orders" / "pvax1_ag_panel" / "ordering_table.csv"
    rows = list(csv.DictReader(table.open()))
    for row in rows:
        if row["name"] == "P4":
            row["sequence"] = row["sequence"][:-1] + ("A" if row["sequence"][-1] != "A" else "C")
    with table.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    assert audit_blocks(sandbox)


# --- nuclear-transcription hazard ------------------------------------------

def test_catches_a_polyadenylation_signal_in_the_orf(sandbox):
    def insert_signal(s):
        i = 20 + 90
        return s[:i] + "AATAAA" + s[i + 6:]
    mutate_fragment(sandbox, "P5", insert_signal)
    assert audit_blocks(sandbox)


# --- the confound that started this ----------------------------------------

def test_catches_a_module_encoded_differently_in_two_constructs(sandbox):
    """The defect that makes routing and codon usage indistinguishable."""
    def resynonymise_antigen_A(s):
        # A silent change inside antigen A, in P3 only: the protein is
        # unchanged, so only the cross-construct encoding check can see it.
        start, _ = epitope_dna_span(s, "SIINFEKL")
        leu = start + 21                       # the L of SIINFEKL
        assert audit_panel.tr(s[leu:leu + 3]) == "L"
        swapped = "CTC" if s[leu:leu + 3] != "CTC" else "CTG"
        return s[:leu] + swapped + s[leu + 3:]
    mutate_fragment(sandbox, "P3", resynonymise_antigen_A)
    assert audit_blocks(sandbox)


# --- and the backbone itself -----------------------------------------------

def test_catches_a_corrupted_backbone(sandbox):
    gb = sandbox / "backbone" / "pvax1_ag_egfp.gb"
    head, body = gb.read_text().split("ORIGIN", 1)
    gb.write_text(head + "ORIGIN" + body.replace("gacattgatt", "gacattgat", 1))
    assert audit_blocks(sandbox)
