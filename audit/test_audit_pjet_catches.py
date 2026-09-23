"""Negative controls for the pJET1.2 audit.

Each test removes or corrupts one of the elements you asked to be able to
verify -- a signal peptide, E5, P2A, an antigen, the LAMP1 tail -- and asserts
the audit blocks. If a mutation stops being caught, this suite fails.
"""

from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_pjet  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "orders" / "pjet12_ivt_units"
pytestmark = pytest.mark.skipif(
    not (PKG / "gblocks.fasta").exists(), reason="pJET package not generated"
)


@pytest.fixture
def sandbox(tmp_path):
    shutil.copytree(PKG, tmp_path / "orders" / "pjet12_ivt_units")
    return tmp_path


def blocks(root: Path) -> bool:
    return audit_pjet.main(root, quiet=True) == 1


def rewrite(root: Path, name: str, fn):
    """Apply fn to one construct's sequence in the fasta, the table and the .gb."""
    d = root / "orders" / "pjet12_ivt_units"

    path = d / "gblocks.fasta"
    out, header, buf = [], None, []

    def flush():
        if header is None:
            return
        seq = "".join(buf)
        if header.startswith(f">{name}|"):
            seq = fn(seq)
        out.append(header)
        out.extend(seq[i:i + 60] for i in range(0, len(seq), 60))

    for line in path.read_text().splitlines():
        if line.startswith(">"):
            flush()
            header, buf = line, []
        else:
            buf.append(line.strip())
    flush()
    path.write_text("\n".join(out) + "\n")

    table = d / "ordering_table.csv"
    rows = list(csv.DictReader(table.open()))
    for row in rows:
        if row["name"] == name:
            row["sequence"] = fn(row["sequence"])
            row["fragment_bp"] = str(len(row["sequence"]))
    with table.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    gb = d / f"{name}.gb"
    head, body = gb.read_text().split("ORIGIN", 1)
    seq = fn("".join(c for c in body.split("//")[0] if c.isalpha()).upper())
    lines = ["ORIGIN"]
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60].lower()
        lines.append(f"{i + 1:>9} " + " ".join(
            chunk[j:j + 10] for j in range(0, len(chunk), 10)))
    gb.write_text(head + "\n".join(lines) + "\n//\n")


def delete_peptide(sequence: str, peptide: str) -> str:
    """Remove a module from the ORF, in frame, leaving everything else."""
    dna = audit_pjet
    orf_start = sequence.index(audit_pjet.LEADER) + len(audit_pjet.LEADER)
    aa = dna.tr(sequence[orf_start:])
    i = aa.index(peptide)
    a = orf_start + i * 3
    b = a + len(peptide) * 3
    return sequence[:a] + sequence[b:]


# --- the control -----------------------------------------------------------

def test_unmutated_package_passes(sandbox):
    assert not blocks(sandbox)


# --- the elements you asked to verify --------------------------------------

def test_catches_a_missing_P2A(sandbox):
    """Without P2A, P7 is one fusion protein instead of two."""
    rewrite(sandbox, "P7", lambda s: delete_peptide(s, audit_pjet.MODULES["P2A"]))
    assert blocks(sandbox)


def test_catches_a_missing_E5_degron(sandbox):
    rewrite(sandbox, "P5", lambda s: delete_peptide(s, audit_pjet.MODULES["E5"]))
    assert blocks(sandbox)


def test_catches_a_missing_signal_peptide(sandbox):
    rewrite(sandbox, "P6", lambda s: delete_peptide(s, audit_pjet.MODULES["LAMP1_SP"]))
    assert blocks(sandbox)


def test_catches_a_missing_LAMP1_tail(sandbox):
    rewrite(sandbox, "P6", lambda s: delete_peptide(s, audit_pjet.MODULES["LAMP1_TMT"]))
    assert blocks(sandbox)


def test_catches_a_missing_antigen(sandbox):
    rewrite(sandbox, "P3", lambda s: delete_peptide(s, audit_pjet.MODULES["B"]))
    assert blocks(sandbox)


def test_catches_a_missing_HA_tag(sandbox):
    rewrite(sandbox, "P1", lambda s: delete_peptide(s, audit_pjet.MODULES["HA"]))
    assert blocks(sandbox)


def test_catches_swapped_antigen_order(sandbox):
    """P3 is A-then-B; swapping it silently turns P3 into P4."""
    def swap(s):
        a, b = audit_pjet.MODULES["A"], audit_pjet.MODULES["B"]
        start = s.index(audit_pjet.LEADER) + len(audit_pjet.LEADER)
        aa = audit_pjet.tr(s[start:])
        ia, ib = aa.index(a), aa.index(b)
        da = s[start + ia * 3: start + (ia + len(a)) * 3]
        db = s[start + ib * 3: start + (ib + len(b)) * 3]
        return (s[:start + ia * 3] + db
                + s[start + (ia + len(a)) * 3: start + ib * 3] + da
                + s[start + (ib + len(b)) * 3:])
    rewrite(sandbox, "P3", swap)
    assert blocks(sandbox)


# --- transcription-unit defects --------------------------------------------

def test_catches_a_broken_T7_promoter(sandbox):
    rewrite(sandbox, "P0", lambda s: s.replace(audit_pjet.T7_CORE,
                                               "TAATACGACTCACTATT"))
    assert blocks(sandbox)


def test_catches_a_lost_AGG_transcription_start(sandbox):
    """CleanCap AG needs the transcript to begin AG."""
    def break_start(s):
        i = s.index(audit_pjet.T7_CORE) + len(audit_pjet.T7_CORE)
        return s[:i] + "CGG" + s[i + 3:]
    rewrite(sandbox, "P0", break_start)
    assert blocks(sandbox)


def test_catches_a_broken_kozak(sandbox):
    rewrite(sandbox, "P2", lambda s: s.replace("GCCACCATG", "GCCAAAATG"))
    assert blocks(sandbox)


def test_catches_a_frameshift(sandbox):
    def shift(s):
        i = s.index(audit_pjet.LEADER) + len(audit_pjet.LEADER) + 30
        return s[:i] + s[i + 1:]
    rewrite(sandbox, "P4", shift)
    assert blocks(sandbox)


def test_catches_a_damaged_epitope(sandbox):
    def break_it(s):
        start = s.index(audit_pjet.LEADER) + len(audit_pjet.LEADER)
        aa = audit_pjet.tr(s[start:])
        i = aa.index("SIINFEKL")
        at = start + (i + 1) * 3
        return s[:at] + "GCG" + s[at + 3:]
    rewrite(sandbox, "P1", break_it)
    assert blocks(sandbox)


def test_catches_a_missing_handle(sandbox):
    rewrite(sandbox, "P5", lambda s: "AAAAAAAAAAAAAAAAAAAA" + s[20:])
    assert blocks(sandbox)


def test_catches_an_encoded_polyA_tract(sandbox):
    """The tail must come from the primer; an A120 here would not synthesise."""
    rewrite(sandbox, "P5", lambda s: s[:-20] + "A" * 40 + s[-20:])
    assert blocks(sandbox)


# --- the paperwork must match the molecule ---------------------------------

def test_catches_a_manifest_that_disagrees_with_the_sequence(sandbox):
    path = sandbox / "orders" / "pjet12_ivt_units" / "element_manifest.csv"
    rows = list(csv.DictReader(path.open()))
    for row in rows:
        if row["construct"] == "P7":
            row["P2A"] = "0"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    assert blocks(sandbox)


def test_catches_a_genbank_that_is_not_the_ordered_sequence(sandbox):
    gb = sandbox / "orders" / "pjet12_ivt_units" / "P3.gb"
    head, body = gb.read_text().split("ORIGIN", 1)
    gb.write_text(head + "ORIGIN" + body.replace("a", "t", 1))
    assert blocks(sandbox)


def test_catches_an_unannotated_module(sandbox):
    """A map that omits an element cannot be used to sanity-check it."""
    gb = sandbox / "orders" / "pjet12_ivt_units" / "P7.gb"
    gb.write_text(gb.read_text().replace('/label="P2A"', '/label="spacer"'))
    assert blocks(sandbox)


def test_catches_the_two_shipped_files_disagreeing(sandbox):
    path = sandbox / "orders" / "pjet12_ivt_units" / "ordering_table.csv"
    rows = list(csv.DictReader(path.open()))
    for row in rows:
        if row["name"] == "P6":
            last = row["sequence"][-1]
            row["sequence"] = row["sequence"][:-1] + ("C" if last != "C" else "A")
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    assert blocks(sandbox)


# --- the requirements added after the first order package ------------------

def test_catches_a_construct_with_a_different_kozak_context(sandbox):
    """The sharpest control here: this mutation is invisible to every other check.

    GCA -> GCC at codon 2 leaves the protein byte-identical, the frame intact,
    every module present and every length unchanged. The only thing it breaks
    is that all nine constructs share one Kozak context -- which is exactly the
    confound the panel cannot afford, and exactly what nothing else would see.
    """
    def repoint(s):
        i = s.index(audit_pjet.LEADER) + len(audit_pjet.LEADER)
        assert s[i:i + 6] == "ATGGCA"
        return s[:i] + "ATGGCC" + s[i + 6:]
    rewrite(sandbox, "P2", repoint)
    assert blocks(sandbox)


def test_the_kozak_mutation_really_is_invisible_to_protein_level_checks(sandbox):
    """Guard against the control above passing for the wrong reason."""
    d = sandbox / "orders" / "pjet12_ivt_units"
    before = audit_pjet.tr(
        read_orf(d / "gblocks.fasta", "P2"))

    def repoint(s):
        i = s.index(audit_pjet.LEADER) + len(audit_pjet.LEADER)
        return s[:i] + "ATGGCC" + s[i + 6:]
    rewrite(sandbox, "P2", repoint)
    after = audit_pjet.tr(read_orf(d / "gblocks.fasta", "P2"))
    assert before == after


def test_catches_a_construct_above_the_gc_ceiling(sandbox, monkeypatch):
    """Drop the ceiling below the panel and the check must fire."""
    monkeypatch.setattr(audit_pjet, "GC_CEILING", 0.30)
    assert blocks(sandbox)


def test_catches_a_panel_whose_gc_spread_is_too_wide(sandbox, monkeypatch, capsys):
    """A spread finding is a warning, so assert it is *reported*, not that it blocks.

    Asserting only "does not block" would pass even if the check had been
    deleted, which is the failure mode this whole file exists to prevent.
    """
    monkeypatch.setattr(audit_pjet, "GC_SPREAD", 0.001)
    audit_pjet.main(sandbox, quiet=False)
    out = capsys.readouterr().out
    assert "wider than the" in out
    assert "WARNING (1)" in out
    assert "DO NOT ORDER" not in out


def test_the_gc_spread_warning_is_absent_at_the_real_threshold(sandbox, capsys):
    audit_pjet.main(sandbox, quiet=False)
    assert "wider than the" not in capsys.readouterr().out


def test_catches_a_lost_kozak_plus_four_G(sandbox):
    def flatten(s):
        i = s.index(audit_pjet.LEADER) + len(audit_pjet.LEADER)
        return s[:i] + "ATGTCA" + s[i + 6:]
    rewrite(sandbox, "P5", flatten)
    assert blocks(sandbox)


def read_orf(path, name: str) -> str:
    """The ORF of one construct, straight out of the fasta."""
    seq, header, buf = None, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header and header.startswith(f">{name}|"):
                seq = "".join(buf)
            header, buf = line, []
        else:
            buf.append(line.strip())
    if header and header.startswith(f">{name}|"):
        seq = "".join(buf)
    start = seq.index(audit_pjet.LEADER) + len(audit_pjet.LEADER)
    return seq[start:]
