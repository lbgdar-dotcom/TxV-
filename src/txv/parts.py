"""Registry of construct parts.

Design rule for this module: **no invented sequences.** Every part carries a
:class:`Provenance` tag saying how much you should trust the literal sequence
stored here.

``CANONICAL``
    Short, universally agreed sequences (T7 promoter, Kozak consensus,
    homopolymeric poly(A)). Safe to use as-is.
``LITERATURE``
    Widely published peptide elements (signal peptides, polyepitope linkers).
    Verify against your own reference before a GMP run, but they are standard.
``PLACEHOLDER``
    Elements where the exact sequence is programme-specific or where you should
    paste in the validated sequence your group actually uses -- notably the
    5'/3' UTRs, which dominate expression and are the part most worth owning.
    A construct that still contains a placeholder fails QC unless you pass
    ``allow_placeholders=True``.

Load your own parts with :func:`load_parts_file` (JSON) and they override the
built-ins by name.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from enum import Enum
from pathlib import Path

from .seqops import CODON_TO_AA, clean


class Provenance(str, Enum):
    CANONICAL = "canonical"
    LITERATURE = "literature"
    PLACEHOLDER = "placeholder"


class PartKind(str, Enum):
    PROMOTER = "promoter"
    UTR5 = "utr5"
    KOZAK = "kozak"
    SIGNAL_PEPTIDE = "signal_peptide"
    TRAFFICKING = "trafficking"
    LINKER = "linker"
    TAG = "tag"
    UTR3 = "utr3"
    POLYA = "polya"
    RESTRICTION = "restriction"


@dataclass(frozen=True)
class Part:
    """A construct element.

    Exactly one of ``dna`` / ``protein`` is the authoritative form. Protein
    parts (signal peptides, linkers) are codon-optimised at assembly time in
    the same pass as the antigen payload, so they get the same uridine and
    motif treatment as the rest of the ORF.
    """

    name: str
    kind: PartKind
    provenance: Provenance
    dna: str | None = None
    protein: str | None = None
    note: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        if (self.dna is None) == (self.protein is None):
            raise ValueError(f"part {self.name!r}: set exactly one of dna/protein")
        if self.dna is not None:
            object.__setattr__(self, "dna", clean(self.dna))
        if self.protein is not None:
            object.__setattr__(self, "protein", self.protein.strip().upper())
            bad = set(self.protein) - set("ACDEFGHIKLMNPQRSTVWY*")
            if bad:
                raise ValueError(f"part {self.name!r}: bad residues {sorted(bad)}")

    @property
    def is_placeholder(self) -> bool:
        return self.provenance is Provenance.PLACEHOLDER


_BUILTIN: list[Part] = [
    # ---- transcription start -------------------------------------------------
    Part(
        "T7_promoter", PartKind.PROMOTER, Provenance.CANONICAL,
        dna="TAATACGACTCACTATAG",
        note="Class III T7 phi6.5 promoter. The final G is +1 of the transcript; "
             "T7 strongly prefers G at +1 and tolerates GG/GGG better than A.",
    ),
    Part(
        "SP6_promoter", PartKind.PROMOTER, Provenance.CANONICAL,
        dna="ATTTAGGTGACACTATAG",
        note="Alternative to T7 if your RNAP or capping chemistry calls for it.",
    ),
    # ---- untranslated regions ------------------------------------------------
    Part(
        "UTR5_hAg", PartKind.UTR5, Provenance.LITERATURE,
        dna="GAGAATAAACTAGTATTCTTCTGGTCCCCACAGACTCAGAGAGAACCC",
        source="BNT162b2 published sequence (5' leader, nt 1-48)",
        note="DEFAULT 5' UTR. The BNT162b2 5' leader up to but not including the "
             "Kozak, which kozak_strong (GCCACC) supplies -- the two concatenated "
             "reproduce the validated 54-nt leader exactly. Human alpha-globin "
             "(HBA1/HBA2-shared) 5' UTR preceded by a 14-nt vector leader "
             "(GAGAATAAACTAGT). No upstream AUG. Two things to know about the "
             "leader: it contains AATAAA at position 3 and an SpeI site (ACTAGT) "
             "at position 8. Both are in the clinically validated sequence and are "
             "inert for a cytoplasmically delivered transcript -- but if you cut "
             "with SpeI, use UTR5_hAg_core instead.",
    ),
    Part(
        "UTR5_hAg_core", PartKind.UTR5, Provenance.LITERATURE,
        dna="GGGATTCTTCTGGTCCCCACAGACTCAGAGAGAACCC",
        source="human alpha-globin 5' UTR from the BNT162b2 leader (nt 15-48), "
               "with a GGG transcription start",
        note="The alpha-globin 5' UTR without the BNT162b2 vector leader, so no "
             "AATAAA and no SpeI site, prefixed with GGG because T7 initiates far "
             "more efficiently from a G-rich +1. Use when you want the validated "
             "UTR but not the cloning scars.",
    ),
    Part(
        "UTR5_placeholder", PartKind.UTR5, Provenance.PLACEHOLDER,
        dna="AGGAAATAAGAGAGAAAAGAAGAGTAAGAAG",
        note="NOT a validated 5' UTR -- an unstructured, A-rich, ATG-free stand-in "
             "so the assembler runs end to end. Replace with your validated 5' UTR "
             "(e.g. an HBA/HBB-derived or screened synthetic leader). Requirements: "
             "no upstream AUG, weak 5'-proximal secondary structure, ~30-60 nt.",
    ),
    Part(
        "UTR3_AES_mtRNR1", PartKind.UTR3, Provenance.LITERATURE,
        dna="CTCGAGCTGGTACTGCATGCACGCAATGCTAGCTGCCCCTTTCCCGTCCTGGGTACCCCGAG"
            "TCTCCCCCGACCTCGGGTCCCAGGTATGCTCCCACCTCCACCTGCCCCACTCACCACCTCTG"
            "CTAGTTCCAGACACCTCCCAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCCACAC"
            "CCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACT"
            "AACCCCAGGGTTGGTCAATTTCGTGCCAGCCACACCCTGGAGCTAGCA",
        source="BNT162b2 published sequence (3' UTR, 296 nt)",
        note="DEFAULT 3' UTR. A composite of two elements selected ex vivo for "
             "RNA stability and total protein output: the amino-terminal enhancer "
             "of split (AES) 3' UTR and the mitochondrially encoded 12S rRNA "
             "(mtRNR1). Verbatim from the validated construct, which means it "
             "carries its cloning scars: XhoI (CTCGAG) at position 0 and NheI "
             "(GCTAGC) at 27 and 289. It also contains AATAAA at 219, inside the "
             "mtRNR1 element. All three are in the clinical sequence; they matter "
             "only if you use those enzymes.",
    ),
    Part(
        "UTR3_placeholder", PartKind.UTR3, Provenance.PLACEHOLDER,
        dna="GCTGCCTTCTGCGGGGCTTGCCTTCTGGCCATGCCCTTCTTCTCTCCCTTGCACCTGT",
        note="NOT a validated 3' UTR -- stand-in only. Replace with your validated "
             "3' UTR. Stability-optimised designs in the field are typically "
             "globin-derived or a screened two-element composite; whatever you use, "
             "keep it free of AU-rich destabilising elements and miR seed matches "
             "for the target tissue.",
    ),
    Part(
        "kozak_strong", PartKind.KOZAK, Provenance.CANONICAL,
        dna="GCCACC",
        note="Placed immediately 5' of the initiator ATG; with the G at +4 this is "
             "the strong Kozak context (gccRccATGG).",
    ),
    # ---- targeting elements (protein-level) ---------------------------------
    Part(
        "SP_tPA", PartKind.SIGNAL_PEPTIDE, Provenance.LITERATURE,
        protein="MDAMKRGLCCVLLLCGAVFVSPS",
        note="Human tissue plasminogen activator signal peptide. Routes the "
             "polypeptide into the secretory pathway; a standard choice for "
             "improving antigen presentation in genetic vaccines.",
    ),
    Part(
        "SP_IgE", PartKind.SIGNAL_PEPTIDE, Provenance.LITERATURE,
        protein="MDWTWILFLVAAATRVHS",
        note="Human IgE heavy-chain leader. Common alternative to tPA.",
    ),
    Part(
        "MITD_placeholder", PartKind.TRAFFICKING, Provenance.PLACEHOLDER,
        protein="GGGGSGGGGS",
        note="Stand-in for an MHC class I trafficking domain (the transmembrane + "
             "cytoplasmic tail of an HLA class I heavy chain) fused C-terminally to "
             "route the payload through the MHC I presentation compartment. Paste in "
             "the exact HLA-derived tail your programme uses.",
    ),
    Part(
        "LAMP1_sorting_motif", PartKind.TRAFFICKING, Provenance.LITERATURE,
        protein="RKRSHAGYQTI",
        note="LAMP1 C-terminal region carrying the YQTI lysosomal sorting signal; "
             "used to divert antigen to the MHC class II loading compartment. Supply "
             "the full LAMP1 TM+tail if you need the complete element.",
    ),
    # ---- polyepitope linkers -------------------------------------------------
    Part(
        "linker_GS10", PartKind.LINKER, Provenance.CANONICAL,
        protein="GGSGGGGSGG",
        note="Flexible glycine-serine spacer. Neutral default between antigen "
             "cassettes; does not itself promote processing.",
    ),
    Part(
        "linker_GS5", PartKind.LINKER, Provenance.CANONICAL,
        protein="GGGGS",
        note="Short flexible spacer for tight packing of long cassettes.",
    ),
    Part(
        "linker_AAY", PartKind.LINKER, Provenance.LITERATURE,
        protein="AAY",
        note="Favours immunoproteasomal cleavage C-terminal to the epitope; a "
             "standard choice between MHC class I (CD8) minimal epitopes.",
    ),
    Part(
        "linker_GPGPG", PartKind.LINKER, Provenance.LITERATURE,
        protein="GPGPG",
        note="Used between MHC class II (CD4/helper) epitopes; disfavours "
             "junctional class II binding.",
    ),
    Part(
        "linker_KK", PartKind.LINKER, Provenance.LITERATURE,
        protein="KK",
        note="Cathepsin-B-favouring dibasic linker, class II contexts.",
    ),
    Part(
        "linker_furin_RAKR", PartKind.LINKER, Provenance.LITERATURE,
        protein="RAKR",
        note="Furin cleavage site: separates cassettes in the secretory pathway so "
             "each is released as a discrete polypeptide.",
    ),
    Part(
        "linker_T2A", PartKind.LINKER, Provenance.LITERATURE,
        protein="GSGEGRGSLLTCGDVEENPGP",
        note="T2A ribosome-skipping element (with the upstream GSG spacer). Use when "
             "cassettes must be genuinely separate proteins, not one polypeptide. "
             "Leaves residual C-terminal residues on the upstream cassette.",
    ),
    # ---- tails ---------------------------------------------------------------
    Part(
        "polyA_120", PartKind.POLYA, Provenance.CANONICAL,
        dna="A" * 120,
        note="DEFAULT tail: homopolymeric 120-nt poly(A), encoded in the template "
             "so every transcript has a defined tail length. Long homopolymers can "
             "contract during propagation in E. coli, so verify the tract length "
             "on each prep rather than assuming it (see segmented_polya for the "
             "alternative that resists contraction).",
    ),
    Part(
        "polyA_100", PartKind.POLYA, Provenance.CANONICAL,
        dna="A" * 100,
        note="Homopolymeric 100-nt poly(A).",
    ),
    # ---- template linearisation ---------------------------------------------
    Part(
        "BspQI_site", PartKind.RESTRICTION, Provenance.CANONICAL,
        dna="GCTCTTC",
        note="Type IIS site for run-off linearisation immediately after the poly(A) "
             "so the transcript ends in A with no vector-derived overhang.",
    ),
    Part(
        "BsaI_site", PartKind.RESTRICTION, Provenance.CANONICAL,
        dna="GGTCTC",
        note="Type IIS alternative for linearisation or Golden Gate assembly.",
    ),
]


def segmented_polya(first: int = 30, linker: str = "GCATATGACT", second: int = 70) -> Part:
    """Segmented poly(A) tail: A(n) - short linker - A(m).

    Splitting the tail keeps it from recombining or contracting during plasmid
    propagation in *E. coli*, which is the usual reason a long homopolymeric
    tail drifts between preps. ``linker`` is programme-specific -- pass the
    spacer your group has validated.
    """
    return Part(
        f"polyA_A{first}L{len(linker)}A{second}",
        PartKind.POLYA,
        Provenance.PLACEHOLDER,
        dna="A" * first + linker + "A" * second,
        note=f"Segmented tail A{first}-linker({len(linker)} nt)-A{second}. The linker "
             "sequence here is a generic spacer: substitute your validated one.",
    )


class PartRegistry:
    """Name -> :class:`Part` lookup with user overrides."""

    def __init__(self, parts: list[Part] | None = None) -> None:
        self._parts: dict[str, Part] = {}
        for part in parts if parts is not None else _BUILTIN:
            self.add(part)

    def add(self, part: Part, overwrite: bool = True) -> Part:
        if part.name in self._parts and not overwrite:
            raise KeyError(f"part {part.name!r} already registered")
        self._parts[part.name] = part
        return part

    def get(self, name: str) -> Part:
        try:
            return self._parts[name]
        except KeyError:
            raise KeyError(
                f"unknown part {name!r}; known: {', '.join(sorted(self._parts))}"
            ) from None

    def of_kind(self, kind: PartKind) -> list[Part]:
        return [p for p in self._parts.values() if p.kind is kind]

    def placeholders(self) -> list[Part]:
        return [p for p in self._parts.values() if p.is_placeholder]

    def names(self) -> list[str]:
        return sorted(self._parts)

    def __contains__(self, name: object) -> bool:
        return name in self._parts

    def __len__(self) -> int:
        return len(self._parts)


def load_parts_file(path: str | Path, registry: PartRegistry | None = None) -> PartRegistry:
    """Merge a JSON parts file into a registry (user parts win).

    Expected shape::

        {"parts": [{"name": "UTR5_house", "kind": "utr5",
                    "provenance": "literature", "dna": "GGGA...",
                    "source": "internal ref 2024-11", "note": "..."}]}
    """
    registry = registry or PartRegistry()
    payload = json.loads(Path(path).read_text())
    for entry in payload.get("parts", []):
        registry.add(
            Part(
                name=entry["name"],
                kind=PartKind(entry["kind"]),
                provenance=Provenance(entry.get("provenance", "literature")),
                dna=entry.get("dna"),
                protein=entry.get("protein"),
                note=entry.get("note", ""),
                source=entry.get("source", ""),
            )
        )
    return registry


def default_registry() -> PartRegistry:
    return PartRegistry()


def replace_part(part: Part, **changes) -> Part:
    return replace(part, **changes)


__all__ = [
    "Part", "PartKind", "PartRegistry", "Provenance",
    "default_registry", "load_parts_file", "segmented_polya", "replace_part",
]
