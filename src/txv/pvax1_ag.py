"""The pVax1_AG / CTLA-4 routing panel: verified modules and the P0-P8 series.

This is the project-specific library, as distinct from the generic parts in
:mod:`txv.parts`. Every sequence here was taken from the audited P0-P8 build
and re-verified: :func:`panel_proteins` reconstructs all nine ORFs from the
twelve modules below, and ``tests/test_pvax1_ag.py`` locks that against the
reference translations. If a module is edited, those tests fail -- which is the
point, because these modules are the part of the design that is already
settled.

Design logic, in one paragraph
------------------------------
MHC-I is loaded in the ER with peptides the proteasome makes in the cytosol;
MHC-II is loaded in the MIIC, a late endosomal compartment. mRNA is translated
on free cytosolic ribosomes, so a bare antigen ORF is CD8-biased by
construction -- the protein is made correctly and is simply in the wrong
compartment. Redirecting it needs a *matched pair* of signals: a signal peptide
(read by SRP, which only works at the extreme N-terminus) to get the chain into
the ER lumen, and a transmembrane domain to stop it being secreted, which
leaves the cassette luminal and the sorting tail cytosolic. The tail is then
the postcode:

* **CTLA-4 tail (YVKM)** -- a YxxPhi motif read by AP-2, so the protein goes
  ER -> Golgi -> plasma membrane -> clathrin-mediated endocytosis -> MIIC.
  ``YVKM`` is internal, with 19 native residues after it, so this module can sit
  upstream of P2A.
* **LAMP1 tail (GYQTI)** -- read by AP-3, a direct TGN -> lysosome route that
  skips the cell surface. ``GYQTI`` must be the **last residues of the
  protein**; append anything and mu3A binding fails. Hence LAMP1 constructs must
  put this module last, and in P7/P8 the LAMP1 arm must sit *downstream* of P2A.

One interpretive consequence is worth repeating: membrane-routed constructs
still show MHC-I presentation, via ERAD retrotranslocation, failed
translocation and DRiPs. **MHC-I signal is therefore not evidence that routing
worked** -- only the class II readout separates the routes.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .parts import Part, PartKind, PartRegistry, Provenance, default_registry

#: Provenance note attached to every module in this library.
SOURCE = "audited P0-P8 build (pVax1_AG panel); re-verified by reconstruction"

# ---------------------------------------------------------------------------
# Protein modules
# ---------------------------------------------------------------------------

MODULES: dict[str, str] = {
    # -- initiation ---------------------------------------------------------
    "MA": "MA",
    # -- routing: CTLA-4 arm ------------------------------------------------
    "CTLA4_SP": "MACLGLRRYKAQLQLPSRTWPFVALLTLLFIPVFS",
    "CTLA4_TMT": (
        "FLLWILVAVSLGLFFYSFLVSAVSLSKMLKKRSPLTTGVYVKMPPTEPECEKQFQPYFIPIN"
    ),
    # -- routing: LAMP1 arm -------------------------------------------------
    "LAMP1_SP": "MAAPGARRPLLLLLLAGLAHGASA",
    "LAMP1_TMT": "MLIPIAVGGALAGLVLIVLIAYLIGRKRSHAGYQTI",
    # -- antigen cassettes (model epitopes in 29-aa native flanks) ----------
    "A": "EVSGLEQLESIINFEKLTEWTSSNVMEER",
    "B": "EEFAKFASFEAQGALANIAVDKANLDVMK",
    # -- detection tags -----------------------------------------------------
    "HA": "YPYDVPDYA",
    "FLAG": "DYKDDDDK",
    # -- other functional modules -------------------------------------------
    "L": "GGGGS",
    "CL1": "ACKNWFSSLSHFVIHL",
    "P2A": "GSGATNFSLLKQAGDVEENPGP",
}

MODULE_NOTES: dict[str, str] = {
    "MA": "Initiator methionine plus the Kozak alanine, giving the G at +4 that "
          "completes the strong Kozak context.",
    "CTLA4_SP": "CTLA-4 signal peptide, N-terminal by necessity: SRP reads the "
                "first thing out of the ribosome, so an SP works nowhere else. "
                "Includes the MA start.",
    "CTLA4_TMT": "CTLA-4 transmembrane domain plus cytoplasmic tail (residues "
                 "~162-223). Carries the YVKM (YxxPhi) motif read by AP-2: rapid "
                 "clathrin-mediated internalisation from the plasma membrane into "
                 "endosomes and on to the MIIC. The extracellular IgV domain is "
                 "deliberately excluded -- that is the part that binds CD80/CD86 "
                 "and does real checkpoint signalling. YVKM is internal, with 19 "
                 "native residues after it, so this module tolerates downstream "
                 "sequence.",
    "LAMP1_SP": "LAMP1 signal peptide; includes the MA start.",
    "LAMP1_TMT": "LAMP1 transmembrane domain plus tail ending in GYQTI, read by "
                 "AP-3 for a direct TGN -> lysosome route. POSITIONAL: GYQTI must "
                 "be the last residues of the protein. Append anything and mu3A "
                 "binding fails, sending the protein to the plasma membrane.",
    "A": "Antigen A: the MHC-I model epitope SIINFEKL in 29 aa of native OVA "
         "flanking sequence. The flanks are not padding -- ERAP1 trims N-terminal "
         "extensions in the ER but nothing trims a C-terminal one, so the "
         "proteasome must make the C-terminus exactly right first time, and "
         "cleavage specificity is set by the surrounding residues.",
    "B": "Antigen B: the MHC-II model epitope (H2-Ea derived, I-A(b) core) in 29 aa "
         "of native flanking sequence, same rationale as A.",
    "HA": "HA epitope tag for anti-HA detection. Answers 'was it made', on a "
          "different reagent from the pMHC antibodies that answer 'was it "
          "presented'. Without it, an absent pMHC signal has three "
          "indistinguishable causes. It is not an assembly overlap, homology arm, "
          "restriction site or selection marker.",
    "FLAG": "FLAG tag on the second cistron of the dual-route constructs, so the "
            "two P2A-separated products can be detected independently.",
    "L": "GGGGS flexible linker.",
    "CL1": "CL1 degron: a hydrophobic degradation signal that drives proteasomal "
           "turnover, boosting the cytosolic (MHC-I) arm.",
    "P2A": "P2A ribosomal skipping element with its GSG spacer. Yields two "
           "separate polypeptides from one ORF, which is what makes a dual-route "
           "construct possible. Skipping is efficient but not complete, so expect "
           "some read-through fusion product.",
}

# ---------------------------------------------------------------------------
# pVax1_AG backbone rails
# ---------------------------------------------------------------------------
# Only the parts of the backbone that were explicitly recorded are reproduced
# here: the primer-annealing rails and the assembly overlaps. The full 5'/3' UTR
# sequences are NOT in this library -- take them from the pvax1_ag_egfp.gb
# record itself and load them with txv.parts.load_parts_file.

BACKBONE_OLIGOS: dict[str, tuple[str, str]] = {
    "backbone_F": (
        "GCTGCCTTCTGCGGGGCTTG",
        "PCR linearisation forward; starts at the retained 3' UTR. Equals the "
        "insert right overlap.",
    ),
    "backbone_R": (
        "GGTGGCTCTTATATTTCTTCTTACT",
        "PCR linearisation reverse; ends at the retained 5' UTR boundary. Its "
        "reverse complement is AGTAAGAAGAAATATAAGAGCCACC, i.e. it anneals across "
        "the 3' end of the 5' UTR *and* the Kozak -- so its last 20 nt are "
        "exactly the insert left overlap.",
    ),
    "insert_left_overlap": (
        "GAAGAAATATAAGAGCCACC",
        "Prepend to each ORF. The final GCCACC is the Kozak; the preceding "
        "GAAGAAATATAAGAG is the 3' end of the pVax1_AG 5' UTR.",
    ),
    "insert_right_overlap": (
        "GCTGCCTTCTGCGGGGCTTG",
        "Append after each ORF stop; the 5' end of the pVax1_AG 3' UTR.",
    ),
}

#: Open decisions carried over from the audit. These are design choices, not
#: defects -- they are recorded so they do not get silently resolved by default.
OPEN_DECISIONS: tuple[str, ...] = (
    "Encoded A120 vs PCR-added A120: the plasmid records carry an encoded A120 "
    "tract, but the IVT reverse primer anneals upstream and adds A120 by PCR. "
    "Removing the encoded tract would improve synthesis and propagation "
    "stability (long homopolymers contract in E. coli). Unresolved.",
    "P1 and P2 are CTLA-4-scaffold single-antigen controls. They are NOT "
    "route-matched controls for P5-P8; do not use them as such.",
    "A HindIII site (AAGCTT) was found at ORF position 251 in the P1 and P2 "
    "builds. Harmless unless HindIII is used for cloning or linearisation -- "
    "re-optimising those ORFs with txv removes it.",
    "Generic commercial pVAX1 is not equivalent to the custom pVax1_AG "
    "sequence. Preserve the exact T7 start context, UTR rails, primer sites, "
    "origin and kanamycin/neomycin cassette.",
)


# ---------------------------------------------------------------------------
# The panel
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PanelConstruct:
    """One member of the P0-P8 series, defined as an ordered module list."""

    name: str
    modules: tuple[str, ...]
    route: str
    question: str
    antigens: tuple[str, ...] = ()

    def protein(self) -> str:
        return "".join(MODULES[m] for m in self.modules)

    def __len__(self) -> int:
        return len(self.protein())


PANEL: tuple[PanelConstruct, ...] = (
    PanelConstruct(
        "P0", ("CTLA4_SP", "HA", "L", "CTLA4_TMT"), "CTLA-4 (AP-2, via surface)",
        "Scaffold-only control: does the routing module alone produce signal?",
        (),
    ),
    PanelConstruct(
        "P1", ("CTLA4_SP", "HA", "L", "A", "L", "CTLA4_TMT"),
        "CTLA-4 (AP-2, via surface)",
        "Single-antigen (A) scaffold control.", ("A",),
    ),
    PanelConstruct(
        "P2", ("CTLA4_SP", "HA", "L", "B", "L", "CTLA4_TMT"),
        "CTLA-4 (AP-2, via surface)",
        "Single-antigen (B) scaffold control.", ("B",),
    ),
    PanelConstruct(
        "P3", ("CTLA4_SP", "HA", "L", "A", "L", "B", "L", "CTLA4_TMT"),
        "CTLA-4 (AP-2, via surface)",
        "Can one transcript drive simultaneous class I and class II "
        "presentation? Order A-then-B.", ("A", "B"),
    ),
    PanelConstruct(
        "P4", ("CTLA4_SP", "HA", "L", "B", "L", "A", "L", "CTLA4_TMT"),
        "CTLA-4 (AP-2, via surface)",
        "Order control for P3: does cassette order change the outcome?",
        ("B", "A"),
    ),
    PanelConstruct(
        "P5", ("MA", "HA", "L", "A", "L", "B", "L", "CL1"), "cytosolic (+CL1 degron)",
        "Cytosolic baseline: proteasome-directed, no membrane routing.",
        ("A", "B"),
    ),
    PanelConstruct(
        "P6", ("LAMP1_SP", "HA", "L", "A", "L", "B", "L", "LAMP1_TMT"),
        "LAMP1 (AP-3, direct to lysosome)",
        "Does the direct TGN -> lysosome route beat the detour via the cell "
        "surface that P0-P4 take?", ("A", "B"),
    ),
    PanelConstruct(
        "P7",
        ("MA", "HA", "L", "A", "L", "CL1", "P2A", "LAMP1_SP", "FLAG", "L", "B",
         "L", "LAMP1_TMT"),
        "dual: A cytosolic + B lysosomal",
        "Dual route from one transcript: A to the proteasome, B to the "
        "endolysosome. LAMP1 arm is downstream of P2A because GYQTI must end "
        "the protein.", ("A", "B"),
    ),
    PanelConstruct(
        "P8",
        ("MA", "HA", "L", "B", "L", "CL1", "P2A", "LAMP1_SP", "FLAG", "L", "A",
         "L", "LAMP1_TMT"),
        "dual: B cytosolic + A lysosomal",
        "Swap control for P7: is the effect about the route or about the "
        "antigen?", ("B", "A"),
    ),
)

PANEL_BY_NAME: dict[str, PanelConstruct] = {p.name: p for p in PANEL}


def panel_proteins() -> dict[str, str]:
    """Reconstruct every panel ORF from the module library."""
    return {p.name: p.protein() for p in PANEL}


def module_registry(base: PartRegistry | None = None) -> PartRegistry:
    """A :class:`~txv.parts.PartRegistry` carrying these modules as parts.

    Modules are registered under a ``pvax1_`` prefix so they never collide with
    the generic built-ins.
    """
    registry = base or default_registry()
    kinds = {
        "CTLA4_SP": PartKind.SIGNAL_PEPTIDE,
        "LAMP1_SP": PartKind.SIGNAL_PEPTIDE,
        "CTLA4_TMT": PartKind.TRAFFICKING,
        "LAMP1_TMT": PartKind.TRAFFICKING,
        "HA": PartKind.TAG,
        "FLAG": PartKind.TAG,
        "L": PartKind.LINKER,
        "P2A": PartKind.LINKER,
        "CL1": PartKind.TAG,
        "MA": PartKind.SIGNAL_PEPTIDE,
        "A": PartKind.TAG,
        "B": PartKind.TAG,
    }
    for name, protein in MODULES.items():
        registry.add(
            Part(
                name=f"pvax1_{name}",
                kind=kinds.get(name, PartKind.LINKER),
                provenance=Provenance.LITERATURE,
                protein=protein,
                note=MODULE_NOTES.get(name, ""),
                source=SOURCE,
            )
        )
    for name, (dna, note) in BACKBONE_OLIGOS.items():
        registry.add(
            Part(
                name=f"pvax1_{name}",
                kind=PartKind.RESTRICTION,
                provenance=Provenance.LITERATURE,
                dna=dna,
                note=note,
                source=SOURCE,
            )
        )
    return registry


# ---------------------------------------------------------------------------
# Building the panel as IVT mRNA constructs
# ---------------------------------------------------------------------------

#: Module -> annotation kind, so each module gets its own coloured feature in
#: the GenBank record and in Benchling rather than one opaque CDS block.
MODULE_KIND: dict[str, str] = {
    "MA": "start",
    "CTLA4_SP": "signal_peptide",
    "LAMP1_SP": "signal_peptide",
    "CTLA4_TMT": "trafficking",
    "LAMP1_TMT": "trafficking",
    "A": "neoepitope",
    "B": "neoepitope",
    "HA": "tag",
    "FLAG": "tag",
    "CL1": "degron",
    "P2A": "linker",
    "L": "linker",
}


def panel_cassette(name: str) -> "Cassette":
    """The panel construct as a :class:`~txv.epitopes.Cassette` of modules.

    Each module becomes its own bead with an empty linker between beads, so the
    assembled construct carries one annotated feature per module. Module order
    is load-bearing and must not be reordered: the signal peptide only works at
    the extreme N-terminus, and ``GYQTI`` only works at the extreme C-terminus.
    """
    from .epitopes import Antigen, Cassette

    panel = PANEL_BY_NAME[name]
    beads = []
    for index, module in enumerate(panel.modules):
        # Module names repeat within a construct (several GGGGS linkers, two
        # tags); suffix them so features stay distinguishable.
        label = f"{module}_{index}" if panel.modules.count(module) > 1 else module
        beads.append(
            Antigen(
                name=label,
                sequence=MODULES[module],
                kind=MODULE_KIND.get(module, "full_length"),
                note=MODULE_NOTES.get(module, ""),
            )
        )
    return Cassette(beads, linker="")


#: Routes available to a generic cassette, built from the verified modules.
#: Each entry is (signal peptide part, trafficking part, note).
ROUTES: dict[str, tuple[str | None, str | None, str]] = {
    "ctla4": ("pvax1_CTLA4_SP", "pvax1_CTLA4_TMT",
              "AP-2 / YVKM: ER -> Golgi -> surface -> endocytosis -> MIIC."),
    "lamp1": ("pvax1_LAMP1_SP", "pvax1_LAMP1_TMT",
              "AP-3 / GYQTI: direct TGN -> lysosome, skipping the surface."),
    "cytosolic": (None, None,
                  "No routing: translated on free ribosomes, proteasome-directed. "
                  "CD8-biased by construction."),
}


def routed_spec(
    route: str = "ctla4",
    name: str = "construct",
    **spec_kwargs,
) -> tuple["ConstructSpec", PartRegistry]:
    """A :class:`~txv.constructs.ConstructSpec` wired to a verified route.

    Returns ``(spec, registry)``. Use this to put a *generic* antigen cassette
    on one of the panel's routes -- the modules are real sequences, so the
    resulting construct contains no placeholders.

    Note the positional constraint the LAMP1 route carries: ``GYQTI`` must end
    the protein, so the spec puts no spacer after the trafficking domain.
    """
    from .constructs import ConstructSpec

    try:
        signal_peptide, trafficking, _ = ROUTES[route]
    except KeyError:
        raise KeyError(
            f"unknown route {route!r}; choose from {sorted(ROUTES)}"
        ) from None

    spec = ConstructSpec(
        name=name,
        signal_peptide=signal_peptide,
        trafficking=trafficking,
        # GYQTI must be the last residues of the protein; a spacer after the
        # trafficking domain would abolish AP-3 binding.
        traffic_spacer=None,
        **spec_kwargs,
    )
    return spec, module_registry()


def build_panel_construct(
    name: str,
    spec: "ConstructSpec | None" = None,
    registry: "PartRegistry | None" = None,
    optimizer_config=None,
):
    """Assemble one panel member as an IVT mRNA template.

    The modules already supply their own initiation, signal peptide and
    trafficking domain, so the generic spec elements for those are switched
    off: adding another signal peptide or a C-terminal domain would break both
    positional constraints at once.
    """
    from .codon import CodonOptimizer
    from .constructs import ConstructBuilder, ConstructSpec

    panel = PANEL_BY_NAME[name]
    spec = spec or ConstructSpec(name=name)
    spec.name = name
    spec.signal_peptide = None
    spec.trafficking = None
    spec.sp_spacer = None
    spec.traffic_spacer = None

    builder = ConstructBuilder(
        registry=registry or module_registry(),
        optimizer=CodonOptimizer(optimizer_config) if optimizer_config else None,
    )
    construct = builder.build(spec, panel_cassette(name))
    construct.design_notes.append(f"{panel.name} route: {panel.route}")
    construct.design_notes.append(f"{panel.name} asks: {panel.question}")
    return construct


def build_panel(
    names: "tuple[str, ...] | None" = None,
    spec_factory=None,
    optimizer_config=None,
) -> dict:
    """Assemble the whole panel. Returns ``{name: Construct}``."""
    from .constructs import ConstructSpec

    selected = names or tuple(p.name for p in PANEL)
    out = {}
    for name in selected:
        spec = spec_factory(name) if spec_factory else ConstructSpec(name=name)
        out[name] = build_panel_construct(
            name, spec=spec, optimizer_config=optimizer_config
        )
    return out


__all__ = [
    "MODULES", "MODULE_NOTES", "MODULE_KIND", "BACKBONE_OLIGOS",
    "OPEN_DECISIONS", "PANEL", "PANEL_BY_NAME", "PanelConstruct", "SOURCE",
    "ROUTES", "routed_spec",
    "panel_proteins", "module_registry", "panel_cassette",
    "build_panel_construct", "build_panel",
]
