"""MHC class II routing panel: one antigen, seven destinations.

A separate study from the P0-P8 panel and deliberately narrower. There the
variable was both antigen and route; here the antigen is fixed to B (the
Ea-bearing cassette, read out by Y-Ae) and the only thing that moves is how
the protein is routed. That makes every pairwise comparison a routing
comparison, which the nine-construct panel cannot claim.

Shared modules are taken from :mod:`txv.pvax1_ag` and keep their canonical
encodings, so a construct here is directly comparable with its counterpart
there rather than merely similar.

Three of the seven need sequence this package does not hold. They are declared
with ``sequence=None`` and the build refuses to emit them rather than guessing:
ordering a gBlock built from a recalled sequence is the one failure that costs
real money and months.
"""

from __future__ import annotations

from dataclasses import dataclass

from .pvax1_ag import MODULES as SHARED_MODULES


@dataclass(frozen=True)
class ExternalModule:
    """A module whose sequence has to come from a database, not from memory."""

    name: str
    protein: str          # what to look up
    organism: str
    region: str           # which residues, and why those boundaries
    why: str              # what the module is for
    sequence: str | None = None

    @property
    def resolved(self) -> bool:
        return self.sequence is not None


#: The three modules this panel needs and does not have. UniProt, NCBI and
#: Ensembl are all unreachable from the build environment, so these are left
#: empty on purpose. Fill `sequence` from the named record and the build
#: completes; leave it and the build skips that construct and says so.
EXTERNAL: dict[str, ExternalModule] = {
    "H2DMB_SP": ExternalModule(
        "H2DMB_SP", "H2-DMb1 (H2-M beta chain)", "Mus musculus",
        "signal peptide only, residues 1 to the signal peptidase site",
        "N-terminal entry into the secretory pathway, matched to the "
        "H2-DM beta tail below so the pair is one protein's own.",
    ),
    "H2DMB_TMT": ExternalModule(
        "H2DMB_TMT", "H2-DMb1 (H2-M beta chain)", "Mus musculus",
        "transmembrane domain through the native C-terminus, which must "
        "retain the tyrosine-based cytoplasmic sorting motif intact and "
        "terminal",
        "The MIIC-targeting determinant. This is the module the construct "
        "exists to test, so a truncation that clips the motif silently "
        "produces a negative result that means nothing.",
    ),
    "II_CLIP_SCAFFOLD": ExternalModule(
        "II_CLIP_SCAFFOLD", "CD74 / invariant chain (Ii p31)", "Mus musculus",
        "full ORF, with the CLIP segment excised and its boundaries recorded "
        "so the epitope can be dropped into that register",
        "Type II membrane protein: the N-terminal cytoplasmic tail carries "
        "the di-leucine endosomal targeting motifs and must not be modified "
        "or tagged.",
    ),
    "LC3B": ExternalModule(
        "LC3B", "MAP1LC3B (LC3B)", "Mus musculus",
        "residues 1-120, ending at the glycine that ATG4 exposes; nothing "
        "may follow it",
        "Autophagosome delivery. The C-terminal glycine must stay free for "
        "lipidation, so this module is terminal by construction and cannot "
        "carry a downstream tag or stop-adjacent element.",
    ),
}


@dataclass(frozen=True)
class RoutingConstruct:
    name: str
    modules: tuple[str, ...]
    route: str
    hypothesis: str
    tag_position: str
    notes: str = ""

    @property
    def external(self) -> tuple[str, ...]:
        return tuple(m for m in self.modules if m in EXTERNAL)

    @property
    def buildable(self) -> bool:
        return all(EXTERNAL[m].resolved for m in self.external)


#: `B` is the Ea-bearing cassette: EEFAKFASFEAQGALANIAVDKANLDVMK, which
#: embeds Ea52-68 (ASFEAQGALANIAVDKA) with six native residues either side so
#: there is something for a protease to cut. Every construct carries the same
#: cassette; only the route changes.
PANEL: tuple[RoutingConstruct, ...] = (
    RoutingConstruct(
        "R0", ("MA", "HA", "L", "B", "L", "E5"),
        "cytosolic", 
        "Baseline. Without it no construct below can be said to have improved "
        "anything, because there is no zero.",
        "N-terminal, cytosolic",
        "Not in the original five. Added because every other construct is a "
        "comparison and this is the thing they are compared to.",
    ),
    RoutingConstruct(
        "R1", ("CTLA4_SP", "HA", "L", "B", "L", "CTLA4_TMT", "E5"),
        "CTLA-4: secretory, to the surface, then AP-2 to the endosome",
        "Benchmark. Identical to P2 of the existing panel, so this study "
        "joins onto that one instead of starting over.",
        "N-terminal, lumenal after signal peptidase",
    ),
    RoutingConstruct(
        "R2", ("LAMP1_SP", "HA", "L", "B", "L", "LAMP1_TMT"),
        "LAMP1: AP-3, direct TGN to lysosome",
        "A direct route to the lysosome presents Ea better than the detour "
        "via the plasma membrane that R1 takes.",
        "N-terminal, lumenal after signal peptidase",
        "No E5: GYQTI must be the last five residues, so nothing may follow "
        "the tail.",
    ),
    RoutingConstruct(
        "R3", ("H2DMB_SP", "HA", "L", "B", "L", "H2DMB_TMT"),
        "H2-M beta tail: to the MHC II loading compartment",
        "Delivery to the compartment where H2-DM catalyses CLIP release "
        "beats delivery to the lysosome generally.",
        "N-terminal, lumenal after signal peptidase",
        "A chimera, not native H2-M. H2-DM is an alpha/beta heterodimer; a "
        "lone beta-derived tail carries the sorting motif but not the "
        "partner, so this tests the targeting determinant, not the enzyme.",
    ),
    RoutingConstruct(
        "R4", ("II_CLIP_SCAFFOLD",),
        "invariant chain, epitope in the CLIP register",
        "Using the class II delivery system itself, rather than arriving in "
        "the same compartment and competing for loading, raises productive "
        "Ea presentation.",
        "C-terminal, lumenal",
        "Type II membrane protein. The N-terminus is CYTOSOLIC and carries "
        "the di-leucine targeting motifs, so it takes no signal peptide and "
        "no N-terminal tag. The epitope replaces CLIP in-register; the B "
        "cassette's flanks are wrong for this one and the Ea core is used "
        "instead.",
    ),
    RoutingConstruct(
        "R5", ("MA", "HA", "L", "B", "L", "LC3B"),
        "macroautophagy: antigen into the autophagosome",
        "Autophagosomes feed the class II loading compartment, giving a "
        "productive route that does not use the secretory pathway at all.",
        "N-terminal, cytosolic",
        "Antigen must be N-terminal to LC3B and nothing may follow LC3B's "
        "C-terminal glycine, so this construct cannot carry E5 or a "
        "C-terminal tag.",
    ),
    RoutingConstruct(
        "R6", ("CTLA4_SP", "HA", "L", "B"),
        "secreted: signal peptide, no membrane anchor",
        "Separates entering the secretory pathway from arriving anywhere in "
        "particular. If R1 and R2 beat R0 but R6 does too, the gain is ER "
        "entry rather than the destination.",
        "N-terminal, lumenal after signal peptidase",
        "Not in the original five. Added because without it, 'routing "
        "works' and 'the secretory pathway works' are the same result.",
    ),
)

PANEL_BY_NAME = {c.name: c for c in PANEL}


def module_protein(name: str) -> str:
    if name in SHARED_MODULES:
        return SHARED_MODULES[name]
    ext = EXTERNAL[name]
    if ext.sequence is None:
        raise KeyError(
            f"{name} has no sequence. It is {ext.protein} ({ext.organism}), "
            f"{ext.region}. Fill EXTERNAL[{name!r}].sequence from that record."
        )
    return ext.sequence


def protein(name: str) -> str:
    """Assembled protein for one construct, or raise naming what is missing."""
    return "".join(module_protein(m) for m in PANEL_BY_NAME[name].modules)
