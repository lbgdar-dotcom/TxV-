# MHC II routing panel -- ORFs for Benchling

Start codon to stop codon. No UTRs, no promoter, no vector.
Codon usage: mouse. Shared modules carry the same canonical
encodings as the P0-P8 panel, so the two studies are directly
comparable.

| construct | bp | aa | route | tag |
|---|---|---|---|---|
| R0 | 168 | 55 | cytosolic | N-terminal, cytosolic |
| R1 | 453 | 150 | CTLA-4: secretory, to the surface, then AP-2 to the endosome | N-terminal, lumenal after signal peptidase |
| R2 | 327 | 108 | LAMP1: AP-3, direct TGN to lysosome | N-terminal, lumenal after signal peptidase |
| R6 | 237 | 78 | secreted: signal peptide, no membrane anchor | N-terminal, lumenal after signal peptidase |

## Not built

These need a sequence this package does not hold; it was
left empty rather than recalled, because a gBlock ordered
from a remembered sequence is an expensive way to find out
it was wrong.

**R3** -- H2-M beta tail: to the MHC II loading compartment
- `H2DMB_SP`: H2-DMb1 (H2-M beta chain) (Mus musculus)
  - region: signal peptide only, residues 1 to the signal peptidase site
  - why: N-terminal entry into the secretory pathway, matched to the H2-DM beta tail below so the pair is one protein's own.
- `H2DMB_TMT`: H2-DMb1 (H2-M beta chain) (Mus musculus)
  - region: transmembrane domain through the native C-terminus, which must retain the tyrosine-based cytoplasmic sorting motif intact and terminal
  - why: The MIIC-targeting determinant. This is the module the construct exists to test, so a truncation that clips the motif silently produces a negative result that means nothing.
- note: A chimera, not native H2-M. H2-DM is an alpha/beta heterodimer; a lone beta-derived tail carries the sorting motif but not the partner, so this tests the targeting determinant, not the enzyme.

**R4** -- invariant chain, epitope in the CLIP register
- `II_CLIP_SCAFFOLD`: CD74 / invariant chain (Ii p31) (Mus musculus)
  - region: full ORF, with the CLIP segment excised and its boundaries recorded so the epitope can be dropped into that register
  - why: Type II membrane protein: the N-terminal cytoplasmic tail carries the di-leucine endosomal targeting motifs and must not be modified or tagged.
- note: Type II membrane protein. The N-terminus is CYTOSOLIC and carries the di-leucine targeting motifs, so it takes no signal peptide and no N-terminal tag. The epitope replaces CLIP in-register; the B cassette's flanks are wrong for this one and the Ea core is used instead.

**R5** -- macroautophagy: antigen into the autophagosome
- `LC3B`: MAP1LC3B (LC3B) (Mus musculus)
  - region: residues 1-120, ending at the glycine that ATG4 exposes; nothing may follow it
  - why: Autophagosome delivery. The C-terminal glycine must stay free for lipidation, so this module is terminal by construction and cannot carry a downstream tag or stop-adjacent element.
- note: Antigen must be N-terminal to LC3B and nothing may follow LC3B's C-terminal glycine, so this construct cannot carry E5 or a C-terminal tag.

