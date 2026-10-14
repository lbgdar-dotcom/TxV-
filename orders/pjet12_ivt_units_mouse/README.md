# P0–P8, codon-optimised for mouse

The same nine constructs as `../pjet12_ivt_units/`, re-optimised against the
*Mus musculus* codon table. **Order from one package or the other, never a
mixture.**

## Why this package exists

The readout is murine: SIINFEKL is presented by H-2K(b), and antigen B's core
is an I-A(b) binder. The original package was optimised against the human
table, which is a mismatch with the cells the panel is read out in.

## What changed, and what did not

**The proteins are byte-identical.** Only synonymous codons moved, so every
length, every element count and every positional constraint is unchanged:

| | fragment | plasmid | nt differing from the human package | identity |
|---|---|---|---|---|
| P0 | 751 bp | 3725 bp | 5 | 99.3% |
| P1, P2 | 853 bp | 3827 bp | 11 | 98.7% |
| P3, P4 | 955 bp | 3929 bp | 22 | 97.7% |
| P5 | 670 bp | 3644 bp | 17 | 97.5% |
| P6 | 829 bp | 3803 bp | 32 | 96.1% |
| P7, P8 | 955 bp | 3929 bp | 38 | 96.0% |

**The primers are the same three oligos**, byte-identical to the human
package. They anneal to the handles, which are not coding sequence, so codon
choice does not touch them. If you already have them, you do not need new ones.

**How far the two codon tables actually are apart:** mean absolute difference
0.012 across 64 codons, largest 0.040 (Lys AAG/AAA and Asn AAC/AAT), and the
single most-used codon changes for exactly one amino acid — proline, human CCC
to mouse CCT. Two mammals are not far apart. This is a correctness fix, not a
large expected gain.

## Honest expectation

Do not expect a dramatic expression difference. Human-optimised constructs
express perfectly well in mouse cells, which is why this was not caught
earlier. The reason to use this package is that optimising for the readout
species is the right default when the readout species is known — not because
the human set would have failed.

The one thing the mismatch could never have done is confound the panel: all
nine constructs shared one table, so any effect applied equally to every arm.
That is still true here.

## What is in this folder

| file | what it is |
|---|---|
| `gblocks.fasta` | the nine sequences to order |
| `ordering_table.csv` | the same, with lengths, GC and synthesis flags |
| `primers.csv` | IVT_F, IVT_R_plain, IVT_R_polyA120 — shared with the human package |
| `element_manifest.csv` | element counts per construct, from the translated protein |
| `ELEMENT_CHECKLIST.md` | every element with its coordinates |
| `P0.gb` … `P8.gb` | annotated fragment maps |
| `benchling_maps/` | the nine finished pJET1.2 plasmids, native vector numbering |
| `IN_SILICO_REPORT.md` | folding, CAI, CpG — CAI scored against the **mouse** table |

## Cloning and IVT

Identical to the human package: `../pjet12_ivt_units/CLONING_PROTOCOL.md`. The
vector, the blunt site at 371/372, the BglII diagnostic and run-off, the
poly(A) PCR and the kit's promoter requirement are all unchanged — only the
codons inside the ORF differ. The protocol is not duplicated here on purpose,
so the two cannot drift apart.

## Rebuilding

```bash
python scripts/build_pjet_package.py mouse
python scripts/build_benchling_maps.py mouse
python scripts/validate_in_silico.py mouse
python audit/audit_pjet.py mouse
```
