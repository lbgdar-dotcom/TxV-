# Where the two UTRs come from

The panel's transcripts are a **hybrid**. The 3' UTR is BNT162b2's; the 5' UTR
is a truncated mRNA-1273 leader. Neither vaccine matches both, which is why
comparing a construct against either one alone looks like a mismatch.

Reference: Jeong, McCoy, Artiles, Ilbay, Fire, Nadeau, Park, Betts, Boyd, Hoh
and Shoura, *Assemblies of putative SARS-CoV2-spike-encoding mRNA sequences for
vaccines BNT-162b2 and mRNA-1273*, v0.21Beta, Stanford University School of
Medicine (github.com/NAalytics). RNA was extracted from vial remnants and
sequenced directly, so these are read off the product rather than copied from a
patent. Figures 1 and 2 of that document are the source for every sequence
below.

## The four reference UTRs, and ours

| | length | identity |
|---|---|---|
| BNT162b2 5' UTR | 54 nt | human alpha-globin (hAg) leader |
| BNT162b2 3' UTR | **296 nt** | AES (136 nt) + mtRNR1 (139 nt), with scars |
| mRNA-1273 5' UTR | 57 nt | Moderna leader |
| mRNA-1273 3' UTR | 119 nt | — |
| **ours, 5'** | **47 nt** | truncated mRNA-1273 leader (see below) |
| **ours, 3'** | **296 nt** | BNT162b2's, byte-identical |

```
BNT162b2  5'  GAGAATAAACTAGTATTCTTCTGGTCCCCACAGACTCAGAGAGAACCCGCCACC
mRNA-1273 5'  GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGACCCCGGCGCCGCCACC
ours      5'  AGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACC
```

## The 3' UTR: BNT162b2's, unmodified

Verified twice, against two independent sources, base for base. Take the
BNT162b2 contig, drop everything through the spike ORF's tandem stop (TGA TGA),
and the remainder is our `UTR3` with no mismatch. It decomposes as:

| position | length | segment |
|---|---|---|
| 1-6 | 6 nt | `CTCGAG` (XhoI) |
| 7-142 | 136 nt | AES / TLE5 3' UTR core |
| 143-281 | 139 nt | mtRNR1 (mitochondrial 12S rRNA) |
| 282-296 | 15 nt | `ACCCTGGAGCTAGCA`, scar ending in NheI |

The two motif cores are **butt-joined** -- no linker, no spacer. That matters
because Orlandini von Niessen et al. (*Mol Ther* 2019;27:824-836) print the
motif core areas separately and never print a composite: the composite is the
concatenation, which is how all 64 of that paper's double-UTRs were built.
AES comes first (the element written F-I in the BioNTech literature, F = AES,
I = mtRNR1), even though the paper's best in-vitro performer is named the other
way round. Both were top of the screen.

Why this one: the screen that found these motifs was run **in human dendritic
cells**, and the composite was then validated in BALB/c mice -- evidence in
both the target cell type and the target species. Nothing here is translated,
so there is no codon usage to adapt for mouse.

## The 5' UTR: mRNA-1273's, truncated

```
mRNA-1273   GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGACCCCGGCGCCGCCACC   57 nt
ours        AGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGA..........GCCACC   47 nt
            ^                                        ^^^^^^^^^^
```

Two differences:

1. **G1 -> A.** Forced by the kit. The VENI manual specifies
   `5'-TAATACGACTCACTATAAGG`, which fixes the cap1 trinucleotide as the
   2'OMe-A form and requires the transcript to begin `AG`. mRNA-1273 begins
   `GGG`. Any leader put here needs this change.
2. **A 10 nt deletion of `CCCCGGCGCC`**, immediately 5' of the Kozak.

The 41 nt body and the `GCCACC` Kozak terminus are identical.

### Why the truncation is kept

Measured head to head over all nine mouse constructs (ViennaRNA, via
`scripts/validate_in_silico.py`):

| leader | GC | CpG (o/e) | cap-proximal MFE | AUG accessibility |
|---|---|---|---|---|
| ours, 47 nt | 34.0% | **0** (0.00) | **-4.4 to -6.6** | 0.513-0.542 |
| mRNA-1273, 57 nt, AGG-headed | 45.6% | 3 (1.04) | -5.4 to -11.1 | 0.454-0.547 |

Both clear every threshold, so this is a preference and not a fix. But the
deleted block carries all three of the element's CpGs -- at o/e 1.04 it would
be the only stretch of the transcript above the 0.80 the panel is held to,
which for a DC-targeted vaccine is a ZAP and innate-sensing lever worth
controlling deliberately -- and it roughly doubles cap-proximal structure,
because a GC-rich block directly 5' of the start codon is the textbook place
for an unwanted hairpin. It costs ~0.02 start accessibility on P5/P7/P8 only,
against a 0.30 floor.

Restoring the full licensed leader would buy provenance at the price of the
properties the leader is there for.

For a methods section: *"a truncated variant of the mRNA-1273 5' UTR,
re-headed AGG for cap1 analog compatibility and with the GC-rich 10 nt block
5' of the Kozak removed."*

## What is pinned in the suite

- `tests/test_ivt_unit.py::test_utr3_segment_decomposition` -- the 3' UTR's
  segment boundaries, so an edit has to say which segment it is editing.
- `tests/test_ivt_unit.py::test_leader_is_a_truncated_mrna_1273_utr5` -- the
  exact relationship above, so nothing silently restores the deleted block.
