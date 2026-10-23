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

### The 3' end straddles a junction

The last 40 nt, `GGGTTGGTCAATTTCGTGCCAGCCACACCCTGGAGCTAGC`, is the stretch
most likely to look unaccountable, because it is not one element:

```
GGGTTGGTCAATTTCGTGCCAGCCAC   ACCCTGGAGCTAGC   A
|---- mtRNR1, 256-281 ----|   |- scar 282-295 -|  296
          26 nt                     14 nt
```

Aligned against human mitochondrial DNA it will match and then stop dead at
`...GCCAGCCAC`, because the native 12S rRNA continues `CGCGGT...` where this
switches to the scar. The divergence point is the junction, not an error.

The scar is cloning debris: BioNTech assembled the composite in a plasmid, so
the element carries the sites used to build it at its outer ends -- XhoI at
1-6, NheI at 290-295. Both are inert here; nothing in this workflow cuts with
either.

The last 61 nt of any fragment in this panel are therefore:

```
GGGTTGGTCAATTTCGTGCCAGCCACACCCTGGAGCTAGCA GTCATGCAGTCGATCACTGA
|---------- end of 3' UTR, 41 nt ---------| |- reverse handle --|
```

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

## Where the 5' UTR's boundary falls

In a Benchling map the leader reads 44 nt, not 47, whenever the T7 promoter is
annotated in the kit's 20 nt form `TAATACGACTCACTATAAGG`: that feature swallows
the `AGG`, and the 5' UTR feature starts after it.

```
TAATACGACTCACTATA AGG AAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACC ATG
|-- T7 core, 17 -| |- transcribed: 47 nt of 5' UTR --------------|
|--- KIT_MIN_PROMOTER, 20 -| |- 5' UTR feature as drawn, 44 nt --|
```

Both are right about the DNA. They disagree about which feature owns the `AGG`,
and the `AGG` is genuinely both things: the last three bases T7 needs to
initiate well, and the first three bases of the transcript. On the RNA there is
no ambiguity -- the 5' UTR is **47 nt**, cap to `A` of `AUG`, and that is the
number for a methods section.

### Expected feature coordinates in the pJET1.2 map

The 5' end is identical on all nine constructs, so these positions hold for
every one of them. Insert goes in blunt at 372, forward handle is 20 nt, so
the T7 core begins at 372 + 20 = 392.

| feature | 1-based, inclusive | length |
|---|---|---|
| forward primer handle | 372-391 | 20 |
| T7 core `TAATACGACTCACTATA` | 392-408 | 17 |
| cap1 `AGG` | 409-411 | 3 |
| 5' UTR body | 412-449 | 38 |
| Kozak `GCCACC` | 450-455 | 6 |
| `ATG` | 456-458 | 3 |

The T7 promoter feature has no canonical length, so the table's 392-408 is one
convention and not the only right answer. Annotations differ in how many
transcribed bases they pull into the promoter:

| span | sequence | convention |
|---|---|---|
| 392-408, 17 nt | `TAATACGACTCACTATA` | core only, `AGG` annotated separately |
| 392-409, 18 nt | `TAATACGACTCACTATAA` | core + 1 |
| 392-410, 19 nt | `TAATACGACTCACTATAAG` | core + 2 -- what **pJET1.2 itself uses** |
| 392-411, 20 nt | `TAATACGACTCACTATAAGG` | the kit manual's minimum template |

The vector's own `T7 promoter` feature is 19 nt (`TAATACGACTCACTATAGG`, 305-323
in stock pJET1.2), so the vendor already counts two transcribed bases as part
of the promoter. Any of these is defensible; pick one and keep it the same
across all nine maps.

Note the manual's wording: *"Minimum T7 promotor sequences:
5'-TAATACGACTCACTATAAGG"*. That constrains the **molecule** -- those 20 bases
must be present on the template -- not the feature boundary on a map. The DNA
satisfies it however the annotation is drawn.

A feature that overlaps the `AGG` by a base or two is not an error. Base 409 is
genuinely both the promoter's last annotated base, under the core+1 convention,
and the first base of the cap1 `AGG`; Benchling allows overlapping features
precisely for this.

### Three numbers for one 5' UTR

Depending on which flanking pieces a map attributes elsewhere, the same DNA
reads as 38, 44 or 47 nt:

| | what it covers |
|---|---|
| 38 nt | body only -- `AGG` and Kozak annotated as their own features |
| 44 nt | body + Kozak -- `AGG` inside the 20 nt T7 promoter feature |
| **47 nt** | `AGG` + body + Kozak -- **cap to AUG, the 5' UTR on the RNA** |

All three are right about the DNA. 47 is the one for a methods section.

## How these compare to the literature

**5' UTR -- Ma et al., *Sci Rep* 2024;14:19845** compared five leaders. As IVT
mRNA the order was beta-globin > albumin > minimal > alpha-globin > CYBA. Ours
was not among them. The two nearest comparators bracket it: their "minimal"
leader (promoter, transcription start and Kozak, nothing else) placed 3rd, and
alpha-globin -- which is BNT162b2's leader, the obvious swap -- placed 4th.
So the licensed alternative is not an upgrade by this paper's data.

Read the paper's own headline result first: plasmid-screened and
IVT-mRNA-delivered rankings did not correlate, and CYBA went from best in
plasmid to worst as mRNA. A ranking that does not survive a change of format
should not be assumed to survive a change of cell type either -- theirs was
HEK293T and BHK-T7, not dendritic cells, and not mouse.

**3' UTR -- Orlandini von Niessen et al., *Mol Ther* 2019;27:824-836** is not a
comparator for ours; it is its source. What ships here is that paper's
AES-mtRNR1 element as BNT162b2 manufactures it, which differs from the paper's
own construct in three small ways: mtRNR1 trimmed 142 -> 139 nt, two C->U
changes in the AES segment, and the XhoI/NheI cloning scars at the outer ends.

## The TMSB10 alternative, evaluated

A candidate pair raised for comparison, from two papers by what appears to be
one group. Both UTRs come from **TMSB10** (thymosin beta-10, NM_021103.4), a
very abundant human transcript -- the same reasoning that puts alpha-globin in
BNT162b2. Not adopted; recorded so the comparison is not redone from scratch.

The screen (MDPI *Vaccines*, Table 1) ranked UTRs from abundant human genes --
TMSB10, ACTB, TMSB4X, B2M, HLA-DRA, HLA-DRB1, FTL, FTH1 and a designed "R27" --
with a tPA-Gaussia luciferase reporter. The vaccine paper (*Sci Rep*, Table S1)
then built SARS-CoV-2 candidates on the TMSB10 pair.

**5' UTR, as actually built** (not as annotated): `AGGG` + the 78 nt TMSB10
5' UTR + `GCCACC` = 88 nt. It begins `AGGG`, so the implied promoter is
`TAATACGACTCACTATAAGGG`, which starts with this kit's minimum
`TAATACGACTCACTATAAGG` -- **compatible with no modification at all**, unlike
the BNT162b2 and mRNA-1273 leaders, which both need re-heading. It ends in the
same `GCCACC` Kozak used here.

**3' UTR, two versions.** The screen used the native 248 nt TMSB10 3' UTR; the
vaccine paper cut it to 58 nt, keeping the native 18 nt head
`GATCCTGGAGGATTTCCT`, the 8 nt tail `TCGTGGCA` and two short internal blocks,
and replacing the rest. The short form has lost the `AATAAA`, which is inert in
an IVT reaction. That a group holding the full native element chose to ship 12%
of it is itself a data point about how much of it was earning its length.

Measured head to head in P0 (mouse ORF, ViennaRNA):

| 5' leader | AUG accessibility | cap-proximal MFE |
|---|---|---|
| ours, 47 nt | **0.529** | **-6.6** |
| TMSB10 as used, 88 nt | 0.519 | -20.1 |

Both clear the floors (0.30 and -25.0). Accessibility is a tie; ours carries a
third of the cap-proximal structure. TMSB10 costs 6 CpGs (o/e 0.87) against
zero. Swapping the 3' UTR moved nothing at the 5' end, as expected 600+ nt
away -- and what a 3' UTR actually does, modulate decay and duration, ViennaRNA
cannot score, so that half of the comparison is not decidable here.

Length, on P0 (751 nt today): 5' +41, 3' -48 with the full form or -238 with
the short one, so net **-7 nt** or **-197 nt**.

Synthesis flags, native 248 nt form only: a 9 nt poly-C run at 152 and a BsaI
site at 138. The 58 nt form has neither.

**Not adopted.** The 3' UTR in use was selected in human dendritic cells,
validated in BALB/c mice and ships in a licensed product; TMSB10's comes from a
reporter screen whose comparator set did not include it. The TMSB10 5' UTR is a
real candidate and the only one so far that needs no re-heading, but it does
not beat the current leader on anything measurable here.

## If the leader is switched to HBB

Human beta-globin's 5' UTR, 50 nt, is
`ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACC` -- the leader Ma et al.
2024 ranked first as IVT mRNA. Two things to get right if it is adopted.

**It does not start `AG`.** It begins `AC`, so the `AGG` has to stay upstream of
it. Kept, the transcript reads `AGGACATTTGCTT...` and the implied promoter is
exactly the kit's minimum. Dropped, the transcript begins `AC` and the cap1
analog cannot incorporate -- a silent, total failure, since transcription still
works.

**Its native end is already a Kozak**, so appending `GCCACC` is optional:
`...ACAGACACC|ATG|G` puts a purine at -3 and G at +4, both satisfied.

Appending it anyway is harmless but not free, and the obvious middle option is
the worst of the three. Measured across all nine mouse constructs:

| leader | nt | CpG | AUG accessibility | cap MFE |
|---|---|---|---|---|
| `AGG` + HBB + `GCCACC` + `GCCACC` | 65 | 2 | 0.511-0.555 | -14.0 |
| `AGG` + HBB + `GCCACC` | 59 | 1 | **0.439-0.541** | -15.7 |
| `AGG` + HBB (native Kozak) | 53 | **0** | **0.542-0.791** | -12.6 to -14.9 |

Deleting one of a duplicated `GCCACC` -- the tidy-looking fix -- lowers start
accessibility on every construct; deleting both raises it. The duplicate
evidently breaks a hairpin the single copy forms. All three clear the 0.30 and
-25.0 floors, the spread is small, and these are ViennaRNA numbers on
unmodified RNA while the transcript is fully m1Psi-substituted, so this ranks
options rather than measuring them. It is still a reason not to assume the
tidy edit is the safe one.

A duplicated Kozak is in itself benign: what would hurt is an upstream AUG, and
`GCCACCGCCACC` contains none. Whatever is chosen has to go on all nine
identically, or initiation context becomes a second variable alongside routing.

**The published element carries one `GCCACC`**, so the middle row is the one
with literature behind it and the table poses no dilemma: a duplicate is an
editing artifact to remove, and dropping the Kozak entirely scores a little
better in the model but matches nothing anyone has measured. The 59 nt form is
the one to build.

## What is pinned in the suite

- `tests/test_ivt_unit.py::test_utr3_segment_decomposition` -- the 3' UTR's
  segment boundaries, so an edit has to say which segment it is editing.
- `tests/test_ivt_unit.py::test_leader_is_a_truncated_mrna_1273_utr5` -- the
  exact relationship above, so nothing silently restores the deleted block.
