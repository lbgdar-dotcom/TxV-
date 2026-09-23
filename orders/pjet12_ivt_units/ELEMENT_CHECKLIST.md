# Element checklist — P0–P8 as ordered for pJET1.2

Every element in every fragment, with its position in the ordered
sequence. Cross-check against `element_manifest.csv` (the same data as a
table) and the `.gb` files (the same data as annotated maps you can open
in Benchling).

Counts below are of the **peptide actually found in the translated**
protein, not of what the design file claims.

## What each tracked element is

| element | what it is | why it is there |
|---|---|---|
| `CTLA4_SP` | CTLA-4 signal peptide | routes into the ER; N-terminal only |
| `LAMP1_SP` | LAMP1 signal peptide | routes into the ER; N-terminal only |
| `M` | bare initiator Met | cytosolic constructs, no signal peptide |
| `HA` | HA tag | detects expression, independently of presentation |
| `FLAG` | FLAG tag | detects the second cistron in the dual-route pair |
| `A` | antigen A (SIINFEKL, 29-aa flanks) | MHC-I readout |
| `B` | antigen B (I-A(b) core, 29-aa flanks) | MHC-II readout |
| `L` | GGGGS linker | flexible spacer; 4 encodings cycled |
| `E5` | E5 acidic C-degron | proteasome turnover; needs a free C-terminus |
| `P2A` | P2A skip peptide | makes one ORF give two proteins |
| `CTLA4_TMT` | CTLA-4 TM + YVKM tail | AP-2, via the plasma membrane |
| `LAMP1_TMT` | LAMP1 TM + GYQTI tail | AP-3, direct to lysosome; must be last |

## Per construct

### P0 — CTLA-4 (AP-2, via surface)

Scaffold-only control: does the routing module alone produce signal?

Fragment **554 bp** · transcript **497 nt** · protein **116 aa** · antigens **none**

| position | element | kind |
|---|---|---|
| 1–20 | `IVT_F_handle` | handle |
| 21–37 | `T7_promoter` | promoter |
| 38–78 | `5'UTR` | utr5 |
| 79–84 | `Kozak` | kozak |
| 85–189 | `CTLA4_SP` | signal_peptide |
| 190–216 | `HA` | tag |
| 217–231 | `L` | linker |
| 232–417 | `CTLA4_TMT` | trafficking |
| 418–432 | `E5` | degron |
| 433–435 | `stop` | stop |
| 436–534 | `3'UTR` | utr3 |
| 535–554 | `IVT_R_handle` | handle |

**Present:** `CTLA4_SP`, `HA`, `L`, `E5`, `CTLA4_TMT`

**Absent (by design):** `LAMP1_SP`, `M`, `FLAG`, `A`, `B`, `P2A`, `LAMP1_TMT`

### P1 — CTLA-4 (AP-2, via surface)

Single-antigen (A) scaffold control.

Fragment **656 bp** · transcript **599 nt** · protein **150 aa** · antigens **A**

| position | element | kind |
|---|---|---|
| 1–20 | `IVT_F_handle` | handle |
| 21–37 | `T7_promoter` | promoter |
| 38–78 | `5'UTR` | utr5 |
| 79–84 | `Kozak` | kozak |
| 85–189 | `CTLA4_SP` | signal_peptide |
| 190–216 | `HA` | tag |
| 217–231 | `L_2` | linker |
| 232–318 | `A` | neoepitope |
| 319–333 | `L_4` | linker |
| 334–519 | `CTLA4_TMT` | trafficking |
| 520–534 | `E5` | degron |
| 535–537 | `stop` | stop |
| 538–636 | `3'UTR` | utr3 |
| 637–656 | `IVT_R_handle` | handle |

**Present:** `CTLA4_SP`, `HA`, `A`, `L`, `E5`, `CTLA4_TMT`

**Absent (by design):** `LAMP1_SP`, `M`, `FLAG`, `B`, `P2A`, `LAMP1_TMT`

### P2 — CTLA-4 (AP-2, via surface)

Single-antigen (B) scaffold control.

Fragment **656 bp** · transcript **599 nt** · protein **150 aa** · antigens **B**

| position | element | kind |
|---|---|---|
| 1–20 | `IVT_F_handle` | handle |
| 21–37 | `T7_promoter` | promoter |
| 38–78 | `5'UTR` | utr5 |
| 79–84 | `Kozak` | kozak |
| 85–189 | `CTLA4_SP` | signal_peptide |
| 190–216 | `HA` | tag |
| 217–231 | `L_2` | linker |
| 232–318 | `B` | neoepitope |
| 319–333 | `L_4` | linker |
| 334–519 | `CTLA4_TMT` | trafficking |
| 520–534 | `E5` | degron |
| 535–537 | `stop` | stop |
| 538–636 | `3'UTR` | utr3 |
| 637–656 | `IVT_R_handle` | handle |

**Present:** `CTLA4_SP`, `HA`, `B`, `L`, `E5`, `CTLA4_TMT`

**Absent (by design):** `LAMP1_SP`, `M`, `FLAG`, `A`, `P2A`, `LAMP1_TMT`

### P3 — CTLA-4 (AP-2, via surface)

Can one transcript drive simultaneous class I and class II presentation? Order A-then-B.

Fragment **758 bp** · transcript **701 nt** · protein **184 aa** · antigens **A+B**

| position | element | kind |
|---|---|---|
| 1–20 | `IVT_F_handle` | handle |
| 21–37 | `T7_promoter` | promoter |
| 38–78 | `5'UTR` | utr5 |
| 79–84 | `Kozak` | kozak |
| 85–189 | `CTLA4_SP` | signal_peptide |
| 190–216 | `HA` | tag |
| 217–231 | `L_2` | linker |
| 232–318 | `A` | neoepitope |
| 319–333 | `L_4` | linker |
| 334–420 | `B` | neoepitope |
| 421–435 | `L_6` | linker |
| 436–621 | `CTLA4_TMT` | trafficking |
| 622–636 | `E5` | degron |
| 637–639 | `stop` | stop |
| 640–738 | `3'UTR` | utr3 |
| 739–758 | `IVT_R_handle` | handle |

**Present:** `CTLA4_SP`, `HA`, `A`, `B`, `L`, `E5`, `CTLA4_TMT`

**Absent (by design):** `LAMP1_SP`, `M`, `FLAG`, `P2A`, `LAMP1_TMT`

### P4 — CTLA-4 (AP-2, via surface)

Order control for P3: does cassette order change the outcome?

Fragment **758 bp** · transcript **701 nt** · protein **184 aa** · antigens **B+A**

| position | element | kind |
|---|---|---|
| 1–20 | `IVT_F_handle` | handle |
| 21–37 | `T7_promoter` | promoter |
| 38–78 | `5'UTR` | utr5 |
| 79–84 | `Kozak` | kozak |
| 85–189 | `CTLA4_SP` | signal_peptide |
| 190–216 | `HA` | tag |
| 217–231 | `L_2` | linker |
| 232–318 | `B` | neoepitope |
| 319–333 | `L_4` | linker |
| 334–420 | `A` | neoepitope |
| 421–435 | `L_6` | linker |
| 436–621 | `CTLA4_TMT` | trafficking |
| 622–636 | `E5` | degron |
| 637–639 | `stop` | stop |
| 640–738 | `3'UTR` | utr3 |
| 739–758 | `IVT_R_handle` | handle |

**Present:** `CTLA4_SP`, `HA`, `A`, `B`, `L`, `E5`, `CTLA4_TMT`

**Absent (by design):** `LAMP1_SP`, `M`, `FLAG`, `P2A`, `LAMP1_TMT`

### P5 — cytosolic (+E5 degron, free C-terminus)

Cytosolic baseline, and the only construct where the degron has a free C-terminus -- so this is the degron-positive comparator.

Fragment **470 bp** · transcript **413 nt** · protein **88 aa** · antigens **A+B**

| position | element | kind |
|---|---|---|
| 1–20 | `IVT_F_handle` | handle |
| 21–37 | `T7_promoter` | promoter |
| 38–78 | `5'UTR` | utr5 |
| 79–84 | `Kozak` | kozak |
| 85–87 | `M` | start |
| 88–114 | `HA` | tag |
| 115–129 | `L_2` | linker |
| 130–216 | `A` | neoepitope |
| 217–231 | `L_4` | linker |
| 232–318 | `B` | neoepitope |
| 319–333 | `L_6` | linker |
| 334–348 | `E5` | degron |
| 349–351 | `stop` | stop |
| 352–450 | `3'UTR` | utr3 |
| 451–470 | `IVT_R_handle` | handle |

**Present:** `M`, `HA`, `A`, `B`, `L`, `E5`

**Absent (by design):** `CTLA4_SP`, `LAMP1_SP`, `FLAG`, `P2A`, `CTLA4_TMT`, `LAMP1_TMT`

### P6 — LAMP1 (AP-3, direct to lysosome)

Does the direct TGN -> lysosome route beat the detour via the cell surface that P0-P4 take?

Fragment **632 bp** · transcript **575 nt** · protein **142 aa** · antigens **A+B**

| position | element | kind |
|---|---|---|
| 1–20 | `IVT_F_handle` | handle |
| 21–37 | `T7_promoter` | promoter |
| 38–78 | `5'UTR` | utr5 |
| 79–84 | `Kozak` | kozak |
| 85–156 | `LAMP1_SP` | signal_peptide |
| 157–183 | `HA` | tag |
| 184–198 | `L_2` | linker |
| 199–285 | `A` | neoepitope |
| 286–300 | `L_4` | linker |
| 301–387 | `B` | neoepitope |
| 388–402 | `L_6` | linker |
| 403–510 | `LAMP1_TMT` | trafficking |
| 511–513 | `stop` | stop |
| 514–612 | `3'UTR` | utr3 |
| 613–632 | `IVT_R_handle` | handle |

**Present:** `LAMP1_SP`, `HA`, `A`, `B`, `L`, `LAMP1_TMT`

**Absent (by design):** `CTLA4_SP`, `M`, `FLAG`, `E5`, `P2A`, `CTLA4_TMT`

### P7 — dual: A cytosolic + B lysosomal

Dual route from one transcript: A to the proteasome, B to the endolysosome. Note E5 is internal here, upstream of P2A, so it has no free C-terminus and is expected to be weak.

Fragment **755 bp** · transcript **698 nt** · protein **183 aa** · antigens **A+B**

| position | element | kind |
|---|---|---|
| 1–20 | `IVT_F_handle` | handle |
| 21–37 | `T7_promoter` | promoter |
| 38–78 | `5'UTR` | utr5 |
| 79–84 | `Kozak` | kozak |
| 85–87 | `M` | start |
| 88–114 | `HA` | tag |
| 115–129 | `L_2` | linker |
| 130–216 | `A` | neoepitope |
| 217–231 | `L_4` | linker |
| 232–246 | `E5` | degron |
| 247–312 | `P2A` | skip_peptide |
| 313–384 | `LAMP1_SP` | signal_peptide |
| 385–408 | `FLAG` | tag |
| 409–423 | `L_9` | linker |
| 424–510 | `B` | neoepitope |
| 511–525 | `L_11` | linker |
| 526–633 | `LAMP1_TMT` | trafficking |
| 634–636 | `stop` | stop |
| 637–735 | `3'UTR` | utr3 |
| 736–755 | `IVT_R_handle` | handle |

**Present:** `LAMP1_SP`, `M`, `HA`, `FLAG`, `A`, `B`, `L`, `E5`, `P2A`, `LAMP1_TMT`

**Absent (by design):** `CTLA4_SP`, `CTLA4_TMT`

### P8 — dual: B cytosolic + A lysosomal

Swap control for P7: is the effect about the route or about the antigen?

Fragment **755 bp** · transcript **698 nt** · protein **183 aa** · antigens **B+A**

| position | element | kind |
|---|---|---|
| 1–20 | `IVT_F_handle` | handle |
| 21–37 | `T7_promoter` | promoter |
| 38–78 | `5'UTR` | utr5 |
| 79–84 | `Kozak` | kozak |
| 85–87 | `M` | start |
| 88–114 | `HA` | tag |
| 115–129 | `L_2` | linker |
| 130–216 | `B` | neoepitope |
| 217–231 | `L_4` | linker |
| 232–246 | `E5` | degron |
| 247–312 | `P2A` | skip_peptide |
| 313–384 | `LAMP1_SP` | signal_peptide |
| 385–408 | `FLAG` | tag |
| 409–423 | `L_9` | linker |
| 424–510 | `A` | neoepitope |
| 511–525 | `L_11` | linker |
| 526–633 | `LAMP1_TMT` | trafficking |
| 634–636 | `stop` | stop |
| 637–735 | `3'UTR` | utr3 |
| 736–755 | `IVT_R_handle` | handle |

**Present:** `LAMP1_SP`, `M`, `HA`, `FLAG`, `A`, `B`, `L`, `E5`, `P2A`, `LAMP1_TMT`

**Absent (by design):** `CTLA4_SP`, `CTLA4_TMT`

## Primers (the same three for all nine)

| name | length | sequence | use |
|---|---|---|---|
| IVT_F | 20 nt | `CAGTCACGTAGCATCGACTG` | forward |
| IVT_R_plain | 20 nt | `TCAGTGATCGACTGCATGAC` | colony PCR, sequencing |
| IVT_R_polyA120 | 140 nt | `T`×120 + `TCAGTGATCGACTGCATGAC` | adds the poly(A) to the IVT template |

