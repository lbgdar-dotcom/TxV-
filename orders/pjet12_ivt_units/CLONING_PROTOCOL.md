# P0–P8 in pJET1.2: ordering and cloning

## The short answer to the orientation question

**You cannot control orientation with pJET1.2, and you don't need to.**

pJET1.2/blunt is a blunt-end vector. Blunt ends carry no sequence information,
so T4 ligase joins the insert in either direction at roughly 50/50. No choice of
restriction enzyme changes that, because **no restriction enzyme is used** —
blunt ligation is just ligation.

These fragments make orientation irrelevant instead. Each one is a **complete
transcription unit**:

```
[IVT_F handle] T7 promoter · 5' leader · Kozak · ORF · 3' UTR [IVT_R handle]
      20 nt         17 nt        47 nt              99 nt         20 nt
```

You recover the IVT template by **PCR with your own primer pair**, not by
cutting the plasmid. A PCR product is double-stranded, and the T7 promoter sits
on a defined strand *of the product*, so it transcribes correctly whichever way
the insert landed. The plasmid is only a way to store and amplify the sequence
in *E. coli* — which is all you asked it to do.

This was verified, not assumed: each fragment was inserted into a mock plasmid
in both orientations and the PCR simulated on each. Same amplicon both ways, all
nine constructs.

## What you are ordering

| | fragment | transcript | protein | GC |
|---|---|---|---|---|
| P0 | 554 bp | 497 nt | 116 aa | 56.1% |
| P1 | 656 bp | 599 nt | 150 aa | 57.2% |
| P2 | 656 bp | 599 nt | 150 aa | 57.3% |
| P3 | 758 bp | 701 nt | 184 aa | 57.8% |
| P4 | 758 bp | 701 nt | 184 aa | 57.8% |
| P5 | 470 bp | 413 nt | 88 aa | 56.4% |
| P6 | 632 bp | 575 nt | 142 aa | 59.5% |
| P7 | 755 bp | 698 nt | 183 aa | 60.0% |
| P8 | 755 bp | 698 nt | 183 aa | 60.0% |

All pass synthesis screening. Order as **dsDNA gene fragments**; no
phosphorylation is needed — pJET1.2/blunt supplies the 5' phosphates and the
remaining nicks are repaired in the cell.

Three primers, in `primers.csv`, and the same three work for all nine:

| | length | use |
|---|---|---|
| `IVT_F` | 20 nt | forward, anneals to the forward handle |
| `IVT_R_plain` | 20 nt | reverse, no tail — colony PCR and sequencing |
| `IVT_R_polyA120` | 140 nt | reverse with a T120 tail — adds the poly(A) to the IVT template |

`IVT_R_polyA120` is an Ultramer-class oligo. One synthesis covers all nine
constructs and every future prep, so the cost is a one-off.

## Cloning

1. **Ligate.** gBlock + pJET1.2/blunt + T4 DNA ligase, per the CloneJET
   protocol. Skip the blunting step — a gene fragment is already blunt; it is
   there for A-tailed PCR products.
2. **Transform and plate on ampicillin.** pJET1.2 is AmpR, not Kan — different
   from pVax1_AG. The lethal `eco47IR` gene means >99% of colonies carry an
   insert, so no blue/white and no no-insert control is needed.
3. **Screen for orientation only if you care.** You don't, for IVT. If you want
   to know anyway, `IVT_F` paired with a pJET1.2 sequencing primer gives a
   product in one orientation and not the other.
4. **Sequence.** The pJET1.2 forward and reverse sequencing primers read into
   the insert from both sides; at 470–758 bp a single read from each end covers
   the fragment. Confirm: the `GCCACC`**`ATG`** junction, no indels in the ORF,
   `GYQTI` ending the protein in P6/P7/P8, and P2A intact in P7/P8.
5. **Make the IVT template.** PCR with `IVT_F` + `IVT_R_polyA120` off the
   miniprep. The product is the transcription unit with a 120-bp A/T tract on
   the end. Column-clean it; that is your template.
6. **Transcribe.** T7, CleanCap AG (the transcript begins `AGG`), m1Ψ if you are
   substituting.

## Why there is no restriction enzyme step

Worth stating because it is the source of the confusion. Three ways to put a
fragment into a vector, and only the first two involve choosing enzymes:

- **Directional restriction cloning** — cut insert and vector with *two
  different* enzymes. The two ends get different, non-complementary overhangs,
  so only one orientation can anneal. This is where enzyme choice matters, and
  the rule is: both sites must be in the vector's MCS, in the right order, and
  **absent from your insert**.
- **Golden Gate** — a Type IIS enzyme cuts *outside* its own recognition site,
  so you choose the 4-base overhang yourself. Directional and scarless, and the
  enzyme site disappears from the product.
- **Blunt cloning (this)** — no overhangs, no enzyme, no direction. You trade
  orientation control for the simplest possible reaction and near-perfect
  positive selection.

We took the third and removed the need for the first.

## Two things this changes from the pVax1_AG design

**No CMV promoter anywhere.** The `AATAAA` and cryptic-splice screens existed
because a nuclear transcript would be spliced and prematurely polyadenylated.
With no eukaryotic promoter there is no nuclear transcript, so those checks
become sequence hygiene rather than functional constraints. They still pass.

**No encoded poly(A).** The tail comes from `IVT_R_polyA120`. That removes the
one element that is genuinely hard to synthesise, and it is why these fragments
order cleanly. If you would rather have it in the plasmid, say so — a segmented
A30–linker–A70 tract is synthesisable where a flat A120 is not.

## Sanity-checking before you order

Three views of the same nine molecules, so a mistake has to survive all three:

| file | what it is for |
|---|---|
| `P0.gb` … `P8.gb` | annotated maps — import into Benchling and look |
| `element_manifest.csv` | one row per construct, one column per element |
| `ELEMENT_CHECKLIST.md` | the same, written out with coordinates |

Counts in the manifest are of the peptide **found in the translated protein**,
not of what the design file claims — so a module that was declared but dropped
shows as absent.

```bash
python audit/audit_pjet.py           # 625 checks; exit 1 blocks ordering
python -m pytest audit/              # 34 negative controls
```

`audit/audit_pjet.py` imports nothing from the design library. It restates the
design from scratch and checks, per construct: the reading frame and translated
protein; every module intact at its expected position; SIINFEKL and the class II
core present; `GYQTI` ending the protein on LAMP1 constructs; E5 terminal where
it has a free C-terminus and internal in P7/P8; P2A present only in P7/P8; the
T7 promoter, the AGG start, the Kozak junction, the 3' UTR and both handles;
that no poly(A) is encoded; and that **the manifest and the GenBank record both
describe the molecule actually being ordered**. Across constructs it checks that
every module has exactly one DNA encoding.

The negative controls delete a signal peptide, E5, P2A, an antigen and the HA
tag; swap the antigen order in P3; break the T7 promoter, the AGG start and the
Kozak; frameshift an ORF; damage SIINFEKL; and make the manifest and the GenBank
disagree with the sequence. Each must block.
