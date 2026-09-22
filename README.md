# txv — IVT mRNA construct design for multi-antigen cancer vaccines

Design, QC, export and Benchling registration for in-vitro-transcribed mRNA
vaccine constructs, with LNP formulation arithmetic.

Two things live here:

1. **A general engine** — codon optimisation, cassette assembly, junctional
   epitope control, construct QC, GenBank export, LNP stoichiometry.
2. **The `pVax1_AG` routing panel** (`txv.pvax1_ag`) — the audited P0–P8 module
   library and construct series, carried over as a regression-locked library.

No third-party dependencies for anything except live Benchling calls.

```bash
pip install -e ".[dev]"            # engine + tests
pip install -e ".[benchling,dev]"  # + live Benchling registration
```

---

## Construct architecture

```
[T7 promoter] [5' UTR] [Kozak] [============ ORF ============] [3' UTR] [poly(A)] [linearisation]
                                 |        |                |
                    signal peptide   antigen cassette   trafficking domain + stop
```

The **whole ORF is codon-optimised in a single pass** — signal peptide, every
bead, every linker, the trafficking domain. That matters because forbidden
motifs and GC excursions are created at element *boundaries*, which an
element-by-element optimiser cannot see. The price is that changing one bead
re-optimises the construct, which is the right trade for a drug substance.

A visible consequence: in a cassette with eight identical `GGGGS` linkers, every
linker gets a **different** codon encoding, because the optimiser penalises long
direct repeats — which otherwise recombine during plasmid propagation.

---

## Quick start

```bash
txv design examples/antigens_example.csv --name TXV-001 --out out/
txv panel --out panel/            # the audited P0-P8 series
txv order --out order/            # synthesis-ready fragments to send a vendor
txv series antigens.csv --axis multiplex --out series/
txv formulate --fasta out/TXV-001.mrna.fasta --mrna-ug 100 --np 6
txv parts                         # what's in the registry, and what's a placeholder
```

```python
from txv.epitopes import Antigen
from txv.pipeline import design_construct, write_outputs

report = design_construct("TXV-001", [
    Antigen("KRAS_G12D", "LVVVGADGVGKSALTIQLIQNHFVD", mutation_offset=6),
    Antigen("NYESO1_157", "SLLMWITQC", kind="tumor_associated"),
    Antigen("PADRE", "AKFVAAWTLKAAA", kind="helper"),
])
print(report.to_text())
write_outputs(report, "out/")
```

---

## The default UTRs

The defaults are the **BNT162b2 UTRs**, extracted from the published sequence
rather than transcribed from memory. The extraction is self-validating: the ORF
between them is exactly 3822 nt in frame, 1273 aa, correct spike N-terminus,
tandem TGA·TGA.

| Part | Length | What it is |
|---|---|---|
| `UTR5_hAg` | 48 nt | Human α-globin (HBA1/HBA2-shared) 5′ UTR behind a 14-nt vector leader. With `kozak_strong` (`GCCACC`) it reconstitutes the validated 54-nt leader exactly — a test asserts this. |
| `UTR5_hAg_core` | 37 nt | The α-globin UTR without the vector leader, prefixed `GGG` for T7. Use if the scars below matter to you. |
| `UTR3_AES_mtRNR1` | 296 nt | AES + mtRNR1 composite, selected *ex vivo* for stability and total protein output. |
| `polyA_120` | 120 nt | Encoded tail. |

**What the validated sequences actually contain**, because you should know
before picking a restriction enzyme: the 5′ leader has `AATAAA` at position 3
and an SpeI site at 8; the 3′ UTR opens with an XhoI scar, carries NheI at 27
and 289, and has `AATAAA` at 219 inside mtRNR1. All of it is in the clinical
sequence. `AAUAAA` only directs cleavage in the *nucleus*, so it is inert in a
cytoplasmically delivered transcript — which is why QC attributes each motif to
its feature and warns **only** for hits in the ORF, where the optimiser could
have avoided them. A hit in a fixed part you chose is reported as a fact, not a
defect.

### Placeholder parts are not real sequences

A default build now contains **no placeholders**. `UTR5_placeholder`,
`UTR3_placeholder` and `MITD_placeholder` remain in the registry as stand-ins,
but nothing selects them by default, and `ConstructSpec.trafficking` defaults to
`None` rather than a placeholder — real routing comes from `txv.pvax1_ag`.

Every part carries a `Provenance` tag (`CANONICAL` / `LITERATURE` /
`PLACEHOLDER`); any construct containing a placeholder reports it, and QC
**fails** on one unless you pass `allow_placeholders=True`.

Supply your own with a JSON parts file:

```json
{"parts": [{"name": "UTR5_house", "kind": "utr5", "provenance": "literature",
            "dna": "GGGA...", "source": "internal ref 2024-11"}]}
```

### The built-in junction scorer is not a binding predictor

`AnchorMotifScorer` matches position-2 and C-terminal anchor residues for six
supertypes. It is there so the pipeline runs without a predictor installed, and
it is only good for *relative ranking*. For anything real, wrap NetMHCpan,
MHCflurry or your in-house model as a `JunctionScorer` (a one-method protocol)
and pass it in.

### Structure checks are proxies

The cap-proximal checks report GC and the longest perfect complementary stem —
not a folding free energy. Nothing here computes a thermodynamic ensemble. Run
ViennaRNA over `construct.transcript` if you need one; QC takes custom checks
via `extra=`.

---

## What each module does

| Module | Responsibility |
|---|---|
| `seqops` | Alphabet, translation, GC, motif search (IUPAC, both strands) |
| `codon` | Beam-search codon optimisation under competing objectives |
| `parts` | Part registry with provenance; JSON overrides |
| `epitopes` | Cassette model, junctional-epitope scanning, linker/order choice |
| `constructs` | Assembly of the annotated DNA template |
| `qc` | Check battery with `PASS` / `WARN` / `FAIL` severities |
| `export` | GenBank and FASTA |
| `lnp` | RNA mass, N/P stoichiometry, lipid masses, mixing volumes, dosing |
| `benchling` | Payload construction, dry run, live registration |
| `pipeline` | The end-to-end flow and file I/O |
| `pvax1_ag` | The audited P0–P8 module library, panel and routes |
| `order` | Synthesis fragments, assembly overlaps, vendor-failure screening |
| `variants` | Matched construct series for the two experimental axes |

### Codon optimisation is a constrained search, not a lookup

Picking the most frequent codon everywhere is wrong for mRNA drug substance,
because the objectives conflict:

* **Translation efficiency** wants common codons.
* **Innate immunogenicity and IVT fidelity** want *low uridine* — uridine drives
  RIG-I/TLR sensing, and long U runs are where T7 slips. m1Ψ reduces but does
  not remove the incentive.
* **Manufacturability** wants a GC band, no linearisation sites, no cryptic
  splice sites, no long homopolymers, no long direct repeats.

So `CodonOptimizer` runs a beam search over synonymous codons, scoring each
partial sequence on a weighted sum of these, with forbidden motifs penalised as
they are created at the codon boundary, then a repair pass that
re-synonymises locally around anything that survived. The weights are in
`OptimizerConfig` so the trade-off is an explicit, reviewable parameter of a
design rather than a property of the tool.

### Junctional epitopes are the multi-antigen failure mode

Fusing bead C-terminus to bead N-terminus creates peptides that exist nowhere in
the patient's proteome or tumour. They compete for presentation and are not a
target. The pipeline settles **linker choice first** (it dominates), then
**bead order** (2-opt; exhaustive at ≤7 beads). On a representative 8-antigen
set the built-in scorer ranks:

```
GGSGGGGSGG   1.67   <- chosen
GGGGS        4.83
GPGPG        9.00
RAKR        22.67
AAY         22.83
```

`fixed_first` / `fixed_last` pin terminal beads when they must abut a signal
peptide or a trafficking domain.

---

## The pVax1_AG routing panel (`txv.pvax1_ag`)

Twelve verified modules that reconstruct all nine audited ORFs exactly;
`tests/test_pvax1_ag.py` locks this against the reference translations, so
editing a module fails the suite by design.

```bash
txv panel --only P3 P7 --out panel/
```

**The mechanism the panel tests.** MHC-I is loaded in the ER with peptides the
proteasome makes in the cytosol; MHC-II is loaded in the MIIC, a late endosome.
mRNA is translated on free cytosolic ribosomes, so a bare antigen ORF is
CD8-biased *by construction* — the protein is made correctly and is simply in
the wrong compartment. That is the mechanistic core of "expressed but not
presented", and it is why a Western blot cannot diagnose it.

Redirecting it takes a **matched pair** of signals, and neither works alone:

* a **signal peptide** to get the chain into the ER lumen — SRP reads the first
  thing out of the ribosome, so this only works at the extreme N-terminus;
* a **transmembrane domain** to stop it being secreted, leaving the cassette
  luminal (ER lumen → endosome lumen → lysosome lumen, where cathepsins and
  MHC-II are) and the sorting tail cytosolic, where sorting machinery reads it.

The tail is the postcode:

| Tail | Motif | Adaptor | Route | Positional constraint |
|---|---|---|---|---|
| CTLA-4 | `YVKM` | AP-2 | ER → Golgi → surface → endocytosis → MIIC | internal, 19 native residues follow — tolerates downstream sequence |
| LAMP1 | `GYQTI` | AP-3 | direct TGN → lysosome | **must end the protein** — append anything and μ3A binding fails |

That constraint is what fixes module order in P7/P8: the LAMP1 arm must sit
*downstream* of P2A, the CTLA-4 arm may sit upstream.

| | Route | Question it answers |
|---|---|---|
| P0 | CTLA-4 | Scaffold-only control |
| P1 / P2 | CTLA-4 | Single-antigen (A / B) scaffold controls |
| P3 / P4 | CTLA-4 | Simultaneous class I + II from one transcript; order control |
| P5 | cytosolic + CL1 degron | Proteasome-directed baseline |
| P6 | LAMP1 | Does the direct route beat the detour via the surface? |
| P7 / P8 | dual, via P2A | One transcript, both routes; antigen-swap control |

**Two elements that look like padding and are not.** The 29-aa native flanks:
ERAP1 trims N-terminal extensions in the ER, but *nothing* trims a C-terminal
one, so the proteasome must make the epitope's C-terminus exactly right first
time, and cleavage specificity is set by surrounding residues. The HA tag: it
answers "was it made" on a different reagent from the pMHC antibodies that
answer "was it presented" — without it, an absent pMHC signal has three
indistinguishable causes.

**The interpretive trap.** Membrane-routed constructs still present on MHC-I,
via ERAD retrotranslocation, failed translocation and DRiPs. **MHC-I signal is
therefore not evidence that routing worked** — only the class II readout
separates the routes. Plan conclusions accordingly.

`OPEN_DECISIONS` carries the unresolved calls from the audit (encoded vs
PCR-added A120; P1/P2 are not route-matched controls; the HindIII site at ORF
position 251 in the old P1/P2 builds, which re-optimisation removes; pVax1_AG
is not interchangeable with commercial pVAX1).

---

## Ordering the plasmid

Two decisions turn a construct into an order, and both are easy to get wrong
silently. `txv order` makes them explicit.

**What are you replacing in the backbone?** `--mode orf` (default) orders only
the coding sequence and drops it between the backbone's existing UTR rails —
minimal change, and the construct's own UTRs are *not* used. `--mode cassette`
orders 5′ UTR + ORF + 3′ UTR and replaces the backbone's UTRs, which is what
you want if upgrading the UTRs is the point — but it **moves the IVT
primer-annealing rails**, so the primers must be re-designed. It is not a
drop-in, and the tool says so on every cassette fragment.

**Is the poly(A) tail encoded or added by PCR?** This one has a concrete
answer. A 120-nt A-tract is genuinely hard to synthesise — most vendors fail,
refuse, or silently deliver a contracted tract — and it contracts again in
*E. coli*. The screen catches it:

```
P3      1053 bp  GC 51.8%  CHECK
         ! 121-nt A-tract at 912: most vendors cannot synthesise this reliably
         ! 81 exact repeat(s) of 40+ nt
```

If your IVT reverse primer already anneals upstream of the encoded tract and
adds the tail by PCR, then the encoded tract does **no work for the transcript**
and only makes synthesis and propagation harder. `--polya-mode pcr_added` (the
default) drops it, and all nine panel fragments then order clean at 346–631 bp
with no flags.

The one thing that assumption rests on is worth measuring once, since it is
cheap: run the IVT product on a denaturing gel or TapeStation against an
untailed control, or do an RNase H / oligo-dT assay. A tail that was never
incorporated looks exactly like a tail that was, until you check.

`txv order --out order/` writes `ordering_table.csv` (vendor-ready, sequence
included), `gblocks.fasta`, and an annotated `*_gblock.gb` per fragment for
Benchling import. The GenBank records deliberately show only what you are
buying — in ORF mode the UTRs come from the backbone, so a record displaying
them would misrepresent the order.

A generated order for the P0-P8 panel, with a step-by-step cloning protocol, is
checked in at [`orders/pvax1_ag_panel/`](orders/pvax1_ag_panel/).

`assemble_into_backbone()` simulates the assembly against the parent plasmid and
returns the full circular sequence. It refuses rather than guessing if either
homology rail is not unique in the backbone.

---

## Two experiments the panel cannot answer

The panel settles *where* an antigen goes. `txv.variants` builds matched series
for the two questions that decide whether a multi-antigen vaccine is potent.
Each moves one variable and holds the rest, because that is the only form in
which the answer is interpretable.

### How many antigens before presentation decays?

Published routing work is mostly one or two antigens. Whether there is a
**presentation budget** — a point past which adding antigens dilutes per-epitope
pMHC instead of adding breadth — is open, and it is the question that decides
how many antigens a cassette should carry.

```bash
txv series antigens.csv --axis multiplex --sizes 1,2,4,6,8 --route ctla4
```

Two properties make the series interpretable, and both are deliberate.
**Nesting**: arm *k* contains every antigen of arm *k−1*, so a drop cannot be
explained by a different antigen set. **An anchor antigen** in position 1 of
every arm, as an internal standard — normalise the others to it and a global
translation drop becomes separable from genuine per-epitope dilution.

### Does uridine depletion trade expression against adjuvanticity?

Depleting uridine raises expression and lowers innate sensing. For a
prophylactic vaccine that is unambiguously good. For a **cancer** vaccine it may
not be: the same sensing you suppress is part of what matures the DC and
licenses it to prime CD8 T cells. So the potency-optimal uridine content may not
be the expression-optimal one.

```bash
txv series antigens.csv --axis uridine --route ctla4
```

Uridine **cannot** be varied independently of codon adaptation — depleting U
forces synonymous choices away from the preferred codon — so the series reports
CAI and GC per arm as covariates rather than pretending the axis is clean. It
also reports `floor_U`, the minimum uridines the encoded protein admits: Phe,
Tyr, Cys, Trp and Ile have no U-free codon, so there is a hard floor, and an
arm sitting on it cannot be pushed lower without changing the protein. An
optimiser that "stops responding" to more weight has simply reached it.

---

## LNP formulation

N/P — moles of ionisable nitrogen over moles of RNA backbone phosphate — is the
quantity that ties mRNA to lipid. One nucleotide is one phosphate, so it follows
from sequence. RNA mass is computed **from the actual base composition**, with
an optional per-uridine offset for m1Ψ, rather than a 330 Da/nt rule of thumb
that drifts by a few percent on a uridine-depleted transcript.

```bash
txv formulate --fasta out/TXV-001.mrna.fasta --mrna-ug 100 --np 6 \
  --stocks '{"SM-102":10,"DSPC":10,"cholesterol":20,"DMG-PEG2000":10}' \
  --doses 0.020,0.025 --mg-per-kg 0.5 --concentration 0.1
```

Gives per-lipid µmol, µg and stock pipetting volumes, the ethanol:aqueous mixing
volumes at your flow rate ratio, and a per-subject dose table. `txv compositions`
lists the shipped lipid compositions.

This is stoichiometry and dilution arithmetic only. It does **not** predict
particle size, PDI or encapsulation efficiency — those are measured, and
`FormulationResult` has fields to record them against the batch.

---

## Benchling

Credentials come from the environment, never from a design file:

```bash
export BENCHLING_TENANT=https://yourcompany.benchling.com
export BENCHLING_API_KEY=...          # or BENCHLING_CLIENT_ID + _SECRET
export BENCHLING_FOLDER_ID=lib_...
export BENCHLING_SCHEMA_ID=ts_...     # optional
export BENCHLING_REGISTRY_ID=src_...  # optional; omit to create unregistered
```

```bash
txv design antigens.csv --name TXV-001 --out out/ --dry-run-register  # review
txv design antigens.csv --name TXV-001 --out out/ --register          # send
```

Two deliberate behaviours:

* **Dry run is the default** whenever credentials are absent, so an accidental
  run in CI produces files rather than registry entries.
* **Registration is refused on a QC failure** (`guard_qc`). A failing construct
  in the registry is worse than no construct at all.

Annotations pass through unchanged — Benchling's coordinates are 0-based,
end-exclusive with an explicit strand, which is this package's internal
convention. Each module gets its own coloured feature rather than one opaque CDS
block. Schema fields are populated with antigen count and names, transcript and
ORF length, GC, uridine fraction, CAI, linker, placeholder list and QC status.

---

## Tests

```bash
python -m pytest        # 218 tests
```

The suite encodes the invariants that matter: codon optimisation never changes
the encoded protein; reordering never increases junctional risk; junction
scanning never flags a peptide contained in one bead; features tile the
template; GenBank coordinates are 1-based inclusive and the ORIGIN block
reproduces the sequence; N/P is exact; signal peptides are N-terminal and
`GYQTI` ends the protein in every panel member; the 5' UTR plus Kozak
reconstitutes the validated BNT162b2 leader exactly; multiplex arms are nested
and share their anchor; uridine arms encode a byte-identical protein.

Sources for the UTR sequences: the [assembled BNT162b2 sequence](https://github.com/NAalytics/Assemblies-of-putative-SARS-CoV2-spike-encoding-mRNA-sequences-for-vaccines-BNT-162b2-and-mRNA-1273),
cross-read against [Xia 2021, *Vaccines*](https://pmc.ncbi.nlm.nih.gov/articles/PMC8310186/)
for the α-globin / AES / mtRNR1 element identities.
