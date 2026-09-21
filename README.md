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

## Read this before you use the output

### Placeholder parts are not real sequences

`UTR5_placeholder`, `UTR3_placeholder` and `MITD_placeholder` are **stand-ins**
so the assembler runs end to end. They are not validated UTRs. Every part
carries a `Provenance` tag (`CANONICAL` / `LITERATURE` / `PLACEHOLDER`), any
construct containing a placeholder reports it, and QC **fails** on one unless
you pass `allow_placeholders=True` (or drop `--strict`).

Supply your own with a JSON parts file:

```bash
txv design antigens.csv --name X --parts my_parts.json --utr5 UTR5_house --strict
```

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
| `pvax1_ag` | The audited P0–P8 module library and panel |

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
python -m pytest        # 153 tests
```

The suite encodes the invariants that matter: codon optimisation never changes
the encoded protein; reordering never increases junctional risk; junction
scanning never flags a peptide contained in one bead; features tile the
template; GenBank coordinates are 1-based inclusive and the ORIGIN block
reproduces the sequence; N/P is exact; signal peptides are N-terminal and
`GYQTI` ends the protein in every panel member.
