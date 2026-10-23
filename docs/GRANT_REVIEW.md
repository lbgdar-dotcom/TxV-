# Committee review — TxV antigen-routing panel

**Scope limitation, stated first.** No proposal was supplied. Section 18 of the
review request was returned with every field unfilled: no funding call, no
duration, no budget, no institutions, no infrastructure, no preliminary data,
no collaborators, no proposal text. What follows reviews the **project as it
exists in this repository and conversation**, not a document. Nothing has been
invented to fill the gaps: where a judgement needs information that does not
exist, it says so.

**Second limitation.** Literature verification from this environment is
severely restricted — PMC, NCBI, ScienceDirect, Nature and Google Patents do
not resolve. Only web search is available. Citations below are marked
**[verified]** where a search returned the specific finding attributed to it,
and **[recalled]** where it comes from training and has not been confirmed.
Phase II's standard — verified DOIs for every central claim — is **not met**,
and that gap must be closed by hand before any submission.

**Third.** This is an AI-generated assessment and does not substitute for
review by qualified human experts.

---

## Deliverable A — Executive assessment

The project asks whether the subcellular destination of an mRNA-encoded
antigen, specified entirely within the transcript, controls the class of MHC
presentation it receives. Nine constructs (P0–P8) hold every module constant
except the routing elements, and two antigens carry the only epitopes in mouse
for which validated pMHC-specific antibodies exist.

**The construct engineering is the strongest part and is genuinely above
normal academic standard.** Module canonicalisation — one DNA sequence per
module across the whole panel — means routing is not confounded with codon
usage, which most comparable studies do not control. The positional logic
(SRP reads an N-terminal signal only; a C-degron needs a free C-terminus;
LAMP1's `GYQTI` must be terminal) is correct and documented at the point of
use. The audit and regression infrastructure is unusual and real.

**The immunological design is largely absent, and one gap is structural.**
The readout cannot test the stated hypothesis: antigen A is measurable only in
class I and antigen B only in class II, so no antigen in the panel can be
observed moving *between* classes. There is no statistical design, no
kinetics, no dose–response, no pathway-blockade control, and no treatment of
cross-presentation — which is precisely the mechanism that would blur the
class I/class II distinction the project rests on.

**It is not, as built, a cancer project.** There is no tumour, no tumour
antigen, no efficacy endpoint and no tumour microenvironment anywhere in the
repository. SIINFEKL and Eα are model epitopes. Submitted to a cancer call in
this form, that is a disqualifying mismatch between title and content.

**Funding readiness: not competitive in its current form** — not because the
science is weak, but because roughly two-thirds of what a proposal must
contain has not been written, and two design problems (readout, normalisation)
would be identified by any competent immunology reviewer.

---

## Deliverable B — Prioritised scientific audit

| # | Problem | Expert | Scientific consequence | Correction | Impact |
|---|---|---|---|---|---|
| **C1** | No tumour, tumour antigen or efficacy endpoint anywhere, in a project titled "cancer vaccines" | 5, 2 | Title/content mismatch; no translational claim is supportable | Either reframe as mechanistic immunology, or add a tumour arm. **B16-F10-OVA or E.G7-OVA** are C57BL/6-compatible and OVA-expressing, so SIINFEKL becomes a true rejection antigen | Decides which call is even eligible |
| **C2** | Readout cannot observe class switching: A readable only in class I (25-D1.16, OT-I), B only in class II (Y-Ae, TEa) | 4 | The central hypothesis is **not testable** by the current panel | Add **OVA323–339** alongside SIINFEKL in one cassette. That antigen is then readable in class I (25-D1.16/OT-I) *and* class II (OT-II, I-A(b)/OVA323–339 tetramer) | Highest-value single change; cheap |
| **C3** | Expression normalisation is confounded with the independent variable | 8, 13 | Routing changes protein half-life *and* HA accessibility (lumenal in P0–P4/P6, cytosolic in P5/P7/P8). "Presentation per unit protein" is circular | Normalise to **delivered mRNA** (RT-qPCR) plus a co-formulated route-invariant reporter; report **absolute pMHC/cell** via calibration beads, not ratios | Determines whether any number means anything |
| **M4** | LNP is undesigned. `lnp.py` is stoichiometry only, by its own docstring | 3, 6, 7 | LNP innate activation alters MHC and costimulation directly — DC preactivation *impairs* class II presentation **[verified]** — so batch variation moves the readout for reasons unrelated to routing | Fix one composition panel-wide; characterise every batch (size, PDI, EE, apparent pKa); include an **empty-LNP arm in every experiment** | Without this the whole panel is uncontrolled |
| **M5** | P7/P8 confound route with ORF position | 13, 4 | For antigen A: P7 = cytosolic AND upstream; P8 = lysosomal AND downstream. 2A skipping is incomplete, so downstream products are systematically scarcer | Either add the position control (costs SP-first, which inhibits 2A separation) or restrict the claim and normalise per product | Already recorded in the panel |
| **M6** | E5 degron inactive in P7/P8 (no free C-terminus behind the 2A remnant) | 4, 8 | P7/P8 are not degron-positive; P5 is the comparator | Already documented; must not be silently treated as E5-containing | Interpretation only |
| **M7** | No statistical design at all | 13 | No n, no power, no effect size, no randomisation or blinding, no multiplicity control across 9 constructs × multiple readouts × timepoints | Pre-specify primary endpoint, power on a stated effect size, randomise, blind the flow analysis, pre-register the comparison set | Mandatory for any serious call |
| **M8** | Single timepoint, no dose–response | 4, 8 | pMHC display is transient and the routes differ in speed (direct TGN→lysosome vs via plasma membrane). Two routes can look identical, or opposite, purely by sampling time | Time course (≥4 points) and ≥3 doses, on at least the key constructs | Without it, results are uninterpretable |
| **M9** | Cross-presentation unaddressed | 4 | In DCs, endosomal/lysosomal antigen also feeds class I. "Lysosomal → class II" is therefore not a clean dichotomy; recent work reports mRNA-LNP CD8 priming that does *not* use the classical cross-presentation pathway **[verified]** | Address explicitly; include cross-presentation-pathway controls | Affects the core interpretation |
| **Mo10** | Signal-peptide cleavage and topology assumed, not measured | 8 | If the CTLA-4 SP cleaves poorly on these chimeras, the routing premise fails silently | Verify by mobility shift ± tunicamycin/EndoH, and surface vs total staining | Cheap, foundational |
| **Mo11** | No pathway-blockade controls | 4 | Correlation with construct design is not mechanism | Bafilomycin/chloroquine (endosomal acidification), proteasome inhibition, TAP- or H2-DM-deficient comparators | Converts association into mechanism |
| **Mo12** | In-silico validation measures the mRNA, not the hypothesis | 1, 13 | Folding/CAI metrics are good QC but no reviewer will credit them as evidence for routing. Folding is computed on unmodified RNA while the product is m1Ψ-substituted | Keep as QC; do not present as mechanistic support | Presentation only |

---

## Deliverable C — Novelty assessment

| Question | Established | Evidence | Contradictory | Remaining gap |
|---|---|---|---|---|
| Does LAMP1 targeting enhance class II presentation? | **Yes — since 1995** | PNAS 1995;92:11671 (LAMP-1 sorting signal fused to HPV-16 E7) **[verified]**; multiple DNA-vaccine studies; granted US patents on LAMP chimeras **[verified]** | A malaria DNA-vaccine study found LAMP fusion gave **no** antibody enhancement and only marginal CD4 proliferation **[verified]** | **Effect size is antigen- and model-dependent and has never been quantified against a matched class I readout** |
| Does mRNA-encoded antigen default to class I? | Yes, conceptually | Endogenous synthesis favours class I **[verified]** | mRNA-LNP CD8 priming reported to require migratory cDC2 but *not* classical cross-presentation **[verified]** | The actual routing of LNP-delivered antigen in DCs is unsettled |
| Has the class I ↔ class II exchange rate been measured directly? | **No evidence found** | — | — | Search returned nothing on direct pMHC quantification across routing variants |

**Verdict on novelty.** "Route antigen to the lysosome and get more CD4 help"
is **thirty years old and patented**. Any proposal claiming that as its
innovation will be rejected on positioning alone.

The defensible novelty is narrower and better:

- **Conceptual** — not *whether* routing works but **what it costs**. Is class I
  and class II presentation a zero-sum exchange, and what is the exchange rate?
- **Methodological** — a module-canonical panel in which routing is the only
  variable, read by **direct pMHC quantification** rather than T-cell proxies.
- **Technological** — dual-route output from a single transcript (P2A), with the
  honest confound stated.

Weak-to-absent: translational novelty (no tumour), chemical novelty (no new
lipid), mechanistic novelty if framed as "LAMP1 works".

---

## Deliverable D — Revised central hypothesis and aims

**Current implicit hypothesis:** routing elements in the transcript control
MHC class. *Problem:* known to be true in part, and untestable with the
current readout.

**Proposed hypothesis:** *MHC class I and class II presentation of a single
mRNA-encoded antigen are competitive, and the subcellular routing encoded in
the transcript sets the exchange rate between them in a quantifiable,
dose- and time-dependent manner.*

Falsifiable, not already answered, and exactly what this panel can measure —
**once one antigen is readable in both classes.**

### Aim 1 — Establish a dual-readable antigen and verify the constructs behave as designed
Add OVA323–339 beside SIINFEKL. Verify SP cleavage, topology, and subcellular
localisation (confocal colocalisation with LAMP1, calnexin, surface stain).
*Endpoint:* each construct localises as designed. **Decision gate:** constructs
that mislocalise are fixed or dropped before any immunology.

### Aim 2 — Quantify the exchange rate
Absolute pMHC molecules/cell (25-D1.16 for Kb/SIINFEKL; Y-Ae for I-Ab/Eα) across
the panel, over ≥4 timepoints and ≥3 doses, normalised to delivered mRNA.
*Primary endpoint:* class I pMHC vs class II pMHC per unit delivered transcript.
*Controls:* empty LNP, untransfected, route-invariant co-reporter.

### Aim 3 — Establish mechanism
Pathway blockade (endosomal acidification, proteasome) and deficient-cell
comparators to show the routes use the pathways claimed, and to quantify the
cross-presentation contribution.

### Aim 4 (only if a cancer call) — Does the exchange rate predict tumour control?
B16-F10-OVA or E.G7-OVA, therapeutic setting. *Endpoint:* tumour growth and
survival against the measured pMHC ratio. This is what converts a mechanistic
result into a translational claim.

---

## Deliverable E — Execution and budget

**Cannot be completed.** No duration, staffing, infrastructure, equipment
access, institution or budget ceiling was supplied, and inventing figures would
breach the review request's own integrity rules. What can be said:

- The **critical path** is Aim 1. Nothing downstream is interpretable until
  localisation is verified, so it must not be run in parallel with Aim 2.
- The **largest cost risk** is the mouse work in Aim 4, and the largest
  *hidden* cost is LNP characterisation, which is currently unbudgeted because
  it is undesigned.
- **Already in hand** (real, not assumed): the construct library, the audit
  suite, the ordering packages, and the in-silico QC.
- **Requires quotation:** gBlock synthesis (9–11 fragments), LNP lipids and
  microfluidics access, 25-D1.16 and Y-Ae antibodies, tetramers (NIH Tetramer
  Core), TCR-transgenic mice (OT-I, OT-II, TEa), tumour lines, animal facility.

---

## Deliverable F — Revised proposal

**Not produced.** A revision requires an original; none was supplied. Writing a
full proposal here would mean inventing aims, budget, institution and
preliminary data — explicitly forbidden by the request. The aims in Deliverable
D are the honest starting skeleton; the remaining sections need input that does
not exist in this repository.

---

## Deliverable G — Scoring

Diagnostic only; not a prediction of funding outcome.

| Criterion | Weight | Score | Basis |
|---|---|---|---|
| Significance & clinical relevance | 15% | **2.5** | Real question; no tumour or clinical link as built |
| Originality | 15% | **3.0** | Low as "LAMP1 targeting"; good as "exchange rate", if repositioned |
| Hypothesis & rationale | 10% | **3.5** | Clear and mechanistic; needs the revision in D |
| Experimental design & rigour | 20% | **2.0** | Construct design excellent; immunological design largely absent |
| Translational potential | 10% | **1.5** | No tumour, no clinical antigen, no path |
| Feasibility & resources | 10% | — | **No information supplied** |
| Project management | 5% | — | **No information supplied** |
| Budget justification | 5% | — | **No information supplied** |
| Writing & compliance | 10% | — | **No document supplied** |

**Weighted score across the 70% that is scoreable: ≈2.5/5.** The unscoreable
30% is unscoreable because the material does not exist, which is itself the
finding.

---

## Deliverable H — The five highest-priority actions

1. **Decide what this project is.** Mechanistic immunology, or cancer vaccine
   with a tumour arm. Everything else — call, aims, budget, framing — follows
   from that one decision, and it cannot be deferred.
2. **Make one antigen readable in both classes** (add OVA323–339). Without it
   the central hypothesis cannot be tested at all.
3. **Fix the normalisation.** Normalise to delivered mRNA and report absolute
   pMHC per cell. The current denominator moves with the independent variable.
4. **Design the LNP, or stop claiming it.** One fixed composition, characterised
   per batch, with an empty-LNP arm in every experiment.
5. **Write a statistical plan.** Primary endpoint, powered n, randomisation,
   blinded analysis, pre-specified comparisons.

**Committee recommendation: major revision required before this is a proposal
at all.** The underlying science is worth doing and the construct work is a
genuine asset; what is missing is the experiment, not the reagent.
