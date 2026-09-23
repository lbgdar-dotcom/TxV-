"""In-silico validation of the nine transcription units.

Run after any design change:  python scripts/validate_in_silico.py

This is a different kind of check from ``audit/``. The auditors verify the
molecule is the one that was designed. This verifies the designed molecule is
one that should work -- folding, codon adaptation, and the hydrophobicity the
routing depends on. It writes IN_SILICO_REPORT.md next to the order package and
exits non-zero if a construct fails a threshold.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from txv.insilico import (  # noqa: E402
    CAP_WINDOW, TM_MIN_LENGTH, cai, cpg_stats, fold_metrics, signal_peptide_call,
    tm_segments,
)
from txv.ivt_unit import LEADER, make_unit  # noqa: E402
from txv.pvax1_ag import MODULES, PANEL, build_panel_construct  # noqa: E402
from txv.seqops import gc_fraction  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "orders" / "pjet12_ivt_units"

# --- thresholds, and why each one is where it is ---------------------------
#: Mean unpaired probability over the AUG window. Below this the start codon is
#: sequestered often enough to matter for scanning. Not a literature constant --
#: it is a panel-consistency bar, set so that a construct which folds its own
#: start codon away stands out from its eight siblings.
MIN_START_ACCESSIBILITY = 0.30
#: MFE of the first CAP_WINDOW nt on its own. A cap-proximal hairpin this
#: stable impedes 43S loading.
MIN_CAP_MFE = -25.0
#: CAI against human usage. The panel should sit well above the ~0.7 that
#: unoptimised human coding sequence averages.
MIN_CAI = 0.75
#: CpG observed/expected. Mammalian transcripts run low; above this the
#: transcript is a better ZAP substrate than it needs to be.
MAX_CPG_OE = 0.80
#: Spread allowed between the best and worst construct on structural measures,
#: as a fraction. The panel compares routes; it should not also be comparing
#: folding energies.
MAX_MFE_PER_NT_SPREAD = 0.25


def main() -> int:
    problems: list[tuple[str, str]] = []
    rows = []

    print("folding nine transcripts (ViennaRNA)...\n")
    for panel in PANEL:
        construct = build_panel_construct(panel.name)
        unit = make_unit(panel.name, construct.orf)
        transcript = unit.transcript
        orf_start = transcript.index(LEADER) + len(LEADER)

        fold = fold_metrics(transcript, orf_start)
        score = cai(construct.orf)
        cpg_n, cpg_oe = cpg_stats(transcript)

        rows.append({
            "construct": panel.name,
            "transcript_nt": len(transcript),
            "gc": f"{gc_fraction(transcript):.4f}",
            "mfe_kcal": f"{fold.mfe:.1f}",
            "mfe_per_nt": f"{fold.mfe_per_nt:.4f}",
            "ensemble_free_energy": f"{fold.ensemble_free_energy:.1f}",
            "ensemble_diversity": f"{fold.ensemble_diversity:.1f}",
            "start_accessibility": f"{fold.start_accessibility:.3f}",
            "cap_proximal_mfe": f"{fold.cap_proximal_mfe:.1f}",
            "cai": f"{score:.4f}",
            "cpg_count": cpg_n,
            "cpg_obs_exp": f"{cpg_oe:.3f}",
        })

        if fold.start_accessibility < MIN_START_ACCESSIBILITY:
            problems.append((panel.name, f"start codon accessibility "
                             f"{fold.start_accessibility:.2f} < "
                             f"{MIN_START_ACCESSIBILITY}"))
        if fold.cap_proximal_mfe < MIN_CAP_MFE:
            problems.append((panel.name, f"cap-proximal {CAP_WINDOW} nt fold to "
                             f"{fold.cap_proximal_mfe:.1f} kcal/mol, below "
                             f"{MIN_CAP_MFE}"))
        if score < MIN_CAI:
            problems.append((panel.name, f"CAI {score:.3f} < {MIN_CAI}"))
        if cpg_oe > MAX_CPG_OE:
            problems.append((panel.name, f"CpG o/e {cpg_oe:.2f} > {MAX_CPG_OE}"))

        print(f"  {panel.name}  {len(transcript):>4} nt  "
              f"MFE {fold.mfe:>8.1f}  ({fold.mfe_per_nt:>6.3f}/nt)  "
              f"AUG access {fold.start_accessibility:>5.2f}  "
              f"cap {fold.cap_proximal_mfe:>6.1f}  "
              f"CAI {score:.3f}  CpG o/e {cpg_oe:.2f}")

    # -- panel consistency ---------------------------------------------------
    per_nt = [abs(float(r["mfe_per_nt"])) for r in rows]
    spread = (max(per_nt) - min(per_nt)) / max(per_nt)
    print(f"\n  MFE/nt spread across the panel: {spread:.1%}")
    if spread > MAX_MFE_PER_NT_SPREAD:
        problems.append(("panel", f"MFE/nt varies {spread:.0%} across the panel, "
                         f"above {MAX_MFE_PER_NT_SPREAD:.0%}; structure is a "
                         "second variable alongside routing"))

    # -- the routing modules have to be what they claim ----------------------
    print("\n  routing modules:")
    for name in ("CTLA4_SP", "LAMP1_SP"):
        call = signal_peptide_call(MODULES[name])
        state = "ok" if call.ok else "; ".join(call.problems)
        print(f"    {name:<10} h-region {call.h_region_length:>2} aa "
              f"(max KD {call.h_region_max:.2f}), n-charge "
              f"{call.n_region_charge:+.0f}, -3/-1 = "
              f"{call.minus_three}/{call.minus_one}  -> {state}")
        for problem in call.problems:
            problems.append((name, problem))

    for name in ("CTLA4_TMT", "LAMP1_TMT"):
        segments = tm_segments(MODULES[name])
        spans = ", ".join(f"{a + 1}-{b + 1} (peak KD {p:.2f})"
                          for a, b, p in segments)
        long_enough = [s for s in segments if s[1] - s[0] + 1 >= TM_MIN_LENGTH]
        print(f"    {name:<10} hydrophobic span(s): {spans or 'none'}")
        if not long_enough:
            problems.append((name, f"no hydrophobic segment >= {TM_MIN_LENGTH} "
                             "residues; this module cannot span a bilayer"))

    # LAMP1's sorting motif only works at the very end of the protein.
    for panel in PANEL:
        protein = build_panel_construct(panel.name).orf_protein
        if "LAMP1_TMT" in panel.modules and not protein.endswith("GYQTI"):
            problems.append((panel.name, "GYQTI is not the C-terminus"))

    # -- report --------------------------------------------------------------
    with (OUT / "in_silico_metrics.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# In-silico validation — P0–P8 transcription units",
        "",
        "Generated by `scripts/validate_in_silico.py`. Folding by ViennaRNA "
        "(Lorenz et al. 2011); CAI per Sharp & Li 1987; hydropathy per Kyte & "
        "Doolittle 1982.",
        "",
        "This answers a different question from `audit/`. The auditors check "
        "the molecule is the one designed. This checks the designed molecule "
        "should behave: that the start codon is reachable, that codon "
        "adaptation and structure are consistent across the panel, and that "
        "the modules the routing depends on really are hydrophobic enough to "
        "do their job.",
        "",
        "## Per construct",
        "",
        "| construct | nt | GC | MFE (kcal/mol) | MFE/nt | AUG access. | "
        "cap-proximal MFE | CAI | CpG o/e |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['construct']} | {r['transcript_nt']} | "
            f"{float(r['gc']):.1%} | {r['mfe_kcal']} | {r['mfe_per_nt']} | "
            f"{r['start_accessibility']} | {r['cap_proximal_mfe']} | "
            f"{r['cai']} | {r['cpg_obs_exp']} |")
    lines += [
        "",
        "**AUG access.** is the mean unpaired probability over a "
        "15-nt window centred on the start codon: 1.0 is a fully single-"
        "stranded start codon, 0.0 a fully base-paired one.",
        "",
        "**cap-proximal MFE** folds the first 70 nt alone, because a "
        "transcript can be relaxed overall and still carry a hairpin on the "
        "cap where the 43S complex has to load.",
        "",
        "## Thresholds",
        "",
        "| measure | bar | basis |",
        "|---|---|---|",
        f"| AUG accessibility | ≥ {MIN_START_ACCESSIBILITY} | panel-"
        "consistency bar, not a literature constant |",
        f"| cap-proximal MFE | ≥ {MIN_CAP_MFE} kcal/mol | a hairpin more "
        "stable than this impedes 43S loading |",
        f"| CAI | ≥ {MIN_CAI} | unoptimised human CDS averages ~0.7 |",
        f"| CpG obs/exp | ≤ {MAX_CPG_OE} | ZAP substrate and innate sensing |",
        f"| MFE/nt spread | ≤ {MAX_MFE_PER_NT_SPREAD:.0%} | structure must not "
        "be a second variable alongside route |",
        "",
        "## What this does not cover",
        "",
        "- **No MHC binding prediction.** The epitopes are literature-"
        "characterised (SIINFEKL/H-2K(b); I-A(b) core); a predicted affinity "
        "would add a number, not information. The auditors check the epitopes "
        "are intact and in frame.",
        "- **No trained signal-peptide or topology predictor.** SignalP and "
        "TMHMM are licence-encumbered and not reproducible from this repo. The "
        "hydropathy architecture check is a weaker stand-in: it will catch a "
        "module that has stopped being a signal peptide, but it does not "
        "predict the cleavage site.",
        "- **No prediction of actual presentation.** Nothing in silico "
        "separates the routes; that is what the panel is for.",
        "",
    ]
    (OUT / "IN_SILICO_REPORT.md").write_text("\n".join(lines) + "\n")

    print(f"\nwrote IN_SILICO_REPORT.md and in_silico_metrics.csv to {OUT}/")
    if problems:
        print(f"\n{len(problems)} FINDING(S):")
        for name, problem in problems:
            print(f"  [{name}] {problem}")
        return 1
    print("\nno findings: all nine pass every in-silico threshold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
