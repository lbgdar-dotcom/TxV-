"""Derive the 5'-terminal codons of the signal-peptide modules.

Why this is a separate step. Local mRNA structure over the start codon is set
by the first ~40 nt of the ORF, and those nt belong to whichever module starts
the construct. The panel has three distinct 5' ends -- CTLA4_SP (P0-P4),
LAMP1_SP (P6) and MA+HA (P5, P7, P8) -- so it has three distinct start-codon
accessibilities whether anyone chose them or not. Left alone they came out at
0.37, 0.50 and 0.52, which makes translation initiation a second variable
running alongside the route the panel is supposed to be testing.

The fix is to re-encode the first few codons of the two signal peptides -- not
to maximise accessibility, but to bring all three groups together. MA+HA is
left alone deliberately: HA is carried by all nine constructs, so tuning it for
the three that start with it would perturb the six that do not.

Accessibility is computed by local folding (RNAplfold; Bernhart et al. 2006,
*Algorithms Mol Biol* 1:3) with a bounded base-pair span. Global MFE folding
was tried first and rejected: over a 700-nt transcript it reported five
different accessibilities for the five constructs whose 5' ends are byte-
identical, which is long-range pairing noise, not signal.

Run this to re-derive; paste the result into FIVE_PRIME_CODONS in pvax1_ag.py.
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import RNA  # noqa: E402

from txv.codon import HUMAN_CODON_USAGE, CodonOptimizer  # noqa: E402
from txv.ivt_unit import LEADER  # noqa: E402
from txv.pvax1_ag import (  # noqa: E402
    MODULES, PANEL, build_panel_construct, canonical_encodings,
    panel_optimizer_config,
)
from txv.seqops import AA_TO_CODONS, gc_fraction  # noqa: E402

#: Codons at the start of the ORF that are free to be re-chosen. Codon 1 is
#: ATG and codon 2 is the pinned Kozak alanine, so tuning starts at codon 3.
FIRST_FREE_CODON = 2
N_FREE = 6

#: RNAplfold parameters. A 40-nt maximum span keeps the measure local, which is
#: both what a scanning 43S sees and the regime where the folding model is
#: trustworthy.
PLFOLD_WINDOW, PLFOLD_MAX_SPAN = 80, 40
#: GC band the first 120 nt must stay inside.
GC_BAND = (0.40, 0.66)
#: How close to the fixed group's accessibility is close enough.
#: Differences below this are well inside the folding model's error.
ACCESS_TOLERANCE = 0.05
#: Resolution at which two accessibilities count as different.
ACCESS_RESOLUTION = 0.02
#: Window centred on the AUG over which mean unpaired probability is taken.
AUG_WINDOW = 15


def accessibility(orf: str) -> float:
    """Mean per-base unpaired probability across the start-codon window."""
    seq = (LEADER + orf)[:160].replace("T", "U")
    start = len(LEADER)
    up = RNA.pfl_fold_up(seq, 1, PLFOLD_WINDOW, PLFOLD_MAX_SPAN)
    flank = (AUG_WINDOW - 3) // 2
    lo, hi = start - flank, start + 3 + flank
    return sum(up[i][1] for i in range(lo, hi)) / (hi - lo)


def search(module: str, target: float, optimizer) -> list[tuple]:
    canon = canonical_encodings()
    current = canon[module]
    context = ""
    for panel in PANEL:
        if panel.modules[0] == module:
            context = build_panel_construct(panel.name).orf
            break

    head = FIRST_FREE_CODON * 3
    aas = MODULES[module][FIRST_FREE_CODON:FIRST_FREE_CODON + N_FREE]
    options = [[c for c in AA_TO_CODONS[a] if HUMAN_CODON_USAGE[c] >= 0.10]
               for a in aas]

    scored = []
    rejected = Counter()
    for combo in product(*options):
        if any(combo[i] == combo[i + 1] for i in range(len(combo) - 1)):
            continue
        variant = current[:head] + "".join(combo) + current[head + N_FREE * 3:]
        candidate = variant + context[len(current):]
        window = candidate[:120]
        if optimizer.find_forbidden(window):
            rejected["forbidden motif"] += 1
            continue
        if re.search(r"(A|C|G|T)\1{3,}", window):
            rejected["homopolymer >=4"] += 1
            continue
        if not GC_BAND[0] <= gc_fraction(window) <= GC_BAND[1]:
            rejected[f"GC outside {GC_BAND}"] += 1
            continue
        access = accessibility(candidate)
        if abs(access - target) > ACCESS_TOLERANCE:
            rejected["accessibility outside the band"] += 1
            continue
        usage = sum(HUMAN_CODON_USAGE[c] for c in combo) / len(combo)
        cpg = sum(1 for i in range(len(window) - 1) if window[i:i + 2] == "CG")
        # Inside the accessibility band all candidates are equivalent for the
        # purpose this exists to serve, so the tie-break is what else we care
        # about: CpG is a ZAP substrate, and rare codons slow elongation.
        # Matching accessibility to three decimal places at the cost of a
        # CpG-dense, rare-codon 5' end would be trading a real liability for a
        # number that is below the resolution of the model producing it.
        # Accessibility still counts inside the band, but only at the
        # resolution the model can actually resolve -- otherwise a candidate
        # drifts to the band edge for a third decimal place of CpG.
        bucket = round(abs(access - target) / ACCESS_RESOLUTION)
        scored.append((cpg, bucket, -usage, access, usage, "".join(combo)))
    scored.sort()
    return scored, rejected


def main() -> int:
    optimizer = CodonOptimizer(panel_optimizer_config())
    canon = canonical_encodings()

    # The group that is not being tuned sets the target the others match.
    fixed = ""
    for panel in PANEL:
        if panel.modules[0] == "MA":
            fixed = build_panel_construct(panel.name).orf
            break
    target = accessibility(fixed)
    print(f"MA+HA start (P5/P7/P8) is fixed at accessibility {target:.3f}")
    print("-- tuning the two signal peptides to match it --\n")

    for module in ("CTLA4_SP", "LAMP1_SP"):
        before = ""
        for panel in PANEL:
            if panel.modules[0] == module:
                before = build_panel_construct(panel.name).orf
                break
        scored, rejected = search(module, target, optimizer)
        if not scored:
            print(f"{module}: no candidate passed the constraints")
            for reason, n in rejected.most_common():
                print(f"    {n:>6} rejected: {reason}")
            return 1
        _, _, _, access, usage, codons = scored[0]
        head = FIRST_FREE_CODON * 3
        print(f"{module}")
        print(f"  current  {canon[module][head:head + N_FREE * 3]}  "
              f"accessibility {accessibility(before):.3f}")
        print(f"  chosen   {codons}  accessibility {access:.3f}  "
              f"mean usage {usage:.2f}  CpG {scored[0][0]} in the first 120 nt")
        print(f"  {len(scored)} candidates inside the accessibility band "
              f"[{target - ACCESS_TOLERANCE:.2f}, {target + ACCESS_TOLERANCE:.2f}]; "
              f"CpG among them {min(x[0] for x in scored)}-"
              f"{max(x[0] for x in scored)}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
