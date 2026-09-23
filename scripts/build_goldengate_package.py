"""Rebuild the BsaI Golden Gate order package.

This package existed on disk with no build step, which meant it silently went
stale the moment the encodings changed. It is regenerated here from the same
modules as every other package so the three cannot disagree.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from txv.goldengate import make_cassette_part, verify_part  # noqa: E402
from txv.order import assess_synthesis  # noqa: E402
from txv.pvax1_ag import PANEL, build_panel_construct  # noqa: E402
from txv.seqops import gc_fraction  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "orders" / "pvax1_ag_panel_goldengate"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    fasta, rows, problems = [], [], []

    for panel in PANEL:
        part = make_cassette_part(panel.name, build_panel_construct(panel.name).orf)
        problems += [(panel.name, p) for p in verify_part(part)]
        seq = part.sequence
        fasta.append(f">{panel.name}|goldengate|{len(seq)}bp gc={gc_fraction(seq):.3f}")
        fasta += [seq[i:i + 60] for i in range(0, len(seq), 60)]
        rows.append({
            "name": panel.name,
            "length_bp": len(seq),
            "gc": f"{gc_fraction(seq):.4f}",
            "assembly": "BsaI Golden Gate",
            "synthesis_flags": "; ".join(assess_synthesis(seq).flags) or "none",
            "sequence": seq,
        })

    (OUT / "gblocks.fasta").write_text("\n".join(fasta) + "\n")
    with (OUT / "ordering_table.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(f"  {row['name']:<3}{row['length_bp']:>6} bp  GC {float(row['gc']):.1%}"
              f"  {row['synthesis_flags']}")
    if problems:
        print("\nPROBLEMS:")
        for name, p in problems:
            print(f"  [{name}] {p}")
        return 1
    print(f"\nwrote gblocks.fasta and ordering_table.csv to {OUT}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
