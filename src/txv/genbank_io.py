"""Minimal GenBank reader.

Enough to read a plasmid map: the LOCUS line, simple feature locations
(``N..M``, ``complement(N..M)``, ``join(...)``), qualifiers, and the ORIGIN
block. It is not a general parser -- it does not handle remote accessions or
fuzzy endpoints -- but it round-trips the records this package writes and the
vector maps that come out of SnapGene, ApE and Benchling.

Coordinates are converted from GenBank's 1-based inclusive convention to this
package's 0-based half-open one on the way in, so a parsed record indexes the
sequence directly.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from .seqops import clean


@dataclass
class GenBankFeature:
    key: str
    start: int          # 0-based, inclusive
    end: int            # 0-based, exclusive
    strand: int         # +1 or -1
    qualifiers: dict[str, str] = field(default_factory=dict)

    @property
    def label(self) -> str:
        for name in ("label", "gene", "product", "note"):
            if name in self.qualifiers:
                return self.qualifiers[name]
        return self.key

    def __len__(self) -> int:
        return self.end - self.start

    def slice(self, sequence: str) -> str:
        from .seqops import revcomp

        body = sequence[self.start : self.end]
        return revcomp(body) if self.strand < 0 else body


@dataclass
class GenBankRecord:
    name: str
    sequence: str
    features: list[GenBankFeature] = field(default_factory=list)
    is_circular: bool = False
    definition: str = ""

    def __len__(self) -> int:
        return len(self.sequence)

    def find(self, label: str) -> GenBankFeature:
        """The single feature with this label; raises if absent or ambiguous."""
        hits = [f for f in self.features if f.label == label]
        if len(hits) != 1:
            raise KeyError(
                f"{len(hits)} feature(s) labelled {label!r} "
                f"(have: {sorted({f.label for f in self.features})})"
            )
        return hits[0]

    def find_all(self, label: str) -> list[GenBankFeature]:
        return [f for f in self.features if f.label == label]


_LOCATION = re.compile(r"(\d+)\.\.(\d+)")


def parse_genbank(text: str) -> GenBankRecord:
    """Parse a single GenBank record."""
    lines = text.splitlines()

    name, circular, definition = "", False, ""
    for line in lines:
        if line.startswith("LOCUS"):
            parts = line.split()
            name = parts[1] if len(parts) > 1 else ""
            circular = "circular" in line.lower()
        elif line.startswith("DEFINITION"):
            definition = line[len("DEFINITION"):].strip()
        if line.startswith("FEATURES"):
            break

    # -- sequence ----------------------------------------------------------
    try:
        origin = text.split("\nORIGIN", 1)[1]
    except IndexError:
        raise ValueError("record has no ORIGIN block") from None
    sequence = clean("".join(re.findall(r"[acgtnACGTN]+",
                                        origin.split("//")[0].replace("ORIGIN", ""))))

    # -- features ----------------------------------------------------------
    features: list[GenBankFeature] = []
    in_features = False
    current: GenBankFeature | None = None
    pending_key: str | None = None
    pending_value: list[str] = []

    def flush_qualifier() -> None:
        if current is not None and pending_key is not None:
            current.qualifiers[pending_key] = (
                "".join(pending_value).strip().strip('"')
            )

    for line in lines:
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if not in_features:
            continue
        if line.startswith("ORIGIN") or line.startswith("//"):
            flush_qualifier()
            break
        if not line.startswith(" "):
            flush_qualifier()
            break

        # A feature header sits at column 5 with the location at column 21.
        if len(line) > 21 and line[5] != " " and not line.lstrip().startswith("/"):
            flush_qualifier()
            pending_key, pending_value = None, []
            key = line[5:21].strip()
            location = line[21:].strip()
            spans = _LOCATION.findall(location)
            if not spans:
                current = None
                continue
            strand = -1 if location.startswith("complement") else 1
            start = min(int(a) for a, _ in spans) - 1
            end = max(int(b) for _, b in spans)
            current = GenBankFeature(key, start, end, strand)
            features.append(current)
        elif line.lstrip().startswith("/"):
            flush_qualifier()
            body = line.strip()[1:]
            if "=" in body:
                pending_key, value = body.split("=", 1)
                pending_value = [value]
            else:
                pending_key, pending_value = body, [""]
        elif current is not None and pending_key is not None:
            pending_value.append(line.strip())

    return GenBankRecord(
        name=name, sequence=sequence, features=features,
        is_circular=circular, definition=definition,
    )


def read_genbank(path: str | Path) -> GenBankRecord:
    return parse_genbank(Path(path).read_text())




def write_genbank(record: GenBankRecord, molecule: str = "ds-DNA") -> str:
    """Render a :class:`GenBankRecord` back to GenBank text."""
    import textwrap
    from datetime import date

    topology = "circular" if record.is_circular else "linear"
    stamp = date.today().strftime("%d-%b-%Y").upper()
    lines = [
        f"LOCUS       {record.name:<24}{len(record.sequence)} bp {molecule}"
        f"     {topology}     {stamp}",
        f"DEFINITION  {record.definition or '.'}",
        "FEATURES             Location/Qualifiers",
    ]
    for feature in sorted(record.features, key=lambda f: (f.start, -len(f))):
        span = f"{feature.start + 1}..{feature.end}"
        location = f"complement({span})" if feature.strand < 0 else span
        lines.append(f"     {feature.key:<16}{location}")
        for key, value in feature.qualifiers.items():
            body = f'/{key}="{str(value).replace(chr(34), chr(39))}"'
            lines += textwrap.wrap(body, width=58, initial_indent=" " * 21,
                                   subsequent_indent=" " * 21) or [" " * 21 + body]
    lines.append("ORIGIN")
    sequence = record.sequence.lower()
    for i in range(0, len(sequence), 60):
        chunk = sequence[i : i + 60]
        groups = " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{i + 1:>9} {groups}")
    lines.append("//")
    return "\n".join(lines) + "\n"


__all__ = [
    "GenBankFeature", "GenBankRecord",
    "parse_genbank", "read_genbank", "write_genbank",
]
