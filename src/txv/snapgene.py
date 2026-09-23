"""Reader for SnapGene ``.dna`` files.

SnapGene is what most people actually have a vector map in, and asking someone
to export to GenBank before the tooling will look at their plasmid is a good
way to have the tooling skipped. This reads the format directly.

The file is a flat sequence of segments, each a one-byte type followed by a
four-byte big-endian length. Only three matter here:

* **0** -- the DNA. One flags byte (bit 0 set means circular) then the bases.
* **5** -- the primers, as XML. Worth reading: a vector's stock sequencing
  primers decide whether an insert can actually be verified.
* **6** -- the notes, as XML, which carry the name and description.
* **10** -- the features, as XML.

Everything else (enzyme sets, history, drawing state) is skipped. This is a
reader, not a writer: the point is to get a vector into
:mod:`txv.genbank_io`'s record type and then work in one representation.
"""

from __future__ import annotations

import re
import struct
from pathlib import Path
from xml.etree import ElementTree

from .genbank_io import GenBankFeature, GenBankRecord

SEG_DNA, SEG_PRIMERS, SEG_NOTES, SEG_FEATURES = 0, 5, 6, 10

#: SnapGene feature types that map onto GenBank keys of a different name.
_KEY_ALIASES = {"misc_feature": "misc_feature", "rep_origin": "rep_origin"}


def _segments(data: bytes):
    offset = 0
    while offset + 5 <= len(data):
        kind = data[offset]
        (length,) = struct.unpack(">I", data[offset + 1 : offset + 5])
        yield kind, data[offset + 5 : offset + 5 + length]
        offset += 5 + length


def _strip_html(text: str) -> str:
    """SnapGene wraps qualifier text in HTML; GenBank wants it plain."""
    text = re.sub(r"<br\s*/?>", " ", text)
    text = re.sub(r"<[^>]+>", "", text)
    return " ".join(text.split())


def read_primers(path: str | Path) -> dict[str, str]:
    """The vector's stored primers, as ``{name: sequence}``."""
    for kind, payload in _segments(Path(path).read_bytes()):
        if kind != SEG_PRIMERS:
            continue
        root = ElementTree.fromstring(payload.rstrip(b"\x00"))
        return {p.get("name", f"primer_{i}"): p.get("sequence", "").upper()
                for i, p in enumerate(root.findall("Primer"))}
    return {}


def read_snapgene(path: str | Path) -> GenBankRecord:
    """Parse a ``.dna`` file into a :class:`~txv.genbank_io.GenBankRecord`."""
    data = Path(path).read_bytes()
    sequence, circular, features, name, definition = "", False, [], "", ""

    for kind, payload in _segments(data):
        if kind == SEG_DNA:
            circular = bool(payload[0] & 1)
            sequence = payload[1:].decode("ascii").upper()
        elif kind == SEG_FEATURES:
            features = _parse_features(payload)
        elif kind == SEG_NOTES:
            root = ElementTree.fromstring(payload.rstrip(b"\x00"))
            label = root.findtext("CustomMapLabel") or ""
            name = label or Path(path).stem
            definition = _strip_html(root.findtext("Description") or "")

    if not sequence:
        raise ValueError(f"{path}: no DNA segment found")

    return GenBankRecord(
        name=name or Path(path).stem,
        sequence=sequence,
        features=features,
        is_circular=circular,
        definition=definition,
    )


def _parse_features(payload: bytes) -> list[GenBankFeature]:
    root = ElementTree.fromstring(payload.rstrip(b"\x00"))
    out: list[GenBankFeature] = []

    for element in root.findall("Feature"):
        segments = element.findall("Segment")
        if not segments:
            continue
        # A feature may be split into segments (a CDS with a signal sequence,
        # a promoter with -10 and -35 boxes). Span them: this reader is for
        # locating things on a vector, not for reproducing SnapGene's drawing.
        starts, ends = [], []
        for segment in segments:
            lo, hi = segment.get("range", "0-0").split("-")
            starts.append(int(lo) - 1)      # 1-based inclusive -> 0-based
            ends.append(int(hi))
        strand = -1 if element.get("directionality") == "2" else 1

        qualifiers = {"label": element.get("name", "")}
        for q in element.findall("Q"):
            value = q.find("V")
            if value is None:
                continue
            text = value.get("text")
            if text is None:
                text = value.get("int", "")
            cleaned = _strip_html(text)
            if cleaned:
                qualifiers[q.get("name", "note")] = cleaned

        key = element.get("type", "misc_feature")
        out.append(GenBankFeature(
            key=_KEY_ALIASES.get(key, key),
            start=min(starts), end=max(ends), strand=strand,
            qualifiers=qualifiers,
        ))
    out.sort(key=lambda f: (f.start, f.end))
    return out
