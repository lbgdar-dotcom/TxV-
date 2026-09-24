"""Read back the generated .docx and print its outline, so the structure can be
checked without a renderer."""

import re
import sys
import zipfile
import xml.etree.ElementTree as ET

#: What the document is expected to contain. These are assertions about the
#: current design, so they change when the design does -- which is the point:
#: the previous set still demanded "BsaI" and "699" long after the project
#: moved to pJET1.2, and would have passed a document describing the wrong
#: vector while failing the right one.
FIGURES = 6
KEY_TERMS = ("GYQTI", "YVKM", "Eco32I", "BglII", "eco47IR", "371/372",
             "P2A", "mtRNR1", "TAATACGACTCACTATAAGG")

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"


def text_of(el) -> str:
    return "".join(t.text or "" for t in el.iter(f"{W}t"))


def main(path: str) -> None:
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml")
        media = [n for n in z.namelist()
                 if n.startswith("word/media/") and not n.endswith("/")]

    body = ET.fromstring(xml).find(f"{W}body")
    headings = tables = images = paras = bullets = page_breaks = 0
    outline = []

    for child in body:
        if child.tag == f"{W}p":
            paras += 1
            style = child.find(f"{W}pPr/{W}pStyle")
            name = style.get(f"{W}val") if style is not None else ""
            has_img = child.find(f".//{A}blip") is not None
            if has_img:
                images += 1
                outline.append("      [FIGURE IMAGE]")
            if child.find(f".//{W}br[@{W}type='page']") is not None:
                page_breaks += 1
                outline.append("      --- page break ---")
            if child.find(f"{W}pPr/{W}numPr") is not None:
                bullets += 1
            body_text = text_of(child).strip()
            if name and name.startswith("Heading"):
                headings += 1
                outline.append(f"{'  ' * int(name[-1])}H{name[-1]}  {body_text}")
            elif body_text.startswith("Figure "):
                outline.append(f"      caption: {body_text[:72]}…")
            elif body_text.isupper() and len(body_text) < 60 and body_text:
                outline.append(f"  eyebrow: {body_text}")
        elif child.tag == f"{W}tbl":
            tables += 1
            rows = child.findall(f"{W}tr")
            cols = len(rows[0].findall(f"{W}tc")) if rows else 0
            first = text_of(rows[0]).strip()[:60] if rows else ""
            outline.append(f"      [TABLE {len(rows)}x{cols}] {first}")

    print("\n".join(outline))
    print()
    print(f"paragraphs   {paras}")
    print(f"headings     {headings}")
    print(f"tables       {tables}")
    print(f"images       {images} (media files in package: {len(media)})")
    print(f"bullets      {bullets}")
    print(f"page breaks  {page_breaks}")

    problems = []
    if images != FIGURES:
        problems.append(f"expected {FIGURES} figures, found {images}")
    if len(media) != FIGURES:
        problems.append(f"expected {FIGURES} media files, found {len(media)}")
    if tables != 3:
        problems.append(f"expected 3 tables (panel, elements, callout), found {tables}")
    full = " ".join(text_of(p) for p in body.iter(f"{W}p"))
    for token in KEY_TERMS:
        if token not in full:
            problems.append(f"missing key term: {token}")
    if re.search(r"\*\*|`[a-z]", full):
        problems.append("unconverted inline markup left in the text")
    print()
    print("PROBLEMS:", problems or "none")


if __name__ == "__main__":
    main(sys.argv[1])
