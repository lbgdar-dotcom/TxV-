"""Generate the design-rationale figures as inline SVG.

Figure 5 in the previous set drew pVax1_AG, BsaI and a CMV promoter long after
the design had moved to pJET1.2 -- because it was hand-typed and nothing tied
it to the design. The panel figure here is built from txv.pvax1_ag, so a change
to the architecture changes the drawing, and a figure cannot quietly describe a
construct that no longer exists.

Writes the <svg> blocks into docs/design_rationale.html between markers, so
docs/build/render_figs.js keeps working unchanged.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from txv.pvax1_ag import MODULES, PANEL, build_panel_construct  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "design_rationale.html"

#: Fill and stroke token per module kind, so one colour means one thing in
#: every figure.
STYLE = {
    "CTLA4_SP":  ("var(--cyto-soft)", "var(--cyto)",  "SP"),
    "LAMP1_SP":  ("var(--endo-soft)", "var(--endo)",  "SP"),
    "MA":        ("var(--surface-2)", "var(--ink-3)", ""),
    "HA":        ("var(--surface-2)", "var(--ink-2)", "HA"),
    "FLAG":      ("var(--surface-2)", "var(--ink-2)", "FLAG"),
    "A":         ("var(--agA-soft)",  "var(--agA)",   "antigen A"),
    "B":         ("var(--agB-soft)",  "var(--agB)",   "antigen B"),
    "L":         ("var(--surface-2)", "var(--rule)",  ""),
    "E5":        ("var(--warn-soft)", "var(--warn)",  ""),
    "P2A":       ("var(--surface)",   "var(--ink)",   "P2A"),
    "CTLA4_TMT": ("var(--cyto-soft)", "var(--cyto)",  "TM + YVKM"),
    "LAMP1_TMT": ("var(--endo-soft)", "var(--endo)",  "TM + GYQTI"),
}

ROUTE_COLOUR = {
    "CTLA-4 (AP-2, via surface)": "var(--cyto)",
    "LAMP1 (AP-3, direct to lysosome)": "var(--endo)",
    "cytosolic (+E5 degron, free C-terminus)": "var(--ink-2)",
    "dual: A cytosolic + B lysosomal": "var(--warn)",
    "dual: B cytosolic + A lysosomal": "var(--warn)",
}

SHORT_ROUTE = {
    "CTLA-4 (AP-2, via surface)": "CTLA-4",
    "LAMP1 (AP-3, direct to lysosome)": "LAMP1",
    "cytosolic (+E5 degron, free C-terminus)": "cytosolic",
    "dual: A cytosolic + B lysosomal": "dual",
    "dual: B cytosolic + A lysosomal": "dual",
}


def esc(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def panel_figure() -> str:
    """Every construct in the panel, to scale, with what each one is for."""
    px_per_aa = 3.4
    x0, row_h, top = 116, 42, 92
    rows = []

    # -- legend ------------------------------------------------------------
    # Colour encodes the routing decision, not the module type: the signal
    # peptide and TM tail of a CTLA-4 construct are both orange, of a LAMP1
    # construct both teal. Labelling the swatch "signal peptide" would have
    # been wrong, since P6-P8's signal peptide is teal.
    legend = [
        ("CTLA4_SP", "CTLA-4 route"),
        ("LAMP1_SP", "LAMP1 route"),
        ("A", "antigen A"),
        ("B", "antigen B"),
        ("L", "linker / tag"),
        ("E5", "E5 degron"),
        ("P2A", "P2A skip"),
    ]
    lx = 8
    for key, label in legend:
        fill, stroke, _ = STYLE[key]
        rows.append(
            f'<rect x="{lx}" y="12" width="13" height="13" rx="3" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>'
            f'<text x="{lx + 18}" y="23" font-size="11.5" '
            f'fill="var(--ink-2)">{esc(label)}</text>')
        lx += 24 + len(label) * 6.3
    rows.append('<text x="8" y="52" font-size="11.5" fill="var(--ink-3)">'
                'colour marks the routing module a construct uses; box width is '
                'proportional to length in amino acids</text>')
    rows.append('<line x1="8" y1="66" x2="852" y2="66" '
                'stroke="var(--rule)" stroke-width="1"/>')

    for i, panel in enumerate(PANEL):
        protein = build_panel_construct(panel.name).orf_protein
        y = top + i * row_h
        colour = ROUTE_COLOUR[panel.route]

        rows.append(
            f'<text x="8" y="{y + 15}" font-size="14" font-weight="600" '
            f'fill="var(--ink)" font-family="var(--mono)">{panel.name}</text>'
            f'<text x="40" y="{y + 15}" font-size="11" fill="{colour}">'
            f'{esc(SHORT_ROUTE[panel.route])}</text>')

        x = x0
        for module in panel.modules:
            width = max(len(MODULES[module]) * px_per_aa, 5.0)
            fill, stroke, label = STYLE[module]
            rows.append(
                f'<rect x="{x:.1f}" y="{y}" width="{width:.1f}" height="21" '
                f'rx="3" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>')
            if label and width > len(label) * 6.2:
                rows.append(
                    f'<text x="{x + width / 2:.1f}" y="{y + 14.5}" '
                    f'font-size="10.5" text-anchor="middle" fill="{stroke}">'
                    f'{esc(label)}</text>')
            x += width + 1.5

        rows.append(
            f'<text x="{x + 8:.1f}" y="{y + 15}" font-size="11" '
            f'fill="var(--ink-3)" font-family="var(--mono)">'
            f'{len(protein)} aa</text>')

    height = top + len(PANEL) * row_h + 8
    body = "\n        ".join(rows)
    return (
        f'<svg viewBox="0 0 860 {height}" role="img" aria-label="All nine '
        f'constructs drawn to scale. P0 to P4 carry the CTLA-4 tail, P5 is '
        f'cytosolic, P6 carries the LAMP1 tail, and P7 and P8 use P2A to make '
        f'two separately routed proteins from one transcript.">\n'
        f'        {body}\n      </svg>')


# --- the figures that are drawn rather than derived -------------------------

def compartment_figure() -> str:
    return '''<svg viewBox="0 0 860 300" role="img" aria-label="A bare mRNA is translated on cytosolic ribosomes, reaches the proteasome and MHC class I, but has no route into the MIIC where class II is loaded.">
        <rect x="8" y="30" width="520" height="250" rx="10" fill="var(--surface-2)" stroke="var(--rule)"/>
        <text x="24" y="54" font-size="12" font-family="var(--mono)" fill="var(--ink-3)">CYTOSOL</text>
        <path d="M30 120 q14 -10 28 0 t28 0 t28 0" fill="none" stroke="var(--ink-2)" stroke-width="2"/>
        <text x="30" y="150" font-size="11.5" fill="var(--ink-2)">IVT mRNA</text>
        <circle cx="132" cy="120" r="14" fill="var(--surface)" stroke="var(--ink)" stroke-width="2"/>
        <text x="132" y="150" font-size="11.5" text-anchor="middle" fill="var(--ink-2)">ribosome</text>
        <line x1="150" y1="120" x2="196" y2="120" stroke="var(--ink)" stroke-width="2" marker-end="url(#ar)"/>
        <rect x="200" y="103" width="132" height="34" rx="5" fill="var(--surface)" stroke="var(--ink)" stroke-width="2"/>
        <text x="266" y="125" font-size="12.5" text-anchor="middle" fill="var(--ink)">antigen protein</text>
        <line x1="266" y1="141" x2="266" y2="186" stroke="var(--cyto)" stroke-width="2" marker-end="url(#arC)"/>
        <text x="276" y="168" font-size="11.5" fill="var(--cyto)">degraded</text>
        <rect x="200" y="190" width="132" height="32" rx="5" fill="var(--cyto-soft)" stroke="var(--cyto)" stroke-width="2"/>
        <text x="266" y="211" font-size="12.5" text-anchor="middle" fill="var(--cyto)">proteasome</text>
        <line x1="334" y1="206" x2="380" y2="206" stroke="var(--cyto)" stroke-width="2" marker-end="url(#arC)"/>
        <text x="357" y="198" font-size="10.5" text-anchor="middle" fill="var(--cyto)">TAP</text>
        <rect x="384" y="190" width="126" height="32" rx="5" fill="var(--surface)" stroke="var(--cyto)" stroke-width="2"/>
        <text x="447" y="211" font-size="12.5" text-anchor="middle" fill="var(--cyto)">ER &#8594; MHC I</text>
        <text x="447" y="246" font-size="12.5" text-anchor="middle" fill="var(--cyto)">CD8 response</text>
        <line x1="447" y1="224" x2="447" y2="234" stroke="var(--cyto)" stroke-width="2"/>
        <line x1="340" y1="120" x2="660" y2="120" stroke="var(--ink-3)" stroke-width="2" stroke-dasharray="7 5" marker-end="url(#arG)"/>
        <text x="500" y="110" font-size="12" text-anchor="middle" fill="var(--warn)">no route from the cytosol</text>
        <line x1="492" y1="112" x2="508" y2="128" stroke="var(--warn)" stroke-width="2.5"/>
        <line x1="508" y1="112" x2="492" y2="128" stroke="var(--warn)" stroke-width="2.5"/>
        <rect x="664" y="60" width="188" height="120" rx="10" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="2"/>
        <text x="758" y="88" font-size="13" font-weight="600" text-anchor="middle" fill="var(--endo)">MIIC</text>
        <text x="758" y="110" font-size="11.5" text-anchor="middle" fill="var(--endo)">late endosome lumen</text>
        <text x="758" y="130" font-size="11.5" text-anchor="middle" fill="var(--endo)">cathepsins cleave here</text>
        <line x1="758" y1="142" x2="758" y2="172" stroke="var(--endo)" stroke-width="2" marker-end="url(#arE)"/>
        <rect x="684" y="190" width="148" height="32" rx="5" fill="var(--surface)" stroke="var(--endo)" stroke-width="2"/>
        <text x="758" y="211" font-size="12.5" text-anchor="middle" fill="var(--endo)">MHC II</text>
        <text x="758" y="246" font-size="12.5" text-anchor="middle" fill="var(--endo)">CD4 help</text>
        <line x1="758" y1="224" x2="758" y2="234" stroke="var(--endo)" stroke-width="2"/>
      </svg>'''


def routes_figure() -> str:
    return '''<svg viewBox="0 0 860 300" role="img" aria-label="Both sorting tails deliver antigen to the MIIC, but the CTLA-4 YVKM tail travels via the plasma membrane and is retrieved by AP-2, while the LAMP1 GYQTI tail goes directly from the trans-Golgi to the lysosome via AP-3.">
        <rect x="8" y="8" width="844" height="128" rx="9" fill="var(--cyto-soft)" stroke="var(--cyto)" stroke-width="1.5"/>
        <text x="24" y="32" font-size="13" font-weight="600" fill="var(--cyto)">CTLA-4 tail &#8212; YVKM, AP-2</text>
        <text x="228" y="32" font-size="11.5" font-family="var(--mono)" fill="var(--cyto)">P0 P1 P2 P3 P4</text>
        <rect x="24" y="52" width="120" height="34" rx="5" fill="var(--surface)" stroke="var(--cyto)" stroke-width="1.5"/>
        <text x="84" y="73" font-size="12" text-anchor="middle" fill="var(--cyto)">ER &#8594; Golgi</text>
        <line x1="146" y1="69" x2="196" y2="69" stroke="var(--cyto)" stroke-width="2" marker-end="url(#arC)"/>
        <rect x="200" y="52" width="150" height="34" rx="5" fill="var(--surface)" stroke="var(--cyto)" stroke-width="1.5"/>
        <text x="275" y="73" font-size="12" text-anchor="middle" fill="var(--cyto)">plasma membrane</text>
        <line x1="352" y1="69" x2="402" y2="69" stroke="var(--cyto)" stroke-width="2" marker-end="url(#arC)"/>
        <text x="377" y="61" font-size="10.5" text-anchor="middle" fill="var(--cyto)">AP-2</text>
        <rect x="406" y="52" width="150" height="34" rx="5" fill="var(--surface)" stroke="var(--cyto)" stroke-width="1.5"/>
        <text x="481" y="73" font-size="12" text-anchor="middle" fill="var(--cyto)">early endosome</text>
        <line x1="558" y1="69" x2="608" y2="69" stroke="var(--cyto)" stroke-width="2" marker-end="url(#arC)"/>
        <rect x="612" y="52" width="130" height="34" rx="5" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="2"/>
        <text x="677" y="73" font-size="12" text-anchor="middle" fill="var(--endo)">MIIC</text>
        <text x="24" y="112" font-size="11.5" fill="var(--cyto)">exposed at the surface on the way &#8212; detectable by flow without permeabilising</text>
        <rect x="8" y="150" width="844" height="128" rx="9" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="1.5"/>
        <text x="24" y="174" font-size="13" font-weight="600" fill="var(--endo)">LAMP1 tail &#8212; GYQTI, AP-3</text>
        <text x="228" y="174" font-size="11.5" font-family="var(--mono)" fill="var(--endo)">P6 P7 P8</text>
        <rect x="24" y="194" width="120" height="34" rx="5" fill="var(--surface)" stroke="var(--endo)" stroke-width="1.5"/>
        <text x="84" y="215" font-size="12" text-anchor="middle" fill="var(--endo)">ER &#8594; Golgi</text>
        <line x1="146" y1="211" x2="196" y2="211" stroke="var(--endo)" stroke-width="2" marker-end="url(#arE)"/>
        <rect x="200" y="194" width="150" height="34" rx="5" fill="var(--surface)" stroke="var(--endo)" stroke-width="1.5"/>
        <text x="275" y="215" font-size="12" text-anchor="middle" fill="var(--endo)">trans-Golgi network</text>
        <line x1="352" y1="211" x2="608" y2="211" stroke="var(--endo)" stroke-width="2" marker-end="url(#arE)"/>
        <text x="480" y="203" font-size="10.5" text-anchor="middle" fill="var(--endo)">AP-3 &#8212; direct, no surface step</text>
        <rect x="612" y="194" width="130" height="34" rx="5" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="2"/>
        <text x="677" y="215" font-size="12" text-anchor="middle" fill="var(--endo)">MIIC</text>
        <text x="24" y="254" font-size="11.5" fill="var(--endo)">never at the surface &#8212; so surface staining separates the two routes cleanly</text>
      </svg>'''


def position_figure() -> str:
    return '''<svg viewBox="0 0 860 250" role="img" aria-label="Three motifs whose function depends on position: the signal peptide must be N-terminal, GYQTI must be the last five residues, and the E5 degron needs a free C-terminus, which it lacks in P7 and P8.">
        <text x="8" y="18" font-size="12.5" font-weight="600" fill="var(--ink)">Motifs that only work in one place</text>
        <rect x="8" y="34" width="180" height="26" rx="4" fill="var(--cyto-soft)" stroke="var(--cyto)" stroke-width="1.5"/>
        <text x="98" y="52" font-size="11.5" text-anchor="middle" fill="var(--cyto)">signal peptide</text>
        <text x="200" y="52" font-size="11.5" fill="var(--ink-2)">must be N-terminal &#8212; SRP reads the chain as it emerges</text>
        <rect x="8" y="72" width="180" height="26" rx="4" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="1.5"/>
        <text x="98" y="90" font-size="11.5" text-anchor="middle" fill="var(--endo)">GYQTI</text>
        <text x="200" y="90" font-size="11.5" fill="var(--ink-2)">must be the last five residues &#8212; otherwise AP-3 fails, silently</text>
        <rect x="8" y="110" width="180" height="26" rx="4" fill="var(--warn-soft)" stroke="var(--warn)" stroke-width="1.5"/>
        <text x="98" y="128" font-size="11.5" text-anchor="middle" fill="var(--warn)">E5 degron</text>
        <text x="200" y="128" font-size="11.5" fill="var(--ink-2)">needs a free C-terminus &#8212; an acidic degron is read from the end</text>
        <line x1="8" y1="152" x2="852" y2="152" stroke="var(--rule)"/>
        <text x="8" y="176" font-size="12" font-weight="600" fill="var(--ink)">The consequence in P7 and P8</text>
        <rect x="8" y="190" width="150" height="24" rx="4" fill="var(--agA-soft)" stroke="var(--agA)" stroke-width="1.5"/>
        <text x="83" y="207" font-size="11" text-anchor="middle" fill="var(--agA)">cytosolic arm</text>
        <rect x="160" y="190" width="54" height="24" rx="4" fill="var(--warn-soft)" stroke="var(--warn)" stroke-width="1.5"/>
        <text x="187" y="207" font-size="11" text-anchor="middle" fill="var(--warn)">E5</text>
        <rect x="216" y="190" width="70" height="24" rx="4" fill="var(--surface)" stroke="var(--ink)" stroke-width="1.5"/>
        <text x="251" y="207" font-size="11" text-anchor="middle" fill="var(--ink)">P2A</text>
        <rect x="288" y="190" width="200" height="24" rx="4" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="1.5"/>
        <text x="388" y="207" font-size="11" text-anchor="middle" fill="var(--endo)">LAMP1 arm</text>
        <rect x="490" y="190" width="76" height="24" rx="4" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="2"/>
        <text x="528" y="207" font-size="11" text-anchor="middle" fill="var(--endo)">GYQTI</text>
        <text x="576" y="207" font-size="11.5" fill="var(--endo)">&#9664; last: works</text>
        <line x1="187" y1="186" x2="187" y2="172" stroke="var(--warn)" stroke-width="1.5"/>
        <text x="187" y="168" font-size="11" text-anchor="middle" fill="var(--warn)">internal: expected weak</text>
        <text x="8" y="238" font-size="11.5" fill="var(--ink-3)">P5 is the degron-positive comparator &#8212; the only construct where E5 ends the protein and the route is cytosolic.</text>
      </svg>'''


def cloning_figure() -> str:
    return '''<svg viewBox="0 0 860 290" role="img" aria-label="The gene fragment is inserted at the single blunt Eco32I site of pJET1.2, which disrupts the lethal eco47IR gene so only recombinants grow, and the two flanking BglII sites both release the insert for checking and produce the IVT template.">
        <text x="8" y="18" font-size="12.5" font-family="var(--mono)" fill="var(--ink-3)">pJET1.2 &#183; 2974 bp &#183; circular &#183; AmpR</text>
        <text x="352" y="70" font-size="11" text-anchor="middle" fill="var(--warn)">Eco32I blunt site 371/372</text>
        <line x1="352" y1="76" x2="352" y2="152" stroke="var(--warn)" stroke-width="2"/>
        <rect x="8" y="96" width="120" height="36" rx="5" fill="var(--surface-2)" stroke="var(--rule)" stroke-width="1.5"/>
        <text x="68" y="114" font-size="10.5" text-anchor="middle" fill="var(--ink-2)">rest of backbone</text>\n        <text x="68" y="127" font-size="9.5" text-anchor="middle" fill="var(--ink-3)">ori, AmpR &#8212; elsewhere on the circle</text>
        <rect x="132" y="96" width="216" height="36" rx="5" fill="var(--warn-soft)" stroke="var(--warn)" stroke-width="1.5"/>
        <text x="240" y="113" font-size="11.5" text-anchor="middle" fill="var(--warn)">eco47IR &#8212; lethal</text>
        <text x="240" y="127" font-size="10" text-anchor="middle" fill="var(--warn)">broken by the insert; that is the selection</text>
        <rect x="356" y="96" width="330" height="36" rx="5" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="2"/>
        <text x="521" y="119" font-size="12.5" font-weight="600" text-anchor="middle" fill="var(--endo)">your gene fragment &#8212; 670 to 955 bp</text>
        <rect x="690" y="96" width="162" height="36" rx="5" fill="var(--surface-2)" stroke="var(--rule)" stroke-width="1.5"/>
        <text x="771" y="119" font-size="11.5" text-anchor="middle" fill="var(--ink-2)">eco47IR continues</text>
        <text x="132" y="88" font-size="10.5" fill="var(--ink-3)">BglII 337</text>
        <line x1="132" y1="92" x2="132" y2="152" stroke="var(--ink-3)" stroke-width="1.5" stroke-dasharray="4 3"/>
        <text x="688" y="88" font-size="10.5" fill="var(--ink-3)">BglII</text>
        <line x1="688" y1="92" x2="688" y2="152" stroke="var(--ink-3)" stroke-width="1.5" stroke-dasharray="4 3"/>
        <line x1="132" y1="162" x2="688" y2="162" stroke="var(--ink-2)" stroke-width="2" marker-start="url(#arS)" marker-end="url(#arG)"/>
        <text x="410" y="182" font-size="11.5" text-anchor="middle" fill="var(--ink-2)">one BglII digest &#8212; the diagnostic band and the IVT template are the same cut</text>
        <line x1="8" y1="200" x2="852" y2="200" stroke="var(--rule)"/>
        <text x="8" y="224" font-size="12" font-weight="600" fill="var(--ink)">Inside the fragment</text>
        <rect x="8" y="236" width="70" height="24" rx="4" fill="var(--surface-2)" stroke="var(--rule)" stroke-width="1.5"/>
        <text x="43" y="253" font-size="10" text-anchor="middle" fill="var(--ink-2)">IVT_F</text>
        <rect x="80" y="236" width="86" height="24" rx="4" fill="var(--surface)" stroke="var(--ink)" stroke-width="1.5"/>
        <text x="123" y="253" font-size="10" text-anchor="middle" fill="var(--ink)">T7 &#8594; AGG</text>
        <rect x="168" y="236" width="86" height="24" rx="4" fill="var(--surface)" stroke="var(--ink-2)" stroke-width="1.5"/>
        <text x="211" y="253" font-size="10" text-anchor="middle" fill="var(--ink-2)">5&#8242;UTR</text>
        <rect x="256" y="236" width="56" height="24" rx="4" fill="var(--surface-2)" stroke="var(--ink-2)" stroke-width="1.5"/>
        <text x="284" y="253" font-size="10" text-anchor="middle" fill="var(--ink-2)">Kozak</text>
        <rect x="314" y="236" width="230" height="24" rx="4" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="2"/>
        <text x="429" y="253" font-size="10.5" text-anchor="middle" fill="var(--endo)">ORF &#8212; the construct</text>
        <rect x="546" y="236" width="150" height="24" rx="4" fill="var(--surface)" stroke="var(--ink-2)" stroke-width="1.5"/>
        <text x="621" y="253" font-size="10" text-anchor="middle" fill="var(--ink-2)">3&#8242;UTR &#8212; AES + mtRNR1</text>
        <rect x="698" y="236" width="70" height="24" rx="4" fill="var(--surface-2)" stroke="var(--rule)" stroke-width="1.5"/>
        <text x="733" y="253" font-size="10" text-anchor="middle" fill="var(--ink-2)">IVT_R</text>
        <text x="776" y="253" font-size="10.5" fill="var(--cyto)">+ A120</text>
        <text x="8" y="280" font-size="11" fill="var(--ink-3)">The poly(A) is added by PCR before cloning &#8212; a 120-nt A-tract cannot be synthesised into a gene fragment.</text>
      </svg>'''


def pair_figure() -> str:
    """Kept from the previous set: why both modules, always."""
    return '''<svg viewBox="0 0 860 240" role="img" aria-label="Three outcomes. With no signal peptide the protein stays cytosolic. With a signal peptide alone it is secreted out of the cell. With both a signal peptide and a transmembrane domain the cassette sits in the lumen and the sorting tail in the cytosol.">
        <text x="8" y="18" font-size="12.5" font-weight="600" fill="var(--ink)">The signal peptide and the TM domain only work as a pair</text>
        <rect x="8" y="34" width="270" height="180" rx="9" fill="var(--surface)" stroke="var(--rule)" stroke-width="1.5"/>
        <text x="143" y="58" font-size="12" font-weight="600" text-anchor="middle" fill="var(--ink-2)">neither</text>
        <rect x="52" y="76" width="182" height="26" rx="4" fill="var(--surface-2)" stroke="var(--ink-3)" stroke-width="1.5"/>
        <text x="143" y="94" font-size="11" text-anchor="middle" fill="var(--ink-2)">antigen cassette</text>
        <line x1="143" y1="110" x2="143" y2="140" stroke="var(--ink-3)" stroke-width="2" marker-end="url(#arG)"/>
        <text x="143" y="162" font-size="11.5" text-anchor="middle" fill="var(--ink-2)">stays cytosolic</text>
        <text x="143" y="182" font-size="11" text-anchor="middle" fill="var(--cyto)">class I only</text>
        <text x="143" y="202" font-size="10.5" text-anchor="middle" fill="var(--ink-3)">this is P5</text>
        <rect x="294" y="34" width="270" height="180" rx="9" fill="var(--surface)" stroke="var(--rule)" stroke-width="1.5"/>
        <text x="429" y="58" font-size="12" font-weight="600" text-anchor="middle" fill="var(--warn)">signal peptide alone</text>
        <rect x="316" y="76" width="50" height="26" rx="4" fill="var(--cyto-soft)" stroke="var(--cyto)" stroke-width="1.5"/>
        <text x="341" y="94" font-size="10" text-anchor="middle" fill="var(--cyto)">SP</text>
        <rect x="370" y="76" width="172" height="26" rx="4" fill="var(--surface-2)" stroke="var(--ink-3)" stroke-width="1.5"/>
        <text x="456" y="94" font-size="11" text-anchor="middle" fill="var(--ink-2)">antigen cassette</text>
        <line x1="429" y1="110" x2="429" y2="140" stroke="var(--warn)" stroke-width="2" marker-end="url(#arW)"/>
        <text x="429" y="162" font-size="11.5" text-anchor="middle" fill="var(--warn)">secreted &#8212; leaves the cell</text>
        <text x="429" y="182" font-size="11" text-anchor="middle" fill="var(--warn)">presented by neither</text>
        <text x="429" y="202" font-size="10.5" text-anchor="middle" fill="var(--ink-3)">not in the panel, by design</text>
        <rect x="580" y="34" width="272" height="180" rx="9" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="2"/>
        <text x="716" y="58" font-size="12" font-weight="600" text-anchor="middle" fill="var(--endo)">both</text>
        <rect x="598" y="76" width="42" height="26" rx="4" fill="var(--cyto-soft)" stroke="var(--cyto)" stroke-width="1.5"/>
        <text x="619" y="94" font-size="10" text-anchor="middle" fill="var(--cyto)">SP</text>
        <rect x="644" y="76" width="124" height="26" rx="4" fill="var(--surface-2)" stroke="var(--ink-3)" stroke-width="1.5"/>
        <text x="706" y="94" font-size="10.5" text-anchor="middle" fill="var(--ink-2)">cassette</text>
        <rect x="772" y="76" width="62" height="26" rx="4" fill="var(--endo-soft)" stroke="var(--endo)" stroke-width="2"/>
        <text x="803" y="94" font-size="10" text-anchor="middle" fill="var(--endo)">TM+tail</text>
        <line x1="716" y1="110" x2="716" y2="140" stroke="var(--endo)" stroke-width="2" marker-end="url(#arE)"/>
        <text x="716" y="162" font-size="11.5" text-anchor="middle" fill="var(--endo)">anchored: cassette luminal,</text>
        <text x="716" y="178" font-size="11.5" text-anchor="middle" fill="var(--endo)">sorting tail cytosolic</text>
        <text x="716" y="202" font-size="10.5" text-anchor="middle" fill="var(--ink-3)">P0&#8211;P4 and P6&#8211;P8</text>
      </svg>'''


#: Arrowhead markers. Injected into every figure rather than shared once, so
#: each <svg> is a self-contained fragment: the renderer pulls figures out of
#: this page one at a time, and a marker defined elsewhere would come out as
#: an arrowless line.
DEFS = '''<defs>
    <marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="var(--ink)"/></marker>
    <marker id="arC" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="var(--cyto)"/></marker>
    <marker id="arE" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="var(--endo)"/></marker>
    <marker id="arG" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="var(--ink-3)"/></marker>
    <marker id="arW" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="var(--warn)"/></marker>
    <marker id="arS" viewBox="0 0 10 10" refX="1" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M10 0 L0 5 L10 10 z" fill="var(--ink-2)"/></marker>
  </defs>'''

FIGURES = [
    (compartment_figure, "Figure 1",
     "The compartments are the whole problem. An mRNA is translated on free "
     "cytosolic ribosomes, so a bare antigen ORF reaches the proteasome and "
     "MHC class I by default and has no route into the MIIC where class II is "
     "loaded. Every construct in this panel is an attempt to build that route."),
    (pair_figure, "Figure 2",
     "Why both modules, always. A signal peptide alone sends the protein out "
     "of the cell, where it is presented by nothing; a transmembrane domain "
     "alone has nothing to insert into. Only the pair leaves the antigen "
     "cassette in the lumen with the sorting tail facing the cytosol, which is "
     "the arrangement the sorting adaptors need. P5 is the deliberate "
     "exception: no signal peptide, cytosolic by design."),
    (panel_figure, "Figure 3",
     "All nine constructs, drawn to scale from the design library rather than "
     "by hand. Reading down: P0 is scaffold only, P1 and P2 add one antigen "
     "each, P3 and P4 carry both in either order, P5 drops the routing modules "
     "entirely, P6 swaps the CTLA-4 tail for LAMP1, and P7 and P8 use P2A to "
     "make two separately routed proteins from one transcript."),
    (routes_figure, "Figure 4",
     "Same destination, different journey. Both tails deliver antigen to the "
     "MIIC, but the CTLA-4 YVKM tail gets there via the plasma membrane and "
     "AP-2 retrieval, while the LAMP1 GYQTI tail goes straight from the "
     "trans-Golgi network by AP-3. The surface step is what makes them "
     "distinguishable: CTLA-4-routed protein is detectable by flow without "
     "permeabilising, LAMP1-routed protein is not."),
    (position_figure, "Figure 5",
     "Three motifs whose function is positional, and what that costs. The "
     "signal peptide must be N-terminal, GYQTI must be the final five "
     "residues, and an acidic degron is read from a free C-terminus. In P7 "
     "and P8 the E5 degron sits upstream of P2A and therefore has no free end, "
     "so those two are expected to show weak degron activity. That is a stated "
     "limitation of the dual-route design, not an oversight."),
    (cloning_figure, "Figure 6",
     "The fragment goes into the single blunt Eco32I site of pJET1.2, which "
     "sits inside the lethal eco47IR gene: breaking it is the selection, so "
     "colonies that grow carry an insert. The two flanking BglII sites do "
     "double duty, releasing the insert for the diagnostic gel and producing "
     "the IVT template in the same digest. No fragment contains a BglII site."),
]


def main() -> int:
    page = PAGE.read_text()
    blocks = []
    for build, label, caption in FIGURES:
        svg = build()
        opening = svg.index(">") + 1        # inject defs after the <svg ...> tag
        svg = svg[:opening] + "\n        " + DEFS + svg[opening:]
        blocks.append(
            f'  <figure>\n      {svg}\n'
            f'    <figcaption><b>{label}.</b> {caption}</figcaption>\n'
            f'  </figure>')

    begin, end = "<!--FIGURES-->", "<!--/FIGURES-->"
    if begin not in page:
        print(f"markers not found in {PAGE}; add {begin} / {end}")
        return 1
    head = page[:page.index(begin) + len(begin)]
    tail = page[page.index(end):]
    PAGE.write_text(head + "\n" + "\n\n".join(blocks) + "\n  " + tail)

    print(f"wrote {len(FIGURES)} figures into {PAGE.relative_to(ROOT)}")
    for _, label, _ in FIGURES:
        print(f"  {label}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
