# Rebuilding the design-rationale document

The figures are authored once, as inline SVG in `docs/design_rationale.html`,
and rendered from there — so the Word document and the web page can never drift
apart.

```bash
cd docs/build
npm install docx playwright        # docx is not preinstalled in every image
node render_figs.js                # SVG -> figures/fig-1..5.png at 3x
node build_docx.js                 # -> P0-P8_design_rationale.docx
python3 inspect_docx.py ../P0-P8_design_rationale.docx   # structural check
```

`render_figs.js` pins the browser to the preinstalled Chromium at
`/opt/pw-browsers/chromium-1194/chrome-linux/chrome`; change that path
elsewhere.

`inspect_docx.py` reads the finished file back and asserts the outline — five
figures each followed by a caption, three tables, and every key term present.
Use it in place of opening the document when no renderer is available.

## Figures are generated, not hand-drawn

`scripts/build_figures.py` writes the `<svg>` blocks into
`docs/design_rationale.html` between the `<!--FIGURES-->` markers. The panel
figure is built from `txv.pvax1_ag`, so changing the architecture changes the
drawing.

That exists because the previous Figure 5 still showed pVax1_AG, BsaI and a
CMV promoter long after the design had moved to pJET1.2 — it was hand-typed,
and nothing tied it to the design.

```bash
python scripts/build_figures.py         # -> SVG into design_rationale.html
cd docs/build
cp ../design_rationale.html panel.html  # render_figs.js reads panel.html
node render_figs.js                     # -> fig-1..6.png at 3x
cp fig-*.png ../figures/
node build_docx.js && cp P0-P8_design_rationale.docx ../
python3 inspect_docx.py ../P0-P8_design_rationale.docx
```

`inspect_docx.py` asserts the figure count and a list of key terms; both are
constants at the top of that file and both describe the *current* design, so
update them when the design changes rather than working around a failure.

**Known gap:** `build_docx.js` carries its own copy of the prose, separate from
`design_rationale.html`. The two must be edited together. Merging them is
worth doing and has not been done.
