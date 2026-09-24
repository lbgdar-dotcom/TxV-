const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, TextRun, ImageRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, BorderStyle, ShadingType, ShadingType: _S,
  LevelFormat, PageBreak, PositionalTab,
} = require('docx');

const dir = __dirname;
const CONTENT = 9360;          // DXA, Letter with 1" margins
const IMG_W = 624;             // px at 96dpi = 6.5"

const INK = '15222C', MUTED = '4A6373', TEAL = '0C6F7A', AMBER = 'B0590F',
      WARN = '9C2C37', RULE = 'D2DBE1', SOFT = 'EFF3F5', WARNSOFT = 'F9ECEE';

const BODY = 'Calibri', HEAD = 'Cambria', MONO = 'Consolas';

function pngSize(file) {
  const b = fs.readFileSync(file);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) };
}

// --- inline text helpers ---------------------------------------------------
// Mini markup: **bold**, _italic_, `mono`
function runs(text, opts = {}) {
  const base = { font: BODY, size: 21, color: INK, ...opts };
  const out = [];
  const re = /(\*\*[^*]+\*\*|_[^_]+_|`[^`]+`)/g;
  let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last) out.push(new TextRun({ ...base, text: text.slice(last, m.index) }));
    const tok = m[0];
    if (tok.startsWith('**')) out.push(new TextRun({ ...base, text: tok.slice(2, -2), bold: true }));
    else if (tok.startsWith('_')) out.push(new TextRun({ ...base, text: tok.slice(1, -1), italics: true }));
    else out.push(new TextRun({ ...base, text: tok.slice(1, -1), font: MONO, size: 19 }));
    last = m.index + tok.length;
  }
  if (last < text.length) out.push(new TextRun({ ...base, text: text.slice(last) }));
  return out;
}

const p = (text, opts = {}) => new Paragraph({
  children: runs(text, opts.run || {}),
  spacing: { after: opts.after ?? 160, line: 276 },
  ...(opts.para || {}),
});

const h1 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_1,
  spacing: { before: 360, after: 140 },
  children: [new TextRun({ text, font: HEAD, size: 30, bold: true, color: INK })],
});

const h2 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_2,
  spacing: { before: 240, after: 120 },
  children: [new TextRun({ text, font: HEAD, size: 24, bold: true, color: INK })],
});

const eyebrow = (text) => new Paragraph({
  spacing: { before: 260, after: 40 },
  children: [new TextRun({
    text: text.toUpperCase(), font: MONO, size: 16, color: MUTED, characterSpacing: 24,
  })],
});

function figure(n, file, caption) {
  const { w, h } = pngSize(file);
  return [
    new Paragraph({
      spacing: { before: 200, after: 80 },
      alignment: AlignmentType.CENTER,
      children: [new ImageRun({
        type: 'png',
        data: fs.readFileSync(file),
        transformation: { width: IMG_W, height: Math.round((h / w) * IMG_W) },
      })],
    }),
    new Paragraph({
      spacing: { after: 240 },
      children: [
        new TextRun({ text: `Figure ${n}. `, font: BODY, size: 19, bold: true, color: INK }),
        ...runs(caption, { size: 19, color: MUTED }),
      ],
    }),
  ];
}

function bullets(items) {
  return items.map((t) => new Paragraph({
    numbering: { reference: 'bul', level: 0 },
    spacing: { after: 110, line: 276 },
    children: runs(t),
  }));
}

function cell(children, width, opts = {}) {
  return new TableCell({
    width: { size: width, type: WidthType.DXA },
    margins: { top: 90, bottom: 90, left: 120, right: 120 },
    shading: opts.shade
      ? { type: ShadingType.CLEAR, fill: opts.shade, color: 'auto' }
      : undefined,
    children,
  });
}

function table(widths, header, rows) {
  const border = { style: BorderStyle.SINGLE, size: 4, color: RULE };
  const mk = (txt, width, o = {}) => cell(
    [new Paragraph({
      spacing: { after: 0, line: 252 },
      children: runs(txt, { size: o.size ?? 19, bold: o.bold, color: o.color,
                            font: o.font }),
    })],
    width,
    { shade: o.shade },
  );
  return new Table({
    width: { size: CONTENT, type: WidthType.DXA },
    columnWidths: widths,
    borders: {
      top: border, bottom: border, left: border, right: border,
      insideHorizontal: border, insideVertical: border,
    },
    rows: [
      new TableRow({
        tableHeader: true,
        children: header.map((t, i) => mk(t, widths[i], {
          bold: true, size: 17, color: MUTED, shade: SOFT,
        })),
      }),
      ...rows.map((r) => new TableRow({
        children: r.map((t, i) => mk(t, widths[i], i === 0 ? { font: MONO, bold: true } : {})),
      })),
    ],
  });
}

function callout(title, body) {
  return new Table({
    width: { size: CONTENT, type: WidthType.DXA },
    columnWidths: [CONTENT],
    borders: {
      top: { style: BorderStyle.NONE }, bottom: { style: BorderStyle.NONE },
      right: { style: BorderStyle.NONE },
      left: { style: BorderStyle.SINGLE, size: 18, color: WARN },
      insideHorizontal: { style: BorderStyle.NONE },
      insideVertical: { style: BorderStyle.NONE },
    },
    rows: [new TableRow({
      children: [cell([
        new Paragraph({
          spacing: { after: 90 },
          children: [new TextRun({ text: title, font: BODY, size: 21, bold: true, color: WARN })],
        }),
        new Paragraph({ spacing: { after: 0, line: 276 }, children: runs(body) }),
      ], CONTENT, { shade: WARNSOFT })],
    })],
  });
}

// ---------------------------------------------------------------------------
const doc = new Document({
  creator: 'TxV',
  title: 'P0-P8 routing panel: design rationale',
  description: 'Why each element of the nine pJET1.2 constructs is there.',
  numbering: {
    config: [{
      reference: 'bul',
      levels: [{
        level: 0, format: LevelFormat.BULLET, text: '•',
        alignment: AlignmentType.LEFT,
        style: { paragraph: { indent: { left: 360, hanging: 240 } } },
      }],
    }],
  },
  styles: {
    default: { document: { run: { font: BODY, size: 21, color: INK } } },
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
      },
    },
    children: [
      // ---- title ------------------------------------------------------
      new Paragraph({
        spacing: { after: 60 },
        children: [new TextRun({
          text: 'pJET1.2  ·  NINE CONSTRUCTS  ·  DESIGN RATIONALE',
          font: MONO, size: 16, color: MUTED, characterSpacing: 24,
        })],
      }),
      new Paragraph({
        spacing: { after: 140 },
        border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: RULE, space: 8 } },
        children: [new TextRun({
          text: 'Why these nine plasmids look the way they do',
          font: HEAD, size: 40, bold: true, color: INK,
        })],
      }),
      p('Every element in these constructs is there to answer a specific question. This is what each one does, and what we lose if it is missing.',
        { run: { size: 23, color: MUTED }, after: 300 }),

      // ---- 1 ----------------------------------------------------------
      eyebrow('01 — the problem'),
      h1('A plain antigen mRNA can only do half the job'),
      p('We want dendritic cells to present our antigens on **both** MHC class I (priming CD8 killers) and class II (recruiting CD4 help). Those two molecules are loaded in _different physical compartments_, and an mRNA delivered by LNP is translated on free ribosomes in the cytosol.'),
      p('So a bare antigen ORF is a CD8-biased vaccine _by construction_. The proteasome is right there in the cytosol; the class II loading compartment is a sealed vesicle the protein has no way into. This is the mechanism behind "expressed but not presented" — and it is why a western blot cannot diagnose it. The protein is made perfectly well. It is in the wrong room.'),
      ...figure(1, path.join(dir, 'fig-1.png'),
        'The compartments are the whole problem. Cytosolic translation feeds MHC I automatically (amber). The MIIC, where MHC II is loaded (teal), is a sealed lumen — a cytosolic protein cannot get in. Every routing module in this panel exists to build that missing arrow.'),

      // ---- 2 ----------------------------------------------------------
      eyebrow('02 — the fix'),
      h1('A signal peptide and a transmembrane domain are a matched pair'),
      p('You cannot push a protein into a membrane compartment; you make the cell’s own machinery carry it, and that needs two signals. Neither does anything useful alone, which is why every routed construct here carries both.'),
      p('The **signal peptide** gets the chain out of the cytosol: SRP grabs its hydrophobic core as it emerges, pauses translation, and threads the chain into the ER lumen. SRP reads the _first_ thing out of the ribosome — so a signal peptide only works at the extreme N-terminus. That single fact fixes module order in P7 and P8.'),
      p('But a signal peptide alone gives you a **secreted** protein, which leaves the cell and is useless to us. The **transmembrane domain** stops that. When the translocon hits ~20 hydrophobic residues, transfer halts and you get type I topology: the antigen cassette ends up in the lumen — and ER lumen becomes endosome lumen becomes lysosome lumen — while the tail stays in the cytosol, where the sorting machinery can read it.'),
      ...figure(2, path.join(dir, 'fig-2.png'),
        'Why both modules, always. A signal peptide alone gives a secreted protein — the one outcome that helps nobody. Only the pair leaves the cassette in the lumen that becomes the MIIC, with the sorting tail in the cytosol where adaptor proteins can read it. **P5 deliberately has neither** — it is the cytosolic baseline.'),

      // ---- 3 ----------------------------------------------------------
      eyebrow('03 — the choice under test'),
      h1('The tail is the postcode, and we are testing two of them'),
      p('With the cassette in the lumen and the tail in the cytosol, the tail’s short sorting motif decides where the whole protein goes. We are comparing two, and the difference is _the route_, not the destination — both end at the MIIC.'),
      ...figure(3, path.join(dir, 'fig-3.png'),
        'All nine constructs, drawn to scale and generated from the design library rather than by hand. Colour marks the routing module: amber for the CTLA-4 tail, teal for LAMP1. **P0 is scaffold only; P5 drops the routing modules entirely; P7 and P8 split into two separately routed proteins at P2A.**'),
      ...figure(4, path.join(dir, 'fig-4.png'),
        'Same destination, different journey. The CTLA-4 tail sends the protein to the cell surface first, where AP-2 and clathrin rapidly pull it back into endosomes; the LAMP1 tail takes an AP-3 shortcut straight from the Golgi. **P6 versus P0–P4 asks whether that detour via the surface helps or hurts.** We keep the CTLA-4 extracellular IgV domain out of the construct — that is the part that would bind CD80/CD86 and do real checkpoint signalling.'),

      h2('Why the module order is not ours to choose'),
      p('`GYQTI` only works as the _last five residues of the protein_. Append anything and the AP-3 adaptor stops binding, and the construct silently reroutes to the plasma membrane — which reads as a negative result, not as a mistake. `YVKM` has no such constraint: it sits internally with 19 native residues after it.'),
      p('That asymmetry is the entire reason P7 and P8 are built the way they are. Each carries two separate proteins from one mRNA, split by a **P2A** ribosomal skipping element — and the LAMP1 arm has to be the _downstream_ one, because only the downstream arm can end the polypeptide.'),
      ...figure(5, path.join(dir, 'fig-5.png'),
        'Three motifs whose function is positional. The failure mode is quiet: a construct with `GYQTI` buried still expresses, still shows an HA signal, and simply routes to the wrong place. The build checks this automatically — every LAMP1-bearing construct is asserted to end in `GYQTI`, and every non-P2A construct carrying E5 is asserted to end in `EEEEE`.'),

      // ---- 4 ----------------------------------------------------------
      new Paragraph({ children: [new PageBreak()] }),
      eyebrow('04 — the panel'),
      h1('What each construct is for'),
      p('Nine constructs, one variable moved at a time. P0 and P1/P2 are the controls that make the rest interpretable — without them, a weak signal from P3 could be the scaffold, the antigen, or the route, and we would not know which.'),
      table([900, 1300, 1500, 5660],
        ['Construct', 'Route', 'Payload', 'The question it answers'],
        [
          ['P0', 'CTLA-4', 'none', 'Does the scaffold alone give signal? Sets the floor.'],
          ['P1', 'CTLA-4', 'A', 'Single-antigen control — is A presented at all on this scaffold?'],
          ['P2', 'CTLA-4', 'B', 'Same for B. Needed before any two-antigen result means anything.'],
          ['P3', 'CTLA-4', 'A + B', 'The core question: can one transcript drive class I and class II presentation at once?'],
          ['P4', 'CTLA-4', 'B + A', 'Order control for P3 — does cassette position change the outcome?'],
          ['P5', 'cytosolic', 'A + B + E5', 'Baseline with no routing, and the only construct where the degron has a free C-terminus — so this is the degron-positive comparator.'],
          ['P6', 'LAMP1', 'A + B', 'Does the direct AP-3 route beat the detour via the cell surface?'],
          ['P7', 'dual', 'A cyto · B lyso', 'Can one mRNA send two antigens to two compartments deliberately? Note E5 is internal here, so it is expected to be weak.'],
          ['P8', 'dual', 'B cyto · A lyso', 'Swap control for P7 — is the effect about the route, or about the antigen?'],
        ]),

      // ---- 5 ----------------------------------------------------------
      eyebrow('05 — the vector'),
      h1('Every element in the plasmid, and what it buys us'),
      p('The backbone is pJET1.2 — a bacterial cloning vector with no eukaryotic machinery, because the plasmid only has to grow in _E. coli_ before being cut for IVT. The fragment goes into its single blunt Eco32I site at 371/372.'),
      ...figure(6, path.join(dir, 'fig-6.png'),
        'The Eco32I site sits inside the lethal `eco47IR` gene, so breaking it is the selection: colonies that grow carry an insert. The two flanking BglII sites do double duty — releasing the insert for the diagnostic gel and producing the IVT template in the same digest. **That is why no fragment may contain a BglII site**; all nine were screened and are clean.'),
      table([2200, 7160], ['Element', 'What it is for'], [
        ['T7, +1 = AGG', 'Not the textbook …CACTATAG. The transcript starts AGG, which is what an AG-initiating cap1 analog needs. The IVT kit manual gives the minimum promoter as 5′-TAATACGACTCACTATAAGG, so this is fixed by the chemistry — and it is the one thing never to disturb if the UTRs are ever changed.'],
        ['5′ UTR, 47 nt', 'Ends in the strong Kozak GCCACC and contains no upstream AUG. An upstream AUG would capture scanning ribosomes and quietly suppress the real ORF.'],
        ['Insert site', 'The blunt Eco32I cut at 371/372, unique in the vector and inside the lethal eco47IR gene. Breaking that gene is the selection: an empty re-ligated vector still kills its host, so colonies that grow carry an insert.'],
        ['3′ UTR, 296 nt', 'The BNT162b2 3′ UTR verbatim — the AES UTR followed by the mtRNR1 element, selected ex vivo for stability and protein output. It replaced a 99-nt sequence carried over from the old backbone, which the parts registry itself marked as a placeholder.'],
        ['poly(A), 120', 'Not in the gene fragment: a flat A120 fails vendor screening, and so does a segmented A30–linker–A70 whose 70-nt run is still too long. It is added by PCR with IVT_R_polyA120 before cloning, so the plasmid carries it encoded and the long primer is spent once rather than every prep.'],
        ['No eukaryotic promoter', 'pJET1.2 has none, which is the point of moving to it. The 809 bp of mammalian-only machinery the previous backbone carried did no work for an IVT template.'],
        ['HA and FLAG tags', 'These answer "was it made", on a different reagent from the pMHC antibodies that answer "was it presented". Without them, an absent pMHC signal has three indistinguishable causes: no translation, no processing, or wrong compartment.'],
        ['29-aa native flanks', 'Not padding. ERAP1 trims N-terminal extensions in the ER, but nothing trims a C-terminal one — so the proteasome has to make the epitope’s C-terminus exactly right, first time, and the surrounding residues are what set that cleavage. Minimal epitopes would gamble on it.'],
        ['E5 degron (EEEEE)', 'An acidic C-degron that accelerates proteasomal turnover, which is what MHC I peptide supply actually depends on. This is the element CVGBM carries and the one that performed best in CureVac\u2019s own routing screen. It is positional: a C-degron needs a free C-terminus, so it works in P0–P5 where it ends the protein, and is expected to be weak or inactive in P7/P8 where it sits upstream of P2A and the skip product has no free C-terminus. Its codons are pinned to GAAGAGGAAGAGGAG — five glutamates in a row are inherently repetitive, and alternating GAA/GAG breaks a periodicity that would otherwise create a 12-nt direct repeat.'],
        ['GGGGS linkers', 'Flexible spacers chosen to minimise junctional epitopes — the novel peptides created where one antigen is fused to the next, which exist nowhere in the tumour and compete for presentation with the real targets.'],
      ]),

      // ---- 6 ----------------------------------------------------------
      new Paragraph({ spacing: { before: 320, after: 0 }, children: [] }),
      callout('The result that will mislead you',
        'Membrane-routed constructs still present on MHC I — via ERAD pulling misfolded protein back into the cytosol, failed translocation, and defective ribosomal products. So **an MHC I signal is not evidence that the routing worked.** Only the class II readout distinguishes the routes. The same care applies to the degron: E5 only has a free C-terminus in P0–P5; in P7 and P8 it sits upstream of the P2A skip and the upstream product never gets one, so **P5 is the degron-positive comparator, not P7/P8.** Plan the conclusions on that basis before the first experiment, not after.'),

      // ---- 7 ----------------------------------------------------------
      eyebrow('06 — not yet tested'),
      h1('What this panel deliberately does not answer'),
      p('This panel settles _where_ an antigen goes. Two questions it cannot reach, both of which decide whether a multi-antigen vaccine is actually potent:'),
      ...bullets([
        '**Is there a presentation budget?** Everything here carries one or two antigens. Whether per-epitope presentation decays as you add a fourth, sixth or eighth — dilution rather than breadth — is open, and it decides how many antigens a cassette should carry.',
        '**Does uridine depletion trade expression against adjuvanticity?** Lowering uridine raises expression and lowers innate sensing. For a prophylactic vaccine that is pure gain. For a cancer vaccine it may not be — the sensing being suppressed is part of what matures the dendritic cell and licenses it to prime CD8 T cells.',
      ]),

      new Paragraph({
        spacing: { before: 360 },
        border: { top: { style: BorderStyle.SINGLE, size: 6, color: RULE, space: 10 } },
        children: runs('Constructs P0–P8 on pJET1.2 (2974 bp backbone; assembled plasmids 3644–3929 bp). Nine gene fragments, 670–955 bp. Every map verified against the parent vector — reading frame, Kozak junction, AGG transcription start, the IVT kit minimum promoter, BglII count and the disruption of eco47IR re-derived from the assembled sequence rather than assumed.',
          { size: 18, color: MUTED }),
      }),
    ],
  }],
});

Packer.toBuffer(doc).then((buf) => {
  const out = path.join(dir, 'P0-P8_design_rationale.docx');
  fs.writeFileSync(out, buf);
  console.log('wrote', out, buf.length, 'bytes');
});
