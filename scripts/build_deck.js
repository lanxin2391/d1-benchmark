// Build the 1118 defence deck from the two content files.
//
//   node scripts/build_deck.js
//
// Design rules taken from the course requirements: 16:9, white/light
// background, Times New Roman, English titles ~42 pt, body ~24 pt
// (secondary text smaller), a left-aligned title, one take-home
// message in the footer of every slide, figures labelled with their
// source, and every slide carrying a visual element.
const fs = require("fs");
const path = require("path");
// pptxgenjs is installed globally (npm i -g pptxgenjs), and Node does not
// resolve global modules from a plain require(), so fall back to the
// global root.
let PptxGenJS;
try {
  PptxGenJS = require("pptxgenjs");
} catch (e) {
  const root = process.env.APPDATA + "\\npm\\node_modules";
  PptxGenJS = require(path.join(root, "pptxgenjs"));
}

const ROOT = path.resolve(__dirname, "..");
const OUT = "D:/Grade3/swxxx/final/D1/1118汇报/D1_defence_deck.pptx";
const FONT = "Times New Roman";
const INK = "2A2723";
const MUTED = "8C877D";
const BLUE = "2A78D6";
const AQUA = "1BAF7A";
const VIOLET = "4A3AA7";
const ORANGE = "EB6834";
const GREEN = "008300";
const RED = "E34948";
const BAND = "F7F5F0";
const HAIR = "D8D3C8";
const M = 0.55;                 // the side margin
const W = 13.333, H = 7.5;
const CW = W - 2 * M;           // the content width

const SECTION = {
  1: GREEN, 2: BLUE, 3: VIOLET, 4: VIOLET, 5: VIOLET, 6: VIOLET,
  7: AQUA, 8: BLUE, 9: BLUE, 10: BLUE, 11: ORANGE, 12: ORANGE,
  13: ORANGE, 14: GREEN, 15: GREEN, 16: GREEN,
};
// the appendix slides are grey: they are evidence, not narrative
const APPENDIX_GREY = "6B675E";

function pngSize(file) {
  const b = fs.readFileSync(file);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) };
}

function titleSize(t) {
  if (t.length <= 42) return 42;
  if (t.length <= 50) return 36;
  if (t.length <= 62) return 32;
  return 28;
}

function readJSON(p) {
  return JSON.parse(fs.readFileSync(p, "utf8"));
}

const part1 = readJSON(path.join(ROOT, "docs/deck/slides_01_08.json"));
const part2 = readJSON(path.join(ROOT, "docs/deck/slides_09_16.json"));
const appx = readJSON(path.join(ROOT, "docs/deck/slides_appendix.json"));
const slides = part1.slides.concat(part2.slides).concat(appx.slides);
const meta = part1.meta;

const pptx = new PptxGenJS();
pptx.defineLayout({ name: "W16x9", width: W, height: H });
pptx.layout = "W16x9";
pptx.title = meta.title;
pptx.subject = meta.subject;
pptx.author = meta.author;
pptx.company = meta.company;

function chrome(s, sl) {
  const colour = sl.kind === "appendix" ? APPENDIX_GREY : (SECTION[sl.n] || BLUE);
  s.background = { color: "FFFFFF" };
  if (sl.tag) {
    s.addText(sl.tag, {
      x: M, y: 0.26, w: CW, h: 0.26, fontFace: FONT, fontSize: 12.5,
      bold: true, color: colour, charSpacing: 2, align: "left",
      margin: 0,
    });
  }
  if (sl.title) {
    s.addText(sl.title, {
      x: M, y: 0.55, w: CW, h: 0.85, fontFace: FONT,
      fontSize: titleSize(sl.title), bold: true, color: INK,
      align: "left", valign: "top", margin: 0,
    });
  }
  // the footer band and the take-home message
  s.addShape(pptx.ShapeType.rect, {
    x: 0, y: 6.82, w: W, h: 0.68, fill: { color: BAND }, line: { type: "none" },
  });
  s.addShape(pptx.ShapeType.rect, {
    x: 0, y: 6.82, w: 0.22, h: 0.68, fill: { color: colour }, line: { type: "none" },
  });
  if (sl.footer) {
    s.addText(sl.footer, {
      x: M, y: 6.90, w: CW - 0.7, h: 0.52, fontFace: FONT, fontSize: 13.5,
      italic: true, color: "4A4741", align: "left", valign: "middle",
      margin: 0,
    });
  }
  s.addText(String(sl.n), {
    x: W - 0.85, y: 6.90, w: 0.5, h: 0.52, fontFace: FONT, fontSize: 13,
    color: MUTED, align: "right", valign: "middle", margin: 0,
  });
}

function image(s, rel, width, top, left) {
  const file = path.join(ROOT, rel);
  const dim = pngSize(file);
  const w = width;
  const h = (w * dim.h) / dim.w;
  const x = left === undefined ? (W - w) / 2 : left;
  s.addImage({ path: file, x, y: top, w, h });
  return h;
}

function bullets(s, items, opt) {
  const o = Object.assign({ x: M, y: 1.5, w: CW, h: 4.6, size: 19 }, opt || {});
  const runs = [];
  items.forEach((t, i) => {
    runs.push({
      text: t,
      options: {
        bullet: { code: "2022" }, breakLine: true, paraSpaceAfter: 10,
      },
    });
  });
  s.addText(runs, {
    x: o.x, y: o.y, w: o.w, h: o.h, fontFace: FONT, fontSize: o.size,
    color: INK, align: "left", valign: "top", lineSpacingMultiple: 1.02,
    margin: 0,
  });
}

// ---------------------------------------------------------------- build
slides.forEach((sl) => {
  const s = pptx.addSlide();

  if (sl.kind === "title") {
    s.background = { color: "FFFFFF" };
    s.addShape(pptx.ShapeType.rect, {
      x: 0, y: 0, w: W, h: 0.30, fill: { color: GREEN }, line: { type: "none" },
    });
    s.addText(sl.tag, {
      x: M, y: 0.85, w: CW, h: 0.3, fontFace: FONT, fontSize: 13,
      bold: true, color: GREEN, charSpacing: 2, margin: 0,
    });
    s.addText(sl.title, {
      x: M, y: 1.25, w: CW, h: 1.5, fontFace: FONT, fontSize: 42, bold: true,
      color: INK, align: "left", valign: "top", margin: 0,
    });
    s.addText(sl.subtitle, {
      x: M, y: 2.85, w: CW - 1.0, h: 0.9, fontFace: FONT, fontSize: 22,
      color: "4A4741", align: "left", valign: "top", margin: 0,
    });
    const rows = sl.members.map((m) => [
      { text: m.split("—")[0].trim(), options: { bold: true, color: GREEN } },
      { text: m.split("—")[1].trim(), options: { color: INK } },
    ]);
    s.addTable(rows, {
      x: M, y: 4.05, w: CW - 1.0, colW: [2.4, 8.8], fontFace: FONT,
      fontSize: 15, border: { type: "none" }, rowH: 0.32, valign: "middle",
      margin: [0, 0.04, 0, 0],
    });
    s.addText(sl.meta, {
      x: M, y: 5.75, w: CW - 1.0, h: 0.6, fontFace: FONT, fontSize: 13,
      color: MUTED, align: "left", valign: "top", margin: 0,
    });
    s.addShape(pptx.ShapeType.rect, {
      x: 0, y: 6.82, w: W, h: 0.68, fill: { color: BAND }, line: { type: "none" },
    });
    s.addShape(pptx.ShapeType.rect, {
      x: 0, y: 6.82, w: 0.22, h: 0.68, fill: { color: GREEN }, line: { type: "none" },
    });
    s.addText(sl.footer, {
      x: M, y: 6.90, w: CW - 0.7, h: 0.52, fontFace: FONT, fontSize: 13.5,
      italic: true, color: "4A4741", align: "left", valign: "middle", margin: 0,
    });
    return;
  }

  if (sl.kind === "image") {
    chrome(s, sl);
    const h = image(s, sl.image, sl.imagew, 1.55);
    return;
  }

  if (sl.kind === "papers") {
    chrome(s, sl);
    s.addText(sl.paper.cite, {
      x: M, y: 1.42, w: CW, h: 0.42, fontFace: FONT, fontSize: 15,
      italic: true, color: VIOLET, align: "left", valign: "top", margin: 0,
    });
    const rows = sl.paper.rows.map((r) => [
      { text: r[0], options: { bold: true, color: VIOLET, fill: { color: "F2F0FA" } } },
      { text: r[1], options: { color: INK } },
    ]);
    s.addTable(rows, {
      x: M, y: 1.95, w: CW, colW: [2.5, 9.73], fontFace: FONT, fontSize: 17,
      border: { type: "solid", color: HAIR, pt: 0.75 },
      rowH: 0.95, valign: "middle", margin: [0.04, 0.08, 0.04, 0.08],
    });
    return;
  }

  if (sl.kind === "hypotheses") {
    chrome(s, sl);
    s.addText(sl.gap, {
      x: M, y: 1.45, w: CW, h: 1.25, fontFace: FONT, fontSize: 18,
      color: "4A4741", align: "left", valign: "top", margin: 0,
    });
    sl.hyp.forEach((h, i) => {
      const y = 2.85 + i * 1.28;
      const colour = [BLUE, AQUA, ORANGE][i];
      s.addShape(pptx.ShapeType.roundRect, {
        x: M, y, w: CW, h: 1.12, fill: { color: "FFFFFF" },
        line: { color: colour, width: 1.25 }, rectRadius: 0.06,
      });
      s.addText(h[0], {
        x: M + 0.18, y: y + 0.08, w: 1.1, h: 0.4, fontFace: FONT,
        fontSize: 19, bold: true, color: colour, margin: 0,
      });
      s.addText(h[1], {
        x: M + 1.35, y: y + 0.08, w: CW - 1.6, h: 0.4, fontFace: FONT,
        fontSize: 18, bold: true, color: INK, margin: 0,
      });
      s.addText(h[2], {
        x: M + 1.35, y: y + 0.52, w: CW - 1.6, h: 0.52, fontFace: FONT,
        fontSize: 15, color: "4A4741", margin: 0,
      });
    });
    return;
  }

  if (sl.kind === "table") {
    chrome(s, sl);
    s.addText(sl.lead, {
      x: M, y: 1.42, w: CW, h: 0.42, fontFace: FONT, fontSize: 18,
      color: "4A4741", align: "left", valign: "top", margin: 0,
    });
    const head = sl.table.head.map((h) => ({
      text: h, options: { bold: true, color: "FFFFFF", fill: { color: AQUA } },
    }));
    const body = sl.table.rows.map((r, i) => r.map((c, j) => ({
      text: c,
      options: {
        color: INK,
        fill: { color: i % 2 ? "FFFFFF" : "F4F8F6" },
        bold: j === 0,
      },
    })));
    s.addTable([head].concat(body), {
      x: M, y: 1.95, w: CW, colW: [2.5, 5.4, 1.7, 1.63], fontFace: FONT,
      fontSize: 15, border: { type: "solid", color: HAIR, pt: 0.75 },
      rowH: 0.42, valign: "middle", margin: [0.03, 0.07, 0.03, 0.07],
    });
    sl.table.rowlabels.forEach((lab, i) => {
      s.addText([
        { text: lab[0] + "   ", options: { bold: true, color: AQUA } },
        { text: lab[1], options: { color: INK } },
      ], {
        x: M, y: 5.05 + i * 0.78, w: CW, h: 0.72, fontFace: FONT,
        fontSize: 15, align: "left", valign: "top", margin: 0,
      });
    });
    return;
  }

  if (sl.kind === "pipeline") {
    chrome(s, sl);
    const cols = 3, rows = 2;
    const bw = (CW - 0.4 * (cols - 1)) / cols;
    const bh = 2.25;
    sl.steps.forEach((st, i) => {
      const c = i % cols, r = Math.floor(i / cols);
      const x = M + c * (bw + 0.4);
      const y = 1.55 + r * (bh + 0.35);
      s.addShape(pptx.ShapeType.roundRect, {
        x, y, w: bw, h: bh, fill: { color: "FFFFFF" },
        line: { color: BLUE, width: 1.1 }, rectRadius: 0.06,
      });
      s.addShape(pptx.ShapeType.ellipse, {
        x: x + 0.16, y: y + 0.16, w: 0.42, h: 0.42, fill: { color: BLUE },
        line: { type: "none" },
      });
      s.addText(String(i + 1), {
        x: x + 0.16, y: y + 0.16, w: 0.42, h: 0.42, fontFace: FONT,
        fontSize: 15, bold: true, color: "FFFFFF", align: "center",
        valign: "middle", margin: 0,
      });
      s.addText(st[0], {
        x: x + 0.7, y: y + 0.18, w: bw - 0.85, h: 0.38, fontFace: FONT,
        fontSize: 18, bold: true, color: INK, margin: 0,
      });
      s.addText(st[1], {
        x: x + 0.18, y: y + 0.72, w: bw - 0.36, h: bh - 0.85, fontFace: FONT,
        fontSize: 13.5, color: "4A4741", align: "left", valign: "top",
        margin: 0, lineSpacingMultiple: 1.0,
      });
    });
    return;
  }

  if (sl.kind === "params") {
    chrome(s, sl);
    const head = ["Parameter", "What it is", "Why this choice", "How we checked it"]
      .map((h) => ({ text: h, options: { bold: true, color: "FFFFFF", fill: { color: BLUE } } }));
    const body = sl.params.map((r, i) => r.map((c, j) => ({
      text: c,
      options: {
        color: INK, bold: j === 0, fill: { color: i % 2 ? "FFFFFF" : "F1F6FC" },
      },
    })));
    s.addTable([head].concat(body), {
      x: M, y: 1.5, w: CW, colW: [2.5, 2.6, 4.0, 3.13], fontFace: FONT,
      fontSize: 13.5, border: { type: "solid", color: HAIR, pt: 0.75 },
      rowH: 0.62, valign: "middle", margin: [0.03, 0.07, 0.03, 0.07],
    });
    s.addShape(pptx.ShapeType.roundRect, {
      x: M, y: 5.35, w: CW, h: 1.15, fill: { color: "FBF6EE" },
      line: { color: ORANGE, width: 1.1 }, rectRadius: 0.05,
    });
    s.addText(sl.honesty, {
      x: M + 0.2, y: 5.45, w: CW - 0.4, h: 0.95, fontFace: FONT, fontSize: 14.5,
      color: INK, align: "left", valign: "middle", margin: 0,
    });
    return;
  }

  if (sl.kind === "controls") {
    chrome(s, sl);
    const head = ["Control", "The artefact it rules out", "What we measured"]
      .map((h) => ({ text: h, options: { bold: true, color: "FFFFFF", fill: { color: BLUE } } }));
    const body = sl.controls.map((r, i) => r.map((c, j) => ({
      text: c,
      options: {
        color: INK, bold: j === 0, fill: { color: i % 2 ? "FFFFFF" : "F1F6FC" },
      },
    })));
    s.addTable([head].concat(body), {
      x: M, y: 1.45, w: CW, colW: [3.6, 3.6, 5.03], fontFace: FONT,
      fontSize: 12.5, border: { type: "solid", color: HAIR, pt: 0.75 },
      rowH: 0.44, valign: "middle", margin: [0.02, 0.07, 0.02, 0.07],
    });
    return;
  }

  if (sl.kind === "heatmap") {
    chrome(s, sl);
    const h = image(s, sl.image, sl.imagew, 1.42);
    const y0 = 1.42 + h + 0.18;
    sl.stats.forEach((st, i) => {
      const x = M + i * 4.05;
      s.addText(st[0], {
        x, y: y0, w: 3.8, h: 0.55, fontFace: FONT, fontSize: 30,
        bold: true, color: ORANGE, align: "left", margin: 0,
      });
      s.addText(st[1], {
        x, y: y0 + 0.53, w: 3.8, h: 0.45, fontFace: FONT, fontSize: 12.5,
        color: MUTED, align: "left", margin: 0,
      });
    });
    return;
  }

  if (sl.kind === "stability") {
    chrome(s, sl);
    const h = image(s, sl.image, 11.4, 1.45);
    bullets(s, sl.bullets, {
      x: M, y: 1.45 + h + 0.12, w: CW, h: 6.7 - (1.45 + h), size: 13,
    });
    return;
  }

  if (sl.kind === "external") {
    chrome(s, sl);
    image(s, sl.image, 6.35, 1.75, M);
    bullets(s, sl.bullets, { x: 7.15, y: 1.55, w: W - 7.15 - M, h: 5.0, size: 12.5 });
    return;
  }

  if (sl.kind === "conclusions") {
    chrome(s, sl);
    sl.points.forEach((p, i) => {
      const y = 1.5 + i * 1.02;
      s.addShape(pptx.ShapeType.ellipse, {
        x: M, y: y + 0.04, w: 0.5, h: 0.5, fill: { color: GREEN },
        line: { type: "none" },
      });
      s.addText(p[0], {
        x: M, y: y + 0.04, w: 0.5, h: 0.5, fontFace: FONT, fontSize: 17,
        bold: true, color: "FFFFFF", align: "center", valign: "middle", margin: 0,
      });
      s.addText(p[1], {
        x: M + 0.68, y, w: CW - 0.68, h: 0.42, fontFace: FONT, fontSize: 19,
        bold: true, color: INK, margin: 0,
      });
      s.addText(p[2], {
        x: M + 0.68, y: y + 0.42, w: CW - 0.68, h: 0.5, fontFace: FONT,
        fontSize: 14.5, color: "4A4741", margin: 0,
      });
    });
    return;
  }

  if (sl.kind === "limits") {
    chrome(s, sl);
    const head = ["Limitation", "Where it bites", "The experiment that would settle it"]
      .map((h) => ({ text: h, options: { bold: true, color: "FFFFFF", fill: { color: GREEN } } }));
    const body = sl.limits.map((r, i) => r.map((c, j) => ({
      text: c,
      options: {
        color: INK, bold: j === 0, fill: { color: i % 2 ? "FFFFFF" : "F1F7F1" },
      },
    })));
    s.addTable([head].concat(body), {
      x: M, y: 1.5, w: CW, colW: [3.9, 3.3, 5.03], fontFace: FONT,
      fontSize: 14, border: { type: "solid", color: HAIR, pt: 0.75 },
      rowH: 0.82, valign: "middle", margin: [0.03, 0.08, 0.03, 0.08],
    });
    return;
  }

  if (sl.kind === "appendix") {
    sl.footer = sl.footer || ("Appendix: the evidence behind the claims -- "
      + "every full table is in the repository "
      + "(results/tables/, with SHA-256 manifests).");
    chrome(s, sl);
    if (sl.note) {
      s.addText(sl.note, {
        x: M, y: 1.38, w: CW, h: 0.5, fontFace: FONT, fontSize: 12,
        italic: true, color: MUTED, align: "left", valign: "top", margin: 0,
      });
    }
    const tb = sl.table;
    const head = tb.head.map((h) => ({
      text: h, options: { bold: true, color: "FFFFFF", fill: { color: APPENDIX_GREY } },
    }));
    const body = tb.rows.map((r, i) => r.map((c, j) => ({
      text: String(c),
      options: {
        color: INK, bold: j === 0, fill: { color: i % 2 ? "FFFFFF" : "F5F3EE" },
      },
    })));
    const ncol = tb.head.length;
    const first = ncol >= 4 ? 3.1 : 3.9;
    const rest = (CW - first) / (ncol - 1);
    const colW = tb.head.map((_, j) => (j === 0 ? first : rest));
    // fit the table into the space between the note and the footer
    const avail = 4.55;
    const rowH = Math.min(tb.rowh || 0.42, avail / (tb.rows.length + 1));
    const fs = Math.min(tb.fontsize || 12, Math.max(8.5, (rowH * 72) / 2.7));
    s.addTable([head].concat(body), {
      x: M, y: 1.95, w: CW, colW, fontFace: FONT, fontSize: fs,
      border: { type: "solid", color: HAIR, pt: 0.75 },
      rowH, valign: "middle", margin: [0.02, 0.06, 0.02, 0.06],
    });
    return;
  }

  if (sl.kind === "closing") {
    chrome(s, sl);
    const rows = sl.links.map((r) => [
      { text: r[0], options: { bold: true, color: GREEN } },
      { text: r[1], options: { color: INK } },
    ]);
    s.addTable(rows, {
      x: M, y: 1.5, w: CW, colW: [3.0, 9.23], fontFace: FONT, fontSize: 14.5,
      border: { type: "none" }, rowH: 0.36, valign: "middle", margin: [0, 0.04, 0, 0],
    });
    s.addShape(pptx.ShapeType.roundRect, {
      x: M, y: 3.5, w: CW, h: 1.62, fill: { color: "F7F5F0" },
      line: { color: HAIR, width: 0.75 }, rectRadius: 0.05,
    });
    s.addText(sl.ai, {
      x: M + 0.22, y: 3.62, w: CW - 0.44, h: 1.4, fontFace: FONT,
      fontSize: 14.5, color: INK, align: "left", valign: "top", margin: 0,
      lineSpacingMultiple: 1.05,
    });
    s.addText(sl.roles, {
      x: M, y: 5.65, w: CW, h: 0.7, fontFace: FONT, fontSize: 14.5,
      color: MUTED, align: "left", valign: "top", margin: 0,
    });
    return;
  }
});

fs.mkdirSync(path.dirname(OUT), { recursive: true });
pptx.writeFile({ fileName: OUT }).then(() => {
  console.log("wrote " + OUT + " (" + slides.length + " slides)");
});
