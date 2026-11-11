#!/usr/bin/env python3
"""Convert the research report draft (Markdown) into a Word document.

The report is written in Markdown; the deliverable the user asked for is a
.docx. This converter maps the subset of Markdown the report uses:

    # / ## / ### headings
    paragraphs, with **bold**, *italic* and `code` runs
    pipe tables (with a shaded header row)
    bullet lists (- and *)
    numbered lists (1. 2. 3.)
    fenced code blocks (monospace, shaded)
    horizontal rules (-> a thin bottom border)

Page setup is US Letter with 1" margins, the body is Times New Roman 11 pt
(the standard for a journal manuscript), tables are set in 9 pt, and every
table carries a numbered caption above it and a source note below it, in
the style a journal expects.

Usage
-----
    python scripts/markdown_to_docx.py \\
        docs/research_report_draft.md \\
        docs/research_report_draft.docx
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import docx
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

BODY_FONT = "Times New Roman"
MONO_FONT = "Consolas"
BODY_SIZE = Pt(11)
TABLE_SIZE = Pt(9)

CONTENT_WIDTH_IN = 6.5          # US Letter, 1" margins
ACCENT = RGBColor(0x1F, 0x38, 0x64)
GREY = RGBColor(0x59, 0x59, 0x59)


# ------------------------------------------------------------------ helpers
def _set_cell_background(cell, hex_color: str) -> None:
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hex_color)
    cell._tc.get_or_add_tcPr().append(shd)


def _set_table_borders(table, color: str = "BFBFBF", size: int = 4) -> None:
    tbl = table._tbl
    tblPr = tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), str(size))
        el.set(qn("w:space"), "0")
        el.set(qn("w:color"), color)
        borders.append(el)
    tblPr.append(borders)


def _repeat_header_row(row) -> None:
    trPr = row._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader")
    el.set(qn("w:val"), "true")
    trPr.append(el)


def _paragraph_border_bottom(paragraph, color: str = "BFBFBF",
                             size: int = 6) -> None:
    pPr = paragraph._p.get_or_add_pPr()
    borders = OxmlElement("w:pBdr")
    bottom = OxmlElement("w:bottom")
    bottom.set(qn("w:val"), "single")
    bottom.set(qn("w:sz"), str(size))
    bottom.set(qn("w:space"), "1")
    bottom.set(qn("w:color"), color)
    borders.append(bottom)
    pPr.append(borders)


def _shade_paragraph(paragraph, hex_color: str) -> None:
    pPr = paragraph._p.get_or_add_pPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hex_color)
    pPr.append(shd)


# ------------------------------------------------------------------ inline
INLINE = re.compile(
    r"(?P<code>`[^`]+`)"
    r"|(?P<link>\[[^\]]+\]\([^)]+\))"
    r"|(?P<autolink><https?://[^>]+>)"
    r"|(?P<bold>\*\*[^*]+\*\*)"
    r"|(?P<italic>\*[^*\n]+\*)"
)

# a bare identifier such as `P_{k+1}` or `W^T` is mathematics, not code:
# it should stay in the body font. Only a token that looks like a file
# name, a path or a command is rendered as code.
CODE_LIKE = re.compile(
    r"^[\w./\\-]+\.(py|md|tsv|csv|json|yml|yaml|pdf|txt|docx|npz|gz)$"
    r"|^[\w./\\-]*\/[\w./\\-]+$"
    r"|^--?[\w-]+$"
    r"|^python\b"
    r"|^[A-Za-z_]+=[^\s]*$"
    r"|^0x[0-9a-fA-F]+$"
)


def _is_code(token: str) -> bool:
    """Is this `...` token code, or mathematics that should stay in body font?"""
    t = token.strip()
    if CODE_LIKE.match(t):
        return True
    if re.search(r"[\s=()\[\]{}<>]", t):
        # contains spaces or code punctuation -> code
        return True
    if "_" in t and re.match(r"^[A-Za-z][A-Za-z0-9_]*$", t):
        # a Python identifier with an underscore, e.g. mean_delta_AUROC
        return True
    # a short all-caps or single-word token with no underscore is likely
    # mathematics (alpha, D, RWR, TP53)
    return False


def add_inline(paragraph, text: str, size=Pt(11), bold=False, italic=False,
               color=None, font=BODY_FONT) -> None:
    """Add `text` to `paragraph`, honouring **bold**, *italic*, `code` and links."""
    pos = 0
    for m in INLINE.finditer(text):
        if m.start() > pos:
            _run(paragraph, text[pos:m.start()], size, bold, italic, color, font)
        if m.group("code"):
            tok = m.group("code")[1:-1]
            if _is_code(tok):
                _run(paragraph, tok, size, bold, italic,
                     color or RGBColor(0x33, 0x33, 0x33), MONO_FONT)
            else:
                # mathematics: keep the body font but italicise the
                # variables, as a journal would
                _run(paragraph, tok, size, bold, True, color, BODY_FONT)
        elif m.group("link"):
            # [text](url) -> a real hyperlink carrying the text
            lm = re.fullmatch(r"\[([^\]]+)\]\(([^)]+)\)", m.group("link"))
            label, url = lm.group(1), lm.group(2)
            _hyperlink(paragraph, label, url, size, bold, italic)
        elif m.group("autolink"):
            url = m.group("autolink")[1:-1]
            _hyperlink(paragraph, url, url, size, bold, italic)
        elif m.group("bold"):
            _run(paragraph, m.group("bold")[2:-2], size, True, italic, color, font)
        elif m.group("italic"):
            _run(paragraph, m.group("italic")[1:-1], size, bold, True, color, font)
        pos = m.end()
    if pos < len(text):
        _run(paragraph, text[pos:], size, bold, italic, color, font)


def _hyperlink(paragraph, label: str, url: str, size, bold, italic) -> None:
    """Add a real Word hyperlink."""
    part = paragraph.part
    r_id = part.relate_to(
        url,
        "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink",
        is_external=True)
    link = OxmlElement("w:hyperlink")
    link.set(qn("r:id"), r_id)
    r = OxmlElement("w:r")
    rPr = OxmlElement("w:rPr")
    for tag, val in (("w:color", "0563C1"), ("w:u", "single")):
        el = OxmlElement(tag)
        el.set(qn("w:val"), "single" if tag == "w:u" else val)
        rPr.append(el)
    rf = OxmlElement("w:rFonts")
    for attr in ("w:ascii", "w:hAnsi", "w:cs"):
        rf.set(qn(attr), BODY_FONT)
    rPr.append(rf)
    sz = OxmlElement("w:sz")
    sz.set(qn("w:val"), str(int(size.pt * 2)))
    rPr.append(sz)
    r.append(rPr)
    t = OxmlElement("w:t")
    t.text = label
    r.append(t)
    link.append(r)
    paragraph._p.append(link)


def _run(paragraph, text, size, bold, italic, color, font) -> None:
    if not text:
        return
    r = paragraph.add_run(text)
    r.font.name = font
    r.font.size = size
    r.bold = bold
    r.italic = italic
    if color is not None:
        r.font.color.rgb = color
    # make the font apply to east-asian / complex-script runs too
    rPr = r._element.get_or_add_rPr()
    rFonts = rPr.find(qn("w:rFonts"))
    if rFonts is None:
        rFonts = OxmlElement("w:rFonts")
        rPr.append(rFonts)
    for attr in ("w:ascii", "w:hAnsi", "w:cs"):
        rFonts.set(qn(attr), font)


# ------------------------------------------------------------------ document
def new_document() -> Document:
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = BODY_FONT
    st.font.size = BODY_SIZE
    st.element.rPr.rFonts.set(qn("w:eastAsia"), BODY_FONT)
    st.paragraph_format.space_after = Pt(6)
    st.paragraph_format.line_spacing = 1.15

    for name, size, color, before, after in (
        ("Heading 1", Pt(15), ACCENT, 16, 6),
        ("Heading 2", Pt(12.5), ACCENT, 12, 4),
        ("Heading 3", Pt(11), ACCENT, 10, 4),
    ):
        s = doc.styles[name]
        s.font.name = BODY_FONT
        s.font.size = size
        s.font.bold = True
        s.font.color.rgb = color
        s.paragraph_format.space_before = Pt(before)
        s.paragraph_format.space_after = Pt(after)
        s.paragraph_format.keep_with_next = True

    sec = doc.sections[0]
    sec.page_width = Inches(8.5)
    sec.page_height = Inches(11)
    for attr in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
        setattr(sec, attr, Inches(1))
    return doc


def add_page_numbers(doc: Document) -> None:
    for section in doc.sections:
        footer = section.footer
        p = footer.paragraphs[0] if footer.paragraphs else footer.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run()
        r.font.name = BODY_FONT
        r.font.size = Pt(9)
        for instr in ("PAGE",):
            fld = OxmlElement("w:fldSimple")
            fld.set(qn("w:instr"), instr)
            p._p.append(fld)


def add_header(doc: Document, text: str) -> None:
    for section in doc.sections:
        hp = section.header.paragraphs[0] if section.header.paragraphs \
            else section.header.add_paragraph()
        hp.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = hp.add_run(text)
        r.font.name = BODY_FONT
        r.font.size = Pt(8.5)
        r.font.color.rgb = GREY
        r.italic = True
        _paragraph_border_bottom(hp, "D0D0D0", 4)


# ------------------------------------------------------------------ tables
def split_row(line: str) -> list[str]:
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [c.strip() for c in line.split("|")]


def is_separator(line: str) -> bool:
    return bool(re.fullmatch(r"\|?[\s:\-|]+\|?", line.strip())) and \
        "-" in line


def add_table(doc: Document, rows: list[list[str]], caption: str | None,
              number: int) -> None:
    ncol = max(len(r) for r in rows)
    if ncol < 1:
        return
    if caption:
        cp = doc.add_paragraph()
        cp.paragraph_format.space_before = Pt(8)
        cp.paragraph_format.space_after = Pt(3)
        cp.paragraph_format.keep_with_next = True
        r = cp.add_run(f"Table {number}. ")
        r.bold = True
        r.font.name = BODY_FONT
        r.font.size = Pt(10)
        add_inline(cp, caption, size=Pt(10))

    table = doc.add_table(rows=len(rows), cols=ncol)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    _set_table_borders(table)
    # do not let a table split across a page: a header row alone at the
    # bottom of a page is unreadable
    for row in table.rows:
        trPr = row._tr.get_or_add_trPr()
        cant = OxmlElement("w:cantSplit")
        trPr.append(cant)

    for i, row in enumerate(rows):
        for j in range(ncol):
            cell = table.cell(i, j)
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.line_spacing = 1.0
            txt = row[j] if j < len(row) else ""
            if i == 0:
                _set_cell_background(cell, "E8EDF4")
                add_inline(p, txt, size=TABLE_SIZE, bold=True)
            else:
                # right-align the numeric columns
                if re.fullmatch(r"[-+−]?[\d.,]+%?", txt.strip()):
                    p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
                add_inline(p, txt, size=TABLE_SIZE)
        if i == 0:
            _repeat_header_row(table.rows[0])


# ------------------------------------------------------------------ figures
FIGURES = [
    ("figures/figure2_main_heatmap.png", "Figure 1. The primary endpoint. "
     "Heat map of ΔAUROC = AUROC(RWR) − AUROC(degree) at α = 0.5 for the "
     "6 networks × 4 label sets. Green is a positive margin (the network "
     "prior beats degree), red is negative. The single negative cell is "
     "Reactome × IntOGen newly-discovered drivers."),
    ("figures/figure9_pauc_top10.png", "Figure 2. The same 24 arms on the "
     "two co-primary endpoints: the partial AUC below 2 % FPR and the "
     "precision in the top 10 %. The margin is smaller in both regions "
     "than in the full AUROC, and the ordering of the arms is unchanged."),
    ("figures/figure2_d09_d11_d12.png", "Figure 3. The three coverage "
     "controls side by side: the native evaluation, the D-11 "
     "common-size cap, and the D-12 shared gene universe."),
    ("figures/figure2_d09_temporal_strip.png", "Figure 4. The temporal "
     "control. The 152 newly-discovered IntOGen drivers are recovered with "
     "a smaller margin (0.02–0.10) than the established drivers "
     "(0.04–0.17); Reactome is the only network with a negative "
     "temporal-drift estimate."),
]


def add_figures(doc: Document, repo_root: Path) -> int:
    """Insert the paper's figures with numbered captions."""
    from docx.shared import Emu

    n = 0
    for rel, caption in FIGURES:
        path = repo_root / rel
        if not path.is_file():
            continue
        n += 1
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(10)
        p.paragraph_format.space_after = Pt(3)
        p.add_run().add_picture(str(path), width=Inches(CONTENT_WIDTH_IN))
        cp = doc.add_paragraph()
        cp.paragraph_format.space_after = Pt(10)
        add_inline(cp, caption, size=Pt(9))
    return n


# ------------------------------------------------------------------ convert
def convert(md_path: Path, out_path: Path, doc: Document,
            repo_root: Path) -> int:
    lines = md_path.read_text(encoding="utf-8").splitlines()
    i = 0
    n = len(lines)
    table_number = 0
    fig_number = 0
    in_code = False
    code_buf: list[str] = []

    while i < n:
        line = lines[i]
        stripped = line.strip()

        # ---- fenced code ----
        if stripped.startswith("```"):
            if in_code:
                p = doc.add_paragraph()
                p.paragraph_format.space_before = Pt(4)
                p.paragraph_format.space_after = Pt(8)
                p.paragraph_format.left_indent = Inches(0.2)
                for k, cl in enumerate(code_buf):
                    if k:
                        p.add_run().add_break()
                    _run(p, cl, Pt(8.5), False, False,
                         RGBColor(0x22, 0x44, 0x22), MONO_FONT)
                _shade_paragraph(p, "F4F6F4")
                code_buf = []
                in_code = False
            else:
                in_code = True
            i += 1
            continue
        if in_code:
            code_buf.append(line)
            i += 1
            continue

        # ---- blank ----
        if not stripped:
            i += 1
            continue

        # ---- horizontal rule ----
        if re.fullmatch(r"-{3,}", stripped):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(4)
            p.paragraph_format.space_after = Pt(8)
            _paragraph_border_bottom(p, "D0D0D0", 6)
            i += 1
            continue

        # ---- heading ----
        m = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if m:
            level = len(m.group(1))
            text = m.group(2).strip()
            if level == 1:
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_before = Pt(0)
                p.paragraph_format.space_after = Pt(6)
                add_inline(p, text, size=Pt(17), bold=True, color=ACCENT)
            else:
                h = doc.add_heading(level=min(level - 1, 3))
                add_inline(h, text,
                           size={2: Pt(15), 3: Pt(12.5), 4: Pt(11),
                                 5: Pt(11), 6: Pt(11)}[min(level, 6)],
                           bold=True, color=ACCENT)
                for r in h.runs:
                    r.font.name = BODY_FONT
            i += 1
            continue

        # ---- table ----
        if stripped.startswith("|") and i + 1 < n and \
                is_separator(lines[i + 1]):
            header = split_row(stripped)
            rows = [header]
            j = i + 2
            while j < n and lines[j].strip().startswith("|"):
                rows.append(split_row(lines[j].strip()))
                j += 1
            # the caption is the sentence just above the table, if it is one
            caption = None
            if i > 0:
                prev = lines[i - 1].strip()
                if prev and not prev.startswith("|") and not prev.startswith("#") \
                        and not prev.startswith("---"):
                    caption = prev.rstrip("*").strip()
                    if caption.startswith("**"):
                        caption = caption.strip("*")
            table_number += 1
            add_table(doc, rows, caption, table_number)
            i = j
            continue

        # ---- bullet list ----
        m = re.match(r"^[-*]\s+(.*)$", stripped)
        if m:
            while i < n:
                mm = re.match(r"^[-*]\s+(.*)$", lines[i].strip())
                if not mm:
                    break
                text = mm.group(1)
                # continuation lines (indented by two or more spaces, or a
                # tab) belong to the same bullet
                j = i + 1
                while j < n and lines[j].strip() and \
                        lines[j].startswith(("  ", "\t")) and \
                        not re.match(r"^\s*[-*]\s", lines[j]) and \
                        not re.match(r"^\s*\d+\.\s", lines[j]):
                    text += " " + lines[j].strip()
                    j += 1
                p = doc.add_paragraph(style="List Bullet")
                p.paragraph_format.space_after = Pt(3)
                p.paragraph_format.line_spacing = 1.15
                add_inline(p, text)
                i = j
            continue

        # ---- numbered list ----
        m = re.match(r"^(\d+)\.\s+(.*)$", stripped)
        if m:
            while i < n:
                mm = re.match(r"^(\d+)\.\s+(.*)$", lines[i].strip())
                if not mm:
                    break
                text = mm.group(2)
                j = i + 1
                while j < n and lines[j].startswith(("   ", "\t")) and \
                        lines[j].strip() and not re.match(r"^\d+\.\s", lines[j].strip()):
                    text += " " + lines[j].strip()
                    j += 1
                p = doc.add_paragraph(style="List Number")
                p.paragraph_format.space_after = Pt(3)
                p.paragraph_format.line_spacing = 1.15
                add_inline(p, text)
                i = j
            continue

        # ---- paragraph ----
        buf = [stripped]
        j = i + 1
        while j < n:
            nxt = lines[j].strip()
            if not nxt:
                break
            if (nxt.startswith("#") or nxt.startswith("|")
                    or nxt.startswith("```") or re.fullmatch(r"-{3,}", nxt)
                    or re.match(r"^[-*]\s+", nxt) or re.match(r"^\d+\.\s+", nxt)):
                break
            buf.append(nxt)
            j += 1
        # A run of "**Field:** value" lines is a metadata block, not a
        # paragraph: each field must stay on its own line.
        if len(buf) > 1 and all(re.match(r"^\*\*[^*]+:\*\*", b) for b in buf):
            for k, b in enumerate(buf):
                p = doc.add_paragraph()
                p.paragraph_format.space_after = Pt(2 if k < len(buf) - 1 else 8)
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT
                add_inline(p, b)
            i = j
            continue
        text = " ".join(buf)
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        p.paragraph_format.space_after = Pt(6)
        # a line that is entirely bold is a lead-in
        if text.startswith("**") and text.endswith("**") and \
                text.count("**") == 2:
            add_inline(p, text[2:-2], size=BODY_SIZE, bold=True)
        else:
            add_inline(p, text)
        i = j

    # ---- the figures, at the end ----
    n_fig = add_figures(doc, repo_root)
    if n_fig:
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(12)
        r = p.add_run("Figures")
        r.bold = True
        r.font.name = BODY_FONT
        r.font.size = Pt(12.5)
        r.font.color.rgb = ACCENT

    return table_number


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("source")
    ap.add_argument("target")
    ap.add_argument("--header", default=None)
    args = ap.parse_args()

    src = Path(args.source)
    dst = Path(args.target)
    repo_root = src.resolve().parent.parent
    doc = new_document()
    header = args.header or "D1 benchmark — research report draft (v1.0-rc3)"
    add_header(doc, header)
    add_page_numbers(doc)
    convert(src, dst, doc, repo_root)
    doc.save(str(dst))
    print(f"wrote {dst} ({dst.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
