"""Write the speaker script INTO the deck, page by page, as the
PowerPoint notes.

Until now the deck carried no notes at all: the notes pane was empty,
and anything typed by hand into the .pptx was lost the next time the
deck was rebuilt from its JSON sources. This generator makes the notes a
build artefact instead: it reads the Chinese 20-minute script, assigns
every paragraph to the slide it belongs to (using the slide cues and the
sub-headings), and writes `docs/deck/notes_zh.json`; `scripts/build_deck.js`
then puts those paragraphs into each slide's notes.

Output:
  docs/deck/notes_zh.json   {"1": ["...", "..."], "2": [...], ...}
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")
ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "docs/report_1118_speaker_script_20min.md"
OUT = ROOT / "docs/deck/notes_zh.json"

SECTION_RE = re.compile(r"^## (\d+)\.\s*(.+)$")
CUE_RE = re.compile(r"【图：([^】]+)】")
SLIDE_NUM_RE = re.compile(r"slide\s*(\d+)")
PAPER_RE = re.compile(r"^\*\*第[一二三四五六七八九十]+篇")
SUBHEAD_RE = re.compile(r"^###")


def parse_allocation(lines: list[str]) -> dict[int, tuple[str, str]]:
    """From the section 0.1 table: slide -> (section, minutes)."""
    out: dict[int, tuple[str, str]] = {}
    for line in lines:
        m = re.match(r"\|\s*\*\*(\d+)\*\*\s*\|[^|]*\|\s*§(\d+(?:\.\d+)?)\s*\|"
                     r"\s*([\d.]+)\s*\|", line)
        if m:
            out[int(m.group(1))] = (m.group(2), m.group(3))
    return out


def slides_of(header: str) -> list[int]:
    """The slides a section covers, from its header.

    Handles 'slide 1', 'slides 8、9', 'slides 3–5' and 'slides 12–13'.
    """
    out: list[int] = []
    for m in re.finditer(r"slides?\s*([0-9]+)\s*([–—-])\s*([0-9]+)", header):
        a, b = int(m.group(1)), int(m.group(3))
        out.extend(range(min(a, b), max(a, b) + 1))
    for m in re.finditer(r"slides?\s*([0-9]+(?:\s*[、,]\s*[0-9]+)*)", header):
        if "–" in m.group(0) or "-" in m.group(0) or "—" in m.group(0):
            continue
        out.extend(int(x) for x in re.findall(r"[0-9]+", m.group(1)))
    return out


def main() -> int:
    lines = SCRIPT.read_text(encoding="utf-8").splitlines()
    alloc = parse_allocation(lines)

    notes: dict[str, list[str]] = {}
    section = None
    title = ""
    body: list[str] = []

    def flush() -> None:
        if section is None:
            return
        # the slides of this section: the section 0.1 allocation table is
        # authoritative; the section header is the fallback (for the
        # opening and the closing, which share slides)
        covers = [s for s, (sec, _) in sorted(alloc.items())
                  if sec == str(section) or sec.startswith(str(section) + ".")]
        if not covers:
            covers = slides_of(title)
        if not covers:
            return
        # the blocks of this section
        blocks: list[list[str]] = []
        cur: list[str] = []
        for ln in body:
            if not ln.strip():
                if cur:
                    blocks.append(cur)
                    cur = []
                continue
            cur.append(ln)
        if cur:
            blocks.append(cur)

        groups: list[tuple[int, str]] = []
        for block in blocks:
            cue_slide = None
            kept = []
            for ln in block:
                c = CUE_RE.search(ln)
                if c:
                    m = SLIDE_NUM_RE.search(c.group(1))
                    if m:
                        cue_slide = int(m.group(1))
                    continue
                kept.append(ln)
            text = "\n".join(kept).strip()
            if not text:
                continue
            # stage directions are not spoken text: the notes carry the
            # words, not the reminders
            if text.startswith("【提示") or text.startswith("---"):
                continue
            if cue_slide is not None and cue_slide in covers:
                slide = cue_slide
            elif (groups and (SUBHEAD_RE.match(text) or PAPER_RE.match(text))
                  and covers and groups[-1][0] != covers[-1]):
                # a new paper or a new sub-section starts the next slide
                idx = covers.index(groups[-1][0]) if groups[-1][0] in covers else -1
                slide = covers[idx + 1] if 0 <= idx < len(covers) - 1 else groups[-1][0]
            elif groups:
                slide = groups[-1][0]
            else:
                slide = covers[0]
            groups.append((slide, text))

        for slide, text in groups:
            key = str(slide)
            sec, mins = alloc.get(slide, (str(section), ""))
            head = (f"【§{sec} · {mins} 分钟】" if mins else f"【§{sec}】")
            notes.setdefault(key, []).append(head)
            notes[key].append(text)

    for line in lines:
        m = SECTION_RE.match(line)
        if m:
            flush()
            section = int(m.group(1))
            title = m.group(2)
            body = []
            continue
        if section is not None and not line.startswith("## "):
            body.append(line)
    flush()

    # the evidence pages: a one-line note, so that no slide has an empty
    # notes pane (the presenter flips to them only when asked)
    POINTERS = {
        16: "A1 · 文献与我们的设计决定一一对应",
        17: "A2 · 每个输入的版本、来源与标识符",
        18: "A3 · 24 个臂的完整结果表",
        19: "A4 · H1 在五种区间构造下的判定",
        20: "A5 · 低 FPR 与 top-10 % 两个共主端点",
        21: "A6 · bait 匹配零模型：可见度吃掉多少",
        22: "A7 · 阴性集审计与清洗后重跑",
        23: "A8 · 等阳性率面板：策展滞后的证据",
        24: "A9 · 标签质量分层",
        25: "A10 · 方差分解：24 臂与 30 臂",
        26: "A11 · 四个设计控制",
        27: "A12 · 度是不是稻草基线（Control 19）",
        28: "A13 · 标签主效应是难度轴（Control 20）",
        29: "A14 · 可见度分解在第二个方法下也成立（Control 21）",
        30: "A15 · 我们如何使用 AI，以及从中提炼的思路（第 10 节主动翻出）",
        31: "A16 · AI 带来的两个设计（备用）",
    }
    for n, pointer in POINTERS.items():
        notes.setdefault(str(n), []).append(
            f"【备用证据页】{pointer}。被问到时翻到这里；"
            f"讲解口径见 20 分钟讲稿的对应小节。")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(notes, ensure_ascii=False, indent=1) + "\n",
                    encoding="utf-8")
    print(f"wrote {OUT.name}: {len(notes)} slides carry notes")
    for k in sorted(notes, key=lambda x: int(x)):
        chars = sum(len(x) for x in notes[k])
        print(f"  slide {k:>2}: {len(notes[k])} blocks, {chars} characters")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
