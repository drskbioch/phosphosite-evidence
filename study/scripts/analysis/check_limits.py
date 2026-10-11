"""Measure the main figures and tables against the limits of the target journal (analysis/submission_limits.yml).

For each main figure (Figures 1 to 5 in manuscript/figures/): width and height of the vector file (PDF page), pixel size and dpi of the PNG and TIFF files, and the depth
of the notes that are printed inside the figure file (the lines below the panels, set in the note font). For each main table (Tables 1 to 3 in manuscript/tables/): number of
rows and columns, and the lines of its notes file. A measured value is compared with a limit only if the limit is not TODO; otherwise the line says "limit not set".
Usage: python check_limits.py        Exit status 1 if a measured value is over a limit that is set.
"""
import csv
import re
import sys
from pathlib import Path

import fitz
from PIL import Image

CUR = Path(__file__).resolve().parents[2]
LIMITS = CUR / "analysis" / "submission_limits.yml"
FIG = CUR / "manuscript" / "figures"
TAB = CUR / "manuscript" / "tables"
MAIN_FIGS = [f"Figure{i}" for i in range(1, 6)]
MAIN_TABLES = ["Table1_data_sources", "Table2_main_results", "Table3_predictors"]
MM_PER_PT = 25.4 / 72
UNIT_MM = {"mm": 1.0, "cm": 10.0, "in": 25.4, "pt": MM_PER_PT}


def read_limits():
    out = {}
    for line in LIMITS.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if ":" in line:
            k, v = line.split(":", 1)
            out[k.strip()] = v.strip().strip('"')
    return out


def number(v):
    """(value, unit) of a limit; None when it is TODO."""
    if v in ("", "TODO"):
        return None
    m = re.fullmatch(r"([0-9.]+)\s*(mm|cm|in|pt)?", v)
    return (float(m.group(1)), m.group(2) or "mm") if m else None


def in_mm(lim):
    return lim[0] * UNIT_MM[lim[1]]


def note_block(page):
    """Top and number of lines of the notes at the bottom of a figure page: the lines of the font size of the lowest line, counted upwards while the spacing stays within 1.6 line heights."""
    lines = []
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            sp = l["spans"][0]
            if "".join(s["text"] for s in l["spans"]).strip():
                lines.append((l["bbox"][1], l["bbox"][3], sp["size"]))
    lines.sort()
    if not lines:
        return None, 0
    size = lines[-1][2]
    block = [lines[-1]]
    for ln in reversed(lines[:-1]):
        if abs(ln[2] - size) < 0.06 and block[0][0] - ln[0] <= 1.6 * size * 1.25 and ln[0] < block[0][0] - 0.5:
            block.insert(0, ln)
        elif abs(ln[2] - size) < 0.06 and abs(ln[0] - block[0][0]) < 0.5:
            continue
        else:
            break
    return block[0][0], len(block)


lim = read_limits()
over = 0
print(f"limits file: {LIMITS.relative_to(CUR.parent)}")
for k in ("main_text_word_limit", "abstract_word_limit", "max_keywords", "reference_style", "preprint_policy"):
    print(f"{k}: {'limit not set' if lim.get(k, 'TODO') == 'TODO' else lim[k] + ' (not measured here: the manuscript text is not final)'}")
n_fig = number(lim.get("max_main_figures", "TODO"))
n_tab = number(lim.get("max_main_tables", "TODO"))
print(f"main figures: {len(MAIN_FIGS)}; max_main_figures: " + ("limit not set" if n_fig is None else f"{n_fig[0]:g} -> {'OK' if len(MAIN_FIGS) <= n_fig[0] else 'OVER'}"))
over += n_fig is not None and len(MAIN_FIGS) > n_fig[0]
print(f"main tables: {len(MAIN_TABLES)}; max_main_tables: " + ("limit not set" if n_tab is None else f"{n_tab[0]:g} -> {'OK' if len(MAIN_TABLES) <= n_tab[0] else 'OVER'}"))
over += n_tab is not None and len(MAIN_TABLES) > n_tab[0]
wmax, dmax, dpimin = (number(lim.get(k, "TODO")) for k in ("figure_width_max", "figure_depth_max", "figure_min_dpi"))
print()
for name in MAIN_FIGS:
    d = fitz.open(FIG / f"{name}.pdf")
    p = d[0]
    w_pt, h_pt = p.rect.width, p.rect.height
    top, n_lines = note_block(p)
    note_pt = h_pt - top if top is not None else 0
    msgs = [f"{name}: PDF {w_pt * MM_PER_PT:.1f} x {h_pt * MM_PER_PT:.1f} mm ({w_pt:.1f} x {h_pt:.1f} pt)"]
    msgs.append(f"notes inside the file: {n_lines} lines, {note_pt * MM_PER_PT:.1f} mm ({note_pt:.1f} pt) deep")
    for ext in ("png", "tif"):
        im = Image.open(FIG / f"{name}.{ext}")
        dpi = tuple(float(x) for x in im.info["dpi"]) if "dpi" in im.info else (None, None)
        msgs.append(f"{ext.upper()} {im.size[0]} x {im.size[1]} px, {dpi[0]:.0f} dpi, mode {im.mode}" if dpi[0] else f"{ext.upper()} {im.size[0]} x {im.size[1]} px, dpi not stored, mode {im.mode}")
        if dpimin is not None and dpi[0] and dpi[0] < dpimin[0]:
            msgs.append(f"** {ext.upper()} dpi {dpi[0]:.0f} below figure_min_dpi {dpimin[0]:g} **")
            over += 1
    if wmax is None:
        msgs.append("figure_width_max: limit not set")
    else:
        ok = w_pt * MM_PER_PT <= in_mm(wmax) + 0.05
        msgs.append(f"width limit {in_mm(wmax):.1f} mm: {'OK' if ok else '** OVER **'}")
        over += not ok
    if dmax is None:
        msgs.append("figure_depth_max: limit not set")
    else:
        ok = h_pt * MM_PER_PT <= in_mm(dmax) + 0.05
        msgs.append(f"depth limit {in_mm(dmax):.1f} mm: {'OK' if ok else '** OVER **'}")
        over += not ok
    if dpimin is None:
        msgs.append("figure_min_dpi: limit not set")
    print("; ".join(msgs))
print()
for stem in MAIN_TABLES:
    rows = list(csv.reader(open(TAB / f"{stem}.csv", encoding="utf-8")))
    notes = TAB / f"{stem}_notes.md"
    n_notes = len([l for l in notes.read_text(encoding="utf-8").splitlines() if l.strip()]) if notes.exists() else 0
    print(f"{stem}: {len(rows) - 1} rows x {len(rows[0])} columns; {n_notes} lines of notes; tables are text (no width, height or dpi); max_main_tables: " + ("limit not set" if n_tab is None else "see above"))
sys.exit(1 if over else 0)
