"""Render Tables 1 to 3 and S1 to S5 to study/analysis/tables.pdf for review.
Reads the CSV files in study/manuscript/tables/ (the same files the tables are built from). Table S1 shows only its first 40 rows.
Nothing is recomputed: each cell is the text of the CSV cell, wrapped to its column."""
import csv
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

CUR = Path(__file__).resolve().parents[2]
TAB = CUR / "manuscript" / "tables"
OUT = CUR / "analysis" / "tables.pdf"
S1_ROWS = 40

TABLES = [("Table1_data_sources", "Table 1. Data sources", "A4"), ("Table2_main_results", "Table 2. Main results", "A4"),
          ("Table3_predictors", "Table 3. Predictors", "A4"), ("TableS1_grid_status", "Table S1. Status of every planned contrast", "A3"),
          ("TableS2_main_results_P-human", "Table S2. Main results, P-human", "A4"), ("TableS3_earlier_to_present_changes", "Table S3. Changes between the two analyses", "A4"),
          ("TableS4_non_primary_pmids", "Table S4. Non-primary PMIDs", "A4"), ("TableS5_uniprot_experimental_vs_inferred", "Table S5. UniProt experimental against inferred evidence", "A4")]
PAGE = {"A4": (11.69, 8.27), "A3": (16.54, 11.69)}
FONT = {"A4": 6.0, "A3": 5.0}
CW = 0.56  # average character width in em


def read(name):
    with open(TAB / f"{name}.csv", encoding="utf-8-sig", newline="") as fh:
        return list(csv.reader(fh))


def layout(rows, width_pt, fs):
    ncol = len(rows[0])
    want = []
    for j in range(ncol):
        lens = sorted(len(r[j]) for r in rows[1:] if j < len(r)) or [1]
        typical = lens[int(0.9 * (len(lens) - 1))]
        want.append(max(min(max(typical, len(rows[0][j]) * 0.6), 48), 5, max(len(w) for w in rows[0][j].split() or [""]) + 1))
    scale = width_pt / (sum(want) * CW * fs + ncol * 4)
    return [w * CW * fs * scale + 4 for w in want]


def wrap(cell, wpt, fs):
    n = max(int((wpt - 4) / (CW * fs)), 3)
    return textwrap.wrap(cell, n, break_long_words=True) or [""]


def render(pdf, rows, title, size, note=None, foot=None):
    pw, ph = PAGE[size]
    fs = FONT[size]
    margin = 0.4 * 72
    width = pw * 72 - 2 * margin
    cols = layout(rows, width, fs)
    scale = width / sum(cols)
    cols = [c * scale for c in cols]
    lh = fs * 1.25
    wrapped = [[wrap(c, cols[j], fs) for j, c in enumerate(r)] for r in rows]
    heights = [max(len(c) for c in r) * lh + 3 for r in wrapped]
    top = ph * 72 - margin - 16
    bottom = margin + (14 if foot or note else 0)
    page_no = 0
    i = 1
    while i < len(rows):
        page_no += 1
        fig = plt.figure(figsize=(pw, ph))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, pw * 72)
        ax.set_ylim(0, ph * 72)
        ax.axis("off")
        ax.text(margin, ph * 72 - margin, title + ("" if page_no == 1 else " (continued)"), fontsize=fs + 3, fontweight="bold", va="top")
        y = top

        def draw(cells, h, bold=False, shade=None):
            nonlocal y
            if shade:
                ax.add_patch(plt.Rectangle((margin, y - h), width, h, fc=shade, ec="none"))
            x = margin
            for j, lines in enumerate(cells):
                ax.text(x + 2, y - 2, "\n".join(lines), fontsize=fs, va="top", ha="left", fontweight="bold" if bold else "normal", linespacing=1.25)
                x += cols[j]
            ax.plot([margin, margin + width], [y - h, y - h], color="#888888", lw=0.3)
            y -= h

        draw(wrapped[0], heights[0], bold=True, shade="#e6e6e6")
        while i < len(rows) and y - heights[i] > bottom:
            draw(wrapped[i], heights[i], shade="#f7f7f7" if i % 2 == 0 else None)
            i += 1
        if i >= len(rows):
            text = "  ".join(t for t in (note, foot) if t)
            if text:
                ax.text(margin, margin, "\n".join(textwrap.wrap(text, int(width / (CW * fs)))), fontsize=fs, va="bottom")
        pdf.savefig(fig)
        plt.close(fig)


def main():
    with PdfPages(OUT) as pdf:
        for name, title, size in TABLES:
            rows = read(name)
            note = None
            foot = None
            if name == "Table3_predictors":
                note = (TAB / "Table3_predictors_note.md").read_text(encoding="utf-8").strip()
            if name == "TableS1_grid_status":
                total = len(rows) - 1
                rows = rows[:S1_ROWS + 1]
                foot = f"First {S1_ROWS} of {total} rows shown; the full table is TableS1_grid_status.csv."
            render(pdf, rows, title, size, note=note, foot=foot)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
