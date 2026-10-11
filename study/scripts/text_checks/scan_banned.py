"""Scan files for banned words (whole-word, case-insensitive).

Usage: python scan_banned.py file [file ...]      (.md, .csv, .txt, .svg, .xlsx, .pdf)
Exit status 1 if any hit is found.

Two exemptions, and no others:
  (a) Table S1 (file name contains "TableS1"): contrast codes of the form G\\d+-\\d+ (for example G7-0123) and the block labels G1 to G23 are
      identifiers, not wording. In other files the pattern G\\d+-\\d+ alone is also exempt (the codes are cited by name in the Supporting Information);
      the bare block labels are exempt only in Table S1 files.
  (b) Table S3 (file name contains "TableS3"): the text of the Title column. Titles are quoted from PubMed and cannot be reworded.
"""
import csv
import io
import re
import sys
from pathlib import Path

WORDS = [w.strip() for w in (Path(__file__).parent / "banned_words.txt").read_text().splitlines() if w.strip()]
pats = [(w, re.compile(r"(?<![A-Za-z0-9])" + re.escape(w).replace(r"\ ", r"\s+") + r"(?![A-Za-z0-9])", re.I)) for w in WORDS]
CODE = re.compile(r"(?<![A-Za-z0-9])G\d+-\d+(?![A-Za-z0-9])")
BLOCK = re.compile(r"(?<![A-Za-z0-9])G\d+(?![A-Za-z0-9])")


def lines_of(path):
    """Text lines of a file; spreadsheets and PDF pages are read as text, one line per row or per text line."""
    suf = path.suffix.lower()
    if suf == ".xlsx":
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        out = []
        for ws in wb.worksheets:
            out.append(ws.title)
            for row in ws.iter_rows(values_only=True):
                out.append(" | ".join("" if c is None else str(c) for c in row))
        return out
    if suf == ".pdf":
        import fitz
        return [ln for page in fitz.open(path) for ln in page.get_text().splitlines()]
    return path.read_text(encoding="utf-8", errors="replace").splitlines()


def mask_title(path, lines):
    """Rule (b): blank the Title column of Table S3 (second cell of a Markdown row, second field of the CSV)."""
    out = []
    if path.suffix.lower() == ".csv":
        for row in csv.reader(io.StringIO("\n".join(lines))):
            if len(row) > 1:
                row[1] = ""
            out.append(" ".join(row))
        return out
    for ln in lines:
        if ln.startswith("|"):
            cells = ln.split("|")
            if len(cells) > 3:
                cells[2] = " "
            ln = "|".join(cells)
        out.append(ln)
    return out


def scan(path):
    path = Path(path)
    lines = lines_of(path)
    if "TableS3" in path.name:
        lines = mask_title(path, lines)
    s1 = "TableS1" in path.name
    hits = []
    for i, line in enumerate(lines, 1):
        text = CODE.sub(" ", line)
        if s1:
            text = BLOCK.sub(" ", text)
        for w, p in pats:
            if p.search(text):
                hits.append((i, w, line.strip()[:120]))
    return hits


def main():
    total = 0
    for name in sys.argv[1:]:
        for i, w, ln in scan(name):
            total += 1
            print(f"{name}:{i}: [{w}] {ln}".encode("ascii", "replace").decode())
    print(f"{total} hit(s)")
    sys.exit(1 if total else 0)


if __name__ == "__main__":
    main()
