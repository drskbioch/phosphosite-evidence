"""Banned-word scan of the Supporting Information package (result: results/si_banned_scan.md).

scan_banned.py exempts the Title column of Table S3 by file name (the name must contain "TableS3"), so a PDF cannot carry that exemption. The package is therefore scanned in
parts: the pages of SupportingInformation.pdf without the Table S3 pages, TableS1_contrast_grid.xlsx, SI_FILE_LIST.md, and Table S3 from its CSV and Markdown files
(the same text as the Table S3 pages of the PDF). The two exemption rules of scan_banned.py apply and no others.
Usage: python scan_si.py        Exit status 1 if any hit remains.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import fitz

CUR = Path(__file__).resolve().parents[2]
SI = CUR / "manuscript" / "supporting_information"
TAB = CUR / "manuscript" / "tables"
SCAN = CUR / "scripts" / "text_checks" / "scan_banned.py"


def scan(files):
    r = subprocess.run([sys.executable, str(SCAN), *map(str, files)], capture_output=True, text=True, encoding="utf-8")
    return r.returncode, (r.stdout + r.stderr).strip()


doc = fitz.open(SI / "SupportingInformation.pdf")
# the Table S3 pages are the pages from the one that starts with "Table S3." up to the page before the next table heading
first = next(i for i, p in enumerate(doc) if p.get_text().lstrip().startswith("Table S3."))
nxt = next(i for i, p in enumerate(doc) if i > first and p.get_text().lstrip().startswith("Table S4."))
s3 = list(range(first, nxt))
parts = []
with tempfile.TemporaryDirectory() as tmp:
    part = Path(tmp) / "SupportingInformation_without_TableS3_pages.pdf"
    d2 = fitz.open(str(SI / "SupportingInformation.pdf"))
    d2.select([i for i in range(len(d2)) if i not in s3])
    d2.save(str(part))
    n_other = len(d2)
    d2.close()
    doc.close()
    results = [("SupportingInformation.pdf, pages other than Table S3 (%d pages)" % n_other, scan([part])),
               ("TableS1_contrast_grid.xlsx", scan([SI / "TableS1_contrast_grid.xlsx"])),
               ("SI_FILE_LIST.md", scan([SI / "SI_FILE_LIST.md"])),
               ("Table S3 (CSV and Markdown; Title column exempt)", scan([TAB / "TableS3_non_primary_pmids.csv", TAB / "TableS3_non_primary_pmids.md", TAB / "TableS3_non_primary_pmids_notes.md"]))]
bad = sum(1 for _, (rc, _) in results if rc)
lines = ["# Banned-word scan of the Supporting Information", "", f"Table S3 pages of the PDF: {[i + 1 for i in s3]} (scanned through the Table S3 files, see `scripts/analysis/scan_si.py`).", "",
         "| Part | Result |", "|---|---|"] + [f"| {n} | {'**' + out.splitlines()[-1] + '**' if rc else out.splitlines()[-1]} |" for n, (rc, out) in results]
(CUR / "analysis" / "results" / "si_banned_scan.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("; ".join(f"{n.split(',')[0].split(' (')[0]}: {out.splitlines()[-1]}" for n, (rc, out) in results))
sys.exit(1 if bad else 0)
