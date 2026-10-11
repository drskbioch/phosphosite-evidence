"""Scan the files listed in RELEASE_SCOPE.md (section 6) for EPSD site-level records and predictor weights.

Usage: python check_release_scope.py [output.json] [--list FILE]
Reads the repository-relative paths from the file list in RELEASE_SCOPE.md and reports, for the files that exist:
 1. file types that store tables of records or model weights (parquet, feather, h5, npz, npy, pkl, pt, pth, onnx, ckpt, safetensors, bin, sqlite, db, gz, zip, tsv, joblib); workbooks (xlsx) are searched by content;
 2. in text (json, csv, md, py, txt, log, yml, cff, svg, toml and text extracted from PDF files when PyMuPDF is installed):
    a UniProt accession next to a residue letter S, T or Y and a position, in either order; EPSD identifiers (EPSD followed by digits);
    amino-acid runs of 20 or more letters; CSV headers with a UniProt column and a position column; accession-like tokens elsewhere;
 3. JSON lists of 500 or more numbers (weights or per-site vectors);
 4. the ten largest files, and any listed path that starts with a prefix of the EXCLUDED list (legacy/, frontend/, legacy/deploy/, study/analysis/internal/, study/_archive/, .claude/, nested repositories,
    study/tools/, study/processed_data/, study/analysis/data_private/) or lies inside a nested repository (a folder with a .git entry).
PNG and TIFF files cannot be searched as text and are counted. Exit status 1 if tests 1 to 3 or the path test find anything other than the known
false positives printed under "Other hits" (the amino-acid alphabet in run_j.py and hex colour codes in fig1_2.py are reported, not hidden).
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LIST = Path(sys.argv[sys.argv.index("--list") + 1]) if "--list" in sys.argv else ROOT / "RELEASE_SCOPE.md"  # --list FILE: check another list (for tests)

BAD_EXT = {".parquet", ".feather", ".h5", ".hdf5", ".npz", ".npy", ".pkl", ".pickle", ".pt", ".pth", ".onnx", ".ckpt", ".safetensors", ".bin", ".sqlite", ".db", ".gz", ".zip", ".tsv", ".joblib"}
ACC = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?"
SEP = r"[ \t_:\-|,;]{1,3}"
SITE = [re.compile(rf"(?<![A-Za-z0-9]){ACC}{SEP}[STY]{SEP}?[0-9]{{1,5}}(?![A-Za-z0-9])"),
        re.compile(rf"(?<![A-Za-z0-9]){ACC}{SEP}[0-9]{{1,5}}{SEP}[STY](?![A-Za-z0-9])"),
        re.compile(rf"(?<![A-Za-z0-9])[STY][0-9]{{1,5}}{SEP}{ACC}(?![A-Za-z0-9])")]
P_EPSD = re.compile(r"EPSD[0-9]{3,}")
P_SEQ = re.compile(r"(?<![A-Za-z])[ACDEFGHIKLMNPQRSTVWY]{20,}(?![A-Za-z])")
P_ACC = re.compile(rf"(?<![A-Za-z0-9]){ACC}(?![A-Za-z0-9])")
TEXT_EXT = {".json", ".csv", ".md", ".py", ".txt", ".log", ".yml", ".cff", ".svg", ".toml"}
EXCLUDED = (  # (path prefix, reason): a listed path that starts with one of these fails the check; RELEASE_SCOPE.md section 3 is written from this list
    ("legacy/", "pre-pivot product: web application, deployment folders and server notes, earlier scripts"),
    ("legacy/deploy/", "deployment folders and server notes"),
    ("legacy/deploy/deepphos-frontend/", "nested repository"),
    ("frontend/", "web application; nested repository"),
    ("study/analysis/internal/", "internal table of the changes between the earlier and the present estimates"),
    ("study/_archive/", "archive of moved files"),
    (".claude/", "local Claude Code settings"),
    ("study/tools/", "predictor code and weights"),
    ("study/processed_data/", "EPSD-derived tables"),
    ("study/analysis/data_private/", "caches and downloads"),
)
PRIVATE = tuple(p for p, _ in EXCLUDED)


def nested_repo(rel):
    """True if a folder above the file (below the repository root) holds a .git entry: the file belongs to a nested repository."""
    d = (ROOT / rel).parent
    while d != ROOT and ROOT in d.parents:
        if (d / ".git").exists():
            return True
        d = d.parent
    return False


def listed_paths():
    paths = []
    in_list = False
    for line in LIST.read_text(encoding="utf-8").splitlines():
        if line.startswith("## 6."):
            in_list = True
            continue
        if in_list:
            m = re.match(r"\| [^|]+ \| `([^`]+)` \| [0-9,]+ \|", line)
            if m:
                paths.append(m.group(1))
    return paths


def numeric_lists(o, path=""):
    out = []
    if isinstance(o, list):
        if len(o) >= 500 and all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in o):
            out.append((path, len(o)))
        for i, x in enumerate(o[:50000]):
            if isinstance(x, (list, dict)):
                out += numeric_lists(x, f"{path}[{i}]")
    elif isinstance(o, dict):
        for k, v in o.items():
            out += numeric_lists(v, f"{path}.{k}" if path else k)
    return out


def xlsx_text(p):
    """Cell values of a workbook as text, one line per row (workbooks are searched by content, not rejected by their extension); also returns the header rows."""
    import openpyxl
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    lines, heads = [], []
    for ws in wb.worksheets:
        for i, row in enumerate(ws.iter_rows(values_only=True)):
            s = ",".join("" if v is None else str(v) for v in row)
            if i == 0:
                heads.append(s.lower())
            lines.append(s)
    return chr(10).join(lines), heads


def pdf_text(p):
    import fitz  # PyMuPDF
    with fitz.open(p) as d:
        return "\n".join(pg.get_text() for pg in d)


def main():
    rels = listed_paths()
    files = [(r, ROOT / r) for r in rels]
    missing = [r for r, p in files if not p.exists()]
    files = [(r, p) for r, p in files if p.exists()]
    res = {"listed": len(rels), "missing": missing, "bytes": sum(p.stat().st_size for _, p in files), "bad_extension": [], "private_path": [r for r in rels if r.startswith(PRIVATE) or nested_repo(r)],
           "site_pattern": [], "epsd_id": [], "csv_header": [], "numeric_lists": [], "sequence_run": [], "accession_like": [], "text_scanned": 0, "pdf_scanned": 0, "pdf_skipped": 0, "images_not_scanned": 0}
    for rel, p in files:
        ext = p.suffix.lower()
        if ext in BAD_EXT:
            res["bad_extension"].append(rel)
        text = None
        if ext in TEXT_EXT or p.name in ("LICENSE", ".zenodo.json"):
            text = p.read_text(encoding="utf-8", errors="replace")
            res["text_scanned"] += 1
        elif ext == ".xlsx":
            try:
                text, heads = xlsx_text(p)
                res["xlsx_scanned"] = res.get("xlsx_scanned", 0) + 1
                if any("uniprot" in h and "position" in h for h in heads):
                    res["csv_header"].append(rel)
            except ImportError:
                res["xlsx_skipped"] = res.get("xlsx_skipped", 0) + 1
        elif ext == ".pdf":
            try:
                text = pdf_text(p)
                res["pdf_scanned"] += 1
            except ImportError:
                res["pdf_skipped"] += 1
        elif ext in (".png", ".tif", ".tiff"):
            res["images_not_scanned"] += 1
        if text is None:
            continue
        for pat in SITE:
            ms = list(pat.finditer(text))
            if ms:
                res["site_pattern"].append({"file": rel, "count": len(ms), "examples": [m.group()[:50] for m in ms[:3]]})
        if P_EPSD.search(text):
            res["epsd_id"].append(rel)
        seq = P_SEQ.findall(text)
        if seq:
            res["sequence_run"].append({"file": rel, "examples": seq[:2]})
        acc = P_ACC.findall(text)
        if acc:
            res["accession_like"].append({"file": rel, "tokens": acc[:3]})
        if ext == ".csv":
            head = text.splitlines()[0].lower() if text else ""
            if "uniprot" in head and "position" in head:
                res["csv_header"].append(rel)
        if ext == ".json":
            try:
                res["numeric_lists"] += [{"file": rel, "path": a, "length": n} for a, n in numeric_lists(json.loads(text))]
            except ValueError:
                pass
    res["largest"] = [{"file": r, "bytes": p.stat().st_size} for r, p in sorted(files, key=lambda t: -t[1].stat().st_size)[:10]]
    failed = bool(missing or res["bad_extension"] or res["private_path"] or res["site_pattern"] or res["epsd_id"] or res["csv_header"] or res["numeric_lists"])
    print(f"files listed {res['listed']}, found {len(files)}, missing {len(missing)}, {res['bytes'] / 1e6:.2f} MB")
    print(f"workbooks scanned by content {res.get('xlsx_scanned', 0)} (skipped {res.get('xlsx_skipped', 0)})")
    print(f"text/SVG scanned {res['text_scanned']}, PDF scanned {res['pdf_scanned']} (skipped {res['pdf_skipped']}), PNG/TIFF not searchable {res['images_not_scanned']}")
    for key in ("bad_extension", "private_path", "site_pattern", "epsd_id", "csv_header", "numeric_lists"):
        print(f"{key}: {len(res[key])}")
    print("Other hits (reported, not counted as failures):")
    for key in ("sequence_run", "accession_like"):
        for x in res[key]:
            print(f"  {key}: {x}")
    if res["private_path"]:
        print("Excluded or nested-repository paths in the list:")
        for r in res["private_path"]:
            print("  " + r)
    if len(sys.argv) > 1 and not sys.argv[1].startswith("--"):
        Path(sys.argv[1]).write_text(json.dumps(res, indent=1), encoding="utf-8")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
