"""Write RELEASE_SCOPE.md: the files of the archive snapshot, their sizes, what is excluded, and the result of the scan for EPSD site-level records and predictor weights.

Sections 1, 2 (untracked files), 3 (excluded paths), 4 (scan result) and 6 (file list) are written from the working tree and from the output of check_release_scope.py; the
prose of section 5 is kept from the existing file. Run from anywhere:  python study/scripts/release/build_release_scope.py
Nothing is uploaded or deposited.
"""
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_release_scope import EXCLUDED as EXCLUDED_LIST  # noqa: E402  (one list for the check and for this document)
DOC = ROOT / "RELEASE_SCOPE.md"
CHECK = ROOT / "study" / "scripts" / "release" / "check_release_scope.py"
NL = "\n"


def files(*patterns, recursive=False, exclude=()):
    out = set()
    for pat in patterns:
        for p in (ROOT.glob(pat) if not recursive else ROOT.rglob(pat)):
            if p.is_file() and "__pycache__" not in p.parts and p.name not in exclude:
                out.add(p.relative_to(ROOT).as_posix())
    return sorted(out)


def old_rows():
    t = DOC.read_text(encoding="utf-8")
    sec = t[t.index("## 6. File list"):]
    return [tuple(c.strip().strip("`") for c in l.strip().strip("|").split("|")) for l in sec.splitlines() if l.startswith("|") and not l.startswith("|---") and not l.startswith("| Group")]


old = old_rows()
keep = lambda grp: [p for g, p, _ in old if g.startswith(grp)]
RESULTS = "study/analysis/results/"
MANIFESTS = ["study/analysis/grid_manifest.csv", RESULTS + "manifest_status.csv", RESULTS + "downloads_manifest.json"]
GROUPS = [  # (group, paths, description of the paths)
    ("Rev 6 scripts", files("study/scripts/analysis/*.py"), "`study/scripts/analysis/*.py` (`__pycache__` is left out), including the table check, the Supporting Information build and scan, and the limits checker"),
    ("Reference copies of the checkers at the freeze", files("study/scripts/analysis/checkers_at_freeze/*"), "`study/scripts/analysis/checkers_at_freeze/`: the three number checkers as they were at commit 3521d5c, for reference; not run"),
    ("Plan, deviation log and wording", ["study/analysis/ANALYSIS_PLAN.md", "study/analysis/table_wording.json", "study/analysis/submission_limits.yml"],
     "`study/analysis/ANALYSIS_PLAN.md` (the deviation log is its section 15), `table_wording.json` (cells whose wording differs from the stored text), `submission_limits.yml` (journal limits, all `TODO`)"),
    ("Manifests", MANIFESTS, "`study/analysis/grid_manifest.csv`, `study/analysis/results/manifest_status.csv`, `study/analysis/results/downloads_manifest.json`"),
    ("Stored results", [p for p in files("study/analysis/results/**/*") if p not in MANIFESTS],
     "everything else in `study/analysis/results/`: result JSON files (including `d2b/`, `prop/`, `f_runs/`), the CSV and Markdown reports of the number checks and of the table and Supporting Information checks, and the run logs (`*.log`)"),
    ("Claims and results files", ["study/analysis/CLAIMS_VS_EVIDENCE.md", "study/analysis/RESULTS_REPORT.md"], "`study/analysis/CLAIMS_VS_EVIDENCE.md`, `study/analysis/RESULTS_REPORT.md`"),
    ("Figure and table code", files("study/scripts/figures/*.py", exclude=("list_old_figures.py",)), "`study/scripts/figures/*.py` except `list_old_figures.py`"),
    ("Composite figures, all formats", files("study/manuscript/figures/Figure*.pdf", "study/manuscript/figures/Figure*.svg", "study/manuscript/figures/Figure*.png", "study/manuscript/figures/Figure*.tif"),
     "`study/manuscript/figures/Figure*.{pdf,svg,png,tif}` (13 figures, 4 formats each)"),
    ("Panel PNGs and panel manifest", files("study/manuscript/figures/panels/Fig*.png", "study/manuscript/figures/panels/PANEL_MANIFEST.md"), "`study/manuscript/figures/panels/Fig*.png` (PNG only) and `PANEL_MANIFEST.md`"),
    ("`phosnet` package", keep("phosnet"), "`study/phosnet/`: source, tests, `LICENSE`, `README.md`, `pyproject.toml`, `CHANGELOG.md` (`validate_engine_modes.py` and `validate_engine_v50.py` import its source)"),
    ("Root documents, requirements, licence and citation files", files("README.md", "PROVENANCE.md", "TITLE.txt", ".gitignore", "scripts/set_title.py", "REPRODUCE.md", "INPUTS.md", "RELEASE_SCOPE.md", "requirements.txt", "LICENSE", "CITATION.cff", ".zenodo.json", "study/scripts/release/*.py"),
     "`README.md`, `PROVENANCE.md`, `TITLE.txt`, `.gitignore`, `scripts/set_title.py`, `REPRODUCE.md`, `INPUTS.md`, `RELEASE_SCOPE.md`, `requirements.txt`, `LICENSE`, `CITATION.cff`, `.zenodo.json`, `study/scripts/release/*.py`"),
    ("Tables written by the table code", files("study/manuscript/tables/Table*.csv", "study/manuscript/tables/Table*.md"), "`study/manuscript/tables/Table*.{csv,md}`: Tables 1 to 3 and S1 to S5, the key to Table S1 and the notes files; the figure and table check reads them"),
    ("Supporting Information package", files("study/manuscript/supporting_information/*"), "`study/manuscript/supporting_information/`: `SupportingInformation.pdf`, `TableS1_contrast_grid.xlsx`, `SI_FILE_LIST.md`"),
    ("Files the number checks read", keep("Files the number checks read"), "`study/analysis/INTERIM_REPORT.md`, `RECONCILIATION.md`, `sources/SOURCES.md`, `nonprimary_pmid_list.md`; `study/audit_reports/phaseH1_step3_predictor_audit.md`; the five result files `study/earlier_results/*.json`"),
    ("Modules a listed script imports", sorted(set(keep("Modules a listed script imports") + files("study/scripts/text_checks/scan_banned.py", "study/scripts/text_checks/banned_words.txt"))), "`study/scripts/earlier_figures/figlib.py` (figure style), `study/scripts/earlier_predictor_audit/step3/score_metrics.py` (score loaders for `run_j.py`), `study/scripts/text_checks/scan_banned.py` and `banned_words.txt` (read by `scan_si.py`)"),
]
EXCLUDED = [p for p, _ in EXCLUDED_LIST]
for g, ps, _ in GROUPS:
    for p in ps:
        assert not any(p.startswith(x) for x in EXCLUDED), p
seen = set()
rows = []
for g, ps, _ in GROUPS:
    for p in ps:
        assert p not in seen, p
        seen.add(p)
        assert (ROOT / p).exists(), p
        rows.append((g, p, (ROOT / p).stat().st_size))
mb = lambda b: f"{b / 1e6:.2f} MB"
tracked = set(subprocess.run(["git", "-C", str(ROOT), "ls-files"], capture_output=True, text=True).stdout.splitlines())
untracked = [p for _, p, _ in rows if p not in tracked]


def sections(scan):
    s1 = ["## 1. Contents and sizes", "", "| Group | Files | Size | Paths |", "|---|---|---|---|"]
    for g, ps, d in GROUPS:
        s1.append(f"| {g} | {len(ps)} | {mb(sum((ROOT / p).stat().st_size for p in ps))} | {d} |")
    s1.append(f"| **Total in the snapshot** | **{len(rows)}** | **{mb(sum(b for _, _, b in rows))}** | |")
    s1 += ["", "The number checks are `check_numbers_sourced.py`, `check_bindings.py`, `check_literals.py` and `check_figure_bindings.py` (in the scripts group) and their reports `number_source_check.md`, `number_binding_check.md`, `builder_numeric_literals.md` and `figure_binding_check.md` (in the stored results group). The table check is `check_table_numbers.py` with its reports `table_number_check.md` and `table_number_binding.csv`. They are part of the groups above, not additional files.", "",
           "The figure code also writes the composite figures to `study/manuscript/earlier_figures/`; that second copy is not in the snapshot.", ""]
    n_logs = len([p for _, p, _ in rows if p.endswith(".log")])
    s2 = ["## 2. Run logs", "", f"- The {n_logs} run logs in `study/analysis/results/` match `.gitignore` (`*.log`) and are tracked by force. `check_literals.py` reads them: without them it lists one more literal as undocumented.",
          "- Redaction: five logs (`f_G8_S1-S2.log`, `f_G8_S3-S4.log`, `run_propensity.log`, `run_propensity2.log`, `run_propensity3.log`) contained Python tracebacks with the local Windows user path. In them the path `C:\\Users\\<name>` is replaced by `<user-home>`. Nothing else in the logs was changed.",
          "- Files in the list that git does not track yet: " + (", ".join(f"`{p}`" for p in untracked) if untracked else "none"), ""]
    s3 = ["## 3. EXCLUDED: not in the snapshot", "",
          "`check_release_scope.py` fails if a path in the file list of section 6 starts with one of these prefixes or lies inside a nested repository (a folder with a `.git` entry). The prefixes are the list `EXCLUDED` in that script.", "",
          "| Excluded path | Reason |", "|---|---|"]
    s3 += [f"| `{p}` | {why} |" for p, why in EXCLUDED_LIST]
    s3 += ["", "Nested repositories (git links or folders with a `.git` entry): `frontend/`, `legacy/deploy/deepphos-frontend/`, and the predictor checkouts inside `study/_archive/`. Their uncommitted changes are not part of any commit of this work.", "",
           "Also not in the snapshot:", "",
           "- `study/manuscript/earlier_figures/final/Figures.pptx` (the author's own file) and everything else under `study/manuscript/` except the figures, panels, tables and Supporting Information package named above.",
           "- `study/scripts/figures/list_old_figures.py` (lists figure files of earlier stages), `study/analysis/PACKAGE_STATUS.md`, `SUBMISSION_STATE.md`, `SUBMISSION_READINESS.md`, `STATEMENT_FACTS.md`, the contact-sheet and review PDFs in `study/analysis/`, `study/analysis/FIGURE_INPUTS.md`, `FIGURE_TABLE_INVENTORY.md`, `OLD_FIGURES_LIST.md`.",
           "- The scripts of earlier stages (`REPRODUCE.md`, section 7) other than the modules named in the last group of section 1, and the web application folders.", ""]
    s4 = ["## 4. Scan for EPSD site-level records and predictor weights", ""]
    if scan is None:
        s4 += ["(scan pending)", ""]
    else:
        r = scan
        s4 += [f"`python study/scripts/release/check_release_scope.py` reads the file list in section 6 of this file and scans every file in it: {r['listed']} files, {r['bytes'] / 1e6:.2f} MB. It searched {r['text_scanned']} text and SVG files, {r.get('xlsx_scanned', 0)} workbook by its cell values, and the text of {r['pdf_scanned']} PDF files. The {r['images_not_scanned']} PNG and TIFF files cannot be searched as text; they are drawn by the same code as the SVG and PDF versions, which were searched.", "",
               "| Test | Result |", "|---|---|",
               f"| File types that hold records or model weights (parquet, feather, h5, npz, npy, pkl, pt, pth, onnx, ckpt, safetensors, bin, sqlite, db, gz, zip, tsv, joblib); workbooks (xlsx) are searched by content | {len(r['bad_extension'])} files |",
               f"| A listed path that starts with an excluded prefix (section 3) or lies inside a nested repository | {len(r['private_path'])} files |",
               f"| A UniProt accession next to a residue letter (S, T or Y) and a position, in either order, with the usual separators | {len(r['site_pattern'])} files |",
               f"| EPSD identifiers (`EPSD` followed by digits) | {len(r['epsd_id'])} files |",
               f"| CSV files and workbook sheets whose header has both a UniProt column and a position column | {len(r['csv_header'])} files |",
               f"| JSON lists of 500 or more numbers (weights or per-site vectors) | {len(r['numeric_lists'])} lists |"]
        seqs = sorted({x["file"] for x in r["sequence_run"]})
        s4.append("| Amino-acid runs of 20 or more letters | " + (f"{len(seqs)} files ({', '.join('`' + f + '`' for f in seqs)}): the string of the 20 amino-acid letters (an encoding table in `run_j.py`, a search pattern in the checker), not a sequence from a protein" if seqs else "0 files") + " |")
        notes = {"study/scripts/figures/fig1_2.py": "hex colour codes", "study/phosnet/tests/test_smoke.py": "not classified",
                 "study/audit_reports/phaseH1_step3_predictor_audit.md": "one accession in an example output line of the DeepPhos test file (accession, position, probability), not an EPSD record"}
        acc = [x["file"] for x in r["accession_like"]]
        s4.append(f"| Accession-like tokens outside the tests above | {len(acc)} files: " + "; ".join(f"`{f}`: {notes.get(f, 'NOT CLASSIFIED')}" for f in acc) + " |")
        s4 += ["", "The largest files are " + ", ".join(f"`{Path(x['file']).name}` {x['bytes'] / 1e6:.2f} MB" for x in r["largest"][:3]) + ".", "",
               "The stored result files hold aggregate statistics: estimates, intervals, counts, and lists of PMIDs and of PubMed titles. Where they name proteins (for example the largest gene clusters) they do so in summary form, never as records of sites. The scan does not read meaning; it shows that no table of site records and no file of weights is present.", ""]
    return s1, s2, s3, s4


def write(scan):
    t = DOC.read_text(encoding="utf-8")
    s5 = t[t.index("## 5. File names"):t.index("## 6. File list")].rstrip()
    s1, s2, s3, s4 = sections(scan)
    head = ["# Release scope", "", "Exactly what the archive contains: the files listed in section 6, with these paths. The archive is deposited from a GitHub release; its version is 1.0.0 (`CITATION.cff`, `.zenodo.json`). Sizes are in bytes or MB (1 MB = 1,000,000 bytes). Nothing has been uploaded or deposited.", ""]
    s6 = ["## 6. File list", "", "| Group | Path | Bytes |", "|---|---|---|"] + [f"| {g} | `{p}` | {b:,} |" for g, p, b in rows]
    DOC.write_text(NL.join(head + s1 + s2 + s3 + s4 + [s5, ""] + s6) + NL, encoding="utf-8")


write(None)
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "scan.json"
    r = subprocess.run([sys.executable, str(CHECK), str(out)], capture_output=True, text=True, encoding="utf-8")
    print(r.stdout)
    scan = json.loads(out.read_text(encoding="utf-8"))
write(scan)
print(f"{len(rows)} files, {sum(b for _, _, b in rows) / 1e6:.2f} MB; scan exit status {r.returncode}")
sys.exit(r.returncode)
