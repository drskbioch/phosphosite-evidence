"""List every numeric literal that the three report builders write into the document text, with its category and where it is documented.

The builders (build_results.py, build_claims.py, build_interim.py) are parsed with `ast`; the constant parts of every string (f-string expressions
are not text) are scanned for numbers, after identifiers, code spans, dates, section references and quoted PubMed titles are removed.
Categories: statistical convention (95% / 90% intervals); design definition (a threshold or setting written in the plan); provenance (publication year, release,
build, taxon or HTTP code); quoted earlier statement (a number from the earlier analysis that a row of the withdrawn-or-changed table or the interim report cites
as it was written); list numbering. A literal is *documented* when the same number appears in the plan, RECONCILIATION.md, sources/SOURCES.md, INTERIM_REPORT.md or a
result file. Anything that is not in a category is listed as UNCLASSIFIED and must be computed from a result file or documented.
Usage: python check_literals.py   -> results/builder_numeric_literals.md
"""
import ast
import glob
import re
import sys
from pathlib import Path

import check_numbers_sourced as N

REV6 = N.REV6
HERE = Path(__file__).resolve().parent
BUILDERS = ["build_results.py", "build_claims.py", "build_interim.py"]
DOCS = {n: (REV6 / n).read_text(encoding="utf-8", errors="ignore") for n in ("ANALYSIS_PLAN.md", "RECONCILIATION.md", "INTERIM_REPORT.md")}
DOCS["scripts/analysis/check_cond_logit.py"] = (HERE / "check_cond_logit.py").read_text(encoding="utf-8", errors="ignore")
DOCS["sources/SOURCES.md"] = (REV6 / "sources" / "SOURCES.md").read_text(encoding="utf-8", errors="ignore") if (REV6 / "sources" / "SOURCES.md").exists() else ""
LOGS = " ".join(open(f, encoding="utf-8", errors="ignore").read() for f in glob.glob(str(REV6 / "results" / "*.log")))

PROV = {"2015", "2020", "2022", "2023", "2025", "2026", "606", "6.0", "9606", "10090", "10116", "403", "1113"}
STAT = re.compile(r"\b(95|90)%|\b(95|90)% ?(CI|interval|equivalence)")
DEF_TOKENS = {"0.2", "20", "50", "7", "35", "0.01", "0.05", "0.1", "10", "3", "5", "2", "31", "40", "0.80", "0.95", "0.99", "0.19", "884", "500", "1,000", "1000", "200", "3,000",
              "5e-04", "1e-04", "1e-9", "-0.05", "1", "0", "-3", "-2", "+1", "+3", "4", "6", "8", "9", "12", "11", "690",
              "2017", "2018"}  # restriction years of the predictor sensitivity rows (plan section 11: 2019 and 2018 for DeepPhos, 2020 and 2017 for MusiteDeep; predictor specification)
ORDINAL = re.compile(r"^\s*\*{0,2}(Limitations to carry forward|Decisions for you|Needed from you)")
QUOTED = re.compile(r"Rev 5|interim|Table 1|text quotes|The study-cluster threshold|inherited")


def documented(tok):
    t = tok.rstrip(",")
    for name, txt in DOCS.items():
        if t in txt:
            return name
    if t in LOGS:
        return "results/*.log"
    return None


rows, unclassified = [], []
SEC7 = next(n.lineno for n in ast.walk(ast.parse((HERE / "build_claims.py").read_text(encoding="utf-8"))) if isinstance(n, ast.Constant) and isinstance(n.value, str) and n.value.startswith("## 7."))
for fn in BUILDERS:
    src = (HERE / fn).read_text(encoding="utf-8")
    for node in ast.walk(ast.parse(src)):
        if not (isinstance(node, ast.Constant) and isinstance(node.value, str)) or len(node.value) < 3:
            continue
        v = node.value
        s = re.sub(r'"[^"]{12,}"', lambda m: " " * len(m.group()), v)
        for p in N.SKIP:
            s = p.sub(lambda m: " " * len(m.group()), s)
        s = re.sub(r"\d+:\d+", lambda m: " " * len(m.group()), s)  # n:n ratios (matching design)
        for m in N.NUM.finditer(s):
            tok = m.group().rstrip(",")
            if tok == "":
                continue
            ctx = v.strip().replace(chr(10), " ")[:110]
            if re.fullmatch(r"0\d+", tok):
                cat = "format spec or identifier fragment"
            elif STAT.search(v) and tok in ("95", "90"):
                cat = "statistical convention"
            elif ORDINAL.search(v) or (tok.isdigit() and len(tok) == 1 and re.search(r"\(%s\)" % tok, v)):
                cat = "list numbering"
            elif tok in PROV:
                cat = "provenance"
            elif fn == "build_interim.py" or QUOTED.search(v) or (fn == "build_claims.py" and node.lineno >= SEC7):
                cat = "quoted earlier statement"
            elif tok in DEF_TOKENS:
                cat = "design definition"
            else:
                cat = "UNCLASSIFIED"
            doc = documented(tok)
            rows.append((fn, node.lineno, tok, cat, doc, ctx))
            if cat == "UNCLASSIFIED" or doc is None:
                unclassified.append((fn, node.lineno, tok, cat, doc, ctx))

from collections import Counter

cnt = Counter(r[3] for r in rows)
out = [f"# Numeric literals in the report builders' document text\n", f"{len(rows)} literals in {len(BUILDERS)} builders. " + "; ".join(f"{k}: {v}" for k, v in sorted(cnt.items())) + ".\n",
       f"**Unclassified or not found in a documenting file: {len(unclassified)}.**\n"]
for fn, ln, tok, cat, doc, ctx in unclassified:
    out.append(f"- {fn}:{ln} `{tok}` ({cat}; documented in: {doc}) in: {ctx}")
out.append("\n## All literals\n\n| Builder | Line | Number | Category | Documented in | Context |\n|---|---|---|---|---|---|")
for fn, ln, tok, cat, doc, ctx in sorted(rows):
    out.append(f"| {fn} | {ln} | {tok} | {cat} | {doc or 'NOT FOUND'} | {ctx.replace('|', '/')} |")
(REV6 / "results" / "builder_numeric_literals.md").write_text(chr(10).join(out) + chr(10), encoding="utf-8")
print(f"{len(rows)} literals; " + "; ".join(f"{k}: {v}" for k, v in sorted(cnt.items())) + f"; unclassified or undocumented: {len(unclassified)}")
