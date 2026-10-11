"""Bind every number in the main tables (1 to 3) and the supplementary tables (S1 to S5), with their notes, to its source.

Each number printed in a table cell or in a notes file is bound, in this order, to
  1. a result value that the table builder formatted (results/figure_bindings.json: file, key path, scale), re-read and recomputed by check_figure_bindings.py;
  2. a cell of text verified against its source text (licences, grades, reasons: results/figure_bindings.json "texts");
  3. a binding of a claims-table row (results/claims_bindings.json);
  4. a count derived from results/manifest_status.csv (Table S1: contrasts per block and task) or a field of earlier_results/dataset_counts.json;
  5. a documented constant (interval level, number of resamples, null values), with the file and line that state it;
  6. a definition or an identifier written in the text (phrases of check_figure_bindings.py, identifiers, dates).
Table S1 (contrast grid) is compared with results/manifest_status.csv: identifiers, estimate, interval limits and attenuation (4 significant digits), and the
second interval, which must be the interval of the paired contrast (the contrast of the source block with the same records).
A number that matches none of these is listed as unbound. If it equals a number of the earlier estimates (Rev 5; results/internal changes table, and the dataset
counts 300,038 sites and 44,112 proteins), it is listed as superseded with the Rev 6 value beside it.
The comparison with the earlier estimates reads analysis/internal/earlier_to_present_changes.csv, a file that is not part of the release snapshot. When it is absent the script prints
"SUPERSEDED-VALUE COMPARISON NOT RUN: internal file not in this tree" (in the output, in the report and in the summary line), runs every other test, and exits 0 only if
those pass. With --require-internal the missing file is an error (exit 2); the repository-level run uses this flag.
Usage: python check_table_numbers.py [--require-internal]        Writes results/table_number_check.md and results/table_number_binding.csv. Exit status 1 if a number is unbound.
"""
import csv
import json
import math
import re
import sys
from pathlib import Path

import check_figure_bindings as F  # runs the figure and table check once and gives its bindings, phrases and helpers
import check_numbers_sourced as N

CUR = F.CUR
R6 = F.R6
TAB = F.TAB
NORM = F.norm
CONSTANTS = {  # (value, decimals) -> source of the definition; value and decimals are those printed
    (95.0, 0): "interval level of the 95% intervals: analysis/ANALYSIS_PLAN.md, task J (line 205: 95% intervals shown)",
    (1000.0, 0): "number of bootstrap resamples: analysis/ANALYSIS_PLAN.md line 7 (1,000 resamples unless stated) and decision D7",
    (0.0, 0): "null value of a difference: definition (GOLD minus BRONZE; no difference is 0)",
    (1.0, 0): "null value of an odds ratio: definition (no difference is 1)",
}
REV5_COUNTS = {(300038.0, 0): "sites of the earlier dataset", (44112.0, 0): "proteins of the earlier dataset"}

# ---- sources of binding
prov = {}  # table -> {(abs value, decimals): (file, key)}
for name, rec in F.BIND.items():
    d = prov.setdefault(name, {})
    for text, rel, path, scale in rec["entries"]:
        try:
            x, dec = N.parse(str(text))
        except ValueError:
            continue
        if F.verify(text, rel, path, scale):
            d.setdefault((abs(x), dec), (rel, path))
claims = {}
for cid, lst in json.load(open(R6 / "claims_bindings.json", encoding="utf-8")).items():
    for text, rel, path, scale in lst:
        try:
            x, dec = N.parse(str(text))
        except ValueError:
            continue
        if F.verify(text, rel, path, scale):
            claims.setdefault((abs(x), dec), (rel, path, cid))
dc = json.load(open(CUR / "earlier_results" / "dataset_counts.json", encoding="utf-8"))
ms_rows = list(csv.DictReader(open(R6 / "manifest_status.csv", encoding="utf-8")))
ms_counts = {}
for r in ms_rows:
    ms_counts[("block", r["block"])] = ms_counts.get(("block", r["block"]), 0) + 1
    ms_counts[("task", r["task"])] = ms_counts.get(("task", r["task"]), 0) + 1
rev5 = {}
INTERNAL = CUR / "analysis" / "internal" / "earlier_to_present_changes.csv"
NOTICE = "SUPERSEDED-VALUE COMPARISON NOT RUN: internal file not in this tree"
REQUIRE_INTERNAL = "--require-internal" in sys.argv
COMPARISON_RUN = INTERNAL.exists()
if not COMPARISON_RUN:
    if REQUIRE_INTERNAL:
        print(f"ERROR: --require-internal given and {INTERNAL.relative_to(CUR.parent)} is missing; the comparison with the earlier estimates cannot run")
        sys.exit(2)
    print(NOTICE)
else:
    for r in csv.DictReader(open(INTERNAL, encoding="utf-8")):
        for tok in F.NUMTOK.finditer(r["Rev 5 value (source file)"]):
            try:
                x, dec = N.parse(tok.group())
            except ValueError:
                continue
            rev5.setdefault((abs(x), dec), (r["ID"], r["Rev 6 value"]))

REC = []  # table, place, token, bound by, source, detail


def match(pset, x, dec):
    for (px, pd_), src in pset.items():
        if abs(px - abs(x)) <= 0.5 * 10 ** (-dec) + 1e-9 and pd_ >= dec:
            return src
    return None


def scan(name, place, cell, texts, pset):
    """Numbers of one cell or one notes line."""
    exact = {NORM(t): (rel, key) for t, rel, key in texts}
    s2 = cell
    if NORM(cell) in exact:
        rel, key = exact[NORM(cell)]
        for m in F.NUMTOK.finditer(cell):
            REC.append((name, place, m.group(), "text cell verified against its source", rel, key))
        return
    for t, rel, key in sorted(texts, key=lambda z: len(z[0]), reverse=True):
        if len(NORM(t)) > 3 and NORM(t) in NORM(s2):
            for m in F.NUMTOK.finditer(t):
                REC.append((name, place, m.group(), "text cell verified against its source", rel, key))
            s2 = s2.replace(t, " " * len(t))
    for pat, why in F.PHRASES:
        def blank(m, why=why):
            for q in F.NUMTOK.finditer(m.group()):
                REC.append((name, place, q.group(), "definition written in the text", "scripts/analysis/check_figure_bindings.py (PHRASES)", why))
            return " " * len(m.group())
        s2 = re.sub(pat, blank, s2)
    for p_ in F.ID_SKIP:
        s2 = p_.sub(lambda m: " " * len(m.group()), s2)
    s2 = re.sub(r"\d{4}-\d{2}-\d{2}|[0-9a-f]{7}", lambda m: " " * len(m.group()), s2)
    for m in F.NUMTOK.finditer(s2):
        tok = m.group()
        try:
            x, dec = N.parse(tok.replace("−", "-"))
        except ValueError:
            continue
        hit = match(pset, x, dec)
        if hit:
            REC.append((name, place, tok, "result value of the table builder", hit[0], hit[1]))
            continue
        hit = match(claims, x, dec)
        if hit:
            REC.append((name, place, tok, "claims-table binding", hit[0], f"{hit[1]} (claim {hit[2]})"))
            continue
        if dec == 0 and abs(x) in (dc.get("n_sites"), dc.get("n_proteins")):
            REC.append((name, place, tok, "dataset count", "earlier_results/dataset_counts.json", "n_sites" if abs(x) == dc.get("n_sites") else "n_proteins"))
            continue
        if (abs(x), dec) in CONSTANTS:
            REC.append((name, place, tok, "documented constant", "see detail", CONSTANTS[(abs(x), dec)]))
            continue
        REC.append((name, place, tok, "UNBOUND", "", ""))


S1_NUM = ("estimate", "ci_low", "ci_high", "second_ci_low", "second_ci_high", "attenuation")


def sig4(x):
    x = float(x)
    return 0.0 if x == 0 else round(x, 3 - int(math.floor(math.log10(abs(x)))))


def table_s1_grid():
    ms = {r["contrast_id"]: r for r in ms_rows}
    rows = list(csv.DictReader(open(TAB / "TableS1_contrast_grid.csv", encoding="utf-8")))
    bad = 0
    for r in rows:
        m = ms.get(r["contrast_id"])
        if m is None:
            bad += 1
            continue
        for c in ("estimate", "ci_low", "ci_high", "attenuation"):
            if r[c] == "":
                bad += m[c] != ""
                continue
            ok = m[c] != "" and abs(float(r[c]) - sig4(m[c])) <= 1e-9
            bad += not ok
            REC.append(("TableS1_contrast_grid", f"{r['contrast_id']} {c}", r[c], "manifest_status.csv" if ok else "UNBOUND", "analysis/results/manifest_status.csv" if ok else "", f"{c} at 4 significant digits" if ok else ""))
    pair_bad = 0
    ci = {r["contrast_id"]: (r["ci_low"], r["ci_high"]) for r in rows}
    sec = {r["contrast_id"]: (r["second_ci_low"], r["second_ci_high"]) for r in rows}
    for r in rows:
        if r["second_ci_low"] == "":
            continue
        # the paired contrast: its interval is this row's second interval and its second interval is this row's interval
        partners = [k for k in ci if k != r["contrast_id"] and ci[k] == sec[r["contrast_id"]] and sec[k] == ci[r["contrast_id"]]]
        ok = bool(partners)
        # and the limits are those of the manifest rows at 4 significant digits
        mp = [ms[k] for k in partners]
        ok = ok and any(abs(float(r["second_ci_low"]) - sig4(x["ci_low"])) <= 1e-9 and abs(float(r["second_ci_high"]) - sig4(x["ci_high"])) <= 1e-9 for x in mp)
        pair_bad += not ok
        for c in ("second_ci_low", "second_ci_high"):
            REC.append(("TableS1_contrast_grid", f"{r['contrast_id']} {c}", r[c], "paired contrast of manifest_status.csv" if ok else "UNBOUND",
                        "analysis/results/manifest_status.csv" if ok else "", f"interval of {', '.join(partners[:1])}" if ok else ""))
    return len(rows), bad, pair_bad


PLAN = (CUR / "analysis" / "ANALYSIS_PLAN.md").read_text(encoding="utf-8").splitlines()


def plan_line(tok, code):
    """First line of the analysis plan that contains the number and the code (or its first word) that the key row explains."""
    word = re.split(r"[ (/]", code)[0]
    if tok in word:  # the number is part of the name of the index (V20, V100, VT2025, VT2026): the plan line that defines the name
        return next((i for i, line in enumerate(PLAN, start=1) if word.lower() in line.lower()), None)
    for i, line in enumerate(PLAN, start=1):
        if re.search(r"(?<![\w.])" + re.escape(tok.replace(",", "")) + r"(?![\w.])", line.replace(",", "")) and word.lower() in line.lower():
            return i
    return None


def table_s1_key():
    grid = list(csv.DictReader(open(TAB / "TableS1_contrast_grid.csv", encoding="utf-8")))
    for i, r in enumerate(csv.DictReader(open(TAB / "TableS1_key.csv", encoding="utf-8"))):
        place = f"{r['Column or code']} {r['Value']}"
        if r["Count"]:
            n = int(r["Count"].split()[0])
            kind = r["Column or code"]
            ok = kind in grid[0] and sum(1 for g in grid if g[kind] == r["Value"]) == n
            REC.append(("TableS1_key", place, r["Count"], "count of rows of the Table S1 grid (equal to manifest_status.csv)" if ok else "UNBOUND",
                        "analysis/results/manifest_status.csv" if ok else "", f"rows with {kind} = {r['Value']}" if ok else ""))
        for c in ("Value", "Meaning"):
            if c == "Value" and r["Count"]:
                continue
            n0 = len(REC)
            scan("TableS1_key", f"{place} / {c}", r[c], [], {})
            for k in range(n0, len(REC)):
                if REC[k][3] == "UNBOUND":
                    ln = plan_line(REC[k][2], r["Value"] if r["Column or code"] in ("index", "tier_rule", "block") else r["Column or code"] + " " + r["Value"])
                    if ln is None:
                        ln = plan_line(REC[k][2], r["Value"])
                    if ln is None and "reference year" in r[c].lower():
                        ln = plan_line(REC[k][2], "reference year")
                    if ln:
                        REC[k] = REC[k][:3] + ("documented definition", "analysis/ANALYSIS_PLAN.md", f"line {ln}")


s1_n, s1_bad, s1_pair_bad = table_s1_grid()
table_s1_key()
for csvf in sorted(TAB.glob("*.csv")):
    name = csvf.name.rsplit(".", 1)[0]
    base = name
    if base in ("TableS1_contrast_grid", "TableS1_key"):
        continue
    # the name used in figure_bindings.json is the part before the first underscore that starts the description
    bname = next((b for b in F.BIND if name == b or name.startswith(b + "_")), None)
    texts = F.BIND[bname].get("texts", []) if bname else []
    pset = prov.get(bname, {})
    rows = list(csv.reader(open(csvf, encoding="utf-8")))
    head = rows[0]
    for ri, row in enumerate(rows[1:], start=1):
        for ci_, cell in enumerate(row):
            scan(name, f"row {ri} / {head[ci_]}", cell, texts, pset)
    for ci_, h in enumerate(head):
        scan(name, f"column header {ci_ + 1}", h, [], {})
    nf = csvf.with_name(name + "_notes.md")
    if nf.exists():
        for li, line in enumerate(nf.read_text(encoding="utf-8").splitlines(), start=1):
            scan(name + "_notes", f"line {li}", line, texts, pset)
    if name.startswith("TableS3") or name.startswith("Table3") or name.startswith("Table1"):
        pass
for nf in sorted(TAB.glob("*_notes.md")):
    nm = nf.name[:-len("_notes.md")]
    if not (TAB / (nm + ".csv")).exists():
        for li, line in enumerate(nf.read_text(encoding="utf-8").splitlines(), start=1):
            scan(nm + "_notes", f"line {li}", line, [], {})

# ---- classify unbound numbers
unb = [r for r in REC if r[3] == "UNBOUND"]
sup = []
oth = []
for r in unb:
    try:
        x, dec = N.parse(r[2].replace("−", "-"))
    except ValueError:
        oth.append(r)
        continue
    k = (abs(x), dec)
    if k in rev5:
        sup.append((r, f"equals the Rev 5 value of {rev5[k][0]}; Rev 6 value: {rev5[k][1]}"))
    elif k in REV5_COUNTS:
        sup.append((r, f"equals the earlier count ({REV5_COUNTS[k]}); Rev 6: {dc.get('n_sites') if 'sites' in REV5_COUNTS[k] else dc.get('n_proteins'):,}"))
    else:
        oth.append(r)

tables = sorted({r[0] for r in REC})
kinds = sorted({r[3] for r in REC})
NL = "\n"
out = ["# Check of the numbers in the tables", "",
       f"Tables checked: {len([t for t in tables if not t.endswith('_notes')])} files (Tables 1 to 3 and S1 to S5) and their notes files. Numbers examined: {len(REC)}. **Unbound: {len(unb)}** (superseded: {len(sup)}; other: {len(oth)}).",
       *([] if COMPARISON_RUN else [f"**{NOTICE}.** The other tests were run; the numbers in the lists below that equal a value of the earlier estimates are reported as unbound without the comparison, not as superseded."]),
       f"Table S1 (contrast grid): {s1_n} rows compared with `results/manifest_status.csv`; fields that differ: {s1_bad}; second intervals without a paired contrast: {s1_pair_bad}.", "",
       "Binding of a number, in this order: result value that the table builder formatted (file and key re-read and recomputed by `check_figure_bindings.py`); text cell verified against its source text; claims-table binding (`claims_bindings.json`); count derived from `manifest_status.csv` or `dataset_counts.json`; documented constant; definition or identifier written in the text. The list of every number with its binding is `table_number_binding.csv` in this folder.", "",
       "## Numbers by table and kind of binding", "", "| Table | " + " | ".join(kinds) + " | total |", "|---|" + "---|" * (len(kinds) + 1)]
for t in tables:
    c = [sum(1 for r in REC if r[0] == t and r[3] == k) for k in kinds]
    out.append(f"| {t} | " + " | ".join(map(str, c)) + f" | {sum(c)} |")
out += ["", "## Unbound numbers", ""]
if not unb:
    out.append("none")
else:
    out += ["| Table | Place | Number | Finding |", "|---|---|---|---|"]
    for r, why in sup:
        out.append(f"| {r[0]} | {r[1]} | {r[2]} | superseded: {why} |")
    for r in oth:
        out.append(f"| {r[0]} | {r[1]} | {r[2]} | no source found |")
out += ["", "## Documented constants used", ""] + [f"- `{k[0]:g}`: {v}" for k, v in CONSTANTS.items() if any(r[3] == "documented constant" and N.parse(r[2])[0] == k[0] for r in REC)]
Path(R6 / "table_number_check.md").write_text(NL.join(out) + NL, encoding="utf-8")
with open(R6 / "table_number_binding.csv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["table", "place", "number", "bound_by", "source_file", "key_or_detail"])
    w.writerows(REC)
print(f"{len(REC)} numbers in {len(tables)} table and notes files; unbound {len(unb)} (superseded {len(sup)}, other {len(oth)}); Table S1: {s1_n} rows, {s1_bad} fields differ, {s1_pair_bad} unpaired second intervals"
      + ("" if COMPARISON_RUN else f"; {NOTICE}"))
sys.exit(1 if unb or s1_bad or s1_pair_bad else 0)
