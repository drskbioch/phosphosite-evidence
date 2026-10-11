"""Extension of check_bindings.py to the Rev 6 figures and tables (imported at the end of check_bindings.py; also runnable alone).

1. Every number that a figure or table script formatted from a tracked result value (results/figure_bindings.json: text, file, key path, scale) is re-read
   from the file and recomputed.
2. results/derived_figure_values.json is recomputed by derive_figure_values.compute() and compared with the stored file.
3. results/figure_counts.json: the input hashes are recomputed and the cross-check against dataset_counts.json must show no difference.
4. Text cells (licences, grades, Rev 5 to Rev 6 table, booleans, reasons) are compared with their source text.
5. The text of every SVG is scanned: every number in it must be a bound number of that figure, a tick label, a number written in a key name of a stored
   file (restriction years), or a definition listed below. Anything else is reported as unsourced.
Writes results/figure_binding_check.md. Target: 0 unsourced.
"""
import hashlib
import html
import json
import re
from pathlib import Path

import check_numbers_sourced as N
import derive_figure_values as DV

CUR = Path(__file__).resolve().parents[2]
R6 = CUR / "analysis" / "results"
FIG = CUR / "manuscript" / "figures"  # the composites; each panel is drawn in one of them, so their text covers every panel
BIND = json.load(open(R6 / "figure_bindings.json", encoding="utf-8"))
_cache = {}
# numbers that are definitions written in the figure text, not results; each with its source
DEFINITIONS = {
    "45,000,000": "identifier above this value is invalid: definition in scripts/earlier_figures/compute_results.py (`num > 45_000_000`)",
}
# phrases of figure text whose numbers are definitions, names or layout constants, with the source of the definition
PHRASES = [
    (r"cluster thresholds 500 and 1,000 sites", "study-cluster thresholds of block G13 (plan; variation cluster_threshold=500 and 1000 in cem_results.json)"),
    (r"EPSD 2\.0|PaxDb 6\.0", "proper names of the resources"),
    (r"\d+:\d+", "matching ratio of the 1:1 caliper matcher (plan, design D2b)"),
    (r"scores? (?:of )?2 or more|score 2 or more", "database agreement score of at least 2 (DAS >= 2 in the plan)"),
    (r"\(?(?:agreement )?score 0 to [23]\)?|\(0 to [23]\)", "range of the agreement score under source set A (0 to 2) and B (0 to 3)"),
    (r"5th to 95th percentile", "whisker definition of the box plots"),
    (r"\b9[05]% (?:gene-cluster |)interval", "interval level: 95% (primary) and 90% (equivalence test), as in the plan"),
    (r"\b95% gene-cluster|\b95% intervals|\b90% interval|\bthe 95%|the 90%", "interval level: 95% (primary) and 90% (equivalence test), as in the plan"),
    (r"deciles 8 to 10 minus deciles 1 to 3", "definition of the predictor gap (top three minus bottom three deciles, plan)"),
    (r"axis begins\s+in 1998", "axis range of the plot (layout); the text may wrap between the words"),
    (r"^in 1998\.$", "axis range of the plot (layout): last line of the wrapped sentence 'the axis begins in 1998'"),
    (r"\b95%\s*$", "interval level 95% (primary), at the end of a wrapped caption line whose next line starts with \"gene-cluster intervals\""),
    (r"per 100,000 sites", "unit of a rate"),
    (r"(?:a score of |score )2", "database agreement score of at least 2"),
]
PCT_SCALE = (1, 100, 1000, 100000)


def navigate(o, path):
    if not path:
        return o
    segs = path.split(".")
    for k in range(1, len(segs) + 1):
        key = ".".join(segs[:k])
        if isinstance(o, dict) and key in o:
            r = navigate(o[key], ".".join(segs[k:]))
            if r is not None:
                return r
        if isinstance(o, list) and k == 1 and segs[0].isdigit() and int(segs[0]) < len(o):
            r = navigate(o[int(segs[0])], ".".join(segs[1:]))
            if r is not None:
                return r
    return None


def field(rel, path):
    if rel not in _cache:
        _cache[rel] = json.load(open(CUR / rel, encoding="utf-8"))
    return navigate(_cache[rel], path)


def verify(text, rel, path, scale):
    try:
        x, dec = N.parse(str(text))
    except ValueError:
        return False
    v = field(rel, path)
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return False
    return abs(abs(v * scale) - abs(x)) <= 0.5 * 10 ** (-dec) + 1e-9


def source_text(rel):
    p = CUR / rel
    if rel.startswith("audit_reports") or rel.startswith("analysis/"):
        return re.sub(r"\s+", " ", p.read_text(encoding="utf-8", errors="replace"))
    return ""


lines = []
ok_n = bad_n = 0
bad = []
printed = {}
for name, rec in BIND.items():
    pset = set()
    for text, rel, path, scale in rec["entries"]:
        if verify(text, rel, path, scale):
            ok_n += 1
            try:
                pset.add(N.parse(str(text)))
            except ValueError:
                pass
        else:
            bad_n += 1
            bad.append(f"- {name}: `{text}` = {rel} : {path} x {scale} does not reproduce")
    printed[name] = pset

# ---- 2. derived values
stored = json.load(open(R6 / "derived_figure_values.json", encoding="utf-8"))
fresh = DV.compute()
der_bad = [k for k in fresh if k not in stored or abs(fresh[k]["value"] - stored[k]["value"]) > 1e-12 or fresh[k]["formula"] != stored[k]["formula"]] + [k for k in stored if k not in fresh]
# ---- 3. counts file
fc = json.load(open(R6 / "figure_counts.json", encoding="utf-8"))
hash_bad = []
hash_missing = []
for rel, info in fc["inputs"].items():
    p = CUR / rel
    if not p.exists():
        hash_missing.append(rel)
        continue
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    if h.hexdigest() != info["sha256"]:
        hash_bad.append(rel)
if hash_missing:
    HASH_TEXT = f"input hashes: skipped: record-level data not present ({len(hash_missing)} of {len(fc['inputs'])} files missing; {len(fc['inputs']) - len(hash_missing)} checked, {len(hash_bad)} changed)"
else:
    HASH_TEXT = f"input hashes recomputed, {len(hash_bad)} changed"
cc = fc["cross_check_summary"]
# ---- 4. text cells
tx_ok = tx_bad = 0
tx_list = []
norm = lambda s: re.sub(r"\s+", " ", s).strip()
WORDING = json.load(open(R6.parent / "table_wording.json", encoding="utf-8"))  # cells whose wording differs from the stored text, with the reason
b4 = json.load(open(R6 / "b4_twoway.json", encoding="utf-8"))
nonp = json.load(open(R6 / "g23_nonprimary_pmids.json", encoding="utf-8"))
cache_txt = {}
for name, rec in BIND.items():
    for cell, rel, key in rec.get("texts", []):
        good = False
        if rel.endswith("b4_twoway.json"):
            cid, fld = key.split(".")
            good = (cell == "yes") == bool(b4[cid][fld]) and cell in ("yes", "no")
        elif rel.endswith("g23_nonprimary_pmids.json"):
            pm = key.split("reasons_b.")[1].split(" ")[0]
            good = (cell == nonp["reasons_b"].get(pm)) or (cell == WORDING.get("analysis/results/g23_nonprimary_pmids.json:reasons_b", {}).get(nonp["reasons_b"].get(pm), {}).get("table")) or (cell == "review-type publication" and any(x["pmid"] == pm and x["criterion_a_review_type"] for x in nonp["nonprimary"])) or cell == ""
        elif rel.endswith("manifest_status.csv"):
            good = norm(cell) in norm((R6 / "manifest_status.csv").read_text(encoding="utf-8"))
        elif rel.endswith("CLAIMS_VS_EVIDENCE.md") and "Strength column" in key:
            cid = key.split(",")[0]
            row = next(l for l in (CUR / "analysis" / "CLAIMS_VS_EVIDENCE.md").read_text(encoding="utf-8").split("\n") if l.startswith(f"| {cid} |"))
            good = norm(row.strip().strip("|").split(" | ")[6]).startswith(norm(cell))
        elif rel.endswith(".json"):
            if rel not in cache_txt:
                cache_txt[rel] = norm(json.dumps(json.load(open(CUR / rel, encoding="utf-8")), ensure_ascii=False))
            good = norm(cell) in cache_txt[rel]
        elif rel.endswith(".md") or rel.endswith(".csv"):
            if rel not in cache_txt:
                cache_txt[rel] = norm((CUR / rel).read_text(encoding="utf-8", errors="replace"))
            good = norm(cell).strip() in cache_txt[rel]
        tx_ok += good
        tx_bad += not good
        if not good:
            tx_list.append(f"- {name}: `{cell[:90]}` from {rel} ({key})")
# ---- 5. SVG scan
NUMTOK = re.compile(r"(?<![\w.])[-+−]?\d[\d,]*(?:\.\d+)?")
ID_SKIP = [p for p in N.SKIP]
key_nums = set()
g11 = json.load(open(R6 / "g11_predictors.json", encoding="utf-8"))
for k, v in g11.items():
    if isinstance(v, dict):
        for kk in v:
            for m in re.finditer(r"\d{4}", kk):
                key_nums.add(float(m.group()))
unb = []
n_checked = 0


def scan_svg(svg, label, ticks, pset):
    """Every number in the text of one SVG must be a bound number of its figure, a tick label, a restriction-year key, or a definition."""
    global n_checked
    t = svg.read_text(encoding="utf-8")
    texts = [html.unescape(re.sub(r"<[^>]+>", "", m)) for m in re.findall(r"<text[^>]*>(.*?)</text>", t, flags=re.S)]
    ticks = set(x.replace("−", "-") for x in ticks)
    for tx_ in texts:
        s = tx_.strip()
        if not s or s.replace("−", "-") in ticks:
            continue
        if re.fullmatch(r"10[\d−-]+", re.sub(r"\s+", "", s)):
            continue  # mathtext tick label of a log axis
        s2 = re.sub(r'"[^"]{12,}"', lambda m: " " * len(m.group()), s)
        for pat, _ in PHRASES:
            s2 = re.sub(pat, lambda m: " " * len(m.group()), s2)
        for p in ID_SKIP:
            s2 = p.sub(lambda m: " " * len(m.group()), s2)
        for m in NUMTOK.finditer(s2):
            tok = m.group()
            if tok in DEFINITIONS:
                continue
            try:
                x, dec = N.parse(tok.replace("−", "-"))
            except ValueError:
                continue
            n_checked += 1
            hit = any(abs(abs(px) - abs(x)) < 1e-9 and pd_ == dec for px, pd_ in pset) or (dec == 0 and float(x) in key_nums)
            if not hit:
                # a printed value may carry more decimals in the label than in the registered text (value x 100 etc.): compare the rounded value
                hit = any(abs(abs(px) - abs(x)) <= 0.5 * 10 ** (-dec) + 1e-9 and pd_ >= dec for px, pd_ in pset)
            if not hit:
                unb.append(f"- {label}: `{tok}` in text: {s[:100]}")


for name in sorted(BIND):
    svg = FIG / f"{name}.svg"
    if svg.exists():
        scan_svg(svg, name, BIND[name]["ticks"], printed[name])
# ---- 6. table scan: every number in a table cell is a bound number of that table, or sits in a verified text cell
import csv as _csv
TAB = CUR / "manuscript" / "tables"
tb_unb, tb_n = [], 0
for name in sorted(BIND):
    if not name.startswith("Table"):
        continue
    fs = sorted(TAB.glob(name + "_*.csv"))
    if not fs:
        continue
    txts = [c for c, _, _ in BIND[name].get("texts", [])]
    pset = printed[name]
    exact = {norm(t) for t in txts}
    for ri_, row in enumerate(_csv.reader(open(fs[0], encoding="utf-8"))):
        if ri_ == 0:
            continue  # column headers (interval level 95% is a definition)
        for cell in row:
            if norm(cell) in exact:
                continue
            s2 = cell
            for t in sorted(txts, key=len, reverse=True):
                if len(norm(t)) > 3 and norm(t) in norm(s2):
                    s2 = s2.replace(t, " " * len(t))
            for pat, _ in PHRASES:
                s2 = re.sub(pat, lambda m: " " * len(m.group()), s2)
            for p_ in ID_SKIP:
                s2 = p_.sub(lambda m: " " * len(m.group()), s2)
            s2 = re.sub(r"\d{4}-\d{2}-\d{2}|[0-9a-f]{7}\b", lambda m: " " * len(m.group()), s2)
            for m in NUMTOK.finditer(s2):
                tok = m.group()
                try:
                    x, dec = N.parse(tok.replace("−", "-"))
                except ValueError:
                    continue
                tb_n += 1
                if not any(abs(abs(px) - abs(x)) <= 0.5 * 10 ** (-dec) + 1e-9 and pd_ >= dec for px, pd_ in pset):
                    tb_unb.append(f"- {name}: `{tok}` in cell: {cell[:90]}")
# Table S1 (contrast grid) against manifest_status.csv: the identifying columns are equal, the estimate, the interval limits and the attenuation equal the stored values at 4 significant digits;
# the wording columns (design, population, index, variation, tier rule, notes) are descriptive text and may be rephrased
import math as _math
_ID_COLS = ("contrast_id", "block", "task", "signal", "signal_label", "decision_bearing", "status", "result_status")
_NUM_COLS = ("estimate", "ci_low", "ci_high", "attenuation")
_ms = {r["contrast_id"]: r for r in _csv.DictReader(open(R6 / "manifest_status.csv", encoding="utf-8"))}
_s1 = list(_csv.DictReader(open(TAB / "TableS1_contrast_grid.csv", encoding="utf-8")))


def _sig4(x):
    x = float(x)
    return 0.0 if x == 0 else round(x, 3 - int(_math.floor(_math.log10(abs(x)))))


s1_bad = 0 if len(_s1) == len(_ms) else 1
for _r in _s1:
    _m = _ms.get(_r["contrast_id"])
    if _m is None:
        s1_bad += 1
        continue
    s1_bad += sum(_r[c] != _m[c] for c in _ID_COLS)
    for c in _NUM_COLS:
        if (_r[c] == "") != (_m[c] == ""):
            s1_bad += 1
        elif _r[c] != "" and abs(float(_r[c]) - _sig4(_m[c])) > 1e-9:
            s1_bad += 1
ms_same = s1_bad == 0
unb = unb + tb_unb
n_checked += tb_n
total_bound = ok_n
head = (f"# Binding of the numbers in the Rev 6 figures and tables{chr(10)}{chr(10)}"
        f"{len(BIND)} figures and tables. **Numbers formatted from tracked result values and recomputed here from the named file and key: {ok_n} verified, {bad_n} not reproduced.** "
        f"Derived values recomputed (`derive_figure_values.py`): {len(fresh)} values, {len(der_bad)} differ. Counts file: {HASH_TEXT}; cross-check against dataset_counts.json: {cc['equal']} of {cc['checked']} equal, differences: {cc['different'] or 'none'}. "
        f"Text cells re-verified against their source text: {tx_ok} verified, {tx_bad} not. Composite SVGs scanned (they contain every panel). SVG text and table scan ({f"Table S1 has the {len(_s1)} contrasts of manifest_status.csv with equal identifiers and values at 4 significant digits" if ms_same else f"TABLE S1 DIFFERS FROM manifest_status.csv in {s1_bad} fields"}): {n_checked} numbers checked, **{len(unb)} unsourced**.{chr(10)}{chr(10)}")
body = ("## Not reproduced" + chr(10) + chr(10) + (chr(10).join(bad) if bad else "none") + chr(10) + chr(10) + "## Derived values that differ" + chr(10) + chr(10) + (", ".join(der_bad) if der_bad else "none") +
        chr(10) + chr(10) + "## Text cells not verified" + chr(10) + chr(10) + (chr(10).join(tx_list) if tx_list else "none") + chr(10) + chr(10) + "## Unsourced numbers in SVG text" + chr(10) + chr(10) + (chr(10).join(unb) if unb else "none") +
        chr(10) + chr(10) + "## Definitions, names and layout constants written in figure text (not results)" + chr(10) + chr(10) + chr(10).join(f"- `{k}`: {v}" for k, v in DEFINITIONS.items()) + chr(10) + chr(10).join(f"- phrase `{p}`: {v}" for p, v in PHRASES) + chr(10))
(R6 / "figure_binding_check.md").write_text(head + body, encoding="utf-8")
print(f"figures and tables: {len(BIND)}; bound numbers verified {ok_n}, not reproduced {bad_n}; derived differ {len(der_bad)}; {HASH_TEXT}; counts cross-check {cc['equal']}/{cc['checked']}; "
      f"text cells {tx_ok} verified, {tx_bad} not; SVG numbers checked {n_checked}, unsourced {len(unb)}")
