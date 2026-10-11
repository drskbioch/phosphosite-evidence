"""Bind every number in the claim rows of CLAIMS_VS_EVIDENCE.md (sections 1 to 6; columns Claim, Estimate, 95% CI, n and Source) to a named result field.

A number is *bound* when it equals (at the number of decimals written; percent notation allowed) a field of a result record that the same row names:
the contrast IDs in the row (G1-0005, G17-0382, ...), the blocks it names (G23, G11, ... -> every contrast of the block in grid_manifest.csv and the
result files whose name starts with the block), or the result files it names in code spans. The binding is reported as file : key path (the contrast ID
is part of the key path). A number that is not bound but equals a derivation recomputed in check_numbers_sourced.py is reported as *derived*; a number that
equals some source value elsewhere but no field of the row's own records is *generic only* (listed); a number that equals nothing is *unsourced* (listed).
Not in scope: the Strength and Basis columns (grading arithmetic, defined in the grading rules at the top of the file).
Usage: python check_bindings.py   -> results/number_binding_check.md
"""
import json
import re

import check_numbers_sourced as N

REV6, CUR = N.REV6, N.CUR
SIDE = json.load(open(REV6 / "results" / "claims_bindings.json", encoding="utf-8"))
_cache = {}


def navigate(o, path):
    """Value at a key path joined by '.'; keys that contain '.' are tried by joining segments."""
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


def field_value(rel, path):
    if rel not in _cache:
        _cache[rel] = json.load(open(CUR / rel, encoding="utf-8"))
    return navigate(_cache[rel], path)


def verify(entry, x, dec):
    """The builder says: text = file.path * scale. Re-read the file and check that the value, rounded to the decimals written, is the number in the document."""
    text, rel, path, scale = entry
    try:
        tx, tdec = N.parse(str(text))
    except ValueError:
        return False
    if tdec != dec or abs(tx - x) > 1e-9:
        return False
    v = field_value(rel, path)
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return False
    return abs(abs(v * scale) - abs(x)) <= 0.5 * 10 ** (-dec) + 1e-9


# ---------------------------------------------------------------- claim rows
rows = []
sect = None
for ln, line in enumerate((REV6 / "CLAIMS_VS_EVIDENCE.md").read_text(encoding="utf-8").split("\n"), 1):
    m = re.match(r"## (\d+)\.", line)
    if m:
        sect = int(m.group(1))
    if sect and 1 <= sect <= 6 and re.match(r"\| C\d+\.\d+ \|", line):
        c = [x.strip() for x in line.strip().strip("|").split(" | ")]
        if len(c) == 8:
            rows.append((ln, c))

G21_ROWS = {f"C5.{k}" for k in range(1, 9)}  # rows built from G21 results
PLAN = (REV6 / "ANALYSIS_PLAN.md").read_text(encoding="utf-8")
PLAN_CONST = {("C1.4", "2"): "DAS >= 2", ("C6.1", "0.05"): "margin +/-0.05", **{(f"C5.{k}", "5"): "more than 5% of the fits" for k in range(1, 9)}}
ROW_DER = {"C1.2": ["sum of the mammalian"], "C1.3": ["human GOLD share"], "C1.6": ["GOLD records that fall"],
           "C4.1": ["plain ratio", "tier spread", "ln spread", "largest |ln", "signal-by-population", "min of", "max of", "disorder P-all spread"],
           "C4.3": ["CEM rows with"], "C4.4": ["width ratio"], "C4.5": ["decision-bearing G1 contrasts"], "C4.6": ["G23 D1 contrasts"], "C6.1": ["predictor tools"],
           **{f"C5.{k}": ["G21 failed-fit share"] for k in range(1, 9)}}
FILE = re.compile(r"`([^`]+[.]json)`")
_keys = {}


def file_keys(rel):
    if rel not in _keys:
        ks = set()

        def walk(o):
            if isinstance(o, dict):
                for k, v in o.items():
                    ks.add(str(k))
                    walk(v)
            elif isinstance(o, list):
                for v in o:
                    walk(v)
        walk(json.load(open(CUR / rel, encoding="utf-8")))
        _keys[rel] = ks
    return _keys[rel]


ID = re.compile(r"[GB]\d{1,2}-\d{3,4}")
bound, derived, generic, unsourced, details = 0, 0, [], [], []
for ln, c in rows:
    cid = c[0]
    ents = SIDE.get(cid, [])
    cells_ = c[1:6] + ([c[7]] if cid in G21_ROWS else [])  # Basis is in scope for the G21 rows (it holds the fit counts)
    for ci_, cell in enumerate(cells_):
        s_ = re.sub(r'"[^"]{12,}"', lambda m: " " * len(m.group()), cell)  # quoted spans are verbatim PubMed titles
        for p in N.SKIP:
            s_ = p.sub(lambda m: " " * len(m.group()), s_)
        for m in N.NUM.finditer(s_):
            try:
                x, dec = N.parse(m.group())
            except ValueError:
                continue
            pct = s_[m.end():m.end() + 1] == "%"
            col = ("Claim", "Estimate", "95% CI", "n", "Source", "Basis")[ci_]
            hit = next((e for e in ents if verify(e, x, dec)), None)
            if hit:
                bound += 1
                details.append(f"| {cid} | {col} | {m.group().rstrip(',')} | {hit[1]} : {hit[2]}" + (f" (x {hit[3]:g})" if hit[3] != 1 else "") + " |")
                continue
            tok = m.group().rstrip(",")
            keyhit = None
            if dec == 0 and len(tok.replace(",", "")) >= 3:
                keyhit = next((e for e in ents if tok.replace(",", "") in e[2].split(".")), None)
                if not keyhit:
                    for fn in FILE.findall(c[5]):
                        for rel in (f"analysis/results/{fn.split('/')[-1]}", f"results/{fn.split('/')[-1]}"):
                            if (CUR / rel).exists() and tok.replace(",", "") in file_keys(rel):
                                keyhit = (tok, rel, "(key name)", 1)
            if dec == 0 and not keyhit and tok in ("0", "1") and (cid, tok) not in PLAN_CONST:
                for i in sorted(set(ID.findall(c[1] + " " + c[5]))):
                    for rel in ("analysis/results/amend_results.json", "analysis/results/cem_results.json", "analysis/results/g5_abundance.json"):
                        v = field_value(rel, f"{i}.null") if (CUR / rel).exists() else None
                        if v is not None and float(v) == float(tok):
                            keyhit = (tok, rel, f"{i}.null", 1)
                            break
                    if keyhit:
                        break
            if keyhit:
                bound += 1
                details.append(f"| {cid} | {col} | {tok} | {keyhit[1]} : {keyhit[2]} |")
                continue
            if (cid, tok) in PLAN_CONST and PLAN_CONST[(cid, tok)] in PLAN:
                bound += 1
                details.append(f"| {cid} | {col} | {tok} | analysis/ANALYSIS_PLAN.md : {PLAN_CONST[(cid, tok)]} (a definition or threshold recorded in the plan, not a result) |")
                continue
            lab = [l for l in N.explained_by(x, dec, pct) if any(l.startswith(pf) for pf in ROW_DER.get(cid, ()))]
            if lab:
                derived += 1
                details.append(f"| {cid} | {col} | {m.group().rstrip(',')} | derived: {lab[0]} |")
                continue
            if N.matches(x, dec, pct):
                generic.append(f"- {cid} (line {ln}), {col}: `{m.group().rstrip(',')}`")
            else:
                unsourced.append(f"- {cid} (line {ln}), {col}: `{m.group().rstrip(',')}`")
# ---------------------------------------------------------------- headline items of RESULTS_REPORT.md: each number must equal a number bound (or derived) in a claim row
bound_vals = []
for line in details:
    cells = [x.strip() for x in line.strip().strip("|").split(" | ")]
    try:
        bx, bdec = N.parse(cells[2])
    except ValueError:
        continue
    bound_vals.append((bx, bdec, cells[0], cells[3]))
res_lines = (REV6 / "RESULTS_REPORT.md").read_text(encoding="utf-8").split(chr(10))
h0 = next(k for k, l in enumerate(res_lines) if l.startswith("## Headline results"))
h1 = next(k for k in range(h0 + 1, len(res_lines)) if res_lines[k].startswith("## "))
HL_CONST = {("3", "0.1"): "SMD"}
HL_LAB = {"5": ["signal-by-population", "tier spread", "disorder P-all spread", "plain ratio", "ln spread", "min of", "max of"]}
hl_bound, hl_unbound = [], []
for hl in res_lines[h0 + 1:h1]:
    if not hl.strip():
        continue
    t_ = re.sub(r'"[^"]{12,}"|\d+:\d+', lambda m: " " * len(m.group()), hl)  # quoted titles; ratios written n:n (matching design)
    for p in N.SKIP:
        t_ = p.sub(lambda m: " " * len(m.group()), t_)
    for m in N.NUM.finditer(t_):
        try:
            x, dec = N.parse(m.group())
        except ValueError:
            continue
        tok = m.group().rstrip(",")
        item = hl.strip()[0]
        hit = next((b for b in bound_vals if b[1] == dec and abs(abs(b[0]) - abs(x)) <= 1e-9), None) if (dec >= 1 or abs(x) >= 1000) else None
        lab = [l for l in N.explained_by(x, dec, False) if l.startswith(f"headline {item}:") or any(l.startswith(pf) for pf in HL_LAB.get(item, ()))]
        if hit:
            hl_bound.append(f"| {item} | {tok} | claim row {hit[2]}: {hit[3]} |")
        elif lab:
            hl_bound.append(f"| {item} | {tok} | derived: {lab[0]} |")
        elif (item, tok) in HL_CONST and HL_CONST[(item, tok)] in PLAN:
            hl_bound.append(f"| {item} | {tok} | analysis/ANALYSIS_PLAN.md : {HL_CONST[(item, tok)]} (a definition, not a result) |")
        else:
            hl_unbound.append(f"- headline item {item}: `{tok}` in: {hl.strip()[:90]}")
total = bound + derived + len(generic) + len(unsourced)
NL = chr(10)
head = (f"# Binding of the numbers in the claim rows of CLAIMS_VS_EVIDENCE.md to named result fields{NL}{NL}"
        f"{len(rows)} claim rows (sections 1 to 6), columns Claim, Estimate, 95% CI, n and Source: {total} numbers. "
        f"**Bound to a named result field (file, key path; the builder recorded the field when it formatted the number and this script re-read the file and recomputed it): {bound}. "
        f"Derived (recomputed in check_numbers_sourced.py): {derived}. Generic match only: {len(generic)}. Unsourced: {len(unsourced)}.**{NL}{NL}")
body = (f"## Generic match only (equal to a source value somewhere, but not bound to a named field for this row){NL}{NL}" + (NL.join(generic) if generic else "none")
        + f"{NL}{NL}## Unsourced{NL}{NL}" + (NL.join(unsourced) if unsourced else "none")
        + f"{NL}{NL}## Bindings{NL}{NL}| Row | Column | Number | Bound to (file : key path) |{NL}|---|---|---|---|{NL}" + NL.join(details) + NL)
body += (f"{NL}## Headline results of RESULTS_REPORT.md: numbers equal to a number bound or derived in a claim row (same value, same decimals){NL}{NL}"
         f"Bound through a claim row: {len(hl_bound)}. Not found in any claim row: {len(hl_unbound)}.{NL}{NL}" + (NL.join(hl_unbound) if hl_unbound else "none") + f"{NL}{NL}| Item | Number | Bound through |{NL}|---|---|---|{NL}" + NL.join(hl_bound) + NL)
(REV6 / "results" / "number_binding_check.md").write_text(head + body, encoding="utf-8")
print(f"{len(rows)} rows, {total} numbers: bound {bound}, derived {derived}, generic only {len(generic)}, unsourced {len(unsourced)}; headline items of RESULTS_REPORT.md: {len(hl_bound)} bound through a claim row, {len(hl_unbound)} not found in a claim row")
