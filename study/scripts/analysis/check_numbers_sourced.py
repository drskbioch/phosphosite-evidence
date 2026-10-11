"""Check that every number in RESULTS_REPORT.md and CLAIMS_VS_EVIDENCE.md has a source in a result file.

Source values: every number in the JSON files under study/analysis/results and study/earlier_results, and every numeric cell of the CSV files in study/analysis
(manifest_status.csv, grid_manifest.csv). A number in a document matches a source value s when, at the number of decimals written, it equals s, abs(s),
100*s (percent notation) or 100*abs(s). Identifiers (C1.6, G21-0491, M1s, K3, D7, V50, ...), dates, section references, file names in code spans, commit
hashes and PMIDs that are written as identifiers are not numbers and are skipped. Numbers that match nothing are listed with their context in
results/number_source_check.md and counted; they are derived values (sums, differences, shares, ratios of source values) or text numbers that need a source.
For those, a second pass tests derivation from two source values (sum, difference, ratio, share, difference of shares) among the values found within
the same line of the document; a number that still matches nothing is reported as unsourced.
Usage: python check_numbers_sourced.py
"""
import bisect
import csv
import glob
import itertools
import json
import re
from pathlib import Path

REV6 = Path(__file__).resolve().parents[2] / "analysis"
CUR = REV6.parent
DOCS = [REV6 / "RESULTS_REPORT.md", REV6 / "CLAIMS_VS_EVIDENCE.md"]
NUM = re.compile(r"(?<![\w.])[-+−]?\d[\d,]*(?:\.\d+)?(?:[eE][-+]?\d+)?")
SKIP = [
    re.compile(r"`[^`]*`"),  # code spans (file names, commands)
    re.compile(r"\b[A-Za-z]+[-_]?\d[\w.-]*"),  # identifiers: C1.6, G21-0491, M1s, K3, D7, V50, S2, P-abund1 ...
    re.compile(r"\b\d{4}-\d{2}-\d{2}\b"),  # dates
    re.compile(r"\b[Ss]ections?\s+\d+[a-z]?(?:\s*(?:and|to|,)\s*\d+[a-z]?)*"),
    re.compile(r"\b(?:Rev|rev|Task|Tasks|limitation|limitations|decision|Decision|item|items|step|batch|Table|Figure)\s+\d+[a-z]?(?:\s*(?:and|to|,)\s*\d+[a-z]?)*"),
    re.compile(r"\b(?:Branch)\s+[AB]\b"),
    re.compile(r"\b[0-9a-f]{7,40}\b"),  # commit hashes
    re.compile(r"^\s*#+\s*\d+[a-z]?\.?", re.M),  # headings
    re.compile(r"^\s*\d+\.\s", re.M),  # numbered list items
    re.compile(r"\[\d+\]"),
]
TOK = re.compile(r"\d")


def values(o, out):
    if isinstance(o, bool):
        return
    if isinstance(o, (int, float)):
        out.append(float(o))
    elif isinstance(o, str):
        for m in NUM.finditer(o):
            try:
                out.append(float(m.group().replace(",", "").replace("−", "-")))
            except ValueError:
                pass
    elif isinstance(o, dict):
        for k, v in o.items():
            values(k, out)
            values(v, out)
    elif isinstance(o, (list, tuple)):
        for v in o:
            values(v, out)


def load_sources():
    vals = []
    files = glob.glob(str(REV6 / "results" / "**" / "*.json"), recursive=True) + glob.glob(str(CUR / "earlier_results" / "*.json"))
    files = [f for f in files if Path(f).name not in ("derived_figure_values.json", "figure_counts.json", "figure_bindings.json")]  # files written for the figures; they are not sources of the claims file
    for f in files:
        try:
            values(json.load(open(f, encoding="utf-8")), vals)
        except Exception:
            pass
    for f in (REV6 / "results" / "manifest_status.csv", REV6 / "grid_manifest.csv"):
        for r in csv.reader(open(f, encoding="utf-8")):
            for c in r:
                values(c, vals)
    return files, vals


files, vals = load_sources()
vals = [v for v in vals if v == v and abs(v) != float("inf")]  # NaN or infinite entries would break the sorted search
absvals = sorted(set(abs(v) for v in vals))
pctvals = sorted(set(100 * v for v in absvals))


def has(arr, x, tol):
    i = bisect.bisect_left(arr, x - tol)
    return i < len(arr) and arr[i] <= x + tol


def matches(x, dec, pct):
    tol = 0.5 * 10 ** (-dec) + 1e-9
    x = abs(x)
    if has(absvals, x, tol):
        return True
    if has(pctvals, x, tol):
        return True
    if pct and has(absvals, x / 100, tol / 100):
        return True
    return False


# ------------------------------------------------------------------ derivations (recomputed from the result files, independent of the report builders)
def _j(name):
    return json.load(open(REV6 / "results" / name, encoding="utf-8"))


def derivations():
    """Labelled values that the documents report but that are not stored in any result file: each is recomputed here from the result files."""
    import math
    import statistics

    import pandas as pd

    out = []  # (label, value)
    man = pd.read_csv(REV6 / "grid_manifest.csv").fillna("")
    cem, g5, am = _j("cem_results.json"), _j("g5_abundance.json"), _j("amend_results.json")
    src_of = {}
    for _, r in man[man["block"] == "G17"].iterrows():
        src_of[re.search(r"same records as (\S+?);", r["note"]).group(1)] = r["contrast_id"]

    def ids(block, sig, pop=None):
        m = man[(man["block"] == block) & (man["signal"] == sig)]
        if pop:
            m = m[m["population"] == pop]
        return {d: i for d, i in zip(m["design"], m["contrast_id"])}

    # (1) width ratio gene-cluster interval / accession-cluster interval, per contrast row; median, min, max over rows
    ratios = []
    for blk, sigs in (("G1", ["S1", "S2", "S3", "S4", "S5"]), ("G2", ["S2", "S3", "S4", "S5"]), ("G3", ["S1"])):
        for s in sigs:
            d = ids(blk, s)
            for dz in ("D0", "D1"):
                a, g = cem[d[dz]]["ci"], am[src_of[d[dz]]]["ci"]
                ratios.append((f"width ratio {blk} {s} {dz}", (g[1] - g[0]) / (a[1] - a[0])))
    for s in ("S1", "S2", "S3", "S4"):
        d = ids("G5", s, "P-abund (human)")
        for dz in ("D0", "D1"):
            a, g = g5[d[dz]]["ci"], am[src_of[d[dz]]]["ci"]
            ratios.append((f"width ratio G5 {s} {dz}", (g[1] - g[0]) / (a[1] - a[0])))
    out += ratios
    w = [v for _, v in ratios]
    out += [("width ratio median", statistics.median(w)), ("width ratio min", min(w)), ("width ratio max", max(w))]

    # (2) tier-definition spread (G19): per population and signal
    tiers = {}
    for cid, v in am.items():
        if cid.startswith("G19"):
            key = (v["signal"], v["population"], v["tier_rule"].split(" (")[0])
            tiers.setdefault(key, {})[v["design"]] = v
    null = lambda s: 1.0 if s == "S1" else 0.0  # noqa: E731
    excl = lambda c, s: c[0] > null(s) or c[1] < null(s)  # noqa: E731
    elig = lambda c, s: excl(c, s) and min(abs(c[0] - null(s)), abs(c[1] - null(s))) >= 0.01  # noqa: E731
    L_adj_all, L_def_all = [], []
    for pop in ("P-human", "P-all"):
        for s in ("S1", "S2", "S3", "S4"):
            ks = [k for k in tiers if k[0] == s and k[1] == pop]
            if not ks:
                continue
            defs = [(tiers[k]["D0"], tiers[k]["D1"]) for k in ks]
            adj = [d1["estimate"] for _, d1 in defs]
            out.append((f"tier spread, all definitions, {pop} {s}", max(adj) / min(adj)))
            el = [d1["estimate"] for _, d1 in defs if elig(d1["ci"], s)]
            if len(el) >= 2:
                out.append((f"tier spread, eligible definitions, {pop} {s}", max(el) / min(el)))
            okc = [(d0, d1) for d0, d1 in defs if excl(d0["ci"], s) and excl(d1["ci"], s) and d0["estimate"] > 0 and d1["estimate"] > 0]
            oka = [d1["estimate"] for _, d1 in defs if excl(d1["ci"], s) and d1["estimate"] > 0]
            La = max(abs(math.log(d1["estimate"] / d0["estimate"])) for d0, d1 in okc) if okc else None
            Ld = math.log(max(oka) / min(oka)) if len(oka) >= 2 else None
            if La is not None:
                out.append((f"largest |ln(adjusted/crude)| within a definition, {pop} {s}", La))
                out.append((f"plain ratio exp(largest within-definition change), {pop} {s}", math.exp(La)))
            if Ld is not None:
                out.append((f"ln spread of adjusted estimates across definitions, {pop} {s}", Ld))
                out.append((f"plain ratio exp(ln spread across definitions), {pop} {s}", math.exp(Ld)))
            if La is not None and Ld is not None:
                L_adj_all.append(La)
                L_def_all.append(Ld)
    out += [("min of within-definition change", min(L_adj_all)), ("max of within-definition change", max(L_adj_all)),
            ("min of ln spread across definitions", min(L_def_all)), ("max of ln spread across definitions", max(L_def_all))]
    out += [("plain ratio exp(min of within-definition change)", math.exp(min(L_adj_all))), ("plain ratio exp(max of within-definition change)", math.exp(max(L_adj_all))),
            ("plain ratio exp(min of ln spread across definitions)", math.exp(min(L_def_all))), ("plain ratio exp(max of ln spread across definitions)", math.exp(max(L_def_all)))]
    s3 = [tiers[k]["D1"]["estimate"] for k in tiers if k[0] == "S3" and k[1] == "P-all"]
    out.append(("disorder P-all spread, largest / smallest adjusted estimate (fold)", max(s3) / min(s3)))

    # (3) matcher results (G18): share of BRONZE records that are paired, where only part of GOLD is matched
    shares = []
    for f in glob.glob(str(REV6 / "results" / "d2b" / "*.json")):
        v = json.load(open(f, encoding="utf-8"))
        if v["share_gold_matched"] < 0.999:
            shares.append(v["n_pairs"] / v["n_bronze"])
    out += [("pairs as share of BRONZE records, minimum (percent)", 100 * min(shares)), ("pairs as share of BRONZE records, maximum (percent)", 100 * max(shares))]

    # (4) PeptideAtlas build 606: peptide instances in total = sum of the counts by number of protein mappings
    out.append(("peptide instances in build 606", float(sum(_j("peptideatlas_606_feasibility.json")["peptide_instance_n_protein_mappings"].values()))))
    # (5) counts and shares written in the claim rows
    dc = json.load(open(CUR / "earlier_results" / "dataset_counts.json", encoding="utf-8"))
    out.append(("sum of the mammalian tier counts", float(sum(dc["tier_counts_mammalian"].values()))))
    nh = dc["human_n_sites"]
    out += [("human GOLD share under the raw PMID rule (percent)", 100 * dc["human_tier_raw_pmid_rule"]["GOLD"] / nh),
            ("human GOLD share under the study-cluster rule (percent)", 100 * dc["human_tier_study_cluster_rule"]["GOLD"] / nh)]
    g23 = _j("g23_tiers_summary.json")
    out.append(("GOLD records that fall to SILVER or BRONZE when non-primary PMIDs are not counted", float(g23["tier_counts_old"]["GOLD"] - g23["tier_counts_new"]["GOLD"])))
    out += [("signal-by-population combinations where the log spread across definitions exceeds the within-definition change", float(sum(1 for la, ld in zip(L_adj_all, L_def_all) if ld > la))),
            ("signal-by-population combinations where both can be computed", float(len(L_adj_all)))]
    g14 = _j("g14_cem.json")
    out.append(("CEM rows with a protein-cluster to Rev 5 interval width ratio (G14)", float(sum(1 for v in g14.values() if isinstance(v, dict) and v.get("ci_width_ratio_protein_cluster_over_rev5") is not None and v.get("status") == "run"))))
    b4 = _j("b4_twoway.json")
    g1d1 = [ids("G1", s)["D1"] for s in ("S1", "S2", "S3", "S4")]
    out += [("decision-bearing G1 contrasts whose two-way interval excludes the null", float(sum(1 for i in g1d1 if b4[i]["survives_two_way"]))),
            ("decision-bearing G1 contrasts with a two-way interval", float(len(g1d1)))]
    g23r = _j("g23_results.json")
    pairs = []
    for cid, v in g23r.items():
        if v["design"] != "D1":
            continue
        src = man[man["contrast_id"] == cid]["note"].iloc[0].split("counterpart of ")[1].split(";")[0]
        p = am[src_of[src]]
        sg = v["signal"]
        pairs.append((excl(v["ci"], sg), (v["estimate"] > null(sg)) == (p["estimate"] > null(sg)), bool(v.get("largely_explained"))))
    out += [("G23 D1 contrasts whose interval excludes the null", float(sum(a for a, _, _ in pairs))), ("G23 D1 contrasts that keep the sign of the primary", float(sum(b for _, b, _ in pairs))),
            ("G23 D1 contrasts largely explained", float(sum(c for _, _, c in pairs))), ("G23 D1 contrasts", float(len(pairs)))]
    g11 = _j("g11_predictors.json")
    ok = 0
    for tool in ("DeepPhos", "MusiteDeep", "PhosF3C"):
        k = [x for x in g11 if x.startswith(tool + "|planned")][0]
        v = g11[k].get("unrestricted")
        if isinstance(v, dict) and "gap" in v and v["gap"] >= 0.05 and (v["gene_ci95"][0] > 0 or v["gene_ci95"][1] < 0):
            ok += 1
    out += [("predictor tools meeting the pre-specified criterion (gap >= 0.05, interval excluding zero)", float(ok)), ("predictor tools audited", 3.0)]

    # (5b) G21: conditional-logit resample fits that did not converge, summed over all contrasts
    g21 = _j("g21_results.json")
    out.append(("G21 conditional-logit resample fits that did not converge, all contrasts", float(sum(v.get("n_boot_unfit", 0) for v in g21.values()))))

    # (5c) G21: share of failed bootstrap fits per contrast (percent), quoted in the C5 rows when above 5%
    for cid, v in g21.items():
        out.append((f"G21 failed-fit share, {cid} (percent)", 100 * v.get("n_boot_unfit", 0) / v["n_boot"]))

    # (6) counts and extremes quoted in the headline results of RESULTS_REPORT.md (labels start with "headline <item>:")
    mid = man.set_index("contrast_id")
    dec_d1 = []
    for blk, sigs in (("G1", ["S1", "S2", "S3", "S4", "S5"]), ("G2", ["S2", "S3", "S4", "S5"]), ("G3", ["S1"])):
        for s in sigs:
            i = ids(blk, s)["D1"]
            if mid.loc[i, "decision_bearing"] == "yes":
                dec_d1.append(i)
    nd = len(dec_d1)
    out += [("headline 1: decision-bearing D1 contrasts in G1 and G2 (and G3)", float(nd)),
            ("headline 1: of them surviving with the gene-cluster interval", float(sum(1 for i in dec_d1 if am[src_of[i]]["survives"]))),
            ("headline 1: of them surviving with the accession-cluster interval", float(sum(1 for i in dec_d1 if cem[i]["survives"]))),
            ("headline 1: of them largely explained", float(sum(1 for i in dec_d1 if am[src_of[i]].get("largely_explained")))),
            ("headline 1: largest attenuation, gene-cluster (percent)", 100 * max(am[src_of[i]]["attenuation"] for i in dec_d1)),
            ("headline 1: decision-bearing G1 contrasts with a two-way interval", float(sum(1 for i in dec_d1 if i in b4))),
            ("headline 1: of them excluding the null under the two-way interval", float(sum(1 for i in dec_d1 if i in b4 and b4[i]["survives_two_way"])))]
    d2b = [json.load(open(f, encoding="utf-8")) for f in glob.glob(str(REV6 / "results" / "d2b" / "*.json"))]
    dd = [v for v in d2b if v["decision_bearing"] == "yes"]
    out += [("headline 4: decision-bearing D2b contrasts", float(len(dd))), ("headline 4: of them surviving", float(sum(1 for v in dd if v["survives"]))),
            ("headline 4: largest attenuation (percent)", 100 * max(v["attenuation"] for v in dd)), ("headline 4: D2b point estimates outside their own interval", float(sum(1 for v in d2b if v["point_outside_own_interval"])))]
    pb = _j("g20_diagnostics.json")["prematching_balance"]
    dbk = [k for k in pb if mid.loc[k, "decision_bearing"] == "yes"]
    pre = [pb[k]["smd_pre"] for k in dbk]
    aft = [pb[k]["smd_after"] for k in dbk]
    out += [("headline 3: decision-bearing contrasts with a before/after balance", float(len(dbk))), ("headline 3: of them with |SMD| of 0.1 or more before matching", float(sum(1 for v in pre if abs(v) >= 0.1))),
            ("headline 3: largest |SMD| before matching", max(abs(v) for v in pre)), ("headline 3: largest |SMD| after matching", max(abs(v) for v in aft))]
    for i in dec_d1:
        if abs(pb[i]["smd_pre"]) >= 0.1:
            out.append((f"headline 3: SMD before matching, {i}", pb[i]["smd_pre"]))
    pr = [v["pearson_r"] for v in _j("g20_diagnostics.json")["index_outcome_association"].values()]
    out += [("headline 3: smallest within-tier correlation of the index with an outcome", min(pr)), ("headline 3: largest within-tier correlation of the index with an outcome", max(pr))]
    prop_ids = []
    for f in glob.glob(str(REV6 / "results" / "prop" / "*.json")):
        v = json.load(open(f, encoding="utf-8"))
        if "+" not in v["contrast_id"] and mid.loc[v["contrast_id"], "block"] in ("G1", "G2", "G3", "G15"):
            prop_ids.append(v)
    out += [("headline 6: planned D2 contrasts", float(len(prop_ids))), ("headline 6: D2 contrasts that could not be fitted", float(sum(1 for v in prop_ids if v["status"] != "run")))]
    labs = [v["gene_label"] for k in g11 if "planned" in k for v in g11[k].values() if isinstance(v, dict) and "gene_label" in v]
    out += [("headline 6: evaluated predictor settings", float(len(labs))), ("headline 6: settings with a gap equivalent to zero", float(sum(1 for lab in labs if lab.startswith("gap absent")))),
            ("headline 6: settings not detectable", float(sum(1 for lab in labs if lab.startswith("not detectable")))), ("headline 6: settings showing a gap", float(sum(1 for lab in labs if lab.startswith("gap present"))))]
    return out


DER = derivations()
DER_VALS = sorted(abs(v) for _, v in DER)
DER_LAB = {}
for lab, v in DER:
    DER_LAB.setdefault(round(abs(v), 6), []).append(lab)


def explained_by(x, dec, pct):
    """Labels of derivations whose value rounds to the number written (percent notation: the value is already in percent)."""
    tol = 0.5 * 10 ** (-dec) + 1e-9
    return sorted({lab for lab, v in DER if abs(abs(v) - abs(x)) <= tol})


def parse(tok):
    t = tok.replace(",", "").replace("−", "-")
    dec = len(t.split(".")[1].split("e")[0].split("E")[0]) if "." in t else 0
    return float(t), dec


def derived(x, dec, line_vals):
    tol = 0.5 * 10 ** (-dec) + 1e-9
    for a, b in itertools.permutations(line_vals, 2):
        for c in (a + b, a - b, abs(a - b), (a / b if b else None), (100 * a / b if b else None), (100 * (a - b) / b if b else None), a * b):
            if c is not None and abs(abs(c) - abs(x)) <= tol:
                return True
    return False


def main():
    report, explained, total, unmatched = [], [], 0, 0
    for doc in DOCS:
        lines = doc.read_text(encoding="utf-8").split("\n")
        for ln, line in enumerate(lines, 1):
            s = line
            for p in SKIP:
                s = p.sub(lambda m: " " * len(m.group()), s)
            toks = [(m.group(), m.end()) for m in NUM.finditer(s)]
            if not toks:
                continue
            lv = []
            for t, _ in toks:
                try:
                    lv.append(parse(t)[0])
                except ValueError:
                    pass
            for t, end in toks:
                try:
                    x, dec = parse(t)
                except ValueError:
                    continue
                total += 1
                pct = s[end:end + 1] == "%"
                if matches(x, dec, pct):
                    continue
                # second pass: line-level values that are themselves sourced
                src_line = [v for v in lv if v != x and matches(v, 0 if float(v).is_integer() else 6, False)]
                if derived(x, dec, src_line):
                    continue
                lab = explained_by(x, dec, pct)
                if lab:
                    explained.append(f"- {doc.name}:{ln}: `{t}` = " + " | ".join(lab))
                    continue
                unmatched += 1
                report.append(f"- {doc.name}:{ln}: `{t}` in: {line.strip()[:220]}")
    out = REV6 / "results" / "number_source_check.md"
    head = (f"# Number-source check of RESULTS_REPORT.md and CLAIMS_VS_EVIDENCE.md\n\nSource values read from {len(files)} JSON files and the two manifest CSV files ({len(vals):,} values). Numbers examined: {total:,}. "
            f"Equal to a source value: {total - len(explained) - unmatched:,}. Equal to a derivation recomputed in this script from the result files: {len(explained)}. **Neither: {unmatched}.**\n\n")
    body = "## Neither a source value nor a derivation\n\n" + ("\n".join(report) if report else "none") + "\n\n## Explained by a derivation\n\n" + "\n".join(explained) + "\n\n## All derivations computed (label: value)\n\n" + "\n".join(f"- {lab}: {v:.6g}" for lab, v in DER) + "\n"
    out.write_text(head + body, encoding="utf-8")
    print(f"{total} numbers examined, {len(explained)} explained by derivation, {unmatched} unexplained; list in {out}")


if __name__ == "__main__":
    main()
