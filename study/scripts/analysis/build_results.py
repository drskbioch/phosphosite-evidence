"""Assemble results/manifest_status.csv (all manifest rows) and study/analysis/RESULTS_REPORT.md from the Rev 6 result files.
Every number in the report is read from a result file or computed from one here; nothing is typed in.
Sections that depend on a file the project lead has not yet placed (PhosphoSitePlus) are reported as not run, with the reason.
"""
import glob
import json
import math
import re

import numpy as np
import pandas as pd

from common import N_BOOT, REV6, RES

J = lambda p: json.load(open(RES / p))  # noqa: E731
man = pd.read_csv(REV6 / "grid_manifest.csv").fillna("")
manid = man.set_index("contrast_id")
cem = J("cem_results.json")
g5 = J("g5_abundance.json")
g14 = J("g14_cem.json")
am = J("amend_results.json")
g20 = J("g20_diagnostics.json")
b4 = J("b4_twoway.json")
g11 = J("g11_predictors.json")
gu = J("gene_unit_summary.json")
g16c = J("g16_coverage.json")
d345 = J("d3_d4_d5_index_checks.json")
pa = J("peptideatlas_606_inspection.json")
val = J("validate_engine.json")
g23_ = J("g23_results.json") if (RES / "g23_results.json").exists() else {}
g22 = J("g22_results.json") if (RES / "g22_results.json").exists() else {}
g21_ = J("g21_results.json") if (RES / "g21_results.json").exists() else {}
g9 = J("g9_measures.json") if (RES / "g9_measures.json").exists() else {}
g9i = J("g9_intersections.json") if (RES / "g9_intersections.json").exists() else {}
prop, d2b, fr = {}, {}, {}
for f_ in glob.glob(str(RES / "prop" / "*.json")):
    r_ = json.load(open(f_))
    prop[r_["contrast_id"]] = r_
for f_ in glob.glob(str(RES / "d2b" / "*.json")):
    r_ = json.load(open(f_))
    d2b[r_["contrast_id"]] = r_
for f_ in glob.glob(str(RES / "f_runs" / "*.json")):
    fr.update(json.load(open(f_)))
PSP_PRESENT = (REV6 / "data_private" / "Phosphorylation_site_dataset.gz").exists()
PSP_MSG = "not run: source file not available (PhosphoSitePlus v7 site provides no bulk download; access requested from the provider)"
SL = {"S1": "Domain overlap", "S2": "Conservation", "S3": "Disorder", "S4": "Motif, sequence-based", "S5": "Motif, database lookup"}
NULL = lambda s: 1.0 if s == "S1" else 0.0  # noqa: E731


def nan(v):
    return v is None or (isinstance(v, float) and math.isnan(v))


def f(v, s, nd=None):
    if nan(v):
        return "n/a"
    nd = nd or (3 if s == "S1" else 4)
    return f"{v:.{nd}f}"


def ci(c, s, nd=None):
    return f"{f(c[0], s, nd)} to {f(c[1], s, nd)}" if c and not nan(c[0]) else "n/a"


def est_ci(e, c, s, nd=None):
    return f"{f(e, s, nd)} ({ci(c, s, nd)})"


def pct(v, nd=1):
    return "n/a" if nan(v) else f"{100 * v:.{nd}f}%"


def excl(c, s):
    return None if (c is None or nan(c[0])) else bool(c[0] > NULL(s) or c[1] < NULL(s))


def yn(b):
    return "n/a" if b is None else ("yes" if b else "no")


def n_(v):
    return f"{int(round(v)):,}"


def side(e, s):
    return "above null" if e > NULL(s) else "below null"


def md(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------ links between blocks
src_of = {}  # source contrast id -> gene-cluster (G17) id
for cid, r in manid.iterrows():
    if r["block"] == "G17":
        src_of[re.search(r"same records as (\S+?);", r["note"]).group(1)] = cid


def d1_d0(block, sig):
    """(D0 id, D1 id) of a G1/G2/G3 signal."""
    m = man[(man["block"] == block) & (man["signal"] == sig)]
    return m[m["design"] == "D0"]["contrast_id"].iloc[0], m[m["design"] == "D1"]["contrast_id"].iloc[0]


# ------------------------------------------------------------------ manifest status
def status_row(r):
    cid, blk = r["contrast_id"], r["block"]
    st, est, c_, att, note = "pending", "", "", "", ""
    if cid in cem:
        x = cem[cid]
        st, est, c_, att = x["status"], x["estimate"], x["ci"], x.get("attenuation", "")
    elif cid in g5:
        x = g5[cid]
        st, est, c_, att, note = x["status"], x.get("estimate", ""), x.get("ci", ""), x.get("attenuation", ""), x.get("reason", "")
    elif cid in g14:
        x = g14[cid]
        st, est, c_ = x["status"], x.get("estimate", ""), x.get("ci_protein_cluster", "")
    elif cid in prop:
        x = prop[cid]
        st, est, c_, att, note = x["status"], x.get("estimate", ""), x.get("ci", ""), x.get("attenuation", ""), x.get("reason", "")
    elif cid in d2b:
        x = d2b[cid]
        st, est, c_, att = x["status"], x["estimate"], x["ci"], x.get("attenuation", "")
    elif cid in am:
        x = am[cid]
        st, est, c_, att, note = x["status"], x.get("estimate", ""), x.get("ci", ""), x.get("attenuation", ""), x.get("reason", "")
    elif cid in g23_:
        x = g23_[cid]
        st, est, c_, att = x["status"], x["estimate"], x["ci"], x.get("attenuation", "")
    elif cid in g22:
        x = g22[cid]
        st, est, c_, att = x["status"], x["estimate"], x["ci"], x.get("attenuation", "")
    elif cid in g21_:
        x = g21_[cid]
        st, est, c_, note = x["status"], x["estimate"], x["ci"], x["estimator"]
    elif cid in fr:
        x = fr[cid]["gene"]
        st, est, c_, att, note = x["status"], x.get("estimate", ""), x.get("ci", ""), x.get("attenuation", ""), x.get("reason", "")
    elif cid in g9:
        x = g9[cid]
        st, note = x["status"], x.get("reason", "")
        y = x.get("gene")
        if y:
            est = y.get("estimate", y.get("adjusted", ""))
            c_ = y.get("ci", y.get("adjusted_ci", ""))
            att = y.get("attenuation", "")
    elif blk == "G9" and r["variation"].startswith("PeptideAtlas"):
        st, note = "not feasible from the available build", "no site-localisation score or site-level evidence category in build 606; see RESULTS_REPORT.md section 9 and results/peptideatlas_606_feasibility.json"
    elif blk == "G9" and r["variation"].startswith("PSP"):
        st, note = (PSP_MSG, "Task H, PhosphoSitePlus measures") if not PSP_PRESENT else ("pending (task H)", "")
    elif blk == "G10":
        st, note = (PSP_MSG, "Task I, within-gene LT-supported against MS-only") if not PSP_PRESENT else ("pending (task I)", "")
    elif blk == "G11":
        tool = r["population"].split(" ")[0]
        key = [k for k in g11 if k.startswith(tool + "|planned")][0]
        x = g11[key][r["variation"]]
        st = "run" if x.get("status") == "evaluated" else "not run"
        est, c_, note = x.get("gap", ""), x.get("gene_ci95", ""), x.get("status", "") + "; " + x.get("gene_label", "")
    elif blk == "G20":
        k = int(cid.split("-")[1])
        if cid in g20["index_outcome_association"]:
            x = g20["index_outcome_association"][cid]
            st, est, c_ = "run", x["pearson_r"], x["ci"]
        elif k == 468:
            st = "run"
        elif k == 487:
            x = g20["index_vs_sites_per_accession"]["human|sites_per_accession"]
            st, est, c_ = "run", x["spearman"], x["ci"]
        elif k == 488:
            st, note = (PSP_MSG, "G20-0488, index against LT_LIT-supported sites per gene") if not PSP_PRESENT else ("pending (task H)", "needs PhosphoSitePlus LT_LIT counts")
        else:
            x = g20["domain_overlap_by_annotation_coverage"]["at least one domain-annotated site" if k == 489 else "no domain-annotated site"]
            st = "run"
            est, c_ = x["crude_or"], x["ci"]
            note = x.get("note", "")
    elif blk in ("G1", "G2", "G3", "G4", "G5", "G12", "G13", "G14", "G15"):
        st = "not run"
    return {**r.to_dict(), "result_status": st, "estimate": "" if nan(est) else est, "ci_low": c_[0] if c_ and not nan(c_[0]) else "", "ci_high": c_[1] if c_ and not nan(c_[1]) else "",
            "attenuation": "" if nan(att) else att, "result_note": note}


ms = pd.DataFrame([status_row(r) for _, r in man.iterrows()])
ms.to_csv(RES / "manifest_status.csv", index=False)
cnt = ms.groupby(["block", "result_status"]).size().unstack(fill_value=0)
cnt = cnt.loc[sorted(cnt.index, key=lambda b: int(b[1:]))]

L = []
A = L.append

# ------------------------------------------------------------------ data used in several sections
# main contrasts (G1, G2, G3) with accession-cluster (cem), gene-cluster (G17) and two-way (B4) results
main = []
for blk, sigs in (("G1", ["S1", "S2", "S3", "S4", "S5"]), ("G2", ["S2", "S3", "S4", "S5"]), ("G3", ["S1"])):
    for s in sigs:
        d0, d1 = d1_d0(blk, s)
        main.append({"block": blk, "s": s, "d0": d0, "d1": d1, "acc": cem[d1], "acc0": cem[d0], "gene": am[src_of[d1]], "gene0": am[src_of[d0]], "tw": b4.get(d1),
                     "db": manid.loc[d1, "decision_bearing"] == "yes", "pop": manid.loc[d1, "population"]})
abund = []
for s in ("S1", "S2", "S3", "S4"):
    m = man[(man["block"] == "G5") & (man["signal"] == s) & (man["population"] == "P-abund (human)")]
    ids = dict(zip(m["design"], m["contrast_id"]))
    abund.append({"s": s, "d0": ids["D0"], "d1": ids["D1"], "acc": g5[ids["D1"]], "acc0": g5[ids["D0"]], "gene": am[src_of[ids["D1"]]], "gene0": am[src_of[ids["D0"]]], "tw": None, "db": True,
                  "pop": "P-abund (human)", "block": "G5"})


def adjusted_estimates(s):
    """Adjusted (D1, D4, D5) point estimates of a signal in every population where it is defined: G1, G2, G3 and P-abund (human, mouse, rat)."""
    out = []
    for x in main:
        if x["s"] == s:
            out.append((f"{x['block']} {x['pop']}", "D1", x["acc"]["estimate"]))
    m = man[(man["block"] == "G5") & (man["signal"] == s)]
    for _, r in m.iterrows():
        if r["design"] in ("D1", "D5", "D4") and r["contrast_id"] in g5 and g5[r["contrast_id"]]["status"] == "run":
            out.append((r["population"], r["design"], g5[r["contrast_id"]]["estimate"]))
    return out


consist = {}
for s in SL:
    e = adjusted_estimates(s)
    sd = {side(v, s) for _, _, v in e}
    consist[s] = {"n": len(e), "consistent": len(sd) == 1, "sides": sd, "rows": e}


def label(s):
    return "direction-consistent" if consist[s]["consistent"] else "population-dependent"


# ------------------------------------------------------------------ report
n_dec = [x for x in main if x["db"]]
surv_acc = sum(1 for x in n_dec if x["acc"]["survives"])
surv_gene = sum(1 for x in n_dec if x["gene"]["survives"])
tw_dec = [x for x in n_dec if x["tw"]]
tw_surv = sum(1 for x in tw_dec if x["tw"]["survives_two_way"])
largest_att = max(x["gene"]["attenuation"] for x in n_dec)
d2b_dec = [v for v in d2b.values() if v["decision_bearing"] == "yes"]
d2b_surv = sum(1 for v in d2b_dec if v["survives"])
d2b_flag = sum(1 for v in d2b.values() if v["point_outside_own_interval"])
d2b_att = max(v["attenuation"] for v in d2b_dec)

d2_ids = [k for k in prop if "+" not in k and manid.loc[k, "block"] in ("G1", "G2", "G3", "G15")]
d2_planned = len(d2_ids)
d2_nofit = sum(1 for k in d2_ids if prop[k]["status"] != "run")
labs = [v["gene_label"] for k in g11 if "planned" in k for v in g11[k].values() if isinstance(v, dict) and "gene_label" in v]
pred_summary = (f"of {len(labs)} evaluated predictor settings {sum(1 for l_ in labs if l_.startswith('gap absent'))} show a gap equivalent to zero, {sum(1 for l_ in labs if l_.startswith('not detectable'))} are not detectable "
                f"and {sum(1 for l_ in labs if l_.startswith('gap present'))} show a gap")

n_h = int(((ms["block"] == "G9") & (ms["result_status"].isin(["pending (task H)", PSP_MSG]))).sum())
n_i = int(((ms["block"] == "G10") & (ms["result_status"].isin(["pending (task I)", PSP_MSG]))).sum())
A("# Results, Rev 6 reanalysis\n")
if PSP_PRESENT:
    status_line = "PhosphoSitePlus is present in `data_private/`."
else:
    status_line = (f"**Final under Branch B (section 14).** Task H (PhosphoSitePlus measures, block G9, {n_h} rows), Task I (block G10, {n_i} rows) and G20-0488 are **{PSP_MSG}**. "
                   "They are listed in Future work (section 15), to be run if the provider supplies the file (possibly during revision); `h_measures.py` is written and has run for the UniProt measure, `cond_logit.py` is validated on real data (G21 check), and the Task I runner is not written. Everything else in the frozen plan has been run.")
A("Branch `rev6-reanalysis`. Written by `scripts/analysis/build_results.py` from the files in `study/analysis/results/`; every number below is read from a result file or computed from one in that script. No manuscript file, reference or legacy file was touched. "
  "The plan is `ANALYSIS_PLAN.md` (deviation log in its section 15, including the amendments of 2026-10-01). `INTERIM_REPORT.md` still holds the interim results (G1 to G5, G12 to G15, index variants, abundance); this report adds the amended and later blocks and re-states the main table with the gene-cluster interval. `results/manifest_status.csv` lists all "
  f"{len(man)} manifest rows with their status.\n")
A(status_line + "\n")

A("## Headline results\n")
HEAD_IDX = len(L)
HL = []
ix = g20["index_outcome_association"]
rmin = min(v["pearson_r"] for v in ix.values())
rmax = max(v["pearson_r"] for v in ix.values())
le = sum(1 for x in n_dec if x["gene"].get("largely_explained"))
d1ids = [x["d1"] for x in n_dec]
smd_dec = [v["smd_pre"] for k, v in g20["prematching_balance"].items() if k in d1ids]
HL.append(f"1. **Matching on the study-effort index (CEM): {surv_gene} of {len(n_dec)} decision-bearing D1 contrasts in G1 and G2 survive with the gene-cluster interval (primary from the amendment) and {surv_acc} of {len(n_dec)} with the accession-cluster interval; {le} are largely explained; the largest attenuation is {pct(largest_att)}.** "
          f"**Under the two-way (protein by study cluster) interval, {tw_surv} of {len(tw_dec)} decision-bearing G1 contrasts exclude the null ({', '.join(SL[x['s']] for x in tw_dec if x['tw']['survives_two_way'])}); "
          f"{', '.join(SL[x['s']] for x in tw_dec if not x['tw']['survives_two_way'])} do not (claim C4.5; approximate interval, sensitivity analysis; section 3).**")
HL.append(f"2. **Domain overlap is {label('S1')}** (it is {'' if consist['S1']['consistent'] else 'not '}on one side of the null in every population): the crude odds ratio is {f(main[0]['gene0']['estimate'], 'S1')} in the full population, {f(abund[0]['gene0']['estimate'], 'S1')} among proteins with an abundance, {f(main[-1]['gene0']['estimate'], 'S1')} among proteins with a domain-annotated site and "
          f"{f(am['G16-0361']['estimate'], 'S1')} among Swiss-Prot reviewed canonical accessions. The other signals: {', '.join(SL[s] + ' ' + label(s) for s in ('S2', 'S3', 'S4', 'S5'))}.")
pbk = g20["prematching_balance"]
dbk = [k for k in pbk if manid.loc[k, "decision_bearing"] == "yes"]
pre = np.array([pbk[k]["smd_pre"] for k in dbk])
aft = np.array([pbk[k]["smd_after"] for k in dbk])
imb = [x for x in n_dec if abs(pbk[x["d1"]]["smd_pre"]) >= 0.1]
HL.append(f"3. **The primary index can carry little confounding.** Before matching GOLD and BRONZE differ in it ({int((np.abs(pre) >= 0.1).sum())} of {len(dbk)} decision-bearing contrasts have |SMD| of 0.1 or more, largest {np.abs(pre).max():.3f}; in G1 and G2: {', '.join(SL[x['s']] + ' ' + x['pop'] + ' ' + f(pbk[x['d1']]['smd_pre'], 'S2', 3) for x in imb) or 'none'}), and matching removes the imbalance (largest SMD after {np.abs(aft).max():.3f}); "
          f"but within a tier the index correlates with an outcome at only {f(rmin, 'S2')} to {f(rmax, 'S2')} (section 4). Adjusting for it therefore moves the estimates little, which is what one expects whether or not site-level study effort confounds the comparison: the design cannot tell the two apart.")
HL.append(f"4. **The matcher the Rev 5 Methods describe (1:1, without replacement, caliper, random order; D2b) gives the same decisions as CEM:** {d2b_surv} of {len(d2b_dec)} decision-bearing D2b contrasts survive, the largest attenuation is {pct(d2b_att)}, {d2b_flag} point estimates lie outside their own interval (section 5).")
HL.append("TIER_HEADLINE")
HL.append(f"6. **Not clean:** the pre-specified propensity design (D2) could not be fitted for {d2_nofit} of {d2_planned} planned contrasts; the pre-specified GO coding measures GOA file coverage and cannot arbitrate between index variants (section 1); the optional UniProt contrast is partly circular for conservation (section 9); {pred_summary} (section 11); PhosphoSitePlus-based measures and Task I are {'run' if PSP_PRESENT else PSP_MSG}.")

# ---------------------------------------------------------------- 1 decisions applied
A("## 1. Decisions applied and fixes\n")
c5 = d345["D5"]
A("- **D5 (index variant), as decided.** The pre-specified rule is applied as written. "
  f"Both pre-specified measures rank {c5['first_by_auroc']} first (Tdark AUROC {f(c5['ranking_auroc_common_set'][c5['first_by_auroc']], 'S2')} against {f(c5['ranking_auroc_common_set']['V50'], 'S2')} for V50; GO Spearman {f(c5['ranking_go_spearman_common_set'][c5['first_by_go']], 'S2')} against {f(c5['ranking_go_spearman_common_set']['V50'], 'S2')}), "
  "so G1 is reported with V50 (primary) and V20 (table below). **The pre-specified GO coding measures GOA file coverage and cannot arbitrate between variants:** "
  f"{n_(d345['info']['human_proteins_without_goa_rows'])} of {n_(d345['info']['human_proteins_with_V0_index'])} human proteins with an index have no row in `goa_human.gaf` and are coded as 0 GO terms, so the correlation is negative for every variant "
  f"({', '.join(f'{k} {f(v, 'S2')}' for k, v in c5['ranking_go_spearman_common_set'].items())}). The GO count restricted to proteins with a GOA row stays as a labelled post-hoc analysis beside it (it ranks {d345['D5_posthoc_GOA_annotated_only']['first_by_go']} first on GO and {d345['D5_posthoc_GOA_annotated_only']['first_by_auroc']} first on the Tdark AUROC; the measures disagree there). No variant changes a decision in G4 (INTERIM_REPORT.md, section 9).\n")
rows = []
for s in ("S1", "S2", "S3", "S4"):
    d1 = main[[x["s"] for x in main].index(s)]["d1"]
    g4 = man[(man["block"] == "G4") & (man["signal"] == s) & (man["index"] == "V20")]
    g4d1 = g4["contrast_id"].iloc[0]
    rows.append([SL[s], est_ci(cem[d1]["estimate"], cem[d1]["ci"], s), est_ci(cem[g4d1]["estimate"], cem[g4d1]["ci"], s), yn(cem[d1]["survives"]), yn(cem[g4d1]["survives"]), pct(cem[d1]["attenuation"]), pct(cem[g4d1]["attenuation"])])
A("G1 D1 with the two indices the rule names (accession-cluster intervals, as in the interim report):\n")
A(md(["Signal", "V50 (G1)", "V20 (G4)", "Survives V50", "Survives V20", "Attenuation V50", "Attenuation V20"], rows))
A(f"- **D2 (propensity matching), as decided.** The {d2_nofit} non-fits of the Rev 5 routine ({', '.join(sorted(k for k in d2_ids if prop[k]['status'] != 'run'))}; `{sorted({prop[k]['reason'] for k in d2_ids if prop[k]['status'] != 'run'})[0]}`) are the result of the pre-specified design (INTERIM_REPORT.md, section 5). The post-hoc S/T/Y-only rows of the interim stage stay beside them; the new matcher is design D2b (section 5).")
bad = []
for fp in sorted(glob.glob(str(REV6 / "*.md"))):
    txt = open(fp, encoding="utf-8").read()
    for ln, line in enumerate(txt.splitlines(), 1):
        if re.search(r"\{(pct|f|ci|DC|rec|val|d2|cov)\b|\{[A-Za-z_]+\[", line):
            bad.append(f"{fp.split(chr(92))[-1]}:{ln}")
A(f"- **Fix, unrendered placeholders.** The interim builder left a template expression unrendered in section 7; fixed and re-rendered (deviation log). A search of every generated report in `study/analysis/` for unrendered template expressions finds {', '.join(bad) if bad else 'none'}. "
  "The same review found one number typed into the interim builder (the count of records with an undefined outcome); it is now read from `results/validate_engine.json`.\n")

# ---------------------------------------------------------------- 2 cluster unit
A("## 2. Cluster unit: gene instead of accession (correction; block G17)\n")
A("**Recorded in the deviation log as a correction made after the accession-level results were seen:** accessions of one gene are not independent, so the primary bootstrap now resamples gene clusters (accessions that share an Entrez gene ID; connected components of the accession-gene graph). All new blocks use the gene cluster. The accession-cluster interval is kept beside it everywhere.\n")
h = gu["human"]
al = gu["all"]
A(f"- **Accessions per gene.** Human: {n_(h['accessions'])} accessions in the site table, {n_(h['accessions_with_gene'])} with a gene mapping in {n_(h['gene_clusters_among_mapped'])} gene clusters (mean {h['accessions_per_gene_cluster']['mean']:.2f}, median {h['accessions_per_gene_cluster']['median']:.0f}, maximum {h['accessions_per_gene_cluster']['max']}; "
  f"{pct(h['accessions_per_gene_cluster']['share_of_clusters_with_more_than_one_accession'])} of clusters have more than one accession). All species: {n_(al['accessions'])} accessions, {n_(al['accessions_with_gene'])} mapped to {n_(al['gene_clusters_among_mapped'])} gene clusters (mean {al['accessions_per_gene_cluster']['mean']:.2f}, maximum {al['accessions_per_gene_cluster']['max']}). "
  f"{n_(h['accessions_without_gene'])} human accessions have no gene mapping and stay single-accession clusters; none of them has an index, so none enters a decision-bearing contrast.")
for key, lab in (("largest_genes_all_records", "all records"), ("largest_genes_human_records", "human records"), ("largest_genes_gold_or_bronze_with_index", "GOLD or BRONZE records with an index")):
    lg = gu[key]
    A(f"- **Records in the 10 largest genes ({lab}).** {n_(lg['records_in_10_largest_genes'])} of {n_(lg['n_records'])} ({pct(lg['share'], 2)}); the 10 largest accessions hold {n_(lg['records_in_10_largest_accessions'])} ({pct(lg['share_10_largest_accessions'], 2)}). "
      f"Largest gene cluster: {lg['top10'][0]['gene_cluster']} with {n_(lg['top10'][0]['records'])} records in {lg['top10'][0]['accessions']} accessions.")
A("")
rows = []
for x in main + abund:
    s = x["s"]
    for tag, a, g in (("D0", x["acc0"], x["gene0"]), ("D1", x["acc"], x["gene"])):
        wa = a["ci"][1] - a["ci"][0]
        wg = g["ci"][1] - g["ci"][0]
        rows.append([x["block"], SL[s] + ("" if x["db"] else " (not decision-bearing)"), x["pop"], tag, f(g["estimate"], s), ci(a["ci"], s), ci(g["ci"], s), f"{wg / wa:.3f}", yn(excl(a["ci"], s)), yn(excl(g["ci"], s))])
A("Accession-cluster against gene-cluster interval (same point estimates; 1,000 resamples each; width ratio = gene / accession):\n")
A(md(["Block", "Signal", "Population", "Design", "Estimate", "Accession-cluster 95% CI", "Gene-cluster 95% CI", "Width ratio", "Excludes null (accession)", "Excludes null (gene)"], rows))
wr = [(g["ci"][1] - g["ci"][0]) / (a["ci"][1] - a["ci"][0]) for x in main + abund for a, g in ((x["acc0"], x["gene0"]), (x["acc"], x["gene"]))]
chg = sum(1 for x in main + abund for a, g in ((x["acc0"], x["gene0"]), (x["acc"], x["gene"])) if excl(a["ci"], x["s"]) != excl(g["ci"], x["s"]))
A(f"The gene-cluster interval is {min(wr):.3f} to {max(wr):.3f} times the width of the accession-cluster interval (median {np.median(wr):.3f}, about {np.median(wr):.2f}; three decimals are shown so that a ratio of about {np.median(wr):.2f} is not rounded without its range); the null is excluded under both intervals in {sum(1 for x in main + abund for a, g in ((x['acc0'], x['gene0']), (x['acc'], x['gene'])) if excl(a['ci'], x['s']) and excl(g['ci'], x['s']))} of {len(wr)} rows and the two intervals disagree on the null in {chg}.\n")

# ---------------------------------------------------------------- 3 main table
A("## 3. Main results: pre-specified rule, direction across populations, two-way interval (G1 to G3, P-abund)\n")
A("Index V50, CEM (D1), set-B tiers. **Survives** is the pre-specified rule (interval excludes the null; shown for the accession-cluster interval of the pre-specified protein-cluster bootstrap and for the gene-cluster interval, which is primary from the amendment). "
  "**Direction-consistent** is a second column: the adjusted estimate lies on the same side of the null in every population where the signal is defined (G1, G2, G3, P-abund in human, mouse and rat under D1, D4 and D5; the estimates used are in the next table). "
  "**Two-way** is the protein-by-study-cluster interval of B4 (G1 D1 only; not run for the other populations); **study share** is V_study / V_two-way. The gene-cluster rule stays primary (decision D4).\n")
rows = []
for x in main + abund:
    s = x["s"]
    tw = x["tw"]
    rows.append([x["block"], SL[s] + ("" if x["db"] else " (not decision-bearing)"), x["pop"], est_ci(x["gene0"]["estimate"], x["gene0"]["ci"], s), est_ci(x["gene"]["estimate"], x["gene"]["ci"], s),
                 ci(x["acc"]["ci"], s), yn(excl(x["acc"]["ci"], s)), yn(x["gene"]["survives"]), pct(x["gene"]["attenuation"]),
                 (ci(tw["ci_two_way"], s) + f" ({yn(tw['survives_two_way'])})") if tw else "not run", pct(tw["V_study"] / tw["V_two_way"]) if tw else "n/a", yn(consist[s]["consistent"]), label(s)])
A(md(["Block", "Signal", "Population", "Crude (gene-cluster CI)", "D1 adjusted (gene-cluster CI)", "D1 accession-cluster CI", "Survives (accession)", "Survives (gene)", "Attenuation (gene)", "Two-way CI (excludes null)", "Study share of variance", "Direction-consistent", "Label"], rows))
tw_ex = [SL[x["s"]] for x in tw_dec if not x["tw"]["survives_two_way"]]
A(f"**Reading.** Excluding the null: {surv_acc} of {len(n_dec)} under the accession-cluster interval, {surv_gene} of {len(n_dec)} under the gene-cluster interval, {tw_surv} of {len(tw_dec)} under the two-way interval (G1 contrasts that do not: {', '.join(tw_ex) if tw_ex else 'none'}). "
  f"The study-cluster component carries {pct(min(x['tw']['V_study'] / x['tw']['V_two_way'] for x in tw_dec))} to {pct(max(x['tw']['V_study'] / x['tw']['V_two_way'] for x in tw_dec))} of the two-way variance. The two-way interval rests on {min(x['tw']['n_study_clusters'] for x in tw_dec)} to {max(x['tw']['n_study_clusters'] for x in tw_dec)} study clusters, "
  "a few of them very large, so it is approximate (INTERIM_REPORT.md, section 7); it is a sensitivity, not the primary rule.\n")
rows = []
for s in SL:
    c = consist[s]
    cells = {f"{p} {d}": v for p, d, v in c["rows"]}
    rows.append([SL[s], len(c["rows"]), "; ".join(f"{k}: {f(v, s)}" for k, v in cells.items()), ", ".join(sorted(c["sides"])), yn(c["consistent"]), label(s)])
A("Adjusted estimates used for the direction column (all from the CEM results; G1 and G2 are D1, P-abund rows are D1, D4 and D5):\n")
A(md(["Signal", "Estimates", "Values", "Side of null", "Direction-consistent", "Label"], rows))

# stratified domain overlap
dc = g20["domain_overlap_by_annotation_coverage"]
A("### Domain overlap by domain-annotation coverage of the protein (G20-0489, G20-0490)\n")
A("An accession is annotated if at least one of its human site records has `in_domain` = 1. Crude, GOLD against BRONZE, set-B tiers, human records of P-all; gene-cluster intervals.\n")
rows = []
for k in ("all proteins", "at least one domain-annotated site", "no domain-annotated site"):
    v = dc[k]
    rows.append([k, n_(v["n_gold"]), n_(v["n_bronze"]), pct(v["share_in_domain_gold"]), pct(v["share_in_domain_bronze"]), est_ci(v["crude_or"], v["ci"], "S1") if not nan(v["crude_or"]) else "undefined (" + v["note"] + ")"])
A(md(["Stratum", "GOLD records", "BRONZE records", "In-domain share, GOLD", "In-domain share, BRONZE", "Crude OR (gene-cluster CI)"], rows))
cp = dc["composition"]
A(f"**The crude odds ratio changes sign across populations:** {f(dc['all proteins']['crude_or'], 'S1')} (all proteins), {f(abund[0]['gene0']['estimate'], 'S1')} (proteins with an abundance, P-abund), {f(dc['at least one domain-annotated site']['crude_or'], 'S1')} (proteins with at least one annotated site; G3 gives the same value), "
  f"{f(am['G16-0361']['estimate'], 'S1')} (canonical reviewed accessions, G16). Where the protein has any domain annotation, GOLD sites are *less* often in a domain than BRONZE sites ({pct(dc['at least one domain-annotated site']['share_in_domain_gold'])} against {pct(dc['at least one domain-annotated site']['share_in_domain_bronze'])}); where it has none the outcome is 0 for every record and no odds ratio exists. "
  f"The pooled odds ratio above 1 comes from composition: {pct(cp['share_gold_records_in_annotated_proteins'])} of GOLD records but {pct(cp['share_bronze_records_in_annotated_proteins'])} of BRONZE records lie in proteins with an annotated site (odds ratio {f(cp['or_gold_vs_bronze_of_being_in_an_annotated_protein'], 'S1')} for being in such a protein, interval {ci(cp['ci'], 'S1')}). "
  "**Domain overlap is therefore labelled population-dependent** (the other signals are direction-consistent, table above); the matched estimate for domain overlap is not a population-stable effect and should not be read as one.\n")

# ---------------------------------------------------------------- 4 G20
A("## 4. Positive control and index diagnostics (block G20)\n")
pb = g20["prematching_balance"]
A("### 4a. Balance before matching\n")
A(f"SMD of the primary index (V50) between GOLD and BRONZE records, before and after CEM, for every D1 and D3 row that was run ({len(pb)} contrasts; decision-bearing flag from the manifest). Balance requirement |SMD| < 0.1 after matching.\n")
rows = []
for k, v in pb.items():
    rows.append([k, v["block"], SL[v["signal"]], v["population"], "yes" if manid.loc[k, "decision_bearing"] == "yes" else "no", f(v["smd_pre"], "S2", 3), f(v["smd_after"], "S2", 3)])
A(md(["Contrast", "Block", "Signal", "Population", "Decision-bearing", "SMD before", "SMD after"], rows))
A(f"Across the {len(dbk)} decision-bearing contrasts the SMD before matching ranges from {pre.min():.3f} to {pre.max():.3f}; {int((np.abs(pre) >= 0.1).sum())} have |SMD| of 0.1 or more before matching and {int((np.abs(aft) >= 0.1).sum())} after (largest after: {np.abs(aft).max():.3f}). "
  f"For the primary domain-overlap contrast (G1-0002) the index is already balanced before matching ({pb['G1-0002']['smd_pre']:.4f}), so matching has nothing to adjust there.\n")
A("### 4b. Crude association of the index with each outcome within a tier\n")
A("Record level, Pearson correlation of the V50 index with the outcome (equal to the point-biserial correlation for binary outcomes), within GOLD records and within BRONZE records, gene-cluster 95% interval.\n")
rows = []
for k, v in ix.items():
    rows.append([k, SL[v["signal"]], v["population"], v["tier"], n_(v["n_records"]), n_(v["n_gene_clusters"]), est_ci(v["pearson_r"], v["ci"], "S2", 3), yn(excl(v["ci"], "S2"))])
A(md(["Contrast", "Signal", "Population", "Tier", "Records", "Gene clusters", "Correlation (CI)", "Excludes 0"], rows))
opp = [k for k in ix if ix[k]["tier"] == "GOLD" and excl(ix[k]["ci"], "S2") and excl(ix[f"G20-{int(k[4:]) + 1:04d}"]["ci"], "S2") and ix[k]["pearson_r"] * ix[f"G20-{int(k[4:]) + 1:04d}"]["pearson_r"] < 0]
A(f"The correlations lie between {rmin:.3f} and {rmax:.3f}. In {len(opp)} of {sum(1 for k in ix if ix[k]['tier'] == 'GOLD')} signal-by-population combinations the correlation excludes zero in both GOLD and BRONZE with opposite signs ({', '.join(SL[ix[k]['signal']] + ' ' + ix[k]['population'] for k in opp) or 'none'}), "
  "which a single index-based stratification cannot represent.\n")
A("### 4c. The index against the number of sites per protein\n")
sp = g20["index_vs_sites_per_accession"]
rows = []
for k, v in sp.items():
    pop, what = k.split("|")
    rows.append([pop, {"sites_per_accession": "all site records per accession", "gold_sites_per_accession": "GOLD records per accession"}[what], n_(v["n_accessions"]), n_(v["n_gene_clusters"]), est_ci(v["spearman"], v["ci"], "S2", 3)])
A(md(["Species set", "Count", "Accessions with an index", "Gene clusters", "Spearman rho (gene-cluster CI)"], rows))
if PSP_PRESENT:
    A("Index against the number of LT_LIT-supported sites per gene: see the Task H run.\n")
else:
    A(f"**Index against the number of LT_LIT-supported sites per gene (G20-0488): {PSP_MSG}.**\n")
A("### 4d. Could matching on the index remove confounding at all?\n")
A("Matching on a variable can remove confounding by that variable only to the extent that (i) the compared groups differ in it and (ii) it is related to the outcome within the groups. Under a linear approximation (omitted-variable formula) the shift it can cause in the GOLD-BRONZE difference, in units of the outcome SD, is about SMD x r, where r is the within-tier correlation. "
  "This is a rough guide to size, not a bias estimate; it is shown for the decision-bearing G1 and G2 contrasts:\n")
rows = []
prods = []
for x in n_dec:
    s, pop = x["s"], x["pop"]
    rg = [v["pearson_r"] for v in ix.values() if v["signal"] == s and v["population"] == pop and v["tier"] == "GOLD"]
    rb = [v["pearson_r"] for v in ix.values() if v["signal"] == s and v["population"] == pop and v["tier"] == "BRONZE"]
    smd0 = pbk[x["d1"]]["smd_pre"]
    rbar = (rg[0] + rb[0]) / 2
    prods.append(abs(smd0 * rbar))
    rows.append([SL[s], pop, f(smd0, "S2", 3), f(rg[0], "S2", 3), f(rb[0], "S2", 3), f(smd0 * rbar, "S2", 4), pct(x["gene"]["attenuation"])])
A(md(["Signal", "Population", "SMD before matching", "r within GOLD", "r within BRONZE", "SMD x mean r", "Observed attenuation (gene)"], rows))
A(f"In plain words: **matching on this index can remove only a small amount of confounding, and a near-zero attenuation does not show that study effort is not a confounder.** In the decision-bearing G1 and G2 contrasts the product is at most {max(prods):.4f} outcome SD, and the within-tier correlations of the index with the outcomes lie between {rmin:.3f} and {rmax:.3f}. "
  f"The index also tracks the number of sites per protein weakly (Spearman {sp['human|sites_per_accession']['spearman']:.3f} among human proteins) and the number of GOLD sites per protein not at all (Spearman {sp['human|gold_sites_per_accession']['spearman']:.3f}, interval {ci(sp['human|gold_sites_per_accession']['ci'], 'S2', 3)}). "
  "The index is a gene-level count of publications; the tier is a property of how a site was reported (study-cluster depth, database agreement). Matching on it does remove an imbalance that exists in many populations, but the imbalance carries little of the outcome, so the results of CEM and D2b show that **adjusting for this index leaves the signals unchanged**; they do not show that site-level study effort cannot explain them.\n")

# ---------------------------------------------------------------- 5 D2b
A("## 5. Matcher D2b: what the Rev 5 Methods describe (block G18, post-hoc)\n")
A("1:1 without replacement, caliper 0.2 SD of the logit (over the records of the resample), serine, threonine and tyrosine records only, random visiting order, 20 order replicates per contrast (min, median, max of the point estimate; the median is the estimate), the logit model and the matching recomputed in every one of 1,000 gene-cluster resamples (one random order per resample). "
  "Rev 5 kept each treated record's 50 nearest controls while its Methods say 1:1 matching: this is recorded in `RECONCILIATION.md`, section 6, and the manuscript Methods must be rewritten to describe what was run (not edited in this cycle). A point estimate outside its own percentile interval is flagged and the contrast is not reported as decision-bearing.\n")
rows = []
for cid in sorted(d2b):
    v = d2b[cid]
    s = v["signal"]
    cp_ = manid.loc[cid, "note"].split("counterpart of ")[1].split(";")[0]
    pr = prop.get(cp_)
    d2s = "not fitted" if pr and pr["status"] != "run" else (est_ci(pr["estimate"], pr["ci"], s) if pr else "n/a")
    rp = v["replicates"]
    rows.append([cid, SL[s], v["population"], est_ci(v["estimate"], v["ci"], s), f"{f(rp['min'], s)} / {f(rp['median'], s)} / {f(rp['max'], s)}", f(v["crude"], s), pct(v["attenuation"]), n_(v["n_pairs"]), pct(v["share_gold_matched"]),
                 f(v["smd"], "S2", 4), "yes" if v["point_outside_own_interval"] else "no", yn(v["survives"]), "yes" if (v["decision_bearing"] == "yes" and not v["point_outside_own_interval"]) else "no", cp_ + ": " + d2s])
A(md(["Contrast", "Signal", "Population", "D2b estimate (95% CI)", "Replicates min / median / max", "Crude", "Attenuation", "Matched pairs", "Share of GOLD matched", "SMD", "Point outside own CI", "Survives", "Decision-bearing", "Pre-specified D2 counterpart"], rows))
A(f"{d2b_surv} of {len(d2b_dec)} decision-bearing D2b contrasts survive; the largest attenuation among them is {pct(d2b_att)}; {d2b_flag} of {len(d2b)} point estimates lie outside their own interval. "
  f"Order variability (max minus min over the 20 replicates) is at most {max(v['replicates']['max'] - v['replicates']['min'] for v in d2b.values() if v['signal'] != 'S1'):.4f} on the difference scale and {max(v['replicates']['max'] - v['replicates']['min'] for v in d2b.values() if v['signal'] == 'S1'):.3f} on the odds-ratio scale, far narrower than the intervals. "
  f"In {sum(1 for v in d2b.values() if v['share_gold_matched'] < 0.999)} of {len(d2b)} contrasts only part of GOLD is matched (share {pct(min(v['share_gold_matched'] for v in d2b.values()))} to {pct(max(v['share_gold_matched'] for v in d2b.values() if v['share_gold_matched'] < 0.999))}); there the pairs number {pct(min(v['n_pairs'] / v['n_bronze'] for v in d2b.values() if v['share_gold_matched'] < 0.999))} to {pct(max(v['n_pairs'] / v['n_bronze'] for v in d2b.values() if v['share_gold_matched'] < 0.999))} of the BRONZE records (the number of controls and the caliper limit the matching under 1:1 without replacement), and the estimate refers to the matched part of GOLD. "
  "The Rev 5 matcher results (including its non-fits) remain in INTERIM_REPORT.md, section 5.\n")

# ---------------------------------------------------------------- 6 G19
A("## 6. Tier-definition sensitivity (block G19, post-hoc)\n")
A("Crude and CEM-adjusted (D1, V50) estimates under four tier definitions: set-B tiers (primary), set-A tiers (human records only), depth-only tiers (GOLD depth >= 10, SILVER 3 to 9, BRONZE <= 2) and an agreement-only tier (GOLD agreement >= 2 against BRONZE agreement 0; records with agreement 1 not compared). Gene-cluster intervals.\n")
tiers = {}
for cid, v in am.items():
    if cid.startswith("G19"):
        key = (v["signal"], v["population"], v["tier_rule"].split(" (")[0].replace(" (primary)", ""))
        tiers.setdefault(key, {})[v["design"]] = v
TNAME = {"set-B tiers": "set-B tiers", "set-A tiers": "set-A tiers", "depth-only tiers": "depth-only tiers", "agreement-only": "agreement-only"}
for pop in ("P-human", "P-all"):
    rows = []
    for s in ("S1", "S2", "S3", "S4"):
        for k in [k for k in tiers if k[0] == s and k[1] == pop]:
            d0, d1 = tiers[k]["D0"], tiers[k]["D1"]
            rows.append([SL[s], k[2], n_(d0["n_sites_gold"]), n_(d0["n_sites_bronze"]), est_ci(d0["estimate"], d0["ci"], s), est_ci(d1["estimate"], d1["ci"], s), pct(d1["attenuation"]), yn(d1["survives"])])
    A(f"**{pop}**\n")
    A(md(["Signal", "Tier definition", "GOLD", "BRONZE", "Crude (CI)", "D1 adjusted (CI)", "Attenuation", "Survives"], rows))


def excl_null(c, s):
    return bool(c[0] > NULL(s) or c[1] < NULL(s))


def eligible_spread(c, s):
    """Rule for the spread column: the interval excludes the null and no limit of it lies within 0.01 of the null."""
    return excl_null(c, s) and min(abs(c[0] - NULL(s)), abs(c[1] - NULL(s))) >= 0.01


A("### Spread of the estimates across tier definitions, compared with the change caused by effort adjustment\n")
A("Two comparisons, both on adjusted estimates from the same four (three for P-all) tier definitions. **(1) Spread, as ratio of the largest to the smallest adjusted estimate:** over all definitions, and excluding any estimate whose interval includes the null or has a limit within 0.01 of the null (the interval reading of \"within 0.01 of zero\"; the point-estimate reading is noted below the table). "
  "**(2) Like-for-like, on the log scale:** the largest relative change from crude to adjusted within one tier definition, max over definitions of |ln(adjusted / crude)|, against the spread of the adjusted estimates across definitions, ln(largest / smallest). Only estimates whose interval excludes the null enter (for the within-definition change both the crude and the adjusted estimate; for the spread the adjusted estimates). "
  "For the odds ratio the log scale is the log odds ratio.\n")
sumrows = []
lf = []
alt_excl = []
int_only = []
for pop in ("P-human", "P-all"):
    for s in ("S1", "S2", "S3", "S4"):
        ks = [k for k in tiers if k[0] == s and k[1] == pop]
        if not ks:
            continue
        defs = [(k[2], tiers[k]["D0"], tiers[k]["D1"]) for k in ks]
        adj = [d1["estimate"] for _, _, d1 in defs]
        spread_all = max(adj) / min(adj)
        el = [(n, d0, d1) for n, d0, d1 in defs if eligible_spread(d1["ci"], s)]
        gone = [n for n, d0, d1 in defs if not eligible_spread(d1["ci"], s)]
        for n, d0, d1 in defs:
            if eligible_spread(d1["ci"], s) and abs(d1["estimate"] - NULL(s)) < 0.01:
                alt_excl.append(f"{SL[s]} {pop} {n} ({f(d1['estimate'], s)})")
            if (not eligible_spread(d1["ci"], s)) and excl_null(d1["ci"], s) and abs(d1["estimate"] - NULL(s)) >= 0.01:
                int_only.append(f"{SL[s]} {pop} {n} (estimate {f(d1['estimate'], s)}, interval {ci(d1['ci'], s)})")
        spread_el = max(d1["estimate"] for _, _, d1 in el) / min(d1["estimate"] for _, _, d1 in el) if len(el) >= 2 else None
        okc = [(n, d0, d1) for n, d0, d1 in defs if excl_null(d0["ci"], s) and excl_null(d1["ci"], s) and d0["estimate"] > 0 and d1["estimate"] > 0]
        L_adj = max(abs(math.log(d1["estimate"] / d0["estimate"])) for _, d0, d1 in okc) if okc else None
        oka = [d1["estimate"] for _, d0, d1 in defs if excl_null(d1["ci"], s) and d1["estimate"] > 0]
        L_def = math.log(max(oka) / min(oka)) if len(oka) >= 2 else None
        lf.append((pop, s, L_adj, L_def))
        sumrows.append([pop, SL[s], len(defs), f"{spread_all:.2f}", ("; ".join(gone) + " excluded") if gone else "none excluded", f"{spread_el:.2f}" if spread_el else "n/a",
                        f"{L_adj:.3f}" if L_adj is not None else "n/a", f"{math.exp(L_adj):.3f}" if L_adj is not None else "n/a", f"{L_def:.3f}" if L_def is not None else "n/a", f"{math.exp(L_def):.2f}" if L_def is not None else "n/a",
                        ("yes" if L_def > L_adj else "no") if (L_adj is not None and L_def is not None) else "n/a"])
A(md(["Population", "Signal", "Tier definitions", "Spread, all adjusted estimates", "Estimates excluded (CI includes the null or is within 0.01 of it)", "Spread, excluding those", "Largest within-definition change, max abs(ln(adjusted / crude)), natural log", "same change as a ratio, exp(change)", "ln spread of adjusted estimates across definitions, natural log", "same spread as a ratio, exp(ln spread)", "Spread larger than the within-definition change"], sumrows))
lfd = [t for t in lf if t[2] is not None and t[3] is not None]
sp_all_ex = [tiers[k]["D1"]["estimate"] for k in tiers if k[0] == "S3" and k[1] == "P-all" and eligible_spread(tiers[k]["D1"]["ci"], "S3")]
tier_bigger = sum(1 for t in lfd if t[3] > t[2])
A(f"Like-for-like, the log-scale spread across tier definitions exceeds the largest log-scale change caused by effort adjustment within a definition in {tier_bigger} of the {len(lfd)} signal-by-population combinations where both can be computed "
  f"(natural-log scale; spread {min(t[3] for t in lfd):.3f} to {max(t[3] for t in lfd):.3f}, as ratios {math.exp(min(t[3] for t in lfd)):.2f} to {math.exp(max(t[3] for t in lfd)):.2f}; within-definition change {min(t[2] for t in lfd):.3f} to {max(t[2] for t in lfd):.3f}, as ratios {math.exp(min(t[2] for t in lfd)):.3f} to {math.exp(max(t[2] for t in lfd)):.3f}). "
  + f"The largest spread ({max(t[3] for t in lfd):.3f}, a {math.exp(max(t[3] for t in lfd)):.2f}-fold ratio; disorder in P-all) rests on a near-zero agreement-only estimate; without it the spread is about {max(sp_all_ex) / min(sp_all_ex):.0f}-fold ({max(sp_all_ex) / min(sp_all_ex):.2f}). "
  + "The interval reading excludes the estimate if either limit of its interval lies within 0.01 of the null"
  + (f"; this excludes, although its point estimate is at least 0.01 from the null (it would stay in under the point-estimate reading): {', '.join(int_only)}" if int_only else "")
  + (f"; under the point-estimate reading the following would be excluded in addition: {', '.join(alt_excl)}" if alt_excl else "; the point-estimate reading would exclude no further estimate")
  + ".\n")
dis_all = tiers[("S3", "P-all", "agreement-only")]["D1"]
dis_h = tiers[("S3", "P-human", "agreement-only")]["D1"]
sp_all = [tiers[k]["D1"]["estimate"] for k in tiers if k[0] == "S3" and k[1] == "P-all"]
sp_h = [tiers[k]["D1"]["estimate"] for k in tiers if k[0] == "S3" and k[1] == "P-human"]
sp_all_ex = [tiers[k]["D1"]["estimate"] for k in tiers if k[0] == "S3" and k[1] == "P-all" and eligible_spread(tiers[k]["D1"]["ci"], "S3")]
sp_h_ex = [tiers[k]["D1"]["estimate"] for k in tiers if k[0] == "S3" and k[1] == "P-human" and eligible_spread(tiers[k]["D1"]["ci"], "S3")]
A(f"**The {max(sp_all) / min(sp_all):.0f}-fold spread for disorder in P-all rests on one near-zero estimate:** the agreement-only adjusted estimate is {f(dis_all['estimate'], 'S3')} (interval {ci(dis_all['ci'], 'S3')}), and the crude agreement-only estimate {f(tiers[('S3', 'P-all', 'agreement-only')]['D0']['estimate'], 'S3')} has an interval that includes zero. "
  f"Without it the spread of the disorder estimates in P-all is about {max(sp_all_ex) / min(sp_all_ex):.0f}-fold ({max(sp_all_ex) / min(sp_all_ex):.2f}). The same holds for P-human ({max(sp_h) / min(sp_h):.1f}-fold with the agreement-only estimate {f(dis_h['estimate'], 'S3')}, about {max(sp_h_ex) / min(sp_h_ex):.0f}-fold ({max(sp_h_ex) / min(sp_h_ex):.2f}) without it).\n")
tn = [t for t in lfd]

# --- 6b pillar decomposition
A("### 6b. Depth-only and agreement-only tiers side by side (no new data; from G19)\n")
A("The tier combines study-cluster depth and database agreement; the two tiers below each use one of them alone. Crude and CEM-adjusted estimates with gene-cluster intervals, same index and design as above. **This is a description of how the estimates differ between the two tier definitions; it carries no causal interpretation.**\n")
pr = []
for pop in ("P-human", "P-all"):
    for s in ("S1", "S2", "S3", "S4"):
        kd, ka = (s, pop, "depth-only tiers"), (s, pop, "agreement-only")
        if kd not in tiers or ka not in tiers:
            continue
        D, Ag = tiers[kd], tiers[ka]
        pr.append([pop, SL[s], f"{n_(D['D0']['n_sites_gold'])} / {n_(D['D0']['n_sites_bronze'])}", est_ci(D["D0"]["estimate"], D["D0"]["ci"], s), est_ci(D["D1"]["estimate"], D["D1"]["ci"], s),
                   f"{n_(Ag['D0']['n_sites_gold'])} / {n_(Ag['D0']['n_sites_bronze'])}", est_ci(Ag["D0"]["estimate"], Ag["D0"]["ci"], s), est_ci(Ag["D1"]["estimate"], Ag["D1"]["ci"], s)])
A(md(["Population", "Signal", "Depth-only: GOLD / BRONZE", "Depth-only crude (CI)", "Depth-only adjusted (CI)", "Agreement-only: GOLD / BRONZE", "Agreement-only crude (CI)", "Agreement-only adjusted (CI)"], pr))
dd, aa = tiers[("S3", "P-all", "depth-only tiers")], tiers[("S3", "P-all", "agreement-only")]
dh, ah = tiers[("S3", "P-human", "depth-only tiers")], tiers[("S3", "P-human", "agreement-only")]
A(f"For disorder, the crude estimate is {f(dd['D0']['estimate'], 'S3')} under the depth-only tier and {f(aa['D0']['estimate'], 'S3')} (interval {ci(aa['D0']['ci'], 'S3')}, includes zero) under the agreement-only tier in P-all; in P-human {f(dh['D0']['estimate'], 'S3')} and {f(ah['D0']['estimate'], 'S3')} (interval {ci(ah['D0']['ci'], 'S3')}, includes zero). "
  f"The adjusted estimates are {f(dd['D1']['estimate'], 'S3')} and {f(aa['D1']['estimate'], 'S3')} (P-all) and {f(dh['D1']['estimate'], 'S3')} and {f(ah['D1']['estimate'], 'S3')} (P-human).\n")

# ---------------------------------------------------------------- 7 G16
A("## 7. Canonical accessions (block G16, post-hoc)\n")
A(f"Restricted to Swiss-Prot reviewed accessions (release 2026_03), one per gene cluster (the accession with the most records where a cluster holds several; {g16c['gene_clusters_with_more_than_one_reviewed_accession']} clusters). **Coverage:** {n_(g16c['accessions_kept'])} accessions kept; human {n_(g16c['human']['accessions_kept'])} of {n_(g16c['human']['accessions'])} accessions and "
  f"{n_(g16c['human']['records_kept'])} of {n_(g16c['human']['records'])} records ({pct(g16c['human']['share_records_kept'])}); mouse {pct(g16c['mouse']['share_records_kept'])} of records; rat {pct(g16c['rat']['share_records_kept'])}. Gene-cluster intervals; the unrestricted gene-cluster result (G17) is shown for comparison.\n")
rows = []
for cid in sorted(k for k in am if k.startswith("G16")):
    v = am[cid]
    if v["design"] != "D1":
        continue
    s = v["signal"]
    unr = [x for x in main if x["s"] == s and (x["pop"] == v["population"].split(",")[0].replace(" (human)", "") or x["pop"].startswith(v["population"].split(",")[0]))]
    unr = unr[0] if unr else None
    c0 = am[cid[:4] + f"{int(cid[4:]) - 1:04d}"]
    rows.append([cid, SL[s], v["population"], n_(v["n_sites_gold"]), n_(v["n_sites_bronze"]), est_ci(c0["estimate"], c0["ci"], s), est_ci(v["estimate"], v["ci"], s), pct(v["attenuation"]), yn(v["survives"]),
                 est_ci(unr["gene"]["estimate"], unr["gene"]["ci"], s) if unr else "n/a", f(v["smd_pre"], "S2", 3) + " / " + f(v["smd"], "S2", 3)])
A(md(["Contrast", "Signal", "Population", "GOLD", "BRONZE", "Crude (CI)", "D1 adjusted (CI)", "Attenuation", "Survives", "Unrestricted D1 (gene CI, G17)", "SMD before / after"], rows))
can = [v for v in am.values() if v["design"] == "D1" and v["variation"] == "canonical accessions only" and v["signal"] != "S1"]
can_same = sum(1 for v in can if v["survives"] and (v["estimate"] > NULL(v["signal"])))
A(f"{can_same} of the {len(can)} contrasts other than domain overlap survive with a positive estimate (as unrestricted). **Domain overlap reverses**: the odds ratio is {f(am['G16-0362']['estimate'], 'S1')} ({ci(am['G16-0362']['ci'], 'S1')}) against {f(main[0]['gene']['estimate'], 'S1')} unrestricted. "
  "The reversal is consistent with the composition account of section 3 (GOLD records sit disproportionately in proteins with a domain-annotated site), but I have not tested that account on the canonical set separately, so it is stated as consistent with the result, not as its demonstrated cause.\n")

# ---------------------------------------------------------------- 8 Task F
A("## 8. Stratification and sensitivity (Task F: G6, G7, G8)\n")
A("Gene-cluster bootstrap (primary) and accession-cluster bootstrap beside it, 1,000 resamples each, in `results/f_runs/`. A sensitivity analysis changes a decision if it moves a contrast across the survives/not or the largely-explained/not boundary relative to the primary design (G1 D1 for the same signal); the count of changes is reported per analysis, for both intervals. A change of sign relative to the primary is shown separately because the decision rule does not detect it.\n")
prim = {x["s"]: x["gene"] for x in main if x["block"] == "G1"}
prim_acc = {x["s"]: x["acc"] for x in main if x["block"] == "G1"}


def dec(x, s):
    return (bool(x["survives"]), bool(x.get("largely_explained", False)))


def f_rows(block, label_of):
    rows, chg_g, chg_a, flips, n_ = [], 0, 0, 0, 0
    for cid, v in fr.items():
        m = manid.loc[cid]
        if m["block"] != block or m["design"] != "D1":
            continue
        s = m["signal"]
        d0 = [i for i in fr if manid.loc[i, "block"] == block and manid.loc[i, "signal"] == s and manid.loc[i, "variation"] == m["variation"] and manid.loc[i, "population"] == m["population"] and manid.loc[i, "design"] == "D0"][0]
        g, a = v["gene"], v["protein"]
        if g["status"] != "run":
            rows.append([label_of(m), SL[s], "not run", g.get("reason", "")] + [""] * 6)
            continue
        n_ += 1
        cg = dec(g, s) != dec(prim[s], s)
        ca = dec(a, s) != dec(prim_acc[s], s)
        fl = (g["estimate"] > NULL(s)) != (prim[s]["estimate"] > NULL(s))
        chg_g += cg
        chg_a += ca
        flips += fl
        rows.append([label_of(m), SL[s], est_ci(fr[d0]["gene"]["estimate"], fr[d0]["gene"]["ci"], s), est_ci(g["estimate"], g["ci"], s), ci(a["ci"], s), pct(g["attenuation"]), yn(g["survives"]), yn(g.get("largely_explained", False)),
                     n_k(g), "changed" if cg else "same", "opposite sign" if fl else "same sign"])
    return rows, n_, chg_g, chg_a, flips


def n_k(g):
    return f"{int(g['n_sites_gold']):,} / {int(g['n_sites_bronze']):,}"


H = ["Analysis", "Signal", "Crude (gene CI)", "D1 adjusted (gene CI)", "D1 accession CI", "Attenuation", "Survives (gene)", "Largely explained", "GOLD / BRONZE", "Decision vs G1 D1 (gene)", "Sign vs G1 D1"]
rows, nn, cg, ca, fl = f_rows("G6", lambda m: m["population"])
A("### G6: residue types separately\n")
A("P-S, P-T and P-Y are the P-all records restricted to serine, threonine and tyrosine.\n")
A(md(H, rows))
A(f"{nn} contrasts run; {cg} change a decision under the gene-cluster interval and {ca} under the accession-cluster interval; {fl} have the opposite sign to the primary ({', '.join(sorted({SL[manid.loc[i, 'signal']] + ' ' + manid.loc[i, 'population'] for i in fr if manid.loc[i, 'block'] == 'G6' and manid.loc[i, 'design'] == 'D1' and (fr[i]['gene']['estimate'] > NULL(manid.loc[i, 'signal'])) != (prim[manid.loc[i, 'signal']]['estimate'] > NULL(manid.loc[i, 'signal']))})) or 'none'}). Reported as found; I have not examined why.\n")
g22 = J("g22_results.json")
A("### Sequence-based motif indicator: its four classes, and a serine/threonine-only contrast (post-hoc, G22)\n")
A("The indicator (`scripts/epsd_phase0/phase1_6_sequence_motif.py`) counts a site as matching if its flanking window satisfies any one of four consensus classes. As written in that source:\n")
A(md(["Class (label in the source)", "Rule, positions relative to the phosphoacceptor", "Pattern as written in the source", "Written as a Ser/Thr motif?"],
     [["Basophilic, PKA-like", "R at -3", "[R]..{S,T}", "yes (acceptor is {S,T})"], ["Basophilic, PKC-like", "R or K at -2", "[RK].{S,T}", "yes (acceptor is {S,T})"],
      ["Proline-directed, CDK/MAPK", "P at +1", "{S,T}P", "yes (acceptor is {S,T})"], ["Acidophilic, CK2-like", "D or E at +3", "{S,T}..[DE]", "yes (acceptor is {S,T})"]]))
A("All four are written as serine/threonine motifs and none is written for tyrosine. The code nevertheless evaluates them for any serine, threonine or tyrosine centre residue, so for tyrosine records the indicator applies a Ser/Thr consensus. "
  "The tyrosine contrast (G6, P-Y) is reported as it was run, and I have not examined why its sign differs.\n")
rows = []
pa_ = {x["s"]: x for x in main if x["block"] == "G1"}["S4"]
rows.append(["P-all (G1; serine, threonine and tyrosine together)", f"{n_(pa_['gene']['n_sites_gold'])} / {n_(pa_['gene']['n_sites_bronze'])}", est_ci(pa_["gene0"]["estimate"], pa_["gene0"]["ci"], "S4"), est_ci(pa_["gene"]["estimate"], pa_["gene"]["ci"], "S4"), pct(pa_["gene"]["attenuation"]), "primary, as reported"])
for popl, lab in (("P-S", "serine only (G6)"), ("P-T", "threonine only (G6)"), ("P-Y", "tyrosine only (G6)")):
    i1 = [i for i in fr if manid.loc[i, "block"] == "G6" and manid.loc[i, "signal"] == "S4" and manid.loc[i, "population"] == popl and manid.loc[i, "design"] == "D1"][0]
    i0 = [i for i in fr if manid.loc[i, "block"] == "G6" and manid.loc[i, "signal"] == "S4" and manid.loc[i, "population"] == popl and manid.loc[i, "design"] == "D0"][0]
    rows.append([lab, f"{n_(fr[i1]['gene']['n_sites_gold'])} / {n_(fr[i1]['gene']['n_sites_bronze'])}", est_ci(fr[i0]["gene"]["estimate"], fr[i0]["gene"]["ci"], "S4"), est_ci(fr[i1]["gene"]["estimate"], fr[i1]["gene"]["ci"], "S4"), pct(fr[i1]["gene"]["attenuation"]), "as reported"])
d0_, d1_ = g22["G22-0515"], g22["G22-0516"]
rows.append(["**serine and threonine together (G22, post-hoc)**", f"{n_(d1_['n_sites_gold'])} / {n_(d1_['n_sites_bronze'])}", est_ci(d0_["estimate"], d0_["ci"], "S4"), est_ci(d1_["estimate"], d1_["ci"], "S4"), pct(d1_["attenuation"]), "post-hoc, not decision-bearing; justified by the definition above"])
A(md(["Population", "GOLD / BRONZE", "Crude (gene CI)", "D1 adjusted (gene CI)", "Attenuation", "Status"], rows))
A(f"The serine/threonine-only difference is {f(d1_['estimate'], 'S4')} (interval {ci(d1_['ci'], 'S4')}), against {f(pa_['gene']['estimate'], 'S4')} for P-all; the interval excludes zero ({yn(d1_['survives'])}) and the attenuation by the index is {pct(d1_['attenuation'])}. The P-all and P-Y results stay as reported.\n")

rows, nn, cg, ca, fl = f_rows("G7", lambda m: m["variation"])
A("### G7: leave one study cluster (or the seven-PMID consortium) out\n")
A("Removed PMIDs are deleted from every site's publication list, study-cluster depth and tier are recomputed, sites with no PMID and no database source flag are dropped. `drop_consortium_7_PMIDs` duplicates `drop_cluster_2` by construction (deviation log) and counts once in the tally below.\n")
A(md(H, rows))
dup = sum(1 for cid in fr if manid.loc[cid, "variation"] == "drop_consortium_7_PMIDs" and manid.loc[cid, "design"] == "D1")
A(f"{nn} D1 contrasts run ({dup} of them the duplicate removal); {cg} change a decision under the gene-cluster interval and {ca} under the accession-cluster interval; {fl} have the opposite sign. "
  "Largest movement of the adjusted estimate relative to G1 D1, per signal: " + "; ".join(
      f"{SL[s]} {manid.loc[i, 'variation']} ({f(fr[i]['gene']['estimate'], s)} against {f(prim[s]['estimate'], s)}, {100 * (fr[i]['gene']['estimate'] / prim[s]['estimate'] - 1):+.0f}%)"
      for s in ("S1", "S2", "S3", "S4")
      for i in [max((i_ for i_ in fr if manid.loc[i_, 'block'] == 'G7' and manid.loc[i_, 'design'] == 'D1' and manid.loc[i_, 'signal'] == s and fr[i_]['gene']['status'] == 'run'), key=lambda i_: abs(fr[i_]['gene']['estimate'] / prim[s]['estimate'] - 1))]) + ".\n")
rows, nn, cg, ca, fl = f_rows("G8", lambda m: m["variation"])
A("### G8: tier thresholds (GOLD depth in {5, 10, 20} x GOLD agreement in {>= 2, >= 3})\n")
A(md(H, rows))
A(f"{nn} D1 contrasts run; {cg} change a decision under the gene-cluster interval and {ca} under the accession-cluster interval; {fl} have the opposite sign. The primary rule (depth 10, agreement 2) reproduces G1 (`results/f_tiers_summary.json`: {'check passed' if J('f_tiers_summary.json')['check_g8_d10_a2_equals_primary_tier'] else 'check failed'}).\n")

# ---------------------------------------------------------------- 9 Task H
A("## 9. Established confidence measures (Task H: G9)\n")
A("### PeptideAtlas build 606: can site-level evidence be reconstructed? (instruction of 2026-10-04)\n")
pf = J("peptideatlas_606_feasibility.json")
tb = pa["tables"]
nm = pf["peptide_instance_n_protein_mappings"]
A(f"The build (`atlas_build_606.tsv.zip`, a phosphopeptide-enrichment build, `is_ptm_build` = Y) holds {len(tb)} tables (names and columns of all of them in `results/peptideatlas_606_inspection.json`). "
  f"`modified_peptide_instance` has {n_(pf['modified_peptide_instances'])} rows, {n_(pf['with_phospho_token'])} with a `[Phospho]` token. What a site-level evidence category would need, and whether the dump has it (`scripts/analysis/peptideatlas_feasibility.py`; descriptive counts only, no category derived):\n")
rows = [
    ["Site identity (accession, residue, position)", "`modified_peptide_instance.modified_peptide_sequence` (residue carrying `[Phospho]`), `.peptide_instance_id` to `peptide_mapping` (`matched_biosequence_id`, `start_in_biosequence`) to `biosequence` (`biosequence_accession`, `biosequence_seq`)",
     f"Derivable. But {n_(nm.get('2 or more', 0))} of {n_(sum(nm.values()))} peptide instances map to two or more protein sequences (`peptide_instance.n_protein_mappings`; isoform and other-sequence handling not examined)."],
    ["Site-localisation probability or score (PTMProphet, Ascore, FLR)", "none", "**No.** Across all tables the only columns that look like confidence are build-level `probability_threshold`, peptide-level `best_probability` and `best_adjusted_probability`, and protein-level `probability` and `confidence` (`columns_matching_localisation_or_site_confidence_pattern` in the feasibility file)."],
    ["Unambiguous placement of the phosphate", "derived: peptide instance with one placement per phospho count",
     f"**Not guaranteed.** {pct(pf['share_groups_with_more_than_one_placement'])} of the {n_(pf['peptide_by_phospho_count_groups'])} (peptide, phospho count) groups appear with more than one placement; {pct(pf['share_on_residues_other_than_STY'])} of the {n_(pf['phospho_tokens_total'])} phospho tokens sit on residues other than S, T or Y ({n_(pf['phospho_tokens_by_residue'].get('A', 0))} on alanine). The alanine tokens look like decoy placements, but nothing in the dump documents that; this is an inference, not verified."],
    ["Number of distinct peptides per site", "count of distinct modified peptide instances per site", "Derivable, but instances with competing placements would count towards several sites."],
    ["Number of PSMs per site", "`modified_peptide_instance.n_observations` (no nulls), `n_adjusted_observations`", "Exists. The dump has no column description, so that these count PSMs is the build's convention and is not verified here."],
    ["Probability", "`modified_peptide_instance.best_probability` (no nulls), `best_adjusted_probability`", f"Peptide-level only, not a localisation probability; `best_adjusted_probability` is empty in {n_(pf['null_counts_in_modified_peptide_instance'].get('best_adjusted_probability', 0))} of {n_(pf['modified_peptide_instances'])} rows."],
    ["Independent samples and search batches", "`sample_ids`, `search_batch_ids`", "Exist."],
]
A(md(["Needed", "Where it would come from", "In this dump?"], rows))
A("**Decision: not feasible from the available build.** The published definitions were read (Kalyuzhnyy et al., J Proteome Res 2022; University of Liverpool repository copy; verbatim text in `sources/SOURCES.md`, section A). For PeptideAtlas the categories are defined by counts of PSM observations of a site at PTMProphet probabilities: "
  "\"The counts of observations with a probability of >0.95 were used as \"positive\" evidence for site phosphorylation. The counts at a probability threshold of ≤0.19 were used as \"negative\" evidence in favor of a site being a nonphosphosite.\" and "
  "\"(\"High\": ≥5 positive observations which accounted for ≥10% of total observations across all probabilities; \"Medium\": ≥5 positive observations, which accounted for <10% of total observations or 2−4 positive observations; \"Low\": 1 positive observation; \"Not phosphorylated\": 0 positive observations and ≥5 negative observations; \"Other\": site did not fall into any described categories).\" "
  "Applying them needs a PTMProphet probability for every PSM of every site. The TSV dump of build 606 has neither a PTMProphet probability nor a per-site PSM table (feasibility table above), so the categories cannot be applied to it, and I have not built a count-based surrogate (instruction of 2026-10-04).\n")
wc_ = J("peptideatlas_606_web_check.json")
A("**What else I checked on the web (2026-10-04; URLs, quotes and outcomes in `sources/SOURCES.md`).** "
  f"No per-site PeptideAtlas file is published on the pages I could read: the build page lists {len(wc_['build_page']['files'])} files (FASTA, TSV and MySQL dumps, coordinate files), none per-site, and shows no licence terms. "
  f"The build-details page shows that the server holds PTMProphet-based evidence (`{'`, `'.join(wc_['build_details_page']['ptmprophet_columns'])}`; \"nP.95-.99: Number of observed PSM with PTMProphet probability >= 0.95 and < 0.99\"), but only as a summary per sequence database and residue (downloaded, {n_(wc_['build_details_page']['saved_bytes'])} bytes) and as a per-protein table ({n_(wc_['protein_level_ptm_summary']['rows_as_displayed'])} proteins, PSM counts aggregated per protein); the bands (>=0.80, >=0.95, >=0.99) do not match the paper's \">0.95\" and \"<=0.19\". "
  "Per-site values would mean scraping a web service, which I did not do. The authors' own site-level tables are listed in the paper's Supporting Information (Table S1, \"Filtered PeptideAtlas build with human STY sites that have at least 1 associated PSM\"; Table S5, \"Summary of all STY sites in our analysis, their conservation data, proximal sites, phosphorylation likelihood, and structural data\"); they describe the 2020 PeptideAtlas build, not build 606, and I could not obtain them (the publisher and figshare pages returned HTTP 403; each fetch was stopped and none was retried). "
  "The 2023 JPR paper on multiple evidence streams could be read only as an abstract (a pipeline from PSM to peptidoform-site level with a decoy-based false localisation rate, shown on rice data); whether it offers site-level categories or data usable here is not determined.\n")
A("**Consequence.** The PeptideAtlas rows of G9 stay \"not feasible from the available build\". If you place Table S1 or S5 in `study/analysis/sources/`, a decision is needed first, because those files describe a different build (2020) from the one downloaded; nothing will be run on them without it.\n")
A("### PhosphoSitePlus measures (LT-supported vs MS-only; category high (n >= 5) vs low (n = 1), n = LT_LIT + MS_LIT + MS_CST)\n")
if PSP_PRESENT:
    A("The file is present; see the g9 results.\n")
else:
    A(f"**{PSP_MSG[0].upper() + PSP_MSG[1:]}.** `h_measures.py psp_lt psp_cat` would run both measures (D0, D1 and D4, gene-cluster and accession-cluster intervals) on the file, and stops if the columns differ from the plan's definition (D3). The set-B extract in `raw_data/` is not a substitute (RECONCILIATION.md section 7) and no surrogate was derived. "
      "**The meaning of `LT_LIT`, `MS_LIT` and `MS_CST` could not be confirmed from a primary page** (the phosphosite.org downloads page returned HTTP 403; the two PhosphoSitePlus papers I could read do not define the column names; `sources/SOURCES.md`, section C). The paper's PhosphoSitePlus rule is \"The sum of observations across all studies was used to predict site phosphorylation likelihood in PSP (\"High\": ≥5 observations; \"Medium\": 2−4 observations; \"Low\": 1 observation).\" It does not name the three columns, so whether n = LT_LIT + MS_LIT + MS_CST is that sum depends on what the download documents. "
      "**When the file arrives I will copy its header lines (licence, column notes), version stamp and any readme into `sources/SOURCES.md`, and I will not assume the meaning of the three columns beyond what the download documents; if it does not document them, Task H stops before it is run and I tell you** (instruction of 2026-10-04, added to decision D3's stop rule). "
      "The old `raw_data/harmonizome/phosphositeplus_gene_attribute_edges.txt.gz` is a gene-level attribute table, not this file, and was not used.\n")
A("### UniProt reviewed MOD_RES, experimental vs inferred (optional item of section 9)\n")
if "uniprot" in g9i:
    iu = g9i["uniprot"]
    A(f"{n_(iu['phospho_features'])} reviewed human phosphosite features (release 2026_03): {n_(iu['experimental_ECO0000269'])} with experimental evidence (ECO:0000269), {n_(iu['inferred_only'])} with inferred evidence only (ECO:0000250, ECO:0000255, ECO:0000305), "
      f"{n_(iu['neither (only ECO:0007744 or mixed or no tag)'])} neither (mass-spectrometry-only ECO:0007744, mixed or untagged; not compared). Matched to the site table on accession, residue and position: {n_(iu['site_records_matched'])} records ({n_(iu['high_class_records'])} experimental, {n_(iu['low_class_records'])} inferred) in {n_(iu['accessions'])} accessions and {n_(iu['gene_clusters'])} gene clusters; "
      f"{pct(iu['share_of_measure_sites_lost'])} of the compared UniProt sites have no record in the site table. Treated group = experimental, control group = inferred; D1 = CEM on the V50 index decile with residue type exact; D4 adds the PaxDb abundance decile (records with an abundance); gene-cluster bootstrap primary.\n")
    A("**Caveat, stated before the results.** The inferred class includes sites annotated by sequence-similarity transfer (ECO:0000250, and the other inferred codes), i.e. partly *because* they are conserved in an ortholog, so conservation is partly circular in this contrast. "
      "**This contrast is not used as evidence for or against the main hypothesis.** It is reported because it was a pre-defined optional item, and it is not compared with the tier results.\n")
    rows = []
    for s in ("S1", "S2", "S3", "S4"):
        m = man[(man["block"] == "G9") & (man["signal"] == s) & (man["variation"].str.startswith("UniProt"))]
        ids = dict(zip(m["design"], m["contrast_id"]))
        d0, d1, d4 = g9[ids["D0"]], g9[ids["D1"]], g9[ids["D4"]]
        if d1["status"] != "run":
            rows.append([SL[s], "not run: " + d1.get("reason", "")] + [""] * 8)
            continue
        g1_, a1 = d1["gene"], d1["protein"]
        rows.append([SL[s], n_(d1["n_high"]) + " / " + n_(d1["n_low"]), est_ci(d0["gene"]["estimate"], d0["gene"]["ci"], s), est_ci(g1_["adjusted"], g1_["adjusted_ci"], s), ci(a1["adjusted_ci"], s), pct(g1_["attenuation"]), yn(g1_["survives"]), yn(g1_["largely_explained"]),
                     est_ci(d4["gene"]["adjusted"], d4["gene"]["adjusted_ci"], s) if d4["status"] == "run" else "not run", f(g1_["smd"], "S2", 3)])
    A(md(["Signal", "Experimental / inferred records", "Crude (gene CI)", "D1 adjusted (gene CI)", "D1 accession CI", "Attenuation", "Survives (gene)", "Largely explained", "D4 adjusted (gene CI)", "SMD after D1"], rows))
    s1, s2, s3, s4 = (g9[[i for i in g9 if manid.loc[i, "signal"] == s and manid.loc[i, "design"] == "D1" and manid.loc[i, "variation"].startswith("UniProt")][0]] for s in ("S1", "S2", "S3", "S4"))
    s4c = g9[[i for i in g9 if manid.loc[i, "signal"] == "S4" and manid.loc[i, "design"] == "D0" and manid.loc[i, "variation"].startswith("UniProt")][0]]["gene"]
    A(f"**Reading (reported as found).** Domain overlap: the crude difference ({f(g9[[i for i in g9 if manid.loc[i, 'signal'] == 'S1' and manid.loc[i, 'design'] == 'D0' and manid.loc[i, 'variation'].startswith('UniProt')][0]]['gene']['estimate'], 'S1')}) disappears after matching on the index ({f(s1['gene']['adjusted'], 'S1')}, interval {ci(s1['gene']['adjusted_ci'], 'S1')} includes 1; attenuation {pct(s1['gene']['attenuation'])}), so this signal does not survive. "
      f"Conservation is lower for experimental than for inferred sites ({f(s2['gene']['adjusted'], 'S2')}). **This contrast is partly circular for conservation:** the inferred class is defined by evidence codes that include sequence-similarity transfer (ECO:0000250), i.e. sites annotated *because* they are conserved in an ortholog. "
      f"Disorder is lower in the experimental class ({f(s3['gene']['adjusted'], 'S3')}) with an attenuation of {pct(s3['gene']['attenuation'])} (largely explained: {yn(s3['gene']['largely_explained'])}). The sequence-based motif difference reverses sign after adjustment (crude {f(s4c['estimate'], 'S4')} with an interval ending at {s4c['ci'][1]:.5f}, adjusted {f(s4['gene']['adjusted'], 'S4')}); under the rule it counts as largely explained (attenuation {pct(s4['gene']['attenuation'])}), but the crude interval excludes zero only marginally, so this row is weak evidence either way. "
      "The population is small and its composition differs from the tier populations.\n")
else:
    A("Not run.\n")

# ---------------------------------------------------------------- 10 Task I
A("## 10. Within-protein selection (Task I: G10)\n")
if PSP_PRESENT:
    A("See the g10 results.\n")
else:
    A(f"**{PSP_MSG[0].upper() + PSP_MSG[1:]}.** It needs the PhosphoSitePlus LT-supported and MS-only classes. What exists: the conditional-logit estimator (`cond_logit.py`), checked against `statsmodels` ConditionalLogit on simulated data (`check_cond_logit.py`; the script asserts agreement of the coefficients within 5e-04 and of the standard errors within 1e-04, a log-likelihood at least as high as statsmodels', and no overflow in a 3,000-record stratum; the assertions pass); the same estimator, with a parallel likelihood, was checked against `statsmodels` on a real subsample for G21 (`results/g21_cond_logit_check.json`). "
      "What does not exist yet: the runner that joins the PhosphoSitePlus classes to the site table, builds the peptide-detectability indicator (tryptic peptide containing the site, cleaved after K or R but not before P, 7 to 35 residues; protein sequences are on disk in `data_private/dl/seq_*.tsv`) and runs the gene-cluster bootstrap with and without that control, plus the within-gene linear model for disorder. It is to be written if the file is supplied, so that it is tested on the real columns.\n")

# ---------------------------------------------------------------- 10b G21
A("## 10b. Within-gene comparison of GOLD and BRONZE records (block G21; amendment of 2026-10-04)\n")
g21 = J("g21_results.json") if (RES / "g21_results.json").exists() else {}
if not g21:
    A("Not run.\n")
else:
    chk = J("g21_cond_logit_check.json")
    A("**Definition (written into the plan and committed before any G21 contrast was run, commit 8807afd; decided after the interim results and the main results table were seen, not pre-specified).** "
      "Records: the GOLD and BRONZE records (set-B tiers) of the G1 populations (P-all, all species; P-human) with serine, threonine or tyrosine residue, outcome and V50 index defined; domain overlap excluded. Stratum: gene cluster (gene fixed effects). Exposure: GOLD against BRONZE. Covariates: residue type (T and Y indicators). "
      "Binary outcomes: conditional logistic regression (odds ratio of GOLD against BRONZE within a gene); disorder: within-gene linear fixed-effect model (difference in means); beside each binary outcome, as a same-scale companion to the CEM differences, the within-gene linear fixed-effect difference in proportions. "
      "**M1 (all records of the population, no detectability control) is the primary model.** M1s (records for which a UniProt reviewed protein sequence is available) and M2 (M1s plus the peptide-detectability indicator) are sensitivity analyses on the reviewed-sequence subset (confirmed by the project lead, 2026-10-04). The indicator is: tryptic peptide containing the site has 7 to 35 residues; cleavage after K or R but not before P. (tryptic peptide containing the site has 7 to 35 residues; cleavage after K or R but not before P). Gene-cluster bootstrap, percentile 95% interval, all models refitted in every resample.\n")
    A(f"**Check before the run (required by the definition).** `cond_logit.py` against `statsmodels` ConditionalLogit on a real subsample ({n_(chk['records'])} records in {n_(chk['gene_clusters'])} gene clusters of P-all, conservation; largest stratum {n_(chk['largest_stratum'])} records): largest coefficient difference {chk['max_abs_diff_coef_serial']:.1e} (tolerance {chk['tolerance_coef']:.0e}), largest standard-error difference {chk['max_abs_diff_se_serial']:.1e} (tolerance {chk['tolerance_se']:.0e}); "
      f"the parallel evaluation used in the bootstrap gives identical coefficients to the serial one (largest difference {chk['max_abs_diff_coef_serial_vs_parallel']:.1e}); the check {'passed' if chk['passes'] else 'FAILED'}. (`results/g21_cond_logit_check.json`.)\n")
    cemref = {}
    for x in main:
        cemref[(x["pop"], x["s"])] = x
    for popl, gl in (("P-all (S/T/Y)", "P-all"), ("P-human (S/T/Y)", "P-human")):
        rows = []
        info_ = None
        for s in ("S2", "S3", "S4", "S5"):
            for mk in ("M1", "M1s", "M2"):
                cid = [i for i in g21 if g21[i]["signal"] == s and g21[i]["population"] == popl and g21[i]["model"] == mk]
                if not cid:
                    continue
                v = g21[cid[0]]
                info_ = v["population_info"]
                if "companion_linear" in v:
                    cmp_ = v["companion_linear"]
                    comp = est_ci(cmp_["estimate"], cmp_["ci"], "S2")
                    est_s = est_ci(v["estimate"], v["ci"], "S1")
                else:
                    comp = "not applicable"
                    est_s = est_ci(v["estimate"], v["ci"], "S3")
                ref = ""
                if mk == "M1":
                    cx = cemref.get((gl, s))
                    if cx:
                        ref = est_ci(cx["gene"]["estimate"], cx["gene"]["ci"], s)
                rows.append([SL[s] + ("" if v["decision_bearing"] == "yes" else " (not decision-bearing)"), mk, n_(v["n_records"]), f"{n_(v['n_gene_clusters_with_both_tiers'])} ({n_(v['gold_in_gene_clusters_with_both_tiers'])} / {n_(v['bronze_in_gene_clusters_with_both_tiers'])})",
                             v["estimator"].split(",")[0].replace("conditional logistic regression", "odds ratio (conditional logit)").replace("within-gene linear fixed effect", "difference (within-gene FE)"), est_s, yn(v["persists"]), comp, ref])
        A(f"**{popl}**\n")
        A(md(["Signal", "Model", "Records", "Gene clusters with both tiers (GOLD / BRONZE records in them)", "Estimator", "Estimate (95% CI)", "Persists within gene (CI excludes null)", "Companion: within-gene linear difference (CI)", "CEM D1 difference, G1/G2 (gene CI), for scale"], rows))
        if info_:
            A(f"Records without a reviewed sequence: {n_(info_['records_without_sequence'])} of {n_(info_['records'])}; with a sequence but a different residue at the stated position: {n_(info_['records_residue_mismatch'])}; records in M1s and M2: {n_(info_['records_with_sequence_and_match'])}. "
              f"Records: {n_(info_['gold'])} GOLD and {n_(info_['bronze'])} BRONZE in {n_(info_['gene_clusters'])} gene clusters, {n_(info_['gene_clusters_with_both_tiers'])} of which have both tiers.\n")
    sc = J("g21_sequence_coverage.json")
    A("**Sequence coverage by tier** (reviewed UniProt sequences only; this decides who enters M1s and M2; `results/g21_sequence_coverage.json`):\n")
    A(md(["Signal", "Population", "GOLD records with a sequence", "BRONZE records with a sequence", "Detectable (7 to 35 residues) among GOLD with sequence", "Detectable among BRONZE with sequence"],
         [[SL[k.split("|")[0]], k.split("|")[1], pct(v["GOLD"]["share_with_sequence"]), pct(v["BRONZE"]["share_with_sequence"]), pct(v["GOLD"]["share_detectable_among_with_sequence"]), pct(v["BRONZE"]["share_detectable_among_with_sequence"])] for k, v in sc.items()]))
    lo = min(v["BRONZE"]["share_with_sequence"] for v in sc.values())
    A(f"The sequence-restricted models M1s and M2 therefore rest on a subset that keeps most GOLD records but only {pct(lo)} to {pct(max(v['BRONZE']['share_with_sequence'] for v in sc.values()))} of BRONZE records, depending on the outcome; their estimates should be compared with M1 with that selection in mind, and the comparison of M2 with M1s (same records, with and without the detectability control) is the clean contrast for the control.\n")
    dbm = [v for v in g21.values() if v["model"] == "M1" and v["decision_bearing"] == "yes"]
    d2m = [v for v in g21.values() if v["model"] == "M2" and v["decision_bearing"] == "yes"]
    unf = sum(v.get("n_boot_unfit", 0) for v in g21.values())
    A(f"Of the {len(dbm)} decision-bearing M1 contrasts, {sum(1 for v in dbm if v['persists'])} persist within gene; of the {len(d2m)} decision-bearing M2 contrasts (with the detectability control), {sum(1 for v in d2m if v['persists'])} persist. "
      f"{'All' if all(v['n_boot'] == N_BOOT for v in g21.values()) else 'Not all'} groups used {N_BOOT:,} resamples (the decision-D7 fallback to 200 {'was not' if not any(v['d7_fallback'] for v in g21.values()) else 'was'} applied); {unf} conditional-logit resample fits did not converge in total.\n")
    A("**In plain words:** a difference between GOLD and BRONZE sites that persists inside the same gene cannot be explained by anything that differs between genes, such as how heavily a gene is studied overall; **it excludes confounding at the protein (gene) level only.** "
      "It does not exclude confounding at the level of the site, such as how much attention a particular site received. Because sequences exist only for reviewed accessions, M1s and M2 describe a subset whose composition differs from M1 (see the shares of records with a sequence above and in `results/g21_sequence_coverage.json`).\n")

# ---------------------------------------------------------------- 10c G23
prim_s1 = am[src_of[d1_d0("G1", "S1")[1]]]
A("## 10c. Non-primary publications counted in depth (block G23; sensitivity only; amendment of 2026-10-04, sixth batch)\n")
g23t = J("g23_tiers_summary.json")
g23r = J("g23_results.json")
A(f"**Why.** PMID 25514926 is the PhosphoSitePlus database paper (PubMed: \"PhosphoSitePlus, 2014: mutations, PTMs and recalibrations\", Nucleic Acids Research 2015; RECONCILIATION.md section 7) and is counted as a publication in {n_([x['records_linked'] for x in g23t['nonprimary'] if x['pmid'] == '25514926'][0])} records. G23 asks how much the results depend on database, resource and review papers being counted as evidence. **The primary analysis is unchanged; G23 is a sensitivity analysis.** "
  f"Rule, written into the plan and the list committed before any depth or contrast was computed (commits 92ab560 and 5650bed): candidate PMIDs = the {n_(J('g23_pubmed_metadata.json')['n_candidates_study_cluster_rule'])} study-cluster candidates (linked to at least 884 records; PMID 25514926 is among them); PubMed metadata fetched by E-utilities; non-primary = PubMed publication type Review or Systematic Review, or a database or resource paper judged from the title alone (list with a reason for each entry in `nonprimary_pmid_list.md`; {len(J('g23_nonprimary_pmids.json')['not_removed_borderline'])} borderline titles are listed there as not removed). "
  "Non-primary PMIDs are deleted from every publication list, the study-cluster assignment of the remaining PMIDs is kept, depth and tier are recomputed (DAS unchanged), and records left with no PMID and no database source flag are dropped.\n")
A("### Descriptive\n")
rows = [[x["pmid"], x["title"], n_(x["records_linked"]), ("Review type" if x["criterion_a_review_type"] else "") + (" + " if x["criterion_a_review_type"] and x["criterion_b_title"] else "") + ("title" if x["criterion_b_title"] else "")] for x in sorted(g23t["nonprimary"], key=lambda x: -x["records_linked"])]
A(md(["PMID", "Title (PubMed)", "Records whose list contains it", "Criterion"], rows))
A(f"**{g23t['n_nonprimary_pmids']} non-primary PMIDs** (of {n_(J('g23_pubmed_metadata.json')['n_candidates_study_cluster_rule'])} candidates; {len(J('g23_nonprimary_pmids.json')['criterion_a_review_type'])} by publication type, {len(J('g23_nonprimary_pmids.json')['criterion_b_title'])} by title; one is both). "
  f"**{n_(g23t['records_with_at_least_one'])} of {n_(g23t['n_records'])} records ({pct(g23t['share_records_with_at_least_one'])}) have at least one in their publication list** (by species: {', '.join(f'{k} {n_(v)}' for k, v in g23t['by_species'].items())}; by tier before correction: {', '.join(f'{k} {n_(v)} ({pct(g23t['share_of_each_old_tier'][k])} of that tier)' for k, v in g23t['by_old_tier'].items())}). "
  f"{n_(g23t['records_whose_only_pmids_are_nonprimary'])} records have only non-primary PMIDs; {n_(g23t['records_dropped_no_pmid_no_source_flag'])} of them have no database source flag either and are dropped.\n")
kp_ = J("g23_kept_nonevidence_pmids.json")
A("**Kept as primary although the title is not phosphorylation evidence (descriptive only; the G23 rule and every tier are unchanged):** "
  + "; ".join(f"{p} (\"{kp_[p]['title']}\") is in the publication lists of {n_(kp_[p]['records'])} records ({100 * kp_[p]['share_of_records']:.2f}%; tiers {', '.join(f'{k} {n_(v)}' for k, v in kp_[p]['by_tier'].items())}; species {', '.join(f'{k} {n_(v)}' for k, v in kp_[p]['by_species'].items())}; only PMID in the list for {n_(kp_[p]['only_pmid_in_list'])})" for p in ("30243273", "37031204", "14702039"))
  + f". Together {n_(kp_['any_of_the_three']['records'])} records ({100 * kp_['any_of_the_three']['share_of_records']:.2f}%) carry at least one. The first has a title on a quantum-chemistry method, which suggests an erroneous PMID in the source; the second is a correction notice that double-counts its original paper; the third is a cDNA sequencing project.\n")
A("### Tier transition under the corrected depth (old tier by new tier)\n")
tt = g23t["transition_table"]
cols_ = ["GOLD", "SILVER", "BRONZE", "dropped"]
A(md(["Old tier \\ new tier"] + cols_ + ["Total"], [[r] + [n_(tt[r].get(c, 0)) for c in cols_] + [n_(sum(tt[r].values()))] for r in ("GOLD", "SILVER", "BRONZE")]))
A(f"GOLD falls from {n_(g23t['tier_counts_old']['GOLD'])} to {n_(g23t['tier_counts_new']['GOLD'])}; {n_(g23t['records_whose_tier_changes_or_dropped'])} records change tier or are dropped ({pct(g23t['records_whose_tier_changes_or_dropped'] / g23t['n_records'])}).\n")
A("### G1 and G2 contrasts under the corrected tiers, beside the primary\n")
A("Same design as the primary (CEM on the V50 decile, gene-cluster bootstrap, 1,000 resamples; D0 and D1 from the same resamples); only the tier column differs. Decision by the pre-specified rule: survives = the gene-cluster interval of the adjusted estimate excludes the null; sign = same side of the null as the primary estimate; largely explained = attenuation of 50% or more with a crude interval excluding the null.\n")
rows = []
nsurv = nsign = nle = ntot = 0
for cid, v in sorted(g23r.items()):
    if v["design"] != "D1":
        continue
    mrow = manid.loc[cid]
    srcid = mrow["note"].split("counterpart of ")[1].split(";")[0]
    s = v["signal"]
    pr = am[src_of[srcid]]
    same = (v["estimate"] > NULL(s)) == (pr["estimate"] > NULL(s))
    d0v = g23r[[k for k in g23r if manid.loc[k, "signal"] == s and manid.loc[k, "population"] == mrow["population"] and manid.loc[k, "design"] == "D0"][0]]
    ntot += 1
    nsurv += int(excl(v["ci"], s))
    nsign += int(same)
    nle += int(bool(v.get("largely_explained")))
    rows.append([("G1" if mrow["variation"].endswith("G1") else "G2"), SL[s] + ("" if mrow["decision_bearing"] == "yes" else " (not decision-bearing)"), mrow["population"], est_ci(pr["estimate"], pr["ci"], s), est_ci(v["estimate"], v["ci"], s),
                 f"{n_(v['n_sites_gold'])} / {n_(v['n_sites_bronze'])}", yn(excl(v["ci"], s)), "same" if same else "OPPOSITE", pct(v["attenuation"]), yn(bool(v.get("largely_explained"))), f"{cid} (counterpart {srcid})"])
A(md(["Block", "Signal", "Population", "Primary D1 (gene CI; G17)", "Corrected tiers D1 (gene CI)", "GOLD / BRONZE (corrected)", "Survives", "Sign vs primary", "Attenuation", "Largely explained", "Contrast"], rows))
A(f"Under the corrected tiers {nsurv} of {ntot} D1 contrasts survive, {nsign} of {ntot} have the same sign as the primary, and {nle} are largely explained; " + ("the decisions are the same as in the primary analysis for every contrast" if (nsurv == ntot and nsign == ntot and nle == 0) else "at least one decision differs from the primary analysis (see the table)") + f". The estimates move modestly (for example domain overlap {f(prim_s1['estimate'], 'S1')} to {f(g23r['G23-0518']['estimate'], 'S1')}, disorder {f(am['G17-0384']['estimate'], 'S3')} to {f(g23r['G23-0522']['estimate'], 'S3')}).\n")
A("### GOLD counts in the held-out designs under the corrected tiers\n")
ho = g23t["held_out"]
nm = {"cluster_disjoint_test": "Cluster-disjoint test set", "cluster_and_protein_disjoint_test": "Cluster- and protein-disjoint test set", "temporal_test_all": "Temporal, all test sites",
      "temporal_protein_disjoint": "Temporal, protein-disjoint", "temporal_study_disjoint": "Temporal, study-disjoint", "temporal_protein_and_study_disjoint": "Temporal, protein- and study-disjoint"}
A(md(["Design", "Test sites", "GOLD (primary tiers)", "GOLD (corrected)", "SILVER (corrected)", "BRONZE (corrected)", "Dropped"],
     [[nm[k], n_(ho[k]["n"]), n_(ho[k]["old"].get("GOLD", 0)), n_(ho[k]["new"].get("GOLD", 0)), n_(ho[k]["new"].get("SILVER", 0)), n_(ho[k]["new"].get("BRONZE", 0)), n_(ho[k]["new"].get("dropped", 0))] for k in nm]
     + [["All sites", n_(g23t["n_records"]), n_(ho["gold_total"]["old"]), n_(ho["gold_total"]["new"]), n_(g23t["tier_counts_new"]["SILVER"]), n_(g23t["tier_counts_new"]["BRONZE"]), n_(g23t["tier_counts_new"]["dropped"])]]))
A("The Rev 5 split assignments (cluster-disjoint, protein-disjoint, temporal) are unchanged; only the tier column is replaced. The held-out designs still retain few GOLD sites (the protein-disjoint test set keeps " + f"{n_(ho['cluster_and_protein_disjoint_test']['new'].get('GOLD', 0))} of {n_(ho['gold_total']['new'])}).\n")

# ---------------------------------------------------------------- 11 Task J
A("## 11. Predictor gaps (Task J: G11, secondary analysis)\n")
A("Top-three minus bottom-three index-decile AUROC gap on the frozen evaluation sample, planned design (DeepPhos and MusiteDeep on serine/threonine sites; PhosF3C on tyrosine sites). Equivalence margin +/-0.05 (90% interval inside the margin = gap absent); gene-cluster resampling within the decile pools (primary) with the Rev 5 site-level interval beside it. Restricted analyses keep positives first reported after the tool's restriction year and need at least 200 positives. No pooled ranking across tools. This analysis is secondary.\n")
rows = []
for k in g11:
    if "planned" not in k:
        continue
    tool = k.split("|")[0]
    for kk, v in g11[k].items():
        if "gap" not in v:
            rows.append([tool, kk, "", n_(v.get("n_positives", 0)), "not evaluated: " + v.get("status", ""), "", "", "", ""])
            continue
        rows.append([tool, kk, n_(v["n_sites"]), n_(v["n_positives"]), f"{f(v['bottom_auroc'], 'S2')} / {f(v['top_auroc'], 'S2')}", f(v["gap"], "S2"), ci(v["gene_ci95"], "S2"), ci(v["gene_ci90"], "S2"), v["gene_label"] + "; site-level interval: " + v["site_label"]])
A(md(["Tool", "Setting", "Sites", "Positives", "AUROC bottom / top", "Gap", "Gene-cluster 95% CI", "Gene-cluster 90% CI", "Label (gene-cluster interval; site-level interval)"], rows))
ov = g11["J4_training_overlap"]["PhosF3C"]
A(f"**Training-set overlap (J4).** Only PhosF3C ships its training file; {n_(ov['n_flagged_ge_40pct'])} of {n_(ov['n_eval_positives_Y'])} tyrosine evaluation positives ({pct(ov['share_flagged'])}) have a 31-residue window with at least 40% ungapped identity to a training positive. "
  f"The DeepPhos and MusiteDeep repositories hold no training data, so the item stops for them. Details in `results/g11_predictors.json`.\n")
agree = [(tool, kk) for tool in ("DeepPhos", "MusiteDeep", "PhosF3C") for kk, v in g11[[k for k in g11 if k.startswith(tool + "|planned")][0]].items() if isinstance(v, dict) and "gene_label" in v]
same = sum(1 for tool, kk in agree if g11[[k for k in g11 if k.startswith(tool + "|planned")][0]][kk]["gene_label"] == g11[[k for k in g11 if k.startswith(tool + "|Rev 5")][0]][kk]["gene_label"])
A(f"The Rev 5 all-residue design is kept beside each planned row in `results/g11_predictors.json` (and the Rev 5 reported intervals in the `rev5_reported` entries); the gene-cluster label of the planned design equals that of the Rev 5 all-residue design in {same} of {len(agree)} settings.\n")
qual = []
for tool in ("DeepPhos", "MusiteDeep", "PhosF3C"):
    for kk, v in g11[[k for k in g11 if k.startswith(tool + "|planned")][0]].items():
        if isinstance(v, dict) and "gap" in v and v["gap"] >= 0.05 and (v["gene_ci95"][0] > 0 or v["gene_ci95"][1] < 0):
            qual.append((tool, kk, v["gap"], v["gene_ci95"]))
qual_tools = sorted({q[0] for q in qual})
A(f"**The criterion written down before any predictor was scored (commit 8e12390, 2026-09-24; recorded in `manuscript/CLAIMS_VS_EVIDENCE.md`) was a gap of at least 0.05 with an interval excluding zero in at least two tools. It is still not met:** only {len(qual)} tool-and-setting"
  f" combination{'s' if len(qual) != 1 else ''} qualif{'y' if len(qual) != 1 else 'ies'}: " + "; ".join(f"{t} {s} (gap {f(g_, 'S2')}; gene-cluster 95% interval {ci(c_, 'S2')})" for t, s, g_, c_ in qual) + f", that is {len(qual_tools)} of 3 tools (criterion: at least 2). This is the criterion set before the results; the 90% equivalence test and the labels above are the additional rules of this plan.\n")

# ---------------------------------------------------------------- 12 departures
A("## 12. Departures, limitations and what I need from you\n")
A("**Logged in the plan (section 15), none changes a result:** all entries of the interim stage; the amendments of 2026-10-01 (decisions D5 and D2, the gene-cluster correction, the added columns, G16 to G20, PeptideAtlas, Task H categories, Task I model); the operational definitions fixed while carrying out the amendment. "
  "Two further items arose in this stage and are logged with this report: (a) the review of the generated reports found a typed number in the interim builder (now computed); (b) the optional UniProt measure was run before the PhosphoSitePlus measures because it does not depend on the file. "
  "The fourth batch of 2026-10-04 added: the definitions read from the paper (PeptideAtlas rows stay not feasible), the definitions of G21 (within-gene GOLD versus BRONZE) and G22 (S/T-only motif), written into the plan and committed before they were run, and the report changes you asked for.\n")
A("**Limitations to carry forward.** (1) The two-way interval rests on a few very large study clusters and is approximate. (2) The index is a gene-level count and is only weakly associated with the outcomes within a tier, so CEM on it cannot adjudicate site-level effort confounding (section 4). "
  "(3) Domain overlap is population-dependent and consistent with a composition effect (section 3), which was not tested further. (4) The D2 non-fits are results of the pre-specified design; D2b is post-hoc. (5) UniProt inferred evidence is partly conservation-derived. (6) Multiplicity is not corrected (plan section 2); the decision-bearing contrasts of G1 to G15 were fixed in the manifest before any contrast was run (commit f199f81, 2026-09-29), and those of the amendment blocks G16 to G20 were fixed in manifest rows committed before they were run (commit 0407605, 2026-10-01), but after the interim results had been seen. "
  "(7) The study-cluster threshold has no effect on any tier (G13). "
  "(8) In the predictor audit the unlabelled negatives on little-studied proteins may include undiscovered phosphosites, which lowers the AUROC there independently of the tool; the gap between well-studied and little-studied deciles is therefore not a clean measure of the tool's performance difference. "
  "(9) The conservation indicator requires the ortholog site to be reported as a phosphosite in the other species, so it depends on how much that species has been studied; a within-gene comparison removes differences between genes but not this dependence, because it is a property of the ortholog sites.\n")
A("**Needed from you (rule 9):**\n")
if not PSP_PRESENT:
    A(f"1. **PhosphoSitePlus.** Nothing is pending from you: Tasks H and I and G20-0488 are {PSP_MSG} and are in Future work.")
A("2. Only if you place Table S1 or S5 of the Kalyuzhnyy paper in `sources/`: a decision on whether files that describe the 2020 PeptideAtlas build may stand in for build 606. The manuscript, references, `legacy/` and `Figures.pptx` are untouched.\n")

# ---------------------------------------------------------------- 13 manifest
A("## 13. Status of the manifest\n")
A(f"{len(man)} rows; every row of `grid_manifest.csv` is in `results/manifest_status.csv` with its status. G7's duplicate removal and G13's identical results count as run.\n")
A(cnt.reset_index().rename(columns={"block": "Block"}).to_markdown(index=False))
A("")
g5nr = ms[(ms["block"] == "G5") & (ms["result_status"] == "not run")]
A(f"**The {len(g5nr)} G5 rows marked \"not run\"** ({', '.join(g5nr['contrast_id'])}) are the {', '.join(sorted(set(g5nr['signal_label'])))} contrasts for {', '.join(sorted(set(g5nr['population'])))} under {', '.join(sorted(set(g5nr['design'])))}. "
  f"Reason recorded in the results file: \"{'; '.join(sorted(set(g5nr['result_note'])))}\" (the domain annotation exists only for the human records). The same signal was run for human (P-abund, human), and every other signal was run for mouse and rat.\n")
# ---------------------------------------------------------------- 14 headline rule, 15 future work
A("## 14. Pre-set headline rule\n")
A("Branch A if Task I shows a within-gene difference between LT-supported and MS-only sites in conservation or disorder (detectability-controlled model, gene-cluster interval excluding zero); Branch B otherwise, or if Task I is not run.\n")
g10r = J("g10_results.json") if (RES / "g10_results.json").exists() else {}
if g10r:
    hits = [v for v in g10r.values() if v.get("signal") in ("S2", "S3") and v.get("model") == "controlled" and v.get("status") == "run" and (v["ci"][0] > 0 or v["ci"][1] < 0)]
    A(f"**Branch {'A' if hits else 'B'} applies** (Task I results in `results/g10_results.json`).\n")
else:
    A(f"**Branch B applies**, because Task I has not been run ({PSP_MSG}). If Task I is run later, the branch is re-evaluated by `build_claims.py` with the same rule.\n")
A("## 15. Future work (ideas outside the frozen scope; not part of any analysis above)\n")
for t in ("Tasks H and I and G20-0488 (PhosphoSitePlus evidence measures; within-gene comparison of LT-supported and MS-only sites, with and without the peptide-detectability control; index against LT_LIT-supported sites per gene), to be run if the provider supplies `Phosphorylation_site_dataset.gz` (possibly during revision); the file's own header and documentation of `LT_LIT`, `MS_LIT` and `MS_CST` are to be recorded first, and the plan's stop rule applies if they are not documented.",
          "Adjusting for study effort at the level of the site (for example the number of publications per site that is not part of the tier definition), which the gene-level index cannot do.",
          "A conservation indicator that accounts for how much the other species has been studied.",
          "A detectability control for every accession (sequences for unreviewed accessions), so that the sequence-restricted models are not limited to reviewed proteins.",
          "A site-level PeptideAtlas evidence measure from per-site PTMProphet probabilities, if a per-site file for the current build becomes available (or a decision is taken on using the authors' 2020-build tables).",
          "Replication of the tier comparisons in a second substrate, which the single-source subsets did not allow.",
          "Reprocessing public phosphoproteomics data at the PSM level with a decoy-based false localisation rate, so that tiers rest on localisation confidence."):
    A("- " + t)
A("")
HL[HL.index("TIER_HEADLINE")] = (f"5. **The tier definition changes the size of the signals more than adjustment for the index does:** like-for-like on the log scale, the spread of the adjusted estimates across tier definitions exceeds the largest change caused by effort adjustment within a definition in {tier_bigger} of {len(lfd)} signal-by-population combinations "
                                 f"(natural-log scale; spread {min(t[3] for t in lfd):.3f} to {max(t[3] for t in lfd):.3f}, as plain ratios {math.exp(min(t[3] for t in lfd)):.2f} to {math.exp(max(t[3] for t in lfd)):.2f}; within-definition change {min(t[2] for t in lfd):.3f} to {max(t[2] for t in lfd):.3f}, as plain ratios {math.exp(min(t[2] for t in lfd)):.3f} to {math.exp(max(t[2] for t in lfd)):.3f}; section 6). "
                                 f"The {max(sp_all) / min(sp_all):.0f}-fold spread for disorder in P-all rests on a near-zero agreement-only estimate; without it the spread is about {max(sp_all_ex) / min(sp_all_ex):.0f}-fold ({max(sp_all_ex) / min(sp_all_ex):.2f}).")
L[HEAD_IDX:HEAD_IDX] = [h + "\n" for h in HL]
open(REV6 / "RESULTS_REPORT.md", "w", encoding="utf-8").write("\n".join(L))
print("written", len(L), "blocks;", "unrendered:", bad)
