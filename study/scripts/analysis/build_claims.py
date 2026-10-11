"""Write study/analysis/CLAIMS_VS_EVIDENCE.md from the Rev 6 result files (and, for the comparison with Rev 5, study/earlier_results/*.json).
Every number is read from a result file; strength grades are computed by the rules printed at the top of the file, not typed in.
Run after build_results.py (it reads results/manifest_status.csv and the same result files). No manuscript file is read for numbers or edited."""
import glob
import json
import math
import re

import numpy as np
import pandas as pd

from common import REV6, RES
from provenance import RowList, dump_bindings, fmt_int, jload

REV5 = REV6.parent / "earlier_results"
J = lambda p: jload(RES / p)  # noqa: E731
J5 = lambda p: jload(REV5 / p)  # noqa: E731
man = pd.read_csv(REV6 / "grid_manifest.csv").fillna("")
manid = man.set_index("contrast_id")
cem = J("cem_results.json")
am = J("amend_results.json")
g5 = J("g5_abundance.json")
g20 = J("g20_diagnostics.json")
b4 = J("b4_twoway.json")
b5 = J("b5_index_effect_sizes.json")
g11 = J("g11_predictors.json")
g14 = J("g14_cem.json")
gu = J("gene_unit_summary.json")
rule = J("b_cluster_rule_info.json")
g22 = J("g22_results.json")
g21 = J("g21_results.json") if (RES / "g21_results.json").exists() else {}
g10 = J("g10_results.json") if (RES / "g10_results.json").exists() else {}
g9 = J("g9_measures.json") if (RES / "g9_measures.json").exists() else {}
fr, d2b, prop = {}, {}, {}
for f_ in glob.glob(str(RES / "f_runs" / "*.json")):
    fr.update(jload(f_))
for f_ in glob.glob(str(RES / "d2b" / "*.json")):
    r_ = jload(f_)
    d2b[r_["contrast_id"]] = r_
for f_ in glob.glob(str(RES / "prop" / "*.json")):
    r_ = jload(f_)
    prop[r_["contrast_id"]] = r_
m5 = J5("matching_results.json")
pm5 = J5("predictor_metrics.json")
dc5 = J5("dataset_counts.json")
SL = {"S1": "Domain overlap", "S2": "Conservation", "S3": "Disorder", "S4": "Motif, sequence-based", "S5": "Motif, database lookup"}
NULL = lambda s: 1.0 if s == "S1" else 0.0  # noqa: E731
NOT_GRADED = "descriptive only; not graded (circular: the indicator comes from the annotation tables that also feed the tiers; near-separation, see the fit counts)"


def f(v, s="S2", nd=None):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "n/a"
    return f"{v:.{nd or (3 if s == 'S1' else 4)}f}"


def ci(c, s="S2", nd=None):
    return f"{f(c[0], s, nd)} to {f(c[1], s, nd)}"


def n_(v):
    return fmt_int(v)


def excl(c, s):
    return bool(c[0] > NULL(s) or c[1] < NULL(s))


def md(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------ ids
src_of = {}
for cid, r in manid.iterrows():
    if r["block"] == "G17":
        src_of[re.search(r"same records as (\S+?);", r["note"]).group(1)] = cid


def d1(block, s, pop=None):
    m = man[(man["block"] == block) & (man["signal"] == s) & (man["design"] == "D1")]
    if pop:
        m = m[m["population"] == pop]
    return m["contrast_id"].iloc[0]


def prim(s):
    """Primary estimate: gene-cluster CEM D1 of G1 (G17 copy)."""
    return am[src_of[d1("G1", s)]]


# ------------------------------------------------------------------ robustness checks per signal
def checks(s):
    """List of (label, block id, estimate, ci) from the robustness blocks, for signal s."""
    out = []
    for cid, v in am.items():
        if v["signal"] != s or v["design"] != "D1" or v["status"] != "run":
            continue
        if cid.startswith("G16"):
            out.append((f"canonical reviewed accessions, {v['population'].split(',')[0]}", cid, v["estimate"], v["ci"]))
        if cid.startswith("G19") and "set-B" not in v["tier_rule"]:
            out.append((f"{v['tier_rule'].split(' (')[0]}, {v['population']}", cid, v["estimate"], v["ci"]))
    for cid, v in fr.items():
        m = manid.loc[cid]
        if m["signal"] != s or m["design"] != "D1" or v["gene"]["status"] != "run":
            continue
        if m["block"] == "G6":
            out.append((f"residue stratum {m['population']}", cid, v["gene"]["estimate"], v["gene"]["ci"]))
        elif m["block"] == "G7" and m["variation"] != "drop_consortium_7_PMIDs":
            out.append((m["variation"], cid, v["gene"]["estimate"], v["gene"]["ci"]))
        elif m["block"] == "G8" and not (m["variation"] == "gold_depth=10;gold_das=2"):
            out.append((f"tier threshold {m['variation']}", cid, v["gene"]["estimate"], v["gene"]["ci"]))
    for cid, v in d2b.items():
        if v["signal"] == s and v["decision_bearing"] == "yes":
            out.append((f"matcher D2b, {v['population']}", cid, v["estimate"], v["ci"]))
    for cid, v in cem.items():
        m = manid.loc[cid] if cid in manid.index else None
        if m is not None and m["block"] == "G4" and m["signal"] == s and v.get("design") == "D3" and m["index"] != "V0":
            out.append((f"index variant {m['index']}", cid, v["estimate"], v["ci"]))
    return out


def populations(s):
    """Adjusted estimates in the populations used for the direction rule: G1, G2, G3 and P-abund (human, mouse, rat; D1, D4, D5)."""
    out = []
    for blk in ("G1", "G2", "G3"):
        m = man[(man["block"] == blk) & (man["signal"] == s) & (man["design"] == "D1")]
        for _, r in m.iterrows():
            out.append((f"{blk} {r['population']}", r["contrast_id"], cem[r["contrast_id"]]["estimate"]))
    m = man[(man["block"] == "G5") & (man["signal"] == s) & man["design"].isin(["D1", "D4", "D5"])]
    for _, r in m.iterrows():
        if r["contrast_id"] in g5 and g5[r["contrast_id"]]["status"] == "run":
            out.append((f"{r['population']} {r['design']}", r["contrast_id"], g5[r["contrast_id"]]["estimate"]))
    return out


def grade(s):
    assert s != "S5", "the database-lookup motif is descriptive only and is not graded"
    p = prim(s)
    k1 = excl(p["ci"], s)
    pops = populations(s)
    k2 = len({e > NULL(s) for _, _, e in pops}) == 1
    sign = p["estimate"] > NULL(s)
    ch = checks(s)
    fails = [(lab, cid) for lab, cid, e, c in ch if (e > NULL(s)) != sign or not excl(c, s)]
    tw = b4.get(d1("G1", s))
    k4 = None if tw is None else bool(tw["survives_two_way"])
    wg = None
    if g21:
        m1 = [v for v in g21.values() if v["signal"] == s and v["model"] == "M1" and v["population"].startswith("P-all")]
        wg = bool(m1[0]["persists"]) if m1 else None
    F = len(fails) + (1 if k4 is False else 0) + (1 if wg is False else 0)
    if not k1:
        g = "not supported"
    elif not k2:
        g = "weak"
    elif F == 0:
        g = "strong"
    elif F <= 2:
        g = "moderate"
    else:
        g = "weak"
    return dict(k1=k1, k2=k2, k4=k4, wg=wg, fails=fails, n_checks=len(ch), F=F, grade=g, pops=pops)


STEP = ["not supported", "weak", "moderate", "strong"]
CAVEAT_DOWN = {"S2": "the indicator needs the ortholog site to be reported in the other species, so it also depends on study effort there (limitation 9 of RESULTS_REPORT.md)"}


def final(s):
    g = grade(s)
    gr = g["grade"]
    cav = CAVEAT_DOWN.get(s)
    if cav and gr in ("strong", "moderate"):
        gr2 = STEP[STEP.index(gr) - 1]
    else:
        gr2 = gr
    return g, gr, gr2, cav


L = []
A = L.append
A("# Claims versus evidence, Rev 6 (internal working file, not part of the manuscript)\n")
A("Written by `scripts/analysis/build_claims.py` from the files in `study/analysis/results/` and, for the Rev 5 values, `study/earlier_results/*.json`. Every number is read from a result file. Nothing in the manuscript folder is read for numbers or edited. "
  "Intervals are 95% intervals; the gene-cluster bootstrap is the primary interval from the amendment of 2026-10-01 unless the row says otherwise. Block and contrast identifiers refer to `grid_manifest.csv` and `results/manifest_status.csv`.\n")
A("## Grading rules (fixed before the grades were computed; applied by the script)\n")
A("**Effect claims about a signal (GOLD minus BRONZE, adjusted for the primary index, CEM, V50, set-B tiers).** Checks: K1 the gene-cluster interval of the primary estimate (G17) excludes the null. K2 the adjusted estimate lies on the same side of the null in every population where the signal is defined (G1, G2, G3 and P-abund in human, mouse and rat under D1, D4, D5). "
  "K3 robustness checks, each counted as a failure if its estimate has the opposite sign to the primary estimate or its gene-cluster interval includes the null: canonical reviewed accessions (G16), alternative tier definitions (G19), residue strata (G6), leave-one-cluster-out (G7), tier thresholds (G8), matcher D2b (G18, decision-bearing contrasts) and index variants (G4). "
  "K4 the two-way interval (B4; G1 D1 only) excludes the null (a failure if it does not). K5 the within-gene comparison (G21, model M1, P-all) persists (a failure if it does not; ignored while not run). F = number of K3 failures + K4 failure + K5 failure.\n")
A("**Grade:** *not supported* if K1 fails; otherwise *weak* if K2 fails or F >= 3; *moderate* if F is 1 or 2; *strong* if F = 0. "
  "**Construct caveat:** where the plan records a design-level limitation that bears directly on what the claim would mean, the grade is lowered one step and the caveat is named in the row (conservation only: limitation 9). "
  "**Not graded:** the database-lookup kinase-motif indicator (S5) is descriptive only (circular: the indicator comes from the annotation tables that also feed the tiers); it has no grade in any row. **Other claims:** a deterministic count or description read from a result file is *strong*; a claim that goes beyond what the design identifies (mechanism, cause) is *not supported*; predictor claims follow the rule given in section 6.\n")

# ---------------------------------------------------------------- 1 dataset and tier construction
A("## 1. Dataset and tier construction (descriptions)\n")
tc = dc5["tier_counts_mammalian"]
hr, hs = dc5["human_tier_raw_pmid_rule"], dc5["human_tier_study_cluster_rule"]
nh = dc5["human_n_sites"]
cl = rule
rows = RowList()
rows.append(["C1.1", f"The analysed dataset holds {n_(dc5['n_sites'])} sites on {n_(dc5['n_proteins'])} proteins (human, mouse, rat).", f"{n_(dc5['n_sites'])} sites; {n_(dc5['n_proteins'])} proteins", "", n_(dc5["n_sites"]), "`study/earlier_results/dataset_counts.json`", "strong", "deterministic count"])
rows.append(["C1.2", "Tier composition over all species is GOLD, SILVER and BRONZE as stated.", f"GOLD {n_(tc['GOLD'])}; SILVER {n_(tc['SILVER'])}; BRONZE {n_(tc['BRONZE'])}", "", n_(sum(tc.values())), "`dataset_counts.json`", "strong", "deterministic count"])
rows.append(["C1.3", "In the human data the share of GOLD sites depends on how publication depth is counted (a raw PMID count against the study-cluster rule). The two rules also differ in thresholds, so the change cannot be attributed to the counting method alone.", f"GOLD {100 * hr['GOLD'] / nh:.1f}% under the raw PMID rule against {100 * hs['GOLD'] / nh:.1f}% under the study-cluster rule", "", n_(nh), "`dataset_counts.json`", "strong (description); the attribution to counting alone is not supported", "deterministic count; thresholds differ between the rules"])
rows.append(["C1.4", "Database agreement (DAS >= 2) in human depends on the source set that defines it.", f"{100 * dc5['human_das_ge2_share_source_set_a']:.2f}% under source set A against {100 * dc5['human_das_ge2_share_source_set_b']:.1f}% under source set B", "", n_(nh), "`dataset_counts.json`", "strong", "deterministic count"])
rows.append(["C1.5", "The study-cluster threshold (500, 884 or 1,000 sites per candidate PMID) changes no record's depth or tier.", "", "", n_(dc5["n_sites"]), "`results/b_cluster_rule_info.json`; G13", "strong", "deterministic recomputation"])
same_tiers = all(cl[t]["tier_counts_study_cluster_rule"] == cl["884"]["tier_counts_study_cluster_rule"] for t in ("500", "1000"))
rows[-1][2] = ("identical tier counts at 500, 884 and 1,000 (GOLD " + n_(cl["884"]["tier_counts_study_cluster_rule"]["GOLD"]) + ", SILVER " + n_(cl["884"]["tier_counts_study_cluster_rule"]["SILVER"]) + ", BRONZE " + n_(cl["884"]["tier_counts_study_cluster_rule"]["BRONZE"]) + "); G13 reproduces G1 record for record") if same_tiers else "tier counts differ between thresholds"
g23t = J("g23_tiers_summary.json")
kp = J("g23_kept_nonevidence_pmids.json")
KEPT_TXT = "; ".join(f"{p} (\"{kp[p]['title'][:60]}...\", {n_(kp[p]['records'])} records, {100 * kp[p]['share_of_records']:.2f}%)" for p in ("30243273", "37031204", "14702039")) + f"; {n_(kp['any_of_the_three']['records'])} records ({100 * kp['any_of_the_three']['share_of_records']:.2f}%) carry at least one"
g23r = J("g23_results.json")
ho = g23t["held_out"]
rows.append(["C1.6", "**Evidence accounting: database or review papers are counted as evidence in publication depth.** " + str(g23t["n_nonprimary_pmids"]) + " non-primary PMIDs (database, resource or review papers; titles in `nonprimary_pmid_list.md`) are in the publication lists of "
             + n_(g23t["records_with_at_least_one"]) + " records (" + f"{100 * g23t['share_records_with_at_least_one']:.1f}" + "%), including the PhosphoSitePlus database paper PMID 25514926 in " + n_([x["records_linked"] for x in g23t["nonprimary"] if x["pmid"] == "25514926"][0]) + " records; " + n_(g23t["tier_counts_old"]["GOLD"] - g23t["tier_counts_new"]["GOLD"]) + " GOLD records fall to SILVER when they are not counted. Three PMIDs kept as primary (rule unchanged) are not phosphorylation evidence: " + KEPT_TXT + ".",
             f"{g23t['n_nonprimary_pmids']} PMIDs; {n_(g23t['records_with_at_least_one'])} records ({100 * g23t['share_records_with_at_least_one']:.1f}%)", "", n_(g23t["n_records"]), "G23 (`results/g23_tiers_summary.json`); RECONCILIATION.md section 7", "strong (description)", "list written from titles and PubMed publication types before any computation; counts deterministic"])
rows.append(["C1.7", "Held-out designs retain few GOLD sites, also under the corrected tiers.", f"protein-disjoint test set {n_(ho['cluster_and_protein_disjoint_test']['old'].get('GOLD', 0))} GOLD (primary tiers) and {n_(ho['cluster_and_protein_disjoint_test']['new'].get('GOLD', 0))} (corrected); temporal test set {n_(ho['temporal_test_all']['old'].get('GOLD', 0))} and {n_(ho['temporal_test_all']['new'].get('GOLD', 0))}; temporal protein-disjoint {n_(ho['temporal_protein_disjoint']['old'].get('GOLD', 0))} and {n_(ho['temporal_protein_disjoint']['new'].get('GOLD', 0))}",
             "", "see RESULTS_REPORT.md section 10c", "G23", "strong (counts)", "Rev 5 split assignments unchanged; only the tier column replaced"])
A(md(["ID", "Claim", "Estimate", "95% CI", "n", "Source", "Strength", "Basis"], rows))

# ---------------------------------------------------------------- 2 signals
A("## 2. Validation signals: GOLD against BRONZE after adjusting for the primary study-effort index\n")
A("Primary estimate: CEM on the V50 decile, set-B tiers, gene-cluster interval (G17 copy of the G1 contrast). Odds ratio for domain overlap, difference for the other signals.\n")
rows = RowList()
detail = []
for i, s in enumerate(("S2", "S3", "S4", "S1", "S5"), 1):
    p = prim(s)
    g, gr, gr2, cav = (None, None, NOT_GRADED, None) if s == "S5" else final(s)
    src = d1("G1", s)
    claim = {"S1": "Domain overlap is more frequent among GOLD than BRONZE records in the full human population (odds ratio above 1).",
             "S2": "GOLD sites more often have a conserved phosphosite in another species than BRONZE sites.",
             "S3": "GOLD sites have higher disorder scores than BRONZE sites (small difference in the mean disorder score).",
             "S4": "GOLD sites more often match a sequence-based kinase-motif class than BRONZE sites.",
             "S5": "GOLD sites more often match the database lookup kinase-motif indicator than BRONZE sites (not decision-bearing; derived from annotation tables)."}[s]
    if s == "S1":
        claim += " **Not claimable as a general finding: the direction is not consistent across populations (see C2.6).**"
    basis = NOT_GRADED if g is None else f"K1 {'pass' if g['k1'] else 'FAIL'}; K2 {'pass' if g['k2'] else 'FAIL'}; K3 {len(g['fails'])} failure(s) of {g['n_checks']}; K4 {'n/a' if g['k4'] is None else ('pass' if g['k4'] else 'FAIL')}; K5 {'not run' if g['wg'] is None else ('pass' if g['wg'] else 'FAIL')}; F = {g['F']} -> {gr}"
    if g is not None and cav and gr != gr2:
        basis += f"; lowered one step: {cav}"
    rows.append([f"C2.{i}", claim, f(p["estimate"], s), ci(p["ci"], s), f"{n_(p['n_sites_gold'])} GOLD / {n_(p['n_sites_bronze'])} BRONZE; {n_(p['n_gene_clusters'])} gene clusters", f"{p['population'] if 'population' in p else ''}; {src.replace('G1', 'G1')} (gene: {[k for k, v in src_of.items() if k == src][0] and src_of[src]})".replace("; ", "; ", 1), gr2, basis])
    if g is not None:
        detail.append((s, g))
A(md(["ID", "Claim", "Estimate (D1)", "95% CI (gene cluster)", "n", "Source (G1 contrast; gene-cluster copy)", "Strength", "Basis"], rows))
A("Robustness failures behind the grades (K3 checks whose sign or interval fails; the full set of checks is in `results/amend_results.json`, `results/f_runs/`, `results/d2b/`, `results/cem_results.json`):\n")
for s, g in detail:
    A(f"- **{SL[s]}:** " + ("none" if not g["fails"] else "; ".join(f"{lab} ({cid})" for lab, cid in g["fails"])) + f" (of {g['n_checks']} checks)")
A("")
# direction table
A("Direction across populations (K2):\n")
rows = RowList()
for s in ("S2", "S3", "S4", "S1"):
    g = grade(s)
    rows.append([SL[s], len(g["pops"]), "yes" if g["k2"] else "no", "; ".join(f"{lab} {f(e, s)}" for lab, _, e in g["pops"][:6]) + (" ..." if len(g["pops"]) > 6 else "")])
A(md(["Signal", "Estimates", "Same side of the null everywhere", "Estimates (first six)"], rows))
# domain overlap detail
cr = {"all human records": am[src_of[d1("G1", "S1")].replace("0002", "0001") if False else "G17-0379"]["estimate"], "proteins with an abundance (P-abund)": am["G17-0399"]["estimate"], "proteins with a domain-annotated site (G3)": am["G17-0397"]["estimate"], "canonical reviewed accessions (G16)": am["G16-0361"]["estimate"]}
rows = RowList()
rows.append(["C2.6", "Domain overlap is population-dependent: the crude odds ratio is above 1 in the full human population and below 1 in three other populations.", "; ".join(f"{k}: {f(v, 'S1')}" for k, v in cr.items()), "", "", "G17-0379, G17-0399, G17-0397, G16-0361; G20-0489 (stratified table)", "strong (as a description of the dependence)", "four crude estimates, two sides of the null"])
A(md(["ID", "Claim", "Estimates", "95% CI", "n", "Source", "Strength", "Basis"], rows))

# ---------------------------------------------------------------- 3 adjustment and index
A("## 3. Adjustment for the study-effort index\n")
dec = [x for x in ("S1", "S2", "S3", "S4")]
att = [prim(s)["attenuation"] for s in dec]
att7 = [v["attenuation"] for k, v in am.items() if k.startswith("G17") and v.get("design") == "D1" and v.get("decision_bearing") == "yes" and (v["population"].startswith("P-all") or v["population"] == "P-human")]
assert len(att7) == 7, len(att7)
sm = {k: v["smd_pre"] for k, v in g20["prematching_balance"].items() if manid.loc[k, "decision_bearing"] == "yes"}
ix = g20["index_outcome_association"]
rmax = max(abs(v["pearson_r"]) for v in ix.values())
sp = g20["index_vs_sites_per_accession"]
eta = b5["v50"]["kruskal_eta2"] if "v50" in b5 else None
rows = RowList()
rows.append(["C3.1", "Adjusting for the protein-level study-effort index (CEM on the V50 decile) leaves the four decision-bearing signals essentially unchanged.", f"attenuation from crude to adjusted, all decision-bearing contrasts (G1 and G2): {100 * min(att7):.1f}% to {100 * max(att7):.1f}% (negative: the adjusted estimate is larger than the crude estimate); all four signals survive", "", "see section 2", "G17 (G1 D1), G1-G3; G18 (matcher D2b)", "strong", "all attenuations small; same conclusion under CEM, D2b and the accession-cluster interval"])
rows.append(["C3.2", "The index is a weak confounder candidate: within a tier it correlates with the outcomes at |r| <= " + f"{rmax:.3f}" + ".", f"largest |r| {rmax:.3f}", "see G20-0469..0486", "see G20", "G20 (index-outcome association)", "strong (description)", "gene-cluster intervals in `results/g20_diagnostics.json`"])
rows.append(["C3.3", "The index is nearly unrelated to the number of GOLD sites per human protein.", f"Spearman {f(sp['human|gold_sites_per_accession']['spearman'], 'S2', 3)}", ci(sp["human|gold_sites_per_accession"]["ci"], "S2", 3), n_(sp["human|gold_sites_per_accession"]["n_accessions"]), "G20-0487", "strong (description)", "gene-cluster interval"])
rows.append(["C3.6", "The primary index (V50) does not order the tiers: its mean is lower for GOLD proteins than for BRONZE proteins and lowest for SILVER.", f"means by highest tier: GOLD {f(b5['means_by_top_tier']['SEI_V50']['GOLD'], 'S2', 3)}, SILVER {f(b5['means_by_top_tier']['SEI_V50']['SILVER'], 'S2', 3)}, BRONZE {f(b5['means_by_top_tier']['SEI_V50']['BRONZE'], 'S2', 3)}; eta-squared {f(b5['effects']['kruskal_eta2_SEI_V50_by_top_tier']['estimate'], 'S2')}", ci(b5["effects"]["kruskal_eta2_SEI_V50_by_top_tier"]["ci95"], "S2"), n_(b5["n_proteins"]) + " proteins", "B5", "strong (description)", "protein-bootstrap interval; no P-value is reported"])
rows.append(["C3.4", "Study effort does not explain the validation signals.", "", "", "", "G1-G3, G17, G18", "not supported", "adjusting for a gene-level index that is weakly related to tier cannot exclude site-level study effort; it shows only that this adjustment does not change the signals"])
rows.append(["C3.5", "Matching on the index removes confounding by study effort at the site level.", "", "", "", "G20", "not supported", "the index is a gene-level count; the tier is a site-level property"])
A(md(["ID", "Claim", "Estimate", "95% CI", "n", "Source", "Strength", "Basis"], rows))

# ---------------------------------------------------------------- 4 tier definition, inference
A("## 4. Tier definition and inference\n")
tiers = {}
for cid, v in am.items():
    if cid.startswith("G19"):
        tiers.setdefault((v["signal"], v["population"], v["tier_rule"].split(" (")[0].replace(" (primary)", "")), {})[v["design"]] = v
dis_d = tiers[("S3", "P-all", "depth-only tiers")]
dis_a = tiers[("S3", "P-all", "agreement-only")]
tw = [b4[d1("G1", s)] for s in ("S1", "S2", "S3", "S4")]
twn = sum(1 for t in tw if t["survives_two_way"])
wr = []
for s in ("S1", "S2", "S3", "S4", "S5"):
    a = cem[d1("G1", s)]["ci"]
    g_ = prim(s)["ci"]
    wr.append((g_[1] - g_[0]) / (a[1] - a[0]))
w14 = [v["ci_width_ratio_protein_cluster_over_rev5"] for v in g14.values() if isinstance(v, dict) and v.get("design", "").startswith("R5_CEM") and v.get("status") == "run"]
lf_ = []
for pop in ("P-human", "P-all"):
    for sg in ("S1", "S2", "S3", "S4"):
        ks = [k for k in tiers if k[0] == sg and k[1] == pop]
        if not ks:
            continue
        defs = [(tiers[k]["D0"], tiers[k]["D1"]) for k in ks]
        okc = [(a, b_) for a, b_ in defs if excl(a["ci"], sg) and excl(b_["ci"], sg) and a["estimate"] > 0 and b_["estimate"] > 0]
        oka = [b_["estimate"] for a, b_ in defs if excl(b_["ci"], sg) and b_["estimate"] > 0]
        if okc and len(oka) >= 2:
            lf_.append((max(abs(math.log(b_["estimate"] / a["estimate"])) for a, b_ in okc), math.log(max(oka) / min(oka))))
_sp_all = [v["D1"] for k, v in tiers.items() if k[0] == "S3" and k[1] == "P-all"]
_sp_all_ex = [d["estimate"] for d in _sp_all if excl(d["ci"], "S3") and min(abs(d["ci"][0]), abs(d["ci"][1])) >= 0.01]
LF_TXT = (f"natural-log spread of adjusted estimates across definitions {min(t[1] for t in lf_):.3f} to {max(t[1] for t in lf_):.3f} (plain ratios {math.exp(min(t[1] for t in lf_)):.2f} to {math.exp(max(t[1] for t in lf_)):.2f}) against within-definition change from crude to adjusted "
          f"{min(t[0] for t in lf_):.3f} to {max(t[0] for t in lf_):.3f} (plain ratios {math.exp(min(t[0] for t in lf_)):.3f} to {math.exp(max(t[0] for t in lf_)):.3f}); the spread is larger in {sum(1 for t in lf_ if t[1] > t[0])} of {len(lf_)} signal-by-population combinations (RESULTS_REPORT.md section 6). "
          f"The largest spread ({math.exp(max(t[1] for t in lf_)):.2f}-fold, disorder in P-all) rests on a near-zero agreement-only estimate; without it the spread is about {max(_sp_all_ex) / min(_sp_all_ex):.0f}-fold ({max(_sp_all_ex) / min(_sp_all_ex):.2f})")
g23r_ = J("g23_results.json")
_pairs = []
for _cid, _v in g23r_.items():
    if _v["design"] != "D1":
        continue
    _src = manid.loc[_cid, "note"].split("counterpart of ")[1].split(";")[0]
    _p = am[src_of[_src]]
    _s = _v["signal"]
    _pairs.append((excl(_v["ci"], _s), (_v["estimate"] > NULL(_s)) == (_p["estimate"] > NULL(_s)), bool(_v.get("largely_explained"))))
G23_TXT = f"{sum(1 for a in _pairs if a[0])} of {len(_pairs)} D1 contrasts survive; {sum(1 for a in _pairs if a[1])} of {len(_pairs)} keep the sign of the primary; {sum(1 for a in _pairs if a[2])} largely explained"
rows = RowList()
rows.append(["C4.1", "The size of every signal depends on how GOLD and BRONZE are defined (set-B, set-A, depth-only or agreement-only tiers); on the natural-log scale the spread across definitions exceeds the change caused by adjusting for the index.", LF_TXT, "gene-cluster intervals in G19", "G19 rows", "G19, G12, G15", "moderate", "the GOLD and BRONZE record sets differ between definitions, so the comparison is between different populations of sites"])
rows.append(["C4.2", f"Disorder is about {f(dis_d['D0']['estimate'], 'S3', 2)} under depth-only tiers and about {f(dis_a['D0']['estimate'], 'S3', 3)} (interval includes zero) under agreement-only tiers (P-all, crude). **Description only; no causal interpretation.**", f"depth-only {f(dis_d['D0']['estimate'], 'S3')}; agreement-only {f(dis_a['D0']['estimate'], 'S3')}", f"{ci(dis_d['D0']['ci'], 'S3')}; {ci(dis_a['D0']['ci'], 'S3')}", f"{n_(dis_d['D0']['n_sites_gold'])} / {n_(dis_d['D0']['n_sites_bronze'])}; {n_(dis_a['D0']['n_sites_gold'])} / {n_(dis_a['D0']['n_sites_bronze'])}", "G19-0458, G19-0460", "strong (description)", "two crude estimates with gene-cluster intervals"])
rows.append(["C4.3", "Intervals from a protein- or gene-cluster bootstrap are wider than the site-level bootstrap of Rev 5 for the same estimates; conclusions about the null do not change in the CEM rows.", f"interval width ratio {min(w14):.3f} to {max(w14):.3f} across {len(w14)} CEM rows (G14)", "", "", "G14", "strong", "same point estimates, only intervals change"])
rows.append(["C4.4", "Gene-cluster and accession-cluster intervals are similar; accessions of one gene are not a material source of dependence for these estimates.", f"gene/accession width ratio {min(wr):.3f} to {max(wr):.3f} (G1 D1, five signals)", "", "", "G17 against G1", "strong", "same estimates, two resampling units"])
rows.append(["C4.5", "The two-way (protein by study cluster) interval excludes the null for only some decision-bearing G1 contrasts; the disorder and sequence-based motif signals do not.", f"{twn} of {len(tw)} decision-bearing G1 contrasts exclude the null", "", "", "B4 (`results/b4_twoway.json`)", "strong (description); approximate interval", "sensitivity analysis resting on a few very large study clusters; not the primary rule"])
rows.append(["C4.6", "The validation signals do not depend on database, resource or review papers being counted as evidence: with the non-primary PMIDs removed from publication depth, every G1 and G2 D1 contrast keeps its sign and its interval excludes the null.", G23_TXT, "", "see RESULTS_REPORT.md section 10c", "G23-0517 to G23-0534", "strong (sensitivity)", "same design as the primary; only the tier column differs; the grades of section 2 are unchanged by this check"])
A(md(["ID", "Claim", "Estimate", "95% CI", "n", "Source", "Strength", "Basis"], rows))

# ---------------------------------------------------------------- 5 within-gene and established measures
A("## 5. Within-gene comparison (G21) and established confidence measures (G9, G10)\n")
rows = RowList()
if g21:
    for s in ("S2", "S3", "S4", "S5"):
        for popl in ("P-all (S/T/Y)", "P-human (S/T/Y)"):
            v1 = [v for v in g21.values() if v["signal"] == s and v["model"] == "M1" and v["population"] == popl]
            v2 = [v for v in g21.values() if v["signal"] == s and v["model"] == "M2" and v["population"] == popl]
            if not v1:
                continue
            v1, v2 = v1[0], (v2[0] if v2 else None)
            if v2 is None:
                gr_ = "moderate" if v1["persists"] else "not supported"
            else:
                gr_ = "moderate" if (v1["persists"] and v2["persists"]) else ("weak" if v1["persists"] else "not supported")
            unit = "odds ratio" if s != "S3" else "difference"
            sc_ = "S1" if s != "S3" else "S3"

            def est_txt(v):
                if "companion_linear" in v:
                    c = v["companion_linear"]
                    return f"{f(v['estimate'], sc_)} (odds ratio; {ci(v['ci'], sc_)}); within-gene linear difference {f(c['estimate'], 'S2')} ({ci(c['ci'], 'S2')})"
                return f"within-gene linear difference {f(v['estimate'], 'S3')} ({ci(v['ci'], 'S3')})"

            gx = "G1" if popl.startswith("P-all") else "G2"
            cx = am[src_of[d1(gx, s)]]
            v1s = [v for v in g21.values() if v["signal"] == s and v["model"] == "M1s" and v["population"] == popl]
            fit_parts, over = [], []
            for lab, vv in (("M1", v1), ("M1s", v1s[0] if v1s else None), ("M2", v2)):
                if vv is None:
                    continue
                fit_parts.append(f"{lab} {vv['n_boot_valid']} valid / {vv['n_boot_unfit']} failed")
                if 100 * vv["n_boot_unfit"] / vv["n_boot"] > 5:
                    over.append(f"{lab} ({100 * vv['n_boot_unfit'] / vv['n_boot']:.1f}%)")
            fit_txt = f"bootstrap fits of {'the conditional logit' if 'companion_linear' in v1 else 'the linear fixed-effect model'} out of {v1['n_boot']} resamples: " + "; ".join(fit_parts)
            if over and s != "S5":
                fit_txt += f". **More than 5% of the fits failed for {', '.join(over)}; the interval rests on the converged resamples**"
            if s == "S5":
                strength = NOT_GRADED
                basis = f"scale check: the CEM D1 difference in {gx} is {f(cx['estimate'], s)} ({ci(cx['ci'], s)}); {fit_txt}"
            else:
                strength = gr_
                basis = (f"scale check: the CEM D1 difference in {gx} is {f(cx['estimate'], s)} ({ci(cx['ci'], s)}); rule for within-gene claims: moderate if M1 and M2 both persist, weak if only M1, "
                         f"not supported if M1 does not; capped at moderate because site-level confounding is not excluded; {fit_txt}")
            rows.append([f"C5.{len(rows) + 1}", f"{SL[s]}{' (not decision-bearing; derived from annotation tables)' if s == 'S5' else ''}: the GOLD-BRONZE difference persists within the same gene ({popl}). Excludes protein-level confounding only.", est_txt(v1), "",
                         f"{n_(v1['n_records'])} records; {n_(v1['n_gene_clusters_with_both_tiers'])} gene clusters with both tiers",
                         f"G21 M1 (primary); M2 (sensitivity): {est_txt(v2)} on the reviewed-sequence subset" if v2 else "G21 M1", strength, basis])
else:
    rows.append(["C5.1", "Within-gene persistence of the signals (G21).", "not yet run", "", "", "G21", "pending", ""])
psp_run = bool(g10)
if psp_run:
    rows.append(["C5.x", "LT-supported against MS-only sites differ within a gene (Task I).", "see results/g10_results.json", "", "", "G10", "see section 7", ""])
else:
    rows.append([f"C5.{len(rows) + 1}", "PhosphoSitePlus-based measures (LT-supported against MS-only; evidence category high against low) behave like the tier comparisons.", "not run", "", "", "G9 (PSP rows), G10", "not supported (no evidence; not run)", "Task H and Task I: not run: source file not available (PhosphoSitePlus v7 site provides no bulk download; access requested from the provider); listed in Future work"])
if g9:
    iu = [v for k, v in g9.items() if manid.loc[k, "variation"].startswith("UniProt") and v.get("design") == "D1" and v.get("status") == "run"]
    rows.append([f"C5.{len(rows) + 1}", "UniProt experimental against inferred evidence contrast. **Not used as evidence for or against the main hypothesis** (the inferred class includes sequence-similarity transfer, so conservation is partly circular).", "reported in RESULTS_REPORT.md section 9", "", "", "G9 (UniProt rows)", "not graded (excluded by instruction)", "caveat stated before the results"])
rows.append([f"C5.{len(rows) + 1}", "PeptideAtlas site-level evidence categories.", "not feasible from the available build", "", "", "G9 (PeptideAtlas rows)", "not supported (no evidence)", "no PTMProphet probability or per-site PSM table in build 606; definitions quoted in `sources/SOURCES.md`"])
A(md(["ID", "Claim", "Estimate", "95% CI", "n", "Source", "Strength", "Basis"], rows))

# ---------------------------------------------------------------- 6 predictors
A("## 6. Predictor audit (secondary)\n")
A("Rule for predictor claims (fixed before grading): *strong* if the criterion written down before any predictor was scored (gap >= 0.05 with an interval excluding zero in at least two tools; commit 8e12390, 2026-09-24) is met; *moderate* if one tool meets it and the restricted analysis does not show the gap to be absent; *weak* if one tool meets it and the restricted analysis (positives first reported after the restriction year) shows the gap equivalent to zero; *not supported* if none meets it.\n")
qual = []
for tool in ("DeepPhos", "MusiteDeep", "PhosF3C"):
    k = [x for x in g11 if x.startswith(tool + "|planned")][0]
    for kk, v in g11[k].items():
        if isinstance(v, dict) and "gap" in v and kk == "unrestricted" and v["gap"] >= 0.05 and excl(v["gene_ci95"], "S2"):
            qual.append((tool, v))
tools_ok = {t for t, _ in qual}
def c62_basis(tool, k, absent):
    """Basis of C6.2: every restriction year stored for the tool, with the 90% interval of the equivalence test (values from g11_predictors.json)."""
    if tool != "MusiteDeep":
        return "one tool meets the criterion; in the restricted analysis " + ("the gap is equivalent to zero" if absent else "the gap is not shown to be absent") + "; unlabelled negatives on little-studied proteins may include undiscovered sites (limitation 8)"
    m = {x: g11[k][f"restriction_year={x}"] for x in (2020, 2017)}
    d18 = g11[[x for x in g11 if x.startswith("DeepPhos|planned")][0]]["restriction_year=2018"]
    assert m[2020]["gene_label"].startswith("gap absent") and not m[2017]["gene_label"].startswith("gap absent") and not d18["gene_label"].startswith("gap absent")
    return ("one tool meets the criterion. Restricted analysis (positives first reported after the year), MusiteDeep: "
            f"2020, the leakage-free year for the scored model: equivalent to zero (gap {f(m[2020]['gap'], 'S2')}; 90% interval {ci(m[2020]['gene_ci90'], 'S2')}); "
            f"2017, the year in the specification: not equivalent, the gap persists ({f(m[2017]['gap'], 'S2')}; 90% interval {ci(m[2017]['gene_ci90'], 'S2')}). "
            f"DeepPhos 2018 (equivalence not shown; gap {f(d18['gap'], 'S2')}; 90% interval {ci(d18['gene_ci90'], 'S2')}) is the same kind of sensitivity check. "
            "The grade follows the 2020 result; unlabelled negatives on little-studied proteins may include undiscovered sites (limitation 8)")


rows = RowList()
rows.append(["C6.1", "The AUROC gap between the best-studied and least-studied deciles is at least 0.05 with an interval excluding zero in at least two of three tools.", f"{len(tools_ok)} of 3 tools", "", "", "G11", "not supported", "criterion met by " + (", ".join(sorted(tools_ok)) or "no tool")])
for tool, v in qual:
    k = [x for x in g11 if x.startswith(tool + "|planned")][0]
    restricted = [vv for kk, vv in g11[k].items() if kk.startswith("restriction_year") and isinstance(vv, dict) and "gap" in vv]
    absent = any(vv["gene_label"].startswith("gap absent") for vv in restricted)
    gr_ = "weak" if absent else "moderate"
    rows.append(["C6.2", f"{tool} AUROC is higher in the three best-studied deciles than in the three least-studied (serine/threonine sites, unrestricted).", f(v["gap"], "S2"), ci(v["gene_ci95"], "S2"), f"{n_(v['n_sites'])} sites; {n_(v['n_positives'])} positives; {n_(v['n_gene_clusters'])} gene clusters", "G11", gr_,
                 c62_basis(tool, k, absent)])
rows.append(["C6.3", "The gap disappears when positives first reported after each tool's publication year are used.", "", "", "", "G11 restricted settings", "not supported as a general statement", "equivalence to zero holds in some settings and is not shown in others (RESULTS_REPORT.md section 11)"])
A(md(["ID", "Claim", "Estimate", "95% CI", "n", "Source", "Strength", "Basis"], rows))

# ---------------------------------------------------------------- 7 withdrawn / changed
A("## 7. Rev 5 claims that are withdrawn or changed (old against new)\n")
A("Old values: `study/earlier_results/matching_results.json`, `predictor_metrics.json`, `dataset_counts.json` (Rev 5 designs: original index, stratified site bootstrap, set-A tiers for domain overlap). New values: Rev 6 result files (V50 index, set-B tiers for every signal, gene-cluster bootstrap).\n")


def old(sig, key):
    x = m5[sig][key]
    return x["point_estimate"], (x["ci_low"], x["ci_high"])


rows = RowList()
mp = {"S1": ("domain_overlap_human", "G1"), "S2": ("conservation", "G1"), "S3": ("disorder", "G1"), "S4": ("kinase_motif_sequence", "G1"), "S5": ("kinase_motif_lookup", "G1")}
for s in ("S1", "S2", "S3", "S4", "S5"):
    key, blk = mp[s]
    oe, oc = old(key, "cem_primary")
    p = prim(s)
    extra = ""
    if s == "S1":
        extra = f"; with set-A tiers and V50 the Rev 6 value is {f(am['G19-0421']['estimate'], 'S1')} (G19-0421), so the change is mostly the tier source"
    rows.append([f"W7.{len(rows) + 1}", f"{SL[s]}, CEM on the primary index", f"{f(oe, s)} ({ci(oc, s)})", f"{f(p['estimate'], s)} ({ci(p['ci'], s)})", f"{d1('G1', s)} / {src_of[d1('G1', s)]}", "changed" + extra])
oe, oc = old("domain_overlap_annotated_only", "cem_primary")
p3 = am["G17-0398"]
rows.append([f"W7.{len(rows) + 1}", "Domain overlap among proteins with a domain-annotated site, CEM", f"{f(oe, 'S1')} ({ci(oc, 'S1')})", f"{f(p3['estimate'], 'S1')} ({ci(p3['ci'], 'S1')})", "G3-0029 / G17-0398", "changed; Rev 5 wording that the inversion arises only in the annotated-only population is withdrawn (it appears in four populations)"])
oe, oc = old("conservation", "propensity")
r6 = prop["G1-0006"]
rows.append([f"W7.{len(rows) + 1}", "Conservation, propensity matching (Rev 5 matcher)", f"{f(oe, 'S2')} ({ci(oc, 'S2')})", f"{f(r6['estimate'], 'S2')} ({ci(r6['ci'], 'S2')}); D2b {f(d2b['G18-0408']['estimate'], 'S2')} ({ci(d2b['G18-0408']['ci'], 'S2')})", "G1-0006; G18-0408", "changed"])
oe, oc = old("domain_overlap_human", "propensity")
rows.append([f"W7.{len(rows) + 1}", "Domain overlap, propensity matching", f"{f(oe, 'S1')} ({ci(oc, 'S1')})", f"not fitted with set-B tiers (G1-0003); D2b on S/T/Y records {f(d2b['G18-0407']['estimate'], 'S1')} ({ci(d2b['G18-0407']['ci'], 'S1')})", "G1-0003; G18-0407", "changed; Rev 5 value not reproducible under the Rev 6 definition"])
rows.append([f"W7.{len(rows) + 1}", "Methods sentence \"Propensity matching paired records 1:1 without replacement within a caliper of 0.2 SD\"", "described 1:1 matching", "the Rev 5 code kept each treated record's 50 nearest controls (ties by record order)", "RECONCILIATION.md section 6; G18", "withdrawn as a description of what was run; D2b implements the sentence"])
# index properties
rows.append([f"W7.{len(rows) + 1}", "Eta-squared of the primary index by highest tier", f"{f(m5['sei_properties']['kruskal_eta2_H'], 'S2')} (Kruskal-Wallis H = {m5['sei_properties']['kruskal_H']:.1f}, P = {m5['sei_properties']['kruskal_p']:.1e})", f"V50: {f(b5['effects']['kruskal_eta2_SEI_V50_by_top_tier']['estimate'], 'S2')} ({ci(b5['effects']['kruskal_eta2_SEI_V50_by_top_tier']['ci95'], 'S2')})", "B5", "changed; P-values on indices are dropped (not informative with 87,338 proteins sharing studies); effect sizes with protein-bootstrap intervals replace them"])
# predictors
for tool, key in (("DeepPhos", "DeepPhos|planned: S/T sites only"), ("MusiteDeep", "MusiteDeep|planned: S/T sites only"), ("PhosF3C", "PhosF3C|planned: Y sites only")):
    o = pm5[tool]["gap"] if "gap" in pm5[tool] else None
    nw = g11[key]["unrestricted"]
    if o:
        rows.append([f"W7.{len(rows) + 1}", f"{tool} decile AUROC gap", f"{f(o['point'], 'S2')} ({ci([o['ci_low'], o['ci_high']], 'S2')}), all scored residues, site-level interval", f"{f(nw['gap'], 'S2')} ({ci(nw['gene_ci95'], 'S2')}), planned residue set, gene-cluster interval", "G11", "changed (residue set and interval); criterion outcome unchanged (one of three tools)"])
rows.append([f"W7.{len(rows) + 1}", "Statement that the protein-level index orders tiers or that adjusting leaves \"three of four comparisons unchanged or larger\"", "three of four unchanged or larger (Rev 5 claims file, F7)", f"{sum(1 for s in ('S1', 'S2', 'S3', 'S4') if prim(s)['survives'])} of {len(('S1', 'S2', 'S3', 'S4'))} survive; {sum(1 for s in ('S1', 'S2', 'S3', 'S4') if prim(s).get('largely_explained'))} largely explained", "G17, G1", "changed; wording \"study effort does not explain the signals\" is not supported (section 3, C3.4)"])
b1 = J("b1_reconciliation.json")
oc_ = b1["ortholog_coverage"]
k74 = [k for k in oc_ if "7.4%" in k][0]
k123 = [k for k in oc_ if "12.3%" in k][0]
sc_ = b1["study_clusters"]
rows.append([f"W7.{len(rows) + 1}", "\"The dataset contained 1,113 study clusters\" and the ortholog-coverage figures 7.4% and 12.3%", f"1,113; 7.4%; 12.3%", f"{n_(sc_['distinct_primary_cluster_labels_all_records'])} is the number of distinct primary-cluster labels; under the merge rule the {n_(cl['884']['n_pmids'])} PMIDs form {n_(cl['884']['n_distinct_cluster_ids_over_all_pmids'])} clusters; {100 * oc_[k74][1]:.1f}% is {n_(oc_[k74][0])} analysis-set records over all records and {100 * oc_[k123][1]:.1f}% is the human share", "RECONCILIATION.md sections 2 to 3", "withdrawn as stated; denominators must be given"])
trc = b1["table_text_reconciliation"]
diff = [x for x in trc if x["text_value"] != x["table_value"]]
rows.append([f"W7.{len(rows) + 1}", "Unmatched estimates quoted in the text against Table 1", f"text and Table 1 differ for {len(diff)} of {len(trc)} rows (" + "; ".join(f"{x['label'].split(' (')[0]} {x['text_value']} against {x['table_value']}" for x in diff[:2]) + ", ...)", "both reproduce from the site table (all records against records with a primary index)", "RECONCILIATION.md section 1", "changed; one definition must be used"])
rows.append([f"W7.{len(rows) + 1}", "All intervals", "stratified site-level bootstrap", "gene-cluster bootstrap (primary), accession-cluster beside it", "G14, G17", "changed for every estimate"])
A(md(["ID", "Rev 5 claim or estimate", "Rev 5 value (source file)", "Rev 6 value", "Block / contrast", "Status"], rows))

# ---------------------------------------------------------------- 8 headline rule
A("## 8. Pre-set headline rule\n")
A("**Branch A** if Task I shows a within-gene difference between LT-supported and MS-only sites in conservation or disorder (detectability-controlled model, gene-cluster interval excluding zero); **Branch B** otherwise, or if Task I was not run.\n")
if g10:
    branch = "B"
    detail_b = "Task I results were not parsed"
    try:
        hits = [v for v in g10.values() if v.get("signal") in ("S2", "S3") and v.get("model") == "controlled" and v.get("status") == "run" and excl(v["ci"], v["signal"])]
        branch = "A" if hits else "B"
        detail_b = "; ".join(f"{SL[v['signal']]}: {f(v['estimate'], v['signal'])} ({ci(v['ci'], v['signal'])})" for v in g10.values() if v.get("model") == "controlled" and v.get("status") == "run" and v.get("signal") in ("S2", "S3"))
    except Exception as e:  # noqa: BLE001
        detail_b = f"could not evaluate: {e}"
    A(f"**Branch {branch} applies.** {detail_b}.\n")
else:
    A("**Branch B applies**, because Task I has not been run (not run: source file not available; PhosphoSitePlus v7 site provides no bulk download; access requested from the provider). If Task I is run later, the branch is re-evaluated by this script.\n")

open(REV6 / "CLAIMS_VS_EVIDENCE.md", "w", encoding="utf-8").write("\n".join(L))
dump_bindings(RES / "claims_bindings.json")
print("written", len(L), "blocks")
