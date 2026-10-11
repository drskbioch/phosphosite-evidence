"""Task B1: reconcile Table 1 and the text of the Rev 5 manuscript with the result files and the site table.

Every estimate is recomputed from `data_private/base.parquet` (the record-level table built from the release files)
under both record sets (all GOLD and BRONZE records with the outcome defined; the subset whose protein has a primary
index) and compared with the value quoted in manuscript.md and in Table 1. Reads manuscript files; edits nothing.
Writes results/b1_reconciliation.json and sections 1 to 5 of study/analysis/RECONCILIATION.md. Sections 6 and 7 of that file are written by hand; if the file exists, everything from its "## 6." heading on is kept unchanged.
"""
import json
import re

import numpy as np
import pandas as pd

import base as B
from common import CUR, REL, REV6, RES

b = B.load()
M = json.load(open(CUR / "earlier_results" / "matching_results.json"))
D = json.load(open(CUR / "earlier_results" / "dataset_counts.json"))
MS = (CUR / "manuscript" / "manuscript.md").read_text(encoding="utf-8")
T1 = (CUR / "manuscript" / "tables" / "_rev5" / "table1_matched_comparisons.md").read_text(encoding="utf-8")
mapped = (b["mapped_to_gene"] == True).to_numpy()  # noqa: E712


def crude(d, y, stat):
    g, c = d["tier"] == "GOLD", d["tier"] == "BRONZE"
    yg, yc = d.loc[g, y].astype(float), d.loc[c, y].astype(float)
    if stat == "or":
        a, bb, cc, dd = yg.sum(), (1 - yg).sum(), yc.sum(), (1 - yc).sum()
        return float(a * dd / (bb * cc)), int(g.sum()), int(c.sum())
    return float(yg.mean() - yc.mean()), int(g.sum()), int(c.sum())


def pop(tier_col, ycol, human=False, protein_filter=None, only_mapped=False):
    m = b[tier_col].isin(["GOLD", "BRONZE"]).to_numpy() & b[ycol].notna().to_numpy()
    if human:
        m &= (b["species"] == "human").to_numpy()
    if protein_filter is not None:
        m &= b["uniprot_id"].isin(protein_filter).to_numpy()
    if only_mapped:
        m &= mapped
    d = b.loc[m, ["tier_setA" if tier_col == "tier_setA" else "tier", ycol]].copy()
    d.columns = ["tier", ycol]
    return d, int(m.sum())


has = b[b["species"] == "human"].groupby("uniprot_id")["in_domain"].max()
annot = set(has[has == 1.0].index)
spec = [  # key, label, tier column, outcome column, human only, statistic, protein filter
    ("domain_overlap_human", "Domain overlap, human (OR)", "tier_setA", "in_domain", True, "or", None),
    ("domain_overlap_annotated_only", "Domain overlap, annotated proteins only, human (OR)", "tier_setA", "in_domain", True, "or", annot),
    ("conservation", "Conservation (difference)", "tier", "conserved", False, "diff", None),
    ("disorder", "Disorder (difference)", "tier", "disorder", False, "diff", None),
    ("kinase_motif_sequence", "Kinase motif, sequence-based (difference)", "tier", "seq_motif_match", False, "diff", None),
    ("kinase_motif_lookup", "Kinase motif, database lookup (difference)", "tier", "kinase_match", False, "diff", None),
]
text_quote = {  # value quoted in the text for the unmatched estimate (regex on manuscript.md)
    "domain_overlap_human": ("2.494", r"odds of domain overlap for GOLD relative to BRONZE records were 2\.494"),
    "domain_overlap_annotated_only": ("0.586", r"the corresponding odds ratios were 0\.586"),
    "conservation": ("0.3760", r"conserved site was 0\.3760"),
    "disorder": ("0.0310", r"mean disorder score was 0\.0310"),
    "kinase_motif_sequence": ("0.0433", r"motif match was 0\.0433"),
    "kinase_motif_lookup": ("0.7780", r"lookup indicator was 0\.7780"),
}
table_row = {"domain_overlap_human": "Domain overlap, human (odds ratio)",
             "domain_overlap_annotated_only": "Domain overlap, annotated proteins only",
             "conservation": "Conservation", "disorder": "Disorder",
             "kinase_motif_sequence": "Kinase motif, sequence-based", "kinase_motif_lookup": "Kinase motif, database lookup"}
rows = []
for key, label, tcol, ycol, human, stat, pf in spec:
    da, na = pop(tcol, ycol, human, pf)
    dm, nm = pop(tcol, ycol, human, pf, only_mapped=True)
    va, ng_a, nb_a = crude(da, ycol, stat)
    vm, ng_m, nb_m = crude(dm, ycol, stat)
    ja, jm = M[key]["unmatched_all_sites"]["point_estimate"], M[key]["unmatched_mapped_sites"]["point_estimate"]
    tq, rx = text_quote[key]
    assert re.search(rx, MS), f"text string for {key} not found in manuscript.md"
    trow = [ln for ln in T1.splitlines() if ln.startswith("| " + table_row[key])][0]
    tval = [c.strip() for c in trow.split("|")[1:-1]][1].split(" ")[0]
    dec = 3 if stat == "or" else 4
    rows.append({
        "key": key, "label": label, "tier_source": "source set A (human file)" if tcol == "tier_setA" else "source set B (all three species)",
        "species": "human" if human else "human, mouse, rat",
        "text_value": tq, "table_value": tval,
        "all_records": {"n": na, "n_gold": ng_a, "n_bronze": nb_a, "recomputed": va, "result_file": ja, "abs_diff": abs(va - ja)},
        "mapped_records": {"n": nm, "n_gold": ng_m, "n_bronze": nb_m, "recomputed": vm, "result_file": jm, "abs_diff": abs(vm - jm)},
        "n_unmapped_excluded_result_file": M[key]["n_unmapped_excluded"], "n_unmapped_recomputed": na - nm,
        "text_matches": "all records" if f"{va:.{dec}f}" == tq else ("records with an index" if f"{vm:.{dec}f}" == tq else "neither"),
        "table_matches": "all records" if f"{va:.{dec}f}" == tval else ("records with an index" if f"{vm:.{dec}f}" == tval else "neither"),
        "n_analysis_matched": M[key]["cem_primary"]["n_analysis"],
    })
    # Rev 5 did not drop undefined outcomes for the kinase-motif rows before counting unmapped records
    if key.startswith("kinase_motif"):
        mm = b[tcol].isin(["GOLD", "BRONZE"]).to_numpy()
        rows[-1]["n_unmapped_recomputed"] = int((mm & ~mapped).sum())
    assert rows[-1]["n_unmapped_excluded_result_file"] == rows[-1]["n_unmapped_recomputed"], key
    assert rows[-1]["all_records"]["abs_diff"] < 1e-9 and rows[-1]["mapped_records"]["abs_diff"] < 1e-9, key

# ------------------------------------------------- ortholog coverage: 7.4 % and 12.3 %
N = len(b)
cn = b["conserved"].notna()
gb = b["tier"].isin(["GOLD", "BRONZE"])
hum = b["species"] == "human"
cov = {
    "n_records": N,
    "defined_all_records": [int(cn.sum()), float(cn.sum() / N)],
    "defined_gold_bronze_mapped_over_all_records (analysis set, text 7.4%)": [int((gb & cn & mapped).sum()), float((gb & cn & mapped).sum() / N)],
    "defined_gold_bronze_over_gold_bronze": [int((gb & cn).sum()), float((gb & cn).sum() / gb.sum())],
    "defined_human_over_human_records (limitations 12.3%)": [int((cn & hum).sum()), float((cn & hum).sum() / hum.sum()), int(hum.sum())],
    "defined_mouse_over_mouse_records": [int((cn & (b.species == 'mouse')).sum()), float((cn & (b.species == 'mouse')).sum() / (b.species == 'mouse').sum())],
    "defined_rat_over_rat_records": [int((cn & (b.species == 'rat')).sum()), float((cn & (b.species == 'rat')).sum() / (b.species == 'rat').sum())],
    "analysis_n_in_matching_results": M["conservation"]["cem_primary"]["n_analysis"],
}
assert abs(cov["defined_gold_bronze_mapped_over_all_records (analysis set, text 7.4%)"][1] - 0.074) < 0.0005
assert abs(cov["defined_human_over_human_records (limitations 12.3%)"][1] - 0.123) < 0.0005

# ------------------------------------------------- study cluster counts: 1,113 and 1,070
ts = pd.read_parquet(REL / "temporal_split.parquet", columns=["study_cluster"])
info = json.load(open(RES / "b_cluster_rule_info.json"))
cs_all = b["primary_cluster"]
K = ["uniprot_id", "residue_type", "position"]
lost_labels = set(cs_all.dropna()) - set(ts["study_cluster"].dropna())
assert not (set(ts["study_cluster"].dropna()) - set(cs_all.dropna()))
clusters = {
    "distinct_primary_cluster_labels_all_records": int(cs_all.nunique()),
    "records_all": int(N), "records_with_label": int(cs_all.notna().sum()),
    "temporal_split_records": int(len(ts)), "distinct_primary_cluster_labels_temporal_split_records": int(ts["study_cluster"].nunique()),
    "records_not_in_temporal_split": int(N - len(ts)),
    "labels_only_on_records_outside_temporal_split": len(lost_labels),
    "distinct_cluster_ids_over_all_50904_pmids_under_rule_884": info["884"]["n_distinct_cluster_ids_over_all_pmids"],
    "candidate_pmids_884": info["884"]["n_candidates"], "multi_pmid_clusters_884": info["884"]["n_multi_pmid_clusters"],
    "candidate_clusters_884": info["884"]["n_candidate_clusters"],
    "years_resolved_not_involved": True,
}
assert clusters["distinct_primary_cluster_labels_all_records"] - clusters["distinct_primary_cluster_labels_temporal_split_records"] == clusters["labels_only_on_records_outside_temporal_split"]

out = {"table_text_reconciliation": rows, "ortholog_coverage": cov, "study_clusters": clusters,
       "cluster_rule_sensitivity": {h: {k: info[h][k] for k in ("n_candidates", "n_candidate_clusters", "n_candidates_in_multi_pmid_cluster", "n_multi_pmid_clusters", "tier_counts_study_cluster_rule", "tier_counts_depth_only_rule")} for h in ("500", "884", "1000")}}
dv = pd.read_parquet(B.PRIV / "depth_variants.parquet", columns=["depth_500", "depth_884", "depth_1000"])
out["cluster_rule_sensitivity"]["records_whose_depth_differs"] = {"500_vs_884": int((dv.depth_500 != dv.depth_884).sum()), "1000_vs_884": int((dv.depth_1000 != dv.depth_884).sum())}
json.dump(out, open(RES / "b1_reconciliation.json", "w"), indent=1, default=str)

# ------------------------------------------------- RECONCILIATION.md
def f(v, key):
    return f"{v:.3f}" if key.startswith("domain") else f"{v:.4f}"

L = []
L.append("# Reconciliation of Table 1 and the text (Rev 5), Task B1\n")
L.append("Written by `scripts/analysis/b1_reconcile.py`. Every value below is recomputed from the site table (`mammalian_sites.parquet`, `human_sites_source_set_a.parquet`, `sei_protein.parquet`) and compared with `study/earlier_results/matching_results.json`, `dataset_counts.json`, Table 1 and the manuscript text. No manuscript file was edited.\n")
L.append("## 1. Unmatched estimates: the text and Table 1 use different record sets\n")
L.append("Two record sets exist in `matching_results.json`: **all records** (`unmatched_all_sites`: every GOLD and BRONZE record with the outcome defined, whether or not its protein has a primary index) and **records with an index** (`unmatched_mapped_sites`: the same records restricted to proteins with a primary index, which is the set the matched designs use). Table 1 states that its unmatched column uses the same records as the matched columns, and it does (records with an index). The text does not do this consistently.\n")
L.append("| Comparison | Tier source; species | Text quotes | Table 1 quotes | Recomputed, all records | Recomputed, records with an index | Text uses | Table 1 uses |")
L.append("|---|---|---|---|---|---|---|---|")
for r in rows:
    a, m2 = r["all_records"], r["mapped_records"]
    L.append(f"| {r['label']} | {r['tier_source']}; {r['species']} | {r['text_value']} | {r['table_value']} | {f(a['recomputed'], r['key'])} (n GOLD {a['n_gold']:,}, BRONZE {a['n_bronze']:,}) | {f(m2['recomputed'], r['key'])} (n GOLD {m2['n_gold']:,}, BRONZE {m2['n_bronze']:,}) | {r['text_matches']} | {r['table_matches']} |")
L.append("")
L.append("All twelve recomputed values equal the result-file values to within 1e-9. Findings:\n")
mism = [r for r in rows if r["text_matches"] != r["table_matches"]]
L.append(f"- The text and Table 1 give **different unmatched values for {len(mism)} of 6 rows** ({', '.join(r['label'].split(' (')[0] for r in mism)}). In each of these the text quotes the all-records value and Table 1 the records-with-an-index value. The text's domain-overlap sentences (full population and annotated-only) quote the records-with-an-index values, so the text mixes both sets. For disorder the two values are 0.0310 (text) and 0.0438 (Table 1), the mismatch noted in the plan.")
L.append("- The matched columns are unaffected (they use records with an index in both places). The records excluded for having no primary index, 127,716 (domain overlap), 5,488 (conservation), 58,400 (disorder) and 213,479 (kinase motif), reproduce exactly (`n_unmapped_recomputed` equals the result file for every row in `results/b1_reconciliation.json`).")
L.append("- **Tier source differs between rows.** Domain overlap (both rows) used tiers computed from source set A (the human file) and is human only; the other four rows used set-B tiers over human, mouse and rat. \"GOLD\" is therefore not the same set of human records in the domain-overlap row and in the other rows. Rev 6 uses set-B tiers for every signal (decision D1); the set-A result is kept as block G15.")
L.append("- The sequence-based kinase-motif row has 46 records with an undefined outcome (non-canonical residues). 43 of them are GOLD or BRONZE records of proteins with an index. Rev 5 counted them in `n_analysis` (803,787) and in the index deciles but never in an estimate; the Rev 6 population excludes them before the deciles are formed (803,744). The effect on the adjusted estimate is 6.8e-6 (0.063056 versus 0.063063); see `results/validate_engine.json`.\n")
L.append("## 2. Ortholog coverage: 7.4% and 12.3% are different quantities\n")
L.append("| Quantity | Count | Denominator | Share |")
L.append("|---|---|---|---|")
c = cov
k1 = "defined_gold_bronze_mapped_over_all_records (analysis set, text 7.4%)"
k2 = "defined_human_over_human_records (limitations 12.3%)"
L.append(f"| Conservation indicator defined, GOLD or BRONZE, protein with an index (the 130,899 analysis records; **text, 7.4%**) | {c[k1][0]:,} | all {N:,} records | {c[k1][1]*100:.2f}% |")
L.append(f"| Conservation indicator defined, human records (**limitations, 12.3%**) | {c[k2][0]:,} | {c[k2][2]:,} human records | {c[k2][1]*100:.2f}% |")
L.append(f"| Conservation indicator defined, all records | {c['defined_all_records'][0]:,} | all {N:,} records | {c['defined_all_records'][1]*100:.2f}% |")
L.append(f"| Conservation indicator defined, GOLD or BRONZE | {c['defined_gold_bronze_over_gold_bronze'][0]:,} | {int(gb.sum()):,} GOLD or BRONZE records | {c['defined_gold_bronze_over_gold_bronze'][1]*100:.2f}% |")
L.append(f"| Defined, mouse records | {c['defined_mouse_over_mouse_records'][0]:,} | mouse records | {c['defined_mouse_over_mouse_records'][1]*100:.2f}% |")
L.append(f"| Defined, rat records | {c['defined_rat_over_rat_records'][0]:,} | rat records | {c['defined_rat_over_rat_records'][1]*100:.2f}% |")
L.append(f"\nThe text's 7.4% is the analysis set as a share of the whole dataset, so it mixes a tier restriction and an index restriction into a coverage figure. The limitations sentence's 12.3% is the human share only, and the sentence calls it \"records\" without saying human. Neither is the dataset-wide coverage of the indicator ({c['defined_all_records'][1]*100:.2f}%). The next manuscript revision should state the denominator for each.\n")
L.append("## 3. Study clusters: 1,113 and 1,070 are the same quantity on different records\n")
cl = clusters
L.append(f"- 1,113 is the number of distinct **primary study cluster** labels over all {cl['records_all']:,} records (a record's primary cluster is the cluster of its globally most frequent PMID; `cluster_split.parquet`). 1,070 is the same count over the {cl['temporal_split_records']:,} records of the temporal split, which excludes {cl['records_not_in_temporal_split']} records (the manuscript text attributes these to 46 non-canonical residues and 2 records with no linked PMID; that split is quoted, not re-derived here). {cl['labels_only_on_records_outside_temporal_split']} labels occur only on those records, and 1,113 minus 43 is 1,070. **Publication-year resolution plays no part.**")
L.append(f"- Neither is the number of study clusters under the merge rule. Under the rule (candidate PMID at >= 884 sites, Jaccard >= 0.90) the 50,904 distinct PMIDs fall into **{cl['distinct_cluster_ids_over_all_50904_pmids_under_rule_884']:,} clusters**: {cl['candidate_pmids_884']} candidate PMIDs, of which 21 sit in {cl['multi_pmid_clusters_884']} multi-PMID clusters (sizes 7, 2, 2, 2, 2, 2, 2, 2), and every other PMID is its own cluster. The sentence \"The dataset contained 1,113 study clusters\" therefore describes the number of distinct primary-cluster labels, not the number of clusters.\n")
L.append("## 4. Study-cluster rule: how 884 arose, and how little it matters\n")
L.append("- **Trace.** The code default of `build_sei_secondary_table(high_freq_min=...)` is 500 (`phosnet/sei.py:126`). The Rev 5 tiers were built with 884. `audit_reports/phase1_mammalian_ingest_report.md` (line 78) records that the human-only stage used a cutoff of 500 sites (500 / 1,071,725 = 0.047%), restated as 0.05% of the pooled count: round(0.0005 x 1,768,068) = 884. It is a proportional rule meant to keep an earlier rate. `phase0_epsd_gate_report.md` (line 216) gives the reason for the 500 cutoff itself as a judgement (a PMID below it \"cannot plausibly be a mega-consortium companion record\"); no data-derived justification for 500 is recorded in the repository. 884 was not derived from a data-quality criterion.")
sens = out["cluster_rule_sensitivity"]
L.append("- **Sensitivity (reported as found).** The three thresholds differ in the number of candidate PMIDs and candidate clusters but **not in any record's study-cluster depth or tier**:\n")
L.append("| Threshold (sites) | Candidate PMIDs | Candidate clusters | Candidates in multi-PMID clusters | Multi-PMID clusters | GOLD / SILVER / BRONZE (study-cluster rule) | GOLD / SILVER / BRONZE (depth-only rule) |")
L.append("|---|---|---|---|---|---|---|")
for h in ("500", "884", "1000"):
    s = sens[h]
    t1, t2 = s["tier_counts_study_cluster_rule"], s["tier_counts_depth_only_rule"]
    L.append(f"| {h} | {s['n_candidates']} | {s['n_candidate_clusters']} | {s['n_candidates_in_multi_pmid_cluster']} | {s['n_multi_pmid_clusters']} | {t1['GOLD']:,} / {t1['SILVER']:,} / {t1['BRONZE']:,} | {t2['GOLD']:,} / {t2['SILVER']:,} / {t2['BRONZE']:,} |")
L.append(f"\nRecords whose depth differs from the 884 rule: {sens['records_whose_depth_differs']['500_vs_884']} at 500 and {sens['records_whose_depth_differs']['1000_vs_884']} at 1,000. The same 21 PMIDs form the same 8 multi-PMID clusters at every threshold (the extra candidates between 500 and 884 sites are singletons). **Block G13 (cluster rule at 500 and 1,000) therefore reproduces G1 record for record and cannot show any dependence on the threshold.** It is kept in the manifest and reported as identical, as planned. This is a null result about the sensitivity of the rule, not evidence that the rule is correct.")
L.append("- Reproduction check: depth under 884 equals the Rev 5 `study_cluster_count` for 100% of records and the tier equals the Rev 5 tier for 100% (`results/b_cluster_rule_info.json`).\n")
L.append("## 5. Where the code departs from the text of the plan\n")
L.append("- The seven consortium PMIDs are exactly the PMID set of the second-largest cluster (`CLUSTER_38359819`); their per-record depth arrays are identical. The planned 11th removal in block G7 (`drop_consortium_7_PMIDs`) therefore duplicates `drop_cluster_2`. Both rows stay in the manifest; the results are identical by construction and the decision count counts the removal once.")
L.append("- Sites with no PMID left after a removal: 130,828 for the consortium removal, of which 13 carry a database source flag; 130,815 have no PMID and no source evidence under either reading (agreement score 0, or no PhosphoSitePlus/Phospho.ELM/dbPTM flag) and are the number dropped.")
_target = REV6 / "RECONCILIATION.md"
_tail = ""
if _target.exists():  # sections 6 and 7 are written by hand: keep everything from the "## 6." heading on
    _old = _target.read_text(encoding="utf-8").replace("\r\n", "\n")
    _i = _old.find("\n## 6.")
    if _i >= 0:
        _tail = _old[_i:]
open(_target, "w", encoding="utf-8").write("\n".join(L) + "\n" + _tail)
print("\n".join(L)[:6000])
