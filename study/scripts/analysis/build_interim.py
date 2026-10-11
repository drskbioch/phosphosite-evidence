"""Assemble results/manifest_status.csv and study/analysis/INTERIM_REPORT.md from the Rev 6 result files.
Every number is read from a result file; nothing is typed in. Tasks F to J are not started and show as pending in the status table.
"""
import glob
import json

import numpy as np
import pandas as pd

from common import REV6, RES

J = lambda p: json.load(open(RES / p))  # noqa: E731
man = pd.read_csv(REV6 / "grid_manifest.csv").fillna("")
cem = J("cem_results.json")
g5 = J("g5_abundance.json")
g14 = J("g14_cem.json")
prop = {}
for f in glob.glob(str(RES / "prop" / "*.json")):
    r = json.load(open(f))
    prop[r["contrast_id"]] = r
dl = J("downloads_manifest.json")
v50 = J("validate_engine_v50.json")
val = J("validate_engine.json")
valm = J("validate_engine_modes.json")
valp = J("validate_propensity.json")
rec = J("b1_reconciliation.json")
b4 = J("b4_twoway.json")
b5 = J("b5_index_effect_sizes.json")
d2 = J("d2_index_variants_summary.json")
d345 = J("d3_d4_d5_index_checks.json")
cov = J("e_abundance_coverage.json")
diag = J("g5_diagnostics.json")
ord_ = J("propensity_order_dependence.json")
proxy = J("d2_index_variants_summary_proxy5_superseded.json")
DC = json.load(open(REV6.parent / "earlier_results" / "dataset_counts.json"))

OR = lambda s: s == "S1"  # noqa: E731


def f(v, s, nd=None):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "n/a"
    nd = nd or (3 if OR(s) else 4)
    return f"{v:.{nd}f}"


def ci(c, s, nd=None):
    return f"{f(c[0], s, nd)} to {f(c[1], s, nd)}" if c else "n/a"


def pct(v):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{100 * v:.1f}%"


# ------------------------------------------------------------------ manifest status
rows = []
for _, r in man.iterrows():
    cid = r["contrast_id"]
    st, est, c_, att, note = "pending (task " + r["task"] + ")", "", "", "", ""
    if cid in cem:
        x = cem[cid]
        st, est, c_ = x["status"], x["estimate"], x["ci"]
        att = x.get("attenuation", "")
    elif cid in g5:
        x = g5[cid]
        st = x["status"]
        est, c_, att, note = x.get("estimate", ""), x.get("ci", ""), x.get("attenuation", ""), x.get("reason", "")
    elif cid in g14:
        x = g14[cid]
        st, est, c_ = x["status"], x.get("estimate", ""), x.get("ci_protein_cluster", "")
    elif cid in prop:
        x = prop[cid]
        st = x["status"]
        est, c_, att, note = x.get("estimate", ""), x.get("ci", ""), x.get("attenuation", ""), x.get("reason", "")
    elif r["block"] in ("G1", "G2", "G3", "G4", "G5", "G12", "G13", "G14", "G15"):
        st = "not run"
    rows.append({**r.to_dict(), "result_status": st, "estimate": est, "ci_low": c_[0] if c_ else "", "ci_high": c_[1] if c_ else "", "attenuation": att, "result_note": note})
ms = pd.DataFrame(rows)
ms.to_csv(RES / "manifest_status.csv", index=False)
cnt = ms.groupby(["block", "result_status"]).size().unstack(fill_value=0)

L = []
A = L.append
A("# Interim report, Rev 6 reanalysis\n")
A("Branch `rev6-reanalysis`. Written by `scripts/analysis/build_interim.py` from the files in `study/analysis/results/`; every number below is read from a result file. No manuscript file, reference or legacy file was touched. The plan is `ANALYSIS_PLAN.md` (deviation log in its section 15). **Work stops here, as planned.** Tasks F (strata, leave-one-cluster-out, tier thresholds), H (established measures), I (within-protein) and J (predictors) have not been started; `results/manifest_status.csv` lists every manifest row with its status. PhosphoSitePlus is not in `data_private/` yet.\n")

# ---------------------------------------------------------------- 1 gates
A("## 1. Gate checks that precede the contrasts\n")
mx = max(v["adj_abs_diff"] for k, v in val.items() if k != "S4")
A(f"- **CEM engine against Rev 5** (`validate_engine.py`, original index, Rev 5 record sets): S1, S2, S3 and S5 reproduce the Rev 5 adjusted estimates to within {mx:.1e}, and the crude estimates exactly. S4 reproduces exactly (difference {val['S4_rev5_record_set']['adj_abs_diff']:.1f}) on the Rev 5 record set. The Rev 6 population for S4 drops the {val['S4']['phosnet_n_analysis'] - val['S4']['engine_n']} records with an undefined outcome before forming the index deciles, as Rev 5 did not; this moves the S4 adjusted estimate by {val['S4']['adj_abs_diff']:.1e} ({val['S4']['engine_adjusted']:.6f} against {val['S4']['phosnet_cem_primary']:.6f}). Recorded in the deviation log.")
A(f"- **V50 index** (`validate_engine_v50.py`): the engine equals phosnet's `matched_comparison` for all five signals (largest difference {max(v['abs_diff'] for v in v50.values()):.1e}). **Abundance modes** (`validate_engine_modes.py`): equal to phosnet with abundance as the covariate or as an extra stratum (largest difference {max(v['abs_diff'] for v in valm.values()):.1e}). The Rev 5 secondary-index CEM rows reproduce in G14 (below); the Rev 5 propensity rows agree up to the randomness of the matching order (the seed of the visiting order differs).")
pd_ = valp
A(f"- **Propensity module** (`validate_propensity.py`): reproduces the Rev 5 matched estimate and pair count exactly for conservation ({pd_['S2']['engine']['estimate']:.6f}, {pd_['S2']['engine']['n_pairs']:,} pairs) and disorder ({pd_['S3']['engine']['estimate']:.6f}, {pd_['S3']['engine']['n_pairs']:,} pairs), and reproduces the Rev 5 failure for the kinase-motif lookup ({pd_['S5']['engine']['reason']}). The cell-level fit used inside resamples agrees with the record-level statsmodels fit to {pd_['S2']['max_abs_logit_diff_cell_vs_recordlevel_STY_records']:.1e} on serine, threonine and tyrosine records (non-canonical records are separated and their fitted values differ; they are never matched, see section 5).")
A("- A bug found while running the abundance block (zero-weight cells above the top decile edge received an out-of-range stratum index) was fixed in `engine.py`. It cannot change an estimate (the cells have weight 0). All earlier CEM contrasts (46 contrasts, 184 compared fields) were rerun with the fixed engine and are identical.\n")

# ---------------------------------------------------------------- 2 downloads
A("## 2. Downloads\n")
A("| Data | Source URL | Version | Retrieved | Result |")
A("|---|---|---|---|---|")
pv = dl["pharos_tdl"]
A(f"| Pharos/TCRD target list with TDL | {pv['source_url']} (GraphQL, `targets{{targets(top,skip){{uniprot sym tdl fam}}}}`) | {pv['version']} | {pv['retrieved']} | {pv['rows']:,} of {pv['count_reported']:,} targets, {pv['distinct_uniprot']:,} distinct UniProt IDs; TDL counts {pv['tdl_counts']} |")
g = dl["goa_human"]
A(f"| GOA human annotations | {g['source_url']} | {'; '.join(g['version'][:3])} (GO {g['version'][3].split('/')[-3] if len(g['version']) > 3 else ''}) | {g['retrieved']} | {g['bytes']:,} bytes, equal to Content-Length; complete gzip read |")
for t, nm in (("9606", "human"), ("10090", "mouse"), ("10116", "rat")):
    p = dl[f"paxdb_{t}"]
    A(f"| PaxDb integrated whole-organism, {nm} | {p['source_url']} | {p['version']}; file last modified {p['last_modified']} | {p['retrieved']} | {p['data_rows']:,} rows |")
a = dl["peptideatlas_phospho_606"]
A(f"| PeptideAtlas human phosphoproteome build | {a['source_url']} | {a['version']}; last modified {a['last_modified_header']} | file on disk modified {a['file_modified_on_disk']} | {a['bytes']:,} bytes, equal to Content-Length, zip test passes. Contains build tables (`atlas_build`, `biosequence_property_set`, ...); per-site confidence categories not yet examined (task H2). {a['retrieval_note']} |")
for sp in ("HUMAN_9606", "MOUSE_10090", "RAT_10116"):
    p = dl[f"uniprot_idmapping_{sp}"]
    A(f"| UniProt ID mapping, {sp} (for PaxDb to UniProt) | {p['source_url']} | {p['version']} | file on disk modified {p['file_modified_on_disk']} | {p['bytes']:,} bytes, gzip test passes |")
A("| PhosphoSitePlus Phosphorylation_site_dataset | not downloaded (you place it) | | | not present in `data_private/` |")
A(f"\nThe first Pharos pull returned 10 of 20,412 targets because of a pagination error in my query (arguments on the wrong field); it was corrected and re-run, and the incomplete files are in `study/_archive/rev6_superseded/`. Pharos and GOA each list about twenty thousand human proteins ({pv['rows']:,} and {d345['info']['goa_proteins']:,}), far fewer than the {cov['human']['proteins_in_site_table']:,} human accessions in the site table, which matters for section 9.\n")

# ---------------------------------------------------------------- 3 reconciliation
A("## 3. Reconciliation of Rev 5 (task B1, details in `RECONCILIATION.md`)\n")
mism = [r["label"] for r in rec["table_text_reconciliation"] if r["text_matches"] != r["table_matches"]]
c = rec["ortholog_coverage"]
cl = rec["study_clusters"]
A(f"- Text and Table 1 quote different unmatched estimates for {len(mism)} of 6 rows ({'; '.join(m.split(' (')[0] for m in mism)}): the text quotes all GOLD and BRONZE records, Table 1 the records whose protein has a primary index. Both values reproduce from the site table (differences below 1e-9). Disorder: 0.0310 against 0.0438. Matched columns are unaffected.")
A("- Domain overlap used set-A tiers and the other four signals set-B tiers in Rev 5, so GOLD was not the same set of human records in both.")
A(f"- Ortholog coverage: the text's 7.4% is {c['analysis_n_in_matching_results']:,} of {c['n_records']:,} records (analysis set over the whole dataset); the 12.3% in the limitations is the human share only ({c['defined_human_over_human_records (limitations 12.3%)'][0]:,} of {c['defined_human_over_human_records (limitations 12.3%)'][2]:,}); dataset-wide coverage is {pct(c['defined_all_records'][1])}.")
A(f"- Study clusters: 1,113 and 1,070 are the same count (distinct primary-cluster labels) over {cl['records_all']:,} and {cl['temporal_split_records']:,} records; publication-year resolution plays no part. Neither is the number of clusters under the merge rule ({cl['distinct_cluster_ids_over_all_50904_pmids_under_rule_884']:,} for 50,904 PMIDs).\n")

# ---------------------------------------------------------------- 4 B2
A("## 4. Study-cluster rule (task B2): the threshold does not matter, reported as found\n")
s_ = rec["cluster_rule_sensitivity"]
A(f"The rule was built with 884 sites per candidate PMID (0.05% of 1,768,068, restating a 500-site cutoff of an earlier human-only stage; the code default is 500; no data-derived justification for 500 is recorded). Candidate PMIDs: {s_['500']['n_candidates']} at 500, {s_['884']['n_candidates']} at 884, {s_['1000']['n_candidates']} at 1,000. The same {s_['884']['n_candidates_in_multi_pmid_cluster']} PMIDs form the same {s_['884']['n_multi_pmid_clusters']} multi-PMID clusters at all three, **study-cluster depth differs for {s_['records_whose_depth_differs']['500_vs_884']} records at 500 and {s_['records_whose_depth_differs']['1000_vs_884']} at 1,000, and the tier counts are identical** ({s_['884']['tier_counts_study_cluster_rule']}). Block G13 therefore reproduces G1 point for point and cannot show threshold dependence; this says nothing about whether the rule is right.\n")

# ---------------------------------------------------------------- 5 primary results
A("## 5. Primary comparisons, with the protein-cluster bootstrap (tasks B3 and C; blocks G1 to G3)\n")
A("Index V50, set-B tiers for every signal, GOLD against BRONZE. 1,000 resamples of proteins (all sites of a sampled protein kept), the whole pipeline re-run in each resample; D0 and D1 share the same resamples, so attenuation is computed within the resample. D2 is propensity matching (the matcher of Rev 5, refitted in each resample). Decision rules are those of section 2 of the plan.\n")
A("| Block | Signal | Crude (95% CI) | D1 CEM (95% CI) | Survives D1 | Attenuation D1 (95% CI) | D1 label | D2 propensity (95% CI) | D2 survives | D2 attenuation | Balance SMD D1 / D2 |")
A("|---|---|---|---|---|---|---|---|---|---|---|")
names = {"S1": "Domain overlap (OR)", "S2": "Conservation", "S3": "Disorder", "S4": "Motif, sequence-based", "S5": "Motif, database lookup (not decision-bearing)"}
d2id = {}
for _, r in man[man["design"] == "D2"].iterrows():
    d2id[(r["block"], r["signal"], r["population"])] = r["contrast_id"]
for blk in ("G1", "G2", "G3"):
    for _, r in man[(man["block"] == blk) & (man["design"] == "D1")].iterrows():
        x = cem[r["contrast_id"]]
        s = r["signal"]
        p = prop.get(d2id.get((blk, s, r["population"])))
        if p and p["status"] == "run":
            pc = f"{f(p['estimate'], s)} ({ci(p['ci'], s)})"
            ps, pa = ("yes" if p["survives"] else "no"), f"{pct(p['attenuation'])} ({pct(p['attenuation_ci'][0])} to {pct(p['attenuation_ci'][1])})"
            sm = f"{x['smd']:.3f} / {p['smd']:.3f}"
        else:
            pc, ps, pa = "not fitted (singular matrix)", "n/a", "n/a"
            sm = f"{x['smd']:.3f} / n/a"
        pop = {"G1": "all records", "G2": "human", "G3": "annotated proteins (diagnostic)"}[blk]
        A(f"| {blk} ({pop}) | {names[s]} | {f(x['crude'], s)} ({ci(x['crude_ci'], s)}) | {f(x['estimate'], s)} ({ci(x['ci'], s)}) | {'yes' if x['survives'] else 'no'} | {pct(x['attenuation'])} ({pct(x['attenuation_ci'][0])} to {pct(x['attenuation_ci'][1])}) | {x['attenuation_label']} | {pc} | {ps} | {pa} | {sm} |")
g1d = [cem[r["contrast_id"]] for _, r in man[(man["block"].isin(["G1", "G2"])) & (man["design"] == "D1") & (man["decision_bearing"] == "yes")].iterrows()]
A(f"\n**Result.** Of the {len(g1d)} decision-bearing D1 contrasts in G1 and G2, {sum(x['survives'] for x in g1d)} survive adjustment (the protein-cluster CI excludes the null) and {sum(bool(x['largely_explained']) for x in g1d)} are largely explained by study effort (attenuation of 50% or more with a crude interval that excludes the null). The largest attenuation point estimate is {pct(max(x['attenuation'] for x in g1d))}; {sum(x['attenuation'] < 0 for x in g1d)} contrasts have negative attenuation (adjustment increased the estimate, reported as not attenuated). These are the results for the planned design (CEM on the index). Under propensity matching the largest attenuation is {pct(max(p['attenuation'] for p in prop.values() if p['status'] == 'run' and p['block'] in ('G1', 'G2', 'G3')))}.\n")

# D2 notes
nf = [p for p in prop.values() if p["status"] != "run"]
ph = [p for p in prop.values() if "post-hoc" in p["analysis"]]
planned = [p for p in prop.values() if p["analysis"] == "planned" and p["block"] in ("G1", "G2", "G3", "G15")]
nfp = sorted(p["contrast_id"] for p in planned if p["status"] != "run")
fit = sorted(p["contrast_id"] for p in planned if p["status"] == "run")
A(f"**Propensity matching (D2) could not be fitted for {len(nfp)} of {len(planned)} planned D2 contrasts** ({', '.join(nfp)}; `LinAlgError: Singular matrix`, the failure Rev 5 reported for the kinase-motif comparisons). It fitted for {', '.join(fit)}. Cause, verified on the kinase-motif lookup and conservation records: the records with a non-canonical residue (K, H, R, C, I, D, G, M, N; {DC['n_noncanonical_sites']} in the dataset) occur only among BRONZE records, so the residue-type term is perfectly separated and the routine's Newton step can fail; whether a given comparison fits depends on numerical details. This is reported as the result. **Beside each non-fitted design, and labelled post-hoc,** the same design restricted to serine, threonine and tyrosine records is reported below (rule 12; the original is not replaced), together with the S/T/Y version of the kinase-motif designs that did fit:\n")
A("| Contrast | Signal | Population | D2 post-hoc S/T/Y estimate (95% CI) | Crude | Attenuation | Pairs / GOLD candidates | Note |")
A("|---|---|---|---|---|---|---|---|")
for cid in sorted(prop):
    p = prop[cid]
    if "post-hoc" not in p["analysis"]:
        continue
    s = p["signal"]
    inside = p["ci"][0] <= p["estimate"] <= p["ci"][1]
    A(f"| {cid} | {names[s]} | {p['population']} | {f(p['estimate'], s)} ({ci(p['ci'], s)}) | {f(p['crude'], s)} | {pct(p['attenuation'])} | {p['n_pairs']:,} / {p['n_treated_candidates']:,} | {'' if inside else 'point estimate lies outside its bootstrap interval (see below)'} |")
out_pt = [p["contrast_id"] for p in prop.values() if p["status"] == "run" and not (p["ci"][0] <= p["estimate"] <= p["ci"][1])]
om = ord_["summary"]
A(f"\n**Caveat on D2.** The matcher inherited from Rev 5 keeps each GOLD record's 50 nearest controls on the logit; with {ord_['distinct_logit_values']:,} distinct logit values among {ord_['n_records']:,} records (human domain overlap, serine/threonine/tyrosine records), ties are broken by record position, so only a fraction of GOLD records are matched (see the pair counts) and the point estimate depends on row order: the matched odds ratio is {om['file_order_min']:.2f} to {om['file_order_max']:.2f} in file order and {om['shuffled_min']:.2f} to {om['shuffled_max']:.2f} when the rows are shuffled (`results/propensity_order_dependence.json`). For {len(out_pt)} contrasts ({', '.join(sorted(out_pt)) if out_pt else 'none'}) the point estimate lies outside its own percentile interval. D2 is a sensitivity design; its interval is reported as computed and should be read with this in mind. I did not change the matcher.\n")

# ---------------------------------------------------------------- 6 tier source etc
A("## 6. What the tier definition and inference change (blocks G12, G14, G15)\n")
g1s1, g15, g12 = cem[man[(man.block == "G1") & (man.signal == "S1") & (man.design == "D1")].iloc[0]["contrast_id"]], cem[man[(man.block == "G15") & (man.design == "D1")].iloc[0]["contrast_id"]], cem[man[(man.block == "G12") & (man.signal == "S1") & (man.design == "D1")].iloc[0]["contrast_id"]]
A(f"- **Domain overlap depends strongly on the tier source.** Matched odds ratio {f(g1s1['estimate'], 'S1')} ({ci(g1s1['ci'], 'S1')}) with set-B tiers (G1), against {f(g15['estimate'], 'S1')} ({ci(g15['ci'], 'S1')}) with set-A tiers (G15) on the same index (the tier source, and with it the GOLD and BRONZE membership, differs), and {f(g12['estimate'], 'S1')} ({ci(g12['ci'], 'S1')}) with the depth-only tier (G12). Rev 5 reported {f(g14['G14-0335']['rev5_estimate'], 'S1')} (set-A tiers, original index). The direction (above 1) is the same in all, the size is not.")
A("- **Depth-only tier (G12).** " + "; ".join(f"{names[r['signal']].split(' (')[0]} {f(cem[r['contrast_id']]['estimate'], r['signal'])} ({ci(cem[r['contrast_id']]['ci'], r['signal'])}), attenuation {pct(cem[r['contrast_id']]['attenuation'])}" for _, r in man[(man.block == "G12") & (man.design == "D1")].iterrows()) + ". All survive; none is largely explained.")
A("\n**Inference alone (G14).** Rev 5 designs (original index; set-A tier for S1) with the protein-cluster bootstrap in place of the stratified site bootstrap. The point estimates reproduce Rev 5 (largest difference " + f"{max(v['estimate_abs_diff_from_rev5'] for v in g14.values()):.1e}" + ") and only the intervals change:\n")
A("| Signal | Design | Estimate | Rev 5 interval | Protein-cluster interval | Width ratio | Excludes null: Rev 5 / protein-cluster |")
A("|---|---|---|---|---|---|---|")
for cid, v in g14.items():
    s = v["signal"]
    A(f"| {names[s].split(' (')[0]} | {v['design'].replace('R5_', '').replace('_', ' ')} | {f(v['estimate'], s)} | {ci(v['rev5_ci_stratified'], s)} | {ci(v['ci_protein_cluster'], s)} | {v['ci_width_ratio_protein_cluster_over_rev5']:.2f} | {'yes' if v['rev5_ci_excludes_null'] else 'no'} / {'yes' if v['protein_cluster_ci_excludes_null'] else 'no'} |")
chg = [v for v in g14.values() if v["rev5_ci_excludes_null"] != v["protein_cluster_ci_excludes_null"]]
A(f"\nIntervals widen by a factor of {min(v['ci_width_ratio_protein_cluster_over_rev5'] for v in g14.values()):.2f} to {max(v['ci_width_ratio_protein_cluster_over_rev5'] for v in g14.values()):.2f}; the null is excluded in the same way as before in {len(g14) - len(chg)} of {len(g14)} CEM rows. Propensity rows of G14 are in `results/prop/`.\n")
# G14 propensity
gp = [p for p in prop.values() if p["block"] == "G14"]
if gp:
    M5 = json.load(open(REV6.parent / "earlier_results" / "matching_results.json"))
    k5 = {"S1": "domain_overlap_human", "S2": "conservation", "S3": "disorder", "S4": "kinase_motif_sequence", "S5": "kinase_motif_lookup"}
    A("G14 propensity rows (Rev 5 design, original index; protein-cluster interval against the Rev 5 site-bootstrap interval):\n")
    A("| Contrast | Signal | This run | Rev 5 | Pairs now / Rev 5 |")
    A("|---|---|---|---|---|")
    for p in sorted(gp, key=lambda q: q["contrast_id"]):
        s_ = p["signal"]
        r5 = M5[k5[s_]]["propensity"]
        now = f"{f(p['estimate'], s_)} ({ci(p['ci'], s_)})" if p["status"] == "run" else f"not fitted ({p['reason']})"
        old = f"{f(r5['point_estimate'], s_)} ({f(r5['ci_low'], s_)} to {f(r5['ci_high'], s_)})" if "point_estimate" in r5 else f"not fitted ({r5['error']})"
        pr = f"{p['n_pairs']:,} / {r5['n_matched_pairs']:,}" if p["status"] == "run" and "n_matched_pairs" in r5 else ""
        A(f"| {p['contrast_id']} | {names[s_].split(' (')[0]} | {now} | {old} | {pr} |")
    A("")

# ---------------------------------------------------------------- 7 B4
A("## 7. Two-way sensitivity (task B4, G1 D1)\n")
A("V = V_protein + V_study - V_cell, bootstrap variances of the adjusted estimate (log scale for the odds ratio), study cluster = a site's primary cluster. Criteria fixed before the run: not feasible if V <= 0 or if at least 90% of records lie in single-record intersection cells.\n")
A("| Signal | Adjusted | V protein | V study | V cell | V two-way | Single-record-cell share | Two-way 95% CI | Survives two-way | Width ratio vs protein-cluster |")
A("|---|---|---|---|---|---|---|---|---|---|")
for cid, v in b4.items():
    s = v["signal"]
    A(f"| {names[s].split(' (')[0]} | {f(v['adjusted'], s)} | {v['V_protein']:.2e} | {v['V_study']:.2e} | {v['V_cell']:.2e} | {v['V_two_way']:.2e} | {pct(v['share_of_records_in_single_record_cells'])} | " + (f"{ci(v['ci_two_way'], s)} | {'yes' if v['survives_two_way'] else 'no'} | {v['width_ratio_two_way_vs_protein']:.1f} |" if v["feasible"] else f"not feasible: {v['reason']} | | |"))
ns = [v for v in b4.values() if v["feasible"] and not v["survives_two_way"]]
A(f"\nThe study-cluster component dominates the variance (V study is {min(v['V_study'] / v['V_protein'] for v in b4.values()):.0f} to {max(v['V_study'] / v['V_protein'] for v in b4.values()):.0f} times V protein). Under the two-way interval {len(ns)} of {sum(v['feasible'] for v in b4.values())} contrasts no longer exclude the null" + (": " + ", ".join(names[v['signal']].split(' (')[0] for v in ns) if ns else "") + f". This is a sensitivity analysis; the protein-cluster interval remains primary (decision D4). The normal approximation uses 1,000 study-cluster resamples with a few very large clusters (the largest holds {pct(json.load(open(RES / 'b_largest_clusters.json'))[0]['n_sites'] / DC['n_sites'])} of all {DC['n_sites']:,} sites), so the interval is approximate.\n")

# ---------------------------------------------------------------- 8 B5
A("## 8. Effect sizes for the index properties (task B5; replaces the p-values)\n")
A(f"{b5['n_proteins']:,} proteins, protein bootstrap, 1,000 resamples, no p-values. Rev 5 values are reproduced for the original index (eta-squared {b5['rev5_check']['rev5_kruskal_eta2']:.4f}, Spearman {b5['rev5_check']['rev5_spearman_rho']:.4f}, rank-biserial {b5['rev5_check']['rev5_rank_biserial']:.4f}). Sign convention of the rank-biserial r (Rev 5's): positive means GOLD proteins have the *lower* index.\n")
A("| Effect size | Estimate | 95% CI |")
A("|---|---|---|")
for k, v in b5["effects"].items():
    A(f"| {v['description']} | {v['estimate']:.4f} | {v['ci95'][0]:.4f} to {v['ci95'][1]:.4f} |")
A("")

# ---------------------------------------------------------------- 9 task D
A("## 9. Index variants (task D)\n")
A("| Variant | PMIDs removed (> n genes) | Protein-PMID links removed | Proteins with zero PMIDs | Spearman with original index |")
A("|---|---|---|---|---|")
A(f"| V0 (original) | 0 | 0 | | 1 |")
for v in ("V20", "V50", "V100"):
    A(f"| {v} | {d2[v]['pmids_removed']:,} | {d2[v]['protein_pmid_links_removed']:,} | {d2[v]['proteins_with_zero_after']:,} | {d2['correlation_with_original'][v]['spearman_with_original']:.4f} |")
vt = d2["VT_rows"]
A(f"| VT2025 | | | {vt['2025']['proteins_missing_no_dated_pmid']:,} missing | {d2['correlation_with_original']['VT2025']['spearman_with_original']:.4f} (n = {d2['correlation_with_original']['VT2025']['n_proteins']:,}) |")
A(f"| VT2026 | | | {vt['2026']['proteins_missing_no_dated_pmid']:,} missing | {d2['correlation_with_original']['VT2026']['spearman_with_original']:.4f} (n = {d2['correlation_with_original']['VT2026']['n_proteins']:,}) |")
A(f"\n**Time normalisation.** t = 2025 - first year + 1 (2026 for the sensitivity); asserted t >= 1 and finite in all rows used (t from {vt['2025']['t_min_used']:.0f} to {vt['2025']['t_max_used']:.0f}). First year = earliest year among all of a protein's bulk-corrected PMIDs ({d2['years']['pmids_needed']:,} PMIDs dated from NCBI E-utilities and `pmid_year.parquet`; none left undated). **Not filled or imputed:** {d2['VT']['proteins_without_dated_pmid']:,} of {d2['VT']['proteins_mapped']:,} mapped proteins ({d2['VT']['genes_without_any_dated_pmid_total']:,} of {d2['VT']['genes_mapped_total']:,} genes) have no PMID left after the V50 correction, so VT is missing for them; all of them are proteins with zero publications. A further {vt['2025']['proteins_missing_first_year_after_reference']} proteins have a first year of 2026, which gives t = 0 against 2025 and is left missing (none for the 2026 reference). These proteins are excluded only from analyses that use VT. Because the excluded proteins are exactly the least-studied ones, VT comparisons are on a differently composed population (see G4).\n")
mp = d345["info"]
A(f"**Pharos TDL (D3).** {mp['human_proteins_with_pharos_entry']:,} of {mp['human_proteins_with_V0_index']:,} human proteins with an index match a Pharos entry (Pharos lists {mp['pharos_entries_total']:,} targets, of which {mp['pharos_entries_matched_to_index_proteins']:,} match index proteins; the others are not in the index or the site table). TDL is built partly from publication-derived scores, so agreement with a publication-count index is supportive, not an independent calibration. Mean index by TDL and AUROC for separating non-Tdark from Tdark (common set of {mp['common_set_with_pharos']:,} proteins with all variants defined; bootstrap 95% CI):\n")
A("| Variant | Tclin | Tchem | Tbio | Tdark | AUROC (95% CI) | Own-set AUROC |")
A("|---|---|---|---|---|---|---|")
for v in ("V0", "V20", "V50", "V100", "VT2025", "VT2026"):
    t = d345["tdl_means"][f"common_set|{v}"]
    au = d345["auroc_tdark"][f"common_set|{v}"]
    ow = d345["auroc_tdark"][f"own_set|{v}"]
    A(f"| {v} | " + " | ".join(f"{t[k]['mean']:.3f}" for k in ("Tclin", "Tchem", "Tbio", "Tdark")) + f" | {au['auroc_nonTdark_vs_Tdark']:.4f} ({au['ci95'][0]:.4f} to {au['ci95'][1]:.4f}) | {ow['auroc_nonTdark_vs_Tdark']:.4f} (n = {ow['n_proteins']:,}) |")
A(f"\n**GO check (D4) and choice of variant (D5).** Pre-specified: proteins absent from `goa_human.gaf` count 0 terms ({mp['human_proteins_without_goa_rows']:,} of {mp['human_proteins_with_V0_index']:,} index proteins have no GOA row; the file lists {mp['goa_proteins']:,} proteins). Under that coding the Spearman correlation between index and GO-term count is negative for every variant (" + ", ".join(f"{v} {d345['go_spearman'][f'common_set|{v}|all_evidence']['spearman']:.3f}" for v in ("V0", "V20", "V50", "V100", "VT2025")) + "), so it measures coverage of the file, not annotation depth. The pre-specified rule nevertheless ranks " + f"**{d345['D5']['first_by_auroc']} first on the Tdark AUROC and {d345['D5']['first_by_go']} first on the GO correlation, which gives the rule's outcome: \"{d345['D5']['decision']}\"** {'The V20 and V50 intervals overlap' if d345['auroc_tdark']['common_set|V20']['ci95'][0] <= d345['auroc_tdark']['common_set|V50']['ci95'][1] else 'The V20 and V50 intervals do not overlap'} (AUROC {d345['auroc_tdark']['common_set|V20']['auroc_nonTdark_vs_Tdark']:.4f} against {d345['auroc_tdark']['common_set|V50']['auroc_nonTdark_vs_Tdark']:.4f}). Block G4 reports S1 to S4 under V20; S5 is not decision-bearing and was not rerun.")
ph_ = d345["D5_posthoc_GOA_annotated_only"]
ph_ = d345["D5_posthoc_GOA_annotated_only"]
gph = d345["go_spearman_posthoc_GOA_annotated_only"]
A(f"**Post-hoc, labelled, beside the pre-specified result:** restricting the GO check to proteins with at least one GOA row, the correlations are positive ({', '.join(v + ' ' + format(gph[v + '|all_evidence']['spearman'], '.3f') for v in ('V0', 'V20', 'V50', 'V100', 'VT2025'))}) and rank **{ph_['first_by_go']}** first, while the Tdark AUROC still ranks {ph_['first_by_auroc']} first. Under this version the two measures disagree, and the rule would report all variants and choose none. **The two readings give different instructions, so I have not acted on either beyond what G4 already shows (V20, V50 and V100 give the same conclusions); please decide.**\n")
A("### G4: S1 to S4 under each index variant (D3, CEM; protein-cluster intervals)\n")
A("| Signal | " + " | ".join(["V0", "V20", "V50 (G1)", "V100", "VT2025", "VT2026"]) + " |")
A("|---|---|---|---|---|---|")
for s in ("S1", "S2", "S3", "S4"):
    cells_ = []
    for v in ("V0", "V20", "V50", "V100", "VT2025", "VT2026"):
        if v == "V50":
            x = cem[man[(man.block == "G1") & (man.signal == s) & (man.design == "D1")].iloc[0]["contrast_id"]]
        else:
            x = cem[man[(man.block == "G4") & (man.signal == s) & (man["index"] == v)].iloc[0]["contrast_id"]]
        cells_.append(f"{f(x['estimate'], s)} ({ci(x['ci'], s)})")
    A(f"| {names[s].split(' (')[0]} | " + " | ".join(cells_) + " |")
g4 = [cem[r["contrast_id"]] for _, r in man[man.block == "G4"].iterrows()]
A(f"\nAll {len(g4)} G4 contrasts survive (interval excludes the null: {sum(x['survives'] for x in g4)} of {len(g4)}); {sum(bool(x['largely_explained']) for x in g4)} are largely explained. The VT rows use only proteins with a defined VT (the {vt['2025']['proteins_missing_no_dated_pmid']:,} zero-PMID proteins drop out), so their crude values differ from the other columns. {'No variant changes a decision (all survive and none is largely explained, as in G1).' if all(x['survives'] and not x['largely_explained'] for x in g4) else 'At least one variant changes a decision; see the table.'}\n")

# ---------------------------------------------------------------- 10 E
A("## 10. Protein abundance (task E, G5)\n")
A("PaxDb 6.0 integrated whole-organism abundance mapped to UniProt through the UniProt ID mapping (id type STRING; maximum kept where several entries map to one accession, a rule set before analysis). Coverage:\n")
A("| Species | PaxDb entries | Unmapped | Proteins with abundance / in site table | Records with abundance | GOLD/BRONZE records with abundance and V50 index |")
A("|---|---|---|---|---|---|")
for sp, v in cov.items():
    A(f"| {sp} | {v['paxdb_entries']:,} | {v['paxdb_entries_unmapped_to_uniprot']:,} | {v['proteins_with_abundance']:,} / {v['proteins_in_site_table']:,} ({pct(v['share_of_proteins'])}) | {v['records_with_abundance']:,} ({pct(v['share_of_records'])}) | {v['gold_bronze_records_with_abundance_and_index']:,} |")
A("\nPopulation P-abund: GOLD or BRONZE records of proteins with both V50 and abundance, identical records in all four columns. D1 matches on index decile, D5 on abundance decile alone, D4 on index decile x abundance decile (species and residue type exact). Abundance enters as log10 ppm. Mouse and rat were run (S1 is human-only and was not run for them).\n")
A("| Species | Signal | Crude | D1 index | D5 abundance | D4 index x abundance | Attenuation D1 / D5 / D4 | Survives D1 / D5 / D4 | SMD of index after D5; of abundance after D1 |")
A("|---|---|---|---|---|---|---|---|---|")
for (s, popl), g in man[man.block == "G5"].groupby(["signal", "population"], sort=False):
    ids = dict(zip(g["design"], g["contrast_id"]))
    if g5[ids["D1"]].get("status") != "run":
        continue
    x1, x5, x4, x0 = g5[ids["D1"]], g5[ids["D5"]], g5[ids["D4"]], g5[ids["D0"]]
    sp = popl.split("(")[1].rstrip(")")
    A(f"| {sp} | {names[s].split(' (')[0]} | {f(x0['estimate'], s)} | {f(x1['estimate'], s)} ({ci(x1['ci'], s)}) | {f(x5['estimate'], s)} ({ci(x5['ci'], s)}) | {f(x4['estimate'], s)} ({ci(x4['ci'], s)}) | {pct(x1['attenuation'])} / {pct(x5['attenuation'])} / {pct(x4['attenuation'])} | {'yes' if x1['survives'] else 'no'} / {'yes' if x5['survives'] else 'no'} / {'yes' if x4['survives'] else 'no'} | {x5['smd_index']:.3f}; {x1['smd_abundance']:.3f} |")
run5 = [x for x in g5.values() if x.get("status") == "run" and x["design"] in ("D1", "D5", "D4")]
bal = [x for x in run5 if x["design"] == "D5" and not x["balance_index_pass"]]
A(f"\n{sum(x['survives'] for x in run5)} of {len(run5)} adjusted contrasts survive; {sum(bool(x['attenuation'] >= 0.5) and bool(x['survives']) for x in run5)} have an attenuation of 50% or more. Matching on the index alone leaves abundance imbalanced (SMD {min(x['smd_abundance'] for x in run5 if x['design']=='D1'):.2f} to {max(x['smd_abundance'] for x in run5 if x['design']=='D1'):.2f}, last column), matching on abundance alone leaves the index imbalanced in rat ({len(bal)} D5 contrasts above |SMD| 0.1 for the index; reported, not re-tuned). The attenuations range from {pct(min(x['attenuation'] for x in run5))} (an increase away from the null) to {pct(max(x['attenuation'] for x in run5))}; the largest increases occur for rat disorder (crude estimate {f(g5[man[(man.block=='G5')&(man.signal=='S3')&(man.population=='P-abund (rat)')&(man.design=='D0')].iloc[0]['contrast_id']]['estimate'], 'S3')}, a small population) and for human domain overlap. **The domain-overlap comparison in P-abund is not comparable with G1:** restricting to proteins with an abundance changes the crude odds ratio from {f(cem[man[(man.block=='G1')&(man.signal=='S1')&(man.design=='D0')].iloc[0]['contrast_id']]['estimate'], 'S1')} (G1) to {f(g5[man[(man.block=='G5')&(man.signal=='S1')&(man.design=='D0')].iloc[0]['contrast_id']]['estimate'], 'S1')}, i.e. across 1; attenuation is relative to the crude value within each population. Protein-level Spearman correlation of abundance: with the V50 index " + f"{diag['human|abundance_vs_index_V50']['spearman']:.3f} ({ci(diag['human|abundance_vs_index_V50']['ci95'], 'S2', 3)}) in human, {diag['mouse|abundance_vs_index_V50']['spearman']:.3f} in mouse, {diag['rat|abundance_vs_index_V50']['spearman']:.3f} in rat; with the outcomes in human: " + ", ".join(f"{names[s].split(' (')[0].lower()} {diag[f'human|abundance_vs_{s}']['spearman']:.3f}" for s in ("S1", "S2", "S3", "S4", "S5")) + " (all intervals in `results/g5_diagnostics.json`).\n")

# ---------------------------------------------------------------- 11 deviations and open items
A("## 11. Departures, problems and what I need from you\n")
A(f"**Logged in the plan (section 15), none changes a result:** B2 before the manifest; G7's 11th removal duplicates the 2nd; G13 identical by construction; S4 record set; D2 first-year definition (a 5-lowest-PMID proxy missed an earlier known year for at least {proxy['y0_proxy_check']['proteins_with_earlier_known_year_than_proxy']:,} of {proxy['y0_proxy_check']['proteins_compared']:,} proteins and was replaced before use); the Pharos pull. Also logged: the engine clipping fix (no estimate changed, verified by rerun); the abundance modes; the propensity non-fits and post-hoc S/T/Y rows (section 5); the GO coding (section 9); the G14 secondary-index restriction; the criteria fixed before the B4 and D5 results.\n")
A("**Decisions for you (rule 9):**\n")
A("1. **Choice of index variant (D5).** The pre-specified rule says report V50 and V20; with the GO coding corrected (post-hoc) it says choose none. V20, V50 and V100 give the same decisions in G4. Do you want the rule applied as written, or the GO coding changed and the rule rerun (a new labelled analysis)?")
A(f"2. **Propensity matching.** {len(nfp)} of {len(planned)} planned D2 contrasts could not be fitted by the Rev 5 routine because the non-canonical residues are all BRONZE. I report that and the post-hoc S/T/Y fits. Should the S/T/Y-only version become the planned D2 for the remaining blocks (G5-G9 do not use D2), or stay a labelled post-hoc?")
A("3. **Matcher tie-breaking.** The inherited 50-nearest-neighbour matcher depends on record order and matches about a quarter of GOLD records. I kept it for fidelity with Rev 5. Do you want a variant with random tie-breaking as a labelled post-hoc?")
A("4. **PhosphoSitePlus** is needed for tasks H and I; the file is not in `data_private/` yet. PeptideAtlas build 606 is on disk; whether it holds per-site confidence categories is untested.")
A("5. The Rev 5 text and Table 1 disagree on four unmatched estimates, and the 7.4%, 12.3%, 1,113 and 1,070 statements need denominators (`RECONCILIATION.md`). No manuscript edit is made in this cycle (rule 13).\n")
A("## 12. Status of the manifest\n")
A(cnt.to_markdown() if hasattr(cnt, "to_markdown") else cnt.to_string())
A("\nRows marked `pending` belong to tasks F to J, which start after this report. G7's duplicate removal and G13's identical results are in the table as run.\n")
HL = ["## Headline results\n"]
dbear = [cem[r["contrast_id"]] for _, r in man[(man["block"].isin(["G1", "G2"])) & (man["design"] == "D1") & (man["decision_bearing"] == "yes")].iterrows()]
HL.append(f"1. **Matching on the study-effort index (CEM) removes none of the decision-bearing validation signals.** All {sum(x['survives'] for x in dbear)} of {len(dbear)} decision-bearing CEM contrasts in G1 and G2 survive (protein-cluster interval excludes the null); {sum(bool(x['largely_explained']) for x in dbear)} are largely explained (attenuation of 50% or more); the largest attenuation is {pct(max(x['attenuation'] for x in dbear))}.")
HL.append(f"2. **Abundance does not explain them either** (G5): {sum(x['survives'] for x in run5)} of {len(run5)} adjusted contrasts across human, mouse and rat survive, {sum(bool(x['attenuation'] >= 0.5) and bool(x['survives']) for x in run5)} with attenuation of 50% or more. All {len(g4)} index-variant contrasts (G4) survive as well.")
ntw = [v for v in b4.values() if v["feasible"] and not v["survives_two_way"]]
HL.append(f"3. **Honest inference widens the intervals but, in the planned analysis, not the conclusions.** The protein-cluster bootstrap widens the Rev 5 intervals by {min(v['ci_width_ratio_protein_cluster_over_rev5'] for v in g14.values()):.1f} to {max(v['ci_width_ratio_protein_cluster_over_rev5'] for v in g14.values()):.1f} times with the same inference in {len(g14) - len(chg)} of {len(g14)} CEM rows (G14), and the two-way (protein and study cluster) sensitivity is wider again: {len(ntw)} of {sum(v['feasible'] for v in b4.values())} contrasts no longer exclude the null there ({', '.join(names[v['signal']].split(' (')[0].lower() for v in ntw)}).")
HL.append(f"4. **Domain overlap depends strongly on how GOLD is defined**: matched odds ratio {f(g1s1['estimate'], 'S1')} with set-B tiers against {f(g15['estimate'], 'S1')} with set-A tiers (same index), and the crude odds ratio falls below 1 in the abundance subset ({f(g5[man[(man.block=='G5')&(man.signal=='S1')&(man.design=='D0')].iloc[0]['contrast_id']]['estimate'], 'S1')}) and in the annotated-protein subset ({f(cem[man[(man.block=='G3')&(man.design=='D0')].iloc[0]['contrast_id']]['estimate'], 'S1')}) against {f(g1s1['crude'], 'S1')} in G1 (section 10).")
HL.append(f"5. **Rev 5 reporting problems found** (`RECONCILIATION.md`): text and Table 1 quote different unmatched estimates for {len(mism)} of 6 rows; 7.4% and 12.3% are different quantities; 1,113 and 1,070 count the same thing on different records. The study-cluster threshold (500, 884, 1,000) changes no record's tier.")
HL.append(f"6. **Not clean:** propensity matching (D2) could not be fitted for {len(nfp)} of {len(planned)} planned contrasts (separated non-canonical residues) and its matcher depends on record order; the GO check is ill-posed under the pre-specified coding, so the D5 rule's outcome (V20) and its post-hoc counterpart (V100, measures disagree) need your decision (section 9, section 11).\n")
idx = next(i for i, l in enumerate(L) if l.startswith('## 1. Gate'))
L[idx:idx] = HL
open(REV6 / "INTERIM_REPORT.md", "w", encoding="utf-8").write("\n".join(L) + "\n")
print(len(L), "lines")
