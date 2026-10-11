"""Descriptive counts for the Rev 6 figures: counts, shares, means, medians, percentiles and histogram bins. No contrast, no bootstrap, no interval.

Writes results/figure_counts.json with the input file names and SHA-256 hashes, and a cross-check against earlier_results/dataset_counts.json (study/earlier_results) for every
count that file also holds. Inputs are record-level data files (data_private/ and processed_data/epsd/release/); none is a result file.
Usage: python count_descriptives.py
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

CUR = Path(__file__).resolve().parents[2]
PRIV = CUR / "analysis" / "data_private"
REL = CUR / "processed_data" / "epsd" / "release"
OUT = CUR / "analysis" / "results" / "figure_counts.json"
COUNTS = json.load(open(CUR / "earlier_results" / "dataset_counts.json", encoding="utf-8"))
INPUTS = [PRIV / "base.parquet", PRIV / "g23_tiers.parquet", PRIV / "depth_variants.parquet", PRIV / "index_variants.parquet",
          REL / "human_sites_source_set_a.parquet", REL / "site_first_year.parquet", REL / "pmid_year.parquet",
          REL / "cluster_split.parquet", REL / "temporal_split.parquet"]
TIERS = ["GOLD", "SILVER", "BRONZE"]
KEY = ["uniprot_id", "residue_type", "position", "species"]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def box(v):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    q = np.percentile(v, [5, 25, 50, 75, 95])
    return {"n": int(len(v)), "mean": float(v.mean()), "p5": float(q[0]), "q1": float(q[1]), "median": float(q[2]), "q3": float(q[3]), "p95": float(q[4])}


def vc(s):
    return {str(k): int(v) for k, v in s.value_counts().items()}


out = {"inputs": {str(p.relative_to(CUR)).replace("\\", "/"): {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in INPUTS},
       "note": "descriptive counts only; see the docstring of scripts/analysis/count_descriptives.py"}
base = pd.read_parquet(PRIV / "base.parquet", columns=KEY + ["DAS", "DAS_setA", "tier", "tier_setA", "pmid_list_str", "primary_cluster", "disorder"])
g23 = pd.read_parquet(PRIV / "g23_tiers.parquet")
dv = pd.read_parquet(PRIV / "depth_variants.parquet", columns=["depth_884"])
assert len(base) == len(g23) == len(dv) == COUNTS["n_sites"]
base["tier_g23"] = g23["tier_g23"].values
base["depth_884"] = dv["depth_884"].values
base["depth_g23"] = g23["depth_g23"].values
human = base[base["species"] == "human"]

# ---- tiers
out["tier_counts_all"] = vc(base["tier"])
out["tier_counts_by_species"] = {s: vc(base.loc[base["species"] == s, "tier"]) for s in ["human", "mouse", "rat"]}
out["n_sites_by_species"] = {s: int((base["species"] == s).sum()) for s in ["human", "mouse", "rat"]}
rel_h = pd.read_parquet(REL / "human_sites_source_set_a.parquet", columns=["tier", "tier_pmid_rule", "DAS"])
out["human_tier_study_cluster_rule_set_b"] = vc(human["tier"])
out["human_tier_study_cluster_rule_set_a"] = vc(rel_h["tier"])
out["human_tier_pmid_count_rule"] = vc(rel_h["tier_pmid_rule"])
out["human_n_sites"] = int(len(human))
out["das_counts_all"] = vc(base["DAS"])
out["human_das_ge2_set_b"] = {"n": int((human["DAS"] >= 2).sum()), "share": float((human["DAS"] >= 2).mean())}
out["human_das_ge2_set_a"] = {"n": int((rel_h["DAS"] >= 2).sum()), "share": float((rel_h["DAS"] >= 2).mean())}
out["human_das_ge2_set_a_from_base"] = {"n": int((human["DAS_setA"] >= 2).sum())}

# ---- depth with and without the non-primary PMIDs (G23)
out["depth_histogram_primary"] = vc(base["depth_884"])
out["depth_histogram_g23"] = vc(base["depth_g23"].dropna())
out["depth_g23_n_missing"] = int(base["depth_g23"].isna().sum())
out["depth_median_primary"] = float(base["depth_884"].median())

# ---- first-reported year by tier (primary and corrected tiers)
fy = pd.read_parquet(REL / "site_first_year.parquet", columns=["site_id"] + KEY + ["evidence_tier", "first_reported_year"])
# site_id is the 1-based row number of the site table; base.parquet has the same rows in the same order (checked below)
assert (base[KEY].values[fy["site_id"].values - 1] == fy[KEY].values).all()
m = base.iloc[fy["site_id"].values - 1][["tier", "tier_g23"]].reset_index(drop=True)
m["first_reported_year"] = fy["first_reported_year"].values
assert (m["tier"].values == fy["evidence_tier"].values).all()
m = m.dropna(subset=["first_reported_year"])
out["first_year_by_tier_primary"] = {t: box(m.loc[m["tier"] == t, "first_reported_year"]) for t in TIERS}
out["first_year_by_tier_g23"] = {t: box(m.loc[m["tier_g23"] == t, "first_reported_year"]) for t in TIERS}
out["first_year_n_sites"] = int(len(m))

# ---- PubMed identifiers and residues
pm = pd.Series(sorted({p for s in base["pmid_list_str"].dropna() for p in s.split(";") if p}))
resolved = set(pd.read_parquet(REL / "pmid_year.parquet")["pmid"].astype(str))
num = pm.astype("int64")
isres = pm.isin(resolved)
out["pmid"] = {"distinct": int(len(pm)), "resolved": int(isres.sum()), "unresolved": int((~isres).sum()),
               "invalid_format_over_45m": int((num[~isres] > 45_000_000).sum()), "valid_format_not_found": int((num[~isres] <= 45_000_000).sum()),
               "unresolved_digit_lengths": vc(pm[~isres].str.len())}
out["noncanonical_by_species"] = vc(base.loc[~base["residue_type"].isin(["S", "T", "Y"]), "species"])
out["noncanonical_n"] = int((~base["residue_type"].isin(["S", "T", "Y"])).sum())

# ---- index by species and by highest tier (V50 and secondary), density bins
iv = pd.read_parquet(PRIV / "index_variants.parquet", columns=["uniprot_id", "species", "mapped_to_gene", "SEI_V50", "SEI_secondary"])
mp = iv[iv["mapped_to_gene"]].dropna(subset=["SEI_V50"])
out["v50_by_species"] = {s: box(mp.loc[mp["species"] == s, "SEI_V50"]) for s in ["human", "mouse", "rat"]}
rank = {"GOLD": 3, "SILVER": 2, "BRONZE": 1}
top = base.assign(r=base["tier"].map(rank)).sort_values("r").groupby("uniprot_id").tail(1)[["uniprot_id", "tier"]]
mt = mp.merge(top, on="uniprot_id")
out["v50_by_highest_tier"] = {t: box(mt.loc[mt["tier"] == t, "SEI_V50"]) for t in TIERS}
ms = mt.dropna(subset=["SEI_secondary"])
out["secondary_by_highest_tier"] = {t: box(ms.loc[ms["tier"] == t, "SEI_secondary"]) for t in TIERS}
both = mp.dropna(subset=["SEI_secondary"])
xe, ye = np.linspace(0, 4.5, 46), np.linspace(0, 3.0, 31)
h2, _, _ = np.histogram2d(both["SEI_V50"], both["SEI_secondary"], bins=[xe, ye])
out["v50_vs_secondary_bins"] = {"x_edges": xe.tolist(), "y_edges": ye.tolist(), "counts": h2.astype(int).tolist(), "n_proteins": int(len(both)),
                                "n_outside_range": int(len(both) - h2.sum())}

# ---- partitions
cs = pd.read_parquet(REL / "cluster_split.parquet", columns=["uniprot_id", "split_cluster_disjoint"])
sp = cs.groupby("uniprot_id")["split_cluster_disjoint"].agg(lambda s: set(s))
out["cluster_split"] = {"train": int((cs["split_cluster_disjoint"] == "train").sum()), "test": int((cs["split_cluster_disjoint"] == "test").sum()),
                        "proteins_total": int(len(sp)), "proteins_in_both": int(sum(1 for v in sp if {"train", "test"} <= v))}
ts = pd.read_parquet(REL / "temporal_split.parquet", columns=["uniprot_id", "temporal_split", "study_cluster"])
spt = ts.groupby("uniprot_id")["temporal_split"].agg(lambda s: set(s))
spc = ts.groupby("study_cluster")["temporal_split"].agg(lambda s: set(s))
out["temporal_split"] = {"train": int((ts["temporal_split"] == "train").sum()), "test": int((ts["temporal_split"] == "test").sum()),
                         "proteins_total": int(len(sp)), "proteins_in_both": int(sum(1 for v in spt if {"train", "test"} <= v)),
                         "study_clusters_total": int(len(spc)), "study_clusters_in_both": int(sum(1 for v in spc if {"train", "test"} <= v))}

# ---- disorder by tier and source set
dis = {"source_set_a_human": {t: {"n": int((human["tier_setA"] == t).sum()), "mean": float(human.loc[human["tier_setA"] == t, "disorder"].mean())} for t in TIERS}}
for s in ["human", "mouse", "rat"]:
    d = base[base["species"] == s]
    dis[f"source_set_b_{s}"] = {t: {"n": int((d["tier"] == t).sum()), "mean": float(d.loc[d["tier"] == t, "disorder"].mean())} for t in TIERS}
out["disorder_mean_by_tier"] = dis

# ---- cross-check against study/earlier_results/dataset_counts.json
cc = []


def chk(name, mine, theirs):
    cc.append({"count": name, "recomputed": mine, "dataset_counts": theirs, "equal": bool(mine == theirs)})


def as_int(d):
    return {k: int(v) for k, v in d.items()}


chk("tier_counts_mammalian", out["tier_counts_all"], as_int(COUNTS["tier_counts_mammalian"]))
chk("das_counts", out["das_counts_all"], {k: int(v) for k, v in COUNTS["das_counts"].items()})
chk("human_tier_raw_pmid_rule", out["human_tier_pmid_count_rule"], as_int(COUNTS["human_tier_raw_pmid_rule"]))
chk("human_tier_study_cluster_rule (source set A human file)", out["human_tier_study_cluster_rule_set_a"], as_int(COUNTS["human_tier_study_cluster_rule"]))
chk("human_das_ge2_share_source_set_a", out["human_das_ge2_set_a"]["share"], COUNTS["human_das_ge2_share_source_set_a"])
chk("human_das_ge2_share_source_set_b", out["human_das_ge2_set_b"]["share"], COUNTS["human_das_ge2_share_source_set_b"])
chk("human_das_ge2_set_a_n_base_vs_release", out["human_das_ge2_set_a_from_base"]["n"], out["human_das_ge2_set_a"]["n"])
chk("sites_by_species", out["n_sites_by_species"], COUNTS["sites_by_species"])
chk("first_year_median_by_tier", {t: out["first_year_by_tier_primary"][t]["median"] for t in TIERS}, {t: COUNTS["first_year_median_by_tier"][t] for t in TIERS})
chk("n_distinct_pmids", out["pmid"]["distinct"], COUNTS["n_distinct_pmids"])
chk("n_pmids_resolved", out["pmid"]["resolved"], COUNTS["n_pmids_resolved"])
chk("n_pmids_over_45m", out["pmid"]["invalid_format_over_45m"], COUNTS["n_pmids_over_45m"])
chk("n_pmids_unresolved_valid_format", out["pmid"]["valid_format_not_found"], COUNTS["n_pmids_unresolved_valid_format"])
chk("unparseable_digit_lengths", out["pmid"]["unresolved_digit_lengths"], {k: int(v) for k, v in COUNTS["unparseable_digit_lengths"].items() if int(v) > 0})
chk("noncanonical_by_species", out["noncanonical_by_species"], COUNTS["noncanonical_by_species"])
chk("n_noncanonical_sites", out["noncanonical_n"], COUNTS["n_noncanonical_sites"])
chk("cluster_split_proteins_in_both", out["cluster_split"]["proteins_in_both"], COUNTS["cluster_split_proteins_in_both"])
chk("cluster_split_proteins_total", out["cluster_split"]["proteins_total"], COUNTS["cluster_split_proteins_total"])
chk("temporal_split_counts", {"train": out["temporal_split"]["train"], "test": out["temporal_split"]["test"]}, COUNTS["temporal_split_counts"])
chk("temporal_proteins_in_both", out["temporal_split"]["proteins_in_both"], COUNTS["temporal_proteins_in_both"])
chk("temporal_clusters_in_both", out["temporal_split"]["study_clusters_in_both"], COUNTS["temporal_clusters_in_both"])
chk("temporal_clusters_total", out["temporal_split"]["study_clusters_total"], COUNTS["temporal_clusters_total"])
chk("disorder_mean_source_set_a_human", {t: round(dis["source_set_a_human"][t]["mean"], 9) for t in TIERS}, {t: round(COUNTS["disorder_mean_source_set_a_human"][t], 9) for t in TIERS})
chk("disorder_mean_source_set_b_human", {t: round(dis["source_set_b_human"][t]["mean"], 9) for t in TIERS}, {t: round(COUNTS["disorder_mean_source_set_b"]["human"][t], 9) for t in TIERS})
b5 = json.load(open(CUR / "analysis" / "results" / "b5_index_effect_sizes.json", encoding="utf-8"))
chk("b5_n_by_top_tier", {t: out["v50_by_highest_tier"][t]["n"] for t in TIERS}, {t: b5["n_by_top_tier"][t] for t in TIERS})
chk("b5_means_V50_by_top_tier", {t: round(out["v50_by_highest_tier"][t]["mean"], 9) for t in TIERS}, {t: round(b5["means_by_top_tier"]["SEI_V50"][t], 9) for t in TIERS})
out["cross_check"] = cc
out["cross_check_summary"] = {"checked": len(cc), "equal": sum(c["equal"] for c in cc), "different": [c["count"] for c in cc if not c["equal"]]}
OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
print(out["cross_check_summary"])
for c in cc:
    if not c["equal"]:
        print("DIFFERENT:", c)
