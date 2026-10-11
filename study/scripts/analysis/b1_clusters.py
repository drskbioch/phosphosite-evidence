"""Task B 2-3 and inputs for G7, G8, G12, G13: study clusters under several rules.

Outputs
  results/b_cluster_rule_info.json     cluster statistics for hf_min 500, 884, 1000 and the 884 reproduction check
  results/b_largest_clusters.json      the ten largest primary clusters and their PMIDs
  data_private/depth_variants.parquet  per-record depth under hf_min 500/884/1000 and with each removal
"""
import json

import numpy as np
import pandas as pd

from common import PRIV, REL, RES, assign_tier, cluster_map, depth_per_site, split_pmids

K = ["uniprot_id", "residue_type", "position"]
m = pd.read_parquet(REL / "mammalian_sites.parquet",
                    columns=K + ["species", "DAS", "study_cluster_count", "tier", "pmid_list_str"])
cs = pd.read_parquet(REL / "cluster_split.parquet", columns=K + ["study_cluster"])
assert (cs[K].values == m[K].values).all()
L = split_pmids(m["pmid_list_str"])
info, dv = {}, pd.DataFrame(index=m.index)

maps = {}
for hf in (500, 884, 1000):
    cid, i, freq = cluster_map(L, hf)
    maps[hf] = cid
    d = depth_per_site(L, cid)
    dv[f"depth_{hf}"] = d
    i = dict(i)
    i["n_distinct_cluster_ids_over_all_pmids"] = len({cid.get(p, p) for p in freq})
    i["size_of_largest_multi_pmid_cluster"] = i["multi_cluster_sizes"][0] if i["multi_cluster_sizes"] else None
    i["multi_cluster_sizes"] = i["multi_cluster_sizes"][:20]
    for rule, kw in (("study_cluster_rule", {}), ("depth_only_rule", {"use_das": False})):
        t = assign_tier(m["DAS"].values, d, **kw)
        dv[f"tier_{rule}_{hf}"] = t
        i[f"tier_counts_{rule}"] = pd.Series(t).value_counts().to_dict()
    info[str(hf)] = i
    print(hf, {k: v for k, v in i.items() if k not in ("multi_cluster_sizes",)}, flush=True)

info["reproduction_of_rev5_at_884"] = {
    "depth_equals_study_cluster_count": float((dv["depth_884"].values == m["study_cluster_count"].values).mean()),
    "tier_equals_rev5_tier": float((dv["tier_study_cluster_rule_884"].values == m["tier"].values).mean())}

# primary cluster of each site (label of the cluster of its globally most frequent PMID), as in cluster_split.parquet
pc = cs["study_cluster"]
vc = pc.value_counts()
info["primary_cluster_labels"] = {"n_distinct": int(pc.nunique()), "top10_share_of_sites": float(vc.head(10).sum() / len(pc)),
                                  "largest_share_of_sites": float(vc.iloc[0] / len(pc))}
cid = maps[884]
members = {}
for p, c in cid.items():
    members.setdefault("CLUSTER_" + c[2:], []).append(p)
top = []
for lab, n in vc.head(10).items():
    pm = members.get(lab, [lab])
    top.append({"label": lab, "n_sites": int(n), "pmids": sorted(pm)})
json.dump(top, open(RES / "b_largest_clusters.json", "w"), indent=1)
info["largest_cluster_pmids_are_labelled_CLUSTER_or_singleton"] = [t["label"][:8] for t in top]

# leave-one-out depth: drop the PMIDs of each of the ten largest primary clusters (884 rule)
for k, t in enumerate(top, 1):
    dv[f"depth_drop_{k}"] = depth_per_site(L, cid, drop=frozenset(t["pmids"]))
    n_lost = int(((dv[f"depth_drop_{k}"] == 0) & (m["DAS"] == 0)).sum())
    info[f"drop_{k}"] = {"label": t["label"], "n_pmids": len(t["pmids"]), "sites_with_no_pmid_left": int((dv[f"depth_drop_{k}"] == 0).sum()),
                         "sites_dropped_no_pmid_and_no_agreement": n_lost}
dv.to_parquet(PRIV / "depth_variants.parquet")
json.dump(info, open(RES / "b_cluster_rule_info.json", "w"), indent=1, default=str)
print("done")
