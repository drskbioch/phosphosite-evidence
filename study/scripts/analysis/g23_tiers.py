"""G23 steps 2, 3 and 5 (definition: ANALYSIS_PLAN.md, sixth batch; list: study/analysis/nonprimary_pmid_list.md, committed before this script was run).
Non-primary PMIDs are deleted from every record's publication list; the study-cluster assignment of the remaining PMIDs is kept; depth is recomputed (depth_per_site);
tiers by the study-cluster rule (GOLD: DAS >= 2 or depth >= 10; SILVER: DAS >= 1 or depth 3-9; BRONZE otherwise); DAS unchanged; records left with no PMID and no
database source flag are dropped (tier missing). First checks that the primary tiers are reproduced when nothing is removed.
Outputs: data_private/g23_tiers.parquet (tier_g23, row-aligned with mammalian_sites.parquet), results/g23_tiers_summary.json."""
import json

import numpy as np
import pandas as pd

from common import PRIV, REL, RES, assign_tier, cluster_map, depth_per_site, split_pmids

K = ["uniprot_id", "residue_type", "position"]
m = pd.read_parquet(REL / "mammalian_sites.parquet", columns=K + ["species", "DAS", "study_cluster_count", "tier", "pmid_list_str", "source_PSP", "source_PhosphoELM", "source_dbPTM"])
L = split_pmids(m["pmid_list_str"])
cid, info, freq = cluster_map(L, 884)
das = m["DAS"].to_numpy()
d0 = depth_per_site(L, cid)
t0 = assign_tier(das, d0)
assert (d0 == m["study_cluster_count"].to_numpy()).all() and (t0 == m["tier"].to_numpy()).all(), "primary tiers not reproduced"
meta = json.load(open(RES / "g23_pubmed_metadata.json"))
lst = json.load(open(RES / "g23_nonprimary_pmids.json"))
npm = frozenset(lst["nonprimary"])
d1 = depth_per_site(L, cid, drop=npm)
t1 = assign_tier(das, d1).astype(object)
any_src = m[["source_PSP", "source_PhosphoELM", "source_dbPTM"]].fillna(0).astype(bool).any(axis=1).to_numpy()
dropped = (d1 == 0) & ~any_src
t1[dropped] = None
pd.DataFrame({"tier_g23": t1, "depth_g23": d1}).to_parquet(PRIV / "g23_tiers.parquet")

# descriptive: records containing at least one non-primary PMID
has = np.array([any(p in npm for p in l) for l in L])
old = m["tier"].to_numpy()
out = {"reproduces_primary": True, "n_nonprimary_pmids": len(npm),
       "nonprimary": [{"pmid": p, "title": meta["records"][p]["title"], "records_linked": int(freq[p]), "criterion_a_review_type": p in lst["criterion_a_review_type"], "criterion_b_title": p in lst["criterion_b_title"]} for p in sorted(npm)],
       "n_records": int(len(m)), "records_with_at_least_one": int(has.sum()), "share_records_with_at_least_one": float(has.mean()),
       "by_species": m.loc[has, "species"].value_counts().to_dict(), "by_old_tier": pd.Series(old[has]).value_counts().to_dict(),
       "share_of_each_old_tier": {t: float(has[old == t].mean()) for t in ("GOLD", "SILVER", "BRONZE")},
       "records_whose_only_pmids_are_nonprimary": int(((d1 == 0) & (d0 > 0)).sum()), "records_dropped_no_pmid_no_source_flag": int(dropped.sum())}
lab = lambda a: np.where(pd.isna(a), "dropped", a)  # noqa: E731
trans = pd.crosstab(pd.Series(old, name="old tier"), pd.Series(lab(t1), name="new tier"))
out["transition_table"] = {r: {c: int(trans.loc[r, c]) for c in trans.columns} for r in trans.index}
out["tier_counts_old"] = pd.Series(old).value_counts().to_dict()
out["tier_counts_new"] = pd.Series(lab(t1)).value_counts().to_dict()
out["records_whose_tier_changes_or_dropped"] = int((old != lab(t1)).sum())

# step 5: held-out designs under the corrected tiers (Rev 5 split assignments unchanged)
cs = pd.read_parquet(REL / "cluster_split.parquet", columns=K + ["split_cluster_disjoint", "split_protein_disjoint"])
assert (cs[K].values == m[K].values).all()
tn = pd.Series(lab(t1), index=m.index)
to = pd.Series(old, index=m.index)


def counts(mask, tcol):
    return pd.Series(np.asarray(tcol)[mask]).value_counts().to_dict()


h = {}
for name, mask in (("cluster_disjoint_test", cs["split_cluster_disjoint"] == "test"), ("cluster_and_protein_disjoint_test", cs["split_protein_disjoint"] == "test")):
    mk = mask.to_numpy()
    h[name] = {"n": int(mk.sum()), "old": counts(mk, to), "new": counts(mk, tn)}
ts = pd.read_parquet(REL / "temporal_split.parquet")
mk_ = m[K].astype(str).agg("|".join, axis=1).to_numpy()
tk = ts[K].astype(str).agg("|".join, axis=1).to_numpy()
idx, j = [], 0
for i, k in enumerate(mk_):
    if j < len(tk) and k == tk[j]:
        idx.append(i)
        j += 1
assert j == len(tk) and (ts["evidence_tier"].to_numpy() == old[idx]).all()
ts = ts.assign(tier_new=tn.to_numpy()[idx])
tr, te = ts[ts["temporal_split"] == "train"], ts[ts["temporal_split"] == "test"]
trp, trc = set(tr["uniprot_id"]), set(tr["study_cluster"].dropna())
designs = {"temporal_test_all": te, "temporal_protein_disjoint": te[~te["uniprot_id"].isin(trp)], "temporal_study_disjoint": te[~te["study_cluster"].isin(trc)],
           "temporal_protein_and_study_disjoint": te[(~te["uniprot_id"].isin(trp)) & (~te["study_cluster"].isin(trc))]}
for name, d in designs.items():
    h[name] = {"n": int(len(d)), "old": d["evidence_tier"].value_counts().to_dict(), "new": d["tier_new"].value_counts().to_dict()}
h["gold_total"] = {"old": int((old == "GOLD").sum()), "new": int((tn == "GOLD").sum())}
out["held_out"] = h
json.dump(out, open(RES / "g23_tiers_summary.json", "w"), indent=1, default=int)
print("records with >=1 non-primary PMID:", out["records_with_at_least_one"], f"({100 * out['share_records_with_at_least_one']:.1f}%)")
print("old counts", out["tier_counts_old"]); print("new counts", out["tier_counts_new"])
print(trans)
print({k: (v["n"], v["old"].get("GOLD", 0), v["new"].get("GOLD", 0)) for k, v in h.items() if "n" in v}, h["gold_total"])
