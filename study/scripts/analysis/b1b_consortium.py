"""G7 addition: remove the seven consortium PMIDs (named in audit_reports/phase0_epsd_gate_report.md) and record
counts under both readings of "no source evidence" (DAS == 0, and no PSP/Phospho.ELM/dbPTM flag)."""
import json
import pandas as pd
from common import PRIV, REL, RES, depth_per_site, split_pmids, cluster_map

PM = ["37582339", "37582357", "37582358", "37582362", "37582371", "37619559", "38359819"]
K = ["uniprot_id", "residue_type", "position"]
m = pd.read_parquet(REL / "mammalian_sites.parquet", columns=K + ["species", "DAS", "pmid_list_str", "source_PSP", "source_PhosphoELM", "source_dbPTM"])
L = split_pmids(m["pmid_list_str"])
cid, _, freq = cluster_map(L, 884)
d = depth_per_site(L, cid, drop=frozenset(PM))
dv = pd.read_parquet(PRIV / "depth_variants.parquet")
dv["depth_drop_consortium7"] = d
dv.to_parquet(PRIV / "depth_variants.parquet")
info = json.load(open(RES / "b_cluster_rule_info.json"))
top = json.load(open(RES / "b_largest_clusters.json"))
any_src = (m[["source_PSP", "source_PhosphoELM", "source_dbPTM"]].fillna(0).astype(bool)).any(axis=1)
nopmid = d == 0
info["drop_consortium7"] = {
    "pmids": PM, "n_pmids_present_in_data": int(sum(p in freq for p in PM)),
    "sites_carrying_any_of_them": int(sum(any(p in PM for p in l) for l in L)),
    "sites_with_no_pmid_left": int(nopmid.sum()),
    "sites_dropped_no_pmid_and_no_database_source_flag": int((nopmid & ~any_src).sum()),
    "sites_dropped_no_pmid_and_DAS_0": int((nopmid & (m["DAS"] == 0)).sum()),
    "identical_pmid_set_to_drop_2": sorted(top[1]["pmids"]) == sorted(PM),
    "drop_2_label_and_pmids": {"label": top[1]["label"], "pmids": sorted(top[1]["pmids"])}}
for k in range(1, 11):  # add the source-flag reading to every removal
    dd = dv[f"depth_drop_{k}"].values == 0
    info[f"drop_{k}"]["sites_dropped_no_pmid_and_no_database_source_flag"] = int((dd & ~any_src.values).sum())
json.dump(info, open(RES / "b_cluster_rule_info.json", "w"), indent=1, default=str)
print(json.dumps(info["drop_consortium7"], indent=1)); print({k: info[f"drop_{k}"]["sites_dropped_no_pmid_and_no_database_source_flag"] for k in range(1, 11)})
print(m[["source_PSP", "source_PhosphoELM", "source_dbPTM"]].dtypes.tolist(), (~any_src).sum())
