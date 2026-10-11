"""Tasks D3 (Pharos/TCRD TDL), D4 (GO annotation count) and D5 (choice of variant) for the index variants.

Decisions fixed before any result:
 - Proteins: human proteins of the site table with a defined index (per variant) and, for D3, a Pharos entry (matched on UniProt accession).
 - Ranking of variants (D5) uses the COMMON set (proteins where V0, V20, V50, V100 and VT2025 are all defined), higher = better supported:
   AUROC for separating non-Tdark (positive) from Tdark with the index as the score; Spearman correlation of the index with the number of distinct GO terms.
   VT2026 is a reference-year sensitivity and is reported but not ranked.
 - GO: GOA human GAF 2.2; distinct GO IDs per protein over rows whose qualifier does not contain NOT; main = all evidence codes, sensitivity = excluding IEA.
   Human proteins without any GOA row count 0 terms (the number is reported).
 - D5 rule from the plan: if both measures rank the same variant first and it is not V50, G1 is reported with V50 and that variant; if they disagree all
   variants are reported and none is chosen.
Protein-bootstrap CIs (1,000 resamples; generator derived from the label).
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

from common import N_BOOT, PRIV, RES, rng_for

iv = pd.read_parquet(PRIV / "index_variants.parquet")
H = iv[(iv["species"] == "human") & (iv["mapped_to_gene"] == True)].copy()  # noqa: E712
VARS = ["V0", "V20", "V50", "V100", "VT2025", "VT2026"]
ICOL = {"V0": "SEI_V0", "V20": "SEI_V20", "V50": "SEI_V50", "V100": "SEI_V100", "VT2025": "SEI_VT2025", "VT2026": "SEI_VT2026"}
RANKED = ["V0", "V20", "V50", "V100", "VT2025"]

# ---- Pharos
ph = pd.read_parquet(PRIV / "pharos_tdl.parquet")[["uniprot", "tdl"]].drop_duplicates("uniprot").rename(columns={"uniprot": "uniprot_id"})
H = H.merge(ph, on="uniprot_id", how="left")
# ---- GO
gaf = pd.read_csv(PRIV / "dl" / "goa_human.gaf.gz", sep="\t", comment="!", header=None, usecols=[1, 3, 4, 6], names=["uniprot_id", "qual", "go", "ev"], dtype=str)
gaf = gaf[~gaf["qual"].fillna("").str.contains("NOT")]
go_all = gaf.groupby("uniprot_id")["go"].nunique().rename("go_all")
go_noiea = gaf[gaf["ev"] != "IEA"].groupby("uniprot_id")["go"].nunique().rename("go_noIEA")
H = H.merge(go_all, on="uniprot_id", how="left").merge(go_noiea, on="uniprot_id", how="left")
n_no_go = int(H["go_all"].isna().sum())
H["go_all"] = H["go_all"].fillna(0)
H["go_noIEA"] = H["go_noIEA"].fillna(0)

info = {"human_proteins_with_V0_index": int(len(H)), "human_proteins_with_pharos_entry": int(H["tdl"].notna().sum()),
        "pharos_entries_total": int(len(ph)), "pharos_entries_matched_to_index_proteins": int(ph["uniprot_id"].isin(H["uniprot_id"]).sum()),
        "tdl_counts_in_index_proteins": H["tdl"].value_counts().to_dict(), "human_proteins_without_goa_rows": n_no_go,
        "goa_rows_after_NOT_filter": int(len(gaf)), "goa_proteins": int(gaf["uniprot_id"].nunique())}
common = H.dropna(subset=[ICOL[v] for v in RANKED]).copy()
info["common_set_proteins"] = int(len(common))
info["common_set_with_pharos"] = int(common["tdl"].notna().sum())
info["proteins_with_index_by_variant"] = {v: int(H[ICOL[v]].notna().sum()) for v in VARS}


def auroc(score, pos):
    r = stats.rankdata(score)
    n1, n0 = pos.sum(), (~pos).sum()
    return (r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def boot(fn, d, label):
    rng = rng_for(label)
    est = fn(d)
    bs = np.array([fn(d.iloc[rng.integers(0, len(d), len(d))]) for _ in range(N_BOOT)])
    return float(est), [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]


res = {"info": info, "tdl_means": {}, "auroc_tdark": {}, "go_spearman": {}}
for setname, S in (("own_set", None), ("common_set", common)):
    for v in VARS:
        d = (H if S is None else S).dropna(subset=[ICOL[v]])
        c = ICOL[v]
        res["tdl_means"][f"{setname}|{v}"] = {t: {"n": int((d["tdl"] == t).sum()), "mean": float(d.loc[d["tdl"] == t, c].mean()), "median": float(d.loc[d["tdl"] == t, c].median())} for t in ("Tclin", "Tchem", "Tbio", "Tdark")}
        dd = d.dropna(subset=["tdl"]).reset_index(drop=True)
        e, ci = boot(lambda x, c=c: auroc(x[c].to_numpy(), (x["tdl"] != "Tdark").to_numpy()), dd, f"D3|{setname}|{v}")
        res["auroc_tdark"][f"{setname}|{v}"] = {"auroc_nonTdark_vs_Tdark": e, "ci95": ci, "n_proteins": int(len(dd)), "n_tdark": int((dd["tdl"] == "Tdark").sum())}
        dd2 = d.reset_index(drop=True)
        for gname, gcol in (("all_evidence", "go_all"), ("excluding_IEA", "go_noIEA")):
            e, ci = boot(lambda x, c=c, gcol=gcol: stats.spearmanr(x[c], x[gcol])[0], dd2, f"D4|{setname}|{v}|{gname}")
            res["go_spearman"][f"{setname}|{v}|{gname}"] = {"spearman": e, "ci95": ci, "n_proteins": int(len(dd2))}
        print(setname, v, round(res["auroc_tdark"][f"{setname}|{v}"]["auroc_nonTdark_vs_Tdark"], 4), round(res["go_spearman"][f"{setname}|{v}|all_evidence"]["spearman"], 4), flush=True)

# ---- post-hoc (labelled): GO correlation restricted to proteins with at least one GOA row. The pre-specified rule counted proteins absent from goa_human.gaf as
# having 0 terms, but most of them are unreviewed accessions that the file does not cover, so the count measures coverage of the file rather than annotation depth.
ann = common[common["go_all"] > 0]
res["go_spearman_posthoc_GOA_annotated_only"] = {}
for v in VARS:
    d = ann.dropna(subset=[ICOL[v]]).reset_index(drop=True) if v != "VT2026" else ann.dropna(subset=[ICOL[v]]).reset_index(drop=True)
    c = ICOL[v]
    for gname, gcol in (("all_evidence", "go_all"), ("excluding_IEA", "go_noIEA")):
        dd = d[d[gcol] > 0].reset_index(drop=True) if gname == "excluding_IEA" else d
        e, ci = boot(lambda x, c=c, gcol=gcol: stats.spearmanr(x[c], x[gcol])[0], dd, f"D4posthoc|{v}|{gname}")
        res["go_spearman_posthoc_GOA_annotated_only"][f"{v}|{gname}"] = {"spearman": e, "ci95": ci, "n_proteins": int(len(dd))}
    print("posthoc", v, round(res["go_spearman_posthoc_GOA_annotated_only"][f"{v}|all_evidence"]["spearman"], 4), flush=True)
gp = {v: res["go_spearman_posthoc_GOA_annotated_only"][f"{v}|all_evidence"]["spearman"] for v in RANKED}
res["D5_posthoc_GOA_annotated_only"] = {"ranking_go_spearman": dict(sorted(gp.items(), key=lambda kv: -kv[1])), "first_by_go": max(gp, key=gp.get),
                                        "first_by_auroc": max({v: res["auroc_tdark"][f"common_set|{v}"]["auroc_nonTdark_vs_Tdark"] for v in RANKED}, key=lambda v: res["auroc_tdark"][f"common_set|{v}"]["auroc_nonTdark_vs_Tdark"]),
                                        "label": "post-hoc, reported beside the pre-specified rule; does not replace it"}
# ---- D5
a = {v: res["auroc_tdark"][f"common_set|{v}"]["auroc_nonTdark_vs_Tdark"] for v in RANKED}
g = {v: res["go_spearman"][f"common_set|{v}|all_evidence"]["spearman"] for v in RANKED}
ga = max(a, key=a.get)
gg = max(g, key=g.get)
res["D5"] = {"ranking_auroc_common_set": dict(sorted(a.items(), key=lambda kv: -kv[1])), "ranking_go_spearman_common_set": dict(sorted(g.items(), key=lambda kv: -kv[1])),
             "first_by_auroc": ga, "first_by_go": gg, "same_first": ga == gg}
if ga == gg and ga != "V50":
    res["D5"]["decision"] = f"Both measures rank {ga} first and it is not V50: report G1 with V50 and with {ga}, and state both."
elif ga == gg:
    res["D5"]["decision"] = "Both measures rank V50 first: V50 remains the only primary index."
else:
    res["D5"]["decision"] = "The measures disagree on the first variant: report all variants and choose none (V50 stays the primary index by decision D2)."
print(json.dumps(res["D5"], indent=1))
json.dump(res, open(RES / "d3_d4_d5_index_checks.json", "w"), indent=1, default=float)
