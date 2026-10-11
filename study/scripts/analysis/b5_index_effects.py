"""Task B5: effect sizes with protein-bootstrap CIs for the index properties (replaces the Rev 5 p-values).

Population as in Rev 5 (scripts/earlier_figures/compute_results.py): proteins mapped to a gene, with the index defined; top tier of a
protein = its highest tier over all records (set-B tiers). Domain-annotation presence = any domain-annotated human site of the
protein (human source-set-A file). Reported for the Rev 5 index (V0) and the Rev 6 primary index (V50); the secondary index is
descriptive. No p-values. Resampling unit: protein (already the unit of observation), 1,000 resamples,
generator np.random.default_rng([20260924, crc32(<statistic name>)]).
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

import base as B
from common import CUR, PRIV, REL, RES, N_BOOT, rng_for

b = B.load()
iv = pd.read_parquet(PRIV / "index_variants.parquet")
rank = {"GOLD": 3, "SILVER": 2, "BRONZE": 1}
top = b.assign(r=b["tier"].map(rank)).sort_values("r").groupby("uniprot_id").tail(1)[["uniprot_id", "tier"]].rename(columns={"tier": "top_tier"})
has = b[b["species"] == "human"].groupby("uniprot_id")["in_domain"].max().rename("has_domain").reset_index()
P = iv[iv["mapped_to_gene"] == True][["uniprot_id", "SEI_V0", "SEI_V50", "SEI_secondary"]].merge(top, on="uniprot_id").merge(has, on="uniprot_id", how="left")  # noqa: E712
P = P.dropna(subset=["SEI_V0"]).reset_index(drop=True)
print("proteins", len(P))


def ci(a):
    a = np.asarray(a, float)
    a = a[~np.isnan(a)]
    return [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]


def kw_eta2(d, col):
    g = [d.loc[d["top_tier"] == t, col].to_numpy() for t in ("GOLD", "SILVER", "BRONZE")]
    H = stats.kruskal(*g).statistic
    n = sum(len(x) for x in g)
    return (H - len(g) + 1) / (n - len(g))


def rank_biserial(d, col):
    g, c = d.loc[d["top_tier"] == "GOLD", col].to_numpy(), d.loc[d["top_tier"] == "BRONZE", col].to_numpy()
    u = stats.mannwhitneyu(g, c).statistic
    return 1 - 2 * u / (len(g) * len(c))


def spearman(d, a, c):
    return stats.spearmanr(d[a], d[c])[0]


def pbis(d, col):
    x = d.dropna(subset=["has_domain", col])
    return np.corrcoef(x["has_domain"].astype(float), x[col])[0, 1]


STATS = {}
for iname, icol in (("V0 (Rev 5 index)", "SEI_V0"), ("V50 (Rev 6 primary index)", "SEI_V50")):
    STATS[f"kruskal_eta2_{icol}_by_top_tier"] = (lambda d, c=icol: kw_eta2(d, c), "Kruskal-Wallis eta-squared of " + iname + " by highest tier")
    STATS[f"rank_biserial_{icol}_GOLD_vs_BRONZE"] = (lambda d, c=icol: rank_biserial(d, c), "rank-biserial r, " + iname + ", GOLD vs BRONZE proteins (highest tier)")
    STATS[f"spearman_{icol}_vs_secondary"] = (lambda d, c=icol: spearman(d, c, "SEI_secondary"), "Spearman rho, " + iname + " vs secondary index")
    STATS[f"pointbiserial_annotation_{icol}"] = (lambda d, c=icol: pbis(d, c), "point-biserial r, domain-annotation presence vs " + iname)
STATS["pointbiserial_annotation_SEI_secondary"] = (lambda d: pbis(d, "SEI_secondary"), "point-biserial r, domain-annotation presence vs secondary index")
STATS["kruskal_eta2_SEI_secondary_by_top_tier"] = (lambda d: kw_eta2(d, "SEI_secondary"), "Kruskal-Wallis eta-squared of secondary index by highest tier")
STATS["rank_biserial_SEI_secondary_GOLD_vs_BRONZE"] = (lambda d: rank_biserial(d, "SEI_secondary"), "rank-biserial r, secondary index, GOLD vs BRONZE proteins")

out = {"n_proteins": int(len(P)), "n_boot": N_BOOT, "resampling_unit": "protein", "p_values": "not reported (dropped in Rev 6; proteins share studies and n is large)",
       "means_by_top_tier": {c: P.groupby("top_tier")[c].mean().to_dict() for c in ("SEI_V0", "SEI_V50", "SEI_secondary")},
       "n_by_top_tier": P["top_tier"].value_counts().to_dict(), "effects": {}}
for name, (fn, desc) in STATS.items():
    rng = rng_for("B5|" + name)
    est = fn(P)
    bs = np.empty(N_BOOT)
    for i in range(N_BOOT):
        bs[i] = fn(P.iloc[rng.integers(0, len(P), len(P))])
    out["effects"][name] = {"description": desc, "estimate": float(est), "ci95": ci(bs), "seed": f"[20260924, crc32(B5|{name})]"}
    print(name, round(est, 4), [round(v, 4) for v in out["effects"][name]["ci95"]], flush=True)
rev5 = json.load(open(CUR / "earlier_results" / "matching_results.json"))["sei_properties"]
out["rev5_check"] = {"n_proteins_rev5": rev5["n_proteins"], "rev5_kruskal_eta2": rev5["kruskal_eta2_H"], "recomputed_kruskal_eta2_V0": out["effects"]["kruskal_eta2_SEI_V0_by_top_tier"]["estimate"],
                     "rev5_spearman_rho": rev5["spearman_rho"], "recomputed_spearman_V0_vs_secondary": out["effects"]["spearman_SEI_V0_vs_secondary"]["estimate"],
                     "rev5_rank_biserial": rev5["gold_vs_bronze_rank_biserial_r"], "recomputed_rank_biserial_V0": out["effects"]["rank_biserial_SEI_V0_GOLD_vs_BRONZE"]["estimate"]}
json.dump(out, open(RES / "b5_index_effect_sizes.json", "w"), indent=1, default=float)
print(out["rev5_check"])
