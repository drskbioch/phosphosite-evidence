"""Tier columns for Task F (G7 leave-one-cluster-out, G8 tier thresholds). Row-aligned with mammalian_sites.parquet (like depth_variants.parquet).

G7: for each removal, the removed PMIDs are deleted from every site's publication list (depth_drop_k in depth_variants.parquet, cluster assignment kept),
the tier is recomputed with the study-cluster rule (GOLD: DAS >= 2 or depth >= 10; SILVER: DAS >= 1 or depth 3-9; BRONZE otherwise), DAS unchanged, and sites left
with no PMID and no database source flag (PSP, Phospho.ELM, dbPTM) are dropped (tier missing).
G8: GOLD if depth >= gd or DAS >= ga; SILVER if not GOLD and (DAS >= 1 or depth >= 3); BRONZE otherwise; depth = study-cluster depth at 884.
"""
import json

import numpy as np
import pandas as pd

from common import PRIV, REL, RES

m = pd.read_parquet(REL / "mammalian_sites.parquet", columns=["DAS", "source_PSP", "source_PhosphoELM", "source_dbPTM"])
dv = pd.read_parquet(PRIV / "depth_variants.parquet")
das = m["DAS"].to_numpy()
any_src = m[["source_PSP", "source_PhosphoELM", "source_dbPTM"]].fillna(0).astype(bool).any(axis=1).to_numpy()
out = pd.DataFrame(index=m.index)
rep = {}


def tier(depth, gd, ga):
    t = np.full(len(depth), "BRONZE", dtype=object)
    t[(das >= 1) | (depth >= 3)] = "SILVER"
    t[(depth >= gd) | (das >= ga)] = "GOLD"
    return t


for k in list(range(1, 11)) + ["consortium7"]:
    d = dv[f"depth_drop_{k}"].to_numpy()
    t = tier(d, 10, 2)
    dropped = (d == 0) & ~any_src
    t[dropped] = None
    out[f"tier_drop_{k}"] = t
    rep[f"drop_{k}"] = {"sites_dropped": int(dropped.sum()), "GOLD": int((t == "GOLD").sum()), "SILVER": int((t == "SILVER").sum()), "BRONZE": int((t == "BRONZE").sum())}
d884 = dv["depth_884"].to_numpy()
for gd in (5, 10, 20):
    for ga in (2, 3):
        t = tier(d884, gd, ga)
        out[f"tier_g8_d{gd}_a{ga}"] = t
        rep[f"g8_d{gd}_a{ga}"] = {"GOLD": int((t == "GOLD").sum()), "SILVER": int((t == "SILVER").sum()), "BRONZE": int((t == "BRONZE").sum())}
# the primary rule must be reproduced by tier_g8_d10_a2
prim = pd.read_parquet(REL / "mammalian_sites.parquet", columns=["tier"])["tier"].to_numpy()
rep["check_g8_d10_a2_equals_primary_tier"] = bool((out["tier_g8_d10_a2"].to_numpy() == prim).all())
out.to_parquet(PRIV / "f_tiers.parquet")
json.dump(rep, open(RES / "f_tiers_summary.json", "w"), indent=1)
print(json.dumps(rep, indent=1)[:3000])
