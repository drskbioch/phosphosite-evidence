"""G23 step 4: D0 and D1 of G1 and G2 re-run with the corrected tiers (tier_g23; g23_tiers.py). CEM on the V50 decile, gene-cluster bootstrap, 1,000 resamples,
D0 and D1 from the same resamples, generator crc32(D1 id). Output results/g23_results.json."""
import json

import pandas as pd

import run_cem_grid as R
from common import PRIV, RES

t = pd.read_parquet(PRIV / "g23_tiers.parquet")
assert len(t) == len(R.b)
R.b["tier_g23"] = t["tier_g23"].values
m = R.man[(R.man["block"] == "G23") & R.man["design"].isin(["D0", "D1"])]
keys = ["signal", "population", "index", "tier_rule", "variation"]
out = {}
for _, g in m.groupby(keys, sort=False):
    res = R.run_group(g, scheme="gene", tier_col="tier_g23")
    out.update(res)
    for cid, r in res.items():
        print(cid, r.get("design"), r.get("status"), r.get("estimate"), r.get("ci"), r.get("attenuation"), r.get("n_sites_gold"), r.get("n_sites_bronze"), flush=True)
    json.dump(out, open(RES / "g23_results.json", "w"), indent=1, default=float)
