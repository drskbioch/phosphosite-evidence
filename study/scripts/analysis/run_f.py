"""Task F (G6 residue strata, G7 leave-one-cluster-out, G8 tier thresholds): D0 and D1 with the gene-cluster bootstrap (primary from the plan amendment of
2026-10-01) and, beside it, the accession-cluster bootstrap of the interim stage. 1,000 resamples each, whole pipeline re-run in every resample.
Seeds: accession cluster rng_for(D1 id), gene cluster rng_for(D1 id + "|gene").
Usage: python run_f.py G7 S1,S2    -> results/f_runs/G7_S1-S2.json  ({contrast id: {"protein": ..., "gene": ...}})
"""
import json
import re
import sys

import pandas as pd

import run_cem_grid as R
from common import PRIV, RES

OUTD = RES / "f_runs"
OUTD.mkdir(exist_ok=True)
b = R.b
ft = pd.read_parquet(PRIV / "f_tiers.parquet")
assert len(ft) == len(b)
for c in ft.columns:
    b[c] = ft[c].values


def tier_col(g):
    block = g["block"].iloc[0]
    v = g["variation"].iloc[0]
    if block == "G6":
        return "tier"
    if block == "G7":
        k = "consortium7" if v == "drop_consortium_7_PMIDs" else v.replace("drop_cluster_", "")
        return f"tier_drop_{k}"
    m = re.match(r"gold_depth=(\d+);gold_das=(\d+)", v)
    return f"tier_g8_d{m.group(1)}_a{m.group(2)}"


if __name__ == "__main__":
    block = sys.argv[1]
    sigs = sys.argv[2].split(",") if len(sys.argv) > 2 else ["S1", "S2", "S3", "S4"]
    m = R.man[(R.man["block"] == block) & R.man["signal"].isin(sigs) & R.man["design"].isin(["D0", "D1"])]
    keys = ["block", "signal", "population", "index", "tier_rule", "variation"]
    out = {}
    fn = OUTD / f"{block}_{'-'.join(sigs)}.json"
    for _, g in m.groupby(keys, sort=False):
        tc = tier_col(g)
        d1id = g[g["design"] == "D1"]["contrast_id"].iloc[0]
        for scheme, seed in (("protein", None), ("gene", d1id + "|gene")):
            res = R.run_group(g, scheme=scheme, tier_col=tc, seed_id=seed)
            for cid, r in res.items():
                out.setdefault(cid, {})[scheme] = r
                r["tier_column"] = tc
            print(d1id, scheme, [(r.get("design"), round(r["estimate"], 4) if "estimate" in r else r.get("reason")) for r in res.values()], flush=True)
        json.dump(out, open(fn, "w"), indent=1, default=float)
