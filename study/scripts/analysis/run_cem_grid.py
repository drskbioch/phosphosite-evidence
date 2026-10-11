"""Run the CEM (D0 crude, D1 CEM) contrasts of grid_manifest.csv with the protein-cluster bootstrap.

D0 and D1 of one contrast come from the same bootstrap: in every resample the whole pipeline (index deciles, strata,
weights, estimate) is re-run on the resampled proteins and the crude estimate is computed on the same records,
so attenuation is computed within the resample. The generator is np.random.default_rng([20260924, crc32(D1 contrast id)]).

Usage: python run_cem_grid.py G1 G2 ...      (blocks; results are merged into results/cem_results.json)
"""
import json
import sys
import time

import numpy as np
import pandas as pd

import base as B
from common import PRIV, REV6, RES, N_BOOT, rng_for
from engine import Cells, summarize
from populations import SIGNALS, population

OUT = RES / "cem_results.json"
man = pd.read_csv(REV6 / "grid_manifest.csv").fillna("")
b = B.load("lean")
dv = pd.read_parquet(PRIV / "depth_variants.parquet")
assert len(dv) == len(b)
b["tier_depth"] = dv["tier_depth_only_rule_884"].values
b["tier_r500"] = dv["tier_study_cluster_rule_500"].values
b["tier_r1000"] = dv["tier_study_cluster_rule_1000"].values

has = b[b["species"] == "human"].groupby("uniprot_id")["in_domain"].max()
ANNOT = set(has[has == 1.0].index)

TIER = {"study-cluster rule (884)": "tier", "study-cluster rule (500)": "tier_r500", "study-cluster rule (1000)": "tier_r1000",
        "set-A tier": "tier_setA", "set-A tier (S1); set-B tier (others)": None}
IDX = {"V50": "SEI_V50", "V0": "SEI_V0", "V0 (original)": "SEI_V0", "V20": "SEI_V20", "V100": "SEI_V100",
       "VT2025": "SEI_VT2025", "VT2026": "SEI_VT2026"}


def tier_col_for(rule, signal):
    if rule.startswith("depth-only"):
        return "tier_depth"
    if rule.startswith("set-A tier (S1)"):
        return "tier_setA" if signal == "S1" else "tier"
    return TIER[rule]


def population_for(row, tier_col, index_col, mask=None):
    s, popl = row["signal"], row["population"]
    kw = {}
    if popl == "P-human":
        kw["species"] = "human"
    if popl.startswith("P-") and popl[2:] in ("S", "T", "Y"):
        kw["residues"] = [popl[2:]]
    d = population(b, s, tier_col=tier_col, index_col=index_col, mask=mask, **kw)
    if popl == "P-annot":
        d = d[d["uniprot_id"].isin(ANNOT)]
    return d


def label(att, crude_ci, null, stat):
    if not ((crude_ci[0] > null) or (crude_ci[1] < null)):
        return "not interpreted (crude interval includes the null)"
    if att is None or np.isnan(att):
        return "undefined"
    if att < 0:
        return "not attenuated"
    return "largely explained (>= 50%)" if att >= 0.5 else "partly attenuated (< 50%)"


def run_group(g, scheme="protein", mask=None, tier_col=None, seed_id=None):
    """scheme "protein": accession clusters (Rev 6 interim); "gene": gene clusters (plan amendment of 2026-10-01).
    mask: optional boolean array over the rows of `b` restricting the population. tier_col overrides the tier column."""
    d0 = g[g["design"] == "D0"].iloc[0] if (g["design"] == "D0").any() else None
    d1 = g[g["design"].isin(["D1", "D3"])].iloc[0]
    s = d1["signal"]
    tcol = tier_col or tier_col_for(d1["tier_rule"], s)
    icol = IDX[d1["index"]]
    d = population_for(d1, tcol, icol, mask=mask)
    ucols = {"protein": ["uniprot_id"]} if scheme == "protein" else {"gene": ["gene_cluster"]}
    c = Cells(d, "y", "x", ucols, strata_cols=SIGNALS[s]["strata"])
    t0 = time.time()
    res = summarize(c, SIGNALS[s]["stat"], scheme, rng_for(seed_id or d1["contrast_id"]), n_boot=N_BOOT)
    if not isinstance(res, tuple):
        return {r["contrast_id"]: {"status": "failed", "reason": res.get("error")} for _, r in g.iterrows()}
    res, bs = res
    null = 1.0 if SIGNALS[s]["stat"] == "or" else 0.0
    common = {"signal": s, "population": d1["population"], "index": d1["index"], "tier_rule": d1["tier_rule"],
              "variation": d1["variation"], "decision_bearing": d1["decision_bearing"], "null": null,
              "n_sites": res["n_sites"], "n_sites_gold": res["n_sites_gold"], "n_sites_bronze": res["n_sites_bronze"],
              "n_proteins": res["n_proteins"], "n_boot": res["n_boot"], "n_boot_valid": res["n_boot_valid"],
              "resampling_unit": "protein (uniprot_id)" if scheme == "protein" else "gene cluster (gene_cluster)", "n_gene_clusters": int(d["gene_cluster"].nunique()), "seed": f"[20260924, crc32({seed_id or d1['contrast_id']})]",
              "seconds": round(time.time() - t0, 1)}
    out = {}
    if d0 is not None:
        out[d0["contrast_id"]] = dict(common, design="D0", status="run", estimate=res["crude"], ci=res["crude_ci"],
                                      ci_excludes_null=bool((res["crude_ci"][0] > null) or (res["crude_ci"][1] < null)),
                                      bootstrap_shared_with=d1["contrast_id"])
    att = res["attenuation"]
    out[d1["contrast_id"]] = dict(common, design=d1["design"], status="run", estimate=res["adjusted"], ci=res["adjusted_ci"],
                                  survives=res["survives"], crude=res["crude"], crude_ci=res["crude_ci"],
                                  attenuation=att, attenuation_ci=res["attenuation_ci"],
                                  attenuation_label=label(att, res["crude_ci"], null, SIGNALS[s]["stat"]),
                                  largely_explained=res["largely_explained"], smd=res["smd"], smd_pre=res["smd_pre"], balance_pass=res["balance_pass"],
                                  n_strata_used=res["n_strata_used"], n_strata_dropped=res["n_strata_dropped"], n_bins=res["n_bins"],
                                  adjusted_boot_sd=res["adjusted_boot_sd"])
    return out


if __name__ == "__main__":
    blocks = sys.argv[1:]
    store = json.load(open(OUT)) if OUT.exists() else {}
    m = man[man["block"].isin(blocks) & man["design"].isin(["D0", "D1", "D3"])]
    keys = ["block", "signal", "population", "index", "tier_rule", "variation"]
    for _, g in m.groupby(keys, sort=False):
        out = run_group(g)
        store.update(out)
        for cid, r in out.items():
            print(cid, r.get("design"), r.get("status"), r.get("estimate"), r.get("ci"), r.get("attenuation"), flush=True)
        json.dump(store, open(OUT, "w"), indent=1, default=float)
