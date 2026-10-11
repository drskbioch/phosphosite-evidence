"""Check that engine.point() reproduces the Rev 5 phosnet CEM point estimates and crude estimates.

Rev 5 populations: original index (SEI_primary), set-B tier for S2 to S5 and set-A tier for S1.
Compares with study/earlier_results/matching_results.json (cem_primary, unmatched_mapped_sites).
"""
import json

import numpy as np

import base as B
from common import CUR, RES
from engine import Cells, point
from populations import SIGNALS, population

M = json.load(open(CUR / "earlier_results" / "matching_results.json"))
names = {"S1": "domain_overlap_human", "S2": "conservation", "S3": "disorder",
         "S4": "kinase_motif_sequence", "S5": "kinase_motif_lookup"}
b = B.load()
rows = {}
for s, key in names.items():
    tier = "tier_setA" if s == "S1" else "tier"
    d = population(b, s, tier_col=tier, index_col="SEI_V0")
    if s == "S4":  # see below: Rev 5 kept sites with an undefined outcome when forming the index deciles
        d5 = population(b, "S5", tier_col=tier, index_col="SEI_V0")
        d5["y"] = b.loc[d5.index, "seq_motif_match"].astype(float).fillna(0.0)  # undefined outcome: set to 0, these sites sit in single-tier strata and never enter an estimate
        c5 = Cells(d5, "y", "x", {"protein": ["uniprot_id"]}, strata_cols=SIGNALS[s]["strata"])
        p5 = point(c5, SIGNALS[s]["stat"])
        rows["S4_rev5_record_set"] = {"engine_adjusted": p5["adj"], "phosnet_cem_primary": M[key]["cem_primary"]["point_estimate"],
                                      "engine_smd": p5["smd"], "phosnet_smd": M[key]["cem_primary"]["smd"], "engine_n": c5.n_sites,
                                      "phosnet_n_analysis": M[key]["cem_primary"]["n_analysis"],
                                      "adj_abs_diff": abs(p5["adj"] - M[key]["cem_primary"]["point_estimate"]),
                                      "note": "46 sites have an undefined sequence-motif outcome (43 in GOLD/BRONZE with an index). Rev 5 kept them in the decile computation; Rev 6 drops them before it."}
        print("S4 (Rev 5 record set)", rows["S4_rev5_record_set"], flush=True)
    c = Cells(d, "y", "x", {"protein": ["uniprot_id"]}, strata_cols=SIGNALS[s]["strata"])
    p = point(c, SIGNALS[s]["stat"])
    ref = M[key]
    rows[s] = {"engine_adjusted": p["adj"], "phosnet_cem_primary": ref["cem_primary"]["point_estimate"],
               "engine_crude": p["crude"], "phosnet_unmatched_mapped": ref["unmatched_mapped_sites"]["point_estimate"],
               "engine_smd": p["smd"], "phosnet_smd": ref["cem_primary"]["smd"],
               "engine_n": c.n_sites, "phosnet_n_analysis": ref["cem_primary"]["n_analysis"],
               "adj_abs_diff": abs(p["adj"] - ref["cem_primary"]["point_estimate"]),
               "crude_abs_diff": abs(p["crude"] - ref["unmatched_mapped_sites"]["point_estimate"])}
    print(s, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in rows[s].items()}, flush=True)
json.dump(rows, open(RES / "validate_engine.json", "w"), indent=2)
