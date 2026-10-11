"""Block G14: the Rev 5 CEM designs (original index; set-A tier for S1, set-B tier for the others; Rev 5 record sets) with
protein-cluster bootstrap CIs, next to the Rev 5 stratified-bootstrap CIs. Shows what changes from the inference alone.
(Propensity rows of G14 are produced by run_propensity.py.)
"""
import json

import numpy as np

import run_cem_grid as R
from common import CUR, N_BOOT, RES, rng_for
from engine import Cells, summarize
from populations import SIGNALS, population

M = json.load(open(CUR / "earlier_results" / "matching_results.json"))
names = {"S1": "domain_overlap_human", "S2": "conservation", "S3": "disorder", "S4": "kinase_motif_sequence", "S5": "kinase_motif_lookup"}
out = {}
for _, row in R.man[(R.man["block"] == "G14") & (R.man["design"].isin(["R5_CEM_primary_index", "R5_CEM_secondary_index"]))].iterrows():
    s = row["signal"]
    sec = row["design"].endswith("secondary_index")
    tcol = "tier_setA" if s == "S1" else "tier"
    icol = "SEI_secondary" if sec else "SEI_V0"
    mapped = (R.b["mapped_to_gene"] == True).to_numpy()  # noqa: E712 - Rev 5 restricted every design to proteins mapped to a gene
    d = population(R.b, s, tier_col=tcol, index_col=icol, mask=mapped)
    note = ""
    if s == "S4":  # Rev 5 record set: the 43 records with an undefined outcome stay in the index deciles (outcome set to 0; they sit in single-tier strata)
        d = population(R.b, "S5", tier_col=tcol, index_col=icol, mask=mapped)
        d["y"] = R.b.loc[d.index, "seq_motif_match"].astype(float).fillna(0.0)
        note = "Rev 5 record set: 43 records with an undefined outcome kept for the index deciles; crude estimate not reported for this row"
    c = Cells(d, "y", "x", {"protein": ["uniprot_id"]}, strata_cols=SIGNALS[s]["strata"])
    r = summarize(c, SIGNALS[s]["stat"], "protein", rng_for(row["contrast_id"]), n_boot=N_BOOT)
    if not isinstance(r, tuple):
        out[row["contrast_id"]] = {"status": "failed", "reason": r.get("error")}
        continue
    r, _ = r
    ref = M[names[s]]["cem_secondary" if sec else "cem_primary"]
    null = 1.0 if SIGNALS[s]["stat"] == "or" else 0.0
    w5, w6 = ref["ci_high"] - ref["ci_low"], r["adjusted_ci"][1] - r["adjusted_ci"][0]
    out[row["contrast_id"]] = {"signal": s, "design": row["design"], "status": "run", "estimate": r["adjusted"], "ci_protein_cluster": r["adjusted_ci"],
                               "rev5_estimate": ref["point_estimate"], "rev5_ci_stratified": [ref["ci_low"], ref["ci_high"]],
                               "estimate_abs_diff_from_rev5": abs(r["adjusted"] - ref["point_estimate"]),
                               "ci_width_ratio_protein_cluster_over_rev5": w6 / w5,
                               "rev5_ci_excludes_null": bool(ref["ci_low"] > null or ref["ci_high"] < null),
                               "protein_cluster_ci_excludes_null": bool(r["adjusted_ci"][0] > null or r["adjusted_ci"][1] < null),
                               "smd": r["smd"], "n_sites": r["n_sites"], "n_proteins": r["n_proteins"], "n_boot": N_BOOT, "note": note,
                               "seed": f"[20260924, crc32({row['contrast_id']})]"}
    print(row["contrast_id"], s, row["design"], round(r["adjusted"], 4), [round(v, 4) for v in r["adjusted_ci"]], "rev5", [round(ref["ci_low"], 4), round(ref["ci_high"], 4)], "ratio", round(w6 / w5, 2), flush=True)
json.dump(out, open(RES / "g14_cem.json", "w"), indent=1, default=float)
