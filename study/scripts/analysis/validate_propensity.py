"""Check propensity.py against the Rev 5 phosnet propensity results (study/earlier_results/matching_results.json), original index,
Rev 5 record sets. Also checks the cell-level Newton fit against the record-level statsmodels fit."""
import json
import numpy as np
import base as B
from common import CUR, RES
from populations import population
import propensity as P

M = json.load(open(CUR / "earlier_results" / "matching_results.json"))
b = B.load()
rows = {}
for s, key in (("S2", "conservation"), ("S3", "disorder"), ("S5", "kinase_motif_lookup")):
    d = population(b, s, tier_col="tier", index_col="SEI_V0")
    pd_, lg, r = P.point(d, "diff", np.random.default_rng(20260924))
    ref = M[key]["propensity"]
    rec = {"engine": r}
    if r["status"] == "fitted":
        eta = pd_.fit_cells(np.ones(pd_.n))
        lc = pd_.logit_from_cells(eta[pd_.cell])
        sty = d["residue_type"].isin(list("STY")).to_numpy()
        rec["max_abs_logit_diff_cell_vs_recordlevel_STY_records"] = float(np.max(np.abs(lc - lg)[sty]))
        rec["max_abs_logit_diff_rare_residue_records_BRONZE_only"] = float(np.max(np.abs(lc - lg)[~sty])) if (~sty).any() else 0.0
        rec["phosnet"] = {k: ref[k] for k in ("point_estimate", "n_matched_pairs", "smd")}
        rec["abs_diff_estimate"] = abs(r["estimate"] - ref["point_estimate"])
        rec["pairs_equal"] = r["n_pairs"] == ref["n_matched_pairs"]
    else:
        rec["phosnet"] = ref
    rows[s] = rec
    print(s, json.dumps(rec, default=float)[:600], flush=True)
json.dump(rows, open(RES / "validate_propensity.json", "w"), indent=1, default=float)
