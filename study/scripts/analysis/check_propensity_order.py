"""Does the inherited propensity matcher depend on record order? Human domain overlap (G1-0003 population, S/T/Y records), cell-level logits,
matched odds ratio under (a) file order with different matching permutations and (b) shuffled row orders. Writes results/propensity_order_dependence.json."""
import json
import numpy as np
import run_cem_grid as R
import propensity as P
from common import RES

row = R.man[R.man.contrast_id == "G1-0003"].iloc[0]
d = R.population_for(row, "tier", "SEI_V50")
d = d[d.residue_type.isin(["S", "T", "Y"])]
out = {"population": "G1-0003 (human domain overlap, set-B tiers, V50), serine/threonine/tyrosine records", "n_records": int(len(d))}
pd0 = P.PropData(d, use_species=False)
eta = pd0.fit_cells(np.ones(pd0.n)); lg = pd0.logit_from_cells(eta[pd0.cell])
out["distinct_logit_values"] = int(len(np.unique(lg)))
out["file_order_different_matching_permutations"] = [float(P.match_estimate(lg, pd0.t, pd0.y, pd0.x, "or", np.random.default_rng(s))["estimate"]) for s in range(8)]
sh = []
for seed in range(8):
    dd = d.sample(frac=1, random_state=seed)
    pdat = P.PropData(dd, use_species=False)
    e = pdat.fit_cells(np.ones(pdat.n)); l2 = pdat.logit_from_cells(e[pdat.cell])
    r = P.match_estimate(l2, pdat.t, pdat.y, pdat.x, "or", np.random.default_rng(100 + seed))
    sh.append({"estimate": float(r["estimate"]), "pairs": int(r["n_pairs"])})
out["shuffled_row_orders"] = sh
fo = out["file_order_different_matching_permutations"]
out["summary"] = {"file_order_min": min(fo), "file_order_max": max(fo), "shuffled_min": min(x["estimate"] for x in sh), "shuffled_max": max(x["estimate"] for x in sh)}
json.dump(out, open(RES / "propensity_order_dependence.json", "w"), indent=1)
print(out["summary"], out["distinct_logit_values"])
