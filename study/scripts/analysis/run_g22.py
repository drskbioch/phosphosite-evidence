"""G22 (post-hoc): sequence-based motif contrast restricted to serine and threonine records. Definition in the plan, deviation log of 2026-10-04, fourth batch.
D0 and D1 from the same gene-cluster bootstrap (1,000 resamples), generator rng_for(D1 id). Output results/g22_results.json."""
import json

import run_cem_grid as R
from common import RES

m = R.man[(R.man["block"] == "G22") & R.man["design"].isin(["D0", "D1"])]
mask = R.b["residue_type"].isin(["S", "T"]).to_numpy()
out = R.run_group(m, scheme="gene", mask=mask)
json.dump(out, open(RES / "g22_results.json", "w"), indent=1, default=float)
for cid, r in out.items():
    print(cid, r.get("design"), r.get("status"), r.get("estimate"), r.get("ci"), r.get("attenuation"), r.get("n_sites_gold"), r.get("n_sites_bronze"), flush=True)
