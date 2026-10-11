"""Check of cond_logit.py against statsmodels ConditionalLogit on a REAL subsample of the G21 records (P-all, conservation, S/T/Y; random 3,000 gene clusters,
seed 20260924), as required by the G21 definition before the block is run. Tolerances: coefficients 5e-04, standard errors 1e-04 (those of check_cond_logit.py).
Output: results/g21_cond_logit_check.json"""
import json
import time
import warnings

import numpy as np
from statsmodels.discrete.conditional_models import ConditionalLogit

import cond_logit as C
import g21_data as G
from common import RES

warnings.filterwarnings("ignore")
d, info = G.prepare("S2", "P-all (S/T/Y)")
rng = np.random.default_rng(20260924)
cl = d["gene_cluster"].unique()
keep = set(rng.choice(cl, size=3000, replace=False))
s = d[d["gene_cluster"].isin(keep)].reset_index(drop=True)
g = s["gene_cluster"].to_numpy()
_, first = np.unique(g, return_index=True)
starts = np.append(first, len(s)).astype(np.int64)
X = s[["gold", "isT", "isY"]].to_numpy(float)
y = s["y"].to_numpy().astype(np.int64)
sizes = np.diff(starts)
t = time.time()
b1, ok1, ll1, cov1 = C.fit(X, y, starts)
t1 = time.time() - t
t = time.time()
b2, ok2, ll2, cov2 = C.fit(X, y, starts, parallel=True)
t2 = time.time() - t
t = time.time()
res = ConditionalLogit(y, X, groups=g).fit(disp=0, method="bfgs", maxiter=500)
t3 = time.time() - t
se1, se2 = np.sqrt(np.diag(cov1)), np.sqrt(np.diag(cov2))
llsm, gsm, _ = C.loglik(X, y, starts, np.ones(len(starts) - 1), np.asarray(res.params))
out = {"records": int(len(s)), "gene_clusters": int(len(starts) - 1), "largest_stratum": int(sizes.max()), "columns": ["gold", "isT", "isY"],
       "own_serial": {"coef": b1.tolist(), "se": se1.tolist(), "converged": bool(ok1), "loglik": float(ll1), "seconds": round(t1, 2)},
       "own_parallel": {"coef": b2.tolist(), "se": se2.tolist(), "converged": bool(ok2), "loglik": float(ll2), "seconds": round(t2, 2)},
       "statsmodels": {"coef": np.asarray(res.params).tolist(), "se": np.asarray(res.bse).tolist(), "seconds": round(t3, 2)},
       "loglik_at_statsmodels_estimate_own_likelihood": float(llsm),
       "max_abs_diff_coef_serial": float(np.abs(b1 - res.params).max()), "max_abs_diff_se_serial": float(np.abs(se1 - res.bse).max()),
       "max_abs_diff_coef_parallel": float(np.abs(b2 - res.params).max()), "max_abs_diff_se_parallel": float(np.abs(se2 - res.bse).max()),
       "max_abs_diff_coef_serial_vs_parallel": float(np.abs(b1 - b2).max()),
       "tolerance_coef": 5e-4, "tolerance_se": 1e-4}
out["passes"] = bool(out["max_abs_diff_coef_serial"] < 5e-4 and out["max_abs_diff_coef_parallel"] < 5e-4 and out["max_abs_diff_se_serial"] < 1e-4 and out["max_abs_diff_se_parallel"] < 1e-4)
json.dump(out, open(RES / "g21_cond_logit_check.json", "w"), indent=1)
print(json.dumps(out, indent=1))
