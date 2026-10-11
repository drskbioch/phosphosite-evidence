"""Check cond_logit.py against statsmodels ConditionalLogit on simulated data (code check only; nothing here enters a result)."""
import numpy as np
from statsmodels.discrete.conditional_models import ConditionalLogit

import cond_logit as C

rng = np.random.default_rng(7)
G = 300
rows = []
for g in range(G):
    n = int(rng.integers(2, 14))
    x1 = rng.integers(0, 2, n)
    x2 = rng.integers(0, 2, n)
    a = rng.normal() * 0.7
    p = 1 / (1 + np.exp(-(a + 0.6 * x1 - 0.4 * x2)))
    y = (rng.random(n) < p).astype(int)
    for i in range(n):
        rows.append((g, x1[i], x2[i], y[i]))
d = np.array(rows)
order = np.argsort(d[:, 0], kind="stable")
d = d[order]
X = d[:, 1:3].astype(float)
y = d[:, 3].astype(np.int64)
_, first = np.unique(d[:, 0], return_index=True)
starts = np.append(first, len(d)).astype(np.int64)
b, ok, ll, cov = C.fit(X, y, starts)
res = ConditionalLogit(y, X, groups=d[:, 0]).fit(disp=0)
print("own", b, np.sqrt(np.diag(cov)), ok, ll)
print("sm ", res.params, res.bse)
ll_sm, g_sm, _ = C.loglik(X, y, starts, np.ones(len(starts) - 1), np.asarray(res.params))
ll_o, g_o, _ = C.loglik(X, y, starts, np.ones(len(starts) - 1), b)
print("log-likelihood own", ll_o, "at statsmodels estimate", ll_sm, "gradient norm own", np.abs(g_o).max(), "at statsmodels estimate", np.abs(g_sm).max())
assert ll_o >= ll_sm - 1e-9 and np.abs(g_o).max() < 1e-6
assert np.allclose(b, res.params, atol=5e-4) and np.allclose(np.sqrt(np.diag(cov)), res.bse, atol=1e-4)
# large stratum (overflow check): 3,000 records in one stratum plus small strata
n = 3000
Xb = np.stack([rng.integers(0, 2, n), rng.integers(0, 2, n)], 1).astype(float)
yb = (rng.random(n) < 0.5).astype(np.int64)
X2 = np.vstack([X, Xb]); y2 = np.concatenate([y, yb]); st2 = np.append(starts[:-1], [len(y), len(y) + n]).astype(np.int64)
b2, ok2, ll2, cov2 = C.fit(X2, y2, st2)
print("with a 3,000-record stratum", b2, ok2, np.isfinite(ll2))
print("ok")
