"""D2b (post-hoc, plan amendment of 2026-10-01): propensity matching as the Rev 5 Methods describe it.

Model: logit(GOLD ~ index + species + residue type) on serine, threonine and tyrosine records. Matching: 1:1 without replacement,
every GOLD record (random visiting order) is paired with the nearest still unused BRONZE record on the logit if the distance is within the
caliper (0.2 x SD of the logit over the analysis records); ties between records with identical logits are broken at random; there is no cap on
the number of candidate controls (the Rev 5 routine looked only at the 50 nearest controls of each GOLD record).
"""
import numpy as np
from numba import njit

from propensity import effect, smd


@njit(cache=True)
def _find(parent, i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


@njit(cache=True)
def _greedy_exact(tl, order, coins, uv, start, cnt0, perm, caliper):
    V = len(uv)
    cnt = cnt0.copy()
    Lp = np.arange(V + 2)
    Rp = np.arange(V + 2)
    pt = np.empty(len(order), np.int64)
    pc = np.empty(len(order), np.int64)
    n = 0
    for oi in order:
        t = tl[oi]
        p = np.searchsorted(uv, t)  # number of values < t
        r = _find(Rp, p + 1)
        l = _find(Lp, p)
        dl = t - uv[l - 1] if l >= 1 else np.inf
        dr = uv[r - 1] - t if r <= V else np.inf
        if dl < dr:
            v, d = l, dl
        elif dr < dl:
            v, d = r, dr
        else:
            if dl == np.inf:
                continue
            v = l if coins[oi] < 0.5 else r
            d = dl
        if d > caliper:
            continue
        taken = cnt0[v - 1] - cnt[v - 1]
        pt[n] = oi
        pc[n] = perm[start[v - 1] + taken]
        n += 1
        cnt[v - 1] -= 1
        if cnt[v - 1] == 0:
            Lp[v] = v - 1
            Rp[v] = v + 1
    return pt[:n], pc[:n]


def match_exact(logit, treated, y, x, stat, rng, caliper_sd=0.2):
    tr = np.flatnonzero(treated == 1)
    co = np.flatnonzero(treated == 0)
    caliper = caliper_sd * np.std(logit, ddof=1)
    lc = logit[co]
    key = rng.random(len(co))
    perm = np.lexsort((key, lc))
    sl = lc[perm]
    uv, start, cnt = np.unique(sl, return_index=True, return_counts=True)
    tl = logit[tr]
    order = rng.permutation(len(tr))
    coins = rng.random(len(tr))
    pt, pc = _greedy_exact(tl, order.astype(np.int64), coins, uv, start.astype(np.int64), cnt.astype(np.int64), co[perm].astype(np.int64), caliper)
    n = len(pt)
    if n == 0:
        return {"n_pairs": 0, "estimate": np.nan, "smd": np.nan, "n_treated": len(tr)}
    it, ic = tr[pt], pc
    return {"n_pairs": n, "n_treated": len(tr), "estimate": float(effect(y[it], y[ic], stat)), "smd": float(smd(x[it], x[ic])), "caliper": float(caliper)}
