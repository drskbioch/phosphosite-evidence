"""Propensity matching (D2) with a protein-cluster bootstrap that refits the model in every resample.

Same design as phosnet.matching.propensity_matched_comparison: logit(GOLD ~ index + species + residue type),
1:1 nearest-neighbour matching without replacement on the logit, caliper 0.2 SD of the logit (SD over all
analysis records), controls limited to the 50 nearest neighbours, treated records visited in random order.

Point estimate: the fit uses the same statsmodels routine as phosnet (Newton, BFGS fallback) on the record-level
data, so a design that phosnet cannot fit is reported as not fitted (cause recorded). Resamples: the covariates are
discrete, so the same maximum-likelihood fit is computed on (index, species, residue) cells with resample
multiplicities (Newton, same tolerance). Checked against phosnet in validate_propensity.py.
"""
import time

import numpy as np
import pandas as pd
from numba import njit
from scipy.spatial import cKDTree

CALIPER_SD = 0.2


@njit(cache=True)
def _greedy(order, dists, nn_idx, caliper, n_c):
    used = np.zeros(n_c, np.bool_)
    pt = np.empty(len(order), np.int64)
    pc = np.empty(len(order), np.int64)
    n = 0
    for ti in order:
        for r in range(dists.shape[1]):
            if dists[ti, r] > caliper:
                break
            ci = nn_idx[ti, r]
            if not used[ci]:
                used[ci] = True
                pt[n] = ti
                pc[n] = ci
                n += 1
                break
    return pt[:n], pc[:n]


def sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))


def effect(yt, yc, stat):
    if stat == "or":
        a, b = yt.sum(), len(yt) - yt.sum()
        c, d = yc.sum(), len(yc) - yc.sum()
        if min(a, b, c, d) <= 0:
            a, b, c, d = a + .5, b + .5, c + .5, d + .5
        return (a * d) / (b * c)
    return yt.mean() - yc.mean()


def smd(xt, xc):
    sd = np.sqrt((np.var(xt) + np.var(xc)) / 2)
    return 0.0 if sd == 0 else (xt.mean() - xc.mean()) / sd


def match_estimate(logit, treated, y, x, stat, rng, caliper_sd=CALIPER_SD):
    """logit, treated (0/1), y, x per record (all analysis records). Returns dict."""
    tr = np.flatnonzero(treated == 1)
    co = np.flatnonzero(treated == 0)
    caliper = caliper_sd * np.std(logit, ddof=1)
    k = min(50, len(co))
    tree = cKDTree(logit[co].reshape(-1, 1))
    dists, nn = tree.query(logit[tr].reshape(-1, 1), k=k)
    if k == 1:
        dists, nn = dists.reshape(-1, 1), nn.reshape(-1, 1)
    order = rng.permutation(len(tr))
    pt, pc = _greedy(order, dists, nn, caliper, len(co))
    n = len(pt)
    if n == 0:
        return {"n_pairs": 0, "estimate": np.nan, "smd": np.nan, "n_treated": len(tr)}
    it, ic = tr[pt], co[pc]
    return {"n_pairs": n, "n_treated": len(tr), "estimate": float(effect(y[it], y[ic], stat)),
            "smd": float(smd(x[it], x[ic])), "caliper": float(caliper)}


class PropData:
    """Record-level arrays of one analysis population (df columns: uniprot_id, species, residue_type, tier, y, x)."""

    def __init__(self, df, use_species=True, unit_col="uniprot_id"):
        self.y = df["y"].to_numpy(float)
        self.x = df["x"].to_numpy(float)
        self.t = (df["tier"] == "GOLD").to_numpy().astype(np.int8)
        self.prot, self.prot_names = pd.factorize(df[unit_col])
        self.U = len(self.prot_names)
        self.use_species = use_species
        sp = df["species"].astype(str).to_numpy() if use_species else np.zeros(len(df), str)
        rs = df["residue_type"].astype(str).to_numpy()
        self.sp_levels, spc = np.unique(sp, return_inverse=True)
        self.rs_levels, rsc = np.unique(rs, return_inverse=True)
        key = np.stack([self.x, spc.astype(float), rsc.astype(float)], axis=1)
        self.cell_keys, self.cell = np.unique(key, axis=0, return_inverse=True)
        self.cell = self.cell.ravel()
        C = len(self.cell_keys)
        cols = [np.ones(C), self.cell_keys[:, 0]]
        self.names = ["const", "x"]
        for j in range(1, len(self.sp_levels)):
            cols.append((self.cell_keys[:, 1] == j).astype(float)); self.names.append("sp:" + self.sp_levels[j])
        for j in range(1, len(self.rs_levels)):
            cols.append((self.cell_keys[:, 2] == j).astype(float)); self.names.append("rs:" + self.rs_levels[j])
        self.Xc = np.stack(cols, axis=1)
        self.n = len(df)

    # ---- point estimate: same routine as phosnet (statsmodels, record level)
    def fit_phosnet(self, df):
        import statsmodels.formula.api as smf
        w = pd.DataFrame({"_treated": self.t, "x": self.x, "residue_type": df["residue_type"].to_numpy(),
                          "species": df["species"].to_numpy()})
        formula = "_treated ~ x + C(species) + C(residue_type)" if self.use_species else "_treated ~ x + C(residue_type)"
        model = smf.logit(formula, data=w).fit(disp=0, method="newton", maxiter=200, tol=1e-8)
        if not model.mle_retvals.get("converged", True):
            model = smf.logit(formula, data=w).fit(disp=0, method="bfgs", maxiter=500)
        p = np.clip(np.asarray(model.predict(w)), 1e-6, 1 - 1e-6)
        return np.log(p / (1 - p))

    # ---- resample fit on cells
    def fit_cells(self, mult_rec, tol=1e-8, maxiter=200):
        """mult_rec: multiplicity of every record. Returns logit per cell, or None if the fit does not converge."""
        C = len(self.Xc)
        nt = np.bincount(self.cell, weights=mult_rec * self.t, minlength=C)
        nn = np.bincount(self.cell, weights=mult_rec, minlength=C)
        keep_rows = nn > 0
        X = self.Xc[keep_rows]
        nt_, nn_ = nt[keep_rows], nn[keep_rows]
        cols = np.flatnonzero((X != 0).any(axis=0) & ((X * nn_[:, None]).sum(axis=0) > 0))
        X = X[:, cols]
        beta = np.zeros(X.shape[1])
        for _ in range(maxiter):
            p = sigmoid(X @ beta)
            g = X.T @ (nt_ - nn_ * p)
            H = (X * (nn_ * p * (1 - p))[:, None]).T @ X
            step = np.linalg.lstsq(H, g, rcond=None)[0]
            beta = beta + step
            if np.max(np.abs(step)) < tol:
                eta = np.zeros(C)
                eta[keep_rows] = X @ beta
                eta[~keep_rows] = np.nan
                return eta
        return None

    def logit_from_cells(self, eta_cells):
        p = np.clip(sigmoid(eta_cells), 1e-6, 1 - 1e-6)
        return np.log(p / (1 - p))


def point(df, stat, rng, use_species=True):
    d = PropData(df, use_species)
    try:
        lg = d.fit_phosnet(df)
    except Exception as e:  # noqa: BLE001 - a design that cannot be fitted is a result
        return d, None, {"status": "not fitted", "reason": f"{type(e).__name__}: {e}"}
    r = match_estimate(lg, d.t, d.y, d.x, stat, rng)
    r.update(status="fitted", n_analysis=d.n)
    return d, lg, r


def crude(mult_rec, d, stat):
    g = mult_rec * d.t
    c = mult_rec * (1 - d.t)
    if stat == "or":
        A, C = (g * d.y).sum(), (c * d.y).sum()
        B, D = g.sum() - A, c.sum() - C
        if min(A, B, C, D) <= 0:
            A, B, C, D = A + .5, B + .5, C + .5, D + .5
        return (A * D) / (B * C)
    return (g * d.y).sum() / g.sum() - (c * d.y).sum() / c.sum()


def bootstrap(d, stat, rng, n_boot, time_log=None):
    """Protein-cluster bootstrap: resample proteins with replacement (all sites kept), refit, re-match."""
    out = {"est": np.full(n_boot, np.nan), "crude": np.full(n_boot, np.nan), "smd": np.full(n_boot, np.nan),
           "pairs": np.full(n_boot, np.nan), "unfit": 0}
    t0 = time.time()
    for b in range(n_boot):
        m = np.bincount(rng.integers(0, d.U, d.U), minlength=d.U)
        mult = m[d.prot].astype(float)
        eta = d.fit_cells(mult)
        if eta is None:
            out["unfit"] += 1
            continue
        rec = np.repeat(np.arange(d.n), m[d.prot])
        lg = d.logit_from_cells(eta[d.cell[rec]])
        r = match_estimate(lg, d.t[rec], d.y[rec], d.x[rec], stat, rng)
        if r["n_pairs"] == 0:
            out["unfit"] += 1
            continue
        out["est"][b], out["smd"][b], out["pairs"][b] = r["estimate"], r["smd"], r["n_pairs"]
        out["crude"][b] = crude(mult, d, stat)
        if time_log is not None and b in (2, 9):
            time_log.append((b + 1, round(time.time() - t0, 1)))
    return out
