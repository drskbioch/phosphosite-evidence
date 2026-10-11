"""Coarsened exact matching with a resampling-unit bootstrap that re-runs the whole pipeline.

Cells are records aggregated over (protein, residue type, tier [, other keys]). For a given resample each
cell has multiplicity m[unit of the cell]. Inside every resample the index deciles are recomputed from the
resampled records (weighted quantile identical to pd.qcut on the expanded data, duplicate edges dropped),
strata (species x residue type x decile) are rebuilt, and the estimate, the crude estimate on the same
records, the balance and the attenuation are recomputed. Point estimates reproduce phosnet's
matched_comparison (checked in validate_engine.py).

stat: "or" (odds ratio, 0.5 correction when a cell count is 0 as in phosnet) | "diff" (difference in means).
"""
import numpy as np
import pandas as pd

from common import N_BOOT


class Cells:
    def __init__(self, df, y, x, unit_cols, strata_cols=("species", "residue_type"), tier_col="tier", z=None):
        """df: record-level frame restricted to the analysis population (GOLD and BRONZE, index and outcome defined).
        y, x: names of the outcome and index columns. unit_cols: dict scheme -> list of columns defining the unit."""
        d = df
        keys = ["uniprot_id", "residue_type", tier_col] + [c for c in unit_cols_all(unit_cols) if c not in ("uniprot_id",)]
        keys = list(dict.fromkeys(keys + list(strata_cols)))
        g = d.groupby(keys, sort=False, observed=True)
        spec = dict(n=(y, "size"), sy=(y, "sum"), x=(x, "first"))
        if z is not None:
            spec["z"] = (z, "first")
        agg = g.agg(**spec).reset_index()
        agg["t"] = (agg[tier_col] == "GOLD").astype(np.int8)
        sp = pd.factorize(agg[list(strata_cols)].astype(str).agg("|".join, axis=1))[0] if len(strata_cols) else np.zeros(len(agg), int)
        self.sp = sp.astype(np.int64)
        self.n = agg["n"].to_numpy(float)
        self.sy = agg["sy"].to_numpy(float)
        self.x = agg["x"].to_numpy(float)
        self.t = agg["t"].to_numpy(np.int8)
        self.units = {}
        for name, cols in unit_cols.items():
            self.units[name] = pd.factorize(agg[cols].astype(str).agg("|".join, axis=1))[0].astype(np.int64)
        self.n_units = {k: int(v.max()) + 1 for k, v in self.units.items()}
        self.order = np.argsort(self.x, kind="stable")
        self.xs = self.x[self.order]
        self.z = agg["z"].to_numpy(float) if z is not None else None
        if z is not None:
            self.order_z = np.argsort(self.z, kind="stable")
            self.zs = self.z[self.order_z]
        self.n_sites = int(self.n.sum())
        self.n_sp = int(self.sp.max()) + 1
        self.cells = agg


def unit_cols_all(unit_cols):
    out = []
    for cols in unit_cols.values():
        out.extend(cols)
    return out


def _edges(xs, cw_sorted, q=10):
    cum = np.cumsum(cw_sorted)
    N = cum[-1]
    hs = (N - 1) * np.linspace(0, 1, q + 1)
    lo, hi = np.floor(hs).astype(np.int64), np.ceil(hs).astype(np.int64)
    vlo = xs[np.minimum(np.searchsorted(cum, lo, side="right"), len(xs) - 1)]
    vhi = xs[np.minimum(np.searchsorted(cum, hi, side="right"), len(xs) - 1)]
    return np.unique(vlo + (hs - lo) * (vhi - vlo))


def _decile(v, vs, order, cw, n_dec):
    edges = _edges(vs, cw[order], n_dec)
    if len(edges) < 2:
        return np.zeros(len(v), dtype=np.int64), 1
    # zero-weight cells (proteins absent from a resample) can lie above the top edge; clipping keeps their stratum index in range (their weight is 0)
    return np.clip(np.searchsorted(edges, v, side="left") - 1, 0, len(edges) - 2), len(edges) - 1


def _estimate(c, w, stat, n_dec=10, mode="x"):
    """mode "x": strata = strata cols x decile of x (the validated Rev 5 design); "z": strata x decile of z; "xz": strata x decile of x x decile of z.
    Deciles are recomputed from the weighted (resampled) records."""
    cw = w * c.n
    if cw.sum() == 0:
        return None
    dec, nb = _decile(c.x, c.xs, c.order, cw, n_dec)
    if mode == "x":
        sid = c.sp * n_dec + dec
        L = c.n_sp * n_dec
    else:
        decz, nbz = _decile(c.z, c.zs, c.order_z, cw, n_dec)
        if mode == "z":
            sid, L, nb = c.sp * n_dec + decz, c.n_sp * n_dec, nbz
        else:
            sid, L, nb = (c.sp * n_dec + dec) * n_dec + decz, c.n_sp * n_dec * n_dec, (nb, nbz)
    t1 = c.t == 1
    Nt = np.bincount(sid[t1], weights=cw[t1], minlength=L)
    Nc = np.bincount(sid[~t1], weights=cw[~t1], minlength=L)
    SYt = np.bincount(sid[t1], weights=(w * c.sy)[t1], minlength=L)
    SYc = np.bincount(sid[~t1], weights=(w * c.sy)[~t1], minlength=L)
    SXt = np.bincount(sid[t1], weights=(cw * c.x)[t1], minlength=L)
    SXc = np.bincount(sid[~t1], weights=(cw * c.x)[~t1], minlength=L)
    SX2t = np.bincount(sid[t1], weights=(cw * c.x * c.x)[t1], minlength=L)
    SX2c = np.bincount(sid[~t1], weights=(cw * c.x * c.x)[~t1], minlength=L)
    if c.z is not None:
        SZt = np.bincount(sid[t1], weights=(cw * c.z)[t1], minlength=L)
        SZc = np.bincount(sid[~t1], weights=(cw * c.z)[~t1], minlength=L)
        SZ2t = np.bincount(sid[t1], weights=(cw * c.z * c.z)[t1], minlength=L)
        SZ2c = np.bincount(sid[~t1], weights=(cw * c.z * c.z)[~t1], minlength=L)
    v = (Nt > 0) & (Nc > 0)
    if not v.any():
        return None
    f = np.where(v, Nt / np.where(Nc > 0, Nc, 1), 0.0)
    T = Nt[v].sum()
    # adjusted
    a = SYt[v].sum()
    c_ = (f * SYc)[v].sum()
    if stat == "or":
        b, d = T - a, T - c_
        if min(a, b, c_, d) <= 0:
            a, b, c_, d = a + .5, b + .5, c_ + .5, d + .5
        adj = (a * d) / (b * c_)
        # crude on the same (resampled) records
        A, C = SYt.sum(), SYc.sum()
        B, D = Nt.sum() - A, Nc.sum() - C
        if min(A, B, C, D) <= 0:
            A, B, C, D = A + .5, B + .5, C + .5, D + .5
        crude = (A * D) / (B * C)
    else:
        adj = a / T - c_ / T
        crude = SYt.sum() / Nt.sum() - SYc.sum() / Nc.sum()
    mt, mc = SXt[v].sum() / T, (f * SXc)[v].sum() / T
    vt = SX2t[v].sum() / T - mt ** 2
    vc = (f * SX2c)[v].sum() / T - mc ** 2
    sd = np.sqrt((vt + vc) / 2)
    smd = (mt - mc) / sd if sd > 0 else 0.0
    # balance before matching: all GOLD and all BRONZE records of the population (record weights)
    m1, m0 = SXt.sum() / Nt.sum(), SXc.sum() / Nc.sum()
    sd_pre = np.sqrt(((SX2t.sum() / Nt.sum() - m1 ** 2) + (SX2c.sum() / Nc.sum() - m0 ** 2)) / 2)
    smd_pre = (m1 - m0) / sd_pre if sd_pre > 0 else 0.0
    smd_z = np.nan
    if c.z is not None:
        zt, zc = SZt[v].sum() / T, (f * SZc)[v].sum() / T
        szd = np.sqrt((SZ2t[v].sum() / T - zt ** 2 + (f * SZ2c)[v].sum() / T - zc ** 2) / 2)
        smd_z = (zt - zc) / szd if szd > 0 else 0.0
    return {"adj": adj, "crude": crude, "smd": smd, "smd_pre": smd_pre, "smd_z": smd_z, "n_strata_used": int(v.sum()),
            "n_strata_dropped": int(((Nt > 0) | (Nc > 0)).sum() - v.sum()),
            "n_gold_used": float(T), "n_bins": nb}


def _att(crude, adj, stat):
    if stat == "or":
        lc, la = np.log(crude), np.log(adj)
        return (lc - la) / lc if lc != 0 else np.nan
    return (crude - adj) / crude if crude != 0 else np.nan


def point(c, stat, mode="x"):
    r = _estimate(c, np.ones(len(c.n)), stat, mode=mode)
    if r is not None:
        r["attenuation"] = _att(r["crude"], r["adj"], stat)
    return r


def bootstrap(c, stat, scheme, rng, n_boot=N_BOOT, mode="x"):
    """Resample units of `scheme` (with replacement, same number as observed) and re-run everything."""
    u = c.units[scheme]
    U = c.n_units[scheme]
    out = {k: np.full(n_boot, np.nan) for k in ("adj", "crude", "smd", "smd_z", "attenuation")}
    for b in range(n_boot):
        m = np.bincount(rng.integers(0, U, U), minlength=U).astype(float)
        r = _estimate(c, m[u], stat, mode=mode)
        if r is None:
            continue
        out["adj"][b], out["crude"][b], out["smd"][b], out["smd_z"][b] = r["adj"], r["crude"], r["smd"], r["smd_z"]
        out["attenuation"][b] = _att(r["crude"], r["adj"], stat)
    return out


def ci(a, lo=2.5, hi=97.5):
    a = a[~np.isnan(a)]
    if len(a) == 0:
        return [np.nan, np.nan]
    return [float(np.percentile(a, lo)), float(np.percentile(a, hi))]


def summarize(c, stat, scheme, rng, n_boot=N_BOOT, mode="x"):
    """Point estimates plus percentile 95% CIs from the resampling-unit bootstrap."""
    p = point(c, stat, mode)
    if p is None:
        return {"error": "no stratum contains both tiers"}
    bs = bootstrap(c, stat, scheme, rng, n_boot, mode)
    null = 1.0 if stat == "or" else 0.0
    adj_ci = ci(bs["adj"])
    res = {"stat": stat, "unit": scheme, "n_boot": n_boot, "n_boot_valid": int((~np.isnan(bs["adj"])).sum()),
           "adjusted": p["adj"], "adjusted_ci": adj_ci, "crude": p["crude"], "crude_ci": ci(bs["crude"]),
           "attenuation": p["attenuation"], "attenuation_ci": ci(bs["attenuation"]),
           "smd": p["smd"], "smd_pre": p["smd_pre"], "smd_z": p["smd_z"], "balance_pass": bool(abs(p["smd"]) < 0.1),
           "n_strata_used": p["n_strata_used"], "n_strata_dropped": p["n_strata_dropped"], "n_bins": p["n_bins"],
           "n_sites": c.n_sites, "n_sites_gold": int(c.n[c.t == 1].sum()), "n_sites_bronze": int(c.n[c.t == 0].sum()),
           "n_proteins": int(len(np.unique(c.cells["uniprot_id"]))), "n_units": c.n_units[scheme],
           "survives": bool(adj_ci[0] > null or adj_ci[1] < null),
           "adjusted_boot_sd": float(np.nanstd(bs["adj"], ddof=1))}
    # decision labels from section 2 of the plan
    res["largely_explained"] = bool((p["attenuation"] is not None) and (not np.isnan(p["attenuation"])) and p["attenuation"] >= 0.5
                                    and ((res["crude_ci"][0] > null) or (res["crude_ci"][1] < null)))
    return res, bs
