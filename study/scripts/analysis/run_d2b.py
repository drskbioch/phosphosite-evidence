"""Block G18: design D2b (post-hoc matcher of the plan amendment of 2026-10-01; see d2b.py).

Point estimate: median of >= 20 order replicates (random visiting order and random tie-breaking; min, median and max reported).
Interval: gene-cluster bootstrap, 1,000 resamples; in each resample the logit model is refitted on the resampled records and the
matching is run once (one random order), so the interval includes the order variability. Crude estimate and attenuation on the same records.
A point estimate outside its own percentile interval is flagged and the contrast is not reported as decision-bearing.
Usage: python run_d2b.py G18-0123 ...  ->  results/d2b/<id>.json
"""
import json
import sys
import time

import numpy as np

import propensity as P
import run_cem_grid as R
from common import N_BOOT, RES, rng_for
from d2b import match_exact
from populations import SIGNALS

OUTD = RES / "d2b"
OUTD.mkdir(exist_ok=True)
N_REP = 20


def ci(a, lo=2.5, hi=97.5):
    a = a[~np.isnan(a)]
    return [float(np.percentile(a, lo)), float(np.percentile(a, hi))] if len(a) else [np.nan, np.nan]


def att(c, a, stat):
    c, a = np.asarray(c, float), np.asarray(a, float)
    if stat == "or":
        lc, la = np.log(c), np.log(a)
        return np.where(lc != 0, (lc - la) / lc, np.nan)
    return np.where(c != 0, (c - a) / c, np.nan)


def run(cid):
    row = R.man[R.man["contrast_id"] == cid].iloc[0]
    s = row["signal"]
    stat = SIGNALS[s]["stat"]
    st = "or" if stat == "or" else "diff"
    null = 1.0 if stat == "or" else 0.0
    tcol = R.tier_col_for(row["tier_rule"], s)
    d = R.population_for(row, tcol, R.IDX[row["index"]])
    n_all = len(d)
    d = d[d["residue_type"].isin(["S", "T", "Y"])]
    pdat = P.PropData(d, "species" in SIGNALS[s]["strata"], unit_col="gene_cluster")
    res = {"contrast_id": cid, "block": "G18", "signal": s, "population": row["population"] + " (S/T/Y only)", "design": "D2b", "index": row["index"],
           "tier_rule": row["tier_rule"], "null": null, "n_population_before_STY_restriction": int(n_all), "n_analysis": int(pdat.n),
           "n_gold": int(pdat.t.sum()), "n_bronze": int((1 - pdat.t).sum()), "n_gene_clusters": int(pdat.U), "n_proteins": int(d["uniprot_id"].nunique()),
           "resampling_unit": "gene cluster (gene_cluster)", "seed": f"[20260924, crc32({cid})]"}
    t0 = time.time()
    eta = pdat.fit_cells(np.ones(pdat.n))
    if eta is None:
        res.update(status="not fitted", reason="logit model did not converge on the S/T/Y records")
        json.dump(res, open(OUTD / f"{cid}.json", "w"), indent=1, default=float)
        print(cid, "not fitted", flush=True)
        return
    lg = pdat.logit_from_cells(eta[pdat.cell])
    rr = rng_for(cid + "|rep")
    reps = [match_exact(lg, pdat.t, pdat.y, pdat.x, st, rr) for _ in range(N_REP)]
    est = np.array([r["estimate"] for r in reps])
    pairs = np.array([r["n_pairs"] for r in reps])
    smds = np.array([r["smd"] for r in reps])
    point = float(np.median(est))
    crude_pt = float(P.crude(np.ones(pdat.n), pdat, st))
    rng = rng_for(cid)
    out = {k: np.full(N_BOOT, np.nan) for k in ("est", "crude", "smd", "pairs", "gold")}
    unfit = 0
    for b in range(N_BOOT):
        m = np.bincount(rng.integers(0, pdat.U, pdat.U), minlength=pdat.U)
        mult = m[pdat.prot].astype(float)
        e = pdat.fit_cells(mult)
        if e is None:
            unfit += 1
            continue
        rec = np.repeat(np.arange(pdat.n), m[pdat.prot])
        l2 = pdat.logit_from_cells(e[pdat.cell[rec]])
        r = match_exact(l2, pdat.t[rec], pdat.y[rec], pdat.x[rec], st, rng)
        if r["n_pairs"] == 0:
            unfit += 1
            continue
        out["est"][b], out["smd"][b], out["pairs"][b], out["gold"][b] = r["estimate"], r["smd"], r["n_pairs"], r["n_treated"]
        out["crude"][b] = P.crude(mult, pdat, st)
    a = att(out["crude"], out["est"], stat)
    cie = ci(out["est"])
    outside = bool(point < cie[0] or point > cie[1])
    res.update(status="run", estimate=point, estimate_rule="median of 20 order replicates", replicates={"n": N_REP, "min": float(est.min()), "median": point, "max": float(est.max()),
               "pairs_min": int(pairs.min()), "pairs_median": float(np.median(pairs)), "pairs_max": int(pairs.max()), "smd_min": float(smds.min()), "smd_max": float(smds.max())},
               ci=cie, crude=crude_pt, crude_ci=ci(out["crude"]), survives=bool(cie[0] > null or cie[1] < null),
               attenuation=float(att(crude_pt, point, stat)), attenuation_ci=ci(a), n_pairs=float(np.median(pairs)), n_gold=int(pdat.t.sum()),
               share_gold_matched=float(np.median(pairs) / pdat.t.sum()), smd=float(np.median(smds)), balance_pass=bool(abs(np.median(smds)) < 0.1),
               point_outside_own_interval=outside, decision_bearing="no" if (outside or row["decision_bearing"] == "no") else "yes",
               n_boot=N_BOOT, n_boot_valid=int((~np.isnan(out["est"])).sum()), n_boot_unfit=unfit,
               boot_pairs_range=[float(np.nanmin(out["pairs"])), float(np.nanmax(out["pairs"]))], seconds=round(time.time() - t0, 1))
    json.dump(res, open(OUTD / f"{cid}.json", "w"), indent=1, default=float)
    print(cid, s, round(point, 4), [round(v, 4) for v in cie], "reps", round(est.min(), 4), round(est.max(), 4), "outside" if outside else "", res["seconds"], flush=True)


if __name__ == "__main__":
    for c in sys.argv[1:]:
        run(c)
