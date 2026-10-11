"""Block G21: within-gene GOLD vs BRONZE (definition in the plan, deviation log of 2026-10-04, fourth batch; manifest rows G21-0491 to G21-0514).

Per signal-by-population group: models M1 (all S/T/Y records), M1s (records with a protein sequence), M2 (M1s plus the detectability indicator).
Binary outcomes: conditional logistic regression with a stratum per gene cluster (exposure coefficient -> odds ratio) and, as a same-scale companion,
the within-gene linear fixed-effect difference in proportions. Disorder: within-gene linear fixed-effect difference.
Gene-cluster bootstrap: every resample draws the multiplicity of each gene cluster of the group's population once, shared by the three models and both estimators;
all models are refitted in every resample (conditional logit warm-started at the full-sample estimate). Generator np.random.default_rng([20260924, crc32(id of the M2 row)]).
Decision D7 applies: if the projected time of a group's bootstrap (from the first resamples) exceeds 2 hours, the group uses 200 resamples, stated in the result.
Usage: python run_g21.py [S2,S3,S4,S5] [P-all,P-human]   -> results/g21_results.json (merged per group)
"""
import json
import sys
import time

import numpy as np
import pandas as pd

import cond_logit as C
import g21_data as G
import run_cem_grid as R
from common import N_BOOT, RES, rng_for

OUT = RES / "g21_results.json"
man = R.man[R.man["block"] == "G21"]
BINARY = {"S2", "S4", "S5"}


class Model:
    def __init__(self, d, cols, mask, cluster_code, K, binary=True):
        s = d.loc[mask].reset_index(drop=True)
        self.n_records = len(s)
        code = cluster_code[mask.to_numpy() if hasattr(mask, "to_numpy") else mask]
        g = s["gene_cluster"].to_numpy()
        first = np.concatenate([[0], np.nonzero(g[1:] != g[:-1])[0] + 1])
        self.starts = np.append(first, len(s)).astype(np.int64)
        self.clus = code[first]  # cluster index of each stratum
        self.X = s[cols].to_numpy(float)
        self.y = s["y"].to_numpy(float)
        self.gold = s["gold"].to_numpy()
        n = np.diff(self.starts)
        self.n_strata = len(n)
        gmin = np.minimum.reduceat(self.gold, first)
        gmax = np.maximum.reduceat(self.gold, first)
        both = (gmin == 0) & (gmax == 1)
        self.n_both = int(both.sum())
        sg = np.add.reduceat(self.gold.astype(int), first)
        self.gold_in_both = int(sg[both].sum())
        self.bronze_in_both = int((n[both] - sg[both]).sum())
        # linear fixed-effect sufficient statistics per stratum (within-stratum demeaned)
        rep = np.repeat(np.arange(len(n)), n)
        Xm = np.add.reduceat(self.X, first, axis=0) / n[:, None]
        ym = np.add.reduceat(self.y, first) / n
        Xd = self.X - Xm[rep]
        yd = self.y - ym[rep]
        p = self.X.shape[1]
        self.Sxx = np.stack([np.add.reduceat(Xd[:, a] * Xd[:, c], first) for a in range(p) for c in range(p)], 1).reshape(len(n), p, p)
        self.Sxy = np.stack([np.add.reduceat(Xd[:, a] * yd, first) for a in range(p)], 1)
        # conditional logit: informative strata only (outcome varies)
        if not binary:
            self.n_inf = 0
            return
        self.yi = self.y.astype(np.int64)
        m = np.add.reduceat(self.yi, first)
        inf = (m > 0) & (m < n)
        self.inf_idx = np.nonzero(inf)[0]
        keep_rows = np.concatenate([np.arange(self.starts[i], self.starts[i + 1]) for i in self.inf_idx]) if len(self.inf_idx) else np.array([], int)
        self.Xc = np.ascontiguousarray(self.X[keep_rows])
        self.yc = self.yi[keep_rows]
        sizes = n[self.inf_idx]
        self.startsc = np.append(0, np.cumsum(sizes)).astype(np.int64)
        self.clus_c = self.clus[self.inf_idx]
        self.n_inf = int(inf.sum())

    def linear(self, mult):
        w = mult[self.clus]
        A = np.tensordot(w, self.Sxx, axes=(0, 0))
        bvec = w @ self.Sxy
        return float(np.linalg.solve(A, bvec)[0])

    def logit(self, mult, b0=None):
        w = mult[self.clus_c].astype(float)
        b, ok, ll, cov = C.fit(self.Xc, self.yc, self.startsc, w, parallel=True, b0=b0)
        return b, ok


def percentile_ci(v):
    v = np.asarray(v)
    v = v[np.isfinite(v)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if len(v) else [float("nan")] * 2


def run_group(signal, popl, store):
    key = f"{popl} ({'S/T/Y'})" if "S/T/Y" not in popl else popl
    rows = man[(man["signal"] == signal) & (man["population"] == key)]
    ids = {r["variation"].split(":")[0]: r["contrast_id"] for _, r in rows.iterrows()}
    d, info = G.prepare(signal, key)
    cols3, cols4 = ["gold", "isT", "isY"], ["gold", "isT", "isY", "det"]
    codes, uniq = pd.factorize(d["gene_cluster"])
    K = len(uniq)
    allm = np.ones(len(d), bool)
    seqm = d["has_seq"].to_numpy()
    binary = signal in BINARY
    models = {"M1": Model(d, cols3, allm, codes, K, binary), "M1s": Model(d, cols3, seqm, codes, K, binary), "M2": Model(d, cols4, seqm, codes, K, binary)}
    one = np.ones(K)
    full = {}
    b0 = {}
    for k, m in models.items():
        r = {"linear": m.linear(one)}
        if binary:
            b, ok = m.logit(one)
            b0[k] = b
            r["logit_coef"] = float(b[0])
            r["logit_ok"] = bool(ok)
        full[k] = r
    rng = rng_for(ids["M2"])
    nb_target = N_BOOT
    boots = {k: {"linear": [], "logit": []} for k in models}
    unfit = {k: 0 for k in models}
    t0 = time.time()
    i = 0
    nb = nb_target
    while i < nb:
        mult = np.bincount(rng.integers(0, K, K), minlength=K).astype(float)
        for k, m in models.items():
            try:
                boots[k]["linear"].append(m.linear(mult))
            except np.linalg.LinAlgError:
                boots[k]["linear"].append(np.nan)
            if binary:
                b, ok = m.logit(mult, b0=b0[k])
                if ok:
                    boots[k]["logit"].append(float(b[0]))
                else:
                    boots[k]["logit"].append(np.nan)
                    unfit[k] += 1
        i += 1
        if i == 5 and nb == nb_target:
            proj = (time.time() - t0) / 5 * nb_target
            if proj > 7200:
                nb = 200
                print(f"{signal} {popl}: projected {proj / 3600:.1f} h for {nb_target} resamples exceeds 2 h; decision D7: {nb} resamples", flush=True)
        if i % 100 == 0:
            print(f"{signal} {popl}: {i}/{nb} resamples, {time.time() - t0:.0f} s", flush=True)
    secs = round(time.time() - t0, 1)
    out = {}
    for k, m in models.items():
        cid = ids[k]
        base = {"signal": signal, "population": key, "model": k, "design": "D6 (within-gene)", "n_records": m.n_records, "n_gene_clusters": m.n_strata,
                "n_gene_clusters_with_both_tiers": m.n_both, "gold_in_gene_clusters_with_both_tiers": m.gold_in_both, "bronze_in_gene_clusters_with_both_tiers": m.bronze_in_both,
                "n_boot": nb, "d7_fallback": bool(nb != nb_target), "resampling_unit": "gene cluster", "seed": f"[20260924, crc32({ids['M2']})]", "seconds_group": secs,
                "population_info": info, "decision_bearing": man[man["contrast_id"] == cid]["decision_bearing"].iloc[0], "status": "run"}
        lin = np.array(boots[k]["linear"])
        if binary:
            lg = np.array(boots[k]["logit"])
            base.update(estimator="conditional logistic regression, odds ratio of GOLD vs BRONZE within gene", null=1.0, estimate=float(np.exp(full[k]["logit_coef"])),
                        ci=[float(np.exp(x)) for x in percentile_ci(lg)], n_informative_gene_clusters=m.n_inf, converged_full=full[k]["logit_ok"],
                        n_boot_valid=int(np.isfinite(lg).sum()), n_boot_unfit=int(unfit[k]),
                        companion_linear=dict(estimator="within-gene linear fixed effect, difference in proportions", null=0.0, estimate=full[k]["linear"], ci=percentile_ci(lin),
                                              n_boot_valid=int(np.isfinite(lin).sum())))
            base["persists"] = bool(base["ci"][0] > 1 or base["ci"][1] < 1)
        else:
            base.update(estimator="within-gene linear fixed effect, difference in means", null=0.0, estimate=full[k]["linear"], ci=percentile_ci(lin), n_boot_valid=int(np.isfinite(lin).sum()), n_boot_unfit=0)
            base["persists"] = bool(base["ci"][0] > 0 or base["ci"][1] < 0)
        out[cid] = base
        print(cid, signal, popl, k, round(base["estimate"], 4), [round(x, 4) for x in base["ci"]], base["persists"], flush=True)
    cur = json.load(open(OUT)) if OUT.exists() else {}
    cur.update(out)
    json.dump(cur, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    sigs = sys.argv[1].split(",") if len(sys.argv) > 1 else ["S2", "S3", "S4", "S5"]
    pops = sys.argv[2].split(",") if len(sys.argv) > 2 else ["P-all", "P-human"]
    chk = json.load(open(RES / "g21_cond_logit_check.json"))
    if not chk["passes"]:
        raise SystemExit("cond_logit check against statsmodels did not pass; G21 is not run (plan definition)")
    done = json.load(open(OUT)) if OUT.exists() else {}
    for popl in pops:
        for s in sigs:
            ids = man[(man["signal"] == s) & (man["population"] == f"{popl} (S/T/Y)")]["contrast_id"].tolist()
            if all(i in done for i in ids):
                print("already run:", s, popl, flush=True)
                continue
            run_group(s, popl, None)
