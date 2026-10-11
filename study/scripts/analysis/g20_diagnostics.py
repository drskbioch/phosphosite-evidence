"""Block G20 (plan amendment of 2026-10-01): positive control and index diagnostics. Gene-cluster bootstrap, 1,000 resamples.
 a. pre-matching balance (SMD of the index, GOLD vs BRONZE, before matching) for every D1 / D3 contrast that has a result
 b. crude association (Pearson r, equal to the point-biserial r for binary outcomes) of the primary index with each outcome within GOLD and within BRONZE records
 c. protein level: Spearman correlation of the primary index with the number of sites per accession (all, GOLD only)
 d. domain overlap by the protein's domain-annotation coverage (no annotated site / at least one annotated site)
Output: results/g20_diagnostics.json
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

import run_cem_grid as R
from common import N_BOOT, RES, rng_for
from engine import Cells, point
from populations import SIGNALS, population

b = R.b
out = {}

# ---- a. pre-matching balance
a = {}
for fn in ("cem_results.json", "amend_results.json"):
    try:
        res = json.load(open(RES / fn))
    except FileNotFoundError:
        continue
    for cid, v in res.items():
        if v.get("design") in ("D1", "D3") and v.get("status") == "run":
            if "smd_pre" in v:
                a[cid] = {"smd_pre": v["smd_pre"], "smd_after": v.get("smd"), "signal": v["signal"], "population": v["population"], "block": cid.split("-")[0]}
                continue
            row = R.man[R.man["contrast_id"] == cid].iloc[0]
            s = row["signal"]
            tcol = R.tier_col_for(row["tier_rule"], s)
            d = R.population_for(row, tcol, R.IDX[row["index"]])
            c = Cells(d, "y", "x", {"gene": ["gene_cluster"]}, strata_cols=SIGNALS[s]["strata"])
            p = point(c, SIGNALS[s]["stat"])
            a[cid] = {"smd_pre": p["smd_pre"], "smd_after": v.get("smd"), "signal": s, "population": v["population"], "block": cid.split("-")[0]}
out["prematching_balance"] = a


# ---- b. association of the index with the outcome within tier
def pearson_boot(x, y, unit, rng):
    codes, uniq = pd.factorize(unit)
    U = len(uniq)
    mom = {k: np.bincount(codes, weights=v, minlength=U) for k, v in {"n": np.ones(len(x)), "x": x, "y": y, "xy": x * y, "xx": x * x, "yy": y * y}.items()}

    def r_of(w):
        n, sx, sy, sxy, sxx, syy = (w @ mom[k] for k in ("n", "x", "y", "xy", "xx", "yy"))
        vx, vy = sxx / n - (sx / n) ** 2, syy / n - (sy / n) ** 2
        return (sxy / n - sx * sy / n ** 2) / np.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan
    pt = r_of(np.ones(U))
    bs = np.array([r_of(np.bincount(rng.integers(0, U, U), minlength=U).astype(float)) for _ in range(N_BOOT)])
    bs = bs[~np.isnan(bs)]
    return float(pt), [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], int(U)


bm = R.man[(R.man["block"] == "G20") & (R.man["design"] == "index-outcome association")]
assoc = {}
for _, row in bm.iterrows():
    s = row["signal"]
    tier = "GOLD" if "GOLD" in row["variation"] else "BRONZE"
    kw = {"species": "human"} if row["population"] == "P-human" else {}
    d = population(b, s, tier_col="tier", index_col="SEI_V50", **kw)
    d = d[d["tier"] == tier]
    pt, ci_, U = pearson_boot(d["x"].to_numpy(float), d["y"].to_numpy(float), d["gene_cluster"].to_numpy(), rng_for(row["contrast_id"]))
    assoc[row["contrast_id"]] = {"signal": s, "population": row["population"], "tier": tier, "n_records": int(len(d)), "n_gene_clusters": U,
                                 "pearson_r": pt, "ci": ci_, "mean_outcome": float(d["y"].mean()), "index_mean": float(d["x"].mean()), "index_sd": float(d["x"].std())}
    print(row["contrast_id"], s, row["population"], tier, round(pt, 4), [round(v, 4) for v in ci_], flush=True)
out["index_outcome_association"] = assoc

# ---- c. index vs number of sites per accession
sub = b[b["SEI_V50"].notna()]
cnt = sub.groupby("uniprot_id").agg(n_sites=("tier", "size"), n_gold=("tier", lambda t: int((t == "GOLD").sum())), index=("SEI_V50", "first"), species=("species", "first"),
                                    gene_cluster=("gene_cluster", "first")).reset_index()


def spearman_boot(x, y, unit, rng):
    codes, uniq = pd.factorize(unit)
    U = len(uniq)
    idx = [np.flatnonzero(codes == u) for u in range(U)]
    pt = stats.spearmanr(x, y)[0]
    bs = []
    for _ in range(N_BOOT):
        sel = np.concatenate([idx[i] for i in rng.integers(0, U, U)])
        bs.append(stats.spearmanr(x[sel], y[sel])[0])
    return float(pt), [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], int(U)


spp = {}
for sp in ("human", "all"):
    cc = cnt if sp == "all" else cnt[cnt["species"] == sp]
    for name, col in (("sites_per_accession", "n_sites"), ("gold_sites_per_accession", "n_gold")):
        pt, ci_, U = spearman_boot(cc["index"].to_numpy(), cc[col].to_numpy(float), cc["gene_cluster"].to_numpy(), rng_for(f"G20|spp|{sp}|{col}"))
        spp[f"{sp}|{name}"] = {"spearman": pt, "ci": ci_, "n_accessions": int(len(cc)), "n_gene_clusters": U}
        print("G20 c", sp, name, round(pt, 3), [round(v, 3) for v in ci_], flush=True)
out["index_vs_sites_per_accession"] = spp

# ---- d. domain overlap by annotation coverage
d = population(b, "S1", tier_col="tier", index_col="SEI_V50")
d["annot"] = d["uniprot_id"].isin(R.ANNOT)


def or_boot(df, rng):
    codes, uniq = pd.factorize(df["gene_cluster"])
    U = len(uniq)
    g = (df["tier"] == "GOLD").to_numpy()
    y = df["y"].to_numpy()
    cnts = np.stack([np.bincount(codes, weights=v, minlength=U) for v in (g & (y == 1), g & (y == 0), ~g & (y == 1), ~g & (y == 0))])

    def o(w):
        A, B, C, D = cnts @ w
        if min(A, B, C, D) <= 0:
            return np.nan
        return A * D / (B * C)
    pt = o(np.ones(U))
    bs = np.array([o(np.bincount(rng.integers(0, U, U), minlength=U).astype(float)) for _ in range(N_BOOT)])
    bs = bs[~np.isnan(bs)]
    return pt, ([float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))] if len(bs) else [np.nan, np.nan]), cnts.sum(axis=1).tolist(), int(U)


dom = {}
for name, sel in (("all proteins", d["annot"] | ~d["annot"]), ("at least one domain-annotated site", d["annot"]), ("no domain-annotated site", ~d["annot"])):
    x = d[sel]
    g = x["tier"] == "GOLD"
    pt, ci_, c4, U = or_boot(x, rng_for(f"G20|dom|{name}"))
    dom[name] = {"n_records": int(len(x)), "n_gold": int(g.sum()), "n_bronze": int((~g).sum()), "n_gene_clusters": U, "share_in_domain_gold": float(x.loc[g, "y"].mean()),
                 "share_in_domain_bronze": float(x.loc[~g, "y"].mean()), "crude_or": float(pt) if pt == pt else None, "ci": ci_,
                 "note": "outcome is constant (0) in this stratum, odds ratio undefined" if pt != pt else ""}
    print("G20 d", name, dom[name]["n_gold"], dom[name]["n_bronze"], pt, ci_, flush=True)
# composition: share of GOLD and of BRONZE records in annotated proteins, odds ratio of being in an annotated protein
g = d["tier"] == "GOLD"
comp = d.assign(y_orig=d["y"], y=d["annot"].astype(float))
pt, ci_, c4, U = or_boot(comp, rng_for("G20|dom|composition"))
dom["composition"] = {"share_gold_records_in_annotated_proteins": float(d.loc[g, "annot"].mean()), "share_bronze_records_in_annotated_proteins": float(d.loc[~g, "annot"].mean()),
                      "or_gold_vs_bronze_of_being_in_an_annotated_protein": float(pt), "ci": ci_}
out["domain_overlap_by_annotation_coverage"] = dom
json.dump(out, open(RES / "g20_diagnostics.json", "w"), indent=1, default=float)
