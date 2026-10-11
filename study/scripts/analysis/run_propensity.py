"""Run propensity-matching (D2 / R5_propensity) contrasts with the protein-cluster bootstrap (refit in every resample).

Usage: python run_propensity.py <contrast_id>[+posthoc_STY] ...   -> results/prop/<id>.json
D7: 1,000 resamples; if the first 10 resamples project to more than 2 hours for 1,000, the run stops at 200 and says so.
A design that phosnet's fitting routine cannot fit is reported as not fitted with the reason (no bootstrap).
`+posthoc_STY` runs the same design on serine, threonine and tyrosine records only (post-hoc, labelled).
"""
import json
import sys
import time

import numpy as np

import propensity as P
import run_cem_grid as R
from common import N_BOOT, RES, rng_for
from populations import SIGNALS

OUTD = RES / "prop"
OUTD.mkdir(exist_ok=True)


def ci(a, lo=2.5, hi=97.5):
    a = a[~np.isnan(a)]
    return [float(np.percentile(a, lo)), float(np.percentile(a, hi))] if len(a) else [np.nan, np.nan]


def att(c, a, stat):
    if stat == "or":
        lc, la = np.log(c), np.log(a)
        return np.where(lc != 0, (lc - la) / lc, np.nan)
    return np.where(c != 0, (c - a) / c, np.nan)


def run(spec):
    cid, _, suffix = spec.partition("+")
    row = R.man[R.man["contrast_id"] == cid].iloc[0]
    s = row["signal"]
    stat = SIGNALS[s]["stat"]
    null = 1.0 if stat == "or" else 0.0
    tcol = R.tier_col_for(row["tier_rule"], s)
    icol = R.IDX[row["index"]]
    d = R.population_for(row, tcol, icol)
    posthoc = suffix == "posthoc_STY"
    if posthoc:
        d = d[d["residue_type"].isin(["S", "T", "Y"])]
    use_species = "species" in SIGNALS[s]["strata"]
    t0 = time.time()
    pdat, lg, pt = P.point(d, "or" if stat == "or" else "diff", rng_for(spec + "|point"), use_species)
    res = {"contrast_id": spec, "block": row["block"], "signal": s, "population": row["population"] + (" (S/T/Y only)" if posthoc else ""),
           "design": "D2", "index": row["index"], "tier_rule": row["tier_rule"], "decision_bearing": "no" if posthoc else row["decision_bearing"],
           "analysis": "post-hoc (labelled): same design restricted to S/T/Y records" if posthoc else "planned", "null": null,
           "n_analysis": int(len(d))}
    if pt["status"] != "fitted":
        res.update(status="not fitted", reason=pt["reason"],
                   note="Design cannot be fitted by the Rev 5 routine (phosnet statsmodels logit). In these records the non-canonical residues occur only among BRONZE records, so the residue-type term is perfectly separated." if "Singular" in pt["reason"] else "")
        json.dump(res, open(OUTD / f"{spec}.json", "w"), indent=1, default=float)
        print(spec, "not fitted", pt["reason"], flush=True)
        return
    crude_pt = P.crude(np.ones(pdat.n), pdat, stat if stat == "or" else "diff")
    rng = rng_for(spec)
    t1 = time.time()
    P.bootstrap(pdat, stat, rng_for(spec + "|timing"), 10)  # timing run on a separate stream, for the D7 projection
    per = (time.time() - t1) / 10
    n_total = N_BOOT if per * N_BOOT <= 7200 else 200
    rest = P.bootstrap(pdat, stat, rng, n_total)
    boot = {k: rest[k] for k in rest}
    est, cr, sm = boot["est"], boot["crude"], boot["smd"]
    a = att(cr, est, stat)
    res.update(status="run", estimate=pt["estimate"], ci=ci(est), crude=float(crude_pt), crude_ci=ci(cr),
               survives=bool(ci(est)[0] > null or ci(est)[1] < null), attenuation=float(att(np.array(crude_pt), np.array(pt["estimate"]), stat)),
               attenuation_ci=ci(a), smd=pt["smd"], balance_pass=bool(abs(pt["smd"]) < 0.1), n_pairs=pt["n_pairs"],
               n_treated_candidates=pt["n_treated"], match_rate=pt["n_pairs"] / pt["n_treated"],
               n_boot=n_total, n_boot_valid=int((~np.isnan(est)).sum()), n_boot_unfit=int(rest["unfit"]),
               boot_pairs_range=[float(np.nanmin(boot["pairs"])), float(np.nanmax(boot["pairs"]))],
               seconds_per_resample=round(per, 2), projected_1000_hours=round(per * 1000 / 3600, 2),
               resample_rule="1,000 resamples" if n_total == N_BOOT else "200 resamples (D7 fallback: projected 1,000-resample run exceeds 2 hours)",
               resampling_unit="protein (uniprot_id)", seed=f"[20260924, crc32({spec})]", seconds=round(time.time() - t0, 1))
    json.dump(res, open(OUTD / f"{spec}.json", "w"), indent=1, default=float)
    print(spec, res["estimate"], res["ci"], res["attenuation"], res["n_boot"], res["seconds_per_resample"], flush=True)


if __name__ == "__main__":
    for sp in sys.argv[1:]:
        run(sp)
