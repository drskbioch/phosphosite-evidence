"""Task E / block G5: protein abundance control. Population P-abund = records (GOLD or BRONZE, set-B tiers) of proteins with both the V50 index and a PaxDb
abundance (log10 ppm), per species. D0 crude, D1 CEM on index decile, D5 CEM on abundance decile alone, D4 CEM on index decile x abundance decile
(species and residue type exact in all). All designs use exactly the same records. Protein-cluster bootstrap, 1,000 resamples, whole pipeline re-run.
Also protein-level diagnostics: Spearman correlation of abundance with the index and with each outcome (protein mean), protein bootstrap CIs.
S1 (domain overlap) is defined for human records only, so mouse and rat rows for S1 are not run.
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

import run_cem_grid as R
from common import N_BOOT, PRIV, RES, rng_for
from engine import Cells, summarize
from populations import SIGNALS, population

b = R.b.merge(pd.read_parquet(PRIV / "abundance.parquet")[["uniprot_id", "abundance_ppm"]].assign(logab=lambda d: np.log10(d["abundance_ppm"])), on="uniprot_id", how="left")
b["tier"] = R.b["tier"].values
assert len(b) == len(R.b)
out = {}
m = R.man[R.man["block"] == "G5"]
MODE = {"D1": "x", "D5": "z", "D4": "xz"}
for (s, popl), g in m.groupby(["signal", "population"], sort=False):
    sp = popl.split("(")[1].rstrip(")")
    ids = dict(zip(g["design"], g["contrast_id"]))
    if s == "S1" and sp != "human":
        for cid in ids.values():
            out[cid] = {"status": "not run", "reason": "domain overlap (S1) is defined for human records only"}
        continue
    d = population(b, s, tier_col="tier", index_col="SEI_V50", species=sp, z_col="logab")
    c = Cells(d, "y", "x", {"protein": ["uniprot_id"]}, strata_cols=SIGNALS[s]["strata"], z="z")
    stat = SIGNALS[s]["stat"]
    null = 1.0 if stat == "or" else 0.0
    for design in ("D1", "D5", "D4"):
        r = summarize(c, stat, "protein", rng_for(ids[design]), n_boot=N_BOOT, mode=MODE[design])
        if not isinstance(r, tuple):
            out[ids[design]] = {"status": "failed", "reason": r.get("error")}
            continue
        r, bs = r
        att = r["attenuation"]
        base = {"signal": s, "population": popl, "index": "V50", "n_sites": r["n_sites"], "n_sites_gold": r["n_sites_gold"], "n_sites_bronze": r["n_sites_bronze"],
                "n_proteins": r["n_proteins"], "n_boot": N_BOOT, "n_boot_valid": r["n_boot_valid"], "resampling_unit": "protein (uniprot_id)",
                "seed": f"[20260924, crc32({ids[design]})]", "null": null}
        out[ids[design]] = dict(base, design=design, status="run", estimate=r["adjusted"], ci=r["adjusted_ci"], survives=r["survives"], crude=r["crude"], crude_ci=r["crude_ci"],
                                attenuation=att, attenuation_ci=r["attenuation_ci"], attenuation_label=R.label(att, r["crude_ci"], null, stat), smd_index=r["smd"], smd_abundance=r["smd_z"],
                                balance_index_pass=bool(abs(r["smd"]) < 0.1), balance_abundance_pass=bool(abs(r["smd_z"]) < 0.1),
                                n_strata_used=r["n_strata_used"], n_strata_dropped=r["n_strata_dropped"], n_bins=r["n_bins"])
        if design == "D1":  # D0 shares the D1 bootstrap
            out[ids["D0"]] = dict(base, design="D0", status="run", estimate=r["crude"], ci=r["crude_ci"], bootstrap_shared_with=ids["D1"])
        print(ids[design], s, sp, design, round(r["adjusted"], 4), [round(v, 4) for v in r["adjusted_ci"]], "att", round(att, 3), "smd_x", round(r["smd"], 3), "smd_z", round(r["smd_z"], 3), flush=True)
json.dump(out, open(RES / "g5_abundance.json", "w"), indent=1, default=float)

# ---- protein-level diagnostics
diag = {}
for sp in ("human", "mouse", "rat"):
    sub = b[(b["species"] == sp) & b["logab"].notna() & b["SEI_V50"].notna() & b["tier"].isin(["GOLD", "BRONZE"])]
    P = sub.drop_duplicates("uniprot_id")[["uniprot_id", "logab", "SEI_V50"]].set_index("uniprot_id")
    for s, sg in SIGNALS.items():
        if s == "S1" and sp != "human":
            continue
        y = b.loc[sub.index, sg["col"]].astype(float)
        P[s] = y.groupby(sub["uniprot_id"].values).mean()
    pairs = [("abundance_vs_index_V50", "logab", "SEI_V50")] + [(f"abundance_vs_{s}", "logab", s) for s in SIGNALS if s in P.columns]
    for name, a, c_ in pairs:
        dd = P[[a, c_]].dropna().reset_index(drop=True)
        rng = rng_for(f"G5diag|{sp}|{name}")
        est = stats.spearmanr(dd[a], dd[c_])[0]
        bs = [stats.spearmanr(*dd.iloc[rng.integers(0, len(dd), len(dd))].to_numpy().T)[0] for _ in range(N_BOOT)]
        diag[f"{sp}|{name}"] = {"spearman": float(est), "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], "n_proteins": int(len(dd))}
        print(sp, name, round(est, 3), flush=True)
json.dump(diag, open(RES / "g5_diagnostics.json", "w"), indent=1)
