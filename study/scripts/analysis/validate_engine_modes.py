"""Check engine modes "z" (abundance deciles only) and "xz" (index decile x abundance decile) against phosnet.matched_comparison on human records
with abundance (point estimates; phosnet given the abundance as the covariate (z) or an abundance-decile column as an extra stratum (xz))."""
import json, sys
sys.path.insert(0, "../../phosnet/src")
import numpy as np, pandas as pd
from phosnet.matching import matched_comparison
from phosnet.sei import sei_deciles
import base as B
from common import PRIV, RES
from engine import Cells, point
from populations import SIGNALS, population
b = B.load("lean")
ab = pd.read_parquet(PRIV / "abundance.parquet")
b = b.merge(ab[["uniprot_id", "abundance_ppm"]].assign(logab=lambda d: np.log10(d["abundance_ppm"])), on="uniprot_id", how="left")
rows = {}
for s in ("S2", "S3"):
    d = population(b, s, tier_col="tier", index_col="SEI_V50", species="human", z_col="logab")
    stat = "proportions_ztest"
    c = Cells(d, "y", "x", {"protein": ["uniprot_id"]}, strata_cols=SIGNALS[s]["strata"], z="z")
    for mode in ("z", "xz"):
        p = point(c, "diff", mode)
        dd = d.rename(columns={"y": "yy", "x": "xx", "z": "zz"}).copy()
        if mode == "z":
            ph = matched_comparison(dd, "tier", "GOLD", "yy", stat, "zz", extra_strata_cols=SIGNALS[s]["strata"], n_boot=2)
        else:
            dec, _ = sei_deciles(dd["zz"])
            dd["zdec"] = dec.to_numpy()
            ph = matched_comparison(dd, "tier", "GOLD", "yy", stat, "xx", extra_strata_cols=tuple(SIGNALS[s]["strata"]) + ("zdec",), n_boot=2)
        rows[f"{s}_{mode}"] = {"engine": p["adj"], "phosnet": ph["point_estimate"], "abs_diff": abs(p["adj"] - ph["point_estimate"]),
                              "smd_engine": p["smd"] if mode == "xz" else p["smd_z"], "smd_phosnet": ph["smd"]}
        print(s, mode, rows[f"{s}_{mode}"], flush=True)
json.dump(rows, open(RES / "validate_engine_modes.json", "w"), indent=1, default=float)
