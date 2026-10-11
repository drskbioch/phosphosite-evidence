"""Check the engine against phosnet's matched_comparison for the Rev 6 primary index (V50), on the Rev 6 populations (n_boot=2: point estimate only)."""
import json, sys
sys.path.insert(0, "../../phosnet/src")
from phosnet.matching import matched_comparison
import base as B
from common import RES
from engine import Cells, point
from populations import SIGNALS, population
b = B.load("lean")
rows = {}
for s in ("S1", "S2", "S3", "S4", "S5"):
    d = population(b, s, tier_col="tier", index_col="SEI_V50")
    c = Cells(d, "y", "x", {"protein": ["uniprot_id"]}, strata_cols=SIGNALS[s]["strata"])
    p = point(c, SIGNALS[s]["stat"])
    ph = matched_comparison(d.rename(columns={"y": "yy", "x": "xx"}), "tier", "GOLD", "yy", "fisher_exact" if s == "S1" else "proportions_ztest", "xx",
                            extra_strata_cols=SIGNALS[s]["strata"], n_boot=2)
    rows[s] = {"engine": p["adj"], "phosnet": ph["point_estimate"], "abs_diff": abs(p["adj"] - ph["point_estimate"]), "smd_engine": p["smd"], "smd_phosnet": ph["smd"]}
    print(s, rows[s], flush=True)
json.dump(rows, open(RES / "validate_engine_v50.json", "w"), indent=1, default=float)
