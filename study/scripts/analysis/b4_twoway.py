"""Task B4: two-way (protein x study cluster) sensitivity for block G1, design D1.

Each record is assigned its primary study cluster (the cluster of its globally most frequent PMID). V_p, V_s and V_ps are
the bootstrap variances of the adjusted estimate when the resampling unit is the protein, the study cluster and the
protein x study-cluster cell (whole pipeline re-run in every resample, 1,000 resamples each), combined as
V = V_p + V_s - V_ps (Cameron-Gelbach-Miller form); CI = estimate +/- 1.96 sqrt(V) (log scale for odds ratios).
Criteria fixed before any result: not feasible if V <= 0, or if at least 90% of records lie in intersection cells that contain one record.
"""
import json

import numpy as np
import pandas as pd

import run_cem_grid as R
from common import N_BOOT, RES, rng_for
from engine import Cells, bootstrap, point
from populations import SIGNALS

res = {}
m = R.man[(R.man["block"] == "G1") & (R.man["design"] == "D1")]
for _, row in m.iterrows():
    s = row["signal"]
    stat = SIGNALS[s]["stat"]
    tcol, icol = R.tier_col_for(row["tier_rule"], s), R.IDX[row["index"]]
    d = R.population_for(row, tcol, icol)
    d["primary_cluster"] = d["primary_cluster"].astype(str)
    c = Cells(d, "y", "x", {"protein": ["uniprot_id"], "study": ["primary_cluster"], "cell": ["uniprot_id", "primary_cluster"]},
              strata_cols=SIGNALS[s]["strata"])
    p = point(c, stat)
    tr = (lambda v: np.log(v)) if stat == "or" else (lambda v: v)
    est = tr(p["adj"])
    V = {}
    for scheme in ("protein", "study", "cell"):
        bs = bootstrap(c, stat, scheme, rng_for(f"B4|{row['contrast_id']}|{scheme}"), N_BOOT)
        a = tr(bs["adj"])
        V[scheme] = float(np.nanvar(a, ddof=1))
    cell_sizes = d.groupby(["uniprot_id", "primary_cluster"]).size()
    single_share = float((cell_sizes == 1).sum() / len(d))
    Vtot = V["protein"] + V["study"] - V["cell"]
    r = {"contrast_id": row["contrast_id"], "signal": s, "scale": "log odds ratio" if stat == "or" else "difference",
         "adjusted": float(p["adj"]), "V_protein": V["protein"], "V_study": V["study"], "V_cell": V["cell"], "V_two_way": float(Vtot),
         "n_records": int(len(d)), "n_proteins": int(d["uniprot_id"].nunique()), "n_study_clusters": int(d["primary_cluster"].nunique()),
         "n_intersection_cells": int(len(cell_sizes)), "share_of_records_in_single_record_cells": single_share,
         "se_protein": float(np.sqrt(V["protein"])), "se_study": float(np.sqrt(V["study"]))}
    if Vtot <= 0:
        r.update(feasible=False, reason="V <= 0")
    elif single_share >= 0.9:
        r.update(feasible=False, reason="intersection cells are essentially singletons (>= 90% of records)")
    else:
        se = float(np.sqrt(Vtot))
        lo, hi = est - 1.96 * se, est + 1.96 * se
        ci = [float(np.exp(lo)), float(np.exp(hi))] if stat == "or" else [float(lo), float(hi)]
        null = 1.0 if stat == "or" else 0.0
        pci = [float(est - 1.96 * np.sqrt(V["protein"])), float(est + 1.96 * np.sqrt(V["protein"]))]
        pci = [float(np.exp(v)) for v in pci] if stat == "or" else pci
        r.update(feasible=True, ci_two_way=ci, ci_protein_normal_approx=pci, se_two_way=se, survives_two_way=bool(ci[0] > null or ci[1] < null),
                 width_ratio_two_way_vs_protein=float((ci[1] - ci[0]) / (pci[1] - pci[0])))
    res[row["contrast_id"]] = r
    print(json.dumps(r, default=float), flush=True)
json.dump(res, open(RES / "b4_twoway.json", "w"), indent=1, default=float)
