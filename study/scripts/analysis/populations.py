"""Analysis populations for each signal (plan section 3)."""
import numpy as np
import pandas as pd

SIGNALS = {
    "S1": {"col": "in_domain", "stat": "or", "label": "domain overlap", "strata": ("residue_type",), "human_only": True},
    "S2": {"col": "conserved", "stat": "diff", "label": "conservation", "strata": ("species", "residue_type")},
    "S3": {"col": "disorder", "stat": "diff", "label": "disorder", "strata": ("species", "residue_type")},
    "S4": {"col": "seq_motif_match", "stat": "diff", "label": "kinase motif (sequence-based)", "strata": ("species", "residue_type")},
    "S5": {"col": "kinase_match", "stat": "diff", "label": "kinase motif (database lookup)", "strata": ("species", "residue_type")},
}


def population(base, signal, tier_col="tier", index_col="SEI_V50", species=None, residues=None,
               mask=None, require_index=True, z_col=None):
    """GOLD and BRONZE records with the outcome and the index defined.
    tier_col may be a column name of `base`. mask: optional boolean array over base rows."""
    s = SIGNALS[signal]
    d = base
    sel = d[tier_col].isin(["GOLD", "BRONZE"]).to_numpy()
    if s.get("human_only"):
        sel &= (d["species"] == "human").to_numpy()
    if species is not None:
        sel &= (d["species"] == species).to_numpy()
    if residues is not None:
        sel &= d["residue_type"].isin(residues).to_numpy()
    if mask is not None:
        sel &= mask
    sel &= d[s["col"]].notna().to_numpy()
    if require_index:
        sel &= d[index_col].notna().to_numpy()
    if z_col is not None:
        sel &= d[z_col].notna().to_numpy()
    cols = ["uniprot_id"] + (["gene_cluster"] if "gene_cluster" in d.columns else []) + ["species", "residue_type", "primary_cluster", s["col"], index_col, tier_col] + ([z_col] if z_col else [])
    out = d.loc[sel, cols].copy()
    out = out.rename(columns={s["col"]: "y", index_col: "x", tier_col: "tier", **({z_col: "z"} if z_col else {})})
    out["y"] = out["y"].astype(float)
    return out
