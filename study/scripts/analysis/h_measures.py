"""Task H / block G9: established confidence measures as contrasts.

For every measure the treated group is the high-confidence class and the control group the low-confidence class:
  psp_lt   : PhosphoSitePlus LT_LIT >= 1 (treated) vs LT_LIT = 0 with MS_LIT + MS_CST >= 1 (control)
  psp_cat  : n = LT_LIT + MS_LIT + MS_CST; high (n >= 5) vs low (n = 1); medium (n = 2 to 4) is not compared (Kalyuzhnyy et al. 2022 categories)
  uniprot  : reviewed UniProt MOD_RES phosphosites; experimental (ECO:0000269 among the evidence tags) vs inferred (only ECO:0000250, ECO:0000255, ECO:0000305)
Records: human records of the site table matching a site of the measure on accession, residue and position. Outcomes S1 to S4 (outcomes of the site table).
D0 crude, D1 CEM on index decile (V50), D4 CEM on index decile x abundance decile (records with a PaxDb abundance; D1 and D0 on those records are reported beside it).
Gene-cluster bootstrap (primary, seed crc32(<id>|gene)) and accession-cluster bootstrap (seed crc32(<id>)), 1,000 resamples, whole pipeline re-run in each resample.

Usage: python h_measures.py psp_lt psp_cat uniprot     (PhosphoSitePlus measures need data_private/Phosphorylation_site_dataset.gz)
Output: results/g9_measures.json, results/g9_intersections.json
"""
import gzip
import io
import json
import re
import sys

import numpy as np
import pandas as pd

import run_cem_grid as R
from common import N_BOOT, PRIV, RES, rng_for
from engine import Cells, summarize
from populations import SIGNALS, population

PSP_FILE = PRIV / "Phosphorylation_site_dataset.gz"
MEASURE_NAMES = {"psp_lt": "PSP LT-supported vs MS-only", "psp_cat": "PSP evidence category high (n>=5) vs low (n=1)",
                 "uniprot": "UniProt MOD_RES experimental vs inferred ECO"}
REQUIRED = ["ACC_ID", "MOD_RSD", "ORGANISM", "LT_LIT", "MS_LIT", "MS_CST"]


def load_psp():
    """Phosphorylation_site_dataset.gz -> one row per (accession, residue, position, organism). Stops (raises) if the columns differ from the plan's definition (D3)."""
    raw = gzip.open(PSP_FILE, "rt", encoding="utf-8", errors="replace").read().splitlines()
    hdr = next(i for i, l in enumerate(raw) if l.startswith("GENE\t"))
    d = pd.read_csv(io.StringIO("\n".join(raw[hdr:])), sep="\t", dtype=str)
    missing = [c for c in REQUIRED if c not in d.columns]
    if missing:
        raise SystemExit(f"PhosphoSitePlus file: columns {missing} not found (found {list(d.columns)}); the column semantics of D3 cannot be applied, stopping")
    info = {"file": PSP_FILE.name, "header_lines_before_columns": hdr, "columns": list(d.columns), "first_lines": raw[:min(hdr, 4)], "rows": int(len(d))}
    d = d[d["MOD_RSD"].str.endswith("-p", na=False)].copy()
    info["phospho_rows"] = int(len(d))
    m = d["MOD_RSD"].str.extract(r"^([A-Z])(\d+)-p$")
    d["residue_type"] = m[0]
    d["position"] = pd.to_numeric(m[1], errors="raise")
    for c in ("LT_LIT", "MS_LIT", "MS_CST"):
        info[f"{c}_missing_set_to_0"] = int(d[c].isna().sum())
        d[c] = pd.to_numeric(d[c], errors="raise").fillna(0)
    d["n"] = d["LT_LIT"] + d["MS_LIT"] + d["MS_CST"]
    d["species"] = d["ORGANISM"].str.lower()
    iso = d["ACC_ID"].str.contains("-", na=False)
    info["isoform_specific_accessions_excluded"] = int(iso.sum())
    d = d[~iso].rename(columns={"ACC_ID": "uniprot_id"})
    dup = d.duplicated(["uniprot_id", "residue_type", "position", "species"], keep=False)
    info["duplicate_site_rows"] = int(dup.sum())
    d = d.sort_values(["uniprot_id", "residue_type", "position", "species", "n"], ascending=[True, True, True, True, False]).drop_duplicates(["uniprot_id", "residue_type", "position", "species"])
    info["sites_by_organism"] = d["species"].value_counts().to_dict()
    return d, info


def psp_measure(kind):
    d, info = load_psp()
    d = d[d["species"] == "human"]
    if kind == "psp_lt":
        g = np.where(d["LT_LIT"] >= 1, 1, np.where((d["LT_LIT"] == 0) & (d["MS_LIT"] + d["MS_CST"] >= 1), 0, -1))
    else:
        g = np.where(d["n"] >= 5, 1, np.where(d["n"] == 1, 0, -1))
    d = d.assign(grp=g)
    info["human_psp_sites"] = int(len(d))
    info["human_sites_high_class"] = int((g == 1).sum())
    info["human_sites_low_class"] = int((g == 0).sum())
    info["human_sites_not_compared"] = int((g == -1).sum())
    return d[d["grp"] >= 0][["uniprot_id", "residue_type", "position", "grp"]], info, int(len(d))


def uniprot_measure():
    d = pd.read_csv(PRIV / "dl" / "seq_9606.tsv", sep="\t")
    rows = []
    for acc, ft in zip(d["Entry"], d["Modified residue"].fillna("")):
        for m in re.finditer(r'MOD_RES (\d+); /note="(Phospho[^"]*)"(?:; /evidence="([^"]*)")?', ft):
            ecos = set(re.findall(r"ECO:\d{7}", m.group(3) or ""))
            rows.append((acc, int(m.group(1)), m.group(2), ecos))
    f = pd.DataFrame(rows, columns=["uniprot_id", "position", "note", "ecos"])
    res = {"Phosphoserine": "S", "Phosphothreonine": "T", "Phosphotyrosine": "Y"}
    f["residue_type"] = f["note"].map(lambda s: next((v for k, v in res.items() if s.startswith(k)), None))
    f = f[f["residue_type"].notna()]
    inferred_codes = {"ECO:0000250", "ECO:0000255", "ECO:0000305"}
    f["grp"] = f["ecos"].map(lambda e: 1 if "ECO:0000269" in e else (0 if e and e <= inferred_codes else -1))
    info = {"source": "UniProt REST stream, reviewed human entries, field ft_mod_res, release 2026_03 (02-September-2026), retrieved 2026-10-01",
            "phospho_features": int(len(f)), "experimental_ECO0000269": int((f["grp"] == 1).sum()), "inferred_only": int((f["grp"] == 0).sum()),
            "neither (only ECO:0007744 or mixed or no tag)": int((f["grp"] == -1).sum())}
    f = f.drop_duplicates(["uniprot_id", "residue_type", "position"])
    return f[f["grp"] >= 0][["uniprot_id", "residue_type", "position", "grp"]], info, int(len(f))


def summarize_pair(c, stat, ids, mode):
    """gene-cluster (primary) and accession-cluster bootstrap on the same cells."""
    out = {}
    for scheme, seed in (("gene", ids + "|gene"), ("protein", ids)):
        r = summarize(c, stat, scheme, rng_for(seed), n_boot=N_BOOT, mode=mode)
        out[scheme] = r[0] if isinstance(r, tuple) else r
    return out


def pack(r, s):
    null = 1.0 if SIGNALS[s]["stat"] == "or" else 0.0
    keys = ("adjusted", "adjusted_ci", "crude", "crude_ci", "attenuation", "attenuation_ci", "survives", "largely_explained", "smd", "smd_pre", "smd_z", "n_strata_used", "n_strata_dropped", "n_bins")
    o = {k: r.get(k) for k in keys}
    o["null"] = null
    o["crude_excludes_null"] = bool(r["crude_ci"][0] > null or r["crude_ci"][1] < null)
    o["balance_pass"] = r["balance_pass"]
    return o


def run_measure(kind, meas, info, n_all):
    name = MEASURE_NAMES[kind]
    b = R.b.merge(meas.assign(species="human"), on=["uniprot_id", "residue_type", "position", "species"], how="left")
    assert len(b) == len(R.b)
    b["mgrp"] = np.where(b["grp"] == 1, "GOLD", np.where(b["grp"] == 0, "BRONZE", None))
    ab = pd.read_parquet(PRIV / "abundance.parquet")[["uniprot_id", "abundance_ppm"]].assign(logab=lambda d: np.log10(d["abundance_ppm"]))
    b = b.merge(ab, on="uniprot_id", how="left")
    m = R.man[(R.man["block"] == "G9") & (R.man["variation"] == name)]
    inter = b[b["mgrp"].notna()]
    meta = dict(info, measure=name, site_records_matched=int(len(inter)), share_of_measure_sites_lost=float(1 - len(inter) / max(len(meas), 1)) if len(meas) else None,
                n_measure_sites_compared=int(len(meas)), n_measure_sites_before_class_filter=int(n_all), high_class_records=int((inter["mgrp"] == "GOLD").sum()),
                low_class_records=int((inter["mgrp"] == "BRONZE").sum()), accessions=int(inter["uniprot_id"].nunique()), gene_clusters=int(inter["gene_cluster"].nunique()))
    out = {}
    for s in ("S1", "S2", "S3", "S4"):
        ids = dict(zip(m[m["signal"] == s]["design"], m[m["signal"] == s]["contrast_id"]))
        stat = SIGNALS[s]["stat"]
        ucols = {"gene": ["gene_cluster"], "protein": ["uniprot_id"]}
        d = population(b, s, tier_col="mgrp", index_col="SEI_V50", species="human")
        base = {"signal": s, "measure": name, "population": "human records matched to the measure", "n_records": int(len(d)), "n_high": int((d["tier"] == "GOLD").sum()),
                "n_low": int((d["tier"] == "BRONZE").sum()), "n_accessions": int(d["uniprot_id"].nunique()), "n_gene_clusters": int(d["gene_cluster"].nunique())}
        if base["n_high"] < 20 or base["n_low"] < 20:
            for dd, cid in ids.items():
                out[cid] = dict(base, design=dd, status="not run", reason="fewer than 20 records in one class")
            continue
        c = Cells(d, "y", "x", ucols, strata_cols=SIGNALS[s]["strata"])
        r = summarize_pair(c, stat, ids["D1"], "x")
        out[ids["D0"]] = dict(base, design="D0", status="run", gene=dict(estimate=r["gene"]["crude"], ci=r["gene"]["crude_ci"]), protein=dict(estimate=r["protein"]["crude"], ci=r["protein"]["crude_ci"]))
        out[ids["D1"]] = dict(base, design="D1", status="run", gene=pack(r["gene"], s), protein=pack(r["protein"], s))
        dz = population(b, s, tier_col="mgrp", index_col="SEI_V50", species="human", z_col="logab")
        bz = dict(base, n_records=int(len(dz)), n_high=int((dz["tier"] == "GOLD").sum()), n_low=int((dz["tier"] == "BRONZE").sum()), n_accessions=int(dz["uniprot_id"].nunique()),
                  n_gene_clusters=int(dz["gene_cluster"].nunique()), population="the same, restricted to accessions with a PaxDb abundance")
        if bz["n_high"] < 20 or bz["n_low"] < 20:
            out[ids["D4"]] = dict(bz, design="D4", status="not run", reason="fewer than 20 records in one class with an abundance")
            continue
        cz = Cells(dz, "y", "x", ucols, strata_cols=SIGNALS[s]["strata"], z="z")
        r4 = summarize_pair(cz, stat, ids["D4"], "xz")
        r1 = summarize_pair(cz, stat, ids["D4"] + "|D1_same_records", "x")
        out[ids["D4"]] = dict(bz, design="D4", status="run", gene=pack(r4["gene"], s), protein=pack(r4["protein"], s),
                              same_records_D1=dict(gene=pack(r1["gene"], s), protein=pack(r1["protein"], s)))
        print(kind, s, "D1", round(r["gene"]["adjusted"], 4), [round(v, 4) for v in r["gene"]["adjusted_ci"]], "D4", round(r4["gene"]["adjusted"], 4), flush=True)
    return meta, out


if __name__ == "__main__":
    kinds = sys.argv[1:]
    store = json.load(open(RES / "g9_measures.json")) if (RES / "g9_measures.json").exists() else {}
    inter = json.load(open(RES / "g9_intersections.json")) if (RES / "g9_intersections.json").exists() else {}
    for k in kinds:
        meas, info, n_all = uniprot_measure() if k == "uniprot" else psp_measure(k)
        meta, out = run_measure(k, meas, info, n_all)
        store.update(out)
        inter[k] = meta
        json.dump(store, open(RES / "g9_measures.json", "w"), indent=1, default=float)
        json.dump(inter, open(RES / "g9_intersections.json", "w"), indent=1, default=float)
