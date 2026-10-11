"""CEM runs of the blocks added by the plan amendment of 2026-10-01: G16 (canonical accessions), G17 (gene-cluster bootstrap), G19 (tier definitions).
Gene-cluster bootstrap (resampling unit = gene_cluster, see gene_unit.py), 1,000 resamples, whole pipeline re-run in each resample, D0 and D1 from the same
resamples, generator np.random.default_rng([20260924, crc32(D1 contrast id)]).
Usage: python run_amend.py G17 [G16] [G19]      (results merged into results/amend_results.json)"""
import json
import sys

import numpy as np
import pandas as pd

import run_cem_grid as R
from common import PRIV, RES, N_BOOT, rng_for
from engine import Cells, summarize
from populations import SIGNALS, population

OUT = RES / "amend_results.json"
b = R.b
b["tier_agree"] = np.where(b["DAS"] >= 2, "GOLD", np.where(b["DAS"] == 0, "BRONZE", "SILVER"))
b.loc[b["DAS"].isna(), "tier_agree"] = None
TIERCOL = {"set-B tiers (primary)": "tier", "set-A tiers": "tier_setA", "depth-only tiers": "tier_depth", "agreement-only (agreement >= 2 vs 0)": "tier_agree"}


def canonical_accessions():
    """Swiss-Prot reviewed accessions (UniProt 2026_03; UniProt lists one canonical accession per entry, isoforms have no accession of their own).
    Where several reviewed accessions belong to one gene cluster (65 clusters), the one with the most records in the site table is kept
    (ties: alphabetical order of the accession); rule fixed before the run."""
    r = pd.read_parquet(PRIV / "swissprot_reviewed.parquet")
    cnt = b["uniprot_id"].value_counts()
    gc = b[["uniprot_id", "gene_cluster"]].drop_duplicates("uniprot_id").set_index("uniprot_id")["gene_cluster"]
    x = pd.DataFrame({"acc": sorted(set(r["Entry"]) & set(cnt.index))})
    x["gc"] = x["acc"].map(gc)
    x["n"] = x["acc"].map(cnt)
    mapped = ~x["gc"].str.startswith("A")
    keep = set(x.loc[~mapped, "acc"])
    y = x[mapped].sort_values(["gc", "n", "acc"], ascending=[True, False, True]).drop_duplicates("gc")
    return keep | set(y["acc"]), {"reviewed_accessions_in_site_table": int(len(x)), "gene_clusters_with_more_than_one_reviewed_accession": int((x[mapped].groupby("gc").size() > 1).sum()),
                                  "accessions_kept": int(len(keep) + len(y))}


def canonical_mask():
    acc, info = canonical_accessions()
    return b["uniprot_id"].isin(acc).to_numpy()


def run_g5_gene(g):
    ab = pd.read_parquet(PRIV / "abundance.parquet")[["uniprot_id", "abundance_ppm"]].assign(logab=lambda d: np.log10(d["abundance_ppm"]))
    bb = b.merge(ab, on="uniprot_id", how="left")
    d0 = g[g["design"] == "D0"].iloc[0]
    d1 = g[g["design"] == "D1"].iloc[0]
    s = d1["signal"]
    d = population(bb, s, tier_col="tier", index_col="SEI_V50", species="human", z_col="logab")
    c = Cells(d, "y", "x", {"gene": ["gene_cluster"]}, strata_cols=SIGNALS[s]["strata"], z="z")
    r, bs = summarize(c, SIGNALS[s]["stat"], "gene", rng_for(d1["contrast_id"]), n_boot=N_BOOT, mode="x")
    return r, d1, d0, d


def pack(r, d1, d0, d, s, extra):
    null = 1.0 if SIGNALS[s]["stat"] == "or" else 0.0
    common = {"signal": s, "population": d1["population"], "index": d1["index"], "tier_rule": d1["tier_rule"], "variation": d1["variation"],
              "decision_bearing": d1["decision_bearing"], "null": null, "n_sites": r["n_sites"], "n_sites_gold": r["n_sites_gold"],
              "n_sites_bronze": r["n_sites_bronze"], "n_proteins": r["n_proteins"], "n_gene_clusters": int(d["gene_cluster"].nunique()),
              "n_boot": r["n_boot"], "n_boot_valid": r["n_boot_valid"], "resampling_unit": "gene cluster (gene_cluster)",
              "seed": f"[20260924, crc32({d1['contrast_id']})]", **extra}
    att = r["attenuation"]
    return {
        d0["contrast_id"]: dict(common, design="D0", status="run", estimate=r["crude"], ci=r["crude_ci"],
                                ci_excludes_null=bool((r["crude_ci"][0] > null) or (r["crude_ci"][1] < null)), bootstrap_shared_with=d1["contrast_id"], smd_pre=r["smd_pre"]),
        d1["contrast_id"]: dict(common, design="D1", status="run", estimate=r["adjusted"], ci=r["adjusted_ci"], survives=r["survives"], crude=r["crude"],
                                crude_ci=r["crude_ci"], attenuation=att, attenuation_ci=r["attenuation_ci"],
                                attenuation_label=R.label(att, r["crude_ci"], null, SIGNALS[s]["stat"]), largely_explained=r["largely_explained"],
                                smd=r["smd"], smd_pre=r["smd_pre"], balance_pass=r["balance_pass"], n_strata_used=r["n_strata_used"],
                                n_strata_dropped=r["n_strata_dropped"], n_bins=r["n_bins"], adjusted_boot_sd=r["adjusted_boot_sd"]),
    }


def run(block, g, store):
    d1 = g[g["design"] == "D1"].iloc[0]
    d0 = g[g["design"] == "D0"].iloc[0]
    s, popl = d1["signal"], d1["population"]
    extra = {}
    if block == "G17" and popl.startswith("P-abund"):
        r, d1, d0, d = run_g5_gene(g)
    else:
        if block == "G16":
            mask = canonical_mask()
            row = d1.copy()
            row["population"] = "P-human" if "P-human" in popl else "P-all"
            res = R.run_group(g.assign(population=row["population"]), scheme="gene", mask=mask)
            for k, v in res.items():
                v["population"] = popl
            return res
        if block == "G19":
            tcol = TIERCOL[d1["tier_rule"]]
            res = R.run_group(g, scheme="gene", tier_col=tcol)
            return res
        res = R.run_group(g, scheme="gene")
        return res
    return pack(r, d1, d0, d, s, extra)


if __name__ == "__main__":
    blocks = sys.argv[1:]
    if "G16" in blocks:
        acc, info = canonical_accessions()
        for sp in ("human", "mouse", "rat"):
            x = b[b["species"] == sp]
            kept = x["uniprot_id"].isin(acc)
            info[sp] = {"records": int(len(x)), "records_kept": int(kept.sum()), "share_records_kept": float(kept.mean()), "accessions": int(x["uniprot_id"].nunique()),
                        "accessions_kept": int(x.loc[kept, "uniprot_id"].nunique())}
        info["source"] = "UniProt REST stream, reviewed:true, organism 9606/10090/10116, release 2026_03 (02-September-2026), retrieved 2026-10-01"
        json.dump(info, open(RES / "g16_coverage.json", "w"), indent=1)
    store = json.load(open(OUT)) if OUT.exists() else {}
    m = R.man[R.man["block"].isin(blocks) & R.man["design"].isin(["D0", "D1"])]
    keys = ["block", "signal", "population", "index", "tier_rule", "variation"]
    for _, g in m.groupby(keys, sort=False):
        blk = g["block"].iloc[0]
        out = run(blk, g, store)
        store.update(out)
        for cid, r in out.items():
            print(cid, r.get("design"), r.get("status"), r.get("estimate"), r.get("ci"), r.get("attenuation"), flush=True)
        # reload before writing so that parallel runs of other blocks are not overwritten
        cur = json.load(open(OUT)) if OUT.exists() else {}
        cur.update(out)
        json.dump(cur, open(OUT, "w"), indent=1, default=float)
