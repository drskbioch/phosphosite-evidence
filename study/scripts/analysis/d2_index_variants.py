"""Task D 1-2: protein-level index variants.

V0   original primary index (reproduced and checked against sei_protein.parquet)
V50  PMIDs linked to > 50 genes (over all organisms in gene2pubmed) excluded  (primary in Rev 6)
V20, V100  sensitivities
VT   log10(1 + n50 / t), t = max(1, ref_year - y0 + 1); y0 = earliest year among the 5 lowest-numbered
     PMIDs of the protein (bulk-corrected set) whose publication year could be resolved.
     ref_year 2025 (primary) and 2026 (sensitivity).

Publication years for the lowest-numbered PMIDs come from NCBI E-utilities esummary (PubDate leading
4-digit year), fetched here and cached in study/analysis/data_private/pmid_year_gene.parquet.
"""
import json
import re
import sys
import time

import numpy as np
import pandas as pd
import requests
from scipy import stats

from common import PRIV, REL, RES

K_LOW = None  # None: y0 from all of a protein's bulk-corrected PMIDs (plan wording). A 5-lowest-PMID proxy was tried first and missed earlier years for >= 4.5% of proteins.
THR = {"V20": 20, "V50": 50, "V100": 100}
URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

sei = pd.read_parquet(REL / "sei_protein.parquet")
pairs = pd.read_parquet(PRIV / "g2p_pairs_3sp.parquet")[["gene", "pmid"]].drop_duplicates()
gpp = pd.read_parquet(PRIV / "g2p_genes_per_pmid.parquet")
gpp["n_genes"] = pd.to_numeric(gpp["n_genes"], errors="raise")

pg = sei.loc[sei["mapped_to_gene"] == True, ["uniprot_id", "gene_ids"]].explode("gene_ids")  # noqa: E712
pg = pg.dropna(subset=["gene_ids"]).rename(columns={"gene_ids": "gene"})
pg["gene"] = pg["gene"].astype(np.int64)
pp = pg.merge(pairs, on="gene").drop_duplicates(["uniprot_id", "pmid"]).merge(gpp, on="pmid", how="left")
print("protein-PMID links", len(pp), "missing gene counts", int(pp["n_genes"].isna().sum()))

out = sei[["uniprot_id", "species", "mapped_to_gene", "SEI_primary", "SEI_secondary", "n_pubmed"]].copy()
n0 = pp.groupby("uniprot_id").size().rename("n_V0")
out = out.merge(n0, on="uniprot_id", how="left")
mapped = out["mapped_to_gene"] == True  # noqa: E712
out.loc[mapped, "n_V0"] = out.loc[mapped, "n_V0"].fillna(0)
chk = out.loc[mapped]
summary = {"reproduces_sei_primary_n_pubmed": bool((chk["n_V0"] == chk["n_pubmed"]).all()),
           "n_mapped_proteins": int(mapped.sum()), "n_protein_pmid_links": int(len(pp)),
           "n_distinct_pmids_3species": int(pp["pmid"].nunique())}
print(summary)

for name, t in THR.items():
    keep = pp[pp["n_genes"] <= t]
    n = keep.groupby("uniprot_id").size().rename(f"n_{name}")
    out = out.merge(n, on="uniprot_id", how="left")
    out.loc[mapped, f"n_{name}"] = out.loc[mapped, f"n_{name}"].fillna(0)
    out[f"SEI_{name}"] = np.log10(1 + out[f"n_{name}"])
    removed_pm = pp.loc[pp["n_genes"] > t, "pmid"].nunique()
    summary[name] = {"threshold_genes_per_pmid": t, "pmids_removed": int(removed_pm),
                     "protein_pmid_links_removed": int((pp["n_genes"] > t).sum()),
                     "proteins_with_changed_n": int((out.loc[mapped, f"n_{name}"] != out.loc[mapped, "n_V0"]).sum()),
                     "proteins_with_zero_after": int((out.loc[mapped, f"n_{name}"] == 0).sum())}
out["SEI_V0"] = out["SEI_primary"]

# ---- publication years of the lowest-numbered PMIDs (bulk-corrected set)
p50 = pp[pp["n_genes"] <= 50][["uniprot_id", "pmid"]].sort_values(["uniprot_id", "pmid"])
low = p50 if K_LOW is None else p50.groupby("uniprot_id").head(K_LOW)
need = set(low["pmid"].unique())
have = pd.read_parquet(REL / "pmid_year.parquet")[["pmid", "year"]]
have["pmid"] = pd.to_numeric(have["pmid"], errors="raise").astype(np.int64)
have["year"] = pd.to_numeric(have["year"], errors="raise")
cache_path = PRIV / "pmid_year_gene.parquet"
cached = pd.read_parquet(cache_path) if cache_path.exists() else pd.DataFrame(columns=["pmid", "year"])
known = set(have["pmid"]) | set(cached["pmid"])
todo = sorted(need - known)
print("PMIDs needing a year", len(need), "to fetch", len(todo), flush=True)
rows = []
failed = 0
t0 = time.time()
B = 500
WORKERS = 8  # NCBI allows 3 requests/s without a key; each request takes several seconds, so 8 in flight stays below it
from concurrent.futures import ThreadPoolExecutor


def fetch_batch(batch):
    for attempt in range(6):
        try:
            r = requests.post(URL, data={"db": "pubmed", "id": ",".join(map(str, batch)), "retmode": "json"}, timeout=180)
            r.raise_for_status()
            res = r.json().get("result", {})
            out_rows = []
            for p in batch:
                rec = res.get(str(p))
                if rec:
                    m = re.match(r"(\d{4})", (rec.get("pubdate") or rec.get("sortpubdate") or "").strip())
                    if m:
                        out_rows.append((p, int(m.group(1))))
            return out_rows, 0
        except Exception:  # noqa: BLE001
            time.sleep(5 * (attempt + 1))
    return [], len(batch)


batches = [todo[i:i + B] for i in range(0, len(todo), B)]
with ThreadPoolExecutor(WORKERS) as ex:
    for k, (r_, f_) in enumerate(ex.map(fetch_batch, batches)):
        rows.extend(r_)
        failed += f_
        if k % 40 == 0:
            print(k, len(batches), round(time.time() - t0), flush=True)
            pd.concat([cached, pd.DataFrame(rows, columns=["pmid", "year"])]).drop_duplicates("pmid").to_parquet(cache_path)
allyears = pd.concat([have, cached, pd.DataFrame(rows, columns=["pmid", "year"])]).drop_duplicates("pmid")
allyears.to_parquet(cache_path)
summary["years"] = {"pmids_needed": len(need), "fetched": len(rows), "fetch_failed_batches_pmids": failed,
                    "pmids_still_without_year": int(len(need - set(allyears["pmid"])))}
low = low.merge(allyears, on="pmid", how="left")
y0 = low.groupby("uniprot_id")["year"].min().rename("y0")
out = out.merge(y0, on="uniprot_id", how="left")
# explicit numeric casts: bad values fail loudly
for c in ("n_V0", "n_V20", "n_V50", "n_V100"):
    out[c] = pd.to_numeric(out[c], errors="raise")
out["y0"] = pd.to_numeric(out["y0"], errors="raise")
allyears["year"] = pd.to_numeric(allyears["year"], errors="raise")
has_y0 = out["y0"].notna()

# VT = log10(1 + n_V50 / t), t = ref_year - y0 + 1 (t >= 1 when y0 <= ref_year).
# Proteins without a dated PMID are NOT filled or imputed: VT stays missing and they are excluded only
# from analyses that use VT. y0 later than the reference year gives t <= 0 (undefined): also left missing, counted.
summary["VT_rows"] = {}
for ref in (2025, 2026):
    t = pd.to_numeric(ref - out["y0"] + 1, errors="raise")
    after = has_y0 & (t < 1)
    ok = has_y0 & (t >= 1)
    assert t[ok].notna().all() and (t[ok] >= 1).all() and np.isfinite(t[ok]).all(), "t must be >= 1 and finite in the rows used"
    assert out.loc[ok, "n_V50"].notna().all()
    out[f"SEI_VT{ref}"] = np.nan
    out.loc[ok, f"SEI_VT{ref}"] = np.log10(1 + out.loc[ok, "n_V50"].astype(float) / t[ok].astype(float))
    assert out.loc[ok, f"SEI_VT{ref}"].notna().all()
    summary["VT_rows"][str(ref)] = {
        "proteins_mapped": int(mapped.sum()), "proteins_with_VT": int((mapped & ok).sum()),
        "proteins_missing_no_dated_pmid": int((mapped & ~has_y0).sum()),
        "proteins_missing_first_year_after_reference": int((mapped & after).sum()),
        "t_min_used": float(t[ok].min()), "t_max_used": float(t[ok].max()),
        "t_nonpositive_or_nan_in_rows_used": int((~(t[ok] >= 1)).sum())}
# who the missing ones are: proteins with no PMID left after the V50 correction have no year to date
nov = mapped & ~has_y0
summary["VT"] = {"proteins_mapped": int(mapped.sum()), "proteins_with_y0": int((mapped & has_y0).sum()),
                 "proteins_without_dated_pmid": int(nov.sum()),
                 "of_which_no_pmid_at_all_after_V50": int((nov & (out["n_V50"] == 0)).sum()),
                 "of_which_have_pmids_but_none_dated": int((nov & (out["n_V50"] > 0)).sum()),
                 "y0_min": float(out["y0"].min()), "y0_max": float(out["y0"].max())}
pgm = pg.merge(pairs.merge(gpp, on="pmid"), on="gene")
pgm = pgm[pgm["n_genes"] <= 50].merge(allyears, on="pmid", how="left")
gen = pgm.groupby("gene")["year"].apply(lambda x: x.notna().any())
all_genes = pg["gene"].drop_duplicates()
genes_with_pm = set(gen.index)
summary["VT"]["genes_mapped_total"] = int(len(all_genes))
summary["VT"]["genes_with_no_pmid_after_V50"] = int(len(all_genes) - len(genes_with_pm))
summary["VT"]["genes_with_pmids_but_none_with_a_resolved_year"] = int((~gen).sum())
summary["VT"]["genes_without_any_dated_pmid_total"] = int(len(all_genes) - int(gen.sum()))

# proxy check: y0 uses the K_LOW lowest-numbered PMIDs. Compare with the earliest year among ALL V50 PMIDs
# of the protein whose year is already known (site-linked PMIDs in pmid_year.parquet): a lower bound on how
# often the proxy misses an earlier publication.
allp = p50.merge(allyears, on="pmid", how="inner").groupby("uniprot_id")["year"].min().rename("y0_known_min")
cmp = out[["uniprot_id", "y0"]].merge(allp, on="uniprot_id")
cmp = cmp[cmp["y0"].notna()]
summary["y0_proxy_check"] = {"K_LOW": K_LOW, "proteins_compared": int(len(cmp)),
                             "proteins_with_earlier_known_year_than_proxy": int((cmp["y0_known_min"] < cmp["y0"]).sum()),
                             "median_years_earlier_when_different": float((cmp["y0"] - cmp["y0_known_min"])[cmp["y0_known_min"] < cmp["y0"]].median()) if (cmp["y0_known_min"] < cmp["y0"]).any() else 0.0}

# ---- agreement with the original index
cor = {}
for v in ["V20", "V50", "V100", "VT2025", "VT2026"]:
    c = f"SEI_{v}"
    ok = mapped & out[c].notna() & out["SEI_V0"].notna()  # VT: proteins with a dated PMID only
    cor[v] = {"spearman_with_original": float(stats.spearmanr(out.loc[ok, c], out.loc[ok, "SEI_V0"])[0]),
              "pearson_with_original": float(np.corrcoef(out.loc[ok, c], out.loc[ok, "SEI_V0"])[0, 1]),
              "n_proteins": int(ok.sum())}
summary["correlation_with_original"] = cor
out.to_parquet(PRIV / "index_variants.parquet")
json.dump(summary, open(RES / "d2_index_variants_summary.json", "w"), indent=2, default=str)
print(json.dumps(summary, indent=1, default=str))
