"""Gene-level resampling unit (plan amendment of 2026-10-01).

The site table has several accessions per gene, so accessions of one gene are not independent resampling units.
gene_cluster = connected component of the accession-gene graph (accessions that share an Entrez gene ID are in the
same cluster), labelled "G<smallest gene ID in the component>". Accessions without a gene mapping are their own cluster
("A<accession>"); none of them has a primary index, so none enters a decision-bearing contrast.

Writes data_private/gene_unit.parquet and results/gene_unit_summary.json.
"""
import json

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from common import PRIV, REL, RES

sei = pd.read_parquet(REL / "sei_protein.parquet", columns=["uniprot_id", "species", "gene_ids", "mapped_to_gene"])
m = sei[sei["mapped_to_gene"] == True].copy()  # noqa: E712
pg = m[["uniprot_id", "species", "gene_ids"]].explode("gene_ids").dropna(subset=["gene_ids"])
pg["gene"] = pg["gene_ids"].astype(np.int64)
pg = pg.drop(columns="gene_ids").drop_duplicates()

acc = pd.factorize(pg["uniprot_id"])
gen = pd.factorize(pg["gene"])
na, ng = len(acc[1]), len(gen[1])
A = coo_matrix((np.ones(len(pg)), (acc[0], na + gen[0])), shape=(na + ng, na + ng))
ncomp, lab = connected_components(A, directed=False)
pg["comp"] = lab[acc[0]]
mn = pg.groupby("comp")["gene"].min()
first = pg.groupby("uniprot_id").agg(comp=("comp", "first"), n_gene_ids=("gene", "nunique"), gene_first=("gene", "min"), species=("species", "first"))
first["gene_cluster"] = "G" + first["comp"].map(mn).astype(str)
first["entrez_gene_id"] = first["gene_first"]
u = sei[["uniprot_id", "species", "mapped_to_gene"]].merge(first[["gene_cluster", "entrez_gene_id", "n_gene_ids"]].reset_index(), on="uniprot_id", how="left")
u["gene_cluster"] = u["gene_cluster"].fillna("A" + u["uniprot_id"])
u.to_parquet(PRIV / "gene_unit.parquet")

summ = {"definition": "entrez_gene_id = smallest Entrez gene ID of the accession; gene_cluster = connected component of accessions sharing a gene ID"}
for sp in ("human", "mouse", "rat", "all"):
    x = u if sp == "all" else u[u["species"] == sp]
    mp = x[x["mapped_to_gene"] == True]  # noqa: E712
    per = mp.groupby("gene_cluster").size()
    summ[sp] = {"accessions": int(len(x)), "accessions_with_gene": int(len(mp)), "accessions_without_gene": int(len(x) - len(mp)),
                "gene_clusters_among_mapped": int(len(per)),
                "accessions_per_gene_cluster": {"mean": float(per.mean()), "median": float(per.median()), "max": int(per.max()),
                                                "share_of_clusters_with_more_than_one_accession": float((per > 1).mean()),
                                                "distribution": {str(k): int(v) for k, v in per.value_counts().sort_index().head(12).items()}},
                "accessions_with_more_than_one_gene_id": int((mp["n_gene_ids"] > 1).sum())}
# records in the largest genes, over the whole site table and over the analysis sites (human, with primary index)
import base as B
b = B.load(["uniprot_id", "species", "tier", "SEI_V50"])
b = b.merge(u[["uniprot_id", "gene_cluster", "entrez_gene_id"]], on="uniprot_id", how="left")
for name, sel in {"all_records": b["uniprot_id"].notna(), "human_records": b["species"] == "human",
                  "gold_or_bronze_with_index": b["tier"].isin(["GOLD", "BRONZE"]) & b["SEI_V50"].notna()}.items():
    x = b[sel]
    cnt = x["gene_cluster"].value_counts()
    acc_cnt = x["uniprot_id"].value_counts()
    top = cnt.head(10)
    summ[f"largest_genes_{name}"] = {"n_records": int(len(x)), "n_gene_clusters": int(len(cnt)), "n_accessions": int(len(acc_cnt)),
                                    "records_in_10_largest_genes": int(top.sum()), "share": float(top.sum() / len(x)),
                                    "records_in_10_largest_accessions": int(acc_cnt.head(10).sum()),
                                    "share_10_largest_accessions": float(acc_cnt.head(10).sum() / len(x)),
                                    "top10": [{"gene_cluster": k, "records": int(v),
                                               "accessions": int(x.loc[x["gene_cluster"] == k, "uniprot_id"].nunique())} for k, v in top.items()]}
json.dump(summ, open(RES / "gene_unit_summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summ.items() if not k.startswith("largest")}, indent=1)[:2500])
for k in summ:
    if k.startswith("largest"):
        print(k, {a: summ[k][a] for a in ("n_records", "records_in_10_largest_genes", "share", "records_in_10_largest_accessions", "share_10_largest_accessions")})
