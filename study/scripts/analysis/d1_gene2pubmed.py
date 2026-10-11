"""Task D step 0: read the full cached gene2pubmed once.

Writes to study/analysis/data_private/ (git-ignored):
  g2p_pairs_3sp.parquet  distinct (tax_id, GeneID, PMID) for human, mouse, rat
  g2p_genes_per_pmid.parquet  PMID -> number of distinct GeneIDs over ALL organisms
"""
import gzip
import json

import numpy as np
import pandas as pd

from common import CUR, PRIV, RES

SRC = CUR / "processed_data" / "epsd" / "sei" / "raw_cache" / "gene2pubmed.gz"
TAX = {9606, 10090, 10116}

pairs, cnt_parts = [], []
n_rows = 0
for ch in pd.read_csv(SRC, sep="\t", usecols=[0, 1, 2], names=["tax", "gene", "pmid"], header=0,
                      dtype={"tax": np.int64, "gene": np.int64, "pmid": np.int64}, chunksize=5_000_000,
                      compression="gzip"):
    n_rows += len(ch)
    ch = ch.drop_duplicates()
    cnt_parts.append(ch[["gene", "pmid"]])
    pairs.append(ch[ch["tax"].isin(TAX)])
    print(n_rows, flush=True)

allp = pd.concat(cnt_parts).drop_duplicates()
gpp = allp.groupby("pmid").size().rename("n_genes").reset_index()
gpp.to_parquet(PRIV / "g2p_genes_per_pmid.parquet")
p3 = pd.concat(pairs).drop_duplicates()
p3.to_parquet(PRIV / "g2p_pairs_3sp.parquet")
summary = {"rows_read": n_rows, "distinct_gene_pmid_pairs_all_taxa": int(len(allp)),
           "distinct_pmids_all_taxa": int(len(gpp)), "pairs_3species": int(len(p3)),
           "pmids_gt_20_genes": int((gpp.n_genes > 20).sum()), "pmids_gt_50_genes": int((gpp.n_genes > 50).sum()),
           "pmids_gt_100_genes": int((gpp.n_genes > 100).sum())}
json.dump(summary, open(RES / "d1_gene2pubmed_summary.json", "w"), indent=2)
print(summary)
