"""Build the record-level analysis table used by every Rev 6 contrast.

One row per record of mammalian_sites.parquet (1,768,068). Adds: seq_motif_match, human in_domain and the
set-A tier (from human_sites_source_set_a.parquet, joined on accession + residue + position), the primary
study cluster of each site (label of the cluster of its globally most frequent PMID; from cluster_split.parquet),
and the protein-level index variants (index_variants.parquet, written by d2_index_variants.py).
Cached in study/analysis/data_private/base.parquet.
"""
import numpy as np
import pandas as pd

from common import PRIV, REL

K = ["uniprot_id", "residue_type", "position"]
CACHE = PRIV / "base.parquet"


def build():
    m = pd.read_parquet(REL / "mammalian_sites.parquet",
                        columns=K + ["species", "DAS", "study_cluster_count", "tier", "pmid_list_str",
                                     "source_PSP", "source_PhosphoELM", "source_dbPTM",
                                     "disorder", "kinase_match", "conserved"])
    sm = pd.read_parquet(REL / "kinase_motif_sequence.parquet", columns=K + ["seq_motif_match", "tier"])
    assert (sm[K].values == m[K].values).all() and (sm["tier"].values == m["tier"].values).all()
    m["seq_motif_match"] = sm["seq_motif_match"].values
    cs = pd.read_parquet(REL / "cluster_split.parquet", columns=K + ["study_cluster"])
    assert (cs[K].values == m[K].values).all()
    m["primary_cluster"] = cs["study_cluster"].values
    ha = pd.read_parquet(REL / "human_sites_source_set_a.parquet", columns=K + ["in_domain", "tier", "DAS"])
    ha = ha.drop_duplicates(K).rename(columns={"tier": "tier_setA", "DAS": "DAS_setA"})
    m = m.merge(ha, on=K, how="left")
    assert len(m) == 1768068
    m["in_domain"] = m["in_domain"].astype(float)
    iv = pd.read_parquet(PRIV / "index_variants.parquet").drop(columns=["species"])
    m = m.merge(iv, on="uniprot_id", how="left")
    assert len(m) == 1768068
    m.to_parquet(CACHE)
    return m


LEAN = ["uniprot_id", "gene_cluster", "species", "residue_type", "position", "DAS", "tier", "tier_setA", "primary_cluster", "mapped_to_gene",
        "in_domain", "conserved", "disorder", "seq_motif_match", "kinase_match",
        "SEI_V0", "SEI_V20", "SEI_V50", "SEI_V100", "SEI_VT2025", "SEI_VT2026", "SEI_secondary"]


def load(columns=None):
    """columns=None: every column (the lean set LEAN omits publication lists and source flags)."""
    if not CACHE.exists():
        build()
    want_gene = columns is None or columns == "lean" or "gene_cluster" in columns
    if columns == "lean":
        import pyarrow.parquet as pq
        names = set(pq.ParquetFile(CACHE).schema.names)
        columns = [c for c in LEAN if c in names]
    if columns is not None:
        columns = [c for c in columns if c != "gene_cluster"]
    d = pd.read_parquet(CACHE, columns=columns)
    if want_gene:
        d["gene_cluster"] = d["uniprot_id"].map(gene_clusters())
    return d


def gene_clusters():
    """accession -> gene cluster label (see gene_unit.py)."""
    u = pd.read_parquet(PRIV / "gene_unit.parquet", columns=["uniprot_id", "gene_cluster"])
    return u.set_index("uniprot_id")["gene_cluster"]


if __name__ == "__main__":
    b = build()
    print(b.shape)
    print(b.columns.tolist())
