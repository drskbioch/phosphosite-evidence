"""Study-Effort Index (SEI) computation.

The primary index is log10(1 + n_pubmed), where n_pubmed is the number of
distinct PubMed IDs linked to a protein's NCBI Gene ID(s) in gene2pubmed.
It is defined from an external literature-to-gene mapping, independent of
any site-level evidence score. A secondary index counts distinct study
clusters instead of PMIDs.
"""
import gzip
from collections import Counter

import numpy as np
import pandas as pd


class _UnionFind:
    """Minimal union-find over hashable items."""

    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def compute_sei_primary(n_pubmed):
    """SEI_primary = log10(1 + n_pubmed). Vectorized; n_pubmed may be a
    scalar, list, or array-like of non-negative publication counts."""
    return np.log10(1 + np.asarray(n_pubmed, dtype=float))


def load_uniprot_to_geneid(path):
    """UniProtKB-AC -> set of GeneID strings, from a by_organism
    idmapping.dat.gz file (UniProt format: AC \\t ID_type \\t value).
    """
    m = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            uid, id_type, val = parts
            if id_type == "GeneID":
                m.setdefault(uid, set()).add(val)
    return m


def load_gene2pubmed_filtered(path, tax_ids):
    """GeneID -> set of PubMed_ID strings, restricted to given tax_ids."""
    tax_strs = {str(t) for t in tax_ids}
    m = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        fh.readline()  # '#tax_id\tGeneID\tPubMed_ID'
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            tax_id, gene_id, pmid = parts
            if tax_id in tax_strs:
                m.setdefault(gene_id, set()).add(pmid)
    return m


def build_sei_table(proteins, gene2pubmed_path, idmap_paths, tax_names=None):
    """Protein-level primary SEI table.

    proteins: DataFrame with columns uniprot_id, taxonomy_id.
    gene2pubmed_path: path to gene2pubmed.gz.
    idmap_paths: dict {taxonomy_id: path to that species' idmapping.dat.gz}.
    tax_names: optional dict {taxonomy_id: species name} for the 'species'
        output column.

    Returns a DataFrame: uniprot_id, taxonomy_id, species, gene_ids,
    n_gene_ids, n_pubmed, SEI_primary, mapped_to_gene, zero_pubmed_flag,
    matched_via_base_accession.

    Isoform accessions are matched by stripping everything from the first
    '-' when the full accession has no gene mapping.
    """
    gene2pmid = load_gene2pubmed_filtered(gene2pubmed_path, set(idmap_paths.keys()))
    uid_to_gene = {tax_id: load_uniprot_to_geneid(path) for tax_id, path in idmap_paths.items()}

    rows = []
    for uid, tax_id in zip(proteins["uniprot_id"], proteins["taxonomy_id"]):
        m = uid_to_gene.get(tax_id, {})
        gene_ids = m.get(uid)
        used_base = False
        if not gene_ids and "-" in uid:
            base = uid.split("-")[0]
            gene_ids = m.get(base)
            used_base = gene_ids is not None
        if not gene_ids:
            rows.append((uid, tax_id, [], 0, np.nan, False, False, used_base))
            continue
        pmids = set()
        any_zero_pmid_gene = False
        for g in gene_ids:
            gp = gene2pmid.get(g)
            if gp:
                pmids |= gp
            else:
                any_zero_pmid_gene = True
        n_pubmed = len(pmids)
        sei = float(compute_sei_primary(n_pubmed))
        rows.append((uid, tax_id, sorted(gene_ids), n_pubmed, sei, True, n_pubmed == 0, used_base))

    out = pd.DataFrame(rows, columns=[
        "uniprot_id", "taxonomy_id", "gene_ids", "n_pubmed", "SEI_primary",
        "mapped_to_gene", "zero_pubmed_flag", "matched_via_base_accession",
    ])
    if tax_names:
        out["species"] = out["taxonomy_id"].map(tax_names)
    out["n_gene_ids"] = out["gene_ids"].apply(len)
    return out


def build_sei_secondary_table(site_publications, site_to_protein, high_freq_min=500,
                               jaccard_threshold=0.90):
    """Secondary SEI: distinct study clusters per protein. PMIDs that
    appear on many sites and have near-identical site co-occurrence
    (Jaccard similarity) are merged into a single study cluster.

    site_publications: DataFrame with columns site_id, pmid.
    site_to_protein: DataFrame with columns site_id, uniprot_id.
    high_freq_min: minimum #sites a PMID must appear on to be clustered
        (PMIDs below this are left as their own singleton cluster).
    jaccard_threshold: co-occurrence Jaccard similarity above which two
        high-frequency PMIDs are merged into the same study-cluster.

    Returns a DataFrame: uniprot_id, n_study_clusters_protein, SEI_secondary.
    """
    pubs = site_publications.copy()
    pubs["pmid"] = pubs["pmid"].astype(str)

    site_pmids = pubs.groupby("site_id")["pmid"].apply(set)

    freq = Counter()
    for s in site_pmids:
        freq.update(s)
    high_freq = {p for p, n in freq.items() if n >= high_freq_min}
    hf_index = {p: i for i, p in enumerate(sorted(high_freq))}
    hf_list = sorted(high_freq)

    pair_counts = Counter()
    for s in site_pmids:
        hf_present = sorted({hf_index[p] for p in s if p in hf_index})
        n = len(hf_present)
        if n < 2:
            continue
        for i in range(n):
            for j in range(i + 1, n):
                pair_counts[(hf_present[i], hf_present[j])] += 1

    uf = _UnionFind()
    for (i, j), inter in pair_counts.items():
        pi, pj = hf_list[i], hf_list[j]
        union_size = freq[pi] + freq[pj] - inter
        jac = inter / union_size if union_size else 0
        if jac >= jaccard_threshold:
            uf.union(pi, pj)

    cluster_id = {p: "CLUSTER_" + uf.find(p) for p in hf_list}

    def get_cluster(p):
        return cluster_id.get(p, p)

    site_to_uid_map = dict(zip(site_to_protein["site_id"], site_to_protein["uniprot_id"]))
    protein_clusters = {}
    for site_id, pmids in site_pmids.items():
        uid = site_to_uid_map.get(site_id)
        if uid is None:
            continue
        clusters = {get_cluster(p) for p in pmids}
        protein_clusters.setdefault(uid, set()).update(clusters)

    rows = [(uid, len(clusters)) for uid, clusters in protein_clusters.items()]
    out = pd.DataFrame(rows, columns=["uniprot_id", "n_study_clusters_protein"])
    # Same log transform as the primary index, applied to cluster counts.
    out["SEI_secondary"] = np.log10(1 + out["n_study_clusters_protein"])
    return out


def sei_deciles(sei_values, n_deciles=10):
    """SEI decile assignment via pd.qcut.

    sei_values: array-like or pandas Series of SEI values (primary or
        secondary; the function is definition-agnostic).
    n_deciles: number of quantile bins.

    Returns (deciles, n_actual_bins). deciles is 0-indexed (0..n-1); add 1
    for the 1..10 coding used by compute_tool_metrics. Duplicate bin edges
    are dropped, so n_actual_bins can be smaller than n_deciles; it is 0
    (all-NaN deciles) when every value is identical.
        """
    series = sei_values if isinstance(sei_values, pd.Series) else pd.Series(sei_values)
    try:
        deciles, bins = pd.qcut(series, q=n_deciles, labels=False, duplicates="drop", retbins=True)
        n_actual_bins = int(deciles.max() + 1) if deciles.notna().any() else 0
        return deciles, n_actual_bins
    except ValueError:
        return pd.Series(0, index=series.index), 1
