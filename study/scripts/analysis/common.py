"""Shared code for the Rev 6 reanalysis (plan: study/analysis/ANALYSIS_PLAN.md).

Everything here is deterministic given SEED. Estimators reproduce the phosnet
coarsened exact matching point estimates exactly (checked in validate_engine.py)
and add a protein-cluster bootstrap that re-runs the whole pipeline (index
deciles, strata, weights, estimate) inside every resample.
"""
import zlib
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CUR = ROOT / "study"
REL = CUR / "processed_data" / "epsd" / "release"
REV6 = CUR / "analysis"
RES = REV6 / "results"
PRIV = REV6 / "data_private"
for d in (RES, PRIV):
    d.mkdir(parents=True, exist_ok=True)

SEED = 20260924
N_BOOT = 1000
JACCARD = 0.90


def rng_for(contrast_id):
    """Generator derived from the base seed and the contrast ID."""
    return np.random.default_rng([SEED, zlib.crc32(contrast_id.encode())])


# ------------------------------------------------------------------ study clusters
class _UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def split_pmids(series):
    return [s.split(";") if s else [] for s in series.fillna("")]


def cluster_map(pmid_lists, hf_min, jaccard=JACCARD):
    """PMID -> cluster label. A PMID linked to >= hf_min sites is a candidate;
    candidates whose site sets have Jaccard similarity >= jaccard are merged
    (union-find, transitive). Other PMIDs are their own cluster."""
    freq = Counter(p for l in pmid_lists for p in set(l))
    hf = sorted(p for p, n in freq.items() if n >= hf_min)
    idx = {p: i for i, p in enumerate(hf)}
    sets = Counter()
    for l in pmid_lists:
        h = tuple(sorted({idx[p] for p in l if p in idx}))
        if len(h) >= 2:
            sets[h] += 1
    pair = Counter()
    for h, n in sets.items():
        for i in range(len(h)):
            for j in range(i + 1, len(h)):
                pair[(h[i], h[j])] += n
    uf = _UF()
    for (i, j), inter in pair.items():
        u = freq[hf[i]] + freq[hf[j]] - inter
        if u and inter / u >= jaccard:
            uf.union(hf[i], hf[j])
    cid = {p: "C_" + uf.find(p) for p in hf}
    info = {"n_pmids": len(freq), "n_candidates": len(hf),
            "n_candidates_in_multi_pmid_cluster": None}
    from collections import defaultdict
    members = defaultdict(list)
    for p, c in cid.items():
        members[c].append(p)
    multi = {c: v for c, v in members.items() if len(v) > 1}
    info["n_candidate_clusters"] = len(members)
    info["n_multi_pmid_clusters"] = len(multi)
    info["n_candidates_in_multi_pmid_cluster"] = int(sum(len(v) for v in multi.values()))
    info["multi_cluster_sizes"] = sorted((len(v) for v in multi.values()), reverse=True)
    return cid, info, freq


def depth_per_site(pmid_lists, cid, drop=frozenset()):
    """Distinct clusters among each site's PMIDs (PMIDs in `drop` removed first)."""
    out = np.empty(len(pmid_lists), dtype=np.int32)
    for k, l in enumerate(pmid_lists):
        out[k] = len({cid.get(p, p) for p in l if p not in drop})
    return out


def assign_tier(das, depth, gold_depth=10, gold_das=2, silver_depth=3, silver_das=1, use_das=True):
    """GOLD if agreement >= gold_das or depth >= gold_depth; SILVER if agreement >= silver_das
    or depth >= silver_depth; else BRONZE. use_das=False gives the depth-only rule."""
    das = np.asarray(das)
    depth = np.asarray(depth)
    if use_das:
        gold = (das >= gold_das) | (depth >= gold_depth)
        silver = (das >= silver_das) | (depth >= silver_depth)
    else:
        gold = depth >= gold_depth
        silver = depth >= silver_depth
    return np.where(gold, "GOLD", np.where(silver, "SILVER", "BRONZE"))
