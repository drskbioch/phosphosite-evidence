"""Data preparation for block G21 (within-gene GOLD vs BRONZE; definition in the plan, deviation log of 2026-10-04, fourth batch).

prepare(signal, population) returns the records of one signal-by-population group, sorted by gene cluster, with
  gold (1 = GOLD), isT, isY (residue indicators; S is the reference), has_seq (protein sequence available and the residue at the stated position
  equals the record's residue type), det (tryptic peptide containing the site has 7 to 35 residues; cleavage after K or R but not before P).
"""
import numpy as np
import pandas as pd

import run_cem_grid as R
from common import PRIV
from populations import SIGNALS, population

STY = ["S", "T", "Y"]
_SEQ = None


def sequences():
    global _SEQ
    if _SEQ is None:
        _SEQ = {}
        for sp in ("9606", "10090", "10116"):
            d = pd.read_csv(PRIV / "dl" / f"seq_{sp}.tsv", sep="\t", usecols=["Entry", "Sequence"])
            _SEQ.update(dict(zip(d["Entry"], d["Sequence"])))
    return _SEQ


def tryptic_length(seq, pos):
    """Length of the tryptic peptide containing 1-based position pos (cleavage after K or R unless the next residue is P); None if pos is outside the sequence."""
    n = len(seq)
    if pos < 1 or pos > n:
        return None
    s = np.frombuffer(seq.encode(), dtype=np.uint8)
    nxt = np.append(s[1:], ord("*"))
    cl = np.nonzero(((s == ord("K")) | (s == ord("R"))) & (nxt != ord("P")))[0] + 1  # 1-based index of the residue after which cleavage occurs
    prev = cl[cl < pos]
    nx = cl[cl >= pos]
    start = (prev[-1] if len(prev) else 0)
    end = (nx[0] if len(nx) else n)
    return int(end - start)


def prepare(signal, popl):
    b = R.b
    d = population(b, signal, tier_col="tier", index_col="SEI_V50", species="human" if popl.startswith("P-human") else None, residues=STY)
    d["position"] = b.loc[d.index, "position"].to_numpy()
    d["gold"] = (d["tier"] == "GOLD").astype(np.int8)
    d["isT"] = (d["residue_type"] == "T").astype(np.int8)
    d["isY"] = (d["residue_type"] == "Y").astype(np.int8)
    seqs = sequences()
    has = np.zeros(len(d), bool)
    det = np.zeros(len(d), np.int8)
    nomatch = 0
    nonseq = 0
    acc = d["uniprot_id"].to_numpy()
    pos = d["position"].to_numpy()
    res = d["residue_type"].to_numpy()
    cache = {}
    for i in range(len(d)):
        a = acc[i]
        s = seqs.get(a)
        if s is None:
            nonseq += 1
            continue
        p = int(pos[i])
        if p < 1 or p > len(s) or s[p - 1] != res[i]:
            nomatch += 1
            continue
        k = (a, p)
        L = tryptic_length(s, p) if k not in cache else cache[k]
        cache[k] = L
        has[i] = True
        det[i] = 1 if 7 <= L <= 35 else 0
    d["has_seq"] = has
    d["det"] = det
    d = d.sort_values("gene_cluster", kind="stable").reset_index(drop=True)
    info = {"records": int(len(d)), "gold": int(d["gold"].sum()), "bronze": int((1 - d["gold"]).sum()), "gene_clusters": int(d["gene_cluster"].nunique()),
            "records_without_sequence": nonseq, "records_residue_mismatch": nomatch, "records_with_sequence_and_match": int(has.sum()),
            "gene_clusters_with_both_tiers": int((d.groupby("gene_cluster")["gold"].agg(["min", "max"]).pipe(lambda x: (x["min"] == 0) & (x["max"] == 1))).sum())}
    return d, info
