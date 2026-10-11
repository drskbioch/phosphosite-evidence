"""Task J / block G11: the predictor audit with a gene-cluster bootstrap and an equivalence test on the AUROC gap.

Sample: processed_data/epsd/predictor_audit/eval_sample_with_sequence.parquet (frozen Rev 5 sample, 10,000 sites; deciles of the original primary index).
Gap = AUROC(deciles 8-10) - AUROC(deciles 1-3), both pooled over the three deciles. Interval: gene-cluster bootstrap (resample gene clusters with replacement
within the bottom pool and within the top pool; 1,000 resamples). The site-level interval of the Rev 5 procedure (resample sites within the pool) is computed
next to it. Equivalence (plan section 2): the 90% percentile interval inside (-0.05, +0.05) -> "gap absent"; otherwise "gap present" if the interval excludes zero and
lies outside the margin, else "not detectable". Restricted analyses keep the positives first reported after the restriction year (all negatives kept) and
need >= 200 positives, as in Rev 5.
Rows: DeepPhos and MusiteDeep on S/T sites, PhosF3C on Y sites; the Rev 5 design (all residues scored by the tool) is reported beside each.
Training overlap (J4): PhosF3C Y only (its training set is in the repository); DeepPhos and MusiteDeep training sets are not available locally.
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

import base as B  # noqa: E402
from common import CUR, N_BOOT, RES, rng_for  # noqa: E402

sys.path.insert(0, str(CUR / "scripts" / "earlier_predictor_audit" / "step3"))
import score_metrics as SM  # noqa: E402

# score_metrics.py sets absolute default locations for the score files; point them at this checkout (the module itself is not edited).
SM.OUT_DIR = str(CUR / "processed_data" / "epsd" / "predictor_audit")
SM.PHOSF3C_RESULT_PATH = str(CUR / "tools" / "phosf3c" / "PhosF3C-main" / "result" / "task_b" / "task_b_sample_Y.json")

SAMPLE = CUR / "processed_data" / "epsd" / "predictor_audit" / "eval_sample_with_sequence.parquet"
MARGIN = 0.05
MIN_POS = 200
BOTTOM, TOP = (1, 2, 3), (8, 9, 10)
sample = pd.read_parquet(SAMPLE)
sample["label"] = (sample["label_type"] == "positive").astype(int)
gc = B.gene_clusters()
sample["gene_cluster"] = sample["uniprot_id"].map(gc)
sample["gene_cluster"] = sample["gene_cluster"].fillna("A" + sample["uniprot_id"])
rev5 = json.load(open(CUR / "earlier_results" / "predictor_metrics.json"))
LOADERS = {"DeepPhos": SM.load_deepphos_scores, "MusiteDeep": SM.load_musitedeep_scores, "PhosF3C": SM.load_phosf3c_scores}


def auc(y, s):
    n1 = y.sum()
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(s)
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def pools(df):
    b = df[df["SEI_decile"].isin(BOTTOM)]
    t = df[df["SEI_decile"].isin(TOP)]
    return b, t


def boot_gap(df, rng, by):
    """by: "gene" (gene-cluster resampling inside each pool) or "site"."""
    b, t = pools(df)
    parts = []
    for p in (b, t):
        y, s = p["label"].to_numpy(), p["score"].to_numpy(float)
        idx = None
        if by == "gene":
            codes, uniq = pd.factorize(p["gene_cluster"])
            order = np.argsort(codes, kind="stable")
            cuts = np.searchsorted(codes[order], np.arange(len(uniq) + 1))
            idx = [order[cuts[u]:cuts[u + 1]] for u in range(len(uniq))]
        parts.append((y, s, idx))
    out = np.full(N_BOOT, np.nan)
    for i in range(N_BOOT):
        a = []
        for y, s, idx in parts:
            if by == "gene":
                sel = np.concatenate([idx[j] for j in rng.integers(0, len(idx), len(idx))])
            else:
                sel = rng.integers(0, len(y), len(y))
            a.append(auc(y[sel], s[sel]))
        out[i] = a[1] - a[0]
    return out[~np.isnan(out)]


def label_of(ci90):
    lo, hi = ci90
    if lo > -MARGIN and hi < MARGIN:
        return "gap absent (equivalent within +/-0.05)"
    if (lo > 0 or hi < 0) and (lo >= MARGIN or hi <= -MARGIN):
        return "gap present (interval excludes zero and lies outside the margin)"
    return "not detectable (interval includes the margin)"


def analyse(df, rng_id):
    b, t = pools(df)
    res = {"n_sites": int(len(df)), "n_positives": int(df["label"].sum()), "n_gene_clusters": int(df["gene_cluster"].nunique()),
           "per_decile": {int(d): {"n": int((df["SEI_decile"] == d).sum()), "positives": int(((df["SEI_decile"] == d) & (df["label"] == 1)).sum())} for d in range(1, 11)}}
    if len(b) == 0 or len(t) == 0 or b["label"].nunique() < 2 or t["label"].nunique() < 2:
        res["status"] = "single class in a pool, gap undefined"
        return res
    ab, at = auc(b["label"].to_numpy(), b["score"].to_numpy(float)), auc(t["label"].to_numpy(), t["score"].to_numpy(float))
    res.update(status="evaluated", bottom_auroc=ab, top_auroc=at, gap=at - ab, n_bottom=int(len(b)), n_top=int(len(t)))
    for by in ("gene", "site"):
        bs = boot_gap(df, rng_for(f"G11|{rng_id}|{by}"), by)
        res[f"{by}_ci95"] = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
        res[f"{by}_ci90"] = [float(np.percentile(bs, 5)), float(np.percentile(bs, 95))]
        res[f"{by}_label"] = label_of(res[f"{by}_ci90"])
        res[f"{by}_n_valid"] = int(len(bs))
    return res


out = {"settings": {"margin": MARGIN, "n_boot": N_BOOT, "min_positives_restricted": MIN_POS, "seed": 20260924,
                    "resampling_unit": "gene cluster within each decile pool (primary); sites within each pool (Rev 5 procedure, shown beside it)"}}
merged = {}
for tool, load in LOADERS.items():
    sc = load()
    m = sample.merge(sc, on=["uniprot_id", "position"], how="inner")
    merged[tool] = m
    resid = ["Y"] if tool == "PhosF3C" else ["S", "T"]
    for design, d in (("planned: " + "/".join(resid) + " sites only", m[m["residue_type"].isin(resid)]), ("Rev 5 design: all scored residues", m)):
        rows = {"unrestricted": analyse(d, f"{tool}|{design}|unrestricted")}
        years = [2019, 2018] if tool == "DeepPhos" else ([2020, 2017] if tool == "MusiteDeep" else [2024])
        for yr in years:
            pos = d[(d["label"] == 1) & (d["first_reported_year"] > yr)]
            neg = d[d["label"] == 0]
            if len(pos) < MIN_POS:
                rows[f"restriction_year={yr}"] = {"status": f"too small to evaluate (n<{MIN_POS})", "n_positives": int(len(pos))}
                continue
            rows[f"restriction_year={yr}"] = analyse(pd.concat([pos, neg]), f"{tool}|{design}|{yr}")
        out[f"{tool}|{design}"] = rows
        for k, r in rows.items():
            print(tool, design, k, r.get("status"), round(r.get("gap", np.nan), 4), r.get("gene_ci95"), r.get("site_ci95"), r.get("gene_label"), flush=True)
    r5 = rev5[tool]
    out[f"{tool}|rev5_reported"] = {"gap": r5["gap_top_minus_bottom"], "restricted": {str(r5["model_pub_year"]): r5["leakage_free"], **r5.get("restriction_year_sensitivity", {})}}

# ---------------------------------------------------------------- J4: training-set overlap (PhosF3C, tyrosine)
TRAIN = CUR / "tools" / "phosf3c" / "PhosF3C-main" / "data" / "merged_train_Y.fasta"
HALF = 15  # hparams window_size: 15 on each side of the site (31-residue window)
tw = []
name, seq = None, []
for line in open(TRAIN):
    line = line.strip()
    if line.startswith(">"):
        if name is not None:
            tw.append((name, "".join(seq)))
        name, seq = line[1:], []
    else:
        seq.append(line)
tw.append((name, "".join(seq)))
AA = "ACDEFGHIKLMNPQRSTVWY"
code = {a: i for i, a in enumerate(AA)}


def enc(window):
    return np.array([code.get(c, -1) for c in window], dtype=np.int8)


tr = []
for nm, s in tw:
    plain, sites = [], []
    for ch in s:
        if ch == "#":
            sites.append(len(plain) - 1)
        else:
            plain.append(ch)
    p = "".join(plain)
    padded = "X" * HALF + p + "X" * HALF
    for i in sites:
        tr.append(enc(padded[i: i + 2 * HALF + 1]))
tr = np.array(tr)
m = merged["PhosF3C"]
pos = m[m["label"] == 1]
ev = np.array([enc(("X" * HALF + sq + "X" * HALF)[p - 1: p - 1 + 2 * HALF + 1]) for sq, p in zip(pos["sequence"], pos["position"])])
ident = np.zeros(len(ev))
for i in range(len(ev)):
    ident[i] = ((tr == ev[i]) & (ev[i] >= 0)).sum(axis=1).max() / (2 * HALF + 1)
flag = ident >= 0.40
flagged_index = pos.index[flag]
keep = m[~m.index.isin(flagged_index)]
j4 = {"training_file": "tools/phosf3c/PhosF3C-main/data/merged_train_Y.fasta (shipped with the PhosF3C repository)", "window": f"{2 * HALF + 1} residues (+/- {HALF}, hparams window_size 15)",
      "n_training_positive_windows": int(len(tr)), "n_eval_positives_Y": int(len(pos)), "n_flagged_ge_40pct": int(flag.sum()), "share_flagged": float(flag.mean()),
      "identity_quantiles_eval_positives": {str(q): float(np.quantile(ident, q)) for q in (0.1, 0.5, 0.9, 0.99)},
      "analysis_without_flagged_positives": analyse(keep, "PhosF3C|J4|overlap_removed"),
      "analysis_all_for_reference": analyse(m, "PhosF3C|J4|all_Y")}
out["J4_training_overlap"] = {"PhosF3C": j4,
                              "DeepPhos": "not run: no training set in the repository (tools/DeepPhos holds weights and a 944-byte test file)",
                              "MusiteDeep": "not run: the repository holds no training data for the deployed phosphorylation models (tool_versions.json: training sources not stated)"}
print(json.dumps({k: v for k, v in j4.items() if not k.startswith("analysis")}, indent=1))
json.dump(out, open(RES / "g11_predictors.json", "w"), indent=1, default=float)
