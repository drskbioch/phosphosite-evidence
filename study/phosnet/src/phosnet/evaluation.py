"""SEI-stratified predictor evaluation.

Computes per-decile AUROC, the difference in AUROC between the
best-studied and least-studied SEI deciles, pooled AUROC, and an optional
secondary analysis restricted to positives first reported after a model's
publication year. Bootstrap confidence intervals use an explicit seed;
each interval derives its own random generator from that seed, so a
result depends only on its data and seed, never on call order.
"""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

N_BOOT = 1000
MIN_LEAKAGE_FREE_POSITIVES = 200
BOTTOM_DECILES = (1, 2, 3)
TOP_DECILES = (8, 9, 10)


def bootstrap_auroc_ci(y_true, y_score, seed, n_boot=N_BOOT):
    """Point AUROC plus a percentile bootstrap 95% CI (n_boot resamples,
    generator seeded from `seed`)."""
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    n = len(y_true)
    if len(np.unique(y_true)) < 2:
        return None, None, None
    point = roc_auc_score(y_true, y_score)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        yt, ys = y_true[idx], y_score[idx]
        boots[i] = roc_auc_score(yt, ys) if len(np.unique(yt)) >= 2 else np.nan
    valid = boots[~np.isnan(boots)]
    if len(valid) == 0:
        return float(point), None, None
    lo, hi = np.percentile(valid, [2.5, 97.5])
    return float(point), float(lo), float(hi)


def bootstrap_gap_ci(y_true_bottom, y_score_bottom, y_true_top, y_score_top, seed, n_boot=N_BOOT):
    """top_auroc - bottom_auroc, plus a percentile bootstrap 95% CI
    (top and bottom groups resampled independently)."""
    yb, sb = np.asarray(y_true_bottom), np.asarray(y_score_bottom)
    yt, st = np.asarray(y_true_top), np.asarray(y_score_top)
    point = roc_auc_score(yt, st) - roc_auc_score(yb, sb)
    rng = np.random.default_rng(seed)
    nb, nt = len(yb), len(yt)
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        idx_b = rng.integers(0, nb, nb)
        idx_t = rng.integers(0, nt, nt)
        yb_i, yt_i = yb[idx_b], yt[idx_t]
        if len(np.unique(yb_i)) < 2 or len(np.unique(yt_i)) < 2:
            diffs[i] = np.nan
            continue
        diffs[i] = roc_auc_score(yt_i, st[idx_t]) - roc_auc_score(yb_i, sb[idx_b])
    valid = diffs[~np.isnan(diffs)]
    lo, hi = np.percentile(valid, [2.5, 97.5])
    return float(point), float(lo), float(hi)


def compute_tool_metrics(merged, seed=20260924, label_col="label_type", score_col="score",
                          decile_col="SEI_decile", positive_label="positive",
                          model_pub_year=None, first_reported_year_col="first_reported_year",
                          min_leakage_free_positives=MIN_LEAKAGE_FREE_POSITIVES, n_boot=N_BOOT):
    """Metrics for one predictor's scored sample: per-decile AUROC, the
    bottom-three vs. top-three decile AUROC difference ("gap"), pooled
    AUROC, the mean of per-decile AUROCs, and (if model_pub_year is given)
    the analysis restricted to positives first reported after that year.

    merged: one row per scored site, with label_col (values including
        positive_label for true positives), score_col, and decile_col
        (SEI decile coded 1..10).
    seed: base seed. The per-decile, gap, pooled and restricted-analysis
        bootstraps each derive their own seed from it.

    Returns a dict with keys n_scored, per_decile, gap_top_minus_bottom,
    pooled_auroc, sei_stratified_mean_auroc and, when requested,
    leakage_free.
    """
    y = (merged[label_col] == positive_label).astype(int).to_numpy()
    s = merged[score_col].to_numpy(dtype=float)
    result = {"n_scored": int(len(merged))}

    deciles_sorted = sorted(merged[decile_col].unique())
    decile_results = {}
    for decile in deciles_sorted:
        sub = merged[merged[decile_col] == decile]
        yy = (sub[label_col] == positive_label).astype(int).to_numpy()
        ss = sub[score_col].to_numpy(dtype=float)
        point, lo, hi = bootstrap_auroc_ci(yy, ss, seed=seed + 1000 + int(decile), n_boot=n_boot)
        decile_results[int(decile)] = {"n": len(sub), "auroc": point, "ci_low": lo, "ci_high": hi}
    result["per_decile"] = decile_results

    bottom3, top3 = deciles_sorted[:3], deciles_sorted[-3:]
    bottom_sub = merged[merged[decile_col].isin(bottom3)]
    top_sub = merged[merged[decile_col].isin(top3)]
    yb = (bottom_sub[label_col] == positive_label).astype(int).to_numpy()
    sb = bottom_sub[score_col].to_numpy(dtype=float)
    yt = (top_sub[label_col] == positive_label).astype(int).to_numpy()
    st = top_sub[score_col].to_numpy(dtype=float)
    gap_point, gap_lo, gap_hi = bootstrap_gap_ci(yb, sb, yt, st, seed=seed + 2000, n_boot=n_boot)
    result["gap_top_minus_bottom"] = {
        "point": gap_point, "ci_low": gap_lo, "ci_high": gap_hi,
        "bottom_deciles": bottom3, "top_deciles": top3,
        "n_bottom": len(bottom_sub), "n_top": len(top_sub),
    }

    pooled_point, pooled_lo, pooled_hi = bootstrap_auroc_ci(y, s, seed=seed, n_boot=n_boot)
    result["pooled_auroc"] = {"point": pooled_point, "ci_low": pooled_lo, "ci_high": pooled_hi}

    decile_aurocs = [v["auroc"] for v in decile_results.values() if v["auroc"] is not None]
    result["sei_stratified_mean_auroc"] = float(np.mean(decile_aurocs)) if decile_aurocs else None

    if model_pub_year is not None:
        positives = merged[merged[label_col] == positive_label]
        negatives = merged[merged[label_col] != positive_label]
        lf_positives = positives[positives[first_reported_year_col] > model_pub_year]
        n_lf_pos = len(lf_positives)
        lf = {"model_pub_year": model_pub_year, "n_positives_after_cutoff": int(n_lf_pos)}
        if n_lf_pos < min_leakage_free_positives:
            lf["status"] = f"too small to evaluate (n<{min_leakage_free_positives})"
        else:
            lf_merged = pd.concat([lf_positives, negatives], ignore_index=True)
            yl = (lf_merged[label_col] == positive_label).astype(int).to_numpy()
            sl = lf_merged[score_col].to_numpy(dtype=float)
            lf_point, lf_lo, lf_hi = bootstrap_auroc_ci(yl, sl, seed=seed + 3000, n_boot=n_boot)
            lf["pooled_auroc"] = {"point": lf_point, "ci_low": lf_lo, "ci_high": lf_hi}
            lf["status"] = "evaluated"

            lf_bottom = lf_merged[lf_merged[decile_col].isin(bottom3)]
            lf_top = lf_merged[lf_merged[decile_col].isin(top3)]
            ylb = (lf_bottom[label_col] == positive_label).astype(int).to_numpy()
            slb = lf_bottom[score_col].to_numpy(dtype=float)
            ylt = (lf_top[label_col] == positive_label).astype(int).to_numpy()
            slt = lf_top[score_col].to_numpy(dtype=float)
            n_pos_bottom, n_pos_top = int(ylb.sum()), int(ylt.sum())
            if n_pos_bottom < 2 or n_pos_top < 2 or len(np.unique(ylb)) < 2 or len(np.unique(ylt)) < 2:
                lf["gap_top_minus_bottom"] = {
                    "status": "too small to evaluate (a decile group has <2 positives or single class)",
                    "n_positives_bottom3": n_pos_bottom, "n_positives_top3": n_pos_top,
                }
            else:
                lf_gap_point, lf_gap_lo, lf_gap_hi = bootstrap_gap_ci(
                    ylb, slb, ylt, slt, seed=seed + 4000, n_boot=n_boot,
                )
                lf["gap_top_minus_bottom"] = {
                    "point": lf_gap_point, "ci_low": lf_gap_lo, "ci_high": lf_gap_hi,
                    "n_positives_bottom3": n_pos_bottom, "n_positives_top3": n_pos_top,
                }
        result["leakage_free"] = lf

    return result
