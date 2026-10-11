"""
Phase H1, Step 3, Task B/C: score parsing + pre-registered metrics.

Pre-registered spec (session brief + DECISIONS_HIGH_TIER_PUSH.md Step 3/G3):
  - Per-SEI-decile AUROC, bootstrap 95% CI, 1000 resamples.
  - Understudied-vs-well-studied gap: bottom 3 vs top 3 deciles, with CI.
  - Pooled AUROC and rank order: standard (raw pooled) vs SEI-stratified
    (mean of the 10 per-decile AUROCs -- an explicit, stated definition,
    since the brief names the comparison but not this exact statistic).
  - Leakage-free secondary: positives restricted to first_reported_year
    strictly after the model's own publication year; negatives unchanged.
    n>=200 positives threshold; below that, report "too small to evaluate".
    Also computes the bottom3-vs-top3 gap (with CI) on this leakage-free
    subset, per Task B.5 -- not just the pooled AUROC.

Gap direction, stated explicitly (close-out clarification, 2026-09-26):
SEI_decile is assigned by build_eval_sample.py via pd.qcut(SEI_primary,
q=10, labels=1..10) -- decile 1 = LOWEST SEI_primary = UNDERSTUDIED
proteins, decile 10 = HIGHEST SEI_primary = WELL-STUDIED proteins (qcut
labels increase with the underlying value). "bottom3" (deciles 1-3) is
therefore the understudied end; "top3" (deciles 8-10) is well-studied.
gap = AUROC(top3, well-studied) - AUROC(bottom3, understudied). A
POSITIVE gap means the tool scores BETTER on well-studied proteins than
on understudied ones -- the direction the confound thesis predicts.

Leakage-free comparability caveat (close-out clarification): the
leakage-free subset changes POSITIVE composition (restricted to sites
first-reported after the model's publication year) while NEGATIVES stay
the full fixed set -- so its pooled AUROC (and gap) is not directly
comparable in absolute terms to the primary analysis's, which uses all
positives. The comparison that IS valid is qualitative: does the
understudied-vs-well-studied pattern (direction, rough magnitude) persist
under leakage-free positives, not whether the two pooled-AUROC numbers
are the same scale.

Seed 20260924 throughout, per standing rule.
"""
import re
import json

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

SEED = 20260924
N_BOOT = 1000
SAMPLE_PATH = "study/processed_data/epsd/predictor_audit/eval_sample_with_sequence.parquet"  # relative to the project root
OUT_DIR = "study/processed_data/epsd/predictor_audit"  # relative to the project root

MODEL_PUB_YEAR = {
    "DeepPhos": 2019,     # Luo et al., Bioinformatics 2019
    "MusiteDeep": 2020,   # Wang et al., NAR 2020 (web-repo's own cited paper)
    "PhosF3C": 2024,      # bioRxiv preprint date 2024-12-25, conservative proxy (no stated cutoff)
}


def load_deepphos_scores():
    st = pd.read_csv(f"{OUT_DIR}/deepphos_scores_ST.txt", sep="\t", header=None,
                      names=["uniprot_id", "position", "score"], quotechar='"')
    y = pd.read_csv(f"{OUT_DIR}/deepphos_scores_Y.txt", sep="\t", header=None,
                     names=["uniprot_id", "position", "score"], quotechar='"')
    scores = pd.concat([st, y], ignore_index=True)
    scores["position"] = scores["position"].astype(int)
    return scores


def load_musitedeep_scores():
    """Windowed-input run (Task B final methodology session): each record in
    musitedeep_windowed_output/output_results.txt is a synthetic 33-mer keyed
    by synth_id ("uniprot_id|position"), fixed-padded so the true target site
    is always at local (1-indexed) position 17 -- see
    build_musitedeep_windows.py's docstring for why this is exact, not an
    approximation. Any other row for the same record is an incidental S/T/Y
    inside the same window and is dropped here."""
    path = f"{OUT_DIR}/musitedeep_windowed_output/output_results.txt"
    rows = []
    current_id = None
    score_re = re.compile(r":([0-9.]+)")
    with open(path, "r") as f:
        next(f)  # header line: Position Residue PTMscores Cutoff=0.5
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].strip()
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            pos = int(parts[0])
            if pos != 17:
                continue  # incidental site elsewhere in the 33-mer, not our target
            m = score_re.search(parts[2])
            if m is None:
                continue
            score = float(m.group(1))
            rows.append((current_id, score))
    hits = pd.DataFrame(rows, columns=["synth_id", "score"])
    key = pd.read_parquet(f"{OUT_DIR}/musitedeep_window_key.parquet")
    merged = key.merge(hits, on="synth_id", how="inner")
    return merged[["uniprot_id", "position", "score"]]


PHOSF3C_RESULT_PATH = ("study/tools/phosf3c/PhosF3C-main"  # relative to the project root
                        "/result/task_b/task_b_sample_Y.json")


def load_phosf3c_scores():
    """Windowed-input run, Y-type checkpoint only (S/T checkpoints incomplete
    this round -- see phaseH1_step3_predictor_audit.md). Each record in
    result/task_b/task_b_sample_Y.json is keyed by synth_id
    ("uniprot_id|position"); the true target's local offset within the
    windowed sub-sequence is recorded per-site in phosf3c_window_key.parquet
    (varies near sequence termini, unlike MusiteDeep's fixed-padding scheme)
    -- see build_phosf3c_windows.py's docstring for the equivalence proof and
    its empirical spot-check (scores matched a full-length-input control run
    to ~1e-6)."""
    with open(PHOSF3C_RESULT_PATH) as f:
        results = json.load(f)
    hits = pd.DataFrame([(r["id"], r["position"], r["output"]) for r in results],
                        columns=["synth_id", "local_offset", "score"])
    key = pd.read_parquet(f"{OUT_DIR}/phosf3c_window_key.parquet")
    merged = key.merge(hits, on="synth_id", how="inner")
    merged = merged[merged["local_offset"] == merged["expected_offset"]]
    return merged[["uniprot_id", "position", "score"]]


def bootstrap_auroc_ci(y_true, y_score, rng, n_boot=N_BOOT):
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    n = len(y_true)
    if len(np.unique(y_true)) < 2:
        return None, None, None
    point = roc_auc_score(y_true, y_score)
    boot_vals = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            boot_vals[i] = np.nan
        else:
            boot_vals[i] = roc_auc_score(yt, ys)
    valid = boot_vals[~np.isnan(boot_vals)]
    ci_low, ci_high = np.percentile(valid, [2.5, 97.5]) if len(valid) > 0 else (np.nan, np.nan)
    return point, ci_low, ci_high


def bootstrap_gap_ci(y_true_bottom, y_score_bottom, y_true_top, y_score_top, rng, n_boot=N_BOOT):
    yb, sb = np.asarray(y_true_bottom), np.asarray(y_score_bottom)
    yt, st_ = np.asarray(y_true_top), np.asarray(y_score_top)
    point = roc_auc_score(yt, st_) - roc_auc_score(yb, sb)
    diffs = np.empty(n_boot)
    nb, nt = len(yb), len(yt)
    for i in range(n_boot):
        idx_b = rng.integers(0, nb, nb)
        idx_t = rng.integers(0, nt, nt)
        yb_i, sb_i = yb[idx_b], sb[idx_b]
        yt_i, st_i = yt[idx_t], st_[idx_t]
        if len(np.unique(yb_i)) < 2 or len(np.unique(yt_i)) < 2:
            diffs[i] = np.nan
            continue
        diffs[i] = roc_auc_score(yt_i, st_i) - roc_auc_score(yb_i, sb_i)
    valid = diffs[~np.isnan(diffs)]
    ci_low, ci_high = np.percentile(valid, [2.5, 97.5])
    return point, ci_low, ci_high


def compute_tool_metrics(tool_name, merged, rng):
    """merged: sample rows with a non-null 'score' column for this tool."""
    y = (merged["label_type"] == "positive").astype(int).values
    s = merged["score"].values
    n_scored = len(merged)
    n_total = 10000
    print(f"\n=== {tool_name}: n_scored={n_scored}/{n_total} "
          f"({n_scored/n_total*100:.2f}%) ===")

    result = {"tool": tool_name, "n_scored": int(n_scored), "n_total": n_total}

    # --- per-decile AUROC ---
    decile_results = {}
    for decile in sorted(merged["SEI_decile"].unique()):
        sub = merged[merged["SEI_decile"] == decile]
        yy = (sub["label_type"] == "positive").astype(int).values
        ss = sub["score"].values
        point, lo, hi = bootstrap_auroc_ci(yy, ss, rng)
        decile_results[int(decile)] = {"n": len(sub), "auroc": point, "ci_low": lo, "ci_high": hi}
        print(f"  decile {decile}: n={len(sub)}, AUROC={point:.4f} [{lo:.4f},{hi:.4f}]" if point else
              f"  decile {decile}: n={len(sub)}, AUROC=undefined (single class)")
    result["per_decile"] = decile_results

    # --- understudied vs well-studied gap (bottom 3 vs top 3 deciles) ---
    deciles_sorted = sorted(merged["SEI_decile"].unique())
    bottom3 = deciles_sorted[:3]
    top3 = deciles_sorted[-3:]
    bottom_sub = merged[merged["SEI_decile"].isin(bottom3)]
    top_sub = merged[merged["SEI_decile"].isin(top3)]
    yb = (bottom_sub["label_type"] == "positive").astype(int).values
    sb = bottom_sub["score"].values
    yt = (top_sub["label_type"] == "positive").astype(int).values
    st_ = top_sub["score"].values
    gap_point, gap_lo, gap_hi = bootstrap_gap_ci(yb, sb, yt, st_, rng)
    result["gap_top_minus_bottom"] = {"point": gap_point, "ci_low": gap_lo, "ci_high": gap_hi,
                                       "bottom_deciles": bottom3, "top_deciles": top3,
                                       "n_bottom": len(bottom_sub), "n_top": len(top_sub)}
    print(f"  GAP (top3 - bottom3 deciles): {gap_point:.4f} [{gap_lo:.4f},{gap_hi:.4f}]")

    # --- pooled AUROC (standard) ---
    pooled_point, pooled_lo, pooled_hi = bootstrap_auroc_ci(y, s, rng)
    result["pooled_auroc"] = {"point": pooled_point, "ci_low": pooled_lo, "ci_high": pooled_hi}
    print(f"  Pooled AUROC (standard): {pooled_point:.4f} [{pooled_lo:.4f},{pooled_hi:.4f}]")

    # --- SEI-stratified summary metric: mean of per-decile AUROCs ---
    decile_aurocs = [v["auroc"] for v in decile_results.values() if v["auroc"] is not None]
    strat_mean = float(np.mean(decile_aurocs))
    result["sei_stratified_mean_auroc"] = strat_mean
    print(f"  SEI-stratified mean-of-deciles AUROC: {strat_mean:.4f}")

    # --- leakage-free secondary analysis ---
    pub_year = MODEL_PUB_YEAR[tool_name]
    positives = merged[merged["label_type"] == "positive"]
    negatives = merged[merged["label_type"] == "unlabelled"]
    lf_positives = positives[positives["first_reported_year"] > pub_year]
    n_lf_pos = len(lf_positives)
    result["leakage_free"] = {"model_pub_year": pub_year, "n_positives_after_cutoff": int(n_lf_pos)}
    print(f"  Leakage-free ({tool_name}, pub_year={pub_year}): "
          f"{n_lf_pos} positives with first_reported_year > {pub_year} "
          f"(of {len(positives)} total positives)")
    if n_lf_pos < 200:
        result["leakage_free"]["status"] = "too small to evaluate (n<200)"
        print(f"    -> too small to evaluate (n<200 threshold)")
    else:
        lf_merged = pd.concat([lf_positives, negatives], ignore_index=True)
        yl = (lf_merged["label_type"] == "positive").astype(int).values
        sl = lf_merged["score"].values
        lf_point, lf_lo, lf_hi = bootstrap_auroc_ci(yl, sl, rng)
        result["leakage_free"]["pooled_auroc"] = {"point": lf_point, "ci_low": lf_lo, "ci_high": lf_hi}
        result["leakage_free"]["status"] = "evaluated"
        result["leakage_free"]["comparability_note"] = (
            "Positive composition differs from the primary analysis (restricted to "
            "first_reported_year > pub_year); negatives are the same fixed set as "
            "primary. Pooled AUROC here is not directly comparable in absolute terms "
            "to the primary pooled AUROC -- only the qualitative pattern (direction, "
            "rough magnitude of the gap) is a valid comparison.")
        print(f"    Leakage-free pooled AUROC: {lf_point:.4f} [{lf_lo:.4f},{lf_hi:.4f}]")
        print(f"    vs primary pooled AUROC {pooled_point:.4f}: "
              f"{'similar' if abs(lf_point - pooled_point) < 0.02 else 'differs'} "
              f"(not directly comparable -- positive set differs, see comparability_note)")

        # --- leakage-free bottom3-vs-top3 gap, with per-decile n (Task B.5) ---
        lf_decile_n = {}
        for decile in bottom3 + top3:
            n_pos_d = len(lf_positives[lf_positives["SEI_decile"] == decile])
            lf_decile_n[int(decile)] = n_pos_d
            flag = " <-- FLAG: n<50 positives" if n_pos_d < 50 else ""
            print(f"    leakage-free decile {decile}: n_positives={n_pos_d}{flag}")

        lf_bottom = lf_merged[lf_merged["SEI_decile"].isin(bottom3)]
        lf_top = lf_merged[lf_merged["SEI_decile"].isin(top3)]
        ylb = (lf_bottom["label_type"] == "positive").astype(int).values
        slb = lf_bottom["score"].values
        ylt = (lf_top["label_type"] == "positive").astype(int).values
        slt = lf_top["score"].values
        n_pos_bottom = int((ylb == 1).sum())
        n_pos_top = int((ylt == 1).sum())
        if n_pos_bottom < 2 or n_pos_top < 2 or len(np.unique(ylb)) < 2 or len(np.unique(ylt)) < 2:
            result["leakage_free"]["gap_top_minus_bottom"] = {
                "status": "too small to evaluate (a decile group has <2 positives or single class)",
                "n_positives_bottom3": n_pos_bottom, "n_positives_top3": n_pos_top}
            print(f"    Leakage-free gap: too small to evaluate "
                  f"(n_pos bottom3={n_pos_bottom}, top3={n_pos_top})")
        else:
            lf_gap_point, lf_gap_lo, lf_gap_hi = bootstrap_gap_ci(ylb, slb, ylt, slt, rng)
            result["leakage_free"]["gap_top_minus_bottom"] = {
                "point": lf_gap_point, "ci_low": lf_gap_lo, "ci_high": lf_gap_hi,
                "n_positives_bottom3": n_pos_bottom, "n_positives_top3": n_pos_top,
                "per_decile_n_positives": lf_decile_n}
            print(f"    Leakage-free GAP (top3-bottom3): {lf_gap_point:.4f} "
                  f"[{lf_gap_lo:.4f},{lf_gap_hi:.4f}] "
                  f"(n_pos bottom3={n_pos_bottom}, top3={n_pos_top})")

    return result


def main():
    rng = np.random.default_rng(SEED)
    sample = pd.read_parquet(SAMPLE_PATH)
    print(f"Sample loaded: {len(sample):,} rows")

    all_results = {}

    # --- DeepPhos ---
    dp_scores = load_deepphos_scores()
    dp_merged = sample.merge(dp_scores, on=["uniprot_id", "position"], how="inner")
    n_missing = len(sample) - len(dp_merged)
    print(f"DeepPhos: {n_missing} sites could not be scored (missing from tool output)")
    all_results["DeepPhos"] = compute_tool_metrics("DeepPhos", dp_merged, rng)
    all_results["DeepPhos"]["n_unscoreable"] = int(n_missing)

    # --- MusiteDeep ---
    try:
        md_scores = load_musitedeep_scores()
        md_merged = sample.merge(md_scores, on=["uniprot_id", "position"], how="inner")
        n_missing = len(sample) - len(md_merged)
        print(f"\nMusiteDeep: {n_missing} sites could not be scored (missing from tool output)")
        all_results["MusiteDeep"] = compute_tool_metrics("MusiteDeep", md_merged, rng)
        all_results["MusiteDeep"]["n_unscoreable"] = int(n_missing)
    except FileNotFoundError:
        print("\nMusiteDeep output not yet available -- skipping this run.")

    # --- PhosF3C (Y-type checkpoint only) ---
    try:
        pf_scores = load_phosf3c_scores()
        pf_merged = sample.merge(pf_scores, on=["uniprot_id", "position"], how="inner")
        n_missing = len(sample) - len(pf_merged)
        print(f"\nPhosF3C: {n_missing} sites could not be scored "
              f"(Y-type checkpoint only this round -- S/T checkpoints incomplete)")
        all_results["PhosF3C"] = compute_tool_metrics("PhosF3C", pf_merged, rng)
        all_results["PhosF3C"]["n_unscoreable"] = int(n_missing)
        all_results["PhosF3C"]["scope_note"] = "Y-type sites only (1,359/10,000); S/T checkpoints incomplete"
    except FileNotFoundError:
        print("\nPhosF3C output not yet available -- skipping this run.")

    # --- rank order check (standard vs SEI-stratified) ---
    if len(all_results) >= 2:
        tools = list(all_results.keys())
        standard_rank = sorted(tools, key=lambda t: -all_results[t]["pooled_auroc"]["point"])
        strat_rank = sorted(tools, key=lambda t: -all_results[t]["sei_stratified_mean_auroc"])
        rank_changed = standard_rank != strat_rank
        print(f"\nStandard-eval rank order: {standard_rank}")
        print(f"SEI-stratified rank order: {strat_rank}")
        print(f"Rank order changed: {rank_changed}")
        all_results["_rank_order"] = {
            "standard": standard_rank, "sei_stratified": strat_rank, "changed": rank_changed}

    with open(f"{OUT_DIR}/task_b_metrics_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nWrote {OUT_DIR}/task_b_metrics_results.json")


if __name__ == "__main__":
    main()
