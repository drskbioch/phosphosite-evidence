"""Matched comparison: coarsened exact matching (CEM, primary design) and
propensity-score matching (sensitivity design).

Both designs compare a treated group with a control group after holding a
study-effort covariate (SEI) constant. Defaults: 10 SEI deciles, a balance
requirement of |SMD| < 0.1, 1000 stratum-preserving bootstrap resamples,
and a propensity caliper of 0.2 SD of the logit.
"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .sei import sei_deciles

SMD_BALANCE_THRESHOLD = 0.1
BOOTSTRAP_N_RESAMPLES = 1000
PROPENSITY_CALIPER_SD_OF_LOGIT = 0.2


def standardized_mean_difference(values_a, values_b, weights_a=None, weights_b=None):
    """Standardized mean difference, weighted if weights are given."""
    if weights_a is not None:
        mean_a = np.average(values_a, weights=weights_a)
        var_a = np.average((values_a - mean_a) ** 2, weights=weights_a)
    else:
        mean_a, var_a = np.mean(values_a), np.var(values_a)
    if weights_b is not None:
        mean_b = np.average(values_b, weights=weights_b)
        var_b = np.average((values_b - mean_b) ** 2, weights=weights_b)
    else:
        mean_b, var_b = np.mean(values_b), np.var(values_b)
    pooled_sd = np.sqrt((var_a + var_b) / 2)
    if pooled_sd == 0:
        return 0.0
    return (mean_a - mean_b) / pooled_sd


def _weighted_effect(treated, control, statistic, treated_weight, control_weight):
    if statistic == "fisher_exact":
        a = (treated * treated_weight).sum()
        b = ((1 - treated) * treated_weight).sum()
        c = (control * control_weight).sum()
        d = ((1 - control) * control_weight).sum()
        ac, bc, cc, dc = (a + 0.5, b + 0.5, c + 0.5, d + 0.5) if min(a, b, c, d) <= 0 else (a, b, c, d)
        return (ac * dc) / (bc * cc)  # odds ratio
    elif statistic in ("proportions_ztest", "mannwhitneyu"):
        treated_mean = np.average(treated, weights=treated_weight)
        control_mean = np.average(control, weights=control_weight)
        return treated_mean - control_mean
    else:
        raise ValueError(statistic)


def _cem_strata_and_weights(df, group_col, treated_value, sei_col, extra_strata_cols, n_deciles=10):
    df = df.copy()
    decile, n_bins = sei_deciles(df[sei_col], n_deciles=n_deciles)
    df["_sei_decile"] = decile
    df["_treated"] = (df[group_col] == treated_value).astype(int)
    df["_stratum"] = list(zip(*([df[c] for c in extra_strata_cols] + [df["_sei_decile"]])))

    strata_counts = df.groupby("_stratum")["_treated"].agg(["sum", "count"])
    strata_counts["n_treat"] = strata_counts["sum"]
    strata_counts["n_control"] = strata_counts["count"] - strata_counts["sum"]
    valid_strata = strata_counts[(strata_counts["n_treat"] > 0) & (strata_counts["n_control"] > 0)].index
    n_dropped = len(strata_counts) - len(valid_strata)

    matched = df[df["_stratum"].isin(valid_strata)].copy()

    weights = []
    for s, sub in matched.groupby("_stratum"):
        n_t = (sub["_treated"] == 1).sum()
        n_c = (sub["_treated"] == 0).sum()
        w = np.where(sub["_treated"] == 1, 1.0, n_t / n_c)
        weights.append(pd.Series(w, index=sub.index))
    matched["_cem_weight"] = pd.concat(weights).reindex(matched.index)

    return matched, n_dropped, n_bins


def _stratified_bootstrap_ci(matched, indicator_col, statistic, seed, n_boot=BOOTSTRAP_N_RESAMPLES):
    rng = np.random.default_rng(seed)
    codes, uniques = pd.factorize(matched["_stratum"].apply(str))
    treated = matched["_treated"].to_numpy()
    indicator = matched[indicator_col].to_numpy(dtype=float)
    weight = matched["_cem_weight"].to_numpy(dtype=float)
    stratum_positions = [np.where(codes == c)[0] for c in range(len(uniques))]

    estimates = np.empty(n_boot)
    for b in range(n_boot):
        idx = np.concatenate([rng.choice(pos, size=len(pos), replace=True) for pos in stratum_positions])
        mask = treated[idx] == 1
        estimates[b] = _weighted_effect(
            indicator[idx][mask], indicator[idx][~mask], statistic,
            weight[idx][mask], weight[idx][~mask],
        )
    lo, hi = np.percentile(estimates, [2.5, 97.5])
    return float(lo), float(hi)


def matched_comparison(df, group_col, treated_value, indicator_col, statistic, sei_col,
                        extra_strata_cols=("residue_type",), n_deciles=10, seed=20260924,
                        n_boot=BOOTSTRAP_N_RESAMPLES):
    """Coarsened exact matching (CEM), primary design.

    df: rows already restricted to the two groups being compared, with
        group_col, sei_col, indicator_col, and every column named in
        extra_strata_cols present.
    group_col / treated_value: e.g. group_col="tier", treated_value="GOLD"
        (the other value present in df is the control).
    indicator_col: the per-row outcome (boolean or numeric).
    statistic: "fisher_exact" (odds ratio), "proportions_ztest" or
        "mannwhitneyu" (weighted difference in means).
    sei_col: the covariate column matched on.
    extra_strata_cols: additional exact-match columns besides the SEI
        decile; pass () for none.

    Returns a dict with point_estimate, ci_low, ci_high, smd, balance_pass,
    n_strata_used and n_strata_dropped. A design that fails the balance
    requirement is reported as failed, not re-tuned.
    """
    matched, n_dropped, n_bins = _cem_strata_and_weights(
        df, group_col, treated_value, sei_col, list(extra_strata_cols), n_deciles=n_deciles,
    )
    treated_mask = matched["_treated"] == 1
    point_est = _weighted_effect(
        matched.loc[treated_mask, indicator_col].to_numpy(dtype=float),
        matched.loc[~treated_mask, indicator_col].to_numpy(dtype=float),
        statistic,
        matched.loc[treated_mask, "_cem_weight"].to_numpy(dtype=float),
        matched.loc[~treated_mask, "_cem_weight"].to_numpy(dtype=float),
    )
    balance = standardized_mean_difference(
        matched.loc[treated_mask, sei_col].to_numpy(dtype=float),
        matched.loc[~treated_mask, sei_col].to_numpy(dtype=float),
        matched.loc[treated_mask, "_cem_weight"].to_numpy(dtype=float),
        matched.loc[~treated_mask, "_cem_weight"].to_numpy(dtype=float),
    )
    balance_pass = abs(balance) < SMD_BALANCE_THRESHOLD
    ci_lo, ci_hi = _stratified_bootstrap_ci(matched, indicator_col, statistic, seed, n_boot=n_boot)
    null = 1 if statistic == "fisher_exact" else 0
    return {
        "design": "CEM",
        "n_analysis": len(df),
        "n_strata_total": matched["_stratum"].nunique() + n_dropped,
        "n_strata_dropped": n_dropped,
        "n_strata_used": matched["_stratum"].nunique(),
        "n_deciles_actual": n_bins,
        "point_estimate": point_est,
        "ci_low": ci_lo,
        "ci_high": ci_hi,
        "smd": balance,
        "balance_pass": bool(balance_pass),
        "ci_excludes_null": bool((ci_lo > null) or (ci_hi < null)),
    }


def propensity_matched_comparison(df, group_col, treated_value, indicator_col, statistic,
                                   sei_col, formula_extra_terms="", seed=20260924,
                                   caliper_sd=PROPENSITY_CALIPER_SD_OF_LOGIT,
                                   n_boot=BOOTSTRAP_N_RESAMPLES):
    """Propensity-score matching (sensitivity design): 1:1 nearest
    neighbour without replacement, within a caliper on the logit.

    formula_extra_terms: additional patsy formula terms appended to the
        propensity model, e.g. "+ C(species) + C(residue_type)".
    """
    import statsmodels.formula.api as smf

    rng = np.random.default_rng(seed)
    work = df.copy().reset_index(drop=True)
    work["_treated"] = (work[group_col] == treated_value).astype(int)

    formula = f"_treated ~ {sei_col} {formula_extra_terms}".strip()
    model = smf.logit(formula, data=work).fit(disp=0, method="newton", maxiter=200, tol=1e-8)
    if not model.mle_retvals.get("converged", True):
        model = smf.logit(formula, data=work).fit(disp=0, method="bfgs", maxiter=500)
    prob = np.clip(model.predict(work), 1e-6, 1 - 1e-6)
    work["_logit"] = np.log(prob / (1 - prob))

    caliper = caliper_sd * work["_logit"].std()
    treated_idx = work.index[work["_treated"] == 1].to_numpy()
    control_idx = work.index[work["_treated"] == 0].to_numpy()
    treated_logit = work.loc[treated_idx, "_logit"].to_numpy()
    control_logit = work.loc[control_idx, "_logit"].to_numpy()

    tree = cKDTree(control_logit.reshape(-1, 1))
    k = min(50, len(control_idx))
    dists, nn_idx = tree.query(treated_logit.reshape(-1, 1), k=k)
    if k == 1:
        dists, nn_idx = dists.reshape(-1, 1), nn_idx.reshape(-1, 1)

    order = rng.permutation(len(treated_idx))
    used_control = np.zeros(len(control_idx), dtype=bool)
    pairs = []
    for ti in order:
        for rank in range(k):
            d = dists[ti, rank]
            ci = nn_idx[ti, rank]
            if d > caliper:
                break
            if not used_control[ci]:
                used_control[ci] = True
                pairs.append((treated_idx[ti], control_idx[ci]))
                break

    n_pairs = len(pairs)
    matched_treated = [t for t, c in pairs]
    matched_control = [c for t, c in pairs]
    matched = pd.concat([work.loc[matched_treated], work.loc[matched_control]], ignore_index=True)
    matched["_cem_weight"] = 1.0

    treated_mask = matched["_treated"] == 1
    point_est = _weighted_effect(
        matched.loc[treated_mask, indicator_col].to_numpy(dtype=float),
        matched.loc[~treated_mask, indicator_col].to_numpy(dtype=float),
        statistic, np.ones(treated_mask.sum()), np.ones((~treated_mask).sum()),
    )
    balance = standardized_mean_difference(
        matched.loc[treated_mask, sei_col].to_numpy(dtype=float),
        matched.loc[~treated_mask, sei_col].to_numpy(dtype=float),
    )
    balance_pass = abs(balance) < SMD_BALANCE_THRESHOLD

    null = 1 if statistic == "fisher_exact" else 0
    if n_pairs > 0:
        indicator = matched[indicator_col].to_numpy(dtype=float)
        gold_slice = np.arange(0, n_pairs)
        bronze_slice = np.arange(n_pairs, 2 * n_pairs)
        estimates = np.empty(n_boot)
        for b in range(n_boot):
            sample = rng.choice(n_pairs, size=n_pairs, replace=True)
            t_idx, c_idx = gold_slice[sample], bronze_slice[sample]
            estimates[b] = _weighted_effect(
                indicator[t_idx], indicator[c_idx], statistic,
                np.ones(n_pairs), np.ones(n_pairs),
            )
        ci_lo, ci_hi = (float(x) for x in np.percentile(estimates, [2.5, 97.5]))
    else:
        ci_lo, ci_hi = float("nan"), float("nan")

    return {
        "design": "propensity",
        "n_analysis": len(df),
        "n_treated_candidates": len(treated_idx),
        "n_matched_pairs": n_pairs,
        "match_rate": n_pairs / len(treated_idx) if len(treated_idx) else 0.0,
        "caliper": float(caliper),
        "point_estimate": point_est,
        "ci_low": ci_lo,
        "ci_high": ci_hi,
        "smd": balance,
        "balance_pass": bool(balance_pass),
        "ci_excludes_null": bool((ci_lo > null) or (ci_hi < null)) if n_pairs > 0 else False,
    }
