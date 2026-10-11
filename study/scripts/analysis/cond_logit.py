"""Conditional logistic regression with one stratum per gene cluster (own implementation, numba).

Conditional likelihood of a stratum with n records, m of them with y = 1:
    L = exp(sum_{y=1} x_i b) / e_m(w_1..w_n),  w_i = exp(x_i b),  e_m = elementary symmetric polynomial of order m.
e_m and its first and second derivatives with respect to b are built by a dynamic programme over the records, with a rescaling
after each record so that large strata do not overflow. A weight per stratum (bootstrap multiplicity) multiplies its log-likelihood.
Newton-Raphson with step halving. Checked against statsmodels ConditionalLogit in check_cond_logit.py.
"""
import numpy as np
from numba import njit, prange


@njit(cache=True)
def _stratum(X, y, b, ll_out, g_out, H_out):
    n, p = X.shape
    m = 0
    for i in range(n):
        m += y[i]
    eta = X @ b
    # numerator
    num = 0.0
    sx = np.zeros(p)
    for i in range(n):
        if y[i] == 1:
            num += eta[i]
            for j in range(p):
                sx[j] += X[i, j]
    mx = eta.max()
    w = np.exp(eta - mx)  # common factor exp(mx) per record: e_m scales by exp(m*mx)
    E0 = np.zeros(m + 1)
    E1 = np.zeros((m + 1, p))
    E2 = np.zeros((m + 1, p, p))
    E0[0] = 1.0
    logscale = 0.0
    for i in range(n):
        top = min(i + 1, m)
        xi = X[i]
        wi = w[i]
        for k in range(top, 0, -1):
            e0 = E0[k - 1]
            for a in range(p):
                for c in range(p):
                    E2[k, a, c] += wi * (E2[k - 1, a, c] + xi[a] * E1[k - 1, c] + E1[k - 1, a] * xi[c] + xi[a] * xi[c] * e0)
            for a in range(p):
                E1[k, a] += wi * (E1[k - 1, a] + xi[a] * e0)
            E0[k] += wi * e0
        s = 0.0
        for k in range(top + 1):
            if E0[k] > s:
                s = E0[k]
        if s > 0:
            for k in range(top + 1):
                E0[k] /= s
                for a in range(p):
                    E1[k, a] /= s
                    for c in range(p):
                        E2[k, a, c] /= s
            logscale += np.log(s)
    e0 = E0[m]
    ll_out[0] = num - (np.log(e0) + logscale + m * mx)
    for a in range(p):
        g_out[a] = sx[a] - E1[m, a] / e0
        for c in range(p):
            H_out[a, c] = -(E2[m, a, c] / e0 - (E1[m, a] / e0) * (E1[m, c] / e0))


@njit(cache=True)
def loglik(X, y, starts, wts, b):
    """X, y sorted by stratum; starts[g]..starts[g+1] rows of stratum g. Returns total log-likelihood, gradient, Hessian."""
    p = X.shape[1]
    ll = 0.0
    g = np.zeros(p)
    H = np.zeros((p, p))
    ll_o = np.zeros(1)
    g_o = np.zeros(p)
    H_o = np.zeros((p, p))
    for s in range(len(starts) - 1):
        if wts[s] == 0:
            continue
        a, c = starts[s], starts[s + 1]
        ys = y[a:c]
        m = ys.sum()
        if m == 0 or m == c - a:
            continue
        _stratum(X[a:c], ys, b, ll_o, g_o, H_o)
        ll += wts[s] * ll_o[0]
        g += wts[s] * g_o
        H += wts[s] * H_o
    return ll, g, H


@njit(parallel=True, cache=True)
def loglik_par(X, y, starts, wts, b):
    """Same quantities as loglik, with the strata evaluated in parallel (each stratum is independent); used by the G21 bootstrap."""
    S = len(starts) - 1
    p = X.shape[1]
    ll_s = np.zeros(S)
    g_s = np.zeros((S, p))
    H_s = np.zeros((S, p, p))
    for s in prange(S):
        if wts[s] == 0:
            continue
        a, c = starts[s], starts[s + 1]
        ys = y[a:c]
        m = ys.sum()
        if m == 0 or m == c - a:
            continue
        ll_o = np.zeros(1)
        g_o = np.zeros(p)
        H_o = np.zeros((p, p))
        _stratum(X[a:c], ys, b, ll_o, g_o, H_o)
        ll_s[s] = wts[s] * ll_o[0]
        g_s[s] = wts[s] * g_o
        H_s[s] = wts[s] * H_o
    return ll_s.sum(), g_s.sum(0), H_s.sum(0)


def fit(X, y, starts, wts=None, maxiter=60, tol=1e-9, parallel=False, b0=None):
    """Returns (beta, converged, loglik, covariance from the Hessian). X must not contain a constant column.
    parallel: evaluate the strata in parallel (same likelihood). b0: starting values (default 0)."""
    loglik_fn = loglik_par if parallel else loglik
    p = X.shape[1]
    if wts is None:
        wts = np.ones(len(starts) - 1)
    b = np.zeros(p) if b0 is None else np.asarray(b0, float).copy()
    ll, g, H = loglik_fn(X, y, starts, wts, b)
    for _ in range(maxiter):
        try:
            step = np.linalg.solve(-H, g)
        except np.linalg.LinAlgError:
            return b, False, ll, None
        t = 1.0
        while True:
            nb = b + t * step
            nll, ng, nH = loglik_fn(X, y, starts, wts, nb)
            if np.isfinite(nll) and nll >= ll - 1e-12:
                break
            t /= 2
            if t < 1e-6:
                return b, False, ll, None
        b, ll, g, H = nb, nll, ng, nH
        if np.max(np.abs(t * step)) < tol:
            try:
                return b, True, ll, np.linalg.inv(-H)
            except np.linalg.LinAlgError:
                return b, True, ll, None
        if np.max(np.abs(b)) > 25:  # separation
            return b, False, ll, None
    return b, False, ll, None
