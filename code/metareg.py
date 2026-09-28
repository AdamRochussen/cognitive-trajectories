"""Minimal but faithful REML random-effects meta-regression (metafor rma.uni equivalent).

Includes: iterative REML tau^2, Knapp-Hartung (t) inference, Q_E, I^2, R^2.
"""
import numpy as np
from scipy import stats


def reml_tau2(y, v, X, tol=1e-10, maxiter=500):
    k, p = X.shape
    tau2 = max(0.0, np.var(y, ddof=1) - np.mean(v))
    for _ in range(maxiter):
        w = 1.0 / (v + tau2)
        W = np.diag(w)
        XtWX = X.T @ W @ X
        P = W - W @ X @ np.linalg.solve(XtWX, X.T @ W)
        num = (y @ P @ P @ y) - np.trace(P)
        den = np.trace(P @ P)
        new = tau2 + num / den
        new = max(0.0, new)
        if abs(new - tau2) < tol:
            tau2 = new
            break
        tau2 = new
    return tau2


def rma(y, v, X, names, knapp=True):
    y = np.asarray(y, float); v = np.asarray(v, float); X = np.asarray(X, float)
    k, p = X.shape
    tau2 = reml_tau2(y, v, X)
    w = 1.0 / (v + tau2)
    W = np.diag(w)
    XtWX = X.T @ W @ X
    XtWXi = np.linalg.inv(XtWX)
    beta = XtWXi @ (X.T @ W @ y)
    resid = y - X @ beta

    if knapp and k > p:
        s2 = float(resid @ W @ resid) / (k - p)
        s2 = max(s2, 1.0)          # metafor's default: do not shrink below 1
        vb = s2 * XtWXi
        df = k - p
        dist = stats.t(df)
    else:
        vb = XtWXi
        df = np.inf
        dist = stats.norm

    se = np.sqrt(np.diag(vb))
    stat = beta / se
    pval = 2 * dist.sf(np.abs(stat))
    crit = dist.ppf(0.975)
    ci_lo, ci_hi = beta - crit * se, beta + crit * se

    # Q_E (residual heterogeneity) using fixed weights 1/v
    wf = 1.0 / v
    Wf = np.diag(wf)
    bf = np.linalg.solve(X.T @ Wf @ X, X.T @ Wf @ y)
    rf = y - X @ bf
    QE = float(rf @ Wf @ rf)
    QE_df = k - p
    QE_p = stats.chi2.sf(QE, QE_df) if QE_df > 0 else np.nan

    # typical within-study variance (Higgins & Thompson)
    P = Wf - Wf @ X @ np.linalg.solve(X.T @ Wf @ X, X.T @ Wf)
    v_typ = (k - p) / np.trace(P)
    I2 = 100 * tau2 / (tau2 + v_typ) if (tau2 + v_typ) > 0 else 0.0

    return dict(k=k, p=p, tau2=tau2, tau=np.sqrt(tau2), beta=beta, se=se, stat=stat,
                pval=pval, ci_lo=ci_lo, ci_hi=ci_hi, names=names, QE=QE, QE_df=QE_df,
                QE_p=QE_p, I2=I2, df=df, v_typ=v_typ)


def report(res, title, tau2_null=None):
    print(f"\n=== {title} ===")
    print(f"k = {res['k']}   tau^2 = {res['tau2']:.4f}  (tau = {res['tau']:.3f})   I^2 = {res['I2']:.1f}%")
    if tau2_null is not None and tau2_null > 0:
        r2 = max(0.0, 100 * (tau2_null - res['tau2']) / tau2_null)
        print(f"R^2 (tau^2 explained vs intercept-only) = {r2:.1f}%")
    lab = "t" if np.isfinite(res['df']) else "z"
    print(f"{'term':<26}{'est':>9}{'se':>9}{'  '+lab:>8}{'p':>9}{'95% CI':>22}")
    for i, nm in enumerate(res['names']):
        ci = f"[{res['ci_lo'][i]:+.3f}, {res['ci_hi'][i]:+.3f}]"
        print(f"{nm:<26}{res['beta'][i]:>+9.4f}{res['se'][i]:>9.4f}{res['stat'][i]:>8.2f}{res['pval'][i]:>9.4f}{ci:>22}")
    print(f"Residual heterogeneity: Q({res['QE_df']}) = {res['QE']:.1f}, p = {res['QE_p']:.2e}")


def var_dz(dz, n, r):
    """Approx sampling variance of a standardised mean change (raw-SD standardised)."""
    return 2 * (1 - r) / n + dz ** 2 / (2 * n)
