"""Trajectory meta-regression.

Each study reports the DEFINITE INTEGRAL of the underlying rate-of-change function
over its own age window:            d_i = INT_{a1}^{a2} rate(a) da
With rate(a) = a + b*a (+ c*a^2) this is exactly linear in the parameters, with
no-intercept design columns
        L = (a2-a1),   Q = (a2^2-a1^2)/2,   C = (a2^3-a1^3)/3
so the rate function can be fitted directly and then integrated to recover the
trajectory.  Note L*midage == Q, so the "change per decade vs midpoint age"
regression is the same model, just read as a derivative.
"""
import numpy as np, pandas as pd
from scipy import stats
from metareg import reml_tau2
from extract import dat

R_ASSUMED = 0.70
ANCHOR = 20.0


def var_d(d, n, r=R_ASSUMED):
    return 2*(1-r)/n + d**2/(2*n)


def design(a1, a2, order):
    cols = [a2-a1, (a2**2-a1**2)/2]
    if order == 2:
        cols.append((a2**3-a1**3)/3)
    return np.column_stack(cols)


def fit(sub, order=1, rho=0.70, clusters=None):
    y = sub["d"].values.astype(float)
    v = np.where(sub["se"].notna(), sub["se"].fillna(0).values**2,
                 var_d(y, sub["n"].values.astype(float), rho))
    X = design(sub["age1"].values.astype(float), sub["age2"].values.astype(float), order)
    k, p = X.shape
    tau2 = reml_tau2(y, v, X)
    w = 1.0/(v+tau2); W = np.diag(w)
    XtWXi = np.linalg.inv(X.T @ W @ X)
    beta = XtWXi @ (X.T @ W @ y)
    e = y - X @ beta
    s2 = max(float(e @ W @ e)/(k-p), 1.0)          # Knapp-Hartung
    vb = s2 * XtWXi
    if clusters is not None:                        # cluster-robust sandwich
        meat = np.zeros((p, p))
        for c in np.unique(clusters):
            idx = clusters == c
            s = (X[idx]*w[idx, None]).T @ e[idx]
            meat += np.outer(s, s)
        m = len(np.unique(clusters))
        adj = m/(m-1) * (k-1)/(k-p)
        vb = XtWXi @ (adj*meat) @ XtWXi
        df = m - p
    else:
        df = k - p
    return dict(beta=beta, vb=vb, df=max(df, 1), k=k, p=p, tau2=tau2, y=y, v=v, X=X)


def rate(a, beta):
    a = np.asarray(a, float)
    r = beta[0] + beta[1]*a
    if len(beta) == 3:
        r = r + beta[2]*a**2
    return r


def traj(a, beta, anchor=ANCHOR):
    """Cumulative change from `anchor` to a."""
    a = np.asarray(a, float)
    t = beta[0]*(a-anchor) + beta[1]*(a**2-anchor**2)/2
    if len(beta) == 3:
        t = t + beta[2]*(a**3-anchor**3)/3
    return t


def traj_se(a, beta, vb, anchor=ANCHOR):
    a = np.asarray(a, float)
    G = [a-anchor, (a**2-anchor**2)/2]
    if len(beta) == 3:
        G.append((a**3-anchor**3)/3)
    G = np.vstack(G)                     # p x len(a)
    return np.sqrt(np.einsum("pi,pq,qi->i", G, vb, G))


def peak_age(res):
    """Age at which rate crosses zero, + Fieller-style CI for the linear model."""
    b = res["beta"]
    if len(b) == 2:
        a0, a1 = b
        pk = -a0/a1
        # delta method on the ratio
        vb = res["vb"]
        g = np.array([-1/a1, a0/a1**2])
        se = np.sqrt(g @ vb @ g)
        t = stats.t.ppf(0.975, res["df"])
        return pk, pk-t*se, pk+t*se
    roots = np.roots([b[2], b[1], b[0]])
    roots = [r.real for r in roots if abs(r.imag) < 1e-9 and 15 < r.real < 95]
    return (roots[0] if roots else np.nan), np.nan, np.nan


def summarise(name, sub, order=1, clusters=None):
    res = fit(sub, order=order, clusters=clusters)
    b, vb = res["beta"], res["vb"]
    t = stats.t.ppf(0.975, res["df"])
    labs = ["alpha (rate at age 0)", "beta (d rate/d age)", "gamma (quadratic)"][:len(b)]
    print(f"\n{'='*78}\n{name}   [k = {res['k']}, tau = {np.sqrt(res['tau2']):.3f}]")
    for i, lab in enumerate(labs):
        se = np.sqrt(vb[i, i]); tv = b[i]/se
        print(f"  {lab:24s} {b[i]:+.5f}  SE {se:.5f}  t {tv:+.2f}  "
              f"p {2*stats.t.sf(abs(tv), res['df']):.4f}")
    pk, lo, hi = peak_age(res)
    if np.isfinite(pk):
        print(f"  PEAK AGE = {pk:.1f}" + (f"   95% CI [{lo:.1f}, {hi:.1f}]" if np.isfinite(lo) else ""))
    print("  rate (SD/decade):  " + "   ".join(
        f"age {A}: {10*rate(A, b):+.3f}" for A in (25, 40, 55, 70, 80)))
    print("  cumulative change from age 20: " + "   ".join(
        f"{A}: {traj(A, b):+.3f}" for A in (30, 40, 50, 60, 70, 80)))
    return res


def fit_all(order=1, include_child=False, rho=R_ASSUMED, cluster_robust=True):
    df = dat if include_child else dat[~dat.child]
    out = {}
    for dom in ("gf", "gc", "gq", "fs"):
        sub = df[df.domain == dom].reset_index(drop=True)
        cl = pd.factorize(sub["cluster"])[0] if cluster_robust else None
        out[dom] = (sub, fit(sub, order=order, rho=rho, clusters=cl))
    return out


if __name__ == "__main__":
    import sys
    from pathlib import Path
    ROOT = Path(__file__).resolve().parents[1]
    (ROOT / "results").mkdir(exist_ok=True)
    adult = dat[~dat.child]
    labels = [("gf", "FLUID (gf)"), ("gc", "CRYSTALLISED (gc)"),
              ("gq", "QUANTITATIVE (gq)"), ("fs", "FULL-SCALE (fs)")]
    rows = []
    with open(ROOT / "results" / "primary_models.txt", "w") as fh:
        _stdout = sys.stdout; sys.stdout = fh
        print("PRIMARY MODELS: adult baselines, raw longitudinal, linear rate, "
              "REML + cluster-robust SEs (cohort clusters)")
        for dom, label in labels:
            sub = adult[adult.domain == dom].reset_index(drop=True)
            res = summarise(label, sub, order=1, clusters=pd.factorize(sub["cluster"])[0])
            pk, lo, hi = peak_age(res)
            rows.append(dict(domain=dom, k=res["k"], n_cohorts=sub.cluster.nunique(),
                             alpha=res["beta"][0], beta=res["beta"][1],
                             tau=float(np.sqrt(res["tau2"])), peak=pk, peak_lo=lo, peak_hi=hi,
                             **{f"change_at_{a}": float(traj(a, res["beta"])) for a in (30, 40, 50, 60, 70, 80)}))
        sys.stdout = _stdout
    pd.DataFrame(rows).round(4).to_csv(ROOT / "results" / "primary_models.csv", index=False)
    print(open(ROOT / "results" / "primary_models.txt").read())
