"""Sensitivity analyses for the primary trajectory models.

Writes results/sensitivity.txt. Each block varies one analytic choice and reports
the fitted peak age (and cumulative change since age 20 at age 70) per domain.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from model import fit, peak_age, traj
from extract import dat

ROOT = Path(__file__).resolve().parents[1]
DOMS = ("gf", "gc", "gq", "fs")
adult = dat[~dat.child].copy()


def line(df, dom, order=1, rho=0.70, robust=True):
    sub = df[df.domain == dom].reset_index(drop=True)
    cl = pd.factorize(sub["cluster"])[0] if robust else None
    r = fit(sub, order=order, rho=rho, clusters=cl)
    pk, lo, hi = peak_age(r)
    ci = f"[{lo:5.1f}, {hi:5.1f}]" if np.isfinite(lo) else " " * 14
    return f"{dom}: peak {pk:5.1f} {ci}  @70 {float(traj(70, r['beta'])):+.2f}"


out = []
p = out.append

p("1. ASSUMED TEST-RETEST CORRELATION r (sampling-variance formula)")
for rho in (0.50, 0.70, 0.85, 0.95):
    p(f"   r = {rho:.2f}   " + "   ".join(line(adult, d, rho=rho) for d in DOMS))

p("\n2. INFERENCE: cluster-robust (primary) vs Knapp-Hartung model-based SEs")
p("   cluster-robust   " + "   ".join(line(adult, d, robust=True) for d in DOMS))
p("   Knapp-Hartung    " + "   ".join(line(adult, d, robust=False) for d in DOMS))
p("   NB: with only 4-5 cohort clusters the sandwich estimator can understate")
p("   uncertainty; the gq interval in particular should be read from the KH row.")

p("\n3. QUADRATIC RATE FUNCTION (rate = a + b*age + c*age^2)")
for d in ("gf", "gc", "fs"):
    p("   " + line(adult, d, order=2))

p("\n4. LEAVE-ONE-SAMPLE-OUT")
for d in DOMS:
    sub = adult[adult.domain == d].reset_index(drop=True)
    pks = []
    for i in sub.index:
        s2 = sub.drop(index=i).reset_index(drop=True)
        if len(s2) < 3:
            continue
        r = fit(s2, 1, clusters=pd.factorize(s2["cluster"])[0])
        pks.append((sub.loc[i, "study"], peak_age(r)[0]))
    lo = min(pks, key=lambda x: x[1]); hi = max(pks, key=lambda x: x[1])
    p(f"   {d}: peak range {lo[1]:.1f} (without {lo[0]}) to {hi[1]:.1f} (without {hi[0]})")

p("\n5. CORLEY et al. (2023) FULL-SCALE: four-domain mean (-0.541, primary) vs -0.31 (Substack table)")
alt = adult.copy()
alt.loc[(alt.study == "Corley et al. (2023)") & (alt.domain == "fs"), "d"] = -0.31
p("   primary  " + line(adult, "fs"))
p("   -0.31    " + line(alt, "fs"))

p("\n6. ADDING THE THREE CHILDHOOD-BASELINE STUDIES (full-scale only)")
p("   adult only        " + line(adult, "fs"))
p("   + child baseline  " + line(dat, "fs"))

p("\n7. DROPPING THE ONE UNVERIFIED VALUE (Dam & Raven 2008, gf)")
p("   with     " + line(adult, "gf"))
p("   without  " + line(adult[adult.study != "Dam & Raven (2008)"], "gf"))

text = "\n".join(out)
(ROOT / "results").mkdir(exist_ok=True)
(ROOT / "results" / "sensitivity.txt").write_text(text + "\n")
print(text)
