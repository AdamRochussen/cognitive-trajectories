"""Figures used in the article.

  figures/fig1_trajectories.(png|pdf|svg)      - fitted trajectories, four domains
  figures/fig2_study_segments.(png|pdf|svg)    - observed study segments vs fitted curves

PDF uses TrueType (fonttype 42) and SVG keeps text as text, so both are editable
in vector editors. Bands are 95% CIs from the primary (cluster-robust) models.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from model import fit, traj, traj_se, peak_age
from extract import dat
from style import fix_svg_fonts

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
# keep text as editable text in PDF/SVG for Illustrator
matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.edgecolor": "#444444",
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
})

COL = {"gf": "#D9541E", "gc": "#1B8A5A", "gq": "#8A6BBE", "fs": "#2A6FB5"}
LAB = {"gf": "Fluid (gf)", "gc": "Crystallised (gc)",
       "gq": "Quantitative (gq)", "fs": "Full-scale composite"}
ORDER = ["gf", "gq", "fs", "gc"]
AGES = np.arange(20, 85.01, 0.5)

adult = dat[~dat.child].copy()


def curves(df):
    out = {}
    for dom in ORDER:
        sub = df[df.domain == dom].reset_index(drop=True)
        r = fit(sub, 1, clusters=pd.factorize(sub["cluster"])[0])
        m = traj(AGES, r["beta"]); s = traj_se(AGES, r["beta"], r["vb"])
        tc = stats.t.ppf(0.975, r["df"])
        out[dom] = dict(m=m, lo=m-tc*s, hi=m+tc*s, peak=peak_age(r)[0],
                        k=r["k"], res=r, sub=sub)
    return out


def draw(ax, cv, bands=True, peaks=True):
    ax.axhline(0, color="#999999", lw=0.7, ls=(0, (3, 3)), zorder=1)
    for dom in ORDER:
        c = cv[dom]
        if bands:
            ax.fill_between(AGES, c["lo"], c["hi"], color=COL[dom], alpha=0.11,
                            lw=0, zorder=2)
    for dom in ORDER:
        c = cv[dom]
        ls = "--" if dom == "fs" else "-"
        ax.plot(AGES, c["m"], color=COL[dom], lw=2.0, ls=ls, zorder=4,
                label=f"{LAB[dom]}  (k={c['k']}, peak {c['peak']:.0f})")
        if peaks and 20 < c["peak"] < 85:
            pk = c["peak"]; pv = traj(np.array([pk]), c["res"]["beta"])[0]
            ax.plot([pk], [pv], "o", ms=4.5, color=COL[dom], mec="white",
                    mew=0.9, zorder=6)
    ax.set_xlim(20, 85); ax.set_ylim(-1.5, 1.0)
    ax.set_xticks(range(20, 86, 10))
    ax.set_yticks(np.arange(-1.5, 1.01, 0.5))
    ax.set_xlabel("Age (years)")
    ax.grid(axis="y", color="#E6E6E6", lw=0.6, zorder=0)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white")
    fix_svg_fonts(OUT / f"{stem}.svg")
    plt.close(fig)
    print("wrote", stem, "(png/pdf/svg)")




cv_main = curves(adult)
fig, ax = plt.subplots(figsize=(6.9, 4.6))
draw(ax, cv_main)
ax.set_ylabel("Cumulative change since age 20 (SD)")
ax.set_title("Adult cognitive trajectories from raw longitudinal studies",
             fontsize=11, loc="left", pad=14, weight="bold")
ax.text(0, 1.015, "Random-effects meta-regression of 9 adult-baseline samples from 7 "
        "independent cohorts; no retest adjustment", transform=ax.transAxes,
        fontsize=7.8, color="#555555")
ax.legend(frameon=False, fontsize=8, loc="lower left", handlelength=2.2)
ax.text(0.99, 0.02, "shaded = 95% CI", transform=ax.transAxes, fontsize=7,
        color="#888888", ha="right")
save(fig, "fig1_trajectories")

# ---------------------------------------------------------------- Figure 2
XL, XR, LABX = 18, 116, 88          # axis left, axis right, label column

SHORT = {
    "Owens (1953)": "Owens 1953",
    "Owens (1966)": "Owens 1966",
    "Schwartzman et al. (1987)": "Schwartzman 1987",
    "Mortensen & Kleven (1993) 50-60": "Mortensen 1993 (50\u201360)",
    "Mortensen & Kleven (1993) 60-70": "Mortensen 1993 (60\u201370)",
    "Deary et al. (1998) HOPE": "Deary 1998 (NART)",
    "Larsen et al. (2008)": "Larsen 2008",
    "Dam & Raven (2008)": "Dam & Raven 2008",
    "Corley et al. (2023)": "Corley 2023",
}


fig, axes = plt.subplots(2, 2, figsize=(11.8, 7.8), sharex=True, sharey=True)
for ax, dom in zip(axes.ravel(), ORDER):
    sub = adult[adult.domain == dom].reset_index(drop=True)
    r = fit(sub, 1, clusters=pd.factorize(sub["cluster"])[0])
    m = traj(AGES, r["beta"]); s = traj_se(AGES, r["beta"], r["vb"])
    tc = stats.t.ppf(0.975, r["df"]); pk = peak_age(r)[0]

    ax.axhline(0, color="#AAAAAA", lw=0.7, ls=(0, (3, 3)), zorder=1)
    ax.fill_between(AGES, m - tc*s, m + tc*s, color=COL[dom], alpha=0.13, lw=0, zorder=2)
    ax.plot(AGES, m, color=COL[dom], lw=2.2, zorder=3)
    if 20 < pk < 85:
        ax.plot([pk], [traj(np.array([pk]), r["beta"])[0]], "o", ms=5.5,
                color=COL[dom], mec="white", mew=1.1, zorder=7)

    ends = []
    for i in sub.index:
        a1, a2, d = sub.loc[i, "age1"], sub.loc[i, "age2"], sub.loc[i, "d"]
        y0 = traj(np.array([a1]), r["beta"])[0]
        ax.plot([a1, a2], [y0, y0 + d], color="#2E2E2E", lw=1.5, alpha=0.9,
                solid_capstyle="round", zorder=5)
        ax.plot([a1], [y0], "|", ms=7, color="#2E2E2E", mew=1.3, zorder=5)
        ax.plot([a2], [y0 + d], "o", ms=3.6, color="#2E2E2E", zorder=6)
        ends.append((a2, y0 + d, SHORT[sub.loc[i, "study"]]))

    # all labels in a right-hand column, spread vertically so none collide
    ends.sort(key=lambda e: e[1])
    ylab = [e[1] for e in ends]
    MIN, LO, HI = 0.145, -1.46, 0.96
    for _ in range(600):
        moved = False
        for j in range(len(ylab) - 1):
            gap = ylab[j + 1] - ylab[j]
            if gap < MIN:
                push = (MIN - gap) / 2 + 1e-4
                ylab[j] -= push; ylab[j + 1] += push; moved = True
        ylab = [min(max(y, LO), HI) for y in ylab]
        if not moved:
            break
    for (x, y, txt), yl in zip(ends, ylab):
        ax.plot([x + 1.2, LABX - 1.0], [y, yl], lw=0.5, color="#B5B5B5",
                zorder=4, solid_capstyle="butt")
        ax.text(LABX, yl, txt, fontsize=6.9, color="#3A3A3A",
                ha="left", va="center", zorder=8)

    ax.set_title(f"{LAB[dom]}      k = {r['k']}      peak {pk:.0f}",
                 fontsize=10.5, loc="left", weight="bold", color=COL[dom], pad=6)
    ax.set_xlim(XL, XR); ax.set_ylim(-1.55, 1.05)
    ax.set_xticks(range(20, 86, 10))
    ax.spines["bottom"].set_bounds(20, 85)
    ax.spines["left"].set_bounds(-1.5, 1.0)
    ax.grid(axis="y", color="#ECECEC", lw=0.6, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

for ax in axes[1]:
    ax.set_xlabel("Age (years)")
for ax in axes[:, 0]:
    ax.set_ylabel("Cumulative change since age 20 (SD)")

fig.suptitle("Observed study changes against the fitted trajectory",
             fontsize=13, weight="bold", x=0.052, ha="left", y=0.985)
fig.text(0.052, 0.950,
         "Each black segment spans one study's age window; its rise or fall is that study's observed change, "
         "anchored to the fitted curve at its baseline age.\n"
         "Adult-baseline raw longitudinal samples only \u2014 9 samples from 7 cohorts, no retest adjustment. "
         "Shaded band = 95% CI. Coloured dot = fitted peak.",
         fontsize=7.8, color="#555555", ha="left", va="top")
fig.tight_layout(rect=[0, 0, 1, 0.925])

for ext in ("png", "pdf", "svg"):
    fig.savefig(OUT / f"fig2_study_segments.{ext}", dpi=300,
                bbox_inches="tight", facecolor="white")
fix_svg_fonts(OUT / "fig2_study_segments.svg")
plt.close(fig)
print("wrote fig2_study_segments (png/pdf/svg)")
