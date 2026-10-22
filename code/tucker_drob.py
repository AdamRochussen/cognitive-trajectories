"""Model-implied mean trajectories from Tucker-Drob et al. (2022), Sci. Adv. 8, eabj2422.

These are the authors' retest-adjusted estimates, reconstructed from the latent
growth-curve means they publish; nothing is re-estimated.

Model (their Eq. 1):   Y[t] = i + A[t]*s + A[t]^2*q + B[t]*r + e[t]
r is a retest (practice) effect: B[t] = 0 at a person's first test and 1 at every
later test, estimated jointly with the growth terms. i, s and q therefore describe
age-related change net of practice. The retest-adjusted mean trajectory is the
curve for a person who has never been retested (their Eq. 5b):
                       Y(age) = i + A*s + A^2*q
  VCAP:   A = (age - 50) / 10     (decades,          Table 1)
  BETULA: A = (age - 50) / 25     (quarter centuries, Table 2)

Plotted on the authors' T metric (their Fig. 3): T = 50 + 10*[Y(age) - Y(35)].
One unit of Y is one SD of the age- and age^2-residualised baseline scores, so
10 T points = one within-age SD.

Full-scale composites (panels C-D). A weighted sum of the domain scores,
C = sum_k w_k Y_k, has mean trajectory sum_k w_k Y_k(age), itself a quadratic,
so its peak is exact: A* = -sum(w s) / (2 sum(w q)). The peak age depends only on
the relative weights. Two weightings are shown:
  WAIS-IV-like: Gc .30, reasoning .30 (Gf and Gv split in VCAP; Gv in Betula),
                memory .20 (episodic memory standing in for working memory),
                speed .20  - the subtest shares of the WAIS-IV FSIQ.
  g-loading:    proportional to each domain's standardised loading on the general
                level factor Fi (Tables 1-2; VCAP Gf's 1.011 capped at 1).
Composites are plotted as the weighted mean of the domain T-scores, i.e. in the
same domain-SD units as panels A-B. A composite re-standardised to its own SD
would show somewhat larger changes; peak ages are unaffected.

Writes data/tucker_drob_parameters.csv, data/tucker_drob_fullscale_weights.csv,
results/tucker_drob.txt and figures/fig3_tucker_drob_trajectories.(png|pdf|svg).
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from style import save_all

ROOT = Path(__file__).resolve().parents[1]
for d in ("data", "results", "figures"):
    (ROOT / d).mkdir(exist_ok=True)

# latent growth-curve means (i, s, q) and their SEs, transcribed from Tables 1-2
PARAMS = [
    # study,   domain, label,                     i,     i_se,  s,      s_se,  q,      q_se,   loading on Fi
    ("VCAP",   "Gs", "Perceptual speed",          0.198, 0.014, -0.315, 0.006, -0.039, 0.003, 0.582),
    ("VCAP",   "Gf", "Abstract reasoning",        0.164, 0.015, -0.186, 0.007, -0.037, 0.003, 1.011),
    ("VCAP",   "Gm", "Memory",                    0.184, 0.015, -0.172, 0.007, -0.039, 0.003, 0.724),
    ("VCAP",   "Gv", "Visuospatial reasoning",    0.090, 0.015, -0.190, 0.007, -0.014, 0.003, 0.817),
    ("VCAP",   "Gc", "Crystallised knowledge",    0.117, 0.010,  0.134, 0.004, -0.038, 0.002, 0.774),
    ("BETULA", "Gs", "Perceptual speed",          0.925, 0.027, -0.936, 0.058, -0.850, 0.043, 0.758),
    ("BETULA", "Gm", "Memory",                    0.577, 0.020, -0.744, 0.041, -0.682, 0.035, 0.711),
    ("BETULA", "Gv", "Visuospatial reasoning",    0.544, 0.022, -0.871, 0.042, -0.525, 0.033, 0.643),
    ("BETULA", "Gc", "Crystallised knowledge",    0.319, 0.016, -0.033, 0.030, -0.580, 0.028, 0.763),
]
P = pd.DataFrame(PARAMS, columns=["study", "domain", "label", "i", "i_se",
                                  "s", "s_se", "q", "q_se", "loading_Fi"])
P["age_scale_years"] = P.study.map({"VCAP": 10, "BETULA": 25})
P["source"] = P.study.map({"VCAP": "Table 1", "BETULA": "Table 2"})

STUDY = {
    # Each curve is drawn over its cohort's observed baseline-age range (VCAP 18-99,
    # Betula 25-86), trimmed to a common 20-90 axis. The original Fig. 3 plots both
    # from 35. Betula is not drawn below 25 because no participant was younger.
    "VCAP":   dict(ages=(20, 90), n=5566, base="18\u201399", title="A.  VCAP (Virginia, US)"),
    "BETULA": dict(ages=(25, 85), n=3320, base="25\u201386", title="B.  Betula (Ume\u00e5, Sweden)"),
}
XAXIS = (20, 90)
ORIGINAL_START = 35
ANCHOR = 35

# Gf and Gc use the fluid/crystallised colours of Figures 1-2
COL = {"Gs": "#A50F15", "Gf": "#D9541E", "Gm": "#F0A30A", "Gv": "#8C6D31", "Gc": "#1B8A5A"}


def y(row, age):
    A = (np.asarray(age, float) - 50.0) / row.age_scale_years
    return row.i + row.s * A + row.q * A ** 2


def T(row, age):
    return 50 + 10 * (y(row, age) - y(row, ANCHOR))


def vertex(row):
    return 50 - row.age_scale_years * row.s / (2 * row.q)


def peak_label(row):
    lo, hi = STUDY[row.study]["ages"]
    v = vertex(row)
    if v < lo:
        return "declines throughout"
    return f"peak {v:.0f}"



def fullscale_weights(study):
    sub = P[P.study == study].set_index("domain")
    doms = list(sub.index)
    reas = [k for k in ("Gf", "Gv") if k in doms]
    wais = {k: 0.0 for k in doms}
    wais.update({k: 0.30 / len(reas) for k in reas}); wais.update(Gc=0.30, Gm=0.20, Gs=0.20)
    lam = sub.loading_Fi.clip(upper=1.0)
    gload = (lam / lam.sum()).to_dict()
    return {"WAIS-IV-like weights": wais, "g-loading weights": gload}


def composite(study, w):
    sub = P[P.study == study].set_index("domain")
    i = sum(w[k] * sub.loc[k, "i"] for k in w)
    s_ = sum(w[k] * sub.loc[k, "s"] for k in w)
    q = sum(w[k] * sub.loc[k, "q"] for k in w)
    return pd.Series(dict(study=study, age_scale_years=sub.age_scale_years.iloc[0],
                          i=i, s=s_, q=q))


# ---------------------------------------------------------------- outputs
P.to_csv(ROOT / "data" / "tucker_drob_parameters.csv", index=False)

lines = ["Tucker-Drob et al. (2022): model-implied mean trajectories (retest-adjusted)",
         "T = 50 + 10*(Y(age) - Y(35)); vertex = 50 - scale*s/(2q)\n"]
for study in ("VCAP", "BETULA"):
    lo, hi = STUDY[study]["ages"]
    lines.append(f"{study}  (plotted {lo}-{hi}; the original Fig. 3 plots from {ORIGINAL_START})")
    for _, r in P[P.study == study].iterrows():
        top = 95 if study == "VCAP" else 85   # upper age plotted in the original Fig. 3
        ages = [a for a in (20, 25, 35, 50, 65, 75, 85, 95) if lo <= a <= top]
        lines.append(f"  {r.domain} {r.label:24s} vertex {vertex(r):6.1f}   " +
                     "  ".join(f"T@{a}={T(r, a):5.1f}" for a in ages))
    lines.append("")
lines.append("Check against the original Fig. 3 (black curves, read by eye): VCAP Gs ~25, "
             "Gv ~36, Gm ~33, Gf ~32, Gc ~51 at 95;")
lines.append("BETULA Gs ~18, Gm ~25, Gv ~24, Gc ~40 at 85.")
lines.append("\nFULL-SCALE COMPOSITES (weighted mean of domain scores)")
wrows = []
for study in ("VCAP", "BETULA"):
    for wname, w in fullscale_weights(study).items():
        c = composite(study, w)
        lines.append(f"  {study:6s} {wname:22s} peak {vertex(c):5.1f}   T@70 {T(c, 70):5.1f}   "
                     "weights " + ", ".join(f"{k} {v:.3f}" for k, v in w.items()))
        wrows += [dict(study=study, weighting=wname, domain=k, weight=round(v, 4)) for k, v in w.items()]
pd.DataFrame(wrows).to_csv(ROOT / "data" / "tucker_drob_fullscale_weights.csv", index=False)
(ROOT / "results" / "tucker_drob.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))

# ---------------------------------------------------------------- figure
FS_COL = "#2A6FB5"          # full-scale colour used in Figures 1-2
FS_LS = {"WAIS-IV-like weights": "-", "g-loading weights": (0, (5, 2.5))}
PANEL = {("VCAP", "fs"): "C.  VCAP full-scale composite",
         ("BETULA", "fs"): "D.  Betula full-scale composite"}

fig, axes = plt.subplots(2, 2, figsize=(11.2, 9.0), sharey=True)


def style_axis(ax, title):
    ax.axhline(50, color="#999999", lw=0.7, ls=(0, (3, 3)), zorder=1)
    ax.set_xlim(*XAXIS); ax.set_ylim(15, 57)
    ax.set_xticks(range(XAXIS[0], XAXIS[1] + 1, 10))
    ax.set_xlabel("Age (years)")
    ax.set_title(title, fontsize=10, loc="left", weight="bold")
    ax.grid(axis="y", color="#E6E6E6", lw=0.6, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=7.6, loc="lower left", handlelength=2.0)


for col, study in enumerate(("VCAP", "BETULA")):
    cfg = STUDY[study]; lo, hi = cfg["ages"]
    ages = np.linspace(lo, hi, 400)

    ax = axes[0, col]                                   # domain panels (unchanged)
    for _, r in P[P.study == study].iterrows():
        ax.plot(ages, T(r, ages), color=COL[r.domain], lw=2.0, zorder=4,
                label=f"{r.label} ({r.domain})  ({peak_label(r)})")
        v = vertex(r)
        if lo <= v <= hi:
            ax.plot([v], [T(r, v)], "o", ms=4.5, color=COL[r.domain], mec="white",
                    mew=0.9, zorder=6)
    style_axis(ax, f"{cfg['title']}   n = {cfg['n']:,}, baseline ages {cfg['base']}")

    ax = axes[1, col]                                   # full-scale panels
    for wname, w in fullscale_weights(study).items():
        c = composite(study, w); v = vertex(c)
        ax.plot(ages, T(c, ages), color=FS_COL, lw=2.0, ls=FS_LS[wname], zorder=4,
                label=f"Full-scale, {wname}  (peak {v:.0f})")
        if lo <= v <= hi:
            ax.plot([v], [T(c, v)], "o", ms=4.5, color=FS_COL, mec="white",
                    mew=0.9, zorder=6)
    style_axis(ax, PANEL[(study, "fs")])

for ax in axes[:, 0]:
    ax.set_ylabel("T score (50 at age 35; 10 points = one SD)")
fig.suptitle("Retest-adjusted longitudinal cognitive trajectories (Tucker-Drob et al. 2022)",
             fontsize=11.5, weight="bold", x=0.075, ha="left", y=1.02)
fig.text(0.075, 0.99, "Model-implied population-mean trajectories from the published "
         "growth-curve parameters (Tables 1\u20132), net of practice effects. "
         "SD = within-age SD at baseline. Dot = peak.\n"
         "Full-scale composites (C\u2013D) are weighted means of the domain scores: "
         "WAIS-IV-like weights (30% Gc, 30% reasoning, 20% memory, 20% speed) or weights "
         "proportional to each domain's g-loading.",
         fontsize=7.8, color="#555555", ha="left", va="top", linespacing=1.5)
fig.tight_layout(rect=[0, 0, 1, 0.955], h_pad=2.2)
save_all(fig, ROOT / "figures", "fig3_tucker_drob_trajectories")
plt.close(fig)
print("wrote figures/fig3_tucker_drob_trajectories (png/pdf/svg)")
