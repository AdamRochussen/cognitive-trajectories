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

Writes data/tucker_drob_parameters.csv, results/tucker_drob.txt and
figures/fig3_tucker_drob_trajectories.(png|pdf|svg).
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
    # study,   domain, label,                     i,     i_se,  s,      s_se,  q,      q_se
    ("VCAP",   "Gs", "Perceptual speed",          0.198, 0.014, -0.315, 0.006, -0.039, 0.003),
    ("VCAP",   "Gf", "Abstract reasoning",        0.164, 0.015, -0.186, 0.007, -0.037, 0.003),
    ("VCAP",   "Gm", "Memory",                    0.184, 0.015, -0.172, 0.007, -0.039, 0.003),
    ("VCAP",   "Gv", "Visuospatial reasoning",    0.090, 0.015, -0.190, 0.007, -0.014, 0.003),
    ("VCAP",   "Gc", "Crystallised knowledge",    0.117, 0.010,  0.134, 0.004, -0.038, 0.002),
    ("BETULA", "Gs", "Perceptual speed",          0.925, 0.027, -0.936, 0.058, -0.850, 0.043),
    ("BETULA", "Gm", "Memory",                    0.577, 0.020, -0.744, 0.041, -0.682, 0.035),
    ("BETULA", "Gv", "Visuospatial reasoning",    0.544, 0.022, -0.871, 0.042, -0.525, 0.033),
    ("BETULA", "Gc", "Crystallised knowledge",    0.319, 0.016, -0.033, 0.030, -0.580, 0.028),
]
P = pd.DataFrame(PARAMS, columns=["study", "domain", "label", "i", "i_se",
                                  "s", "s_se", "q", "q_se"])
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
(ROOT / "results" / "tucker_drob.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), sharey=True)
for ax, study in zip(axes, ("VCAP", "BETULA")):
    cfg = STUDY[study]; lo, hi = cfg["ages"]
    ages = np.linspace(lo, hi, 400)
    ax.axhline(50, color="#999999", lw=0.7, ls=(0, (3, 3)), zorder=1)
    for _, r in P[P.study == study].iterrows():
        ax.plot(ages, T(r, ages), color=COL[r.domain], lw=2.0, zorder=4,
                label=f"{r.label} ({r.domain})  ({peak_label(r)})")
        v = vertex(r)
        if lo <= v <= hi:
            ax.plot([v], [T(r, v)], "o", ms=4.5, color=COL[r.domain], mec="white",
                    mew=0.9, zorder=6)
    ax.set_xlim(*XAXIS); ax.set_ylim(15, 57)
    ax.set_xticks(range(XAXIS[0], XAXIS[1] + 1, 10))
    ax.set_xlabel("Age (years)")
    ax.set_title(f"{cfg['title']}   n = {cfg['n']:,}, baseline ages {cfg['base']}",
                 fontsize=10, loc="left", weight="bold")
    ax.grid(axis="y", color="#E6E6E6", lw=0.6, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=7.6, loc="lower left", handlelength=2.0)
axes[0].set_ylabel("T score (50 at age 35; 10 points = one SD)")
fig.suptitle("Retest-adjusted longitudinal cognitive trajectories (Tucker-Drob et al. 2022)",
             fontsize=11.5, weight="bold", x=0.075, ha="left", y=1.035)
fig.text(0.075, 0.975, "Model-implied population-mean trajectories from the published "
         "growth-curve parameters (Tables 1\u20132), net of practice effects. "
         "SD = within-age SD at baseline. Dot = peak.",
         fontsize=7.8, color="#555555", ha="left")
fig.tight_layout()
save_all(fig, ROOT / "figures", "fig3_tucker_drob_trajectories")
plt.close(fig)
print("wrote figures/fig3_tucker_drob_trajectories (png/pdf/svg)")
