"""Domain-coded longitudinal dataset, extracted from the uploaded primary sources.

Standardiser: BASELINE-wave SD throughout (the convention the Substack table states
in its footnote 3), so that differential decline in the SD does not inflate d.

Domains (following Corley et al. 2023's own partition of the LBC1936 battery):
  gf  = visuospatial / reasoning / speed / working memory  ("fluid")
  gc  = vocabulary, information, reading, knowledge        ("crystallised")
  gq  = arithmetic / numerical  (reported separately, in NEITHER curve)
  fs  = the study's own full-scale / total composite
"""
import numpy as np, pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# ----------------------------------------------------------------- OWENS
# Owens (1966) Table 9: component difference scores 1950-1919, normalised standard
# scores with norm SD = 1, for the 96 retested (D_T) and 31 not-retested (D_U).
# Owens 1953's N=127 figure = n-weighted average of the two.
w96, w31 = 96, 31
ow53 = {k: (w96*t + w31*u)/(w96+w31) for k, (t, u) in {
    "gc": (0.631, 0.679),    # Verbal component
    "gq": (-0.059, -0.151),  # Numerical component
    "gf": (0.150, -0.035),   # Reasoning component
    "fs": (0.583, 0.447),    # Total
}.items()}
# Owens (1966) Table 3: 1961 minus 1950, same metric (norm SD = 1).
ow66 = {"gc": -0.032, "gq": -0.203, "gf": 0.015, "fs": -0.047}

# ----------------------------------------------- SCHWARTZMAN ET AL. 1987 (Table 3)
# M-Test, N=259, WWII (age ~25) -> current (M age 64.7), normal time limits.
# (mean_t1, sd_t1, mean_t2)
schwartz = {
    "Picture Completion": (15.74, 2.92, 14.36, "gf"),
    "Picture Anomalies":  (14.09, 3.30, 11.54, "gf"),
    "Paper Formboard":    (13.09, 4.42,  9.70, "gf"),
    "Verbal Analogies":   (26.49, 6.34, 24.96, "gf"),
    "Tool Recognition":   (21.66, 4.67, 22.33, "gc"),
    "Mechanical Info":    (21.57, 6.06, 22.92, "gc"),
    "Vocabulary":         (18.08, 6.14, 21.08, "gc"),
    "Arithmetic":         (10.67, 3.10,  9.69, "gq"),
    "TOTAL":              (142.21, 27.86, 137.01, "fs"),
}
sch_d = {k: (v[2]-v[0])/v[1] for k, v in schwartz.items()}
sch = {}
for dom in ("gf", "gc", "gq"):
    vals = [sch_d[k] for k, v in schwartz.items() if v[3] == dom]
    sch[dom] = float(np.mean(vals))
sch["fs"] = sch_d["TOTAL"]

# ------------------------------------------------ MORTENSEN & KLEVEN 1993 (Table 1)
# N=141, born 1914, WAIS at 50 / 60 / 70. (sd_at_50, delta_60_50, delta_70_60)
mort = {
    "Information":         (3.53,  0.26,  0.11, "gc"),
    "Comprehension":       (3.05,  0.32, -0.67, "gc"),
    "Similarities":        (3.95, -0.57,  0.07, "gc"),
    "Vocabulary":         (10.12,  0.46, -2.04, "gc"),
    "Digit Span":          (1.61, -0.68, -0.09, "gf"),
    "Digit Symbol":       (10.94, -0.01, -4.86, "gf"),
    "Picture Completion":  (3.13, -0.26, -0.34, "gf"),
    "Block Design":        (7.87, -0.95, -3.80, "gf"),
    "Picture Arrangement": (5.85, -1.26, -0.79, "gf"),
    "Object Assembly":     (6.79, -0.66, -2.46, "gf"),
    "Arithmetic":          (2.99, -0.35, -0.43, "gq"),
    "FSIQ":               (13.55, -1.88, -3.57, "fs"),
}
mo1, mo2 = {}, {}
for dom in ("gf", "gc", "gq"):
    mo1[dom] = float(np.mean([v[1]/v[0] for v in mort.values() if v[3] == dom]))
    mo2[dom] = float(np.mean([v[2]/v[0] for v in mort.values() if v[3] == dom]))
mo1["fs"] = mort["FSIQ"][1]/mort["FSIQ"][0]
mo2["fs"] = mort["FSIQ"][2]/mort["FSIQ"][0]

# ------------------------------------------- DEARY, STARR & MacLENNAN 1998 (Table 1)
# HOPE study, N=387, age 75.1 -> 79.3. The ONLY test given at BOTH waves is the NART
# (115.8, SD 7.5 -> 113.7). Raven's and Logical Memory were wave-2 only ("NT").
hope_gc = (113.7 - 115.8) / 7.5           # = -0.280  -> matches the Substack table

# --------------------------------------------------------- LARSEN ET AL. 2008 (Table 1)
lar_gc = (116.52 - 107.16) / 22.26        # ACB verbal
lar_gq = (104.61 - 104.43) / 22.01        # ACB arithmetic

# ------------------------------------------------------- CORLEY ET AL. 2023 (Table 3)
# LGC slope means/SEs and intercept variances; SD change/yr = slope / SD(intercept).
YRS = 12.53
corley_raw = {   # domain: (slope_mean, slope_se, intercept_variance)
    "visuospatial": (-0.201, 0.059, 13.711),
    "speed":        (-0.413, 0.046, 21.971),
    "memory":       (-0.361, 0.066, 171.046),
    "verbal":       (-0.022, 0.010, 59.634),
}
cor = {}
for k, (m, se, ivar) in corley_raw.items():
    sd = np.sqrt(ivar)
    cor[k] = (m/sd*YRS, se/sd*YRS)        # (total z change, SE)

def combine(keys, rho=0.5):
    """Mean of correlated domain changes + SE assuming exchangeable correlation rho."""
    ds  = np.array([cor[k][0] for k in keys])
    ses = np.array([cor[k][1] for k in keys])
    m = len(keys)
    var = (ses**2).sum() + 2*rho*sum(ses[i]*ses[j] for i in range(m) for j in range(i+1, m))
    return float(ds.mean()), float(np.sqrt(var)/m)

cor_gf = combine(["visuospatial", "speed", "memory"])
cor_gc = (cor["verbal"][0], cor["verbal"][1])
cor_fs = combine(["visuospatial", "speed", "memory", "verbal"])

# ============================================================ ASSEMBLE LONG FORMAT
rows = []
def add(study, cluster, n, a1, a2, dom, d, se=None, child=False, note=""):
    rows.append(dict(study=study, cluster=cluster, n=n, age1=a1, age2=a2,
                     domain=dom, d=d, se=se, child=child, note=note))

for dom, d in ow53.items():
    add("Owens (1953)", "Owens/Iowa", 127, 19.0, 50.0, dom, d, note="Owens 1966 Tab.9")
for dom, d in ow66.items():
    add("Owens (1966)", "Owens/Iowa", 96, 50.0, 61.0, dom, d, note="Owens 1966 Tab.3")
for dom in ("gf", "gc", "gq", "fs"):
    add("Schwartzman et al. (1987)", "Concordia", 259, 25.0, 64.7, dom, sch[dom])
for dom in ("gf", "gc", "gq", "fs"):
    add("Mortensen & Kleven (1993) 50-60", "Glostrup", 141, 50.0, 60.0, dom, mo1[dom])
for dom in ("gf", "gc", "gq", "fs"):
    add("Mortensen & Kleven (1993) 60-70", "Glostrup", 141, 60.0, 70.0, dom, mo2[dom])
add("Deary et al. (1998) HOPE", "HOPE", 387, 75.1, 79.3, "gc", hope_gc, note="NART only")
add("Larsen et al. (2008)", "VES", 4384, 19.92, 38.35, "gc", lar_gc)
add("Larsen et al. (2008)", "VES", 4385, 19.92, 38.35, "gq", lar_gq)
add("Dam & Raven (2008)", "DamRaven", 99, 20.0, 70.0, "gf", -0.610, note="APM")
add("Corley et al. (2023)", "LBC1936", 1091, 69.5, 82.03, "gf", cor_gf[0], se=cor_gf[1])
add("Corley et al. (2023)", "LBC1936", 1091, 69.5, 82.03, "gc", cor_gc[0], se=cor_gc[1])
add("Corley et al. (2023)", "LBC1936", 1091, 69.5, 82.03, "fs", cor_fs[0], se=cor_fs[1])
# child-baseline (excluded from primary adult models)
add("Bradway & Thompson (1962)", "Bradway", 111, 13.6, 29.5, "fs", 0.690, child=True)
add("Deary et al. (2000)", "SMS1932", 101, 11.0, 77.0, "fs", 0.920, child=True)
add("Deary et al. (2004)", "SMS1932/47", 493, 11.0, 80.0, "fs", 1.080, child=True)

dat = pd.DataFrame(rows)
dat["interval"] = dat.age2 - dat.age1
dat["midage"]   = (dat.age1 + dat.age2) / 2
dat["per_dec"]  = dat.d / dat.interval * 10

def source_values_table():
    """Every primary-source number used, with its table reference."""
    r = []
    for k, (t, u) in {"Verbal": (0.631, 0.679), "Numerical": (-0.059, -0.151),
                      "Reasoning": (0.150, -0.035), "Total": (0.583, 0.447)}.items():
        r.append(("Owens (1966)", "Table 9", f"{k} component, 1950-1919, retested N=96", t))
        r.append(("Owens (1966)", "Table 9", f"{k} component, 1950-1919, not retested N=31", u))
    for k, v in {"Verbal": -0.032, "Numerical": -0.203, "Reasoning": 0.015, "Total": -0.047}.items():
        r.append(("Owens (1966)", "Table 3", f"{k} component, 1961-1950", v))
    for k, (m1, sd1, m2, _) in schwartz.items():
        r.append(("Schwartzman et al. (1987)", "Table 3", f"{k} WWII mean", m1))
        r.append(("Schwartzman et al. (1987)", "Table 3", f"{k} WWII SD", sd1))
        r.append(("Schwartzman et al. (1987)", "Table 3", f"{k} current mean (normal time)", m2))
    for k, (sd, d1, d2, _) in mort.items():
        r.append(("Mortensen & Kleven (1993)", "Table 1", f"{k} SD at age 50", sd))
        r.append(("Mortensen & Kleven (1993)", "Table 1", f"{k} mean difference 60-50", d1))
        r.append(("Mortensen & Kleven (1993)", "Table 1", f"{k} mean difference 70-60", d2))
    r += [("Deary, Starr & MacLennan (1998)", "Table 1", "NART wave 1 mean", 115.8),
          ("Deary, Starr & MacLennan (1998)", "Table 1", "NART wave 1 SD", 7.5),
          ("Deary, Starr & MacLennan (1998)", "Table 1", "NART wave 2 mean", 113.7),
          ("Larsen et al. (2008)", "Table 1", "ACB verbal t1 mean", 107.16),
          ("Larsen et al. (2008)", "Table 1", "ACB verbal t1 SD", 22.26),
          ("Larsen et al. (2008)", "Table 1", "ACB verbal t2 mean", 116.52),
          ("Larsen et al. (2008)", "Table 1", "ACB arithmetic t1 mean", 104.43),
          ("Larsen et al. (2008)", "Table 1", "ACB arithmetic t1 SD", 22.01),
          ("Larsen et al. (2008)", "Table 1", "ACB arithmetic t2 mean", 104.61),
          ("Dam & Raven (2008)", "secondary: Substack table", "APM standardised change 20->70", -0.610)]
    for k, (m, se, iv) in corley_raw.items():
        r.append(("Corley et al. (2023)", "Table 3", f"{k} slope mean", m))
        r.append(("Corley et al. (2023)", "Table 3", f"{k} slope SE", se))
        r.append(("Corley et al. (2023)", "Table 3", f"{k} intercept variance", iv))
    return pd.DataFrame(r, columns=["study", "source_table", "quantity", "value"])


if __name__ == "__main__":
    (ROOT / "data").mkdir(exist_ok=True)
    cols = ["study", "cluster", "n", "age1", "age2", "interval", "midage",
            "domain", "d", "se", "child", "note"]
    dat[cols].round(4).to_csv(ROOT / "data" / "effect_sizes.csv", index=False)
    source_values_table().to_csv(ROOT / "data" / "source_values.csv", index=False)
    print(f"wrote data/effect_sizes.csv ({len(dat)} rows) and data/source_values.csv")
    pd.set_option("display.width", 200)
    print(dat.pivot_table(index=["study", "age1", "age2"], columns="domain",
                          values="d").round(3).to_string())
