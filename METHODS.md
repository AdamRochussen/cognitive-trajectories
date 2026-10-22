# Methods

This document describes the original analyses relating to  *What age does intelligence
really peak?* published by Aporia and *At what age do we start to get stupider?* published by Fruit of Eden, by Adam Rochussen. There is enough detail here to replicate the analyses
without the code. Every input number is also listed, with its source table, in
`data/source_values.csv`.

There are three original analyses:

1. **A trajectory meta-regression of raw longitudinal studies** (Figures 1 and 2).
2. **Arithmetic on the published estimates of McArdle et al. (2002)** — the
   practice-effect comparison in the section "Practice effects can't be dismissed".
3. **A reconstruction of the retest-adjusted trajectories of Tucker-Drob et al.
   (2022)** from their published parameters (Figure 3).

---

## 1. Trajectory meta-regression

### 1.1 Design constraint

The analysis deliberately uses **only raw longitudinal data**: within-person
test–retest change with no adjustment for practice/retest effects. This is not
an endorsement of that choice. It tests whether the studies in the Heretical
Insights table support that article's conclusions *on the article's own
premises*. Consequences:

* Any per-occasion practice gain remains in the data, biasing every trajectory
  upward, and biasing short-interval studies more than long ones.
* Selective attrition is not corrected (except where authors used FIML — Corley
  et al., 2023). Attrition generally also biases trajectories upward.
* Cohort (Flynn) effects do not affect within-person change.

### 1.2 Study selection

The starting set was the 13 rows of the table in the Heretical Insights article.

**Included (9 samples, 7 independent cohorts):**

| Sample | Cohort | N | Ages | Numbers taken from |
|---|---|---|---|---|
| Owens (1953) | Iowa State | 127 | 19 → 50 | Owens (1966) Table 9 |
| Owens (1966) | Iowa State | 96 | 50 → 61 | Owens (1966) Table 3 |
| Schwartzman et al. (1987) | Canadian WWII veterans | 259 | 25 → 64.7 | Table 3 (normal time limits) |
| Mortensen & Kleven (1993), 50–60 | Glostrup 1914 | 141 | 50 → 60 | Table 1 |
| Mortensen & Kleven (1993), 60–70 | Glostrup 1914 | 141 | 60 → 70 | Table 1 |
| Deary et al. (1998) | HOPE | 387 | 75.1 → 79.3 | Deary, Starr & MacLennan (1998) Table 1 |
| Larsen et al. (2008) | Vietnam Experience Study | 4,384 | 19.9 → 38.4 | Table 1 |
| Dam & Raven (2008) | university sample | 99 | ~20 → ~70 | Heretical Insights table (not verified) |
| Corley et al. (2023) | LBC1936 | 1,091 | 69.5 → 82.0 | Table 3 |

**Excluded:**

* **Bradway & Thompson (1962)** (baseline 13.6), **Deary et al. (2000)** and
  **Deary et al. (2004)** (baseline 11): childhood baselines confound adult
  change with maturation. Bradway & Thompson also report age-normed IQ rather
  than raw scores. These three appear only in one sensitivity analysis.
* **Tucker-Drob et al. (2022)**: its estimates are adjusted for retest effects,
  which the design constraint excludes.

**Raw-data check.** Each included source was checked for any practice
adjustment. None applies one. Notes: Owens reports normalised standard scores
against a fixed 1,000-person norm (M = 5, SD = 1), so differences are already in
SD units; Mortensen & Kleven refer all three waves to a single fixed norm set, so
there is no age-norming drift; in the HOPE study only the NART was administered
at both waves (Raven's and Logical Memory were wave-2 only); Corley et al. use
FIML for missing data, which is not a retest adjustment.

### 1.3 Effect size

For each sample and domain:

```
d = (mean at follow-up − mean at baseline) / SD at baseline
```

The baseline-wave SD is used (the convention stated in footnote 3 of the
Heretical Insights article). Using pooled SDs instead changes no `d` by more
than 0.02. `d` is in within-sample SD units; samples differ in range
restriction, and no correction is applied.

### 1.4 Domain coding

Following the partition Corley et al. (2023) apply to the LBC1936 battery
(visuospatial, speed and memory as fluid; verbal as crystallised), extended to
the older batteries by content:

* **gf (fluid):** reasoning, visuospatial, processing speed, working memory.
* **gc (crystallised):** vocabulary, information, comprehension, word reading,
  mechanical knowledge.
* **gq (quantitative):** arithmetic/numerical tests — a separate curve, in
  neither gf nor gc.
* **fs (full-scale):** the study's own total or composite score.

| Study | gf | gc | gq | fs |
|---|---|---|---|---|
| Owens | Reasoning component | Verbal component | Numerical component | Total |
| Schwartzman | Picture Completion, Picture Anomalies, Paper Formboard, Verbal Analogies | Tool Recognition, Mechanical Information, Vocabulary | Arithmetic | Total M score |
| Mortensen & Kleven | Digit Span, Digit Symbol, Picture Completion, Block Design, Picture Arrangement, Object Assembly | Information, Comprehension, Similarities, Vocabulary | Arithmetic | FSIQ |
| HOPE | — | NART | — | — |
| Larsen | — | ACB verbal | ACB arithmetic | — |
| Dam & Raven | APM | — | — | — |
| Corley | visuospatial, speed, memory | verbal | — | mean of all four |

Where a domain has several tests, its `d` is the **unweighted mean** of the
test-level `d` values.

### 1.5 Extraction, worked

**Owens (1953).** Owens (1966) Table 9 gives 1950 − 1919 differences separately
for the 96 men later retested (D_T) and the 31 not retested (D_U). The N = 127
value is `(96·D_T + 31·D_U)/127`: Verbal +0.643, Numerical −0.081, Reasoning
+0.105, Total +0.550. (Owens states the total gain for the 127 as .55σ, which
this reproduces.)

**Owens (1966).** Table 3, 1961 − 1950: Verbal −0.032, Numerical −0.203,
Reasoning +0.015, Total −0.047.

**Schwartzman et al. (1987).** Table 3, `(current normal-time mean − WWII mean)
/ WWII SD`:

| Test | d | domain |
|---|---|---|
| Picture Completion | −0.473 | gf |
| Picture Anomalies | −0.773 | gf |
| Paper Formboard | −0.767 | gf |
| Verbal Analogies | −0.241 | gf |
| Tool Recognition | +0.143 | gc |
| Mechanical Information | +0.223 | gc |
| Vocabulary | +0.489 | gc |
| Arithmetic | −0.316 | gq |
| Total | −0.187 | fs |

Domain means: gf −0.563, gc +0.285.

**Mortensen & Kleven (1993).** Table 1 gives age-50 SDs and the mean
differences 60−50 and 70−60 for each raw subtest; each difference is divided by
the age-50 SD. Domain means: 50→60 gf −0.157, gc +0.020, gq −0.117,
fs −0.139 (FSIQ −1.88 / 13.55); 60→70 gf −0.265, gc −0.093, gq −0.144,
fs −0.263 (−3.57 / 13.55).

**HOPE.** NART 115.8 (SD 7.5) → 113.7: d = −0.280.

**Larsen et al. (2008).** Verbal (116.52 − 107.16)/22.26 = +0.420; arithmetic
(104.61 − 104.43)/22.01 = +0.008. No absolute full-scale estimate exists.

**Corley et al. (2023).** Table 3 gives each latent domain's slope mean, slope
SE and intercept variance. SD change per year = slope / √(intercept variance);
total change = that × 12.53 years:

| Domain | total d (SE) |
|---|---|
| Visuospatial | −0.680 (0.200) |
| Processing speed | −1.104 (0.123) |
| Memory | −0.346 (0.063) |
| Verbal | −0.036 (0.016) |

gf = mean of the first three = −0.710; gc = −0.036; fs = mean of all four =
−0.541. Composite SEs assume correlation 0.5 between domain slopes. Table 3
reports no general-factor slope, so fs is the four-domain mean; the Heretical
Insights table's −0.31 could not be derived from Table 3 (see sensitivity 5).

### 1.6 Sampling variances

Corley et al.'s SEs are used directly. For all other samples:

```
v(d) = 2(1 − r)/n + d²/(2n),    r = 0.70
```

`r` is the assumed test–retest correlation (reported values in these samples
span roughly .54–.91). Results are insensitive to it (sensitivity 1).

### 1.7 Model

Each sample reports change over its own age window, from 3 to 50 years long.
Each `d` is therefore treated as the **definite integral** of an underlying
rate-of-change function over that window:

```
d_i = ∫[a1_i, a2_i] rate(a) da
```

With a linear rate `rate(a) = α + β·a`, this is linear in the parameters with
**no intercept**:

```
d_i = α·L_i + β·Q_i,     L_i = a2_i − a1_i,     Q_i = (a2_i² − a1_i²)/2
```

**Estimation.** Random-effects meta-regression: τ² by REML (iterative
estimating equation, truncated at zero); coefficients by weighted least squares
with weights `1/(v_i + τ²)`. Primary inference uses **cluster-robust
(sandwich) standard errors clustered on cohort** (so the two Owens samples and
the two Mortensen waves are not treated as independent), with small-sample
correction `m/(m−1)·(k−1)/(k−p)` and `m − p` degrees of freedom (m clusters,
k samples, p = 2 parameters). Knapp–Hartung model-based intervals are reported
as a sensitivity.

**Trajectory.** The antiderivative of the fitted rate, anchored at age 20:

```
T(a) = α(a − 20) + β(a² − 400)/2
```

95% bands by the delta method with gradient `[(a − 20), (a² − 400)/2]` and
t critical values on the model's df. Bands are zero at age 20 by construction:
they describe uncertainty in change since 20.

**Peak age.** `a* = −α/β`, the age at which the fitted rate crosses zero;
CI by the delta method on the ratio.

### 1.8 Results

| Domain | k (cohorts) | Peak (95% CI) | Change since 20, at 50 | at 70 | at 80 |
|---|---|---|---|---|---|
| Fluid | 7 (5) | 35.2 (30.1–40.3) | +0.01 | −0.58 | −1.05 |
| Quantitative | 6 (4) | 29.9 (29.4–30.4) | −0.08 | −0.39 | −0.63 |
| Full-scale | 6 (4) | 43.9 (34.2–53.5) | +0.39 | −0.08 | −0.53 |
| Crystallised | 8 (6) | 56.7 (42.5–70.8) | +0.52 | +0.46 | +0.32 |

Full output: `results/primary_models.txt`.

### 1.9 Sensitivity analyses

Full output: `results/sensitivity.txt`.

1. **Assumed r = 0.50 / 0.85 / 0.95:** no peak moves by more than 0.5 years.
2. **Knapp–Hartung instead of cluster-robust intervals:** peaks identical;
   intervals become gf 27.9–42.6, gq 26.8–33.0, fs 39.8–48.0, gc 46.0–67.3.
   With only 4–6 clusters the sandwich estimator can understate uncertainty,
   and the very narrow cluster-robust gq interval should not be taken at face
   value.
3. **Quadratic rate function:** gf 35.2, gc 52.0, fs 42.5; the quadratic term
   is not identified at these k.
4. **Leave-one-sample-out:** gf 32.7–36.3, gc 51.4–59.8, gq 28.4–30.0,
   fs 40.9–47.0.
5. **Corley full-scale = −0.31 (Heretical Insights value):** fs peak
   43.9 → 44.6.
6. **Adding the three childhood-baseline studies to full-scale:** fs peak
   43.9 → 52.7 (44.0–61.3).
7. **Dropping the unverified Dam & Raven value:** gf peak 35.2 → 35.6.

### 1.10 Figures

* **Figure 1** plots `T(a)` with 95% bands for each domain, ages 20–85.
* **Figure 2** plots, for each domain, each sample as a line segment from
  `(a1, T(a1))` to `(a2, T(a1) + d)` — i.e., its observed change drawn from the
  fitted curve at its baseline age — over the fitted curve and band.

### 1.11 Limitations

* Few samples per curve (6–8) and substantial between-study heterogeneity.
* Within-sample standardisers differ in range restriction.
* A linear rate function is a simplification; true decline is probably
  accelerating in late life.
* One value is unverified (Dam & Raven) and one reconstructed (Corley
  full-scale).
* The raw-data constraint leaves practice and attrition biases in place; both
  push the trajectories upward.

---

## 2. McArdle et al. (2002) practice-effect arithmetic

All inputs are published values; nothing is re-estimated.
Script: `code/mcardle_practice.py`; output: `results/mcardle_practice.txt`.

* **Retest intervals** (Table 1, all ages): 20.5% under one year; 48.5% under
  two years; median 2.2 years. Adult median lags: 3.0 years (ages 20–49) and
  3.1 years (ages 50–95).
* **Model fit** (text, p. 130): dual-exponential model without practice,
  χ²(5) = 1,838; with practice, χ²(7) = 2,070 (both relative to the same
  baseline; the implied gain is 232 on 2 df). The bilinear model gives a
  directly stated comparison with versus without practice, χ²(2) = 319 (p. 129).
* **Size of the practice effect for BCA:** practice mean 2.7 W (Table 8);
  average rate of change ages 20–75, −0.3 W per year (Table 9).
  2.7 / 0.3 = **9 years** of average decline per retest.
* **Expected observed change for a typical adult:** 2.7 + (−0.3 × 3.0) =
  **+1.8 W** despite underlying decline.
* **Robustness:** the bilinear model (Table 7) gives 3.2 W practice against
  −0.39 W per year, i.e. 8.2 years.

---

## 3. Tucker-Drob et al. (2022) trajectory reconstruction

Script: `code/tucker_drob.py`; inputs: `data/tucker_drob_parameters.csv`;
output: `results/tucker_drob.txt`, `figures/fig3_tucker_drob_trajectories.*`.
Nothing is re-estimated: the curves are computed from the published latent
growth-curve means.

**Source model.** Tucker-Drob et al. fit multivariate age-based growth curves
(their Eq. 1), `Y[t] = i + A[t]·s + A[t]²·q + B[t]·r + e[t]`, in two cohorts:
VCAP (Virginia, N = 5,566, baseline ages 18–99) with age basis
`A = (age − 50)/10`, and Betula (Umeå, N = 3,320) with `A = (age − 50)/25`.
`r` is a retest effect (B = 0 at a person's first occasion, 1 thereafter).
Outcomes are scaled so that one unit is one SD of the age- and
age²-residualised baseline scores.

**Inputs.** The means of i, s and q for each domain: Table 1 (VCAP: Gs, Gf, Gm,
Gv, Gc) and Table 2 (Betula: Gs, Gm, Gv, Gc; Betula has no abstract-reasoning
measure). For Gc these are the second-step (constrained) estimates, as reported.

**Retest adjustment.** Because the retest term `r` is estimated jointly with
the growth terms, `i`, `s` and `q` describe age-related change net of practice.
The retest-adjusted mean trajectory is the model's prediction for a person who
has never been retested, i.e. with `B[t] = 0` — which is how the authors derive
their plotted trajectories (their Eq. 5b):

```
Y(age) = i + A·s + A²·q
T(age) = 50 + 10·[Y(age) − Y(35)]
```

This is the metric of their Figure 3 (T = 50 at age 35; 10 T points = one
within-age SD). Peak age is the vertex, `50 − scale·s/(2q)`, with scale 10
(VCAP) or 25 (Betula).

**Validation.** Computed T-scores at the right-hand end of each curve match the
black mean curves in the original Figure 3 as read by eye (VCAP at 95: Gs 24.1,
Gv 36.1, Gm 32.7, Gf 32.2, Gc 51.2; Betula at 85: Gs 17.7, Gm 24.2, Gv 24.2,
Gc 40.1).

**Plotted age ranges.** The original Figure 3 plots both cohorts from age 35
(to 95 for VCAP, 85 for Betula). Here each curve is drawn over its cohort's
observed baseline-age range, on a common 20–90 axis: VCAP 20–90 (baseline ages
18–99) and Betula 25–85 (baseline ages 25–86). Betula is not drawn below 25
because no participant was younger, so values there would be extrapolation
beyond the data rather than a model-implied mean. VCAP data are dense below
about 85, so its young end is well supported; Betula data are dense only from
about 40, so its curve below 40 (including the visuospatial peak at 29) is the
least certain part of the figure.

**Resulting peaks.** VCAP: abstract reasoning 24.9, memory 27.9, crystallised
67.6; perceptual speed and visuospatial reasoning have vertices below the youngest
observed age (9.6 and −17.9) and so decline throughout. Betula: visuospatial
29.3, memory 36.4, perceptual speed 36.2, crystallised 49.3. Peaks are marked
with a dot in Figure 3.

**Full-scale composites (Figure 3, panels C–D).** A full-scale score is a
weighted sum of domain scores, `C = Σ w_k Y_k`, so its mean trajectory is the
same weighted sum of the domain trajectories. Because each is quadratic, so is
the composite, and its peak is exact: `A* = −Σ w_k s_k / (2 Σ w_k q_k)`. The peak
age depends only on the relative weights; the correlations among domains affect
only how large the composite's changes are in its own SD units. Two weightings
are used (`data/tucker_drob_fullscale_weights.csv`):

* **WAIS-IV-like:** the subtest shares of the WAIS-IV FSIQ — 30% crystallised
  (VCI), 30% reasoning (PRI; split equally between Gf and Gv in VCAP, all Gv in
  Betula, which has no abstract-reasoning measure), 20% memory (WMI; episodic
  memory is the closest available domain, as working memory was not measured)
  and 20% processing speed (PSI).
* **g-loading:** proportional to each domain's standardised loading on the
  general level factor `Fi` (Tables 1–2), normalised to sum to 1. VCAP Gf's
  loading of 1.011 is capped at 1.

Resulting peaks: VCAP 33.6 (WAIS-IV-like) and 28.8 (g-loading); Betula 38.1
under both. Composites are plotted as the weighted mean of the domain T-scores,
i.e. in the same domain-SD units as panels A–B. A composite re-standardised to
its own SD would show somewhat larger changes, but the same peak ages. Giving
crystallised ability more weight moves the peak later: in VCAP by roughly five
years per extra 10% of weight (44 at a 50/50 crystallised/fluid split, as in
older Wechsler scales); in Betula by about 1.5 years per 10%.

**Caveats.** These are population-mean curves only; the original figure also
shows trajectories stratified by the general factor of fluid decline, which are
not reproduced. The full-scale composites are constructed here, not reported by
the original authors. Quadratic curves are least reliable at the edges of each
cohort's age range. The paper gives Betula's baseline age range as 25–86 in the
Results and 25–95 in the Methods; the figure uses 25–86.

---

## 4. Software

Python 3.12, NumPy 2.4, SciPy 1.17, pandas 3.0, Matplotlib 3.10. All figures
use DejaVu Sans (bundled with Matplotlib), embedded in the PDFs; the SVGs keep
text editable and declare a sans-serif fallback chain (`code/style.py`). The
meta-regression is implemented directly (`code/metareg.py`, `code/model.py`)
and follows the REML estimator used by `metafor::rma.uni`.

Data extraction, analysis code and this document were produced with the
assistance of Claude Opus 5.5 (Anthropic), working from the author's
specification. Every extracted input is listed with its source table in
`data/source_values.csv` so it can be checked against the original papers.
