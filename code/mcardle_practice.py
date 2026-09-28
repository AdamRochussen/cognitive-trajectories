"""McArdle, Ferrer-Caja, Hamagami & Woodcock (2002), Developmental Psychology 38,
115-142. Reproduces every number the article derives from this paper.

All inputs are transcribed from the published tables; nothing is re-estimated.
Writes results/mcardle_practice.txt.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Table 1: retest-interval distribution (all ages, N = 1,193)
PCT_LAG_UNDER_1Y = 20.5
PCT_LAG_1_TO_2Y = 28.0
MEDIAN_LAG_ALL_AGES = 2.2
MEDIAN_LAG_20_49 = 3.0
MEDIAN_LAG_50_95 = 3.1

# Broad Cognitive Ability (BCA), dual-exponential with-practice model
RATE_20_75 = -0.3        # Table 9: average rate of change, W units per year, ages 20-75
PRACTICE_NONLIN = 2.7    # Table 8: practice mean mu_p, W units
PEAK_BCA = 26.2          # Table 9
PEAK_GF, PEAK_GC = 22.8, 35.6

# Model-fit statistics quoted in the text (p. 130; all relative to a common baseline)
CHI2_NO_PRACTICE = (5, 1838)     # dual exponential without practice
CHI2_WITH_PRACTICE = (7, 2070)   # dual exponential with practice (= Table 8 BCA)
CHI2_BILINEAR_VS_NOPRACTICE = (2, 319)   # p. 129, bilinear with vs without practice

# Bilinear model (Table 7), for the robustness comparison
RATE_BILINEAR = -0.39
PRACTICE_BILINEAR = 3.2

lines = []
p = lines.append
p("McArdle et al. (2002): quantities quoted in the article\n")
p(f"Retest intervals (Table 1, all ages): {PCT_LAG_UNDER_1Y}% < 1 y; "
  f"{PCT_LAG_UNDER_1Y + PCT_LAG_1_TO_2Y:.1f}% < 2 y; median {MEDIAN_LAG_ALL_AGES} y")
p(f"Adult median lags: {MEDIAN_LAG_20_49} y (ages 20-49), {MEDIAN_LAG_50_95} y (ages 50-95)\n")
p(f"Peak ages (Table 9): BCA {PEAK_BCA}, Gf {PEAK_GF}, Gc {PEAK_GC}\n")

df0, c0 = CHI2_NO_PRACTICE; df1, c1 = CHI2_WITH_PRACTICE
p(f"Fit: without practice chi2({df0}) = {c0:,}; with practice chi2({df1}) = {c1:,}")
p(f"     implied gain from practice term: chi2({df1 - df0}) = {c1 - c0} (derived)")
p(f"     bilinear, with vs without practice: chi2({CHI2_BILINEAR_VS_NOPRACTICE[0]}) = "
  f"{CHI2_BILINEAR_VS_NOPRACTICE[1]} (stated, p. 129)\n")

yrs = PRACTICE_NONLIN / abs(RATE_20_75)
p(f"One retest = {PRACTICE_NONLIN} W / {abs(RATE_20_75)} W per year = {yrs:.1f} years of decline")
for lag, grp in ((MEDIAN_LAG_20_49, "20-49"), (MEDIAN_LAG_50_95, "50-95")):
    net = PRACTICE_NONLIN + RATE_20_75 * lag
    p(f"Expected observed change at median lag {lag} y (ages {grp}): "
      f"{PRACTICE_NONLIN} + ({RATE_20_75} x {lag}) = {net:+.2f} W")
p(f"\nRobustness, bilinear model (Table 7): {PRACTICE_BILINEAR} / {abs(RATE_BILINEAR)} = "
  f"{PRACTICE_BILINEAR / abs(RATE_BILINEAR):.1f} years of decline per retest")

text = "\n".join(lines)
(ROOT / "results").mkdir(exist_ok=True)
(ROOT / "results" / "mcardle_practice.txt").write_text(text + "\n")
print(text)
