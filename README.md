# At what age does intelligence really peak? — data and code

Data, code and methods for the original analyses in *At what age does
intelligence really peak? A rebuttal to Aaron Dymarskiy* (Adam Rochussen), a
response to *Age and Cognitive Ability* (Heretical Insights).

## What is here

| Path | Contents |
|---|---|
| `METHODS.md` | Full methods, sufficient to replicate without the code |
| `data/source_values.csv` | Every number transcribed from a source paper, with its table |
| `data/effect_sizes.csv` | Derived standardised changes, one row per sample × domain |
| `code/extract.py` | Transcription and effect-size computation |
| `code/metareg.py`, `code/model.py` | Random-effects meta-regression and trajectory model |
| `code/sensitivity.py` | Sensitivity analyses |
| `code/mcardle_practice.py` | Arithmetic on McArdle et al. (2002)'s published estimates |
| `code/figures.py` | Figures 1 and 2 (PNG, plus PDF/SVG with editable text) |
| `code/style.py` | Shared figure style (font, SVG font fallback) |
| `code/tucker_drob.py` | Figure 3: Tucker-Drob et al. (2022) trajectories reconstructed from their published parameters |
| `data/tucker_drob_parameters.csv` | Growth-curve means transcribed from Tucker-Drob et al. Tables 1–2 |
| `results/` | Model output as produced by the code |
| `figures/` | Figures exactly as used in the article |

## Reproduce

```bash
pip install -r requirements.txt
./run_all.sh
```

Runs in a few seconds and regenerates everything in `data/`, `results/` and
`figures/`.

## Summary of the analysis

Nine adult-baseline samples (seven cohorts) from the Heretical Insights table,
using raw test–retest change only (no practice adjustment), were meta-analysed
with a random-effects model in which each study's observed change is the
integral of an age-dependent rate of change over its own age window. Integrating
the fitted rate gives the trajectory. Peak ages: quantitative 30, fluid 35,
full-scale 44, crystallised 57. Separately, the retest-adjusted trajectories
from Tucker-Drob et al. (2022) are reconstructed from their published
growth-curve parameters (Figure 3). See `METHODS.md` for details.

## Notes

* The Dam & Raven (2008) value is taken from the Heretical Insights table and has
  not been checked against the original.
* Code and documentation were prepared with the assistance of Claude Opus 5.5
  (Anthropic).
