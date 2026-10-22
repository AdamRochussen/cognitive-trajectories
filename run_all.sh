#!/usr/bin/env bash
# Reproduce every original analysis in the article. Takes a few seconds.
set -euo pipefail
cd "$(dirname "$0")/code"
python3 extract.py          > /dev/null && echo "data/        effect sizes + source values"
python3 model.py            > /dev/null && echo "results/     primary_models.txt / .csv"
python3 sensitivity.py      > /dev/null && echo "results/     sensitivity.txt"
python3 mcardle_practice.py > /dev/null && echo "results/     mcardle_practice.txt"
python3 figures.py          > /dev/null && echo "figures/     fig1, fig2 (png/pdf/svg)"
python3 tucker_drob.py      > /dev/null && echo "figures/     fig3 (incl. full-scale panels) + results/tucker_drob.txt"
