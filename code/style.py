"""Shared figure style for every figure in the repository.

All figures use DejaVu Sans (matplotlib's bundled font). PDFs embed it. SVGs keep
text as editable text, so each SVG is given a sans-serif fallback chain: without
it, a viewer that lacks DejaVu Sans (e.g. most macOS apps) substitutes its
default font, which is usually a serif.
"""
import re
from pathlib import Path
import matplotlib

RC = {
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.edgecolor": "#444444",
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
}
matplotlib.rcParams.update(RC)

SANS_STACK = "'DejaVu Sans', 'Helvetica Neue', Helvetica, Arial, sans-serif"


def fix_svg_fonts(path):
    p = Path(path)
    s = p.read_text()
    s = re.sub(r"font-family:\s*'DejaVu Sans'", f"font-family: {SANS_STACK}", s)
    p.write_text(s)


def save_all(fig, outdir, stem, dpi=300):
    outdir = Path(outdir)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(outdir / f"{stem}.{ext}", dpi=dpi, bbox_inches="tight",
                    facecolor="white")
    fix_svg_fonts(outdir / f"{stem}.svg")
