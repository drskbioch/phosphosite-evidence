"""Shared figure style and helpers.

Journal widths: single column 85 mm, double column 180 mm. One sans-serif
font (Arial), all text at least 7 pt at final size, bold panel letters at a
constant offset, one colour-blind-safe palette for tiers and a separate one for
species. Figures are written as vector PDF and SVG plus 600 dpi PNG and TIFF.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "study" / "earlier_results"
OUTDIR = ROOT / "study" / "manuscript" / "earlier_figures"

MM = 1 / 25.4
W_SINGLE = 85 * MM
W_DOUBLE = 180 * MM

FS = 7.0          # base text size (pt)
FS_TITLE = 7.5
FS_LETTER = 9.0

INK = "#1a1a1a"
GREY = "#6b6b6b"
LIGHT = "#d9d9d9"

TIER = {"GOLD": "#E69F00", "SILVER": "#0072B2", "BRONZE": "#CC79A7"}
SPECIES = {"human": "#332288", "mouse": "#44AA99", "rat": "#DDCC77"}
TOOL = {"DeepPhos": "#0077BB", "MusiteDeep": "#EE7733", "PhosF3C": "#009988"}
NEUTRAL = "#4d4d4d"


def load(name):
    with open(RESULTS / name) as f:
        return json.load(f)


def setup():
    matplotlib.rcParams.update({
        "font.family": "Arial",
        "font.size": FS,
        "axes.labelsize": FS,
        "axes.titlesize": FS_TITLE,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.titlepad": 4,
        "xtick.labelsize": FS,
        "ytick.labelsize": FS,
        "legend.fontsize": FS,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "lines.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.grid": False,
    })


def letter(ax, text, dx=-24, dy=6):
    """Bold panel letter at a constant offset (points) from the axes' top-left."""
    ax.annotate(text, xy=(0, 1), xycoords="axes fraction", xytext=(dx, dy),
                textcoords="offset points", fontsize=FS_LETTER, fontweight="bold",
                ha="left", va="bottom", annotation_clip=False)


def title(ax, text):
    ax.set_title(text, fontsize=FS_TITLE, fontweight="bold", loc="left", pad=4)


def fmt_int(v):
    return f"{int(round(v)):,}"


def save(fig, name):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / f"{name}.pdf")
    fig.savefig(OUTDIR / f"{name}.svg")
    fig.savefig(OUTDIR / f"{name}.png", dpi=600)
    with Image.open(OUTDIR / f"{name}.png") as im:
        im.convert("RGB").save(OUTDIR / f"{name}.tif", dpi=(600, 600), compression="tiff_lzw")
    plt.close(fig)


def hbar_stack(ax, parts, colors, labels, total, min_label=8):
    """One horizontal 100% stacked bar. parts: list of counts."""
    left = 0.0
    for p, c, lab in zip(parts, colors, labels):
        w = p / total * 100
        ax.barh(0, w, left=left, color=c, edgecolor="white", linewidth=0.6, height=0.6)
        left += w
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.7, 0.7)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Percent")


def lum(hex_):
    h = hex_.lstrip("#")
    r, g, b = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


def text_on(color):
    return "white" if lum(color) < 0.25 else INK
