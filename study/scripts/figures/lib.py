"""Shared style, loaders and binding log for the Rev 6 figures and tables.

Every number printed in a figure or table comes from a result file read with `provenance.jload`, which returns tracked numbers: formatting one writes
(text, file, key path, scale) to `provenance.LOG`. `finish()` stores the entries of one figure or table in results/figure_bindings.json together with the
tick-label strings of the figure (axis ticks are scales, not data); `check_bindings.py` re-reads the files and verifies every entry, and scans the SVG text.
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

CUR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CUR / "scripts" / "analysis"))
sys.path.insert(0, str(CUR / "scripts" / "earlier_figures"))
import provenance as P  # noqa: E402
from provenance import fmt_int, jload  # noqa: E402,F401
import figlib as F  # noqa: E402

Image.MAX_IMAGE_PIXELS = None
R6 = CUR / "analysis" / "results"
OUTDIR = CUR / "manuscript" / "earlier_figures"
FIGDIR6 = CUR / "manuscript" / "figures"
PANELDIR = FIGDIR6 / "panels"
PANEL_BIND = R6 / "panel_bindings.json"
CAPTIONS_JSON = R6 / "figure_captions.json"
TABDIR = CUR / "manuscript" / "tables"
BIND = R6 / "figure_bindings.json"
W = 504 / 72  # composite width in inches: exactly 504 pt (177.8 mm), the width limit of the journal
FS, FS_TITLE, INK, GREY, LIGHT = F.FS, F.FS_TITLE, F.INK, F.GREY, F.LIGHT
TIER = F.TIER
TIERS = ["GOLD", "SILVER", "BRONZE"]
# populations and designs: none of these is a tier colour
POP = {"P-all": "#1a1a1a", "P-human": "#5E4FA2", "P-abund": "#00796B"}
TOOL = {"DeepPhos": "#117733", "MusiteDeep": "#882255", "PhosF3C": "#6e6e6e"}
SIGNAL = {"S2": "Conservation", "S3": "Disorder", "S4": "Sequence motif", "S1": "Domain overlap", "S5": "Lookup motif"}
_cache = {}
TEXTLOG = []


def tx(cell, rel, key):
    """A table cell copied from a text source (claims row, SOURCES.md row, boolean field); logged so that check_bindings.py can re-verify it."""
    TEXTLOG.append((str(cell), rel, key))
    return str(cell)




def load(name):
    """Tracked Rev 6 result file (name relative to study/analysis/results)."""
    if name not in _cache:
        _cache[name] = jload(R6 / name)
    return _cache[name]


def load_dir(sub):
    """All JSON files of a result sub-folder as {file stem: tracked record}."""
    return {p.stem: load(f"{sub}/{p.name}") for p in sorted((R6 / sub).glob("*.json"))}


def setup():
    F.setup()
    plt.rcParams.update({"axes.titlesize": FS_TITLE, "legend.frameon": False})


# ---------------------------------------------------------------------------- panel mode
# A figure function draws every panel of its composite. In panel mode the same function is run once per panel: the selected panel is drawn into its own
# figure (MODE["pfig"], 85 or 180 mm wide), every other panel into a throwaway axes of a scratch figure, panel letters and the caption are left out, and
# `save` writes the panel files. Standalone-only additions (a legend that the composite shares, a title the composite lacks, row labels that the composite
# prints only in its first column) sit behind `standalone()`.
MODE = {"letter": None, "pfig": None, "name": None, "rect": None, "fscale": 1.0, "notight": False}
CAPTIONS = {}


def standalone():
    return MODE["letter"] is not None


def newfig(figsize):
    return plt.figure(figsize=figsize)


def pax(fig, spec, pid, multi=False):
    """Axes of panel `pid`. Composite: in `fig` at `spec` (a SubplotSpec or a [l, b, w, h] rectangle). Panel mode: the selected panel gets axes in the standalone
    figure (at the composite's grid position when `multi`, e.g. the cells of one panel), the others a throwaway axes."""
    if not standalone():
        return fig.add_axes(spec) if isinstance(spec, (list, tuple)) else fig.add_subplot(spec)
    if pid != MODE["letter"]:
        return fig.add_axes([0.1, 0.1, 0.8, 0.8])
    if multi:
        return MODE["pfig"].add_subplot(spec)
    return MODE["pfig"].add_subplot(111)


def letter(ax, t, dx=-26, dy=6):
    if standalone():
        return
    F.letter(ax, t, dx, dy)


def title(ax, t):
    F.title(ax, t)


# definitions of the two evidence sources behind the agreement score (columns source_PSP, source_PhosphoELM, source_dbPTM of the site files)
SRC_A = "Source set A (human file): PhosphoSitePlus and Phospho.ELM flags; agreement score 0 to 2."
SRC_B = "Source set B (all species, primary): adds the dbPTM flag; agreement score 0 to 3."
SRC_A_SHORT = "Source set A (human): PhosphoSitePlus and Phospho.ELM flags, score 0 to 2."
SRC_B_SHORT = "Source set B (all species): adds dbPTM, score 0 to 3."
SRC_B_FULL = "Source set B (all species, primary): PhosphoSitePlus, Phospho.ELM and dbPTM flags; agreement score 0 to 3."


def src_box(ax, xy, text, width, ha="left", va="top", fs=FS - 0.8):
    """Boxed explanation of the source sets, wrapped to `width` characters, placed in axes coordinates."""
    import textwrap
    ax.annotate(textwrap.fill(text, width), xy=xy, xycoords="axes fraction", fontsize=fs, ha=ha, va=va, linespacing=1.2, annotation_clip=False, zorder=6,
                bbox=dict(boxstyle="round,pad=0.35", fc="#f3f3f3", ec=GREY, lw=0.5))


class Rec:
    """Binding log of one figure or table: marks the start of LOG when created."""

    def __init__(self, name):
        self.name, self.start, self.tstart = name, len(P.LOG), len(TEXTLOG)
        MODE["name"] = name

    def finish(self, fig=None):
        ticks = []
        if fig is not None:
            fig.canvas.draw()
            allax = list(fig.axes) + [c for a in fig.axes for c in a.child_axes]
            for ax in allax:
                for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
                    ticks.append(t.get_text())
                for t in list(ax.get_xticklabels(minor=True)) + list(ax.get_yticklabels(minor=True)):
                    ticks.append(t.get_text())
        d = json.load(open(BIND, encoding="utf-8")) if BIND.exists() else {}
        d[self.name] = {"entries": [list(e) for e in P.LOG[self.start:]], "texts": [list(e) for e in TEXTLOG[self.tstart:]], "ticks": sorted(set(ticks))}
        BIND.write_text(json.dumps(d, indent=0), encoding="utf-8")


def scale_fonts(fig, k):
    """Multiply every text size of a panel figure by k (tick labels through tick_params, so that regenerated ticks keep the size)."""
    import matplotlib.text as mtext
    axes = list(fig.axes) + [c for a in fig.axes for c in a.child_axes]
    ticks = set()
    for ax in axes:
        for axis in (ax.xaxis, ax.yaxis):
            ticks.update(axis.get_ticklabels(which="both"))
            ticks.add(axis.get_offset_text())
    for t in fig.findobj(mtext.Text):
        if t not in ticks:
            t.set_fontsize(t.get_fontsize() * k)
    for ax in axes:
        for name, axis in (("x", ax.xaxis), ("y", ax.yaxis)):
            labs = axis.get_ticklabels()
            ax.tick_params(axis=name, which="both", labelsize=(labs[0].get_fontsize() if labs else FS) * k)


def save_rgb_png(fig, path):
    """600 dpi PNG in RGB mode (no alpha channel); the page is white, so the conversion changes no pixel."""
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=600)
    buf.seek(0)
    with Image.open(buf) as im:
        im.convert("RGB").save(path, dpi=(600, 600))


def write_all(fig, directory, name):
    """PDF, SVG (text as text), 600 dpi RGB PNG and LZW TIFF of one figure."""
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / f"{name}.pdf")
    fig.savefig(directory / f"{name}.svg")
    save_rgb_png(fig, directory / f"{name}.png")
    with Image.open(directory / f"{name}.png") as im:
        im.convert("RGB").save(directory / f"{name}.tif", dpi=(600, 600), compression="tiff_lzw")


def write_png(fig, directory, name):
    """600 dpi RGB PNG of one panel (panels are released as PNG only)."""
    directory.mkdir(parents=True, exist_ok=True)
    save_rgb_png(fig, directory / f"{name}.png")


def tick_strings(fig):
    fig.canvas.draw()
    ticks = []
    for ax in list(fig.axes) + [c for a in fig.axes for c in a.child_axes]:
        for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()) + list(ax.get_xticklabels(minor=True)) + list(ax.get_yticklabels(minor=True)):
            ticks.append(t.get_text())
    return sorted(set(ticks))


def save(fig, name, rec):
    """Composite: written to figures/. Panel mode: the selected panel is written to figures/panels/ as a PNG and its tick strings stored."""
    if standalone():
        pf, pid = MODE["pfig"], f"{name.replace('Figure', 'Fig')}{MODE['letter']}"
        if MODE["fscale"] != 1:
            scale_fonts(pf, MODE["fscale"])
        if not MODE["notight"]:
            pf.tight_layout(pad=0.3, rect=MODE["rect"]) if MODE["rect"] else pf.tight_layout(pad=0.3)
        d = json.load(open(PANEL_BIND, encoding="utf-8")) if PANEL_BIND.exists() else {}
        d[pid] = {"parent": name, "fscale": MODE["fscale"], "ticks": tick_strings(pf), "width_mm": round(pf.get_figwidth() * 25.4, 1), "height_mm": round(pf.get_figheight() * 25.4, 1)}
        PANEL_BIND.write_text(json.dumps(d, indent=0), encoding="utf-8")
        write_png(pf, PANELDIR, pid)
        plt.close("all")
        return
    rec.finish(fig)
    write_all(fig, FIGDIR6, name)
    plt.close(fig)


def legend_line(fig, text, y=0.012, size=FS - 0.5):
    """Wrapped note under the figure (n, intervals, formulas); every number in it is built from tracked values.
    `text` is a string, or a dict {"": shared text, "A": text of panel A, ...}: the composite prints the shared text and then "A: ..." for each panel; the
    dict is stored in results/figure_captions.json so that the panel manifest gives each panel its own caption."""
    if isinstance(text, dict):
        CAPTIONS[MODE["name"]] = text
        text = " ".join(([text[""]] if text.get("") else []) + [f"{k}: {v}" for k, v in text.items() if k])
    else:
        CAPTIONS[MODE["name"]] = {"": text}
    if standalone():
        return
    fig.text(0.012, y, wrap_to_figure(fig, text, size), fontsize=size, va="bottom", ha="left", color=INK)


def wrap_to_figure(fig, text, size, margin=0.012):
    """Line breaks so that no line is wider than the figure minus a margin on both sides (the wrapping of matplotlib runs to the figure edge)."""
    r = fig.canvas.get_renderer()
    limit = fig.get_figwidth() * fig.dpi * (1 - 2 * margin)
    t = fig.text(0, 0, "", fontsize=size)
    lines, cur = [], ""
    for word in text.split(" "):
        trial = (cur + " " + word).strip()
        t.set_text(trial)
        if cur and t.get_window_extent(r).width > limit:
            lines.append(cur)
            cur = word
        else:
            cur = trial
    lines.append(cur)
    t.remove()
    return chr(10).join(lines)


def pct(v, nd=1):
    """Tracked share (0 to 1) as a percentage string."""
    return f"{v * 100:.{nd}f}%"


def box_draw(ax, x, stats, color, width=0.5, lw=0.6, alpha=1.0, hatch=None):
    """Box with whiskers at the 5th and 95th percentile from stored summary statistics."""
    ax.plot([x, x], [stats["p5"], stats["p95"]], color=INK, lw=0.7, zorder=1)
    ax.add_artist(plt.Rectangle((x - width / 2, stats["q1"]), width, stats["q3"] - stats["q1"], fc=color, ec=INK, lw=lw, alpha=alpha, hatch=hatch, zorder=2))
    ax.plot([x - width / 2, x + width / 2], [stats["median"]] * 2, color=INK, lw=1.0, zorder=3)


def ci_marker(ax, y, est, ci, color=INK, marker="o", filled=True, ms=3.8, lw=0.8, cap=1.5, orient="h", zorder=3):
    face = color if filled else "white"
    if orient == "h":
        ax.errorbar(est, y, xerr=[[est - ci[0]], [ci[1] - est]], fmt=marker, ms=ms, mfc=face, mec=color, mew=0.8, ecolor=color, elinewidth=lw, capsize=cap, zorder=zorder)
    else:
        ax.errorbar(y, est, yerr=[[est - ci[0]], [ci[1] - est]], fmt=marker, ms=ms, mfc=face, mec=color, mew=0.8, ecolor=color, elinewidth=lw, capsize=cap, zorder=zorder)
