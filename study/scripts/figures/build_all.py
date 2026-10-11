"""Build every Rev 6 figure: the composites (study/manuscript/figures/) and every panel as its own file (figures/panels/).

    python build_all.py                 composites, then all panels
    python build_all.py composites      composites only (also writes results/figure_captions.json)
    python build_all.py panels [ID ..]  all panels, or only the named ones (e.g. Fig3A FigS2G)

A panel is drawn by the same function as its composite, run in panel mode (lib.py): the selected panel goes into its own figure of the width below, the other
panels into a throwaway axes, and panel letters and captions are left out. Legends the composite shares across panels, and titles it lacks, are added in
panel mode. Single column is 85 mm, wide panels 180 mm; Arial 7 pt as in the composites.
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

import lib  # puts scripts/earlier_figures (the older figure scripts, which also have a fig3.py) first on the path ...
sys.path.insert(0, str(Path(__file__).resolve().parent))  # ... so the Rev 6 scripts must come before it
import fig1_2  # noqa: E402
import fig3  # noqa: E402
import fig4_5  # noqa: E402
import fig_supp  # noqa: E402
from lib import CAPTIONS, CAPTIONS_JSON, MODE, setup  # noqa: E402

assert Path(fig3.__file__).resolve().parent == Path(__file__).resolve().parent, fig3.__file__

MM = 1 / 25.4
# (figure name, script, function name, function, panel letters)
FIGURES = [("Figure1", "fig1_2.py", "figure1", fig1_2.figure1, [""]), ("Figure2", "fig1_2.py", "figure2", fig1_2.figure2, list("ABCDE")),
           ("Figure3", "fig3.py", "figure3", fig3.figure3, list("ABC")), ("Figure4", "fig4_5.py", "figure4", fig4_5.figure4, list("ABC")),
           ("Figure5", "fig4_5.py", "figure5", fig4_5.figure5, list("AB")), ("FigureS1", "fig_supp.py", "figS1", fig_supp.figS1, list("ABCDE")),
           ("FigureS2", "fig_supp.py", "figS2", fig_supp.figS2, list("ABCDEFGH")), ("FigureS3", "fig_supp.py", "figS3", fig_supp.figS3, ["A"]),
           ("FigureS4", "fig_supp.py", "figS4", fig_supp.figS4, list("ABCD")), ("FigureS5", "fig_supp.py", "figS5", fig_supp.figS5, list("ABCD")),
           ("FigureS6", "fig_supp.py", "figS6", fig_supp.figS6, list("ABC")), ("FigureS7", "fig_supp.py", "figS7", fig_supp.figS7, ["A"]),
           ("FigureS8", "fig_supp.py", "figS8", fig_supp.figS8, ["A"])]
# panel id -> (width mm, height mm)
# panels whose text is enlarged (factor on every text size; 7 pt base)
FONT_SCALE = {pid: 1.35 for pid in ("Fig3A", "Fig3B", "Fig3C", "Fig4B", "Fig4C", "FigS2G", "FigS3A", "FigS4D", "FigS6A", "FigS6B", "FigS6C", "FigS7A", "FigS8A")}
FONT_SCALE["Fig4A"] = 1.2
SIZES = {"Fig1": (180, 105), "Fig2A": (85, 80), "Fig2B": (85, 100), "Fig2C": (85, 65), "Fig2D": (180, 65), "Fig2E": (85, 75),
         "Fig3A": (180, 135), "Fig3B": (180, 135), "Fig3C": (180, 135), "Fig4A": (180, 150), "Fig4B": (180, 115), "Fig4C": (180, 100), "Fig5A": (180, 100), "Fig5B": (180, 80),
         "FigS1A": (85, 60), "FigS1B": (85, 65), "FigS1C": (85, 60), "FigS1D": (85, 60), "FigS1E": (85, 60),
         "FigS2A": (85, 60), "FigS2B": (85, 60), "FigS2C": (85, 60), "FigS2D": (85, 65), "FigS2E": (85, 55), "FigS2F": (85, 60), "FigS2G": (180, 105), "FigS2H": (85, 50),
         "FigS3A": (180, 100), "FigS4A": (85, 70), "FigS4B": (180, 85), "FigS4C": (85, 70), "FigS4D": (180, 160),
         "FigS5A": (85, 55), "FigS5B": (180, 45), "FigS5C": (180, 60), "FigS5D": (85, 50),
         "FigS6A": (180, 150), "FigS6B": (180, 150), "FigS6C": (180, 150), "FigS7A": (180, 85), "FigS8A": (180, 90)}


def panel_id(name, letter_):
    return name.replace("Figure", "Fig") + letter_


def composites():
    for name, _, _, fn, _ in FIGURES:
        fn()
        print(name, "written")
    CAPTIONS_JSON.write_text(json.dumps(CAPTIONS, indent=1), encoding="utf-8")


def panels(only=None):
    if not only and lib.PANEL_BIND.exists():
        lib.PANEL_BIND.unlink()
    for name, _, _, fn, letters in FIGURES:
        for L in letters:
            pid = panel_id(name, L)
            if only and pid not in only:
                continue
            w, h = SIZES[pid]
            MODE["letter"], MODE["pfig"], MODE["fscale"] = L, plt.figure(figsize=(w * MM, h * MM)), FONT_SCALE.get(pid, 1.0)
            try:
                fn()
            finally:
                MODE["letter"], MODE["pfig"], MODE["rect"], MODE["fscale"], MODE["notight"] = None, None, None, 1.0, False
                plt.close("all")
            print(pid, "written")


def manifest():
    """panels/PANEL_MANIFEST.md: panel, file, size, source, caption."""
    meta = json.loads(lib.PANEL_BIND.read_text(encoding="utf-8"))
    caps = json.loads(CAPTIONS_JSON.read_text(encoding="utf-8"))
    rows = []
    for name, script, fn_name, _, letters in FIGURES:
        for L in letters:
            pid = panel_id(name, L)
            m = meta[pid]
            c = caps[name]
            cap = " ".join(t for t in (c.get(""), c.get(L)) if t)
            cap = cap.replace("|", "/").replace(chr(10), " ")
            rows.append(f"| {pid} | `{pid}.png` | {m['width_mm']:g} x {m['height_mm']:g} | {"Arial 7 pt" if m.get("fscale", 1) == 1 else "Arial " + format(7 * m["fscale"], ".1f") + " pt (all text x " + format(m["fscale"], "g") + ")"} | `scripts/figures/{script}`, `{fn_name}()` (panel {L or 'whole figure'}) | {cap} |")
    head = ("# Panel files of Figures 1 to 5 and S1 to S8" + chr(10) + chr(10) +
            f"{len(rows)} panels, each as a 600 dpi PNG file, in this folder. The composite figures (PDF, SVG, PNG and TIFF) are in the folder above (`figures/`). "
            "Each panel file holds its axes, its legend (a legend that the composite shares across panels is repeated in every panel that uses it) and its title; the panel letter and the caption are left out. "
            "Width is 85 mm (single column) or 180 mm (wide panels); Arial 7 pt as in the composites, except the panels with a larger text size given in the table (all text scaled by one factor). Figure 1 and the figures with a single panel (S3, S7, S8) are one file each." + chr(10) + chr(10) +
            "Written by `scripts/figures/build_all.py`: each panel is drawn by the function of its composite, run in panel mode. The caption column is the caption text that belongs to the panel (the text shared by the whole figure, then the text of the panel). "
            "Sizes are those stored when the panel was drawn. Every panel is drawn from one of the composites, and the numbers in the SVG text of the composites are checked by `scripts/analysis/check_figure_bindings.py` against the bound values; the panels are not checked separately." + chr(10) + chr(10) +
            "| Panel | File | Size (mm, width x height) | Text size | Source script and function | Caption |" + chr(10) + "|---|---|---|---|---|---|" + chr(10))
    (lib.PANELDIR / "PANEL_MANIFEST.md").write_text(head + chr(10).join(rows) + chr(10), encoding="utf-8")
    print("manifest:", len(rows), "panels")


if __name__ == "__main__":
    setup()
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "composites"):
        composites()
    if what in ("all", "panels"):
        panels(set(sys.argv[2:]) or None)
        manifest()
