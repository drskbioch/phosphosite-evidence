"""Contact sheet of the Rev 6 figures: one figure per page at full width (study/analysis/figures_contact_sheet.pdf), and a greyscale render of each main figure
(study/analysis/figures_greyscale.pdf). Reads the 600 dpi PNG files in study/manuscript/figures/.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
CUR = Path(__file__).resolve().parents[2]
FIG = CUR / "manuscript" / "figures"
MAIN = [f"Figure{i}" for i in range(1, 6)]
SUPP = [f"FigureS{i}" for i in range(1, 9)]
TITLE = {"Figure1": "Figure 1. Study design", "Figure2": "Figure 2. What a confidence tier counts", "Figure3": "Figure 3. The signals under every adjustment",
         "Figure4": "Figure 4. What the signals depend on", "Figure5": "Figure 5. Predictor audit (secondary)", "FigureS1": "Figure S1. Data quality and provenance",
         "FigureS2": "Figure S2. Study-effort index", "FigureS3": "Figure S3. Covariate balance", "FigureS4": "Figure S4. Inference", "FigureS5": "Figure S5. Held-out partitions",
         "FigureS6": "Figure S6. Robustness", "FigureS7": "Figure S7. Lookup motif, descriptive only", "FigureS8": "Figure S8. Disorder by tier and source set"}


def pages(names, out, grey=False):
    with PdfPages(out) as pdf:
        for n in names:
            im = Image.open(FIG / f"{n}.png").convert("RGB")
            if grey:
                im = im.convert("L").convert("RGB")
            w, h = im.size
            im = im.resize((2400, int(2400 * h / w)), Image.LANCZOS)
            PW = 8.27
            iw = PW - 0.6
            ih = iw * im.size[1] / im.size[0]
            fig = plt.figure(figsize=(PW, ih + 0.9))
            fig.text(0.3 / PW, 1 - 0.3 / (ih + 0.9), TITLE[n] + (" (greyscale)" if grey else ""), fontsize=10, fontweight="bold", va="top")
            ax = fig.add_axes([0.3 / PW, 0.3 / (ih + 0.9), iw / PW, ih / (ih + 0.9)])
            ax.imshow(im)
            ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)


if __name__ == "__main__":
    pages(MAIN + SUPP, CUR / "analysis" / "figures_contact_sheet.pdf")
    pages(MAIN, CUR / "analysis" / "figures_greyscale.pdf", grey=True)
    print("contact sheet and greyscale render written")
