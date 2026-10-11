"""Figure 1 (study design) and Figure 2 (what a confidence tier counts)."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from lib import (FS, FS_TITLE, GREY, INK, LIGHT, MODE, SRC_A, SRC_A_SHORT, SRC_B, SRC_B_SHORT, TIER, TIERS, W, Rec, box_draw, src_box, fmt_int, legend_line, letter, load, newfig, pax, save, setup, standalone, title)


def figure1():
    rec = Rec("Figure1")
    g23, gu, am = load("g23_tiers_summary.json"), load("gene_unit_summary.json"), load("amend_results.json")
    fig = newfig((W, 4.7))
    ax = pax(fig, [0.125, 0.11, 0.775, 0.77], "")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 66)
    ax.axis("off")
    nl = chr(10)

    # one colour per stage (fill, edge); none of these is a tier colour
    stages = [("#D9EEFA", "#1F6E9C"), ("#D3EEE6", "#00795A"), ("#F8E0D0", "#B04A00"), ("#E4E0F0", "#5E4FA2")]

    def stage_of(x):
        return 0 if x < 20 else 1 if x < 45 else 2 if x < 77 else 3

    def box(x, y, w, h, text, fc=None, ls="-", bold_first=False):
        fill, edge = stages[stage_of(x)]
        ax.add_artist(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.8", fc=fc or fill, ec=edge, lw=0.9, ls=ls))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=FS, linespacing=1.3)

    def arrow(x0, y0, x1, y1, ls="-"):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.7, color=INK, shrinkA=0, shrinkB=0, mutation_scale=7, ls=ls))

    for x, t in [(9.5, "Inputs"), (33, "Site table, covariates"), (60.5, "Designs"), (89, "Inference and signals")]:
        ax.text(x, 64.3, t, ha="center", va="center", fontsize=FS, fontweight="bold", color=stages[stage_of(x)][1])
    box(0, 46, 19, 16, "EPSD 2.0" + nl + "site records" + nl + "(human, mouse, rat)" + nl + "(PhosphoSitePlus," + nl + "Phospho.ELM, dbPTM" + nl + "source flags)")
    box(0, 34, 19, 9, "gene2pubmed links," + nl + "UniProt ID mapping")
    box(0, 22, 19, 9, "PaxDb 6.0" + nl + "protein abundance")
    box(0, 8, 19, 11, "Pharos" + nl + "target development" + nl + "level")
    box(23, 39, 20, 22, "Site table" + nl + f"{fmt_int(g23['n_records'])} sites" + nl + f"{fmt_int(gu['all']['accessions'])} proteins" + nl + nl + "agreement score," + nl + "study-cluster depth," + nl + "set-B tiers:" + nl + "GOLD, SILVER," + nl + "BRONZE")
    box(23, 26, 20, 10.5, "Check: " + f"{g23['n_nonprimary_pmids']} database or" + nl + "review PMIDs" + nl + "not counted", ls=(0, (3, 2)))
    box(23, 6, 20, 16, "Covariates" + nl + nl + "V50 study-effort" + nl + "index (per gene)," + nl + "abundance")
    arrow(19, 54, 23, 55)
    arrow(19, 38.5, 23, 47)
    arrow(19, 26.5, 23, 16)
    arrow(19, 13.5, 23, 10)
    arrow(33, 39, 33, 36.5, ls=(0, (3, 2)))
    box(47, 6, 28, 55, "GOLD against BRONZE" + nl + nl + "crude" + nl + "CEM on the V50 decile (primary)" + nl + "index x abundance deciles" + nl + "abundance deciles alone" + nl + "1:1 caliper matcher" + nl + "within-gene model (M1, M2)" + nl + "non-primary PMIDs removed" + nl + "tier-definition variants" + nl + "residue strata" + nl + "leave-one-cluster-out" + nl + "tier thresholds")
    arrow(43, 50, 47, 50)
    arrow(43, 14, 47, 14)
    nb = am["G17-0382"]["n_boot"]
    box(79, 44, 21, 17, "Gene-cluster bootstrap," + nl + f"{fmt_int(nb)} resamples" + nl + "(primary); accession-" + nl + "cluster and two-way" + nl + "(protein x study)" + nl + "intervals beside it")
    box(79, 22, 21, 20, "Signals" + nl + "decision-bearing:" + nl + "conservation, disorder," + nl + "sequence motif" + nl + "population-dependent:" + nl + "domain overlap" + nl + "descriptive only:" + nl + "lookup motif")
    box(79, 6, 21, 14, "Predictor audit" + nl + "(secondary): DeepPhos," + nl + "MusiteDeep, PhosF3C", ls=(0, (3, 2)))
    arrow(75, 52.5, 79, 52.5)
    arrow(75, 33, 79, 33)
    arrow(75, 13, 79, 13)
    if standalone():
        title(ax, "Study design")
        ax.text(0, -3, "Dashed outline: check or secondary analysis", fontsize=FS - 0.5, va="center", ha="left")
    legend_line(fig, "Boxes with dashed outlines are checks or secondary analyses.", y=0.005)
    save(fig, "Figure1", rec)


HATCH = {"GOLD": "", "SILVER": "", "BRONZE": "////"}


def figure2():
    rec = Rec("Figure2")
    c, d, g23 = load("figure_counts.json"), load("derived_figure_values.json"), load("g23_tiers_summary.json")
    fig = newfig((W, 5.6))
    gs = fig.add_gridspec(1, 6, left=0.075, right=0.99, top=0.935, bottom=0.62, wspace=1.9)
    gs2 = fig.add_gridspec(1, 2, left=0.33, right=0.99, top=0.46, bottom=0.215, wspace=0.55, width_ratios=[1.0, 0.9])

    def stack(ax, x, shares, counts_n=None, wd=0.62, minlab=0.04):
        bottom = 0.0
        for t in ["BRONZE", "SILVER", "GOLD"]:
            v = shares.get(t)
            if v is None:
                continue
            ax.bar(x, float(v) * 100, bottom=bottom, color=TIER[t], width=wd, edgecolor="white", lw=0.5, hatch=HATCH[t])
            if float(v) >= minlab:
                ax.text(x, bottom + float(v) * 50, f"{v * 100:.1f}%", ha="center", va="center", fontsize=FS - 0.5, color="white" if t == "SILVER" else INK)
            bottom += float(v) * 100

    # A tiers by species
    ax = pax(fig, gs[0, 0:2], "A"); letter(ax, "A")
    order = [("all", "All"), ("human", "Human"), ("mouse", "Mouse"), ("rat", "Rat")]
    for i, (k, lab) in enumerate(order):
        stack(ax, i, {t: d[f"fig2a|{k}|{t}"]["value"] for t in TIERS})
    ax.set_xticks(range(4)); ax.set_xticklabels([lab for k, lab in order])
    ax.set_ylabel("Sites (%)"); ax.set_ylim(0, 100)
    title(ax, "Tier composition (study-cluster rule)")
    ax.legend([plt.Rectangle((0, 0), 1, 1, fc=TIER[t], ec="white", hatch=HATCH[t]) for t in TIERS], TIERS, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.27), handlelength=0.9, columnspacing=0.8, handletextpad=0.3)
    # B two rules
    ax = pax(fig, gs[0, 2:4], "B"); letter(ax, "B")
    bars = [("fig2b|pmid_rule_set_a", "PubMed-ID" + chr(10) + "count" + chr(10) + "(set A)"), ("fig2b|cluster_rule_set_a", "Study-" + chr(10) + "cluster" + chr(10) + "rule" + chr(10) + "(set A)"),
            ("fig2b|cluster_rule_set_b", "Study-" + chr(10) + "cluster" + chr(10) + "rule" + chr(10) + "(set B," + chr(10) + "primary)")]
    for i, (k, lab) in enumerate(bars):
        stack(ax, i, {t: d[f"{k}|{t}"]["value"] for t in TIERS})
    ax.set_xticks(range(3)); ax.set_xticklabels([b[1] for b in bars], fontsize=FS - 0.7)
    ax.set_ylabel(f"Human sites (%), n = {fmt_int(c['human_n_sites'])}"); ax.set_ylim(0, 100)
    sv = c["human_tier_pmid_count_rule"]["SILVER"]
    ax.annotate("SILVER:" + chr(10) + f"{fmt_int(sv)} sites" + chr(10) + f"({d['fig2b|pmid_rule_set_a|SILVER']['value'] * 100:.1f}%)", xy=(0, float(d['fig2b|pmid_rule_set_a|BRONZE']['value']) * 100 + 0.3), xytext=(0, 65), fontsize=FS - 1.0, ha="center", va="center", bbox=dict(fc="white", ec="none", pad=0.6),
                arrowprops=dict(arrowstyle="-", lw=0.5, color=INK))
    ax.set_xlim(-0.6, 2.6)
    title(ax, "Two rules (thresholds differ)")
    if standalone():
        MODE["notight"] = True
        ax.set_position([0.2, 0.42, 0.76, 0.50])
        ax.legend([plt.Rectangle((0, 0), 1, 1, fc=TIER[t], ec="white", hatch=HATCH[t]) for t in TIERS], TIERS, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.46), handlelength=0.9, columnspacing=0.8, handletextpad=0.3)
        src_box(ax, (0.5, -0.8), SRC_A + " " + SRC_B, 48, ha="center", va="bottom")
    # C DAS >= 2
    ax = pax(fig, gs[0, 4:6], "C"); letter(ax, "C")
    a_, b_ = c["human_das_ge2_set_a"], c["human_das_ge2_set_b"]
    ax.bar([0, 1], [float(a_["share"]) * 100, float(b_["share"]) * 100], color=GREY, width=0.55)
    ax.text(0, float(a_["share"]) * 100 + 0.5, f"{a_['share'] * 100:.2f}%" + chr(10) + f"n = {fmt_int(a_['n'])}", ha="center", va="bottom", fontsize=FS - 0.5)
    ax.text(1, float(b_["share"]) * 100 + 0.5, f"{b_['share'] * 100:.1f}%" + chr(10) + f"n = {fmt_int(b_['n'])}", ha="center", va="bottom", fontsize=FS - 0.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Source set A", "Source set B"])
    ax.set_ylabel("Human sites, score 2 or more (%)"); ax.set_ylim(0, 34); ax.set_xlim(-0.6, 1.6)
    title(ax, "Database agreement")
    src_box(ax, (0.5, 0.98), SRC_A_SHORT + " " + SRC_B_SHORT, 36 if not standalone() else 46, ha="center")
    # D evidence accounting
    ax = pax(fig, gs2[0, 0], "D")
    if not standalone():
        fig.text(0.012, 0.465, "D", fontsize=9, fontweight="bold", va="bottom")
    items = [(g23["records_with_at_least_one"], f"Records with at least one non-primary PMID" + chr(10) + f"({g23['n_nonprimary_pmids']} database or review PMIDs)", g23["share_records_with_at_least_one"]),
             (next(x for x in g23["nonprimary"] if x["pmid"] == "25514926")["records_linked"], "Records with the PhosphoSitePlus database paper" + chr(10) + "(PMID 25514926)", d["fig2d|psp_paper_share_of_records"]["value"]),
             (g23["transition_table"]["GOLD"]["SILVER"], "GOLD records that fall to SILVER when" + chr(10) + "these PMIDs are not counted", None)]
    ys = [2, 1, 0]
    ax.barh(ys, [float(v[0]) / 1000 for v in items], color=[GREY, GREY, TIER["GOLD"]], height=0.55, edgecolor=INK, lw=0.4)
    for y, (v, lab, sh) in zip(ys, items):
        ax.text(float(v) / 1000 + 12, y, fmt_int(v) + (f" ({sh * 100:.1f}% of records)" if sh is not None else f" ({d['fig2d|gold_to_silver_share_of_gold']['value'] * 100:.1f}% of GOLD)"), va="center", fontsize=FS - 0.3)
    ax.set_yticks(ys); ax.set_yticklabels([i[1] for i in items], fontsize=FS - 0.5)
    ax.set_xlabel(f"Records (thousands); {fmt_int(g23['n_records'])} records in all"); ax.set_xlim(0, 1500)
    title(ax, "Database and review papers counted as evidence")
    # E first-reported year
    ax = pax(fig, gs2[0, 1], "E"); letter(ax, "E")
    for i, t in enumerate(TIERS):
        for dx, key, al in ((-0.2, "first_year_by_tier_primary", 1.0), (0.2, "first_year_by_tier_g23", 0.45)):
            st = c[key][t]
            sts = {k: float(v) for k, v in st.items()}
            box_draw(ax, i + dx, sts, TIER[t], width=0.34, alpha=al)
            ax.text(i + dx, sts["p95"] + 0.7, f"{st['median']:.0f}", ha="center", va="bottom", fontsize=FS - 0.6)
    ax.set_xticks(range(3)); ax.set_xticklabels(TIERS); ax.set_xlim(-0.6, 2.6); ax.set_ylim(1998, 2036)
    ax.set_ylabel("First-reported year")
    ax.legend([plt.Rectangle((0, 0), 1, 1, fc=GREY, ec=INK), plt.Rectangle((0, 0), 1, 1, fc=GREY, ec=INK, alpha=0.45)], ["Study-cluster depth", "Non-primary PMIDs not counted"], loc="upper left", fontsize=FS - 0.5, handlelength=0.9)
    title(ax, "First report by tier")
    n1, n2 = c["first_year_by_tier_primary"], c["first_year_by_tier_g23"]
    legend_line(fig, {
        "A": "n = " + ", ".join(f"{k} {fmt_int(c['n_sites_by_species'][k])}" for k in ('human', 'mouse', 'rat')) + f", all {fmt_int(g23['n_records'])} sites; set-B tiers, study-cluster rule; BRONZE is hatched.",
        "B": "the PubMed-ID count rule and the study-cluster rule differ in thresholds, so the change cannot be attributed to the counting method alone; set A tiers come from the human file with two source flags (agreement score 0 to 2); set B adds the dbPTM flag (0 to 3); BRONZE is hatched.",
        "C": "same sources; n is the number of sites with a score of 2 or more.",
        "D": "records with at least one non-primary PMID, records linked to the PhosphoSitePlus database paper, and GOLD records that fall to SILVER when the non-primary PMIDs are not counted.",
        "E": "set-B tiers, study-cluster rule; left box primary tiers, right box tiers with the non-primary PMIDs not counted (sites with a year: " + ", ".join(f"{t} {fmt_int(n1[t]['n'])} and {fmt_int(n2[t]['n'])}" for t in TIERS) +
             "); boxes span the interquartile range, whiskers the 5th to 95th percentile, printed values are medians; the axis begins in 1998."}, y=0.004)
    save(fig, "Figure2", rec)


if __name__ == "__main__":
    setup()
    figure1()
    figure2()
    print("Figures 1 and 2 written")
