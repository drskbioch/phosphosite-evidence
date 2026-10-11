"""Supplementary figures S1 to S8 of the Rev 6 rebuild."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from lib import (FS, FS_TITLE, GREY, INK, LIGHT, POP, SIGNAL, SRC_A, SRC_B, SRC_B_FULL, TIER, TIERS, W, Rec, src_box, box_draw, ci_marker, fmt_int, legend_line, letter, load, load_dir, newfig, pax, save, setup, standalone, title)

NL = chr(10)


def pick(d, prefix="", **cond):
    hits = [(k, v) for k, v in d.items() if isinstance(v, dict) and k.startswith(prefix) and all(v.get(c) == x for c, x in cond.items())]
    assert len(hits) == 1, (cond, [k for k, _ in hits])
    return hits[0]


def fci(r, key="ci"):
    return float(r["estimate"]), [float(x) for x in r[key]]


# ---------------------------------------------------------------------------- S1
def figS1():
    rec = Rec("FigureS1")
    c, d = load("figure_counts.json"), load("derived_figure_values.json")
    fig = newfig((W, 6.2))
    gs = fig.add_gridspec(2, 6, left=0.075, right=0.975, top=0.94, bottom=0.2, wspace=2.2, hspace=0.75)
    ax = pax(fig, gs[0, 0:2], "A"); letter(ax, "A")
    xs = ["0", "1", "2", "3"]
    ax.bar(range(4), [float(c["das_counts_all"][k]) / 1000 for k in xs], color=GREY, width=0.65)
    for i, k in enumerate(xs):
        ax.text(i, float(c["das_counts_all"][k]) / 1000 + 15, f"{d['figS1|das|' + k]['value'] * 100:.1f}%", ha="center", fontsize=FS - 0.5)
    ax.set_xticks(range(4)); ax.set_xticklabels(xs); ax.set_xlabel("Database agreement score (source set B)"); ax.set_ylabel("Sites (thousands)"); ax.set_ylim(0, 1000)
    title(ax, "Database agreement")
    src_box(ax, (0.97, 0.97), SRC_B_FULL + " Source set A (human file): the first two flags; score 0 to 2.", 27, ha="right")
    ax = pax(fig, gs[0, 2:4], "B"); letter(ax, "B")
    for key, col, lab, ls in (("depth_histogram_primary", INK, "study-cluster depth", "-"), ("depth_histogram_g23", TIER["GOLD"], "non-primary PMIDs not counted", "--")):
        h = c[key]
        xs_ = sorted(int(k) for k in h)
        ax.step(xs_, [float(h[str(k)]) / 1000 for k in xs_], where="mid", color=col, lw=0.9, ls=ls, label=lab)
    ax.axvline(float(c["depth_median_primary"]), color=GREY, lw=0.6, ls=":")
    ax.text(float(c["depth_median_primary"]) * 1.15, 400, f"median {c['depth_median_primary']:.0f}", fontsize=FS - 0.5)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Study clusters per site"); ax.set_ylabel("Sites (thousands)"); ax.legend(loc="lower left", fontsize=FS - 0.6)
    title(ax, "Publication depth")
    ax = pax(fig, gs[0, 4:6], "C"); letter(ax, "C")
    parts = [("resolved", "Year found", GREY), ("valid_format_not_found", "Valid format, not found", LIGHT), ("invalid_format_over_45m", "Invalid format", "white")]
    left = 0
    for k, lab, col in parts:
        w_ = float(d[f"figS1|pmid|{k}"]["value"]) * 100
        ax.barh(0, w_, left=left, color=col, edgecolor=INK, lw=0.5, height=0.5)
        left += w_
    ax.text(float(d["figS1|pmid|resolved"]["value"]) * 50, 0, f"Year found" + NL + f"{fmt_int(c['pmid']['resolved'])} ({d['figS1|pmid|resolved']['value'] * 100:.1f}%)", ha="center", va="center", color="white", fontsize=FS - 0.5)
    ax.annotate(f"Valid format, not found" + NL + f"{fmt_int(c['pmid']['valid_format_not_found'])} ({d['figS1|pmid|valid_format_not_found']['value'] * 100:.1f}%)", xy=(87.5, -0.25), xytext=(62, -0.95), fontsize=FS - 0.6, ha="center", arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.annotate(f"Invalid format" + NL + f"{fmt_int(c['pmid']['invalid_format_over_45m'])} ({d['figS1|pmid|invalid_format_over_45m']['value'] * 100:.1f}%)", xy=(93, 0.25), xytext=(80, 0.95), fontsize=FS - 0.6, ha="center", arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.set_xlim(0, 100); ax.set_ylim(-1.5, 1.6); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlabel(f"Percent of {fmt_int(c['pmid']['distinct'])} distinct PubMed IDs")
    title(ax, "Publication-year lookup")
    ax = pax(fig, gs[1, 0:3], "D"); letter(ax, "D")
    dl = c["pmid"]["unresolved_digit_lengths"]
    ks = sorted(int(k) for k in dl)
    ax.bar(range(len(ks)), [float(dl[str(k)]) for k in ks], color=GREY, width=0.6)
    for i, k in enumerate(ks):
        ax.text(i, float(dl[str(k)]) * 1.2, fmt_int(dl[str(k)]), ha="center", fontsize=FS - 0.5)
    ax.set_yscale("log"); ax.set_ylim(1, 3e4); ax.set_xticks(range(len(ks))); ax.set_xticklabels([str(k) for k in ks])
    ax.set_xlabel("Digits in the PubMed ID"); ax.set_ylabel(f"Identifiers with no publication year")
    title(ax, f"Identifiers with no year ({fmt_int(c['pmid']['unresolved'])})")
    ax = pax(fig, gs[1, 3:6], "E"); letter(ax, "E")
    sp = ["human", "mouse", "rat"]
    cols = {"human": "#332288", "mouse": "#44AA99", "rat": "#DDCC77"}
    for i, s in enumerate(sp):
        r = d[f"figS1|noncanonical_rate|{s}"]["value"]
        ax.bar(i, float(r) * 1e5, color=cols[s], edgecolor=INK, lw=0.4, width=0.6)
        ax.text(i, float(r) * 1e5 + 0.15, f"{fmt_int(c['noncanonical_by_species'][s])} of {fmt_int(c['n_sites_by_species'][s])}", ha="center", fontsize=FS - 0.7)
    ax.set_xticks(range(3)); ax.set_xticklabels([s.capitalize() for s in sp]); ax.set_ylabel("Non-canonical residue sites" + NL + "per 100,000 sites"); ax.set_ylim(0, 6)
    title(ax, "Non-canonical residue type")
    legend_line(fig, {"A": f"all {fmt_int(load('g23_tiers_summary.json')['n_records'])} sites; source set B counts the PhosphoSitePlus, Phospho.ELM and dbPTM flags (score 0 to 3).",
                      "B": "sites with the study-cluster depth, and with the non-primary PMIDs not counted; dotted line: median depth.",
                      "C": "identifiers above 45,000,000 are treated as invalid; identifiers of valid format without a year were not found.",
                      "D": "identifiers above 45,000,000 are treated as invalid; identifiers of valid format without a year were not found.",
                      "E": f"{fmt_int(c['noncanonical_n'])} sites in all."}, y=0.004)
    save(fig, "FigureS1", rec)


# ---------------------------------------------------------------------------- S2
def figS2():
    rec = Rec("FigureS2")
    c, b5, d3, g20, gu = load("figure_counts.json"), load("b5_index_effect_sizes.json"), load("d3_d4_d5_index_checks.json"), load("g20_diagnostics.json"), load("gene_unit_summary.json")
    H = 6.9  # inches (496.8 pt); positions in inches from the top: grid 0.41 to 3.64 (rows 1.24, gap 0.75), lower row 4.36 to 5.62, then the note
    fig = newfig((W, H))
    gs = fig.add_gridspec(2, 6, left=0.075, right=0.99, top=1 - 0.41 / H, bottom=1 - 3.64 / H, wspace=2.3, hspace=0.605)
    gs3 = fig.add_gridspec(1, 2, left=0.27, right=0.99, top=1 - 4.36 / H, bottom=1 - 5.62 / H, wspace=0.55, width_ratios=[1.6, 1])
    ax = pax(fig, gs[0, 0:2], "A"); letter(ax, "A")
    cols = {"human": "#332288", "mouse": "#44AA99", "rat": "#DDCC77"}
    for i, s in enumerate(["human", "mouse", "rat"]):
        st = {k: float(v) for k, v in c["v50_by_species"][s].items()}
        box_draw(ax, i, st, cols[s], width=0.5)
        ax.plot(i, st["mean"], "o", mfc="white", mec=INK, ms=3, mew=0.6, zorder=4)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Human", "Mouse", "Rat"]); ax.set_ylabel("V50 index"); title(ax, "V50 by species")
    if standalone():
        ax.legend([plt.Line2D([], [], marker="o", mfc="white", mec=INK, ls="none", ms=3, mew=0.6)], ["mean"], loc="upper center", bbox_to_anchor=(0.5, -0.1), fontsize=FS - 0.8, handletextpad=0.2)
    ax = pax(fig, gs[0, 2:4], "B"); letter(ax, "B")
    labs = []
    for i, t in enumerate(TIERS):
        st = c["v50_by_highest_tier"][t]
        box_draw(ax, i, {k: float(v) for k, v in st.items()}, TIER[t], width=0.5)
        ax.plot(i, float(st["mean"]), "o", mfc="white", mec=INK, ms=3, mew=0.6, zorder=4)
        labs.append(f"{t}")
    ax.set_xticks(range(3)); ax.set_xticklabels(labs); ax.set_ylabel("V50 index")
    e = b5["effects"]["kruskal_eta2_SEI_V50_by_top_tier"]
    ax.text(0.03, 0.97, f"eta-squared {e['estimate']:.4f}" + NL + f"({e['ci95'][0]:.4f} to {e['ci95'][1]:.4f})", transform=ax.transAxes, va="top", fontsize=FS - 0.5)
    title(ax, "V50 by highest tier")
    if standalone():
        ax.legend([plt.Line2D([], [], marker="o", mfc="white", mec=INK, ls="none", ms=3, mew=0.6)], ["mean"], loc="upper center", bbox_to_anchor=(0.5, -0.1), fontsize=FS - 0.8, handletextpad=0.2)
    ax = pax(fig, gs[0, 4:6], "C"); letter(ax, "C")
    for i, t in enumerate(TIERS):
        st = c["secondary_by_highest_tier"][t]
        box_draw(ax, i, {k: float(v) for k, v in st.items()}, TIER[t], width=0.5)
        ax.plot(i, float(st["mean"]), "o", mfc="white", mec=INK, ms=3, mew=0.6, zorder=4)
    rb = b5["effects"]["rank_biserial_SEI_secondary_GOLD_vs_BRONZE"]
    ax.set_xticks(range(3)); ax.set_xticklabels(TIERS); ax.set_ylabel("Secondary index")
    ax.text(0.97, 0.97, f"rank-biserial, GOLD against BRONZE:" + NL + f"{rb['estimate']:.4f} ({rb['ci95'][0]:.4f} to {rb['ci95'][1]:.4f})", transform=ax.transAxes, va="top", ha="right", fontsize=FS - 0.8)
    title(ax, "Secondary index by tier")
    if standalone():
        ax.legend([plt.Line2D([], [], marker="o", mfc="white", mec=INK, ls="none", ms=3, mew=0.6)], ["mean"], loc="upper center", bbox_to_anchor=(0.5, -0.1), fontsize=FS - 0.8, handletextpad=0.2)
    ax = pax(fig, gs[1, 0:2], "D"); letter(ax, "D")
    hb = c["v50_vs_secondary_bins"]
    cnt = np.array(hb["counts"], float).T
    pm = ax.pcolormesh(np.array([float(x) for x in hb["x_edges"]]), np.array([float(x) for x in hb["y_edges"]]), np.ma.masked_equal(cnt, 0), cmap="Greys", norm=LogNorm())
    cb = ax.figure.colorbar(pm, ax=ax, fraction=0.05, pad=0.03); cb.set_label("Proteins per bin", fontsize=FS - 0.5); cb.ax.tick_params(labelsize=FS - 1)
    rho = b5["effects"]["spearman_SEI_V50_vs_secondary"]
    ax.text(0.97, 0.04, f"Spearman {rho['estimate']:.4f}" + NL + f"({rho['ci95'][0]:.4f} to {rho['ci95'][1]:.4f})" + NL + f"n = {fmt_int(b5['n_proteins'])} proteins", transform=ax.transAxes, va="bottom", ha="right", fontsize=FS - 0.8, bbox=dict(fc="white", ec="none", pad=1.0))
    ax.set_xlabel("V50 index"); ax.set_ylabel("Secondary index"); title(ax, "V50 against secondary index")
    ax = pax(fig, gs[1, 2:4], "E"); letter(ax, "E")
    pb = b5["effects"]["pointbiserial_annotation_SEI_V50"]
    ci_marker(ax, 0, float(pb["estimate"]), [float(x) for x in pb["ci95"]], INK, ms=4)
    ax.axvline(0, color=GREY, lw=0.5)
    ax.set_yticks([0]); ax.set_yticklabels(["V50 index"]); ax.set_ylim(-1, 1); ax.set_xlim(-0.3, 0.1)
    ax.text(float(pb["ci95"][1]) + 0.01, 0, f"{pb['estimate']:.3f}", va="center", fontsize=FS - 0.6)
    ax.set_xlabel("Point-biserial correlation with" + NL + "domain-annotation presence")
    title(ax, "Annotation presence")
    ax = pax(fig, gs[1, 4:6], "F"); letter(ax, "F")
    tdl = d3["tdl_means"]["own_set|V50"]
    for i, k in enumerate(["Tclin", "Tchem", "Tbio", "Tdark"]):
        ax.bar(i, float(tdl[k]["mean"]), color=GREY, width=0.6)
        ax.text(i, float(tdl[k]["mean"]) + 0.05, f"n = {fmt_int(tdl[k]['n'])}", ha="center", fontsize=FS - 0.9)
    ax.set_xticks(range(4)); ax.set_xticklabels(["Tclin", "Tchem", "Tbio", "Tdark"]); ax.set_ylabel("Mean V50 index"); ax.set_ylim(0, 3.6)
    au = d3["auroc_tdark"]["own_set|V50"]
    ax.text(0.02, 0.99, f"AUROC non-Tdark against Tdark:" + NL + f"{au['auroc_nonTdark_vs_Tdark']:.4f} ({au['ci95'][0]:.4f} to {au['ci95'][1]:.4f})", transform=ax.transAxes, va="top", ha="left", fontsize=FS - 0.8)
    title(ax, "Pharos target development level")
    ax = pax(fig, gs3[0, 0], "G"); letter(ax, "G", dx=-132)
    ioa = g20["index_outcome_association"]
    ys, labs = [], []
    y = 0
    recs = [(k, v) for k, v in sorted(ioa.items()) if v.get("tier") in ("GOLD", "BRONZE") and v["signal"] in ("S2", "S3", "S4")]
    for k, v in recs[::-1]:
        ci_marker(ax, y, float(v["pearson_r"]), [float(x) for x in v["ci"]], TIER[v["tier"]], ms=3.2, cap=1.2)
        ys.append(y); labs.append(f"{SIGNAL[v['signal']]}, {v['population'].replace('P-', '')}, {v['tier']}")
        y += 1
    ax.axvline(0, color=GREY, lw=0.5)
    ax.set_yticks(ys); ax.set_yticklabels(labs, fontsize=FS - 1.1)
    ax.set_xlabel("Pearson correlation of the V50 index with the outcome," + NL + "within a tier (95% gene-cluster interval)" if standalone() else "Pearson correlation of the V50 index with the outcome, within a tier (95% gene-cluster interval)")
    title(ax, "Index and outcome within a tier")
    if standalone():
        ax.legend([plt.Line2D([], [], marker="o", color=TIER["GOLD"], ls="-", ms=3.2, lw=0.8), plt.Line2D([], [], marker="o", color=TIER["BRONZE"], ls="-", ms=3.2, lw=0.8)], ["GOLD sites", "BRONZE sites"], loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=2, fontsize=FS - 0.8, handletextpad=0.3)
    ax = pax(fig, gs3[0, 1], "H"); letter(ax, "H", dx=-8)
    sv = g20["index_vs_sites_per_accession"]
    items = [("human|sites_per_accession", "All sites per protein"), ("human|gold_sites_per_accession", "GOLD sites per protein")]
    for i, (k, lab) in enumerate(items[::-1]):
        ci_marker(ax, i, float(sv[k]["spearman"]), [float(x) for x in sv[k]["ci"]], INK, ms=4)
        ax.text(float(sv[k]["ci"][1]) + 0.004, i, f"{sv[k]['spearman']:.3f}", va="center", fontsize=FS - 0.6)
    ax.axvline(0, color=GREY, lw=0.5)
    ax.set_yticks([0, 1]); ax.set_yticklabels([i[1] for i in items[::-1]], fontsize=FS - 0.8); ax.set_ylim(-0.6, 1.6); ax.set_xlim(-0.02, 0.14)
    ax.set_xlabel("Spearman correlation with V50")
    title(ax, "Index and sites per protein")
    nb = f"{fmt_int(sv['human|sites_per_accession']['n_accessions'])} human proteins, {fmt_int(sv['human|sites_per_accession']['n_gene_clusters'])} gene clusters"
    tiers_n = ", ".join(f"{t} {fmt_int(c['v50_by_highest_tier'][t]['n'])}" for t in TIERS)
    legend_line(fig, {"A": "boxes span the interquartile range, whiskers the 5th to 95th percentile, open circles are means; " + ", ".join(f"{s} {fmt_int(c['v50_by_species'][s]['n'])}" for s in ('human', 'mouse', 'rat')) + " proteins.",
                      "B": "boxes span the interquartile range, whiskers the 5th to 95th percentile, open circles are means; proteins by highest tier: " + tiers_n + "; no interval is drawn on the means; the eta-squared interval is the stored protein-bootstrap interval.",
                      "C": "boxes span the interquartile range, whiskers the 5th to 95th percentile, open circles are means; proteins by highest tier: " + tiers_n + "; no interval is drawn on the means; the secondary index is built from the tiers' own publications, which is why it was not used for matching.",
                      "D": "counts per bin.",
                      "E": fmt_int(gu['human']['accessions_with_gene']) + " human proteins; annotation presence is a proxy for the outcome, not an independent variable.",
                      "F": nb + "; Pharos target classes, V50 index.",
                      "G": "GOLD and BRONZE sites in P-all and P-human; 95% gene-cluster intervals.",
                      "H": nb + "."}, y=0.004)
    save(fig, "FigureS2", rec)


# ---------------------------------------------------------------------------- S3
def figS3():
    rec = Rec("FigureS3")
    am, g20 = load("amend_results.json"), load("g20_diagnostics.json")
    pb = g20["prematching_balance"]
    rows = []
    for k, v in sorted(pb.items()):
        r = am.get(k)
        if k.startswith("G17") and r is not None and r["design"] == "D1":
            rows.append((f"{SIGNAL[r['signal']]}, {r['population']}", v))
    fig = newfig((W, 4.2))
    ax = pax(fig, [0.28, 0.2, 0.69, 0.72], "A"); letter(ax, "A", dx=-170)
    ax.axvspan(-0.1, 0.1, color="#ececec", zorder=0)
    ys = list(range(len(rows)))[::-1]
    for y, (lab, v) in zip(ys, rows):
        ax.plot([float(v["smd_pre"]), float(v["smd_after"])], [y, y], color=LIGHT, lw=1.2, zorder=1)
        ax.plot(float(v["smd_pre"]), y, "o", mfc="white", mec=INK, ms=4, zorder=3)
        ax.plot(float(v["smd_after"]), y, "o", color=INK, ms=4, zorder=3)
    ax.axvline(0, color=GREY, lw=0.5)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows], fontsize=FS - 0.5)
    ax.set_xlabel("Standardised mean difference of the V50 index, GOLD against BRONZE"); title(ax, "Covariate balance before (open) and after (filled) matching on the V50 decile")
    if standalone():
        ax.legend([plt.Line2D([], [], marker="o", mfc="white", mec=INK, ls="none", ms=4), plt.Line2D([], [], marker="o", color=INK, ls="none", ms=4), plt.Rectangle((0, 0), 1, 1, fc="#ececec", ec="none")],
                  ["before matching", "after matching", "balance criterion of the plan"], loc="upper center", bbox_to_anchor=(0.4, -0.2), ncol=3, fontsize=FS - 0.8, handletextpad=0.3, columnspacing=1.0)
    legend_line(fig, {"": "Shaded band: the balance criterion of the plan. CEM on the V50 decile, set-B tiers; the secondary-index design is not used."}, y=0.01)
    save(fig, "FigureS3", rec)


# ---------------------------------------------------------------------------- S4
def figS4():
    rec = Rec("FigureS4")
    g14, d, b4, g21 = load("g14_cem.json"), load("derived_figure_values.json"), load("b4_twoway.json"), load("g21_results.json")
    H = 7.2  # inches (518 pt); the rows are given in inches: top margin 0.35, rows 2.2 and 2.45, gap 0.95, bottom 1.2 (rotated tick labels and the note)
    fig = newfig((W, H))
    gs = fig.add_gridspec(2, 2, left=0.2, right=0.965, top=1 - 0.35 / H, bottom=1.2 / H, wspace=0.75, hspace=0.41, height_ratios=[2.2, 2.45])
    ax = pax(fig, gs[0, 0], "A"); letter(ax, "A", dx=-60)
    sig = {"S1": "Domain overlap", "S2": "Conservation", "S3": "Disorder", "S4": "Sequence motif", "S5": "Lookup motif"}
    ys = list(range(5))[::-1]
    for y, s in zip(ys, ["S1", "S2", "S3", "S4", "S5"]):
        for key, mk, filled in (("R5_CEM_primary_index", "o", True), ("R5_CEM_secondary_index", "D", False)):
            r = next(v for v in g14.values() if v["signal"] == s and v["design"] == key)
            ax.plot(float(r["ci_width_ratio_protein_cluster_over_rev5"]), y + (0.12 if filled else -0.12), mk, mfc=INK if filled else "white", mec=INK, ms=4)
    ax.axvline(1, color=GREY, lw=0.5, ls=":")
    ax.set_yticks(ys); ax.set_yticklabels([sig[s] for s in ["S1", "S2", "S3", "S4", "S5"]], fontsize=FS - 0.6)
    ax.set_xlabel("Protein-cluster interval width / site-level interval width"); title(ax, "Cluster against site-level interval")
    ax.legend([plt.Line2D([], [], marker="o", color=INK, ls="none", ms=4), plt.Line2D([], [], marker="D", mfc="white", mec=INK, ls="none", ms=4)], ["CEM on the original index", "CEM on the secondary index"], numpoints=1, loc="lower right", fontsize=FS - 0.8)
    ax = pax(fig, gs[0, 1], "B"); letter(ax, "B", dx=-8)
    keys = [k for k in d if k.startswith("figS4|")]
    keys = sorted(keys)
    ys = list(range(len(keys)))[::-1]
    for y, k in zip(ys, keys):
        ax.plot(float(d[k]["value"]), y, "o", color=POP["P-human"] if "P-human" in k else INK, ms=4)
    ax.axvline(1, color=GREY, lw=0.5, ls=":")
    ax.set_yticks(ys); ax.set_yticklabels([k.split("|")[1] + ", " + sig.get(k.split("|")[2], k.split("|")[2]) for k in keys], fontsize=FS - 1.2)
    ax.set_xlabel("Gene-cluster interval width /" + NL + "accession-cluster interval width"); title(ax, "Gene against accession cluster")
    if standalone():
        ax.legend([plt.Line2D([], [], marker="o", color=INK, ls="none", ms=4), plt.Line2D([], [], marker="o", color=POP["P-human"], ls="none", ms=4)], ["P-all", "P-human"], loc="lower right", fontsize=FS - 0.8, handletextpad=0.2)
    ax = pax(fig, gs[1, 0], "C"); letter(ax, "C", dx=-60)
    ks = ["G1-0002", "G1-0005", "G1-0008", "G1-0011", "G1-0014"]
    names = ["Domain overlap", "Conservation", "Disorder", "Sequence motif", "Lookup motif (descriptive)"]
    for i, (k, nm) in enumerate(zip(ks, names)):
        ax.bar(i, float(b4[k]["width_ratio_two_way_vs_protein"]), color=GREY, width=0.6)
        ax.text(i, float(b4[k]["width_ratio_two_way_vs_protein"]) * 1.1, f"{b4[k]['width_ratio_two_way_vs_protein']:.1f}", ha="center", fontsize=FS - 0.6)
    ax.set_yscale("log"); ax.set_ylim(1, 90); ax.axhline(1, color=GREY, lw=0.5)
    ax.set_xticks(range(5)); ax.set_xticklabels([n.replace(" (descriptive)", "") for n in names], rotation=30, ha="right", fontsize=FS - 0.8)
    ax.set_ylabel("Two-way interval width /" + NL + "protein interval width"); title(ax, "Two-way interval (P-all, primary contrast)")
    ax = pax(fig, gs[1, 1], "D"); letter(ax, "D", dx=-8)
    rows = [(k, v) for k, v in sorted(g21.items()) if isinstance(v, dict)]
    ys = list(range(len(rows)))[::-1]
    for y, (k, v) in zip(ys, rows):
        valid = float(v["n_boot_valid"])
        ax.barh(y, valid, color=GREY if float(v["n_boot_unfit"]) == 0 else TIER["GOLD"], height=0.7)
        if float(v["n_boot_unfit"]) > 0:
            ax.text(8, y, f"{fmt_int(v['n_boot_unfit'])} failed", va="center", ha="left", fontsize=FS - 1.4, color=INK)
    ax.set_yticks(ys); ax.set_yticklabels([f"{sig.get(v['signal'], v['signal'])}, {v['population'].split(' ')[0].replace('P-', '')}, {v['model']}" for _, v in rows], fontsize=FS - 1.6)
    ax.set_xlim(0, 1000); ax.set_xticks([0, 250, 500, 750, 1000]); ax.set_xlabel("Valid bootstrap fits"); title(ax, "Within-gene models")
    if standalone():
        ax.legend([plt.Rectangle((0, 0), 1, 1, fc=GREY), plt.Rectangle((0, 0), 1, 1, fc=TIER["GOLD"])], ["no failed fit", "some fits failed"], loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=FS - 0.8, handletextpad=0.3)
    nb = g21[rows[0][0]]["n_boot"]
    legend_line(fig, {"A": "ratio of the protein-cluster interval to the site-level interval for the same estimate (estimates from the site-level analysis; compares interval methods only).",
                      "B": "gene-cluster against accession-cluster interval, D1 contrasts, populations P-all (black) and P-human (purple).",
                      "C": "ratio of the width of the two-way (protein x study) interval to the width of the protein-cluster interval, primary contrast in P-all.",
                      "D": f"out of {fmt_int(nb)} resamples; orange bars mark rows with failed fits."}, y=0.008)
    save(fig, "FigureS4", rec)


# ---------------------------------------------------------------------------- S5
def figS5():
    rec = Rec("FigureS5")
    c, d, g23, b1 = load("figure_counts.json"), load("derived_figure_values.json"), load("g23_tiers_summary.json"), load("b1_reconciliation.json")
    fig = newfig((W, 5.3))
    cs, ts = c["cluster_split"], c["temporal_split"]
    ax = pax(fig, [0.115, 0.63, 0.30, 0.27], "A"); letter(ax, "A", dx=-55)
    for i, (lab, dd) in enumerate((("Cluster split", cs), ("Temporal split", ts))):
        ax.barh(1 - i, float(dd["train"]) / 1e6, color=GREY, height=0.5)
        ax.barh(1 - i, float(dd["test"]) / 1e6, left=float(dd["train"]) / 1e6, color=LIGHT, height=0.5)
        ax.text(float(dd["train"]) / 2e6, 1 - i, f"{fmt_int(dd['train'])}", ha="center", va="center", color="white", fontsize=FS - 0.8)
        ax.text(float(dd["train"]) / 1e6 + float(dd["test"]) / 2e6, 1 - i, f"{fmt_int(dd['test'])}", ha="center", va="center", fontsize=FS - 0.8)
    ax.set_yticks([1, 0]); ax.set_yticklabels(["Cluster split", "Temporal split"]); ax.set_xlabel("Sites (millions): training, test"); title(ax, "Split sizes")
    if standalone():
        ax.legend([plt.Rectangle((0, 0), 1, 1, fc=GREY), plt.Rectangle((0, 0), 1, 1, fc=LIGHT)], ["training", "test"], loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2, fontsize=FS - 0.8, handletextpad=0.3)
    ax = pax(fig, [0.665, 0.63, 0.31, 0.27], "B"); letter(ax, "B", dx=-118)
    items = [("Proteins, cluster split", "figS5|cluster_split|proteins_in_both", cs["proteins_in_both"], cs["proteins_total"]),
             ("Proteins, temporal split", "figS5|temporal_split|proteins_in_both", ts["proteins_in_both"], ts["proteins_total"]),
             ("Study clusters, temporal split", "figS5|temporal_split|clusters_in_both", ts["study_clusters_in_both"], ts["study_clusters_total"])]
    for i, (lab, k, a, b) in enumerate(items[::-1]):
        ax.barh(i, float(d[k]["value"]) * 100, color=GREY, height=0.5)
        ax.text(float(d[k]["value"]) * 100 + 0.8, i, f"{d[k]['value'] * 100:.1f}% ({fmt_int(a)} of {fmt_int(b)})", va="center", fontsize=FS - 0.6)
    ax.set_yticks(range(3)); ax.set_yticklabels([x[0] for x in items[::-1]], fontsize=FS - 0.6); ax.set_xlim(0, 100); ax.set_xlabel("Present in both partitions (%)"); title(ax, "Overlap between training and test")
    ax = pax(fig, [0.285, 0.2, 0.30, 0.27], "C"); letter(ax, "C", dx=-130)
    ho = g23["held_out"]
    names = [("cluster_disjoint_test", "Cluster-disjoint test set"), ("cluster_and_protein_disjoint_test", "Cluster- and protein-disjoint"), ("temporal_test_all", "Temporal test set"),
             ("temporal_protein_disjoint", "Temporal, protein-disjoint"), ("temporal_study_disjoint", "Temporal, study-disjoint"), ("temporal_protein_and_study_disjoint", "Temporal, protein- and study-disjoint")]
    for i, (k, lab) in enumerate(names[::-1]):
        ax.barh(i, float(ho[k]["n"]) / 1e3, color=GREY, height=0.5)
        ax.text(float(ho[k]["n"]) / 1e3 * 1.15, i, fmt_int(ho[k]["n"]), va="center", fontsize=FS - 0.6)
    ax.set_xscale("log"); ax.set_xlim(0.02, 8000); ax.set_yticks(range(6)); ax.set_yticklabels([n[1] for n in names[::-1]], fontsize=FS - 0.6)
    ax.set_xlabel("Sites in the test set (thousands)"); title(ax, "Held-out test sets")
    ax = pax(fig, [0.785, 0.2, 0.2, 0.27], "D"); letter(ax, "D", dx=-85)
    sc = b1["study_clusters"]
    vals = [("Labels, all sites", sc["distinct_primary_cluster_labels_all_records"]), ("Labels, temporal", sc["distinct_primary_cluster_labels_temporal_split_records"]),
            ("Clusters, all PMIDs", sc["distinct_cluster_ids_over_all_50904_pmids_under_rule_884"])]
    for i, (lab, v) in enumerate(vals[::-1]):
        ax.barh(i, float(v), color=GREY, height=0.5)
        ax.text(float(v) * 1.15, i, fmt_int(v), va="center", fontsize=FS - 0.6)
    ax.set_xscale("log"); ax.set_xlim(100, 3e6); ax.set_yticks(range(3)); ax.set_yticklabels([v[0] for v in vals[::-1]], fontsize=FS - 0.8)
    ax.set_xlabel("Study clusters"); title(ax, "Cluster denominators")
    split_txt = f"cluster split and temporal split of the {fmt_int(c['cluster_split']['proteins_total'])} proteins."
    legend_line(fig, {"A": split_txt, "B": split_txt, "C": "test sets as assigned before scoring.",
                      "D": "the number of study clusters depends on what is counted: distinct primary cluster labels (all sites, and sites in the temporal split) "
                           f"or clusters formed over all {fmt_int(c['pmid']['distinct'])} distinct PubMed IDs under the merge rule."}, y=0.008)
    save(fig, "FigureS5", rec)


# ---------------------------------------------------------------------------- S6
def figS6():
    rec = Rec("FigureS6")
    am, cem, g23, g18 = load("amend_results.json"), load("cem_results.json"), load("g23_results.json"), load_dir("d2b")
    g6, g7a, g7b, g8 = (load("f_runs/G6_S1-S2-S3-S4.json"), load("f_runs/G7_S1-S2.json"), load("f_runs/G7_S3-S4.json"), load("f_runs/G8_S1-S2-S3-S4.json"))
    g7 = {**g7a, **g7b}
    fig = newfig((W, 6.6))
    gs = fig.add_gridspec(1, 3, left=0.27, right=0.985, top=0.9, bottom=0.2, wspace=0.12)
    rows = ["Primary (CEM on the V50 decile)", "Reviewed canonical accessions", "Residue S", "Residue T", "Residue Y", "Leave one study cluster out", "Tier thresholds",
            "Depth-only tiers", "Agreement-only tiers", "1:1 caliper matcher", "Non-primary PMIDs not counted", "Cluster threshold (protein-cluster interval)", "Index variants (protein-cluster interval)"]
    ys = {r: len(rows) - 1 - i for i, r in enumerate(rows)}
    for ci_, s in enumerate(["S2", "S3", "S4"]):
        ax = pax(fig, gs[0, ci_], "ABC"[ci_]); letter(ax, "ABC"[ci_], dx=0 if ci_ else -130, dy=20)
        ax.set_title(SIGNAL[s], fontsize=FS_TITLE, fontweight="bold", loc="left", pad=4)
        r = pick(am, "G17", signal=s, population="P-all", design="D1")[1]
        e, ci = fci(r); ci_marker(ax, ys[rows[0]], e, ci, INK, ms=3.6)
        ax.axvline(e, color=INK, lw=0.4, ls=(0, (1.5, 1.5)), zorder=0); ax.axvline(0, color=GREY, lw=0.5, zorder=0)
        r = pick(am, "G16", signal=s, population="P-all, Swiss-Prot reviewed canonical", design="D1")[1]; ci_marker(ax, ys[rows[1]], *fci(r), INK, ms=3.6)
        for pop, lab in (("P-S", "Residue S"), ("P-T", "Residue T"), ("P-Y", "Residue Y")):
            r = next(v["gene"] for v in g6.values() if v["gene"]["signal"] == s and v["gene"]["population"] == pop and v["gene"]["design"] == "D1")
            ci_marker(ax, ys[lab], *fci(r), INK, ms=3.6)
        for lab, src, fam in (("Leave one study cluster out", g7, None), ("Tier thresholds", g8, None)):
            vs = [v["gene"] for v in src.values() if v["gene"]["signal"] == s and v["gene"]["design"] == "D1"]
            lo, hi = min(float(v["ci"][0]) for v in vs), max(float(v["ci"][1]) for v in vs)
            ax.plot([lo, hi], [ys[lab]] * 2, color=LIGHT, lw=3.5, solid_capstyle="butt", zorder=1)
            ax.plot([float(v["estimate"]) for v in vs], [ys[lab]] * len(vs), "o", color=INK, ms=2.6, zorder=3)
        for lab, tl in (("Depth-only tiers", "depth-only tiers"), ("Agreement-only tiers", "agreement-only")):
            r = next(v for k, v in am.items() if k.startswith("G19") and v["signal"] == s and v["population"] == "P-all" and v["design"] == "D1" and v["tier_rule"].startswith(tl[:5]))
            ci_marker(ax, ys[lab], *fci(r), INK, ms=3.6)
        r = pick(g18, signal=s, population="P-all (S/T/Y only)")[1]; ci_marker(ax, ys[rows[9]], *fci(r), INK, ms=3.6)
        r = pick(g23, signal=s, population="P-all", design="D1")[1]; ci_marker(ax, ys[rows[10]], *fci(r), INK, ms=3.6)
        for lab, fam, extra in ((rows[11], "G13", lambda v: v["design"] == "D1"), (rows[12], "G4", lambda v: True)):
            vs = [v for k, v in cem.items() if k.startswith(fam) and v["signal"] == s and v["population"] == "P-all" and extra(v)] if fam == "G13" else [v for k, v in cem.items() if k.startswith(fam) and v["signal"] == s]
            ax.plot([min(float(v["ci"][0]) for v in vs), max(float(v["ci"][1]) for v in vs)], [ys[lab]] * 2, color=LIGHT, lw=3.5, solid_capstyle="butt", zorder=1)
            ax.plot([float(v["estimate"]) for v in vs], [ys[lab]] * len(vs), "o", mfc="white", mec=INK, ms=2.8, zorder=3)
        first = ci_ == 0 or standalone()
        ax.set_yticks([ys[r_] for r_ in rows]); ax.set_yticklabels(rows if first else [], fontsize=FS - 0.6)
        if not first:
            ax.tick_params(axis="y", length=0)
        ax.set_ylim(-0.7, len(rows) - 0.3)
        ax.set_xlabel("Difference in mean disorder score" if s == "S3" else "Difference in proportion"); ax.margins(x=0.08)
        if standalone():
            ax.legend([plt.Line2D([], [], marker="o", color=INK, ls="-", ms=3.6, lw=0.8), plt.Line2D([], [], marker="o", color=INK, ls="none", ms=2.6), plt.Line2D([], [], color=LIGHT, lw=3.5),
                       plt.Line2D([], [], marker="o", mfc="white", mec=INK, ls="none", ms=2.8), plt.Line2D([], [], color=INK, lw=0.4, ls=(0, (1.5, 1.5)))],
                      ["estimate with 95% gene-cluster interval", "variant estimate", "range of the interval limits over variants", "estimate with stored protein-cluster interval", "primary estimate"],
                      loc="upper center", bbox_to_anchor=(0.35, -0.12), ncol=2, fontsize=FS - 0.8, handletextpad=0.3, columnspacing=1.0)
    legend_line(fig, {"": "GOLD against BRONZE in P-all, set-B tiers unless named. Intervals are 95% gene-cluster intervals, except the last two rows (stored protein-cluster intervals, open circles; the cluster thresholds 500 and 1,000 sites change no record, so their estimates equal the primary estimate). "
                "Rows with several variants (leave one study cluster out, tier thresholds, cluster threshold, index variants) show every variant estimate as a dot and the range of the interval limits as a grey bar. Dotted vertical line: primary estimate."}, y=0.008)
    save(fig, "FigureS6", rec)


# ---------------------------------------------------------------------------- S7
def figS7():
    rec = Rec("FigureS7")
    am, g21 = load("amend_results.json"), load("g21_results.json")
    fig = newfig((W, 3.6))
    ax = pax(fig, [0.3, 0.27, 0.68, 0.58], "A"); letter(ax, "A", dx=-170)
    rows = []
    for pop in ("P-all", "P-human"):
        r0 = pick(am, "G17", signal="S5", population=pop, design="D0")[1]
        r1 = pick(am, "G17", signal="S5", population=pop, design="D1")[1]
        m1 = pick(g21, signal="S5", population=lambda p: False, model="M1")[1] if False else next(v for v in g21.values() if isinstance(v, dict) and v["signal"] == "S5" and v["population"].startswith(pop + " (") and v["model"] == "M1")
        m2 = next(v for v in g21.values() if isinstance(v, dict) and v["signal"] == "S5" and v["population"].startswith(pop + " (") and v["model"] == "M2")
        rows += [(f"{pop}, crude", fci(r0)), (f"{pop}, CEM on the V50 decile", fci(r1)), (f"{pop}, within gene M1", fci(m1["companion_linear"])), (f"{pop}, within gene M2", fci(m2["companion_linear"]))]
    ys = list(range(len(rows)))[::-1]
    for y, (lab, (e, ci)) in zip(ys, rows):
        ci_marker(ax, y, e, ci, POP["P-human"] if "P-human" in lab else INK, ms=3.8)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows], fontsize=FS - 0.4)
    ax.set_xlabel("Difference in proportion of sites with a database-lookup motif match"); ax.margins(x=0.08)
    title(ax, "Lookup motif: descriptive only")
    r = pick(am, "G17", signal="S5", population="P-all", design="D1")[1]
    if standalone():
        ax.legend([plt.Line2D([], [], marker="o", color=INK, ls="-", ms=3.8, lw=0.8), plt.Line2D([], [], marker="o", color=POP["P-human"], ls="-", ms=3.8, lw=0.8)], ["P-all", "P-human"], loc="upper center", bbox_to_anchor=(0.4, -0.22), ncol=2, fontsize=FS - 0.8, handletextpad=0.3)
    legend_line(fig, {"": "Descriptive only, not graded: the indicator comes from the annotation tables that also feed the tiers (circular) and separates the tiers almost completely. GOLD against BRONZE, set-B tiers, 95% gene-cluster intervals; "
                f"P-all: {fmt_int(r['n_sites_gold'])} / {fmt_int(r['n_sites_bronze'])} sites, {fmt_int(r['n_gene_clusters'])} gene clusters."}, y=0.01)
    save(fig, "FigureS7", rec)


# ---------------------------------------------------------------------------- S8
def figS8():
    rec = Rec("FigureS8")
    c = load("figure_counts.json")
    dis = c["disorder_mean_by_tier"]
    fig = newfig((W, 3.6))
    ax = pax(fig, [0.085, 0.27, 0.895, 0.6], "A"); letter(ax, "A", dx=-34)
    groups = [("Source set A, human", "source_set_a_human"), ("Source set B, human", "source_set_b_human"), ("Source set B, mouse", "source_set_b_mouse"), ("Source set B, rat", "source_set_b_rat")]
    for i, (lab, k) in enumerate(groups):
        for j, t in enumerate(TIERS):
            x = i + (j - 1) * 0.25
            ax.bar(x, float(dis[k][t]["mean"]), width=0.23, color=TIER[t], edgecolor=INK, lw=0.4)
            ax.text(x, float(dis[k][t]["mean"]) + 0.01, f"{dis[k][t]['mean']:.3f}", ha="center", fontsize=FS - 0.8)
    ax.set_xticks(range(4)); ax.set_xticklabels([g[0] for g in groups]); ax.set_ylabel("Mean disorder score"); ax.set_ylim(0, 0.8)
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=TIER[t]) for t in TIERS], TIERS, ncol=3, loc="upper right")
    title(ax, "Disorder by tier and source set")
    src_box(ax, (0.02, 0.98), SRC_A + " " + SRC_B, 64)
    legend_line(fig, {"": "Means, no interval. Sites GOLD / SILVER / BRONZE: " + "; ".join(f"{lab[0].lower() + lab[1:]} " + " / ".join(fmt_int(dis[k][t]["n"]) for t in TIERS) for lab, k in groups) +
                ". Source set A: tiers from the human file with two source flags; source set B: three source flags."}, y=0.01)
    save(fig, "FigureS8", rec)


if __name__ == "__main__":
    setup()
    for f in (figS1, figS2, figS3, figS4, figS5, figS6, figS7, figS8):
        f()
        print(f.__name__, "written")
