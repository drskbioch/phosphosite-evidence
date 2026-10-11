"""Figure 4 (what the signals depend on) and Figure 5 (predictor audit, secondary)."""
import sys

import matplotlib.pyplot as plt
import numpy as np

from lib import (CUR, FS, FS_TITLE, GREY, INK, LIGHT, MODE, POP, SIGNAL, TIER, TOOL, W, Rec, ci_marker, fmt_int, legend_line, letter, load, newfig, pax, save, setup, standalone, title)

sys.path.insert(0, str(CUR / "scripts" / "analysis"))
from derive_figure_values import tier_label  # noqa: E402

DEFS = ["set-B tiers", "set-A tiers", "depth-only tiers", "agreement-only"]
DEF_LAB = {"set-B tiers": "set-B (primary)", "set-A tiers": "set-A", "depth-only tiers": "depth only", "agreement-only": "agreement only"}
NL = chr(10)


def ci_(r):
    return [float(x) for x in r["ci"]]


def figure4():
    rec = Rec("Figure4")
    am, d = load("amend_results.json"), load("derived_figure_values.json")
    tiers = {}
    for k, v in am.items():
        if k.startswith("G19") and v["design"] == "D1":
            tiers.setdefault((v["signal"], v["population"]), {})[tier_label(v["tier_rule"])] = v
    H = 6.8  # composite height in inches (489.6 pt); positions below are given in inches from the top and divided by H
    fig = newfig((W, H))
    gs = fig.add_gridspec(2, 4, **(dict(left=0.17, right=0.93, top=0.89, bottom=0.14, wspace=0.75, hspace=0.8) if standalone() else dict(left=0.115, right=0.99, top=1 - 0.30 / H, bottom=1 - 3.30 / H, wspace=0.28, hspace=0.46)))
    gb = fig.add_gridspec(1, 1, left=0.085, right=0.40, top=1 - 3.95 / H, bottom=1 - 5.35 / H)
    gc = fig.add_gridspec(1, 1, left=0.665, right=0.985, top=1 - 3.95 / H, bottom=1 - 5.35 / H)
    cells = [("P-human", ["S1", "S2", "S3", "S4"]), ("P-all", [None, "S2", "S3", "S4"])]
    for ri, (pop, sigs) in enumerate(cells):
        for ci_idx, s in enumerate(sigs):
            ax = pax(fig, gs[ri, ci_idx], "A", multi=True)
            if s is None:
                ax.axis("off")
                if standalone() and MODE["letter"] == "A":
                    MODE["pfig"].legend([plt.Line2D([], [], marker="o", color=POP["P-human"], mfc="white", ls="none", ms=3.4), plt.Line2D([], [], marker="o", color=POP["P-human"], ls="none", ms=3.4),
                                         plt.Line2D([], [], marker="o", color=POP["P-all"], mfc="white", ls="none", ms=3.4), plt.Line2D([], [], marker="o", color=POP["P-all"], ls="none", ms=3.4)],
                                        ["P-human, crude", "P-human, CEM on the V50 decile", "P-all, crude", "P-all, CEM on the V50 decile"], loc="lower center", ncol=4, fontsize=FS - 0.8, handletextpad=0.2, columnspacing=1.2)
                    MODE["notight"] = True
                continue
            if ri == 0 and ci_idx == 0:
                letter(ax, "A", dx=-54, dy=12)
            defs = [t for t in DEFS if t in tiers[(s, pop)]]
            ys = list(range(len(defs)))[::-1]
            col = POP[pop]
            for y, t in zip(ys, defs):
                r = tiers[(s, pop)][t]
                ci_marker(ax, y + 0.17, float(r["crude"]), [float(x) for x in r["crude_ci"]], col, filled=False, ms=3.4, cap=1.2)
                ci_marker(ax, y - 0.17, float(r["estimate"]), ci_(r), col, filled=True, ms=3.4, cap=1.2)
                if s == "S3" and t == "agreement-only":
                    ax.annotate(f"crude {r['crude']:.4f}" + NL + f"({r['crude_ci'][0]:.4f} to {r['crude_ci'][1]:.4f})" + NL + f"adjusted {r['estimate']:.4f}" + NL + f"({r['ci'][0]:.4f} to {r['ci'][1]:.4f})",
                                xy=(0.30, y), xycoords=("axes fraction", "data"), fontsize=FS - 1.4, ha="left", va="center", linespacing=1.1)
            if s == "S1":
                ax.set_xscale("log")
                ax.axvline(1, color=GREY, lw=0.5, ls=":")
            else:
                ax.axvline(0, color=GREY, lw=0.5, ls=":")
            ax.set_ylim(-0.7, len(defs) - 0.1 + 0.55)
            ax.set_yticks(ys)
            first_in_row = (ri == 0 and ci_idx == 0) or (ri == 1 and ci_idx == 1)
            ax.set_yticklabels([DEF_LAB[t] for t in defs] if first_in_row else [])
            if not first_in_row:
                ax.tick_params(axis="y", length=0)
            ax.set_title(f"{SIGNAL[s]}, {pop}", fontsize=FS - 0.2, fontweight="bold", loc="left", pad=3)
            ax.set_xlabel("Odds ratio" if s == "S1" else ("Difference in mean score" if s == "S3" else "Difference in proportion"), fontsize=FS - 0.5)
            ax.margins(x=0.1)
    if standalone() and MODE["letter"] == "A":
        MODE["pfig"].suptitle("GOLD against BRONZE by tier definition", fontsize=FS_TITLE, fontweight="bold", x=0.01, ha="left")
    # B scatter
    axb = pax(fig, gb[0, 0], "B"); letter(axb, "B", dx=-38, dy=16)
    if standalone():
        title(axb, "Spread across tier definitions against change from adjustment")
    pts = []
    for pop in ("P-human", "P-all"):
        for s in ("S1", "S2", "S3", "S4"):
            k = f"fig4b|{pop}|{SIGNAL[s].lower()}|"
            if k + "ln_change" in d:
                pts.append((pop, s, d[k + "ln_change"]["value"], d[k + "ln_spread"]["value"], d[k + "spread_all"]["value"], d.get(k + "spread_eligible", {}).get("value")))
    lim = 3.15
    axb.plot([0, lim], [0, lim], color=GREY, lw=0.6, ls="--", zorder=0)
    axb.text(0.40, 0.46, "y = x", fontsize=FS - 0.5, color=GREY)
    for pop, s, x, y, sp_all, sp_el in pts:
        axb.plot(float(x), float(y), "o", color=POP[pop], ms=4.3, zorder=3)
        lab = f"{SIGNAL[s]}, {pop}"
        off = {("P-human", "S2"): (0.012, 0.08, "left", "bottom"), ("P-all", "S2"): (0.012, 0.0, "left", "center"), ("P-human", "S1"): (0.012, 0.08, "left", "bottom"),
               ("P-all", "S4"): (0.012, 0.0, "left", "center"), ("P-human", "S4"): (-0.012, 0.08, "right", "bottom"), ("P-human", "S3"): (0.012, -0.05, "left", "top")}
        dx, dy, ha, va = off.get((pop, s), (0.012, 0.075, "left", "bottom"))
        if pop == "P-all" and s == "S3":
            lab = f"{SIGNAL[s]}, {pop}:" + NL + f"{sp_all:.0f}-fold; {sp_el:.0f}-fold without the near-zero" + NL + "agreement-only estimate"
            dx, dy, ha, va = 0.012, 0.1, "left", "bottom"
        axb.text(float(x) + dx, float(y) + dy, lab, fontsize=FS - 1.0, ha=ha, va=va, linespacing=1.15)
    axb.set_xlim(0, 0.46); axb.set_ylim(0, 3.95)
    axb.set_xlabel("Change from crude to adjusted within a tier definition (natural log)", fontsize=FS - 0.3)
    axb.set_ylabel("Spread across tier definitions (natural log)", fontsize=FS - 0.3)
    sx = axb.secondary_xaxis("top", functions=(np.exp, np.log)); sx.set_xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5]); sx.tick_params(labelsize=FS - 1)
    sy = axb.secondary_yaxis("right", functions=(np.exp, np.log)); sy.set_yticks([1, 2, 3, 5, 10, 20]); sy.tick_params(labelsize=FS - 1)
    axb.legend([plt.Line2D([], [], marker="o", color=POP["P-all"], ls="none", ms=4.3), plt.Line2D([], [], marker="o", color=POP["P-human"], ls="none", ms=4.3)], ["P-all", "P-human"],
               loc="center right", bbox_to_anchor=(1.0, 0.55), fontsize=FS - 0.8, handletextpad=0.2, frameon=False)
    # C domain overlap across populations
    axc = pax(fig, gc[0, 0], "C"); letter(axc, "C", dx=-122, dy=16)
    pop_rows = [("All human records", "G17-0379", "G17-0380", INK), ("Proteins with an abundance value", "G17-0399", "G17-0400", INK),
                ("Proteins with a domain-annotated site", "G17-0397", "G17-0398", INK), ("Reviewed canonical accessions", "G16-0361", "G16-0362", INK),
                ("All human records, set-A tiers", "G19-0420", "G19-0421", INK)]
    ys = list(range(len(pop_rows)))[::-1]
    notes = []
    for y, (lab, k0, k1, col) in zip(ys, pop_rows):
        r0, r1 = am[k0], am[k1]
        ci_marker(axc, y + 0.17, float(r0["estimate"]), ci_(r0), col, filled=False, ms=3.6, cap=1.2)
        ci_marker(axc, y - 0.17, float(r1["estimate"]), ci_(r1), col, filled=True, ms=3.6, cap=1.2)
        axc.text(max(float(r0["ci"][1]), float(r1["ci"][1])) * 1.06, y, f"crude {r0['estimate']:.2f}, adjusted {r1['estimate']:.2f}", fontsize=FS - 0.8, va="center")
        notes.append(f"{lab[0].lower() + lab[1:]} {fmt_int(r1['n_sites_gold'])} / {fmt_int(r1['n_sites_bronze'])} ({fmt_int(r1['n_gene_clusters'])} gene clusters)")
    axc.set_xscale("log"); axc.axvline(1, color=GREY, lw=0.5, ls=":")
    axc.set_yticks(ys); axc.set_yticklabels([r[0] for r in pop_rows], fontsize=FS - 0.5)
    axc.set_xlim(0.45, 70); axc.set_ylim(-0.7, len(pop_rows) - 0.3)
    axc.set_xlabel("Domain overlap, GOLD against BRONZE (odds ratio)", fontsize=FS - 0.3)
    title(axc, "Domain overlap across populations")
    if standalone():
        axc.legend([plt.Line2D([], [], marker="o", color=INK, mfc="white", ls="none", ms=3.6), plt.Line2D([], [], marker="o", color=INK, ls="none", ms=3.6)], ["crude", "adjusted (CEM on the V50 decile)"], loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2, fontsize=FS - 0.8, handletextpad=0.2)
    # legend line
    fm = {k: d[k]["formula"] for k in ("fig4b|P-all|disorder|ln_change", "fig4b|P-all|disorder|ln_spread", "fig4b|P-all|disorder|spread_eligible")}
    ns = []
    for pop, sigs in cells:
        for s in sigs:
            if s is None:
                continue
            r = tiers[(s, pop)]["set-B tiers"]
            ns.append(f"{pop} {SIGNAL[s].lower()} {fmt_int(r['n_sites_gold'])} / {fmt_int(r['n_sites_bronze'])} ({fmt_int(r['n_gene_clusters'])})")
    legend_line(fig, {"A": "GOLD against BRONZE; open marker crude, filled marker CEM on the V50 decile, 95% gene-cluster intervals. GOLD / BRONZE sites (gene clusters) under set-B tiers: " + "; ".join(ns) + "; all other cells are listed in Table S5.",
                      "B": "x = " + fm["fig4b|P-all|disorder|ln_change"] + "; y = " + fm["fig4b|P-all|disorder|ln_spread"] + "; second value of the labelled point: " + fm["fig4b|P-all|disorder|spread_eligible"] + ". Secondary axes: the same values as plain ratios.",
                      "C": "open marker crude, filled adjusted; GOLD / BRONZE sites (gene clusters): " + "; ".join(notes) + "; set-B tiers except where stated."}, y=0.003, size=FS - 0.8)
    save(fig, "Figure4", rec)


def figure5():
    rec = Rec("Figure5")
    g11, g23 = load("g11_predictors.json"), load("g23_tiers_summary.json")
    fig = newfig((W, 6.6))
    gs = fig.add_gridspec(2, 1, left=0.34, right=0.97, top=0.93, bottom=0.255, hspace=0.62, height_ratios=[1.1, 1])
    # A: gap by study effort
    ax = pax(fig, gs[0, 0], "A"); letter(ax, "A", dx=-165, dy=14)
    margin = float(g11["settings"]["margin"])
    ax.axvspan(-margin, margin, color="#ececec", zorder=0)
    ax.axvline(0, color=GREY, lw=0.5)
    ax.axvline(margin, color=INK, lw=0.6, ls="--")
    tools = [("DeepPhos", "DeepPhos|planned: S/T sites only", ["unrestricted", "restriction_year=2019", "restriction_year=2018"]),
             ("MusiteDeep", "MusiteDeep|planned: S/T sites only", ["unrestricted", "restriction_year=2020", "restriction_year=2017"]),
             ("PhosF3C", "PhosF3C|planned: Y sites only", ["unrestricted", "restriction_year=2024"])]
    y = 9.0
    ylabs, yt, notes = [], [], []
    for tool, key, sets in tools:
        blk = g11[key]
        for st in sets:
            r = blk[st]
            if r["status"] != "evaluated":
                ax.text(-0.095, y, "not evaluated: too few positives", fontsize=FS - 0.8, va="center", ha="left", color=GREY)
                ylabs.append(f"{tool}, positives after {st.split('=')[1]}")
            else:
                restricted = st != "unrestricted"
                ci = [float(x) for x in (r["gene_ci90"] if restricted else r["gene_ci95"])]
                ci_marker(ax, y, float(r["gap"]), ci, TOOL[tool], marker="o", filled=not restricted, ms=3.8, cap=1.5)
                ax.text(0.132, y, f"{r['gap']:.4f}", fontsize=FS - 0.8, va="center")
                ylabs.append(f"{tool}, " + ("all positives" if not restricted else f"positives after {st.split('=')[1]}"))
                notes.append(f"{tool} {'all' if not restricted else st.split('=')[1]}: {fmt_int(r['n_sites'])} sites, {fmt_int(r['n_positives'])} positives, {fmt_int(r['n_gene_clusters'])} gene clusters")
            yt.append(y)
            y -= 1
        y -= 0.6
    ax.set_yticks(yt); ax.set_yticklabels(ylabs, fontsize=FS - 0.6)
    ax.set_ylim(min(yt) - 0.8, 9.6)
    ax.set_xlim(-0.1, 0.17)
    ax.set_xlabel("AUROC gap: deciles 8 to 10 minus deciles 1 to 3" + NL + "(deciles of the original index; evaluation sample set before scoring)", fontsize=FS - 0.5)
    ax.text(margin + 0.003, min(yt) - 0.45, f"{margin:.2f}: criterion set before scoring", fontsize=FS - 0.8, ha="left", va="center")
    title(ax, "Predictor gap by study effort")
    if standalone():
        ax.legend([plt.Line2D([], [], marker="o", color=GREY, ls="-", ms=3.8, lw=0.8), plt.Line2D([], [], marker="o", color=GREY, mfc="white", ls="-", ms=3.8, lw=0.8), plt.Rectangle((0, 0), 1, 1, fc="#ececec", ec="none"), plt.Line2D([], [], color=INK, lw=0.6, ls="--")],
                  ["all positives (95% gene-cluster interval)", "positives first reported after the year (90% interval)", "equivalence margin", "criterion set before scoring"], loc="upper center", bbox_to_anchor=(0.35, -0.3), ncol=2, fontsize=FS - 0.8, handletextpad=0.3, columnspacing=1.0)
    ax.title.set_position((0.0, 1.03))
    # B: GOLD sites kept
    axb = pax(fig, gs[1, 0], "B"); letter(axb, "B", dx=-165, dy=14)
    ho = g23["held_out"]
    designs = [("gold_total", "All records"), ("cluster_disjoint_test", "Cluster-disjoint test set"), ("cluster_and_protein_disjoint_test", "Cluster- and protein-disjoint test set"),
               ("temporal_test_all", "Temporal test set"), ("temporal_protein_disjoint", "Temporal, protein-disjoint"), ("temporal_study_disjoint", "Temporal, study-disjoint"),
               ("temporal_protein_and_study_disjoint", "Temporal, protein- and study-disjoint")]
    gold = lambda k, which: ho[k][which] if k == "gold_total" else ho[k][which].get("GOLD")
    omitted = [(lab, ho[k]["n"]) for k, lab in designs if gold(k, "old") is None and gold(k, "new") is None]
    designs = [(k, lab) for k, lab in designs if not (gold(k, "old") is None and gold(k, "new") is None)]
    yy = list(range(len(designs)))[::-1]
    for y_, (k, lab) in zip(yy, designs):
        blk = ho[k]
        old, new = gold(k, "old"), gold(k, "new")
        for off, v, al, hat in ((0.2, old, 1.0, None), (-0.2, new, 0.45, "///")):
            if v is None:
                continue
            axb.barh(y_ + off, float(v), height=0.34, color=TIER["GOLD"], edgecolor=INK, lw=0.4, alpha=al, hatch=hat)
            axb.text(float(v) * 1.15 + 0.3, y_ + off, fmt_int(v), fontsize=FS - 0.8, va="center")
    axb.set_xscale("symlog", linthresh=1); axb.set_xlim(0, 3e6)
    axb.set_yticks(yy); axb.set_yticklabels([d[1] for d in designs], fontsize=FS - 0.6); axb.set_ylim(-0.75, len(designs) - 0.3)
    axb.set_xlabel("GOLD sites (log scale)", fontsize=FS - 0.5)
    axb.legend([plt.Rectangle((0, 0), 1, 1, fc=TIER["GOLD"], ec=INK), plt.Rectangle((0, 0), 1, 1, fc=TIER["GOLD"], ec=INK, alpha=0.45, hatch="///")], ["Primary tiers", "Non-primary PMIDs not counted"],
               loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=FS - 0.7, handlelength=1.0)
    title(axb, "GOLD sites kept by each held-out design")
    prim_note = "primary restriction year per tool: " + ", ".join(f"{t} {sets[1].split('=')[1]}" for t, _, sets in tools) + " (year of the first publication describing the version scored); the other restriction year shown for DeepPhos and MusiteDeep is a sensitivity analysis"
    legend_line(fig, {"": "Secondary analysis.",
                      "A": "serine/threonine sites for DeepPhos and MusiteDeep, tyrosine sites for PhosF3C; filled marker with the 95% gene-cluster interval for all positives, open marker with the 90% interval for positives first reported after the stated year; "
                           + prim_note + "; grey band: equivalence margin; sites, positives and gene clusters: " + "; ".join(notes) + ".",
                      "B": "set-B tiers from the study-cluster rule and tiers with the non-primary PMIDs not counted; test sets as assigned before scoring. Designs with no GOLD site are not drawn: " + "; ".join(f"{lab.lower()} (test set of {fmt_int(n)} sites)" for lab, n in omitted) + "."}, y=0.004, size=FS - 1.1)
    save(fig, "Figure5", rec)


if __name__ == "__main__":
    setup()
    figure4()
    figure5()
    print("Figures 4 and 5 written")
