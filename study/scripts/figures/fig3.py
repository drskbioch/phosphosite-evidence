"""Figure 3: the three decision-bearing signals under every adjustment (C2.1 to C2.3, C3.1, C4.5, C4.6, C5.1 to C5.6)."""
import matplotlib.pyplot as plt
import numpy as np

from lib import (FS, FS_TITLE, GREY, INK, POP, SIGNAL, W, Rec, ci_marker, fmt_int, legend_line, letter, load, load_dir, newfig, pax, save, setup, standalone, MODE)


def pick(d, prefix="", **cond):
    hits = [(k, v) for k, v in d.items() if isinstance(v, dict) and k.startswith(prefix) and all((v.get(c) == x) if not callable(x) else x(v.get(c)) for c, x in cond.items())]
    assert len(hits) == 1, (cond, [k for k, _ in hits])
    return hits[0]


def est_ci(rec):
    return float(rec["estimate"]), [float(x) for x in rec["ci"]]


def figure3():
    rec = Rec("Figure3")
    am, g23, g21, b4, g5, g18 = (load("amend_results.json"), load("g23_results.json"), load("g21_results.json"), load("b4_twoway.json"), load("g5_abundance.json"), load_dir("d2b"))
    sigs = ["S2", "S3", "S4"]
    rows = [("Crude", "crude"), ("CEM on the V50 decile (primary)", "primary"), ("Two-way interval (protein x study), P-all only", "twoway"),
            ("1:1 caliper matcher (S/T/Y records)", "d2b"), ("Non-primary PMIDs not counted", "g23"), ("Within gene, M1", "m1"), ("Within gene, M2 (reviewed sequences)", "m2")]
    ab_rows = [("Crude (D0)", "D0"), ("CEM on the V50 decile (D1)", "D1"), ("Index x abundance deciles (D4)", "D4"), ("Abundance deciles alone (D5)", "D5")]
    ny = len(rows) + len(ab_rows) + 1.6
    fig = newfig((W, 5.5))
    gs = fig.add_gridspec(1, 3, left=0.30, right=0.985, top=0.855, bottom=0.265, wspace=0.12)
    ys = {}
    y = ny
    for lab, key in rows:
        ys[key] = y
        y -= 1
    y -= 1.5
    for lab, key in ab_rows:
        ys["ab" + key] = y
        y -= 1
    # marker key
    h = [plt.Line2D([], [], marker="o", color=POP["P-all"], ls="-", ms=3.4, lw=0.8, label="P-all (filled circle)"),
         plt.Line2D([], [], marker="s", color=POP["P-human"], ls="-", ms=3.4, lw=0.8, label="P-human (filled square)"),
         plt.Line2D([], [], marker="o", color=POP["P-all"], mfc="white", ls="-", ms=3.4, lw=0.8, label="two-way interval (open circle, P-all only)"),
         plt.Line2D([], [], marker="D", color=POP["P-abund"], ls="-", ms=3.2, lw=0.8, label="P-abund (diamond)")]
    prim_ref = {}
    notes = []
    for ci_, s in enumerate(sigs):
        ax = pax(fig, gs[0, ci_], "ABC"[ci_])
        letter(ax, "ABC"[ci_], dx=0 if ci_ else -130, dy=33)
        ax.set_title(SIGNAL[s], fontsize=FS_TITLE, fontweight="bold", loc="left", pad=17 * MODE["fscale"])
        data = {}
        for pop in ("P-all", "P-human"):
            d0 = pick(am, "G17", signal=s, population=pop, design="D0")[1]
            d1k, d1 = pick(am, "G17", signal=s, population=pop, design="D1")
            data[pop] = {"crude": est_ci(d0), "primary": est_ci(d1)}
            data[pop]["d2b"] = est_ci(pick(g18, signal=s, population=f"{pop} (S/T/Y only)")[1])
            data[pop]["g23"] = est_ci(pick(g23, signal=s, population=pop, design="D1")[1])
            for m, key in (("M1", "m1"), ("M2", "m2")):
                r = pick(g21, signal=s, population=lambda p, pop=pop: p is not None and p.startswith(pop + " ("), model=m)[1]
                data[pop][key] = est_ci(r["companion_linear"] if "companion_linear" in r else r)
            data[pop]["att"] = d1["attenuation"]
            data[pop]["rec"] = d1
        tw = {"S2": "G1-0005", "S3": "G1-0008", "S4": "G1-0011"}[s]
        data["P-all"]["twoway"] = (float(b4[tw]["adjusted"]), [float(x) for x in b4[tw]["ci_two_way"]])
        for key in ("crude", "primary", "d2b", "g23", "m1", "m2", "twoway"):
            for pop, dy, mk in (("P-all", 0.16, "o"), ("P-human", -0.16, "s")):
                if key == "twoway" and pop == "P-human":
                    continue
                e, ci = data[pop][key]
                off = 0 if key == "twoway" else dy
                ci_marker(ax, ys[key] + off, e, ci, POP[pop], marker=mk, filled=(key != "twoway"), ms=3.4, cap=1.2)
        # abundance block: stored protein-cluster intervals
        for lab, d in ab_rows:
            r = pick(g5, signal=s, population="P-abund (human)", design=d)[1]
            assert r["resampling_unit"].startswith("protein")
            ci_marker(ax, ys["ab" + d], float(r["estimate"]), [float(x) for x in r["ci"]], POP["P-abund"], marker="D", ms=3.2, cap=1.2)
        ax.axvline(0, color=GREY, lw=0.5, zorder=0)
        ref = data["P-all"]["primary"][0]
        ax.axvline(ref, color=INK, lw=0.7, ls=(0, (1, 2)), zorder=0)
        ax.set_ylim(min(ys.values()) - 0.7, ny + 0.7)
        ax.set_yticks([ys[k] for _, k in rows] + [ys["ab" + k] for _, k in ab_rows])
        first = ci_ == 0 or standalone()
        ax.set_yticklabels([l for l, _ in rows] + [l for l, _ in ab_rows] if first else [])
        if not first:
            ax.tick_params(axis="y", length=0)
        if first:
            ax.annotate(xy=(0, ys["abD0"] + 1.05), xycoords=("axes fraction", "data"), xytext=(-6, 0), textcoords="offset points", text="Human proteins with an" + chr(10) + "abundance value (P-abund);" + chr(10) + "protein-cluster intervals",
                        fontsize=FS - 0.7, fontweight="bold", va="center", ha="right", linespacing=1.2, annotation_clip=False)
        ax.set_xlabel("Difference in mean disorder score" if s == "S3" else "Difference in proportion")
        ax.margins(x=0.06)
        if standalone():
            h2 = h + [plt.Line2D([], [], color=INK, lw=0.7, ls=(0, (1, 2)), label="primary P-all estimate (dotted)"), plt.Line2D([], [], color=GREY, lw=0.5, label="no difference (solid)")]
            ax.legend(handles=h2, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=FS - 0.7, handletextpad=0.3, columnspacing=1.0)
        att_all, att_h = data["P-all"]["att"], data["P-human"]["att"]
        ax.text(0.0, 1.012, f"Attenuation, crude to primary:"+chr(10)+f"P-all {att_all * 100:.1f}%, P-human {att_h * 100:.1f}%", transform=ax.transAxes, fontsize=FS - 0.7, va="bottom", ha="left", linespacing=1.2)
        r = data["P-all"]["rec"]
        notes.append(f"{SIGNAL[s].lower()} {fmt_int(r['n_sites_gold'])} / {fmt_int(r['n_sites_bronze'])} sites, {fmt_int(r['n_gene_clusters'])} gene clusters (P-human {fmt_int(data['P-human']['rec']['n_sites_gold'])} / {fmt_int(data['P-human']['rec']['n_sites_bronze'])}, {fmt_int(data['P-human']['rec']['n_gene_clusters'])})")
    if not standalone():
        fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.26, 0.175), ncol=4, fontsize=FS - 0.7, handletextpad=0.3, columnspacing=0.9)
    legend_line(fig, {"": f"Differences between GOLD and BRONZE sites, set-B tiers, 95% gene-cluster intervals ({fmt_int(data['P-all']['rec']['n_boot'])} resamples) unless marked. Dotted vertical line: primary P-all estimate; solid grey vertical line: no difference. "
                "Primary row, GOLD / BRONZE sites and gene clusters. No two-way interval is stored for P-human, so that marker is absent. "
                "Within-gene rows use the linear within-gene difference (S/T/Y records). The P-abund block (human proteins with an abundance value) shows the crude estimate and the adjusted designs with the stored protein-cluster intervals.",
                      "A": notes[0], "B": notes[1], "C": notes[2]}, y=0.005)
    save(fig, "Figure3", rec)


if __name__ == "__main__":
    setup()
    figure3()
    print("Figure 3 written")
