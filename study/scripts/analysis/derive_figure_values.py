"""Values the Rev 6 figures print that are plain arithmetic on stored result values (ratios, shares, logarithms). No new estimate and no interval.

`compute()` returns {name: {"value", "formula", "inputs"}}; `inputs` lists the result file and key (or contrast ID and field) each term comes from.
Running the file writes results/derived_figure_values.json; `check_bindings.py` calls `compute()` again and compares it with the stored file.
Rules reproduced from build_results.py (section 6 of RESULTS_REPORT.md): an adjusted estimate enters the spread only if its interval excludes the null; for the
within-definition change both the crude and the adjusted estimate must exclude the null; the 3-fold ("without the near-zero estimate") value drops every
estimate whose interval includes the null or has a limit within 0.01 of the null.
Usage: python derive_figure_values.py
"""
import json
import math
from pathlib import Path

CUR = Path(__file__).resolve().parents[2]
R6 = CUR / "analysis" / "results"
OUT = R6 / "derived_figure_values.json"
SIGNALS = {"S1": "domain overlap", "S2": "conservation", "S3": "disorder", "S4": "sequence motif"}
NULL = {"S1": 1.0, "S2": 0.0, "S3": 0.0, "S4": 0.0}


def jl(name):
    return json.load(open(R6 / name, encoding="utf-8"))


def excl_null(c, s):
    return bool(c[0] > NULL[s] or c[1] < NULL[s])


def eligible_spread(c, s):
    return excl_null(c, s) and min(abs(c[0] - NULL[s]), abs(c[1] - NULL[s])) >= 0.01


def tier_label(rule):
    for k, v in (("set-B", "set-B tiers"), ("set-A", "set-A tiers"), ("depth-only", "depth-only tiers"), ("agreement-only", "agreement-only")):
        if rule.startswith(k):
            return v
    raise ValueError(rule)


def compute():
    amend, cem, counts = jl("amend_results.json"), jl("cem_results.json"), jl("figure_counts.json")
    o = {}
    # ---- Fig 4B: tier-definition spread against the change caused by effort adjustment (G19; D1 records carry `estimate`, `crude` and both intervals)
    tiers = {}
    for k, v in amend.items():
        if k.startswith("G19") and v["design"] == "D1":
            tiers.setdefault((v["signal"], v["population"]), {})[tier_label(v["tier_rule"])] = (k, v)
    for (s, pop), defs in sorted(tiers.items()):
        adj = [v["estimate"] for _, v in defs.values()]
        okc = [(n, v) for n, (_, v) in defs.items() if excl_null(v["crude_ci"], s) and excl_null(v["ci"], s) and v["estimate"] > 0 and v["crude"] > 0]
        oka = [v["estimate"] for _, (_, v) in defs.items() if excl_null(v["ci"], s) and v["estimate"] > 0]
        el = [v["estimate"] for _, (_, v) in defs.items() if eligible_spread(v["ci"], s)]
        ids = sorted(k for k, _ in defs.values())
        base = f"{pop}|{SIGNALS[s]}"
        inp = [["analysis/results/amend_results.json", f"{i}.estimate, .crude, .ci, .crude_ci"] for i in ids]
        if len(oka) >= 2 and okc:
            o[f"fig4b|{base}|ln_change"] = {"value": max(abs(math.log(v["estimate"] / v["crude"])) for _, v in okc), "inputs": inp,
                                            "formula": "max over tier definitions of abs(ln(adjusted / crude)); definitions where the crude and the adjusted interval both exclude the null"}
            o[f"fig4b|{base}|ln_spread"] = {"value": math.log(max(oka) / min(oka)), "inputs": inp,
                                            "formula": "ln(largest adjusted / smallest adjusted) over tier definitions whose adjusted interval excludes the null"}
            o[f"fig4b|{base}|n_definitions"] = {"value": len(oka), "inputs": inp, "formula": "number of tier definitions whose adjusted interval excludes the null"}
        o[f"fig4b|{base}|spread_all"] = {"value": max(adj) / min(adj), "inputs": inp, "formula": "largest adjusted / smallest adjusted over all tier definitions"}
        if len(el) >= 2:
            o[f"fig4b|{base}|spread_eligible"] = {"value": max(el) / min(el), "inputs": inp,
                                                  "formula": "largest / smallest adjusted estimate after dropping any estimate whose interval includes the null or has a limit within 0.01 of the null"}
    # ---- Fig S4: gene-cluster against accession-cluster interval width (G17 against G1/G2 and G16-free rows)
    for pop_g17, blk in (("P-all", "G1"), ("P-human", "G2")):
        for s in ("S2", "S3", "S4", "S5"):
            g17 = next(k for k, v in amend.items() if k.startswith("G17") and v["signal"] == s and v["population"] == pop_g17 and v["design"] == "D1")
            g1 = next(k for k, v in cem.items() if k.startswith(blk) and v["signal"] == s and v["design"] == "D1" and v["population"] == pop_g17)
            wg, wp = amend[g17]["ci"][1] - amend[g17]["ci"][0], cem[g1]["ci"][1] - cem[g1]["ci"][0]
            o[f"figS4|{pop_g17}|{s}|width_ratio_gene_over_accession"] = {
                "value": wg / wp, "formula": "(upper - lower of the gene-cluster interval) / (upper - lower of the accession-cluster interval), same contrast, D1",
                "inputs": [["analysis/results/amend_results.json", f"{g17}.ci"], ["analysis/results/cem_results.json", f"{g1}.ci"]]}
    g1 = next(k for k, v in cem.items() if k.startswith("G1") and v["signal"] == "S1" and v["design"] == "D1")
    g17 = next(k for k, v in amend.items() if k.startswith("G17") and v["signal"] == "S1" and v["population"] == "P-all (human)" and v["design"] == "D1")
    o["figS4|P-all (human)|S1|width_ratio_gene_over_accession"] = {
        "value": (amend[g17]["ci"][1] - amend[g17]["ci"][0]) / (cem[g1]["ci"][1] - cem[g1]["ci"][0]),
        "formula": "(upper - lower of the gene-cluster interval) / (upper - lower of the accession-cluster interval), same contrast, D1",
        "inputs": [["analysis/results/amend_results.json", f"{g17}.ci"], ["analysis/results/cem_results.json", f"{g1}.ci"]]}
    # ---- shares of counts (figure_counts.json)
    f = "analysis/results/figure_counts.json"

    def share(name, num, den, nk, dk, what):
        o[name] = {"value": num / den, "formula": what, "inputs": [[f, nk], [f, dk]]}

    for sp in ("human", "mouse", "rat"):
        n = counts["n_sites_by_species"][sp]
        for t, v in counts["tier_counts_by_species"][sp].items():
            share(f"fig2a|{sp}|{t}", v, n, f"tier_counts_by_species.{sp}.{t}", f"n_sites_by_species.{sp}", "tier count / sites of the species")
    n_all = sum(counts["n_sites_by_species"].values())
    for t, v in counts["tier_counts_all"].items():
        o[f"fig2a|all|{t}"] = {"value": v / n_all, "formula": "tier count / all sites", "inputs": [[f, f"tier_counts_all.{t}"], [f, "n_sites_by_species (sum)"]]}
    nh = counts["human_n_sites"]
    for key, nm in (("human_tier_pmid_count_rule", "pmid_rule_set_a"), ("human_tier_study_cluster_rule_set_a", "cluster_rule_set_a"), ("human_tier_study_cluster_rule_set_b", "cluster_rule_set_b")):
        for t, v in counts[key].items():
            share(f"fig2b|{nm}|{t}", v, nh, f"{key}.{t}", "human_n_sites", "tier count / human sites")
    for t, v in counts["das_counts_all"].items():
        o[f"figS1|das|{t}"] = {"value": v / n_all, "formula": "sites with agreement score / all sites", "inputs": [[f, f"das_counts_all.{t}"], [f, "n_sites_by_species (sum)"]]}
    p = counts["pmid"]
    for k in ("resolved", "valid_format_not_found", "invalid_format_over_45m"):
        share(f"figS1|pmid|{k}", p[k], p["distinct"], f"pmid.{k}", "pmid.distinct", "identifiers in the class / distinct identifiers")
    for sp, v in counts["noncanonical_by_species"].items():
        share(f"figS1|noncanonical_rate|{sp}", v, counts["n_sites_by_species"][sp], f"noncanonical_by_species.{sp}", f"n_sites_by_species.{sp}", "non-canonical residue sites / sites of the species")
    share("figS1|noncanonical_rate|all", counts["noncanonical_n"], n_all, "noncanonical_n", "n_sites_by_species (sum)", "non-canonical residue sites / all sites")
    cs, ts = counts["cluster_split"], counts["temporal_split"]
    share("figS5|cluster_split|proteins_in_both", cs["proteins_in_both"], cs["proteins_total"], "cluster_split.proteins_in_both", "cluster_split.proteins_total", "proteins in both partitions / proteins")
    share("figS5|temporal_split|proteins_in_both", ts["proteins_in_both"], ts["proteins_total"], "temporal_split.proteins_in_both", "temporal_split.proteins_total", "proteins in both partitions / proteins")
    share("figS5|temporal_split|clusters_in_both", ts["study_clusters_in_both"], ts["study_clusters_total"], "temporal_split.study_clusters_in_both", "temporal_split.study_clusters_total", "study clusters in both partitions / study clusters")
    # ---- Fig 2D: GOLD share lost to SILVER and the share of records with a non-primary PMID are stored; the PSP paper share is a ratio
    g23 = jl("g23_tiers_summary.json")
    psp = next(x for x in g23["nonprimary"] if x["pmid"] == "25514926")
    o["fig2d|psp_paper_share_of_records"] = {"value": psp["records_linked"] / g23["n_records"], "formula": "records linked to PMID 25514926 / all records",
                                             "inputs": [["analysis/results/g23_tiers_summary.json", "nonprimary[pmid=25514926].records_linked"], ["analysis/results/g23_tiers_summary.json", "n_records"]]}
    o["fig2d|gold_to_silver_share_of_gold"] = {"value": g23["transition_table"]["GOLD"]["SILVER"] / g23["tier_counts_old"]["GOLD"], "formula": "GOLD records that fall to SILVER / GOLD records",
                                               "inputs": [["analysis/results/g23_tiers_summary.json", "transition_table.GOLD.SILVER"], ["analysis/results/g23_tiers_summary.json", "tier_counts_old.GOLD"]]}
    return o


if __name__ == "__main__":
    d = compute()
    OUT.write_text(json.dumps(d, indent=1), encoding="utf-8")
    print(len(d), "values written to", OUT)
    for k in sorted(d):
        if k.startswith("fig4b"):
            print(f"{k:60s} {d[k]['value']:.4f}")
