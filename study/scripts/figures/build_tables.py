"""Tables 1 to 3 and Tables S1 to S6 of the Rev 6 rebuild. Every number is read from a Rev 6 result file through `provenance.jload`; text cells come from
analysis/sources/SOURCES.md (licences), CLAIMS_VS_EVIDENCE.md (grades, Rev 5 to Rev 6 table) or a boolean field, and are logged for check_bindings.py.
Writes CSV and Markdown files to study/manuscript/tables/; the table of changes between the two analyses is internal and goes to study/analysis/internal/.
Footnotes are written below the table in the Markdown file and, for the tables that have one, to <table>_notes.md.
"""
import csv
import json
import re
import shutil

import pandas as pd

from lib import BIND, CUR, R6, TABDIR, Rec, fmt_int, load, load_dir, tx

INTERNAL = R6.parent / "internal"
WORDING = json.load(open(R6.parent / "table_wording.json", encoding="utf-8")) if (R6.parent / "table_wording.json").exists() else {}

NL = chr(10)
SOURCES = (R6.parent / "sources" / "SOURCES.md").read_text(encoding="utf-8")
CLAIMS = (R6.parent / "CLAIMS_VS_EVIDENCE.md").read_text(encoding="utf-8")


def md_table_rows(text, header_start):
    """Rows of the markdown table that follows the first line starting with header_start."""
    lines = text.split(chr(10))
    i = next(k for k, l in enumerate(lines) if l.startswith(header_start))
    rows = []
    for l in lines[i + 2:]:
        if not l.startswith("|"):
            break
        rows.append([c.strip() for c in l.strip().rstrip(chr(13)).strip("|").split(" | ")])
    return rows


LIC = {r[0]: r for r in md_table_rows(SOURCES, "| key | URL")}


def jtx(val, rel, key):
    return tx(str(val), rel, key)


def licence(key):
    r = LIC[key]
    return tx(r[4], "analysis/sources/SOURCES.md", f"section F, key {key}, cell text")


def claim_row(cid):
    for l in CLAIMS.split(chr(10)):
        if l.startswith(f"| {cid} |"):
            return [c.strip() for c in l.rstrip(chr(13)).strip("|").split(" | ")]
    raise KeyError(cid)


def grade(cid):
    g = claim_row(cid)[6]
    g = g.split(" (")[0]
    return tx(g, "analysis/CLAIMS_VS_EVIDENCE.md", f"{cid}, Strength column (text before the first parenthesis)")


def write(df, name, notes=None, directory=None):
    """CSV and Markdown file of a table; `notes` (list of strings) are the footnotes: appended to the Markdown file and written to <name>_notes.md."""
    directory = directory or TABDIR
    directory.mkdir(parents=True, exist_ok=True)
    df.to_csv(directory / f"{name}.csv", index=False, encoding="utf-8")
    md = df.to_markdown(index=False)
    if notes:
        text = NL.join(notes)
        md += NL + NL + "Notes. " + text.replace(NL, NL + NL)
        (directory / f"{name}_notes.md").write_text(text + NL, encoding="utf-8")
    (directory / f"{name}.md").write_text(md, encoding="utf-8")


def fd(x, nd):
    """x at nd decimals; a value that would print as zero but is not zero gets decimals until a digit shows (no '-0.0000')."""
    s = format(x, f".{nd}f")
    k = nd
    while float(s) == 0 and float(x) != 0 and k < 8:
        k += 1
        s = format(x, f".{k}f")
    return s


def sig4(x):
    return float(f"{float(x):.4g}")


def pick(d, prefix="", **cond):
    hits = [(k, v) for k, v in d.items() if isinstance(v, dict) and k.startswith(prefix) and all(v.get(c) == x for c, x in cond.items())]
    assert len(hits) == 1, (cond, [k for k, _ in hits])
    return hits[0]


def ec(r, nd, ci=None):
    ci = ci if ci is not None else r["ci"]
    return f"{r['estimate']:.{nd}f} ({ci[0]:.{nd}f} to {ci[1]:.{nd}f})"


def sites(r):
    return f"{fmt_int(r['n_sites_gold'])} / {fmt_int(r['n_sites_bronze'])}; {fmt_int(r['n_gene_clusters'])}"


# ------------------------------------------------------------------ Table 1
def table1():
    rec = Rec("Table1")
    g23, gu, ds, ab, g11 = load("g23_tiers_summary.json"), load("gene_unit_summary.json"), load("downloads_manifest.json"), load("e_abundance_coverage.json"), load("g11_predictors.json")
    d1, b5, d3, cov = load("d1_gene2pubmed_summary.json"), load("b5_index_effect_sizes.json"), load("d3_d4_d5_index_checks.json"), load("figure_counts.json")
    tv = load("../../earlier_results/tool_versions.json") if False else None
    import provenance as P
    tv = P.jload(CUR / "earlier_results" / "tool_versions.json")
    ms = list(csv.DictReader(open(R6 / "manifest_status.csv", encoding="utf-8")))
    psp_note = next(r["result_status"] for r in ms if r["block"] == "G9" and r["result_status"].startswith("not run"))
    pa_note = next(r["result_status"] for r in ms if r["block"] == "G9" and r["result_status"].startswith("not feasible"))
    rows = [
        ["EPSD 2.0 (site records)", tx("2.0", "analysis/sources/SOURCES.md", "section F, EPSD provenance, version"), tx("on or before 2026-09-13", "analysis/sources/SOURCES.md", "section F, EPSD provenance, download date"), licence("epsd"),
         f"{fmt_int(g23['n_records'])} sites on {fmt_int(gu['all']['accessions'])} proteins",
         f"human {fmt_int(cov['n_sites_by_species']['human'])}, mouse {fmt_int(cov['n_sites_by_species']['mouse'])}, rat {fmt_int(cov['n_sites_by_species']['rat'])} sites"],
        ["NCBI gene2pubmed", tx("gene2pubmed.gz from the NCBI Gene FTP", "analysis/sources/SOURCES.md", "section F, gene2pubmed provenance"), tx("2026-09-24", "analysis/sources/SOURCES.md", "section F, gene2pubmed provenance, retrieved date"), licence("gene2pubmed"),
         f"{fmt_int(d1['pairs_3species'])} gene-PMID pairs for the three species",
         f"V50 index defined for {fmt_int(b5['n_proteins'])} of {fmt_int(gu['all']['accessions'])} proteins"],
        ["UniProt ID mapping", jtx(ds["uniprot_idmapping_HUMAN_9606"]["version"], "analysis/results/downloads_manifest.json", "uniprot_idmapping_HUMAN_9606.version"), jtx(ds["uniprot_idmapping_HUMAN_9606"]["retrieved"], "analysis/results/downloads_manifest.json", "uniprot_idmapping_HUMAN_9606.retrieved"), licence("uniprot"),
         f"human, mouse and rat ID-mapping files",
         f"human accessions with a gene: {fmt_int(gu['human']['accessions_with_gene'])} of {fmt_int(gu['human']['accessions'])}"],
        ["PaxDb", jtx(ds["paxdb_9606"]["version"], "analysis/results/downloads_manifest.json", "paxdb_9606.version"), jtx(ds["paxdb_9606"]["retrieved"], "analysis/results/downloads_manifest.json", "paxdb_9606.retrieved"), licence("paxdb"),
         f"{fmt_int(ab['human']['paxdb_entries'])} human entries",
         f"human proteins with abundance {fmt_int(ab['human']['proteins_with_abundance'])} of {fmt_int(ab['human']['proteins_in_site_table'])}; records {fmt_int(ab['human']['records_with_abundance'])} of {fmt_int(ab['human']['records_in_site_table'])} ({ab['human']['share_of_records'] * 100:.1f}%)"],
        ["Pharos (TCRD)", jtx(ds["pharos_tdl"]["version"], "analysis/results/downloads_manifest.json", "pharos_tdl.version"), jtx(ds["pharos_tdl"]["retrieved"], "analysis/results/downloads_manifest.json", "pharos_tdl.retrieved"), licence("pharos"),
         f"{fmt_int(ds['pharos_tdl']['rows'])} targets",
         f"Tclin {fmt_int(ds['pharos_tdl']['tdl_counts']['Tclin'])}, Tchem {fmt_int(ds['pharos_tdl']['tdl_counts']['Tchem'])}, Tbio {fmt_int(ds['pharos_tdl']['tdl_counts']['Tbio'])}, Tdark {fmt_int(ds['pharos_tdl']['tdl_counts']['Tdark'])}"],
    ]
    for tool, key, lic in (("DeepPhos", "DeepPhos|planned: S/T sites only", "deepphos"), ("MusiteDeep", "MusiteDeep|planned: S/T sites only", "musitedeep"), ("PhosF3C", "PhosF3C|planned: Y sites only", "phosf3c")):
        u = g11[key]["unrestricted"]
        ver = jtx(tv[tool]["scored_checkout"], "earlier_results/tool_versions.json", f"{tool}.scored_checkout") if tool != "PhosF3C" else jtx(tv[tool]["repository"], "earlier_results/tool_versions.json", "PhosF3C.repository") + "; " + jtx(tv[tool]["scored_weights"], "earlier_results/tool_versions.json", "PhosF3C.scored_weights")
        rows.append([f"{tool} (predictor)", ver, "not recorded", licence(lic),
                     f"{fmt_int(u['n_sites'])} sites scored", f"{fmt_int(u['n_positives'])} positives, {fmt_int(u['n_gene_clusters'])} gene clusters"])
    FLAGS = "not used as a file; flags and publication lists carried inside EPSD records"
    rows.append(["PhosphoSitePlus", "not used as a file", "not used as a file", "not used as a file", FLAGS,
                 "separate PhosphoSitePlus analyses not run: " + tx(psp_note.replace("not run: ", ""), "analysis/results/manifest_status.csv", "G9 PSP rows, result_status")])
    rows.append(["Phospho.ELM", "not used as a file", "not used as a file", "not used as a file", FLAGS, "not used as a file"])
    rows.append(["dbPTM", "not used as a file", "not used as a file", "not used as a file", FLAGS, "not used as a file"])
    rows.append(["PeptideAtlas (human phospho build)", jtx(ds["peptideatlas_phospho_606"]["version"], "analysis/results/downloads_manifest.json", "peptideatlas_phospho_606.version"), jtx(ds["peptideatlas_phospho_606"]["retrieved"], "analysis/results/downloads_manifest.json", "peptideatlas_phospho_606.retrieved"), "not used",
                 "not used: " + tx(pa_note, "analysis/results/manifest_status.csv", "G9 PeptideAtlas rows, result_status"), "not used"])
    df = pd.DataFrame(rows, columns=["Source", "Version or build", "Download date", "Licence", "Records used", "Coverage"])
    write(df, "Table1_data_sources", notes=["Licence texts are quoted from the page of each resource on the date given in the table. The three database rows (PhosphoSitePlus, Phospho.ELM, dbPTM) are the evidence sources of the agreement score: their flags and publication lists are carried inside the EPSD 2.0 records."])
    rec.finish()
    return df


# ------------------------------------------------------------------ Table 2 and S2
def results_rows(pop):
    am, g21, b4, g19 = load("amend_results.json"), load("g21_results.json"), load("b4_twoway.json"), load("g19_unused.json") if False else None
    rows = []
    spec = [("S2", "Conservation", 4, "C2.1"), ("S3", "Disorder", 4, "C2.2"), ("S4", "Sequence motif", 4, "C2.3"), ("S1", "Domain overlap (population-dependent)", 3, "C2.4"), ("S5", "Lookup motif (descriptive only)", 4, "C2.5")]
    two = {"S2": "G1-0005", "S3": "G1-0008", "S4": "G1-0011", "S1": "G1-0002", "S5": "G1-0014"}
    for s, lab, nd, cid in spec:
        if s == "S1" and pop == "P-human":
            continue
        if s == "S1" and pop == "P-human":
            d0 = pick(am, "G19", signal="S1", population="P-human", design="D0", tier_rule="set-B tiers (primary)")[1] if False else None
            keys = [(k, v) for k, v in am.items() if k.startswith("G19") and v["signal"] == "S1" and v["population"] == "P-human" and v["tier_rule"].startswith("set-B")]
            r0 = next(v for k, v in keys if v["design"] == "D0")
            r1 = next(v for k, v in keys if v["design"] == "D1")
        else:
            ppop = "P-all (human)" if s == "S1" else pop
            r0 = pick(am, "G17", signal=s, population=ppop, design="D0")[1]
            r1 = pick(am, "G17", signal=s, population=ppop, design="D1")[1]
        if s == "S1":
            m1 = "not run"
        else:
            g = pick(g21, signal=s, population=f"{pop} (S/T/Y)", model="M1")[1]
            m1 = ec(g["companion_linear"], nd) if "companion_linear" in g else ec(g, nd)
        if pop == "P-all":
            tw = tx("yes" if b4[two[s]]["survives_two_way"] else "no", "analysis/results/b4_twoway.json", f"{two[s]}.survives_two_way")
        else:
            tw = "not computed"
        gr = grade(cid) if pop == "P-all" or s != "S1" else grade(cid)
        rows.append([lab, ec(r0, nd, r0["crude_ci"] if "crude_ci" in r0 and r0.get("design") == "D1" else r0["ci"]) if False else ec(r0, nd), ec(r1, nd), m1, tw, f"{fmt_int(r1['n_sites_gold'])} / {fmt_int(r1['n_sites_bronze'])}; {fmt_int(r1['n_gene_clusters'])}", gr])
    return rows


def table2_notes(nboot, pop):
    n = ["Each estimate is GOLD minus BRONZE with set-B tiers: a difference in proportion (conservation, sequence motif, lookup motif) or in mean score (disorder), and an odds ratio for domain overlap (no difference is 0, and 1 for the odds ratio).",
         f"Intervals are 95% gene-cluster bootstrap intervals from {nboot} resamples. Crude: no adjustment, the same records as the primary column. Primary: coarsened exact matching (CEM) on deciles of the V50 study-effort index. Within gene, M1: GOLD against BRONZE within the same gene, all records of the population; not run for domain overlap.",
         "Two-way interval excludes the null: whether the interval that allows for dependence among proteins and among studies (protein x study) excludes the null value; not computed for the human-only population. GOLD / BRONZE sites; gene clusters: records and gene clusters in the primary estimate.",
         "Grade: strength of the claim under the grading rules of the claims table (strong, moderate, weak, not supported); the database-lookup motif is descriptive only because its indicator comes from the annotation tables that also feed the tiers, and is not graded."]
    if pop == "P-human":
        n.append("Domain overlap is defined for human records only; its estimate is the one in Table 2.")
    return n


def table2():
    rec = Rec("Table2")
    am = load("amend_results.json")
    nb = pick(am, "G17", signal="S2", population="P-all", design="D1")[1]["n_boot"]
    rows = results_rows("P-all")
    cols = ["Signal", "Crude difference (95% CI)", "Primary: CEM on the V50 decile (95% CI)", "Within gene, M1 (95% CI)", "Two-way interval excludes the null", "GOLD / BRONZE sites; gene clusters", "Grade"]
    df = pd.DataFrame(rows, columns=cols)
    write(df, "Table2_main_results", notes=table2_notes(fmt_int(nb), "P-all"))
    rec.finish()
    rec = Rec("TableS2")
    rows = results_rows("P-human")
    nb = pick(am, "G17", signal="S2", population="P-human", design="D1")[1]["n_boot"]
    write(pd.DataFrame(rows, columns=cols), "TableS2_main_results_P-human", notes=table2_notes(fmt_int(nb), "P-human"))
    rec.finish()


# ------------------------------------------------------------------ Table 3
def table3():
    rec = Rec("Table3")
    import provenance as P
    g11 = load("g11_predictors.json")
    tv = P.jload(CUR / "earlier_results" / "tool_versions.json")
    rows = []
    for tool, key, lic in (("DeepPhos", "DeepPhos|planned: S/T sites only", "deepphos"), ("MusiteDeep", "MusiteDeep|planned: S/T sites only", "musitedeep"), ("PhosF3C", "PhosF3C|planned: Y sites only", "phosf3c")):
        u = g11[key]["unrestricted"]
        ver = jtx(tv[tool]["scored_checkout"], "earlier_results/tool_versions.json", f"{tool}.scored_checkout") if tool != "PhosF3C" else jtx(tv[tool]["repository"], "earlier_results/tool_versions.json", "PhosF3C.repository") + "; " + jtx(tv[tool]["scored_weights"], "earlier_results/tool_versions.json", "PhosF3C.scored_weights")
        rows.append([tool, ver, jtx(tv[tool]["training_data_sources"], "earlier_results/tool_versions.json", f"{tool}.training_data_sources"), licence(lic), "serine/threonine" if tool != "PhosF3C" else "tyrosine", fmt_int(u["n_sites"]), fmt_int(u["n_positives"]), f"{tv[tool]['year_used_for_restriction']:d}"])
    sens = []
    for tool, key in (("DeepPhos", "DeepPhos|planned: S/T sites only"), ("MusiteDeep", "MusiteDeep|planned: S/T sites only"), ("PhosF3C", "PhosF3C|planned: Y sites only")):
        years = [k.split("=")[1] for k in g11[key] if k.startswith("restriction_year=")]
        main = f"{tv[tool]['year_used_for_restriction']:d}"
        other = [y for y in years if y != main]
        sens.append(tx(", ".join(other), "analysis/results/g11_predictors.json", f"{key}: keys restriction_year=...") if other else "none")
    for r, sy in zip(rows, sens):
        r.append(sy)
    df = pd.DataFrame(rows, columns=["Tool", "Version scored", "Training data as documented", "Licence", "Residues scored", "Sites scored", "Positives", "Year used for restriction", "Sensitivity year"])
    note = ("LMPhosSite: not scored: " + tx("timebox exceeded before a test prediction was attempted", "audit_reports/phaseH1_step3_predictor_audit.md", "section 'LMPhosSite - NOT COMPLETED, timebox exceeded'") +
            " (there is no predictor result for it).")
    write(df, "Table3_predictors", notes=[
        note,
        "Year used for restriction: the publication year of the version scored; the secondary analysis keeps positives first reported after that year. Sensitivity year: an earlier restriction year that was also run (none for PhosF3C). No tool states when its training data collection ended.",
        "Sites scored and positives refer to the evaluation sample of human sites (serine and threonine for DeepPhos and MusiteDeep, tyrosine for PhosF3C)."])
    rec.finish()


# ------------------------------------------------------------------ supplementary tables
BLOCK_LABEL = {
    "G1": "primary comparison (all species)", "G2": "human-only comparison", "G3": "domain overlap among proteins with a domain-annotated site (diagnostic)", "G4": "study-effort index variants",
    "G5": "protein abundance (PaxDb) as a covariate", "G6": "residue-specific populations (S, T, Y)", "G7": "leave one study cluster out", "G8": "tier thresholds",
    "G9": "established confidence measures", "G10": "within protein (PhosphoSitePlus sites)", "G11": "predictors: AUROC gap by study effort", "G12": "depth-only tier",
    "G13": "study-cluster threshold (500 and 1,000 sites)", "G14": "index without the exclusion and set-A tiers, with protein-cluster intervals", "G15": "set-A tier (human)",
    "G16": "Swiss-Prot reviewed canonical accessions", "G17": "gene-cluster bootstrap of the comparisons in G1, G2, G3 and G5", "G18": "1:1 caliper matcher (D2b)",
    "G19": "tier-definition sensitivity", "G20": "diagnostics: balance and index-outcome association", "G21": "within-gene comparison", "G22": "serine/threonine-only sequence motif",
    "G23": "non-primary PMIDs removed from publication depth"}
TASK_LABEL = {"B,C": "comparison of GOLD and BRONZE with a study-effort index (tasks B and C)", "D": "index variants", "E": "abundance", "F": "residue strata, leave-one-cluster-out, tier thresholds",
              "H": "established measures (PhosphoSitePlus, UniProt, PeptideAtlas)", "I": "within protein (PhosphoSitePlus)", "B": "reconciliation and inference (task B)", "J": "predictors",
              "amend": "not pre-specified: added to the grid after the first results"}
SIGNAL_LABEL = {"S1": "domain overlap", "S2": "conservation", "S3": "disorder", "S4": "kinase motif, sequence-based", "S5": "kinase motif, database lookup", "AUROC gap": "AUROC in deciles 8 to 10 minus deciles 1 to 3 of the index"}
DESIGN_LABEL = {"D0": "crude: no adjustment, the same records as the adjusted design", "D1": "coarsened exact matching (CEM) on deciles of the index (primary)",
                "D2": "propensity matching on the index (sensitivity)", "D2b": "1:1 caliper matcher without replacement (caliper 0.2 SD of the logit)", "D3": "CEM on an index variant",
                "D4": "CEM on index decile x abundance decile", "D5": "CEM on abundance decile alone", "D6 (within-gene)": "within-gene model (one intercept per gene)"}
INDEX_LABEL = {"V50": "primary index: PMIDs linked to more than 50 genes excluded", "V0 (original)": "index without the exclusion", "V0": "index without the exclusion", "V20": "PMIDs linked to more than 20 genes excluded",
               "V100": "PMIDs linked to more than 100 genes excluded", "VT2025": "time-normalised index, reference year 2025", "VT2026": "time-normalised index, reference year 2026"}
# wording of the manifest that the banned-word list does not allow in a released table (the contrast codes and block labels are exempt in scan_banned.py; these are not)
RELABEL = {"as in Rev 5": "index without the exclusion; set-A tiers for domain overlap", "sensitivity: corrected tiers": "sensitivity: tiers with non-primary PMIDs removed",
           "V0 (original)": "V0 (no exclusion)", "within-protein fixed effects": "within-protein model (one intercept per protein)",
           "gene fixed effects": "one intercept per gene", "within-gene linear fixed effect": "within-gene linear model (one intercept per gene)",
           "added by the plan amendment of 2026-10-01 (after the interim report)": "not pre-specified: added to the grid after the first results",
           "added by the plan amendment of 2026-10-04 (after the results report was reviewed)": "not pre-specified: added to the grid after the full results were reviewed",
           "added by the plan amendment of 2026-10-04 (sixth batch; scope freeze lifted for G23 only)": "not pre-specified: added to the grid after the full results were reviewed",
           "Task H, PhosphoSitePlus measures": "PhosphoSitePlus measures", "Task I, within-gene": "within-gene", "see RESULTS_REPORT.md section 9 and ": "see "}


def sources_s1():
    """Per contrast: unit of the interval that manifest_status.csv reports, from the result record the interval was read from (same order as build_results.py)."""
    import glob

    def J(p):
        f = R6 / p
        return json.load(open(f, encoding="utf-8")) if f.exists() else {}

    cem, g5, g14, am, g23_, g22, g21_, g9 = (J(p) for p in ("cem_results.json", "g5_abundance.json", "g14_cem.json", "amend_results.json", "g23_results.json", "g22_results.json", "g21_results.json", "g9_measures.json"))
    prop, d2b, fr = {}, {}, {}
    for f_ in glob.glob(str(R6 / "prop" / "*.json")):
        r_ = json.load(open(f_, encoding="utf-8"))
        prop[r_["contrast_id"]] = r_
    for f_ in glob.glob(str(R6 / "d2b" / "*.json")):
        r_ = json.load(open(f_, encoding="utf-8"))
        d2b[r_["contrast_id"]] = r_
    for f_ in glob.glob(str(R6 / "f_runs" / "*.json")):
        fr.update(json.load(open(f_, encoding="utf-8")))

    def norm(u):
        u = (u or "").lower()
        return "gene cluster" if "gene" in u else ("protein cluster" if "protein" in u else "")

    def unit(cid, blk):
        if cid in cem:
            return norm(cem[cid].get("resampling_unit"))
        if cid in g5:
            return norm(g5[cid].get("resampling_unit"))
        if cid in g14:
            return "protein cluster"
        if cid in prop:
            return norm(prop[cid].get("resampling_unit"))
        if cid in d2b:
            return norm(d2b[cid].get("resampling_unit"))
        if cid in am:
            return norm(am[cid].get("resampling_unit"))
        if cid in g23_:
            return norm(g23_[cid].get("resampling_unit"))
        if cid in g22:
            return norm(g22[cid].get("resampling_unit"))
        if cid in g21_:
            return "gene cluster"
        if cid in fr:
            return norm(fr[cid].get("gene", {}).get("resampling_unit", "gene cluster"))
        if cid in g9:
            return "gene cluster" if g9[cid].get("gene") else ""
        if blk == "G11":
            return "gene cluster"
        return ""
    return unit


def table_s1():
    ms = pd.read_csv(R6 / "manifest_status.csv", encoding="utf-8").fillna("")
    unit = sources_s1()
    ms["interval_unit"] = [unit(c, b) if (lo != "" and hi != "") else "" for c, b, lo, hi in zip(ms["contrast_id"], ms["block"], ms["ci_low"], ms["ci_high"])]
    ms["second_interval_unit"], ms["second_ci_low"], ms["second_ci_high"] = "", "", ""
    key = {(r.block, r.signal, r.design, r.population): i for i, r in ms.iterrows()}
    paired = 0
    for i, r in ms.iterrows():
        m = re.match(r"gene-cluster bootstrap; source block (G\d+)", str(r["variation"]))
        if not m:
            continue
        j = key.get((m.group(1), r["signal"], r["design"], r["population"]))
        if j is None:
            cand = ms[(ms["block"] == m.group(1)) & (ms["signal"] == r["signal"]) & (ms["design"] == r["design"])]
            cand = cand[cand["tier_rule"] == r["tier_rule"]] if len(cand) > 1 else cand
            j = cand.index[0] if len(cand) == 1 else None
        if j is None or ms.at[j, "ci_low"] == "" or r["ci_low"] == "":
            continue
        for a, b in ((i, j), (j, i)):
            ms.at[a, "second_interval_unit"], ms.at[a, "second_ci_low"], ms.at[a, "second_ci_high"] = ms.at[b, "interval_unit"], ms.at[b, "ci_low"], ms.at[b, "ci_high"]
        paired += 1
    for col in ("estimate", "ci_low", "ci_high", "second_ci_low", "second_ci_high", "attenuation"):
        ms[col] = [sig4(v) if v != "" else "" for v in ms[col]]
    for col in ms.columns:
        if ms[col].dtype == object:
            for a, b in RELABEL.items():
                ms[col] = ms[col].astype(str).str.replace(a, b, regex=False)
    cols = ["contrast_id", "block", "task", "signal", "signal_label", "population", "design", "index", "tier_rule", "variation", "decision_bearing", "status", "note", "result_status",
            "estimate", "ci_low", "ci_high", "interval_unit", "second_interval_unit", "second_ci_low", "second_ci_high", "attenuation", "result_note"]
    ms = ms[cols]
    # key
    man = pd.read_csv(R6 / "manifest_status.csv", encoding="utf-8").fillna("")
    rows = [["contrast_id", "G<block>-<row>", "block number and row number in the planned grid of contrasts", ""],
            ["estimate, ci_low, ci_high", "4 significant digits", "estimate of the contrast and its 95% interval; interval_unit says which resampling unit the interval comes from (gene cluster is primary)", ""],
            ["second_interval_unit, second_ci_low, second_ci_high", "", "the interval of the same estimate from the other resampling unit, where the comparison was run with both (protein cluster and gene cluster); empty otherwise", ""],
            ["attenuation", "share", "(crude - adjusted) / crude on the same records, for odds ratios on the log scale; empty where it was not computed", ""]]
    for b in sorted(man["block"].unique(), key=lambda x: int(x[1:])):
        rows.append(["block", b, BLOCK_LABEL.get(b, ""), f"{int((man['block'] == b).sum())} contrasts"])
    for t, lab in TASK_LABEL.items():
        rows.append(["task", t, lab, f"{int((man['task'] == t).sum())} contrasts"])
    for s_, lab in SIGNAL_LABEL.items():
        rows.append(["signal", s_, lab, f"{int((man['signal'] == s_).sum())} contrasts"])
    for d_, n in man["design"].value_counts().items():
        rows.append(["design", d_, DESIGN_LABEL.get(d_, "named as in the grid of contrasts"), f"{int(n)} contrasts"])
    for x_, n in man["index"].value_counts().items():
        rows.append(["index", x_, INDEX_LABEL.get(x_, "as named in the grid of contrasts"), f"{int(n)} contrasts"])
    rows.append(["tier_rule", "study-cluster rule (884)", "tiers from agreement score and study-cluster depth (a PMID linked to at least 884 sites is a cluster candidate); the other values name the variant tested", ""])
    rows.append(["decision_bearing", "yes / no", "whether the contrast enters the decision on a claim (no: diagnostic or sensitivity)", ""])
    rows.append(["status", "planned / conditional", "planned: in the plan; conditional: run only if the data exist", ""])
    for rs, n in man["result_status"].value_counts().items():
        rows.append(["result_status", rs, "outcome of the contrast", f"{int(n)} contrasts"])
    key_df = pd.DataFrame(rows, columns=["Column or code", "Value", "Meaning", "Count"])
    for col in key_df.columns:
        for x, y in RELABEL.items():
            key_df[col] = key_df[col].astype(str).str.replace(x, y, regex=False)
    TABDIR.mkdir(parents=True, exist_ok=True)
    ms.to_csv(TABDIR / "TableS1_contrast_grid.csv", index=False, encoding="utf-8")
    key_df.to_csv(TABDIR / "TableS1_key.csv", index=False, encoding="utf-8")
    (TABDIR / "TableS1_key.md").write_text(key_df.to_markdown(index=False), encoding="utf-8")
    (TABDIR / "TableS1_contrast_grid.md").write_text("Key to the codes: TableS1_key.md" + NL + NL + ms.to_markdown(index=False), encoding="utf-8")
    print("S1:", len(ms), "rows;", paired, "pairs with two interval types")
    return ms, key_df


def tables_s():
    table_s1()
    # internal: section 7 of the claims file, verbatim (not a supplementary table)
    sect = CLAIMS[CLAIMS.index("## 7."):]
    rows = md_table_rows(sect, "| ID | Rev 5 claim")
    hdr = [c.strip() for c in next(l for l in sect.split(chr(10)) if l.startswith("| ID | Rev 5 claim")).strip("|").split(" | ")]
    rec = Rec("InternalChanges")
    for r in rows:
        for c in r:
            tx(c, "analysis/CLAIMS_VS_EVIDENCE.md", f"section 7, row {r[0]}")
    write(pd.DataFrame(rows, columns=hdr), "earlier_to_present_changes", directory=INTERNAL)
    rec.finish()
    # S3: non-primary PMIDs and kept PMIDs
    rec = Rec("TableS3")
    g23, npm, kept = load("g23_tiers_summary.json"), load("g23_nonprimary_pmids.json"), load("g23_kept_nonevidence_pmids.json")
    rows = []
    for x in g23["nonprimary"]:
        why = npm["reasons_b"].get(x["pmid"]) if x["pmid"] in npm["reasons_b"] else ("review-type publication" if x["criterion_a_review_type"] else "")
        stored_why = why
        why = WORDING.get("analysis/results/g23_nonprimary_pmids.json:reasons_b", {}).get(why, {}).get("table", why)  # reworded cell, listed with its reason in results/table_wording.json
        rows.append([x["pmid"], jtx(x["title"], "analysis/results/g23_tiers_summary.json", f"nonprimary[pmid={x['pmid']}].title"), fmt_int(x["records_linked"]), "removed from depth", tx(why, "analysis/results/g23_nonprimary_pmids.json", f"reasons_b.{x['pmid']} or criterion_a_review_type")])
    for pm in ("30243273", "37031204", "14702039"):
        k = kept[pm]
        rows.append([pm, jtx(k["title"], "analysis/results/g23_kept_nonevidence_pmids.json", f"{pm}.title"), fmt_int(k["records"]), "kept (not phosphorylation evidence, rule unchanged)", f"{k['share_of_records'] * 100:.2f}% of records"])
    write(pd.DataFrame(rows, columns=["PMID", "Title", "Records", "Treatment", "Reason or share"]), "TableS3_non_primary_pmids", notes=[
        "Records: site records whose publication list contains the PMID. Treatment: removed from depth = the PMID is not counted in study-cluster depth or in the tiers (sensitivity analysis); kept = the PMID is not phosphorylation evidence but stays in the count, because the rule is unchanged.",
        f"Reason or share: the reason a PMID was classed as a database, resource or review publication, or the share of all {fmt_int(g23['n_records'])} records for a kept PMID. Titles are quoted from PubMed."])
    rec.finish()
    # S4: UniProt experimental against inferred
    rec = Rec("TableS4")
    g9 = load("g9_measures.json")
    rows = []
    for k, v in sorted(g9.items()):
        if isinstance(v, dict) and v.get("status") == "run" and "UniProt" in v.get("measure", ""):
            nd = 3 if v["signal"] == "S1" else 4
            g = v["gene"]
            est, cii = (g["estimate"], g["ci"]) if "estimate" in g else (g["adjusted"], g["adjusted_ci"])
            dn = {"D0": "crude", "D1": "CEM on the V50 decile", "D4": "index x abundance deciles"}[v["design"]]
            rows.append([{"S1": "Domain overlap", "S2": "Conservation", "S3": "Disorder", "S4": "Sequence motif"}[v["signal"]], dn, f"{fd(est, nd)} ({fd(cii[0], nd)} to {fd(cii[1], nd)})",
                         f"{fmt_int(v['n_high'])} / {fmt_int(v['n_low'])}; {fmt_int(v['n_gene_clusters'])}"])
    write(pd.DataFrame(rows, columns=["Signal", "Design", "Experimental minus inferred (95% gene-cluster CI)", "Experimental / inferred records; gene clusters"]), "TableS4_uniprot_experimental_vs_inferred", notes=[
        "Experimental and inferred refer to the evidence code (ECO) of the UniProt modified-residue annotation of the site. Estimates are experimental minus inferred, as a difference in proportion or mean score, and as an odds ratio for domain overlap; intervals are 95% gene-cluster bootstrap intervals.",
        "A limit that would print as zero at the usual number of decimals is given with more decimals so that its sign shows."])
    rec.finish()
    # S5: sample sizes behind Figure 4A
    rec = Rec("TableS5")
    import sys
    sys.path.insert(0, str(CUR / "scripts" / "analysis"))
    from derive_figure_values import tier_label
    am = load("amend_results.json")
    rows = []
    for k, v in sorted(am.items()):
        if k.startswith("G19") and v["design"] == "D1":
            nd = 3 if v["signal"] == "S1" else 4
            rows.append([v["population"], {"S1": "Domain overlap", "S2": "Conservation", "S3": "Disorder", "S4": "Sequence motif"}[v["signal"]], tier_label(v["tier_rule"]), fmt_int(v["n_sites_gold"]), fmt_int(v["n_sites_bronze"]), fmt_int(v["n_gene_clusters"]),
                         f"{v['crude']:.{nd}f} ({v['crude_ci'][0]:.{nd}f} to {v['crude_ci'][1]:.{nd}f})", ec(v, nd)])
    write(pd.DataFrame(rows, columns=["Population", "Signal", "Tier definition", "GOLD sites", "BRONZE sites", "Gene clusters", "Crude (95% CI)", "CEM on the V50 decile (95% CI)"]), "TableS5_figure4_sample_sizes", notes=[
        "The cells behind Figure 4A: GOLD against BRONZE under each tier definition, with the numbers of sites and gene clusters. Intervals are 95% gene-cluster bootstrap intervals; domain overlap is an odds ratio, the other signals are differences."])
    rec.finish()


def purge_bindings():
    """Remove the entries of earlier tables from results/figure_bindings.json so that every table is bound under its present name."""
    d = json.load(open(BIND, encoding="utf-8")) if BIND.exists() else {}
    for k in [k for k in d if k.startswith(("Table", "Internal"))]:
        del d[k]
    BIND.write_text(json.dumps(d, indent=0), encoding="utf-8")


if __name__ == "__main__":
    purge_bindings()
    table1()
    table2()
    table3()
    tables_s()
    print("tables written")
