"""G23: write the list of non-primary PMIDs (criterion b, from titles alone) with the reason for each entry, and the borderline titles that are NOT removed.
Reasons were written by me from the titles after fetching them (results/g23_pubmed_metadata.json) and before any depth, tier or contrast was computed.
Outputs: study/analysis/nonprimary_pmid_list.md and results/g23_nonprimary_pmids.json. Criterion (a), PubMed publication type Review or Systematic Review, is applied mechanically from the metadata."""
import json

from common import REV6, RES

meta = json.load(open(RES / "g23_pubmed_metadata.json"))
R = meta["records"]
T = lambda p: R[p]["title"]  # noqa: E731

# (b) database or resource paper, judged from the title alone. Principle: the title uses the word database, knowledgebase / knowledge base, repository, resource,
# data set(s) / dataset or collection as the subject of the paper, or is a named database followed by an update, year or release marker.
B = {
    "15174125": "title calls it a bioinformatics resource (PhosphoSite)",
    "15489334": "title is a project report on a gene collection (Mammalian Gene Collection), a resource",
    "16381945": "title: an information repository (dbPTM)",
    "18988627": "named database with an update marker (Human Protein Reference Database, 2009 update)",
    "19060867": "title: a database (PhosphoPep)",
    "19366988": "title: a systematic resource (SysPTM)",
    "20428315": "title: a database (PhosphoGRID)",
    "21062810": "title: a database, update 2011 (Phospho.ELM)",
    "21081558": "title: the posttranslational modification database (PHOSIDA 2011)",
    "21082442": "title: a phosphoproteome resource",
    "22135298": "title: a comprehensive resource (PhosphoSitePlus)",
    "23193284": "title: a database (PTMcode)",
    "23193290": "named database with a release marker, title also says informative resource (DbPTM 3.0)",
    "24174276": "title: proteomic data sets as the subject",
    "24705204": "title: an updated systematic resource (SysPTM 2.0)",
    "25514926": "named database with a year marker and update topics (PhosphoSitePlus, 2014: mutations, PTMs and recalibrations)",
    "26425664": "title: a dataset paper (Dataset from the global phosphoproteomic mapping ...)",
    "29145615": "title: an integrated resource (iPTMnet)",
    "30087585": "title: tissue repositories as a resource",
    "30395287": "title: a hub of protein knowledge (UniProt), a knowledgebase",
    "30418626": "named database with a year marker (dbPTM in 2019)",
    "30445427": "named database with an anniversary marker (15 years of PhosphoSitePlus)",
    "30476227": "title: an interaction database, 2019 update (BioGRID)",
    "32508104": "title: a comprehensive resource (Scop3P)",
    "32618424": "title: a phosphoproteomic resource (LymphoAtlas)",
    "35608653": "title: a database (iKiP-DB)",
    "37582339": "title: data and resources for pan-cancer analysis (also PubMed type Review, criterion a)",
    "37619559": "title: a data-driven knowledge base of human cancer",
}
# borderline or odd titles that are NOT removed, with the reason (so that the boundary is visible)
N = {
    "14702039": "sequencing and characterization of cDNAs; no database or resource word in the title",
    "18669648": "an atlas as the result of a primary study; no database or resource word",
    "21183079": "an atlas as the result of a primary study; no database or resource word",
    "26074081": "a phosphopeptide atlas as the result of a primary multi-protease study; no database or resource word",
    "36631611": "an atlas of substrate specificities as the result of a primary study; no database or resource word",
    "21630457": "a bioinformatics platform (software), not named as a database or resource",
    "28150246": "iPTMnet title says integrative bioinformatics, not database or resource (its 2018 title, 29145615, does and is removed)",
    "23312004": "a primary identification study that mentions databases in its title; the subject is not a database",
    "22673903": "quantitative maps of phosphorylation sites; a primary study",
    "37582357": "pan-cancer primary analysis paper from the same consortium release",
    "37582358": "pan-cancer primary analysis paper from the same consortium release",
    "37582362": "pan-cancer primary analysis paper from the same consortium release",
    "37582371": "pan-cancer primary analysis paper from the same consortium release",
    "38359819": "pan-cancer primary analysis paper from the same consortium release",
    "37031204": "a correction notice to a methods paper; not covered by criterion (a) or (b)",
    "30243273": "title on a quantum-chemistry method, unrelated to phosphorylation; not covered by criterion (a) or (b)",
}
rev = sorted(p for p, v in R.items() if any(t in ("Review", "Systematic Review") for t in (v.get("pubtype") or [])))
assert all(p in R for p in B) and all(p in R for p in N)
nonprimary = sorted(set(rev) | set(B))
lines = ["# G23: non-primary publications (written before any depth, tier or contrast was computed)\n",
         f"Candidates: {len(meta['candidates'])} PMIDs (the study-cluster candidates linked to at least 884 records, which include PMID 25514926); PubMed metadata retrieved {meta['retrieved']} (`results/g23_pubmed_metadata.json`). "
         "Rule (ANALYSIS_PLAN.md, sixth batch): a PMID is non-primary if (a) its PubMed publication types include Review or Systematic Review, or (b) it is a database or resource paper, judged from the title alone. "
         "The reasons below were written from the titles after fetching them. Principle for (b): the title uses database, knowledgebase or knowledge base, repository, resource, data set(s) or dataset, or collection as the subject of the paper, or is a named database followed by an update, year, release or anniversary marker.\n",
         f"## (a) PubMed publication type Review or Systematic Review: {len(rev)} PMID(s)\n"]
for p in rev:
    lines.append(f"- {p}: {T(p)} (types: {', '.join(R[p]['pubtype'])})")
lines.append(f"\n## (b) Database or resource papers, from the title alone: {len(B)} PMIDs\n")
lines.append("| PMID | Title | Reason |\n|---|---|---|")
for p in sorted(B):
    lines.append(f"| {p} | {T(p)} | {B[p]} |")
lines.append(f"\n## Not removed: borderline or odd titles ({len(N)} PMIDs)\n")
lines.append("| PMID | Title | Reason for not removing |\n|---|---|---|")
for p in sorted(N):
    lines.append(f"| {p} | {T(p)} | {N[p]} |")
lines.append(f"\nTotal non-primary PMIDs removed in G23: {len(nonprimary)} (union of (a) and (b)). Every other candidate is treated as primary.\n")
open(REV6 / "nonprimary_pmid_list.md", "w", encoding="utf-8").write("\n".join(lines))
json.dump({"nonprimary": nonprimary, "criterion_a_review_type": rev, "criterion_b_title": sorted(B), "reasons_b": B, "not_removed_borderline": N}, open(RES / "g23_nonprimary_pmids.json", "w"), indent=1)
print("non-primary:", len(nonprimary), "| (a):", rev, "| (b):", len(B), "| not removed listed:", len(N))
