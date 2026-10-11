# Sources consulted on the web (Rev 6; for the plan only, not for the manuscript or the reference list)

Retrieved 2026-10-04. Rule 13 stays in force: nothing here goes into the manuscript or the reference list. Quotes are in quotation marks. **How the quotes were obtained:** the fetch tool returns an extraction made by a small model from the page, which was instructed to quote verbatim; I could not check the quotes character by character against the page, so they are "as returned by the fetch tool". **Section A is the exception: its quotes were extracted locally from the PDF of the paper (details there).** Where a page could not be read, that is stated and nothing is filled in from memory. Nothing below states a threshold that was not read on a page.

## A. Evidence-category definitions in Kalyuzhnyy et al., J Proteome Res 2022

**Definitions obtained (2026-10-04, from the University of Liverpool repository copy, approved by the project lead).** Blocked fetches earlier the same day are listed in section E and were not retried.

- **Source.** https://livrepository.liverpool.ac.uk/3165686/ , file "Kalyuzhnyy FLR global phospho acs.jproteome.2c00131.pdf" (5,413,399 bytes, SHA-256 508158ec880afbed83c3afea3c363992dea0134c69ec0bcbd6ff15ead64e81e6), downloaded to `data_private/dl/kalyuzhnyy_liverpool.pdf` (git-ignored; not committed or redistributed). The repository page lists the file as "Published version" while its header line says "Preprint (2021)" and states no licence; the PDF itself carries "(c) 2022 The Authors. Published by American Chemical Society" and no licence text. The PDF is the journal version: `Cite This: J. Proteome Res. 2022, 21, 1510−1524`, and `Received: March 4, 2022 Published: May 9, 2022 Article pubs.acs.org/jpr © 2022 The Authors. Published by American Chemical Society`.
- **How quoted.** The text was extracted from the PDF locally (PyMuPDF), not through a summariser. The blocks below are the extracted character spans, with line breaks and end-of-line hyphenation of the PDF kept as extracted, the typographic ligatures (fi, fl, ff) written out as plain letters, and nothing else changed.

### A1. PeptideAtlas categories (Methods, "Processing and Categorizing Phosphorylation Data in PSP and PA")

```
Phosphorylation data in PeptideAtlas (PA) (2020 build)26
was filtered to only include human Ser/Thr/Tyr sites from
canonical UniProt protein sequences with at least one PSM
observation (1 069 709 sites across 63 616 sequences) (Table
S1). The sites were categorized according to the number of
[page break in the PDF]
PSM observations with a certain phosphorylation probability
score assigned by PTMProphet.21 The counts of observations
with a probability of >0.95 were used as “positive” evidence
for site phosphorylation. The counts at a probability threshold
of ≤0.19 were used as “negative” evidence in favor of a site
being a nonphosphosite. The total number of PSM
observations per site was considered to distinguish sites for
which ≥10% of all associated PSMs had a PTM probability
>0.95, from sites where a small minority (<10%) of associated
PSMs had this probability. Based on this, selected confidence
categories were applied to predict site phosphorylation
likelihood in PA (“High”: ≥5 positive observations which
accounted for ≥10% of total observations across all
probabilities; “Medium”: ≥5 positive observations, which
accounted for <10% of total observations or 2−4 positive
observations; “Low”: 1 positive observation; “Not phosphory-
lated”: 0 positive observations and ≥5 negative observations;
“Other”: site did not fall into any described categories).
```

### A2. PhosphoSitePlus categories (Methods, same section)

```
PhosphoSitePlus (PSP) data (11/03/20 build; Phosphoryla-
tion_site_dataset.gz)29 was filtered to only include human
Ser/Thr/Tyr sites from canonical protein sequences labeled
by UniProt identifiers (231 607 sites across 17 830 sequences)
(Table S2). The sites were ranked based on the number of
times they have been characterized in low/high-throughput
studies. The sum of observations across all studies was used to
predict site phosphorylation likelihood in PSP (“High”: ≥5
observations; “Medium”: 2−4 observations; “Low”: 1 observa-
tion).
```

**What this does and does not say about columns.** The paper says the sites "were ranked based on the number of times they have been characterized in low/high-throughput studies" and that "The sum of observations across all studies" was used. It does **not** name the columns LT_LIT, MS_LIT or MS_CST (the strings do not occur in the PDF), and it does not say which columns were summed. Whether n = LT_LIT + MS_LIT + MS_CST reproduces "the sum of observations across all studies" is therefore not documented by this paper; it depends on what the PhosphoSitePlus download documents (see section C and the project lead's instruction of 2026-10-04 to stop before Task H if the download does not document the columns). The paper used the "11/03/20 build; Phosphorylation_site_dataset.gz".

### A3. Table 1 as extracted

```
Table 1. Categorizing Serine (Ser), Threonine (Thr), and Tyrosine (Tyr) Sites from UniProt’s Reference Human Proteome
into Phosphorylation Likelihood Sets Based on Available Phosphorylation Evidence in PhosphoSitePlus (PSP) and
PeptideAtlas (PA)
Phosphorylation
likelihood set
Phosphorylation evidence per site
Ser
count
Thr
count
Tyr
count
High in PSP
5+ pieces of evidence
32 306
8161
9763
Medium in PSP
2−4 pieces of evidence
34 154
12 197
7228
Low in PSP
1 piece of evidence
66 777
35 173
20 191
High in PA
5+ observations at PTM score >0.95 which is ≥10% of total observations
26 186
4204
1169
Medium in PA
5+ observations at PTM score >0.95 which is <10% of total observations OR 2−4 observations at
PTM score >0.95
20 517
5297
1460
Low in PA
1 observation at PTM score >0.95
12 950
4895
1324
Not phosphorylated
0 observations at PTM score >0.19 AND 5+ observations at PTM score ≤0.19 AND no evidence in
PSP
13 892
8462
2184
Other sites
At least 1 observation in PA but does not fall into any other PA categories AND no evidence in PSP
60 221
35 000
10 009
```

**Discrepancy inside the paper (reported, not resolved).** The Methods define "Not phosphorylated" as "0 positive observations and >=5 negative observations" (positive = probability >0.95, negative = probability <=0.19). Table 1 defines it as "0 observations at PTM score >0.19 AND 5+ observations at PTM score <=0.19 AND no evidence in PSP". The zero-count condition differs (no observation above 0.95 versus no observation above 0.19), and the Table adds the PSP condition (which the text states separately: "a site was only categorized as a nonphosphosite if it had no evidence in PSP"). Not relevant for the high-versus-low contrasts used in Task H, which use only the "High" and "Low" definitions.

### A4. Context sentences

```
The latest PA builds incorporate the use of the PTMProphet
algorithm for phosphosite localization where each potential
phosphosite within an observed PSM is assigned a probability
score between 0 and 1 of being phosphorylated.21 As with
PSP, researchers should also be careful when accepting sites in
PA with only a single piece of identification evidence (i.e., a
single associated PSM) as positively identified phosphosites.
Instead, phosphosites that not only have multiple PSM
observations in PA but also have high phosphorylation
probability scores assigned within the majority of those
PSMs are most likely to be true-positive identifications.
```

```
The evidence for phosphorylation comes from
manually curated reviews of literature primarily describing
tandem MS studies and also low-throughput experiments, or
from in-house MS studies.29 Interestingly, the majority of
phosphosites in PSP only have a single piece of evidence
associated with their identification (i.e., there is only one
study identifying the phosphosite). As mentioned in the PSP
documentation itself, researchers should be cautious when
accepting such sites as true positives.
```

```
To eliminate
potential false assignments when considering evidence in both
PSP and PA, a site was only categorized as a nonphosphosite
if it had no evidence in PSP in addition to having “negative”
phosphorylation evidence in PA (Table 1).
```

### A5. What the authors publish as supplementary files (Supporting Information, "ASSOCIATED CONTENT")

```
Filtered PeptideAtlas build with human STY sites that
have at least 1 associated PSM (Table S1) (ZIP)
Filtered PhosphoSitePlus build with human STY sites
from canonical protein sequences (Table S2) (ZIP)
...
Summary of all STY sites in our analysis, their
conservation data, proximal sites, phosphorylation
likelihood, and structural data (Table S5) (ZIP)
...
Counts of STY sites in phosphorylation likelihood sets
based on evidence in PeptideAtlas before considering
the evidence in PhosphoSitePlus (Table S8) (XLSX)
Cross-referencing sets of sites between PhosphoSitePlus
and PeptideAtlas (Table S9) (XLSX)
STY sites in human proteome with plenty of
phosphorylation evidence in both PhosphoSitePlus
and PeptideAtlas (Table S10) (XLSX)
```

The PDF gives the location as "https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00131" (blocked for me, HTTP 403 on the sister ACS page); file names, sizes and licence of the supplementary files are not stated in the PDF. **Table S1 ("Filtered PeptideAtlas build with human STY sites that have at least 1 associated PSM") and Table S5 are site-level files from the authors**; they are not available to me. The project lead will place supplementary files in `study/analysis/sources/` if obtainable.

### A6. Why the PeptideAtlas categories cannot be applied to build 606 as downloaded

Every PeptideAtlas category above is defined by counts of PSM observations of a site at a PTMProphet probability ("5+ observations at PTM score >0.95", "<=0.19", ">0.19", and the share of total observations). The TSV dump of build 606 has no PTMProphet probability and no per-site PSM table (`results/peptideatlas_606_feasibility.json`); the server's PTM Coverage table uses probability bands (>=0.80, >=0.95, >=0.99) that do not match ">0.95" and "<=0.19" and is aggregated per protein, not per site (section B). So the categories cannot be applied faithfully to build 606 from what we hold.

## B. Site-level PeptideAtlas file for the human phosphoproteome

**Result: no per-site file found.**

1. Build page, https://peptideatlas.org/builds/human/phospho/ (fetched; content below as returned):
   - Statistics: "3,091,839 distinct phosphopeptides" across "11,852 canonical proteins" from "138 datasets" and "240 experiments".
   - Files listed (size as shown on the page): Biosequence Set in FASTA format [662MB] `202512/Homo_sapiens.fasta`; Database tables exported as TSV dump file [523MB] `202512/atlas_build_606.tsv.zip`; Database tables exported as mysql dump file [600MB] `202512/atlas_build_606.mysql.gz`; Peptide CDS and chromosomal coordinates [1GB] `202512/coordinate_mapping.txt`; Peptide CDS coordinates [765MB] `202512/APD_ensembl_hits.tsv`; Peptide sequences in FASTA format [11MB] `202512/APD_Hs_all.fasta`.
   - Description: "a compendium of results from uniformly reprocessed mass spectrometry proteomics datasets", processed with "the Trans-Proteomic Pipeline suite of tools".
   - The returned content mentions no site-level file, no PTMProphet localisation file and no licence or terms of use. Licence terms: **not found**.
   - Directory listing https://peptideatlas.org/builds/human/phospho/202512/: HTTP 403 Forbidden. Fetch stopped.
2. Build details page, https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/buildDetails?atlas_build_id=606 (fetched):
   - "PSM FDR threshold: 0.0000500"; "Probability threshold: n/a".
   - It has a "PTM Coverage" table with columns `nP.80-.95`, `nP.95-.99`, `nP.99-1` and `no-choice`, explained on the page as: "nP.95-.99: Number of observed PSM with PTMProphet probability >= 0.95 and < 0.99"; "nP.99-1: PTMProphet probability >= 0.99"; "no-choice: Number of PSMs covering this site for which there was no choice in the localization of the PTM. Only one residue was an S, T, or Y".
   - The page does not define "total_S_sites", "obs_protein_S_sites" or "covered_S_sites".
   - Its "Download as: TSV" link gives the PTM Coverage table itself. I downloaded it (2,277 bytes, saved as `data_private/dl/pa606_ptm_coverage.tsv`). It is a **summary** of site counts per sequence database (CoreProteome, GENCODE, UP000005640, UniProt TrEMBL, Ensembl, ...) and residue (A, S, T, Y), not a per-site file.
3. Per-protein PTM summary, https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPTMSummary?display_options=proteinView&action=QUERY&redundancy_constraint=on&proteome_component=Description,coreProteome&ptm_type=ASTY:Phospho&atlas_build_id=606&apply_action=QUERY (fetched, first page only): column headers "Biosequence Name | Presence Level | Redundancy Relationship | Redundant With Respect To | PTM Type | N_PTM_Sites nP.8-1 | nP.80-.95 | nP.95-.99 | nP.99-1 | no-choice | Biosequence Gene Name | Biosequence Description"; "Displayed rows 1 - 50 of 13884"; each row is one protein. The PTMProphet columns count PSMs (per the page's legend text as extracted: "Number of observed PSM with PTMProphet probability"), aggregated **per protein**, not per site. The page showed no "Download as" link with a file address; a generic "Download ResultSet in Format:" option exists. I did not download or scrape this table (13,884 proteins; per-protein pages would be needed for any per-site view, which would be scraping a web service, not a published file).

So: the PeptideAtlas server holds PTMProphet-based evidence, but the TSV dump of build 606 (what I inspected) does not, and no per-site file is published on the pages above.

## C. Meaning of LT_LIT, MS_LIT, MS_CST (PhosphoSitePlus)

**Not confirmed from a primary page.**

- phosphosite.org downloads page https://www.phosphosite.org/staticDownloads.action: HTTP 403 Forbidden. Fetch stopped; other phosphosite.org pages not tried.
- PhosphoSitePlus 2014 paper, https://academic.oup.com/nar/article/43/D1/D512/2439467 (fetched): it does **not** define the column names. Sentences returned: "only 4.5% of the sites in PSP are derived from LTP studies"; "Over 95% of the sites are from mass spectrometry (MS) experiments"; "Since 2003, over 340 000 non-redundant MS2 PTM sites have been curated into PSP from over 8000 separate experiments performed at CST and from 202 publications"; and a reference to downloads containing "the number of HTP and LTP records".
- PhosphoSitePlus 2012 paper, https://academic.oup.com/nar/article/40/D1/D261/2903142 (fetched): it does not define the column names either. Sentence returned: "7168 (67%) of these LTP sites are associated with only 1 record, 2745 (26%) with 2-5 records, 615 (6%) with 6-25 records and 173 (2%) with >25 records."
- Several web-search summaries attributed meanings to the three columns, but they were not quotations from a page I read and name no readable source; **they are not used here** and no definition is asserted.
- **File version.** I cannot confirm the version of the file you will place. On arrival I will record the file name, size, modification date and the header lines inside the file (which normally carry the version stamp and the licence line) here, and apply the stop rule of decision D3 (the task stops if the columns do not match the plan's definition). The licence terms of PhosphoSitePlus could not be read (403); local-only use is rule 15 of the project instructions.

## D. Camacho et al., Journal of Proteome Research 2023, "Assessing Multiple Evidence Streams to Decide on Confidence for Identification of Post-Translational Modifications, within and Across Data Sets"

- Located: bioRxiv preprint https://www.biorxiv.org/content/10.1101/2022.12.15.520504v1.full (not fetched, bioRxiv returned 403 for the other paper); journal page https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00823 (fetched): HTTP 403 Forbidden. Fetch stopped. A search-tool summary gave the DOI 10.1021/acs.jproteome.2c00823 and "Volume 22, Issue 6, pages 1828-1842" (not read from a page).
- Only the abstract (Consensus connector, earlier in this cycle) was read: "Here, we describe a simple pipeline aiming to maximize the information extracted from these studies by objectively collapsing from peptide-spectrum match to the peptidoform-site level, as well as combining findings from multiple studies while maintaining track of false localization rates." and "In our case study using eight rice phosphoproteomics data sets, 6368 unique sites were confidently identified using our decoy approach compared to 4687 using traditional thresholding in which false localization rates are unknown."
- **Whether it provides site-level categories or data usable here: not determined.** The abstract describes a pipeline and a rice case study; it says nothing about human data, categories or a downloadable file. The full text, the data availability statement and the supplementary material were not readable.

## E. Summary of blocked fetches (also in the deviation log)

PMC (reCAPTCHA); bioRxiv full text (403); figshare record (403); phosphosite.org downloads page (403); ACS journal page (403); PeptideAtlas build directory listing (403). None was retried or approached by another route.


## F. Licences and provenance for Tables 1 and 3 (fetched 2026-10-07)

Each row: the page or file fetched, the access date, the sentence quoted from it, and the cell text the tables use. A cell never comes from memory. "not stated" means the page or repository was reached and states no licence. `scripts/figures/build_tables.py` reads the table below (key, cell text) and the quoted text is kept here for the check.

| key | URL | accessed | quoted sentence | cell text |
|---|---|---|---|---|
| epsd | http://epsd.biocuckoo.cn/ | 2026-10-07 | "all data sets and annotations are freely accessed for all users"; footer: "Copyright © 2010-2026.The CUCKOO Workgroup. All Rights Reserved" | not stated (checked 2026-10-07); site text: data sets and annotations freely accessed |
| gene2pubmed | https://www.ncbi.nlm.nih.gov/home/about/policies/ | 2026-10-07 | "NCBI itself places no restrictions on the use or distribution of the data contained therein." (the page adds that some submitters of the original data may claim rights in it) | NCBI places no restrictions on use or distribution (checked 2026-10-07) |
| uniprot | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/LICENSE | 2026-10-07 | "We have chosen to apply the Creative Commons Attribution 4.0 International (CC BY 4.0) License (https://creativecommons.org/licenses/by/4.0/) to all copyrightable parts of our databases." | CC BY 4.0 (checked 2026-10-07) |
| paxdb | https://pax-db.org/ | 2026-10-07 | structured data on the home page: "license": "https://creativecommons.org/licenses/by/4.0/" | CC BY 4.0, declared in the page's structured data (checked 2026-10-07) |
| pharos | https://pharos.nih.gov/about | 2026-10-07 | "Data accessed from Pharos and TCRD is publicly available from the primary sources listed above. Please respect their individual licenses regarding proper use and redistribution." | no single licence; the primary sources' licences apply (checked 2026-10-07) |
| deepphos | https://github.com/USTC-HIlab/DeepPhos (contents listing via the GitHub API; README.md and readme.txt read) | 2026-10-07 | no LICENSE file (files: README.md, methods, models, predict.py, readme.txt, requirement.txt, sh.exe.stackdump, test data.csv, train.py); neither readme contains the word "licen" | not stated (checked 2026-10-07) |
| musitedeep | https://raw.githubusercontent.com/duolinwang/MusiteDeep_web/master/LICENSE | 2026-10-07 | "MIT License ... Copyright (c) 2019 duolinwang ... Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction" | MIT (checked 2026-10-07) |
| phosf3c | https://raw.githubusercontent.com/SkywalkerLuke/PhosF3C/main/LICENSE | 2026-10-07 | "MIT License ... Copyright (c) 2024 SkywalkerLuke ... Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction" | MIT (checked 2026-10-07) |

Not fetched: PhosphoSitePlus and PeptideAtlas are not used (Table 1 gives the reason, section A to E above). LMPhosSite was not scored in Rev 6 (no G11 result), so it has no row in Table 3.

**EPSD provenance.** Version: the home page text reads "Now in the 2.0 version" (fetched 2026-10-07). The human file `Homo sapiens.zip` was downloaded from `http://epsd.biocuckoo.cn/Download/Homo%20sapiens.zip` (289 MB, SHA-256 beginning `1a760f6f`); the Phase 0 report recording it is dated 2026-09-13, so the download date is **on or before 2026-09-13** and the exact date was not recorded. **gene2pubmed:** `processed_data/epsd/sei/raw_cache/gene2pubmed.gz.sidecar.json` records `retrieved_date` 2026-09-24 for the file gene2pubmed.gz from the NCBI Gene FTP ( 287,184,886 bytes, SHA-256 beginning `6b8153bd`).
