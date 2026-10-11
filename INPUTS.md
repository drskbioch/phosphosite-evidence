# Inputs

Every input that the scripts in `study/scripts/analysis/` read, with its source, version, retrieval date, SHA-256 and licence. Licence text is quoted from `study/analysis/sources/SOURCES.md`, section F (pages fetched 2026-10-07). Where section F has no row for an input, the entry says so; no licence is stated from memory. SHA-256 values come from `study/analysis/results/downloads_manifest.json` unless the entry says otherwise.

**Reports, figures, tables and the number checks need none of these inputs.** They run from the stored results in `study/analysis/results/` (see `REPRODUCE.md`, steps R, F and C). The inputs below are needed only to recompute the results.

Default locations: downloads go to `study/analysis/data_private/dl/` (git-ignored; project rule 15), except `gene2pubmed.gz`, which the first script reads from `study/processed_data/epsd/sei/raw_cache/`.

**Not redistributed.** EPSD-derived site tables (section 1); PhosphoSitePlus, Phospho.ELM and dbPTM (section 2); predictor code and weights (section 5). Nothing in sections 1 and 2 is part of the archive, which is deposited from a GitHub release (`RELEASE_SCOPE.md`).

## 1. EPSD and the tables derived from it

### 1.1 EPSD 2.0 site files (human, mouse, rat)
- **Source:** http://epsd.biocuckoo.cn/ , download page `Download.php`. Files: `Homo sapiens.zip`, `Mus musculus.zip`, `Rattus norvegicus.zip` at `http://epsd.biocuckoo.cn/Download/<name with %20 for the space>.zip` (human 303,083,518 bytes; mouse 84,745,295; rat 23,821,878).
- **Version:** the home page text reads "Now in the 2.0 version" (fetched 2026-10-07).
- **Retrieved:** the files used for the stored results: human, on or before 2026-09-13 (the exact date was not recorded; the report that records the download is dated 2026-09-13); mouse and rat: not recorded at download. The mouse and rat URLs above were not recorded at download either; they are the URLs on the download page on 2026-10-08.
- **SHA-256.** The original download files are no longer on the machine. The values below are those of a fresh download on 2026-10-08. The server reports Last-Modified 2024-10-08 for all three files, which is earlier than their use. The human hash begins `1a760f6f` and ends `d3c8e0`, as the 2026-09-13 report recorded for the file that was used, so it is the same file. The mouse and rat hashes cannot be compared with a recorded value (none was kept), but the row counts of their text files (534,246 and 162,097) equal the site counts of the tables built from them, as do the human rows (1,071,725).
  - `Homo sapiens.zip` `1a760f6f85f6c4cdd125e88a67ae8636c820b10d60a4053131ddcde419d3c8e0`
  - `Mus musculus.zip` `489bc622dde34bdacfb536ccf11b724d6b38a1a9d878f40bc2078acee378632e`
  - `Rattus norvegicus.zip` `296ee68186a23e96cb02a79d24f49298f7def99e290055d7a075df00f66675ca`
- **Licence (section F, row `epsd`):** "not stated (checked 2026-10-07); site text: data sets and annotations freely accessed". Quoted from the site: "all data sets and annotations are freely accessed for all users"; footer: "Copyright © 2010-2026.The CUCKOO Workgroup. All Rights Reserved".
- **Redistribution:** not redistributed: licence not stated by the provider (permission requested). The fresh copies are held in `study/analysis/data_private/dl/` (git-ignored).
- **You download it:** yes, from EPSD, if you want to run the earlier-stage extraction. That extraction is not a Rev 6 step, and the repository states it is not reproducible from the repository alone (`README.md`, "Note on reproducibility").

### 1.2 EPSD-derived tables used by the Rev 6 scripts
Copied under their release names and column names by `study/scripts/internal/build_release_data.py` from tables that the earlier-stage scripts (`study/ingest/`, `study/scripts/epsd_phase0/`, `study/scripts/earlier_predictor_audit/`) built from the EPSD files. They hold one row per site, protein or PMID. Location: `study/processed_data/epsd/release/`. SHA-256 below is of the local copy, computed on 2026-10-08 (the download manifest does not cover files built from EPSD).

| File | Bytes | SHA-256 (local copy) |
|---|---|---|
| `mammalian_sites.parquet` | 40,185,391 | `37db691c82da28ec476512cddd53f7b4038b598ee9cc2fd5cacc495514554976` |
| `cluster_split.parquet` | 42,365,333 | `d8fd63904472c52e99c9a0d51c48c34a564295905400f493b0c92b371db6d4a3` |
| `kinase_motif_sequence.parquet` | 38,496,289 | `77c77df4882211fdeb1e8a43c6689dd4d6b50463163e8c7237833677429d21ef` |
| `human_sites_source_set_a.parquet` | 28,388,522 | `64daed74999dec025355f79a6b636124f1bfbc80138267d154aaa2eb1eb1ffc1` |
| `temporal_split.parquet` | 17,799,768 | `f033690c4d0e37c4accda008af65ba7d6560cd0980da04815d4c7e4259c6ffca` |
| `site_first_year.parquet` | 15,911,106 | `0979a2baa268c37a44ae2297d5dd0bb541061408770736f9c18ff1d4485326f6` |
| `sei_protein.parquet` | 2,170,157 | `a07fb1e0d4703b1f3f7607747836b10aa9a1e5287a071d1b666cb3b99af04cf0` |
| `pmid_year.parquet` | 418,534 | `2c0575c58205d67b47531c4f8365c326f88933fb7ba15e1ff25206ebd29fe7c9` |

- **Licence and redistribution:** not redistributed: licence not stated by the provider (permission requested). Same row as 1.1.
- **You download it:** not downloadable. Obtain EPSD (1.1) and run the earlier-stage pipeline, or ask the authors of this repository for the tables once EPSD permission is settled.
- **Also EPSD-derived and not redistributed:** `study/processed_data/epsd/predictor_audit/eval_sample_with_sequence.parquet` (10,000-site evaluation sample with sequences; 12,867,923 bytes; SHA-256 `00b2b14cb2cd7b113d2c20cf45a4ffc5c70dd378ca027e8a6c360d746a9a15c6`, local copy, computed 2026-10-08) and the predictor score files in the same folder (section 5).
- **Caches written by the Rev 6 scripts** (`study/analysis/data_private/*.parquet`, 1.1 GB with the downloads): derived from the tables above and the downloads below; they inherit the restrictions of their inputs and are not redistributed.

## 2. PhosphoSitePlus, Phospho.ELM, dbPTM

- **How they enter:** only as membership flags and publication lists inside the EPSD-derived tables (columns `source_PSP`, `source_PhosphoELM`, `source_dbPTM`; EPSD annotation exports `PhosphoSitePlus.txt` and `Phospho.ELM.txt`). No Rev 6 script reads a file downloaded from these databases.
- **PhosphoSitePlus bulk file (`Phosphorylation_site_dataset.gz`):** not obtained ("PhosphoSitePlus v7 site provides no bulk download; access requested from the provider"). Tasks H and I and G20-0488 were not run (`study/analysis/ANALYSIS_PLAN.md`, section 15).
- **Licence:** section F has no row for these three databases.
- **Redistribution:** not redistributed (project rule 15: used locally only).
- **You download it:** not needed for any stored result.

## 3. Public downloads (you download these to recompute the results)

### 3.1 NCBI `gene2pubmed.gz`
- **Source:** https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz (NCBI Gene FTP). The file changes over time; identify the file used by its hash.
- **Version:** the file as retrieved; 287,184,886 bytes.
- **Retrieved:** 2026-09-24 (sidecar `study/processed_data/epsd/sei/raw_cache/gene2pubmed.gz.sidecar.json`).
- **SHA-256:** `6b8153bdaa908849c754b59154f2ff3c2bb30b7de56706c057c7a2a7d5b28dee` (sidecar; not in `downloads_manifest.json`; matches the local file, checked 2026-10-08).
- **Licence (section F, row `gene2pubmed`):** "NCBI itself places no restrictions on the use or distribution of the data contained therein." (https://www.ncbi.nlm.nih.gov/home/about/policies/; the page adds that some submitters of the original data may claim rights in it.)
- **You download it:** yes. Read by `d1_gene2pubmed.py`.

### 3.2 UniProt ID mapping, human, mouse, rat
- **Source:** `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{HUMAN_9606,MOUSE_10090,RAT_10116}_idmapping.dat.gz`. The path serves the current release; a later download may return a newer release, so check the hash.
- **Version:** UniProt Knowledgebase release 2026_03 (02-Sep-2026); Last-Modified 03 Sep 2026.
- **Retrieved:** 2026-10-01 (manifest). The same three files (same hashes) were recorded on 2026-09-24 in the sidecars next to `gene2pubmed.gz`.
- **SHA-256:**
  - HUMAN_9606: `5fb356e81b07a938ac4d8e8a4b5e0faa89001bcff2aef59422ff30eaa8e72996` (38,464,542 bytes)
  - MOUSE_10090: `c9cb0e15fb417b550ecd7151e3eed9f42eecb1e334cb35d174fca90fd292a2ea` (15,536,610 bytes)
  - RAT_10116: `45ba24ca5c95c19afaa3a6790fc317795a84fc5c3ccbf951b6104abee19d0a7f` (9,415,486 bytes)
- **Licence (section F, row `uniprot`):** "We have chosen to apply the Creative Commons Attribution 4.0 International (CC BY 4.0) License (https://creativecommons.org/licenses/by/4.0/) to all copyrightable parts of our databases." Cell text: CC BY 4.0.
- **You download it:** yes (`python d_downloads.py idmap`). Read by `e_abundance_map.py`.

### 3.3 UniProt reviewed-entry lists and sequences, human, mouse, rat
- **Files:** `rev_{9606,10090,10116}.tsv` (columns Entry, Reviewed, Entry Name, Gene Names (primary); 20,431, 17,283 and 8,233 rows) and `seq_{9606,10090,10116}.tsv` (columns Entry, Sequence, Modified residue; same rows).
- **Source:** the UniProt REST service, endpoint `https://rest.uniprot.org/uniprotkb/stream`, format `tsv`, query `(reviewed:true) AND (organism_id:<taxon>)`; fields `accession,reviewed,id,gene_primary` for the lists and `accession,sequence,ft_mod_res` for the sequences. **Script:** `study/scripts/analysis/d_uniprot_reviewed.py` (it also writes `swissprot_reviewed.parquet`, the three lists stacked in the order 9606, 10090, 10116 with a `taxon` column, 45,947 rows).
- **Version:** UniProt release 2026_03 (release date 02-September-2026), read from the response headers.
- **Retrieved:** the files used for the stored results on 2026-10-01 (HTTP `Date` header: `rev_9606` 16:20:20 GMT, `seq_9606` 16:47:26 GMT). The script was written afterwards, with the query reconstructed from the column headers of those files, and run on 2026-10-08.
- **Check of the script.** On 2026-10-08 the service still served release 2026_03, and a fresh download with the script gave files whose SHA-256 values equal the ones below, for all six files. `swissprot_reviewed.parquet` written by the script has the same content as the stored one (the stored file kept the separate index of each list; the script numbers the rows 0 to 45,946; `run_amend.py` reads only the `Entry` column). When UniProt publishes a later release the service will return different files; compare the release number and the hashes.
- **SHA-256** (recorded by the script in `downloads_manifest.json`, keys `uniprot_reviewed_list_<taxon>` and `uniprot_reviewed_sequences_<taxon>`):
  - `rev_9606.tsv` `cdc8c18163a8c027f2c38bc5ea486ee41d63988cdfd3629f86fcd7ab2956fda4`
  - `rev_10090.tsv` `889a3108d4607a0b254e8e4fe18ae7d075dbbc3972550e90b9de0605746f6b3f`
  - `rev_10116.tsv` `9f29b2ee7be0ddb515dc9146f88050f68323d80359d22cecaa6f1e2510ce0910`
  - `seq_9606.tsv` `01c8cb12ab66b991e521e17918fd4058cfbcb777ad9e561dd5e10fb1f824a1a1`
  - `seq_10090.tsv` `92365f51124cefaee54787efc7c222c5dbc47a1ce4f9a84aa155f6cd738fbc65`
  - `seq_10116.tsv` `ad5ee08c64a1cdd6db6b7835a1d5003ff33b5334f0e19044e9a32dca3e1ccb51`
- **Licence:** UniProt row of section F (3.2): CC BY 4.0.
- **You download it:** yes (`python d_uniprot_reviewed.py`). Read by `g21_data.py` and `h_measures.py` (`seq_*.tsv`) and, through `swissprot_reviewed.parquet`, by `run_amend.py` (`rev_*.tsv`).

### 3.4 PaxDb 6.0 integrated whole-organism abundance, human, mouse, rat
- **Source:** `https://pax-db.org/downloads/6.0/datasets/{9606,10090,10116}/{taxid}-WHOLE_ORGANISM-integrated.txt`
- **Version:** PaxDb 6.0 (v6.0 download path); Last-Modified 01 Oct 2025. Dataset ids in the file headers: 1825491670 (human), 3905858159 (mouse), 4094364361 (rat).
- **Retrieved:** 2026-10-01.
- **SHA-256:**
  - 9606: `e019aaf12f0520775674f44ccfdfd16769de3f442de483147556f25205222373` (644,763 bytes, 19,483 data rows)
  - 10090: `5985202b190952c8ed8ad7f10f85749226bec330aa38e8d13c3f20b0ae64be66` (736,156 bytes, 19,774 data rows)
  - 10116: `35ef575fb510be6c6f92d6863b893ebb21bc9a01c78ff9504dc149be08eac9d1` (605,415 bytes, 16,184 data rows)
- **Licence (section F, row `paxdb`):** structured data on the home page: "license": "https://creativecommons.org/licenses/by/4.0/". Cell text: CC BY 4.0, declared in the page's structured data.
- **You download it:** yes (`python d_downloads.py paxdb`). Saved as `px_{taxid}.txt`. Read by `e_abundance_map.py`.

### 3.5 Pharos / TCRD target development levels
- **Source:** https://pharos-api.ncats.io/graphql (GraphQL query `targets{targets(top,skip){uniprot sym tdl fam}}`).
- **Version:** `pharos319`; 20,412 targets (Tbio 12,303; Tdark 5,501; Tchem 1,904; Tclin 704).
- **Retrieved:** 2026-10-01.
- **SHA-256:** none; the response is not a file. The stored `pharos_tdl.parquet` (334,859 bytes) has SHA-256 `6700082294e68ca1708848d3839ef3bfb82971f30316269c6566cf366ee5687f` (local copy, computed 2026-10-08; a parquet file written by another pandas or pyarrow version will have a different hash).
- **Licence (section F, row `pharos`):** "Data accessed from Pharos and TCRD is publicly available from the primary sources listed above. Please respect their individual licenses regarding proper use and redistribution." Cell text: no single licence; the primary sources' licences apply.
- **You download it:** yes (`python d_downloads.py pharos`). Read by `d3_d4_d5_index_checks.py`.

### 3.6 Gene Ontology annotation, human (GOA)
- **Source:** https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz
- **Version:** GAF 2.2, generated 2026-07-28 16:23 by UniProt; go-version http://purl.obolibrary.org/obo/go/releases/2026-07-26/extensions/go-plus.ofn; Last-Modified 28 Jul 2026.
- **Retrieved:** 2026-10-01.
- **SHA-256:** `92504b0208fc37d1b652beaab883298c169b0be342888b91924975a49778cfdb` (11,020,779 bytes). No `.md5` file is published at the expected path; integrity rests on the size equal to Content-Length and a complete gzip read.
- **Licence:** section F has no row for GOA.
- **You download it:** yes (`python d_downloads.py goa`). Read by `d3_d4_d5_index_checks.py`. The path serves the current file, so check the hash.

### 3.7 PubMed record data through NCBI E-utilities
- **Source:** https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi (db=pubmed). Used for (a) publication years of PMIDs not in `pmid_year.parquet` (`d2_index_variants.py`, stored as `pmid_year_gene.parquet`) and (b) title, journal, date and publication type of 690 candidate PMIDs (`g23_fetch.py`, stored in `study/analysis/results/g23_pubmed_metadata.json`).
- **Version:** PubMed as served on the retrieval date.
- **Retrieved:** (a) 2026-10-01 (date of the file write; the file records no retrieval date); (b) 2026-10-04 (field `retrieved` of the result file).
- **SHA-256:** none; API responses.
- **Licence:** section F has no row for PubMed or E-utilities; the NCBI policy page quoted in 3.1 is the only NCBI statement recorded.
- **You download it:** yes, through the two scripts (network access needed; `d2_index_variants.py` posts batches of 500 identifiers with 8 parallel requests for every PMID that is in neither `pmid_year.parquet` nor the cache; the first full fetch took about one hour according to a note made at the time, and the duration is not in a log). The stored result `g23_pubmed_metadata.json` is in the snapshot.

### 3.8 PeptideAtlas human phosphoproteome build 606 (feasibility check only)
- **Source:** https://peptideatlas.org/builds/human/phospho/202512/atlas_build_606.tsv.zip
- **Version:** human phospho build 606 (directory 202512); Last-Modified 16 Dec 2025.
- **Retrieved:** 2026-10-01.
- **SHA-256:** `d49b8af4e865b887dfc25ff8fcd9670a5d78cd0007dc1f69326983ef2e02b53e` (549,143,491 bytes).
- **Licence:** section F has no row for PeptideAtlas ("PhosphoSitePlus and PeptideAtlas are not used" in the tables; the build holds peptide-level evidence only, so no contrast was run).
- **You download it:** only to repeat `peptideatlas_feasibility.py` (549 MB; `python d_downloads.py atlas`). No result depends on it other than `peptideatlas_606_feasibility.json`.

## 4. Results of earlier stages that the Rev 6 scripts read (in the repository)

`study/earlier_results/` holds five small JSON files (`dataset_counts.json`, `matching_results.json`, `predictor_metrics.json`, `tool_versions.json`, `consequence_experiment_results.json`). They are aggregate statistics and are in the archive, which is deposited from a GitHub release. `b1_reconcile.py` also reads the earlier manuscript text and tables in `study/manuscript/` (not in the snapshot); its stored result is `study/analysis/results/b1_reconciliation.json`.

## 5. Predictors (block G11)

The predictor audit (`run_j.py`) reads score files that were produced by running three tools on the evaluation sample. The score files (`study/processed_data/epsd/predictor_audit/`) are site-level and EPSD-derived: not redistributed. The tools, with their weights, are not part of the snapshot either; download them from their repositories. Versions are from `study/earlier_results/tool_versions.json`; retrieval dates and hashes are not recorded.

| Tool | Repository and version scored | Licence (section F) |
|---|---|---|
| DeepPhos | https://github.com/USTC-HIlab/DeepPhos , checkout 6f94394 (2019-03-16) | "no LICENSE file (files: README.md, methods, models, predict.py, readme.txt, requirement.txt, sh.exe.stackdump, test data.csv, train.py); neither readme contains the word "licen"". Cell text: not stated (checked 2026-10-07). Do not redistribute. |
| MusiteDeep | https://github.com/duolinwang/MusiteDeep_web , checkout cca9f0c (2020-12-04); models Phosphoserine_Phosphothreonine and Phosphotyrosine | "MIT License ... Copyright (c) 2019 duolinwang ... Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction". Cell text: MIT. |
| PhosF3C | https://github.com/SkywalkerLuke/PhosF3C ; the tyrosine checkpoint linked from the repository README; checkout not recorded | "MIT License ... Copyright (c) 2024 SkywalkerLuke ... Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction". Cell text: MIT. |

- **Training-set overlap (J4):** `run_j.py` also reads `PhosF3C-main/data/merged_train_Y.fasta` from the PhosF3C repository (expected at `study/tools/phosf3c/PhosF3C-main/data/`).
- **You download them:** yes, only to recompute block G11. The stored result is `study/analysis/results/g11_predictors.json`.
