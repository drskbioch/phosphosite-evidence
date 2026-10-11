# Reproducing the analysis

This file gives the run order of the scripts in `study/scripts/analysis/` (the Rev 6 scripts), from the raw inputs to the stored results, and then the reports, the figures and the number checks. Inputs are listed in `INPUTS.md`; the contents of the archive, which is deposited from a GitHub release are in `RELEASE_SCOPE.md`.

There are two ways to use it.

| Level | What you need | What you get | Time |
|---|---|---|---|
| **Stored results** (sections 2 to 4) | the snapshot and the Python packages in `requirements.txt` | the two reports, the 13 composite figures, the tables, the number checks, the table check, the Supporting Information package and the 41 panel PNGs | reports and checks under 10 s; composite figures and tables 42 s; panels 15 s; Supporting Information 2 s (measured) |
| **Raw inputs** (sections 5 and 6) | the inputs in `INPUTS.md`, several of which are not in the repository | every file in `study/analysis/results/` | about 15 hours of computation if the steps run one after another (14.6 hours by the times in section 5, which include about 1 hour for the first fetch of publication years from PubMed and leave out the downloads), of which `run_propensity.py` is 6.2 hours and **`run_g21.py` 4.9 hours** |

**Steps that need inputs that are not in the repository** are marked in the "Needs" column of the stage tables in section 5 and collected in section 6. The EPSD-derived site tables cannot be downloaded; they are built by the earlier stages (section 7) from EPSD files and are not redistributed.

## 1. Environment

- Python 3.13.7 on Windows 10. Commands below are bash (Git Bash on Windows). Paths are computed from each script's own location, so nothing needs editing.
- Install: `pip install -r requirements.txt`. This is the only environment file. It was installed into a new Python 3.13.7 virtual environment on 2026-10-08 and every step in sections 2 to 4 was run in it. A conda environment file is not provided: conda was not available to test one. `requirements.txt` ends with two packages (PyMuPDF, openpyxl) that only the Supporting Information scripts, the limits checker and the release-scope scan use.
- **Versions.** The repository holds no lock file, no environment export and no `pip freeze` output (searched on 2026-10-11), and no result file records package versions. The installed set was compared with `requirements.txt` on 2026-10-11 and is the same (14 packages). The versions in `requirements.txt` are those installed in the Python 3.13.7 interpreter that produced the results; every one of them was installed before the first result file was written (2026-09-29), judged from the installation dates of the packages.
- Figures use the font Arial (`font.family` in `study/scripts/earlier_figures/figlib.py`). The stored PNG and TIFF files were drawn with Arial installed; the byte-for-byte match in section 8 was obtained on a machine that has it.
- Run times in this file are either **recorded** (the `seconds` fields written into the result files by the original runs, summed here) or **measured** on 2026-10-08 on an Intel Core i5-1235U (12 threads) with 24 GB of memory, with up to three jobs running together. Each process that loads the record-level table (`base.parquet`) needs about 1 GB of memory.
- Names used below: `PRIV` = `study/analysis/data_private/` (git-ignored caches and downloads), `RES` = `study/analysis/results/`, `REL` = `study/processed_data/epsd/release/`.

## 2. Reports (from stored results)

Run from `study/scripts/analysis`.

| Step | Command | Reads | Writes | Time |
|---|---|---|---|---|
| R1 | `python build_results.py` | `RES/*.json` (including `d2b/`, `prop/`, `f_runs/`), `grid_manifest.csv`, the Markdown files in `study/analysis/` | `study/analysis/RESULTS_REPORT.md`, `RES/manifest_status.csv` | 3 s (measured) |
| R2 | `python build_claims.py` | the same result files, `study/earlier_results/*.json` | `study/analysis/CLAIMS_VS_EVIDENCE.md`, `RES/claims_bindings.json` | 1 s (measured) |

`build_interim.py` also lives in this folder. It writes `INTERIM_REPORT.md`, a record of an intermediate stage, and is not part of the run order. Do not run it on top of the stored report.

## 3. Number checks (from stored results)

Run from `study/scripts/analysis`, after section 2.

| Step | Command | Checks | Writes | Time |
|---|---|---|---|---|
| C1 | `python check_numbers_sourced.py` | every number in `RESULTS_REPORT.md` and `CLAIMS_VS_EVIDENCE.md` has a source value in a result file | `RES/number_source_check.md` | 1 s (measured) |
| C2 | `python check_bindings.py` | every number in the claim rows is bound to a named field of a result file; then it runs the figure and table check (C4) | `RES/number_binding_check.md`, and the file of C4 | 1 s (measured) |
| C3 | `python check_literals.py` | every numeric literal that the three report builders write into text is classified and documented | `RES/builder_numeric_literals.md` | 1 s (measured) |
| C4 | `python check_figure_bindings.py` (also run at the end of C2) | every number printed in a figure or table is bound to a result field; derived values and counts are calculated again from the result files; the numbers in the SVG text of the 13 composite figures (which contain every panel) are bound | `RES/figure_binding_check.md` | 1 s (measured) |
| C5 | `python check_table_numbers.py --require-internal` | every number in Tables 1 to 3 and S1 to S5 and in their notes files is bound to a result value, a verified text cell, a claims-table binding, a count or a documented constant; Table S1 agrees with `manifest_status.csv` for all 534 contrasts (identifiers, estimates and interval limits at 4 significant digits, second intervals from the paired contrast); it runs C4 first | `RES/table_number_check.md`, `RES/table_number_binding.csv` | 1 s (measured) |
| C6 | `python scan_si.py` | banned-word scan of the Supporting Information (the two exemptions of `study/scripts/text_checks/scan_banned.py` apply) | `RES/si_banned_scan.md` | 1 s (measured) |

C5 compares the numbers that are not bound with the values of the earlier estimates in `study/analysis/internal/earlier_to_present_changes.csv`, a file that is not in the snapshot (it is excluded, `RELEASE_SCOPE.md` section 3; the table builder F3 writes it again). In the repository run C5 uses `--require-internal`, which makes a missing file an error (exit status 2). In the snapshot run C5 without the flag prints "SUPERSEDED-VALUE COMPARISON NOT RUN: internal file not in this tree" in its output, its report and its summary line, and runs every other test; it exits 0 only if those pass.

C1, C2 and C3 are the three number checks. Every one of them finishes with 0 unsourced and 0 unclassified numbers, from the snapshot alone and from the full working tree (section 8). Points about running them from the snapshot:

- C4 re-hashes the record-level input files listed in `RES/figure_counts.json` (`PRIV/base.parquet`, `PRIV/g23_tiers.parquet`, `PRIV/depth_variants.parquet`, `PRIV/index_variants.parquet` and four files in `REL`). When `PRIV` and `REL` are absent, as in the snapshot, the step is skipped: the report and the summary line say "input hashes: skipped: record-level data not present (9 of 9 files missing)", and every other part of C4 and of C2 runs. With the files present the line gives the number of files whose hash changed (0). C4 also reads `study/audit_reports/phaseH1_step3_predictor_audit.md` (one table cell cites it); the snapshot includes that file.
- C3 reads the run logs `RES/*.log` (the literal `0001` in `build_claims.py` is documented only there). The logs are tracked in the snapshot.
- C4 scans the SVG text of the composite figures in `study/manuscript/figures/`. Each panel is drawn from one of the composites, so the panels need no separate scan and C4 does not depend on panel files.

The files of the snapshot are checked for EPSD site-level records and predictor weights with `python study/scripts/release/check_release_scope.py`, run from anywhere; it reads the file list in `RELEASE_SCOPE.md` (section 4 of that file gives the result). PDF text is searched only if PyMuPDF is installed.

## 4. Figures and tables (from stored results)

Run from `study/scripts/analysis` for F1, and from `study/scripts/figures` for F2 to F4.

| Step | Command | Reads | Writes | Time |
|---|---|---|---|---|
| F1 | `python derive_figure_values.py` | `RES/*.json` | `RES/derived_figure_values.json` (ratios, shares and logarithms printed in the figures) | under 1 s (measured) |
| F2 | `python build_all.py composites` | `RES/*.json`, `study/earlier_results/*.json` | the 13 composite figures as PDF, SVG, PNG and TIFF in `study/manuscript/figures/` and `study/manuscript/earlier_figures/`; `RES/figure_captions.json`, `RES/figure_bindings.json` | 40 s (measured) |
| F3 | `python build_tables.py` | `RES/*.json`, `study/analysis/CLAIMS_VS_EVIDENCE.md`, `study/analysis/sources/SOURCES.md`, `study/analysis/table_wording.json`, `study/earlier_results/tool_versions.json` | Tables 1 to 3 and S1 to S5 (CSV and Markdown, with notes files and the key to Table S1) in `study/manuscript/tables/`; the table of changes between the two analyses in `study/analysis/internal/` (excluded from the release) | 2 s (measured) |
| F4 | `python build_all.py panels` | as F2 | the 41 panels as 600 dpi PNG files in `study/manuscript/figures/panels/` (PNG only), `panels/PANEL_MANIFEST.md`, `RES/panel_bindings.json` | 15 s (measured) |
| F5 | `python build_si.py` (run from `study/scripts/analysis`) | the composite PDFs `FigureS1` to `FigureS8` and the CSV files of Tables S1 to S5 | `study/manuscript/supporting_information/SupportingInformation.pdf`, `TableS1_contrast_grid.xlsx` (the PDF and the workbook carry a creation date, so they are not byte-identical between runs; their text and cell values are) | 1 s (measured) |

C4, C5 and C6 (section 3) belong after F2, F3 and F5, because they read the composites, the tables and the package. `SI_FILE_LIST.md` in the same folder is written by hand.

The limits of the target journal are in `study/analysis/submission_limits.yml` (all `TODO`); `python check_limits.py` (from `study/scripts/analysis`) measures the five main figures and three main tables against them and prints "limit not set" for every limit that is `TODO`.

**Expected result of sections 2 to 4, checked on 2026-10-11 in a copy of the snapshot** (details in section 8.4): the two reports, `claims_bindings.json`, `manifest_status.csv`, `derived_figure_values.json`, `figure_captions.json`, `figure_bindings.json`, all 25 table files, the number-check reports, `table_number_check.md`, `table_number_binding.csv`, `si_banned_scan.md`, `PANEL_MANIFEST.md`, `panel_bindings.json` and all 41 panel PNGs are byte-identical to the stored files. The composite PNG and TIFF files are byte-identical. The composite PDF and SVG files differ in the creation date and in the random element identifiers that matplotlib writes, and are identical once those are removed. `SupportingInformation.pdf` and `TableS1_contrast_grid.xlsx` differ in the creation date only (page text, bookmarks and cell values are identical). `figure_binding_check.md` differs in one line: the stored copy was written with the record-level files present and gives 0 changed files; from the snapshot alone the line reports the step as skipped.

## 5. From raw inputs to the stored results

All commands run from `study/scripts/analysis`. They overwrite the files they write in `RES` and `PRIV`. Order matters where a stage reads what an earlier stage wrote. "Needs": **E** = the EPSD-derived tables in `REL` (INPUTS.md section 1); **G** = `gene2pubmed.gz` (3.1); **D** = downloads in `PRIV/dl/` (3.2 to 3.8); **N** = network access; **P** = predictor score files and the PhosF3C training file (section 5 of INPUTS.md); **C** = caches written by an earlier stage; **M** = earlier manuscript files.

Scripts that are not run directly (imported by the others): `common.py` (paths, seed 20260924, 1,000 resamples, generators derived from the seed and the contrast ID), `base.py`, `engine.py`, `populations.py`, `propensity.py`, `d2b.py`, `cond_logit.py`, `g21_data.py`, `provenance.py`.

### Stage 0. Manifest (written once, stored)

`make_manifest.py`, `make_manifest_amendment.py`, `make_manifest_g21_g22.py` and `make_manifest_g23.py` wrote `study/analysis/grid_manifest.csv` (534 contrasts) from the plan, in that order, before any contrast of their blocks was run. The stored file is their output. Do not run `make_manifest.py` on top of it: it writes the file from scratch. The later manifest scripts append their blocks and do nothing, or refuse, if their rows already exist.

### Stage 1. Downloads

| Step | Command | Needs | Writes | Time |
|---|---|---|---|---|
| 1.1 | `python d_downloads.py pharos goa paxdb atlas idmap` | N | `PRIV/dl/*`, `PRIV/pharos_tdl.parquet`, `RES/downloads_manifest.json` (URL, version, date, size, SHA-256) | network-bound; the PeptideAtlas file is 549 MB |
| 1.2 | `python verify_downloads.py` | D | the final state of each download in `RES/downloads_manifest.json` | seconds |
| 1.3 | `python d_uniprot_reviewed.py` (`--out <folder>` makes a test download and leaves the stored files alone) | N | `PRIV/dl/rev_*.tsv`, `PRIV/dl/seq_*.tsv`, `PRIV/swissprot_reviewed.parquet`; one record per file (URL, release, date, SHA-256) in `RES/downloads_manifest.json` | network-bound |

### Stage 2. Gene, index, study-cluster and tier tables

| Step | Command | Needs | Reads | Writes | Time |
|---|---|---|---|---|---|
| 2.1 | `python d1_gene2pubmed.py` | G | `gene2pubmed.gz` | `PRIV/g2p_pairs_3sp.parquet`, `PRIV/g2p_genes_per_pmid.parquet`, `RES/d1_gene2pubmed_summary.json` | 60 s (measured) |
| 2.2 | `python d2_index_variants.py` | E, C, N | 2.1 outputs, `REL/pmid_year.parquet`, `REL/sei_protein.parquet`, PubMed through E-utilities | `PRIV/index_variants.parquet` (index variants V0, V20, V50, V100, VT2025, VT2026), `PRIV/pmid_year_gene.parquet`, `RES/d2_index_variants_summary.json` | network-bound; the first fetch of publication years took about 1 hour (note made at the time) |
| 2.3 | `python gene_unit.py` | E | `REL/sei_protein.parquet` | `PRIV/gene_unit.parquet`, `RES/gene_unit_summary.json` | 5 s (measured) |
| 2.4 | `python b1_clusters.py` | E | `REL/mammalian_sites.parquet`, `REL/cluster_split.parquet` | `PRIV/depth_variants.parquet`, `RES/b_cluster_rule_info.json`, `RES/b_largest_clusters.json` | 271 s (measured) |
| 2.5 | `python b1b_consortium.py` (after 2.4, before 2.8 and stage 4) | E, C | `PRIV/depth_variants.parquet`, `REL/mammalian_sites.parquet` | adds the consortium-removal columns to `PRIV/depth_variants.parquet`; adds the counts to `RES/b_cluster_rule_info.json` | 93 s (measured) |
| 2.6 | `python base.py` (also built on first use by any later script) | E, C | `REL/mammalian_sites.parquet`, `kinase_motif_sequence.parquet`, `cluster_split.parquet`, `human_sites_source_set_a.parquet`; `PRIV/index_variants.parquet` | `PRIV/base.parquet` (1,768,068 records) | 7 s (measured) |
| 2.7 | `python e_abundance_map.py` | D, C | `PRIV/dl/px_*.txt`, `PRIV/dl/*_idmapping.dat.gz`, `PRIV/base.parquet` | `PRIV/abundance.parquet`, `RES/e_abundance_coverage.json` | 7 s (measured) |
| 2.8 | `python f_tiers.py` | E, C | `PRIV/depth_variants.parquet`, `REL/mammalian_sites.parquet` | `PRIV/f_tiers.parquet`, `RES/f_tiers_summary.json` | 8 s (measured) |
| 2.9 | `python g23_fetch.py`, then `python g23_list.py`, then `python g23_tiers.py`, then `python g23_kept_pmids.py` | E, C, N (g23_fetch only) | `REL/*`, PubMed through E-utilities (g23_fetch), then each step reads the one before | `RES/g23_pubmed_metadata.json`; `study/analysis/nonprimary_pmid_list.md` and `RES/g23_nonprimary_pmids.json`; `PRIV/g23_tiers.parquet`, `RES/g23_tiers_summary.json`; `RES/g23_kept_nonevidence_pmids.json` | g23_fetch: network (690 records); g23_list: 1 s; g23_tiers: 100 s; g23_kept_pmids: 29 s (measured) |
| 2.10 | `python b1_reconcile.py` | E, C, M | `PRIV/depth_variants.parquet`, `REL/temporal_split.parquet`, `RES/b_cluster_rule_info.json`, `study/manuscript/manuscript.md`, `study/manuscript/tables/_rev5/table1_matched_comparisons.md` | `RES/b1_reconciliation.json`, `study/analysis/RECONCILIATION.md` (sections 1 to 5) | 2 s (measured) |

`b1_reconcile.py` writes sections 1 to 5 of `RECONCILIATION.md` and keeps everything from the heading "## 6." on unchanged. Sections 6 and 7 are written by hand; the file names the result files each is drawn from, and says which counts come from direct queries that are not stored in a result file.

`g23_list.py` writes the list of non-primary PMIDs that was reviewed before `g23_tiers.py` was run (plan, deviation log of 2026-10-04).

### Stage 3. Validation of the estimators

| Step | Command | Needs | Writes | Time |
|---|---|---|---|---|
| 3.1 | `python validate_engine.py` | E, C | `RES/validate_engine.json` | 10 s (measured) |
| 3.2 | `python validate_engine_modes.py` | C | `RES/validate_engine_modes.json` | 5 s (measured) |
| 3.3 | `python validate_engine_v50.py` | C | `RES/validate_engine_v50.json` | 15 s (measured) |
| 3.4 | `python validate_propensity.py` | C | `RES/validate_propensity.json` | 101 s (measured) |
| 3.5 | `python check_cond_logit.py` | none | prints the agreement of the own conditional logit with statsmodels (simulated data) | 7 s (measured) |
| 3.6 | `python check_cond_logit_real.py` | C | `RES/g21_cond_logit_check.json` (required by stage 4, step 4.14) | 31 s (measured) |
| 3.7 | `python check_propensity_order.py` | C | `RES/propensity_order_dependence.json` | 29 s (measured) |

3.1 to 3.3 import the package source from `study/phosnet/src` (relative path `../../phosnet/src`) to compare the estimates.

### Stage 4. Contrasts

| Step | Command | Needs | Reads | Writes | Time |
|---|---|---|---|---|---|
| 4.1 | `for b in 1 2 3 4 12 13 15; do python run_cem_grid.py "G$b"; done` | C | `PRIV/base.parquet`, `PRIV/depth_variants.parquet`, `grid_manifest.csv` | `RES/cem_results.json` (blocks merged in) | 14 min recorded (626 s measured) |
| 4.2 | `python run_g14.py` | C | `PRIV/base.parquet`, `study/earlier_results/matching_results.json` | `RES/g14_cem.json` | 159 s (measured) |
| 4.3 | `python run_propensity.py $(ls ../../analysis/results/prop \| sed 's/\.json$//')` | C | `PRIV/base.parquet` | `RES/prop/<id>.json` (24 files, among them the runs labelled post-hoc) | 6.2 hours recorded for the 18 files that were fitted; 442 s measured for one file (G14-0337, recorded 675 s) |
| 4.4 | `python b4_twoway.py`; `python b5_index_effects.py` | C | `PRIV/base.parquet`, `PRIV/index_variants.parquet`, `study/earlier_results/matching_results.json` | `RES/b4_twoway.json`; `RES/b5_index_effect_sizes.json` | 495 s; 206 s (measured) |
| 4.5 | `python d3_d4_d5_index_checks.py` | D, C | `PRIV/base.parquet`, `PRIV/index_variants.parquet`, `PRIV/pharos_tdl.parquet`, `PRIV/dl/goa_human.gaf.gz` | `RES/d3_d4_d5_index_checks.json` | 391 s (measured) |
| 4.6 | `python run_g5.py` | C | `PRIV/base.parquet`, `PRIV/abundance.parquet` | `RES/g5_abundance.json`, `RES/g5_diagnostics.json` | 210 s (measured) |
| 4.7 | `for b in 16 17 19; do python run_amend.py "G$b"; done` | C | `PRIV/base.parquet`, `PRIV/abundance.parquet`, `PRIV/swissprot_reviewed.parquet` | `RES/amend_results.json` (blocks merged in), `RES/g16_coverage.json` | 10.5 min recorded (585 s measured) |
| 4.8 | `python run_d2b.py $(ls ../../analysis/results/d2b \| sed 's/\.json$//')` | C | `PRIV/base.parquet` | `RES/d2b/<id>.json` (11 files) | 0.56 hours recorded; 151 s measured for G18-0407 (recorded 185 s) |
| 4.9 | `python g20_diagnostics.py` (after 4.1 and 4.7) | C | `RES/cem_results.json`, `RES/amend_results.json`, `PRIV/base.parquet` | `RES/g20_diagnostics.json` | 202 s (measured) |
| 4.10 | `python run_f.py G6 S1,S2,S3,S4`; `python run_f.py G7 S1,S2`; `python run_f.py G7 S3,S4`; `python run_f.py G8 S1,S2,S3,S4` | C | `PRIV/f_tiers.parquet`, `PRIV/base.parquet` | `RES/f_runs/<block>_<signals>.json` | 44 min recorded in all (G6 151 s, G7 576 s and 996 s, G8 903 s); 36 min measured (123 s, 400 s, 818 s, 798 s) |
| 4.11 | `python run_j.py` | E, P, C | `study/processed_data/epsd/predictor_audit/eval_sample_with_sequence.parquet` and the predictor score files in the same folder, `study/tools/phosf3c/PhosF3C-main/data/merged_train_Y.fasta`, `study/earlier_results/predictor_metrics.json`, `PRIV/gene_unit.parquet`, and `study/scripts/earlier_predictor_audit/step3/score_metrics.py` (imported) | `RES/g11_predictors.json` | 38 s (measured) |
| 4.12 | `python h_measures.py uniprot` | D, C | `PRIV/dl/seq_9606.tsv`, `PRIV/abundance.parquet`, `PRIV/base.parquet` | `RES/g9_measures.json`, `RES/g9_intersections.json` | 20 s (measured) |
| 4.13 | `python run_g22.py` | C | `PRIV/base.parquet` | `RES/g22_results.json` | 28 s (measured) |
| 4.14 | `python run_g21.py` (**many hours**) | D, C | `PRIV/base.parquet`, `PRIV/dl/seq_*.tsv`, `RES/g21_cond_logit_check.json` | `RES/g21_results.json` (24 contrasts, written group by group) | **4.9 hours recorded** (8 groups; two groups took 4,847 s and 4,867 s). A group that is interrupted restarts from zero; groups already in `g21_results.json` are skipped. Use `nohup` or a background job. |
| 4.15 | `python g21_sequence_coverage.py` | D, C | `PRIV/dl/seq_*.tsv`, `PRIV/base.parquet` | `RES/g21_sequence_coverage.json` | 69 s (measured) |
| 4.16 | `python run_g23.py` | C | `PRIV/g23_tiers.parquet`, `PRIV/base.parquet` | `RES/g23_results.json` | 2.6 min recorded (130 s measured) |
| 4.17 | `python peptideatlas_feasibility.py` | D | `PRIV/dl/atlas_build_606.tsv.zip`, `RES/peptideatlas_606_inspection.json` | `RES/peptideatlas_606_feasibility.json` | not timed (549 MB archive) |

The Task H and I contrasts that need PhosphoSitePlus (`h_measures.py psp_lt psp_cat`, the Task I runner) were not run: the source file was not available (INPUTS.md section 2). `RES/peptideatlas_606_inspection.json` and `RES/peptideatlas_606_web_check.json` have no script in the repository; they are stored results. `RES/prop_ids.txt` and `RES/prop_ids2.txt` list the propensity contrast IDs of two runs.

**Outputs that differ from the stored files when the current scripts run.** `run_cem_grid.py` (4.1) writes two fields that the stored `cem_results.json` does not have: `n_gene_clusters` for all 66 contrasts and `smd_pre` for 43. `g20_diagnostics.py` (4.9) writes four more entries than the stored `g20_diagnostics.json` in `prematching_balance`: `G19-0461`, `G19-0463`, `G19-0465` and `G19-0467`. Every value that is in the stored files is identical.

### Stage 5. Counts for the figures

| Step | Command | Needs | Reads | Writes | Time |
|---|---|---|---|---|---|
| 5.1 | `python count_descriptives.py` | E, C | `PRIV/base.parquet`, `PRIV/depth_variants.parquet`, `PRIV/g23_tiers.parquet`, `PRIV/index_variants.parquet`, `REL/*.parquet`, `RES/b5_index_effect_sizes.json`, `study/earlier_results/dataset_counts.json` | `RES/figure_counts.json` (counts and SHA-256 of each input file, cross-check against `dataset_counts.json`) | 12 s (measured) |

Then the reports (section 2), the number checks (section 3) and the figures and tables (section 4).

## 6. Which steps need inputs that are not in the repository

| Input | Needed by | Obtainable |
|---|---|---|
| EPSD-derived tables in `REL` and the evaluation sample | steps 2.2 to 2.6, 2.8, 2.9, 2.10, 3.1, 4.11 and 5.1 | Not downloadable; built from the EPSD files by the earlier stages, which the repository states are not reproducible from it. Not redistributed (licence not stated by the provider; permission requested). |
| `gene2pubmed.gz` | 2.1 | Download (INPUTS.md 3.1). |
| UniProt ID mapping, PaxDb, Pharos, GOA, PeptideAtlas | 1.1, 2.7, 4.5 and 4.17 | Download with `d_downloads.py`. |
| UniProt reviewed lists and sequences | 4.7 (through `swissprot_reviewed.parquet`), 4.12, 4.14 and 4.15 | Download with `d_uniprot_reviewed.py` (a fresh download on 2026-10-08 gave the stored hashes; release 2026_03). |
| PubMed through E-utilities | 2.2 and 2.9 | Network access. |
| Predictor score files, predictor repositories, PhosF3C training file | 4.11 | Run the three tools on the evaluation sample; DeepPhos has no stated licence. |
| Earlier manuscript files | 2.10 (`b1_reconcile.py`) | `study/manuscript/manuscript.md` and `study/manuscript/tables/_rev5/table1_matched_comparisons.md` in the repository; not in the snapshot. |
| PhosphoSitePlus file | not run | Not obtained. |

Everything in sections 2 to 4 needs none of these, except that C2 and C4 need the record-level input files named in section 3.

## 7. Earlier stages, not needed for Rev 6 results

These folders hold the scripts of earlier stages. They are not edited by the Rev 6 work and are not needed to produce, check or draw any Rev 6 result:

- `study/ingest/` (EPSD parsing, study clustering, site database)
- `study/scripts/epsd_phase0/` (the four validation signals, pilots and diagnostics)
- `study/scripts/earlier_predictor_audit/` (index, temporal split, predictor audit)
- `study/scripts/reproduction/`, `study/scripts/internal/`, `study/scripts/text_checks/`
- `study/scripts/earlier_figures/` (earlier figure and table scripts)

The web application and its deployment folders at the repository root are unrelated to the results and are not part of the snapshot.

Four exceptions, because a Rev 6 script imports or reads them (they are in the snapshot, unchanged):

- `study/scripts/earlier_figures/figlib.py` (figure style, imported by the Rev 6 figure code),
- `study/scripts/earlier_predictor_audit/step3/score_metrics.py` (score loaders, imported by `run_j.py`). Its three default paths (the evaluation sample, the output folder of the predictor scores and the PhosF3C result file) are strings relative to the project root; `run_j.py` sets the output folder and the PhosF3C path from its own location, and the evaluation sample is read by the earlier-stage `main()` of this module only. In this archive these defaults replace the absolute paths of the machine on which the earlier stages were run (see `PROVENANCE.md`, Redactions),
- `study/phosnet/src/` (imported by `validate_engine_modes.py` and `validate_engine_v50.py` only; the whole package is in the snapshot).
- `study/scripts/text_checks/scan_banned.py` and `banned_words.txt` (read by `scan_si.py`, step C6).

## 8. What was checked, and what differs

Sections 8.1 to 8.3 are the checks of 2026-10-08. They were made before the figures and tables were regenerated and renumbered; the counts of figure and table numbers in them (795 verified, 809 in the SVG text, 19 table files) are those of the earlier files. Section 8.4 is the check of the present files and gives their counts.

### 8.1 Stored-results steps from the snapshot alone

The 339 files of the snapshot (`RELEASE_SCOPE.md`) were copied to a scratch folder without the data folders. `requirements.txt` was installed into a new virtual environment. The steps ran in the order R1, R2, C1, F1, F2, F3, F4, C2, C3, C4, followed by `check_release_scope.py`. Every step exited normally. The outputs were compared with the stored files.

| Step | Result |
|---|---|
| R1, R2 | `RESULTS_REPORT.md`, `manifest_status.csv`, `CLAIMS_VS_EVIDENCE.md`, `claims_bindings.json`: byte-identical to the stored files. |
| C1 | 4,060 numbers examined from 87 result files: 4,013 equal to a source value, 47 explained by derivation, **0 unexplained**. `number_source_check.md` byte-identical. |
| C2 | 39 claim rows, 331 numbers: 300 bound, 31 derived, **0 unsourced**; 44 headline items bound. `number_binding_check.md` byte-identical. The input-hash step of the figure check was skipped (record-level data not present) and the run finished. |
| C3 | 269 literals: 88 design definitions, 5 format fragments, 13 list numbers, 35 provenance, 95 quoted statements, 33 statistical conventions; **0 unclassified or undocumented**. `builder_numeric_literals.md` byte-identical. |
| C4 | 795 bound numbers verified, 0 not reproduced; derived values 0 differ; counts cross-check 26 of 26; text cells 195 verified, 0 not; 809 numbers in the SVG text of the composites, **0 unsourced**; input-hash step skipped. `figure_binding_check.md` differs from the stored copy in the input-hash statement only (the stored copy gives 0 changed files). |
| F1, F3 | `derived_figure_values.json` and all 19 table files: byte-identical. |
| F2 | `figure_captions.json` and `figure_bindings.json`: byte-identical. The 13 composite PNG and TIFF files: byte-identical; the PDF and SVG files: identical after removing the creation date and the element identifiers. |
| F4 | All 41 panel PNG files, `PANEL_MANIFEST.md` and `panel_bindings.json`: byte-identical. |

### 8.2 The same steps in the full working tree

R1, R2, C1, C2 (with C4), C3 were run in the repository with the data folders present. Reports unchanged. C1: 0 unexplained. C2: 0 unsourced. C3: 0 unclassified. C4: 795 verified, 0 not reproduced, input hashes 0 changed, 809 numbers in the SVG text, 0 unsourced.

### 8.3 Steps from raw inputs (section 5)

Scripts of stages 2 to 5 were run in a scratch copy of the scripts, with copies of the input files and of the stored results, and the outputs compared with the stored files.

- **Identical to the stored files** (timing fields such as `seconds` excluded; NaN compared as equal): the outputs of steps 2.1, 2.3, 2.4, 2.5, 2.7, 2.8, 2.9 (all but `g23_fetch.py`), 2.10 (`b1_reconciliation.json`; sections 1 to 5 of `RECONCILIATION.md`, which the script now writes while keeping sections 6 and 7), 3.1 to 3.7, 4.1 (all seven blocks), 4.2, 4.4 to 4.13, 4.15, 4.16 and 5.1. `g11_predictors.json` from 4.11 is byte-identical. The caches written by 2.4, 2.5, 2.6 and 2.9 (`depth_variants.parquet`, `base.parquet`, `g23_tiers.parquet`) have the SHA-256 values recorded in `RES/figure_counts.json`.
- **Downloads checked:** the six UniProt files (1.3) from a fresh run of `d_uniprot_reviewed.py` have the stored SHA-256 values (release 2026_03 still served); the three EPSD archives were downloaded again (INPUTS.md 1.1), the human one has the recorded partial hash, and the row counts of the three text files equal the site counts.
- **Run in part:** 4.3 (1 of 24 files, `G14-0337`: identical except the projected run time), 4.8 (1 of 11 files, `G18-0407`) and 4.14 (1 of 8 groups, signal S2 in P-human: 227 s measured, 220 s recorded; identical).
- **Differences:** the two extra fields in `cem_results.json` and the four extra entries in `g20_diagnostics.json` (section 5, stage 4); every stored value is identical.
- **Not run:** 1.1 and 1.2 (downloads), 2.2 and `g23_fetch.py` (network), 4.3 (23 files), 4.8 (10 files), 4.14 (7 groups), 4.17 (archive not copied), and the `make_manifest*.py` scripts (the manifest is stored).

### 8.4 The present snapshot, from the stored results (2026-10-11)

The 364 files of the snapshot (`RELEASE_SCOPE.md`) were copied to a scratch folder without the data folders. With the installed set of `requirements.txt` (section 1), these steps ran in this order and all exited normally: R1, R2, C1, F1, F2, F3, F4, F5, C2 (with C4), C3, C5, C6, `check_limits.py`, `check_release_scope.py`. The outputs were compared with the stored files.

| Step | Result |
|---|---|
| R1, R2, F1, F3, F4 | byte-identical: `RESULTS_REPORT.md`, `manifest_status.csv`, `CLAIMS_VS_EVIDENCE.md`, `claims_bindings.json`, `derived_figure_values.json`, all 25 table files, all 41 panel PNGs, `PANEL_MANIFEST.md`, `panel_bindings.json`. The table of changes between the two analyses (`study/analysis/internal/`, excluded from the snapshot) was written again and is identical to the stored copy. |
| F2 | `figure_captions.json`, `figure_bindings.json` and the 13 composite PNG and TIFF files byte-identical; the PDF and SVG files identical in page text, page size and number of drawing objects (PDF) and in the SVG text after removing the creation date and the element identifiers. |
| F5 | `SupportingInformation.pdf`: 14 pages, page text and bookmarks identical; `TableS1_contrast_grid.xlsx`: both sheets, all cell values identical; the files differ in bytes by the creation date. |
| C1 | 4,060 numbers, 47 explained by derivation, **0 unexplained**; `number_source_check.md` byte-identical. |
| C2, C4 | 39 claim rows, 331 numbers: 300 bound, 31 derived, **0 unsourced**; figure and table step: 793 bound numbers verified, 0 not reproduced, derived values 0 differ, counts cross-check 26 of 26, 196 text cells verified, 800 numbers in the SVG text, **0 unsourced**; input hashes skipped (record-level data not present). `number_binding_check.md` byte-identical; `figure_binding_check.md` differs in the input-hash statement only. |
| C3 | 269 literals, **0 unclassified or undocumented**; `builder_numeric_literals.md` byte-identical. |
| C5 | 2,478 numbers in 15 table and notes files, **0 unbound**; Table S1: 534 rows, 0 fields differ, 0 second intervals without a paired contrast; reports byte-identical (this run followed F3, which writes `study/analysis/internal/` again; run on the unzipped snapshot without F3 and without the flag, C5 prints "SUPERSEDED-VALUE COMPARISON NOT RUN: internal file not in this tree" and its report carries that line). |
| C6 | 0 hits in the four parts of the Supporting Information; `si_banned_scan.md` byte-identical. |
| `check_limits.py` | runs; every limit prints "limit not set". |
| `check_release_scope.py` | 364 files found, 0 missing, no bad file type, no private path, no site-level pattern, no EPSD identifier, no CSV or workbook header with a UniProt and a position column, no JSON list of 500 or more numbers. |

Everything else in the snapshot is unchanged by these steps (29 files differ in bytes: 26 composite PDF and SVG files, `figure_binding_check.md`, `SupportingInformation.pdf` and `TableS1_contrast_grid.xlsx`, for the reasons given above).

