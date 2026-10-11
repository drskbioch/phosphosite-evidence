# Release scope

Exactly what the archive contains: the files listed in section 6, with these paths. The archive is deposited from a GitHub release; its version is 1.0.0 (`CITATION.cff`, `.zenodo.json`). Sizes are in bytes or MB (1 MB = 1,000,000 bytes). Nothing has been uploaded or deposited.

## 1. Contents and sizes

| Group | Files | Size | Paths |
|---|---|---|---|
| Rev 6 scripts | 67 | 0.54 MB | `study/scripts/analysis/*.py` (`__pycache__` is left out), including the table check, the Supporting Information build and scan, and the limits checker |
| Reference copies of the checkers at the freeze | 4 | 0.04 MB | `study/scripts/analysis/checkers_at_freeze/`: the three number checkers as they were at commit 3521d5c, for reference; not run |
| Plan, deviation log and wording | 3 | 0.08 MB | `study/analysis/ANALYSIS_PLAN.md` (the deviation log is its section 15), `table_wording.json` (cells whose wording differs from the stored text), `submission_limits.yml` (journal limits, all `TODO`) |
| Manifests | 3 | 0.24 MB | `study/analysis/grid_manifest.csv`, `study/analysis/results/manifest_status.csv`, `study/analysis/results/downloads_manifest.json` |
| Stored results | 122 | 1.79 MB | everything else in `study/analysis/results/`: result JSON files (including `d2b/`, `prop/`, `f_runs/`), the CSV and Markdown reports of the number checks and of the table and Supporting Information checks, and the run logs (`*.log`) |
| Claims and results files | 2 | 0.14 MB | `study/analysis/CLAIMS_VS_EVIDENCE.md`, `study/analysis/RESULTS_REPORT.md` |
| Figure and table code | 9 | 0.13 MB | `study/scripts/figures/*.py` except `list_old_figures.py` |
| Composite figures, all formats | 52 | 24.09 MB | `study/manuscript/figures/Figure*.{pdf,svg,png,tif}` (13 figures, 4 formats each) |
| Panel PNGs and panel manifest | 42 | 10.04 MB | `study/manuscript/figures/panels/Fig*.png` (PNG only) and `PANEL_MANIFEST.md` |
| `phosnet` package | 10 | 0.04 MB | `study/phosnet/`: source, tests, `LICENSE`, `README.md`, `pyproject.toml`, `CHANGELOG.md` (`validate_engine_modes.py` and `validate_engine_v50.py` import its source) |
| Root documents, requirements, licence and citation files | 14 | 0.13 MB | `README.md`, `PROVENANCE.md`, `TITLE.txt`, `.gitignore`, `scripts/set_title.py`, `REPRODUCE.md`, `INPUTS.md`, `RELEASE_SCOPE.md`, `requirements.txt`, `LICENSE`, `CITATION.cff`, `.zenodo.json`, `study/scripts/release/*.py` |
| Tables written by the table code | 25 | 0.77 MB | `study/manuscript/tables/Table*.{csv,md}`: Tables 1 to 3 and S1 to S5, the key to Table S1 and the notes files; the figure and table check reads them |
| Supporting Information package | 3 | 0.62 MB | `study/manuscript/supporting_information/`: `SupportingInformation.pdf`, `TableS1_contrast_grid.xlsx`, `SI_FILE_LIST.md` |
| Files the number checks read | 10 | 0.17 MB | `study/analysis/INTERIM_REPORT.md`, `RECONCILIATION.md`, `sources/SOURCES.md`, `nonprimary_pmid_list.md`; `study/audit_reports/phaseH1_step3_predictor_audit.md`; the five result files `study/earlier_results/*.json` |
| Modules a listed script imports | 4 | 0.02 MB | `study/scripts/earlier_figures/figlib.py` (figure style), `study/scripts/earlier_predictor_audit/step3/score_metrics.py` (score loaders for `run_j.py`), `study/scripts/text_checks/scan_banned.py` and `banned_words.txt` (read by `scan_si.py`) |
| **Total in the snapshot** | **370** | **38.85 MB** | |

The number checks are `check_numbers_sourced.py`, `check_bindings.py`, `check_literals.py` and `check_figure_bindings.py` (in the scripts group) and their reports `number_source_check.md`, `number_binding_check.md`, `builder_numeric_literals.md` and `figure_binding_check.md` (in the stored results group). The table check is `check_table_numbers.py` with its reports `table_number_check.md` and `table_number_binding.csv`. They are part of the groups above, not additional files.

The figure code also writes the composite figures to `study/manuscript/earlier_figures/`; that second copy is not in the snapshot.

## 2. Run logs

- The 29 run logs in `study/analysis/results/` match `.gitignore` (`*.log`) and are tracked by force. `check_literals.py` reads them: without them it lists one more literal as undocumented.
- Redaction: five logs (`f_G8_S1-S2.log`, `f_G8_S3-S4.log`, `run_propensity.log`, `run_propensity2.log`, `run_propensity3.log`) contained Python tracebacks with the local Windows user path. In them the path `C:\Users\<name>` is replaced by `<user-home>`. Nothing else in the logs was changed.
- Files in the list that git does not track yet: none

## 3. EXCLUDED: not in the snapshot

`check_release_scope.py` fails if a path in the file list of section 6 starts with one of these prefixes or lies inside a nested repository (a folder with a `.git` entry). The prefixes are the list `EXCLUDED` in that script.

| Excluded path | Reason |
|---|---|
| `legacy/` | pre-pivot product: web application, deployment folders and server notes, earlier scripts |
| `legacy/deploy/` | deployment folders and server notes |
| `legacy/deploy/deepphos-frontend/` | nested repository |
| `frontend/` | web application; nested repository |
| `study/analysis/internal/` | internal table of the changes between the earlier and the present estimates |
| `study/_archive/` | archive of moved files |
| `.claude/` | local Claude Code settings |
| `study/tools/` | predictor code and weights |
| `study/processed_data/` | EPSD-derived tables |
| `study/analysis/data_private/` | caches and downloads |

Nested repositories (git links or folders with a `.git` entry): `frontend/`, `legacy/deploy/deepphos-frontend/`, and the predictor checkouts inside `study/_archive/`. Their uncommitted changes are not part of any commit of this work.

Also not in the snapshot:

- `study/manuscript/earlier_figures/final/Figures.pptx` (the author's own file) and everything else under `study/manuscript/` except the figures, panels, tables and Supporting Information package named above.
- `study/scripts/figures/list_old_figures.py` (lists figure files of earlier stages), `study/analysis/PACKAGE_STATUS.md`, `SUBMISSION_STATE.md`, `SUBMISSION_READINESS.md`, `STATEMENT_FACTS.md`, the contact-sheet and review PDFs in `study/analysis/`, `study/analysis/FIGURE_INPUTS.md`, `FIGURE_TABLE_INVENTORY.md`, `OLD_FIGURES_LIST.md`.
- The scripts of earlier stages (`REPRODUCE.md`, section 7) other than the modules named in the last group of section 1, and the web application folders.

## 4. Scan for EPSD site-level records and predictor weights

`python study/scripts/release/check_release_scope.py` reads the file list in section 6 of this file and scans every file in it: 370 files, 38.85 MB. It searched 286 text and SVG files, 1 workbook by its cell values, and the text of 14 PDF files. The 67 PNG and TIFF files cannot be searched as text; they are drawn by the same code as the SVG and PDF versions, which were searched.

| Test | Result |
|---|---|
| File types that hold records or model weights (parquet, feather, h5, npz, npy, pkl, pt, pth, onnx, ckpt, safetensors, bin, sqlite, db, gz, zip, tsv, joblib); workbooks (xlsx) are searched by content | 0 files |
| A listed path that starts with an excluded prefix (section 3) or lies inside a nested repository | 0 files |
| A UniProt accession next to a residue letter (S, T or Y) and a position, in either order, with the usual separators | 0 files |
| EPSD identifiers (`EPSD` followed by digits) | 0 files |
| CSV files and workbook sheets whose header has both a UniProt column and a position column | 0 files |
| JSON lists of 500 or more numbers (weights or per-site vectors) | 0 lists |
| Amino-acid runs of 20 or more letters | 2 files (`study/scripts/analysis/run_j.py`, `study/scripts/release/check_release_scope.py`): the string of the 20 amino-acid letters (an encoding table in `run_j.py`, a search pattern in the checker), not a sequence from a protein |
| Accession-like tokens outside the tests above | 3 files: `study/scripts/figures/fig1_2.py`: hex colour codes; `study/phosnet/tests/test_smoke.py`: not classified; `study/audit_reports/phaseH1_step3_predictor_audit.md`: one accession in an example output line of the DeepPhos test file (accession, position, probability), not an EPSD record |

The largest files are `FigureS2.tif` 1.46 MB, `Figure4.tif` 1.38 MB, `Figure2.tif` 1.36 MB.

The stored result files hold aggregate statistics: estimates, intervals, counts, and lists of PMIDs and of PubMed titles. Where they name proteins (for example the largest gene clusters) they do so in summary form, never as records of sites. The scan does not read meaning; it shows that no table of site records and no file of weights is present.

## 5. File names

Some file names in the list contain words that the figure and manuscript text avoids: the contrast identifiers of the first three blocks in `study/analysis/results/prop/`, `d2_index_variants_summary_proxy5_superseded.json`, and two run logs. The names are the keys used in the manifests and the scripts and are left as they are.

## 6. File list

| Group | Path | Bytes |
|---|---|---|
| Rev 6 scripts | `study/scripts/analysis/b1_clusters.py` | 3,451 |
| Rev 6 scripts | `study/scripts/analysis/b1_reconcile.py` | 20,022 |
| Rev 6 scripts | `study/scripts/analysis/b1b_consortium.py` | 2,230 |
| Rev 6 scripts | `study/scripts/analysis/b4_twoway.py` | 3,504 |
| Rev 6 scripts | `study/scripts/analysis/b5_index_effects.py` | 4,949 |
| Rev 6 scripts | `study/scripts/analysis/base.py` | 3,297 |
| Rev 6 scripts | `study/scripts/analysis/build_claims.py` | 44,997 |
| Rev 6 scripts | `study/scripts/analysis/build_interim.py` | 40,778 |
| Rev 6 scripts | `study/scripts/analysis/build_release_snapshot.py` | 3,683 |
| Rev 6 scripts | `study/scripts/analysis/build_results.py` | 100,878 |
| Rev 6 scripts | `study/scripts/analysis/build_si.py` | 8,221 |
| Rev 6 scripts | `study/scripts/analysis/check_bindings.py` | 11,790 |
| Rev 6 scripts | `study/scripts/analysis/check_cond_logit.py` | 1,926 |
| Rev 6 scripts | `study/scripts/analysis/check_cond_logit_real.py` | 2,798 |
| Rev 6 scripts | `study/scripts/analysis/check_figure_bindings.py` | 15,635 |
| Rev 6 scripts | `study/scripts/analysis/check_limits.py` | 5,802 |
| Rev 6 scripts | `study/scripts/analysis/check_literals.py` | 5,811 |
| Rev 6 scripts | `study/scripts/analysis/check_numbers_sourced.py` | 20,565 |
| Rev 6 scripts | `study/scripts/analysis/check_propensity_order.py` | 1,774 |
| Rev 6 scripts | `study/scripts/analysis/check_table_numbers.py` | 16,667 |
| Rev 6 scripts | `study/scripts/analysis/common.py` | 4,017 |
| Rev 6 scripts | `study/scripts/analysis/cond_logit.py` | 5,180 |
| Rev 6 scripts | `study/scripts/analysis/count_descriptives.py` | 12,141 |
| Rev 6 scripts | `study/scripts/analysis/d1_gene2pubmed.py` | 1,627 |
| Rev 6 scripts | `study/scripts/analysis/d2_index_variants.py` | 10,140 |
| Rev 6 scripts | `study/scripts/analysis/d2b.py` | 2,788 |
| Rev 6 scripts | `study/scripts/analysis/d3_d4_d5_index_checks.py` | 7,719 |
| Rev 6 scripts | `study/scripts/analysis/d_downloads.py` | 6,099 |
| Rev 6 scripts | `study/scripts/analysis/d_uniprot_reviewed.py` | 3,717 |
| Rev 6 scripts | `study/scripts/analysis/derive_figure_values.py` | 9,831 |
| Rev 6 scripts | `study/scripts/analysis/e_abundance_map.py` | 3,085 |
| Rev 6 scripts | `study/scripts/analysis/engine.py` | 9,825 |
| Rev 6 scripts | `study/scripts/analysis/f_tiers.py` | 2,388 |
| Rev 6 scripts | `study/scripts/analysis/g20_diagnostics.py` | 7,488 |
| Rev 6 scripts | `study/scripts/analysis/g21_data.py` | 3,325 |
| Rev 6 scripts | `study/scripts/analysis/g21_sequence_coverage.py` | 979 |
| Rev 6 scripts | `study/scripts/analysis/g23_fetch.py` | 2,162 |
| Rev 6 scripts | `study/scripts/analysis/g23_kept_pmids.py` | 1,387 |
| Rev 6 scripts | `study/scripts/analysis/g23_list.py` | 7,015 |
| Rev 6 scripts | `study/scripts/analysis/g23_tiers.py` | 5,547 |
| Rev 6 scripts | `study/scripts/analysis/gene_unit.py` | 4,680 |
| Rev 6 scripts | `study/scripts/analysis/h_measures.py` | 11,015 |
| Rev 6 scripts | `study/scripts/analysis/make_manifest.py` | 5,621 |
| Rev 6 scripts | `study/scripts/analysis/make_manifest_amendment.py` | 5,258 |
| Rev 6 scripts | `study/scripts/analysis/make_manifest_g21_g22.py` | 2,366 |
| Rev 6 scripts | `study/scripts/analysis/make_manifest_g23.py` | 1,377 |
| Rev 6 scripts | `study/scripts/analysis/peptideatlas_feasibility.py` | 3,493 |
| Rev 6 scripts | `study/scripts/analysis/populations.py` | 2,010 |
| Rev 6 scripts | `study/scripts/analysis/propensity.py` | 8,284 |
| Rev 6 scripts | `study/scripts/analysis/provenance.py` | 3,978 |
| Rev 6 scripts | `study/scripts/analysis/run_amend.py` | 7,042 |
| Rev 6 scripts | `study/scripts/analysis/run_cem_grid.py` | 6,265 |
| Rev 6 scripts | `study/scripts/analysis/run_d2b.py` | 5,442 |
| Rev 6 scripts | `study/scripts/analysis/run_f.py` | 2,175 |
| Rev 6 scripts | `study/scripts/analysis/run_g14.py` | 3,439 |
| Rev 6 scripts | `study/scripts/analysis/run_g21.py` | 9,459 |
| Rev 6 scripts | `study/scripts/analysis/run_g22.py` | 793 |
| Rev 6 scripts | `study/scripts/analysis/run_g23.py` | 1,008 |
| Rev 6 scripts | `study/scripts/analysis/run_g5.py` | 4,865 |
| Rev 6 scripts | `study/scripts/analysis/run_j.py` | 9,527 |
| Rev 6 scripts | `study/scripts/analysis/run_propensity.py` | 4,516 |
| Rev 6 scripts | `study/scripts/analysis/scan_si.py` | 2,900 |
| Rev 6 scripts | `study/scripts/analysis/validate_engine.py` | 2,978 |
| Rev 6 scripts | `study/scripts/analysis/validate_engine_modes.py` | 1,953 |
| Rev 6 scripts | `study/scripts/analysis/validate_engine_v50.py` | 1,152 |
| Rev 6 scripts | `study/scripts/analysis/validate_propensity.py` | 1,603 |
| Rev 6 scripts | `study/scripts/analysis/verify_downloads.py` | 3,271 |
| Reference copies of the checkers at the freeze | `study/scripts/analysis/checkers_at_freeze/README.md` | 1,235 |
| Reference copies of the checkers at the freeze | `study/scripts/analysis/checkers_at_freeze/check_bindings.py` | 11,615 |
| Reference copies of the checkers at the freeze | `study/scripts/analysis/checkers_at_freeze/check_literals.py` | 5,625 |
| Reference copies of the checkers at the freeze | `study/scripts/analysis/checkers_at_freeze/check_numbers_sourced.py` | 20,360 |
| Plan, deviation log and wording | `study/analysis/ANALYSIS_PLAN.md` | 79,028 |
| Plan, deviation log and wording | `study/analysis/table_wording.json` | 759 |
| Plan, deviation log and wording | `study/analysis/submission_limits.yml` | 677 |
| Manifests | `study/analysis/grid_manifest.csv` | 90,350 |
| Manifests | `study/analysis/results/manifest_status.csv` | 136,002 |
| Manifests | `study/analysis/results/downloads_manifest.json` | 11,566 |
| Stored results | `study/analysis/results/amend_results.json` | 88,449 |
| Stored results | `study/analysis/results/b1.log` | 1,224 |
| Stored results | `study/analysis/results/b1_reconciliation.json` | 7,186 |
| Stored results | `study/analysis/results/b4.log` | 3,649 |
| Stored results | `study/analysis/results/b4_twoway.json` | 4,041 |
| Stored results | `study/analysis/results/b5.log` | 987 |
| Stored results | `study/analysis/results/b5_index_effect_sizes.json` | 4,376 |
| Stored results | `study/analysis/results/b_cluster_rule_info.json` | 4,779 |
| Stored results | `study/analysis/results/b_largest_clusters.json` | 1,000 |
| Stored results | `study/analysis/results/builder_numeric_literals.md` | 41,691 |
| Stored results | `study/analysis/results/cem_results.json` | 59,690 |
| Stored results | `study/analysis/results/claims_bindings.json` | 33,352 |
| Stored results | `study/analysis/results/d1_g2p.log` | 371 |
| Stored results | `study/analysis/results/d1_gene2pubmed_summary.json` | 234 |
| Stored results | `study/analysis/results/d2.log` | 2,911 |
| Stored results | `study/analysis/results/d2_index_variants_summary.json` | 2,872 |
| Stored results | `study/analysis/results/d2_index_variants_summary_proxy5_superseded.json` | 2,867 |
| Stored results | `study/analysis/results/d2b/G18-0407.json` | 1,388 |
| Stored results | `study/analysis/results/d2b/G18-0408.json` | 1,395 |
| Stored results | `study/analysis/results/d2b/G18-0409.json` | 1,411 |
| Stored results | `study/analysis/results/d2b/G18-0410.json` | 1,407 |
| Stored results | `study/analysis/results/d2b/G18-0411.json` | 1,394 |
| Stored results | `study/analysis/results/d2b/G18-0412.json` | 1,388 |
| Stored results | `study/analysis/results/d2b/G18-0413.json` | 1,414 |
| Stored results | `study/analysis/results/d2b/G18-0414.json` | 1,400 |
| Stored results | `study/analysis/results/d2b/G18-0415.json` | 1,390 |
| Stored results | `study/analysis/results/d2b/G18-0416.json` | 1,394 |
| Stored results | `study/analysis/results/d2b/G18-0417.json` | 1,368 |
| Stored results | `study/analysis/results/d2b_a.log` | 731 |
| Stored results | `study/analysis/results/d2b_b.log` | 246 |
| Stored results | `study/analysis/results/d2b_c.log` | 253 |
| Stored results | `study/analysis/results/d3_d4_d5_index_checks.json` | 14,966 |
| Stored results | `study/analysis/results/derived_figure_values.json` | 34,895 |
| Stored results | `study/analysis/results/downloads.log` | 5,329 |
| Stored results | `study/analysis/results/downloads2.log` | 2,423 |
| Stored results | `study/analysis/results/downloads3.log` | 0 |
| Stored results | `study/analysis/results/e_abundance_coverage.json` | 1,674 |
| Stored results | `study/analysis/results/f_G6_S1-S4.log` | 1,723 |
| Stored results | `study/analysis/results/f_G7_a.log` | 3,136 |
| Stored results | `study/analysis/results/f_G7_b.log` | 3,140 |
| Stored results | `study/analysis/results/f_G8.log` | 3,425 |
| Stored results | `study/analysis/results/f_G8_S1-S2.log` | 3,356 |
| Stored results | `study/analysis/results/f_G8_S3-S4.log` | 3,357 |
| Stored results | `study/analysis/results/f_runs/G6_S1-S2-S3-S4.json` | 45,432 |
| Stored results | `study/analysis/results/f_runs/G7_S1-S2.json` | 85,430 |
| Stored results | `study/analysis/results/f_runs/G7_S3-S4.json` | 85,992 |
| Stored results | `study/analysis/results/f_runs/G8_S1-S2-S3-S4.json` | 96,885 |
| Stored results | `study/analysis/results/f_tiers_summary.json` | 1,611 |
| Stored results | `study/analysis/results/figure_binding_check.md` | 2,569 |
| Stored results | `study/analysis/results/figure_bindings.json` | 102,102 |
| Stored results | `study/analysis/results/figure_captions.json` | 8,840 |
| Stored results | `study/analysis/results/figure_counts.json` | 34,201 |
| Stored results | `study/analysis/results/g11.log` | 3,545 |
| Stored results | `study/analysis/results/g11_predictors.json` | 25,675 |
| Stored results | `study/analysis/results/g14_cem.json` | 6,817 |
| Stored results | `study/analysis/results/g16.log` | 1,625 |
| Stored results | `study/analysis/results/g16_coverage.json` | 737 |
| Stored results | `study/analysis/results/g19.log` | 4,537 |
| Stored results | `study/analysis/results/g20.log` | 1,413 |
| Stored results | `study/analysis/results/g20_diagnostics.json` | 22,736 |
| Stored results | `study/analysis/results/g21.log` | 4,088 |
| Stored results | `study/analysis/results/g21_cond_logit_check.json` | 1,257 |
| Stored results | `study/analysis/results/g21_results.json` | 30,816 |
| Stored results | `study/analysis/results/g21_sequence_coverage.json` | 3,338 |
| Stored results | `study/analysis/results/g22_results.json` | 1,927 |
| Stored results | `study/analysis/results/g23_kept_nonevidence_pmids.json` | 1,144 |
| Stored results | `study/analysis/results/g23_nonprimary_pmids.json` | 4,513 |
| Stored results | `study/analysis/results/g23_pubmed_metadata.json` | 261,956 |
| Stored results | `study/analysis/results/g23_results.json` | 17,669 |
| Stored results | `study/analysis/results/g23_tiers_summary.json` | 8,796 |
| Stored results | `study/analysis/results/g5.log` | 2,787 |
| Stored results | `study/analysis/results/g5_abundance.json` | 34,258 |
| Stored results | `study/analysis/results/g5_diagnostics.json` | 2,453 |
| Stored results | `study/analysis/results/g9_intersections.json` | 642 |
| Stored results | `study/analysis/results/g9_measures.json` | 20,381 |
| Stored results | `study/analysis/results/g9_uniprot.log` | 200 |
| Stored results | `study/analysis/results/gene_unit_summary.json` | 5,753 |
| Stored results | `study/analysis/results/number_binding_check.md` | 34,417 |
| Stored results | `study/analysis/results/number_source_check.md` | 14,461 |
| Stored results | `study/analysis/results/panel_bindings.json` | 14,071 |
| Stored results | `study/analysis/results/peptideatlas_606_feasibility.json` | 1,452 |
| Stored results | `study/analysis/results/peptideatlas_606_inspection.json` | 11,795 |
| Stored results | `study/analysis/results/peptideatlas_606_web_check.json` | 1,414 |
| Stored results | `study/analysis/results/prop/G1-0003+posthoc_STY.json` | 1,124 |
| Stored results | `study/analysis/results/prop/G1-0003.json` | 539 |
| Stored results | `study/analysis/results/prop/G1-0006.json` | 1,027 |
| Stored results | `study/analysis/results/prop/G1-0009+posthoc_STY.json` | 1,126 |
| Stored results | `study/analysis/results/prop/G1-0009.json` | 531 |
| Stored results | `study/analysis/results/prop/G1-0012+posthoc_STY.json` | 1,116 |
| Stored results | `study/analysis/results/prop/G1-0012.json` | 1,028 |
| Stored results | `study/analysis/results/prop/G1-0015+posthoc_STY.json` | 1,115 |
| Stored results | `study/analysis/results/prop/G1-0015.json` | 1,026 |
| Stored results | `study/analysis/results/prop/G14-0337.json` | 1,058 |
| Stored results | `study/analysis/results/prop/G14-0340.json` | 1,057 |
| Stored results | `study/analysis/results/prop/G14-0343.json` | 1,066 |
| Stored results | `study/analysis/results/prop/G14-0346.json` | 1,058 |
| Stored results | `study/analysis/results/prop/G14-0349.json` | 560 |
| Stored results | `study/analysis/results/prop/G15-0352+posthoc_STY.json` | 1,101 |
| Stored results | `study/analysis/results/prop/G15-0352.json` | 518 |
| Stored results | `study/analysis/results/prop/G2-0018.json` | 1,026 |
| Stored results | `study/analysis/results/prop/G2-0021+posthoc_STY.json` | 1,129 |
| Stored results | `study/analysis/results/prop/G2-0021.json` | 533 |
| Stored results | `study/analysis/results/prop/G2-0024+posthoc_STY.json` | 1,131 |
| Stored results | `study/analysis/results/prop/G2-0024.json` | 1,038 |
| Stored results | `study/analysis/results/prop/G2-0027+posthoc_STY.json` | 1,119 |
| Stored results | `study/analysis/results/prop/G2-0027.json` | 532 |
| Stored results | `study/analysis/results/prop/G3-0030.json` | 1,028 |
| Stored results | `study/analysis/results/prop_ids.txt` | 161 |
| Stored results | `study/analysis/results/prop_ids2.txt` | 193 |
| Stored results | `study/analysis/results/propensity_order_dependence.json` | 1,064 |
| Stored results | `study/analysis/results/run_cem_g1_3_12_13_15.log` | 4,116 |
| Stored results | `study/analysis/results/run_cem_rerun_after_fix.log` | 6,074 |
| Stored results | `study/analysis/results/run_propensity.log` | 20,396 |
| Stored results | `study/analysis/results/run_propensity2.log` | 1,894 |
| Stored results | `study/analysis/results/run_propensity3.log` | 1,606 |
| Stored results | `study/analysis/results/si_banned_scan.md` | 404 |
| Stored results | `study/analysis/results/table_number_binding.csv` | 356,406 |
| Stored results | `study/analysis/results/table_number_check.md` | 2,484 |
| Stored results | `study/analysis/results/validate_engine.json` | 2,432 |
| Stored results | `study/analysis/results/validate_engine_modes.json` | 746 |
| Stored results | `study/analysis/results/validate_engine_v50.json` | 895 |
| Stored results | `study/analysis/results/validate_propensity.json` | 1,275 |
| Claims and results files | `study/analysis/CLAIMS_VS_EVIDENCE.md` | 29,432 |
| Claims and results files | `study/analysis/RESULTS_REPORT.md` | 113,404 |
| Figure and table code | `study/scripts/figures/build_all.py` | 7,252 |
| Figure and table code | `study/scripts/figures/build_tables.py` | 33,707 |
| Figure and table code | `study/scripts/figures/fig1_2.py` | 12,055 |
| Figure and table code | `study/scripts/figures/fig3.py` | 7,422 |
| Figure and table code | `study/scripts/figures/fig4_5.py` | 16,415 |
| Figure and table code | `study/scripts/figures/fig_supp.py` | 36,171 |
| Figure and table code | `study/scripts/figures/lib.py` | 12,116 |
| Figure and table code | `study/scripts/figures/make_contact_sheet.py` | 2,330 |
| Figure and table code | `study/scripts/figures/render_tables_pdf.py` | 4,690 |
| Composite figures, all formats | `study/manuscript/figures/Figure1.pdf` | 51,404 |
| Composite figures, all formats | `study/manuscript/figures/Figure1.png` | 740,715 |
| Composite figures, all formats | `study/manuscript/figures/Figure1.svg` | 18,856 |
| Composite figures, all formats | `study/manuscript/figures/Figure1.tif` | 1,346,978 |
| Composite figures, all formats | `study/manuscript/figures/Figure2.pdf` | 57,999 |
| Composite figures, all formats | `study/manuscript/figures/Figure2.png` | 1,073,432 |
| Composite figures, all formats | `study/manuscript/figures/Figure2.svg` | 52,528 |
| Composite figures, all formats | `study/manuscript/figures/Figure2.tif` | 1,357,586 |
| Composite figures, all formats | `study/manuscript/figures/Figure3.pdf` | 62,619 |
| Composite figures, all formats | `study/manuscript/figures/Figure3.png` | 690,552 |
| Composite figures, all formats | `study/manuscript/figures/Figure3.svg` | 67,962 |
| Composite figures, all formats | `study/manuscript/figures/Figure3.tif` | 925,728 |
| Composite figures, all formats | `study/manuscript/figures/Figure4.pdf` | 70,216 |
| Composite figures, all formats | `study/manuscript/figures/Figure4.png` | 1,208,880 |
| Composite figures, all formats | `study/manuscript/figures/Figure4.svg` | 110,683 |
| Composite figures, all formats | `study/manuscript/figures/Figure4.tif` | 1,379,034 |
| Composite figures, all formats | `study/manuscript/figures/Figure5.pdf` | 56,118 |
| Composite figures, all formats | `study/manuscript/figures/Figure5.png` | 849,407 |
| Composite figures, all formats | `study/manuscript/figures/Figure5.svg` | 37,877 |
| Composite figures, all formats | `study/manuscript/figures/Figure5.tif` | 1,182,024 |
| Composite figures, all formats | `study/manuscript/figures/FigureS1.pdf` | 63,711 |
| Composite figures, all formats | `study/manuscript/figures/FigureS1.png` | 722,127 |
| Composite figures, all formats | `study/manuscript/figures/FigureS1.svg` | 95,757 |
| Composite figures, all formats | `study/manuscript/figures/FigureS1.tif` | 990,052 |
| Composite figures, all formats | `study/manuscript/figures/FigureS2.pdf` | 86,770 |
| Composite figures, all formats | `study/manuscript/figures/FigureS2.png` | 1,105,315 |
| Composite figures, all formats | `study/manuscript/figures/FigureS2.svg` | 334,542 |
| Composite figures, all formats | `study/manuscript/figures/FigureS2.tif` | 1,458,550 |
| Composite figures, all formats | `study/manuscript/figures/FigureS3.pdf` | 50,208 |
| Composite figures, all formats | `study/manuscript/figures/FigureS3.png` | 358,986 |
| Composite figures, all formats | `study/manuscript/figures/FigureS3.svg` | 21,792 |
| Composite figures, all formats | `study/manuscript/figures/FigureS3.tif` | 559,444 |
| Composite figures, all formats | `study/manuscript/figures/FigureS4.pdf` | 57,365 |
| Composite figures, all formats | `study/manuscript/figures/FigureS4.png` | 939,101 |
| Composite figures, all formats | `study/manuscript/figures/FigureS4.svg` | 54,411 |
| Composite figures, all formats | `study/manuscript/figures/FigureS4.tif` | 1,218,652 |
| Composite figures, all formats | `study/manuscript/figures/FigureS5.pdf` | 50,726 |
| Composite figures, all formats | `study/manuscript/figures/FigureS5.png` | 539,386 |
| Composite figures, all formats | `study/manuscript/figures/FigureS5.svg` | 40,123 |
| Composite figures, all formats | `study/manuscript/figures/FigureS5.tif` | 737,092 |
| Composite figures, all formats | `study/manuscript/figures/FigureS6.pdf` | 60,046 |
| Composite figures, all formats | `study/manuscript/figures/FigureS6.png` | 467,737 |
| Composite figures, all formats | `study/manuscript/figures/FigureS6.svg` | 54,579 |
| Composite figures, all formats | `study/manuscript/figures/FigureS6.tif` | 791,328 |
| Composite figures, all formats | `study/manuscript/figures/FigureS7.pdf` | 48,818 |
| Composite figures, all formats | `study/manuscript/figures/FigureS7.png` | 271,385 |
| Composite figures, all formats | `study/manuscript/figures/FigureS7.svg` | 17,081 |
| Composite figures, all formats | `study/manuscript/figures/FigureS7.tif` | 397,662 |
| Composite figures, all formats | `study/manuscript/figures/FigureS8.pdf` | 47,847 |
| Composite figures, all formats | `study/manuscript/figures/FigureS8.png` | 364,729 |
| Composite figures, all formats | `study/manuscript/figures/FigureS8.svg` | 16,464 |
| Composite figures, all formats | `study/manuscript/figures/FigureS8.tif` | 724,334 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig1.png` | 745,306 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig2A.png` | 141,243 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig2B.png` | 259,668 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig2C.png` | 147,640 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig2D.png` | 190,109 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig2E.png` | 118,538 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig3A.png` | 482,386 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig3B.png` | 493,292 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig3C.png` | 483,536 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig4A.png` | 442,426 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig4B.png` | 368,112 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig4C.png` | 328,095 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig5A.png` | 369,117 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/Fig5B.png` | 207,776 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS1A.png` | 173,840 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS1B.png` | 127,812 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS1C.png` | 112,749 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS1D.png` | 98,896 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS1E.png` | 100,001 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS2A.png` | 62,916 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS2B.png` | 85,867 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS2C.png` | 101,702 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS2D.png` | 125,384 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS2E.png` | 74,729 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS2F.png` | 114,887 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS2G.png` | 363,556 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS2H.png` | 82,823 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS3A.png` | 426,431 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS4A.png` | 128,128 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS4B.png` | 194,408 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS4C.png` | 160,053 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS4D.png` | 504,422 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS5A.png` | 89,762 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS5B.png` | 143,584 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS5C.png` | 148,499 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS5D.png` | 73,087 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS6A.png` | 398,678 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS6B.png` | 405,504 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS6C.png` | 409,903 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS7A.png` | 236,078 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/FigS8A.png` | 307,535 |
| Panel PNGs and panel manifest | `study/manuscript/figures/panels/PANEL_MANIFEST.md` | 16,261 |
| `phosnet` package | `study/phosnet/.gitignore` | 189 |
| `phosnet` package | `study/phosnet/CHANGELOG.md` | 42 |
| `phosnet` package | `study/phosnet/LICENSE` | 1,072 |
| `phosnet` package | `study/phosnet/pyproject.toml` | 1,304 |
| `phosnet` package | `study/phosnet/README.md` | 4,241 |
| `phosnet` package | `study/phosnet/src/phosnet/__init__.py` | 762 |
| `phosnet` package | `study/phosnet/src/phosnet/evaluation.py` | 7,670 |
| `phosnet` package | `study/phosnet/src/phosnet/matching.py` | 11,170 |
| `phosnet` package | `study/phosnet/src/phosnet/sei.py` | 7,925 |
| `phosnet` package | `study/phosnet/tests/test_smoke.py` | 5,043 |
| Root documents, requirements, licence and citation files | `.gitignore` | 320 |
| Root documents, requirements, licence and citation files | `.zenodo.json` | 583 |
| Root documents, requirements, licence and citation files | `CITATION.cff` | 422 |
| Root documents, requirements, licence and citation files | `INPUTS.md` | 18,769 |
| Root documents, requirements, licence and citation files | `LICENSE` | 1,072 |
| Root documents, requirements, licence and citation files | `PROVENANCE.md` | 6,072 |
| Root documents, requirements, licence and citation files | `README.md` | 3,564 |
| Root documents, requirements, licence and citation files | `RELEASE_SCOPE.md` | 39,424 |
| Root documents, requirements, licence and citation files | `REPRODUCE.md` | 33,927 |
| Root documents, requirements, licence and citation files | `TITLE.txt` | 109 |
| Root documents, requirements, licence and citation files | `requirements.txt` | 901 |
| Root documents, requirements, licence and citation files | `scripts/set_title.py` | 2,075 |
| Root documents, requirements, licence and citation files | `study/scripts/release/build_release_scope.py` | 14,023 |
| Root documents, requirements, licence and citation files | `study/scripts/release/check_release_scope.py` | 9,432 |
| Tables written by the table code | `study/manuscript/tables/Table1_data_sources.csv` | 2,610 |
| Tables written by the table code | `study/manuscript/tables/Table1_data_sources.md` | 8,135 |
| Tables written by the table code | `study/manuscript/tables/Table1_data_sources_notes.md` | 274 |
| Tables written by the table code | `study/manuscript/tables/Table2_main_results.csv` | 866 |
| Tables written by the table code | `study/manuscript/tables/Table2_main_results.md` | 2,913 |
| Tables written by the table code | `study/manuscript/tables/Table2_main_results_notes.md` | 1,125 |
| Tables written by the table code | `study/manuscript/tables/Table3_predictors.csv` | 1,051 |
| Tables written by the table code | `study/manuscript/tables/Table3_predictors.md` | 3,082 |
| Tables written by the table code | `study/manuscript/tables/Table3_predictors_notes.md` | 548 |
| Tables written by the table code | `study/manuscript/tables/TableS1_contrast_grid.csv` | 124,115 |
| Tables written by the table code | `study/manuscript/tables/TableS1_contrast_grid.md` | 552,697 |
| Tables written by the table code | `study/manuscript/tables/TableS1_key.csv` | 5,846 |
| Tables written by the table code | `study/manuscript/tables/TableS1_key.md` | 28,479 |
| Tables written by the table code | `study/manuscript/tables/TableS2_main_results_P-human.csv` | 774 |
| Tables written by the table code | `study/manuscript/tables/TableS2_main_results_P-human.md` | 2,710 |
| Tables written by the table code | `study/manuscript/tables/TableS2_main_results_P-human_notes.md` | 1,211 |
| Tables written by the table code | `study/manuscript/tables/TableS3_non_primary_pmids.csv` | 5,788 |
| Tables written by the table code | `study/manuscript/tables/TableS3_non_primary_pmids.md` | 12,421 |
| Tables written by the table code | `study/manuscript/tables/TableS3_non_primary_pmids_notes.md` | 466 |
| Tables written by the table code | `study/manuscript/tables/TableS4_uniprot_experimental_vs_inferred.csv` | 1,084 |
| Tables written by the table code | `study/manuscript/tables/TableS4_uniprot_experimental_vs_inferred.md` | 2,549 |
| Tables written by the table code | `study/manuscript/tables/TableS4_uniprot_experimental_vs_inferred_notes.md` | 413 |
| Tables written by the table code | `study/manuscript/tables/TableS5_figure4_sample_sizes.csv` | 2,961 |
| Tables written by the table code | `study/manuscript/tables/TableS5_figure4_sample_sizes.md` | 4,784 |
| Tables written by the table code | `study/manuscript/tables/TableS5_figure4_sample_sizes_notes.md` | 241 |
| Supporting Information package | `study/manuscript/supporting_information/SI_FILE_LIST.md` | 1,344 |
| Supporting Information package | `study/manuscript/supporting_information/SupportingInformation.pdf` | 553,748 |
| Supporting Information package | `study/manuscript/supporting_information/TableS1_contrast_grid.xlsx` | 68,296 |
| Files the number checks read | `study/analysis/INTERIM_REPORT.md` | 35,973 |
| Files the number checks read | `study/analysis/RECONCILIATION.md` | 15,971 |
| Files the number checks read | `study/analysis/sources/SOURCES.md` | 21,337 |
| Files the number checks read | `study/analysis/nonprimary_pmid_list.md` | 8,989 |
| Files the number checks read | `study/audit_reports/phaseH1_step3_predictor_audit.md` | 61,290 |
| Files the number checks read | `study/earlier_results/consequence_experiment_results.json` | 1,234 |
| Files the number checks read | `study/earlier_results/dataset_counts.json` | 6,068 |
| Files the number checks read | `study/earlier_results/matching_results.json` | 11,819 |
| Files the number checks read | `study/earlier_results/predictor_metrics.json` | 9,124 |
| Files the number checks read | `study/earlier_results/tool_versions.json` | 2,734 |
| Modules a listed script imports | `study/scripts/earlier_figures/figlib.py` | 3,827 |
| Modules a listed script imports | `study/scripts/earlier_predictor_audit/step3/score_metrics.py` | 17,121 |
| Modules a listed script imports | `study/scripts/text_checks/banned_words.txt` | 376 |
| Modules a listed script imports | `study/scripts/text_checks/scan_banned.py` | 3,196 |
