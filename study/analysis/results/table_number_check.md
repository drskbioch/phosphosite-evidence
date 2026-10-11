# Check of the numbers in the tables

Tables checked: 9 files (Tables 1 to 3 and S1 to S5) and their notes files. Numbers examined: 2478. **Unbound: 0** (superseded: 0; other: 0).
Table S1 (contrast grid): 534 rows compared with `results/manifest_status.csv`; fields that differ: 0; second intervals without a paired contrast: 0.

Binding of a number, in this order: result value that the table builder formatted (file and key re-read and recomputed by `check_figure_bindings.py`); text cell verified against its source text; claims-table binding (`claims_bindings.json`); count derived from `manifest_status.csv` or `dataset_counts.json`; documented constant; definition or identifier written in the text. The list of every number with its binding is `table_number_binding.csv` in this folder.

## Numbers by table and kind of binding

| Table | claims-table binding | count of rows of the Table S1 grid (equal to manifest_status.csv) | definition written in the text | documented constant | documented definition | manifest_status.csv | paired contrast of manifest_status.csv | result value of the table builder | text cell verified against its source | total |
|---|---|---|---|---|---|---|---|---|---|---|
| Table1_data_sources | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 30 | 64 | 95 |
| Table1_data_sources_notes | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| Table2_main_results | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 57 | 0 | 60 |
| Table2_main_results_notes | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 3 | 0 | 4 |
| Table3_predictors | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 9 | 27 | 36 |
| TableS1_contrast_grid | 0 | 0 | 0 | 0 | 0 | 1668 | 112 | 0 | 0 | 1780 |
| TableS1_key | 4 | 71 | 9 | 0 | 7 | 0 | 0 | 0 | 0 | 91 |
| TableS2_main_results_P-human | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 48 | 0 | 51 |
| TableS2_main_results_P-human_notes | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 3 | 0 | 4 |
| TableS3_non_primary_pmids | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 34 | 19 | 53 |
| TableS3_non_primary_pmids_notes | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 |
| TableS4_uniprot_experimental_vs_inferred | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 72 | 0 | 73 |
| TableS4_uniprot_experimental_vs_inferred_notes | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| TableS5_figure4_sample_sizes | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 225 | 0 | 227 |
| TableS5_figure4_sample_sizes_notes | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |

## Unbound numbers

none

## Documented constants used

- `95`: interval level of the 95% intervals: analysis/ANALYSIS_PLAN.md, task J (line 205: 95% intervals shown)
