# Supporting Information: file list

| File name | File type | Description |
|---|---|---|
| SupportingInformation.pdf | PDF, 14 pages (A4) | Contents page; Figures S1 to S8 (one page each, vector graphics, notes below each figure); Tables S2 to S5 with their notes. |
| TableS1_contrast_grid.xlsx | Excel workbook, 2 sheets | Sheet "Table S1": the 534 contrasts of the planned grid (identifier, block, task, signal, population, design, index, tier rule, status, estimate, 95% interval and its resampling unit, second interval where the comparison was run with both units, attenuation, result status and note). Sheet "Key": the meaning of the columns and codes, with the number of contrasts for each. |
| SI_FILE_LIST.md | Markdown | This list. |

Each figure is also available as a stand-alone file in PDF, SVG, PNG (600 dpi) and TIFF (600 dpi) in `study/manuscript/figures/` (FigureS1 to FigureS8), and each panel as a PNG file in `study/manuscript/figures/panels/`; these are not part of the package.
Tables S2 to S5 are also available as CSV files in `study/manuscript/tables/`; the numbers in them are bound to result files by `study/scripts/analysis/check_table_numbers.py` (report: `study/analysis/results/table_number_check.md`).
The package is built by `study/scripts/analysis/build_si.py` from the figure and table files named above.
