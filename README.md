# <!--title-->Counting rules and study effort in phosphosite evidence tiers: an audit of the signals used to validate them<!--/title-->

This archive holds the analysis code, the stored results, the figures and tables, the Supporting Information and the number checks for the study "<!--title-->Counting rules and study effort in phosphosite evidence tiers: an audit of the signals used to validate them<!--/title-->". The code compares the evidence tiers GOLD, SILVER and BRONZE of phosphosite records in the EPSD 2.0 resource (human, mouse and rat), measures how the tier assignment depends on the counting rule and on the amount of study effort behind each protein, and tests four validation signals (domain overlap, conservation, disorder and sequence motif) and three phosphorylation-site predictors against these quantities.

## What is in the archive

- `study/scripts/analysis/`: the analysis scripts, from the inputs to the stored results, the report builders and the number checks. `REPRODUCE.md` gives the run order.
- `study/analysis/`: the analysis plan with its deviation log (`ANALYSIS_PLAN.md`), the grid of planned contrasts (`grid_manifest.csv`), the stored results (`results/`), the claims table (`CLAIMS_VS_EVIDENCE.md`) and the results report (`RESULTS_REPORT.md`).
- `study/scripts/figures/`, `study/manuscript/figures/`, `study/manuscript/tables/` and `study/manuscript/supporting_information/`: the code and the files of the figures, the panels, the tables and the Supporting Information.
- `study/phosnet/`: the `phosnet` Python package. phosnet 0.1.0 implements earlier estimators and is not used for the results; the results were produced by the scripts in `study/scripts/analysis/`.
- `REPRODUCE.md`, `INPUTS.md`, `RELEASE_SCOPE.md` (the list of the files taken from the working repository), `PROVENANCE.md`, `requirements.txt`, `LICENSE`, `CITATION.cff`, `.zenodo.json`, `TITLE.txt` and `scripts/set_title.py` (writes the title into this file, `CITATION.cff` and `.zenodo.json`).

## What is not included, and why

- **EPSD site-level data.** The EPSD files and the tables derived from them (one row per site, protein or PMID) are not redistributed: the provider does not state a licence, and permission has been requested.
- **PhosphoSitePlus, Phospho.ELM and dbPTM.** Site-level data from these databases are not redistributed. They enter the analysis only as evidence flags and publication lists inside the EPSD records.
- **Predictor code, weights and site-level scores.** The three predictors (DeepPhos, MusiteDeep, PhosF3C) are obtained from their own repositories (`INPUTS.md`, section 5); their scores for the evaluation sample are EPSD-derived and are not redistributed.
- **Public downloads** (NCBI, UniProt, PaxDb, Pharos, GOA, PeptideAtlas) are not included; `INPUTS.md` gives the source, version, retrieval date and SHA-256 of each.

## How to reproduce

`REPRODUCE.md` has two levels. The reports, the figures, the tables, the Supporting Information and the number checks can be run from the stored results alone (under two minutes, measured). Recomputing the stored results needs the inputs in `INPUTS.md`, several of which are not in the archive, and about 15 hours of computation. Install the packages with `pip install -r requirements.txt`.

## How to cite

See `CITATION.cff`. DOI: TO BE SET AFTER RELEASE.

## Author and licence

Dr Shahid Ullah, drskbioch@gmail.com, https://www.habdsk.org. The code and the files of this archive are released under the MIT licence (`LICENSE`): Copyright (c) 2026 Dr Shahid Ullah.
