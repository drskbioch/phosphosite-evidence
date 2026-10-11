# Provenance

- The analysis was frozen at the tag `rev6-final` (commit `3521d5c`) in the private working repository.
- Changes made after the freeze are listed in section 15 of `study/analysis/ANALYSIS_PLAN.md`, in the entries headed "Post-freeze extension". The result files that hold estimates are byte-identical between the tag and this folder. Among the stored results, the files that differ from the tag are `claims_bindings.json`, `downloads_manifest.json` and the reports `number_source_check.md`, `number_binding_check.md` and `builder_numeric_literals.md`, together with files added after the tag (new bindings, check reports, run logs); those entries describe each.
- This folder was built from commit `3a481f81662db3863a02633d4a855eb038e059d5` of the private working repository: the files listed in `RELEASE_SCOPE.md`, taken from that commit, with the folder and file names changed as listed under Renamed paths, and with these changes and additions: `CITATION.cff` and `.zenodo.json` (title, version 1.0.0, author), a new `README.md`, this file, `TITLE.txt`, `scripts/set_title.py` and `.gitignore`.
- The private working repository, with its history, is not part of the archive.

## Redactions

The private working repository is unchanged. The files of this folder differ from it only in the ways listed here, in the renamed paths below and in the additions named above.

| File | What was replaced | Why |
|---|---|---|
| `study/analysis/results/run_propensity.log` | 16 occurrences of the absolute path of the project folder on the build machine, replaced by `<project-root>` | the path appears in Python tracebacks and names the machine's folders; nothing else in the log changed |
| `study/analysis/results/f_G8_S1-S2.log` | 4 occurrences, same replacement | same |
| `study/analysis/results/f_G8_S3-S4.log` | 4 occurrences, same replacement | same |
| `study/analysis/results/d2b_a.log` | 3 occurrences, same replacement | same |
| `study/scripts/earlier_predictor_audit/step3/score_metrics.py` | the three default paths (evaluation sample, output folder, PhosF3C result file) are strings relative to the project root | the defaults were absolute paths of the build machine; `run_j.py` sets the output folder and the PhosF3C path itself (`REPRODUCE.md`, section 7) |
| `study/phosnet/pyproject.toml` | the `[project.urls]` table (three links to the private working repository) is removed | the links point to a repository that is not part of the archive |
| `REPRODUCE.md`, `INPUTS.md` | wording that named a Zenodo snapshot now says the archive is deposited from a GitHub release; `REPRODUCE.md` section 7 describes the `score_metrics.py` change | the archive is deposited from a GitHub release |
| `study/scripts/release/build_release_scope.py`, `RELEASE_SCOPE.md` | the builder lists the five files added in this folder and the header names the GitHub release; `RELEASE_SCOPE.md` is written again by the builder and lists all files with their sizes; the list of excluded paths names `study/_archive/` once instead of the predictor checkouts inside it | the file list must name every file of the folder |

No estimate, figure or table was changed. The names of the folders and files below and the path strings inside the stored files were changed; every number in every stored file is identical.

## Renamed paths

The folders and files of the archive carry the names below; the private working repository (and the tag `rev6-final`) use the names on the left. Every path written in the documents, the scripts and the stored files of this folder (for example in `claims_bindings.json` and `figure_bindings.json`) uses the new names.

| Name in the working repository | Name in the archive | Holds |
|---|---|---|
| `current/` | `study/` | the top folder of the working repository that holds the study |
| `current/rev6/` | `study/analysis/` | analysis plan and deviation log, manifests, claims table, results report, `results/`, `sources/` |
| `current/scripts/rev6/` | `study/scripts/analysis/` | analysis, report and check scripts |
| `current/scripts/rev6/_frozen_at_rev6-final/` | `study/scripts/analysis/checkers_at_freeze/` | the three number checkers as they were at the freeze, for reference |
| `current/scripts/figures_rev6/` | `study/scripts/figures/` | figure and table code |
| `current/scripts/figures/` | `study/scripts/earlier_figures/` | figure style module of the earlier stages |
| `current/scripts/high_tier/` | `study/scripts/earlier_predictor_audit/` | score loaders of the earlier predictor audit |
| `current/scripts/manuscript_checks/` | `study/scripts/text_checks/` | banned-word scan |
| `current/manuscript/figures_rev6/` | `study/manuscript/figures/` | composite figures and panels |
| `current/manuscript/supplementary/` | `study/manuscript/supporting_information/` | Supporting Information package |
| `current/results/` | `study/earlier_results/` | the five result files of the earlier stages that the scripts read |
| `current/phosnet/, current/audit_reports/, current/manuscript/tables/, current/scripts/release/` | `study/phosnet/, study/audit_reports/, study/manuscript/tables/, study/scripts/release/` | only the top folder changed |

| File name in the working repository | File name in the archive |
|---|---|
| `ANALYSIS_PLAN_REV6.md` | `ANALYSIS_PLAN.md` |
| `RESULTS_REV6.md` | `RESULTS_REPORT.md` |
| `build_claims_rev6.py` | `build_claims.py` |
| `build_tables_rev6.py` | `build_tables.py` |
| `count_descriptives_rev6.py` | `count_descriptives.py` |
| `f_G6_S1,S2,S3,S4.log` | `f_G6_S1-S4.log` |
| `f_G8_S1,S2.log` | `f_G8_S1-S2.log` |
| `f_G8_S3,S4.log` | `f_G8_S3-S4.log` |
| `g23_nonprimary_list.md` | `nonprimary_pmid_list.md` |

Not renamed: the script and result files whose names carry a task, block or contrast identifier (for example `run_g21.py`, `g23_results.json`, `prop/G1-0003+posthoc_STY.json`), because the identifier is the key used in the plan, the manifest and the stored results; the names `phosnet` and `rev6-final`; and the names of the identifiers `REV6` and `R6` inside the scripts.
