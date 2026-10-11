# Reconciliation of Table 1 and the text (Rev 5), Task B1

Written by `scripts/analysis/b1_reconcile.py`. Every value below is recomputed from the site table (`mammalian_sites.parquet`, `human_sites_source_set_a.parquet`, `sei_protein.parquet`) and compared with `study/earlier_results/matching_results.json`, `dataset_counts.json`, Table 1 and the manuscript text. No manuscript file was edited.

## 1. Unmatched estimates: the text and Table 1 use different record sets

Two record sets exist in `matching_results.json`: **all records** (`unmatched_all_sites`: every GOLD and BRONZE record with the outcome defined, whether or not its protein has a primary index) and **records with an index** (`unmatched_mapped_sites`: the same records restricted to proteins with a primary index, which is the set the matched designs use). Table 1 states that its unmatched column uses the same records as the matched columns, and it does (records with an index). The text does not do this consistently.

| Comparison | Tier source; species | Text quotes | Table 1 quotes | Recomputed, all records | Recomputed, records with an index | Text uses | Table 1 uses |
|---|---|---|---|---|---|---|---|
| Domain overlap, human (OR) | source set A (human file); human | 2.494 | 2.494 | 2.734 (n GOLD 157,785, BRONZE 607,096) | 2.494 (n GOLD 153,738, BRONZE 483,427) | records with an index | records with an index |
| Domain overlap, annotated proteins only, human (OR) | source set A (human file); human | 0.586 | 0.586 | 0.570 (n GOLD 147,185, BRONZE 284,373) | 0.586 (n GOLD 143,906, BRONZE 242,308) | records with an index | records with an index |
| Conservation (difference) | source set B (all three species); human, mouse, rat | 0.3760 | 0.3986 | 0.3760 (n GOLD 95,973, BRONZE 40,414) | 0.3986 (n GOLD 93,308, BRONZE 37,591) | all records | records with an index |
| Disorder (difference) | source set B (all three species); human, mouse, rat | 0.0310 | 0.0438 | 0.0310 (n GOLD 337,716, BRONZE 257,209) | 0.0438 (n GOLD 320,780, BRONZE 215,745) | all records | records with an index |
| Kinase motif, sequence-based (difference) | source set B (all three species); human, mouse, rat | 0.0433 | 0.0617 | 0.0433 (n GOLD 342,439, BRONZE 674,782) | 0.0617 (n GOLD 324,907, BRONZE 478,837) | all records | records with an index |
| Kinase motif, database lookup (difference) | source set B (all three species); human, mouse, rat | 0.7780 | 0.7934 | 0.7780 (n GOLD 342,439, BRONZE 674,827) | 0.7934 (n GOLD 324,907, BRONZE 478,880) | all records | records with an index |

All twelve recomputed values equal the result-file values to within 1e-9. Findings:

- The text and Table 1 give **different unmatched values for 4 of 6 rows** (Conservation, Disorder, Kinase motif, sequence-based, Kinase motif, database lookup). In each of these the text quotes the all-records value and Table 1 the records-with-an-index value. The text's domain-overlap sentences (full population and annotated-only) quote the records-with-an-index values, so the text mixes both sets. For disorder the two values are 0.0310 (text) and 0.0438 (Table 1), the mismatch noted in the plan.
- The matched columns are unaffected (they use records with an index in both places). The records excluded for having no primary index, 127,716 (domain overlap), 5,488 (conservation), 58,400 (disorder) and 213,479 (kinase motif), reproduce exactly (`n_unmapped_recomputed` equals the result file for every row in `results/b1_reconciliation.json`).
- **Tier source differs between rows.** Domain overlap (both rows) used tiers computed from source set A (the human file) and is human only; the other four rows used set-B tiers over human, mouse and rat. "GOLD" is therefore not the same set of human records in the domain-overlap row and in the other rows. Rev 6 uses set-B tiers for every signal (decision D1); the set-A result is kept as block G15.
- The sequence-based kinase-motif row has 46 records with an undefined outcome (non-canonical residues). 43 of them are GOLD or BRONZE records of proteins with an index. Rev 5 counted them in `n_analysis` (803,787) and in the index deciles but never in an estimate; the Rev 6 population excludes them before the deciles are formed (803,744). The effect on the adjusted estimate is 6.8e-6 (0.063056 versus 0.063063); see `results/validate_engine.json`.

## 2. Ortholog coverage: 7.4% and 12.3% are different quantities

| Quantity | Count | Denominator | Share |
|---|---|---|---|
| Conservation indicator defined, GOLD or BRONZE, protein with an index (the 130,899 analysis records; **text, 7.4%**) | 130,899 | all 1,768,068 records | 7.40% |
| Conservation indicator defined, human records (**limitations, 12.3%**) | 131,687 | 1,071,725 human records | 12.29% |
| Conservation indicator defined, all records | 284,413 | all 1,768,068 records | 16.09% |
| Conservation indicator defined, GOLD or BRONZE | 136,387 | 1,017,266 GOLD or BRONZE records | 13.41% |
| Defined, mouse records | 103,465 | mouse records | 19.37% |
| Defined, rat records | 49,261 | rat records | 30.39% |

The text's 7.4% is the analysis set as a share of the whole dataset, so it mixes a tier restriction and an index restriction into a coverage figure. The limitations sentence's 12.3% is the human share only, and the sentence calls it "records" without saying human. Neither is the dataset-wide coverage of the indicator (16.09%). The next manuscript revision should state the denominator for each.

## 3. Study clusters: 1,113 and 1,070 are the same quantity on different records

- 1,113 is the number of distinct **primary study cluster** labels over all 1,768,068 records (a record's primary cluster is the cluster of its globally most frequent PMID; `cluster_split.parquet`). 1,070 is the same count over the 1,768,020 records of the temporal split, which excludes 48 records (the manuscript text attributes these to 46 non-canonical residues and 2 records with no linked PMID; that split is quoted, not re-derived here). 43 labels occur only on those records, and 1,113 minus 43 is 1,070. **Publication-year resolution plays no part.**
- Neither is the number of study clusters under the merge rule. Under the rule (candidate PMID at >= 884 sites, Jaccard >= 0.90) the 50,904 distinct PMIDs fall into **50,891 clusters**: 690 candidate PMIDs, of which 21 sit in 8 multi-PMID clusters (sizes 7, 2, 2, 2, 2, 2, 2, 2), and every other PMID is its own cluster. The sentence "The dataset contained 1,113 study clusters" therefore describes the number of distinct primary-cluster labels, not the number of clusters.

## 4. Study-cluster rule: how 884 arose, and how little it matters

- **Trace.** The code default of `build_sei_secondary_table(high_freq_min=...)` is 500 (`phosnet/sei.py:126`). The Rev 5 tiers were built with 884. `audit_reports/phase1_mammalian_ingest_report.md` (line 78) records that the human-only stage used a cutoff of 500 sites (500 / 1,071,725 = 0.047%), restated as 0.05% of the pooled count: round(0.0005 x 1,768,068) = 884. It is a proportional rule meant to keep an earlier rate. `phase0_epsd_gate_report.md` (line 216) gives the reason for the 500 cutoff itself as a judgement (a PMID below it "cannot plausibly be a mega-consortium companion record"); no data-derived justification for 500 is recorded in the repository. 884 was not derived from a data-quality criterion.
- **Sensitivity (reported as found).** The three thresholds differ in the number of candidate PMIDs and candidate clusters but **not in any record's study-cluster depth or tier**:

| Threshold (sites) | Candidate PMIDs | Candidate clusters | Candidates in multi-PMID clusters | Multi-PMID clusters | GOLD / SILVER / BRONZE (study-cluster rule) | GOLD / SILVER / BRONZE (depth-only rule) |
|---|---|---|---|---|---|---|
| 500 | 757 | 744 | 21 | 8 | 342,439 / 750,802 / 674,827 | 246,111 / 500,437 / 1,021,520 |
| 884 | 690 | 677 | 21 | 8 | 342,439 / 750,802 / 674,827 | 246,111 / 500,437 / 1,021,520 |
| 1000 | 674 | 661 | 21 | 8 | 342,439 / 750,802 / 674,827 | 246,111 / 500,437 / 1,021,520 |

Records whose depth differs from the 884 rule: 0 at 500 and 0 at 1,000. The same 21 PMIDs form the same 8 multi-PMID clusters at every threshold (the extra candidates between 500 and 884 sites are singletons). **Block G13 (cluster rule at 500 and 1,000) therefore reproduces G1 record for record and cannot show any dependence on the threshold.** It is kept in the manifest and reported as identical, as planned. This is a null result about the sensitivity of the rule, not evidence that the rule is correct.
- Reproduction check: depth under 884 equals the Rev 5 `study_cluster_count` for 100% of records and the tier equals the Rev 5 tier for 100% (`results/b_cluster_rule_info.json`).

## 5. Where the code departs from the text of the plan

- The seven consortium PMIDs are exactly the PMID set of the second-largest cluster (`CLUSTER_38359819`); their per-record depth arrays are identical. The planned 11th removal in block G7 (`drop_consortium_7_PMIDs`) therefore duplicates `drop_cluster_2`. Both rows stay in the manifest; the results are identical by construction and the decision count counts the removal once.
- Sites with no PMID left after a removal: 130,828 for the consortium removal, of which 13 carry a database source flag; 130,815 have no PMID and no source evidence under either reading (agreement score 0, or no PhosphoSitePlus/Phospho.ELM/dbPTM flag) and are the number dropped.

## 6. The Rev 5 propensity matcher versus its Methods description

*Written by hand, not by `b1_reconcile.py`, from these files: the source of `phosnet/src/phosnet/matching.py` (lines 158 to 225); the Methods sentence quoted below, from the earlier manuscript text; `results/propensity_order_dependence.json`; the not-fitted records in `results/prop/`; and `INTERIM_REPORT.md`.*

- **What the Rev 5 code does** (`phosnet/src/phosnet/matching.py`, `propensity_matched_comparison`, lines 158 to 225). The logit of the propensity score is fitted by statsmodels on all analysis records. For each GOLD record, visited in random order, the code looks at that record's **50 nearest BRONZE records on the logit** (`k = min(50, len(control_idx))`, k-d tree order, so ties on the logit are resolved by the order in which records appear in the table), takes the first one that lies within the caliper (0.2 x SD of the logit over all records) and has not been used, and moves on. A GOLD record whose 50 nearest controls are all used or lie outside the caliper stays unmatched even when a free control within the caliper exists further away. The covariates are discrete, so many records share one logit value and the 50-record window is often filled by records with an identical logit; the result then depends on row order (`results/propensity_order_dependence.json`).
- **What the Rev 5 Methods say.** "Propensity matching paired records 1:1 without replacement within a caliper of 0.2 standard deviations of the logit." The pairing is indeed 1:1 and without replacement and uses that caliper. The Methods do not describe the 50-nearest-controls window or the resolution of ties by record position, and they do not say that the fit fails when the residue-type term separates the groups perfectly (non-canonical residues occur only among BRONZE records; 5 of 11 planned D2 contrasts could not be fitted for that reason, `results/prop/`).
- **Consequence.** The manuscript Methods must be rewritten to describe what was run (1:1 without replacement, caliper, a window of the 50 nearest controls, ties by record order, no fit for contrasts that include non-canonical residues). Rev 6 adds the matcher the Methods sentence suggests as a labelled post-hoc design, D2b (block G18; `scripts/analysis/d2b.py`): nearest unused control with no window, ties between identical logits broken at random, serine, threonine and tyrosine records only, at least 20 order replicates per contrast. The Rev 5 matcher results are kept beside it in `INTERIM_REPORT.md` and in the results file. The manuscript is not edited in this cycle.


## 7. The set-B PhosphoSitePlus extract and PMID 25514926 (checked 2026-10-04)

*Written by hand, not by `b1_reconcile.py`, from these files: `results/g23_pubmed_metadata.json` (the PubMed record of PMID 25514926 and its 225,187 linked records) and `results/g23_tiers_summary.json` (the same record count); the source of `study/ingest/phase1_5_das_repair.py`, `epsd_parser.py`, `phase1_ingest.py` and `epsd_pds_study_clustering.py`. The other counts in this section (rows and sources of the workbook `raw_data/new data of psp (1).xlsx`; species, tier and flag counts of the 225,187 records; the tier changes) come from direct queries of that workbook and of `mammalian_sites.parquet` on 2026-10-04 and are not stored in a result file.*

**The extract.** `raw_data/new data of psp (1).xlsx` (37,262,190 bytes; workbook metadata: last modified 2020-06-09) has three sheets; `Sheet3` holds 655,976 rows without a header in six columns (UniProt accession, peptide window with `p` marking the phosphorylated residue, residue and position, species, PMID, `Source`). `Source` is `PhosphoSitePlus` in 371,633 rows, `Phospho.ELM` in 45,802 and blank in 238,541. It has no per-site columns equivalent to `LT_LIT`, `MS_LIT` or `MS_CST`, no version stamp and no licence text, so it cannot stand in for `Phosphorylation_site_dataset.gz`; no PMID-count surrogate was derived from it (decision of the project lead, 2026-10-04).

**What PMID 25514926 is.** PubMed (E-utilities esummary, 2026-10-04): "PhosphoSitePlus, 2014: mutations, PTMs and recalibrations.", Nucleic Acids Research, 2015 Jan (electronic 2014 Dec 16), volume 43, pages D512-20, doi 10.1093/nar/gku1267, Hornbeck PV and colleagues; publication types Journal Article and Research Support. It is a database paper, not a primary study. In the workbook it is the PMID of 383,429 of the 655,976 rows (58%).

**Did PMIDs from the workbook enter publication lists, study clusters or depth?**
- *By the code:* no. `study/ingest/phase1_5_das_repair.py` reads `Sheet3`, assigns column names including `PMID`, and then uses only the (accession, residue, position) triplets of the rows whose `Source` is `PhosphoSitePlus` to set the set-B `source_PSP` flag (lines 31 to 43 and 110); the `PMID` column is never used. The publication lists of the current dataset are built from the EPSD download's `reference_raw` field (`epsd_parser.py` line 83, `phase1_ingest.py` line 73) and the study clusters and `PDS_v2` from those lists (`epsd_pds_study_clustering.py`). Three legacy scripts (`legacy/scripts/phase1/`) also load the workbook; they belong to the earlier project and were not read for the current dataset.
- *By the data:* **the same PMID is nevertheless in the dataset's publication lists.** 225,187 of 1,768,068 records (12.7%; human 146,964, mouse 57,339, rat 20,884) have 25514926 in their publication list (`mammalian_sites.parquet`, `pmid_list_str`); 84.3% of them carry the `source_PSP` flag (11.4% of the other records); 163,303 are GOLD, 58,466 SILVER and 3,418 BRONZE; for 3,208 of them it is the only PMID. The route is the EPSD reference field, not the workbook (the code does not read the workbook's PMIDs). It is therefore counted as a publication in study-cluster depth for those records.
- *Size of the effect (descriptive, not an analysis result and not used anywhere):* if the PMID were not counted (depth minus 1, same database agreement score, assuming it forms its own study unit; whether it was merged into one of the seven two-PMID clusters was not checked), the tier of 6,149 records would change: 984 from GOLD to SILVER and 5,165 from SILVER to BRONZE.
- *Status:* this is a finding about the publication-depth input (a database paper counted as a study), not about the workbook. It is reported to the project lead (the stop rule of the instruction of 2026-10-04 applies) and no tier, estimate or result has been changed.
