# Interim report, Rev 6 reanalysis

Branch `rev6-reanalysis`. Written by `scripts/analysis/build_interim.py` from the files in `study/analysis/results/`; every number below is read from a result file. No manuscript file, reference or legacy file was touched. The plan is `ANALYSIS_PLAN.md` (deviation log in its section 15). **Work stops here, as planned.** Tasks F (strata, leave-one-cluster-out, tier thresholds), H (established measures), I (within-protein) and J (predictors) have not been started; `results/manifest_status.csv` lists every manifest row with its status. PhosphoSitePlus is not in `data_private/` yet.

## Headline results

1. **Matching on the study-effort index (CEM) removes none of the decision-bearing validation signals.** All 7 of 7 decision-bearing CEM contrasts in G1 and G2 survive (protein-cluster interval excludes the null); 0 are largely explained (attenuation of 50% or more); the largest attenuation is 1.4%.
2. **Abundance does not explain them either** (G5): 30 of 30 adjusted contrasts across human, mouse and rat survive, 0 with attenuation of 50% or more. All 20 index-variant contrasts (G4) survive as well.
3. **Honest inference widens the intervals but, in the planned analysis, not the conclusions.** The protein-cluster bootstrap widens the Rev 5 intervals by 1.2 to 5.6 times with the same inference in 10 of 10 CEM rows (G14), and the two-way (protein and study cluster) sensitivity is wider again: 2 of 5 contrasts no longer exclude the null there (disorder, motif, sequence-based).
4. **Domain overlap depends strongly on how GOLD is defined**: matched odds ratio 5.912 with set-B tiers against 2.646 with set-A tiers (same index), and the crude odds ratio falls below 1 in the abundance subset (0.741) and in the annotated-protein subset (0.692) against 6.065 in G1 (section 10).
5. **Rev 5 reporting problems found** (`RECONCILIATION.md`): text and Table 1 quote different unmatched estimates for 4 of 6 rows; 7.4% and 12.3% are different quantities; 1,113 and 1,070 count the same thing on different records. The study-cluster threshold (500, 884, 1,000) changes no record's tier.
6. **Not clean:** propensity matching (D2) could not be fitted for 5 of 11 planned contrasts (separated non-canonical residues) and its matcher depends on record order; the GO check is ill-posed under the pre-specified coding, so the D5 rule's outcome (V20) and its post-hoc counterpart (V100, measures disagree) need your decision (section 9, section 11).

## 1. Gate checks that precede the contrasts

- **CEM engine against Rev 5** (`validate_engine.py`, original index, Rev 5 record sets): S1, S2, S3 and S5 reproduce the Rev 5 adjusted estimates to within 8.9e-16, and the crude estimates exactly. S4 reproduces exactly (difference 0.0) on the Rev 5 record set. The Rev 6 population for S4 drops the 43 records with an undefined outcome before forming the index deciles, as Rev 5 did not; this moves the S4 adjusted estimate by 6.8e-06 (0.063056 against 0.063063). Recorded in the deviation log.
- **V50 index** (`validate_engine_v50.py`): the engine equals phosnet's `matched_comparison` for all five signals (largest difference 1.1e-16). **Abundance modes** (`validate_engine_modes.py`): equal to phosnet with abundance as the covariate or as an extra stratum (largest difference 1.1e-16). The Rev 5 secondary-index CEM rows reproduce in G14 (below); the Rev 5 propensity rows agree up to the randomness of the matching order (the seed of the visiting order differs).
- **Propensity module** (`validate_propensity.py`): reproduces the Rev 5 matched estimate and pair count exactly for conservation (0.390794, 31,001 pairs) and disorder (0.046485, 81,575 pairs), and reproduces the Rev 5 failure for the kinase-motif lookup (LinAlgError: Singular matrix). The cell-level fit used inside resamples agrees with the record-level statsmodels fit to 1.7e-04 on serine, threonine and tyrosine records (non-canonical records are separated and their fitted values differ; they are never matched, see section 5).
- A bug found while running the abundance block (zero-weight cells above the top decile edge received an out-of-range stratum index) was fixed in `engine.py`. It cannot change an estimate (the cells have weight 0). All earlier CEM contrasts (46 contrasts, 184 compared fields) were rerun with the fixed engine and are identical.

## 2. Downloads

| Data | Source URL | Version | Retrieved | Result |
|---|---|---|---|---|
| Pharos/TCRD target list with TDL | https://pharos-api.ncats.io/graphql (GraphQL, `targets{targets(top,skip){uniprot sym tdl fam}}`) | pharos319 | 2026-10-01 | 20,412 of 20,412 targets, 20,412 distinct UniProt IDs; TDL counts {'Tbio': 12303, 'Tdark': 5501, 'Tchem': 1904, 'Tclin': 704} |
| GOA human annotations | https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz | !gaf-version: 2.2; !date-generated: 2026-07-28 16:23; !generated-by: UniProt (GO 2026-07-26) | 2026-10-01 | 11,020,779 bytes, equal to Content-Length; complete gzip read |
| PaxDb integrated whole-organism, human | https://pax-db.org/downloads/6.0/datasets/9606/9606-WHOLE_ORGANISM-integrated.txt | PaxDb 6.0 (v6.0 download path); file last modified Wed, 01 Oct 2025 09:00:12 GMT | 2026-10-01 | 19,483 rows |
| PaxDb integrated whole-organism, mouse | https://pax-db.org/downloads/6.0/datasets/10090/10090-WHOLE_ORGANISM-integrated.txt | PaxDb 6.0 (v6.0 download path); file last modified Wed, 01 Oct 2025 08:58:51 GMT | 2026-10-01 | 19,774 rows |
| PaxDb integrated whole-organism, rat | https://pax-db.org/downloads/6.0/datasets/10116/10116-WHOLE_ORGANISM-integrated.txt | PaxDb 6.0 (v6.0 download path); file last modified Wed, 01 Oct 2025 08:59:08 GMT | 2026-10-01 | 16,184 rows |
| PeptideAtlas human phosphoproteome build | https://peptideatlas.org/builds/human/phospho/202512/atlas_build_606.tsv.zip | PeptideAtlas human phospho build 606 (directory 202512), file atlas_build_606.tsv.zip; last modified Tue, 16 Dec 2025 00:26:21 GMT | file on disk modified 2026-10-01T14:40:31 | 549,143,491 bytes, equal to Content-Length, zip test passes. Contains build tables (`atlas_build`, `biosequence_property_set`, ...); per-site confidence categories not yet examined (task H2). Fully retrieved and verified (size equals the Content-Length header; archive test passes). Later re-download attempts on 2026-10-01 failed on the network (DNS or read time-out) and were discarded; the verified file was not replaced. |
| UniProt ID mapping, HUMAN_9606 (for PaxDb to UniProt) | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz | UniProt Knowledgebase release 2026_03 (02-Sep-2026) | file on disk modified 2026-10-01T14:43:45 | 38,464,542 bytes, gzip test passes |
| UniProt ID mapping, MOUSE_10090 (for PaxDb to UniProt) | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz | UniProt Knowledgebase release 2026_03 (02-Sep-2026) | file on disk modified 2026-10-01T14:44:53 | 15,536,610 bytes, gzip test passes |
| UniProt ID mapping, RAT_10116 (for PaxDb to UniProt) | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/RAT_10116_idmapping.dat.gz | UniProt Knowledgebase release 2026_03 (02-Sep-2026) | file on disk modified 2026-10-01T14:45:12 | 9,415,486 bytes, gzip test passes |
| PhosphoSitePlus Phosphorylation_site_dataset | not downloaded (you place it) | | | not present in `data_private/` |

The first Pharos pull returned 10 of 20,412 targets because of a pagination error in my query (arguments on the wrong field); it was corrected and re-run, and the incomplete files are in `study/_archive/rev6_superseded/`. Pharos and GOA each list about twenty thousand human proteins (20,412 and 19,677), far fewer than the 75,248 human accessions in the site table, which matters for section 9.

## 3. Reconciliation of Rev 5 (task B1, details in `RECONCILIATION.md`)

- Text and Table 1 quote different unmatched estimates for 4 of 6 rows (Conservation; Disorder; Kinase motif, sequence-based; Kinase motif, database lookup): the text quotes all GOLD and BRONZE records, Table 1 the records whose protein has a primary index. Both values reproduce from the site table (differences below 1e-9). Disorder: 0.0310 against 0.0438. Matched columns are unaffected.
- Domain overlap used set-A tiers and the other four signals set-B tiers in Rev 5, so GOLD was not the same set of human records in both.
- Ortholog coverage: the text's 7.4% is 130,899 of 1,768,068 records (analysis set over the whole dataset); the 12.3% in the limitations is the human share only (131,687 of 1,071,725); dataset-wide coverage is 16.1%.
- Study clusters: 1,113 and 1,070 are the same count (distinct primary-cluster labels) over 1,768,068 and 1,768,020 records; publication-year resolution plays no part. Neither is the number of clusters under the merge rule (50,891 for 50,904 PMIDs).

## 4. Study-cluster rule (task B2): the threshold does not matter, reported as found

The rule was built with 884 sites per candidate PMID (0.05% of 1,768,068, restating a 500-site cutoff of an earlier human-only stage; the code default is 500; no data-derived justification for 500 is recorded). Candidate PMIDs: 757 at 500, 690 at 884, 674 at 1,000. The same 21 PMIDs form the same 8 multi-PMID clusters at all three, **study-cluster depth differs for 0 records at 500 and 0 at 1,000, and the tier counts are identical** ({'SILVER': 750802, 'BRONZE': 674827, 'GOLD': 342439}). Block G13 therefore reproduces G1 point for point and cannot show threshold dependence; this says nothing about whether the rule is right.

## 5. Primary comparisons, with the protein-cluster bootstrap (tasks B3 and C; blocks G1 to G3)

Index V50, set-B tiers for every signal, GOLD against BRONZE. 1,000 resamples of proteins (all sites of a sampled protein kept), the whole pipeline re-run in each resample; D0 and D1 share the same resamples, so attenuation is computed within the resample. D2 is propensity matching (the matcher of Rev 5, refitted in each resample). Decision rules are those of section 2 of the plan.

| Block | Signal | Crude (95% CI) | D1 CEM (95% CI) | Survives D1 | Attenuation D1 (95% CI) | D1 label | D2 propensity (95% CI) | D2 survives | D2 attenuation | Balance SMD D1 / D2 |
|---|---|---|---|---|---|---|---|---|---|---|
| G1 (all records) | Domain overlap (OR) | 6.065 (5.734 to 6.435) | 5.912 (5.588 to 6.254) | yes | 1.4% (1.1% to 1.9%) | partly attenuated (< 50%) | not fitted (singular matrix) | n/a | n/a | -0.003 / n/a |
| G1 (all records) | Conservation | 0.3986 (0.3878 to 0.4089) | 0.3993 (0.3897 to 0.4094) | yes | -0.2% (-1.9% to 1.3%) | not attenuated | 0.3846 (0.3719 to 0.3946) | yes | 3.5% (1.5% to 6.0%) | 0.004 / 0.015 |
| G1 (all records) | Disorder | 0.0438 (0.0390 to 0.0483) | 0.0457 (0.0413 to 0.0500) | yes | -4.5% (-8.2% to -1.7%) | not attenuated | not fitted (singular matrix) | n/a | n/a | 0.005 / n/a |
| G1 (all records) | Motif, sequence-based | 0.0617 (0.0583 to 0.0653) | 0.0618 (0.0586 to 0.0652) | yes | -0.1% (-1.5% to 1.3%) | not attenuated | 0.0450 (0.0408 to 0.0519) | yes | 27.2% (16.3% to 33.8%) | 0.004 / 0.009 |
| G1 (all records) | Motif, database lookup (not decision-bearing) | 0.7934 (0.7902 to 0.7965) | 0.7919 (0.7888 to 0.7951) | yes | 0.2% (0.1% to 0.3%) | partly attenuated (< 50%) | 0.7738 (0.7683 to 0.7769) | yes | 2.5% (2.1% to 3.1%) | 0.004 / 0.008 |
| G2 (human) | Conservation | 0.4560 (0.4424 to 0.4695) | 0.4501 (0.4367 to 0.4646) | yes | 1.3% (0.6% to 1.9%) | partly attenuated (< 50%) | 0.4313 (0.4166 to 0.4457) | yes | 5.4% (3.4% to 7.2%) | 0.007 / 0.009 |
| G2 (human) | Disorder | 0.0542 (0.0479 to 0.0606) | 0.0534 (0.0478 to 0.0599) | yes | 1.4% (-1.4% to 3.1%) | partly attenuated (< 50%) | not fitted (singular matrix) | n/a | n/a | 0.003 / n/a |
| G2 (human) | Motif, sequence-based | 0.0331 (0.0289 to 0.0371) | 0.0368 (0.0328 to 0.0406) | yes | -11.1% (-14.8% to -7.8%) | not attenuated | 0.0278 (0.0230 to 0.0365) | yes | 16.2% (-9.2% to 29.8%) | -0.003 / -0.000 |
| G2 (human) | Motif, database lookup (not decision-bearing) | 0.7983 (0.7937 to 0.8026) | 0.7964 (0.7917 to 0.8006) | yes | 0.2% (0.2% to 0.3%) | partly attenuated (< 50%) | not fitted (singular matrix) | n/a | n/a | -0.003 / n/a |
| G3 (annotated proteins (diagnostic)) | Domain overlap (OR) | 0.692 (0.647 to 0.739) | 0.679 (0.634 to 0.721) | yes | -5.5% (-10.2% to -2.8%) | not attenuated | 0.686 (0.642 to 0.737) | yes | -2.4% (-21.4% to 15.9%) | 0.005 / 0.003 |

**Result.** Of the 7 decision-bearing D1 contrasts in G1 and G2, 7 survive adjustment (the protein-cluster CI excludes the null) and 0 are largely explained by study effort (attenuation of 50% or more with a crude interval that excludes the null). The largest attenuation point estimate is 1.4%; 4 contrasts have negative attenuation (adjustment increased the estimate, reported as not attenuated). These are the results for the planned design (CEM on the index). Under propensity matching the largest attenuation is 27.8%.

**Propensity matching (D2) could not be fitted for 5 of 11 planned D2 contrasts** (G1-0003, G1-0009, G15-0352, G2-0021, G2-0027; `LinAlgError: Singular matrix`, the failure Rev 5 reported for the kinase-motif comparisons). It fitted for G1-0006, G1-0012, G1-0015, G2-0018, G2-0024, G3-0030. Cause, verified on the kinase-motif lookup and conservation records: the records with a non-canonical residue (K, H, R, C, I, D, G, M, N; 46 in the dataset) occur only among BRONZE records, so the residue-type term is perfectly separated and the routine's Newton step can fail; whether a given comparison fits depends on numerical details. This is reported as the result. **Beside each non-fitted design, and labelled post-hoc,** the same design restricted to serine, threonine and tyrosine records is reported below (rule 12; the original is not replaced), together with the S/T/Y version of the kinase-motif designs that did fit:

| Contrast | Signal | Population | D2 post-hoc S/T/Y estimate (95% CI) | Crude | Attenuation | Pairs / GOLD candidates | Note |
|---|---|---|---|---|---|---|---|
| G1-0003+posthoc_STY | Domain overlap (OR) | P-all (human) (S/T/Y only) | 7.041 (5.993 to 6.834) | 6.067 | -8.3% | 56,283 / 211,848 | point estimate lies outside its bootstrap interval (see below) |
| G1-0009+posthoc_STY | Disorder | P-all (S/T/Y only) | 0.0416 (0.0366 to 0.0470) | 0.0438 | 4.8% | 66,857 / 320,780 |  |
| G1-0012+posthoc_STY | Motif, sequence-based | P-all (S/T/Y only) | 0.0446 (0.0408 to 0.0514) | 0.0617 | 27.8% | 87,200 / 324,907 |  |
| G1-0015+posthoc_STY | Motif, database lookup (not decision-bearing) | P-all (S/T/Y only) | 0.7753 (0.7682 to 0.7768) | 0.7934 | 2.3% | 87,208 / 324,907 |  |
| G15-0352+posthoc_STY | Domain overlap (OR) | human (S/T/Y only) | 3.054 (2.655 to 3.045) | 2.495 | -22.1% | 49,030 / 153,738 | point estimate lies outside its bootstrap interval (see below) |
| G2-0021+posthoc_STY | Disorder | P-human (S/T/Y only) | 0.0565 (0.0480 to 0.0616) | 0.0542 | -4.4% | 40,991 / 208,322 |  |
| G2-0024+posthoc_STY | Motif, sequence-based | P-human (S/T/Y only) | 0.0281 (0.0234 to 0.0367) | 0.0331 | 15.2% | 56,278 / 211,848 |  |
| G2-0027+posthoc_STY | Motif, database lookup (not decision-bearing) | P-human (S/T/Y only) | 0.7785 (0.7716 to 0.7823) | 0.7983 | 2.5% | 56,303 / 211,848 |  |

**Caveat on D2.** The matcher inherited from Rev 5 keeps each GOLD record's 50 nearest controls on the logit; with 3,169 distinct logit values among 535,848 records (human domain overlap, serine/threonine/tyrosine records), ties are broken by record position, so only a fraction of GOLD records are matched (see the pair counts) and the point estimate depends on row order: the matched odds ratio is 7.01 to 7.08 in file order and 6.87 to 6.97 when the rows are shuffled (`results/propensity_order_dependence.json`). For 3 contrasts (G1-0003+posthoc_STY, G14-0337, G15-0352+posthoc_STY) the point estimate lies outside its own percentile interval. D2 is a sensitivity design; its interval is reported as computed and should be read with this in mind. I did not change the matcher.

## 6. What the tier definition and inference change (blocks G12, G14, G15)

- **Domain overlap depends strongly on the tier source.** Matched odds ratio 5.912 (5.588 to 6.254) with set-B tiers (G1), against 2.646 (2.491 to 2.817) with set-A tiers (G15) on the same index (the tier source, and with it the GOLD and BRONZE membership, differs), and 2.644 (2.482 to 2.807) with the depth-only tier (G12). Rev 5 reported 2.775 (set-A tiers, original index). The direction (above 1) is the same in all, the size is not.
- **Depth-only tier (G12).** Domain overlap 2.644 (2.482 to 2.807), attenuation -6.5%; Conservation 0.5308 (0.5227 to 0.5384), attenuation -12.1%; Disorder 0.1369 (0.1330 to 0.1416), attenuation 1.8%; Motif, sequence-based 0.1270 (0.1235 to 0.1304), attenuation 3.4%. All survive; none is largely explained.

**Inference alone (G14).** Rev 5 designs (original index; set-A tier for S1) with the protein-cluster bootstrap in place of the stratified site bootstrap. The point estimates reproduce Rev 5 (largest difference 8.9e-16) and only the intervals change:

| Signal | Design | Estimate | Rev 5 interval | Protein-cluster interval | Width ratio | Excludes null: Rev 5 / protein-cluster |
|---|---|---|---|---|---|---|
| Domain overlap | CEM primary index | 2.775 | 2.741 to 2.807 | 2.601 to 2.970 | 5.60 | yes / yes |
| Domain overlap | CEM secondary index | 0.767 | 0.754 to 0.779 | 0.722 to 0.808 | 3.51 | yes / yes |
| Conservation | CEM primary index | 0.3976 | 0.3912 to 0.4038 | 0.3878 to 0.4081 | 1.61 | yes / yes |
| Conservation | CEM secondary index | 0.4124 | 0.4060 to 0.4195 | 0.4023 to 0.4214 | 1.41 | yes / yes |
| Disorder | CEM primary index | 0.0422 | 0.0406 to 0.0437 | 0.0377 to 0.0471 | 3.10 | yes / yes |
| Disorder | CEM secondary index | 0.0540 | 0.0518 to 0.0560 | 0.0494 to 0.0588 | 2.27 | yes / yes |
| Motif, sequence-based | CEM primary index | 0.0631 | 0.0607 to 0.0653 | 0.0597 to 0.0664 | 1.47 | yes / yes |
| Motif, sequence-based | CEM secondary index | 0.1387 | 0.1352 to 0.1423 | 0.1344 to 0.1428 | 1.20 | yes / yes |
| Motif, database lookup | CEM primary index | 0.7920 | 0.7907 to 0.7934 | 0.7889 to 0.7951 | 2.27 | yes / yes |
| Motif, database lookup | CEM secondary index | 0.8253 | 0.8239 to 0.8269 | 0.8217 to 0.8287 | 2.33 | yes / yes |

Intervals widen by a factor of 1.20 to 5.60; the null is excluded in the same way as before in 10 of 10 CEM rows. Propensity rows of G14 are in `results/prop/`.

G14 propensity rows (Rev 5 design, original index; protein-cluster interval against the Rev 5 site-bootstrap interval):

| Contrast | Signal | This run | Rev 5 | Pairs now / Rev 5 |
|---|---|---|---|---|
| G14-0337 | Domain overlap | 3.043 (2.591 to 2.934) | 3.033 (2.966 to 3.108) | 62,490 / 62,493 |
| G14-0340 | Conservation | 0.3896 (0.3753 to 0.3972) | 0.3908 (0.3839 to 0.3983) | 30,995 / 31,001 |
| G14-0343 | Disorder | 0.0462 (0.0406 to 0.0495) | 0.0465 (0.0440 to 0.0487) | 81,549 / 81,575 |
| G14-0346 | Motif, sequence-based | 0.0452 (0.0428 to 0.0531) | not fitted (LinAlgError: Singular matrix) |  |
| G14-0349 | Motif, database lookup | not fitted (LinAlgError: Singular matrix) | not fitted (LinAlgError: Singular matrix) |  |

## 7. Two-way sensitivity (task B4, G1 D1)

V = V_protein + V_study - V_cell, bootstrap variances of the adjusted estimate (log scale for the odds ratio), study cluster = a site's primary cluster. Criteria fixed before the run: not feasible if V <= 0 or if at least 90% of records lie in single-record intersection cells.

| Signal | Adjusted | V protein | V study | V cell | V two-way | Single-record-cell share | Two-way 95% CI | Survives two-way | Width ratio vs protein-cluster |
|---|---|---|---|---|---|---|---|---|---|
| Domain overlap | 5.912 | 8.55e-04 | 1.83e-01 | 7.79e-04 | 1.83e-01 | 9.3% | 2.556 to 13.678 | yes | 16.4 |
| Conservation | 0.3993 | 2.76e-05 | 5.19e-03 | 1.96e-05 | 5.20e-03 | 17.6% | 0.2580 to 0.5407 | yes | 13.7 |
| Disorder | 0.0457 | 5.37e-06 | 1.52e-03 | 4.37e-06 | 1.52e-03 | 11.5% | -0.0308 to 0.1223 | no | 16.8 |
| Motif, sequence-based | 0.0618 | 2.82e-06 | 4.09e-03 | 2.50e-06 | 4.09e-03 | 10.1% | -0.0636 to 0.1872 | no | 38.1 |
| Motif, database lookup | 0.7919 | 2.55e-06 | 1.68e-03 | 2.48e-06 | 1.68e-03 | 10.1% | 0.7116 to 0.8723 | yes | 25.7 |

The study-cluster component dominates the variance (V study is 188 to 1451 times V protein). Under the two-way interval 2 of 5 contrasts no longer exclude the null: Disorder, Motif, sequence-based. This is a sensitivity analysis; the protein-cluster interval remains primary (decision D4). The normal approximation uses 1,000 study-cluster resamples with a few very large clusters (the largest holds 25.4% of all 1,768,068 sites), so the interval is approximate.

## 8. Effect sizes for the index properties (task B5; replaces the p-values)

87,338 proteins, protein bootstrap, 1,000 resamples, no p-values. Rev 5 values are reproduced for the original index (eta-squared 0.0067, Spearman 0.1206, rank-biserial 0.0458). Sign convention of the rank-biserial r (Rev 5's): positive means GOLD proteins have the *lower* index.

| Effect size | Estimate | 95% CI |
|---|---|---|
| Kruskal-Wallis eta-squared of V0 (Rev 5 index) by highest tier | 0.0067 | 0.0057 to 0.0078 |
| rank-biserial r, V0 (Rev 5 index), GOLD vs BRONZE proteins (highest tier) | 0.0458 | 0.0357 to 0.0556 |
| Spearman rho, V0 (Rev 5 index) vs secondary index | 0.1206 | 0.1146 to 0.1270 |
| point-biserial r, domain-annotation presence vs V0 (Rev 5 index) | -0.1530 | -0.1612 to -0.1435 |
| Kruskal-Wallis eta-squared of V50 (Rev 6 primary index) by highest tier | 0.0047 | 0.0038 to 0.0057 |
| rank-biserial r, V50 (Rev 6 primary index), GOLD vs BRONZE proteins (highest tier) | 0.0475 | 0.0380 to 0.0573 |
| Spearman rho, V50 (Rev 6 primary index) vs secondary index | 0.0756 | 0.0689 to 0.0824 |
| point-biserial r, domain-annotation presence vs V50 (Rev 6 primary index) | -0.1280 | -0.1369 to -0.1193 |
| point-biserial r, domain-annotation presence vs secondary index | 0.6486 | 0.6438 to 0.6536 |
| Kruskal-Wallis eta-squared of secondary index by highest tier | 0.7660 | 0.7636 to 0.7682 |
| rank-biserial r, secondary index, GOLD vs BRONZE proteins | -0.9986 | -0.9988 to -0.9984 |

## 9. Index variants (task D)

| Variant | PMIDs removed (> n genes) | Protein-PMID links removed | Proteins with zero PMIDs | Spearman with original index |
|---|---|---|---|---|
| V0 (original) | 0 | 0 | | 1 |
| V20 | 10,240 | 3,300,460 | 6,781 | 0.9201 |
| V50 | 3,527 | 2,956,031 | 5,320 | 0.9290 |
| V100 | 2,227 | 2,771,305 | 4,638 | 0.9388 |
| VT2025 | | | 5,320 missing | 0.9096 (n = 81,967) |
| VT2026 | | | 5,320 missing | 0.9117 (n = 82,018) |

**Time normalisation.** t = 2025 - first year + 1 (2026 for the sensitivity); asserted t >= 1 and finite in all rows used (t from 1 to 106). First year = earliest year among all of a protein's bulk-corrected PMIDs (1,979,539 PMIDs dated from NCBI E-utilities and `pmid_year.parquet`; none left undated). **Not filled or imputed:** 5,320 of 87,338 mapped proteins (4,527 of 45,871 genes) have no PMID left after the V50 correction, so VT is missing for them; all of them are proteins with zero publications. A further 51 proteins have a first year of 2026, which gives t = 0 against 2025 and is left missing (none for the 2026 reference). These proteins are excluded only from analyses that use VT. Because the excluded proteins are exactly the least-studied ones, VT comparisons are on a differently composed population (see G4).

**Pharos TDL (D3).** 18,624 of 48,162 human proteins with an index match a Pharos entry (Pharos lists 20,412 targets, of which 18,624 match index proteins; the others are not in the index or the site table). TDL is built partly from publication-derived scores, so agreement with a publication-count index is supportive, not an independent calibration. Mean index by TDL and AUROC for separating non-Tdark from Tdark (common set of 17,819 proteins with all variants defined; bootstrap 95% CI):

| Variant | Tclin | Tchem | Tbio | Tdark | AUROC (95% CI) | Own-set AUROC |
|---|---|---|---|---|---|---|
| V0 | 2.418 | 2.254 | 1.876 | 1.205 | 0.9104 (0.9049 to 0.9152) | 0.9242 (n = 18,624) |
| V20 | 2.286 | 2.063 | 1.505 | 0.513 | 0.9540 (0.9506 to 0.9575) | 0.9619 (n = 18,624) |
| V50 | 2.300 | 2.083 | 1.540 | 0.612 | 0.9509 (0.9474 to 0.9544) | 0.9595 (n = 18,624) |
| V100 | 2.315 | 2.098 | 1.568 | 0.659 | 0.9496 (0.9459 to 0.9533) | 0.9584 (n = 18,624) |
| VT2025 | 0.876 | 0.728 | 0.404 | 0.086 | 0.9288 (0.9247 to 0.9329) | 0.9288 (n = 17,819) |
| VT2026 | 0.867 | 0.719 | 0.396 | 0.080 | 0.9335 (0.9295 to 0.9375) | 0.9328 (n = 17,825) |

**GO check (D4) and choice of variant (D5).** Pre-specified: proteins absent from `goa_human.gaf` count 0 terms (29,828 of 48,162 index proteins have no GOA row; the file lists 19,677 proteins). Under that coding the Spearman correlation between index and GO-term count is negative for every variant (V0 -0.062, V20 -0.020, V50 -0.020, V100 -0.023, VT2025 -0.027), so it measures coverage of the file, not annotation depth. The pre-specified rule nevertheless ranks **V20 first on the Tdark AUROC and V20 first on the GO correlation, which gives the rule's outcome: "Both measures rank V20 first and it is not V50: report G1 with V50 and with V20, and state both."** The V20 and V50 intervals overlap (AUROC 0.9540 against 0.9509). Block G4 reports S1 to S4 under V20; S5 is not decision-bearing and was not rerun.
**Post-hoc, labelled, beside the pre-specified result:** restricting the GO check to proteins with at least one GOA row, the correlations are positive (V0 0.705, V20 0.719, V50 0.726, V100 0.730, VT2025 0.703) and rank **V100** first, while the Tdark AUROC still ranks V20 first. Under this version the two measures disagree, and the rule would report all variants and choose none. **The two readings give different instructions, so I have not acted on either beyond what G4 already shows (V20, V50 and V100 give the same conclusions); please decide.**

### G4: S1 to S4 under each index variant (D3, CEM; protein-cluster intervals)

| Signal | V0 | V20 | V50 (G1) | V100 | VT2025 | VT2026 |
|---|---|---|---|---|---|
| Domain overlap | 6.042 (5.692 to 6.398) | 5.898 (5.558 to 6.242) | 5.912 (5.588 to 6.254) | 5.921 (5.601 to 6.272) | 5.969 (5.623 to 6.321) | 5.970 (5.630 to 6.325) |
| Conservation | 0.3976 (0.3884 to 0.4084) | 0.3994 (0.3896 to 0.4107) | 0.3993 (0.3897 to 0.4094) | 0.3988 (0.3889 to 0.4094) | 0.3988 (0.3897 to 0.4094) | 0.3986 (0.3890 to 0.4098) |
| Disorder | 0.0422 (0.0383 to 0.0472) | 0.0458 (0.0415 to 0.0506) | 0.0457 (0.0413 to 0.0500) | 0.0456 (0.0414 to 0.0503) | 0.0446 (0.0403 to 0.0496) | 0.0447 (0.0403 to 0.0495) |
| Motif, sequence-based | 0.0631 (0.0600 to 0.0665) | 0.0617 (0.0586 to 0.0651) | 0.0618 (0.0586 to 0.0652) | 0.0619 (0.0587 to 0.0654) | 0.0598 (0.0565 to 0.0630) | 0.0598 (0.0564 to 0.0630) |

All 20 G4 contrasts survive (interval excludes the null: 20 of 20); 0 are largely explained. The VT rows use only proteins with a defined VT (the 5,320 zero-PMID proteins drop out), so their crude values differ from the other columns. No variant changes a decision (all survive and none is largely explained, as in G1).

## 10. Protein abundance (task E, G5)

PaxDb 6.0 integrated whole-organism abundance mapped to UniProt through the UniProt ID mapping (id type STRING; maximum kept where several entries map to one accession, a rule set before analysis). Coverage:

| Species | PaxDb entries | Unmapped | Proteins with abundance / in site table | Records with abundance | GOLD/BRONZE records with abundance and V50 index |
|---|---|---|---|---|---|
| human | 19,483 | 545 | 18,524 / 75,248 (24.6%) | 555,053 (51.8%) | 276,157 |
| mouse | 19,774 | 1,231 | 16,335 / 45,003 (36.3%) | 281,362 (52.7%) | 135,525 |
| rat | 16,184 | 2,363 | 9,680 / 16,217 (59.7%) | 124,030 (76.5%) | 24,268 |

Population P-abund: GOLD or BRONZE records of proteins with both V50 and abundance, identical records in all four columns. D1 matches on index decile, D5 on abundance decile alone, D4 on index decile x abundance decile (species and residue type exact). Abundance enters as log10 ppm. Mouse and rat were run (S1 is human-only and was not run for them).

| Species | Signal | Crude | D1 index | D5 abundance | D4 index x abundance | Attenuation D1 / D5 / D4 | Survives D1 / D5 / D4 | SMD of index after D5; of abundance after D1 |
|---|---|---|---|---|---|---|---|---|
| human | Domain overlap | 0.741 | 0.711 (0.670 to 0.753) | 0.652 (0.617 to 0.691) | 0.644 (0.609 to 0.676) | -13.9% / -42.9% / -46.8% | yes / yes / yes | 0.059; 0.428 |
| human | Conservation | 0.4755 | 0.4720 (0.4610 to 0.4842) | 0.4737 (0.4622 to 0.4855) | 0.4753 (0.4644 to 0.4878) | 0.7% / 0.4% / 0.0% | yes / yes / yes | 0.069; 0.427 |
| mouse | Conservation | 0.4352 | 0.4319 (0.4186 to 0.4446) | 0.4161 (0.4033 to 0.4298) | 0.4144 (0.4004 to 0.4277) | 0.8% / 4.4% / 4.8% | yes / yes / yes | 0.082; 0.311 |
| rat | Conservation | 0.1473 | 0.1394 (0.1088 to 0.1718) | 0.1103 (0.0863 to 0.1479) | 0.1235 (0.0954 to 0.1566) | 5.3% / 25.1% / 16.2% | yes / yes / yes | 0.193; 0.436 |
| human | Disorder | 0.1362 | 0.1352 (0.1279 to 0.1413) | 0.1419 (0.1353 to 0.1485) | 0.1426 (0.1371 to 0.1487) | 0.7% / -4.2% / -4.7% | yes / yes / yes | 0.054; 0.424 |
| mouse | Disorder | 0.0969 | 0.0901 (0.0835 to 0.0968) | 0.0896 (0.0834 to 0.0957) | 0.0897 (0.0837 to 0.0950) | 7.0% / 7.5% / 7.4% | yes / yes / yes | 0.076; 0.384 |
| rat | Disorder | 0.0203 | 0.0320 (0.0163 to 0.0517) | 0.0357 (0.0177 to 0.0569) | 0.0595 (0.0449 to 0.0799) | -57.4% / -75.3% / -192.7% | yes / yes / yes | 0.282; 0.511 |
| human | Motif, sequence-based | 0.1006 | 0.1041 (0.0991 to 0.1091) | 0.1109 (0.1060 to 0.1166) | 0.1120 (0.1074 to 0.1174) | -3.5% / -10.2% / -11.4% | yes / yes / yes | 0.059; 0.428 |
| mouse | Motif, sequence-based | 0.1803 | 0.1784 (0.1721 to 0.1856) | 0.1844 (0.1774 to 0.1913) | 0.1848 (0.1772 to 0.1917) | 1.1% / -2.3% / -2.5% | yes / yes / yes | 0.087; 0.406 |
| rat | Motif, sequence-based | 0.2005 | 0.2009 (0.1755 to 0.2193) | 0.2208 (0.1997 to 0.2421) | 0.2067 (0.1857 to 0.2340) | -0.2% / -10.1% / -3.1% | yes / yes / yes | 0.295; 0.528 |

30 of 30 adjusted contrasts survive; 0 have an attenuation of 50% or more. Matching on the index alone leaves abundance imbalanced (SMD 0.31 to 0.53, last column), matching on abundance alone leaves the index imbalanced in rat (3 D5 contrasts above |SMD| 0.1 for the index; reported, not re-tuned). The attenuations range from -192.7% (an increase away from the null) to 25.1%; the largest increases occur for rat disorder (crude estimate 0.0203, a small population) and for human domain overlap. **The domain-overlap comparison in P-abund is not comparable with G1:** restricting to proteins with an abundance changes the crude odds ratio from 6.065 (G1) to 0.741, i.e. across 1; attenuation is relative to the crude value within each population. Protein-level Spearman correlation of abundance: with the V50 index 0.370 (0.355 to 0.383) in human, 0.195 in mouse, 0.231 in rat; with the outcomes in human: domain overlap 0.078, conservation 0.213, disorder 0.066, motif, sequence-based -0.059, motif, database lookup 0.303 (all intervals in `results/g5_diagnostics.json`).

## 11. Departures, problems and what I need from you

**Logged in the plan (section 15), none changes a result:** B2 before the manifest; G7's 11th removal duplicates the 2nd; G13 identical by construction; S4 record set; D2 first-year definition (a 5-lowest-PMID proxy missed an earlier known year for at least 3,660 of 82,018 proteins and was replaced before use); the Pharos pull. Also logged: the engine clipping fix (no estimate changed, verified by rerun); the abundance modes; the propensity non-fits and post-hoc S/T/Y rows (section 5); the GO coding (section 9); the G14 secondary-index restriction; the criteria fixed before the B4 and D5 results.

**Decisions for you (rule 9):**

1. **Choice of index variant (D5).** The pre-specified rule says report V50 and V20; with the GO coding corrected (post-hoc) it says choose none. V20, V50 and V100 give the same decisions in G4. Do you want the rule applied as written, or the GO coding changed and the rule rerun (a new labelled analysis)?
2. **Propensity matching.** 5 of 11 planned D2 contrasts could not be fitted by the Rev 5 routine because the non-canonical residues are all BRONZE. I report that and the post-hoc S/T/Y fits. Should the S/T/Y-only version become the planned D2 for the remaining blocks (G5-G9 do not use D2), or stay a labelled post-hoc?
3. **Matcher tie-breaking.** The inherited 50-nearest-neighbour matcher depends on record order and matches about a quarter of GOLD records. I kept it for fidelity with Rev 5. Do you want a variant with random tie-breaking as a labelled post-hoc?
4. **PhosphoSitePlus** is needed for tasks H and I; the file is not in `data_private/` yet. PeptideAtlas build 606 is on disk; whether it holds per-site confidence categories is untested.
5. The Rev 5 text and Table 1 disagree on four unmatched estimates, and the 7.4%, 12.3%, 1,113 and 1,070 statements need denominators (`RECONCILIATION.md`). No manuscript edit is made in this cycle (rule 13).

## 12. Status of the manifest

| block   |   not fitted |   not run |   pending (task F) |   pending (task H) |   pending (task I) |   pending (task J) |   run |
|:--------|-------------:|----------:|-------------------:|-------------------:|-------------------:|-------------------:|------:|
| G1      |            2 |         0 |                  0 |                  0 |                  0 |                  0 |    13 |
| G10     |            0 |         0 |                  0 |                  0 |                  4 |                  0 |     0 |
| G11     |            0 |         0 |                  0 |                  0 |                  0 |                  8 |     0 |
| G12     |            0 |         0 |                  0 |                  0 |                  0 |                  0 |     8 |
| G13     |            0 |         0 |                  0 |                  0 |                  0 |                  0 |    16 |
| G14     |            1 |         0 |                  0 |                  0 |                  0 |                  0 |    14 |
| G15     |            1 |         0 |                  0 |                  0 |                  0 |                  0 |     2 |
| G2      |            2 |         0 |                  0 |                  0 |                  0 |                  0 |    10 |
| G3      |            0 |         0 |                  0 |                  0 |                  0 |                  0 |     3 |
| G4      |            0 |         0 |                  0 |                  0 |                  0 |                  0 |    20 |
| G5      |            0 |         8 |                  0 |                  0 |                  0 |                  0 |    40 |
| G6      |            0 |         0 |                 24 |                  0 |                  0 |                  0 |     0 |
| G7      |            0 |         0 |                 88 |                  0 |                  0 |                  0 |     0 |
| G8      |            0 |         0 |                 48 |                  0 |                  0 |                  0 |     0 |
| G9      |            0 |         0 |                  0 |                 48 |                  0 |                  0 |     0 |

Rows marked `pending` belong to tasks F to J, which start after this report. G7's duplicate removal and G13's identical results are in the table as run.

