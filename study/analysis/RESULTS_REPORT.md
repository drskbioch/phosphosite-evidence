# Results, Rev 6 reanalysis

Branch `rev6-reanalysis`. Written by `scripts/analysis/build_results.py` from the files in `study/analysis/results/`; every number below is read from a result file or computed from one in that script. No manuscript file, reference or legacy file was touched. The plan is `ANALYSIS_PLAN.md` (deviation log in its section 15, including the amendments of 2026-10-01). `INTERIM_REPORT.md` still holds the interim results (G1 to G5, G12 to G15, index variants, abundance); this report adds the amended and later blocks and re-states the main table with the gene-cluster interval. `results/manifest_status.csv` lists all 534 manifest rows with their status.

**Final under Branch B (section 14).** Task H (PhosphoSitePlus measures, block G9, 24 rows), Task I (block G10, 4 rows) and G20-0488 are **not run: source file not available (PhosphoSitePlus v7 site provides no bulk download; access requested from the provider)**. They are listed in Future work (section 15), to be run if the provider supplies the file (possibly during revision); `h_measures.py` is written and has run for the UniProt measure, `cond_logit.py` is validated on real data (G21 check), and the Task I runner is not written. Everything else in the frozen plan has been run.

## Headline results

1. **Matching on the study-effort index (CEM): 7 of 7 decision-bearing D1 contrasts in G1 and G2 survive with the gene-cluster interval (primary from the amendment) and 7 of 7 with the accession-cluster interval; 0 are largely explained; the largest attenuation is 1.4%.** **Under the two-way (protein by study cluster) interval, 2 of 4 decision-bearing G1 contrasts exclude the null (Domain overlap, Conservation); Disorder, Motif, sequence-based do not (claim C4.5; approximate interval, sensitivity analysis; section 3).**

2. **Domain overlap is population-dependent** (it is not on one side of the null in every population): the crude odds ratio is 6.065 in the full population, 0.741 among proteins with an abundance, 0.692 among proteins with a domain-annotated site and 0.746 among Swiss-Prot reviewed canonical accessions. The other signals: Conservation direction-consistent, Disorder direction-consistent, Motif, sequence-based direction-consistent, Motif, database lookup direction-consistent.

3. **The primary index can carry little confounding.** Before matching GOLD and BRONZE differ in it (41 of 78 decision-bearing contrasts have |SMD| of 0.1 or more, largest 0.270; in G1 and G2: Conservation P-all 0.160, Conservation P-human 0.187), and matching removes the imbalance (largest SMD after 0.017); but within a tier the index correlates with an outcome at only -0.0881 to 0.0671 (section 4). Adjusting for it therefore moves the estimates little, which is what one expects whether or not site-level study effort confounds the comparison: the design cannot tell the two apart.

4. **The matcher the Rev 5 Methods describe (1:1, without replacement, caliper, random order; D2b) gives the same decisions as CEM:** 7 of 7 decision-bearing D2b contrasts survive, the largest attenuation is 2.9%, 0 point estimates lie outside their own interval (section 5).

5. **The tier definition changes the size of the signals more than adjustment for the index does:** like-for-like on the log scale, the spread of the adjusted estimates across tier definitions exceeds the largest change caused by effort adjustment within a definition in 7 of 7 signal-by-population combinations (natural-log scale; spread 0.432 to 2.840, as plain ratios 1.54 to 17.11; within-definition change 0.015 to 0.399, as plain ratios 1.015 to 1.490; section 6). The 17-fold spread for disorder in P-all rests on a near-zero agreement-only estimate; without it the spread is about 3-fold (2.99).

6. **Not clean:** the pre-specified propensity design (D2) could not be fitted for 5 of 11 planned contrasts; the pre-specified GO coding measures GOA file coverage and cannot arbitrate between index variants (section 1); the optional UniProt contrast is partly circular for conservation (section 9); of 7 evaluated predictor settings 2 show a gap equivalent to zero, 5 are not detectable and 0 show a gap (section 11); PhosphoSitePlus-based measures and Task I are not run: source file not available (PhosphoSitePlus v7 site provides no bulk download; access requested from the provider).

## 1. Decisions applied and fixes

- **D5 (index variant), as decided.** The pre-specified rule is applied as written. Both pre-specified measures rank V20 first (Tdark AUROC 0.9540 against 0.9509 for V50; GO Spearman -0.0203 against -0.0204), so G1 is reported with V50 (primary) and V20 (table below). **The pre-specified GO coding measures GOA file coverage and cannot arbitrate between variants:** 29,828 of 48,162 human proteins with an index have no row in `goa_human.gaf` and are coded as 0 GO terms, so the correlation is negative for every variant (V20 -0.0203, V50 -0.0204, V100 -0.0227, VT2025 -0.0273, V0 -0.0618). The GO count restricted to proteins with a GOA row stays as a labelled post-hoc analysis beside it (it ranks V100 first on GO and V20 first on the Tdark AUROC; the measures disagree there). No variant changes a decision in G4 (INTERIM_REPORT.md, section 9).

G1 D1 with the two indices the rule names (accession-cluster intervals, as in the interim report):

| Signal | V50 (G1) | V20 (G4) | Survives V50 | Survives V20 | Attenuation V50 | Attenuation V20 |
|---|---|---|---|---|---|---|
| Domain overlap | 5.912 (5.588 to 6.254) | 5.898 (5.558 to 6.242) | yes | yes | 1.4% | 1.5% |
| Conservation | 0.3993 (0.3897 to 0.4094) | 0.3994 (0.3896 to 0.4107) | yes | yes | -0.2% | -0.2% |
| Disorder | 0.0457 (0.0413 to 0.0500) | 0.0458 (0.0415 to 0.0506) | yes | yes | -4.5% | -4.6% |
| Motif, sequence-based | 0.0618 (0.0586 to 0.0652) | 0.0617 (0.0586 to 0.0651) | yes | yes | -0.1% | 0.1% |

- **D2 (propensity matching), as decided.** The 5 non-fits of the Rev 5 routine (G1-0003, G1-0009, G15-0352, G2-0021, G2-0027; `LinAlgError: Singular matrix`) are the result of the pre-specified design (INTERIM_REPORT.md, section 5). The post-hoc S/T/Y-only rows of the interim stage stay beside them; the new matcher is design D2b (section 5).
- **Fix, unrendered placeholders.** The interim builder left a template expression unrendered in section 7; fixed and re-rendered (deviation log). A search of every generated report in `study/analysis/` for unrendered template expressions finds none. The same review found one number typed into the interim builder (the count of records with an undefined outcome); it is now read from `results/validate_engine.json`.

## 2. Cluster unit: gene instead of accession (correction; block G17)

**Recorded in the deviation log as a correction made after the accession-level results were seen:** accessions of one gene are not independent, so the primary bootstrap now resamples gene clusters (accessions that share an Entrez gene ID; connected components of the accession-gene graph). All new blocks use the gene cluster. The accession-cluster interval is kept beside it everywhere.

- **Accessions per gene.** Human: 75,248 accessions in the site table, 48,162 with a gene mapping in 18,757 gene clusters (mean 2.57, median 2, maximum 50; 56.7% of clusters have more than one accession). All species: 136,468 accessions, 87,338 mapped to 45,545 gene clusters (mean 1.92, maximum 50). 27,086 human accessions have no gene mapping and stay single-accession clusters; none of them has an index, so none enters a decision-bearing contrast.
- **Records in the 10 largest genes (all records).** 22,074 of 1,768,068 (1.25%); the 10 largest accessions hold 12,320 (0.70%). Largest gene cluster: G7273 with 5,174 records in 15 accessions.
- **Records in the 10 largest genes (human records).** 20,100 of 1,071,725 (1.88%); the 10 largest accessions hold 10,367 (0.97%). Largest gene cluster: G7273 with 5,174 records in 15 accessions.
- **Records in the 10 largest genes (GOLD or BRONZE records with an index).** 13,109 of 803,787 (1.63%); the 10 largest accessions hold 6,344 (0.79%). Largest gene cluster: G7273 with 3,247 records in 15 accessions.

Accession-cluster against gene-cluster interval (same point estimates; 1,000 resamples each; width ratio = gene / accession):

| Block | Signal | Population | Design | Estimate | Accession-cluster 95% CI | Gene-cluster 95% CI | Width ratio | Excludes null (accession) | Excludes null (gene) |
|---|---|---|---|---|---|---|---|---|---|
| G1 | Domain overlap | P-all (human) | D0 | 6.065 | 5.734 to 6.435 | 5.716 to 6.403 | 0.981 | yes | yes |
| G1 | Domain overlap | P-all (human) | D1 | 5.912 | 5.588 to 6.254 | 5.571 to 6.241 | 1.007 | yes | yes |
| G1 | Conservation | P-all | D0 | 0.3986 | 0.3878 to 0.4089 | 0.3876 to 0.4103 | 1.075 | yes | yes |
| G1 | Conservation | P-all | D1 | 0.3993 | 0.3897 to 0.4094 | 0.3892 to 0.4099 | 1.045 | yes | yes |
| G1 | Disorder | P-all | D0 | 0.0438 | 0.0390 to 0.0483 | 0.0388 to 0.0487 | 1.062 | yes | yes |
| G1 | Disorder | P-all | D1 | 0.0457 | 0.0413 to 0.0500 | 0.0411 to 0.0512 | 1.161 | yes | yes |
| G1 | Motif, sequence-based | P-all | D0 | 0.0617 | 0.0583 to 0.0653 | 0.0584 to 0.0649 | 0.936 | yes | yes |
| G1 | Motif, sequence-based | P-all | D1 | 0.0618 | 0.0586 to 0.0652 | 0.0587 to 0.0652 | 0.993 | yes | yes |
| G1 | Motif, database lookup (not decision-bearing) | P-all | D0 | 0.7934 | 0.7902 to 0.7965 | 0.7902 to 0.7966 | 1.022 | yes | yes |
| G1 | Motif, database lookup (not decision-bearing) | P-all | D1 | 0.7919 | 0.7888 to 0.7951 | 0.7887 to 0.7951 | 1.026 | yes | yes |
| G2 | Conservation | P-human | D0 | 0.4560 | 0.4424 to 0.4695 | 0.4418 to 0.4691 | 1.006 | yes | yes |
| G2 | Conservation | P-human | D1 | 0.4501 | 0.4367 to 0.4646 | 0.4349 to 0.4642 | 1.051 | yes | yes |
| G2 | Disorder | P-human | D0 | 0.0542 | 0.0479 to 0.0606 | 0.0468 to 0.0618 | 1.183 | yes | yes |
| G2 | Disorder | P-human | D1 | 0.0534 | 0.0478 to 0.0599 | 0.0474 to 0.0611 | 1.133 | yes | yes |
| G2 | Motif, sequence-based | P-human | D0 | 0.0331 | 0.0289 to 0.0371 | 0.0293 to 0.0370 | 0.927 | yes | yes |
| G2 | Motif, sequence-based | P-human | D1 | 0.0368 | 0.0328 to 0.0406 | 0.0331 to 0.0404 | 0.939 | yes | yes |
| G2 | Motif, database lookup (not decision-bearing) | P-human | D0 | 0.7983 | 0.7937 to 0.8026 | 0.7937 to 0.8026 | 1.002 | yes | yes |
| G2 | Motif, database lookup (not decision-bearing) | P-human | D1 | 0.7964 | 0.7917 to 0.8006 | 0.7918 to 0.8007 | 0.996 | yes | yes |
| G3 | Domain overlap (not decision-bearing) | P-annot | D0 | 0.692 | 0.647 to 0.739 | 0.655 to 0.733 | 0.838 | yes | yes |
| G3 | Domain overlap (not decision-bearing) | P-annot | D1 | 0.679 | 0.634 to 0.721 | 0.642 to 0.717 | 0.851 | yes | yes |
| G5 | Domain overlap | P-abund (human) | D0 | 0.741 | 0.700 to 0.786 | 0.701 to 0.783 | 0.958 | yes | yes |
| G5 | Domain overlap | P-abund (human) | D1 | 0.711 | 0.670 to 0.753 | 0.672 to 0.752 | 0.962 | yes | yes |
| G5 | Conservation | P-abund (human) | D0 | 0.4755 | 0.4643 to 0.4875 | 0.4628 to 0.4869 | 1.034 | yes | yes |
| G5 | Conservation | P-abund (human) | D1 | 0.4720 | 0.4610 to 0.4842 | 0.4597 to 0.4830 | 1.005 | yes | yes |
| G5 | Disorder | P-abund (human) | D0 | 0.1362 | 0.1294 to 0.1422 | 0.1298 to 0.1428 | 1.021 | yes | yes |
| G5 | Disorder | P-abund (human) | D1 | 0.1352 | 0.1279 to 0.1413 | 0.1284 to 0.1421 | 1.016 | yes | yes |
| G5 | Motif, sequence-based | P-abund (human) | D0 | 0.1006 | 0.0957 to 0.1054 | 0.0954 to 0.1052 | 1.002 | yes | yes |
| G5 | Motif, sequence-based | P-abund (human) | D1 | 0.1041 | 0.0991 to 0.1091 | 0.0988 to 0.1090 | 1.015 | yes | yes |

The gene-cluster interval is 0.838 to 1.183 times the width of the accession-cluster interval (median 1.006, about 1.01; three decimals are shown so that a ratio of about 1.01 is not rounded without its range); the null is excluded under both intervals in 28 of 28 rows and the two intervals disagree on the null in 0.

## 3. Main results: pre-specified rule, direction across populations, two-way interval (G1 to G3, P-abund)

Index V50, CEM (D1), set-B tiers. **Survives** is the pre-specified rule (interval excludes the null; shown for the accession-cluster interval of the pre-specified protein-cluster bootstrap and for the gene-cluster interval, which is primary from the amendment). **Direction-consistent** is a second column: the adjusted estimate lies on the same side of the null in every population where the signal is defined (G1, G2, G3, P-abund in human, mouse and rat under D1, D4 and D5; the estimates used are in the next table). **Two-way** is the protein-by-study-cluster interval of B4 (G1 D1 only; not run for the other populations); **study share** is V_study / V_two-way. The gene-cluster rule stays primary (decision D4).

| Block | Signal | Population | Crude (gene-cluster CI) | D1 adjusted (gene-cluster CI) | D1 accession-cluster CI | Survives (accession) | Survives (gene) | Attenuation (gene) | Two-way CI (excludes null) | Study share of variance | Direction-consistent | Label |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| G1 | Domain overlap | P-all (human) | 6.065 (5.716 to 6.403) | 5.912 (5.571 to 6.241) | 5.588 to 6.254 | yes | yes | 1.4% | 2.556 to 13.678 (yes) | 100.0% | no | population-dependent |
| G1 | Conservation | P-all | 0.3986 (0.3876 to 0.4103) | 0.3993 (0.3892 to 0.4099) | 0.3897 to 0.4094 | yes | yes | -0.2% | 0.2580 to 0.5407 (yes) | 99.8% | yes | direction-consistent |
| G1 | Disorder | P-all | 0.0438 (0.0388 to 0.0487) | 0.0457 (0.0411 to 0.0512) | 0.0413 to 0.0500 | yes | yes | -4.5% | -0.0308 to 0.1223 (no) | 99.9% | yes | direction-consistent |
| G1 | Motif, sequence-based | P-all | 0.0617 (0.0584 to 0.0649) | 0.0618 (0.0587 to 0.0652) | 0.0586 to 0.0652 | yes | yes | -0.1% | -0.0636 to 0.1872 (no) | 100.0% | yes | direction-consistent |
| G1 | Motif, database lookup (not decision-bearing) | P-all | 0.7934 (0.7902 to 0.7966) | 0.7919 (0.7887 to 0.7951) | 0.7888 to 0.7951 | yes | yes | 0.2% | 0.7116 to 0.8723 (yes) | 100.0% | yes | direction-consistent |
| G2 | Conservation | P-human | 0.4560 (0.4418 to 0.4691) | 0.4501 (0.4349 to 0.4642) | 0.4367 to 0.4646 | yes | yes | 1.3% | not run | n/a | yes | direction-consistent |
| G2 | Disorder | P-human | 0.0542 (0.0468 to 0.0618) | 0.0534 (0.0474 to 0.0611) | 0.0478 to 0.0599 | yes | yes | 1.4% | not run | n/a | yes | direction-consistent |
| G2 | Motif, sequence-based | P-human | 0.0331 (0.0293 to 0.0370) | 0.0368 (0.0331 to 0.0404) | 0.0328 to 0.0406 | yes | yes | -11.1% | not run | n/a | yes | direction-consistent |
| G2 | Motif, database lookup (not decision-bearing) | P-human | 0.7983 (0.7937 to 0.8026) | 0.7964 (0.7918 to 0.8007) | 0.7917 to 0.8006 | yes | yes | 0.2% | not run | n/a | yes | direction-consistent |
| G3 | Domain overlap (not decision-bearing) | P-annot | 0.692 (0.655 to 0.733) | 0.679 (0.642 to 0.717) | 0.634 to 0.721 | yes | yes | -5.5% | not run | n/a | no | population-dependent |
| G5 | Domain overlap | P-abund (human) | 0.741 (0.701 to 0.783) | 0.711 (0.672 to 0.752) | 0.670 to 0.753 | yes | yes | -13.9% | not run | n/a | no | population-dependent |
| G5 | Conservation | P-abund (human) | 0.4755 (0.4628 to 0.4869) | 0.4720 (0.4597 to 0.4830) | 0.4610 to 0.4842 | yes | yes | 0.7% | not run | n/a | yes | direction-consistent |
| G5 | Disorder | P-abund (human) | 0.1362 (0.1298 to 0.1428) | 0.1352 (0.1284 to 0.1421) | 0.1279 to 0.1413 | yes | yes | 0.7% | not run | n/a | yes | direction-consistent |
| G5 | Motif, sequence-based | P-abund (human) | 0.1006 (0.0954 to 0.1052) | 0.1041 (0.0988 to 0.1090) | 0.0991 to 0.1091 | yes | yes | -3.5% | not run | n/a | yes | direction-consistent |

**Reading.** Excluding the null: 7 of 7 under the accession-cluster interval, 7 of 7 under the gene-cluster interval, 2 of 4 under the two-way interval (G1 contrasts that do not: Disorder, Motif, sequence-based). The study-cluster component carries 99.8% to 100.0% of the two-way variance. The two-way interval rests on 489 to 871 study clusters, a few of them very large, so it is approximate (INTERIM_REPORT.md, section 7); it is a sensitivity, not the primary rule.

Adjusted estimates used for the direction column (all from the CEM results; G1 and G2 are D1, P-abund rows are D1, D4 and D5):

| Signal | Estimates | Values | Side of null | Direction-consistent | Label |
|---|---|---|---|---|---|
| Domain overlap | 5 | G1 P-all (human) D1: 5.912; G3 P-annot D1: 0.679; P-abund (human) D1: 0.711; P-abund (human) D5: 0.652; P-abund (human) D4: 0.644 | above null, below null | no | population-dependent |
| Conservation | 11 | G1 P-all D1: 0.3993; G2 P-human D1: 0.4501; P-abund (human) D1: 0.4720; P-abund (human) D5: 0.4737; P-abund (human) D4: 0.4753; P-abund (mouse) D1: 0.4319; P-abund (mouse) D5: 0.4161; P-abund (mouse) D4: 0.4144; P-abund (rat) D1: 0.1394; P-abund (rat) D5: 0.1103; P-abund (rat) D4: 0.1235 | above null | yes | direction-consistent |
| Disorder | 11 | G1 P-all D1: 0.0457; G2 P-human D1: 0.0534; P-abund (human) D1: 0.1352; P-abund (human) D5: 0.1419; P-abund (human) D4: 0.1426; P-abund (mouse) D1: 0.0901; P-abund (mouse) D5: 0.0896; P-abund (mouse) D4: 0.0897; P-abund (rat) D1: 0.0320; P-abund (rat) D5: 0.0357; P-abund (rat) D4: 0.0595 | above null | yes | direction-consistent |
| Motif, sequence-based | 11 | G1 P-all D1: 0.0618; G2 P-human D1: 0.0368; P-abund (human) D1: 0.1041; P-abund (human) D5: 0.1109; P-abund (human) D4: 0.1120; P-abund (mouse) D1: 0.1784; P-abund (mouse) D5: 0.1844; P-abund (mouse) D4: 0.1848; P-abund (rat) D1: 0.2009; P-abund (rat) D5: 0.2208; P-abund (rat) D4: 0.2067 | above null | yes | direction-consistent |
| Motif, database lookup | 2 | G1 P-all D1: 0.7919; G2 P-human D1: 0.7964 | above null | yes | direction-consistent |

### Domain overlap by domain-annotation coverage of the protein (G20-0489, G20-0490)

An accession is annotated if at least one of its human site records has `in_domain` = 1. Crude, GOLD against BRONZE, set-B tiers, human records of P-all; gene-cluster intervals.

| Stratum | GOLD records | BRONZE records | In-domain share, GOLD | In-domain share, BRONZE | Crude OR (gene-cluster CI) |
|---|---|---|---|---|---|
| all proteins | 211,848 | 324,037 | 63.9% | 22.6% | 6.065 (5.751 to 6.429) |
| at least one domain-annotated site | 199,302 | 97,200 | 68.0% | 75.4% | 0.692 (0.657 to 0.730) |
| no domain-annotated site | 12,546 | 226,837 | 0.0% | 0.0% | undefined (outcome is constant (0) in this stratum, odds ratio undefined) |

**The crude odds ratio changes sign across populations:** 6.065 (all proteins), 0.741 (proteins with an abundance, P-abund), 0.692 (proteins with at least one annotated site; G3 gives the same value), 0.746 (canonical reviewed accessions, G16). Where the protein has any domain annotation, GOLD sites are *less* often in a domain than BRONZE sites (68.0% against 75.4%); where it has none the outcome is 0 for every record and no odds ratio exists. The pooled odds ratio above 1 comes from composition: 94.1% of GOLD records but 30.0% of BRONZE records lie in proteins with an annotated site (odds ratio 37.073 for being in such a protein, interval 34.010 to 40.699). **Domain overlap is therefore labelled population-dependent** (the other signals are direction-consistent, table above); the matched estimate for domain overlap is not a population-stable effect and should not be read as one.

## 4. Positive control and index diagnostics (block G20)

### 4a. Balance before matching

SMD of the primary index (V50) between GOLD and BRONZE records, before and after CEM, for every D1 and D3 row that was run (87 contrasts; decision-bearing flag from the manifest). Balance requirement |SMD| < 0.1 after matching.

| Contrast | Block | Signal | Population | Decision-bearing | SMD before | SMD after |
|---|---|---|---|---|---|---|
| G1-0002 | G1 | Domain overlap | P-all (human) | yes | -0.002 | -0.003 |
| G1-0005 | G1 | Conservation | P-all | yes | 0.160 | 0.004 |
| G1-0008 | G1 | Disorder | P-all | yes | 0.078 | 0.005 |
| G1-0011 | G1 | Motif, sequence-based | P-all | yes | -0.016 | 0.004 |
| G1-0014 | G1 | Motif, database lookup | P-all | no | -0.016 | 0.004 |
| G2-0017 | G2 | Conservation | P-human | yes | 0.187 | 0.007 |
| G2-0020 | G2 | Disorder | P-human | yes | 0.075 | 0.003 |
| G2-0023 | G2 | Motif, sequence-based | P-human | yes | -0.001 | -0.003 |
| G2-0026 | G2 | Motif, database lookup | P-human | no | -0.002 | -0.003 |
| G3-0029 | G3 | Domain overlap | P-annot | no | 0.151 | 0.005 |
| G12-0102 | G12 | Domain overlap | P-all | yes | 0.141 | 0.005 |
| G13-0104 | G13 | Domain overlap | P-all | yes | -0.002 | -0.003 |
| G13-0106 | G13 | Domain overlap | P-all | yes | -0.002 | -0.003 |
| G12-0178 | G12 | Conservation | P-all | yes | 0.270 | 0.009 |
| G13-0180 | G13 | Conservation | P-all | yes | 0.160 | 0.004 |
| G13-0182 | G13 | Conservation | P-all | yes | 0.160 | 0.004 |
| G12-0254 | G12 | Disorder | P-all | yes | 0.258 | 0.009 |
| G13-0256 | G13 | Disorder | P-all | yes | 0.078 | 0.005 |
| G13-0258 | G13 | Disorder | P-all | yes | 0.078 | 0.005 |
| G12-0330 | G12 | Motif, sequence-based | P-all | yes | 0.155 | 0.011 |
| G13-0332 | G13 | Motif, sequence-based | P-all | yes | -0.016 | 0.004 |
| G13-0334 | G13 | Motif, sequence-based | P-all | yes | -0.016 | 0.004 |
| G15-0351 | G15 | Domain overlap | human | no | 0.142 | 0.005 |
| G4-0031 | G4 | Domain overlap | P-all | yes | 0.070 | -0.000 |
| G4-0032 | G4 | Domain overlap | P-all | yes | -0.008 | -0.002 |
| G4-0033 | G4 | Domain overlap | P-all | yes | 0.005 | -0.002 |
| G4-0034 | G4 | Domain overlap | P-all | yes | 0.006 | 0.002 |
| G4-0035 | G4 | Domain overlap | P-all | yes | 0.006 | 0.002 |
| G4-0107 | G4 | Conservation | P-all | yes | 0.223 | 0.017 |
| G4-0108 | G4 | Conservation | P-all | yes | 0.156 | 0.004 |
| G4-0109 | G4 | Conservation | P-all | yes | 0.170 | 0.006 |
| G4-0110 | G4 | Conservation | P-all | yes | 0.124 | 0.003 |
| G4-0111 | G4 | Conservation | P-all | yes | 0.125 | 0.002 |
| G4-0183 | G4 | Disorder | P-all | yes | 0.136 | 0.013 |
| G4-0184 | G4 | Disorder | P-all | yes | 0.073 | 0.004 |
| G4-0185 | G4 | Disorder | P-all | yes | 0.085 | 0.005 |
| G4-0186 | G4 | Disorder | P-all | yes | 0.054 | -0.002 |
| G4-0187 | G4 | Disorder | P-all | yes | 0.054 | -0.002 |
| G4-0259 | G4 | Motif, sequence-based | P-all | yes | 0.013 | 0.010 |
| G4-0260 | G4 | Motif, sequence-based | P-all | yes | -0.019 | 0.004 |
| G4-0261 | G4 | Motif, sequence-based | P-all | yes | -0.012 | 0.005 |
| G4-0262 | G4 | Motif, sequence-based | P-all | yes | -0.013 | 0.001 |
| G4-0263 | G4 | Motif, sequence-based | P-all | yes | -0.013 | 0.001 |
| G17-0380 | G17 | Domain overlap | P-all (human) | yes | -0.002 | -0.003 |
| G17-0382 | G17 | Conservation | P-all | yes | 0.160 | 0.004 |
| G17-0384 | G17 | Disorder | P-all | yes | 0.078 | 0.005 |
| G17-0386 | G17 | Motif, sequence-based | P-all | yes | -0.016 | 0.004 |
| G17-0388 | G17 | Motif, database lookup | P-all | no | -0.016 | 0.004 |
| G17-0390 | G17 | Conservation | P-human | yes | 0.187 | 0.007 |
| G17-0392 | G17 | Disorder | P-human | yes | 0.075 | 0.003 |
| G17-0394 | G17 | Motif, sequence-based | P-human | yes | -0.001 | -0.003 |
| G17-0396 | G17 | Motif, database lookup | P-human | no | -0.002 | -0.003 |
| G17-0398 | G17 | Domain overlap | P-annot | no | 0.151 | 0.005 |
| G17-0400 | G17 | Domain overlap | P-abund (human) | yes | 0.246 | 0.010 |
| G17-0402 | G17 | Conservation | P-abund (human) | yes | 0.204 | 0.008 |
| G17-0404 | G17 | Disorder | P-abund (human) | yes | 0.235 | 0.010 |
| G17-0406 | G17 | Motif, sequence-based | P-abund (human) | yes | 0.247 | 0.010 |
| G16-0362 | G16 | Domain overlap | P-all (human), Swiss-Prot reviewed canonical | yes | 0.252 | 0.011 |
| G16-0364 | G16 | Conservation | P-all, Swiss-Prot reviewed canonical | yes | 0.164 | 0.006 |
| G16-0366 | G16 | Disorder | P-all, Swiss-Prot reviewed canonical | yes | 0.212 | 0.010 |
| G16-0368 | G16 | Motif, sequence-based | P-all, Swiss-Prot reviewed canonical | yes | 0.228 | 0.011 |
| G16-0370 | G16 | Motif, database lookup | P-all, Swiss-Prot reviewed canonical | no | 0.228 | 0.011 |
| G16-0372 | G16 | Conservation | P-human, Swiss-Prot reviewed canonical | yes | 0.205 | 0.008 |
| G16-0374 | G16 | Disorder | P-human, Swiss-Prot reviewed canonical | yes | 0.240 | 0.011 |
| G16-0376 | G16 | Motif, sequence-based | P-human, Swiss-Prot reviewed canonical | yes | 0.253 | 0.011 |
| G16-0378 | G16 | Motif, database lookup | P-human, Swiss-Prot reviewed canonical | no | 0.252 | 0.011 |
| G19-0419 | G19 | Domain overlap | P-human | yes | -0.002 | -0.003 |
| G19-0421 | G19 | Domain overlap | P-human | yes | 0.142 | 0.005 |
| G19-0423 | G19 | Domain overlap | P-human | yes | 0.141 | 0.005 |
| G19-0425 | G19 | Domain overlap | P-human | yes | -0.036 | -0.005 |
| G19-0427 | G19 | Conservation | P-human | yes | 0.187 | 0.007 |
| G19-0429 | G19 | Conservation | P-human | yes | 0.239 | 0.008 |
| G19-0431 | G19 | Conservation | P-human | yes | 0.237 | 0.008 |
| G19-0433 | G19 | Conservation | P-human | yes | 0.195 | 0.007 |
| G19-0435 | G19 | Disorder | P-human | yes | 0.075 | 0.003 |
| G19-0437 | G19 | Disorder | P-human | yes | 0.236 | 0.014 |
| G19-0439 | G19 | Disorder | P-human | yes | 0.235 | 0.014 |
| G19-0441 | G19 | Disorder | P-human | yes | 0.034 | -0.000 |
| G19-0443 | G19 | Motif, sequence-based | P-human | yes | -0.001 | -0.003 |
| G19-0445 | G19 | Motif, sequence-based | P-human | yes | 0.143 | 0.005 |
| G19-0447 | G19 | Motif, sequence-based | P-human | yes | 0.141 | 0.005 |
| G19-0449 | G19 | Motif, sequence-based | P-human | yes | -0.036 | -0.005 |
| G19-0451 | G19 | Conservation | P-all | yes | 0.160 | 0.004 |
| G19-0453 | G19 | Conservation | P-all | yes | 0.270 | 0.009 |
| G19-0455 | G19 | Conservation | P-all | yes | 0.172 | 0.000 |
| G19-0457 | G19 | Disorder | P-all | yes | 0.078 | 0.005 |
| G19-0459 | G19 | Disorder | P-all | yes | 0.258 | 0.009 |

Across the 78 decision-bearing contrasts the SMD before matching ranges from -0.036 to 0.270; 41 have |SMD| of 0.1 or more before matching and 0 after (largest after: 0.017). For the primary domain-overlap contrast (G1-0002) the index is already balanced before matching (-0.0015), so matching has nothing to adjust there.

### 4b. Crude association of the index with each outcome within a tier

Record level, Pearson correlation of the V50 index with the outcome (equal to the point-biserial correlation for binary outcomes), within GOLD records and within BRONZE records, gene-cluster 95% interval.

| Contrast | Signal | Population | Tier | Records | Gene clusters | Correlation (CI) | Excludes 0 |
|---|---|---|---|---|---|---|---|
| G20-0469 | Domain overlap | P-all (human) | GOLD | 211,848 | 15,021 | 0.067 (0.044 to 0.090) | yes |
| G20-0470 | Domain overlap | P-all (human) | BRONZE | 324,037 | 16,421 | -0.063 (-0.078 to -0.048) | yes |
| G20-0471 | Conservation | P-all | GOLD | 93,308 | 14,880 | -0.088 (-0.103 to -0.074) | yes |
| G20-0472 | Conservation | P-all | BRONZE | 37,591 | 11,289 | -0.081 (-0.106 to -0.058) | yes |
| G20-0473 | Conservation | P-human | GOLD | 49,939 | 6,744 | -0.012 (-0.031 to 0.009) | no |
| G20-0474 | Conservation | P-human | BRONZE | 17,477 | 5,553 | 0.010 (-0.022 to 0.046) | no |
| G20-0475 | Disorder | P-all | GOLD | 320,780 | 30,233 | -0.074 (-0.091 to -0.058) | yes |
| G20-0476 | Disorder | P-all | BRONZE | 215,745 | 27,094 | 0.006 (-0.015 to 0.026) | no |
| G20-0477 | Disorder | P-human | GOLD | 208,322 | 15,008 | -0.028 (-0.047 to -0.007) | yes |
| G20-0478 | Disorder | P-human | BRONZE | 129,438 | 14,922 | 0.046 (0.019 to 0.072) | yes |
| G20-0479 | Motif, sequence-based | P-all | GOLD | 324,907 | 30,327 | -0.073 (-0.078 to -0.068) | yes |
| G20-0480 | Motif, sequence-based | P-all | BRONZE | 478,837 | 31,279 | -0.023 (-0.029 to -0.017) | yes |
| G20-0481 | Motif, sequence-based | P-human | GOLD | 211,848 | 15,021 | -0.044 (-0.050 to -0.037) | yes |
| G20-0482 | Motif, sequence-based | P-human | BRONZE | 324,000 | 16,420 | -0.033 (-0.041 to -0.026) | yes |
| G20-0483 | Motif, database lookup | P-all | GOLD | 324,907 | 30,327 | -0.032 (-0.039 to -0.025) | yes |
| G20-0484 | Motif, database lookup | P-all | BRONZE | 478,880 | 31,282 | -0.033 (-0.041 to -0.026) | yes |
| G20-0485 | Motif, database lookup | P-human | GOLD | 211,848 | 15,021 | -0.018 (-0.028 to -0.010) | yes |
| G20-0486 | Motif, database lookup | P-human | BRONZE | 324,037 | 16,421 | 0.006 (-0.002 to 0.013) | no |

The correlations lie between -0.088 and 0.067. In 2 of 9 signal-by-population combinations the correlation excludes zero in both GOLD and BRONZE with opposite signs (Domain overlap P-all (human), Disorder P-human), which a single index-based stratification cannot represent.

### 4c. The index against the number of sites per protein

| Species set | Count | Accessions with an index | Gene clusters | Spearman rho (gene-cluster CI) |
|---|---|---|---|---|
| human | all site records per accession | 48,162 | 18,757 | 0.100 (0.087 to 0.113) |
| human | GOLD records per accession | 48,162 | 18,757 | 0.008 (-0.001 to 0.017) |
| all | all site records per accession | 87,338 | 45,545 | 0.155 (0.146 to 0.164) |
| all | GOLD records per accession | 87,338 | 45,545 | 0.043 (0.036 to 0.050) |

**Index against the number of LT_LIT-supported sites per gene (G20-0488): not run: source file not available (PhosphoSitePlus v7 site provides no bulk download; access requested from the provider).**

### 4d. Could matching on the index remove confounding at all?

Matching on a variable can remove confounding by that variable only to the extent that (i) the compared groups differ in it and (ii) it is related to the outcome within the groups. Under a linear approximation (omitted-variable formula) the shift it can cause in the GOLD-BRONZE difference, in units of the outcome SD, is about SMD x r, where r is the within-tier correlation. This is a rough guide to size, not a bias estimate; it is shown for the decision-bearing G1 and G2 contrasts:

| Signal | Population | SMD before matching | r within GOLD | r within BRONZE | SMD x mean r | Observed attenuation (gene) |
|---|---|---|---|---|---|---|
| Domain overlap | P-all (human) | -0.002 | 0.067 | -0.063 | -0.0000 | 1.4% |
| Conservation | P-all | 0.160 | -0.088 | -0.081 | -0.0135 | -0.2% |
| Disorder | P-all | 0.078 | -0.074 | 0.006 | -0.0026 | -4.5% |
| Motif, sequence-based | P-all | -0.016 | -0.073 | -0.023 | 0.0008 | -0.1% |
| Conservation | P-human | 0.187 | -0.012 | 0.010 | -0.0002 | 1.3% |
| Disorder | P-human | 0.075 | -0.028 | 0.046 | 0.0007 | 1.4% |
| Motif, sequence-based | P-human | -0.001 | -0.044 | -0.033 | 0.0001 | -11.1% |

In plain words: **matching on this index can remove only a small amount of confounding, and a near-zero attenuation does not show that study effort is not a confounder.** In the decision-bearing G1 and G2 contrasts the product is at most 0.0135 outcome SD, and the within-tier correlations of the index with the outcomes lie between -0.088 and 0.067. The index also tracks the number of sites per protein weakly (Spearman 0.100 among human proteins) and the number of GOLD sites per protein not at all (Spearman 0.008, interval -0.001 to 0.017). The index is a gene-level count of publications; the tier is a property of how a site was reported (study-cluster depth, database agreement). Matching on it does remove an imbalance that exists in many populations, but the imbalance carries little of the outcome, so the results of CEM and D2b show that **adjusting for this index leaves the signals unchanged**; they do not show that site-level study effort cannot explain them.

## 5. Matcher D2b: what the Rev 5 Methods describe (block G18, post-hoc)

1:1 without replacement, caliper 0.2 SD of the logit (over the records of the resample), serine, threonine and tyrosine records only, random visiting order, 20 order replicates per contrast (min, median, max of the point estimate; the median is the estimate), the logit model and the matching recomputed in every one of 1,000 gene-cluster resamples (one random order per resample). Rev 5 kept each treated record's 50 nearest controls while its Methods say 1:1 matching: this is recorded in `RECONCILIATION.md`, section 6, and the manuscript Methods must be rewritten to describe what was run (not edited in this cycle). A point estimate outside its own percentile interval is flagged and the contrast is not reported as decision-bearing.

| Contrast | Signal | Population | D2b estimate (95% CI) | Replicates min / median / max | Crude | Attenuation | Matched pairs | Share of GOLD matched | SMD | Point outside own CI | Survives | Decision-bearing | Pre-specified D2 counterpart |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| G18-0407 | Domain overlap | P-all (human) (S/T/Y only) | 5.758 (5.423 to 6.050) | 5.715 / 5.758 / 5.794 | 6.067 | 2.9% | 211,848 | 100.0% | 0.0001 | no | yes | yes | G1-0003: not fitted |
| G18-0408 | Conservation | P-all (S/T/Y only) | 0.3981 (0.3875 to 0.4086) | 0.3951 / 0.3981 / 0.4009 | 0.3985 | 0.1% | 37,576 | 40.3% | 0.0089 | no | yes | yes | G1-0006: 0.3846 (0.3719 to 0.3946) |
| G18-0409 | Disorder | P-all (S/T/Y only) | 0.0460 (0.0411 to 0.0506) | 0.0456 / 0.0460 / 0.0464 | 0.0438 | -5.1% | 213,732 | 66.6% | 0.0183 | no | yes | yes | G1-0009: not fitted |
| G18-0410 | Motif, sequence-based | P-all (S/T/Y only) | 0.0608 (0.0578 to 0.0639) | 0.0597 / 0.0608 / 0.0614 | 0.0617 | 1.5% | 316,649 | 97.5% | 0.0193 | no | yes | yes | G1-0012: 0.0450 (0.0408 to 0.0519) |
| G18-0411 | Motif, database lookup | P-all (S/T/Y only) | 0.7894 (0.7860 to 0.7928) | 0.7890 / 0.7894 / 0.7896 | 0.7934 | 0.5% | 316,649 | 97.5% | 0.0194 | no | yes | no | G1-0015: 0.7738 (0.7683 to 0.7769) |
| G18-0412 | Conservation | P-human (S/T/Y only) | 0.4481 (0.4340 to 0.4623) | 0.4443 / 0.4481 / 0.4550 | 0.4559 | 1.7% | 17,462 | 35.0% | 0.0060 | no | yes | yes | G2-0018: 0.4313 (0.4166 to 0.4457) |
| G18-0413 | Disorder | P-human (S/T/Y only) | 0.0554 (0.0480 to 0.0625) | 0.0546 / 0.0554 / 0.0564 | 0.0542 | -2.3% | 129,402 | 62.1% | 0.0180 | no | yes | yes | G2-0021: not fitted |
| G18-0414 | Motif, sequence-based | P-human (S/T/Y only) | 0.0376 (0.0344 to 0.0418) | 0.0362 / 0.0376 / 0.0390 | 0.0331 | -13.4% | 211,848 | 100.0% | 0.0001 | no | yes | yes | G2-0024: 0.0278 (0.0230 to 0.0365) |
| G18-0415 | Motif, database lookup | P-human (S/T/Y only) | 0.7957 (0.7912 to 0.7997) | 0.7951 / 0.7957 / 0.7961 | 0.7983 | 0.3% | 211,848 | 100.0% | 0.0001 | no | yes | no | G2-0027: not fitted |
| G18-0416 | Domain overlap | P-annot (S/T/Y only) | 0.687 (0.652 to 0.726) | 0.681 / 0.687 / 0.691 | 0.693 | -2.4% | 97,166 | 48.8% | 0.0200 | no | yes | no | G3-0030: 0.686 (0.642 to 0.737) |
| G18-0417 | Domain overlap | human (S/T/Y only) | 2.651 (2.510 to 2.812) | 2.639 / 2.651 / 2.671 | 2.495 | -6.6% | 153,738 | 100.0% | 0.0023 | no | yes | no | G15-0352: not fitted |

7 of 7 decision-bearing D2b contrasts survive; the largest attenuation among them is 2.9%; 0 of 11 point estimates lie outside their own interval. Order variability (max minus min over the 20 replicates) is at most 0.0107 on the difference scale and 0.079 on the odds-ratio scale, far narrower than the intervals. In 7 of 11 contrasts only part of GOLD is matched (share 35.0% to 97.5%); there the pairs number 66.1% to 100.0% of the BRONZE records (the number of controls and the caliper limit the matching under 1:1 without replacement), and the estimate refers to the matched part of GOLD. The Rev 5 matcher results (including its non-fits) remain in INTERIM_REPORT.md, section 5.

## 6. Tier-definition sensitivity (block G19, post-hoc)

Crude and CEM-adjusted (D1, V50) estimates under four tier definitions: set-B tiers (primary), set-A tiers (human records only), depth-only tiers (GOLD depth >= 10, SILVER 3 to 9, BRONZE <= 2) and an agreement-only tier (GOLD agreement >= 2 against BRONZE agreement 0; records with agreement 1 not compared). Gene-cluster intervals.

**P-human**

| Signal | Tier definition | GOLD | BRONZE | Crude (CI) | D1 adjusted (CI) | Attenuation | Survives |
|---|---|---|---|---|---|---|---|
| Domain overlap | set-B tiers | 211,848 | 324,037 | 6.065 (5.731 to 6.387) | 5.912 (5.588 to 6.222) | 1.4% | yes |
| Domain overlap | set-A tiers | 153,738 | 483,427 | 2.494 (2.346 to 2.644) | 2.646 (2.491 to 2.803) | -6.5% | yes |
| Domain overlap | depth-only tiers | 153,571 | 483,453 | 2.492 (2.343 to 2.641) | 2.644 (2.487 to 2.799) | -6.5% | yes |
| Domain overlap | agreement-only | 174,814 | 412,834 | 8.099 (7.620 to 8.546) | 7.581 (7.122 to 8.008) | 3.2% | yes |
| Conservation | set-B tiers | 49,939 | 17,477 | 0.4560 (0.4422 to 0.4692) | 0.4501 (0.4360 to 0.4645) | 1.3% | yes |
| Conservation | set-A tiers | 35,760 | 49,705 | 0.6241 (0.6142 to 0.6333) | 0.6162 (0.6064 to 0.6255) | 1.3% | yes |
| Conservation | depth-only tiers | 35,718 | 49,709 | 0.6245 (0.6152 to 0.6331) | 0.6165 (0.6066 to 0.6257) | 1.3% | yes |
| Conservation | agreement-only | 43,671 | 19,757 | 0.3978 (0.3803 to 0.4142) | 0.3920 (0.3756 to 0.4090) | 1.5% | yes |
| Disorder | set-B tiers | 208,322 | 129,438 | 0.0542 (0.0466 to 0.0610) | 0.0534 (0.0466 to 0.0606) | 1.4% | yes |
| Disorder | set-A tiers | 150,249 | 288,135 | 0.1614 (0.1557 to 0.1670) | 0.1507 (0.1451 to 0.1566) | 6.6% | yes |
| Disorder | depth-only tiers | 150,082 | 288,149 | 0.1616 (0.1555 to 0.1678) | 0.1508 (0.1448 to 0.1570) | 6.6% | yes |
| Disorder | agreement-only | 174,728 | 164,082 | 0.0052 (-0.0034 to 0.0141) | 0.0107 (0.0044 to 0.0191) | -104.2% | yes |
| Motif, sequence-based | set-B tiers | 211,848 | 324,000 | 0.0331 (0.0295 to 0.0368) | 0.0368 (0.0334 to 0.0403) | -11.1% | yes |
| Motif, sequence-based | set-A tiers | 153,738 | 483,390 | 0.1035 (0.0999 to 0.1073) | 0.1008 (0.0973 to 0.1045) | 2.6% | yes |
| Motif, sequence-based | depth-only tiers | 153,571 | 483,416 | 0.1035 (0.1000 to 0.1073) | 0.1008 (0.0974 to 0.1047) | 2.6% | yes |
| Motif, sequence-based | agreement-only | 174,814 | 412,796 | 0.0177 (0.0139 to 0.0215) | 0.0264 (0.0231 to 0.0299) | -49.0% | yes |

**P-all**

| Signal | Tier definition | GOLD | BRONZE | Crude (CI) | D1 adjusted (CI) | Attenuation | Survives |
|---|---|---|---|---|---|---|---|
| Conservation | set-B tiers | 93,308 | 37,591 | 0.3986 (0.3884 to 0.4094) | 0.3993 (0.3897 to 0.4097) | -0.2% | yes |
| Conservation | depth-only tiers | 66,155 | 113,047 | 0.4734 (0.4654 to 0.4815) | 0.5308 (0.5229 to 0.5385) | -12.1% | yes |
| Conservation | agreement-only | 81,550 | 43,466 | 0.3348 (0.3213 to 0.3486) | 0.3445 (0.3336 to 0.3573) | -2.9% | yes |
| Disorder | set-B tiers | 320,780 | 215,745 | 0.0438 (0.0384 to 0.0493) | 0.0457 (0.0412 to 0.0514) | -4.5% | yes |
| Disorder | depth-only tiers | 229,448 | 506,153 | 0.1393 (0.1351 to 0.1438) | 0.1369 (0.1327 to 0.1415) | 1.8% | yes |
| Disorder | agreement-only | 270,355 | 270,815 | 0.0018 (-0.0039 to 0.0070) | 0.0080 (0.0036 to 0.0138) | -350.0% | yes |
| Motif, sequence-based | set-B tiers | 324,907 | 478,837 | 0.0617 (0.0583 to 0.0652) | 0.0618 (0.0585 to 0.0650) | -0.1% | yes |
| Motif, sequence-based | depth-only tiers | 233,219 | 771,537 | 0.1314 (0.1281 to 0.1349) | 0.1270 (0.1237 to 0.1304) | 3.4% | yes |
| Motif, sequence-based | agreement-only | 270,840 | 604,860 | 0.0469 (0.0437 to 0.0505) | 0.0518 (0.0484 to 0.0549) | -10.4% | yes |

### Spread of the estimates across tier definitions, compared with the change caused by effort adjustment

Two comparisons, both on adjusted estimates from the same four (three for P-all) tier definitions. **(1) Spread, as ratio of the largest to the smallest adjusted estimate:** over all definitions, and excluding any estimate whose interval includes the null or has a limit within 0.01 of the null (the interval reading of "within 0.01 of zero"; the point-estimate reading is noted below the table). **(2) Like-for-like, on the log scale:** the largest relative change from crude to adjusted within one tier definition, max over definitions of |ln(adjusted / crude)|, against the spread of the adjusted estimates across definitions, ln(largest / smallest). Only estimates whose interval excludes the null enter (for the within-definition change both the crude and the adjusted estimate; for the spread the adjusted estimates). For the odds ratio the log scale is the log odds ratio.

| Population | Signal | Tier definitions | Spread, all adjusted estimates | Estimates excluded (CI includes the null or is within 0.01 of it) | Spread, excluding those | Largest within-definition change, max abs(ln(adjusted / crude)), natural log | same change as a ratio, exp(change) | ln spread of adjusted estimates across definitions, natural log | same spread as a ratio, exp(ln spread) | Spread larger than the within-definition change |
|---|---|---|---|---|---|---|---|---|---|---|
| P-human | Domain overlap | 4 | 2.87 | none excluded | 2.87 | 0.066 | 1.068 | 1.053 | 2.87 | yes |
| P-human | Conservation | 4 | 1.57 | none excluded | 1.57 | 0.015 | 1.015 | 0.453 | 1.57 | yes |
| P-human | Disorder | 4 | 14.08 | agreement-only excluded | 2.82 | 0.069 | 1.071 | 2.645 | 14.08 | yes |
| P-human | Motif, sequence-based | 4 | 3.81 | none excluded | 3.81 | 0.399 | 1.490 | 1.339 | 3.81 | yes |
| P-all | Conservation | 3 | 1.54 | none excluded | 1.54 | 0.114 | 1.121 | 0.432 | 1.54 | yes |
| P-all | Disorder | 3 | 17.11 | agreement-only excluded | 2.99 | 0.044 | 1.045 | 2.840 | 17.11 | yes |
| P-all | Motif, sequence-based | 3 | 2.45 | none excluded | 2.45 | 0.099 | 1.104 | 0.897 | 2.45 | yes |

Like-for-like, the log-scale spread across tier definitions exceeds the largest log-scale change caused by effort adjustment within a definition in 7 of the 7 signal-by-population combinations where both can be computed (natural-log scale; spread 0.432 to 2.840, as ratios 1.54 to 17.11; within-definition change 0.015 to 0.399, as ratios 1.015 to 1.490). The largest spread (2.840, a 17.11-fold ratio; disorder in P-all) rests on a near-zero agreement-only estimate; without it the spread is about 3-fold (2.99). The interval reading excludes the estimate if either limit of its interval lies within 0.01 of the null; this excludes, although its point estimate is at least 0.01 from the null (it would stay in under the point-estimate reading): Disorder P-human agreement-only (estimate 0.0107, interval 0.0044 to 0.0191); the point-estimate reading would exclude no further estimate.

**The 17-fold spread for disorder in P-all rests on one near-zero estimate:** the agreement-only adjusted estimate is 0.0080 (interval 0.0036 to 0.0138), and the crude agreement-only estimate 0.0018 has an interval that includes zero. Without it the spread of the disorder estimates in P-all is about 3-fold (2.99). The same holds for P-human (14.1-fold with the agreement-only estimate 0.0107, about 3-fold (2.82) without it).

### 6b. Depth-only and agreement-only tiers side by side (no new data; from G19)

The tier combines study-cluster depth and database agreement; the two tiers below each use one of them alone. Crude and CEM-adjusted estimates with gene-cluster intervals, same index and design as above. **This is a description of how the estimates differ between the two tier definitions; it carries no causal interpretation.**

| Population | Signal | Depth-only: GOLD / BRONZE | Depth-only crude (CI) | Depth-only adjusted (CI) | Agreement-only: GOLD / BRONZE | Agreement-only crude (CI) | Agreement-only adjusted (CI) |
|---|---|---|---|---|---|---|---|
| P-human | Domain overlap | 153,571 / 483,453 | 2.492 (2.343 to 2.641) | 2.644 (2.487 to 2.799) | 174,814 / 412,834 | 8.099 (7.620 to 8.546) | 7.581 (7.122 to 8.008) |
| P-human | Conservation | 35,718 / 49,709 | 0.6245 (0.6152 to 0.6331) | 0.6165 (0.6066 to 0.6257) | 43,671 / 19,757 | 0.3978 (0.3803 to 0.4142) | 0.3920 (0.3756 to 0.4090) |
| P-human | Disorder | 150,082 / 288,149 | 0.1616 (0.1555 to 0.1678) | 0.1508 (0.1448 to 0.1570) | 174,728 / 164,082 | 0.0052 (-0.0034 to 0.0141) | 0.0107 (0.0044 to 0.0191) |
| P-human | Motif, sequence-based | 153,571 / 483,416 | 0.1035 (0.1000 to 0.1073) | 0.1008 (0.0974 to 0.1047) | 174,814 / 412,796 | 0.0177 (0.0139 to 0.0215) | 0.0264 (0.0231 to 0.0299) |
| P-all | Conservation | 66,155 / 113,047 | 0.4734 (0.4654 to 0.4815) | 0.5308 (0.5229 to 0.5385) | 81,550 / 43,466 | 0.3348 (0.3213 to 0.3486) | 0.3445 (0.3336 to 0.3573) |
| P-all | Disorder | 229,448 / 506,153 | 0.1393 (0.1351 to 0.1438) | 0.1369 (0.1327 to 0.1415) | 270,355 / 270,815 | 0.0018 (-0.0039 to 0.0070) | 0.0080 (0.0036 to 0.0138) |
| P-all | Motif, sequence-based | 233,219 / 771,537 | 0.1314 (0.1281 to 0.1349) | 0.1270 (0.1237 to 0.1304) | 270,840 / 604,860 | 0.0469 (0.0437 to 0.0505) | 0.0518 (0.0484 to 0.0549) |

For disorder, the crude estimate is 0.1393 under the depth-only tier and 0.0018 (interval -0.0039 to 0.0070, includes zero) under the agreement-only tier in P-all; in P-human 0.1616 and 0.0052 (interval -0.0034 to 0.0141, includes zero). The adjusted estimates are 0.1369 and 0.0080 (P-all) and 0.1508 and 0.0107 (P-human).

## 7. Canonical accessions (block G16, post-hoc)

Restricted to Swiss-Prot reviewed accessions (release 2026_03), one per gene cluster (the accession with the most records where a cluster holds several; 65 clusters). **Coverage:** 40,583 accessions kept; human 19,612 of 75,248 accessions and 563,610 of 1,071,725 records (52.6%); mouse 50.1% of records; rat 40.4%. Gene-cluster intervals; the unrestricted gene-cluster result (G17) is shown for comparison.

| Contrast | Signal | Population | GOLD | BRONZE | Crude (CI) | D1 adjusted (CI) | Attenuation | Survives | Unrestricted D1 (gene CI, G17) | SMD before / after |
|---|---|---|---|---|---|---|---|---|---|---|
| G16-0362 | Domain overlap | P-all (human), Swiss-Prot reviewed canonical | 199,886 | 77,656 | 0.746 (0.707 to 0.793) | 0.714 (0.675 to 0.757) | -14.8% | yes | 5.912 (5.571 to 6.241) | 0.252 / 0.011 |
| G16-0364 | Conservation | P-all, Swiss-Prot reviewed canonical | 89,922 | 32,040 | 0.4489 (0.4404 to 0.4575) | 0.4364 (0.4280 to 0.4457) | 2.8% | yes | 0.3993 (0.3892 to 0.4099) | 0.164 / 0.006 |
| G16-0366 | Disorder | P-all, Swiss-Prot reviewed canonical | 300,067 | 119,893 | 0.1179 (0.1135 to 0.1220) | 0.1182 (0.1137 to 0.1226) | -0.2% | yes | 0.0457 (0.0411 to 0.0512) | 0.212 / 0.010 |
| G16-0368 | Motif, sequence-based | P-all, Swiss-Prot reviewed canonical | 300,448 | 123,125 | 0.1277 (0.1237 to 0.1318) | 0.1266 (0.1226 to 0.1309) | 0.9% | yes | 0.0618 (0.0587 to 0.0652) | 0.228 / 0.011 |
| G16-0370 | Motif, database lookup | P-all, Swiss-Prot reviewed canonical | 300,448 | 123,164 | 0.8336 (0.8301 to 0.8369) | 0.8340 (0.8304 to 0.8374) | -0.0% | yes | 0.7919 (0.7887 to 0.7951) | 0.228 / 0.011 |
| G16-0372 | Conservation | P-human, Swiss-Prot reviewed canonical | 49,799 | 16,797 | 0.4754 (0.4635 to 0.4867) | 0.4718 (0.4601 to 0.4840) | 0.8% | yes | 0.4501 (0.4349 to 0.4642) | 0.205 / 0.008 |
| G16-0374 | Disorder | P-human, Swiss-Prot reviewed canonical | 199,633 | 76,400 | 0.1362 (0.1301 to 0.1424) | 0.1351 (0.1287 to 0.1414) | 0.8% | yes | 0.0534 (0.0474 to 0.0611) | 0.240 / 0.011 |
| G16-0376 | Motif, sequence-based | P-human, Swiss-Prot reviewed canonical | 199,886 | 77,622 | 0.1007 (0.0963 to 0.1051) | 0.1044 (0.0996 to 0.1092) | -3.6% | yes | 0.0368 (0.0331 to 0.0404) | 0.253 / 0.011 |
| G16-0378 | Motif, database lookup | P-human, Swiss-Prot reviewed canonical | 199,886 | 77,656 | 0.8292 (0.8243 to 0.8333) | 0.8301 (0.8253 to 0.8345) | -0.1% | yes | 0.7964 (0.7918 to 0.8007) | 0.252 / 0.011 |

8 of the 8 contrasts other than domain overlap survive with a positive estimate (as unrestricted). **Domain overlap reverses**: the odds ratio is 0.714 (0.675 to 0.757) against 5.912 unrestricted. The reversal is consistent with the composition account of section 3 (GOLD records sit disproportionately in proteins with a domain-annotated site), but I have not tested that account on the canonical set separately, so it is stated as consistent with the result, not as its demonstrated cause.

## 8. Stratification and sensitivity (Task F: G6, G7, G8)

Gene-cluster bootstrap (primary) and accession-cluster bootstrap beside it, 1,000 resamples each, in `results/f_runs/`. A sensitivity analysis changes a decision if it moves a contrast across the survives/not or the largely-explained/not boundary relative to the primary design (G1 D1 for the same signal); the count of changes is reported per analysis, for both intervals. A change of sign relative to the primary is shown separately because the decision rule does not detect it.

### G6: residue types separately

P-S, P-T and P-Y are the P-all records restricted to serine, threonine and tyrosine.

| Analysis | Signal | Crude (gene CI) | D1 adjusted (gene CI) | D1 accession CI | Attenuation | Survives (gene) | Largely explained | GOLD / BRONZE | Decision vs G1 D1 (gene) | Sign vs G1 D1 |
|---|---|---|---|---|---|---|---|---|---|---|
| P-S | Domain overlap | 6.226 (5.864 to 6.611) | 6.218 (5.851 to 6.591) | 5.838 to 6.609 | 0.1% | yes | no | 137,274 / 215,518 | same | same sign |
| P-T | Domain overlap | 5.576 (5.243 to 5.942) | 5.504 (5.164 to 5.853) | 5.145 to 5.900 | 0.8% | yes | no | 49,247 / 81,649 | same | same sign |
| P-Y | Domain overlap | 6.821 (6.269 to 7.398) | 6.844 (6.309 to 7.353) | 6.358 to 7.362 | -0.2% | yes | no | 25,327 / 26,833 | same | same sign |
| P-S | Conservation | 0.4121 (0.4006 to 0.4244) | 0.4145 (0.4040 to 0.4258) | 0.4038 to 0.4262 | -0.6% | yes | no | 60,810 / 21,174 | same | same sign |
| P-T | Conservation | 0.3914 (0.3771 to 0.4063) | 0.4044 (0.3898 to 0.4182) | 0.3902 to 0.4179 | -3.3% | yes | no | 20,036 / 11,151 | same | same sign |
| P-Y | Conservation | 0.2874 (0.2688 to 0.3049) | 0.3126 (0.2962 to 0.3318) | 0.2958 to 0.3311 | -8.8% | yes | no | 12,462 / 5,251 | same | same sign |
| P-S | Disorder | 0.0442 (0.0391 to 0.0493) | 0.0462 (0.0414 to 0.0520) | 0.0420 to 0.0516 | -4.6% | yes | no | 216,980 / 134,167 | same | same sign |
| P-T | Disorder | 0.0474 (0.0410 to 0.0536) | 0.0533 (0.0472 to 0.0603) | 0.0477 to 0.0591 | -12.5% | yes | no | 70,548 / 60,562 | same | same sign |
| P-Y | Disorder | 0.0257 (0.0187 to 0.0333) | 0.0262 (0.0190 to 0.0346) | 0.0199 to 0.0328 | -1.7% | yes | no | 33,252 / 20,974 | same | same sign |
| P-S | Motif, sequence-based | 0.0779 (0.0741 to 0.0820) | 0.0759 (0.0720 to 0.0798) | 0.0720 to 0.0799 | 2.5% | yes | no | 220,476 / 322,061 | same | same sign |
| P-T | Motif, sequence-based | 0.0646 (0.0596 to 0.0699) | 0.0639 (0.0587 to 0.0690) | 0.0584 to 0.0695 | 1.0% | yes | no | 71,075 / 119,641 | same | same sign |
| P-Y | Motif, sequence-based | -0.0362 (-0.0438 to -0.0280) | -0.0353 (-0.0432 to -0.0270) | -0.0422 to -0.0283 | 2.6% | yes | no | 33,356 / 37,135 | same | opposite sign |

12 contrasts run; 0 change a decision under the gene-cluster interval and 0 under the accession-cluster interval; 1 have the opposite sign to the primary (Motif, sequence-based P-Y). Reported as found; I have not examined why.

### Sequence-based motif indicator: its four classes, and a serine/threonine-only contrast (post-hoc, G22)

The indicator (`scripts/epsd_phase0/phase1_6_sequence_motif.py`) counts a site as matching if its flanking window satisfies any one of four consensus classes. As written in that source:

| Class (label in the source) | Rule, positions relative to the phosphoacceptor | Pattern as written in the source | Written as a Ser/Thr motif? |
|---|---|---|---|
| Basophilic, PKA-like | R at -3 | [R]..{S,T} | yes (acceptor is {S,T}) |
| Basophilic, PKC-like | R or K at -2 | [RK].{S,T} | yes (acceptor is {S,T}) |
| Proline-directed, CDK/MAPK | P at +1 | {S,T}P | yes (acceptor is {S,T}) |
| Acidophilic, CK2-like | D or E at +3 | {S,T}..[DE] | yes (acceptor is {S,T}) |

All four are written as serine/threonine motifs and none is written for tyrosine. The code nevertheless evaluates them for any serine, threonine or tyrosine centre residue, so for tyrosine records the indicator applies a Ser/Thr consensus. The tyrosine contrast (G6, P-Y) is reported as it was run, and I have not examined why its sign differs.

| Population | GOLD / BRONZE | Crude (gene CI) | D1 adjusted (gene CI) | Attenuation | Status |
|---|---|---|---|---|---|
| P-all (G1; serine, threonine and tyrosine together) | 324,907 / 478,837 | 0.0617 (0.0584 to 0.0649) | 0.0618 (0.0587 to 0.0652) | -0.1% | primary, as reported |
| serine only (G6) | 220,476 / 322,061 | 0.0779 (0.0741 to 0.0820) | 0.0759 (0.0720 to 0.0798) | 2.5% | as reported |
| threonine only (G6) | 71,075 / 119,641 | 0.0646 (0.0596 to 0.0699) | 0.0639 (0.0587 to 0.0690) | 1.0% | as reported |
| tyrosine only (G6) | 33,356 / 37,135 | -0.0362 (-0.0438 to -0.0280) | -0.0353 (-0.0432 to -0.0270) | 2.6% | as reported |
| **serine and threonine together (G22, post-hoc)** | 291,551 / 441,702 | 0.0761 (0.0727 to 0.0797) | 0.0729 (0.0692 to 0.0765) | 4.3% | post-hoc, not decision-bearing; justified by the definition above |

The serine/threonine-only difference is 0.0729 (interval 0.0692 to 0.0765), against 0.0618 for P-all; the interval excludes zero (yes) and the attenuation by the index is 4.3%. The P-all and P-Y results stay as reported.

### G7: leave one study cluster (or the seven-PMID consortium) out

Removed PMIDs are deleted from every site's publication list, study-cluster depth and tier are recomputed, sites with no PMID and no database source flag are dropped. `drop_consortium_7_PMIDs` duplicates `drop_cluster_2` by construction (deviation log) and counts once in the tally below.

| Analysis | Signal | Crude (gene CI) | D1 adjusted (gene CI) | D1 accession CI | Attenuation | Survives (gene) | Largely explained | GOLD / BRONZE | Decision vs G1 D1 (gene) | Sign vs G1 D1 |
|---|---|---|---|---|---|---|---|---|---|---|
| drop_cluster_1 | Domain overlap | 6.002 (5.704 to 6.353) | 5.856 (5.564 to 6.197) | 5.510 to 6.180 | 1.4% | yes | no | 210,052 / 326,160 | same | same sign |
| drop_cluster_2 | Domain overlap | 4.423 (4.164 to 4.704) | 4.418 (4.151 to 4.678) | 4.152 to 4.680 | 0.1% | yes | no | 210,538 / 219,278 | same | same sign |
| drop_cluster_3 | Domain overlap | 6.381 (6.046 to 6.728) | 6.192 (5.864 to 6.519) | 5.833 to 6.558 | 1.6% | yes | no | 210,610 / 335,695 | same | same sign |
| drop_cluster_4 | Domain overlap | 6.065 (5.752 to 6.417) | 5.912 (5.602 to 6.234) | 5.557 to 6.274 | 1.4% | yes | no | 211,848 / 324,037 | same | same sign |
| drop_cluster_5 | Domain overlap | 6.091 (5.802 to 6.442) | 5.937 (5.652 to 6.274) | 5.564 to 6.269 | 1.4% | yes | no | 211,737 / 323,072 | same | same sign |
| drop_cluster_6 | Domain overlap | 6.074 (5.735 to 6.444) | 5.884 (5.552 to 6.223) | 5.554 to 6.222 | 1.8% | yes | no | 210,494 / 319,403 | same | same sign |
| drop_cluster_7 | Domain overlap | 6.065 (5.733 to 6.405) | 5.912 (5.591 to 6.231) | 5.558 to 6.256 | 1.4% | yes | no | 211,848 / 324,037 | same | same sign |
| drop_cluster_8 | Domain overlap | 6.089 (5.760 to 6.423) | 5.935 (5.621 to 6.251) | 5.591 to 6.294 | 1.4% | yes | no | 209,450 / 324,128 | same | same sign |
| drop_cluster_9 | Domain overlap | 6.065 (5.742 to 6.409) | 5.912 (5.596 to 6.228) | 5.580 to 6.251 | 1.4% | yes | no | 211,848 / 324,037 | same | same sign |
| drop_cluster_10 | Domain overlap | 6.108 (5.773 to 6.469) | 5.925 (5.606 to 6.272) | 5.586 to 6.304 | 1.7% | yes | no | 210,782 / 321,013 | same | same sign |
| drop_consortium_7_PMIDs | Domain overlap | 4.423 (4.168 to 4.713) | 4.418 (4.161 to 4.686) | 4.179 to 4.664 | 0.1% | yes | no | 210,538 / 219,278 | same | same sign |
| drop_cluster_1 | Conservation | 0.3980 (0.3875 to 0.4087) | 0.4009 (0.3918 to 0.4114) | 0.3916 to 0.4113 | -0.7% | yes | no | 92,346 / 38,339 | same | same sign |
| drop_cluster_2 | Conservation | 0.3950 (0.3837 to 0.4060) | 0.4001 (0.3894 to 0.4108) | 0.3905 to 0.4102 | -1.3% | yes | no | 93,226 / 36,410 | same | same sign |
| drop_cluster_3 | Conservation | 0.3986 (0.3876 to 0.4090) | 0.3993 (0.3889 to 0.4092) | 0.3894 to 0.4095 | -0.2% | yes | no | 93,209 / 37,853 | same | same sign |
| drop_cluster_4 | Conservation | 0.3980 (0.3875 to 0.4086) | 0.3963 (0.3861 to 0.4062) | 0.3858 to 0.4067 | 0.4% | yes | no | 93,015 / 36,164 | same | same sign |
| drop_cluster_5 | Conservation | 0.3987 (0.3867 to 0.4100) | 0.3968 (0.3861 to 0.4070) | 0.3864 to 0.4067 | 0.5% | yes | no | 93,184 / 36,413 | same | same sign |
| drop_cluster_6 | Conservation | 0.3981 (0.3873 to 0.4092) | 0.3991 (0.3893 to 0.4098) | 0.3899 to 0.4082 | -0.3% | yes | no | 93,242 / 37,475 | same | same sign |
| drop_cluster_7 | Conservation | 0.4006 (0.3889 to 0.4121) | 0.3999 (0.3897 to 0.4104) | 0.3902 to 0.4108 | 0.2% | yes | no | 93,181 / 37,007 | same | same sign |
| drop_cluster_8 | Conservation | 0.3986 (0.3879 to 0.4100) | 0.3985 (0.3889 to 0.4093) | 0.3881 to 0.4092 | 0.0% | yes | no | 92,887 / 37,625 | same | same sign |
| drop_cluster_9 | Conservation | 0.3983 (0.3867 to 0.4092) | 0.3991 (0.3887 to 0.4090) | 0.3895 to 0.4099 | -0.2% | yes | no | 93,111 / 37,538 | same | same sign |
| drop_cluster_10 | Conservation | 0.3988 (0.3876 to 0.4098) | 0.3995 (0.3888 to 0.4097) | 0.3902 to 0.4097 | -0.2% | yes | no | 93,200 / 37,406 | same | same sign |
| drop_consortium_7_PMIDs | Conservation | 0.3950 (0.3835 to 0.4058) | 0.4001 (0.3909 to 0.4104) | 0.3906 to 0.4097 | -1.3% | yes | no | 93,226 / 36,410 | same | same sign |
| drop_cluster_1 | Disorder | 0.0448 (0.0388 to 0.0501) | 0.0461 (0.0411 to 0.0518) | 0.0419 to 0.0508 | -3.0% | yes | no | 317,443 / 219,635 | same | same sign |
| drop_cluster_2 | Disorder | 0.0429 (0.0376 to 0.0481) | 0.0454 (0.0409 to 0.0511) | 0.0414 to 0.0502 | -6.0% | yes | no | 320,035 / 201,116 | same | same sign |
| drop_cluster_3 | Disorder | 0.0435 (0.0385 to 0.0489) | 0.0454 (0.0410 to 0.0509) | 0.0410 to 0.0500 | -4.3% | yes | no | 320,024 / 217,536 | same | same sign |
| drop_cluster_4 | Disorder | 0.0434 (0.0382 to 0.0483) | 0.0447 (0.0403 to 0.0498) | 0.0401 to 0.0500 | -3.0% | yes | no | 319,919 / 209,702 | same | same sign |
| drop_cluster_5 | Disorder | 0.0448 (0.0395 to 0.0503) | 0.0459 (0.0416 to 0.0516) | 0.0417 to 0.0507 | -2.6% | yes | no | 320,359 / 212,232 | same | same sign |
| drop_cluster_6 | Disorder | 0.0419 (0.0366 to 0.0470) | 0.0444 (0.0400 to 0.0497) | 0.0398 to 0.0493 | -5.8% | yes | no | 319,820 / 214,646 | same | same sign |
| drop_cluster_7 | Disorder | 0.0429 (0.0381 to 0.0482) | 0.0441 (0.0396 to 0.0495) | 0.0400 to 0.0490 | -2.8% | yes | no | 320,390 / 211,102 | same | same sign |
| drop_cluster_8 | Disorder | 0.0433 (0.0378 to 0.0485) | 0.0452 (0.0406 to 0.0506) | 0.0407 to 0.0497 | -4.3% | yes | no | 318,472 / 215,830 | same | same sign |
| drop_cluster_9 | Disorder | 0.0436 (0.0383 to 0.0487) | 0.0456 (0.0411 to 0.0510) | 0.0413 to 0.0500 | -4.5% | yes | no | 320,256 / 215,473 | same | same sign |
| drop_cluster_10 | Disorder | 0.0417 (0.0365 to 0.0471) | 0.0440 (0.0395 to 0.0493) | 0.0399 to 0.0487 | -5.5% | yes | no | 319,874 / 214,976 | same | same sign |
| drop_consortium_7_PMIDs | Disorder | 0.0429 (0.0372 to 0.0477) | 0.0454 (0.0405 to 0.0509) | 0.0414 to 0.0504 | -6.0% | yes | no | 320,035 / 201,116 | same | same sign |
| drop_cluster_1 | Motif, sequence-based | 0.0628 (0.0597 to 0.0659) | 0.0628 (0.0597 to 0.0660) | 0.0593 to 0.0661 | -0.0% | yes | no | 321,403 / 482,891 | same | same sign |
| drop_cluster_2 | Motif, sequence-based | 0.0755 (0.0717 to 0.0791) | 0.0743 (0.0706 to 0.0780) | 0.0710 to 0.0777 | 1.5% | yes | no | 323,597 / 374,078 | same | same sign |
| drop_cluster_3 | Motif, sequence-based | 0.0590 (0.0557 to 0.0621) | 0.0592 (0.0560 to 0.0623) | 0.0561 to 0.0624 | -0.5% | yes | no | 323,559 / 490,769 | same | same sign |
| drop_cluster_4 | Motif, sequence-based | 0.0561 (0.0528 to 0.0592) | 0.0564 (0.0533 to 0.0594) | 0.0532 to 0.0598 | -0.6% | yes | no | 323,999 / 468,503 | same | same sign |
| drop_cluster_5 | Motif, sequence-based | 0.0641 (0.0610 to 0.0673) | 0.0640 (0.0609 to 0.0673) | 0.0605 to 0.0674 | 0.2% | yes | no | 324,393 / 467,816 | same | same sign |
| drop_cluster_6 | Motif, sequence-based | 0.0588 (0.0554 to 0.0620) | 0.0594 (0.0560 to 0.0624) | 0.0563 to 0.0625 | -1.0% | yes | no | 323,553 / 474,203 | same | same sign |
| drop_cluster_7 | Motif, sequence-based | 0.0653 (0.0618 to 0.0683) | 0.0652 (0.0617 to 0.0681) | 0.0619 to 0.0685 | 0.1% | yes | no | 324,466 / 465,501 | same | same sign |
| drop_cluster_8 | Motif, sequence-based | 0.0632 (0.0600 to 0.0662) | 0.0632 (0.0601 to 0.0661) | 0.0600 to 0.0665 | -0.1% | yes | no | 322,509 / 478,928 | same | same sign |
| drop_cluster_9 | Motif, sequence-based | 0.0616 (0.0582 to 0.0646) | 0.0616 (0.0585 to 0.0644) | 0.0584 to 0.0650 | 0.1% | yes | no | 324,382 / 478,540 | same | same sign |
| drop_cluster_10 | Motif, sequence-based | 0.0589 (0.0553 to 0.0623) | 0.0594 (0.0559 to 0.0627) | 0.0561 to 0.0627 | -0.8% | yes | no | 323,641 / 476,039 | same | same sign |
| drop_consortium_7_PMIDs | Motif, sequence-based | 0.0755 (0.0717 to 0.0789) | 0.0743 (0.0708 to 0.0779) | 0.0708 to 0.0777 | 1.5% | yes | no | 323,597 / 374,078 | same | same sign |

44 D1 contrasts run (4 of them the duplicate removal); 0 change a decision under the gene-cluster interval and 0 under the accession-cluster interval; 0 have the opposite sign. Largest movement of the adjusted estimate relative to G1 D1, per signal: Domain overlap drop_cluster_2 (4.418 against 5.912, -25%); Conservation drop_cluster_4 (0.3963 against 0.3993, -1%); Disorder drop_cluster_10 (0.0440 against 0.0457, -4%); Motif, sequence-based drop_cluster_2 (0.0743 against 0.0618, +20%).

### G8: tier thresholds (GOLD depth in {5, 10, 20} x GOLD agreement in {>= 2, >= 3})

| Analysis | Signal | Crude (gene CI) | D1 adjusted (gene CI) | D1 accession CI | Attenuation | Survives (gene) | Largely explained | GOLD / BRONZE | Decision vs G1 D1 (gene) | Sign vs G1 D1 |
|---|---|---|---|---|---|---|---|---|---|---|
| gold_depth=5;gold_das=2 | Domain overlap | 5.257 (5.025 to 5.510) | 5.011 (4.800 to 5.247) | 4.791 to 5.239 | 2.9% | yes | no | 330,364 / 324,037 | same | same sign |
| gold_depth=5;gold_das=3 | Domain overlap | 5.179 (4.957 to 5.426) | 4.948 (4.733 to 5.180) | 4.711 to 5.177 | 2.8% | yes | no | 320,907 / 324,037 | same | same sign |
| gold_depth=10;gold_das=2 | Domain overlap | 6.065 (5.758 to 6.403) | 5.912 (5.619 to 6.236) | 5.564 to 6.245 | 1.4% | yes | no | 211,848 / 324,037 | same | same sign |
| gold_depth=10;gold_das=3 | Domain overlap | 5.224 (4.907 to 5.541) | 5.317 (4.996 to 5.641) | 4.979 to 5.677 | -1.1% | yes | no | 154,515 / 324,037 | same | same sign |
| gold_depth=20;gold_das=2 | Domain overlap | 6.709 (6.322 to 7.143) | 6.487 (6.098 to 6.889) | 6.157 to 6.909 | 1.8% | yes | no | 186,879 / 324,037 | same | same sign |
| gold_depth=20;gold_das=3 | Domain overlap | 5.169 (4.784 to 5.577) | 5.436 (5.033 to 5.849) | 5.037 to 5.880 | -3.1% | yes | no | 86,505 / 324,037 | same | same sign |
| gold_depth=5;gold_das=2 | Conservation | 0.3400 (0.3301 to 0.3509) | 0.3430 (0.3340 to 0.3527) | 0.3338 to 0.3523 | -0.9% | yes | no | 125,851 / 37,591 | same | same sign |
| gold_depth=5;gold_das=3 | Conservation | 0.3508 (0.3401 to 0.3610) | 0.3529 (0.3439 to 0.3628) | 0.3431 to 0.3617 | -0.6% | yes | no | 121,809 / 37,591 | same | same sign |
| gold_depth=10;gold_das=2 | Conservation | 0.3986 (0.3877 to 0.4091) | 0.3993 (0.3892 to 0.4091) | 0.3894 to 0.4095 | -0.2% | yes | no | 93,308 / 37,591 | same | same sign |
| gold_depth=10;gold_das=3 | Conservation | 0.5106 (0.5000 to 0.5214) | 0.5150 (0.5053 to 0.5255) | 0.5047 to 0.5252 | -0.8% | yes | no | 66,723 / 37,591 | same | same sign |
| gold_depth=20;gold_das=2 | Conservation | 0.3874 (0.3765 to 0.3975) | 0.3914 (0.3821 to 0.4018) | 0.3814 to 0.4017 | -1.0% | yes | no | 84,882 / 37,591 | same | same sign |
| gold_depth=20;gold_das=3 | Conservation | 0.5929 (0.5820 to 0.6041) | 0.6184 (0.6070 to 0.6302) | 0.6076 to 0.6294 | -4.3% | yes | no | 36,166 / 37,591 | same | same sign |
| gold_depth=5;gold_das=2 | Disorder | 0.0289 (0.0245 to 0.0334) | 0.0345 (0.0304 to 0.0395) | 0.0307 to 0.0386 | -19.6% | yes | no | 459,002 / 215,745 | same | same sign |
| gold_depth=5;gold_das=3 | Disorder | 0.0345 (0.0296 to 0.0392) | 0.0397 (0.0352 to 0.0447) | 0.0362 to 0.0437 | -15.2% | yes | no | 444,381 / 215,745 | same | same sign |
| gold_depth=10;gold_das=2 | Disorder | 0.0438 (0.0389 to 0.0486) | 0.0457 (0.0413 to 0.0508) | 0.0413 to 0.0502 | -4.5% | yes | no | 320,780 / 215,745 | same | same sign |
| gold_depth=10;gold_das=3 | Disorder | 0.0932 (0.0878 to 0.0987) | 0.0885 (0.0835 to 0.0940) | 0.0840 to 0.0936 | 5.1% | yes | no | 231,137 / 215,745 | same | same sign |
| gold_depth=20;gold_das=2 | Disorder | 0.0328 (0.0273 to 0.0383) | 0.0364 (0.0320 to 0.0418) | 0.0318 to 0.0415 | -11.0% | yes | no | 285,496 / 215,745 | same | same sign |
| gold_depth=20;gold_das=3 | Disorder | 0.1243 (0.1184 to 0.1303) | 0.1172 (0.1116 to 0.1239) | 0.1119 to 0.1233 | 5.7% | yes | no | 125,090 / 215,745 | same | same sign |
| gold_depth=5;gold_das=2 | Motif, sequence-based | 0.0434 (0.0403 to 0.0468) | 0.0440 (0.0409 to 0.0474) | 0.0413 to 0.0470 | -1.3% | yes | no | 488,413 / 478,837 | same | same sign |
| gold_depth=5;gold_das=3 | Motif, sequence-based | 0.0457 (0.0424 to 0.0491) | 0.0459 (0.0428 to 0.0493) | 0.0430 to 0.0490 | -0.6% | yes | no | 473,757 / 478,837 | same | same sign |
| gold_depth=10;gold_das=2 | Motif, sequence-based | 0.0617 (0.0584 to 0.0650) | 0.0618 (0.0584 to 0.0647) | 0.0586 to 0.0650 | -0.1% | yes | no | 324,907 / 478,837 | same | same sign |
| gold_depth=10;gold_das=3 | Motif, sequence-based | 0.0905 (0.0870 to 0.0939) | 0.0872 (0.0838 to 0.0905) | 0.0836 to 0.0910 | 3.7% | yes | no | 234,911 / 478,837 | same | same sign |
| gold_depth=20;gold_das=2 | Motif, sequence-based | 0.0698 (0.0665 to 0.0734) | 0.0707 (0.0675 to 0.0741) | 0.0673 to 0.0741 | -1.2% | yes | no | 286,367 / 478,837 | same | same sign |
| gold_depth=20;gold_das=3 | Motif, sequence-based | 0.1312 (0.1270 to 0.1351) | 0.1261 (0.1221 to 0.1299) | 0.1216 to 0.1306 | 3.9% | yes | no | 125,517 / 478,837 | same | same sign |

24 D1 contrasts run; 0 change a decision under the gene-cluster interval and 0 under the accession-cluster interval; 0 have the opposite sign. The primary rule (depth 10, agreement 2) reproduces G1 (`results/f_tiers_summary.json`: check passed).

## 9. Established confidence measures (Task H: G9)

### PeptideAtlas build 606: can site-level evidence be reconstructed? (instruction of 2026-10-04)

The build (`atlas_build_606.tsv.zip`, a phosphopeptide-enrichment build, `is_ptm_build` = Y) holds 30 tables (names and columns of all of them in `results/peptideatlas_606_inspection.json`). `modified_peptide_instance` has 3,097,865 rows, 2,908,370 with a `[Phospho]` token. What a site-level evidence category would need, and whether the dump has it (`scripts/analysis/peptideatlas_feasibility.py`; descriptive counts only, no category derived):

| Needed | Where it would come from | In this dump? |
|---|---|---|
| Site identity (accession, residue, position) | `modified_peptide_instance.modified_peptide_sequence` (residue carrying `[Phospho]`), `.peptide_instance_id` to `peptide_mapping` (`matched_biosequence_id`, `start_in_biosequence`) to `biosequence` (`biosequence_accession`, `biosequence_seq`) | Derivable. But 370,940 of 381,252 peptide instances map to two or more protein sequences (`peptide_instance.n_protein_mappings`; isoform and other-sequence handling not examined). |
| Site-localisation probability or score (PTMProphet, Ascore, FLR) | none | **No.** Across all tables the only columns that look like confidence are build-level `probability_threshold`, peptide-level `best_probability` and `best_adjusted_probability`, and protein-level `probability` and `confidence` (`columns_matching_localisation_or_site_confidence_pattern` in the feasibility file). |
| Unambiguous placement of the phosphate | derived: peptide instance with one placement per phospho count | **Not guaranteed.** 51.3% of the 384,770 (peptide, phospho count) groups appear with more than one placement; 12.7% of the 3,538,479 phospho tokens sit on residues other than S, T or Y (447,926 on alanine). The alanine tokens look like decoy placements, but nothing in the dump documents that; this is an inference, not verified. |
| Number of distinct peptides per site | count of distinct modified peptide instances per site | Derivable, but instances with competing placements would count towards several sites. |
| Number of PSMs per site | `modified_peptide_instance.n_observations` (no nulls), `n_adjusted_observations` | Exists. The dump has no column description, so that these count PSMs is the build's convention and is not verified here. |
| Probability | `modified_peptide_instance.best_probability` (no nulls), `best_adjusted_probability` | Peptide-level only, not a localisation probability; `best_adjusted_probability` is empty in 3,097,865 of 3,097,865 rows. |
| Independent samples and search batches | `sample_ids`, `search_batch_ids` | Exist. |

**Decision: not feasible from the available build.** The published definitions were read (Kalyuzhnyy et al., J Proteome Res 2022; University of Liverpool repository copy; verbatim text in `sources/SOURCES.md`, section A). For PeptideAtlas the categories are defined by counts of PSM observations of a site at PTMProphet probabilities: "The counts of observations with a probability of >0.95 were used as "positive" evidence for site phosphorylation. The counts at a probability threshold of ≤0.19 were used as "negative" evidence in favor of a site being a nonphosphosite." and "("High": ≥5 positive observations which accounted for ≥10% of total observations across all probabilities; "Medium": ≥5 positive observations, which accounted for <10% of total observations or 2−4 positive observations; "Low": 1 positive observation; "Not phosphorylated": 0 positive observations and ≥5 negative observations; "Other": site did not fall into any described categories)." Applying them needs a PTMProphet probability for every PSM of every site. The TSV dump of build 606 has neither a PTMProphet probability nor a per-site PSM table (feasibility table above), so the categories cannot be applied to it, and I have not built a count-based surrogate (instruction of 2026-10-04).

**What else I checked on the web (2026-10-04; URLs, quotes and outcomes in `sources/SOURCES.md`).** No per-site PeptideAtlas file is published on the pages I could read: the build page lists 6 files (FASTA, TSV and MySQL dumps, coordinate files), none per-site, and shows no licence terms. The build-details page shows that the server holds PTMProphet-based evidence (`nP.80-.95`, `nP.95-.99`, `nP.99-1`, `no-choice`; "nP.95-.99: Number of observed PSM with PTMProphet probability >= 0.95 and < 0.99"), but only as a summary per sequence database and residue (downloaded, 2,277 bytes) and as a per-protein table (13,884 proteins, PSM counts aggregated per protein); the bands (>=0.80, >=0.95, >=0.99) do not match the paper's ">0.95" and "<=0.19". Per-site values would mean scraping a web service, which I did not do. The authors' own site-level tables are listed in the paper's Supporting Information (Table S1, "Filtered PeptideAtlas build with human STY sites that have at least 1 associated PSM"; Table S5, "Summary of all STY sites in our analysis, their conservation data, proximal sites, phosphorylation likelihood, and structural data"); they describe the 2020 PeptideAtlas build, not build 606, and I could not obtain them (the publisher and figshare pages returned HTTP 403; each fetch was stopped and none was retried). The 2023 JPR paper on multiple evidence streams could be read only as an abstract (a pipeline from PSM to peptidoform-site level with a decoy-based false localisation rate, shown on rice data); whether it offers site-level categories or data usable here is not determined.

**Consequence.** The PeptideAtlas rows of G9 stay "not feasible from the available build". If you place Table S1 or S5 in `study/analysis/sources/`, a decision is needed first, because those files describe a different build (2020) from the one downloaded; nothing will be run on them without it.

### PhosphoSitePlus measures (LT-supported vs MS-only; category high (n >= 5) vs low (n = 1), n = LT_LIT + MS_LIT + MS_CST)

**Not run: source file not available (PhosphoSitePlus v7 site provides no bulk download; access requested from the provider).** `h_measures.py psp_lt psp_cat` would run both measures (D0, D1 and D4, gene-cluster and accession-cluster intervals) on the file, and stops if the columns differ from the plan's definition (D3). The set-B extract in `raw_data/` is not a substitute (RECONCILIATION.md section 7) and no surrogate was derived. **The meaning of `LT_LIT`, `MS_LIT` and `MS_CST` could not be confirmed from a primary page** (the phosphosite.org downloads page returned HTTP 403; the two PhosphoSitePlus papers I could read do not define the column names; `sources/SOURCES.md`, section C). The paper's PhosphoSitePlus rule is "The sum of observations across all studies was used to predict site phosphorylation likelihood in PSP ("High": ≥5 observations; "Medium": 2−4 observations; "Low": 1 observation)." It does not name the three columns, so whether n = LT_LIT + MS_LIT + MS_CST is that sum depends on what the download documents. **When the file arrives I will copy its header lines (licence, column notes), version stamp and any readme into `sources/SOURCES.md`, and I will not assume the meaning of the three columns beyond what the download documents; if it does not document them, Task H stops before it is run and I tell you** (instruction of 2026-10-04, added to decision D3's stop rule). The old `raw_data/harmonizome/phosphositeplus_gene_attribute_edges.txt.gz` is a gene-level attribute table, not this file, and was not used.

### UniProt reviewed MOD_RES, experimental vs inferred (optional item of section 9)

40,916 reviewed human phosphosite features (release 2026_03): 5,082 with experimental evidence (ECO:0000269), 10,372 with inferred evidence only (ECO:0000250, ECO:0000255, ECO:0000305), 25,462 neither (mass-spectrometry-only ECO:0007744, mixed or untagged; not compared). Matched to the site table on accession, residue and position: 14,185 records (5,023 experimental, 9,162 inferred) in 4,547 accessions and 4,542 gene clusters; 8.1% of the compared UniProt sites have no record in the site table. Treated group = experimental, control group = inferred; D1 = CEM on the V50 index decile with residue type exact; D4 adds the PaxDb abundance decile (records with an abundance); gene-cluster bootstrap primary.

**Caveat, stated before the results.** The inferred class includes sites annotated by sequence-similarity transfer (ECO:0000250, and the other inferred codes), i.e. partly *because* they are conserved in an ortholog, so conservation is partly circular in this contrast. **This contrast is not used as evidence for or against the main hypothesis.** It is reported because it was a pre-defined optional item, and it is not compared with the tier results.

| Signal | Experimental / inferred records | Crude (gene CI) | D1 adjusted (gene CI) | D1 accession CI | Attenuation | Survives (gene) | Largely explained | D4 adjusted (gene CI) | SMD after D1 |
|---|---|---|---|---|---|---|---|---|---|
| Domain overlap | 5,021 / 9,106 | 1.336 (1.138 to 1.583) | 1.026 (0.874 to 1.229) | 0.869 to 1.222 | 91.0% | no | yes | 1.068 (0.887 to 1.243) | 0.040 |
| Conservation | 1,467 / 2,586 | -0.1757 (-0.2166 to -0.1412) | -0.1297 (-0.1881 to -0.0873) | -0.1813 to -0.0863 | 26.2% | yes | no | -0.1561 (-0.2014 to -0.1098) | 0.017 |
| Disorder | 5,014 / 9,091 | -0.0637 (-0.0797 to -0.0483) | -0.0318 (-0.0489 to -0.0147) | -0.0494 to -0.0143 | 50.1% | yes | yes | -0.0346 (-0.0465 to -0.0123) | 0.041 |
| Motif, sequence-based | 5,021 / 9,106 | -0.0217 (-0.0446 to -0.0000) | 0.0327 (0.0087 to 0.0545) | 0.0076 to 0.0562 | 250.7% | yes | yes | 0.0383 (0.0134 to 0.0674) | 0.040 |

**Reading (reported as found).** Domain overlap: the crude difference (1.336) disappears after matching on the index (1.026, interval 0.874 to 1.229 includes 1; attenuation 91.0%), so this signal does not survive. Conservation is lower for experimental than for inferred sites (-0.1297). **This contrast is partly circular for conservation:** the inferred class is defined by evidence codes that include sequence-similarity transfer (ECO:0000250), i.e. sites annotated *because* they are conserved in an ortholog. Disorder is lower in the experimental class (-0.0318) with an attenuation of 50.1% (largely explained: yes). The sequence-based motif difference reverses sign after adjustment (crude -0.0217 with an interval ending at -0.00004, adjusted 0.0327); under the rule it counts as largely explained (attenuation 250.7%), but the crude interval excludes zero only marginally, so this row is weak evidence either way. The population is small and its composition differs from the tier populations.

## 10. Within-protein selection (Task I: G10)

**Not run: source file not available (PhosphoSitePlus v7 site provides no bulk download; access requested from the provider).** It needs the PhosphoSitePlus LT-supported and MS-only classes. What exists: the conditional-logit estimator (`cond_logit.py`), checked against `statsmodels` ConditionalLogit on simulated data (`check_cond_logit.py`; the script asserts agreement of the coefficients within 5e-04 and of the standard errors within 1e-04, a log-likelihood at least as high as statsmodels', and no overflow in a 3,000-record stratum; the assertions pass); the same estimator, with a parallel likelihood, was checked against `statsmodels` on a real subsample for G21 (`results/g21_cond_logit_check.json`). What does not exist yet: the runner that joins the PhosphoSitePlus classes to the site table, builds the peptide-detectability indicator (tryptic peptide containing the site, cleaved after K or R but not before P, 7 to 35 residues; protein sequences are on disk in `data_private/dl/seq_*.tsv`) and runs the gene-cluster bootstrap with and without that control, plus the within-gene linear model for disorder. It is to be written if the file is supplied, so that it is tested on the real columns.

## 10b. Within-gene comparison of GOLD and BRONZE records (block G21; amendment of 2026-10-04)

**Definition (written into the plan and committed before any G21 contrast was run, commit 8807afd; decided after the interim results and the main results table were seen, not pre-specified).** Records: the GOLD and BRONZE records (set-B tiers) of the G1 populations (P-all, all species; P-human) with serine, threonine or tyrosine residue, outcome and V50 index defined; domain overlap excluded. Stratum: gene cluster (gene fixed effects). Exposure: GOLD against BRONZE. Covariates: residue type (T and Y indicators). Binary outcomes: conditional logistic regression (odds ratio of GOLD against BRONZE within a gene); disorder: within-gene linear fixed-effect model (difference in means); beside each binary outcome, as a same-scale companion to the CEM differences, the within-gene linear fixed-effect difference in proportions. **M1 (all records of the population, no detectability control) is the primary model.** M1s (records for which a UniProt reviewed protein sequence is available) and M2 (M1s plus the peptide-detectability indicator) are sensitivity analyses on the reviewed-sequence subset (confirmed by the project lead, 2026-10-04). The indicator is: tryptic peptide containing the site has 7 to 35 residues; cleavage after K or R but not before P. (tryptic peptide containing the site has 7 to 35 residues; cleavage after K or R but not before P). Gene-cluster bootstrap, percentile 95% interval, all models refitted in every resample.

**Check before the run (required by the definition).** `cond_logit.py` against `statsmodels` ConditionalLogit on a real subsample (22,504 records in 3,000 gene clusters of P-all, conservation; largest stratum 214 records): largest coefficient difference 5.1e-05 (tolerance 5e-04), largest standard-error difference 1.4e-07 (tolerance 1e-04); the parallel evaluation used in the bootstrap gives identical coefficients to the serial one (largest difference 0.0e+00); the check passed. (`results/g21_cond_logit_check.json`.)

**P-all (S/T/Y)**

| Signal | Model | Records | Gene clusters with both tiers (GOLD / BRONZE records in them) | Estimator | Estimate (95% CI) | Persists within gene (CI excludes null) | Companion: within-gene linear difference (CI) | CEM D1 difference, G1/G2 (gene CI), for scale |
|---|---|---|---|---|---|---|---|---|
| Conservation | M1 | 130,884 | 8,390 (71,464 / 30,847) | odds ratio (conditional logit) | 9.144 (8.596 to 9.698) | yes | 0.3763 (0.3675 to 0.3846) | 0.3993 (0.3892 to 0.4099) |
| Conservation | M1s | 121,829 | 8,054 (69,418 / 27,401) | odds ratio (conditional logit) | 9.981 (9.405 to 10.561) | yes | 0.3943 (0.3854 to 0.4019) |  |
| Conservation | M2 | 121,829 | 8,054 (69,418 / 27,401) | odds ratio (conditional logit) | 9.512 (8.959 to 10.041) | yes | 0.3833 (0.3747 to 0.3908) |  |
| Disorder | M1 | 536,483 | 21,800 (287,865 / 198,523) | difference (within-gene FE) | 0.0453 (0.0431 to 0.0475) | yes | not applicable | 0.0457 (0.0411 to 0.0512) |
| Disorder | M1s | 419,142 | 20,536 (266,425 / 107,918) | difference (within-gene FE) | 0.0688 (0.0662 to 0.0714) | yes | not applicable |  |
| Disorder | M2 | 419,142 | 20,536 (266,425 / 107,918) | difference (within-gene FE) | 0.0713 (0.0686 to 0.0738) | yes | not applicable |  |
| Motif, sequence-based | M1 | 803,744 | 23,705 (302,590 / 450,332) | odds ratio (conditional logit) | 1.343 (1.326 to 1.359) | yes | 0.0692 (0.0665 to 0.0720) | 0.0618 (0.0587 to 0.0652) |
| Motif, sequence-based | M1s | 422,852 | 20,561 (266,819 / 108,570) | odds ratio (conditional logit) | 1.845 (1.809 to 1.882) | yes | 0.1435 (0.1389 to 0.1480) |  |
| Motif, sequence-based | M2 | 422,852 | 20,561 (266,819 / 108,570) | odds ratio (conditional logit) | 1.992 (1.953 to 2.033) | yes | 0.1591 (0.1547 to 0.1636) |  |
| Motif, database lookup (not decision-bearing) | M1 | 803,744 | 23,705 (302,590 / 450,332) | odds ratio (conditional logit) | 160.195 (152.282 to 167.980) | yes | 0.7872 (0.7831 to 0.7910) | 0.7919 (0.7887 to 0.7951) |
| Motif, database lookup (not decision-bearing) | M1s | 422,852 | 20,561 (266,819 / 108,570) | odds ratio (conditional logit) | 485.168 (436.454 to 536.950) | yes | 0.8157 (0.8109 to 0.8204) |  |
| Motif, database lookup (not decision-bearing) | M2 | 422,852 | 20,561 (266,819 / 108,570) | odds ratio (conditional logit) | 509.494 (457.686 to 563.220) | yes | 0.8180 (0.8132 to 0.8227) |  |

Records without a reviewed sequence: 379,629 of 803,744; with a sequence but a different residue at the stated position: 1,263; records in M1s and M2: 422,852. Records: 324,907 GOLD and 478,837 BRONZE in 37,901 gene clusters, 23,705 of which have both tiers.

**P-human (S/T/Y)**

| Signal | Model | Records | Gene clusters with both tiers (GOLD / BRONZE records in them) | Estimator | Estimate (95% CI) | Persists within gene (CI excludes null) | Companion: within-gene linear difference (CI) | CEM D1 difference, G1/G2 (gene CI), for scale |
|---|---|---|---|---|---|---|---|---|
| Conservation | M1 | 67,404 | 4,331 (40,361 / 15,054) | odds ratio (conditional logit) | 10.127 (9.425 to 10.904) | yes | 0.4070 (0.3955 to 0.4181) | 0.4501 (0.4349 to 0.4642) |
| Conservation | M1s | 66,499 | 4,300 (40,076 / 14,570) | odds ratio (conditional logit) | 10.678 (9.880 to 11.536) | yes | 0.4145 (0.4025 to 0.4254) |  |
| Conservation | M2 | 66,499 | 4,300 (40,076 / 14,570) | odds ratio (conditional logit) | 10.288 (9.520 to 11.132) | yes | 0.4059 (0.3941 to 0.4174) |  |
| Disorder | M1 | 337,724 | 12,784 (196,610 / 122,851) | difference (within-gene FE) | 0.0473 (0.0443 to 0.0502) | yes | not applicable | 0.0534 (0.0474 to 0.0611) |
| Disorder | M1s | 275,381 | 12,475 (184,462 / 69,842) | difference (within-gene FE) | 0.0766 (0.0732 to 0.0798) | yes | not applicable |  |
| Disorder | M2 | 275,381 | 12,475 (184,462 / 69,842) | difference (within-gene FE) | 0.0791 (0.0758 to 0.0824) | yes | not applicable |  |
| Motif, sequence-based | M1 | 535,848 | 13,801 (206,974 / 314,556) | odds ratio (conditional logit) | 1.226 (1.209 to 1.243) | yes | 0.0479 (0.0447 to 0.0511) | 0.0368 (0.0331 to 0.0404) |
| Motif, sequence-based | M1s | 276,943 | 12,487 (184,725 / 70,288) | odds ratio (conditional logit) | 1.635 (1.599 to 1.673) | yes | 0.1152 (0.1102 to 0.1204) |  |
| Motif, sequence-based | M2 | 276,943 | 12,487 (184,725 / 70,288) | odds ratio (conditional logit) | 1.765 (1.726 to 1.805) | yes | 0.1311 (0.1262 to 0.1360) |  |
| Motif, database lookup (not decision-bearing) | M1 | 535,848 | 13,801 (206,974 / 314,556) | odds ratio (conditional logit) | 192.318 (180.580 to 204.799) | yes | 0.7903 (0.7850 to 0.7955) | 0.7964 (0.7918 to 0.8007) |
| Motif, database lookup (not decision-bearing) | M1s | 276,943 | 12,487 (184,725 / 70,288) | odds ratio (conditional logit) | 492.682 (430.476 to 564.578) | yes | 0.8028 (0.7964 to 0.8091) |  |
| Motif, database lookup (not decision-bearing) | M2 | 276,943 | 12,487 (184,725 / 70,288) | odds ratio (conditional logit) | 518.874 (451.143 to 593.605) | yes | 0.8055 (0.7991 to 0.8117) |  |

Records without a reviewed sequence: 257,926 of 535,848; with a sequence but a different residue at the stated position: 979; records in M1s and M2: 276,943. Records: 211,848 GOLD and 324,000 BRONZE in 17,640 gene clusters, 13,801 of which have both tiers.

**Sequence coverage by tier** (reviewed UniProt sequences only; this decides who enters M1s and M2; `results/g21_sequence_coverage.json`):

| Signal | Population | GOLD records with a sequence | BRONZE records with a sequence | Detectable (7 to 35 residues) among GOLD with sequence | Detectable among BRONZE with sequence |
|---|---|---|---|---|---|
| Conservation | P-all (S/T/Y) | 96.3% | 85.2% | 83.9% | 75.6% |
| Disorder | P-all (S/T/Y) | 93.4% | 55.4% | 84.3% | 75.8% |
| Motif, sequence-based | P-all (S/T/Y) | 92.3% | 25.7% | 84.3% | 75.8% |
| Motif, database lookup | P-all (S/T/Y) | 92.3% | 25.7% | 84.3% | 75.8% |
| Conservation | P-human (S/T/Y) | 99.6% | 96.0% | 83.2% | 75.5% |
| Disorder | P-human (S/T/Y) | 95.6% | 58.8% | 83.8% | 75.6% |
| Motif, sequence-based | P-human (S/T/Y) | 94.2% | 23.9% | 83.8% | 75.6% |
| Motif, database lookup | P-human (S/T/Y) | 94.2% | 23.9% | 83.8% | 75.6% |

The sequence-restricted models M1s and M2 therefore rest on a subset that keeps most GOLD records but only 23.9% to 96.0% of BRONZE records, depending on the outcome; their estimates should be compared with M1 with that selection in mind, and the comparison of M2 with M1s (same records, with and without the detectability control) is the clean contrast for the control.

Of the 6 decision-bearing M1 contrasts, 6 persist within gene; of the 6 decision-bearing M2 contrasts (with the detectability control), 6 persist. All groups used 1,000 resamples (the decision-D7 fallback to 200 was not applied); 856 conditional-logit resample fits did not converge in total.

**In plain words:** a difference between GOLD and BRONZE sites that persists inside the same gene cannot be explained by anything that differs between genes, such as how heavily a gene is studied overall; **it excludes confounding at the protein (gene) level only.** It does not exclude confounding at the level of the site, such as how much attention a particular site received. Because sequences exist only for reviewed accessions, M1s and M2 describe a subset whose composition differs from M1 (see the shares of records with a sequence above and in `results/g21_sequence_coverage.json`).

## 10c. Non-primary publications counted in depth (block G23; sensitivity only; amendment of 2026-10-04, sixth batch)

**Why.** PMID 25514926 is the PhosphoSitePlus database paper (PubMed: "PhosphoSitePlus, 2014: mutations, PTMs and recalibrations", Nucleic Acids Research 2015; RECONCILIATION.md section 7) and is counted as a publication in 225,187 records. G23 asks how much the results depend on database, resource and review papers being counted as evidence. **The primary analysis is unchanged; G23 is a sensitivity analysis.** Rule, written into the plan and the list committed before any depth or contrast was computed (commits 92ab560 and 5650bed): candidate PMIDs = the 690 study-cluster candidates (linked to at least 884 records; PMID 25514926 is among them); PubMed metadata fetched by E-utilities; non-primary = PubMed publication type Review or Systematic Review, or a database or resource paper judged from the title alone (list with a reason for each entry in `nonprimary_pmid_list.md`; 16 borderline titles are listed there as not removed). Non-primary PMIDs are deleted from every publication list, the study-cluster assignment of the remaining PMIDs is kept, depth and tier are recomputed (DAS unchanged), and records left with no PMID and no database source flag are dropped.

### Descriptive

| PMID | Title (PubMed) | Records whose list contains it | Criterion |
|---|---|---|---|
| 30418626 | dbPTM in 2019: exploring disease association and cross-talk of post-translational modifications. | 448,403 | title |
| 29145615 | iPTMnet: an integrated resource for protein post-translational modification network discovery. | 421,351 | title |
| 22135298 | PhosphoSitePlus: a comprehensive resource for investigating the structure and function of experimentally determined post-translational modifications in man and mouse. | 245,532 | title |
| 25514926 | PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. | 225,187 | title |
| 37582339 | Proteogenomic data and resources for pan-cancer analysis. | 222,623 | Review type + title |
| 37619559 | A proteogenomics data-driven knowledge base of human cancer. | 222,623 | title |
| 15174125 | PhosphoSite: A bioinformatics resource dedicated to physiological protein phosphorylation. | 219,383 | title |
| 30445427 | 15 years of PhosphoSitePlus®: integrating post-translationally modified sites, disease variants and isoforms. | 219,383 | title |
| 23193284 | PTMcode: a database of known and predicted functional associations between post-translational modifications in proteins. | 127,502 | title |
| 32508104 | Scop3P: A Comprehensive Resource of Human Phosphosites within Their Full Context. | 93,664 | title |
| 16381945 | dbPTM: an information repository of protein post-translational modification. | 82,357 | title |
| 23193290 | DbPTM 3.0: an informative resource for investigating substrate site specificity and functional association of protein post-translational modifications. | 56,670 | title |
| 18988627 | Human Protein Reference Database--2009 update. | 55,853 | title |
| 21082442 | Phosphoproteome resource for systems biology research. | 52,111 | title |
| 30395287 | UniProt: a worldwide hub of protein knowledge. | 52,022 | title |
| 21081558 | PHOSIDA 2011: the posttranslational modification database. | 49,545 | title |
| 24705204 | SysPTM 2.0: an updated systematic resource for post-translational modification. | 40,010 | title |
| 20428315 | PhosphoGRID: a database of experimentally verified in vivo protein phosphorylation sites from the budding yeast Saccharomyces cerevisiae. | 35,520 | title |
| 30476227 | The BioGRID interaction database: 2019 update. | 35,520 | title |
| 21062810 | Phospho.ELM: a database of phosphorylation sites--update 2011. | 34,495 | title |
| 19060867 | PhosphoPep--a database of protein phosphorylation sites in model organisms. | 20,397 | title |
| 26425664 | Dataset from the global phosphoproteomic mapping of early mitotic exit in human cells. | 17,813 | title |
| 35608653 | In Vitro Kinase-to-Phosphosite Database (iKiP-DB) Predicts Kinase Activity in Phosphoproteomic Datasets. | 12,725 | title |
| 19366988 | SysPTM: a systematic resource for proteomic research on post-translational modifications. | 7,974 | title |
| 32618424 | LymphoAtlas: a dynamic and integrated phosphoproteomic resource of TCR signaling in primary T cells reveals ITSN2 as a regulator of effector functions. | 5,532 | title |
| 24174276 | Comprehensive proteomic data sets for studying adipocyte-macrophage cell-cell communication. | 5,343 | title |
| 30087585 | Residual tissue repositories as a resource for population-based cancer proteomic studies. | 4,120 | title |
| 15489334 | The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC). | 2,626 | title |

**28 non-primary PMIDs** (of 690 candidates; 1 by publication type, 28 by title; one is both). **776,218 of 1,768,068 records (43.9%) have at least one in their publication list** (by species: human 576,167, mouse 147,362, rat 52,689; by tier before correction: GOLD 327,155 (95.5% of that tier), SILVER 250,397 (33.4% of that tier), BRONZE 198,666 (29.4% of that tier)). 32,921 records have only non-primary PMIDs; 28,436 of them have no database source flag either and are dropped.

**Kept as primary although the title is not phosphorylation evidence (descriptive only; the G23 rule and every tier are unchanged):** 30243273 ("An efficient pair natural orbital based configuration interaction scheme for the calculation of open-shell ionization potentials.") is in the publication lists of 3,786 records (0.21%; tiers GOLD 3,685, SILVER 100, BRONZE 1; species human 3,786; only PMID in the list for 0); 37031204 ("Correction to: Phospho-heavy-labeled-spiketide FAIMS stepped-CV DDA (pHASED) provides real-time phosphoproteomics data to aid in cancer drug selection.") is in the publication lists of 3,893 records (0.22%; tiers GOLD 3,743, SILVER 103, BRONZE 47; species mouse 3,890, human 3; only PMID in the list for 37); 14702039 ("Complete sequencing and characterization of 21,243 full-length human cDNAs.") is in the publication lists of 1,864 records (0.11%; tiers GOLD 1,851, SILVER 13; species human 1,864; only PMID in the list for 0). Together 9,388 records (0.53%) carry at least one. The first has a title on a quantum-chemistry method, which suggests an erroneous PMID in the source; the second is a correction notice that double-counts its original paper; the third is a cDNA sequencing project.

### Tier transition under the corrected depth (old tier by new tier)

| Old tier \ new tier | GOLD | SILVER | BRONZE | dropped | Total |
|---|---|---|---|---|---|
| GOLD | 332,882 | 9,557 | 0 | 0 | 342,439 |
| SILVER | 0 | 733,367 | 12,294 | 5,141 | 750,802 |
| BRONZE | 0 | 0 | 651,532 | 23,295 | 674,827 |

GOLD falls from 342,439 to 332,882; 50,287 records change tier or are dropped (2.8%).

### G1 and G2 contrasts under the corrected tiers, beside the primary

Same design as the primary (CEM on the V50 decile, gene-cluster bootstrap, 1,000 resamples; D0 and D1 from the same resamples); only the tier column differs. Decision by the pre-specified rule: survives = the gene-cluster interval of the adjusted estimate excludes the null; sign = same side of the null as the primary estimate; largely explained = attenuation of 50% or more with a crude interval excluding the null.

| Block | Signal | Population | Primary D1 (gene CI; G17) | Corrected tiers D1 (gene CI) | GOLD / BRONZE (corrected) | Survives | Sign vs primary | Attenuation | Largely explained | Contrast |
|---|---|---|---|---|---|---|---|---|---|---|
| G1 | Domain overlap | P-all (human) | 5.912 (5.571 to 6.241) | 6.676 (6.302 to 7.082) | 207,256 / 318,105 | yes | same | 1.9% | no | G23-0518 (counterpart G1-0002) |
| G1 | Conservation | P-all | 0.3993 (0.3892 to 0.4099) | 0.3899 (0.3800 to 0.3997) | 91,504 / 35,222 | yes | same | -2.5% | no | G23-0520 (counterpart G1-0005) |
| G1 | Disorder | P-all | 0.0457 (0.0411 to 0.0512) | 0.0388 (0.0339 to 0.0445) | 313,881 / 204,181 | yes | same | -13.8% | no | G23-0522 (counterpart G1-0008) |
| G1 | Motif, sequence-based | P-all | 0.0618 (0.0587 to 0.0652) | 0.0582 (0.0548 to 0.0613) | 317,166 / 471,558 | yes | same | -2.1% | no | G23-0524 (counterpart G1-0011) |
| G1 | Motif, database lookup (not decision-bearing) | P-all | 0.7919 (0.7887 to 0.7951) | 0.7976 (0.7944 to 0.8007) | 317,166 / 471,601 | yes | same | 0.1% | no | G23-0526 (counterpart G1-0014) |
| G2 | Conservation | P-human | 0.4501 (0.4349 to 0.4642) | 0.4349 (0.4190 to 0.4495) | 49,198 / 15,527 | yes | same | 1.3% | no | G23-0528 (counterpart G2-0017) |
| G2 | Disorder | P-human | 0.0534 (0.0474 to 0.0611) | 0.0456 (0.0389 to 0.0533) | 204,532 / 119,319 | yes | same | -2.0% | no | G23-0530 (counterpart G2-0020) |
| G2 | Motif, sequence-based | P-human | 0.0368 (0.0331 to 0.0404) | 0.0322 (0.0285 to 0.0358) | 207,256 / 318,068 | yes | same | -17.4% | no | G23-0532 (counterpart G2-0023) |
| G2 | Motif, database lookup (not decision-bearing) | P-human | 0.7964 (0.7918 to 0.8007) | 0.7992 (0.7949 to 0.8030) | 207,256 / 318,105 | yes | same | 0.3% | no | G23-0534 (counterpart G2-0026) |

Under the corrected tiers 9 of 9 D1 contrasts survive, 9 of 9 have the same sign as the primary, and 0 are largely explained; the decisions are the same as in the primary analysis for every contrast. The estimates move modestly (for example domain overlap 5.912 to 6.676, disorder 0.0457 to 0.0388).

### GOLD counts in the held-out designs under the corrected tiers

| Design | Test sites | GOLD (primary tiers) | GOLD (corrected) | SILVER (corrected) | BRONZE (corrected) | Dropped |
|---|---|---|---|---|---|---|
| Cluster-disjoint test set | 361,908 | 37,959 | 37,277 | 171,206 | 150,276 | 3,149 |
| Cluster- and protein-disjoint test set | 29,727 | 537 | 525 | 11,153 | 17,556 | 493 |
| Temporal, all test sites | 536,037 | 28 | 23 | 36,521 | 485,589 | 13,904 |
| Temporal, protein-disjoint | 168,154 | 4 | 1 | 8,177 | 159,546 | 430 |
| Temporal, study-disjoint | 411 | 0 | 0 | 0 | 411 | 0 |
| Temporal, protein- and study-disjoint | 47 | 0 | 0 | 0 | 47 | 0 |
| All sites | 1,768,068 | 342,439 | 332,882 | 742,924 | 663,826 | 28,436 |

The Rev 5 split assignments (cluster-disjoint, protein-disjoint, temporal) are unchanged; only the tier column is replaced. The held-out designs still retain few GOLD sites (the protein-disjoint test set keeps 525 of 332,882).

## 11. Predictor gaps (Task J: G11, secondary analysis)

Top-three minus bottom-three index-decile AUROC gap on the frozen evaluation sample, planned design (DeepPhos and MusiteDeep on serine/threonine sites; PhosF3C on tyrosine sites). Equivalence margin +/-0.05 (90% interval inside the margin = gap absent); gene-cluster resampling within the decile pools (primary) with the Rev 5 site-level interval beside it. Restricted analyses keep positives first reported after the tool's restriction year and need at least 200 positives. No pooled ranking across tools. This analysis is secondary.

| Tool | Setting | Sites | Positives | AUROC bottom / top | Gap | Gene-cluster 95% CI | Gene-cluster 90% CI | Label (gene-cluster interval; site-level interval) |
|---|---|---|---|---|---|---|---|---|
| DeepPhos | unrestricted | 8,641 | 4,485 | 0.7277 / 0.7684 | 0.0406 | 0.0124 to 0.0673 | 0.0184 to 0.0639 | not detectable (interval includes the margin); site-level interval: not detectable (interval includes the margin) |
| DeepPhos | restriction_year=2019 | 5,383 | 1,227 | 0.7938 / 0.7924 | -0.0013 | -0.0386 to 0.0393 | -0.0322 to 0.0320 | gap absent (equivalent within +/-0.05); site-level interval: gap absent (equivalent within +/-0.05) |
| DeepPhos | restriction_year=2018 | 5,780 | 1,624 | 0.7872 / 0.8124 | 0.0252 | -0.0102 to 0.0600 | -0.0039 to 0.0536 | not detectable (interval includes the margin); site-level interval: not detectable (interval includes the margin) |
| MusiteDeep | unrestricted | 8,641 | 4,485 | 0.6955 / 0.7567 | 0.0611 | 0.0324 to 0.0894 | 0.0372 to 0.0858 | not detectable (interval includes the margin); site-level interval: not detectable (interval includes the margin) |
| MusiteDeep | restriction_year=2020 | 5,124 | 968 | 0.8153 / 0.8027 | -0.0125 | -0.0504 to 0.0265 | -0.0450 to 0.0215 | gap absent (equivalent within +/-0.05); site-level interval: gap absent (equivalent within +/-0.05) |
| MusiteDeep | restriction_year=2017 | 5,890 | 1,734 | 0.7634 / 0.8002 | 0.0368 | 0.0024 to 0.0714 | 0.0064 to 0.0678 | not detectable (interval includes the margin); site-level interval: not detectable (interval includes the margin) |
| PhosF3C | unrestricted | 1,359 | 515 | 0.6658 / 0.7199 | 0.0541 | -0.0221 to 0.1248 | -0.0063 to 0.1142 | not detectable (interval includes the margin); site-level interval: not detectable (interval includes the margin) |
| PhosF3C | restriction_year=2024 |  | 0 | not evaluated: too small to evaluate (n<200) |  |  |  |  |

**Training-set overlap (J4).** Only PhosF3C ships its training file; 293 of 515 tyrosine evaluation positives (56.9%) have a 31-residue window with at least 40% ungapped identity to a training positive. The DeepPhos and MusiteDeep repositories hold no training data, so the item stops for them. Details in `results/g11_predictors.json`.

The Rev 5 all-residue design is kept beside each planned row in `results/g11_predictors.json` (and the Rev 5 reported intervals in the `rev5_reported` entries); the gene-cluster label of the planned design equals that of the Rev 5 all-residue design in 7 of 7 settings.

**The criterion written down before any predictor was scored (commit 8e12390, 2026-09-24; recorded in `manuscript/CLAIMS_VS_EVIDENCE.md`) was a gap of at least 0.05 with an interval excluding zero in at least two tools. It is still not met:** only 1 tool-and-setting combination qualifies: MusiteDeep unrestricted (gap 0.0611; gene-cluster 95% interval 0.0324 to 0.0894), that is 1 of 3 tools (criterion: at least 2). This is the criterion set before the results; the 90% equivalence test and the labels above are the additional rules of this plan.

## 12. Departures, limitations and what I need from you

**Logged in the plan (section 15), none changes a result:** all entries of the interim stage; the amendments of 2026-10-01 (decisions D5 and D2, the gene-cluster correction, the added columns, G16 to G20, PeptideAtlas, Task H categories, Task I model); the operational definitions fixed while carrying out the amendment. Two further items arose in this stage and are logged with this report: (a) the review of the generated reports found a typed number in the interim builder (now computed); (b) the optional UniProt measure was run before the PhosphoSitePlus measures because it does not depend on the file. The fourth batch of 2026-10-04 added: the definitions read from the paper (PeptideAtlas rows stay not feasible), the definitions of G21 (within-gene GOLD versus BRONZE) and G22 (S/T-only motif), written into the plan and committed before they were run, and the report changes you asked for.

**Limitations to carry forward.** (1) The two-way interval rests on a few very large study clusters and is approximate. (2) The index is a gene-level count and is only weakly associated with the outcomes within a tier, so CEM on it cannot adjudicate site-level effort confounding (section 4). (3) Domain overlap is population-dependent and consistent with a composition effect (section 3), which was not tested further. (4) The D2 non-fits are results of the pre-specified design; D2b is post-hoc. (5) UniProt inferred evidence is partly conservation-derived. (6) Multiplicity is not corrected (plan section 2); the decision-bearing contrasts of G1 to G15 were fixed in the manifest before any contrast was run (commit f199f81, 2026-09-29), and those of the amendment blocks G16 to G20 were fixed in manifest rows committed before they were run (commit 0407605, 2026-10-01), but after the interim results had been seen. (7) The study-cluster threshold has no effect on any tier (G13). (8) In the predictor audit the unlabelled negatives on little-studied proteins may include undiscovered phosphosites, which lowers the AUROC there independently of the tool; the gap between well-studied and little-studied deciles is therefore not a clean measure of the tool's performance difference. (9) The conservation indicator requires the ortholog site to be reported as a phosphosite in the other species, so it depends on how much that species has been studied; a within-gene comparison removes differences between genes but not this dependence, because it is a property of the ortholog sites.

**Needed from you (rule 9):**

1. **PhosphoSitePlus.** Nothing is pending from you: Tasks H and I and G20-0488 are not run: source file not available (PhosphoSitePlus v7 site provides no bulk download; access requested from the provider) and are in Future work.
2. Only if you place Table S1 or S5 of the Kalyuzhnyy paper in `sources/`: a decision on whether files that describe the 2020 PeptideAtlas build may stand in for build 606. The manuscript, references, `legacy/` and `Figures.pptx` are untouched.

## 13. Status of the manifest

534 rows; every row of `grid_manifest.csv` is in `results/manifest_status.csv` with its status. G7's duplicate removal and G13's identical results count as run.

| Block   |   not feasible from the available build |   not fitted |   not run |   not run: source file not available (PhosphoSitePlus v7 site provides no bulk download; access requested from the provider) |   run |
|:--------|----------------------------------------:|-------------:|----------:|-----------------------------------------------------------------------------------------------------------------------------:|------:|
| G1      |                                       0 |            2 |         0 |                                                                                                                            0 |    13 |
| G2      |                                       0 |            2 |         0 |                                                                                                                            0 |    10 |
| G3      |                                       0 |            0 |         0 |                                                                                                                            0 |     3 |
| G4      |                                       0 |            0 |         0 |                                                                                                                            0 |    20 |
| G5      |                                       0 |            0 |         8 |                                                                                                                            0 |    40 |
| G6      |                                       0 |            0 |         0 |                                                                                                                            0 |    24 |
| G7      |                                       0 |            0 |         0 |                                                                                                                            0 |    88 |
| G8      |                                       0 |            0 |         0 |                                                                                                                            0 |    48 |
| G9      |                                      12 |            0 |         0 |                                                                                                                           24 |    12 |
| G10     |                                       0 |            0 |         0 |                                                                                                                            4 |     0 |
| G11     |                                       0 |            0 |         1 |                                                                                                                            0 |     7 |
| G12     |                                       0 |            0 |         0 |                                                                                                                            0 |     8 |
| G13     |                                       0 |            0 |         0 |                                                                                                                            0 |    16 |
| G14     |                                       0 |            1 |         0 |                                                                                                                            0 |    14 |
| G15     |                                       0 |            1 |         0 |                                                                                                                            0 |     2 |
| G16     |                                       0 |            0 |         0 |                                                                                                                            0 |    18 |
| G17     |                                       0 |            0 |         0 |                                                                                                                            0 |    28 |
| G18     |                                       0 |            0 |         0 |                                                                                                                            0 |    11 |
| G19     |                                       0 |            0 |         0 |                                                                                                                            0 |    50 |
| G20     |                                       0 |            0 |         0 |                                                                                                                            1 |    22 |
| G21     |                                       0 |            0 |         0 |                                                                                                                            0 |    24 |
| G22     |                                       0 |            0 |         0 |                                                                                                                            0 |     2 |
| G23     |                                       0 |            0 |         0 |                                                                                                                            0 |    18 |

**The 8 G5 rows marked "not run"** (G5-0040, G5-0041, G5-0042, G5-0043, G5-0044, G5-0045, G5-0046, G5-0047) are the domain overlap contrasts for P-abund (mouse), P-abund (rat) under D0, D1, D4, D5. Reason recorded in the results file: "domain overlap (S1) is defined for human records only" (the domain annotation exists only for the human records). The same signal was run for human (P-abund, human), and every other signal was run for mouse and rat.

## 14. Pre-set headline rule

Branch A if Task I shows a within-gene difference between LT-supported and MS-only sites in conservation or disorder (detectability-controlled model, gene-cluster interval excluding zero); Branch B otherwise, or if Task I is not run.

**Branch B applies**, because Task I has not been run (not run: source file not available (PhosphoSitePlus v7 site provides no bulk download; access requested from the provider)). If Task I is run later, the branch is re-evaluated by `build_claims.py` with the same rule.

## 15. Future work (ideas outside the frozen scope; not part of any analysis above)

- Tasks H and I and G20-0488 (PhosphoSitePlus evidence measures; within-gene comparison of LT-supported and MS-only sites, with and without the peptide-detectability control; index against LT_LIT-supported sites per gene), to be run if the provider supplies `Phosphorylation_site_dataset.gz` (possibly during revision); the file's own header and documentation of `LT_LIT`, `MS_LIT` and `MS_CST` are to be recorded first, and the plan's stop rule applies if they are not documented.
- Adjusting for study effort at the level of the site (for example the number of publications per site that is not part of the tier definition), which the gene-level index cannot do.
- A conservation indicator that accounts for how much the other species has been studied.
- A detectability control for every accession (sequences for unreviewed accessions), so that the sequence-restricted models are not limited to reviewed proteins.
- A site-level PeptideAtlas evidence measure from per-site PTMProphet probabilities, if a per-site file for the current build becomes available (or a decision is taken on using the authors' 2020-build tables).
- Replication of the tier comparisons in a second substrate, which the single-source subsets did not allow.
- Reprocessing public phosphoproteomics data at the PSM level with a decoy-based false localisation rate, so that tiers rest on localisation confidence.
