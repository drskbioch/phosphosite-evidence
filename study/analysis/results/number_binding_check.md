# Binding of the numbers in the claim rows of CLAIMS_VS_EVIDENCE.md to named result fields

39 claim rows (sections 1 to 6), columns Claim, Estimate, 95% CI, n and Source: 331 numbers. **Bound to a named result field (file, key path; the builder recorded the field when it formatted the number and this script re-read the file and recomputed it): 300. Derived (recomputed in check_numbers_sourced.py): 31. Generic match only: 0. Unsourced: 0.**

## Generic match only (equal to a source value somewhere, but not bound to a named field for this row)

none

## Unsourced

none

## Bindings

| Row | Column | Number | Bound to (file : key path) |
|---|---|---|---|
| C1.1 | Claim | 1,768,068 | earlier_results/dataset_counts.json : n_sites |
| C1.1 | Claim | 136,468 | earlier_results/dataset_counts.json : n_proteins |
| C1.1 | Estimate | 1,768,068 | earlier_results/dataset_counts.json : n_sites |
| C1.1 | Estimate | 136,468 | earlier_results/dataset_counts.json : n_proteins |
| C1.1 | n | 1,768,068 | earlier_results/dataset_counts.json : n_sites |
| C1.2 | Estimate | 342,439 | earlier_results/dataset_counts.json : tier_counts_mammalian.GOLD |
| C1.2 | Estimate | 750,802 | earlier_results/dataset_counts.json : tier_counts_mammalian.SILVER |
| C1.2 | Estimate | 674,827 | earlier_results/dataset_counts.json : tier_counts_mammalian.BRONZE |
| C1.2 | n | 1,768,068 | derived: sum of the mammalian tier counts |
| C1.3 | Estimate | 48.4 | derived: human GOLD share under the raw PMID rule (percent) |
| C1.3 | Estimate | 14.7 | derived: human GOLD share under the study-cluster rule (percent) |
| C1.3 | n | 1,071,725 | earlier_results/dataset_counts.json : human_n_sites |
| C1.4 | Claim | 2 | analysis/ANALYSIS_PLAN.md : DAS >= 2 (a definition or threshold recorded in the plan, not a result) |
| C1.4 | Estimate | 0.16 | earlier_results/dataset_counts.json : human_das_ge2_share_source_set_a (x 100) |
| C1.4 | Estimate | 16.4 | earlier_results/dataset_counts.json : human_das_ge2_share_source_set_b (x 100) |
| C1.4 | n | 1,071,725 | earlier_results/dataset_counts.json : human_n_sites |
| C1.5 | Claim | 500 | analysis/results/b_cluster_rule_info.json : (key name) |
| C1.5 | Claim | 884 | analysis/results/b_cluster_rule_info.json : 884.tier_counts_study_cluster_rule.GOLD |
| C1.5 | Claim | 1,000 | analysis/results/b_cluster_rule_info.json : (key name) |
| C1.5 | Estimate | 500 | analysis/results/b_cluster_rule_info.json : (key name) |
| C1.5 | Estimate | 884 | analysis/results/b_cluster_rule_info.json : 884.tier_counts_study_cluster_rule.GOLD |
| C1.5 | Estimate | 1,000 | analysis/results/b_cluster_rule_info.json : (key name) |
| C1.5 | Estimate | 342,439 | analysis/results/b_cluster_rule_info.json : 884.tier_counts_study_cluster_rule.GOLD |
| C1.5 | Estimate | 750,802 | analysis/results/b_cluster_rule_info.json : 884.tier_counts_study_cluster_rule.SILVER |
| C1.5 | Estimate | 674,827 | analysis/results/b_cluster_rule_info.json : 884.tier_counts_study_cluster_rule.BRONZE |
| C1.5 | n | 1,768,068 | earlier_results/dataset_counts.json : n_sites |
| C1.6 | Claim | 28 | analysis/results/g23_tiers_summary.json : n_nonprimary_pmids |
| C1.6 | Claim | 776,218 | analysis/results/g23_tiers_summary.json : records_with_at_least_one |
| C1.6 | Claim | 43.9 | analysis/results/g23_tiers_summary.json : share_records_with_at_least_one (x 100) |
| C1.6 | Claim | 225,187 | analysis/results/g23_tiers_summary.json : nonprimary.15.records_linked |
| C1.6 | Claim | 9,557 | derived: GOLD records that fall to SILVER or BRONZE when non-primary PMIDs are not counted |
| C1.6 | Claim | 3,786 | analysis/results/g23_kept_nonevidence_pmids.json : 30243273.records |
| C1.6 | Claim | 0.21 | analysis/results/g23_kept_nonevidence_pmids.json : 30243273.share_of_records (x 100) |
| C1.6 | Claim | 3,893 | analysis/results/g23_kept_nonevidence_pmids.json : 37031204.records |
| C1.6 | Claim | 0.22 | analysis/results/g23_kept_nonevidence_pmids.json : 37031204.share_of_records (x 100) |
| C1.6 | Claim | 1,864 | analysis/results/g23_kept_nonevidence_pmids.json : 14702039.records |
| C1.6 | Claim | 0.11 | analysis/results/g23_kept_nonevidence_pmids.json : 14702039.share_of_records (x 100) |
| C1.6 | Claim | 9,388 | analysis/results/g23_kept_nonevidence_pmids.json : any_of_the_three.records |
| C1.6 | Claim | 0.53 | analysis/results/g23_kept_nonevidence_pmids.json : any_of_the_three.share_of_records (x 100) |
| C1.6 | Estimate | 28 | analysis/results/g23_tiers_summary.json : n_nonprimary_pmids |
| C1.6 | Estimate | 776,218 | analysis/results/g23_tiers_summary.json : records_with_at_least_one |
| C1.6 | Estimate | 43.9 | analysis/results/g23_tiers_summary.json : share_records_with_at_least_one (x 100) |
| C1.6 | n | 1,768,068 | analysis/results/g23_tiers_summary.json : n_records |
| C1.7 | Estimate | 537 | analysis/results/g23_tiers_summary.json : held_out.cluster_and_protein_disjoint_test.old.GOLD |
| C1.7 | Estimate | 525 | analysis/results/g23_tiers_summary.json : held_out.cluster_and_protein_disjoint_test.new.GOLD |
| C1.7 | Estimate | 28 | analysis/results/g23_tiers_summary.json : held_out.temporal_test_all.old.GOLD |
| C1.7 | Estimate | 23 | analysis/results/g23_tiers_summary.json : held_out.temporal_test_all.new.GOLD |
| C1.7 | Estimate | 4 | analysis/results/g23_tiers_summary.json : held_out.temporal_protein_disjoint.old.GOLD |
| C1.7 | Estimate | 1 | analysis/results/g23_tiers_summary.json : held_out.temporal_protein_disjoint.new.GOLD |
| C2.1 | Estimate | 0.3993 | analysis/results/amend_results.json : G17-0382.estimate |
| C2.1 | 95% CI | 0.3892 | analysis/results/amend_results.json : G17-0382.ci.0 |
| C2.1 | 95% CI | 0.4099 | analysis/results/amend_results.json : G17-0382.ci.1 |
| C2.1 | n | 93,308 | analysis/results/amend_results.json : G17-0382.n_sites_gold |
| C2.1 | n | 37,591 | analysis/results/amend_results.json : G17-0382.n_sites_bronze |
| C2.1 | n | 17,777 | analysis/results/amend_results.json : G17-0382.n_gene_clusters |
| C2.2 | Estimate | 0.0457 | analysis/results/amend_results.json : G17-0384.estimate |
| C2.2 | 95% CI | 0.0411 | analysis/results/amend_results.json : G17-0384.ci.0 |
| C2.2 | 95% CI | 0.0512 | analysis/results/amend_results.json : G17-0384.ci.1 |
| C2.2 | n | 320,780 | analysis/results/amend_results.json : G17-0384.n_sites_gold |
| C2.2 | n | 215,745 | analysis/results/amend_results.json : G17-0384.n_sites_bronze |
| C2.2 | n | 35,524 | analysis/results/amend_results.json : G17-0384.n_gene_clusters |
| C2.3 | Estimate | 0.0618 | analysis/results/amend_results.json : G17-0386.estimate |
| C2.3 | 95% CI | 0.0587 | analysis/results/amend_results.json : G17-0386.ci.0 |
| C2.3 | 95% CI | 0.0652 | analysis/results/amend_results.json : G17-0386.ci.1 |
| C2.3 | n | 324,907 | analysis/results/amend_results.json : G17-0386.n_sites_gold |
| C2.3 | n | 478,837 | analysis/results/amend_results.json : G17-0386.n_sites_bronze |
| C2.3 | n | 37,901 | analysis/results/amend_results.json : G17-0386.n_gene_clusters |
| C2.4 | Claim | 1 | analysis/results/cem_results.json : G1-0002.null |
| C2.4 | Estimate | 5.912 | analysis/results/amend_results.json : G17-0380.estimate |
| C2.4 | 95% CI | 5.571 | analysis/results/amend_results.json : G17-0380.ci.0 |
| C2.4 | 95% CI | 6.241 | analysis/results/amend_results.json : G17-0380.ci.1 |
| C2.4 | n | 211,848 | analysis/results/amend_results.json : G17-0380.n_sites_gold |
| C2.4 | n | 324,037 | analysis/results/amend_results.json : G17-0380.n_sites_bronze |
| C2.4 | n | 17,641 | analysis/results/amend_results.json : G17-0380.n_gene_clusters |
| C2.5 | Estimate | 0.7919 | analysis/results/amend_results.json : G17-0388.estimate |
| C2.5 | 95% CI | 0.7887 | analysis/results/amend_results.json : G17-0388.ci.0 |
| C2.5 | 95% CI | 0.7951 | analysis/results/amend_results.json : G17-0388.ci.1 |
| C2.5 | n | 324,907 | analysis/results/amend_results.json : G17-0388.n_sites_gold |
| C2.5 | n | 478,880 | analysis/results/amend_results.json : G17-0388.n_sites_bronze |
| C2.5 | n | 37,903 | analysis/results/amend_results.json : G17-0388.n_gene_clusters |
| C2.6 | Claim | 1 | analysis/results/amend_results.json : G16-0361.null |
| C2.6 | Claim | 1 | analysis/results/amend_results.json : G16-0361.null |
| C2.6 | Estimate | 6.065 | analysis/results/amend_results.json : G17-0379.estimate |
| C2.6 | Estimate | 0.741 | analysis/results/amend_results.json : G17-0399.estimate |
| C2.6 | Estimate | 0.692 | analysis/results/amend_results.json : G17-0397.estimate |
| C2.6 | Estimate | 0.746 | analysis/results/amend_results.json : G16-0361.estimate |
| C3.1 | Estimate | -11.1 | analysis/results/amend_results.json : G17-0394.attenuation (x 100) |
| C3.1 | Estimate | 1.4 | analysis/results/amend_results.json : G17-0392.attenuation (x 100) |
| C3.2 | Claim | 0.088 | analysis/results/g20_diagnostics.json : index_outcome_association.G20-0471.pearson_r |
| C3.2 | Estimate | 0.088 | analysis/results/g20_diagnostics.json : index_outcome_association.G20-0471.pearson_r |
| C3.3 | Estimate | 0.008 | analysis/results/g20_diagnostics.json : index_vs_sites_per_accession.human|gold_sites_per_accession.spearman |
| C3.3 | 95% CI | -0.001 | analysis/results/g20_diagnostics.json : index_vs_sites_per_accession.human|gold_sites_per_accession.ci.0 |
| C3.3 | 95% CI | 0.017 | analysis/results/g20_diagnostics.json : index_vs_sites_per_accession.human|gold_sites_per_accession.ci.1 |
| C3.3 | n | 48,162 | analysis/results/g20_diagnostics.json : index_vs_sites_per_accession.human|gold_sites_per_accession.n_accessions |
| C3.6 | Estimate | 1.372 | analysis/results/b5_index_effect_sizes.json : means_by_top_tier.SEI_V50.GOLD |
| C3.6 | Estimate | 1.298 | analysis/results/b5_index_effect_sizes.json : means_by_top_tier.SEI_V50.SILVER |
| C3.6 | Estimate | 1.435 | analysis/results/b5_index_effect_sizes.json : means_by_top_tier.SEI_V50.BRONZE |
| C3.6 | Estimate | 0.0047 | analysis/results/b5_index_effect_sizes.json : effects.kruskal_eta2_SEI_V50_by_top_tier.estimate |
| C3.6 | 95% CI | 0.0038 | analysis/results/b5_index_effect_sizes.json : effects.kruskal_eta2_SEI_V50_by_top_tier.ci95.0 |
| C3.6 | 95% CI | 0.0057 | analysis/results/b5_index_effect_sizes.json : effects.kruskal_eta2_SEI_V50_by_top_tier.ci95.1 |
| C3.6 | n | 87,338 | analysis/results/b5_index_effect_sizes.json : n_proteins |
| C4.1 | Estimate | 0.432 | derived: ln spread of adjusted estimates across definitions, P-all S2 |
| C4.1 | Estimate | 2.840 | derived: ln spread of adjusted estimates across definitions, P-all S3 |
| C4.1 | Estimate | 1.54 | derived: plain ratio exp(ln spread across definitions), P-all S2 |
| C4.1 | Estimate | 17.11 | derived: disorder P-all spread, largest / smallest adjusted estimate (fold) |
| C4.1 | Estimate | 0.015 | derived: largest |ln(adjusted/crude)| within a definition, P-human S2 |
| C4.1 | Estimate | 0.399 | derived: largest |ln(adjusted/crude)| within a definition, P-human S4 |
| C4.1 | Estimate | 1.015 | derived: plain ratio exp(largest within-definition change), P-human S2 |
| C4.1 | Estimate | 1.490 | derived: plain ratio exp(largest within-definition change), P-human S4 |
| C4.1 | Estimate | 7 | derived: signal-by-population combinations where both can be computed |
| C4.1 | Estimate | 7 | derived: signal-by-population combinations where both can be computed |
| C4.1 | Estimate | 17.11 | derived: disorder P-all spread, largest / smallest adjusted estimate (fold) |
| C4.1 | Estimate | 3 | derived: ln spread of adjusted estimates across definitions, P-all S3 |
| C4.1 | Estimate | 2.99 | derived: tier spread, eligible definitions, P-all S3 |
| C4.2 | Claim | 0.14 | analysis/results/amend_results.json : G19-0458.estimate |
| C4.2 | Claim | 0.002 | analysis/results/amend_results.json : G19-0460.estimate |
| C4.2 | Estimate | 0.1393 | analysis/results/amend_results.json : G19-0458.estimate |
| C4.2 | Estimate | 0.0018 | analysis/results/amend_results.json : G19-0460.estimate |
| C4.2 | 95% CI | 0.1351 | analysis/results/amend_results.json : G19-0458.ci.0 |
| C4.2 | 95% CI | 0.1438 | analysis/results/amend_results.json : G19-0458.ci.1 |
| C4.2 | 95% CI | -0.0039 | analysis/results/amend_results.json : G19-0460.ci.0 |
| C4.2 | 95% CI | 0.0070 | analysis/results/amend_results.json : G19-0460.ci.1 |
| C4.2 | n | 229,448 | analysis/results/amend_results.json : G19-0458.n_sites_gold |
| C4.2 | n | 506,153 | analysis/results/amend_results.json : G19-0458.n_sites_bronze |
| C4.2 | n | 270,355 | analysis/results/amend_results.json : G19-0460.n_sites_gold |
| C4.2 | n | 270,815 | analysis/results/amend_results.json : G19-0460.n_sites_bronze |
| C4.3 | Estimate | 1.198 | analysis/results/g14_cem.json : G14-0345.ci_width_ratio_protein_cluster_over_rev5 |
| C4.3 | Estimate | 5.600 | analysis/results/g14_cem.json : G14-0335.ci_width_ratio_protein_cluster_over_rev5 |
| C4.3 | Estimate | 10 | derived: CEM rows with a protein-cluster to Rev 5 interval width ratio (G14) |
| C4.4 | Estimate | 0.993 | derived: width ratio G1 S4 D1 |
| C4.4 | Estimate | 1.161 | derived: width ratio G1 S3 D1 |
| C4.5 | Estimate | 2 | derived: decision-bearing G1 contrasts whose two-way interval excludes the null |
| C4.5 | Estimate | 4 | derived: decision-bearing G1 contrasts with a two-way interval |
| C4.6 | Estimate | 9 | derived: G23 D1 contrasts |
| C4.6 | Estimate | 9 | derived: G23 D1 contrasts |
| C4.6 | Estimate | 9 | derived: G23 D1 contrasts |
| C4.6 | Estimate | 9 | derived: G23 D1 contrasts |
| C4.6 | Estimate | 0 | derived: G23 D1 contrasts largely explained |
| C5.1 | Estimate | 9.144 | analysis/results/g21_results.json : G21-0491.estimate |
| C5.1 | Estimate | 8.596 | analysis/results/g21_results.json : G21-0491.ci.0 |
| C5.1 | Estimate | 9.698 | analysis/results/g21_results.json : G21-0491.ci.1 |
| C5.1 | Estimate | 0.3763 | analysis/results/g21_results.json : G21-0491.companion_linear.estimate |
| C5.1 | Estimate | 0.3675 | analysis/results/g21_results.json : G21-0491.companion_linear.ci.0 |
| C5.1 | Estimate | 0.3846 | analysis/results/g21_results.json : G21-0491.companion_linear.ci.1 |
| C5.1 | n | 130,884 | analysis/results/g21_results.json : G21-0491.n_records |
| C5.1 | n | 8,390 | analysis/results/g21_results.json : G21-0491.n_gene_clusters_with_both_tiers |
| C5.1 | Source | 9.512 | analysis/results/g21_results.json : G21-0493.estimate |
| C5.1 | Source | 8.959 | analysis/results/g21_results.json : G21-0493.ci.0 |
| C5.1 | Source | 10.041 | analysis/results/g21_results.json : G21-0493.ci.1 |
| C5.1 | Source | 0.3833 | analysis/results/g21_results.json : G21-0493.companion_linear.estimate |
| C5.1 | Source | 0.3747 | analysis/results/g21_results.json : G21-0493.companion_linear.ci.0 |
| C5.1 | Source | 0.3908 | analysis/results/g21_results.json : G21-0493.companion_linear.ci.1 |
| C5.1 | Basis | 0.3993 | analysis/results/amend_results.json : G17-0382.estimate |
| C5.1 | Basis | 0.3892 | analysis/results/amend_results.json : G17-0382.ci.0 |
| C5.1 | Basis | 0.4099 | analysis/results/amend_results.json : G17-0382.ci.1 |
| C5.1 | Basis | 1000 | analysis/results/g21_results.json : G21-0491.n_boot |
| C5.1 | Basis | 966 | analysis/results/g21_results.json : G21-0491.n_boot_valid |
| C5.1 | Basis | 34 | analysis/results/g21_results.json : G21-0491.n_boot_unfit |
| C5.1 | Basis | 959 | analysis/results/g21_results.json : G21-0492.n_boot_valid |
| C5.1 | Basis | 41 | analysis/results/g21_results.json : G21-0492.n_boot_unfit |
| C5.1 | Basis | 939 | analysis/results/g21_results.json : G21-0493.n_boot_valid |
| C5.1 | Basis | 61 | analysis/results/g21_results.json : G21-0493.n_boot_unfit |
| C5.1 | Basis | 5 | analysis/ANALYSIS_PLAN.md : more than 5% of the fits (a definition or threshold recorded in the plan, not a result) |
| C5.1 | Basis | 6.1 | derived: G21 failed-fit share, G21-0493 (percent) |
| C5.2 | Estimate | 10.127 | analysis/results/g21_results.json : G21-0503.estimate |
| C5.2 | Estimate | 9.425 | analysis/results/g21_results.json : G21-0503.ci.0 |
| C5.2 | Estimate | 10.904 | analysis/results/g21_results.json : G21-0503.ci.1 |
| C5.2 | Estimate | 0.4070 | analysis/results/g21_results.json : G21-0503.companion_linear.estimate |
| C5.2 | Estimate | 0.3955 | analysis/results/g21_results.json : G21-0503.companion_linear.ci.0 |
| C5.2 | Estimate | 0.4181 | analysis/results/g21_results.json : G21-0503.companion_linear.ci.1 |
| C5.2 | n | 67,404 | analysis/results/g21_results.json : G21-0503.n_records |
| C5.2 | n | 4,331 | analysis/results/g21_results.json : G21-0503.n_gene_clusters_with_both_tiers |
| C5.2 | Source | 10.288 | analysis/results/g21_results.json : G21-0505.estimate |
| C5.2 | Source | 9.520 | analysis/results/g21_results.json : G21-0505.ci.0 |
| C5.2 | Source | 11.132 | analysis/results/g21_results.json : G21-0505.ci.1 |
| C5.2 | Source | 0.4059 | analysis/results/g21_results.json : G21-0505.companion_linear.estimate |
| C5.2 | Source | 0.3941 | analysis/results/g21_results.json : G21-0505.companion_linear.ci.0 |
| C5.2 | Source | 0.4174 | analysis/results/g21_results.json : G21-0505.companion_linear.ci.1 |
| C5.2 | Basis | 0.4501 | analysis/results/amend_results.json : G17-0390.estimate |
| C5.2 | Basis | 0.4349 | analysis/results/amend_results.json : G17-0390.ci.0 |
| C5.2 | Basis | 0.4642 | analysis/results/amend_results.json : G17-0390.ci.1 |
| C5.2 | Basis | 1000 | analysis/results/g21_results.json : G21-0503.n_boot |
| C5.2 | Basis | 950 | analysis/results/g21_results.json : G21-0503.n_boot_valid |
| C5.2 | Basis | 50 | analysis/results/g21_results.json : G21-0503.n_boot_unfit |
| C5.2 | Basis | 946 | analysis/results/g21_results.json : G21-0504.n_boot_valid |
| C5.2 | Basis | 54 | analysis/results/g21_results.json : G21-0504.n_boot_unfit |
| C5.2 | Basis | 952 | analysis/results/g21_results.json : G21-0505.n_boot_valid |
| C5.2 | Basis | 48 | analysis/results/g21_results.json : G21-0505.n_boot_unfit |
| C5.2 | Basis | 5 | analysis/ANALYSIS_PLAN.md : more than 5% of the fits (a definition or threshold recorded in the plan, not a result) |
| C5.2 | Basis | 5.4 | derived: G21 failed-fit share, G21-0504 (percent) |
| C5.3 | Estimate | 0.0453 | analysis/results/g21_results.json : G21-0494.estimate |
| C5.3 | Estimate | 0.0431 | analysis/results/g21_results.json : G21-0494.ci.0 |
| C5.3 | Estimate | 0.0475 | analysis/results/g21_results.json : G21-0494.ci.1 |
| C5.3 | n | 536,483 | analysis/results/g21_results.json : G21-0494.n_records |
| C5.3 | n | 21,800 | analysis/results/g21_results.json : G21-0494.n_gene_clusters_with_both_tiers |
| C5.3 | Source | 0.0713 | analysis/results/g21_results.json : G21-0496.estimate |
| C5.3 | Source | 0.0686 | analysis/results/g21_results.json : G21-0496.ci.0 |
| C5.3 | Source | 0.0738 | analysis/results/g21_results.json : G21-0496.ci.1 |
| C5.3 | Basis | 0.0457 | analysis/results/amend_results.json : G17-0384.estimate |
| C5.3 | Basis | 0.0411 | analysis/results/amend_results.json : G17-0384.ci.0 |
| C5.3 | Basis | 0.0512 | analysis/results/amend_results.json : G17-0384.ci.1 |
| C5.3 | Basis | 1000 | analysis/results/g21_results.json : G21-0494.n_boot_valid |
| C5.3 | Basis | 1000 | analysis/results/g21_results.json : G21-0494.n_boot_valid |
| C5.3 | Basis | 0 | analysis/results/g21_results.json : G21-0494.n_boot_unfit |
| C5.3 | Basis | 1000 | analysis/results/g21_results.json : G21-0494.n_boot_valid |
| C5.3 | Basis | 0 | analysis/results/g21_results.json : G21-0494.n_boot_unfit |
| C5.3 | Basis | 1000 | analysis/results/g21_results.json : G21-0494.n_boot_valid |
| C5.3 | Basis | 0 | analysis/results/g21_results.json : G21-0494.n_boot_unfit |
| C5.4 | Estimate | 0.0473 | analysis/results/g21_results.json : G21-0506.estimate |
| C5.4 | Estimate | 0.0443 | analysis/results/g21_results.json : G21-0506.ci.0 |
| C5.4 | Estimate | 0.0502 | analysis/results/g21_results.json : G21-0506.ci.1 |
| C5.4 | n | 337,724 | analysis/results/g21_results.json : G21-0506.n_records |
| C5.4 | n | 12,784 | analysis/results/g21_results.json : G21-0506.n_gene_clusters_with_both_tiers |
| C5.4 | Source | 0.0791 | analysis/results/g21_results.json : G21-0508.estimate |
| C5.4 | Source | 0.0758 | analysis/results/g21_results.json : G21-0508.ci.0 |
| C5.4 | Source | 0.0824 | analysis/results/g21_results.json : G21-0508.ci.1 |
| C5.4 | Basis | 0.0534 | analysis/results/amend_results.json : G17-0392.estimate |
| C5.4 | Basis | 0.0474 | analysis/results/amend_results.json : G17-0392.ci.0 |
| C5.4 | Basis | 0.0611 | analysis/results/amend_results.json : G17-0392.ci.1 |
| C5.4 | Basis | 1000 | analysis/results/g21_results.json : G21-0506.n_boot_valid |
| C5.4 | Basis | 1000 | analysis/results/g21_results.json : G21-0506.n_boot_valid |
| C5.4 | Basis | 0 | analysis/results/g21_results.json : G21-0506.n_boot_unfit |
| C5.4 | Basis | 1000 | analysis/results/g21_results.json : G21-0506.n_boot_valid |
| C5.4 | Basis | 0 | analysis/results/g21_results.json : G21-0506.n_boot_unfit |
| C5.4 | Basis | 1000 | analysis/results/g21_results.json : G21-0506.n_boot_valid |
| C5.4 | Basis | 0 | analysis/results/g21_results.json : G21-0506.n_boot_unfit |
| C5.5 | Estimate | 1.343 | analysis/results/g21_results.json : G21-0497.estimate |
| C5.5 | Estimate | 1.326 | analysis/results/g21_results.json : G21-0497.ci.0 |
| C5.5 | Estimate | 1.359 | analysis/results/g21_results.json : G21-0497.ci.1 |
| C5.5 | Estimate | 0.0692 | analysis/results/g21_results.json : G21-0497.companion_linear.estimate |
| C5.5 | Estimate | 0.0665 | analysis/results/g21_results.json : G21-0497.companion_linear.ci.0 |
| C5.5 | Estimate | 0.0720 | analysis/results/g21_results.json : G21-0497.companion_linear.ci.1 |
| C5.5 | n | 803,744 | analysis/results/g21_results.json : G21-0497.n_records |
| C5.5 | n | 23,705 | analysis/results/g21_results.json : G21-0497.n_gene_clusters_with_both_tiers |
| C5.5 | Source | 1.992 | analysis/results/g21_results.json : G21-0499.estimate |
| C5.5 | Source | 1.953 | analysis/results/g21_results.json : G21-0499.ci.0 |
| C5.5 | Source | 2.033 | analysis/results/g21_results.json : G21-0499.ci.1 |
| C5.5 | Source | 0.1591 | analysis/results/g21_results.json : G21-0499.companion_linear.estimate |
| C5.5 | Source | 0.1547 | analysis/results/g21_results.json : G21-0499.companion_linear.ci.0 |
| C5.5 | Source | 0.1636 | analysis/results/g21_results.json : G21-0499.companion_linear.ci.1 |
| C5.5 | Basis | 0.0618 | analysis/results/amend_results.json : G17-0386.estimate |
| C5.5 | Basis | 0.0587 | analysis/results/amend_results.json : G17-0386.ci.0 |
| C5.5 | Basis | 0.0652 | analysis/results/amend_results.json : G17-0386.ci.1 |
| C5.5 | Basis | 1000 | analysis/results/g21_results.json : G21-0497.n_boot |
| C5.5 | Basis | 998 | analysis/results/g21_results.json : G21-0497.n_boot_valid |
| C5.5 | Basis | 2 | analysis/results/g21_results.json : G21-0497.n_boot_unfit |
| C5.5 | Basis | 997 | analysis/results/g21_results.json : G21-0498.n_boot_valid |
| C5.5 | Basis | 3 | analysis/results/g21_results.json : G21-0498.n_boot_unfit |
| C5.5 | Basis | 998 | analysis/results/g21_results.json : G21-0497.n_boot_valid |
| C5.5 | Basis | 2 | analysis/results/g21_results.json : G21-0497.n_boot_unfit |
| C5.6 | Estimate | 1.226 | analysis/results/g21_results.json : G21-0509.estimate |
| C5.6 | Estimate | 1.209 | analysis/results/g21_results.json : G21-0509.ci.0 |
| C5.6 | Estimate | 1.243 | analysis/results/g21_results.json : G21-0509.ci.1 |
| C5.6 | Estimate | 0.0479 | analysis/results/g21_results.json : G21-0509.companion_linear.estimate |
| C5.6 | Estimate | 0.0447 | analysis/results/g21_results.json : G21-0509.companion_linear.ci.0 |
| C5.6 | Estimate | 0.0511 | analysis/results/g21_results.json : G21-0509.companion_linear.ci.1 |
| C5.6 | n | 535,848 | analysis/results/g21_results.json : G21-0509.n_records |
| C5.6 | n | 13,801 | analysis/results/g21_results.json : G21-0509.n_gene_clusters_with_both_tiers |
| C5.6 | Source | 1.765 | analysis/results/g21_results.json : G21-0511.estimate |
| C5.6 | Source | 1.726 | analysis/results/g21_results.json : G21-0511.ci.0 |
| C5.6 | Source | 1.805 | analysis/results/g21_results.json : G21-0511.ci.1 |
| C5.6 | Source | 0.1311 | analysis/results/g21_results.json : G21-0511.companion_linear.estimate |
| C5.6 | Source | 0.1262 | analysis/results/g21_results.json : G21-0511.companion_linear.ci.0 |
| C5.6 | Source | 0.1360 | analysis/results/g21_results.json : G21-0511.companion_linear.ci.1 |
| C5.6 | Basis | 0.0368 | analysis/results/amend_results.json : G17-0394.estimate |
| C5.6 | Basis | 0.0331 | analysis/results/amend_results.json : G17-0394.ci.0 |
| C5.6 | Basis | 0.0404 | analysis/results/amend_results.json : G17-0394.ci.1 |
| C5.6 | Basis | 1000 | analysis/results/g21_results.json : G21-0510.n_boot_valid |
| C5.6 | Basis | 999 | analysis/results/g21_results.json : G21-0509.n_boot_valid |
| C5.6 | Basis | 1 | analysis/results/g21_results.json : G21-0509.n_boot_unfit |
| C5.6 | Basis | 1000 | analysis/results/g21_results.json : G21-0510.n_boot_valid |
| C5.6 | Basis | 0 | analysis/results/g21_results.json : G21-0510.n_boot_unfit |
| C5.6 | Basis | 999 | analysis/results/g21_results.json : G21-0509.n_boot_valid |
| C5.6 | Basis | 1 | analysis/results/g21_results.json : G21-0509.n_boot_unfit |
| C5.7 | Estimate | 160.195 | analysis/results/g21_results.json : G21-0500.estimate |
| C5.7 | Estimate | 152.282 | analysis/results/g21_results.json : G21-0500.ci.0 |
| C5.7 | Estimate | 167.980 | analysis/results/g21_results.json : G21-0500.ci.1 |
| C5.7 | Estimate | 0.7872 | analysis/results/g21_results.json : G21-0500.companion_linear.estimate |
| C5.7 | Estimate | 0.7831 | analysis/results/g21_results.json : G21-0500.companion_linear.ci.0 |
| C5.7 | Estimate | 0.7910 | analysis/results/g21_results.json : G21-0500.companion_linear.ci.1 |
| C5.7 | n | 803,744 | analysis/results/g21_results.json : G21-0500.n_records |
| C5.7 | n | 23,705 | analysis/results/g21_results.json : G21-0500.n_gene_clusters_with_both_tiers |
| C5.7 | Source | 509.494 | analysis/results/g21_results.json : G21-0502.estimate |
| C5.7 | Source | 457.686 | analysis/results/g21_results.json : G21-0502.ci.0 |
| C5.7 | Source | 563.220 | analysis/results/g21_results.json : G21-0502.ci.1 |
| C5.7 | Source | 0.8180 | analysis/results/g21_results.json : G21-0502.companion_linear.estimate |
| C5.7 | Source | 0.8132 | analysis/results/g21_results.json : G21-0502.companion_linear.ci.0 |
| C5.7 | Source | 0.8227 | analysis/results/g21_results.json : G21-0502.companion_linear.ci.1 |
| C5.7 | Basis | 0.7919 | analysis/results/amend_results.json : G17-0388.estimate |
| C5.7 | Basis | 0.7887 | analysis/results/amend_results.json : G17-0388.ci.0 |
| C5.7 | Basis | 0.7951 | analysis/results/amend_results.json : G17-0388.ci.1 |
| C5.7 | Basis | 1000 | analysis/results/g21_results.json : G21-0500.n_boot |
| C5.7 | Basis | 927 | analysis/results/g21_results.json : G21-0500.n_boot_valid |
| C5.7 | Basis | 73 | analysis/results/g21_results.json : G21-0500.n_boot_unfit |
| C5.7 | Basis | 889 | analysis/results/g21_results.json : G21-0501.n_boot_valid |
| C5.7 | Basis | 111 | analysis/results/g21_results.json : G21-0501.n_boot_unfit |
| C5.7 | Basis | 892 | analysis/results/g21_results.json : G21-0502.n_boot_valid |
| C5.7 | Basis | 108 | analysis/results/g21_results.json : G21-0502.n_boot_unfit |
| C5.8 | Estimate | 192.318 | analysis/results/g21_results.json : G21-0512.estimate |
| C5.8 | Estimate | 180.580 | analysis/results/g21_results.json : G21-0512.ci.0 |
| C5.8 | Estimate | 204.799 | analysis/results/g21_results.json : G21-0512.ci.1 |
| C5.8 | Estimate | 0.7903 | analysis/results/g21_results.json : G21-0512.companion_linear.estimate |
| C5.8 | Estimate | 0.7850 | analysis/results/g21_results.json : G21-0512.companion_linear.ci.0 |
| C5.8 | Estimate | 0.7955 | analysis/results/g21_results.json : G21-0512.companion_linear.ci.1 |
| C5.8 | n | 535,848 | analysis/results/g21_results.json : G21-0512.n_records |
| C5.8 | n | 13,801 | analysis/results/g21_results.json : G21-0512.n_gene_clusters_with_both_tiers |
| C5.8 | Source | 518.874 | analysis/results/g21_results.json : G21-0514.estimate |
| C5.8 | Source | 451.143 | analysis/results/g21_results.json : G21-0514.ci.0 |
| C5.8 | Source | 593.605 | analysis/results/g21_results.json : G21-0514.ci.1 |
| C5.8 | Source | 0.8055 | analysis/results/g21_results.json : G21-0514.companion_linear.estimate |
| C5.8 | Source | 0.7991 | analysis/results/g21_results.json : G21-0514.companion_linear.ci.0 |
| C5.8 | Source | 0.8117 | analysis/results/g21_results.json : G21-0514.companion_linear.ci.1 |
| C5.8 | Basis | 0.7964 | analysis/results/amend_results.json : G17-0396.estimate |
| C5.8 | Basis | 0.7918 | analysis/results/amend_results.json : G17-0396.ci.0 |
| C5.8 | Basis | 0.8007 | analysis/results/amend_results.json : G17-0396.ci.1 |
| C5.8 | Basis | 1000 | analysis/results/g21_results.json : G21-0512.n_boot |
| C5.8 | Basis | 908 | analysis/results/g21_results.json : G21-0512.n_boot_valid |
| C5.8 | Basis | 92 | analysis/results/g21_results.json : G21-0512.n_boot_unfit |
| C5.8 | Basis | 919 | analysis/results/g21_results.json : G21-0513.n_boot_valid |
| C5.8 | Basis | 81 | analysis/results/g21_results.json : G21-0513.n_boot_unfit |
| C5.8 | Basis | 906 | analysis/results/g21_results.json : G21-0514.n_boot_valid |
| C5.8 | Basis | 94 | analysis/results/g21_results.json : G21-0514.n_boot_unfit |
| C6.1 | Claim | 0.05 | analysis/ANALYSIS_PLAN.md : margin +/-0.05 (a definition or threshold recorded in the plan, not a result) |
| C6.1 | Estimate | 1 | derived: predictor tools meeting the pre-specified criterion (gap >= 0.05, interval excluding zero) |
| C6.1 | Estimate | 3 | derived: predictor tools audited |
| C6.2 | Estimate | 0.0611 | analysis/results/g11_predictors.json : MusiteDeep|planned: S/T sites only.unrestricted.gap |
| C6.2 | 95% CI | 0.0324 | analysis/results/g11_predictors.json : MusiteDeep|planned: S/T sites only.unrestricted.gene_ci95.0 |
| C6.2 | 95% CI | 0.0894 | analysis/results/g11_predictors.json : MusiteDeep|planned: S/T sites only.unrestricted.gene_ci95.1 |
| C6.2 | n | 8,641 | analysis/results/g11_predictors.json : MusiteDeep|planned: S/T sites only.unrestricted.n_sites |
| C6.2 | n | 4,485 | analysis/results/g11_predictors.json : MusiteDeep|planned: S/T sites only.unrestricted.n_positives |
| C6.2 | n | 5,209 | analysis/results/g11_predictors.json : MusiteDeep|planned: S/T sites only.unrestricted.n_gene_clusters |

## Headline results of RESULTS_REPORT.md: numbers equal to a number bound or derived in a claim row (same value, same decimals)

Bound through a claim row: 44. Not found in any claim row: 0.

none

| Item | Number | Bound through |
|---|---|---|
| 1 | 7 | derived: headline 1: decision-bearing D1 contrasts in G1 and G2 (and G3) |
| 1 | 7 | derived: headline 1: decision-bearing D1 contrasts in G1 and G2 (and G3) |
| 1 | 7 | derived: headline 1: decision-bearing D1 contrasts in G1 and G2 (and G3) |
| 1 | 7 | derived: headline 1: decision-bearing D1 contrasts in G1 and G2 (and G3) |
| 1 | 0 | derived: headline 1: of them largely explained |
| 1 | 1.4 | claim row C3.1: analysis/results/amend_results.json : G17-0392.attenuation (x 100) |
| 1 | 2 | derived: headline 1: of them excluding the null under the two-way interval |
| 1 | 4 | derived: headline 1: decision-bearing G1 contrasts with a two-way interval |
| 2 | 6.065 | claim row C2.6: analysis/results/amend_results.json : G17-0379.estimate |
| 2 | 0.741 | claim row C2.6: analysis/results/amend_results.json : G17-0399.estimate |
| 2 | 0.692 | claim row C2.6: analysis/results/amend_results.json : G17-0397.estimate |
| 2 | 0.746 | claim row C2.6: analysis/results/amend_results.json : G16-0361.estimate |
| 3 | 41 | derived: headline 3: of them with |SMD| of 0.1 or more before matching |
| 3 | 78 | derived: headline 3: decision-bearing contrasts with a before/after balance |
| 3 | 0.1 | derived: headline 3: largest within-tier correlation of the index with an outcome |
| 3 | 0.270 | derived: headline 3: largest |SMD| before matching |
| 3 | 0.160 | derived: headline 3: SMD before matching, G1-0005 |
| 3 | 0.187 | derived: headline 3: SMD before matching, G2-0017 |
| 3 | 0.017 | claim row C3.3: analysis/results/g20_diagnostics.json : index_vs_sites_per_accession.human|gold_sites_per_accession.ci.1 |
| 3 | -0.0881 | derived: headline 3: smallest within-tier correlation of the index with an outcome |
| 3 | 0.0671 | derived: headline 3: largest within-tier correlation of the index with an outcome |
| 4 | 7 | derived: headline 4: decision-bearing D2b contrasts |
| 4 | 7 | derived: headline 4: decision-bearing D2b contrasts |
| 4 | 2.9 | derived: headline 4: largest attenuation (percent) |
| 4 | 0 | derived: headline 4: D2b point estimates outside their own interval |
| 5 | 7 | derived: signal-by-population combinations where both can be computed |
| 5 | 7 | derived: signal-by-population combinations where both can be computed |
| 5 | 0.432 | claim row C4.1: derived: ln spread of adjusted estimates across definitions, P-all S2 |
| 5 | 2.840 | claim row C4.1: derived: ln spread of adjusted estimates across definitions, P-all S3 |
| 5 | 1.54 | claim row C4.1: derived: plain ratio exp(ln spread across definitions), P-all S2 |
| 5 | 17.11 | claim row C4.1: derived: disorder P-all spread, largest / smallest adjusted estimate (fold) |
| 5 | 0.015 | claim row C4.1: derived: largest |ln(adjusted/crude)| within a definition, P-human S2 |
| 5 | 0.399 | claim row C4.1: derived: largest |ln(adjusted/crude)| within a definition, P-human S4 |
| 5 | 1.015 | claim row C4.1: derived: plain ratio exp(largest within-definition change), P-human S2 |
| 5 | 1.490 | claim row C4.1: derived: plain ratio exp(largest within-definition change), P-human S4 |
| 5 | 17 | derived: disorder P-all spread, largest / smallest adjusted estimate (fold) |
| 5 | 3 | derived: ln spread of adjusted estimates across definitions, P-all S3 |
| 5 | 2.99 | claim row C4.1: derived: tier spread, eligible definitions, P-all S3 |
| 6 | 5 | derived: headline 6: D2 contrasts that could not be fitted |
| 6 | 11 | derived: headline 6: planned D2 contrasts |
| 6 | 7 | derived: headline 6: evaluated predictor settings |
| 6 | 2 | derived: headline 6: settings with a gap equivalent to zero |
| 6 | 5 | derived: headline 6: D2 contrasts that could not be fitted |
| 6 | 0 | derived: headline 6: settings showing a gap |
