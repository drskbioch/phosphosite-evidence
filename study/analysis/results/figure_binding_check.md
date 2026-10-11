# Binding of the numbers in the Rev 6 figures and tables

21 figures and tables. **Numbers formatted from tracked result values and recomputed here from the named file and key: 793 verified, 0 not reproduced.** Derived values recomputed (`derive_figure_values.py`): 81 values, 0 differ. Counts file: input hashes recomputed, 0 changed; cross-check against dataset_counts.json: 26 of 26 equal, differences: none. Text cells re-verified against their source text: 196 verified, 0 not. Composite SVGs scanned (they contain every panel). SVG text and table scan (Table S1 has the 534 contrasts of manifest_status.csv with equal identifiers and values at 4 significant digits): 800 numbers checked, **0 unsourced**.

## Not reproduced

none

## Derived values that differ

none

## Text cells not verified

none

## Unsourced numbers in SVG text

none

## Definitions, names and layout constants written in figure text (not results)

- `45,000,000`: identifier above this value is invalid: definition in scripts/earlier_figures/compute_results.py (`num > 45_000_000`)
- phrase `cluster thresholds 500 and 1,000 sites`: study-cluster thresholds of block G13 (plan; variation cluster_threshold=500 and 1000 in cem_results.json)
- phrase `EPSD 2\.0|PaxDb 6\.0`: proper names of the resources
- phrase `\d+:\d+`: matching ratio of the 1:1 caliper matcher (plan, design D2b)
- phrase `scores? (?:of )?2 or more|score 2 or more`: database agreement score of at least 2 (DAS >= 2 in the plan)
- phrase `\(?(?:agreement )?score 0 to [23]\)?|\(0 to [23]\)`: range of the agreement score under source set A (0 to 2) and B (0 to 3)
- phrase `5th to 95th percentile`: whisker definition of the box plots
- phrase `\b9[05]% (?:gene-cluster |)interval`: interval level: 95% (primary) and 90% (equivalence test), as in the plan
- phrase `\b95% gene-cluster|\b95% intervals|\b90% interval|\bthe 95%|the 90%`: interval level: 95% (primary) and 90% (equivalence test), as in the plan
- phrase `deciles 8 to 10 minus deciles 1 to 3`: definition of the predictor gap (top three minus bottom three deciles, plan)
- phrase `axis begins\s+in 1998`: axis range of the plot (layout); the text may wrap between the words
- phrase `^in 1998\.$`: axis range of the plot (layout): last line of the wrapped sentence 'the axis begins in 1998'
- phrase `\b95%\s*$`: interval level 95% (primary), at the end of a wrapped caption line whose next line starts with "gene-cluster intervals"
- phrase `per 100,000 sites`: unit of a rate
- phrase `(?:a score of |score )2`: database agreement score of at least 2
