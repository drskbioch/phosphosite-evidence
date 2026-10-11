# phosnet

Phosphosite databases and the predictors validated against them are
usually evaluated as though every protein had been studied equally.
`phosnet` provides the pieces needed to check that assumption:

- **Study-Effort Index (SEI)**: a per-protein covariate,
  `SEI = log10(1 + n_pubmed)`, where `n_pubmed` is the number of distinct
  PubMed IDs linked to the protein's NCBI Gene ID(s) in `gene2pubmed`. It
  comes from an external literature-to-gene mapping, so it is independent
  of any site-level evidence score.
- **Matched comparison**: coarsened exact matching (primary) and
  propensity-score matching (sensitivity) that compare two groups of sites,
  such as high- and low-confidence tiers, while holding SEI constant.
- **SEI-stratified evaluation**: per-SEI-decile AUROC for any predictor, the
  difference in AUROC between the best-studied and least-studied deciles,
  pooled AUROC, and an optional analysis restricted to positives first
  reported after a model's publication year.

phosnet 0.1.0 implements earlier estimators; the results of the manuscript were produced by the scripts in study/scripts/analysis/ (see REPRODUCE.md).

Bootstrap intervals take an explicit `seed`; every interval derives its own
random generator from it, so results do not depend on call order.

## Installation

```bash
pip install phosnet
```

Python 3.9 or later. Dependencies: numpy, pandas, scipy, scikit-learn,
statsmodels. For a source checkout: `pip install -e ".[test]"`, then
`pytest`.

## Example

The example below is self-contained and uses simulated data.

```python
import numpy as np
import pandas as pd
from phosnet import compute_sei_primary, sei_deciles, compute_tool_metrics, matched_comparison

rng = np.random.default_rng(0)

# 1. SEI from publication counts
n_pubmed = rng.negative_binomial(1, 0.02, size=5000)
sei = compute_sei_primary(n_pubmed)
deciles, n_bins = sei_deciles(pd.Series(sei), n_deciles=10)  # 0-based

# 2. SEI-stratified evaluation of a predictor
label = rng.random(5000) < 0.5
sites = pd.DataFrame({
    "label_type": np.where(label, "positive", "negative"),
    "score": label * 1.0 + rng.normal(size=5000),
    "SEI_decile": deciles.to_numpy() + 1,          # coded 1..10
    "first_reported_year": rng.integers(1995, 2025, 5000),
})
res = compute_tool_metrics(sites, seed=1, model_pub_year=2015, n_boot=200)
gap = res["gap_top_minus_bottom"]
print(f"pooled AUROC {res['pooled_auroc']['point']:.3f}; "
      f"top-3 minus bottom-3 decile AUROC {gap['point']:+.3f} "
      f"[{gap['ci_low']:+.3f}, {gap['ci_high']:+.3f}]")

# 3. Matched comparison of two groups at equal SEI
group = rng.random(5000) < 0.5
df = pd.DataFrame({
    "tier": np.where(group, "GOLD", "BRONZE"),
    "SEI": sei,
    "disorder": rng.normal(0.45 + 0.04 * group, 0.2),
    "residue_type": rng.choice(list("STY"), 5000),
})
out = matched_comparison(
    df, group_col="tier", treated_value="GOLD", indicator_col="disorder",
    statistic="mannwhitneyu", sei_col="SEI",
    extra_strata_cols=("residue_type",), seed=1, n_boot=200,
)
print(f"matched difference {out['point_estimate']:.3f} "
      f"[{out['ci_low']:.3f}, {out['ci_high']:.3f}]; "
      f"|SMD| < 0.1: {out['balance_pass']}")
```

To build SEI for real proteins, use `phosnet.build_sei_table(proteins,
gene2pubmed_path, idmap_paths)` with NCBI `gene2pubmed.gz` and UniProt
`idmapping.dat.gz` files.

## API summary

| Function | Purpose |
|---|---|
| `compute_sei_primary(n_pubmed)` | `log10(1 + n_pubmed)` |
| `build_sei_table(...)` | protein-level primary SEI from `gene2pubmed` and UniProt ID mappings |
| `build_sei_secondary_table(...)` | SEI based on distinct study clusters instead of PMIDs |
| `sei_deciles(values, n_deciles=10)` | quantile bins of an SEI vector |
| `matched_comparison(...)` | coarsened exact matching on SEI and optional strata |
| `standardized_mean_difference(a, b)` | balance diagnostic |
| `compute_tool_metrics(...)` | SEI-stratified evaluation of one predictor |
| `bootstrap_auroc_ci`, `bootstrap_gap_ci` | percentile bootstrap intervals |

## Citation

If you use phosnet, please cite the accompanying paper (reference to be
added on publication) and the software version:
`phosnet` 0.1.0.

## License

MIT. See `LICENSE`.
