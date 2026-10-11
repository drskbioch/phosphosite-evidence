"""Smoke and unit tests on tiny in-memory fixtures (no external data)."""
import gzip

import numpy as np
import pandas as pd
import pytest

from phosnet import (
    bootstrap_auroc_ci,
    build_sei_table,
    compute_sei_primary,
    compute_tool_metrics,
    matched_comparison,
    sei_deciles,
    standardized_mean_difference,
)


@pytest.fixture
def sei_inputs(tmp_path):
    """Three proteins in one species, with a gene2pubmed table."""
    idmap = tmp_path / "idmapping_9606.dat.gz"
    with gzip.open(idmap, "wt") as fh:
        fh.write("P00001\tGeneID\t100\n")
        fh.write("P00002\tGeneID\t200\n")
        fh.write("P00003\tGeneID\t300\n")
    g2p = tmp_path / "gene2pubmed.gz"
    with gzip.open(g2p, "wt") as fh:
        fh.write("#tax_id\tGeneID\tPubMed_ID\n")
        for pmid in range(1, 10):  # gene 100: 9 PMIDs
            fh.write(f"9606\t100\t{pmid}\n")
        fh.write("9606\t200\t1\n")  # gene 200: 1 PMID
        fh.write("10090\t300\t5\n")  # other species, ignored
    proteins = pd.DataFrame({
        "uniprot_id": ["P00001", "P00002", "P00003", "P99999", "P00001-2"],
        "taxonomy_id": [9606] * 5,
    })
    return proteins, g2p, {9606: idmap}


def test_compute_sei_primary_formula():
    assert compute_sei_primary(0) == 0.0
    assert compute_sei_primary(9) == pytest.approx(1.0)
    assert np.allclose(compute_sei_primary([0, 99]), [0.0, 2.0])


def test_build_sei_table_smoke(sei_inputs):
    proteins, g2p, idmaps = sei_inputs
    table = build_sei_table(proteins, g2p, idmaps, tax_names={9606: "human"})
    by_id = table.set_index("uniprot_id")

    assert by_id.loc["P00001", "n_pubmed"] == 9
    assert by_id.loc["P00001", "SEI_primary"] == pytest.approx(1.0)
    assert by_id.loc["P00002", "n_pubmed"] == 1
    assert by_id.loc["P00002", "SEI_primary"] == pytest.approx(np.log10(2))
    # mapped to a gene, but no PubMed records for that species
    assert by_id.loc["P00003", "mapped_to_gene"]
    assert by_id.loc["P00003", "n_pubmed"] == 0
    assert by_id.loc["P00003", "zero_pubmed_flag"]
    # unknown accession
    assert not by_id.loc["P99999", "mapped_to_gene"]
    assert np.isnan(by_id.loc["P99999", "SEI_primary"])
    # isoform falls back to base accession
    assert by_id.loc["P00001-2", "matched_via_base_accession"]
    assert by_id.loc["P00001-2", "n_pubmed"] == 9
    assert (table["species"] == "human").all()


def test_sei_deciles_are_balanced():
    values = pd.Series(np.arange(1000, dtype=float))
    deciles, n_bins = sei_deciles(values, n_deciles=10)
    assert n_bins == 10
    assert deciles.min() == 0 and deciles.max() == 9
    assert (deciles.value_counts() == 100).all()


def test_sei_deciles_constant_input_forms_no_bins():
    deciles, n_bins = sei_deciles(pd.Series([1.0] * 50), n_deciles=10)
    assert n_bins == 0
    assert deciles.isna().all()


def test_standardized_mean_difference():
    a = np.array([1.0, 2.0, 3.0])
    assert standardized_mean_difference(a, a) == 0.0
    assert standardized_mean_difference(a + 1.0, a) > 0
    assert standardized_mean_difference(np.ones(4), np.ones(4)) == 0.0


def test_matched_comparison_recovers_group_effect():
    rng = np.random.default_rng(0)
    n = 4000
    sei = rng.uniform(0, 3, n)
    treated = rng.random(n) < 0.5
    indicator = rng.normal(0, 1, n) + 0.5 * treated  # true difference 0.5
    df = pd.DataFrame({
        "tier": np.where(treated, "GOLD", "BRONZE"),
        "sei": sei,
        "x": indicator,
        "species": "human",
    })
    result = matched_comparison(
        df, group_col="tier", treated_value="GOLD", indicator_col="x",
        statistic="proportions_ztest", sei_col="sei",
        extra_strata_cols=("species",), seed=1, n_boot=200,
    )
    assert result["balance_pass"]
    assert result["ci_excludes_null"]
    assert result["ci_low"] < 0.5 < result["ci_high"]


def test_bootstrap_auroc_ci_is_seed_deterministic():
    rng = np.random.default_rng(3)
    y = rng.integers(0, 2, 300)
    s = y * 0.8 + rng.normal(size=300)
    first = bootstrap_auroc_ci(y, s, seed=11, n_boot=100)
    second = bootstrap_auroc_ci(y, s, seed=11, n_boot=100)
    assert first == second
    assert first[1] < first[0] < first[2]


def test_compute_tool_metrics_smoke():
    rng = np.random.default_rng(5)
    n = 3000
    label = rng.integers(0, 2, n)
    df = pd.DataFrame({
        "label_type": np.where(label == 1, "positive", "negative"),
        "score": label * 1.0 + rng.normal(size=n),
        "SEI_decile": rng.integers(1, 11, n),
        "first_reported_year": rng.integers(1995, 2025, n),
    })
    result = compute_tool_metrics(df, seed=7, model_pub_year=2010, n_boot=100)

    assert result["n_scored"] == n
    assert set(result["per_decile"]) == set(range(1, 11))
    assert 0.6 < result["pooled_auroc"]["point"] < 0.9
    gap = result["gap_top_minus_bottom"]
    assert gap["bottom_deciles"] == [1, 2, 3] and gap["top_deciles"] == [8, 9, 10]
    assert result["leakage_free"]["model_pub_year"] == 2010
    assert result["leakage_free"]["status"] == "evaluated"
