import numpy as np
import pandas as pd
import pytest

from src.analysis import (
    calculate_library_sizes,
    calculate_sample_correlations,
    normalize_log_cpm,
    run_pca,
)


def make_counts():
    return pd.DataFrame(
        {
            "control_1": [100, 10, 50, 5],
            "control_2": [110, 12, 55, 6],
            "treated_1": [20, 90, 50, 10],
            "treated_2": [25, 100, 52, 11],
        },
        index=["gene_1", "gene_2", "gene_3", "gene_4"],
    )


def make_metadata():
    return pd.DataFrame(
        {
            "condition": [
                "control",
                "control",
                "treated",
                "treated",
            ]
        },
        index=[
            "control_1",
            "control_2",
            "treated_1",
            "treated_2",
        ],
    )


def test_library_sizes_are_correct():
    result = calculate_library_sizes(make_counts())

    assert result["total_counts"].tolist() == [
        165,
        183,
        170,
        188,
    ]


def test_log_cpm_is_finite():
    normalized = normalize_log_cpm(make_counts())

    assert normalized.shape == make_counts().shape
    assert np.isfinite(normalized.to_numpy()).all()


def test_pca_returns_all_samples():
    result = run_pca(make_counts(), make_metadata())

    assert len(result.coordinates) == 4
    assert {"sample_id", "PC1", "PC2", "condition"}.issubset(
        result.coordinates.columns
    )
    assert result.variable_gene_count >= 2
    assert result.explained_variance[0] > 0


def test_sample_correlation_diagonal_is_one():
    correlations = calculate_sample_correlations(
        make_counts()
    )

    assert np.allclose(
        np.diag(correlations),
        np.ones(4),
    )


def test_pca_rejects_mismatched_samples():
    metadata = make_metadata().drop(index="treated_2")

    with pytest.raises(ValueError):
        run_pca(make_counts(), metadata)