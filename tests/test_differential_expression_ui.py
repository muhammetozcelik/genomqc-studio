import pandas as pd
import pytest

from src.differential_expression_ui import (
    _apply_significance_thresholds,
)


def make_results() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_id": [
                "gene_up",
                "gene_down",
                "gene_high_padj",
                "gene_small_change",
                "gene_missing_padj",
            ],
            "padj": [
                0.01,
                0.02,
                0.20,
                0.01,
                float("nan"),
            ],
            "log2FoldChange": [
                2.0,
                -1.5,
                3.0,
                0.5,
                -2.0,
            ],
        }
    )


def test_significance_thresholds_classify_genes():
    results = make_results()

    classified = _apply_significance_thresholds(
        results=results,
        adjusted_pvalue_threshold=0.05,
        log2_fold_change_threshold=1.0,
    )

    assert classified["significant"].tolist() == [
        True,
        True,
        False,
        False,
        False,
    ]

    assert classified["direction"].tolist() == [
        "Upregulated",
        "Downregulated",
        "Not significant",
        "Not significant",
        "Not significant",
    ]


def test_lower_fold_change_threshold_changes_classification():
    results = make_results()

    classified = _apply_significance_thresholds(
        results=results,
        adjusted_pvalue_threshold=0.05,
        log2_fold_change_threshold=0.5,
    )

    assert classified.loc[
        classified["gene_id"] == "gene_small_change",
        "direction",
    ].iloc[0] == "Upregulated"


def test_original_results_are_not_modified():
    results = make_results()

    _apply_significance_thresholds(
        results=results,
        adjusted_pvalue_threshold=0.05,
        log2_fold_change_threshold=1.0,
    )

    assert "significant" not in results.columns
    assert "direction" not in results.columns


@pytest.mark.parametrize(
    "threshold",
    [-0.01, 1.01],
)
def test_invalid_adjusted_pvalue_threshold_is_rejected(
    threshold: float,
):
    with pytest.raises(
        ValueError,
        match="Adjusted p-value threshold",
    ):
        _apply_significance_thresholds(
            results=make_results(),
            adjusted_pvalue_threshold=threshold,
            log2_fold_change_threshold=1.0,
        )


def test_negative_fold_change_threshold_is_rejected():
    with pytest.raises(
        ValueError,
        match="Log2 fold-change threshold",
    ):
        _apply_significance_thresholds(
            results=make_results(),
            adjusted_pvalue_threshold=0.05,
            log2_fold_change_threshold=-1.0,
        )