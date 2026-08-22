from pathlib import Path
import shutil

import pytest

from src.data_loader import (
    load_count_matrix,
    load_metadata,
)
from src.differential_expression import run_deseq2


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEMO_COUNTS = PROJECT_ROOT / "demo_data" / "counts.csv"
DEMO_METADATA = PROJECT_ROOT / "demo_data" / "metadata.csv"


@pytest.mark.integration
@pytest.mark.skipif(
    shutil.which("Rscript") is None,
    reason="Rscript is required for the DESeq2 integration test.",
)
def test_demo_dataset_runs_through_deseq2():
    counts = load_count_matrix(DEMO_COUNTS)
    metadata = load_metadata(DEMO_METADATA)

    result = run_deseq2(
        counts=counts,
        metadata=metadata,
        condition_column="condition",
        reference_condition="control",
        comparison_condition="treated",
        timeout_seconds=180,
    )

    expected_result_columns = {
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
        "significant",
        "direction",
    }

    assert not result.results.empty
    assert expected_result_columns.issubset(
        result.results.columns
    )
    assert set(result.results["gene_id"]).issubset(
        set(counts.index)
    )

    expected_count_columns = {
        "gene_id",
        *counts.columns,
    }

    assert (
        set(result.normalized_counts.columns)
        == expected_count_columns
    )
    assert (
        set(result.vst_counts.columns)
        == expected_count_columns
    )

    adjusted_pvalues = result.results["padj"].dropna()
    assert adjusted_pvalues.between(0, 1).all()

    assert result.reference_condition == "control"
    assert result.comparison_condition == "treated"
    assert result.summary["input_genes"] == str(
        counts.shape[0]
    )
    assert result.summary["samples"] == str(
        counts.shape[1]
    )