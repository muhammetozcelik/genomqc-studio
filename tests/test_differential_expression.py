import pandas as pd
import pytest

from src.differential_expression import (
    DESeq2Error,
    run_deseq2,
    validate_comparison,
)


SAMPLES = [
    "control_1",
    "control_2",
    "control_3",
    "treated_1",
    "treated_2",
    "treated_3",
]


def make_metadata():
    return pd.DataFrame(
        {
            "condition": [
                "control",
                "control",
                "control",
                "treated",
                "treated",
                "treated",
            ]
        },
        index=SAMPLES,
    )


def make_counts():
    return pd.DataFrame(
        {
            sample: [10, 20, 30]
            for sample in SAMPLES
        },
        index=["gene_1", "gene_2", "gene_3"],
    )


def test_valid_comparison_passes():
    validate_comparison(
        metadata=make_metadata(),
        condition_column="condition",
        reference_condition="control",
        comparison_condition="treated",
    )


def test_identical_conditions_are_rejected():
    with pytest.raises(DESeq2Error):
        validate_comparison(
            metadata=make_metadata(),
            condition_column="condition",
            reference_condition="control",
            comparison_condition="control",
        )


def test_unknown_condition_column_is_rejected():
    with pytest.raises(DESeq2Error):
        validate_comparison(
            metadata=make_metadata(),
            condition_column="group",
            reference_condition="control",
            comparison_condition="treated",
        )


def test_insufficient_replication_is_rejected():
    metadata = make_metadata().drop(
        index=["treated_2", "treated_3"]
    )

    with pytest.raises(DESeq2Error):
        validate_comparison(
            metadata=metadata,
            condition_column="condition",
            reference_condition="control",
            comparison_condition="treated",
        )


def test_mismatched_samples_are_rejected():
    counts = make_counts().drop(
        columns="treated_3"
    )

    with pytest.raises(DESeq2Error):
        run_deseq2(
            counts=counts,
            metadata=make_metadata(),
            condition_column="condition",
            reference_condition="control",
            comparison_condition="treated",
        )