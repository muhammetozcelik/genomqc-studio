from io import StringIO

import pytest

from src.data_loader import (
    DataFormatError,
    load_count_matrix,
    load_metadata,
)


def test_load_valid_count_matrix():
    content = StringIO(
        "gene_id,control_1,treated_1\n"
        "gene_a,10,20\n"
        "gene_b,5,15\n"
    )

    counts = load_count_matrix(content)

    assert list(counts.columns) == ["control_1", "treated_1"]
    assert list(counts.index) == ["gene_a", "gene_b"]
    assert counts.loc["gene_a", "treated_1"] == 20


def test_load_valid_metadata():
    content = StringIO(
        "sample_id\tcondition\treplicate\n"
        "control_1\tcontrol\t1\n"
        "treated_1\ttreated\t1\n"
    )

    metadata = load_metadata(content)

    assert list(metadata.index) == ["control_1", "treated_1"]
    assert metadata.loc["treated_1", "condition"] == "treated"


def test_negative_count_is_rejected():
    content = StringIO(
        "gene_id,control_1,treated_1\n"
        "gene_a,10,-1\n"
    )

    with pytest.raises(DataFormatError):
        load_count_matrix(content)


def test_decimal_count_is_rejected():
    content = StringIO(
        "gene_id,control_1,treated_1\n"
        "gene_a,10,4.5\n"
    )

    with pytest.raises(DataFormatError):
        load_count_matrix(content)


def test_missing_metadata_column_is_rejected():
    content = StringIO(
        "sample_id,group\n"
        "sample_1,control\n"
    )

    with pytest.raises(DataFormatError):
        load_metadata(content)