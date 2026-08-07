import pandas as pd

from src.validator import validate_dataset


SAMPLES = [
    "control_1",
    "control_2",
    "control_3",
    "treated_1",
    "treated_2",
    "treated_3",
]


def make_counts():
    genes = [f"gene_{number}" for number in range(120)]

    return pd.DataFrame(
        {
            sample: [100] * len(genes)
            for sample in SAMPLES
        },
        index=genes,
    )


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


def issue_codes(issues):
    return {issue.code for issue in issues}


def test_valid_dataset_passes():
    report = validate_dataset(make_counts(), make_metadata())

    assert report.is_valid
    assert report.score == 100
    assert "VALIDATION_PASSED" in issue_codes(report.issues)


def test_missing_metadata_sample_is_error():
    metadata = make_metadata().drop(index="treated_3")

    report = validate_dataset(make_counts(), metadata)

    assert not report.is_valid
    assert "SAMPLES_MISSING_FROM_METADATA" in issue_codes(report.errors)


def test_single_replicate_is_error():
    counts = make_counts().drop(
        columns=["treated_2", "treated_3"]
    )
    metadata = make_metadata().drop(
        index=["treated_2", "treated_3"]
    )

    report = validate_dataset(counts, metadata)

    assert not report.is_valid
    assert "INSUFFICIENT_REPLICATION" in issue_codes(report.errors)


def test_low_library_size_is_warning():
    counts = make_counts()
    counts["control_1"] = 1

    report = validate_dataset(counts, make_metadata())

    assert "LOW_LIBRARY_SIZE" in issue_codes(report.warnings)


def test_empty_library_is_error():
    counts = make_counts()
    counts["control_1"] = 0

    report = validate_dataset(counts, make_metadata())

    assert not report.is_valid
    assert "EMPTY_LIBRARIES" in issue_codes(report.errors)