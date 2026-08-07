from dataclasses import dataclass, field
from typing import Literal

import pandas as pd


Severity = Literal["error", "warning", "info"]


@dataclass(frozen=True)
class ValidationIssue:
    severity: Severity
    code: str
    message: str


@dataclass
class ValidationReport:
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def errors(self) -> list[ValidationIssue]:
        return [
            issue for issue in self.issues
            if issue.severity == "error"
        ]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [
            issue for issue in self.issues
            if issue.severity == "warning"
        ]

    @property
    def is_valid(self) -> bool:
        return not self.errors

    @property
    def score(self) -> int:
        penalty = len(self.errors) * 25 + len(self.warnings) * 8
        return max(0, 100 - penalty)

    def add(
        self,
        severity: Severity,
        code: str,
        message: str,
    ) -> None:
        self.issues.append(
            ValidationIssue(
                severity=severity,
                code=code,
                message=message,
            )
        )


def validate_dataset(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
) -> ValidationReport:
    """Validate count-matrix structure and experimental design."""
    report = ValidationReport()

    count_samples = set(counts.columns)
    metadata_samples = set(metadata.index)

    missing_metadata = sorted(count_samples - metadata_samples)
    if missing_metadata:
        report.add(
            "error",
            "SAMPLES_MISSING_FROM_METADATA",
            "Count-matrix samples missing from metadata: "
            + ", ".join(missing_metadata),
        )

    missing_counts = sorted(metadata_samples - count_samples)
    if missing_counts:
        report.add(
            "error",
            "SAMPLES_MISSING_FROM_COUNTS",
            "Metadata samples missing from the count matrix: "
            + ", ".join(missing_counts),
        )

    if counts.shape[0] < 100:
        report.add(
            "warning",
            "LOW_GENE_COUNT",
            (
                f"The matrix contains only {counts.shape[0]} genes. "
                "Confirm that this is an unfiltered gene-level matrix."
            ),
        )

    all_zero_genes = int((counts.sum(axis=1) == 0).sum())
    if all_zero_genes:
        percentage = all_zero_genes / counts.shape[0] * 100
        report.add(
            "info",
            "ALL_ZERO_GENES",
            (
                f"{all_zero_genes:,} genes ({percentage:.1f}%) "
                "have zero counts across every sample."
            ),
        )

    library_sizes = counts.sum(axis=0)

    zero_libraries = library_sizes[library_sizes == 0].index.tolist()
    if zero_libraries:
        report.add(
            "error",
            "EMPTY_LIBRARIES",
            "Samples with no assigned counts: "
            + ", ".join(zero_libraries),
        )

    positive_libraries = library_sizes[library_sizes > 0]
    if not positive_libraries.empty:
        median_library = float(positive_libraries.median())
        low_threshold = median_library * 0.25

        low_libraries = library_sizes[
            (library_sizes > 0) & (library_sizes < low_threshold)
        ]

        for sample, size in low_libraries.items():
            report.add(
                "warning",
                "LOW_LIBRARY_SIZE",
                (
                    f"Sample {sample} has {int(size):,} counts, "
                    f"below 25% of the median library size "
                    f"({median_library:,.0f})."
                ),
            )

    shared_samples = [
        sample
        for sample in counts.columns
        if sample in metadata.index
    ]

    if shared_samples:
        aligned_metadata = metadata.loc[shared_samples]

        condition_count = aligned_metadata["condition"].nunique()
        if condition_count < 2:
            report.add(
                "error",
                "SINGLE_CONDITION",
                "Differential analysis requires at least two conditions.",
            )

        replicate_counts = (
            aligned_metadata["condition"]
            .value_counts()
            .sort_index()
        )

        for condition, replicate_count in replicate_counts.items():
            if replicate_count < 2:
                report.add(
                    "error",
                    "INSUFFICIENT_REPLICATION",
                    (
                        f"Condition '{condition}' has only "
                        f"{replicate_count} sample. At least two "
                        "biological replicates are required."
                    ),
                )
            elif replicate_count < 3:
                report.add(
                    "warning",
                    "LIMITED_REPLICATION",
                    (
                        f"Condition '{condition}' has only "
                        f"{replicate_count} replicates. Three or more "
                        "are recommended for stronger inference."
                    ),
                )

    if report.is_valid:
        report.add(
            "info",
            "VALIDATION_PASSED",
            "The dataset passed all blocking validation checks.",
        )

    return report