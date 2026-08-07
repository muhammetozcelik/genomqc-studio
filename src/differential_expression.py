from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
import tempfile

import pandas as pd


class DESeq2Error(RuntimeError):
    """Raised when a DESeq2 analysis cannot be completed."""


@dataclass(frozen=True)
class DifferentialExpressionResult:
    results: pd.DataFrame
    normalized_counts: pd.DataFrame
    vst_counts: pd.DataFrame
    summary: dict[str, str]
    reference_condition: str
    comparison_condition: str


def validate_comparison(
    metadata: pd.DataFrame,
    condition_column: str,
    reference_condition: str,
    comparison_condition: str,
) -> None:
    if condition_column not in metadata.columns:
        raise DESeq2Error(
            f"Metadata column not found: {condition_column}"
        )

    condition_values = (
        metadata[condition_column]
        .astype(str)
        .str.strip()
    )

    available_conditions = set(condition_values)

    if reference_condition not in available_conditions:
        raise DESeq2Error(
            "Reference condition not found: "
            f"{reference_condition}"
        )

    if comparison_condition not in available_conditions:
        raise DESeq2Error(
            "Comparison condition not found: "
            f"{comparison_condition}"
        )

    if reference_condition == comparison_condition:
        raise DESeq2Error(
            "Reference and comparison conditions must differ."
        )

    selected_counts = condition_values.value_counts()

    for condition in [
        reference_condition,
        comparison_condition,
    ]:
        replicate_count = int(
            selected_counts.get(condition, 0)
        )

        if replicate_count < 2:
            raise DESeq2Error(
                f"Condition '{condition}' has only "
                f"{replicate_count} sample(s). "
                "At least two replicates are required."
            )


def _read_output_table(
    output_directory: Path,
    filename: str,
) -> pd.DataFrame:
    output_path = output_directory / filename

    if not output_path.exists():
        raise DESeq2Error(
            f"Expected DESeq2 output was not created: "
            f"{filename}"
        )

    return pd.read_csv(output_path)


def run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    condition_column: str,
    reference_condition: str,
    comparison_condition: str,
    timeout_seconds: int = 300,
) -> DifferentialExpressionResult:
    validate_comparison(
        metadata=metadata,
        condition_column=condition_column,
        reference_condition=reference_condition,
        comparison_condition=comparison_condition,
    )

    if set(counts.columns) != set(metadata.index):
        raise DESeq2Error(
            "Count matrix and metadata samples do not match."
        )

    rscript_executable = shutil.which("Rscript")

    if rscript_executable is None:
        raise DESeq2Error(
            "Rscript was not found in the active environment."
        )

    runner_path = Path(__file__).with_name(
        "run_deseq2.R"
    )

    if not runner_path.exists():
        raise DESeq2Error(
            f"DESeq2 runner was not found: {runner_path}"
        )

    ordered_metadata = metadata.loc[counts.columns].copy()

    with tempfile.TemporaryDirectory(
        prefix="genomqc_deseq2_"
    ) as temporary_directory:
        working_directory = Path(temporary_directory)

        counts_path = working_directory / "counts.csv"
        metadata_path = working_directory / "metadata.csv"
        output_directory = working_directory / "output"

        output_directory.mkdir(parents=True)

        counts_export = counts.copy()
        counts_export.index.name = "gene_id"
        counts_export.reset_index().to_csv(
            counts_path,
            index=False,
        )

        metadata_export = ordered_metadata.copy()
        metadata_export.index.name = "sample_id"
        metadata_export.reset_index().to_csv(
            metadata_path,
            index=False,
        )

        command = [
            rscript_executable,
            str(runner_path),
            str(counts_path),
            str(metadata_path),
            condition_column,
            reference_condition,
            comparison_condition,
            str(output_directory),
        ]

        try:
            completed_process = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            raise DESeq2Error(
                "DESeq2 analysis exceeded the time limit."
            ) from error

        if completed_process.returncode != 0:
            error_message = (
                completed_process.stderr.strip()
                or completed_process.stdout.strip()
                or "Unknown DESeq2 error."
            )

            raise DESeq2Error(error_message)

        results = _read_output_table(
            output_directory,
            "deseq2_results.csv",
        )
        normalized_counts = _read_output_table(
            output_directory,
            "normalized_counts.csv",
        )
        vst_counts = _read_output_table(
            output_directory,
            "vst_counts.csv",
        )

        summary_path = output_directory / "run_summary.tsv"

        if not summary_path.exists():
            raise DESeq2Error(
                "DESeq2 run summary was not created."
            )

        summary_table = pd.read_csv(
            summary_path,
            sep="\t",
            dtype=str,
        )

        summary = dict(
            zip(
                summary_table["metric"],
                summary_table["value"],
            )
        )

    results["significant"] = (
        results["padj"].notna()
        & (results["padj"] < 0.05)
        & (results["log2FoldChange"].abs() >= 1)
    )

    results["direction"] = "Not significant"
    results.loc[
        results["significant"]
        & (results["log2FoldChange"] >= 1),
        "direction",
    ] = "Upregulated"
    results.loc[
        results["significant"]
        & (results["log2FoldChange"] <= -1),
        "direction",
    ] = "Downregulated"

    return DifferentialExpressionResult(
        results=results,
        normalized_counts=normalized_counts,
        vst_counts=vst_counts,
        summary=summary,
        reference_condition=reference_condition,
        comparison_condition=comparison_condition,
    )