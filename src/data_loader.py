from pathlib import Path
from typing import BinaryIO, TextIO

import pandas as pd


TableSource = str | Path | BinaryIO | TextIO


class DataFormatError(ValueError):
    """Raised when an uploaded table has an invalid structure."""


def read_table(source: TableSource) -> pd.DataFrame:
    """Read a CSV or TSV table by automatically detecting its separator."""
    if hasattr(source, "seek"):
        source.seek(0)

    try:
        table = pd.read_csv(source, sep=None, engine="python")
    except Exception as error:
        raise DataFormatError(f"Could not read the table: {error}") from error

    if table.empty:
        raise DataFormatError("The uploaded table is empty.")

    table.columns = [str(column).strip() for column in table.columns]
    return table


def load_count_matrix(source: TableSource) -> pd.DataFrame:
    """
    Load a gene-count matrix.

    Expected structure:
    gene_id, sample_1, sample_2, ...
    """
    table = read_table(source)

    if table.shape[1] < 3:
        raise DataFormatError(
            "The count matrix must contain a gene column "
            "and at least two sample columns."
        )

    gene_column = table.columns[0]
    table = table.rename(columns={gene_column: "gene_id"})

    table["gene_id"] = table["gene_id"].astype(str).str.strip()

    if table["gene_id"].eq("").any():
        raise DataFormatError("Gene identifiers cannot be empty.")

    if table["gene_id"].duplicated().any():
        duplicated = table.loc[
            table["gene_id"].duplicated(), "gene_id"
        ].iloc[0]
        raise DataFormatError(
            f"Gene identifiers must be unique. Duplicate: {duplicated}"
        )

    sample_columns = list(table.columns[1:])

    if len(sample_columns) != len(set(sample_columns)):
        raise DataFormatError("Sample names must be unique.")

    for sample in sample_columns:
        table[sample] = pd.to_numeric(table[sample], errors="coerce")

    if table[sample_columns].isna().any().any():
        raise DataFormatError(
            "All count values must be numeric and cannot be missing."
        )

    if (table[sample_columns] < 0).any().any():
        raise DataFormatError("Count values cannot be negative.")

    non_integer = (table[sample_columns] % 1 != 0).any().any()
    if non_integer:
        raise DataFormatError("Count values must be whole numbers.")

    table[sample_columns] = table[sample_columns].astype("int64")
    return table.set_index("gene_id")


def load_metadata(source: TableSource) -> pd.DataFrame:
    """
    Load sample metadata.

    Required columns:
    sample_id, condition
    """
    table = read_table(source)

    required_columns = {"sample_id", "condition"}
    missing_columns = required_columns.difference(table.columns)

    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise DataFormatError(
            f"Metadata is missing required columns: {missing}"
        )

    table["sample_id"] = table["sample_id"].astype(str).str.strip()
    table["condition"] = table["condition"].astype(str).str.strip()

    if table["sample_id"].eq("").any():
        raise DataFormatError("Sample identifiers cannot be empty.")

    if table["condition"].eq("").any():
        raise DataFormatError("Condition values cannot be empty.")

    if table["sample_id"].duplicated().any():
        duplicated = table.loc[
            table["sample_id"].duplicated(), "sample_id"
        ].iloc[0]
        raise DataFormatError(
            f"Sample identifiers must be unique. Duplicate: {duplicated}"
        )

    return table.set_index("sample_id")