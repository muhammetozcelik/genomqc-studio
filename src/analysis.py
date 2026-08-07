from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


@dataclass(frozen=True)
class PCAResult:
    coordinates: pd.DataFrame
    explained_variance: tuple[float, float]
    variable_gene_count: int


def calculate_library_sizes(
    counts: pd.DataFrame,
) -> pd.DataFrame:
    library_sizes = counts.sum(axis=0)

    return pd.DataFrame(
        {
            "sample_id": library_sizes.index,
            "total_counts": library_sizes.values,
        }
    )


def normalize_log_cpm(
    counts: pd.DataFrame,
    pseudocount: float = 1.0,
) -> pd.DataFrame:
    library_sizes = counts.sum(axis=0)

    if (library_sizes <= 0).any():
        empty_samples = library_sizes[
            library_sizes <= 0
        ].index.tolist()

        raise ValueError(
            "Cannot normalize empty libraries: "
            + ", ".join(empty_samples)
        )

    counts_per_million = counts.divide(
        library_sizes,
        axis="columns",
    ) * 1_000_000

    return np.log2(counts_per_million + pseudocount)


def run_pca(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    top_variable_genes: int = 500,
) -> PCAResult:
    if set(counts.columns) != set(metadata.index):
        raise ValueError(
            "Count matrix and metadata samples do not match."
        )

    normalized = normalize_log_cpm(counts)

    gene_variances = normalized.var(axis=1)
    variable_genes = (
        gene_variances[gene_variances > 0]
        .sort_values(ascending=False)
        .head(top_variable_genes)
        .index
    )

    if len(variable_genes) < 2:
        raise ValueError(
            "At least two variable genes are required for PCA."
        )

    analysis_matrix = normalized.loc[
        variable_genes,
        counts.columns,
    ].transpose()

    model = PCA(n_components=2)
    principal_components = model.fit_transform(
        analysis_matrix
    )

    coordinates = pd.DataFrame(
        principal_components,
        columns=["PC1", "PC2"],
        index=counts.columns,
    )
    coordinates.index.name = "sample_id"

    ordered_metadata = metadata.loc[counts.columns]
    coordinates = coordinates.join(ordered_metadata)
    coordinates = coordinates.reset_index()

    explained_variance = tuple(
        model.explained_variance_ratio_[:2] * 100
    )

    return PCAResult(
        coordinates=coordinates,
        explained_variance=explained_variance,
        variable_gene_count=len(variable_genes),
    )


def calculate_sample_correlations(
    counts: pd.DataFrame,
) -> pd.DataFrame:
    normalized = normalize_log_cpm(counts)
    return normalized.corr(method="pearson")