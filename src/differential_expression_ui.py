import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from src.differential_expression import (
    DESeq2Error,
    run_deseq2,
)


DIRECTION_COLORS = {
    "Upregulated": "#35D0A6",
    "Downregulated": "#FF6B81",
    "Not significant": "#647A74",
}


def _dataset_fingerprint(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
) -> str:
    count_hash = pd.util.hash_pandas_object(
        counts,
        index=True,
    ).sum()

    metadata_hash = pd.util.hash_pandas_object(
        metadata,
        index=True,
    ).sum()

    return f"{int(count_hash)}:{int(metadata_hash)}"


def _prepare_volcano_data(
    results: pd.DataFrame,
) -> pd.DataFrame:
    volcano_data = results.copy()

    smallest_value = np.finfo(float).tiny
    adjusted_pvalues = volcano_data["padj"].clip(
        lower=smallest_value
    )

    volcano_data["minus_log10_padj"] = -np.log10(
        adjusted_pvalues
    )

    return volcano_data


def _apply_significance_thresholds(
    results: pd.DataFrame,
    adjusted_pvalue_threshold: float,
    log2_fold_change_threshold: float,
) -> pd.DataFrame:
    if not 0 <= adjusted_pvalue_threshold <= 1:
        raise ValueError(
            "Adjusted p-value threshold must be between 0 and 1."
        )

    if log2_fold_change_threshold < 0:
        raise ValueError(
            "Log2 fold-change threshold cannot be negative."
        )

    classified_results = results.copy()

    adjusted_pvalues = classified_results["padj"].fillna(1.0)
    fold_changes = classified_results["log2FoldChange"]

    passes_pvalue_threshold = (
        adjusted_pvalues <= adjusted_pvalue_threshold
    )
    passes_fold_change_threshold = (
        fold_changes.abs() >= log2_fold_change_threshold
    )
    has_direction = fold_changes != 0

    significant = (
        passes_pvalue_threshold
        & passes_fold_change_threshold
        & has_direction
    )

    classified_results["significant"] = significant
    classified_results["direction"] = "Not significant"

    classified_results.loc[
        significant & (fold_changes > 0),
        "direction",
    ] = "Upregulated"

    classified_results.loc[
        significant & (fold_changes < 0),
        "direction",
    ] = "Downregulated"

    return classified_results


def _render_volcano_plot(
    results: pd.DataFrame,
    log2_fold_change_threshold: float,
    adjusted_pvalue_threshold: float,
):
    volcano_data = _prepare_volcano_data(results)

    figure = px.scatter(
        volcano_data,
        x="log2FoldChange",
        y="minus_log10_padj",
        color="direction",
        color_discrete_map=DIRECTION_COLORS,
        hover_name="gene_id",
        hover_data={
            "baseMean": ":.2f",
            "log2FoldChange": ":.3f",
            "padj": ":.3e",
            "minus_log10_padj": False,
        },
        labels={
            "log2FoldChange": "Log2 fold change",
            "minus_log10_padj": "−Log10 adjusted p-value",
            "direction": "Result",
        },
        title="Differential-expression volcano plot",
    )

    if log2_fold_change_threshold > 0:
        figure.add_vline(
            x=-log2_fold_change_threshold,
            line_dash="dash",
            line_color="#A9BBB6",
        )
        figure.add_vline(
            x=log2_fold_change_threshold,
            line_dash="dash",
            line_color="#A9BBB6",
        )

    if adjusted_pvalue_threshold > 0:
        figure.add_hline(
            y=-np.log10(adjusted_pvalue_threshold),
            line_dash="dash",
            line_color="#A9BBB6",
        )

    figure.update_traces(
        marker={
            "size": 7,
            "opacity": 0.78,
        }
    )
    figure.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    st.plotly_chart(
        figure,
        use_container_width=True,
    )


def _render_ma_plot(
    results: pd.DataFrame,
    log2_fold_change_threshold: float,
):
    plot_data = results.copy()
    plot_data["log10_base_mean"] = np.log10(
        plot_data["baseMean"] + 1
    )

    figure = px.scatter(
        plot_data,
        x="log10_base_mean",
        y="log2FoldChange",
        color="direction",
        color_discrete_map=DIRECTION_COLORS,
        hover_name="gene_id",
        hover_data={
            "baseMean": ":.2f",
            "log2FoldChange": ":.3f",
            "padj": ":.3e",
            "log10_base_mean": False,
        },
        labels={
            "log10_base_mean": "Log10 mean expression",
            "log2FoldChange": "Log2 fold change",
            "direction": "Result",
        },
        title="MA plot",
    )

    figure.add_hline(
        y=0,
        line_color="#A9BBB6",
    )

    if log2_fold_change_threshold > 0:
        figure.add_hline(
            y=-log2_fold_change_threshold,
            line_dash="dash",
            line_color="#A9BBB6",
        )
        figure.add_hline(
            y=log2_fold_change_threshold,
            line_dash="dash",
            line_color="#A9BBB6",
        )

    figure.update_traces(
        marker={
            "size": 7,
            "opacity": 0.78,
        }
    )
    figure.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    st.plotly_chart(
        figure,
        use_container_width=True,
    )


def _render_expression_heatmap(
    results: pd.DataFrame,
    vst_counts: pd.DataFrame,
):
    top_genes = (
        results.loc[results["significant"]]
        .sort_values("padj")
        .head(30)["gene_id"]
        .tolist()
    )

    if len(top_genes) < 2:
        st.info(
            "At least two significant genes are required "
            "to create the expression heatmap."
        )
        return

    expression_matrix = (
        vst_counts.set_index("gene_id")
        .loc[top_genes]
        .astype(float)
    )

    row_means = expression_matrix.mean(axis=1)
    row_standard_deviations = (
        expression_matrix.std(axis=1).replace(0, 1)
    )

    scaled_expression = expression_matrix.sub(
        row_means,
        axis=0,
    ).div(
        row_standard_deviations,
        axis=0,
    )

    figure = px.imshow(
        scaled_expression,
        aspect="auto",
        color_continuous_scale=[
            "#315CFF",
            "#F4FAF8",
            "#FF5F6D",
        ],
        color_continuous_midpoint=0,
        labels={
            "x": "Sample",
            "y": "Gene",
            "color": "Row Z-score",
        },
        title="Top differential genes",
    )

    figure.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        height=max(500, len(top_genes) * 20),
    )

    st.plotly_chart(
        figure,
        use_container_width=True,
    )


def render_differential_expression(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    validation_report,
):
    st.subheader("Differential expression")
    st.caption(
        "DESeq2 analysis with Benjamini–Hochberg "
        "multiple-testing correction."
    )

    if not validation_report.is_valid:
        st.warning(
            "Differential expression is disabled until "
            "blocking validation errors are corrected."
        )
        return

    conditions = (
        metadata["condition"]
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    if len(conditions) < 2:
        st.error(
            "At least two experimental conditions are required."
        )
        return

    selection_columns = st.columns(2)

    reference_condition = selection_columns[0].selectbox(
        "Reference condition",
        options=conditions,
    )

    comparison_options = [
        condition
        for condition in conditions
        if condition != reference_condition
    ]

    comparison_condition = selection_columns[1].selectbox(
        "Comparison condition",
        options=comparison_options,
    )

    threshold_columns = st.columns(2)

    adjusted_pvalue_threshold = threshold_columns[0].slider(
        "Adjusted p-value (FDR) threshold",
        min_value=0.01,
        max_value=0.10,
        value=0.05,
        step=0.01,
        help=(
            "Genes must have an adjusted p-value at or below "
            "this threshold."
        ),
    )

    log2_fold_change_threshold = threshold_columns[1].slider(
        "Absolute log2 fold-change threshold",
        min_value=0.0,
        max_value=5.0,
        value=1.0,
        step=0.1,
        help=(
            "Genes must meet or exceed this absolute "
            "log2 fold-change threshold."
        ),
    )

    st.caption(
        "Positive fold changes indicate higher expression in "
        f"**{comparison_condition}** relative to "
        f"**{reference_condition}**."
    )

    fingerprint = _dataset_fingerprint(
        counts,
        metadata,
    )
    analysis_key = (
        f"{fingerprint}:"
        f"{reference_condition}:"
        f"{comparison_condition}"
    )

    run_analysis = st.button(
        "Run DESeq2 analysis",
        type="primary",
        use_container_width=True,
    )

    if run_analysis:
        try:
            with st.spinner(
                "Running DESeq2 differential-expression analysis..."
            ):
                result = run_deseq2(
                    counts=counts,
                    metadata=metadata,
                    condition_column="condition",
                    reference_condition=reference_condition,
                    comparison_condition=comparison_condition,
                )

            st.session_state["deseq2_analysis"] = {
                "key": analysis_key,
                "result": result,
            }

        except DESeq2Error as error:
            st.error(f"DESeq2 analysis failed: {error}")
            return

    stored_analysis = st.session_state.get(
        "deseq2_analysis"
    )

    if (
        stored_analysis is None
        or stored_analysis["key"] != analysis_key
    ):
        st.info(
            "Choose the comparison and run DESeq2 to "
            "generate statistical results."
        )
        return

    result = stored_analysis["result"]

    result_table = _apply_significance_thresholds(
        results=result.results,
        adjusted_pvalue_threshold=adjusted_pvalue_threshold,
        log2_fold_change_threshold=log2_fold_change_threshold,
    )

    upregulated_count = int(
        (
            result_table["direction"]
            == "Upregulated"
        ).sum()
    )
    downregulated_count = int(
        (
            result_table["direction"]
            == "Downregulated"
        ).sum()
    )
    significant_count = (
        upregulated_count + downregulated_count
    )

    metric_columns = st.columns(4)
    metric_columns[0].metric(
        "Analyzed genes",
        f"{len(result_table):,}",
    )
    metric_columns[1].metric(
        "Significant",
        f"{significant_count:,}",
    )
    metric_columns[2].metric(
        "Upregulated",
        f"{upregulated_count:,}",
    )
    metric_columns[3].metric(
        "Downregulated",
        f"{downregulated_count:,}",
    )

    volcano_tab, ma_tab, heatmap_tab, results_tab = st.tabs(
        [
            "Volcano plot",
            "MA plot",
            "Expression heatmap",
            "Results table",
        ]
    )

    with volcano_tab:
        _render_volcano_plot(
            results=result_table,
            log2_fold_change_threshold=(
                log2_fold_change_threshold
            ),
            adjusted_pvalue_threshold=(
                adjusted_pvalue_threshold
            ),
        )

    with ma_tab:
        _render_ma_plot(
            results=result_table,
            log2_fold_change_threshold=(
                log2_fold_change_threshold
            ),
        )

    with heatmap_tab:
        _render_expression_heatmap(
            result_table,
            result.vst_counts,
        )

    with results_tab:
        significant_results = (
            result_table.loc[result_table["significant"]]
            .sort_values("padj")
        )

        st.dataframe(
            significant_results,
            use_container_width=True,
            hide_index=True,
            column_config={
                "baseMean": st.column_config.NumberColumn(
                    format="%.2f"
                ),
                "log2FoldChange": (
                    st.column_config.NumberColumn(
                        format="%.3f"
                    )
                ),
                "padj": st.column_config.NumberColumn(
                    format="%.3e"
                ),
            },
        )

        download_columns = [
            "gene_id",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
            "direction",
        ]

        st.download_button(
            label="Download complete DESeq2 results",
            data=result_table[
                download_columns
            ].to_csv(index=False),
            file_name=(
                f"deseq2_{comparison_condition}"
                f"_vs_{reference_condition}.csv"
            ),
            mime="text/csv",
        )