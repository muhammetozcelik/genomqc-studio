import json
from pathlib import Path

import plotly.express as px
import streamlit as st

from src.analysis import (
    calculate_library_sizes,
    calculate_sample_correlations,
    run_pca,
)
from src.data_loader import (
    DataFormatError,
    load_count_matrix,
    load_metadata,
)
from src.validator import validate_dataset


ROOT_DIRECTORY = Path(__file__).resolve().parent
DEMO_COUNTS = ROOT_DIRECTORY / "demo_data" / "counts.csv"
DEMO_METADATA = ROOT_DIRECTORY / "demo_data" / "metadata.csv"

COLORS = [
    "#35D0A6",
    "#FFB86B",
    "#6EA8FE",
    "#D98BFF",
    "#FF6B81",
]


st.set_page_config(
    page_title="GenomQC Studio",
    page_icon="🧬",
    layout="wide",
)

st.markdown(
    """
    <style>
        .stApp {
            background:
                radial-gradient(
                    circle at top right,
                    rgba(53, 208, 166, 0.10),
                    transparent 30%
                ),
                #071916;
        }

        [data-testid="stSidebar"] {
            background-color: #0B2420;
        }

        h1, h2, h3 {
            letter-spacing: -0.03em;
        }

        .hero-label {
            color: #35D0A6;
            font-size: 0.82rem;
            font-weight: 800;
            letter-spacing: 0.14em;
            text-transform: uppercase;
        }

        .hero-title {
            color: #F4FAF8;
            font-size: clamp(2.6rem, 6vw, 5.6rem);
            font-weight: 800;
            line-height: 0.95;
            letter-spacing: -0.065em;
            margin: 0.4rem 0 1rem 0;
        }

        .hero-copy {
            color: #B8CCC7;
            font-size: 1.12rem;
            max-width: 760px;
            line-height: 1.65;
        }

        div[data-testid="stMetric"] {
            background: rgba(16, 54, 47, 0.75);
            border: 1px solid rgba(108, 195, 172, 0.22);
            border-radius: 16px;
            padding: 1rem;
        }

        div[data-testid="stMetric"] label {
            color: #9DB8B1;
        }

        .source-badge {
            display: inline-block;
            margin-top: 0.8rem;
            padding: 0.35rem 0.7rem;
            color: #71E1C3;
            background: rgba(53, 208, 166, 0.12);
            border: 1px solid rgba(53, 208, 166, 0.28);
            border-radius: 999px;
            font-size: 0.82rem;
            font-weight: 700;
        }
    </style>
    """,
    unsafe_allow_html=True,
)


def display_issue(issue):
    text = f"**{issue.code}** — {issue.message}"

    if issue.severity == "error":
        st.error(text)
    elif issue.severity == "warning":
        st.warning(text)
    else:
        st.info(text)


def prepare_passport(
    counts,
    metadata,
    report,
    source_name,
):
    library_sizes = calculate_library_sizes(counts)

    return {
        "application": "GenomQC Studio",
        "dataset_source": source_name,
        "validation": {
            "passed": report.is_valid,
            "score": report.score,
            "errors": len(report.errors),
            "warnings": len(report.warnings),
        },
        "dataset": {
            "genes": int(counts.shape[0]),
            "samples": int(counts.shape[1]),
            "conditions": int(
                metadata["condition"].nunique()
            ),
        },
        "library_sizes": {
            row["sample_id"]: int(row["total_counts"])
            for row in library_sizes.to_dict("records")
        },
        "issues": [
            {
                "severity": issue.severity,
                "code": issue.code,
                "message": issue.message,
            }
            for issue in report.issues
        ],
    }


with st.sidebar:
    st.markdown("## GenomQC Studio")
    st.caption("Transcriptomics dataset intelligence")

    use_demo = st.toggle(
        "Use simulated demo dataset",
        value=True,
    )

    count_file = None
    metadata_file = None

    if not use_demo:
        count_file = st.file_uploader(
            "Count matrix",
            type=["csv", "tsv", "txt"],
        )
        metadata_file = st.file_uploader(
            "Sample metadata",
            type=["csv", "tsv", "txt"],
        )

    st.divider()

    st.markdown("### Required structure")

    st.code(
        "gene_id,control_1,treated_1\n"
        "GENE_0001,120,340",
        language="text",
    )

    st.code(
        "sample_id,condition\n"
        "control_1,control\n"
        "treated_1,treated",
        language="text",
    )

    st.caption(
        "Uploaded files remain in the current application "
        "session and are not written to the project directory."
    )


st.markdown(
    """
    <div class="hero-label">RNA-seq dataset intelligence</div>
    <div class="hero-title">
        Know your data<br>
        before you trust the result.
    </div>
    <div class="hero-copy">
        Validate count matrices and experimental metadata,
        identify design risks, inspect library composition and
        explore sample relationships before differential
        expression analysis.
    </div>
    """,
    unsafe_allow_html=True,
)


if use_demo:
    source_name = "Simulated demonstration dataset"
    st.markdown(
        '<div class="source-badge">'
        'SIMULATED DEMO · 3 CONTROL · 3 TREATED'
        "</div>",
        unsafe_allow_html=True,
    )

    counts_source = DEMO_COUNTS
    metadata_source = DEMO_METADATA
else:
    source_name = "User-uploaded dataset"

    if count_file is None or metadata_file is None:
        st.info(
            "Upload both a count matrix and a metadata file "
            "from the sidebar to begin."
        )
        st.stop()

    count_file.seek(0)
    metadata_file.seek(0)

    counts_source = count_file
    metadata_source = metadata_file


try:
    counts = load_count_matrix(counts_source)
    metadata = load_metadata(metadata_source)
except DataFormatError as error:
    st.error(f"Dataset could not be loaded: {error}")
    st.stop()
except Exception as error:
    st.error(f"Unexpected loading error: {error}")
    st.stop()


report = validate_dataset(counts, metadata)

st.write("")

metric_columns = st.columns(4)

metric_columns[0].metric(
    "Genes",
    f"{counts.shape[0]:,}",
)
metric_columns[1].metric(
    "Samples",
    counts.shape[1],
)
metric_columns[2].metric(
    "Conditions",
    metadata["condition"].nunique(),
)
metric_columns[3].metric(
    "Validation score",
    f"{report.score}/100",
)

overview_tab, audit_tab, exploration_tab, data_tab = st.tabs(
    [
        "Overview",
        "Quality audit",
        "Exploration",
        "Data preview",
    ]
)


with overview_tab:
    st.subheader("Dataset overview")

    library_sizes = calculate_library_sizes(counts)
    metadata_table = metadata.reset_index()

    library_table = library_sizes.merge(
        metadata_table,
        on="sample_id",
        how="left",
    )

    library_figure = px.bar(
        library_table,
        x="sample_id",
        y="total_counts",
        color="condition",
        color_discrete_sequence=COLORS,
        labels={
            "sample_id": "Sample",
            "total_counts": "Total counts",
            "condition": "Condition",
        },
        title="Library sizes",
    )

    library_figure.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend_title_text="Condition",
    )

    st.plotly_chart(
        library_figure,
        use_container_width=True,
    )

    condition_counts = (
        metadata["condition"]
        .value_counts()
        .rename_axis("condition")
        .reset_index(name="sample_count")
    )

    st.markdown("#### Experimental groups")
    st.dataframe(
        condition_counts,
        use_container_width=True,
        hide_index=True,
    )


with audit_tab:
    st.subheader("Pre-analysis quality audit")
    st.progress(report.score / 100)

    if report.is_valid:
        st.success(
            "The dataset passed all blocking validation checks."
        )
    else:
        st.error(
            "Blocking issues must be corrected before analysis."
        )

    for issue in report.issues:
        display_issue(issue)

    passport = prepare_passport(
        counts,
        metadata,
        report,
        source_name,
    )

    st.download_button(
        label="Download analysis passport",
        data=json.dumps(passport, indent=2),
        file_name="genomqc_analysis_passport.json",
        mime="application/json",
    )


with exploration_tab:
    st.subheader("Sample-level exploration")

    if not report.is_valid:
        st.warning(
            "Exploratory analysis is disabled until blocking "
            "validation errors are corrected."
        )
    else:
        try:
            pca_result = run_pca(counts, metadata)

            pc1_variance = pca_result.explained_variance[0]
            pc2_variance = pca_result.explained_variance[1]

            scatter_arguments = {
                "data_frame": pca_result.coordinates,
                "x": "PC1",
                "y": "PC2",
                "color": "condition",
                "hover_name": "sample_id",
                "color_discrete_sequence": COLORS,
                "labels": {
                    "PC1": f"PC1 ({pc1_variance:.1f}%)",
                    "PC2": f"PC2 ({pc2_variance:.1f}%)",
                    "condition": "Condition",
                },
                "title": (
                    "PCA of the most variable genes "
                    f"(n={pca_result.variable_gene_count})"
                ),
            }

            if "batch" in pca_result.coordinates.columns:
                scatter_arguments["symbol"] = "batch"

            pca_figure = px.scatter(**scatter_arguments)
            pca_figure.update_traces(
                marker={"size": 14, "line": {"width": 1}}
            )
            pca_figure.update_layout(
                template="plotly_dark",
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
            )

            st.plotly_chart(
                pca_figure,
                use_container_width=True,
            )

            correlations = calculate_sample_correlations(
                counts
            )

            correlation_figure = px.imshow(
                correlations,
                text_auto=".2f",
                aspect="auto",
                color_continuous_scale=[
                    "#102E29",
                    "#35D0A6",
                    "#F4FAF8",
                ],
                zmin=0,
                zmax=1,
                title="Pearson sample correlation",
            )
            correlation_figure.update_layout(
                template="plotly_dark",
                paper_bgcolor="rgba(0,0,0,0)",
            )

            st.plotly_chart(
                correlation_figure,
                use_container_width=True,
            )

        except Exception as error:
            st.error(f"Exploration failed: {error}")


with data_tab:
    st.subheader("Count matrix")
    st.caption("First 30 genes")
    st.dataframe(
        counts.head(30),
        use_container_width=True,
    )

    st.subheader("Sample metadata")
    st.dataframe(
        metadata.reset_index(),
        use_container_width=True,
        hide_index=True,
    )