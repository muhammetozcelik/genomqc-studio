# GenomQC Studio

[![CI](https://github.com/muhammetozcelik/genomqc-studio/actions/workflows/ci.yml/badge.svg)](https://github.com/muhammetozcelik/genomqc-studio/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.12-3776AB)
![Streamlit](https://img.shields.io/badge/Streamlit-1.61-FF4B4B)

A local-first application for auditing, validating and exploring RNA-seq count datasets before differential expression analysis.

GenomQC Studio checks whether a count matrix and its experimental metadata are structurally and statistically ready for downstream analysis. It turns common dataset problems into clear, actionable findings and produces a downloadable analysis passport.

## Why this project exists

RNA-seq workflows can produce technically valid results even when the underlying experimental design contains sample mismatches, insufficient replication or unusable libraries.

GenomQC Studio adds a quality gate before differential expression analysis:

1. Load the count matrix and sample metadata.
2. Validate file structure and experimental design.
3. Identify blocking errors and statistical risks.
4. Explore library sizes, PCA and sample correlations.
5. Export a machine-readable analysis passport.

## Current capabilities

- CSV and TSV count-matrix loading
- Experimental metadata loading
- Sample identity matching
- Duplicate identifier detection
- Negative, decimal and missing-count detection
- Empty and unusually small library detection
- Experimental-condition and replication checks
- Dataset readiness score
- Library-size visualization
- Log-CPM normalization
- Principal component analysis
- Pearson sample-correlation heatmap
- Downloadable JSON analysis passport
- Deterministic simulated demonstration dataset
- Automated test suite

## Validation checks

| Check | Severity | Purpose |
|---|---|---|
| Missing metadata samples | Error | Prevents incorrect sample annotation |
| Missing count-matrix samples | Error | Detects incomplete input datasets |
| Empty libraries | Error | Blocks unusable sequencing libraries |
| Single experimental condition | Error | Prevents invalid group comparison |
| Insufficient replication | Error | Detects conditions with fewer than two samples |
| Limited replication | Warning | Highlights two-replicate experimental groups |
| Low library size | Warning | Identifies potential sequencing-depth outliers |
| Low gene count | Warning | Flags unexpectedly small count matrices |
| All-zero genes | Information | Reports uninformative features |

## Quick start

Create the reproducible environment:

```bash
micromamba create -f environment.yml
micromamba activate genomqc-studio
```

Run the automated tests:

```bash
pytest -q
```

Launch the application:

```bash
streamlit run app.py
```

Then open:

```text
http://localhost:8501
```

## Input files

### Count matrix

The first column contains unique gene identifiers. Remaining columns contain non-negative integer counts for each sample.

```csv
gene_id,control_1,control_2,treated_1,treated_2
GENE_0001,120,135,340,312
GENE_0002,52,49,18,21
```

### Sample metadata

The metadata file must contain `sample_id` and `condition` columns. Additional columns such as `replicate` and `batch` are retained for exploration.

```csv
sample_id,condition,replicate,batch
control_1,control,1,batch_1
control_2,control,2,batch_2
treated_1,treated,1,batch_1
treated_2,treated,2,batch_2
```

Sample identifiers must match the count-matrix column names exactly.

## Demonstration dataset

The repository includes a deterministic simulated dataset containing:

- 500 genes
- 3 control samples
- 3 treated samples
- 60 simulated differentially expressed genes

The dataset is intended only for demonstration and software testing. It is not presented as biological evidence.

Regenerate it with:

```bash
python scripts/generate_demo_data.py
```

## Project structure

```text
genomqc-studio/
├── app.py
├── demo_data/
├── reports/
├── scripts/
│   └── generate_demo_data.py
├── src/
│   ├── analysis.py
│   ├── data_loader.py
│   └── validator.py
├── tests/
├── environment.yml
├── pytest.ini
└── README.md
```

## Reproducibility and privacy

- The software environment is defined in `environment.yml`.
- Demo-data generation uses a fixed random seed.
- Core validation and analysis functions are covered by automated tests.
- Uploaded datasets are processed within the active application session.
- Uploaded files are not automatically written to the project directory.

## Development roadmap

- DESeq2 differential expression analysis
- Design-formula and confounding audit
- Volcano and MA plots
- Differential-expression result tables
- Expression heatmaps
- Exportable HTML analysis report
- Containerized deployment
- Continuous integration

## Limitations

GenomQC Studio is a research and educational tool. Its output must be interpreted alongside the experimental design and relevant biological context. It is not intended for clinical diagnosis.

## Author

**Muhammet Ozcelik**  
Computational genomics and reproducible bioinformatics workflows  
GitHub: [muhammetozcelik](https://github.com/muhammetozcelik)