from pathlib import Path

import numpy as np
import pandas as pd


RANDOM_SEED = 42
GENE_COUNT = 500
REPLICATES_PER_CONDITION = 3

OUTPUT_DIRECTORY = Path("demo_data")
OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(RANDOM_SEED)

gene_ids = [
    f"GENE_{number:04d}"
    for number in range(1, GENE_COUNT + 1)
]

sample_ids = [
    "control_1",
    "control_2",
    "control_3",
    "treated_1",
    "treated_2",
    "treated_3",
]

conditions = [
    "control",
    "control",
    "control",
    "treated",
    "treated",
    "treated",
]

base_means = 10 ** rng.uniform(0.5, 3.2, GENE_COUNT)

fold_changes = np.ones(GENE_COUNT)
fold_changes[:30] = 4.0
fold_changes[30:60] = 0.25

count_data = {"gene_id": gene_ids}

for sample_id, condition in zip(sample_ids, conditions):
    library_scale = rng.lognormal(mean=0, sigma=0.08)
    expected_counts = base_means * library_scale

    if condition == "treated":
        expected_counts = expected_counts * fold_changes

    dispersion_size = 20
    probability = dispersion_size / (
        dispersion_size + expected_counts
    )

    count_data[sample_id] = rng.negative_binomial(
        dispersion_size,
        probability,
    )

counts = pd.DataFrame(count_data)

metadata = pd.DataFrame(
    {
        "sample_id": sample_ids,
        "condition": conditions,
        "replicate": [1, 2, 3, 1, 2, 3],
        "batch": [
            "batch_1",
            "batch_2",
            "batch_1",
            "batch_1",
            "batch_2",
            "batch_1",
        ],
    }
)

truth = pd.DataFrame(
    {
        "gene_id": gene_ids,
        "expected_log2_fold_change": np.log2(fold_changes),
        "is_differential": fold_changes != 1,
    }
)

counts.to_csv(OUTPUT_DIRECTORY / "counts.csv", index=False)
metadata.to_csv(OUTPUT_DIRECTORY / "metadata.csv", index=False)
truth.to_csv(
    OUTPUT_DIRECTORY / "simulation_truth.csv",
    index=False,
)

print("Demo dataset generated successfully.")
print(f"Genes: {GENE_COUNT}")
print(f"Samples: {len(sample_ids)}")
print("Differential genes: 60")