#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
})


fail <- function(message_text) {
    message("ERROR: ", message_text)
    quit(status = 1)
}


arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) != 6) {
    fail(
        paste(
            "Usage: Rscript run_deseq2.R",
            "<counts.csv> <metadata.csv>",
            "<condition_column> <reference>",
            "<comparison> <output_directory>"
        )
    )
}

counts_path <- arguments[[1]]
metadata_path <- arguments[[2]]
condition_column <- arguments[[3]]
reference_condition <- arguments[[4]]
comparison_condition <- arguments[[5]]
output_directory <- arguments[[6]]

dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
)

counts_table <- read.csv(
    counts_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

metadata <- read.csv(
    metadata_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (!"gene_id" %in% colnames(counts_table)) {
    fail("Count matrix must contain a gene_id column.")
}

if (!"sample_id" %in% colnames(metadata)) {
    fail("Metadata must contain a sample_id column.")
}

if (!condition_column %in% colnames(metadata)) {
    fail(
        paste(
            "Metadata does not contain condition column:",
            condition_column
        )
    )
}

if (anyDuplicated(counts_table$gene_id)) {
    fail("Gene identifiers must be unique.")
}

if (anyDuplicated(metadata$sample_id)) {
    fail("Sample identifiers must be unique.")
}

gene_ids <- counts_table$gene_id
counts_table$gene_id <- NULL

count_matrix <- as.matrix(counts_table)
rownames(count_matrix) <- gene_ids

if (!is.numeric(count_matrix)) {
    fail("All count values must be numeric.")
}

if (any(is.na(count_matrix))) {
    fail("Count matrix contains missing values.")
}

if (any(count_matrix < 0)) {
    fail("Count matrix contains negative values.")
}

if (any(count_matrix != floor(count_matrix))) {
    fail("DESeq2 requires whole-number counts.")
}

metadata_sample_ids <- metadata$sample_id

if (!setequal(colnames(count_matrix), metadata_sample_ids)) {
    fail(
        paste(
            "Count matrix and metadata sample identifiers",
            "do not match."
        )
    )
}

rownames(metadata) <- metadata$sample_id
metadata$sample_id <- NULL
metadata <- metadata[colnames(count_matrix), , drop = FALSE]

condition_values <- as.character(
    metadata[[condition_column]]
)

if (!reference_condition %in% condition_values) {
    fail(
        paste(
            "Reference condition was not found:",
            reference_condition
        )
    )
}

if (!comparison_condition %in% condition_values) {
    fail(
        paste(
            "Comparison condition was not found:",
            comparison_condition
        )
    )
}

if (reference_condition == comparison_condition) {
    fail(
        "Reference and comparison conditions must differ."
    )
}

metadata[[condition_column]] <- relevel(
    factor(metadata[[condition_column]]),
    ref = reference_condition
)

keep_genes <- rowSums(count_matrix) >= 10

if (sum(keep_genes) < 2) {
    fail(
        "Fewer than two genes passed the count filter."
    )
}

filtered_counts <- count_matrix[
    keep_genes,
    ,
    drop = FALSE
]

design_formula <- reformulate(condition_column)

dataset <- DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = metadata,
    design = design_formula
)

dataset <- DESeq(
    dataset,
    quiet = TRUE
)

result <- results(
    dataset,
    contrast = c(
        condition_column,
        comparison_condition,
        reference_condition
    ),
    alpha = 0.05
)

result_table <- as.data.frame(result)
result_table$gene_id <- rownames(result_table)

result_table <- result_table[
    ,
    c(
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj"
    )
]

result_table <- result_table[
    order(result_table$padj, na.last = TRUE),
    ,
    drop = FALSE
]

write.csv(
    result_table,
    file.path(
        output_directory,
        "deseq2_results.csv"
    ),
    row.names = FALSE,
    na = ""
)

normalized_matrix <- counts(
    dataset,
    normalized = TRUE
)

normalized_table <- data.frame(
    gene_id = rownames(normalized_matrix),
    normalized_matrix,
    check.names = FALSE
)

write.csv(
    normalized_table,
    file.path(
        output_directory,
        "normalized_counts.csv"
    ),
    row.names = FALSE
)

transformed_dataset <- varianceStabilizingTransformation(
    dataset,
    blind = FALSE
)

vst_matrix <- assay(transformed_dataset)

vst_table <- data.frame(
    gene_id = rownames(vst_matrix),
    vst_matrix,
    check.names = FALSE
)

write.csv(
    vst_table,
    file.path(
        output_directory,
        "vst_counts.csv"
    ),
    row.names = FALSE
)

significant <- (
    !is.na(result_table$padj)
    & result_table$padj < 0.05
)

significant_up <- (
    significant
    & result_table$log2FoldChange >= 1
)

significant_down <- (
    significant
    & result_table$log2FoldChange <= -1
)

run_summary <- data.frame(
    metric = c(
        "reference_condition",
        "comparison_condition",
        "input_genes",
        "analyzed_genes",
        "samples",
        "significant_genes",
        "significant_up",
        "significant_down"
    ),
    value = c(
        reference_condition,
        comparison_condition,
        nrow(count_matrix),
        nrow(filtered_counts),
        ncol(filtered_counts),
        sum(significant),
        sum(significant_up),
        sum(significant_down)
    )
)

write.table(
    run_summary,
    file.path(
        output_directory,
        "run_summary.tsv"
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

message(
    "DESeq2 analysis completed successfully: ",
    comparison_condition,
    " vs ",
    reference_condition
)