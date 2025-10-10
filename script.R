library(GEOquery)
library(dplyr)
library(hgu133plus2.db)
library(tibble)
library(msigdbr)
library(ggplot2)
library(tidyr)
library(GENIE3)

# download data
gse <- tryCatch(
    getGEO("GSE45827", GSEMatrix = TRUE),
    error = function(e) {
        stop("Failed to download GEO dataset")
    }
)

# extract expression data and metadata
if (is.list(gse)) gse <- gse[[1]]
expression_data <- exprs(gse)
metadata <- pData(gse)

# map probe IDs to gene symbols
probe_ids <- rownames(expression_data)
gene_symbols <- mapIds(
    hgu133plus2.db,
    keys = probe_ids,
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
)

# remove NA gene symbols
valid_idx <- !is.na(gene_symbols)
expression_data <- expression_data[valid_idx, ]
gene_symbols <- gene_symbols[valid_idx]

# summarize duplicate genes by taking the mean expression
expression_data <- aggregate(expression_data, by = list(Gene = gene_symbols), FUN = mean)
rownames(expression_data) <- expression_data$Gene
expression_data$Gene <- NULL

# histogram of the expression data
expression_data_long <- gather(expression_data, key = "Sample", value = "Expression")
ggplot(expression_data_long, aes(x = Expression)) +
    geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Histogram of Gene Expression", x = "Expression Level", y = "Frequency") +
    theme_minimal()

# get Hallmark gene sets for ESTROGEN_RESPONSE_EARLY and ESTROGEN_RESPONSE_LATE
hallmark_genes <- msigdbr(species = "Homo sapiens", collection = "H")
estrogen_genes <- hallmark_genes %>%
    filter(gs_name %in% c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE")) %>%
    pull(gene_symbol)

# get unique estrogen gene symbols
estrogen_genes <- unique(estrogen_genes)

# sanity check: are all estrogen_genes are in the expression_data
missing_genes <- setdiff(estrogen_genes, rownames(expression_data))
if (length(missing_genes) > 0) {
    message("The following estrogen genes are missing from the expression data: ", paste(missing_genes, collapse = ", "))
} else {
    message("All estrogen genes are present in the expression data.")
}

# subset expression_data to include only estrogen genes
estrogen_expression <- expression_data[rownames(expression_data) %in% estrogen_genes, ]

# histogram of the expression data
estrogen_expression_long <- gather(estrogen_expression, key = "Sample", value = "Expression")
ggplot(estrogen_expression_long, aes(x = Expression)) +
    geom_histogram(binwidth = 0.5, fill = "red", color = "black", alpha = 0.7) +
    labs(title = "Histogram of Estrogen Gene Expression", x = "Expression Level", y = "Frequency") +
    theme_minimal()
