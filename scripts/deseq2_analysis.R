# ============================================================================
# DESeq2 Differential Expression Analysis Pipeline
# ============================================================================

# Load required libraries
library(DESeq2)
library(tidyverse)
library(readxl)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files
count_file <- "count_data.csv"           #count matrix file
metadata_file <- "metadata.xlsx"          #sample metadata file

# Output file prefix (will be added to all output files)
output_prefix <- "analysis"              

# Analysis parameters
reference_condition <- "Healthy"          # Reference level for comparison
comparison_condition <- "Diseased"        # Condition to compare against reference
log2fc_threshold <- 1.5                   # Log2 fold change threshold
padj_threshold <- 0.05                    # Adjusted p-value threshold

# Output directory
output_dir <- "results"                   

# ============================================================================
# SETUP
# ============================================================================

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

# ============================================================================
# Step 1: Read and prepare count data
# ============================================================================

counts_data <- read.csv(count_file)

# Remove duplicates based on gene_symbol and set row names
counts_data <- counts_data %>% 
  distinct(gene_symbol, .keep_all = TRUE) %>%
  remove_rownames() %>% 
  column_to_rownames(var = "gene_symbol")

# Check for non-numeric columns
sapply(counts_data, class)

# ============================================================================
# Step 2: Read and prepare metadata
# ============================================================================

colData <- read_xlsx(metadata_file)

# Remove duplicates and set row names
colData <- colData %>% 
  distinct(Gene_id, .keep_all = TRUE) %>%
  remove_rownames() %>% 
  column_to_rownames(var = "Gene_id")

# Ensure column names match row names in column data
if (!all(colnames(counts_data) %in% rownames(colData))) {
  stop("Column names in count data do not match row names in column data.")
}

# ============================================================================
# Step 3: Construct DESeqDataSet object
# ============================================================================

# Ensure the row names in colData match the column names in counts_data
count_data <- counts_data
colData <- colData[colnames(count_data), , drop = FALSE]

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = round(count_data),
  colData = colData,
  design = ~ Condition
)

# ============================================================================
# Step 4: Pre-filtering and setting factor level
# ============================================================================

# Set the reference level for comparison
dds$Condition <- relevel(dds$Condition, ref = reference_condition)

# ============================================================================
# Step 5: Run DESeq analysis
# ============================================================================
dds <- DESeq(dds)
res <- results(dds)

# Summary of results
summary(res)

# Results with adjusted p-value threshold
res_filtered <- results(dds, alpha = padj_threshold)
summary(res_filtered)

# ============================================================================
# Step 6: Perform specific contrast
# ============================================================================

# List available contrasts
print(resultsNames(dds))

# Perform contrast
contrast_results <- results(
  dds, 
  contrast = c("Condition", comparison_condition, reference_condition)
)

# ============================================================================
# Step 7: Filter for upregulated and downregulated genes
# ============================================================================

# Upregulated genes
upregulated_genes <- subset(
  contrast_results, 
  log2FoldChange >= log2fc_threshold & padj <= padj_threshold
)

# Downregulated genes
downregulated_genes <- subset(
  contrast_results, 
  log2FoldChange <= -log2fc_threshold & padj <= padj_threshold
)

# All differentially expressed genes (up or down)
de_genes <- contrast_results[
  !is.na(contrast_results$padj) & 
  ((contrast_results$log2FoldChange >= log2fc_threshold | 
    contrast_results$log2FoldChange <= -log2fc_threshold) & 
   contrast_results$padj <= padj_threshold), 
]

# Print summary statistics
message("  Upregulated genes: ", nrow(upregulated_genes))
message("  Downregulated genes: ", nrow(downregulated_genes))
message("  Total DE genes: ", nrow(de_genes))

# ============================================================================
# Step 8: Save results to files
# ============================================================================

# Save all results
write.csv(
  as.data.frame(contrast_results),
  file = file.path(output_dir, paste0(output_prefix, "_all_results.csv")),
  row.names = TRUE
)

# Save upregulated genes
write.csv(
  as.data.frame(upregulated_genes),
  file = file.path(output_dir, paste0(output_prefix, "_upregulated_genes.csv")),
  row.names = TRUE
)

# Save downregulated genes
write.csv(
  as.data.frame(downregulated_genes),
  file = file.path(output_dir, paste0(output_prefix, "_downregulated_genes.csv")),
  row.names = TRUE
)

# Save all DE genes
write.csv(
  as.data.frame(de_genes),
  file = file.path(output_dir, paste0(output_prefix, "_differentially_expressed_genes.csv")),
  row.names = TRUE
)


# ============================================================================
# Save DESeq2 object for later use
# ============================================================================

saveRDS(dds, file = file.path(output_dir, paste0(output_prefix, "_dds_object.rds")))