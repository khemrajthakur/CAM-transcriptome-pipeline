# ============================================================================
# Transcriptomics Visualization Pipeline
# Volcano Plot and Heatmap Generation
# ============================================================================

# Load required libraries
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(gplots)
library(openxlsx)
library(grid)
library(DESeq2)  

set.seed(123)

# ============================================================================
# CONFIGURATION - Update these variables for your analysis
# ============================================================================

# Input files
count_file <- "count_data.csv"                    # Raw count data
deseq_results_file <- "analysis_all_results.csv"  # DESeq2 results from previous analysis

# Column patterns for identifying groups (for heatmap)
healthy_pattern <- "^H"      # Pattern to identify healthy/control samples (e.g., "^C", "^Control")
diseased_pattern <- "^R"     # Pattern to identify diseased/patient samples (e.g., "^P", "^Patient")

# Analysis parameters
padj_threshold <- 0.05       # Adjusted p-value threshold
log2fc_threshold <- 1.5      # Log2 fold change threshold
heatmap_diff_threshold <- 4  # Difference threshold for heatmap gene selection
n_genes_per_group <- 25      # Number of top genes per group for heatmap

# Output settings
output_prefix <- "analysis"   # Prefix for all output files
output_dir <- "results"       # Output directory

# Plot titles
volcano_title <- "Volcano Plot: Diseased vs Healthy"
heatmap_title <- "Global Z-score Heatmap\n(Diseased vs Healthy Reference)"

# ============================================================================
# SETUP
# ============================================================================

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

# ============================================================================
# PART 1: VOLCANO PLOT
# ============================================================================

# Step 1: Load DESeq2 results
res_df <- read.csv(deseq_results_file, row.names = 1)

# Remove rows with NA values
res_df <- na.omit(res_df)

# Add gene names
res_df$gene <- rownames(res_df)

# Step 2: Define significance categories
res_df$significance <- "Not Significant"
res_df$significance[res_df$log2FoldChange > log2fc_threshold & res_df$padj < padj_threshold] <- "Upregulated"
res_df$significance[res_df$log2FoldChange < -log2fc_threshold & res_df$padj < padj_threshold] <- "Downregulated"
res_df$significance <- factor(res_df$significance, 
                              levels = c("Upregulated", "Downregulated", "Not Significant"))

# Step 3: Create volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  
  # Color scheme
  scale_color_manual(values = c("Upregulated" = "#d62728",      # Red
                                "Downregulated" = "#1f77b4",    # Blue
                                "Not Significant" = "grey70")) +
  
  # Add threshold lines
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), 
             linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(padj_threshold), 
             linetype = "dashed", color = "black", alpha = 0.5) +
  
  # Labels and theme
  labs(title = volcano_title,
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Regulation") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank())

# Step 4: Add gene labels for top significant genes
top_genes <- res_df %>%
  filter(significance != "Not Significant") %>%
  arrange(padj) %>%
  head(20)

volcano_plot_labeled <- volcano_plot +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 20,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.2)

# Step 6: Save volcano plots

ggsave(
  file.path(output_dir, paste0(output_prefix, "_volcano_plot_labeled.png")),
  plot = volcano_plot_labeled,
  width = 10,
  height = 7,
  dpi = 300
)

# ============================================================================
# PART 2: HEATMAP
# ============================================================================

# Step 1: Load count data
df <- fread(count_file, check.names = FALSE)

# Ensure gene symbol column exists
possible_gene_cols <- c("gene_symbol", "Gene.symbol", "Gene", "gene", "SYMBOL", "symbol")
gene_col <- intersect(possible_gene_cols, colnames(df))
if (length(gene_col) == 0) {
  stop("No gene symbol column found. Please ensure a column named 'gene_symbol' is present.")
}
gene_col <- gene_col[1]
setnames(df, gene_col, "Gene_symbol")

# Step 2: Detect healthy and diseased columns
healthy_cols  <- grep(healthy_pattern, colnames(df), ignore.case = FALSE, value = TRUE)
diseased_cols <- grep(diseased_pattern, colnames(df), ignore.case = FALSE, value = TRUE)

if (length(healthy_cols) == 0) {
  stop("No healthy/control columns found with pattern: ", healthy_pattern)
}
if (length(diseased_cols) == 0) {
  stop("No diseased/patient columns found with pattern: ", diseased_pattern)
}

message("Found ", length(healthy_cols), " healthy samples and ", length(diseased_cols), " diseased samples")

# Step 3: Compute means and log2 transform
df[, Healthy_Mean  := rowMeans(.SD, na.rm = TRUE), .SDcols = healthy_cols]
df[, Diseased_Mean := rowMeans(.SD, na.rm = TRUE), .SDcols = diseased_cols]

# Save mean expressions
fwrite(
  df[, .(Gene_symbol, Healthy_Mean, Diseased_Mean)],
  file.path(output_dir, paste0(output_prefix, "_mean_expression.csv"))
)

# Log2 transform means
df <- df %>%
  mutate(
    Healthy  = log2(Healthy_Mean  + 1),
    Diseased = log2(Diseased_Mean + 1)
  )

fwrite(
  df[, .(Gene_symbol, Healthy, Diseased)],
  file.path(output_dir, paste0(output_prefix, "_log2_mean_expression.csv"))
)

# Step 4: Calculate fold change (Diseased - Healthy)
df <- df %>% mutate(Diff_DH = Diseased - Healthy)

# Step 5: Select top up/downregulated genes
upregulated_in_diseased <- df %>%
  filter(is.finite(Diff_DH) & Diff_DH >= heatmap_diff_threshold) %>%
  arrange(desc(Diff_DH))

downregulated_in_diseased <- df %>%
  filter(is.finite(Diff_DH) & Diff_DH <= -heatmap_diff_threshold) %>%
  arrange(Diff_DH)

# Fallbacks if not enough genes pass threshold
if (nrow(upregulated_in_diseased) < n_genes_per_group) {
  message("Not enough genes above threshold, selecting top ", n_genes_per_group, " by fold change")
  upregulated_in_diseased <- df %>% slice_max(order_by = Diff_DH, n = n_genes_per_group)
}
if (nrow(downregulated_in_diseased) < n_genes_per_group) {
  message("Not enough genes below threshold, selecting top ", n_genes_per_group, " by fold change")
  downregulated_in_diseased <- df %>% slice_min(order_by = Diff_DH, n = n_genes_per_group)
}

upregulated_top <- upregulated_in_diseased %>% 
  slice_head(n = min(n_genes_per_group, nrow(upregulated_in_diseased)))

downregulated_top <- downregulated_in_diseased %>% 
  slice_head(n = min(n_genes_per_group, nrow(downregulated_in_diseased)))

# Combine and shuffle
selected_genes <- bind_rows(upregulated_top, downregulated_top)
shuffled_genes <- selected_genes %>% slice_sample(n = nrow(selected_genes))

# Save selected genes
write.xlsx(
  as.data.frame(shuffled_genes),
  file.path(output_dir, paste0(output_prefix, "_heatmap_selected_genes.xlsx")),
  rowNames = FALSE
)

# Clean gene symbols
shuffled_genes <- shuffled_genes %>%
  filter(!is.na(Gene_symbol) & Gene_symbol != "" & Gene_symbol != "NA") %>%
  mutate(Gene_symbol = make.unique(as.character(Gene_symbol), sep = "_"))

# Step 6: Prepare matrix for heatmap
heatmat <- shuffled_genes %>%
  dplyr::select(Gene_symbol, Healthy, Diseased) %>%
  column_to_rownames("Gene_symbol") %>%
  as.matrix()

# Replace non-finite values
if (any(!is.finite(heatmat))) {
  finite_min <- min(heatmat[is.finite(heatmat)], na.rm = TRUE)
  heatmat[!is.finite(heatmat)] <- finite_min - 1
}

# Step 7: Global Z-score normalization
z_vals <- scale(as.vector(heatmat))
z_heatmat <- matrix(
  z_vals,
  nrow = nrow(heatmat),
  ncol = ncol(heatmat),
  byrow = FALSE
)

rownames(z_heatmat) <- rownames(heatmat)
colnames(z_heatmat) <- colnames(heatmat)

# Step 8: Define color palette (green = low, black = mid, red = high)
my_colors <- colorRampPalette(c(
  "green",  # low expression
  "black",  # mid
  "red"     # high expression (upregulated in diseased)
))(200)

# Step 9: Generate heatmap with pheatmap
message("Generating pheatmap...")
p <- pheatmap(
  z_heatmat,
  color         = my_colors,
  cluster_rows  = TRUE,
  cluster_cols  = TRUE,
  show_rownames = TRUE,
  main          = heatmap_title,
  legend        = TRUE,
  fontsize_row  = 5,
  fontsize_col  = 10,
  border_color  = NA,
  filename      = file.path(output_dir, paste0(output_prefix, "_heatmap_pheatmap.png")),
  width         = 8,
  height        = 10
)

# Step 10: Generate heatmap with histogram using heatmap.2
png(
  file.path(output_dir, paste0(output_prefix, "_heatmap_with_histogram.png")),
  width = 1200,
  height = 1000,
  res = 120
)

heatmap.2(
  z_heatmat,
  col          = my_colors,
  scale        = "none",
  trace        = "none",
  density.info = "histogram",
  key          = TRUE,
  key.title    = "Color Key & Histogram",
  key.xlab     = "Global Z-score",
  dendrogram   = "both",
  Rowv         = TRUE,
  Colv         = TRUE,
  margins      = c(6, 8),
  lhei         = c(1, 4),
  cexCol       = 1.2,
  srtCol       = 0,
  adjCol       = c(0.5, 1)
)

dev.off()