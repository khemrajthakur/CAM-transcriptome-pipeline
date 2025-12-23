# ============================================================
# LOAD LIBRARIES
# ============================================================
library(DESeq2)
library(gprofiler2)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(openxlsx)

# ============================================================
# LOAD DESEQ2 RESULTS
# ============================================================
deseq_results <- read.csv("upregulated_genes_neutro.csv")

# ============================================================
# FILTER SIGNIFICANT GENES
# ============================================================
significant_genes <- deseq_results %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) %>%
  pull(gene_symbol)

# ============================================================
# MAP GENE SYMBOLS TO ENTREZ IDs
# ============================================================
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = significant_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

entrez_ids <- na.omit(entrez_ids)

# ============================================================
# GO ENRICHMENT ANALYSIS (g:Profiler)
# ============================================================
go_results <- gost(
  query = names(entrez_ids),
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "GO:CC")
)

# ============================================================
# PROCESS GO RESULTS
# ============================================================
go_results_df <- go_results$result %>%
  mutate(across(where(is.list), ~ sapply(.x, toString))) %>%
  mutate(fold_enrichment = intersection_size / effective_domain_size)

# Save GO results
write.csv(go_results_df, "GO_results_up_neutro.csv", row.names = FALSE)

# ============================================================
# LOAD GO RESULTS FOR PLOTTING
# ============================================================
new_go_data <- read.csv("GO_results_up_neutro.csv")
new_go_data <- as.data.frame(new_go_data)

# ============================================================
# SELECT TOP GO TERMS
# ============================================================
top_go <- new_go_data %>%
  arrange(desc(fold_enrichment)) %>%
  head(20)

# ============================================================
# CAPITALIZE GO TERM NAMES (KEY FIX)
# ============================================================
top_go <- top_go %>%
  mutate(term_name = tools::toTitleCase(tolower(term_name)))

# ============================================================
# BUBBLE PLOT
# ============================================================
bubble_plot <- ggplot(
  top_go,
  aes(
    x = term_name,
    y = fold_enrichment,
    size = intersection_size,
    fill = source,
    color = -log10(p_value)
  )
) +
  geom_point(alpha = 1.0, shape = 21, stroke = 1.5) +
  coord_flip() +
  theme_minimal() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "GO Term",
    y = "Fold Enrichment",
    size = "Gene Count",
    fill = "Process Type",
    color = "-log10(p-value)"
  ) +
  theme(
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5), order = 1),
    size = guide_legend(order = 2),
    color = guide_colorbar(order = 3)
  )

# ============================================================
# DISPLAY PLOT
# ============================================================
print(bubble_plot)

# ============================================================
# SAVE HIGH-RESOLUTION FIGURE
# ============================================================
ggsave(
  filename = "GO_bubble_plot_up_neutro.tiff",
  plot = bubble_plot,
  width = 12,
  height = 9,
  dpi = 600,
  bg = "white"
)
