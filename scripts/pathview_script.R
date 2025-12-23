# =============================================================================
# Pathview Pathway Visualization for DEGs
# CVID – Monocyte / Neutrophil Immune Pathways
# =============================================================================
# REQUIREMENTS:
# 1. A DESeq2 results object named `res` must exist in the environment
# 2. Row names of `res` must be gene SYMBOLS
# =============================================================================

# ==============================
# USER-DEFINED PARAMETERS
# ==============================
padj_cutoff   <- 0.05
log2fc_cutoff <- 1.5
species_id    <- "hsa"
output_suffix <- "cvid"

# ==============================
# LOAD REQUIRED PACKAGES
# ==============================
required_pkgs <- c("pathview", "org.Hs.eg.db", "png", "grid", "gridExtra")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
}

library(pathview)
library(org.Hs.eg.db)
library(png)
library(grid)
library(gridExtra)

# ==============================
# CHECK INPUT OBJECT
# ==============================
stopifnot(exists("res"))

# ==============================
# STEP 1: Prepare DEG data
# ==============================
deg_data <- as.data.frame(res)
deg_data$gene_symbol <- rownames(deg_data)

deg_data <- deg_data[
  !is.na(deg_data$log2FoldChange) &
  !is.na(deg_data$padj),
]

cat("Total genes:", nrow(deg_data), "\n")
cat("Significant genes (padj <", padj_cutoff, "):",
    sum(deg_data$padj < padj_cutoff), "\n")

# ==============================
# STEP 2: Filter DEGs
# ==============================
filtered_deg <- deg_data[
  deg_data$padj < padj_cutoff &
  abs(deg_data$log2FoldChange) > log2fc_cutoff,
]

cat("Genes passing filters:", nrow(filtered_deg), "\n")
cat("Upregulated:", sum(filtered_deg$log2FoldChange > log2fc_cutoff), "\n")
cat("Downregulated:", sum(filtered_deg$log2FoldChange < -log2fc_cutoff), "\n")

if (nrow(filtered_deg) == 0) {
  stop("No genes meet filtering criteria.")
}

# ==============================
# STEP 3: SYMBOL → ENTREZ
# ==============================
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = filtered_deg$gene_symbol,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

gene_data <- filtered_deg$log2FoldChange
names(gene_data) <- entrez_ids
gene_data <- gene_data[!is.na(names(gene_data))]

cat("Mapped genes:", length(gene_data), "\n")
cat("Unmapped genes:", sum(is.na(entrez_ids)), "\n")

# ==============================
# STEP 4: Immune Pathways
# ==============================
immune_pathways <- list(

  "hsa04630" = "JAK-STAT signaling pathway",
  "hsa04620" = "Toll-like receptor signaling pathway",
  "hsa04621" = "NOD-like receptor signaling pathway",
  "hsa04622" = "RIG-I-like receptor signaling pathway",
  "hsa04623" = "Cytosolic DNA-sensing pathway",
  "hsa04625" = "C-type lectin receptor signaling pathway",
  "hsa04210" = "Apoptosis pathway",

  "hsa04062" = "Chemokine signaling pathway",
  "hsa04670" = "Leukocyte transendothelial migration",
  "hsa04666" = "Fc gamma R-mediated phagocytosis",
  "hsa04145" = "Phagosome",
  "hsa05140" = "Leishmaniasis",
  "hsa04380" = "Osteoclast differentiation",

  "hsa04060" = "Cytokine-cytokine receptor interaction",
  "hsa04064" = "NF-kappa B signaling pathway",
  "hsa04668" = "TNF signaling pathway",
  "hsa04657" = "IL-17 signaling pathway",
  "hsa04010" = "MAPK signaling pathway",
  "hsa04151" = "PI3K-Akt signaling pathway",

  "hsa05208" = "Chemical carcinogenesis - ROS",
  "hsa04146" = "Peroxisome",
  "hsa04142" = "Lysosome"
)

# ==============================
# STEP 5: Pathview Visualization
# ==============================
output_dir <- "pathview_immune_pathways"
dir.create(output_dir, showWarnings = FALSE)

successful <- character()
skipped    <- character()

for (pid in names(immune_pathways)) {

  cat("\nProcessing:", immune_pathways[[pid]], "\n")

  tryCatch({
    pathview(
      gene.data   = gene_data,
      pathway.id  = pid,
      species     = species_id,
      gene.idtype = "ENTREZ",
      out.suffix  = output_suffix,
      kegg.native = TRUE,
      same.layer  = TRUE,
      limit       = list(gene = c(-3, 3)),
      out.dir     = output_dir
    )
    successful <- c(successful, immune_pathways[[pid]])
    cat("  ✓ Success\n")
  }, error = function(e) {
    skipped <- c(skipped, immune_pathways[[pid]])
    cat("  ✗ Skipped:", e$message, "\n")
  })
}

# ==============================
# STEP 6: Save DEG Table
# ==============================
write.csv(
  filtered_deg,
  file.path(output_dir, "filtered_DEGs_log2FC_gt_1.5.csv"),
  row.names = FALSE
)

# ==============================
# STEP 7: Combine PNGs → PDF
# ==============================
png_files <- list.files(
  output_dir,
  pattern = paste0(output_suffix, ".*\\.png$"),
  full.names = TRUE
)

png_files <- sort(png_files)

if (length(png_files) == 0) {
  stop("No PNG files found.")
}

output_pdf <- file.path(output_dir, "All_Immune_Pathways_CVID.pdf")

pdf(output_pdf, width = 11, height = 8.5)

for (i in seq_along(png_files)) {
  img <- readPNG(png_files[i])
  grid.newpage()

  pathway_name <- gsub(
    paste0("\\.", output_suffix, "\\.png"),
    "",
    basename(png_files[i])
  )

  grid.text(
    pathway_name,
    x = 0.5, y = 0.97,
    gp = gpar(fontsize = 14, fontface = "bold")
  )

  grid.raster(img, x = 0.5, y = 0.48, width = 0.95, height = 0.9)
}

dev.off()

cat("\nPDF created:", output_pdf, "\n")
cat("========================================\n")
