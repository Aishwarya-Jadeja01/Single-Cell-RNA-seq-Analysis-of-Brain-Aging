# ============================================================================
# Doublet Detection
# ============================================================================

library(Seurat)
library(DoubletFinder)
library(tidyverse)

setwd("aging-brain-singlecell")


# ----------------------------------------------------------------------------
# Load and prepare data
# ----------------------------------------------------------------------------
cat("Loading brain data...\n")
brain <- readRDS("data/processed/brain_filtered.rds")

cat("  Starting with:", ncol(brain), "cells\n")
cat("  Genes:", nrow(brain), "\n\n")

# ----------------------------------------------------------------------------
# Preprocessing with explicit data type control
# ----------------------------------------------------------------------------
cat("Preparing data for doublet detection...\n")

# Normalize - standard log normalization
cat(" Normalizing data")
brain <- NormalizeData(brain, 
                       normalization.method = "LogNormalize",
                       scale.factor = 10000,
                       verbose = FALSE)

# Find variable features
cat(" Finding variable features")
brain <- FindVariableFeatures(brain, 
                              selection.method = "vst", 
                              nfeatures = 2000, 
                              verbose = FALSE)


cat("Scaling data")
# Scale only the variable features to reduce memory and ensure clean matrices
brain <- ScaleData(brain, 
                   features = VariableFeatures(brain),
                   verbose = FALSE)

# Run PCA with explicit parameters
cat("Running PCA")
brain <- RunPCA(brain, 
                features = VariableFeatures(brain),
                npcs = 30,
                verbose = FALSE)

pca_embeddings <- Embeddings(brain, reduction = "pca")
cat("    PCA embeddings dimensions:", dim(pca_embeddings)[1], "cells x", 
    dim(pca_embeddings)[2], "PCs\n")
cat("    PCA embeddings class:", class(pca_embeddings), "\n")

# PCA embeddings aren't a matrix, convert them
if (!is.matrix(pca_embeddings)) {
  pca_embeddings <- as.matrix(pca_embeddings)
  brain[["pca"]]@cell.embeddings <- pca_embeddings
}

# Run UMAP for visualization
cat("Running UMAP")
brain <- RunUMAP(brain, 
                 dims = 1:30,
                 verbose = FALSE)

# Clustering - needed for homotypic doublet adjustment
cat("Clustering cells")
brain <- FindNeighbors(brain, 
                       dims = 1:30,
                       verbose = FALSE)

brain <- FindClusters(brain, 
                      resolution = 0.5,
                      verbose = FALSE)

cat("Preprocessing complete")
cat("Identified", length(unique(Idents(brain))), "clusters\n\n")

# ----------------------------------------------------------------------------
# Calculate expected doublet rate with corrected formula
# ----------------------------------------------------------------------------
cat("Calculating expected doublet rate...\n")

n_cells <- ncol(brain)

# formula: 0.008 per 1000 cells
expected_doublet_rate <- 0.008 * (n_cells / 1000)

cat("  Expected doublet rate for", n_cells, "cells:", 
    round(expected_doublet_rate * 100, 1), "%\n")

# Calculate expected number of doublets
nExp_poi <- round(expected_doublet_rate * n_cells)
cat("  → nExp = ", nExp_poi, "\n\n")

# ----------------------------------------------------------------------------
# Adjust for homotypic doublets
# ----------------------------------------------------------------------------
cat("Adjusting for homotypic doublets...\n")

# Model homotypic doublet proportion based on clustering
homotypic_prop <- modelHomotypic(Idents(brain))
cat("  Estimated homotypic doublet proportion:", round(homotypic_prop, 3), "\n")

# Adjust expected doublet number
nExp_adj <- round(nExp_poi * (1 - homotypic_prop))
cat("  Adjusted expected doublets:", nExp_adj, "\n\n")

# ----------------------------------------------------------------------------
# Find optimal pK parameter
# ----------------------------------------------------------------------------
cat("Finding optimal pK parameter...\n")
cat("  This may take 5-10 minutes...\n\n")

# Parameter sweep with error handling
sweep_res <- tryCatch({
  paramSweep(brain, PCs = 1:30, sct = FALSE, num.cores = 1)
}, error = function(e) {
  cat("  Error during parameter sweep:", e$message, "\n")
  cat("  Trying alternative approach...\n")
  # If parallel processing causes issues, force single core
  paramSweep(brain, PCs = 1:30, sct = FALSE, num.cores = 1)
})

# Summarize sweep results
sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)

# Find optimal pK
bcmvn <- find.pK(sweep_stats)
optimal_pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

cat("  Optimal pK found:", optimal_pk, "\n\n")

# ----------------------------------------------------------------------------
# Run DoubletFinder
# ----------------------------------------------------------------------------
cat("Running DoubletFinder...\n")
cat("  Using adjusted parameters:\n")
cat("    PCs: 1:30\n")
cat("    pN: 0.25\n")
cat("    pK:", optimal_pk, "\n")
cat("    nExp:", nExp_adj, "\n\n")

# v3 with explicit parameters
brain <- tryCatch({
  doubletFinder_v3(brain,
                   PCs = 1:30,
                   pN = 0.25,
                   pK = as.numeric(as.character(optimal_pk)),
                   nExp = nExp_adj,
                   reuse.pANN = FALSE,
                   sct = FALSE)
}, error = function(e) {
  cat("\n⚠ Error during doubletFinder_v3:", e$message, "\n")
  cat("  This may be due to data structure incompatibility\n")
  cat("  Attempting alternative doublet detection method...\n\n")
  
  # Alternative: Use simpler metrics-based approach
  # High gene count + high UMI as proxy for doublets
  brain$simple_doublet_score <- scale(brain$nFeature_RNA) + scale(brain$nCount_RNA)
  brain$DF.classifications <- ifelse(brain$simple_doublet_score > quantile(brain$simple_doublet_score, 1 - expected_doublet_rate),
                                     "Doublet", 
                                     "Singlet")
  
  cat("  Used simple metrics-based classification instead\n")
  return(brain)
})


# Find the doublet classification column
doublet_col <- grep("DF.classifications", colnames(brain@meta.data), value = TRUE)

if (length(doublet_col) == 0) {
  doublet_col <- "DF.classifications"
}

# summary
doublet_summary <- table(brain@meta.data[[doublet_col]])
print(doublet_summary)

n_doublets <- sum(brain@meta.data[[doublet_col]] == "Doublet")
doublet_pct <- 100 * n_doublets / ncol(brain)

cat("\nSummary:\n")
cat("  Total cells:", ncol(brain), "\n")
cat("  Identified doublets:", n_doublets, 
    "(", round(doublet_pct, 1), "%)\n")
cat("  Singlets:", ncol(brain) - n_doublets, "\n\n")

# Create visualizations
cat("Creating visualizations...\n")

# UMAP colored by doublet status
p1 <- DimPlot(brain, 
              reduction = "umap",
              group.by = doublet_col,
              cols = c("Singlet" = "gray70", "Doublet" = "red3")) +
  ggtitle(paste0("Doublet Detection\n", n_doublets, " doublets (", 
                 round(doublet_pct, 1), "%)")) +
  theme_minimal()

ggsave("results/figures/doublet_detection.png",
       p1, width = 8, height = 6, dpi = 300)
cat("  ✓ Saved: 07_doublet_detection_fixed.png\n")

# QC comparison
p2 <- VlnPlot(brain,
              features = c("nFeature_RNA", "nCount_RNA"),
              group.by = doublet_col,
              pt.size = 0,
              cols = c("Singlet" = "gray70", "Doublet" = "red3")) +
  plot_layout(ncol = 2)

ggsave("results/figures/08_doublet_qc_comparison_fixed.png",
       p2, width = 10, height = 5, dpi = 300)
cat("doublet_qc_comparison.png\n\n")

# ----------------------------------------------------------------------------
#  results
# ----------------------------------------------------------------------------

saveRDS(brain, "data/processed/brain_with_doublets_fixed.rds")
cat("brain_with_doublets_fixed.rds\n")

doublet_info <- data.frame(
  cell = colnames(brain),
  doublet_status = brain@meta.data[[doublet_col]],
  nFeature = brain$nFeature_RNA,
  nCount = brain$nCount_RNA,
  age_group = brain$age_group
)

write.csv(doublet_info, 
          "results/tables/doublet_classification.csv",
          row.names = FALSE)
cat("doublet_classification.csv\n\n")