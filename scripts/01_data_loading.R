library(Seurat)
library(readxl)
library(tidyverse)

setwd("aging-brain-singlecell")

# ======================================================================
# Load Excel files
# ======================================================================

# Load counts (genes as rows, cells as columns)
counts <- readxl::read_excel("data/raw/Brain_FACS_counts.xlsx")
counts <- as.data.frame(counts)

# Expect first column = gene names
rownames(counts) <- counts[,1]
counts <- counts[,-1]

# Convert to numeric
counts <- apply(counts, 2, as.numeric)
rownames(counts) <- rownames(readxl::read_excel("data/raw/Brain_FACS_counts.xlsx", col_names = TRUE)[,1])

cat("Dimensions of counts matrix:", nrow(counts), "genes x", ncol(counts), "cells\n\n")

# Load metadata
cat("Loading metadata...\n")
meta <- readxl::read_excel("data/raw/Brain_FACS_counts.xlsx")
meta <- as.data.frame(meta)

# metadata row names match cell IDs in counts
rownames(meta) <- meta$cell_id  
meta <- meta[colnames(counts), ]

cat("Metadata columns:\n")
print(colnames(meta))
cat("\n")

# ======================================================================
#  Create Seurat object
# ======================================================================

brain <- CreateSeuratObject(counts = counts, meta.data = meta, project = "AgingBrain")

cat("Seurat object created:\n")
print(brain)
cat("\n")

# ======================================================================
# STEP 3: Save Seurat object
# ======================================================================
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(brain, "data/processed/brain_seurat_full.rds")

cat("Saved Seurat object to data/processed/brain_seurat_full.rds\n\n")
