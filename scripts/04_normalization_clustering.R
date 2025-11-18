# ============================================================================
# Normalization, Dimensionality Reduction, and Clustering
# ============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)

setwd("aging-brain-singlecell")


# ----------------------------------------------------------------------------
# Load Quality-Controlled Data
# ----------------------------------------------------------------------------
cat("Loading filtered brain data...\n")

# Load the data that passed quality control
brain <- readRDS("data/processed/brain_final_qc.rds")

cat("  Loaded:", ncol(brain), "cells\n")
cat("  Genes:", nrow(brain), "genes\n")
cat("  Age groups:", unique(brain$age_group), "\n\n")

# ----------------------------------------------------------------------------
# Normalization
# ----------------------------------------------------------------------------
# LogNormalize method:
#  Log transform: log1p(x) = log(1 + x)

cat("Normalizing data...\n")
cat("  Method: LogNormalize\n")
cat("  Scale factor: 10,000\n")
cat("  Rationale: Account for different sequencing depths between cells\n\n")

brain <- NormalizeData(brain, 
                       normalization.method = "LogNormalize",
                       scale.factor = 10000,
                       verbose = FALSE)

cat("Normalization complete\n")
cat("Data is now in log-normalized space\n\n")

# ----------------------------------------------------------------------------
# Identify Highly Variable Features
# ----------------------------------------------------------------------------
cat("highly variable genes...\n")

brain <- FindVariableFeatures(brain,
                              selection.method = "vst",
                              nfeatures = 2000,
                              verbose = FALSE)

# Get the top variable genes to see what they are
top10_variable <- head(VariableFeatures(brain), 10)
cat("  Top most variable genes:\n")
for(i in 1:length(top10_variable)) {
  cat("    ", i, ". ", top10_variable[i], "\n", sep = "")
}
cat("\n")

# Visualize variable features
p_variable <- VariableFeaturePlot(brain)
p_variable <- LabelPoints(plot = p_variable, 
                          points = top10_variable,
                          repel = TRUE,
                          max.overlaps = 10)

ggsave("results/figures/09_variable_features.png",
       p_variable, width = 10, height = 6, dpi = 300)

cat("variable_features.png\n")

# ----------------------------------------------------------------------------
# Scaling
# ----------------------------------------------------------------------------

# Scale all genes (this can take a few minutes with 20k+ cells)
cat("Scaling all genes")
all_genes <- rownames(brain)

brain <- ScaleData(brain,
                   features = all_genes,
                   vars.to.regress = "nCount_RNA",
                   verbose = FALSE)

# ----------------------------------------------------------------------------
# Principal Component Analysis (PCA)
# ----------------------------------------------------------------------------

brain <- RunPCA(brain,
                features = VariableFeatures(brain),
                npcs = 50,
                verbose = FALSE)

# Look at genes that define each PC
cat("Top genes contributing to first 3 PCs:\n\n")
print(brain[["pca"]], dims = 1:3, nfeatures = 5)
cat("\n")

# Visualize PCA results colored by age
p_pca <- DimPlot(brain,
                 reduction = "pca",
                 group.by = "age_group",
                 pt.size = 0.5) +
  ggtitle("PCA: Young vs Old Brain Cells")

ggsave("results/figures/PCA_by_age.png",
       p_pca, width = 8, height = 6, dpi = 300)

cat("PCA_by_age.png\n")


# The elbow plot shows variance explained by each PC

p_elbow <- ElbowPlot(brain, ndims = 50) +
  ggtitle("Elbow Plot: Variance Explained by Each PC") +
  geom_vline(xintercept = 30, linetype = "dashed", color = "red") +
  annotate("text", x = 35, y = max(brain[["pca"]]@stdev[1:50]),
           label = "Using 30 PCs\n(after elbow)",
           color = "red", hjust = 0)

ggsave("results/figures/elbow_plot.png",
       p_elbow, width = 8, height = 6, dpi = 300)

cat("elbow_plot.png")

# Based on elbow plot, we'll use 30 PCs for downstream analysis
n_pcs_to_use <- 30
cat("  Decision: Using", n_pcs_to_use, "PCs for clustering\n")

# ----------------------------------------------------------------------------
# UMAP for Visualization
# ----------------------------------------------------------------------------

brain <- RunUMAP(brain,
                 dims = 1:n_pcs_to_use,
                 verbose = FALSE)

# Visualize UMAP colored by age
p_umap_age <- DimPlot(brain,
                      reduction = "umap",
                      group.by = "age_group",
                      pt.size = 0.5) +
  ggtitle("UMAP: Brain Cells Colored by Age") +
  theme_minimal()

ggsave("results/figures/UMAP_by_age.png",
       p_umap_age, width = 8, height = 6, dpi = 300)


# ----------------------------------------------------------------------------
# Clustering
# ----------------------------------------------------------------------------

cat("Clustering cells")

# Build neighbor graph
cat("  Building neighbor graph...\n")
brain <- FindNeighbors(brain,
                       dims = 1:n_pcs_to_use,
                       verbose = FALSE)

# Try multiple resolutions
cat("  Testing multiple clustering resolutions...\n")
resolutions_to_try <- c(0.2, 0.4, 0.6, 0.8)

for(res in resolutions_to_try) {
  cat("    Resolution", res, "... ")
  brain <- FindClusters(brain,
                        resolution = res,
                        verbose = FALSE)
  n_clusters <- length(unique(Idents(brain)))
  cat(n_clusters, "clusters\n")
}

cat("\n")

# Visualize different resolutions
p1 <- DimPlot(brain, reduction = "umap", group.by = "RNA_snn_res.0.2", 
              label = TRUE, label.size = 6) + 
  ggtitle("Resolution 0.2") + NoLegend()

p2 <- DimPlot(brain, reduction = "umap", group.by = "RNA_snn_res.0.4",
              label = TRUE, label.size = 6) +
  ggtitle("Resolution 0.4") + NoLegend()

p3 <- DimPlot(brain, reduction = "umap", group.by = "RNA_snn_res.0.6",
              label = TRUE, label.size = 6) +
  ggtitle("Resolution 0.6") + NoLegend()

p4 <- DimPlot(brain, reduction = "umap", group.by = "RNA_snn_res.0.8",
              label = TRUE, label.size = 6) +
  ggtitle("Resolution 0.8") + NoLegend()

p_resolutions <- (p1 + p2) / (p3 + p4)

ggsave("clustering_resolutions.png",
       p_resolutions, width = 14, height = 12, dpi = 300)

cat("  ✓ Saved: 13_clustering_resolutions.png\n\n")

# choice of resolution 0.4 as default 
# You can adjust this based on what looks biological
Idents(brain) <- "RNA_snn_res.0.4"
chosen_resolution <- 0.4
n_final_clusters <- length(unique(Idents(brain)))

cat("  Final choice: Resolution", chosen_resolution, "\n")
cat("  Number of clusters:", n_final_clusters, "\n")
cat("  These clusters likely represent major brain cell types\n\n")

# Create final clustering visualization
p_final_clusters <- DimPlot(brain,
                            reduction = "umap",
                            label = TRUE,
                            label.size = 6,
                            pt.size = 0.5) +
  ggtitle(paste0("Brain Cell Clusters (n=", n_final_clusters, ")")) +
  theme_minimal()

ggsave("results/figures/final_clusters.png",
       p_final_clusters, width = 10, height = 8, dpi = 300)

cat("final_clusters.png")

# Show clusters split by age to see if clustering is driven by age
p_clusters_by_age <- DimPlot(brain,
                             reduction = "umap",
                             split.by = "age_group",
                             label = TRUE,
                             label.size = 4) +
  ggtitle("Clusters in Young vs Old Brain")

ggsave("results/figures/clusters_split_by_age.png",
       p_clusters_by_age, width = 14, height = 6, dpi = 300)

cat("clusters_split_by_age.png")

# ----------------------------------------------------------------------------
# Cluster Composition Analysis
# ----------------------------------------------------------------------------
cluster_sizes <- table(Idents(brain))
cat("  Cluster sizes:\n")
print(cluster_sizes)
cat("\n")

# Cluster composition by age
cluster_age_table <- table(Idents(brain), brain$age_group)
cat("  Cells per cluster by age group:\n")
print(cluster_age_table)
cat("\n")

# Calculate proportions
cluster_age_prop <- prop.table(cluster_age_table, margin = 1)
cat("  Proportion of young vs old in each cluster:\n")
print(round(cluster_age_prop, 3))
cat("\n")

# saving tables
write.csv(cluster_sizes, "results/tables/cluster_sizes.csv")
write.csv(cluster_age_table, "results/tables/cluster_composition_by_age.csv")
write.csv(cluster_age_prop, "results/tables/cluster_age_proportions.csv")


# Visualize cluster composition
composition_df <- as.data.frame(cluster_age_table)
colnames(composition_df) <- c("Cluster", "Age", "Count")

p_composition <- ggplot(composition_df, aes(x = Cluster, y = Count, fill = Age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Cluster Composition by Age",
       y = "Proportion",
       x = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("results/figures/cluster_composition.png",
       p_composition, width = 10, height = 6, dpi = 300)

cat("cluster_composition.png\n\n")