# ============================================================================
# Cell Type Annotation & Aging Analysis
# ============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)

setwd("aging-brain-singlecell")


# ----------------------------------------------------------------------------
# Load clustered data 
# ----------------------------------------------------------------------------
brain <- readRDS("data/processed/brain_clustered.rds")

cat("  Cells:", ncol(brain), "\n")
cat("  Clusters:", length(unique(Idents(brain))), "\n\n")

# Set to use resolution 0.4 clusters
Idents(brain) <- "RNA_snn_res.0.4"


cat("canonical cell type markers...\n\n")

# Define markers
markers <- c(
  # Astrocytes
  "Gfap", "Aldh1l1", "Aqp4",
  # Endothelial
  "Cldn5", "Kdr", "Esam", "Flt1",
  # Pericytes/SMC
  "Pdgfrb", "Rgs5", "Acta2",
  # Microglia
  "P2ry12", "Cx3cr1", "Aif1", "Tmem119",
  # Oligodendrocytes
  "Mbp", "Cldn11", "Mog", "Plp1",
  # OPCs
  "Sox10", "Pdgfra", "Cspg4",
  # Neurons
  "Slc17a7", "Rbfox3", "Snap25"
)

# Check which markers are present
markers_present <- markers[markers %in% rownames(brain)]
cat("  Found", length(markers_present), "out of", length(markers), "markers\n\n")

# Create improved DotPlot
p_dotplot <- DotPlot(brain, 
                     features = markers_present,
                     cols = c("lightgrey", "red"),
                     dot.scale = 8) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 11),
        legend.position = "right") +
  labs(title = "Canonical Cell Type Markers Across Clusters",
       x = "Gene", y = "Cluster")

ggsave("results/figures/canonical_markers_dotplot.png",
       p_dotplot, width = 14, height = 8, dpi = 300)

cat("canonical_markers_dotplot.png")

# Feature plots for key markers
key_markers <- c("Cx3cr1", "Mbp", "Gfap", "Cldn5")
key_markers_present <- key_markers[key_markers %in% rownames(brain)]

if(length(key_markers_present) > 0) {
  p_features <- FeaturePlot(brain, 
                            features = key_markers_present,
                            ncol = 2,
                            pt.size = 0.1) &
    theme(legend.position = "right")
  
  ggsave("results/figures/key_markers_featureplot.png",
         p_features, width = 12, height = 10, dpi = 300)
  
  cat("key_markers_featureplot.png")
}

# ----------------------------------------------------------------------------
# FIND MARKERS FOR CLUSTER
# ----------------------------------------------------------------------------

all.markers <- FindAllMarkers(brain, 
                              only.pos = TRUE, 
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              verbose = FALSE)

cat("  Total markers:", nrow(all.markers), "\n\n")

# Save all markers
write.csv(all.markers, 
          "results/tables/markers_per_cluster.csv", 
          row.names = FALSE)

# Get top per cluster
top10_per_cluster <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(top10_per_cluster,
          "results/tables/top10_markers_per_cluster.csv",
          row.names = FALSE)

for(clust in sort(unique(all.markers$cluster))) {
  top3 <- all.markers %>%
    filter(cluster == clust) %>%
    arrange(desc(avg_log2FC)) %>%
    head(3)
  
  cat("  Cluster", clust, ":", 
      paste(top3$gene, collapse = ", "), "\n")
}
cat("\n")

cat("markers_per_cluster.csv\n")
cat("top10_markers_per_cluster.csv\n\n")

# Calculate proportions
cluster_counts <- table(Idents(brain), brain$age_group)
cluster_props <- prop.table(cluster_counts, margin = 2)

cat("  Cluster proportions by age:\n")
print(round(cluster_props, 3))
cat("\n")

# Chi-square test per cluster
chi_results <- data.frame(
  cluster = character(),
  chi_sq = numeric(),
  p_value = numeric(),
  old_prop = numeric(),
  young_prop = numeric(),
  fold_change = numeric(),
  stringsAsFactors = FALSE
)

for(clust in rownames(cluster_counts)) {
  # Create contingency table for this cluster vs all others
  this_cluster <- cluster_counts[clust, ]
  other_clusters <- colSums(cluster_counts) - this_cluster
  
  contingency <- rbind(this_cluster, other_clusters)
  
  # Chi-square test
  chi_test <- chisq.test(contingency)
  
  # Calculate fold change (old/young proportion)
  old_prop <- cluster_props[clust, "Old_18-24m"]
  young_prop <- cluster_props[clust, "Young_3m"]
  fc <- old_prop / young_prop
  
  chi_results <- rbind(chi_results, data.frame(
    cluster = clust,
    chi_sq = chi_test$statistic,
    p_value = chi_test$p.value,
    old_prop = old_prop,
    young_prop = young_prop,
    fold_change = fc
  ))
}

# Adjust p-values for multiple testing
chi_results$p_adj <- p.adjust(chi_results$p_value, method = "BH")

# Sort by significance
chi_results <- chi_results %>% arrange(p_adj)

write.csv(chi_results,
          "results/tables/cluster_abundance_age_test.csv",
          row.names = FALSE)

cat("Age-biased clusters (p.adj < 0.05):\n\n")
sig_clusters <- chi_results %>% filter(p_adj < 0.05)
if(nrow(sig_clusters) > 0) {
  print(sig_clusters[, c("cluster", "fold_change", "p_adj")])
} else {
  cat("No significant changes detected\n")
}
cat("\n")

chi_results$direction <- ifelse(chi_results$fold_change > 1, "Old-enriched", "Young-enriched")
chi_results$significant <- ifelse(chi_results$p_adj < 0.05, "Significant", "Not significant")

p_abundance <- ggplot(chi_results, 
                      aes(x = reorder(cluster, fold_change), 
                          y = log2(fold_change),
                          fill = significant)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "red") +
  scale_fill_manual(values = c("Significant" = "red3", 
                               "Not significant" = "grey70")) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Age-Biased Cluster Abundance",
       subtitle = "Log2(Old proportion / Young proportion)",
       x = "Cluster",
       y = "Log2 Fold Change (Old/Young)",
       fill = "FDR < 0.05") +
  theme(legend.position = "top",
        axis.text = element_text(size = 11),
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/age_biased_abundance.png",
       p_abundance, width = 8, height = 10, dpi = 300)

cat("age_biased_abundance.png\n")
cat("cluster_abundance_age_test.csv\n\n")


# Focus on largest and most age-biased clusters
clusters_to_test <- c(0, 1, 3, 5)  

age_de_results <- list()

for(clust in clusters_to_test) {
  cat("  Analyzing cluster", clust, "...\n")
  
  # Subset to this cluster
  cells_in_cluster <- colnames(brain)[Idents(brain) == clust]
  
  if(length(cells_in_cluster) < 50) {
    cat("    Skipping - too few cells\n")
    next
  }
  
  brain_subset <- subset(brain, cells = cells_in_cluster)
  
  # Find age DE genes
  Idents(brain_subset) <- "age_group"
  
  tryCatch({
    de_genes <- FindMarkers(brain_subset,
                            ident.1 = "Old_18-24m",
                            ident.2 = "Young_3m",
                            min.pct = 0.1,
                            logfc.threshold = 0.25,
                            verbose = FALSE)
    
    de_genes$gene <- rownames(de_genes)
    de_genes$cluster <- clust
    
    age_de_results[[as.character(clust)]] <- de_genes
    
    n_up <- sum(de_genes$avg_log2FC > 0 & de_genes$p_val_adj < 0.05)
    n_down <- sum(de_genes$avg_log2FC < 0 & de_genes$p_val_adj < 0.05)
    
    cat("    Upregulated in old:", n_up, "genes\n")
    cat("    Downregulated in old:", n_down, "genes\n")
    
    # Save individual cluster results
    write.csv(de_genes,
              paste0("results/tables/age_DE_cluster_", clust, ".csv"),
              row.names = FALSE)
    
  }, error = function(e) {
    cat("    Error:", e$message, "\n")
  })
}

# Combine all age DE results
if(length(age_de_results) > 0) {
  age_de_combined <- bind_rows(age_de_results)
  write.csv(age_de_combined,
            "results/tables/age_DE_all_tested_clusters.csv",
            row.names = FALSE)
  
  cat(length(age_de_results), "clusters\n\n")
}


# Define aging gene signatures
aging_signatures <- list(
  Interferon = c("Ifit1", "Ifit3", "Irf7", "Stat1", "Mx1", "Oas1a"),
  Inflammation = c("Il1b", "Tnf", "Il6", "Cxcl10", "Ccl2", "Ccl5"),
  MHC_II = c("H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74"),
  Senescence = c("Cdkn1a", "Cdkn2a", "Il6", "Cxcl1"),
  Stress = c("Hspa1a", "Hspa1b", "Hsp90aa1", "Hspb1")
)

# Check which genes are present
for(sig_name in names(aging_signatures)) {
  genes_present <- aging_signatures[[sig_name]][
    aging_signatures[[sig_name]] %in% rownames(brain)
  ]
  cat("  ", sig_name, ":", length(genes_present), "of", 
      length(aging_signatures[[sig_name]]), "genes present\n")
}
cat("\n")

# Score cells for each signature
for(sig_name in names(aging_signatures)) {
  genes_present <- aging_signatures[[sig_name]][
    aging_signatures[[sig_name]] %in% rownames(brain)
  ]
  
  if(length(genes_present) > 2) {
    brain <- AddModuleScore(brain,
                            features = list(genes_present),
                            name = paste0(sig_name, "_Score"))
  }
}

# Visualize signature scores by age and cluster
signature_cols <- grep("_Score1$", colnames(brain@meta.data), value = TRUE)

if(length(signature_cols) > 0) {
  for(i in 1:length(signature_cols)) {
    new_name <- gsub("_Score1$", "_Score", signature_cols[i])
    brain@meta.data[[new_name]] <- brain@meta.data[[signature_cols[i]]]
  }
  
  signature_cols_clean <- grep("_Score$", colnames(brain@meta.data), value = TRUE)
  
  # violin plots
  plot_list <- list()
  for(sig in signature_cols_clean) {
    p <- VlnPlot(brain, features = sig, group.by = "age_group", pt.size = 0) +
      ggtitle(gsub("_Score", "", sig)) +
      theme(legend.position = "none")
    plot_list[[sig]] <- p
  }
  
  p_signatures <- wrap_plots(plot_list, ncol = 2)
  
  ggsave("results/figures/aging_signature_scores.png",
         p_signatures, width = 10, height = 12, dpi = 300)
  
  cat("aging_signature_scores.png\n\n")
}



#Variable Features
top20_variable <- head(VariableFeatures(brain), 20)

p_variable_improved <- VariableFeaturePlot(brain) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(title = "Highly Variable Genes in Aging Brain",
       subtitle = paste("Top 2000 genes |", 
                        "Red = Selected for downstream analysis"))

p_variable_improved <- LabelPoints(plot = p_variable_improved,
                                   points = top20_variable,
                                   repel = TRUE,
                                   max.overlaps = 20,
                                   size = 3)

ggsave("results/figures/variable_features_IMPROVED.png",
       p_variable_improved, width = 12, height = 8, dpi = 300)

cat("variable_features_IMPROVED.png\n")

# Elbow Plot 
pca_variance <- brain[["pca"]]@stdev^2
pca_var_explained <- pca_variance / sum(pca_variance) * 100

elbow_data <- data.frame(
  PC = 1:50,
  StDev = brain[["pca"]]@stdev[1:50],
  VarExplained = pca_var_explained[1:50],
  CumVar = cumsum(pca_var_explained[1:50])
)

p_elbow_improved <- ggplot(elbow_data, aes(x = PC, y = StDev)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_vline(xintercept = 30, linetype = "dashed", 
             color = "red", size = 1) +
  annotate("text", x = 35, y = max(elbow_data$StDev) * 0.9,
           label = "Using 30 PCs\n(captures major variation)",
           color = "red", hjust = 0, size = 4, fontface = "bold") +
  annotate("text", x = 30, y = max(elbow_data$StDev) * 0.7,
           label = paste0(round(elbow_data$CumVar[30], 1), 
                          "% variance\nexplained"),
           hjust = 1.1, size = 3.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "PCA Variance Explained (Elbow Plot)",
       subtitle = "Each PC captures decreasing amounts of variation",
       x = "Principal Component",
       y = "Standard Deviation")

ggsave("results/figures/elbow_plot_IMPROVED.png",
       p_elbow_improved, width = 10, height = 6, dpi = 300)

cat("elbow_plot_IMPROVED.png\n\n")

saveRDS(brain, "data/processed/brain_annotated.rds")
