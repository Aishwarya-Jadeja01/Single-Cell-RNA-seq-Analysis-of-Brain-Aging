# ============================================================================
# Aging Analysis 
# ============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(scales)

setwd("aging-brain-singlecell")


brain <- readRDS("data/processed/brain_annotated.rds")

cat("  Cells:", ncol(brain), "\n")
cat("  Clusters:", length(unique(Idents(brain))), "\n\n")

# Make sure we're using the right clustering
Idents(brain) <- "RNA_snn_res.0.4"

# ----------------------------------------------------------------------------
# Annotate Clusters with Cell Type Names
# ----------------------------------------------------------------------------
cat("cluster annotation based on markers")

# Based on  marker results,

cluster_annotations <- c(
  "0" = "Microglia_Homeostatic",
  "1" = "Microglia_Activated", 
  "2" = "Oligodendrocytes",
  "3" = "OPC_Progenitors",
  "4" = "Astrocytes",
  "5" = "Microglia_DAM",
  "6" = "Oligodendrocytes_Myelinating",
  "7" = "Astrocytes_Reactive",
  "8" = "Endothelial",
  "9" = "Pericytes_SMC",
  "10" = "Macrophages_BAM",
  "11" = "Mixed_Population",
  "12" = "Ependymal",
  "13" = "Rare_CellType1",
  "14" = "Rare_CellType2",
  "15" = "Senescent_Cells"
)

cat("  Cluster annotations:\n")
for(i in names(cluster_annotations)) {
  cat("    Cluster", i, "→", cluster_annotations[i], "\n")
}
cat("\n")

# current cluster identities
current_idents <- as.character(Idents(brain))

# column by mapping
cell_types <- cluster_annotations[current_idents]
names(cell_types) <- colnames(brain) 

# metadata using AddMetaData 
brain <- AddMetaData(brain, metadata = cell_types, col.name = "cell_type")

cat(" Cell type annotation complete\n")
cat("  Cell types assigned:\n")
print(table(brain$cell_type))


#  major analyses
brain$cell_class <- case_when(
  grepl("Microglia|Macrophage", brain$cell_type) ~ "Microglia/Macrophages",
  grepl("Oligodendrocyte|OPC", brain$cell_type) ~ "Oligodendrocyte Lineage",
  grepl("Astrocyte", brain$cell_type) ~ "Astrocytes",
  grepl("Endothelial|Pericyte", brain$cell_type) ~ "Vascular",
  TRUE ~ "Other"
)
print(table(brain$cell_class))


# ----------------------------------------------------------------------------
# Annotated UMAP Figures
# ----------------------------------------------------------------------------

# Main UMAP with cell type labels
p_umap_annotated <- DimPlot(brain, 
                            reduction = "umap",
                            group.by = "cell_type",
                            label = TRUE,
                            label.size = 3,
                            repel = TRUE,
                            pt.size = 0.3) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Brain Cell Types in Aging",
       subtitle = paste(ncol(brain), "cells across 16 populations"))

ggsave("results/figures/UMAP_annotated_celltypes.png",
       p_umap_annotated, width = 12, height = 8, dpi = 300)

cat("UMAP_annotated_celltypes.png\n")

# UMAP split by age with annotations
p_umap_age_split <- DimPlot(brain,
                            reduction = "umap",
                            group.by = "cell_type",
                            split.by = "age_group",
                            label = TRUE,
                            label.size = 2.5,
                            repel = TRUE,
                            pt.size = 0.2,
                            ncol = 2) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Cell Types: Young vs Old Brain")

ggsave("results/figures/UMAP_age_split_annotated.png",
       p_umap_age_split, width = 14, height = 6, dpi = 300)

cat("UMAP_age_split_annotated.png\n")

# UMAP by broad cell class
p_umap_class <- DimPlot(brain,
                        reduction = "umap",
                        group.by = "cell_class",
                        pt.size = 0.5,
                        cols = brewer.pal(5, "Set2")) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Major Brain Cell Classes")

ggsave("results/figures/UMAP_cell_classes.png",
       p_umap_class, width = 10, height = 7, dpi = 300)

cat("UMAP_cell_classes.png\n\n")

# ----------------------------------------------------------------------------
# Cell Type Composition Analysis
# ----------------------------------------------------------------------------

#  cell type proportions
celltype_counts <- table(brain$cell_type, brain$age_group)
celltype_props <- prop.table(celltype_counts, margin = 2)

#  dataframe for visualization
comp_df <- as.data.frame(celltype_counts)
colnames(comp_df) <- c("CellType", "Age", "Count")

comp_df <- comp_df %>%
  group_by(Age) %>%
  mutate(Proportion = Count / sum(Count))

# Stacked bar plot
p_composition <- ggplot(comp_df, aes(x = Age, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(16)) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 11)) +
  labs(title = "Brain Cell Type Composition: Young vs Old",
       y = "Proportion",
       x = "Age Group",
       fill = "Cell Type") +
  scale_y_continuous(labels = percent)

ggsave("results/figures/cell_composition_by_age.png",
       p_composition, width = 10, height = 8, dpi = 300)


fc_df <- as.data.frame(celltype_props)
colnames(fc_df) <- c("CellType", "Age", "Proportion")


fc_wide <- fc_df %>%
  pivot_wider(names_from = Age, values_from = Proportion)


cat("\n  Debug - Column names after pivot:\n")
print(colnames(fc_wide))

# Extract the two age columns dynamically
age_columns <- colnames(fc_wide)[colnames(fc_wide) != "CellType"]
cat("  Age columns found:", paste(age_columns, collapse = ", "), "\n\n")

# fold change using column positions
if(length(age_columns) == 2) {
  # Determine which is old and which is young
  old_idx <- which(grepl("Old|old|18|24", age_columns))
  young_idx <- which(grepl("Young|young|3m", age_columns))
  
  if(length(old_idx) == 0) old_idx <- 1
  if(length(young_idx) == 0) young_idx <- 2
  
  cat("  Using column", age_columns[old_idx], "as Old\n")
  cat("  Using column", age_columns[young_idx], "as Young\n\n")
  
  fc_wide$Old_Proportion <- fc_wide[[age_columns[old_idx]]]
  fc_wide$Young_Proportion <- fc_wide[[age_columns[young_idx]]]
  fc_wide$FoldChange <- fc_wide$Old_Proportion / fc_wide$Young_Proportion
  fc_wide$Log2FC <- log2(fc_wide$FoldChange)
}

# Bar plot of fold changes
p_fc <- ggplot(fc_wide, aes(x = reorder(CellType, Log2FC), y = Log2FC)) +
  geom_bar(stat = "identity", 
           aes(fill = Log2FC > 0),
           show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("TRUE" = "red3", "FALSE" = "steelblue")) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 9)) +
  labs(title = "Cell Type Abundance Changes with Age",
       subtitle = "Log2 Fold Change (Old/Young proportion)",
       x = "Cell Type",
       y = "Log2(Old/Young)")

ggsave("results/figures/celltype_abundance_foldchange.png",
       p_fc, width = 10, height = 8, dpi = 300)

cat("celltype_abundance_foldchange.png\n\n")

# composition data
write.csv(fc_wide, "results/tables/celltype_abundance_changes.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
#  Aging Genes Across Cell Types (Heatmap)
# ----------------------------------------------------------------------------

age_de_files <- list.files("results/tables", 
                           pattern = "age_DE_cluster_.*\\.csv",
                           full.names = TRUE)

if(length(age_de_files) > 0) {
  # Combine all age DE results
  all_age_de <- lapply(age_de_files, read.csv) %>%
    bind_rows()
  
  # Get top aging genes (both up and down)
  top_aging_genes <- all_age_de %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    slice_head(n = 5) %>%
    pull(gene) %>%
    unique()
  
  if(length(top_aging_genes) > 0 & length(top_aging_genes) <= 100) {
    # Calculate average expression per cell type and age
    avg_exp <- AverageExpression(brain,
                                 features = top_aging_genes,
                                 group.by = c("cell_type", "age_group"),
                                 slot = "data")
    
    # Extract matrix
    mat <- avg_exp$RNA
    
    # Create heatmap
    png("results/figures/aging_genes_heatmap.png", 
        width = 12, height = 14, units = "in", res = 300)
    
    pheatmap(mat,
             scale = "row",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             main = "Top Aging-Associated Genes Across Cell Types",
             fontsize_row = 6,
             fontsize_col = 8,
             border_color = NA)
    
    dev.off()
    
    cat("aging_genes_heatmap.png\n\n")
  }
}

# ----------------------------------------------------------------------------
# Key Aging Genes (H2 genes, inflammatory markers)
# ----------------------------------------------------------------------------

# Define key aging genes
key_aging_genes <- c(
  # MHC-II / Antigen presentation
  "H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74",
  # Inflammation
  "Il1b", "Tnf", "Ccl2", "Ccl5", "Cxcl10",
  # Aging markers
  "Apoe", "Cst7", "Lpl", "Axl",
  # Interferon response
  "Ifit3", "Stat1", "Irf7",
  # Stress/senescence
  "Cdkn1a", "Cdkn2a"
)

# Check which are present
key_genes_present <- key_aging_genes[key_aging_genes %in% rownames(brain)]

cat("  Key aging genes present:", length(key_genes_present), 
    "out of", length(key_aging_genes), "\n\n")

if(length(key_genes_present) > 0) {
  # Violin plots by age for each gene
  for(gene in head(key_genes_present, 8)) {  # Top 8 for space
    p <- VlnPlot(brain,
                 features = gene,
                 group.by = "cell_type",
                 split.by = "age_group",
                 pt.size = 0) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "top") +
      labs(title = paste("Expression of", gene, "in Young vs Old"))
    
    assign(paste0("p_", gene), p)
  }
  
  # Combine plots
  if(length(key_genes_present) >= 4) {
    p_aging_genes <- wrap_plots(
      mget(paste0("p_", head(key_genes_present, 4))),
      ncol = 2
    )
    
    ggsave("results/figures/key_aging_genes_violin.png",
           p_aging_genes, width = 14, height = 12, dpi = 300)
    
    cat("key_aging_genes_violin.png\n")
  }
  
  # Feature plots for top aging markers
  if(length(key_genes_present) >= 4) {
    p_features <- FeaturePlot(brain,
                              features = head(key_genes_present, 4),
                              split.by = "age_group",
                              ncol = 4,
                              pt.size = 0.1)
    
    ggsave("results/figures/aging_genes_featureplot.png",
           p_features, width = 16, height = 8, dpi = 300)
    
    cat("aging_genes_featureplot.png\n\n")
  }
}

# ----------------------------------------------------------------------------
# Cell-Type-Specific Aging Signatures 
# ----------------------------------------------------------------------------

# Load all age DE results
if(length(age_de_files) > 0) {
  aging_summary <- data.frame(
    CellType = character(),
    Cluster = character(),
    Upregulated_Genes = integer(),
    Downregulated_Genes = integer(),
    Top_Upregulated = character(),
    Top_Downregulated = character(),
    stringsAsFactors = FALSE
  )
  
  for(file in age_de_files) {
    de_data <- read.csv(file)
    cluster_id <- gsub(".*cluster_(\\d+)\\.csv", "\\1", file)
    
    # Get cell type name
    celltype <- cluster_annotations[cluster_id]
    
    # Count up/down genes
    n_up <- sum(de_data$avg_log2FC > 0 & de_data$p_val_adj < 0.05)
    n_down <- sum(de_data$avg_log2FC < 0 & de_data$p_val_adj < 0.05)
    
    # Get top genes
    top_up <- de_data %>%
      filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
      arrange(desc(avg_log2FC)) %>%
      head(3) %>%
      pull(gene) %>%
      paste(collapse = ", ")
    
    top_down <- de_data %>%
      filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
      arrange(avg_log2FC) %>%
      head(3) %>%
      pull(gene) %>%
      paste(collapse = ", ")
    
    aging_summary <- rbind(aging_summary, data.frame(
      CellType = celltype,
      Cluster = cluster_id,
      Upregulated_Genes = n_up,
      Downregulated_Genes = n_down,
      Top_Upregulated = top_up,
      Top_Downregulated = top_down
    ))
  }
  
  write.csv(aging_summary, 
            "results/tables/aging_signatures_summary.csv",
            row.names = FALSE)
  
  cat("  Aging signatures by cell type:\n\n")
  print(aging_summary[, 1:4])
  cat("\n")
  
  cat("aging_signatures_summary.csv\n\n")
}

saveRDS(brain, "data/processed/brain_final_annotated.rds")

