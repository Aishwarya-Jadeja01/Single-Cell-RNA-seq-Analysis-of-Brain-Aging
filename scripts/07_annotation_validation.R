# ============================================================================
#  Cell Type Annotations Based on Canonical Marker Validation
# ============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)

setwd("aging-brain-singlecell")


# Load final annotated object
brain <- readRDS("data/processed/brain_FINAL_annotated.rds")

# Set identities to cluster numbers
Idents(brain) <- "RNA_snn_res.0.4"

# annotations based on canonical marker validation
cluster_annotations_corrected <- c(
  "0" = "Microglia",                    # P2ry12+, P2ry13+ (homeostatic)
  "1" = "Microglia_Activated",          # Likely activated microglia
  "2" = "Endothelial",                  # Cldn5+ (CORRECTED)
  "3" = "Proliferating_Cells",          # Histones+ (cycling cells)
  "4" = "Oligodendrocytes",             # Mog+, Cldn11+, Opalin+ (CORRECTED)
  "5" = "Microglia_Inflammatory",       # Likely inflammatory microglia
  "6" = "Oligodendrocytes_Immature",    # Oligodendrocyte subset
  "7" = "Neurons",                      # Rbfox3+ (CORRECTED - rare!)
  "8" = "Astrocytes",                   # Aqp4+ (CORRECTED)
  "9" = "Pericytes",                    # Vascular mural cells
  "10" = "Microglia_Border",            # Border-associated macrophages
  "11" = "Smooth_Muscle_Cells",         # Acta2+ (vascular SMC)
  "12" = "Ependymal_Choroid",           # Ependymal/choroid plexus
  "13" = "Rare_Population_1",           # Small rare cluster
  "14" = "OPC",                         # Pdgfra+ (CORRECTED - OPCs!)
  "15" = "Age_Associated_Cells"         # 95% old-specific
)


# corrected annotations
current_clusters <- as.character(Idents(brain))
cell_types_corrected <- cluster_annotations_corrected[current_clusters]
names(cell_types_corrected) <- colnames(brain)

brain <- AddMetaData(brain, metadata = cell_types_corrected, col.name = "cell_type_corrected")

# corrections
cat("=== CORRECTED ANNOTATIONS ===\n\n")
for(i in names(cluster_annotations_corrected)) {
  cat("Cluster", i, ":", cluster_annotations_corrected[i], "\n")
}
cat("\n")

# cell counts
cat("Cell type distribution:\n")
print(table(brain$cell_type_corrected))
cat("\n")

# Create broad categories
brain$cell_class_corrected <- case_when(
  grepl("Microglia", brain$cell_type_corrected) ~ "Microglia",
  grepl("Oligodendrocyte|OPC", brain$cell_type_corrected) ~ "Oligodendrocyte_Lineage",
  grepl("Astrocyte", brain$cell_type_corrected) ~ "Astrocytes",
  grepl("Endothelial|Pericyte|Smooth_Muscle", brain$cell_type_corrected) ~ "Vascular",
  grepl("Neuron", brain$cell_type_corrected) ~ "Neurons",
  TRUE ~ "Other"
)

cat("Broad cell classes:\n")
print(table(brain$cell_class_corrected))
cat("\n")

# UMAP
p_corrected <- DimPlot(brain,
                       reduction = "umap",
                       group.by = "cell_type_corrected",
                       label = TRUE,
                       label.size = 3,
                       repel = TRUE,
                       pt.size = 0.3) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Cell Type Annotations",
       subtitle = "Based on canonical marker validation")

ggsave("results/figures/UMAP_annotations.png",
       p_corrected, width = 12, height = 8, dpi = 300)


# Split by age
p_corrected_split <- DimPlot(brain,
                             reduction = "umap",
                             group.by = "cell_type_corrected",
                             split.by = "age_group",
                             label = TRUE,
                             label.size = 2,
                             repel = TRUE,
                             pt.size = 0.2,
                             ncol = 2) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "CORRECTED Annotations: Young vs Old")

ggsave("results/figures/UMAP_CORRECTED_by_age.png",
       p_corrected_split, width = 14, height = 6, dpi = 300)


# validation DotPlot with corrected labels
canonical_markers_all <- c(
  # Microglia
  "P2ry12", "P2ry13", "Cx3cr1", "Tmem119",
  # Oligodendrocytes
  "Mbp", "Mog", "Plp1", "Cldn11", "Opalin",
  # OPCs
  "Pdgfra", "Cspg4", "Sox10",
  # Astrocytes
  "Gfap", "Aldh1l1", "Aqp4", "Aldoc",
  # Endothelial
  "Cldn5", "Pecam1", "Flt1",
  # Pericytes/SMC
  "Pdgfrb", "Rgs5", "Acta2",
  # Neurons
  "Snap25", "Rbfox3", "Tubb3"
)

markers_present <- canonical_markers_all[canonical_markers_all %in% rownames(brain)]

# cell types on y-axis
Idents(brain) <- factor(
  brain$cell_type_corrected,
  levels = sort(unique(brain$cell_type_corrected))
)

dot_theme <- theme_bw(base_size = 18) +
  theme(
    plot.title      = element_text(size = 20, face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y     = element_text(size = 14),
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    panel.grid      = element_blank(),
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 14)
  )

p_dotplot <- DotPlot(
  brain,
  features  = markers_present,
  cols      = c("grey90", "red"),
  dot.scale = 7
) +
  RotatedAxis() +
  scale_color_gradientn(
    colours = c("#f7fbff", "#fc9272", "#cb181d")
  ) +
  scale_size(range = c(1.5, 8)) +
  guides(
    colour = guide_colourbar(title = "Scaled expression"),
    size   = guide_legend(title = "Fraction of cells")
  ) +
  labs(title = "Canonical markers validating cell types") +
  dot_theme

ggsave(
  "results/figures/dotplot_CORRECTED_validation.png",
  p_dotplot,
  width  = 18,         
  height = 10,
  dpi    = 600,
  type   = "cairo"      
)

# Vector PDF (for manuscripts)
ggsave(
  "results/figures/dotplot_CORRECTED_validation.pdf",
  p_dotplot,
  width  = 18,
  height = 10,
  device = cairo_pdf
)


# Recalculate composition with corrected annotations
comp_corrected <- as.data.frame(table(brain$cell_type_corrected, brain$age_group))
colnames(comp_corrected) <- c("CellType", "Age", "Count")

comp_corrected <- comp_corrected %>%
  group_by(Age) %>%
  mutate(Proportion = Count / sum(Count))

p_comp <- ggplot(comp_corrected, aes(x = Age, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(16)) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +
  labs(title = "Cell Composition - Annotations",
       y = "Proportion") +
  scale_y_continuous(labels = scales::percent)

ggsave("results/figures/composition_CORRECTED.png",
       p_comp, width = 10, height = 7, dpi = 300)

# Save corrected object
saveRDS(brain, "data/processed/brain_CORRECTED_FINAL.rds")