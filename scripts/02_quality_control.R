# ============================================================================
# Quality Control for Aging Brain Analysis
# ============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)


# DIRECTORIES
cat("Creating output directories...\n")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed",  recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Load Seurat object
# ----------------------------------------------------------------------------
cat("Step 1: Loading Seurat object...\n")
brain <- readRDS("data/processed/brain_seurat_full.rds")
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^MT-")

cat("  Loaded:", ncol(brain), "cells ×", nrow(brain), "genes\n")
cat("  Current gene IDs (first 5):\n")
print(head(rownames(brain), 5))
cat("\n")

cat("Starting dataset:\n")
cat("  Cells:", ncol(brain), "\n")
cat("  Genes:", nrow(brain), "\n\n")
#anndata::read_h5ad()
# ============================================================================
#  Filter to BRAIN tissues only
# ============================================================================

brain_only <- subset(brain, subset = tissue %in% c("Brain_Myeloid", "Brain_Non-Myeloid"))

cat("  Brain cells:", ncol(brain_only), "\n")
cat("  Removed:", ncol(brain) - ncol(brain_only), "non-brain cells\n")
cat("  Retention:", round(100 * ncol(brain_only) / ncol(brain), 1), "%\n\n")

# ============================================================================
#  Simplify age groups 
# ============================================================================

brain_only$age_group <- ifelse(brain_only$age == "3m", "Young_3m", "Old_18-24m")

cat("  Age distribution:\n")
age_dist <- table(brain_only$age_group)
print(age_dist)
cat("\n")

# ============================================================================
#  Calculate QC metrics
# ============================================================================

brain_only[["percent.mt"]] <- PercentageFeatureSet(brain_only, pattern = "^mt-")
brain_only[["percent.ribo"]] <- PercentageFeatureSet(brain_only, pattern = "^Rp[sl]")

cat("QC metrics calculated\n\n")

# ============================================================================
# Visualize QC before filtering
# ============================================================================
pub_theme <- theme_bw(base_size = 16) +
  theme(
    plot.title   = element_text(face = "bold", size = 18),
    axis.text    = element_text(size = 14),
    axis.title.x = element_text(size = 0),    # effectively hide x-title
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.position = "none"
  )

age_labels <- c(
  "Young_3m"    = "Young (3m)",
  "Old_18-24m"  = "Old (18–24m)"
)

# ---- Genes per cell (nFeature_RNA) ----
p_genes_before <- VlnPlot(
  brain_only,
  features = "nFeature_RNA",
  group.by = "age_group",
  pt.size  = 0
) +
  geom_hline(yintercept = c(200, 6000),
             linetype = "dashed",
             color = "red",
             linewidth = 0.5) +
  scale_x_discrete(labels = age_labels) +
  labs(
    title = "QC BEFORE Filtering",
    x     = NULL,
    y     = "Number of Genes Detected"
  ) +
  pub_theme

# ---- UMI counts per cell (nCount_RNA) ----
p_umis_before <- VlnPlot(
  brain_only,
  features = "nCount_RNA",
  group.by = "age_group",
  pt.size  = 0
) +
  scale_x_discrete(labels = age_labels) +
  labs(
    title = NULL,
    x     = NULL,
    y     = "UMI Counts per Cell"
  ) +
  pub_theme

# Combine side by side (left = genes, right = UMIs)
p_before_combined <- p_genes_before + p_umis_before + plot_layout(ncol = 2)

ggsave(
  "results/figures/01_QC_before_filtering.png",
  p_before_combined,
  width  = 12,
  height = 5,
  dpi    = 600
)

cat("01_QC_before_filtering.png")

# ============================================================================
# QC statistics
# ============================================================================
cat("Total cells:", ncol(brain_only), "\n")
cat("Median genes per cell:", median(brain_only$nFeature_RNA), "\n")
cat("Median UMIs per cell:", median(brain_only$nCount_RNA), "\n")
cat("Median % mitochondrial:", round(median(brain_only$percent.mt), 2), "%\n")
cat("Median % ribosomal:", round(median(brain_only$percent.ribo), 2), "%\n\n")

# ============================================================================
# Filter cells based on QC metrics
# ============================================================================
cat("  Filtering rationale:\n")
cat("    - nFeature < 200: Empty droplets/debris\n")
cat("    - nFeature > 6000: Likely doublets (2 cells stuck together)\n")
cat("    - percent.mt > 10%: Dying/stressed cells\n\n")

cat("  Applying filters:\n")
cat("    - Genes per cell: 200 - 6000\n")
cat("    - Mitochondrial content: < 10%\n\n")

brain_filtered <- subset(brain_only, 
                         subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 6000 & 
                           percent.mt < 10)

cat("Remaining cells:", ncol(brain_filtered), "\n")
cat("Cells removed:", ncol(brain_only) - ncol(brain_filtered), "\n")
cat("Percent retained:", round(100 * ncol(brain_filtered) / ncol(brain_only), 1), "%\n\n")

# Age distribution after filtering
cat("Age distribution (after filtering):\n")
age_dist_filtered <- table(brain_filtered$age_group)
print(age_dist_filtered)
cat("\n")

# Retention by age group
cat("Retention by age group:\n")
young_before <- sum(brain_only$age_group == "Young_3m")
young_after <- sum(brain_filtered$age_group == "Young_3m")
cat("  Young: ", young_after, "/", young_before, 
    " (", round(100*young_after/young_before, 1), "%)\n", sep="")

old_before <- sum(brain_only$age_group == "Old_18-24m")
old_after <- sum(brain_filtered$age_group == "Old_18-24m")
cat("  Old: ", old_after, "/", old_before,
    " (", round(100*old_after/old_before, 1), "%)\n\n", sep="")

# ============================================================================
# Visualize after filtering
# ============================================================================

# Genes per cell AFTER filtering
p_genes_after <- VlnPlot(
  brain_filtered,
  features = "nFeature_RNA",
  group.by = "age_group",
  pt.size  = 0
) +
  geom_hline(yintercept = c(200, 6000),
             linetype = "dashed",
             color = "red",
             linewidth = 0.5) +
  scale_x_discrete(labels = age_labels) +
  labs(
    title = "QC AFTER Filtering",
    x     = NULL,
    y     = "Number of Genes Detected"
  ) +
  pub_theme

# UMI counts per cell AFTER filtering
p_umis_after <- VlnPlot(
  brain_filtered,
  features = "nCount_RNA",
  group.by = "age_group",
  pt.size  = 0
) +
  scale_x_discrete(labels = age_labels) +
  labs(
    title = NULL,
    x     = NULL,
    y     = "UMI Counts per Cell"
  ) +
  pub_theme

p_after_combined <- p_genes_after + p_umis_after + plot_layout(ncol = 2)

ggsave(
  "results/figures/03_QC_after_filtering.png",
  p_after_combined,
  width  = 12,
  height = 5,
  dpi    = 600
)

cat("03_QC_after_filtering.png")


qc_summary <- data.frame(
  Metric = c("Total cells", 
             "Young cells", 
             "Old cells", 
             "Median genes/cell", 
             "Median UMIs/cell", 
             "Median mt%",
             "Median ribo%"),
  Before_Filtering = c(ncol(brain_only), 
                       sum(brain_only$age_group == "Young_3m"),
                       sum(brain_only$age_group == "Old_18-24m"),
                       median(brain_only$nFeature_RNA),
                       median(brain_only$nCount_RNA),
                       round(median(brain_only$percent.mt), 2),
                       round(median(brain_only$percent.ribo), 2)),
  After_Filtering = c(ncol(brain_filtered),
                      sum(brain_filtered$age_group == "Young_3m"),
                      sum(brain_filtered$age_group == "Old_18-24m"),
                      median(brain_filtered$nFeature_RNA),
                      median(brain_filtered$nCount_RNA),
                      round(median(brain_filtered$percent.mt), 2),
                      round(median(brain_filtered$percent.ribo), 2))
)

write.csv(qc_summary, "results/tables/qc_summary.csv", row.names = FALSE)
cat("results/tables/qc_summary.csv")

print(qc_summary)

saveRDS(brain_filtered, "data/processed/brain_filtered.rds")
