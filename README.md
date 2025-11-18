# Single-Cell RNA-seq Analysis of Brain Aging

A comprehensive single-cell transcriptomic analysis investigating cellular composition changes in the aging mouse brain using the Tabula Muris Senis dataset.

[![R](https://img.shields.io/badge/R-%3E%3D4.3.0-blue)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.0.0-green)](https://satijalab.org/seurat/)

## 🔬 Project Overview

This project analyzes ~17,000 high-quality brain cells from young (3 months) and aged (18-24 months) mice to understand how brain cellular composition changes with aging. The analysis reveals dramatic cell type-specific alterations, including severe oligodendrocyte precursor cell (OPC) depletion and increased microglial activation in aged brain tissue.

### Key Findings

- **Dramatic OPC Depletion:** ~5-fold reduction in oligodendrocyte precursor cells in aged brain, potentially explaining myelin deterioration during aging
- **Microglial Activation:** Increased inflammatory and disease-associated microglia (DAM) populations in aged tissue
- **Senescent Cell Accumulation:** 4-fold increase in senescent cells with aging
- **16 Distinct Cell Populations:** Successfully identified and validated using canonical markers including neurons, astrocytes, oligodendrocytes, microglia, and vascular cells

## 📊 Results Summary

### Quality Control
![image](https://github.com/Aishwarya-Jadeja01/Single-Cell-RNA-seq-Analysis-of-Brain-Aging/blob/main/results/figures/fig1_qc_before_filtering.png)



*High-quality cells retained after stringent filtering (200-6000 genes, <10% mitochondrial content)*

 ![image](https://github.com/Aishwarya-Jadeja01/Single-Cell-RNA-seq-Analysis-of-Brain-Aging/blob/main/results/figures/fig2_qc_after_filtering.png.png)


*Final dataset: 17,000+ cells with tight quality distributions*

### Cell Type Identification
![image](https://github.com/Aishwarya-Jadeja01/Single-Cell-RNA-seq-Analysis-of-Brain-Aging/blob/main/results/figures/fig3_umap_cell_types.png.png)

*16 distinct brain cell populations identified and validated with canonical markers*

![image](https://github.com/Aishwarya-Jadeja01/Single-Cell-RNA-seq-Analysis-of-Brain-Aging/blob/main/results/figures/fig4_umap_by_age.png.png)

*Side-by-side comparison showing cell type distributions in young vs aged brain*

### Aging-Associated Changes
![image](https://github.com/Aishwarya-Jadeja01/Single-Cell-RNA-seq-Analysis-of-Brain-Aging/blob/main/results/figures/fig5_cell_composition.png)

*Stacked bar plot showing proportional changes in cell type composition with age*

![image](https://github.com/Aishwarya-Jadeja01/Single-Cell-RNA-seq-Analysis-of-Brain-Aging/blob/main/results/figures/fig6_celltype_abundance.png)

*Log2 fold-change analysis reveals dramatic OPC depletion and senescent cell accumulation*

### Validation
![image](https://github.com/Aishwarya-Jadeja01/Single-Cell-RNA-seq-Analysis-of-Brain-Aging/blob/main/results/figures/fig10_dotplot_validation.png)

*Canonical marker expression confirms accurate cell type annotations (e.g., P2ry12 for microglia, Pdgfra for OPCs, Gfap for astrocytes)*

## 🧬 Biological Interpretation

### Why OPC Depletion Matters
Oligodendrocyte precursor cells are essential for:
- Myelin maintenance and repair
- Neuronal support and metabolism
- Cognitive function preservation

The 5-fold reduction in aged brain suggests:
- Impaired remyelination capacity
- Potential contribution to cognitive decline
- Novel therapeutic target for age-related neurodegeneration

### Microglial Activation Patterns
Increased inflammatory microglia indicate:
- Chronic neuroinflammation
- Potential contribution to neurodegenerative processes
- Link between aging and Alzheimer's disease pathology

## 📁 Repository Structure

```
aging-brain-singlecell/
│
├── README.md                          # This file
├── requirements.txt                   # R package dependencies
│
├── scripts/                           # Analysis pipeline
│   ├── 01_data_loading.R             # Load raw data, create Seurat object
│   ├── 02_quality_control.R          # QC filtering and visualization
│   ├── 03_doublet_detection.R        # Remove doublets with DoubletFinder
│   ├── 04_normalization_clustering.R # Normalize, PCA, UMAP, clustering
│   ├── 05_cell_type_annotation.R     # Annotate clusters with cell types
│   ├── 06_annotation_validation.R    # Validate with canonical markers
│   ├── 07_aging_analysis.R           # Differential abundance analysis
│   └── 08_pathway_enrichment.R       # Functional enrichment analysis
│
├── figures/                           # Publication-quality figures
│   ├── fig1_qc_before_filtering.png
│   ├── fig2_qc_after_filtering.png
│   ├── fig3_umap_cell_types.png
│   ├── fig4_umap_by_age.png
│   ├── fig5_cell_composition.png
│   ├── fig6_celltype_abundance.png
│   └── fig7_aging_genes_heatmap.png
│   └── fig8_key_aging_genes_violin.png
│   └── fig9_aging_genes_featureplot.png
│   └── fig10_dotplot_validation.png
│
├── tables/                           # Analysis outputs
│   └── summary_statistics.csv
```

## 🚀 Installation & Usage

### Prerequisites

- R >= 4.3.0
- RStudio (recommended)
- Required packages listed in `requirements.txt`

### Installation

1. **Clone this repository:**
```bash
git clone https://github.com/yourusername/aging-brain-singlecell.git
cd aging-brain-singlecell
```

2. **Install R packages:**
```r
# Install CRAN packages
install.packages(c("Seurat", "tidyverse", "patchwork", "readxl", "viridis", 
                   "pheatmap", "scales", "cowplot"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "enrichplot", "org.Mm.eg.db", 
                       "DESeq2", "limma", "ComplexHeatmap"))

# Install DoubletFinder from GitHub
devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
```

3. **Download data:**
   
The analysis uses the Tabula Muris Senis Brain FACS dataset. Download from:
- **Source:** [Tabula Muris Senis](https://tabula-muris-senis.ds.czbiohub.org/)
- **Direct link:** [Brain FACS data](https://figshare.com/projects/Tabula_Muris_Senis/64982)

Place the downloaded files in:
```
data/raw/Brain_FACS_counts.xlsx
data/raw/Brain_FACS_metadata.xlsx
```

4. **Run analysis pipeline:**
```r
# Run scripts sequentially from project root directory
source("scripts/01_data_loading.R")
source("scripts/02_quality_control.R")
source("scripts/03_doublet_detection.R")
source("scripts/04_normalization_clustering.R")
source("scripts/05_cell_type_annotation.R")
source("scripts/06_annotation_validation.R")
source("scripts/07_aging_analysis.R")
source("scripts/08_pathway_enrichment.R")
```

**Important:** Always run scripts from the project root directory to ensure relative paths work correctly.

## 🔍 Methods Overview

### Quality Control
- **Filtering criteria:** 200-6000 genes per cell, <10% mitochondrial content
- **Doublet removal:** DoubletFinder algorithm to remove technical artifacts
- **Final dataset:** ~17,000 high-quality cells

### Normalization & Clustering
- **Normalization:** Log-normalization with scaling
- **Feature selection:** 2000 highly variable genes
- **Dimensionality reduction:** PCA (30 components) + UMAP
- **Clustering:** Louvain algorithm (resolution 0.8)

### Cell Type Annotation
- **Strategy:** Systematic validation using canonical markers
- **Key markers used:**
  - Microglia: P2ry12, Cx3cr1, Tmem119
  - Astrocytes: Gfap, Aqp4, Aldh1l1
  - Oligodendrocytes: Mobp, Plp1, Mog
  - OPCs: Pdgfra, Cspg4
  - Neurons: Snap25, Rbfox3, Syp
  - Endothelial: Cldn5, Pecam1

### Statistical Analysis
- **Differential abundance:** Log2 fold-change (Old/Young proportions)
- **Threshold:** |Log2FC| > 1 considered biologically significant
- **Pathway enrichment:** ClusterProfiler with GO and KEGG databases

## 📈 Key Statistics

| Metric | Value |
|--------|-------|
| Total cells analyzed | 17,436 |
| Young (3m) cells | 8,947 |
| Old (18-24m) cells | 8,489 |
| Cell types identified | 16 |
| Median genes/cell | 2,487 |
| Median UMIs/cell | 7,891 |

## 🎯 Critical Methodological Decisions

### Why This Analysis Is Rigorous

1. **Canonical Marker Validation**
   - Systematically validated each cluster with known markers
   - Corrected misannotations (e.g., initial algorithm misclassified some microglia)

2. **Biological Interpretation First**
   - Every analytical decision driven by understanding and study material help
   - Statistical thresholds informed by literature
   - Focus on reproducibility and transparency

3. **Conservative Quality Control**
   - Stringent filtering removes low-quality cells
   - Doublet detection prevents technical artifacts
   - High-quality final dataset suitable for downstream analysis

## 📚 Data Source

This analysis uses data from:

**Tabula Muris Senis Consortium** (2020). A single-cell transcriptomic atlas characterizes ageing tissues in the mouse. *Nature*, 583(7817), 590-595.

- **Dataset:** Brain FACS (Fluorescence-Activated Cell Sorting)
- **Species:** *Mus musculus*
- **Ages:** 3 months (young adult) vs 18-24 months (aged)
- **Sequencing platform:** Smart-seq2
- **Repository:** [Tabula Muris Senis](https://tabula-muris-senis.ds.czbiohub.org/)

## 🤝 Contributing

This is a portfolio project, but suggestions and feedback are welcome! Please open an issue for:
- Bug reports
- Methodological suggestions
- Analysis questions

## 📄 License

This project is licensed under the MIT License - see LICENSE file for details.

## 👤 Author

**Aishwarya**
- [GitHub](https://github.com/Aishwarya-Jadeja01) |

## 🙏 Acknowledgments

- **Tabula Muris Senis Consortium** for providing open-access aging datasets
- **Satija Lab** for developing Seurat

---

**Last Updated:** November 2024
