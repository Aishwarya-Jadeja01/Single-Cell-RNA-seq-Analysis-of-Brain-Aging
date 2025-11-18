# ============================================================================
# Pathway Enrichment Cell Types 
# ============================================================================

suppressPackageStartupMessages{
  library(Seurat)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Mm.eg.db)  
  library(enrichplot)
  library(ggplot2)



obj_path <- "data/processed/brain_CORRECTED_FINAL.rds"
Idents(brain) <- "age_group"

brain <- readRDS(obj_path)

# Metadata sanity
age_col <- if ("age_group" %in% colnames(brain@meta.data)) "age_group" else if ("age" %in% colnames(brain@meta.data)) "age" else stop("No age column found.")
stopifnot(all(c("Young_3m","Old_18-24m") %in% unique(brain$age_group)))
ct_col <- "cell_type_corrected"
stopifnot(ct_col %in% colnames(brain@meta.data))

#Differential expression per corrected cell type 
message("Running DE per corrected cell type (Old_18-24m vs Young_3m) ==")
cell_types <- sort(unique(brain$cell_type_corrected))
de_list <- list()

for (ct in cell_types) {
  cells_ct <- subset(brain, subset = cell_type_corrected == ct)
  age_tab <- table(cells_ct$age_group)
  
  # skip rare or missing-age cell types
  if (!all(c("Young_3m", "Old_18-24m") %in% names(age_tab)) || sum(age_tab) < 40) {
    message("Skipping ", ct, " (missing ages or too few cells: ", paste(names(age_tab), age_tab, collapse = ", "), ")")
    next
  }
  
  # set identities to age_group inside this subset
  Idents(cells_ct) <- "age_group"
  
  de_res <- tryCatch({
    FindMarkers(
      cells_ct,
      ident.1 = "Old_18-24m",
      ident.2 = "Young_3m",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      test.use = "wilcox"
    ) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::arrange(dplyr::desc(abs(avg_log2FC))) %>%
      dplyr::mutate(cell_type = ct)
  }, error = function(e) {
    message("DE error in ", ct, ": ", e$message)
    NULL
  })
  
  if (!is.null(de_res)) {
    out_csv <- file.path("results/tables", paste0("DE_aging_", gsub("[^A-Za-z0-9]+","_", ct), ".csv"))
    write.csv(de_res, out_csv, row.names = FALSE)
    de_list[[ct]] <- de_res
  }
}

de_all <- dplyr::bind_rows(de_list)
if (nrow(de_all) == 0) stop("No DE results produced.")
write.csv(de_all, "results/tables/DE_aging_all_celltypes_CORRECTED.csv", row.names = FALSE)

# symbol->ENTREZ mapping 
sym2entrez <- function(genes_char) {
  genes_char <- unique(genes_char)
  suppressMessages({
    bitr(genes_char, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  }) %>% distinct(SYMBOL, .keep_all = TRUE)
}

# Enrichment runner per cell type
run_enrichment_for_ct <- function(df_ct, ct_label, suffix="all") {
  #  significant genes (padj < 0.05), then split by log2FC
  sig <- df_ct %>% filter(!is.na(p_val_adj), p_val_adj < 0.05)
  up  <- sig %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
  dn  <- sig %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()
  
  do_enrich <- function(genes, direction) {
    if (length(genes) < 5) return(NULL)
    mapdf <- sym2entrez(genes)
    if (nrow(mapdf) < 5) return(NULL)
    
    # GO 
    ego <- tryCatch({
      enrichGO(
        gene          = mapdf$ENTREZID,
        OrgDb         = org.Mm.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        readable      = TRUE
      )
    }, error = function(e) NULL)
    
    # KEGG
    ekegg <- tryCatch({
      enrichKEGG(
        gene          = mapdf$ENTREZID,
        organism      = "mmu",
        pvalueCutoff  = 0.05
      )
    }, error = function(e) NULL)
    
      
      # symbol input
      ehall <- tryCatch({
        enricher(
          gene          = genes,
          TERM2GENE     = stack(hallmark_list) %>%
            select(values, ind) %>%
            rename(gene=values, set=ind),
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.05
        )
      }, error = function(e) NULL)
    }
    
    list(GO=ego, KEGG=ekegg, HALL=ehall, n_genes=length(genes), direction=direction)
  }
  
  out_up <- do_enrich(up, "Up_in_24m")
  out_dn <- do_enrich(dn, "Down_in_24m")
  
  # Save tables
  save_enr <- function(obj, kind) {
    if (!is.null(obj) && !is.null(obj[[kind]]) && nrow(as.data.frame(obj[[kind]]))>0) {
      df <- as.data.frame(obj[[kind]])
      fn <- paste0("results/tables/ENR_", kind, "_", gsub("[^A-Za-z0-9]+","_", ct_label), "_", suffix, "_", obj$direction, ".csv")
      write.csv(df, fn, row.names = FALSE)
      return(fn)
    }
    return(NULL)
  }
  files <- list(
    GO_up   = save_enr(out_up, "GO"),
    KEGG_up = save_enr(out_up, "KEGG"),
    H_up    = save_enr(out_up, "HALL"),
    GO_dn   = save_enr(out_dn, "GO"),
    KEGG_dn = save_enr(out_dn, "KEGG"),
    H_dn    = save_enr(out_dn, "HALL")
  )
  
  # Dotplots
  plot_enr <- function(obj, kind, topN=15) {
    if (!is.null(obj) && !is.null(obj[[kind]]) && nrow(as.data.frame(obj[[kind]]))>0) {
      p <- suppressWarnings(dotplot(obj[[kind]], showCategory = topN) +
                              ggtitle(paste0(ct_label, " • ", obj$direction, " • ", kind)))
      fp <- paste0("results/figures/1", gsub("[^A-Za-z0-9]+","_", ct_label), "_", suffix, "_", obj$direction, "_", kind, "_dotplot.png")
      ggsave(fp, p, width = 10, height = 7, dpi = 300)
      return(fp)
    }
    return(NULL)
  }
  figs <- list(
    GO_up   = plot_enr(out_up, "GO"),
    KEGG_up = plot_enr(out_up, "KEGG"),
    H_up    = plot_enr(out_up, "HALL"),
    GO_dn   = plot_enr(out_dn, "GO"),
    KEGG_dn = plot_enr(out_dn, "KEGG"),
    H_dn    = plot_enr(out_dn, "HALL")
  )
  
  list(files = files, figs = figs,
       counts = tibble(cell_type = ct_label,
                       up_genes = length(unique(up)),
                       down_genes = length(unique(dn))))
}

#  Run enrichment per cell type
for (ct in cell_types) {
  df_ct <- de_all %>% filter(cell_type == ct)
  if (nrow(df_ct) == 0) { message("No DE for ", ct); next }
  res <- run_enrichment_for_ct(df_ct, ct, "sig")
  if (!is.null(res)) enr_summary[[ct]] <- res$counts
}
enr_counts <- bind_rows(enr_summary)
write.csv(enr_counts, "results/tables/ENR_counts_up_down_by_celltype.csv", row.names = FALSE)
