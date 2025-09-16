## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile")
library(Seurat)
packageVersion("Seurat")
library(parallel)
library(purrr)
library(tibble)
library(presto)
library(dplyr)
library(patchwork)
library(plyr)
library(RColorBrewer)
library(multtest)
library(metap) 
library(ggprism)
library(glmGamPoi)
library(msigdbr)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(forcats)
library(EnhancedVolcano)
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(SingleCellExperiment))
#BiocManager::install("scran")
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))


## Load the object ----
T1D_Timepoints<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v6.rds")

DimPlot(T1D_Timepoints,group.by = "CellSubType")
#T1D_Timepoints <- UpdateSeuratObject(T1D_Timepoints)

table(T1D_Timepoints$sample,T1D_Timepoints$group)
# # Change 'Week6_6' group from 'Non-Progressor' to 'Progressor'
# T1D_Timepoints$group[T1D_Timepoints$sample == "Week6_6"] <- "Progressor"
# # Check that the change has been made
# table(T1D_Timepoints$sample, T1D_Timepoints$group)
# saveRDS(T1D_Timepoints , file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v6.rds")
# Convert to character to allow new labels
T1D_Timepoints@meta.data$CellType <- as.character(T1D_Timepoints@meta.data$CellSubType)

# Modify labels
T1D_Timepoints@meta.data$CellType[T1D_Timepoints@meta.data$CellSubType %in% c("CD8 memory", "CD8 exhausted effector-like")] <- "CD8 T Cell"

T1D_Timepoints@meta.data$CellType[T1D_Timepoints@meta.data$CellSubType %in% c("Tcon memory", "Tcon exhausted effector-like", "Tregs", "Tcon activated ", "Tcon Interferon Sensing")] <- "CD4 T Cell"

T1D_Timepoints@meta.data$CellType[T1D_Timepoints@meta.data$CellSubType %in% c("ILC2", "ILC3")] <- "ILC"

# Convert back to factor if needed
T1D_Timepoints@meta.data$CellType <- factor(T1D_Timepoints@meta.data$CellType)

unique(T1D_Timepoints$CellType)

# Convert time to character to avoid factor level issues
T1D_Timepoints@meta.data$time <- as.character(T1D_Timepoints@meta.data$time)

# Rename Week6 to Early and Week12 to Intermediate
T1D_Timepoints@meta.data$time[T1D_Timepoints@meta.data$time == "Week6"] <- "Early"
T1D_Timepoints@meta.data$time[T1D_Timepoints@meta.data$time == "Week12"] <- "Intermediate"

# Convert back to factor if necessary
T1D_Timepoints@meta.data$time <- factor(T1D_Timepoints@meta.data$time)
unique(T1D_Timepoints$time)

# keep <- colnames(T1D_Timepoints)[!(T1D_Timepoints$CellType %in% unwanted_cells)]
# length(keep)                                  # 55459
# length(intersect(keep, Cells(T1D_Timepoints[["SCT"]])))  # <-- very likely 0
# DefaultAssay(T1D_Timepoints)                  # "SCT"


# Remove unwanted cell types
unwanted_cells <- c("Gamma Delta T Cell", "Mesenchymal-like", "Acinar Cell")



# Keep only the desired cells
T1D_Timepoints <- subset(T1D_Timepoints, subset = !(CellSubType %in% unwanted_cells))

# Check remaining cell types
unique(T1D_Timepoints$CellType)

## DEG Analysis----



# Get unique cell types and time points
cell_types <- unique(T1D_Timepoints@meta.data$CellType)
time_points <- unique(T1D_Timepoints@meta.data$time)
out_dir <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/"
# Define safe filename helper
safe <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}
# Create an empty list to store results
DEG_results <- list()
DEG_results_stats <- list()   ### NEW: list to hold DEG tables with stats
markers$avg_log2FC
# Loop over each cell type and time point
for (cell in cell_types) {
  for (tp in time_points) {
    cell_time_subset <- subset(T1D_Timepoints, subset = (CellType == cell & time == tp))
    Idents(cell_time_subset) <- cell_time_subset$group
    
    # size check
    grp_counts <- table(cell_time_subset$group)
    n_prog <- as.integer(grp_counts["Progressor"])
    n_non  <- as.integer(grp_counts["Non-Progressor"])
    
    if (is.na(n_prog)) n_prog <- 0
    if (is.na(n_non))  n_non  <- 0
    
    if (n_prog < 3 || n_non < 3) {
      warning(sprintf("Skipping %s | %s: Progressor=%d, Non-Progressor=%d (<3).",
                      cell, tp, n_prog, n_non))
      next
    }
    
    markers <- FindMarkers(
      cell_time_subset,
      ident.1 = "Progressor", ident.2 = "Non-Progressor",
      assay = "SCT", recorrect_umi = FALSE
    )
    markers$p_val_adj <- p.adjust(markers$p_val, method = "BH")
    
    upregulated   <- rownames(markers[markers$avg_log2FC >=  1 & markers$p_val_adj <= 0.1, , drop=FALSE])
    downregulated <- rownames(markers[markers$avg_log2FC <= -1 & markers$p_val_adj <= 0.1, , drop=FALSE])
    
    DEG_results[[paste0(cell, "_Up_", tp)]]   <- upregulated
    DEG_results[[paste0(cell, "_Down_", tp)]] <- downregulated
    ### NEW: store stats for significant DEGs
    # NEW: store stats for significant DEGs
    sig_genes <- c(upregulated, downregulated)
    if (length(sig_genes) > 0) {
      sig_markers <- markers[sig_genes, c("avg_log2FC","p_val","p_val_adj"), drop=FALSE]
      sig_markers <- tibble::rownames_to_column(sig_markers, "Gene")
      sig_markers$CellType  <- cell
      sig_markers$Timepoint <- tp
      DEG_results_stats[[paste0(cell, "_", tp)]] <- sig_markers
    }
  }
}

library(tibble)
library(dplyr)

# Convert list to data frame
DEG_table <- tibble(
  Comparison = names(DEG_results),
  Genes = sapply(DEG_results, function(x) paste(x, collapse = "; "))
)

# Save as CSV for Excel
write.csv(DEG_table, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/Pancreas_CellType_DEGs.csv", row.names = FALSE)
### NEW: write combined table with FC + pvals
if (length(DEG_results_stats) > 0) {
  DEG_table_stats <- bind_rows(DEG_results_stats)
  write.csv(DEG_table_stats,
            file.path(out_dir, "WholePancreas", "Pancreas_CellType_DEGs_with_stats.csv"),
            row.names = FALSE)
}

## Plot DEGs----
library(Seurat)
library(dplyr)
library(ggplot2)

library(Seurat)
library(ggplot2)


## DEG extraction for a given cell & time ------
library(Seurat)
library(dplyr)
library(tibble)

get_markers_one <- function(obj, cell, tp) {
  # subset exactly like in your loop
  cell_time_subset <- subset(obj, subset = (CellType == cell & time == tp))
  Idents(cell_time_subset) <- cell_time_subset$group
  
  # size check (same thresholds)
  grp_counts <- table(cell_time_subset$group)
  n_prog <- as.integer(grp_counts["Progressor"]); if (is.na(n_prog)) n_prog <- 0
  n_non  <- as.integer(grp_counts["Non-Progressor"]); if (is.na(n_non)) n_non <- 0
  
  if (n_prog < 3 || n_non < 3) {
    stop(sprintf("Skipping %s | %s: Progressor=%d, Non-Progressor=%d (<3).",
                 cell, tp, n_prog, n_non), call. = FALSE)
  }
  
  # FindMarkers with identical args
  markers <- FindMarkers(
    cell_time_subset,
    ident.1 = "Progressor", ident.2 = "Non-Progressor",
    assay = "SCT", recorrect_umi = FALSE
  )
  
  # adjust p-values (BH), exactly as before
  markers$p_val_adj <- p.adjust(markers$p_val, method = "BH")
  
  # sanity check: require the same columns used in your loop
  if (!all(c("avg_log2FC","p_val","p_val_adj") %in% colnames(markers))) {
    stop("Expected columns 'avg_log2FC', 'p_val', and 'p_val_adj' not found in FindMarkers output.", call. = FALSE)
  }
  
  # identical thresholds
  up_genes   <- rownames(markers[markers$avg_log2FC >=  1 & markers$p_val_adj <= 0.1, , drop = FALSE])
  down_genes <- rownames(markers[markers$avg_log2FC <= -1 & markers$p_val_adj <= 0.1, , drop = FALSE])
  
  # build a stats table exactly like you stored in DEG_results_stats
  sig_genes <- c(up_genes, down_genes)
  sig_table <- NULL
  if (length(sig_genes) > 0) {
    sig_table <- markers[sig_genes, c("avg_log2FC","p_val","p_val_adj"), drop = FALSE]
    sig_table <- tibble::rownames_to_column(sig_table, "Gene")
    sig_table$CellType  <- cell
    sig_table$Timepoint <- tp
  }
  
  # return a structured list
  return(list(
    markers_table = markers,        # full table
    up_genes = up_genes,            # vector
    down_genes = down_genes,        # vector
    sig_table = sig_table,          # stats for sig genes (like your DEG_results_stats entries)
    n_progressor = n_prog,
    n_nonprogressor = n_non
  ))
}

## ---- Volcano for a single result table --------------------------------
# install.packages("BiocManager"); BiocManager::install("EnhancedVolcano")  # if needed
library(EnhancedVolcano)

plot_volcano_one <- function(res_one,
                             title = "",
                             genes_to_label = NULL,
                             fc_cutoff = 1,
                             padj_col = "p_val_adj") {
  
  stopifnot(is.list(res_one), !is.null(res_one$markers_table))
  df <- res_one$markers_table
  
  # Ensure required columns & numeric types
  req <- c("avg_log2FC", "p_val", padj_col)
  miss <- setdiff(req, colnames(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  df$avg_log2FC <- as.numeric(df$avg_log2FC)
  df$p_val      <- as.numeric(df$p_val)
  df[[padj_col]] <- as.numeric(df[[padj_col]])
  
  # Handle NA/zero padj so plotting works and those points remain "not significant"
  df[[padj_col]][is.na(df[[padj_col]])] <- 1
  df[[padj_col]][df[[padj_col]] <= .Machine$double.xmin] <- .Machine$double.xmin
  
  # Colors: only two “sig + big effect” classes; everything else grey
  keyvals <- ifelse(df[[padj_col]] < 0.1 & df$avg_log2FC >=  fc_cutoff, "#7E0000",
                    ifelse(df[[padj_col]] < 0.1 & df$avg_log2FC <= -fc_cutoff, "#002E66",
                           "gray70"))
  names(keyvals) <- ifelse(df[[padj_col]] < 0.1 & df$avg_log2FC >=  fc_cutoff, "Progressor-Upregulated",
                           ifelse(df[[padj_col]] < 0.1 & df$avg_log2FC <= -fc_cutoff, "Non-Progressor-Upregulated",
                                  "Not significant"))
  
  # Labels: only points that satisfy FDR<0.1 & |log2FC|≥1; if genes_to_label is
  # provided, keep only those that also meet the significance/effect criteria
  sig_big <- rownames(df)[df[[padj_col]] < 0.1 & abs(df$avg_log2FC) >= fc_cutoff]
  if (is.null(genes_to_label)) {
    sel_labels <- sig_big
  } else {
    sel_labels <- intersect(genes_to_label, sig_big)
  }
  
  EnhancedVolcano(
    df,
    lab = rownames(df),
    selectLab = sel_labels,
    x = "avg_log2FC",
    y = padj_col,                     # use adjusted p; y-axis is -log10(padj)
    title = title,
    pCutoff = 0.1,                    # horizontal line at -log10(0.1) = 1
    FCcutoff = fc_cutoff,             # vertical lines at -1 and +1
    cutoffLineType = "dashed",
    cutoffLineWidth = 0.6,
    cutoffLineCol = "grey30",
    colCustom = keyvals,
    pointSize = 2.2,
    labSize = 5.2,
    colAlpha = 0.85,
    legendPosition = "right",
    drawConnectors = TRUE,
    widthConnectors = 0.6,
    boxedLabels = TRUE,
    xlab = expression("log"[2]~"Fold Change (Progressor / Non-Progressor)"),
    ylab = expression("-log"[10]~"(adj p)"),
    #lim = c(-8, 8),
    ylim = c(0, 4),
    max.overlaps = 20,
    caption = NULL
  )
}

save_deg_csv <- function(res_one, filename) {
  stopifnot(is.list(res_one), !is.null(res_one$markers_table))
  df <- res_one$markers_table
  
  out <- data.frame(
    Gene      = rownames(df),
    log2FC    = df$avg_log2FC,
    pvalue    = df$p_val,
    padj      = df$p_val_adj,
    stringsAsFactors = FALSE
  )
  
  write.csv(out, filename, row.names = FALSE)
}



### Macrophage Early----
res_mac_early <- get_markers_one(T1D_Timepoints, cell = "Macrophage", tp = "Early")
save_deg_csv(res_mac_early,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/Macrophage_Early_DEGs.csv")
plot_volcano_one(
  res_mac_early,
  title = "Macrophage - Early",
  genes_to_label = c("Spp1","Crp","Fn1","Socs6","Zfp36l2","S100a16","Scara3","Cd34","Flt1","Fstl1","Igfbp5","Igfbp7","Tsc22d3","Tnfrsf19","Tnfrsf10b",
                     "Fos","Junb","Jdp2","Il6ra","Cd163","Cd33","Cd209f","Ltc4s","Pf4","Mertk","Siglec1","Timd4","Slc40a1")

)


### Macrophage Intermediate----
res_mac_int <- get_markers_one(T1D_Timepoints, cell = "Macrophage", tp = "Intermediate")
save_deg_csv(res_mac_int,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/Macrophage_Intermediate_DEGs.csv")
plot_volcano_one(
  res_mac_int,
  title = "Macrophage - Intermediate",
  genes_to_label = c(c("Mertk","Eepd1","Dusp1","Irf2bp2","Sirpa","Ctss","Cd300a","Ifitm2","Ifngr1","Rnd3","Picalm","Sh3pxd2b","Rgl2",
                       "Hspa1b","Ighd","Igkc","Rac2","Trac","Cd79a","Itk","Hspa1a","Lat","Gimap1",
                       "Tnfrsf1b")
)
)

### CD4 T Cell Early----
res_CD4_early <- get_markers_one(T1D_Timepoints, cell = "CD4 T Cell", tp = "Early")
save_deg_csv(res_CD4_early,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/CD4TCell_Early_DEGs.csv")
plot_volcano_one(
  res_CD4_early,
  title = "CD4 T Cell - Early",
  genes_to_label = c("Tnfaip3","Cxcr4","Cd74","H2-Eb1","H2-Ab1","H2-Aa","Dusp1","Dusp2","Zfp36l2","Ddit4",
                     "Tnfrsf9","Ifit1","Oas2","Plscr1","Alms1","Krt10","Ar","Rab2b","Extl2","Cacna2d4")
)

### CD4 T Cell Intermediate----
res_CD4_int <- get_markers_one(T1D_Timepoints, cell = "CD4 T Cell", tp = "Intermediate")
save_deg_csv(res_CD4_int,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/CD4TCell_Intermediate_DEGs.csv")
plot_volcano_one(
  res_CD4_int,
  title = "CD4 T Cell - Intermediate",
  genes_to_label = c(c("Zfp36l2","Creb3l2","Cd40lg","Ifngr1",
                       "Hspa1b","Hspa1a","Hsp90ab1","Hsph1","Hspe1","Hspa4l",
                       "Fos","Fosb","Isg15","Irf7","H2-T24","Id3")
  )
)

### CD8 T Cell Early----
res_CD8_early <- get_markers_one(T1D_Timepoints, cell = "CD8 T Cell", tp = "Early")
save_deg_csv(res_CD8_early,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/CD8TCell_Early_DEGs.csv")
plot_volcano_one(
  res_CD8_early,
  title = "CD8 T Cell - Early",
  c("Zfp36l2","Tsc22d3","Cxcr4","Tnfaip3","Cd74","H2-Eb1","H2-Ab1","H2-Aa","Dusp2","Fos",
    "Usp18","Ifit3","Trim30c","Fut7","Ptprf","Mgst3","Gstk1","Ahcy","Rbbp9")
)

## CD8 T Cell Intermediate----
res_CD8_int <- get_markers_one(T1D_Timepoints, cell = "CD8 T Cell", tp = "Intermediate")
save_deg_csv(res_CD8_int,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/CD8TCell_Intermediate_DEGs.csv")
plot_volcano_one(
  res_CD8_int,
  title = "CD8 T Cell - Intermediate",
  c("Zfp36l2","Tyrobp","Fasl","Gzmb","Tcrg-C2","Trdc","Lgals3","Rgs1","Tnfrsf4","Egr2",
    "Hspa1a","Hspa1b","Hsph1","Hspa4l","Dnajb1","Dnajb4","Fos","Fosb","Isg15","Nfkbid","Id3")
)

### B Cell Early----
res_BCell_early <- get_markers_one(T1D_Timepoints, cell = "B Cell", tp = "Early")
save_deg_csv(res_BCell_early,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/BCell_Early_DEGs.csv")
plot_volcano_one(
  res_BCell_early,
  title = "B Cell  - Early",
  c("Zfp36l2","Tsc22d3","Ddit4","Dusp2","Itk","Lef1","Ptger4","Irs2","Igfbp4",
    "Ackr2","C5ar2","Sema7a","Spic","Pstpip2","Plxna1","Klf4","Anxa2","Hmga2","Cd300lf")
)

### B Cell Intermediate----
res_BCell_int <- get_markers_one(T1D_Timepoints, cell = "B Cell", tp = "Intermediate")
save_deg_csv(res_BCell_int,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/BCell_Intermediate_DEGs.csv")
plot_volcano_one(
  res_BCell_int,
  title = "B Cell - Intermediate",
  c("Zfp36l2","Tsc22d3","Igfbp4","Rgs1","Rgs10","Cd40lg","Cd27","Cd5","Prkcq","Camk4",
    "Fos","Fosb","Hspa1a","Hspa1b","Hsph1","Hspa4l","Hspe1","Dnajb1","Lifr","Nfkbid")
)

### ILC Early----
res_ILC_early <- get_markers_one(T1D_Timepoints, cell = "ILC", tp = "Early")
save_deg_csv(res_ILC_early,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/ILC_Early_DEGs.csv")
plot_volcano_one(
  res_ILC_early,
  title = "ILC - Early",
  c("Zfp36l2","Tsc22d3","Txnip",
    "Gzmb","Vdr","Gpr15","Tead2","Arhgef26","Lclat1","Derl3","Pdia2","Cbs","Nes")
)

### ILC Intermediate----
res_ILC_int <- get_markers_one(T1D_Timepoints, cell = "ILC", tp = "Intermediate")
save_deg_csv(res_ILC_int,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/ILC_Intermediate_DEGs.csv")
plot_volcano_one(
  res_ILC_int,
  title = "ILC - Intermediate",
  c("Hspa1b")
)

### DC Early----
res_DC_early <- get_markers_one(T1D_Timepoints, cell = "Dendritic Cell", tp = "Early")
save_deg_csv(res_DC_early,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/DC_Early_DEGs.csv")

### DC Early----
res_DC_int <- get_markers_one(T1D_Timepoints, cell = "Dendritic Cell", tp = "Intermediate")
save_deg_csv(res_DC_int,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/DC_Intermediate_DEGs.csv")


### Plasma Early----
res_Plasma_early <- get_markers_one(T1D_Timepoints, cell = "Plasma Cell", tp = "Early")
save_deg_csv(res_Plasma_early,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/Plasma_Early_DEGs.csv")

### Plasma Intermediate----
res_Plasma_int <- get_markers_one(T1D_Timepoints, cell = "Plasma Cell", tp = "Intermediate")
save_deg_csv(res_Plasma_int,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/Plasma_Intermediate_DEGs.csv")

### NK Early----
res_NK_early <- get_markers_one(T1D_Timepoints, cell = "NK Cell", tp = "Early")
save_deg_csv(res_NK_early,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/NK_Early_DEGs.csv")
plot_volcano_one(
  res_NK_early,
  title = "NK - Early",

)

### NK Intermediate----
res_NK_int <- get_markers_one(T1D_Timepoints, cell = "NK Cell", tp = "Intermediate")
save_deg_csv(res_NK_int,
             "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/NK_Intermediate_DEGs.csv")
plot_volcano_one(
  res_NK_int,
  title = "NK - Intermediate",
  
)
## Pathway Analysis----


# Now DEG_results contains all up/downregulated genes for each cell type and time point
library(clusterProfiler)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)  # Use org.Hs.eg.db for human genes

# Create an empty list to store GO enrichment results
GO_results <- list()

# Loop through DEG_results
for (name in names(DEG_results)) {
  genes <- DEG_results[[name]]
  
  # Skip if there are no DEGs
  if (length(genes) == 0) next  
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)  # Remove NAs
  # Define universe (all genes in your dataset or a background gene set)

  # Perform GO enrichment
  ego_all <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = "ALL",  # Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  pool = TRUE)

  
  # Filter out CC terms
  # Ensure ego_all is still an enrichResult object (not a data.frame)
  ego <- ego_all
  
  # Filter the @result slot inside the S4 object
  ego@result <- ego_all@result[ego_all@result$ONTOLOGY %in% c("BP", "MF"), ]
  # ego <- enrichKEGG(
  #   gene= entrez_ids,
  #   organism = "mmu",
  #   keyType = "kegg",
  #   pvalueCutoff = 1,
  #   pAdjustMethod = "BH",
  #   minGSSize = 10,
  #   maxGSSize = 500,
  #   qvalueCutoff = 1,
  #  # use_internal_data = FALSE
  # )
  
  # Store results if enrichment was found
  if (!is.null(ego) && nrow(ego@result) > 0) {
    GO_results[[name]] <- ego
  }
}

library(dplyr)
library(pheatmap)
library(purrr)

# Initialize an empty list to store p-values
go_pvalues <- list()

for (name in names(GO_results)) {
  df <- GO_results[[name]]@result  # Extract the result table
  
  if (!is.null(df) && nrow(df) > 0) {
    # Ensure df is a data.frame (Base R)
    df <- as.data.frame(df)  
    
    # Explicitly use dplyr::select to avoid conflicts
    go_pvalues[[name]] <- dplyr::select(df, Description, p.adjust)  # Select columns with dplyr
  }
}

# Merge all data into one dataframe using bind_rows (alternative to reduce)
go_pval_df <- bind_rows(go_pvalues, .id = "CellType_Time")
library(tidyr)  # Make sure tidyr is loaded

go_pval_df_wide <- tidyr::pivot_wider(go_pval_df, names_from = CellType_Time, values_from = p.adjust)
go_pval_df_wide <- as.data.frame(go_pval_df_wide)   # drop tibble class
rownames(go_pval_df_wide) <- go_pval_df_wide$Description
go_pval_df_wide <- dplyr::select(go_pval_df_wide, -Description)


go_pval_matrix <- as.matrix(go_pval_df_wide)
go_pval_matrix[is.na(go_pval_matrix)] <- 1
# Filter pathways where p-value < 0.05 for at least one group
filtered_matrix <- go_pval_matrix[rowSums(go_pval_matrix <=0.1, na.rm = TRUE) >= 1, ]

# Convert p-values to -log10(p-value) for better visualization
filtered_matrix[!is.na(filtered_matrix)] <- -log10(filtered_matrix[!is.na(filtered_matrix)])
# filtered_matrix is already the matrix with -log10(p.adjust)

# Convert matrix to data frame with rownames preserved
filtered_df <- as.data.frame(filtered_matrix)
filtered_df <- tibble::rownames_to_column(filtered_df, var = "Pathway")

# Write to CSV
write.csv(filtered_df,
          file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/WholePancreas/Filtered_enrichedpathways.csv",
          row.names = FALSE)


# Initialize an empty vector to collect pathways
selected_pathways <- c()
# 
# # Loop through each column
# for (col in colnames(filtered_matrix)) {
#   
#   # For each column, get the top 1- pathways with highest -log10(p)
#   top10 <- rownames(filtered_matrix)[order(filtered_matrix[, col], decreasing = TRUE)][1:20]
#   
#   # Add to the selected list
#   selected_pathways <- c(selected_pathways, top10)
# }
#Select biologucally revelant pathways with FDR<=0.1
selected_pathways<-c(#Macrophage-Up Early
                      "interleukin-6 production",
                     "regulation of interleukin-6 production",
                     "positive regulation of inflammatory response",
                     "regulation of inflammatory response",
                     "chemokine-mediated signaling pathway",
                     "cytokine-mediated signaling pathway",
                     "positive regulation of pattern recognition receptor signaling pathway",
                     "alpha-beta T cell activation",
                     "regulation of T-helper 17 type immune response",
                     "complement activation",
                     #Macrophage-Down Early
                     "positive regulation of chemotaxis",
                     "cellular response to transforming growth factor beta stimulus",
                     "humoral immune response",
                     #Macrophage-Up Intermediate
                     "antigen processing and presentation of exogenous peptide antigen via MHC class II",
                     "positive regulation of T cell activation",
                     "toll-like receptor signaling pathway",
                     "canonical NF-kappaB signal transduction",
                     "tumor necrosis factor production",
                     "interleukin-1 beta production",
                     "cytokine production involved in immune response",
                     "macrophage activation",
                     "cellular response to lipopolysaccharide",
                     "type I interferon production",
                     #Macrophage-Down Intermediate
                     "antigen processing and presentation",
                     "positive regulation of T cell activation",
                     "immune response-activating signaling pathway",
                     "canonical NF-kappaB signal transduction",
                     "toll-like receptor signaling pathway",
                     "tumor necrosis factor production",
                     "interleukin-1 beta production",
                     "cellular response to lipopolysaccharide",
                     "type I interferon production",
                     "macrophage activation",
                     #CD8 T Cells- Up Early
                     "positive regulation of cell killing",
                     "positive regulation of leukocyte mediated cytotoxicity",
                     "positive regulation of natural killer cell activation",
                     "regulation of natural killer cell activation",
                     "immune response-activating signaling pathway",
                     "toll-like receptor signaling pathway",
                     "canonical NF-kappaB signal transduction",
                     "tumor necrosis factor production",
                     "interleukin-1 beta production",
                     "macrophage activation involved in immune response",
                     #CD8 T Cells- Down Early
                     "positive regulation of lymphocyte activation",
                     "positive regulation of B cell activation",
                     "positive regulation of T cell migration",
                     "leukocyte activation involved in immune response",
                     "cell activation involved in immune response",
                     "response to lipopolysaccharide",
                     "chronic inflammatory response",
                     "positive regulation of neutrophil migration",
                     "acute inflammatory response to antigenic stimulus",
                     "T cell mediated immunity",
                     #CD8 T Cells- Up Intermediate
                     "antigen processing and presentation of exogenous peptide antigen via MHC class II",
                     "antigen processing and presentation of peptide antigen via MHC class II",
                     "positive regulation of lymphocyte activation",
                     "positive regulation of leukocyte activation",
                     "B cell activation",
                     "lymphocyte proliferation",
                     "mononuclear cell proliferation",
                     "leukocyte proliferation",
                     "regulation of T cell activation",
                     "macrophage activation involved in immune response",
                     "neutrophil activation involved in immune response",
                     "positive regulation of immune effector process",
                     "T cell proliferation",
                     "positive regulation of natural killer cell activation",
                     "positive regulation of interleukin-8 production",
                     "leukocyte activation involved in inflammatory response",
                     "cell activation involved in immune response",
                     "macrophage chemotaxis",
                     "cell killing",
                     "leukocyte mediated cytotoxicity",
                     "immunological synapse formation",
                     #CD8 T Cells- Down Intermediate
                     "type I interferon production",
                     "regulation of type I interferon production",
                     "interferon-beta production",
                     "regulation of interferon-beta production",
                     "negative regulation of NF-kappaB transcription factor activity",
                     "NF-kappaB binding",
                     "double-stranded RNA binding",
                     "cellular response to unfolded protein",
                     "Hsp70 protein binding",
                     #CD4 T Cells- Up Early
                     "antigen processing and presentation of exogenous peptide antigen via MHC class II",
                   "antigen processing and presentation of peptide or polysaccharide antigen via MHC class II",
                     "antigen processing and presentation of exogenous antigen",
                     "positive regulation of lymphocyte activation",
                     "positive regulation of leukocyte activation",
                     "positive regulation of cell activation",
                     "B cell activation",
                     "positive regulation of T cell activation",
                     "immune receptor activity",
                     #CD4 T Cells- Down Early
                   "CD4 receptor binding",
                   "T cell differentiation in thymus",
                   "mast cell activation",
                   "acute inflammatory response",
                   "positive regulation of tumor necrosis factor production",
                   "positive regulation of tumor necrosis factor superfamily cytokine production",
                   "positive regulation of NF-kappaB transcription factor activity",
                   "type I interferon-mediated signaling pathway",
                   #CD4 T Cells- Up Intermediate
                   "B cell activation",
                   "regulation of B cell activation",
                   "B cell differentiation",
                   "regulation of B cell differentiation",
                   "regulation of hemopoiesis",
                   "regulation of lymphocyte differentiation",
                   "T cell differentiation",
                   "T cell differentiation in thymus",
                   #CD4 T Cells- Down Intermediate
                   "type I interferon production",
                   "positive regulation of type I interferon production",
                   "type I interferon-mediated signaling pathway",
                   "cellular response to type I interferon",
                   "interferon-mediated signaling pathway",
                   "defense response to virus",
                   #B Cell-Up Early
                   "B cell activation",
                   "regulation of B cell activation",
                   "B cell differentiation",
                   "lymphocyte proliferation",
                   "positive regulation of lymphocyte activation",
                   "interleukin-6 production",
                   #B Cell- Down Early
                   "negative regulation of chemotaxis",
                   "S100 protein binding",            
                   #B Cell-Up Intermediate
                   "positive regulation of lymphocyte activation",
                   "B cell activation",
                   "T cell differentiation",
                   "positive regulation of T cell proliferation",
                   "positive regulation of B cell activation",
                   "positive regulation of CD4-positive, alpha-beta T cell activation",
                   "T cell receptor signaling pathway",
                   "T cell costimulation",
                   "lymphocyte costimulation",
                   "positive regulation of interleukin-2 production",
                   "positive regulation of adaptive immune response",
                   "immunoglobulin mediated immune response",
                   "B cell mediated immunity",
                   "positive regulation of tolerance induction",
                   "dendritic cell differentiation",
                   #B Cell-Down Intermediate- 
                   "cytokine-mediated signaling pathway",
                   "response to reactive oxygen species",
                   #NK Cells- Up Early
                   "regulation of leukocyte apoptotic process",
                   "lymphocyte homeostasis",
                   "negative regulation of lymphocyte apoptotic process",
                   "T cell apoptotic process",
                   "cellular response to tumor necrosis factor",
                   "response to transforming growth factor beta",
                  "activation-induced cell death of T cells",
                   "T cell homeostasis",
                   #NK Cells- Down Late
                  "lipopolysaccharide-mediated signaling pathway",
                  "cellular response to lipopolysaccharide",
                  "p38MAPK cascade",
                  "positive regulation of JUN kinase activity",
                  #NK Cells- Up Intermediate
                  #"lipid transport",
                  #"cholesterol transport",
                  #"sterol transport",
                  #"lipid catabolic process",
                  #NK Cells- Down Intermediate
                  "NF-kappaB binding",
                  "double-stranded RNA binding",
                  "positive regulation of proteasomal ubiquitin-dependent protein catabolic process",
                  "regulation of ubiquitin-dependent protein catabolic process",
                  "positive regulation of proteolysis",
                  "regulation of cysteine-type endopeptidase activity involved in apoptotic process",
                  #ILC- Up Early
                  "B cell activation",
                  "T cell differentiation",
                  "regulation of B cell activation",
                  "lymphocyte differentiation",
                  "T cell differentiation in thymus",
                  "cellular response to tumor necrosis factor",
                  "response to tumor necrosis factor",
                  "cellular response to transforming growth factor beta stimulus",
                  "response to granulocyte macrophage colony-stimulating factor",
                  "activation-induced cell death of T cells",
                  #ILC- Down Early
                  "homotypic cell-cell adhesion",
                  #ILC- Up Intermediate- Nothing
                  #ILC-Down Intermediate
                  "regulation of cysteine-type endopeptidase activity involved in apoptotic process",
                  "positive regulation of protein catabolic process",
                  "chaperone-mediated protein folding",
                  "unfolded protein binding",
                  #"heat shock protein binding",
                  "NF-kappaB binding"
                     )
# Get unique pathways (since some pathways may be selected by multiple columns)
selected_pathways <- unique(selected_pathways)



# Now subset the filtered_matrix to keep only the selected pathways
final_matrix_to_plot <- filtered_matrix[selected_pathways, ]
# Better, aesthetic color palette for heatmaps (white -> light blue -> deep blue)
#new_colors <- colorRampPalette(c("white", "skyblue", "royalblue", "navy"))(1000)
# Soft pink -> medium red -> deep dark red
# Light yellow -> orange -> red -> dark red
new_colors <- colorRampPalette(c("white","lightyellow","lemonchiffon", "lightgoldenrodyellow", ,"yellow", "orange", "red", "darkred"))(10000)

# Cap the values at 10 for better visualization
mat_capped <- pmin(final_matrix_to_plot, 8)

pheatmap(mat_capped,
         cluster_rows = TRUE,
         cluster_cols = FALSE,   # no column clustering
         color = new_colors,     # <<< Use the better color palette
         main = "",
         fontsize_row = 18,
         fontsize_col = 18,
         border_color = NA)      # <<< no gridlines, cleaner for publication



