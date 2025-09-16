## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
BiocManager::install("dittoSeq")
BiocManager::install("escape")
suppressPackageStartupMessages(library(escape))
library(escape)

# List all functions in the escape package
ls("package:escape")
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))
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
T1D_Timepoints<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v4.rds")
DimPlot(T1D_Timepoints)
unique(T1D_Timepoints$orig.ident)

# Convert to character to allow new labels
T1D_Timepoints@meta.data$CellType <- as.character(T1D_Timepoints@meta.data$CellSubType)

# Modify labels
T1D_Timepoints@meta.data$CellType[T1D_Timepoints@meta.data$CellSubType %in% c("CD8 memory", "CD8 exhausted effector-like")] <- "CD8 T Cell"

T1D_Timepoints@meta.data$CellType[T1D_Timepoints@meta.data$CellSubType %in% c("Tcon memory", "Tcon exhausted effector-like", "Tregs", "Tcon activated ", "Tcon Interferon Sensing")] <- "CD4 T Cell"

T1D_Timepoints@meta.data$CellType[T1D_Timepoints@meta.data$CellSubType %in% c("ILC2", "ILC3")] <- "ILC"

# Convert back to factor if needed
T1D_Timepoints@meta.data$CellType <- factor(T1D_Timepoints@meta.data$CellType)

unique(T1D_Timepoints$CellType)

# BiocManager is needed for pseudobulk
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("muscat")
remotes::install_github("immunogenomics/pseudobulk")

library(Seurat)
library(SingleCellExperiment)
library(muscat)
library(GSVA)


library(org.Mm.eg.db)
library(GO.db)
library(AnnotationDbi)

get_GO_genesets_symbols <- function(ontology = c("BP", "MF"), organism_db = org.Mm.eg.db) {
  ontology <- match.arg(ontology)

  # Get gene-to-GO mappings
  gene2go <- AnnotationDbi::select(
    organism_db,
    keys = keys(organism_db, keytype = "ENTREZID"),
    columns = c("GO", "ONTOLOGY", "SYMBOL"),
    keytype = "ENTREZID"
  )

  # Filter for the desired ontology
  gene2go <- subset(gene2go, ONTOLOGY == ontology & !is.na(SYMBOL))

  # Create GOID -> list of gene SYMBOLs
  go_list <- split(gene2go$SYMBOL, gene2go$GO)

  # Remove small gene sets
  go_list <- go_list[lengths(go_list) > 5]

  # Add GO term names for clarity
  go_terms <- Term(GOTERM[names(go_list)])
  names(go_list) <- paste(names(go_list), go_terms, sep = "_")

  return(go_list)
}


# Example usage
GS.BP.symbol <- get_GO_genesets_symbols("BP")  # Biological Process
GS.MF.symbol <- get_GO_genesets_symbols("MF")  # Molecular Function



#### CD8 T Cells----
cd8_seurat <- subset(T1D_Timepoints, subset = CellType == "CD8 T Cell")
sce_cd8 <- as.SingleCellExperiment(cd8_seurat)

# Aggregate logcounts by sample
pb <- muscat::aggregateData(
  sce_cd8,
  assay = "logcounts",
  by = "sample",
  fun = "mean"
)

# Extract the pseudobulk matrix
pseudobulk_matrix <- assay(pb)


# Assuming GS.MF.symbol is your list of gene sets
# Create parameter object using GSVAParams
params <- gsvaParam(
  as.matrix(pseudobulk_matrix),
   GS.MF.symbol
)
gsva_res <- gsva(params)

library(pheatmap)

top_var <- apply(gsva_res, 1, sd)
top_pathways <- names(sort(top_var, decreasing = TRUE))[1:25]

pheatmap::pheatmap(
  gsva_res[top_pathways, ],
  scale = "row",
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  main = "CD8 T Cell Pathway Enrichment (Pseudobulk + GSVA)"
)

# Define your gene list (Mouse gene symbols, adjust if needed)
cd8_cytotoxic_genes <- c("Gzma", "Gzmb", "Prf1", "Ifng", "Fasl", "Cd8a", "Cd8b1")
tnf_alpha_genes <- c("Tnf", "Tnfrsf1a", "Tnfrsf1b")
nfkb_pathway_genes <- c("Nfkb1", "Nfkb2", "Rel", "Rela", "Chuk", "Ikbkb", "Ikbkg")

# Combine gene sets
selected_genes <- unique(c(cd8_cytotoxic_genes, tnf_alpha_genes, nfkb_pathway_genes))

# Check which genes are in the pseudobulk matrix
genes_present <- selected_genes[selected_genes %in% rownames(pseudobulk_matrix)]

# Subset expression matrix
expr_subset <- pseudobulk_matrix[genes_present, ]
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Plot the heatmap
pheatmap::pheatmap(
  as.matrix(expr_subset),
  scale = "row",
  cluster_cols = F,
  cluster_rows = F,
  main = "Pancreas CD8 T Cell Enrichment",
  color = custom_colors
)

### Macrophages----
unique(T1D_Timepoints$CellType)
Mac_seurat <- subset(T1D_Timepoints, subset = CellType == "Macrophage")
sce_Mac <- as.SingleCellExperiment(Mac_seurat)

# Aggregate logcounts by sample
pb_Mac <- muscat::aggregateData(
  sce_Mac,
  assay = "logcounts",
  by = "sample",
  fun = "mean"
)

# Extract the pseudobulk matrix
pseudobulk_matrix_Mac <- assay(pb_Mac)

tnf_alpha_genes <- c("Tnf", "Tnfrsf1a", "Tnfrsf1b")
nfkb_pathway_genes <- c("Nfkb1", "Nfkb2", "Rel", "Rela", "Chuk", "Ikbkb", "Ikbkg")

# Combine gene sets
selected_genes <- unique(c(tnf_alpha_genes, nfkb_pathway_genes))

# Check which genes are in the pseudobulk matrix
genes_present <- selected_genes[selected_genes %in% rownames(pseudobulk_matrix_Mac)]

# Subset expression matrix
expr_subset_Mac <- pseudobulk_matrix_Mac[genes_present, ]
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Plot the heatmap
pheatmap::pheatmap(
  as.matrix(expr_subset_Mac),
  scale = "row",
  cluster_cols = F,
  cluster_rows = F,
  main = "Pancreas Macrophage Enrichment",
  color = custom_colors
)

# sce <- as.SingleCellExperiment(T1D_Timepoints)
# 
# # Load required libraries

# # Preview
# names(GS.BP.symbol)[1:3]
# GS.BP.symbol[[1]][1:5]
# 
# ES.seurat <- runEscape(
#   sce,                      # Seurat object
#   gene.sets = GS.MF.symbol,           # List of gene sets with gene symbols
#   min.size = 5,                       # Optional: gene set size filter
#   ncores = 10                         # Parallel threads
# )
# 
# 
