# 1. Load Libraries ----
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

# 2. Load Single Cell Object ----
getwd()
T1D_Timepoints = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v4.rds")
unique(T1D_Timepoints$CellSubType)
Idents(T1D_Timepoints)<-T1D_Timepoints$CellSubType
week12_data <- subset(T1D_Timepoints,subset =time  %in% c("Week12"))
week12_data$sample_id <- ifelse(week12_data$group == "Non-Progressor",
                               paste0("NonProgressor", gsub("Week12_", "", week6_data$sample)),
                               paste0("Progressor", gsub("Week12_", "", week6_data$sample)))

# Extract raw counts and metadata to create SingleCellExperiment object
counts_w12 <- week12_data@assays$RNA$counts 
metadata_w12 <- week12_data@meta.data
# Set up metadata as desired for aggregation and DE analysis
metadata_w12$cluster_id <- factor(week12_data@active.ident)

# Create single cell experiment object
sce_w12 <- SingleCellExperiment(assays = list(counts = counts_w12), 
                               colData = metadata_w12)

dim(colData(sce_w12))
head(colData(sce_w12))

# 3. Expression of Intermediate Scaffold Genes in CellTypes ----
ScafGenes <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/ScaffoldDESEQ/DEGAnalysis_Volcano_Intermediate.csv", sep=",", header=T) # Metadata file
rownames(ScafGenes)<-ScafGenes$X
ScafGenes$X<-NULL
colnames(ScafGenes)
# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 2, ])
# Ensure selected genes exist in the dataset
selected_genes <- selected_genes[selected_genes %in% rownames(week12_data)]
library(Seurat)
# Generate the heatmap
DoHeatmap(
  week12_data,
  features = selected_genes,
  #group.by = "cell_type",  # Replace with your cell type metadata column
  slot = "scale.data"  # Uses scaled expression data
) + scale_fill_gradientn(colors = c("blue", "white", "red"))


dotplot <- DotPlot(
  week12_data,
  features = selected_genes,
  #group.by = "cell_type",
  dot.scale = 8
) + RotatedAxis()  # Rotate x-axis labels for better visibility

# Adjust color gradient manually
dotplot + scale_color_gradientn(colors = c("blue", "white", "red"))+  
  scale_size_continuous(range = c(1, 8))  +  # Min size = 1, Max size = 8
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Increase x-axis text size & rotate
    axis.text.y = element_text(size = 16),  # Increase y-axis text size
    axis.title.x = element_text(size = 18, face = "bold"),  # Increase x-axis title size
    axis.title.y = element_text(size = 18, face = "bold")   # Increase y-axis title size
  )+
  labs(
    x = "Immunological Niche Based Differentially Expressed Genes-Intermediate Stage",  # Custom x-axis label
    y = "Pancreas Cell Types"  # Custom y-axis label
  )

# 5. DE Genes in Cell Types ----

rm(week6_data)
## I. Macrophage ----
Macrophage_Week12 <- subset(week12_data, idents = "Macrophage")
Idents(Macrophage_Week12) <- Macrophage_Week12$group
Macrophage_Week12 = PrepSCTFindMarkers(Macrophage_Week12)

Macrophage_Week12_PvNP = FindMarkers(Macrophage_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                    recorrect_umi = FALSE)

colnames(ScafGenes)
# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])
ScafDEGenes<-ScafGenes[selected_genes,]
# Filter for DE genes in the Macrophage_Week12_PvNP dataset
Macrophage_Week12_PvNP_filtered <- Macrophage_Week12_PvNP[
  (Macrophage_Week12_PvNP$p_val <= 0.05) & 
    (abs(Macrophage_Week12_PvNP$avg_log2FC) >= 1), 
]

# Ensure that the trends (direction) of fold change are the same (both ≥1 or both ≤-1)
same_trend_genes <- intersect(
  rownames(ScafDEGenes[ScafDEGenes$log2FoldChange >= 1,]),
  rownames(Macrophage_Week12_PvNP_filtered[Macrophage_Week6_PvNP_filtered$avg_log2FC >= 1,])
)

# Combine genes with the same trend for both datasets (log2FoldChange <= -1)
same_trend_genes <- union(
  same_trend_genes,
  intersect(
    rownames(ScafDEGenes[ScafDEGenes$log2FoldChange <= -1,]),
    rownames(Macrophage_Week12_PvNP_filtered[Macrophage_Week6_PvNP_filtered$avg_log2FC <= -1,])
  )
)
# Display the common DE genes
same_trend_genes
VlnPlot(Macrophage_Week12, features = c("Hpgd","Asb2","Cox6a2","Fasn"))

options(future.globals.maxSize = 5000 * 1024^2)  # Increase to 5GB

## II. CD8 exhausted effector-like ----
CD8Exhaus_Week12 <- subset(week12_data, idents = "CD8 exhausted effector-like")
Idents(CD8Exhaus_Week12) <- CD8Exhaus_Week12$group
CD8Exhaus_Week12 = PrepSCTFindMarkers(CD8Exhaus_Week12)

CD8Exhaus_Week12_PvNP = FindMarkers(CD8Exhaus_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                     recorrect_umi = FALSE)

colnames(ScafGenes)
# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])
ScafDEGenes<-ScafGenes[selected_genes,]
# Filter for DE genes in the Macrophage_Week12_PvNP dataset
Macrophage_Week12_PvNP_filtered <- Macrophage_Week12_PvNP[
  (Macrophage_Week12_PvNP$p_val <= 0.05) & 
    (abs(Macrophage_Week12_PvNP$avg_log2FC) >= 1), 
]

# Ensure that the trends (direction) of fold change are the same (both ≥1 or both ≤-1)
same_trend_genes <- intersect(
  rownames(ScafDEGenes[ScafDEGenes$log2FoldChange >= 1,]),
  rownames(Macrophage_Week12_PvNP_filtered[Macrophage_Week6_PvNP_filtered$avg_log2FC >= 1,])
)

# Combine genes with the same trend for both datasets (log2FoldChange <= -1)
same_trend_genes <- union(
  same_trend_genes,
  intersect(
    rownames(ScafDEGenes[ScafDEGenes$log2FoldChange <= -1,]),
    rownames(Macrophage_Week12_PvNP_filtered[Macrophage_Week6_PvNP_filtered$avg_log2FC <= -1,])
  )
)
# Display the common DE genes
same_trend_genes
VlnPlot(Macrophage_Week12, features = c("Hpgd","Asb2","Cox6a2","Fasn"))

options(future.globals.maxSize = 5000 * 1024^2)  # Increase to 5GB
