setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/CellChat")
library(Seurat)
library(SoupX)
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
library(ComplexHeatmap)
library(circlize)

## ----load Object and subset objects-----
T1D_Timepoints<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v6.rds")
Idents(T1D_Timepoints)<-T1D_Timepoints$CellSubType
DimPlot(T1D_Timepoints)
W6_T1D <- subset(T1D_Timepoints, subset = time == "Week6")

table(T1D_Timepoints$sample,T1D_Timepoints$group)
T1DGene_top100<- read.csv("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered_PLSDATop100.csv", sep=",", header=T,check.names = FALSE)
T1DGene_top100<-as.data.frame(T1DGene_top100)

Top100Genes<-T1DGene_top100[,1]

# Average expression per CellSubType
avg_exp <- AverageExpression(W6_T1D, features = Top100Genes, return.seurat = FALSE)$RNA
W6_T1D$group
# Scale the data per gene (row-wise z-score)
zscore_matrix <- t(scale(t(avg_exp)))

# Cap extreme values for visualization clarity (optional but helps)
zscore_matrix[zscore_matrix > 2] <- 2
zscore_matrix[zscore_matrix < -2] <- -2

# Choose a color palette


custom_colors <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# Remove rows (genes) with any NA, NaN, or Inf values
zscore_matrix_clean <- zscore_matrix[complete.cases(zscore_matrix) & 
                                       apply(zscore_matrix, 1, function(x) all(is.finite(x))), ]

# Check how many genes remain
cat("Remaining genes after cleaning:", nrow(zscore_matrix_clean), "\n")

# Plot high-quality heatmap
pheatmap(as.matrix(zscore_matrix_clean),
         color = custom_colors,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 7,
         fontsize_col = 10,
         main = "Scaffold Gene Signature Expresssion- Pancreas Microenvironment",
         border_color = NA)



