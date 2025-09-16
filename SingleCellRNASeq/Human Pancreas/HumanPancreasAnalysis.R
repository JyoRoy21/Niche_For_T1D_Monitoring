## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/HumanSingleCellRNASeq/")
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
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
HumanPancreas<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/HumanSingleCellRNASeq/T1D_T2D_20220428.rds")
unique(Idents(HumanPancreas))
#Idents(Tcells)<-Tcells$TCellSubType
HumanPancreasImmune<-subset(HumanPancreas, idents='Immune')
rm(HumanPancreas)
HumanPancreasImmune


ElbowPlot(HumanPancreasImmune, ndims = 50)
HumanPancreasImmune <- RunUMAP(HumanPancreasImmune, dims = 1:30)
HumanPancreasImmune <- FindNeighbors(HumanPancreasImmune, dims = 1:30, verbose = F)
HumanPancreasImmune <- FindClusters(HumanPancreasImmune, res = 0.5)#graph.name = "SCT_nn"


DimPlot(HumanPancreasImmune, reduction = "umap", combine = F, label = T)
unique(HumanPancreasImmune$sample_id)


HumanPancreasImmune <- PrepSCTFindMarkers(HumanPancreasImmune)
all.markers <- FindAllMarkers(HumanPancreasImmune , only.pos = TRUE,recorrect_umi = FALSE)

top20_sct=all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)

# Save the top 20 markers grouped by cluster as a CSV file
#write.csv(top20_sct, file = "top20_sct_Tcells_T1D_Timepoints.csv", row.names = FALSE)


DotPlot(HumanPancreasImmune, assay = "SCT", features = c(
  "CD3D", "CD3E", "CD3G",  # Pan-T Cell
  "CD4", "IL7R", "CCR7", "FOXP3",  # CD4
  "CD8A", "CD8B", "GZMK", "GZMB",  # CD8
  "CD19", "MS4A1", "CD79A", "CD79B", "IGKC",  # B Cells
  "CD68", "CD14", "C1QA", "C1QB", "C1QC", "CSF1R", "MRC1",  # Macrophage
  "IL1B", "NOS2", "TNF", "HLA-DRA", "CD80","CD86",  # M1-like Macrophages
  "CD163", "MSR1", "IL10", "TGFB1",  # M2-like Macrophages
  "CD1C", "FLT3", "CLEC10A",  # cDC (Removed duplicate "HLA-DRA")
  "L3RA", "TCF4", "LILRA4", "NRP1",  # pDC
  "NCAM1", "KLRD1", "KLRF1", "NKG7", "FCGR3A",  # NK
  "FCGR3B", "S100A8", "S100A9", "MPO", "ELANE",  # Neutrophils
  "TMOD4", "ADAMTS5", "PDGFRL","ACTA2", #Mesenchymal
  "INS","GCG","SST","PTPRC"
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

# Define the T cell subtypes for each cluster
ImmuneSubtype <- c("Macrophage",  # 0
                  "Macrophage",  # 1
                  "Monocytes",  # 2
                  "Macrophage",  # 3
                  "Non-Immune",  # 4
                  "Non-Immune",  # 5
                  "Macrophage", # 6
                  "CD8 T Cell"
)        

# Ensure the Idents are in numeric format
ident_integers <- as.integer(Idents(HumanPancreasImmune))

# Check if there are any unexpected levels
unique(ident_integers)

# Add the T cell subtype information to the metadata of the Seurat object
HumanPancreasImmune$CellType <- ImmuneSubtype[as.integer(Idents(HumanPancreasImmune))]

# Check the updated metadata
head(HumanPancreasImmune@meta.data)
# Subset the object by removing Non-Immune cells
HumanPancreasImmune <- subset(HumanPancreasImmune, subset = CellType != "Non-Immune")

# Plot the UMAP with the updated object
DimPlot(HumanPancreasImmune, reduction = 'umap', group.by = "CellType", label = TRUE)
saveRDS(HumanPancreasImmune,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/HumanSingleCellRNASeq/HumanPancreasImmuneCells.rds")

unique(HumanPancreasImmune$disease_state)

## Macrophage ----

#### DEG Analysis---- 
Idents(HumanPancreasImmune)<-HumanPancreasImmune$CellType
Macrophage_Human_T1D_AAb <- subset(HumanPancreasImmune, idents = "Macrophage", subset = disease_state %in% c("T1D","AAB"))
Macrophage_Human_T1D_Control <- subset(HumanPancreasImmune, idents = "Macrophage", subset = disease_state %in% c("T1D","Control"))
Macrophage_Human_AAb_Control <- subset(HumanPancreasImmune, idents = "Macrophage", subset = disease_state %in% c("AAB","Control"))

table(Macrophage_Human_T1D_AAb$disease_state)
table(Macrophage_Human_T1D_Control$disease_state)
table(Macrophage_Human_AAb_Control$disease_state)

Idents(Macrophage_Human_T1D_AAb) <- Macrophage_Human_T1D_AAb$disease_state
Idents(Macrophage_Human_T1D_Control) <- Macrophage_Human_T1D_Control$disease_state
Idents(Macrophage_Human_AAb_Control) <- Macrophage_Human_AAb_Control$disease_state


Macrophage_Human_T1D_AAb = PrepSCTFindMarkers(Macrophage_Human_T1D_AAb)
Macrophage_Human_T1D_Control = PrepSCTFindMarkers(Macrophage_Human_T1D_Control)
Macrophage_Human_AAb_Control<-PrepSCTFindMarkers(Macrophage_Human_AAb_Control)

##### T1D Vs AAB----

Macrophage_Human_T1DVsAAb = FindMarkers(Macrophage_Human_T1D_AAb, 
                            ident.1 = "T1D", 
                            ident.2 = "AAB", 
                            assay = "SCT",
                            recorrect_umi = FALSE)
Macrophage_Human_T1DVsAAb$p_val_adj <- p.adjust(Macrophage_Human_T1DVsAAb$p_val, method = "BH")

##### T1D Vs Control----
Macrophage_Human_T1DVsControl = FindMarkers(Macrophage_Human_T1D_Control,
                             ident.1 = "T1D",
                             ident.2 = "Control",
                             assay= "SCT",
                             recorrect_umi = FALSE)
Macrophage_Human_T1DVsControl$p_val_adj <- p.adjust(Macrophage_Human_T1DVsControl$p_val, method = "BH")

colnames(Macrophage_Human_T1DVsControl)


##### AAB Vs Control----
Macrophage_Human_AAbVsControl = FindMarkers(Macrophage_Human_AAb_Control,
                                            ident.1 = "AAB",
                                            ident.2 = "Control",
                                            assay= "SCT",
                                            recorrect_umi = FALSE)
Macrophage_Human_AAbVsControl$p_val_adj <- p.adjust(Macrophage_Human_AAbVsControl$p_val, method = "BH")



#### Volcano Plot ----

##### T1D Vs AAb----
Macrophage_Human_T1DVsAAb_data = Macrophage_Human_T1DVsAAb

keyvals <- ifelse(
  (Macrophage_Human_T1DVsAAb_data$avg_log2FC < -1 & Macrophage_Human_T1DVsAAb_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Macrophage_Human_T1DVsAAb_data$avg_log2FC > 1 & Macrophage_Human_T1DVsAAb_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non-Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated-T1D'
names(keyvals)[keyvals == '#1465AC'] <- 'Upregulated-AAB'


#png(file = "AllergyScaf_NPPBSCh4_VolcanoPlot.png",units = "in",width = 10, height = 10, res = 400)
EnhancedVolcano(Macrophage_Human_T1DVsAAb_data,
                lab = rownames(Macrophage_Human_T1DVsAAb_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                ylab = bquote(~-Log[10] ~ italic(P[adj])),
                pCutoff = 0.05,
                FCcutoff = 1,
                colCustom = keyvals,
                selectLab = rownames(Macrophage_Human_T1DVsAAb_data)[which(names(keyvals) %in% c('Upregulated-T1D', 'Upregulated-AAB'))],
                title = "Macrophage: T1D vs AAB",
                colAlpha = 0.75,,
                max.overlaps = 10,
                drawConnectors = T,
                widthConnectors = 0.75,
                ylim = c(0,15),
                xlim=c(-8,8),
                boxedLabels = F)
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/HumanSingleCellRNASeq/Volcanoplot_Macrophage_T1DvsAAB.png", width = 10, height = 8, dpi = 600)

#dev.off()
write.csv(Macrophage_Human_T1DVsAAb_data,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/HumanSingleCellRNASeq/Macrophage_DEG_T1DVsAAB.csv",row.names = TRUE)

##### T1D Vs Control----
Macrophage_Human_T1DVsControl_data = Macrophage_Human_T1DVsControl

keyvals <- ifelse(
  (Macrophage_Human_T1DVsControl_data$avg_log2FC < -1 & Macrophage_Human_T1DVsControl_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Macrophage_Human_T1DVsControl_data$avg_log2FC > 1 & Macrophage_Human_T1DVsControl_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non-Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated-T1D'
names(keyvals)[keyvals == '#1465AC'] <- 'Upregulated-Control'


#png(file = "AllergyScaf_NPPBSCh4_VolcanoPlot.png",units = "in",width = 10, height = 10, res = 400)
EnhancedVolcano(Macrophage_Human_T1DVsControl_data,
                lab = rownames(Macrophage_Human_T1DVsControl_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                ylab = bquote(~-Log[10] ~ italic(P[adj])),
                pCutoff = 0.05,
                FCcutoff = 1,
                colCustom = keyvals,
                selectLab = rownames(Macrophage_Human_T1DVsControl_data)[which(names(keyvals) %in% c('Upregulated-T1D', 'Upregulated-AAB'))],
                title = "Macrophage: T1D vs Control",
                colAlpha = 0.75,,
                max.overlaps = 10,
                drawConnectors = T,
                widthConnectors = 0.75,
                ylim = c(0,30),
                xlim=c(-8,8),
                boxedLabels = F)
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/HumanSingleCellRNASeq/Volcanoplot_Macrophage_T1DvsControl.png", width = 10, height = 8, dpi = 600)

write.csv(Macrophage_Human_T1DVsControl_data,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/HumanSingleCellRNASeq/Macrophage_DEG_T1DVsControl.csv",row.names = TRUE)


##### AAB Vs Control----
Macrophage_Human_AAbVsControl_data = Macrophage_Human_AAbVsControl

keyvals <- ifelse(
  (Macrophage_Human_AAbVsControl_data$avg_log2FC < -1 & Macrophage_Human_AAbVsControl_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Macrophage_Human_AAbVsControl_data$avg_log2FC > 1 & Macrophage_Human_AAbVsControl_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non-Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated-AAB'
names(keyvals)[keyvals == '#1465AC'] <- 'Upregulated-Control'


#png(file = "AllergyScaf_NPPBSCh4_VolcanoPlot.png",units = "in",width = 10, height = 10, res = 400)
EnhancedVolcano(Macrophage_Human_AAbVsControl_data,
                lab = rownames(Macrophage_Human_AAbVsControl_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                ylab = bquote(~-Log[10] ~ italic(P[adj])),
                pCutoff = 0.05,
                FCcutoff = 1,
                colCustom = keyvals,
                selectLab = rownames(Macrophage_Human_AAbVsControl_data)[which(names(keyvals) %in% c('Upregulated-T1D', 'Upregulated-AAB'))],
                title = "Macrophage: AAB vs Control",
                colAlpha = 0.75,,
                max.overlaps = 10,
                drawConnectors = T,
                widthConnectors = 0.75,
                ylim = c(0,10),
                boxedLabels = F)
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/HumanSingleCellRNASeq/Volcanoplot_Macrophage_AABvsControl.png", width = 10, height = 8, dpi = 600)

#dev.off()
write.csv(Macrophage_Human_AAbVsControl_data,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/HumanSingleCellRNASeq/Macrophage_DEG_AABVsControl.csv",row.names = TRUE)


## CD8 T Cells ----

#### DEG Analysis---- 
Idents(HumanPancreasImmune)<-HumanPancreasImmune$CellType
Macrophage_Human_T1D_AAb <- subset(HumanPancreasImmune, idents = "Macrophage", subset = disease_state %in% c("T1D","AAB"))
