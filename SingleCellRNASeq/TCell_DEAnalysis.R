## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
library(Seurat)
packageVersion("Seurat")
library(parallel)
library(purrr)
library(tibble)
library(presto)
library(dplyr)
library(patchwork)
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
BiocManager::install("scran")
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))


## Subset of T Cells ----
TCells<-readRDS("./annotated_TCells_T1D_Timepoints_v1.rds")

DimPlot(TCells)

CD4TCells <- subset(TCells, idents = "CD4 T Cell")
CD8TCells <- subset(TCells, idents = "CD8 T Cell")

#### CD4 T Cells ----

##### All CD4 T Cells ----

CD4TCells_Week6 <- subset(CD4TCells, subset = time %in% "Week6")
CD4TCells_Week12 <- subset(CD4TCells, subset = time %in% "Week12")

Idents(CD4TCells_Week6) <- CD4TCells_Week6$group
Idents(CD4TCells_Week12) <- CD4TCells_Week12$group

CD4TCells_Week6 = PrepSCTFindMarkers(CD4TCells_Week6)
CD4TCells_Week12 = PrepSCTFindMarkers(CD4TCells_Week12)

CD4TCells_Week6_PvNP = FindMarkers(CD4TCells_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                    recorrect_umi = TRUE)
CD4TCells_Week6_PvNP$p_val_adj<-p.adjust(CD4TCells_Week6_PvNP$p_val, method = "BH")


CD4TCells_Week12_PvNP = FindMarkers(CD4TCells_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                     recorrect_umi = FALSE)
CD4TCells_Week12_PvNP$p_val_adj<-p.adjust(CD4TCells_Week12_PvNP$p_val, method = "BH")

CD4TCells_Week6_PvNP$symbol = rownames(CD4TCells_Week6_PvNP)
CD4TCells_Week12_PvNP$symbol = rownames(CD4TCells_Week12_PvNP)



###### Week 6----
CD4TCells_Week6_PvNP_data = CD4TCells_Week6_PvNP

keyvals <- ifelse(
  (CD4TCells_Week6_PvNP_data$avg_log2FC < -1 & CD4TCells_Week6_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((CD4TCells_Week6_PvNP_data$avg_log2FC > 1 & CD4TCells_Week6_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'

CD4TCells_Week6_PvNP_data$p_val_adj
EnhancedVolcano(CD4TCells_Week6_PvNP_data,
                lab = rownames(CD4TCells_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(CD4TCells_Week6_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "CD4TCells: Progressor vs Non-Progressor at Week 6",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_CD4TCells_T1DTimepoints_Week6.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(CD4TCells_Week6_PvNP_data,"CD4TCells_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)

###### Week 12----
CD4TCells_Week12_PvNP_data = CD4TCells_Week12_PvNP

keyvals <- ifelse(
  (CD4TCells_Week12_PvNP_data$avg_log2FC < -1 & CD4TCells_Week12_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((CD4TCells_Week12_PvNP_data$avg_log2FC > 1 & CD4TCells_Week12_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'


EnhancedVolcano(CD4TCells_Week12_PvNP_data,
                lab = rownames(CD4TCells_Week12_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(CD4TCells_Week12_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "CD4TCells: Progressor vs Non-Progressor at Week 12",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_CD4TCells_T1DTimepoints_Week12.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(CD4TCells_Week12_PvNP_data,"CD4TCells_T1DTimepoints_DEG_Week12.csv",row.names = TRUE)

##### Tcon Memory Cells ----

unique(CD4TCells$TCellSubType)
TconMemory <- subset(CD4TCells, subset = TCellSubType %in% "Tcon memory")

TconMemory_Week6 <- subset(TconMemory, subset = time %in% "Week6")
TconMemory_Week12 <- subset(TconMemory, subset = time %in% "Week12")

Idents(TconMemory_Week6) <- TconMemory_Week6$group
Idents(TconMemory_Week12) <- TconMemory_Week12$group

TconMemory_Week6 = PrepSCTFindMarkers(TconMemory_Week6)
TconMemory_Week12 = PrepSCTFindMarkers(TconMemory_Week12)

TconMemory_Week6_PvNP = FindMarkers(TconMemory_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                   recorrect_umi = TRUE)
TconMemory_Week6_PvNP$p_val_adj<-p.adjust(TconMemory_Week6_PvNP$p_val, method = "BH")


TconMemory_Week12_PvNP = FindMarkers(TconMemory_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                    recorrect_umi = FALSE)
TconMemory_Week12_PvNP$p_val_adj<-p.adjust(TconMemory_Week12_PvNP$p_val, method = "BH")

TconMemory_Week6_PvNP$symbol = rownames(TconMemory_Week6_PvNP)
TconMemory_Week12_PvNP$symbol = rownames(TconMemory_Week12_PvNP)


###### Week 6----
TconMemory_Week6_PvNP_data = TconMemory_Week6_PvNP

keyvals <- ifelse(
  (TconMemory_Week6_PvNP_data$avg_log2FC < -1 & TconMemory_Week6_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((TconMemory_Week6_PvNP_data$avg_log2FC > 1 & TconMemory_Week6_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'

EnhancedVolcano(TconMemory_Week6_PvNP_data,
                lab = rownames(TconMemory_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(TconMemory_Week6_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "TconMemory: Progressor vs Non-Progressor at Week 6",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_TconMemory_T1DTimepoints_Week6.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(TconMemory_Week6_PvNP_data,"TconMemory_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)


###### Week 12----
TconMemory_Week12_PvNP_data = TconMemory_Week12_PvNP

keyvals <- ifelse(
  (TconMemory_Week12_PvNP_data$avg_log2FC < -1 & TconMemory_Week12_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((TconMemory_Week12_PvNP_data$avg_log2FC > 1 & TconMemory_Week12_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'


EnhancedVolcano(TconMemory_Week12_PvNP_data,
                lab = rownames(TconMemory_Week12_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(TconMemory_Week12_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "TconMemory: Progressor vs Non-Progressor at Week 12",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_TconMemory_T1DTimepoints_Week12.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(TconMemory_Week12_PvNP_data,"TconMemory_T1DTimepoints_DEG_Week12.csv",row.names = TRUE)


##### Tcon effector exhausted Cells ----

unique(CD4TCells$TCellSubType)
Tcon_exhausted_effector <- subset(CD4TCells, subset = TCellSubType %in% "Tcon exhausted effector-like")

Tcon_exhausted_effector_Week6 <- subset(Tcon_exhausted_effector, subset = time %in% "Week6")
Tcon_exhausted_effector_Week12 <- subset(Tcon_exhausted_effector, subset = time %in% "Week12")

Idents(Tcon_exhausted_effector_Week6) <- Tcon_exhausted_effector_Week6$group
Idents(Tcon_exhausted_effector_Week12) <- Tcon_exhausted_effector_Week12$group

Tcon_exhausted_effector_Week6 = PrepSCTFindMarkers(Tcon_exhausted_effector_Week6)
Tcon_exhausted_effector_Week12 = PrepSCTFindMarkers(Tcon_exhausted_effector_Week12)

Tcon_exhausted_effector_Week6_PvNP = FindMarkers(Tcon_exhausted_effector_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                    recorrect_umi = TRUE)
Tcon_exhausted_effector_Week6_PvNP$p_val_adj<-p.adjust(Tcon_exhausted_effector_Week6_PvNP$p_val, method = "BH")


Tcon_exhausted_effector_Week12_PvNP = FindMarkers(Tcon_exhausted_effector_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                     recorrect_umi = TRUE)
Tcon_exhausted_effector_Week12_PvNP$p_val_adj<-p.adjust(Tcon_exhausted_effector_Week12_PvNP$p_val, method = "BH")

Tcon_exhausted_effector_Week6_PvNP$symbol = rownames(Tcon_exhausted_effector_Week6_PvNP)
Tcon_exhausted_effector_Week12_PvNP$symbol = rownames(Tcon_exhausted_effector_Week12_PvNP)


###### Week 6----
Tcon_exhausted_effector_Week6_PvNP_data = Tcon_exhausted_effector_Week6_PvNP

keyvals <- ifelse(
  (Tcon_exhausted_effector_Week6_PvNP_data$avg_log2FC < -1 & Tcon_exhausted_effector_Week6_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Tcon_exhausted_effector_Week6_PvNP_data$avg_log2FC > 1 & Tcon_exhausted_effector_Week6_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'

EnhancedVolcano(Tcon_exhausted_effector_Week6_PvNP_data,
                lab = rownames(Tcon_exhausted_effector_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(Tcon_exhausted_effector_Week6_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "Tcon_exhausted_effector: Progressor vs Non-Progressor at Week 6",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_Tcon_exhausted_effector_T1DTimepoints_Week6.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(Tcon_exhausted_effector_Week6_PvNP_data,"Tcon_exhausted_effector_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)


###### Week 12----
Tcon_exhausted_effector_Week12_PvNP_data = Tcon_exhausted_effector_Week12_PvNP

keyvals <- ifelse(
  (Tcon_exhausted_effector_Week12_PvNP_data$avg_log2FC < -1 & Tcon_exhausted_effector_Week12_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Tcon_exhausted_effector_Week12_PvNP_data$avg_log2FC > 1 & Tcon_exhausted_effector_Week12_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'


EnhancedVolcano(Tcon_exhausted_effector_Week12_PvNP_data,
                lab = rownames(Tcon_exhausted_effector_Week12_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(Tcon_exhausted_effector_Week12_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "Tcon_exhausted_effector: Progressor vs Non-Progressor at Week 12",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_Tcon_exhausted_effector_T1DTimepoints_Week12.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(Tcon_exhausted_effector_Week12_PvNP_data,"Tcon_exhausted_effector_T1DTimepoints_DEG_Week12.csv",row.names = TRUE)


##### Tregs Cells ----

unique(CD4TCells$TCellSubType)
Tregs <- subset(CD4TCells, subset = TCellSubType %in% "Tregs")
table(FetchData(Tregs, vars = c("time", "group")))
Tregs_Week6 <- subset(Tregs, subset = time %in% "Week6")
Tregs_Week12 <- subset(Tregs, subset = time %in% "Week12")

Idents(Tregs_Week6) <- Tregs_Week6$group
Idents(Tregs_Week12) <- Tregs_Week12$group

Tregs_Week6 = PrepSCTFindMarkers(Tregs_Week6)
Tregs_Week12 = PrepSCTFindMarkers(Tregs_Week12)

Tregs_Week6_PvNP = FindMarkers(Tregs_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                                 recorrect_umi = TRUE)
Tregs_Week6_PvNP$p_val_adj<-p.adjust(Tregs_Week6_PvNP$p_val, method = "BH")


Tregs_Week12_PvNP = FindMarkers(Tregs_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                                  recorrect_umi = FALSE)
Tregs_Week12_PvNP$p_val_adj<-p.adjust(Tregs_Week12_PvNP$p_val, method = "BH")

Tregs_Week6_PvNP$symbol = rownames(Tregs_Week6_PvNP)
Tregs_Week12_PvNP$symbol = rownames(Tregs_Week12_PvNP)


###### Week 6----
Tregs_Week6_PvNP_data = Tregs_Week6_PvNP

keyvals <- ifelse(
  (Tregs_Week6_PvNP_data$avg_log2FC < -1 & Tregs_Week6_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Tregs_Week6_PvNP_data$avg_log2FC > 1 & Tregs_Week6_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'

EnhancedVolcano(Tregs_Week6_PvNP_data,
                lab = rownames(Tregs_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(Tregs_Week6_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "Tregs: Progressor vs Non-Progressor at Week 6",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_Tregs_T1DTimepoints_Week6.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(Tregs_Week6_PvNP_data,"Tregs_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)


###### Week 12----
Tregs_Week12_PvNP_data = Tregs_Week12_PvNP

keyvals <- ifelse(
  (Tregs_Week12_PvNP_data$avg_log2FC < -1 & Tregs_Week12_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Tregs_Week12_PvNP_data$avg_log2FC > 1 & Tregs_Week12_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'


EnhancedVolcano(Tregs_Week12_PvNP_data,
                lab = rownames(Tregs_Week12_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(Tregs_Week12_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "Tregs: Progressor vs Non-Progressor at Week 12",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_Tregs_T1DTimepoints_Week12.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(Tregs_Week12_PvNP_data,"Tregs_T1DTimepoints_DEG_Week12.csv",row.names = TRUE)


##### Tcon activated ----

unique(CD4TCells$TCellSubType)
Tcon_activated <- subset(CD4TCells, subset = TCellSubType %in% "Tcon activated ")
table(FetchData(Tcon_activated, vars = c("time", "group")))
Tcon_activated_Week6 <- subset(Tcon_activated, subset = time %in% "Week6")
Tcon_activated_Week12 <- subset(Tcon_activated, subset = time %in% "Week12")

Idents(Tcon_activated_Week6) <- Tcon_activated_Week6$group
Idents(Tcon_activated_Week12) <- Tcon_activated_Week12$group

Tcon_activated_Week6 = PrepSCTFindMarkers(Tcon_activated_Week6)
Tcon_activated_Week12 = PrepSCTFindMarkers(Tcon_activated_Week12)

Tcon_activated_Week6_PvNP = FindMarkers(Tcon_activated_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                               recorrect_umi = TRUE)
Tcon_activated_Week6_PvNP$p_val_adj<-p.adjust(Tcon_activated_Week6_PvNP$p_val, method = "BH")


Tcon_activated_Week12_PvNP = FindMarkers(Tcon_activated_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                recorrect_umi = TRUE)
Tcon_activated_Week12_PvNP$p_val_adj<-p.adjust(Tcon_activated_Week12_PvNP$p_val, method = "BH")

Tcon_activated_Week6_PvNP$symbol = rownames(Tcon_activated_Week6_PvNP)
Tcon_activated_Week12_PvNP$symbol = rownames(Tcon_activated_Week12_PvNP)


###### Week 6----
Tcon_activated_Week6_PvNP_data = Tcon_activated_Week6_PvNP

keyvals <- ifelse(
  (Tcon_activated_Week6_PvNP_data$avg_log2FC < -1 & Tcon_activated_Week6_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Tcon_activated_Week6_PvNP_data$avg_log2FC > 1 & Tcon_activated_Week6_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'

EnhancedVolcano(Tcon_activated_Week6_PvNP_data,
                lab = rownames(Tcon_activated_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(Tcon_activated_Week6_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "Tcon_activated: Progressor vs Non-Progressor at Week 6",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_Tcon_activated_T1DTimepoints_Week6.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(Tcon_activated_Week6_PvNP_data,"Tcon_activated_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)


###### Week 12----
Tcon_activated_Week12_PvNP_data = Tcon_activated_Week12_PvNP

keyvals <- ifelse(
  (Tcon_activated_Week12_PvNP_data$avg_log2FC < -1 & Tcon_activated_Week12_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Tcon_activated_Week12_PvNP_data$avg_log2FC > 1 & Tcon_activated_Week12_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'


EnhancedVolcano(Tcon_activated_Week12_PvNP_data,
                lab = rownames(Tcon_activated_Week12_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(Tcon_activated_Week12_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "Tcon_activated: Progressor vs Non-Progressor at Week 12",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_Tcon_activated_T1DTimepoints_Week12.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(Tcon_activated_Week12_PvNP_data,"Tcon_activated_T1DTimepoints_DEG_Week12.csv",row.names = TRUE)

##### Tcon Interferon Sensing activated ----

unique(CD4TCells$TCellSubType)
TconInterferonSensing <- subset(CD4TCells, subset = TCellSubType %in% "Tcon Interferon Sensing")

table(FetchData(TconInterferonSensing, vars = c("time", "group")))
TconInterferonSensing_Week6 <- subset(TconInterferonSensing, subset = time %in% "Week6")
TconInterferonSensing_Week12 <- subset(TconInterferonSensing, subset = time %in% "Week12")

Idents(TconInterferonSensing_Week6) <- TconInterferonSensing_Week6$group
Idents(TconInterferonSensing_Week12) <- TconInterferonSensing_Week12$group

TconInterferonSensing_Week6 = PrepSCTFindMarkers(TconInterferonSensing_Week6)
TconInterferonSensing_Week12 = PrepSCTFindMarkers(TconInterferonSensing_Week12)

TconInterferonSensing_Week6_PvNP = FindMarkers(TconInterferonSensing_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                        recorrect_umi = TRUE)
TconInterferonSensing_Week6_PvNP$p_val_adj<-p.adjust(TconInterferonSensing_Week6_PvNP$p_val, method = "BH")


TconInterferonSensing_Week12_PvNP = FindMarkers(TconInterferonSensing_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                         recorrect_umi = TRUE)
TconInterferonSensing_Week12_PvNP$p_val_adj<-p.adjust(TconInterferonSensing_Week12_PvNP$p_val, method = "BH")

TconInterferonSensing_Week6_PvNP$symbol = rownames(TconInterferonSensing_Week6_PvNP)
TconInterferonSensing_Week12_PvNP$symbol = rownames(TconInterferonSensing_Week12_PvNP)


###### Week 6----
TconInterferonSensing_Week6_PvNP_data = TconInterferonSensing_Week6_PvNP

keyvals <- ifelse(
  (TconInterferonSensing_Week6_PvNP_data$avg_log2FC < -1 & TconInterferonSensing_Week6_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((TconInterferonSensing_Week6_PvNP_data$avg_log2FC > 1 & TconInterferonSensing_Week6_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'

EnhancedVolcano(TconInterferonSensing_Week6_PvNP_data,
                lab = rownames(TconInterferonSensing_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(TconInterferonSensing_Week6_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "TconInterferonSensing: Progressor vs Non-Progressor at Week 6",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_TconInterferonSensing_T1DTimepoints_Week6.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(TconInterferonSensing_Week6_PvNP_data,"TconInterferonSensing_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)


###### Week 12----
TconInterferonSensing_Week12_PvNP_data = TconInterferonSensing_Week12_PvNP

keyvals <- ifelse(
  (TconInterferonSensing_Week12_PvNP_data$avg_log2FC < -1 & TconInterferonSensing_Week12_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((TconInterferonSensing_Week12_PvNP_data$avg_log2FC > 1 & TconInterferonSensing_Week12_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'


EnhancedVolcano(TconInterferonSensing_Week12_PvNP_data,
                lab = rownames(TconInterferonSensing_Week12_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(TconInterferonSensing_Week12_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "TconInterferonSensing: Progressor vs Non-Progressor at Week 12",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)

# Save the plot
ggsave("Volcanoplot_TconInterferonSensing_T1DTimepoints_Week12.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(TconInterferonSensing_Week12_PvNP_data,"TconInterferonSensing_T1DTimepoints_DEG_Week12.csv",row.names = TRUE)


### CD8 T Cells ----

##### All CD8 T Cells ----

CD8TCells_Week6 <- subset(CD8TCells, subset = time %in% "Week6")
CD8TCells_Week12 <- subset(CD8TCells, subset = time %in% "Week12")

Idents(CD8TCells_Week6) <- CD8TCells_Week6$group
Idents(CD8TCells_Week12) <- CD8TCells_Week12$group

CD8TCells_Week6 = PrepSCTFindMarkers(CD8TCells_Week6)
CD8TCells_Week12 = PrepSCTFindMarkers(CD8TCells_Week12)

CD8TCells_Week6_PvNP = FindMarkers(CD8TCells_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                   recorrect_umi = TRUE)
CD8TCells_Week6_PvNP$p_val_adj<-p.adjust(CD8TCells_Week6_PvNP$p_val, method = "BH")


CD8TCells_Week12_PvNP = FindMarkers(CD8TCells_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                    recorrect_umi = TRUE)
CD8TCells_Week12_PvNP$p_val_adj<-p.adjust(CD8TCells_Week12_PvNP$p_val, method = "BH")

CD8TCells_Week6_PvNP$symbol = rownames(CD8TCells_Week6_PvNP)
CD8TCells_Week12_PvNP$symbol = rownames(CD8TCells_Week12_PvNP)



###### Week 6----
CD8TCells_Week6_PvNP_data = CD8TCells_Week6_PvNP

keyvals <- ifelse(
  (CD8TCells_Week6_PvNP_data$avg_log2FC < -1 & CD8TCells_Week6_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((CD8TCells_Week6_PvNP_data$avg_log2FC > 1 & CD8TCells_Week6_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'

CD8TCells_Week6_PvNP_data$p_val_adj
EnhancedVolcano(CD8TCells_Week6_PvNP_data,
                lab = rownames(CD8TCells_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(CD8TCells_Week6_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "CD8TCells: Progressor vs Non-Progressor at Week 6",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_CD8TCells_T1DTimepoints_Week6.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(CD8TCells_Week6_PvNP_data,"CD8TCells_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)

###### Week 12----
CD8TCells_Week12_PvNP_data = CD8TCells_Week12_PvNP

keyvals <- ifelse(
  (CD8TCells_Week12_PvNP_data$avg_log2FC < -1 & CD8TCells_Week12_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((CD8TCells_Week12_PvNP_data$avg_log2FC > 1 & CD8TCells_Week12_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'


EnhancedVolcano(CD8TCells_Week12_PvNP_data,
                lab = rownames(CD8TCells_Week12_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(CD8TCells_Week12_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "CD8TCells: Progressor vs Non-Progressor at Week 12",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_CD8TCells_T1DTimepoints_Week12.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(CD8TCells_Week12_PvNP_data,"CD8TCells_T1DTimepoints_DEG_Week12.csv",row.names = TRUE)

##### CD8 Memory T Cells ----

unique(CD8TCells$TCellSubType)
CD8Memory <- subset(CD8TCells, subset = TCellSubType %in% "CD8 memory")

table(FetchData(CD8Memory, vars = c("time", "group")))
CD8Memory_Week6 <- subset(CD8Memory, subset = time %in% "Week6")
CD8Memory_Week12 <- subset(CD8Memory, subset = time %in% "Week12")

Idents(CD8Memory_Week6) <- CD8Memory_Week6$group
Idents(CD8Memory_Week12) <- CD8Memory_Week12$group

CD8Memory_Week6 = PrepSCTFindMarkers(CD8Memory_Week6)
CD8Memory_Week12 = PrepSCTFindMarkers(CD8Memory_Week12)


CD8Memory_Week6_PvNP = FindMarkers(CD8Memory_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                   recorrect_umi = TRUE)
CD8Memory_Week6_PvNP$p_val_adj<-p.adjust(CD8Memory_Week6_PvNP$p_val, method = "BH")


CD8Memory_Week12_PvNP = FindMarkers(CD8Memory_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                    recorrect_umi = TRUE)
CD8Memory_Week12_PvNP$p_val_adj<-p.adjust(CD8Memory_Week12_PvNP$p_val, method = "BH")

CD8Memory_Week6_PvNP$symbol = rownames(CD8Memory_Week6_PvNP)
CD8Memory_Week12_PvNP$symbol = rownames(CD8Memory_Week12_PvNP)



###### Week 6----
CD8Memory_Week6_PvNP_data = CD8Memory_Week6_PvNP

keyvals <- ifelse(
  (CD8Memory_Week6_PvNP_data$avg_log2FC < -1 & CD8Memory_Week6_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((CD8Memory_Week6_PvNP_data$avg_log2FC > 1 & CD8Memory_Week6_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'

CD8Memory_Week6_PvNP_data$p_val_adj
EnhancedVolcano(CD8Memory_Week6_PvNP_data,
                lab = rownames(CD8Memory_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(CD8Memory_Week6_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "CD8Memory: Progressor vs Non-Progressor at Week 6",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_CD8Memory_T1DTimepoints_Week6.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(CD8Memory_Week6_PvNP_data,"CD8Memory_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)

###### Week 12----
CD8Memory_Week12_PvNP_data = CD8Memory_Week12_PvNP

keyvals <- ifelse(
  (CD8Memory_Week12_PvNP_data$avg_log2FC < -1 & CD8Memory_Week12_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((CD8Memory_Week12_PvNP_data$avg_log2FC > 1 & CD8Memory_Week12_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'


EnhancedVolcano(CD8Memory_Week12_PvNP_data,
                lab = rownames(CD8Memory_Week12_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(CD8Memory_Week12_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "CD8Memory: Progressor vs Non-Progressor at Week 12",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,6),
                boxedLabels = F,
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_CD8Memory_T1DTimepoints_Week12.png", width = 12, height = 8, dpi = 600)
#dev.off()
write.csv(CD8Memory_Week12_PvNP_data,"CD8Memory_T1DTimepoints_DEG_Week12.csv",row.names = TRUE)




## SCPA Analysis-Week 6 ----

devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org")
devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
BiocManager::install("singscore")
devtools::install_github("jackbibby1/SCPA")
library(SCPA)
library(msigdbr)
library(magrittr)
library(dplyr)
library(Seurat)
library(cowplot)
library(plyr)
library(ggplot2)
library(openxlsx)
library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(SeuratWrappers)
library(colorRamp2)
library(ggplot2)

# Define Gene Sets for Enrichment Analysis
all_gene_sets = msigdbr(species = "Mus musculus")
unique(all_gene_sets$gs_subcat)
unique(all_gene_sets$gs_cat)
# kegg_gene_sets_msig = msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG")
# biocarta_gene_sets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:BIOCARTA")
# Reactome_gene_sets = msigdbr(species = "mouse", subcategory = "CP:REACTOME")
# Wiki_gene_sets = msigdbr(species = "mouse", subcategory = "CP:WIKIPATHWAYS")
# HM_gene_sets = msigdbr(species = "mouse", category = "H")
# GO_gene_sets = msigdbr(species = "mouse", subcategory = "GO:BP")
# combined_data_set = rbind(kegg_gene_sets_msig, 
#                           HM_gene_sets)

library(msigdbr)
library(dplyr)

# Load MSigDB H category pathways
msig_h <- msigdbr("Mus musculus", "H") %>%
  format_pathways()
# Load KEGG pathways from C2 category
kegg_gene_sets_msig <- msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG") %>%
  format_pathways()
# Load Gene Ontology Biological Process (GO:BP) pathways
GO_gene_sets <- msigdbr(species = "mouse", subcategory = "GO:BP") %>%
  format_pathways()

# Combine the pathways into a single data frame
# Combine the pathways, ensuring the format remains the same
pathways <- c(msig_h,kegg_gene_sets_msig)


Macrophage_W6_Progressor <- seurat_extract(Macrophage_Week6,
                                           meta1 = "group", value_meta1 = "Progressor")
Macrophage_W6_NonProgressor <- seurat_extract(Macrophage_Week6,
                                              meta1 = "group", value_meta1 = "Non-Progressor")
scpa_out_Mac_W6 <- compare_pathways(samples = list(Macrophage_W6_Progressor, Macrophage_W6_NonProgressor), 
                                    pathways = pathways)
head(scpa_out_Mac_W6)




library(ggplot2)

# Filter for significant pathways
filtered_data <- scpa_out_Mac_W6 %>%
  filter(Pval < 0.05, adjPval < 0.9,FC>1) %>%
  filter(!Pathway %in% c("KEGG_ALZHEIMERS_DISEASE", "KEGG_ACUTE_MYELOID_LEUKEMIA","HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN" )) %>% 
  arrange(adjPval)  # Sort by adjusted p-value for better visualization

# Create the plot
# Create the plot with improved aesthetics
p <- ggplot(filtered_data, aes(x = -log10(adjPval), y = reorder(Pathway, -log10(adjPval)))) +
  geom_bar(stat = "identity", fill = "#B2182B", color = "black") +  # Red bars with black border
  labs(
    x = expression(bold("-log"[10] ~ "(p"[adj] * "-value)")),
    y = expression(bold("Pathway")),
    title = "Significant Pathways- Macrophages (Early)"
  ) +
  theme_classic(base_size = 16) +  # Remove gridlines, use classic theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  )

# Print the plot
print(p)

# Save the plot in high quality
ggsave("SCPA_UpRegulated_Barplot_Macrophage_W6.png", plot = p, width = 12, height = 6, dpi = 600)



write.csv(scpa_out_Mac_W6,"Macrophage_T1DTimepoints_SCPA_Week6.csv",row.names = TRUE)

