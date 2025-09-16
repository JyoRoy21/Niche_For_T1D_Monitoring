## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/")
library(Seurat)
packageVersion("Seurat")
library(SoupX)
library(scRNAseq)
#devtools::install_github("ayshwaryas/ddqc_R")
library(ddqcR)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(parallel)
library(purrr)
library(tibble)
#install_github('immunogenomics/presto')
library(presto)
library(dplyr)
library(patchwork)
library(plyr)
library(RColorBrewer)
#BiocManager::install("multtest")
library(multtest)
library(metap) 
library(ggprism)
#BiocManager::install('glmGamPoi')
library(glmGamPoi)

library(dplyr)
library(Seurat)
library(patchwork)

# Annotation- the NOD Rag ----
NODRag.data <- Read10X(data.dir = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/")
# Initialize the Seurat object with the raw (non-normalized data).
NODRag <- CreateSeuratObject(counts = NODRag.data, project = "NODRag", min.cells = 3, min.features = 200)
NODRag
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
NODRag[["percent.mt"]] <- PercentageFeatureSet(NODRag, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(NODRag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(NODRag, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NODRag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

NODRag <- subset(NODRag, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
NODRag  <- SCTransform(NODRag, vars.to.regress = "percent.mt", verbose=F)

NODRag <- RunPCA(NODRag)
NODRag <- RunUMAP(NODRag, dims = 1:30)
NODRag[[]]

ElbowPlot(NODRag, ndims = 50)
NODRag <- FindNeighbors(NODRag, dims = 1:30, verbose = F)
NODRag <- FindClusters(NODRag, verbose =F, resolution = 0.4)
DimPlot(NODRag, label=T)
NODRag$sample<-NODRag$orig.ident

NODRag <- PrepSCTFindMarkers(NODRag)
all.markers <- FindAllMarkers(NODRag , only.pos = TRUE)

top20_sct=all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)

# Plot the heatmap for the top 20 markers
DoHeatmap(
  object = NODRag, 
  features = top20_sct$gene,  # Use the top 20 marker genes
  size = 4,                   # Adjust text size for readability
  group.by = "ident",         # Group cells by cluster identity
  label = TRUE                # Show cluster labels
) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) # Adjust color scale








# Extract top genes for cluster 0
cluster_0_genes <- top20_sct %>% 
  filter(cluster == 0) %>% 
  pull(gene)

# Extract top genes for cluster 110
cluster_110_genes <- top20_sct %>% 
  filter(cluster == 110) %>% 
  pull(gene)

# Print genes for cluster 0
cat("Top 20 genes for Cluster 0:\n")
print(cluster_0_genes)

# Print genes for cluster 110
cat("\nTop 20 genes for Cluster 110:\n")
print(cluster_110_genes)


DotPlot(NODRag, features = c( "Ptprc", #CD45
                                      "Trbc1","Trbc2","Cd3g","Cd3e","Cd3d","Lat", #T Cells
                                      "C1qa","C1qb","C1qc","Lyz2","Adgre1", #Macrophages
                                      "Cd79a", "Cd79b","Ms4a1","Cd19",#B Cells
                                      "S100a4","S100a11","Mdh2","Gm2a","Pglyrp1","H2-Oa",#cDc
                                      "Siglech","Cd209a","Cox6a2","Rnase6",#pDC
                                      "Igj","Ighg2b","Ighg2c","Igha","Jchain",#Plasma
                                      #"Ctrb1","Cela3b","Cpa1","Prss1",#Acinar
                                      "Cd34","Eng","Col15a1","Tie1","Kdr","Pecam",
                                      "Il34","Pdgfrb","Pdgfa","Col3a1" #Mesenchymal
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")



NODRag = RenameIdents(NODRag, 
                              "0" = "Endothelial Cells",
                              "1" = "Macrophage",
                              "2" = "Mesenchymal-like",
                              "3" = "Endothelial Cells",
                              "4" = "Macrophage",
                              "5" = "Macrophage",
                              "6" = "T Cell",
                              "7" = "Dendritic Cells",
                              "8" = "Mesenchymal-like",
                              "9" = "Dendritic Cells",
                              "10" = "Mesenchymal-like",
                              "11" = "Endothelial Cells",
                              "12" = "T Cell",
                              "13" = "Macrophage",
                              "14" = "Macrophage"
)



DimPlot(NODRag,label = T,label.size = 8)
NODRag$CellType<-Idents(NODRag)
saveRDS(NODRag,file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/NODRag.rds")

intersect(rownames(T1D_Timepoints), rownames(NODRag)) %>% length()
T1D_Timepoints_NODRag <- merge(T1D_Timepoints, NODRag, merge.data = TRUE, assay = "RNA")

## ----Integration of NOD Rah and T1DTimepoints-------------------------------------------------------------------------------------------------------------
T1D_Timepoints = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v4.rds")
NODRag_anno<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/NODRag.rds")
unique(colnames(T1D_Timepoints@meta.data))
unique(colnames(NODRag@meta.data))

library(Seurat)
library(dplyr)


NODRag.data <- Read10X(data.dir = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/")
# Initialize the Seurat object with the raw (non-normalized data).
NODRag <- CreateSeuratObject(counts = NODRag.data, project = "NODRag", min.cells = 3, min.features = 200)
NODRag
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
NODRag[["percent.mt"]] <- PercentageFeatureSet(NODRag, pattern = "^mt-")
NODRag <- subset(NODRag, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
# Ensure both Seurat objects have the same cell names
common_cells <- intersect(Cells(NODRag), Cells(NODRag_anno))
# Subset the CellType metadata from NODRag_anno based on common cells
CellType_metadata <- NODRag_anno@meta.data[common_cells, "CellType", drop = FALSE]
# Add the CellType metadata to NODRag
NODRag@meta.data <- cbind(NODRag@meta.data, CellType_metadata)
# Check that the new metadata has been added correctly
head(NODRag@meta.data)
NODRag$sample <-NODRag$orig.ident
NODRag$group<-NODRag$sample
NODRag$time <- NODRag$sample
T1D_Timepoints$CellType<-Idents(T1D_Timepoints)
# Define columns to keep
columns_to_keep <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt",
                     "sample","time","group","CellType")

# Subset metadata
T1D_Timepoints@meta.data <- T1D_Timepoints@meta.data %>% select(all_of(columns_to_keep))
NODRag@meta.data <- NODRag@meta.data %>% select(all_of(columns_to_keep))
unique(colnames(T1D_Timepoints@meta.data))
unique(colnames(NODRag@meta.data))
#Change Assay to RNA
DefaultAssay(NODRag)<-"RNA"
DefaultAssay(T1D_Timepoints)<-"RNA"
NODRag<- NormalizeData(NODRag)
NODRag <- ScaleData(NODRag, verbose = F)


options(future.globals.maxSize = 50 * 1024^3)  # Set to 2 GB, for example
# Load the future package if not already loaded
library(future)

# Set up a multi-core plan for parallel processing
#plan("multicore", workers = 4)  # Adjust the number of workers as needed

T1D_Timepoints_NODRag = merge(x = T1D_Timepoints, y = NODRag, add.cell.ids = c("T1D_Timepoints", "NODRag"))

T1D_Timepoints_NODRag  <- SCTransform(T1D_Timepoints_NODRag , vars.to.regress = "percent.mt", verbose=T)

# Increase the max size limit for future globals
saveRDS(T1D_Timepoints_NODRag  , file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/T1D_Timepoints_NODRag_SCT_merged.rds")
T1D_Timepoints_NODRag<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/T1D_Timepoints_NODRag_SCT_merged.rds")


T1D_Timepoints_NODRag <- RunPCA(T1D_Timepoints_NODRag)
T1D_Timepoints_NODRag <- RunUMAP(T1D_Timepoints_NODRag, dims = 1:30)
T1D_Timepoints_NODRag[[]]
DimPlot(T1D_Timepoints_NODRag, reduction = "umap", group.by = c("sample","CellType"))

ElbowPlot(T1D_Timepoints_NODRag, ndims = 50)
T1D_Timepoints_NODRag <- FindNeighbors(T1D_Timepoints_NODRag, dims = 1:30, verbose = F)
T1D_Timepoints_NODRag <- FindClusters(T1D_Timepoints_NODRag, verbose =F, resolution = 0.4)
DimPlot(T1D_Timepoints_NODRag, group.by = c("sample","CellType"), label=F)

Macrophage_T1D_Timepoints_NODRag <- subset(T1D_Timepoints_NODRag, CellType %in% c("Macrophage"))
DimPlot(Macrophage_T1D_Timepoints_NODRag, group.by = c("sample","CellType"), label=F)

#' # Integration
#' Seurat integration vignette: https://satijalab.org/seurat/articles/integration_introduction.html
#' Use integration method for SCtransform, not the regular one 
## ----integration clustering---------------------------------------------------------------------------------------------------------------------------

int.T1D_Timepoints_NODRag<- IntegrateLayers(object = T1D_Timepoints_NODRag, method = CCAIntegration, normalization.method = "SCT", verbose =T)
int.T1D_Timepoints_NODRag <- FindNeighbors(int.T1D_Timepoints_NODRag, reduction = "integrated.dr", dims = 1:30)
int.T1D_Timepoints_NODRag <- FindClusters(int.T1D_Timepoints_NODRag, res = 0.8)

#Mac
int.Macrophage_T1D_Timepoints_NODRag<- IntegrateLayers(object = Macrophage_T1D_Timepoints_NODRag, method = CCAIntegration, normalization.method = "SCT", verbose =T)
int.Macrophage_T1D_Timepoints_NODRag <- FindNeighbors(int.Macrophage_T1D_Timepoints_NODRag, reduction = "integrated.dr", dims = 1:30)
int.Macrophage_T1D_Timepoints_NODRag <- FindClusters(int.Macrophage_T1D_Timepoints_NODRag, res = 0.4)


## ----cluster integrated DimPlot-----------------------------------------------------------------------------------------------------------------------
int.T1D_Timepoints_NODRag <- RunUMAP(int.T1D_Timepoints_NODRag, dims = 1:30, reduction = "integrated.dr")
DimPlot(int.T1D_Timepoints_NODRag, reduction = "umap", group.by = c( "sample","CellType"), combine = F)

int.Macrophage_T1D_Timepoints_NODRag <- RunUMAP(int.Macrophage_T1D_Timepoints_NODRag, dims = 1:30, reduction = "integrated.dr")
DimPlot(int.Macrophage_T1D_Timepoints_NODRag, reduction = "umap", group.by = c( "sample","CellType"), combine = F)

table(int.Macrophage_T1D_Timepoints_NODRag$group)

DefaultAssay(int.Macrophage_T1D_Timepoints_NODRag)
getPalette = colorRampPalette(brewer.pal(17, "Set1"))


## ----save integrated object---------------------------------------------------------------------------------------------------------------------------
saveRDS(int.Macrophage_T1D_Timepoints_NODRag, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/int.Macrophage_T1D_Timepoints_NODRag.rds")
getwd()



# Subclustering of Macrophage Subtypes----
int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag <-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/int.Macrophage_T1D_Timepoints_NODRag.rds")
DimPlot(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag,group.by = "group")
ElbowPlot(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag, ndims = 50)
int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag <- RunUMAP(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag, dims = 1:30)
int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag <- FindNeighbors(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag, dims = 1:30, verbose = F)
int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag <- FindClusters(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag, res = 0.1)#graph.name = "SCT_nn"
#table(Macrophage$SCT_snn_res.0.1)

DimPlot(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag, reduction = "umap", combine = F, label = T,repel = TRUE,label.box =TRUE,
        label.size = 8,alpha=1,pt.size = 2)

DimPlot(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag, reduction = "umap", group.by = "group",combine = F, label = T,repel = TRUE,label.box =TRUE,
        label.size = 8,alpha=1,pt.size = 2)

table(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag$group)

DefaultAssay(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag) = "SCT"
int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag <- PrepSCTFindMarkers(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag)
all.markers <- FindAllMarkers(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag , only.pos = TRUE,recorrect_umi = FALSE)

top20_sct <- all.markers %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.5) %>% 
  top_n(n = 20, wt = avg_log2FC)




DotPlot(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag,assay = "SCT",features = c("Apoe","Selenop","Maf","Mrc1","F13a1","Csf1r",#Cluster 0-Mac1 --Anti-Inflammatory like
                                              "H2-Ab1", "H2-Aa","Cd74","Tyrobp","Itgb2","Cd52","Gm2a", #Cluster 2-Mac 2-Inflammatory, Antigen Presenting
                                              "Clec9a","Tap1","Psmb8","Igkc","Shisa5","Ets1","Limd2","Lgals3" #Cluster 1-Mac-3- Tissue Resident
                                              
                                              
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")


DotPlot(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag, assay = "SCT", features = c(
  # Inflammatory markers
  "Ccl2", "Ccl3", "Ccl4", "Nos2", "Il1b", "Tnf", "Fcgr1",  
  # Regulatory (M2-like) and Tissue-Resident Markers
  "Arg1", "Chil3", "Retnla", "Il10", "Tgfb1", "Cd163", "Marco", 
  # Dendritic-like Macrophage Markers
  "Batf3", "Irf8", "Flt3", "Zbtb46", "Cd80", "Cd86",  
  # Monocyte-Derived/Recruited Macrophage Markers
  "Cx3cr1", "Itgam", "Tgfbi", "Spp1"  
)) +  
  RotatedAxis() +  
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

# Define the Mac subtypes for each cluster
MacSubtype <- c("Mac-1", 
                "Mac-2",
                "Mac-3",
                "Mac-4",
                "Mac-5")        

# Ensure the Idents are in numeric format
ident_integers <- as.integer(Idents(bar_colors))

# Check if there are any unexpected levels
unique(ident_integers)

# Add the Macrophage subtype information to the metadata of the Seurat object
int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag$MacType <- MacSubtype[as.integer(Idents(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag))]

all.markers <- FindAllMarkers(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag , only.pos = TRUE,recorrect_umi = FALSE)
top20_sct <- all.markers %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.5) %>% 
  top_n(n = 20, wt = avg_log2FC)

## Plot the barplot ----

Idents(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag)<-int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag$MacType
# Create frequency tables for Macrophage subtypes of Week 6 using active.ident and group
freq_Mac <- prop.table(x = table(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag@active.ident,int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag@meta.data[, "group"]), margin =2)

library(ggsci)
# Extract the colors from the plot
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=5)
#color_list <- pal_npg("nrc", alpha = 0.8)(5)
int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag$MacType<-as.factor(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag$MacType)
cell_types<-levels(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag$MacType)
# Create frequency tables for Week 18 using active.ident and group
# Create a vector of colors for each bar
named_colors<-setNames(color_list,cell_types)
bar_colors <- named_colors[rownames(freq_Mac)]

# Plot for Week 6
par(mar = c(5, 6, 4, 6))  # Increasing the right margin to 6

# Plot the bar plot with modified x-axis label rotation
barplot(
  freq_Mac, 
  xlim = c(0, ncol(freq_Mac) + 3),
  col = bar_colors, 
  legend.text = TRUE, 
  args.legend = list(
    x = ncol(freq_Mac) + 3, 
    y = max(colSums(freq_Mac)), 
    bty = "n", 
    cex = 1.5
  ), 
  width = 1.25, 
  ylab = "Fraction of Macrophages",
  cex.lab = 2,   # Increase font size of axis titles
  cex.names = 2,  # Increase font size of x-axis labels
  cex.axis = 1.5     # Increase font size of axis tick labels
)

table(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag$group,int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag$MacType)



## ----Non-Progressor vs NOD.Rag----

int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag <- subset(int.Macrophage_T1D_Timepoints_NODRag, subset = group %in% c("NODRag","Non-Progressor"))
int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag = PrepSCTFindMarkers(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag)
table(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag$group)
Idents(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag)<-int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag$group
Macrophage_NPvsNODRag = FindMarkers(int.Macrophage_T1D_Timepoints_NonProgressorVsNODRag, 
                            ident.1 = "Non-Progressor", 
                            ident.2 = "NODRag", 
                            assay = "SCT",
                            recorrect_umi = FALSE)
# Apply Benjamini-Hochberg correction
Macrophage_NPvsNODRag$p_val_adj <- p.adjust(Macrophage_NPvsNODRag$p_val, method = "BH")



keyvals <- ifelse(
  (Macrophage_NPvsNODRag$avg_log2FC < -1 & Macrophage_NPvsNODRag$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Macrophage_NPvsNODRag$avg_log2FC > 1 & Macrophage_NPvsNODRag$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non-Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Non-Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Upregulated- NOD.Rag'


#png(file = "AllergyScaf_NPPBSCh4_VolcanoPlot.png",units = "in",width = 10, height = 10, res = 400)
EnhancedVolcano(Macrophage_NPvsNODRag,
                lab = rownames(Macrophage_NPvsNODRag),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colCustom = keyvals,
                selectLab = rownames(Macrophage_NPvsNODRag)[which(names(keyvals) %in% c('Upregulated- Non-Progressor', 'Upregulated- NOD.Rag'))],
                title = "Macrophage: Non-Progressor vs NOD.Rag",
                drawConnectors = F,
                widthConnectors = 0.75,
                ylab = bquote(~-Log[10] ~ italic(P[adj])),
                boxedLabels = F)
dev.off()
write.csv(Macrophage_NPvsNODRag,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/Macrophage_NPvsNODRag.csv",row.names = TRUE)


## ----Progressor vs NOD.Rag---------------------------------------------------------------------------------------------------------------------------

int.Macrophage_T1D_Timepoints_ProgressorVsNODRag <- subset(int.Macrophage_T1D_Timepoints_NODRag, subset = group %in% c("NODRag","Progressor"))
int.Macrophage_T1D_Timepoints_ProgressorVsNODRag = PrepSCTFindMarkers(int.Macrophage_T1D_Timepoints_ProgressorVsNODRag)
table(int.Macrophage_T1D_Timepoints_ProgressorVsNODRag$group)
Idents(int.Macrophage_T1D_Timepoints_ProgressorVsNODRag)<-int.Macrophage_T1D_Timepoints_ProgressorVsNODRag$group
Macrophage_PvsNODRag = FindMarkers(int.Macrophage_T1D_Timepoints_ProgressorVsNODRag, 
                                    ident.1 = "Progressor", 
                                    ident.2 = "NODRag", 
                                    assay = "SCT",
                                    recorrect_umi = FALSE)
# Apply Benjamini-Hochberg correction
Macrophage_PvsNODRag$p_val_adj <- p.adjust(Macrophage_PvsNODRag$p_val, method = "BH")



keyvals <- ifelse(
  (Macrophage_PvsNODRag$avg_log2FC < -1 & Macrophage_PvsNODRag$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Macrophage_PvsNODRag$avg_log2FC > 1 & Macrophage_PvsNODRag$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non-Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Upregulated- NOD.Rag'


#png(file = "AllergyScaf_NPPBSCh4_VolcanoPlot.png",units = "in",width = 10, height = 10, res = 400)
EnhancedVolcano(Macrophage_PvsNODRag,
                lab = rownames(Macrophage_PvsNODRag),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colCustom = keyvals,
                selectLab = rownames(Macrophage_PvsNODRag)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Upregulated- NOD.Rag'))],
                title = "Macrophage: Progressor vs NOD.Rag",
                drawConnectors = F,
                widthConnectors = 0.75,
                ylab = bquote(~-Log[10] ~ italic(P[adj])),
                boxedLabels = F)
dev.off()
write.csv(Macrophage_PvsNODRag,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/Macrophage_PvsNODRag.csv",row.names = TRUE)

## Venn Diagram ----

# Install and load required packages
install.packages("VennDiagram")
library(VennDiagram)
library(grid)

# Define upregulated and downregulated genes for each comparison
Upregulated_NPvsNODRag <- rownames(Macrophage_NPvsNODRag)[Macrophage_NPvsNODRag$avg_log2FC > 1.5 & Macrophage_NPvsNODRag$p_val_adj < 0.05]
Downregulated_NPvsNODRag <- rownames(Macrophage_NPvsNODRag)[Macrophage_NPvsNODRag$avg_log2FC < -1.5 & Macrophage_NPvsNODRag$p_val_adj < 0.05]

Upregulated_PvsNODRag <- rownames(Macrophage_PvsNODRag)[Macrophage_PvsNODRag$avg_log2FC > 1.5 & Macrophage_PvsNODRag$p_val_adj < 0.05]
Downregulated_PvsNODRag <- rownames(Macrophage_PvsNODRag)[Macrophage_PvsNODRag$avg_log2FC < -1.5 & Macrophage_PvsNODRag$p_val_adj < 0.05]

# Create a high-quality Venn diagram for upregulated genes with consistent circle positioning
venn_up <- venn.diagram(
  x = list(
    "Non-Progressor vs NOD.Rag" = Upregulated_NPvsNODRag, 
    "Progressor vs NOD.Rag" = Upregulated_PvsNODRag
  ),
  filename = NULL,
  fill = c("#3B9AB2", "#E41A1C"), # Soft blue and red shades
  alpha = 0.6,  # Adjust transparency
  cex = 1.6,  # Increase text size
  cat.cex = 1.5, # Category text size
  cat.pos = c(0, 0), # Set category labels to be on the same side
  cat.dist = c(0.1, 0.1), # Distance of category labels from the circles
  main = "Upregulated Genes", # Title
  main.cex = 2, # Title font size
  main.fontface = "bold", # Bold title
  sub.cex = 1.4, # Subtitle font size
  sub.fontface = "italic", # Subtitle style
  lwd = 2, # Line width
  lty = 1, # Line type (solid)
  margin = 0.05, # Adjust diagram margin
  col = "black", # Outline color
  x.pos = c(0.5, 1.5), # Adjust position of circles to be aligned consistently with downregulated diagram
  y.pos = c(0.5, 0.5) # Keep y position consistent
)
dev.off()

grid.draw(venn_up)

dev.off()
# Create a high-quality Venn diagram for downregulated genes with circles swapped
venn_down <- venn.diagram(
  x = list(
    "Non-Progressor vs NOD.Rag" = Downregulated_PvsNODRag,  # Swap the two sets here
    "Progressor vs NOD.Rag" = Downregulated_NPvsNODRag   # and here to switch the circles
  ),
  filename = NULL,
  fill = c("#3B9AB2", "#E41A1C"), # Soft blue and red shades
  alpha = 0.6,  # Adjust transparency
  cex = 1.6,  # Increase text size
  cat.cex = 1.5, # Category text size
  cat.dist = c(0.1, 0.1), # Distance of category labels from the circles
  main = "Downregulated Genes", # Title
  main.cex = 2, # Title font size
  main.fontface = "bold", # Bold title
  sub.cex = 1.4, # Subtitle font size
  sub.fontface = "italic", # Subtitle style
  lwd = 2, # Line width
  lty = 1, # Line type (solid)
  margin = 0.05, # Adjust diagram margin
  col = "black", # Outline color
  x.pos = c(0.5, 1.5), # Same positioning for alignment
  y.pos = c(0.5, 0.5) # Keep vertical position consistent
)
dev.off()
# Display the Venn diagram
grid.draw(venn_down)

# Load necessary library
library(dplyr)

# Function to categorize genes
generate_gene_csv <- function(set1, set2, label1, label2, filename) {
  # Create a data frame with categories
  df <- data.frame(
    Gene = unique(c(set1, set2)),  # All genes from both sets
    Category = sapply(unique(c(set1, set2)), function(gene) {
      if (gene %in% set1 & gene %in% set2) {
        return("Common")
      } else if (gene %in% set1) {
        return(label1)
      } else {
        return(label2)
      }
    })
  )
  
  # Save as CSV
  write.csv(df, file = filename, row.names = FALSE)
}

# Create CSV for Upregulated genes
generate_gene_csv(
  Upregulated_PvsNODRag, Upregulated_NPvsNODRag,
  "Unique to Progressor vs NOD.Rag", "Unique to Non-Progressor vs NOD.Rag",
  "Upregulated_Genes_Categorized.csv"
)

# Create CSV for Downregulated genes
generate_gene_csv(
  Downregulated_PvsNODRag, Downregulated_NPvsNODRag,
  "Unique to Progressor vs NOD.Rag", "Unique to Non-Progressor vs NOD.Rag",
  "Downregulated_Genes_Categorized.csv"
)




## ----Progressor vs Non-Progressor---------------------------------------------------------------------------------------------------------------------------

int.Macrophage_T1D_Timepoints_ProgressorVsNonProgressor <- subset(int.Macrophage_T1D_Timepoints_NODRag, subset = group %in% c("Progressor","Non-Progressor"))
int.Macrophage_T1D_Timepoints_ProgressorVsNonProgressor = PrepSCTFindMarkers(int.Macrophage_T1D_Timepoints_ProgressorVsNonProgressor)
table(int.Macrophage_T1D_Timepoints_ProgressorVsNonProgressor$group)
Idents(int.Macrophage_T1D_Timepoints_ProgressorVsNonProgressor)<-int.Macrophage_T1D_Timepoints_ProgressorVsNonProgressor$group
Macrophage_PvsNP = FindMarkers(int.Macrophage_T1D_Timepoints_ProgressorVsNonProgressor, 
                                   ident.1 = "Progressor", 
                                   ident.2 = "Non-Progressor", 
                                   assay = "SCT",
                                   recorrect_umi = FALSE)
# Apply Benjamini-Hochberg correction
Macrophage_PvsNP$p_val_adj <- p.adjust(Macrophage_PvsNP$p_val, method = "BH")



keyvals <- ifelse(
  (Macrophage_PvsNP$avg_log2FC < -1.5 & Macrophage_PvsNP$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Macrophage_PvsNP$avg_log2FC > 1.5 & Macrophage_PvsNP$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non-Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Upregulated- Non-Progressor'


#png(file = "AllergyScaf_NPPBSCh4_VolcanoPlot.png",units = "in",width = 10, height = 10, res = 400)
EnhancedVolcano(Macrophage_PvsNP,
                lab = rownames(Macrophage_PvsNP),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoffCol="p_val_adj",
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colCustom = keyvals,
                selectLab = rownames(Macrophage_PvsNP)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Upregulated- Non-Progressor'))],
                title = "Macrophage: Progressor vs Non-Progressor",
                drawConnectors = F,
                widthConnectors = 0.75,
                ylab = bquote(~-Log[10] ~ italic(P[adj])),
                boxedLabels = F)
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Volcanoplot_Macrophage_T1DTimepoints_integrated.png", width = 12, height = 8, dpi = 600)

dev.off()
write.csv(Macrophage_PvsNP,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RagNOD/Macrophage_PvsNP_intgreated.csv",row.names = TRUE)





## Seurat to Scanpy ----
getwd()
library(SeuratData)
library(SeuratDisk)
int.Macrophage_T1D_Timepoints_NODRag[["RNA_counts"]] <- CreateAssayObject(counts = GetAssayData(int.Macrophage_T1D_Timepoints_NODRag, slot = "counts"))
int.Macrophage_T1D_Timepoints_NODRag[["RNA_counts"]]$data
SaveH5Seurat(int.Macrophage_T1D_Timepoints_NODRag, filename = "Macrophage_NODRag_int.h5Seurat")
Convert("Macrophage_NODRag_int.h5Seurat", dest = "h5ad",assay = "RNA_counts")
