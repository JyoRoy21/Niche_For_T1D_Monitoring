#' ---
#' title: "Pre-Diabetes to Diabetes Notebook"
#' author: "Jyotirmoy Roy adapted from Kate Griffin"
#' Date: "`r Sys.Date()`"
#' output:
#'   pdf_document: default
#'   html_notebook: default
#' editor_options: 
#'   markdown: 
#'     wrap: 72
#' ---

## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq")
library(Seurat)
packageVersion("Seurat")
library(SoupX)
devtools::install_github("ayshwaryas/ddqc_R")
library(ddqcR)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(parallel)
library(purrr)
library(tibble)
install_github('immunogenomics/presto')
library(presto)
library(dplyr)
library(patchwork)
library(plyr)
library(RColorBrewer)
BiocManager::install("multtest")
library(multtest)
library(metap) 
library(ggprism)
BiocManager::install('glmGamPoi')
library(glmGamPoi)
#' 
## ----save merged and cluster unintegrated-------------------------------------------------------------------------------------------------------------
T1D_Week18 = readRDS("./T1D_Week18_sct_filtered.rds")


T1D_Week18 <- RunPCA(T1D_Week18)
T1D_Week18 <- RunUMAP(T1D_Week18, dims = 1:30)
T1D_Week18[[]]
DimPlot(T1D_Week18, reduction = "umap", group.by = c("sample", "seurat_clusters"))
DimPlot(T1D_Week18, reduction = "umap", split.by = c("sample"))

ElbowPlot(T1D_Week18, ndims = 50)
T1D_Week18 <- FindNeighbors(T1D_Week18, dims = 1:30, verbose = F)
T1D_Week18 <- FindClusters(T1D_Week18, verbose =F, resolution = 0.4)
DimPlot(T1D_Week18, group.by = c("sample"), label=F)


table(Idents(T1D_Week18),T1D_Week18@meta.data$sample)
#' 
#' # Integration
#' Seurat integration vignette: https://satijalab.org/seurat/articles/integration_introduction.html
#' Use integration method for SCtransform, not the regular one 
## ----integration clustering---------------------------------------------------------------------------------------------------------------------------
int.T1D_Week18 <- IntegrateLayers(object = T1D_Week18, method = CCAIntegration, normalization.method = "SCT", verbose =T)
int.T1D_Week18 <- FindNeighbors(int.T1D_Week18, reduction = "integrated.dr", dims = 1:30)
int.T1D_Week18 <- FindClusters(int.T1D_Week18, res = 0.8)

 
## ----cluster integrated DimPlot-----------------------------------------------------------------------------------------------------------------------
int.T1D_Week18 <- RunUMAP(int.T1D_Week18, dims = 1:30, reduction = "integrated.dr")
DimPlot(int.T1D_Week18, reduction = "umap", group.by = c( "seurat_clusters"), split.by = c( "time"), combine = F)


getPalette = colorRampPalette(brewer.pal(17, "Set1"))

#png(file = "UMAP_31clusters.png",    units = "in",width = 11, height = 11, res = 400)
# DimPlot(int.T1D_Week18, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T, 
#         group.by = "seurat_clusters",
#         cols = getPalette(31)) +
#   theme_prism(base_size = 18) + theme(legend.text = element_text(size = 28)) +
#   labs( title = "UMAP of Cell Types") 
# ggsave(file = "UMAP_31clusters.png",
#     units = "in",width = 11, height = 11, dpi = 400)
# dev.off()
# 
# dev.set(dev.prev())



# visualize

DimPlot(int.T1D_Week18, reduction = "umap", group.by = c("sample"))
DimPlot(int.T1D_Week18, reduction = "umap", split.by = c("sample"))
DimPlot(int.T1D_Week18, reduction = "umap", group.by = c( "seurat_clusters"), split.by = c( "time"), combine = F)
DimPlot(T1D_Week18, reduction = "umap", group.by = c( "seurat_clusters"), split.by = c( "time"), combine = F)




#' 
## ----table of integrated clusters---------------------------------------------------------------------------------------------------------------------
table(Idents(int.T1D_Week18),int.T1D_Week18@meta.data$sample)

#' 
#' 
#' Save integrated object
## ----save integrated object---------------------------------------------------------------------------------------------------------------------------
saveRDS(int.T1D_Week18, file = "./sct_integrated_T1D_Week18_v2.rds")
getwd()

#' 
#' 
## ----save workspace-----------------------------------------------------------------------------------------------------------------------------------
#save.image("./JR_KGpipeline_integrated.T1D_Week18.RData")

#' 
## ----rejoin layers and find conserved markers---------------------------------------------------------------------------------------------------------
# now that integration is complete, rejoin layers
Layers(int.T1D_Week18)
int.T1D_Week18
DefaultAssay(int.T1D_Week18) = "RNA"
int.T1D_Week18 <- JoinLayers(int.T1D_Week18)




get_conserved <- function(cluster){
  Seurat::FindConservedMarkers(int.T1D_Week18,
                               ident.1 = cluster,
                               grouping.var = "sample",
                               only.pos = TRUE) %>%
    rownames_to_column(var = "gene")  %>%
    #left_join(y = unique(annotations[, c("gene_name", "description")]),
    #          by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
#devtools::install_github('immunogenomics/presto')
FindConservedMarkers(int.T1D_Week18,
                     ident.1 = 0,
                     grouping.var = "sample",
                     only.pos = TRUE)

#conserved_markers <- map_dfr(unique(new.cluster.ids), get_conserved)
install.packages("tibble")
library(purrr)
library(tibble)
DimPlot(int.T1D_Week18,reduction = "umap")
conserved_markers <- map_dfr(0:20, get_conserved)


head(conserved_markers)
conserved_markers$Week18_4_avg_log2FC
# Extract top 10 markers per cluster
top20 <- conserved_markers %>% 
  mutate(avg_fc = (`Week18_1_avg_log2FC` + 
                     `Week18_2_avg_log2FC` + 
                     `Week18_3_avg_log2FC` +
                     `Week18_4_avg_log2FC` +
                     `Week18_5_avg_log2FC` +
                     `Week18_6_avg_log2FC`) /6) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)

DefaultAssay(int.T1D_Week18) = "SCT"


## ----cluster ID---------------------------------------------------------------------------------------------------------------------------------------



int.T1D_Week18 <- PrepSCTFindMarkers(int.T1D_Week18)
all.markers <- FindAllMarkers(int.T1D_Week18 , only.pos = TRUE)

top20_sct=all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)

# Plot the heatmap for the top 20 markers
DoHeatmap(
  object = int.T1D_Week18, 
  features = top20_sct$gene,  # Use the top 20 marker genes
  size = 4,                   # Adjust text size for readability
  group.by = "ident",         # Group cells by cluster identity
  label = TRUE                # Show cluster labels
) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) # Adjust color scale


# Extract top genes for cluster 6
cluster_6_genes <- top20_sct %>% 
  filter(cluster == 6) %>% 
  pull(gene)


# Print genes for cluster 6
cat("Top 20 genes for Cluster 6:\n")
print(cluster_6_genes)




DotPlot(int.T1D_Week18, features = c( "Ptprc", #CD45
                                      "Trbc1","Trbc2","Cd3g","Cd3e","Cd3d","Lat", #T Cells
                                      "C1qa","C1qb","C1qc","Lyz2","Adgre1", #Macrophages
                                      "Cd79a", "Cd79b","Ms4a1","Cd19",#B Cells
                                      "S100a4","S100a11","Itgax","Zbtb46","Flt3","H2-Aa","Cd74","Mdh2","Gm2a","Pglyrp1","H2-Oa",#cDc
                                      "Siglech","Cd209a","Cox6a2","Rnase6",#pDC
                                      "H2-M2","Cd209c",#Monocyte Derived DC
                                      "Ighg2b","Ighg2c","Igha","Igkc","Jchain",#Plasma
                                      "Ctrb1","Cela3b","Cpa1","Prss1","Fcgr1","Itgae"#Acinar
                                      # "Il34","Pdgfrb","Pdgfa","Col3a1","Foxs1","Mapk10" #Mesenchymal
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")



int.T1D_Week18 = RenameIdents(int.T1D_Week18, 
                              "0" = "T Cell",
                              "1" = "B Cell",
                              "2" = "T Cell",
                              "3" = "Macrophage",
                              "4" = "B Cell",
                              "5" = "T Cell",
                              "6" = "Monocyte Derived DC",
                              "7" = "Macrophage",
                              "8" = "T Cell",
                              "9" = "T Cell",
                              "10" = "B Cell",
                              "11" = "T Cell",
                              "12" = "Plasma Cell",
                              "13" = "T Cell",
                              "14" = "T Cell",
                              "15" = "T Cell",
                              "16" = "T Cell",
                              "17" = "T Cell",
                              "18" = "Acinar Cell",
                              "19" = "Dendritic Cell",
                              "20" = "T Cell"
)




DimPlot(int.T1D_Week18, reduction = "umap", 
        #split.by = c( "Tx"), 
        combine = F, label = T)

library(ggpubr)
DimPlot(int.T1D_Week18, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 7, repel = T) +
  theme_prism(base_size = 16, base_fontface = "bold") + 
  theme(legend.text = element_text(size = 18)) +
  labs( title = "UMAP of Cell Types") 
ggsave(file = "UMAP_CellTypeclusters_Week18.png",
       units = "in",width = 11, height = 11, dpi = 600)


saveRDS(int.T1D_Week18, file = "./annotated_int.T1D_Week18_v3.rds")
int.T1D_Week18<-readRDS("./annotated_int.T1D_Week18_v3.rds")
DimPlot(int.T1D_Week18)
# Extract the colors from the plot
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=7)
# Create frequency tables for Week 18 using active.ident and group
freq_w18<- prop.table(x = table(int.T1D_Week18@active.ident,int.T1D_Week18@meta.data[, "group"]), margin =2)
# Create a vector of colors for each bar
bar_colors <- rep(color_list, length.out = nrow(freq_w18))
# Set the resolution to 300 DPI
# Plot for Week 18 T Cells
par(mar = c(5, 6, 4, 6))  # Increasing the right margin to 6

# Plot the bar plot with modified x-axis label rotation
barplot(
  freq_w18, 
  xlim = c(0, ncol(freq_w18) + 3),
  col = bar_colors, 
  legend.text = TRUE, 
  args.legend = list(
    x = ncol(freq_w18) + 4.2, 
    y = max(colSums(freq_w18)), 
    bty = "n", 
    cex = 1.5
  ), 
  width = 1.25, 
  ylab = "Fraction of Cell Population",
  cex.lab = 2,   # Increase font size of axis titles
  cex.names = 2,  # Increase font size of x-axis labels
  cex.axis = 1.5     # Increase font size of axis tick labels
)

## ----T Cell subclustering---------------------------------------------------------------------------------------------------------------------------------------

int.T1D_Week18 = FindSubCluster(
  object = int.T1D_Week18,
  cluster = "T Cell",
  graph.name = "SCT_nn",
  resolution = 0.5,
  subcluster.name = "T Cell Subtypes",
  algorithm = 2
)

unique(T1D_Timepoints$`T Cell Subtypes`)

#Idents(T1D_Timepoints) = T1D_Timepoints$`T Cell Subtypes`

#DimPlot(T1D_Timepoints, reduction = "umap", group.by = , combine = F, label = T)



Tcells_w18 = subset(x = int.T1D_Week18, idents = "T Cell")
ElbowPlot(Tcells_w18, ndims = 50)

Tcells_w18 <- FindNeighbors(Tcells_w18, dims = 1:30, verbose = F,reduction = "integrated.dr")
Tcells_w18 <- FindClusters(Tcells_w18,graph.name = "SCT_nn", res = 0.7)#graph.name = "SCT_nn"
Tcells_w18 <- RunUMAP(Tcells_w18, dims = 1:30,reduction = "integrated.dr")

DimPlot(Tcells_w18, reduction = "umap", combine = F, label = T)
DimPlot(Tcells_w18, reduction = "umap", combine = F, label = T,group.by = ("sample"))

get_conserved <- function(cluster){
  Seurat::FindConservedMarkers(Tcells_w18,
                               ident.1 = cluster,
                               grouping.var = "sample",
                               only.pos = TRUE) %>%
    rownames_to_column(var = "gene")  %>%
    #left_join(y = unique(annotations[, c("gene_name", "description")]),
    #          by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
unique(Idents(Tcells_w18))

conserved_markers.Tcells_w18 <- map_dfr(0:10, get_conserved)
conserved_markers.Tcells_w18[is.na(conserved_markers.Tcells_w18)] <- 0
head(conserved_markers.Tcells_w18)

DefaultAssay(Tcells_w18) = "SCT"

# Extract top 10 markers per cluster
top20_TCells_w18 <- conserved_markers.Tcells_w18 %>% 
  mutate(avg_fc = (`Week18_1_avg_log2FC` + 
                     `Week18_2_avg_log2FC` + 
                     `Week18_3_avg_log2FC` +
                     `Week18_4_avg_log2FC` +
                     `Week18_5_avg_log2FC` +
                     `Week18_6_avg_log2FC`) /10) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)


Tcells_w18 <- PrepSCTFindMarkers(Tcells_w18)
all.markers <- FindAllMarkers(Tcells_w18 , only.pos = TRUE,recorrect_umi = FALSE)

top20_sct_Tcell_w18=all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)

# Save the top 20 markers grouped by cluster as a CSV file
write.csv(top20_sct_Tcell_w18, file = "top20_sct_Tcells_Week18.csv", row.names = FALSE)


DotPlot(Tcells_w18,assay = "SCT", features = c( "Tnfrsf4","Tbc1d4","Cd4", #CD4
                                            "Cd8a","Cd8b1","Blk",#CD8
                                            "Ccr7","Il7r","Lef1","Sell",#Memory
                                            "Gzma","Gzmb","Klrg1","Tbx21","Zeb2",#SLEC
                                            "Ccl5", "Gzmk", "Klrc1", "Nkg7", "Cd38", "Cxcr3", "Cxcr6", "Fasl", "Ifng",#CD8 Effector
                                            "Cd40lg","Il6ra",#Tcon effector
                                            "Nt5e","Izumo1r",#Anergy
                                            "Slamf6","Tcf7",#Progenitor
                                            "Cd200","Il21","Tnfsf8", #Th21
                                            "Rora","Stat3", "Ccr6",#Th17
                                            "Bcl6","Cxcr5","Tox2",#Tfh
                                            "Ctla4","Foxp3","Icos","Ikzf2","Il2ra","Tigit",#Tregs
                                            "Cd69","Egr1","Egr2","Nr4a1",#Activation
                                            "Ifit1","Ifit3","Isg15","Stat1",#Interferon sensing
                                            "Mki67","Stmn1",#Proliferation
                                            "Tox","Pdcd1","Lag3",#Exhaustion
                                            "Klrb1b","Klra8","Klrb1a","Eomes","Ncr1","Klre1", #NK
                                            "Tcrg-C2",#Gamm_Delta-,"Trac","Cd3e"
                                            "Hs3st1","Gata3",#ILC2
                                            "Hebp1","Ncoa7","Cd74","Cited4","H2-Ab1","Rorc"#ILC3
                                            
                                            
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

# Check average expression of Nkg7 per cluster
AverageExpression(Tcells_w18, features = "Cd8b1", group.by = "seurat_clusters")


# Extract top genes for cluster 4
cluster_4_genes <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 4) %>% 
  pull(gene)
# Print genes for cluster 0
cat("Top 20 genes for Cluster 4:\n")
print(cluster_4_genes)
cluster_4 <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 4)
cluster_4

# Extract top genes for cluster 5
cluster_5_genes <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 5) %>% 
  pull(gene)
# Print genes for cluster 0
cat("Top 20 genes for Cluster 5:\n")
print(cluster_5_genes)
cluster_5 <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 5)

# Extract top genes for cluster 6
cluster_6_genes <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 6) %>% 
  pull(gene)
# Print genes for cluster 6
cat("Top 20 genes for Cluster 6:\n")
print(cluster_6_genes)
cluster_6 <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 6)
cluster_6 

# Extract top genes for cluster 7
cluster_7_genes <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 7) %>% 
  pull(gene)
# Print genes for cluster 7
cat("Top 20 genes for Cluster 7:\n")
print(cluster_7_genes)
cluster_7 <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 7)
cluster_7 

# Extract top genes for cluster 7
cluster_7_genes <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 7) %>% 
  pull(gene)
# Print genes for cluster 7
cat("Top 20 genes for Cluster 7:\n")
print(cluster_7_genes)
cluster_7 <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 7)
cluster_7 

# Extract top genes for cluster 9
cluster_9_genes <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 9) %>% 
  pull(gene)
# Print genes for cluster 9
cat("Top 20 genes for Cluster 9:\n")
print(cluster_9_genes)
cluster_9 <- top20_sct_Tcell_w18  %>% 
  filter(cluster == 9)
cluster_9 

# Define the T cell subtypes for each cluster
TCelltype <- c("CD8 T Cell",  # 0
                  "CD4 T Cell",  #1 
                  "CD4 T Cell",  #2 
                  "CD8 T Cell",  #3 
                  "CD4 T Cell",  #4 
                  "NK Cell",  #5 
                  "CD4 T Cell",  #6 Tregs
                  "CD8 T Cell",  #7 
                  "CD4 T Cell",  #8 
                  "CD4 T Cell",  #9 TBD
                  "NK Cell"  #10 
)  

# Ensure the Idents are in numeric format
ident_integers <- as.integer(Idents(Tcells_w18))

# Check if there are any unexpected levels
unique(ident_integers)

# Add the T cell subtype information to the metadata of the Seurat object
Tcells_w18$TCellType <- TCelltype[as.integer(Idents(Tcells_w18))]

DimPlot(Tcells_w18,reduction = 'umap',group.by = "TCellType" ,label = T)


# Define the T cell subtypes for each cluster
TCellDetailedtype <- c("CD8 memory",  # 0
               "Tcon memory",  #1 
               "Tcon memory",  #2 
               "CD8 exhausted effector-like",  #3 
               "Th21-like",  #4 
               "NK Cell",  #5 
               "Tregs",  #6 Tregs
               "CD8 memory",  #7 
               "Tcon activated",  #8 
               "Tcon Interferon Sensing",  #9 TBD
               "NK Cell"  #10 
)  

# Ensure the Idents are in numeric format
ident_integers <- as.integer(Idents(Tcells_w18))

# Check if there are any unexpected levels
unique(ident_integers)

# Add the T cell subtype information to the metadata of the Seurat object
Tcells_w18$TCellSubType <- TCellDetailedtype[as.integer(Idents(Tcells_w18))]

DimPlot(Tcells_w18,reduction = 'umap',group.by = "TCellSubType" ,label = T)
Tcells_w18$Idents<-Idents(Tcells_w18)
Idents(Tcells_w18)<-Tcells_w18$TCellType
DimPlot(Tcells_w18,reduction = 'umap',label = T)
saveRDS(Tcells_w18, file = "./annotated_Tcells_w18_int_T1D_Week18_v1.rds")
Tcells_w18<-readRDS("./annotated_Tcells_w18_int_T1D_Week18_v1.rds")

Idents(Tcells_w18)<-Tcells_w18$TCellSubType
DimPlot(Tcells_w18)


# Extract the colors from the plot
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=8)
# Create frequency tables for Week 6 and Week 12 using active.ident and group
freq_Tcells_w18<- prop.table(x = table(Tcells_w18@active.ident,Tcells_w18@meta.data[, "group"]), margin =2)
# Create a vector of colors for each bar
bar_colors <- rep(color_list, length.out = nrow(freq_Tcells_w18))
# Set the resolution to 300 DPI
# Plot for Week 18 T Cells
par(mar = c(5, 6, 4, 6))  # Increasing the right margin to 6

# Plot the bar plot with modified x-axis label rotation
barplot(
  freq_Tcells_w18, 
  xlim = c(0, ncol(freq_Tcells_w18) + 3),
  col = bar_colors, 
  legend.text = TRUE, 
  args.legend = list(
    x = ncol(freq_Tcells_w18) + 4.2, 
    y = max(colSums(freq_Tcells_w18)), 
    bty = "n", 
    cex = 1.5
  ), 
  width = 1.25, 
  ylab = "Fraction of T Cell Population",
  cex.lab = 2,   # Increase font size of axis titles
  cex.names = 2,  # Increase font size of x-axis labels
  cex.axis = 1.5     # Increase font size of axis tick labels
)


dev.off()

## Rename T Cell Clusters in integrated.T1D_Week18 ----

setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/CellChat")
T1D_Week18<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_int.T1D_Week18_v3.rds")
T1D_TCells_Week18<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_Tcells_w18_int_T1D_Week18_v1.rds")
# Step 1: Ensure the cell names in T1D_TCells_Week18 exist in T1D_Week18
cells_to_update <- colnames(T1D_TCells_Week18)

# Step 2: Initialize CellSubType with the original cell identities
T1D_Week18$CellSubType <- Idents(T1D_Week18)

# Step 3: Visualize before updating
DimPlot(T1D_Week18, group.by = "CellSubType")

# Step 4: Ensure T1D_TCells_Week18$TCellSubType is a character vector
T1D_TCells_Week18$TCellSubType <- as.character(T1D_TCells_Week18$TCellSubType)

# Step 5: Expand levels in T1D_Week18$CellSubType to include new ones from TCellSubType
T1D_Week18$CellSubType <- factor(T1D_Week18$CellSubType, 
                                 levels = unique(c(levels(T1D_Week18$CellSubType), unique(T1D_TCells_Week18$TCellSubType))))

# Step 6: Assign T cell subtypes without generating invalid factor levels
T1D_Week18$CellSubType[cells_to_update] <- T1D_TCells_Week18$TCellSubType

# Step 7: Verify the update
levels(T1D_Week18$CellSubType)

# Step 8: Visualize after updating
DimPlot(T1D_Week18, group.by = "CellSubType")

saveRDS(T1D_Week18,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_int.T1D_Week18_v4.rds")

# ----Macrophage subclustering---------------------------------------------------------------------------------------------------------------------------------------

int.T1D_Week18 = FindSubCluster(
  object = int.T1D_Week18,
  cluster = "Macrophage",
  graph.name = "SCT_nn",
  resolution = 0.5,
  subcluster.name = "Macrophage Subtypes",
  algorithm = 2
)

unique(int.T1D_Week18$`Macrophage Subtypes`)

#Idents(T1D_Timepoints) = T1D_Timepoints$`T Cell Subtypes`

#DimPlot(T1D_Timepoints, reduction = "umap", group.by = , combine = F, label = T)



Mac_w18 = subset(x = int.T1D_Week18, idents = "Macrophage")
ElbowPlot(Mac_w18, ndims = 50)

Mac_w18 <- FindNeighbors(Mac_w18, dims = 1:30, verbose = F,reduction = "integrated.dr")
Mac_w18 <- FindClusters(Mac_w18,graph.name = "SCT_nn", res = 0.3)#graph.name = "SCT_nn"
Mac_w18 <- RunUMAP(Mac_w18, dims = 1:30,reduction = "integrated.dr")

DimPlot(Mac_w18, reduction = "umap", combine = F, label = T)
DimPlot(Mac_w18, reduction = "umap", combine = F, label = T,group.by = ("sample"))

get_conserved <- function(cluster){
  Seurat::FindConservedMarkers(Mac_w18,
                               ident.1 = cluster,
                               grouping.var = "sample",
                               only.pos = TRUE) %>%
    rownames_to_column(var = "gene")  %>%
    #left_join(y = unique(annotations[, c("gene_name", "description")]),
    #          by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
unique(Idents(Mac_w18))
DefaultAssay(Mac_w18) = "RNA"
conserved_markers.Mac_w18 <- map_dfr(0:4, get_conserved)
conserved_markers.Mac_w18[is.na(conserved_markers.Mac_w18)] <- 0
head(conserved_markers.Mac_w18)



# Extract top 10 markers per cluster
top20_Mac_w18 <- conserved_markers.Mac_w18 %>% 
  mutate(avg_fc = (`Week18_1_avg_log2FC` + 
                     `Week18_2_avg_log2FC` + 
                     `Week18_3_avg_log2FC` +
                     `Week18_4_avg_log2FC` +
                     `Week18_5_avg_log2FC` +
                     `Week18_6_avg_log2FC`) /10) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)

DefaultAssay(Mac_w18) = "SCT"
Mac_w18 <- PrepSCTFindMarkers(Mac_w18)
all.markers <- FindAllMarkers(Mac_w18 , only.pos = TRUE,recorrect_umi = FALSE)

top20_sct_Mac_w18=all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)

# Save the top 20 markers grouped by cluster as a CSV file
write.csv(top20_sct_Mac_w18, file = "top20_sct_Mac_Week18.csv", row.names = FALSE)


# Extract top genes for cluster 0
cluster_0_genes <- top20_sct_Mac_w18  %>% 
  filter(cluster == 0) %>% 
  pull(gene)
# Print genes for cluster 0
cat("Top 20 genes for Cluster 0:\n")
print(cluster_0_genes)
cluster_0 <- top20_sct_Mac_w18  %>% 
  filter(cluster == 0)
cluster_0

# Extract top genes for cluster 1
cluster_1_genes <- top20_sct_Mac_w18  %>% 
  filter(cluster == 1) %>% 
  pull(gene)
# Print genes for cluster 1
cat("Top 20 genes for Cluster 1:\n")
print(cluster_1_genes)
cluster_1 <- top20_sct_Mac_w18  %>% 
  filter(cluster == 1)

# Extract top genes for cluster 2
cluster_2_genes <- top20_sct_Mac_w18  %>% 
  filter(cluster == 2) %>% 
  pull(gene)
# Print genes for cluster 2
cat("Top 20 genes for Cluster 2:\n")
print(cluster_2_genes)
cluster_2 <- top20_sct_Mac_w18  %>% 
  filter(cluster == 2)
cluster_2 

# Extract top genes for cluster 3
cluster_3_genes <- top20_sct_Mac_w18  %>% 
  filter(cluster == 3) %>% 
  pull(gene)
# Print genes for cluster 3
cat("Top 20 genes for Cluster 3:\n")
print(cluster_3_genes)
cluster_3 <- top20_sct_Mac_w18  %>% 
  filter(cluster == 3)
cluster_3


# Extract top genes for cluster 4
cluster_4_genes <- top20_sct_Mac_w18  %>% 
  filter(cluster == 4) %>% 
  pull(gene)
# Print genes for cluster 4
cat("Top 20 genes for Cluster 4:\n")
print(cluster_4_genes)
cluster_4 <- top20_sct_Mac_w18  %>% 
  filter(cluster == 4)
cluster_4


DotPlot(Mac_w18,assay = "SCT",group.by = "MacType", features = c(  "Jun","Atf3","Fos","Ccl3","Egr1","Junb", #NF-KB 
                                             "Apoe","Trem2",#Not Activated
                                             "Cxcl9","Stat1","Ccl5","Tapbp","Tap1","Cd40","Il12b","Tap2","Psmb8",#IFN-Y and NF-KB activated
                                             "Prdx1","Il1rn","Lgals1","Lgals3","Cd36","Anxa1","Anxa4","Anxa5","Anxa2","Cxcl14",#Anti-Inflammatory, cell death
                                             "Birc5","Stmn1","Cdca3","Mki67","Ccna2","Cdk1","Cdca8"#Cell-cycle
                                                
                                                
                                                
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

# Check average expression of Nkg7 per cluster
AverageExpression(Mac_w18, features = "Cd8b1", group.by = "seurat_clusters")


# Define the T cell subtypes for each cluster
Mactype <- c("Mac-1",  # 0
               "Mac-2",  #1 -NFKB Activated
               "Mac-3",  #2 
               "Mac-4",  #3 
               "Mac-2"  #4 
)  

# Ensure the Idents are in numeric format
ident_integers <- as.integer(Idents(Mac_w18))

# Check if there are any unexpected levels
unique(ident_integers)

# Add the T cell subtype information to the metadata of the Seurat object
Mac_w18$MacType <- Mactype[as.integer(Idents(Mac_w18))]

DimPlot(Mac_w18,reduction = 'umap',group.by = "MacType" ,label = T)

DefaultAssay(Mac_w18) = "SCT"
Mac_w18 <- PrepSCTFindMarkers(Mac_w18)
all.markers <- FindAllMarkers(Mac_w18 , only.pos = TRUE,recorrect_umi = FALSE,group.by = "MacType")

top20_sct_Mac_w18=all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)

# Save the top 20 markers grouped by cluster as a CSV file
write.csv(top20_sct_Mac_w18, file = "top20_sct_Mac_Week18_anno.csv", row.names = FALSE)


# Extract top genes for cluster Mac-1
cluster_Mac1_genes <- top20_sct_Mac_w18  %>% 
  filter(cluster == 'Mac-1') %>% 
  pull(gene)
# Print genes for cluster Mac-1
cat("Top 20 genes for Cluster Mac-1:\n")
print(cluster_Mac1_genes)
cluster_Mac1 <- top20_sct_Mac_w18  %>% 
  filter(cluster == 'Mac-1')
cluster_0

# Extract top genes for cluster Mac-2
cluster_Mac2_genes <- top20_sct_Mac_w18  %>% 
  filter(cluster == "Mac-2") %>% 
  pull(gene)
# Print genes for cluster Mac-2
cat("Top 20 genes for Cluster Mac-2:\n")
print(cluster_Mac2_genes)
cluster_Mac2 <- top20_sct_Mac_w18  %>% 
  filter(cluster == "Mac-2")

# Extract top genes for cluster Mac-3
cluster_Mac3_genes <- top20_sct_Mac_w18  %>% 
  filter(cluster == "Mac-3") %>% 
  pull(gene)
# Print genes for cluster Mac-3
cat("Top 20 genes for Cluster Mac-3:\n")
print(cluster_Mac3_genes)
cluster_Mac3 <- top20_sct_Mac_w18  %>% 
  filter(cluster == "Mac-3")
cluster_Mac3 

# Extract top genes for cluster Mac-4
cluster_Mac4_genes <- top20_sct_Mac_w18  %>% 
  filter(cluster == Mac4) %>% 
  pull(gene)
# Print genes for cluster 3
cat("Top 20 genes for Cluster Mac-4:\n")
print(cluster_Mac4_genes)
cluster_Mac4 <- top20_sct_Mac_w18  %>% 
  filter(cluster == "Mac-4")
cluster_Mac4


# Extract top genes for cluster 4
cluster_4_genes <- top20_sct_Mac_w18  %>% 
  filter(cluster == 4) %>% 
  pull(gene)
# Print genes for cluster 4
cat("Top 20 genes for Cluster 4:\n")
print(cluster_4_genes)
cluster_4 <- top20_sct_Mac_w18  %>% 
  filter(cluster == 4)
cluster_4



Mac_w18$Idents<-Idents(Mac_w18)

Idents(Mac_w18)<-Mac_w18$MacType
DimPlot(Mac_w18,reduction = 'umap',label = T)


color_list <- ggplotColours(n=4)
# Create frequency tables for Mac Week 18
freq_Mac_w18<- prop.table(x = table(Mac_w18@active.ident,Mac_w18@meta.data[, "group"]), margin =2)
# Create a vector of colors for each bar
bar_colors <- rep(color_list, length.out = nrow(freq_Mac_w18))
# Set the resolution to 300 DPI
# Plot for Week 18 Mac
par(mar = c(5, 6, 4, 6))  # Increasing the right margin to 6

# Plot the bar plot with modified x-axis label rotation
barplot(
  freq_Mac_w18, 
  xlim = c(0, ncol(freq_Mac_w18) + 3),
  col = bar_colors, 
  legend.text = TRUE, 
  args.legend = list(
    x = ncol(freq_Mac_w18) + 2, 
    y = max(colSums(freq_Mac_w18)), 
    bty = "n", 
    cex = 1.5
  ), 
  width = 1.25, 
  ylab = "Fraction of Macrophage Population",
  cex.lab = 2,   # Increase font size of axis titles
  cex.names = 2,  # Increase font size of x-axis labels
  cex.axis = 1.5     # Increase font size of axis tick labels
)
saveRDS(Mac_w18, file = "./annotated_Mac_w18_int_T1D_Week18_v1.rds")
