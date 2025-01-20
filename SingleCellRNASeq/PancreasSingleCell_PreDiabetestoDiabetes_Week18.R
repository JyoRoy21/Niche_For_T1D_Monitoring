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

