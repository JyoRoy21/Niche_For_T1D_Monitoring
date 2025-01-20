#' ---
#' title: "R Notebook"
#' author: "Jyotirmoy Roy adapted from Kate Griffin"
#' Date: "`r Sys.Date()`"
#' output:
#'   pdf_document: default
#'   html_notebook: default
#' editor_options: 
#'   markdown: 
#'     wrap: 72
#' ---
#' 
#' This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you
#' execute code within the notebook, the results appear beneath the code.
#' 
#' Try executing this chunk by clicking the *Run* button within the chunk
#' or by placing your cursor inside it and pressing *Cmd+Shift+Enter*.
#' 
## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' Add a new chunk by clicking the *Insert Chunk* button on the toolbar or
#' by pressing *Cmd+Option+I*.
#' 
#' When you save the notebook, an HTML file containing the code and output
#' will be saved alongside it (click the *Preview* button or press
#' *Cmd+Shift+K* to preview the HTML file).
#' 
#' The preview shows you a rendered HTML copy of the contents of the
#' editor. Consequently, unlike *Knit*, *Preview* does not run any R code
#' chunks. Instead, the output of the chunk when it was last run in the
#' editor is displayed.
#' 
#' # Load libraries
#' 
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
## ----source files-------------------------------------------------------------------------------------------------------------------------------------
source("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/SingleCellRNASeq/remove_soup.R")
source("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/SingleCellRNASeq/find_doublets.R")
source("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/SingleCellRNASeq/run_ddqcR.R")

#' 
#' # Read in 10X data
#' 
#' Some meta data info:
#' 
#' Sample2 "PBS PreCh1" Sample3 "NP Ch4" Sample4 "PBS Ch4" Sample5 "NP Ch7"
#' Sample6 "PBS Ch7"
#' 
#' 
#' 
## ----read in 10x data---------------------------------------------------------------------------------------------------------------------------------
data_dir = "11329_JR_pool/"

sample.names = list.files(path = '11329_JR_pool/')
sample.names

folder.names = paste0(sample.names, "/")

sample_dir_list = folder.names



#' 
#' # SoupX
#' 
#' run SoupX from function file "CellRanger_SoupX_Processing.R" params:
#' sample, data_dir, sample_dir
#' 
## ----run soupX----------------------------------------------------------------------------------------------------------------------------------------
#Install glmGamPoi before running soupX
#BiocManager::install('glmGamPoi')
#library("glmGamPoi")

#rm(sample_list)

#in loop format
#for (x in 1:5) {
#  sample_list[x] = remove_soup(sample_list[x], data_dir, sample_dir_list[x])
#}

Week18_1 = remove_soup(Week18_1, data_dir, sample_dir_list[1])
Week18_2 =  remove_soup(Week18_2, data_dir, sample_dir_list[9])
Week18_3 =  remove_soup(Week18_3, data_dir, sample_dir_list[10])
Week18_4 =  remove_soup(NWeek18_4, data_dir, sample_dir_list[11])
Week18_5 =  remove_soup(Week18_5, data_dir, sample_dir_list[12])
Week18_6 =  remove_soup(Week18_6, data_dir, sample_dir_list[13])

Week12_1 = remove_soup(Week12_1, data_dir, sample_dir_list[14])
Week12_2 =  remove_soup(Week12_2, data_dir, sample_dir_list[15])
Week12_3 =  remove_soup(Week12_3, data_dir, sample_dir_list[16])
Week12_4 =  remove_soup(NWeek12_4, data_dir, sample_dir_list[2])

Week6_1 = remove_soup(Week6_1, data_dir, sample_dir_list[3])
Week6_2 =  remove_soup(Week6_2, data_dir, sample_dir_list[4])
Week6_3 =  remove_soup(Week6_3, data_dir, sample_dir_list[5])
Week6_4 =  remove_soup(NWeek6_4, data_dir, sample_dir_list[6])
Week6_5 =  remove_soup(Week6_5, data_dir, sample_dir_list[7])
Week6_6 =  remove_soup(Week6_6, data_dir, sample_dir_list[8])

#save.image(file = "Allergy_Katepiepeline_v1.RData")


#' 
#' # ddqcR
#' 
#' run ddqcR with "filter_ddqcr" function params: seurat object, "sample
#' name", data_dir, sample_dir
#' 
## ----run ddqcR----------------------------------------------------------------------------------------------------------------------------------------

sample_list<- c("Week18_1", "Week18_2", "Week18_3", "Week18_4", "Week18_5","Week18_6",
                "Week12_1","Week12_2","Week12_3","Week12_4",
                "Week6_1","Week6_2","Week6_3","Week6_4","Week6_5","Week6_6")

Week18_1 = filter_ddqcr(Week18_1, sample_list[1], sample_dir_list[1])
Week18_2 =  filter_ddqcr(Week18_2, sample_list[2], sample_dir_list[9])
Week18_3 =  filter_ddqcr(Week18_3, sample_list[3], sample_dir_list[10])
Week18_4 =  filter_ddqcr(Week18_4, sample_list[4], sample_dir_list[11])
Week18_5 =  filter_ddqcr(Week18_5, sample_list[5], sample_dir_list[12])
Week18_6 =  filter_ddqcr(Week18_6, sample_list[6], sample_dir_list[13])

Week12_1 = filter_ddqcr(Week12_1, sample_list[7], sample_dir_list[14])
Week12_2 =  filter_ddqcr(Week12_2, sample_list[8], sample_dir_list[15])
Week12_3 =  filter_ddqcr(Week12_3, sample_list[9], sample_dir_list[16])
Week12_4 =  filter_ddqcr(Week12_4, sample_list[10], sample_dir_list[2])

Week6_1 = filter_ddqcr(Week6_1, sample_list, sample_dir_list[3])
Week6_2 =  filter_ddqcr(Week6_2, sample_list, sample_dir_list[4])
Week6_3 =  filter_ddqcr(Week6_3, sample_list, sample_dir_list[5])
Week6_4 =  filter_ddqcr(Week6_4, sample_list, sample_dir_list[6])
Week6_5 =  filter_ddqcr(Week6_5, sample_list, sample_dir_list[7])
Week6_6 =  filter_ddqcr(Week6_6, sample_list, sample_dir_list[8])



#' 
#' # Run Doublet Finder on Samples
#' 
#' run doublet finder with "find_doublets" function params: seu_object,
#' "sample name", data_dir, sample_dir
#' 
#' 
## ----Run Doublet Finder on Samples--------------------------------------------------------------------------------------------------------------------

Week18_1 = find_doublets(Week18_1, sample_list[1], sample_dir_list[1])
Week18_2 =  find_doublets(Week18_2, sample_list[2], sample_dir_list[9])
Week18_3 =  find_doublets(Week18_3, sample_list[3], sample_dir_list[10])
Week18_4 =  find_doublets(Week18_4, sample_list[4], sample_dir_list[11])
Week18_5 =  find_doublets(Week18_5, sample_list[5], sample_dir_list[12])
Week18_6 =  find_doublets(Week18_6, sample_list[6], sample_dir_list[13])

Week12_1 = find_doublets(Week12_1, sample_list[7], sample_dir_list[14])
Week12_2 =  find_doublets(Week12_2, sample_list[8], sample_dir_list[15])
Week12_3 =  find_doublets(Week12_3, sample_list[9], sample_dir_list[16])
Week12_4 =  find_doublets(Week12_4, sample_list[10], sample_dir_list[2])

Week6_1 = find_doublets(Week6_1, sample_list, sample_dir_list[3])
Week6_2 =  find_doublets(Week6_2, sample_list, sample_dir_list[4])
Week6_3 =  find_doublets(Week6_3, sample_list, sample_dir_list[5])
Week6_4 =  find_doublets(Week6_4, sample_list, sample_dir_list[6])
Week6_5 =  find_doublets(Week6_5, sample_list, sample_dir_list[7])
Week6_6 =  find_doublets(Week6_6, sample_list, sample_dir_list[8])




#' 
## ----save filtered seurat objects---------------------------------------------------------------------------------------------------------------------
saveRDS(Week18_1, file = "./filtered_Week18_1.rds")
saveRDS(Week18_2, file = "./filtered_Week18_2.rds")
saveRDS(Week18_3, file = "./filtered_Week18_3.rds")
saveRDS(Week18_4, file = "./filtered_Week18_4.rds")
saveRDS(Week18_5, file = "./filtered_Week18_5.rds")
saveRDS(Week18_6, file = "./filtered_Week18_6.rds")

saveRDS(Week12_1, file = "./filtered_Week12_1.rds")
saveRDS(Week12_2, file = "./filtered_Week12_2.rds")
saveRDS(Week12_3, file = "./filtered_Week12_3.rds")
saveRDS(Week12_4, file = "./filtered_Week12_4.rds")

saveRDS(Week6_1, file = "./filtered_Week6_1.rds")
saveRDS(Week6_2, file = "./filtered_Week6_2.rds")
saveRDS(Week6_3, file = "./filtered_Week6_3.rds")
saveRDS(Week6_4, file = "./filtered_Week6_4.rds")
saveRDS(Week6_5, file = "./filtered_Week6_5.rds")
saveRDS(Week6_6, file = "./filtered_Week6_6.rds")








#' 
## ----add to metadata of each seurat object------------------------------------------------------------------------------------------------------------
#Assign an identity to the cells - this is only necessary when you are combining more than one file/Seurat object
## Sample Name
#Week 18
Week18_1@meta.data$sample<-"Week18_1"
Week18_2@meta.data$sample<-"Week18_2"
Week18_3@meta.data$sample<-"Week18_3"
Week18_4@meta.data$sample<-"Week18_4"
Week18_5@meta.data$sample<-"Week18_5"
Week18_6@meta.data$sample<-"Week18_6"
#Week 12
Week12_1@meta.data$sample<-"Week12_1"
Week12_2@meta.data$sample<-"Week12_2"
Week12_3@meta.data$sample<-"Week12_3"
Week12_4@meta.data$sample<-"Week12_4"
#Week 6
Week6_1@meta.data$sample<-"Week6_1"
Week6_2@meta.data$sample<-"Week6_2"
Week6_3@meta.data$sample<-"Week6_3"
Week6_4@meta.data$sample<-"Week6_4"
Week6_5@meta.data$sample<-"Week6_5"
Week6_6@meta.data$sample<-"Week6_6"

## Time
#Week 18
Week18_1@meta.data$time<-"Week18"
Week18_2@meta.data$time<-"Week18"
Week18_3@meta.data$time<-"Week18"
Week18_4@meta.data$time<-"Week18"
Week18_5@meta.data$time<-"Week18"
Week18_6@meta.data$time<-"Week18"
#Week 12
Week12_1@meta.data$time<-"Week12"
Week12_2@meta.data$time<-"Week12"
Week12_3@meta.data$time<-"Week12"
Week12_4@meta.data$time<-"Week12"
#Week 6
Week6_1@meta.data$time<-"Week6"
Week6_2@meta.data$time<-"Week6"
Week6_3@meta.data$time<-"Week6"
Week6_4@meta.data$time<-"Week6"
Week6_5@meta.data$time<-"Week6"
Week6_6@meta.data$time<-"Week6"

## Group
#Week 1
Week18_1@meta.data$group<-"Pre-Progressor"
Week18_2@meta.data$group<-"Pre-Progressor"
Week18_3@meta.data$group<-"Pre-Progressor"
Week18_4@meta.data$group<-"Progressor"
Week18_5@meta.data$group<-"Progressor"
Week18_6@meta.data$group<-"Progressor"
#Week 12
Week12_1@meta.data$group<-"Progressor"
Week12_2@meta.data$group<-"Progressor"
Week12_3@meta.data$group<-"Non-Progressor"
Week12_4@meta.data$group<-"Non-Progressor"
#Week 6
Week6_1@meta.data$group<-"Non-Progressor"
Week6_2@meta.data$group<-"Progressor"
Week6_3@meta.data$group<-"Progressor"
Week6_4@meta.data$group<-"Progressor"
Week6_5@meta.data$group<-"Progressor"
Week6_6@meta.data$group<-"Progressor"

#' 
#' # Merge THEN scTransform
#' 
#' Do not run if using scTransform method
#' 
## ----Merge THEN scTransform method--------------------------------------------------------------------------------------------------------------------
#Create Seurat objects
# use count matrix to create seurat object, which has data (ie count matrix)
# and analysis (ie PCA/clustering)
# initialized with raw data

Week18_1 <- NormalizeData(Week18_1) # normalize the data with log norm
Week18_1 <- ScaleData(Week18_1, verbose = F) # linear transformation prior to
# dimensional reduction like PCA. Shifts expression so mean exp is 0 across all cells,
# scales exp, so var across cells is 1
# verbose calls progress bar
Week18_2<- NormalizeData(Week18_2)
Week18_2 <- ScaleData(Week18_2, verbose = F)
Week18_3<- NormalizeData(Week18_3)
Week18_3 <- ScaleData(Week18_3, verbose = F)
Week18_4<- NormalizeData(Week18_4)
Week18_4 <- ScaleData(Week18_4, verbose = F)
Week18_5<- NormalizeData(Week18_5)
Week18_5 <- ScaleData(Week18_5, verbose = F)
Week18_6<- NormalizeData(Week18_6)
Week18_6 <- ScaleData(Week18_6, verbose = F)

Week12_1 <- NormalizeData(Week12_1) 
Week12_1 <- ScaleData(Week12_1, verbose = F) 
Week12_2<- NormalizeData(Week12_2)
Week12_2 <- ScaleData(Week12_2, verbose = F)
Week12_3<- NormalizeData(Week12_3)
Week12_3 <- ScaleData(Week12_3, verbose = F)
Week12_4<- NormalizeData(Week12_4)
Week12_4 <- ScaleData(Week12_4, verbose = F)


Week6_1 <- NormalizeData(Week6_1) # normalize the data with log norm
Week6_1 <- ScaleData(Week6_1, verbose = F) 
Week6_2<- NormalizeData(Week6_2)
Week6_2 <- ScaleData(Week6_2, verbose = F)
Week6_3<- NormalizeData(Week6_3)
Week6_3 <- ScaleData(Week6_3, verbose = F)
Week6_4<- NormalizeData(Week6_4)
Week6_4 <- ScaleData(Week6_4, verbose = F)
Week6_5<- NormalizeData(Week6_5)
Week6_5 <- ScaleData(Week6_5, verbose = F)
Week6_6<- NormalizeData(Week6_6)
Week6_6 <- ScaleData(Week6_6, verbose = F)


#' 
#' ## merge
#' 
#' merge tutorial: <https://satijalab.org/seurat/archive/v4.3/merge>
#' 
## ----merge Week 6 and Week 12 timepoints--------------------------------------------------------------------------------------------------------------------------------------------
T1D_Timepoints = merge(x = Week6_1, y = list( Week6_2, Week6_3,Week6_4, Week6_5,Week6_6,
                                              Week12_1,Week12_2,Week12_3,Week12_4))


#' 
#' ## scTransform
#' 
#' scTransform vignette:
#' <https://satijalab.org/seurat/articles/sctransform_vignette.html>
#' install glmGamPoi before using, significantly improves speed
#' BiocManager::install("glmGamPoi")
#' 
## ----scTransform on Week 6 and Week 12 timepoints --------------------------------------------------------------------------------------------------------------------------------------
T1D_Timepoints  <- SCTransform(T1D_Timepoints, vars.to.regress = "percent.mt", verbose=F)
saveRDS(T1D_Timepoints , file = "./T1D_Timepoints_sct_filtered.rds")


## ----merge Week 18 timepoints--------------------------------------------------------------------------------------------------------------------------------------------
T1D_Week18 = merge(x = Week18_1, y = list( Week18_2, Week18_3,Week18_4, Week18_5,Week18_6))


#' 
#' ## scTransform
#' 
#' scTransform vignette:
#' <https://satijalab.org/seurat/articles/sctransform_vignette.html>
#' install glmGamPoi before using, significantly improves speed
#' BiocManager::install("glmGamPoi")
#' 
## ----scTransform--------------------------------------------------------------------------------------------------------------------------------------
T1D_Week18  <- SCTransform(T1D_Week18 , vars.to.regress = "percent.mt", verbose=T)
saveRDS(T1D_Week18  , file = "./T1D_Week18_sct_filtered.rds")


#' 
## ----save merged and cluster unintegrated-------------------------------------------------------------------------------------------------------------
T1D_Timepoints = readRDS("./T1D_Timepoints_sct_filtered.rds")
#T1D_Timepoints <- SCTransform(T1D_Timepoints)
T1D_Timepoints <- RunPCA(T1D_Timepoints)
T1D_Timepoints <- RunUMAP(T1D_Timepoints, dims = 1:30)
T1D_Timepoints[[]]
DimPlot(T1D_Timepoints, reduction = "umap", group.by = c("sample", "seurat_clusters"))
DimPlot(T1D_Timepoints, reduction = "umap", split.by = c("sample"))

ElbowPlot(T1D_Timepoints, ndims = 50)
T1D_Timepoints <- FindNeighbors(T1D_Timepoints, dims = 1:30, verbose = F)
T1D_Timepoints <- FindClusters(T1D_Timepoints, verbose =F, resolution = 0.4)
DimPlot(T1D_Timepoints, group.by = c("time"), label=F)

DimPlot(T1D_Timepoints, reduction = "umap", group.by = c("sample"))
DimPlot(int.T1D_Timepoints, reduction = "umap", group.by = c("sample"))

table(Idents(T1D_Timepoints),T1D_Timepoints@meta.data$sample)Dim
#' 
#' # Integration
#' Seurat integration vignette: https://satijalab.org/seurat/articles/integration_introduction.html
#' Use integration method for SCtransform, not the regular one 
## ----integration clustering---------------------------------------------------------------------------------------------------------------------------
int.T1D_Timepoints <- IntegrateLayers(object = T1D_Timepoints, method = CCAIntegration, normalization.method = "SCT", verbose =T)
int.T1D_Timepoints <- FindNeighbors(int.T1D_Timepoints, reduction = "integrated.dr", dims = 1:30)
int.T1D_Timepoints <- FindClusters(int.T1D_Timepoints, res = 0.6)

#' 
## ----cluster integrated DimPlot-----------------------------------------------------------------------------------------------------------------------
int.T1D_Timepoints <- RunUMAP(int.T1D_Timepoints, dims = 1:30, reduction = "integrated.dr")
DimPlot(int.T1D_Timepoints, reduction = "umap", group.by = c( "seurat_clusters"), split.by = c( "time"), combine = F)
DimPlot(int.T1D_Timepoints, reduction = "umap", group.by = c( "SCT_snn_res.0.3"), split.by = c( "time"), combine = F)
int.T1D_Timepoints$SCT_snn_res.0.5

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#png(file = "UMAP_31clusters.png",    units = "in",width = 11, height = 11, res = 400)
# DimPlot(int.T1D_Timepoints, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T, 
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
DimPlot(int.T1D_Timepoints, reduction = "umap", group.by = c("sample", "seurat_clusters"))
DimPlot(int.T1D_Timepoints, reduction = "umap", split.by = c("sample"))
DimPlot(int.T1D_Timepoints, reduction = "umap", group.by = c( "seurat_clusters"), split.by = c( "time"), combine = F)
DimPlot(T1D_Timepoints, reduction = "umap", group.by = c( "seurat_clusters"), split.by = c( "time"), combine = F)

int.T1D_Timepoints$SCT_snn_res.0.3



#' 
## ----table of integrated clusters---------------------------------------------------------------------------------------------------------------------
table(Idents(int.T1D_Timepoints),int.T1D_Timepoints@meta.data$sample)

#' 
#' 
#' Save integrated object
## ----save integrated object---------------------------------------------------------------------------------------------------------------------------
saveRDS(int.T1D_Timepoints, file = "./sct_integrated_T1D_Timepoints_v2.rds")
getwd()

#' 
#' 
## ----save workspace-----------------------------------------------------------------------------------------------------------------------------------
save.image("./JR_KGpipeline_integrated.T1D_Timepoints.RData")

#' 
## ----rejoin layers and find conserved markers---------------------------------------------------------------------------------------------------------
# now that integration is complete, rejoin layers
Layers(T1D_Timepoints)
T1D_Timepoints
DefaultAssay(T1D_Timepoints) = "RNA"
T1D_Timepoints <- JoinLayers(T1D_Timepoints)




get_conserved <- function(cluster){
  Seurat::FindConservedMarkers(T1D_Timepoints,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene")  %>%
    #left_join(y = unique(annotations[, c("gene_name", "description")]),
    #          by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
#devtools::install_github('immunogenomics/presto')
FindConservedMarkers(T1D_Timepoints,
                     ident.1 = 0,
                     grouping.var = "sample",
                     only.pos = TRUE)

#conserved_markers <- map_dfr(unique(new.cluster.ids), get_conserved)
install.packages("tibble")
library(purrr)
library(tibble)
DimPlot(T1D_Timepoints,reduction = "umap")
conserved_markers <- map_dfr(0:22, get_conserved)

DefaultAssay(T1D_Timepoints) = "SCT"
head(conserved_markers)
conserved_markers$Week6_2_avg_log2FC
# Extract top 10 markers per cluster
top20 <- conserved_markers %>% 
  mutate(avg_fc = (`Week6_1_avg_log2FC` + 
                     `Week6_2_avg_log2FC` + 
                     `Week6_3_avg_log2FC` +
                     `Week6_4_avg_log2FC` +
                     `Week6_5_avg_log2FC` +
                     `Week6_6_avg_log2FC` +
                     `Week12_1_avg_log2FC` + 
                     `Week12_2_avg_log2FC` +
                     `Week12_3_avg_log2FC` +
                     `Week12_4_avg_log2FC` ) /10) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)




## ----cluster ID---------------------------------------------------------------------------------------------------------------------------------------



T1D_Timepoints <- PrepSCTFindMarkers(T1D_Timepoints)
all.markers <- FindAllMarkers(T1D_Timepoints , only.pos = TRUE)

top20_sct=all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)

# Plot the heatmap for the top 20 markers
DoHeatmap(
  object = T1D_Timepoints, 
  features = top20_sct$gene,  # Use the top 20 marker genes
  size = 4,                   # Adjust text size for readability
  group.by = "ident",         # Group cells by cluster identity
  label = TRUE                # Show cluster labels
) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) # Adjust color scale







#' 
#' 
#' "Cd33",
# "Cd3e", "Cd4", "Cd8a","Foxp3","Ctla4", 
# "Ncr1", "Klrb1c", "Klre1", 
# "Cd79a", "Cd19" ,"Ms4a1", "Ighm", 
# "Cd300a", "Ly6c2", "Cd14", "Ccr2", 
# "Smox", "Apoe", "Ms4a7", "Ccr5", "Fcgr1", "Adgre1", "Cd68",
# "Ftl1", "Fth1", 
# "Tmem119", "Trem2", "Aif1", "Csf1r",
# "S100a9", "Camp", "Retnlg", "Ly6g",
# "Ifitm1", "Siglech", "Klk1", 
# "Zbtb46", "Itgax", "Cd86", "Cd209a",
# "Bst2", "Cmah", "Ly6a", "Nrp1", "Clec4g", 
# "Cacnb3", "Fscn1", "Syn3", "Tmem150c", 
# "Snca", "Ube2o"
#"Ncr1",  "Klre1",#NK
# Extract top genes for cluster 18
cluster_18_genes <- top20_sct %>% 
  filter(cluster == 18) %>% 
  pull(gene)

# Extract top genes for cluster 19
cluster_19_genes <- top20_sct %>% 
  filter(cluster == 19) %>% 
  pull(gene)

# Print genes for cluster 18
cat("Top 20 genes for Cluster 18:\n")
print(cluster_18_genes)

# Print genes for cluster 19
cat("\nTop 20 genes for Cluster 19:\n")
print(cluster_19_genes)



DotPlot(T1D_Timepoints, features = c( "Ptprc", #CD45
                                      "Trbc1","Trbc2","Cd3g","Cd3e","Cd3d","Lat", #T Cells
                                      "C1qa","C1qb","C1qc","Lyz2","Adgre1", #Macrophages
                                      "Cd79a", "Cd79b","Ms4a1","Cd19",#B Cells
                                      "S100a4","S100a11","Itgax","Zbtb46","Flt3","H2-Aa","Cd74","Mdh2","Gm2a","Pglyrp1","H2-Oa",#cDc
                                      "Siglech","Cd209a","Cox6a2","Rnase6",#pDC
                                      "Ighg2b","Ighg2c","Igha","Igkc","Jchain",#Plasma
                                      "Ctrb1","Cela3b","Cpa1","Prss1",#Acinar
                                      "Fam184b","Chad","Kcnk1","Pdgfrb"#Mesenchymal Progenitors-Lingo4
                                     # "Il34","Pdgfrb","Pdgfa","Col3a1","Foxs1","Mapk10" #Mesenchymal
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")



T1D_Timepoints = RenameIdents(T1D_Timepoints, 
                                  "0" = "B Cell",
                                  "1" = "T Cell",
                                  "2" = "T Cell",
                                  "3" = "T Cell",
                                  "4" = "T Cell",
                                  "5" = "T Cell",
                                  "6" = "T Cell",
                                  "7" = "B Cell",
                                  "8" = "B Cell",
                                  "9" = "T Cell",
                                  "10" = "T Cell",
                                  "11" = "Macrophage",
                                  "12" = "T Cell",
                                  "13" = "T Cell",
                                  "14" = "T Cell",
                                  "15" = "T Cell",
                                  "16" = "T Cell",
                                  "17" = "Acinar Cell",
                                  "18" = "Mesenchymal-like", #TBD
                                  "19" = "Dendritic Cell", #TBD
                                  "20" = "Plasma Cell",
                                  "21" = "T Cell",
                                  "22" = "Acinar Cell"
)




DimPlot(T1D_Timepoints, reduction = "umap", 
        #split.by = c( "Tx"), 
        combine = F, label = T)

library(ggpubr)
DimPlot(T1D_Timepoints, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 7, repel = T) +
  theme_prism(base_size = 16, base_fontface = "bold") + 
  theme(legend.text = element_text(size = 18)) +
  labs( title = "UMAP of Cell Types") 
ggsave(file = "UMAP_CellTypeclusters.png",
    units = "in",width = 11, height = 11, dpi = 600)


saveRDS(T1D_Timepoints, file = "./annotated_T1D_Timepoints_v3.rds")


table(T1D_Timepoints$time,T1D_Timepoints$group)


# Subset the data for Week 6 and Week 12
week6_data <- subset(T1D_Timepoints,subset =time  %in% c("Week6"))
week12_data <- subset(T1D_Timepoints,subset =time  %in% c("Week12"))

# Create frequency tables for Week 6 and Week 12 using active.ident and group
freq_week6 <- prop.table(x = table(week6_data@active.ident,week6_data@meta.data[, "group"]), margin =2)

freq_week12 <- prop.table(x = table(week12_data@active.ident,week12_data@meta.data[, "group"]), margin =2)

         # Set the resolution to 300 DPI
# Plot for Week 6
par(mar = c(5, 6, 4, 6))  # Increasing the right margin to 6

# Plot the bar plot with modified x-axis label rotation
barplot(
  freq_week6, 
  xlim = c(0, ncol(freq_week6) + 3),
  col = bar_colors, 
  legend.text = TRUE, 
  args.legend = list(
    x = ncol(freq_week6) + 3, 
    y = max(colSums(freq_week6)), 
    bty = "n", 
    cex = 1.5
  ), 
  width = 1.25, 
  ylab = "Fraction of Cell Population",
  cex.lab = 2,   # Increase font size of axis titles
  cex.names = 2,  # Increase font size of x-axis labels
  cex.axis = 1.5     # Increase font size of axis tick labels
)


# Plot for Week 12
par(mar = c(5, 6, 4, 6))  # Increasing the right margin to 6

# Plot the bar plot with modified x-axis label rotation
barplot(
  freq_week12, 
  xlim = c(0, ncol(freq_week12) + 3),
  col = bar_colors, 
  legend.text = TRUE, 
  args.legend = list(
    x = ncol(freq_week12) + 3, 
    y = max(colSums(freq_week12)), 
    bty = "n", 
    cex = 1.5
  ), 
  width = 1.25, 
  ylab = "Fraction of Cell Population",
  cex.lab = 2,   # Increase font size of axis titles
  cex.names = 2,  # Increase font size of x-axis labels
  cex.axis = 1.5     # Increase font size of axis tick labels
)

ggsave(file = "CellProportions_9clusters.png",
    units = "in",width = 11, height = 11, dpi = 400)
dev.off()



## ----table of total cells per sample------------------------------------------------------------------------------------------------------------------
table(T1D_Timepoints@meta.data$sample,T1D_Timepoints@meta.data$sample)

#' 
## ----more cluster ids---------------------------------------------------------------------------------------------------------------------------------

immune.markers.paper = read.csv("mmc2.csv",
                                check.names=FALSE,header=T, sep=",")
unique(immune.markers.paper$Column2)

immune.markers.paper[which(immune.markers.paper$Column2 == "Macrophage" | immune.markers.paper$Column2 == "Macrophages"),3]

plots <- VlnPlot(int.T1D_Timepoints, features = c("Cd3e", "Cd4", "Cd8a","Cd8b1"), 
                 split.by = "Tx",
                 group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)


plots <- VlnPlot(int.T1D_Timepoints, features = c("Cd3e", "Cd4", "Cd8a","Cd8b1", "Cd19", "Ms4a1", "Sdc1", "Cd27", "Slamf7"), 
                 split.by = "Tx",
                 group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)


#macs
c("Marco", "Itgam", "Adgre1",
                               "Apoe","Cd14", "Lyz",
                               "Ms4a7", "Pf4", "Fapb4", "Gpnmb", "Gpr84"
)

#stromal
c("Ptprc","Col1a1", "Col6a1", "Dcn", "Lum", "Postn", "Eef1akmt3"
                               
)


#mast cell
c(
  "Mcpt1", "Fcer1a", "Ighe", "Kit", "Itgb7"
)

DotPlot(int.T1D_Timepoints, features = c("Cd3e","Marco", "Itgam", "Adgre1", "Adgre4",
                               "Apoe","Cd14", "Lyz",
                               "Ms4a7", "Pf4", "Fapb4", "Gpnmb", "Gpr84", "Cd163", "Cd68",
                               "Mrc1", "Msr1", "Fcgr3", "Ly6c2", "Arg1"
)
, group.by = "seurat_clusters") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")


DotPlot(int.T1D_Timepoints, features = c("Jchain","S100a9", "S100a8",
                               "Ccr3", "Epx", "Edn1", "Siglecf",
                               "Il1a", "Lin28a", "Fcgr3",
                               "Mki67", "Tuba1b", "Prg3", "Ear1", "Ear2", "Ear6", "Alox15", "S100a6", "S100a10", "Aldh2", "Il5",
                            "Itgam", "Il5ra","Siglec5"

                               ), group.by = "seurat_clusters") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

DotPlot(int.T1D_Timepoints, features = 
          unique(immune.markers.paper[which(immune.markers.paper$Column2 == "Macrophage" | immune.markers.paper$Column2 == "Macrophages"),3])
          , group.by = "seurat_clusters") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

DotPlot(int.T1D_Timepoints, features = 
          unique(immune.markers.paper[which(immune.markers.paper$Column2 == "monocyte" | immune.markers.paper$Column2 == "Monocyte" 
                                            | immune.markers.paper$Column2 == "Monocytes"),3])
          , group.by = "seurat_clusters") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

DotPlot(int.T1D_Timepoints, features = 
          unique(immune.markers.paper[which(immune.markers.paper$Column2 == "NK_cells" |immune.markers.paper$Column2 == "NK_cell"|
                                              immune.markers.paper$Column2 == "NK"),3])
          , group.by = "seurat_clusters") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

DotPlot(int.T1D_Timepoints, features = c("Itgam", "Itgax")
          , group.by = "seurat_clusters") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")


DotPlot(int.T1D_Timepoints, features =c("Tbx21", "Gzmb", "Gzmc","Il7r",
                              "Cd3e","Cd5","Cx3cl1","Rorc",
                              "Gata3","Ifng","Il4","Il22","Il17a",
                              "Foxs1","Gpx1","Tcf7", "Cxcl2","Podnl1","Socs3",
                              "Atf3", "Maff", "Ostf1","Rcor1","Smox",
                              "Ccl20", "Tnfsf11","Il23r","Lta","Tnf",
                              "Hk2", "Hk3", "Pfki", "Pfkp","Eno1","Ldhb",
                              "Bpgm","Tns4","Tns3","Arg1","Pdk1","Mgll",
                              "Mtmr12","Plcg2")
          , group.by = "seurat_clusters") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")


DotPlot(int.T1D_Timepoints, features = top20.int[which(top20.int$cluster_id == "11"),2]
          , group.by = "test.ids") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")


plots <- VlnPlot(int.T1D_Timepoints, features = c("Cd19", "Ms4a1", "Sdc1", "Tnfrsf13b", "Slamf7", "Ifngr1", "Cd27"), 
                 split.by = "sample",
                 group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

FeaturePlot(int.T1D_Timepoints, 
            reduction = "umap", 
            features = c("Cd4", "Cd8a", "Cd8b1"), 
            #sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

plots <- VlnPlot(int.T1D_Timepoints, features = c("Ptprc", "Flt1", "Ly6a",
                                    "Id3", "Fhl2", "Gng11", "Cd34", "Acvrl1",
                                    "Itgb3", "Vcam1", "Thbd", "Cd55", "Vegfa", "Vegfd",
                                    "Ece1", "Plec", "Ecscr","Tgfbr2", "Icam1"), split.by = "Tx",
                 group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)

plots <- VlnPlot(int.T1D_Timepoints, features = c("Cd14", "Lyz1"), 
                 split.by = "Tx",
                 group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)

#DCs
plots <- VlnPlot(int.T1D_Timepoints, features = c("Cd209a", "Btla", "Flt3", "Itgax", #cd11c
                                    "Clec9a", "Ly75", "Sirpa", "Bst2", "Itgae", #CD103
                                    "Itgam", #Cd11b
                                    "Cd80", "Cd86", "Il4i1"
), 
split.by = "Tx",
group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

plots <- VlnPlot(int.T1D_Timepoints, features = c("Flt3","H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa",
                                    "Cd80", "Cd86", "Cd83"
), 
split.by = "sample",
group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)
Idents(TNBC) = "seurat_clusters"
DCs_types_markers <- FindConservedMarkers(int.T1D_Timepoints, ident.1 = "14", ident.2 =  c("7", "12"),grouping.var = "Antigen",
                                          verbose = FALSE)

data[[]]
plots <- VlnPlot(int.T1D_Timepoints, features = c(
  "Mcpt1", "Fcer1a", "Ighe", "Kit", "Itgb7"
), 
split.by = "Tx",
group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

FeaturePlot(int.T1D_Timepoints, 
            reduction = "umap", 
            features = c("Mcpt1", "Fcer1a", "Ighe", "Kit", "Itgb7"), 
            #sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

plots <- VlnPlot(int.T1D_Timepoints, features = c(
  "S100a8", "S100a9"
), 
split.by = "Tx",
group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

plots <- VlnPlot(int.T1D_Timepoints, features = c(
  "Ncr1", "Klrb1a", "Klre1"
), 
split.by = "Tx",
group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)




DotPlot(int.T1D_Timepoints, features = c("Cd4","Marco", "Itgam","Itgax", "Adgre1", "Adgre4",
                               "Apoe","Cd14", "Lyz",
                               "Ms4a7", "Pf4", "Fapb4", "Gpnmb", "Gpr84", "Cd163", "Cd68",
                               "Mrc1", "Msr1", "Fcgr3", "Ly6c2", "Arg1", "Mx1", "Ccr3",
                               "Timd4", "Spp1", "Slamf7", "Cxcl2", "Ccr2", "Cx3cr1", "Il1a", 
                               "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa"
)
, group.by = "seurat_clusters") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

Idents(int.T1D_Timepoints) = int.T1D_Timepoints$SCT_snn_res.0.6
My11v24 = FindMarkers(object = int.T1D_Timepoints, ident.1 = "11", ident.2 = "24")

VlnPlot(int.T1D_Timepoints, features = c(
  "Tmsb4x"
), 
split.by = "Tx",
group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)


## ----plots--------------------------------------------------------------------------------------------------------------------------------------------


DimPlot(int.T1D_Timepoints, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T, 
        group.by = "test.ids",
        cols = getPalette(9)) +
  theme_prism(base_size = 16, base_fontface = "bold") + 
  theme(legend.text = element_text(size = 18)) +
  labs( title = "UMAP of Cell Types") 
ggsave(file = "UMAP_9clusters.png",
    units = "in",width = 11, height = 11, dpi = 400)


int.T1D_Timepoints$Cell.types = factor(int.T1D_Timepoints$Cell.types, levels =  c("CD4 T Cells", "CD8 T Cells","Th17 ILC","Th2 ILC",  "B Cells", "Plasma Cells", "NK", "Neutrophil", "Mast Cells", 
    "DCs", "CD103+ DCs", "Inflam. Mac", "MHCII+ Mac", "Stromal Cells"  ) )


DimPlot(int.T1D_Timepoints, reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 4, repel = T, 
        group.by = "Cell.types",
        cols = getPalette(14)) +
  theme_prism(base_size = 12, base_fontface = "bold") + 
  theme(legend.text = element_text(size = 14)) +
  labs( title = "UMAP of Cell Types") 
ggsave(file = "UMAP_14clusters.png",
    units = "in",width = 6, height = 4, dpi = 400)

DimPlot(int.T1D_Timepoints, reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 5, repel = T,
        group.by = "Cell.types",
        split.by = c( "Tx"),
        cols = getPalette(14)) +
  theme_prism(base_size = 16, base_fontface = "bold") + 
  theme(legend.text = element_text(size = 18)) +
  labs( title = "UMAP of Cell Types") 

ggsave(file = "UMAP_14clusters_SplitTx.png",
    units = "in",width = 10, height = 5, dpi = 400)

Idents(int.T1D_Timepoints) = "Tx"

DimPlot(subset(int.T1D_Timepoints, idents = "PBS"), reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 5, repel = T,
        group.by = "Cell.types",
        #split.by = c( "Tx"),
        cols = getPalette(14)) +
  theme_prism(base_size = 16, base_fontface = "bold") + 
  theme(legend.text = element_text(size = 18)
       ) +
  labs( title = "UMAP of Cell Types") 

ggsave(file = "UMAP_14clusters_PBS.png",
    units = "in",width = 7, height = 5, dpi = 400)



DimPlot(int.T1D_Timepoints, reduction = "umap", group.by = c("Cell.types"), split.by = c( "Tx"), combine = F, label = T)

table(int.T1D_Timepoints$test.ids,int.T1D_Timepoints$sample)

CellTypeTable = as.data.frame(proportions(table(int.T1D_Timepoints$test.ids,
                                                int.T1D_Timepoints$sample), 
                                          margin = 2))
CellTypeTable$Freq = CellTypeTable$Freq*100
levels(CellTypeTable$Var1)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#display.brewer.all()
CellTypeTable[which(CellTypeTable$Var1 == "DCs"),]


CellTypeTable$Var2 = factor(CellTypeTable$Var2, levels = c("PBS PreCh1", "PBS Ch4", "PBS Ch7", "NP Ch4", "NP Ch7"))
#png(file = "CellProportions_24clusters.png",units = "in",width = 14, height = 13, res = 400)
ggplot(CellTypeTable,aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_col(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism( base_size = 16)+
  labs(x = "Condition", y= "Percent of Cells") + 
  scale_fill_manual(values = getPalette(length(unique(int.T1D_Timepoints$test.ids))), name = "Cell Type") + 
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  )
ggsave(file = "CellProportions_9clusters.png",
    units = "in",width = 11, height = 11, dpi = 400)
dev.off()


CellTypeTable = as.data.frame(proportions(table(int.T1D_Timepoints$Cell.types,
                                                int.T1D_Timepoints$sample), 
                                          margin = 2))
CellTypeTable$Freq = CellTypeTable$Freq*100
levels(CellTypeTable$Var1)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#display.brewer.all()


CellTypeTable$Var2 = factor(CellTypeTable$Var2, levels = c("PBS PreCh1", "PBS Ch4", "PBS Ch7", "NP Ch4", "NP Ch7"))
#png(file = "CellProportions_24clusters.png",units = "in",width = 14, height = 13, res = 400)
ggplot(CellTypeTable, aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_col(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism( base_size = 10)+
  labs(x = "Condition", y= "Percent of Cells") + 
  scale_fill_manual(values = getPalette(length(unique(int.T1D_Timepoints$Cell.types))), name = "Cell Type") + 
  theme(legend.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  )
ggsave(file = "CellProportions_14clusters.png",
    units = "in",width = 5, height = 4.5, dpi = 400)


CellTypeTable$TP = gsub("PBS ", "", CellTypeTable$Var2)
ggplot(CellTypeTable[
  which(CellTypeTable$Var2 == "PBS PreCh1" |
        CellTypeTable$Var2 == "PBS Ch4" |
          CellTypeTable$Var2 == "PBS Ch7"
          ),
],aes(x=TP,y=Freq,fill=Var1)) + 
  geom_col(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism( base_size = 14)+
  labs( y= "Percent of Cells") + 
  scale_fill_manual(values = getPalette(length(unique(int.T1D_Timepoints$Cell.types))), name = "Cell Type") + 
  theme(legend.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        #plot.title = element_text(hjust = 0.5)
        axis.title.x = element_blank()
  )
ggsave(file = "CellProportions_14clusters_PBS.png",
    units = "in",width = 5, height = 4.5, dpi = 400)
dev.off()

#Cluster IDs
DotPlot(int.T1D_Timepoints, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
  "Col4a1","Flt1","Nkx2-3", "Tie1", #Stromal
  "Itgam", "Itgax", "Adgre4","Adgre1","Cd80", "Clec9a",
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", "H2-K1", #Myeloid
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe"  #Plasma 
  )
          , group.by = "test.ids") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism(base_family = "Arial", base_size = 20) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1))

ggsave(file = "DotplotMarkergenes_9clusters.png",
    units = "in",width = 17, height = 8, dpi = 400)


#Cluster IDs
DotPlot(int.T1D_Timepoints, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
  "Col4a1","Flt1","Nkx2-3", "Tie1", #Stromal
  "S100a8", "S100a9", #Neut
  "Clec9a",  "Itgax", "Cd209a", # DCs /CD11c+ SiglecH DC-sign(Cd209a/d)
  "Itgae",#CD103
  "Itgam", "Adgre1", "Adgre4","Apoe","Cd14","Ms4a7", "Cd68", #Mac
   "Mrc1", "Msr1", "Fcgr3",  "Ccr3", "Cxcl2",  "Cx3cr1", "Il1a", "Il1b","Ccr2", #Mac
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe"  #Plasma 
  )
          ) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism(base_family = "Arial", base_size = 12) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1, size = 8, face = "plain"))

ggsave(file = "DotplotMarkergenes_14clusters.png",
    units = "in",width = 10, height = 5, dpi = 400)

CellTypeData = CellTypeTable
# load stringr library
library(stringr)
library(ggpubr)

# Split name column into firstname and last name
CellTypeData[c('Tx', 'TP')] <- str_split_fixed(CellTypeData$Var2, ' ', 2)

plotCellType = function(data, celltype){
  ggplot(data[which(data$Var1 == celltype),], aes(x = Tx, y= Freq, fill = TP)) + 
    geom_bar(position = position_dodge(),stat = 'identity') +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    theme_prism(base_family = "Arial", base_size = 12)+
    labs( y= "Percent of Cells", title = paste(celltype)) + 
    theme(legend.text = element_text(size = 12) 
    ) + xlab(NULL)
}

CellTypeData$Tx = factor(CellTypeData$Tx, levels = c("PBS", "NP"))
myplots = lapply(unique(CellTypeData$Var1), FUN = plotCellType, data = CellTypeData)

#png(file = "CellProportions_24clusters_individual_grid.png", units = "in",width = 20, height = 16, res = 400)
ggarrange(plotlist = myplots, ncol=6, nrow=4, 
          common.legend = T, legend = "top")
dev.off()

plotCellTypePresentation = function(data,celltype){
  ggplot(data[which(data$Var1 == celltype),], aes(x = Tx, y= Freq, fill = TP)) + 
    geom_bar(position = position_dodge(),stat = 'identity', width = 0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    theme_prism(base_family = "Arial", base_size = 16)+
    labs(x = "Condition", y= "Percent of Cells", title = paste(celltype)) + 
    theme(legend.text = element_text(size = 16)
    )
}

CellTypeData$TP = factor(CellTypeData$TP, levels = c("PreCh1", "Ch4", "Ch7"))


plotCellTypePresentation(CellTypeData, "ILC")


myplots = lapply(unique(CellTypeData$Var1), FUN = plotCellTypePresentation, data = CellTypeData)

png(file = "ImmuneCell_prop.png",
    units = "in",width = 17, height = 8, res = 400)
ggarrange(plotlist = myplots, ncol=5, nrow=3, 
          common.legend = T, legend = "right")

ggsave(file = "ImmuneCell_prop_14cluster.png",
    units = "in",width = 17, height = 15
    , dpi = 400)
dev.off()





#' 
## ----Cluster ID dotplot-------------------------------------------------------------------------------------------------------------------------------
#Cluster IDs
Idents(int.T1D_Timepoints) ="Cell.types"
DotPlot(int.T1D_Timepoints, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
  "Col4a1","Flt1","Nkx2-3", "Tie1", #Stromal
  "S100a8", "S100a9", #Neut
  "Clec9a",  "Itgax", "Cd209a", # DCs /CD11c+ SiglecH DC-sign(Cd209a/d)
  "Itgae",#CD103
  "Itgam", "Adgre1", "Adgre4","Apoe","Cd14","Ms4a7", "Cd68", #Mac
   "Mrc1", "Msr1", "Fcgr3",  "Ccr3", "Cxcl2",  "Cx3cr1", "Il1a", "Il1b","Ccr2", #Mac
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe"  #Plasma 
  ), group.by = "Cell.types"
          ) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism( base_size = 12) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1, size = 8, face = "plain"))

#' 
## ----Heatmap of cluster ids---------------------------------------------------------------------------------------------------------------------------

DoHeatmap(int.T1D_Timepoints, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
  "Col4a1","Flt1","Nkx2-3", "Tie1", #Stromal
  "S100a8", "S100a9", #Neut
  "Clec9a",  "Itgax", "Cd209a", # DCs /CD11c+ SiglecH DC-sign(Cd209a/d)
  "Itgae",#CD103
  "Itgam", "Adgre1", "Adgre4","Apoe","Cd14","Ms4a7", "Cd68", #Mac
   "Mrc1", "Msr1", "Fcgr3",  "Ccr3", "Cxcl2",  "Cx3cr1", "Il1a", "Il1b","Ccr2", #Mac
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe"  #Plasma 
  ), group.by = "Cell.types", assay = "SCT", slot = "scale.data"
          ) 

#' 
#' #Subcluster B cells
#' 
## ----subcluster B cells-------------------------------------------------------------------------------------------------------------------------------
Idents(int.T1D_Timepoints) = int.T1D_Timepoints$test.ids

names(int.T1D_Timepoints@graphs)

int.T1D_Timepoints = FindSubCluster(
  object = int.T1D_Timepoints,
  cluster = "B Cells",
  graph.name = "SCT_nn",
  resolution = 0.5,
  subcluster.name = "B Cells Subtypes",
  algorithm = 2
)

unique(int.T1D_Timepoints$`B Cells Subtypes`)

Idents(int.T1D_Timepoints) = int.T1D_Timepoints$`B Cells Subtypes`

DimPlot(int.T1D_Timepoints, reduction = "umap", split.by = c( "Tx"), combine = F, label = T)

table(Idents(Bcells),Bcells@meta.data$sample)



Idents(int.T1D_Timepoints) = int.T1D_Timepoints$test.ids
Bcells = subset(x = int.T1D_Timepoints, idents = "B Cells")
Bcells <- FindNeighbors(Bcells, reduction = "integrated.dr", dims = 1:30)
Bcells <- FindClusters(Bcells, res = 0.3, graph.name = "SCT_nn")
Bcells <- RunUMAP(Bcells, dims = 1:20, reduction = "integrated.dr")

DimPlot(Bcells, reduction = "umap", split.by = c( "Tx"), combine = F, label = T)


get_conserved <- function(cluster){
  Seurat::FindConservedMarkers(Bcells,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene")  %>%
    #left_join(y = unique(annotations[, c("gene_name", "description")]),
    #          by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
unique(Idents(Bcells))

conserved_markers.Bcells <- map_dfr(0:5, get_conserved)
conserved_markers.Bcells[is.na(conserved_markers.Bcells)] <- 0

DefaultAssay(Bcells) = "SCT"
head(conserved_markers.Bcells)
conserved_markers$`PBS PreCh1_avg_log2FC`
# Extract top 10 markers per cluster
conserved_markers.Bcells <- conserved_markers.Bcells %>% 
  mutate(avg_fc = (`PBS PreCh1_avg_log2FC` + 
                     `PBS Ch4_avg_log2FC` + 
                     `NP Ch4_avg_log2FC` +
                     `PBS Ch7_avg_log2FC` +
                     `NP Ch7_avg_log2FC`) /5)



top20.Bcells = top_n(
    conserved_markers.Bcells[which(conserved_markers.Bcells$cluster_id == 0),],
        n = 20, 
        wt = avg_fc)

for(i in 1:5){
  
 
   subset = top_n(
  conserved_markers.Bcells[which(conserved_markers.Bcells$cluster_id == 5),],
  n = 20, 
        wt = avg_fc)
  
  top20.Bcells = rbind(top20.Bcells, subset)
}

DotPlot(Bcells, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #Myeloid
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3"  #Plasma 
  ), scale = F
          ) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") 
  
colnames(Bcells[[]])

Idents(Bcells)

B0 = subset(x = Bcells, idents = "0")
B1 = subset(x = Bcells, idents = "1")
B2 = subset(x = Bcells, idents = "2")
Idents(B0) = B0$sample
Idents(B1) = B1$sample
Idents(B2) = B2$sample

B0 = PrepSCTFindMarkers(B0)
B1 = PrepSCTFindMarkers(B1)
B2 = PrepSCTFindMarkers(B2)
B0_ch7_PBSvNP = FindMarkers(B0, ident.1 = "PBS Ch7", ident.2 = "NP Ch7", assay= "SCT",
                            recorrect_umi = FALSE)
B1_ch7_PBSvNP = FindMarkers(B1, ident.1 = "PBS Ch7", ident.2 = "NP Ch7", assay= "SCT",
                            recorrect_umi = FALSE)
B2_ch7_PBSvNP = FindMarkers(B2, ident.1 = "PBS Ch7", ident.2 = "NP Ch7", assay= "SCT",
                            recorrect_umi = FALSE)

B0_ch7_PBSvNP$symbol = rownames(B0_ch7_PBSvNP)
B1_ch7_PBSvNP$symbol = rownames(B1_ch7_PBSvNP)
B2_ch7_PBSvNP$symbol = rownames(B2_ch7_PBSvNP)

library(msigdbr)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(forcats)
library(EnhancedVolcano)

all_gene_sets = msigdbr(species = "Mus musculus")
unique(all_gene_sets$gs_subcat)
unique(all_gene_sets$gs_cat)

kegg_gene_sets_msig = msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG")


biocarta_gene_sets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:BIOCARTA")

Reactome_gene_sets = msigdbr(species = "mouse", subcategory = "CP:REACTOME")
Wiki_gene_sets = msigdbr(species = "mouse", subcategory = "CP:WIKIPATHWAYS")

HM_gene_sets = msigdbr(species = "mouse", category = "H")

GO_gene_sets = msigdbr(species = "mouse", subcategory = "GO:BP")


B0_ch7_PBSvNP_geneList =  B0_ch7_PBSvNP %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(symbol, avg_log2FC)

B0_ch7_PBSvNP_ranks<- deframe(B0_ch7_PBSvNP_geneList)

B1_ch7_PBSvNP_geneList =  B1_ch7_PBSvNP %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(symbol, avg_log2FC)

B1_ch7_PBSvNP_ranks<- deframe(B1_ch7_PBSvNP_geneList)


B2_ch7_PBSvNP_geneList =  B2_ch7_PBSvNP %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(symbol, avg_log2FC)

B2_ch7_PBSvNP_ranks<- deframe(B2_ch7_PBSvNP_geneList)


combined_data_set = rbind(kegg_gene_sets_msig, 
                          biocarta_gene_sets,
                          Reactome_gene_sets,
                          Wiki_gene_sets,
                          HM_gene_sets,
                          GO_gene_sets)

pathway_gene_set = combined_data_set[,c("gs_name", "gene_symbol")]




y = GSEA(B0_ch7_PBSvNP_ranks,
         pvalueCutoff = 0.5,
         pAdjustMethod = "fdr",
         #OrgDb = org.Mm.eg.db, 
         #TERM2GENE = CTD_gene_sets,
         #TERM2GENE = kegg_gene_sets_msig[,c("gs_description", "entrez_gene")],
         #TERM2GENE = biocarta_gene_sets[,c("gs_description", "entrez_gene")],
         #TERM2GENE = Reactome_gene_sets[,c("gs_description", "entrez_gene")],
         TERM2GENE = pathway_gene_set,
         minGSSize = 1
         
)

y1 = GSEA(B1_ch7_PBSvNP_ranks,
         pvalueCutoff = 0.5,
         pAdjustMethod = "fdr",
         #OrgDb = org.Mm.eg.db, 
         #TERM2GENE = CTD_gene_sets,
         #TERM2GENE = kegg_gene_sets_msig[,c("gs_description", "entrez_gene")],
         #TERM2GENE = biocarta_gene_sets[,c("gs_description", "entrez_gene")],
         #TERM2GENE = Reactome_gene_sets[,c("gs_description", "entrez_gene")],
         TERM2GENE = pathway_gene_set,
         minGSSize = 1
         
)

y2 = GSEA(B2_ch7_PBSvNP_ranks,
         pvalueCutoff = 0.5,
         pAdjustMethod = "fdr",
         #OrgDb = org.Mm.eg.db, 
         #TERM2GENE = CTD_gene_sets,
         #TERM2GENE = kegg_gene_sets_msig[,c("gs_description", "entrez_gene")],
         #TERM2GENE = biocarta_gene_sets[,c("gs_description", "entrez_gene")],
         #TERM2GENE = Reactome_gene_sets[,c("gs_description", "entrez_gene")],
         TERM2GENE = pathway_gene_set,
         minGSSize = 1
         
)


head(y)

y_filt = y[which(abs(y$qvalue)<0.25),]
#y_filt = as.data.frame(y)

y_filt$Condition = "PBS"
y_filt[which((y_filt$NES)<0),]$Condition = "NP" 

y_filt$Day = "Challenge 7"
y_filt$Subset = "Cluster 0 Bcells"

y_filt1 = y1[which(abs(y1$qvalue)<0.25),]
#y_filt = as.data.frame(y)

y_filt1$Condition = "PBS"
y_filt1[which((y_filt1$NES)<0),]$Condition = "NP" 

y_filt1$Day = "Challenge 7"
y_filt1$Subset = "Cluster 1 Bcells"

y_filt2 = y2[which(abs(y2$qvalue)<0.25),]
#y_filt = as.data.frame(y)

y_filt2$Condition = "PBS"
y_filt2[which((y_filt2$NES)<0),]$Condition = "NP" 

y_filt2$Day = "Challenge 7"
y_filt2$Subset = "Cluster 2 Bcells"

y_filt3 = rbind(y_filt, y_filt1, y_filt2)

y_filt3 = y_filt3[which(abs(y_filt3$NES)>1.9),]

ggplot(y_filt3, aes(x = Subset, y = fct_reorder(Description, abs(NES)))) + 
  geom_point(aes(color = abs(NES), size = qvalue) ) +
  theme_bw(base_size = 14) +
  scale_colour_gradient2( low = "blue",high="red", mid = "white",na.value = "black") +
  ylab(NULL) +
  ggtitle("Biocarta Pathways GSEA at Challenge 7 in T1D_Timepoints") +
  facet_grid(~Condition) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 16),
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x =element_text(size = 16, face = "bold"),
        axis.title.x = element_blank(),
        strip.text.x  = element_text(size = 16, face = "bold")) 
ggsave(file = "BcellsubsetGSEA.png",
    units = "in",width = 25, height = 20, dpi = 400)

data = B0_ch7_PBSvNP
keyvals <- ifelse(
  (data$avg_log2FC < -1.5 & data$p_val_adj <0.05), 'royalblue',
  ifelse((data$avg_log2FC > 1.5 & data$p_val_adj < 0.05), 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'black'] <- 'Not DEG'
names(keyvals)[keyvals == 'red'] <- 'PBS'
names(keyvals)[keyvals == 'royalblue'] <- 'NP'

#png(file = "AllergyScaf_NPPBSCh4_VolcanoPlot.png",units = "in",width = 10, height = 10, res = 400)
EnhancedVolcano(data,
                lab = rownames(data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colCustom = keyvals,
                selectLab = rownames(data)[which(names(keyvals) %in% c('NP', 'PBS'))],
                title = "DEGs b/w NP and PBS at Challenge 7",
                drawConnectors = F,
                widthConnectors = 0.75,
                boxedLabels = F)
dev.off()



# pseudobulk ----
Idents(int.T1D_Timepoints) = "sample"


# pseudobulk cells by stimulation condition AND cell type AND donor
bulk <- AggregateExpression(int.T1D_Timepoints, group.by = c("Tx", "Cell.types", "TP"),
                            assays = "RNA",
                            features = NULL,
                            return.seurat = FALSE,
                            #group.by = "Cell.types",
                            add.ident = NULL,
                            slot = "counts",
                            verbose = TRUE)

pseudobulk_PBS1 = AggregateExpression(
  subset(int.T1D_Timepoints, idents = c("PBS PreCh1")),
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "Cell.types",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)

Idents(TNBC) = "Day"
pseudobulkD9 = AggregateExpression(
  subset(TNBC,idents = c("D9")),
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "new.cluster.ids",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)
Idents(TNBC) = "new.cluster.ids"
head(pseudobulkD9)

Idents(TNBC) = "Day"
pseudobulkD7 = AggregateExpression(
  subset(TNBC,idents = c("D7")),
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "new.cluster.ids",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)
Idents(TNBC) = "new.cluster.ids"
head(pseudobulkD7)


signature_genes1 = read.csv("/Users/lailarad/Documents/BI_EAE/aaron_pEAE/Signature_genes1.csv",row.names = NULL)
signature_genesTx = read.csv("/Users/lailarad/Documents/BI_EAE/Cohort1_RNA_seq/signature_genes_VLA4tx.csv",row.names = NULL)
signature_genes_day9_pEAE_bulk = read.csv("/Users/lailarad/Documents/BI_EAE/aaron_pEAE/aaronpEAE_signature_genes_day9.csv",row.names = NULL)




Features(int.T1D_Timepoints)
Features(int.T1D_Timepoints)[grepl("^H2-", Features(int.T1D_Timepoints))]
int.T1D_Timepoints[["RNA"]]$counts[grepl("^H2-", Features(int.T1D_Timepoints)),]
mat[grepl("^H2-", Features(int.T1D_Timepoints)),]

pseudobulkset = pseudobulk_PBS1
mat = pseudobulkset$RNA
mat=log(mat + 1, base = 2)

head(mat)

pseudobulkset = bulk
mat = pseudobulkset$RNA
mat=log(mat + 1, base = 2)
mat = mat - rowMeans(mat)



max(mat)



sigGenes = (intersect(rownames(mat), signature_genes1$x))
sigGenes = (intersect(rownames(mat), signature_genes_day9_pEAE_bulk$x))
sigGenes = (intersect(rownames(mat), signature_genesTx$genes))

pheatmap::pheatmap(t(mat),
                   #annotation_col = anno,
                   show_colnames = T,
                   fontsize_row = 18,
                   fontsize_col = 12,
                   cluster_cols = T,
                   cluster_rows = T,
                   #annotation_colors = ann_colors,
                   #filename = paste(path, "signatureCelltypesD79_vla4tx_log2.png",sep=""),
                   width = 8,
                   height = 5
                   
                   
)
dev.off()
getwd()
plots <- VlnPlot(TNBC, features = sigGenes[(25:31)], 
                 split.by = "sample",
                 group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)

DotPlot(TNBC, features = sigGenes) + 
  RotatedAxis()


## ----Top Gene changing over time per Cell Type in PBS mice--------------------------------------------------------------------------------------------
unique(int.T1D_Timepoints$sample)
Idents(int.T1D_Timepoints) = int.T1D_Timepoints$Cell.Type.TP.Tx

Idents(int.T1D_Timepoints) = "sample"

int.T1D_Timepoints$Cell.Type.TP.Tx = paste0(int.T1D_Timepoints$Cell.types, "_", int.T1D_Timepoints$TP, "_", int.T1D_Timepoints$Tx)
unique(int.T1D_Timepoints$Cell.Type.TP.Tx)

PBS1_7 = subset(int.T1D_Timepoints, idents = c("PBS PreCh1", "PBS Ch7") )

TH17ILC_1_4 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 = "Th17 ILC_Ch4_PBS",
  ident.2 = "Th17 ILC_PreCh1_PBS", 
  only.pos = T
)
TH17ILC_4_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "Th17 ILC_Ch7_PBS",
  ident.2 =  "Th17 ILC_Ch4_PBS",
  only.pos = T
)

TH17ILC_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "Th17 ILC_Ch7_PBS",
  ident.2 =  "Th17 ILC_PreCh1_PBS",
  only.pos = T
)

TH17ILC_1_4 = TH17ILC_1_4[which(TH17ILC_1_4$p_val_adj <0.05),]
TH17ILC_4_7 = TH17ILC_4_7[which(TH17ILC_4_7$p_val_adj <0.05),]
TH17ILC_1_7 = TH17ILC_1_7[which(TH17ILC_1_7$p_val_adj <0.05& TH17ILC_1_7$pct.1 > 0.5),]
TH17ILC_4_7[intersect(rownames(TH17ILC_1_4), rownames(TH17ILC_4_7)),]
TH17ILC_1_4[intersect(rownames(TH17ILC_1_4), rownames(TH17ILC_4_7)),]

int.T1D_Timepoints$TP = factor(int.T1D_Timepoints$TP,levels = c("PreCh1", "Ch4", "Ch7") )
Idents(PBS1_7) = "Cell.types"
unique(Idents(PBS1_7))
TH17 = subset(PBS1_7, idents = "Th17 ILC")
VlnPlot(TH17, features = "Bcl2",
        group.by = "TP",assay = "RNA")


Idents(int.T1D_Timepoints) = int.T1D_Timepoints$Cell.Type.TP.Tx
unique(Idents(int.T1D_Timepoints))

CD8T_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "CD8 T Cells_Ch7_PBS",
  ident.2 =  "CD8 T Cells_PreCh1_PBS",
  only.pos = F
)

CD8T_1_7 = CD8T_1_7[which(CD8T_1_7$p_val_adj <0.05),]

unique(Idents(PBS1_7))
CD8Tsub = subset(PBS1_7, idents = "CD8 T Cells")

CD8_vp = VlnPlot(CD8Tsub, features = "Ccl5",
        group.by = "TP",assay = "RNA", pt.size = 1) + 
  labs(x = "CD8 T Cells") +ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

CD8_vp$layers[[1]]$aes_params$size = 1.2
CD8_vp + NoLegend()

CD4T_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "CD4 T Cells_Ch7_PBS",
  ident.2 =  "CD4 T Cells_PreCh1_PBS",
  only.pos = T
)

CD4T_1_7 = CD4T_1_7[which(CD4T_1_7$p_val_adj <0.05 & CD4T_1_7$pct.1 > 0.5),]

unique(Idents(PBS1_7))
CD4Tsub = subset(PBS1_7, idents = "CD4 T Cells")
CD4_vp = VlnPlot(CD4Tsub, features = "Ctla4",
        group.by = "TP",assay = "RNA", pt.size = 1) + labs(x = "CD4 T Cells") +
   ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

CD4_vp$layers[[1]]$aes_params$size = 1.2
CD4_vp + NoLegend()

Bcell_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "B Cells_Ch7_PBS",
  ident.2 =  "B Cells_PreCh1_PBS",
  only.pos = T
)

Idents(int.T1D_Timepoints) = int.T1D_Timepoints$Cell.Type.TP.Tx
Bcell_1_7df = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "B Cells_Ch7_PBS",
  ident.2 =  "B Cells_PreCh1_PBS",
  only.pos = F
)


write.csv(Bcell_1_7df, file = "BcellT1D_Timepoints_1_7.csv")

Bcell_1_7 = Bcell_1_7[which(Bcell_1_7$p_val_adj <0.05 & Bcell_1_7$pct.1 > 0.5),]

unique(Idents(PBS1_7))
Bcellsub = subset(PBS1_7, idents = "B Cells")
gene = "Egr1"
gene = "Apoe"
gene = "Scd1"
B_vp = VlnPlot(Bcellsub, features = gene,
        group.by = "TP",assay = "RNA", pt.size = 1) + labs(x = "B Cells") +
 # ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

B_vp$layers[[1]]$aes_params$size = 1.2
B_vp + NoLegend()

#Early Growth Response-1 Plays a Non-redundant Role in the Differentiation of B Cells into Plasma Cells
#https://doi.org/10.4110%2Fin.2015.15.3.161 


MHCMac_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "MHCII+ Mac_Ch7_PBS",
  ident.2 =  "MHCII+ Mac_PreCh1_PBS",
  only.pos = F
)


MHCMac_1_7df = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "MHCII+ Mac_Ch7_PBS",
  ident.2 =  "MHCII+ Mac_PreCh1_PBS",
  only.pos = F
)

write.csv(MHCMac_1_7df, file = "MHCMac_1_7df.csv")

MHCMac_1_7 = MHCMac_1_7[which(MHCMac_1_7$p_val_adj <0.05),]
MHCMacsub = subset(PBS1_7, idents = "MHCII+ Mac")
gene =  "H2-Aa"
gene = "Il1b"
gene = "Hcar2"
gene = "Selenop"
MHCMac_vp = VlnPlot(MHCMacsub, features = gene,
        group.by = "TP",assay = "RNA", pt.size = 1) + 
  labs(x = "MHC-II Mac") +ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

MHCMac_vp$layers[[1]]$aes_params$size = 1.2
MHCMac_vp + NoLegend()

InfMac_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "Inflam. Mac_Ch7_PBS",
  ident.2 =  "Inflam. Mac_PreCh1_PBS",
  only.pos = F
)


InfMac_1_7df = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "Inflam. Mac_Ch7_PBS",
  ident.2 =  "Inflam. Mac_PreCh1_PBS",
  only.pos = F
)


write.csv(InfMac_1_7df, file = "InfMac_1_7df.csv")

InfMac_1_7 = InfMac_1_7[which(InfMac_1_7$p_val_adj <0.05 & InfMac_1_7$pct.2 > 0.3),]

unique(Idents(PBS1_7))
InfMacsub = subset(PBS1_7, idents = "Inflam. Mac")
gene = "Tff2"
gene= "Il1b"
gene = "Cxcl2"
InfMac_vp = VlnPlot(InfMacsub, features = gene,
        group.by = "TP",assay = "RNA", pt.size = 1) + 
  labs(x = "Inflam. Mac") +ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

InfMac_vp$layers[[1]]$aes_params$size = 1.2
InfMac_vp + NoLegend()

DC_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "DCs_Ch7_PBS",
  ident.2 =  "DCs_PreCh1_PBS",
  only.pos = F
)

DC_1_7df = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "DCs_Ch7_PBS",
  ident.2 =  "DCs_PreCh1_PBS",
  only.pos = F
)


write.csv(DC_1_7df, file = "DC_1_7df.csv")

DC_1_7 = DC_1_7[which(DC_1_7$p_val_adj <0.05 & DC_1_7$pct.2 > 0.3),]

unique(Idents(PBS1_7))
DCsub = subset(PBS1_7, idents = "DCs")
gene = "Irf8"
gene = "Klra17"
gene = "Zdhhc14"
DC_vp = VlnPlot(DCsub, features = gene,
        group.by = "TP",assay = "RNA", pt.size = 1) + 
  labs(x = "DCs") +ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

DC_vp$layers[[1]]$aes_params$size = 1.2
DC_vp + NoLegend() 

CD103DC_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "CD103+ DCs_Ch7_PBS",
  ident.2 =  "CD103+ DCs_PreCh1_PBS",
  only.pos = F
)

CD103DC_1_7df = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "CD103+ DCs_Ch7_PBS",
  ident.2 =  "CD103+ DCs_PreCh1_PBS",
  only.pos = F
)

write.csv(CD103DC_1_7df, file = "CD103DC_1_7df.csv")

CD103DC_1_7 = CD103DC_1_7[which(CD103DC_1_7$p_val_adj <0.05),]

unique(Idents(PBS1_7))
CD103DCsub = subset(PBS1_7, idents = "CD103+ DCs")
gene = "Slpi"
gene = "Zdhhc14"
CD103DC_vp = VlnPlot(CD103DCsub, features = gene,
        group.by = "TP",assay = "RNA", pt.size = 1) +
  labs(x = "CD103 DCs") + ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

CD103DC_vp$layers[[1]]$aes_params$size = 1.2
CD103DC_vp + NoLegend() 
#Secretory leukoprotease inhibitor in mucosal lymph node dendritic cells regulates the threshold for mucosal tolerance

Mast_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "Mast Cells_Ch7_PBS",
  ident.2 =  "Mast Cells_PreCh1_PBS",
  only.pos = F
)


Mast_1_7df = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "Mast Cells_Ch7_PBS",
  ident.2 =  "Mast Cells_PreCh1_PBS",
  only.pos = F
)

write.csv(Mast_1_7df, file = "Mast_1_7df.csv")

Mast_1_7 = Mast_1_7[which(Mast_1_7$p_val_adj <0.05 & Mast_1_7$pct.1>0.5),]

unique(Idents(PBS1_7))
Mastsub = subset(PBS1_7, idents = "Mast Cells")
gene = "Serpina3g"
gene = "Muc13"
gene = "Hdc"
mast_vp = VlnPlot(Mastsub, features = gene,
        group.by = "TP",assay = "RNA", , pt.size = 1) + 
  labs(x = "Mast Cells") + ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

mast_vp$layers[[1]]$aes_params$size = 1.2
mast_vp + NoLegend() 

#high in BMMC and low perotineal 
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7564378/


Neut_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "Neutrophil_Ch7_PBS",
  ident.2 =  "Neutrophil_PreCh1_PBS",
  only.pos = F
)

Neut_1_7df = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "Neutrophil_Ch7_PBS",
  ident.2 =  "Neutrophil_PreCh1_PBS",
  only.pos = F
)

write.csv(Neut_1_7df, file = "Neut_1_7df.csv")

Neut_1_7 = Neut_1_7[which(Neut_1_7$p_val_adj <0.05),]

unique(Idents(PBS1_7))
Neutsub = subset(PBS1_7, idents = "Neutrophil")
gene = "Ppp1r16b"
gene = "Cxcr2"
Neut_vp = VlnPlot(Neutsub, features = gene ,
        group.by = "TP",assay = "RNA", pt.size = 1, layer = "data") +
  labs(x = "Neut") + ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

Neut_vp$layers[[1]]$aes_params$size = 1.2
Neut_vp + NoLegend() 

NK_1_7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "NK_Ch7_PBS",
  ident.2 =  "NK_PreCh1_PBS",
  only.pos = F
)

NK_1_7df = FindMarkers(
  int.T1D_Timepoints,
  ident.1 =  "NK_Ch7_PBS",
  ident.2 =  "NK_PreCh1_PBS",
  only.pos = F
)

write.csv(NK_1_7df, file = "NK_1_7df.csv")

NK_1_7 = NK_1_7[which(NK_1_7$p_val_adj <0.05 & NK_1_7$pct.1 > 0.5),]

unique(Idents(PBS1_7))
NKsub = subset(PBS1_7, idents = "NK")
gene = "Gzma"
gene = "Acaca"
gene = "Nr4a2"
NK_vp = VlnPlot(NKsub, features = gene,
        group.by = "TP",assay = "RNA", pt.size = 1) +
   labs(x = "NKs") + ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

NK_vp$layers[[1]]$aes_params$size = 1.2
NK_vp + NoLegend() 
library(ggpubr)

png(file = "CellType_VlnPlts_Ch17.png", units = "in",width = 5, height = 20, res = 400)
ggarrange(Neut_vp,NK_vp,mast_vp,CD103DC_vp , DC_vp , InfMac_vp , MHCMac_vp , B_vp , CD8_vp , CD4_vp, ncol = 10, nrow = 1, legend = F) 
ggsave(file = "CellType_VlnPlts_Ch17.png",
    units = "px",width = 10000, height = 1200, dpi = 400)



patchwork = Neut_vp + NK_vp + mast_vp + CD103DC_vp + DC_vp + InfMac_vp + MHCMac_vp + B_vp + CD8_vp + CD4_vp


int.T1D_Timepoints$cell.types.combT = int.T1D_Timepoints$Cell.types

Idents(int.T1D_Timepoints) = int.T1D_Timepoints$cell.types.combT
unique(Idents(int.T1D_Timepoints))

Idents(int.T1D_Timepoints) = "sample"
PBS1_7 = subset(int.T1D_Timepoints, idents = c("PBS PreCh1", "PBS Ch7") )
Idents(PBS1_7) = PBS1_7$Cell.types
PBS1_7 = RenameIdents(PBS1_7, 

             "CD4 T Cells" = "T Cells",
             "CD8 T Cells" = "T Cells"
             )

PBS1_7$cell.types.combT = Idents(PBS1_7)

PBS1_7$Cell.Type.TP.Tx_combT = paste0(PBS1_7$cell.types.combT, "_", PBS1_7$TP, "_", PBS1_7$Tx)
unique(PBS1_7$Cell.Type.TP.Tx_combT)
Idents(PBS1_7) = "Cell.Type.TP.Tx_combT"

T_C_1_7 = FindMarkers(
  PBS1_7,
  ident.1 = "T Cells_Ch7_PBS",
  ident.2 = "T Cells_PreCh1_PBS", 
  only.pos = F
)


T_C_1_7df = FindMarkers(
  PBS1_7,
  ident.1 = "T Cells_Ch7_PBS",
  ident.2 = "T Cells_PreCh1_PBS", 
  only.pos = F
)
write.csv(T_C_1_7df, file = "T_C_1_7df.csv")

T_C_1_7 = T_C_1_7[which(T_C_1_7$p_val_adj <0.05 & T_C_1_7$pct.1 > 0.4),]
Idents(PBS1_7) = "cell.types.combT"
unique(Idents(PBS1_7))
Tsub = subset(PBS1_7, idents = "T Cells")
gene = "Il17rb"
gene = "Dtx3"
gene = "Sh2d1a"

T_vp = VlnPlot(Tsub, features = gene,
        group.by = "TP",assay = "RNA", pt.size = 1) + labs(x = "T Cells") +
   ylim(0,8) +
  theme_prism(base_size = 14,base_line_size = 1.5) +
  theme(
    plot.title = element_text(angle = 0, vjust = 0.5),
    axis.title.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(linewidth = 2)
  ) 

T_vp$layers[[1]]$aes_params$size = 1.2
T_vp + NoLegend()
    
cairo_pdf("CellType_VlnPlts_combT_Ch17_DEGs.pdf", width = 10, height =  12
  
)
ggarrange(Neut_vp,
          NK_vp,
          mast_vp,
          CD103DC_vp , 
          DC_vp , InfMac_vp , MHCMac_vp , B_vp , T_vp, ncol = 10, nrow = 1, legend = F) 

ggarrange( T_vp,
           B_vp,
          mast_vp, ncol = 10, nrow = 1, legend = F) 
ggsave(file = "CellType_VlnPlts_combT_Ch17_DEGs.pdf",
    units = "px",width = 7000, height = 1200, dpi = 400)


#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 

#Redo Cluster IDs B cells and Plasma Cells -----
#Int.T1D_Timepoints clusters: 
# Cluster 0 Plasma Cells? Jchain Sdc1 Tnfrsf13b Slamf7
# Cluster 1 ILC type 3 Rorc+ Th17 https://doi.org/10.1016/j.cell.2016.07.043 NK? Ncr1 No CD45/Ptprc
# Cluster 2 DCs CD4+ Cd8b1 in NP? DC? Itgax/CD11c+ SiglecH DC-sign(Cd209a/d) Clec9a
# Cluster 3 Naive? B cell Cd19 Ighd Cd22 Cd40 Cd20 Mhc
# Cluster 4 Tgd? cells CD3 CD4 Tcrg-V4 
# Cluster 5 ?? Ighm No CD45 Jchain
# Cluster 6 ILC type 3 Rorc+ Th17 https://doi.org/10.1016/j.cell.2016.07.043
# Cluster 7 ILC type 3 Gata3+ Th2 https://doi.org/10.1016/j.cell.2016.07.043
# Cluster 8 NK Ncr1 Klre1 Klri1 Gmza
# Cluster 9 Mast Cells CD3 MCPT1 Fcer1a No CD45? Kit/CD117+
# Cluster 10 Plasma?? Ighg1 Igha Jchain Sdc1 No CD45
# Cluster 11 Mac?? or Eosinophil? f4/80/Adgre1 (apprently eosinophil have F4/80 in mice) Metalloendopeptidase ADAM-like Decysin 1 (ADAMDEC1) Marker for Inflammatory Bowel Disease https://doi.org/10.3390%2Fijms23095007
# Cluster 12 Treg Cd3 CD4 Foxp3 Il10 Ctla4
# Cluster 13 ?? No CD45 low Jchain
# Cluster 14 Different kind of Treg? Cd3 Cd4 Foxp3 Gzma Gzmb Il10 Cd5
# Cluster 15 ?? No CD45 low Jchain
# Cluster 16 Cd103+ DC Cd86 Clec9a Itgae/Cd103 Mhc
# Cluster 17 Stromal Collagen DCn No CD45
# Cluster 18 Th2 T cells Cd3 Cd4 Il17rb Il4 Gata3 Il13 Il5
# Cluster 19 ? low exp Cd3 Cd4 low Jchain? No CD45
# Cluster 20 ? low exp Cd3 Np Cd4? Jchain No CD45
# Cluster 21 Plasma? Jchain No CD45
# Cluster 22 Plasma? Jchain No CD45
# Cluster 23 Neut?  S100a8/9 has CD45
# Cluster 24 Mac/Mon/eosinophil>??? Mhc Cd14 Apoe Lyz2 Adgre1/f4/80 Itgam Fcgr3/Cd16 
# Cluster 25 Plasma? Jchain Igha No CD45
# Cluster 26 Plasma? Jchain Igha Ly6a/m? No CD45
# Cluster 27 CD8 T cell Cd3 Cd8a Cd8b Tcrg-V7 Ifng
# Cluster 28 B cells Cd19 Cd20 Mhc
# Cluster 29 ILC type 3 Rorc+ Th17 https://doi.org/10.1016/j.cell.2016.07.043
# Cluster 30 Acinar


#old IDs
DotPlot(int.T1D_Timepoints, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
  "Col4a1","Flt1","Nkx2-3", "Tie1", #Stromal
  "S100a8", "S100a9", #Neut
  "Clec9a",  "Itgax", "Cd209a", # DCs /CD11c+ SiglecH DC-sign(Cd209a/d)
  "Itgae",#CD103
  "Itgam", "Adgre1", "Adgre4","Apoe","Cd14","Ms4a7", "Cd68", #Mac
  "Mrc1", "Msr1", "Fcgr3",  "Ccr3", "Cxcl2",  "Cx3cr1", "Il1a", "Il1b","Ccr2", #Mac
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe"  #Plasma 
), group.by = "Cell.types", scale = F
) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism( base_size = 12) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1, size = 8, face = "plain"))


DotPlot(int.T1D_Timepoints, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
  "Col4a1","Flt1","Nkx2-3", "Tie1", #Stromal
  "S100a8", "S100a9", #Neut
  "Clec9a",  "Itgax", "Cd209a", # DCs /CD11c+ SiglecH DC-sign(Cd209a/d)
  "Itgae",#CD103
  "Itgam", "Adgre1", "Adgre4","Apoe","Cd14","Ms4a7", "Cd68", #Mac
  "Mrc1", "Msr1", "Fcgr3",  "Ccr3", "Cxcl2",  "Cx3cr1", "Il1a", "Il1b","Ccr2", #Mac
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Tnfrsf17", "Slamf7" , "Slc3a2", "Slc7a5",
  "Ly6a", "Ly6c1", "Ly6c2", "Ly6k", "Cd28",
  "Igha", "Ighe", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3", "Ighmbp2",  #Plasma
  "Ctrc", "Pnliprp1", "Rnase1", #Acinar
  "Muc2", "Spink4", "ClCa1", "Krt7","Krt8", "Prom1" #Goblet/Gastric Secreting Parietal cells
  
),  scale = F, group.by = "Cell.types2"
) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism(base_family = "Arial", base_size = 12) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1, size = 8, face = "plain"))

ggsave(file = "DotplotMarkergenes_17clusters.png",
       units = "in",width = 14, height = 5, dpi = 400)


Idents(int.T1D_Timepoints) = "seurat_clusters"

Features(int.T1D_Timepoints)[grepl("^Igh", Features(int.T1D_Timepoints))]

int.T1D_Timepoints = RenameIdents(int.T1D_Timepoints, 
                        "0" = "Plasma Cells 1",
                        "1" = "Th17 ILC",
                        "2" = "DCs",
                        "3" = "B Cells",
                        "4" = "CD4 T Cells 1",
                        "5" = "Plasma Cells 2",
                        "6" = "Th17 ILC",
                        "7" = "Th2 ILC",
                        "8" = "NK",
                        "9" = "Mast Cells",
                        "10" = "IgG1+ Plasma Cells",
                        "11" = "Inflam. Mac",
                        "12" = "CD4 T Cells 2",
                        "13" = "Stromal Cells",
                        "14" = "CD4 T Cells 3",
                        "15" = "Stromal Cells",
                        "16" = "CD103+ DCs",
                        "17" = "Stromal Cells",
                        "18" = "CD4 T Cells 4",
                        "19" = "Plasma Cells 4",
                        "20" = "Plasma Cells 5",
                        "21" = "Plasma Cells 6",
                        "22" = "Plasma Cells 7",
                        "23" = "Neutrophil",
                        "24" = "MHCII+ Mac",
                        "25" = "Plasma Cells 8",
                        "26" = "Plasma Cells 9",
                        "27" = "CD8 T Cells",
                        "28" = "IgG1 B Cells",
                        "29" = "Th17 ILC",
                        "30" = "Acinar"
)

int.T1D_Timepoints = RenameIdents(int.T1D_Timepoints, 
                        "0" = "Plasma Cells",
                        "1" = "Th17 ILC",
                        "2" = "DCs",
                        "3" = "B Cells",
                        "4" = "CD4 T Cells",
                        "5" = "Plasma Cells",
                        "6" = "Th17 ILC",
                        "7" = "Th2 ILC",
                        "8" = "NK",
                        "9" = "Mast Cells",
                        "10" = "IgG1+ Plasma Cells",
                        "11" = "Inflam. Mac",
                        "12" = "CD4 T Cells",
                        "13" = "Stromal Cells",
                        "14" = "CD4 T Cells",
                        "15" = "Stromal Cells",
                        "16" = "CD103+ DCs",
                        "17" = "Stromal Cells",
                        "18" = "CD4 T Cells",
                        "19" = "Plasma Cells",
                        "20" = "Plasma Cells",
                        "21" = "Plasma Cells",
                        "22" = "Plasma Cells",
                        "23" = "Neutrophil",
                        "24" = "MHCII+ Mac",
                        "25" = "Plasma Cells",
                        "26" = "Plasma Cells",
                        "27" = "CD8 T Cells",
                        "28" = "IgG1+ B Cells",
                        "29" = "Th17 ILC",
                        "30" = "Acinar"
)

int.T1D_Timepoints[[]]
DimPlot(int.T1D_Timepoints, reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 4, repel = T, 
        group.by = "Cell.types2",
        cols = getPalette(17)) +
  theme_prism(base_size = 12, base_fontface = "bold") + 
  theme(legend.text = element_text(size = 14)) +
  labs( title = "UMAP of Cell Types") 
ggsave(file = "UMAP_17clusters.png",
       units = "in",width = 6, height = 4, dpi = 400)

saveRDS(int.T1D_Timepoints, file = "./updated_17cluster_int.T1D_Timepoints.rds")

CellTypeTable_redo = as.data.frame(proportions(table(Idents(int.T1D_Timepoints),
                                                int.T1D_Timepoints$sample), 
                                          margin = 2))
CellTypeTable_redo$Freq = CellTypeTable_redo$Freq*100
levels(CellTypeTable_redo$Var1)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#display.brewer.all()


CellTypeTable_redo$Var2 = factor(CellTypeTable_redo$Var2, levels = c("PBS PreCh1", "PBS Ch4", "PBS Ch7", "NP Ch4", "NP Ch7"))
#png(file = "CellProportions_24clusters.png",units = "in",width = 14, height = 13, res = 400)
ggplot(CellTypeTable_redo, aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_col(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism( base_size = 10)+
  labs(x = "Condition", y= "Percent of Cells") + 
  scale_fill_manual(values = getPalette(17), name = "Cell Type") + 
  theme(legend.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  )
ggsave(file = "CellProportions_17clusters.png",
       units = "in",width = 5, height = 5.5, dpi = 400)




CellTypeData_redo = CellTypeTable_redo
# load stringr library
library(stringr)
library(ggpubr)

# Split name column into firstname and last name
CellTypeData_redo[c('Tx', 'TP')] <- str_split_fixed(CellTypeData_redo$Var2, ' ', 2)


CellTypeData_redo$Tx = factor(CellTypeData_redo$Tx, levels = c("PBS", "NP"))
CellTypeData_redo$TP = factor(CellTypeData_redo$TP, levels = c("PreCh1", "Ch4", "Ch7"))

myplots = lapply(unique(CellTypeData_redo$Var1), FUN = plotCellType, data = CellTypeData_redo)

png(file = "CellProportions_17clusters_individual_grid.png", units = "in",width = 28, height = 9, res = 400)
ggarrange(plotlist = myplots, ncol=9, nrow=2, 
          common.legend = T, legend = "top")
dev.off()





plotCellTypePresentation(CellTypeData_redo, "B Cells")


myplots = lapply(unique(CellTypeData_redo$Var1), FUN = plotCellTypePresentation, data = CellTypeData_redo)

png(file = "ImmuneCell_prop.png",
    units = "in",width = 17, height = 8, res = 400)
ggarrange(plotlist = myplots, ncol=6, nrow=3, 
          common.legend = T, legend = "right")

ggsave(file = "ImmuneCell_prop_17cluster.png",
       units = "in",width = 17, height = 15
       , dpi = 400)
dev.off()



ggarrange(plotlist = myplots[c(1,9,4,16)], ncol=4, nrow=1, 
          common.legend = T, legend = "right")

ggsave(file = "ImmuneCell_prop_bcells.png",
       units = "in",width = 25, height = 7
       , dpi = 400)

head(int.T1D_Timepoints@meta.data)

int.T1D_Timepoints$Cell.types2 = Idents(int.T1D_Timepoints)
unique(int.T1D_Timepoints$Cell.Type2.TP.Tx)
int.T1D_Timepoints$Cell.Type2.TP.Tx = paste0(int.T1D_Timepoints$Cell.types2, "_",
                                   int.T1D_Timepoints$TP, "_", int.T1D_Timepoints$Tx)

Idents(int.T1D_Timepoints) = "Cell.Type2.TP.Tx"
int.T1D_Timepoints=PrepSCTFindMarkers(int.T1D_Timepoints)

BcellsNPvPBSCh7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 = "B Cells_Ch7_NP",
  ident.2 = "B Cells_Ch7_PBS",
  assay = "SCT",
  test.use = "wilcox"
  
)

IgG1BcellsNPvPBSCh7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 = "IgG1+ B Cells_Ch7_NP",
  ident.2 = "IgG1+ B Cells_Ch7_PBS",
  assay = "SCT",
  test.use = "wilcox"
  
)

PlasmaNPvPBSCh7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 = "Plasma Cells_Ch7_NP",
  ident.2 = "Plasma Cells_Ch7_PBS",
  assay = "SCT",
  test.use = "wilcox"
  
)


IgG1PlasmaNPvPBSCh7 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 = "IgG1+ Plasma Cells_Ch7_NP",
  ident.2 = "IgG1+ Plasma Cells_Ch7_PBS",
  assay = "SCT",
  test.use = "wilcox"
  
)

BcellsNPvPBSCh4 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 = "B Cells_Ch4_NP",
  ident.2 = "B Cells_Ch4_PBS",
  assay = "SCT",
  test.use = "wilcox"
  
)

IgG1BcellsNPvPBSCh4 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 = "IgG1+ B Cells_Ch4_NP",
  ident.2 = "IgG1+ B Cells_Ch4_PBS",
  assay = "SCT",
  test.use = "wilcox"
  
)

PlasmaNPvPBSCh4 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 = "Plasma Cells_Ch4_NP",
  ident.2 = "Plasma Cells_Ch4_PBS",
  assay = "SCT",
  test.use = "wilcox"
  
)


IgG1PlasmaNPvPBSCh4 = FindMarkers(
  int.T1D_Timepoints,
  ident.1 = "IgG1+ Plasma Cells_Ch4_NP",
  ident.2 = "IgG1+ Plasma Cells_Ch4_PBS",
  assay = "SCT",
  test.use = "wilcox"
  
)


#GSEA of B cells Plasma Cells ----
B2_ch7_PBSvNP$symbol = rownames(B2_ch7_PBSvNP)

library(msigdbr)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(forcats)
library(EnhancedVolcano)

all_gene_sets = msigdbr(species = "Mus musculus")
unique(all_gene_sets$gs_subcat)
unique(all_gene_sets$gs_cat)

kegg_gene_sets_msig = msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG")


biocarta_gene_sets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:BIOCARTA")

Reactome_gene_sets = msigdbr(species = "mouse", subcategory = "CP:REACTOME")
Wiki_gene_sets = msigdbr(species = "mouse", subcategory = "CP:WIKIPATHWAYS")

HM_gene_sets = msigdbr(species = "mouse", category = "H")

GO_gene_sets = msigdbr(species = "mouse", subcategory = "GO:BP")

PlasmaNPvPBSCh7$symbol = rownames(PlasmaNPvPBSCh7)

Plasma_ch7_PBSvNP_geneList =  PlasmaNPvPBSCh7 %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(symbol, avg_log2FC)

Plasma_ch7_PBSvNP_PBSvNP_ranks<- deframe(Plasma_ch7_PBSvNP_geneList)


IgG1PlasmaNPvPBSCh7$symbol = rownames(IgG1PlasmaNPvPBSCh7)

G1Plasma_ch7_PBSvNP_geneList =  IgG1PlasmaNPvPBSCh7 %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(symbol, avg_log2FC)

G1Plasma_ch7_PBSvNP_ranks<- deframe(G1Plasma_ch7_PBSvNP_geneList)


BcellsNPvPBSCh7$symbol = rownames(BcellsNPvPBSCh7)

B_ch7_PBSvNP_geneList =  BcellsNPvPBSCh7 %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(symbol, avg_log2FC)

B_ch7_PBSvNP_PBSvNP_ranks<- deframe(B_ch7_PBSvNP_geneList)


IgG1BcellsNPvPBSCh7$symbol = rownames(IgG1BcellsNPvPBSCh7)

G1B_ch7_PBSvNP_geneList =  IgG1BcellsNPvPBSCh7 %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(symbol, avg_log2FC)

G1B_ch7_PBSvNP_PBSvNP_ranks<- deframe(G1B_ch7_PBSvNP_geneList)


combined_data_set = rbind(kegg_gene_sets_msig, 
                          biocarta_gene_sets,
                          Reactome_gene_sets,
                          Wiki_gene_sets,
                          HM_gene_sets,
                          GO_gene_sets)

pathway_gene_set = combined_data_set[,c("gs_name", "gene_symbol")]




y = GSEA(Plasma_ch7_PBSvNP_PBSvNP_ranks,
         pvalueCutoff = 0.5,
         pAdjustMethod = "fdr",
         #OrgDb = org.Mm.eg.db, 
         #TERM2GENE = CTD_gene_sets,
         #TERM2GENE = kegg_gene_sets_msig[,c("gs_description", "entrez_gene")],
         #TERM2GENE = biocarta_gene_sets[,c("gs_description", "entrez_gene")],
         #TERM2GENE = Reactome_gene_sets[,c("gs_description", "entrez_gene")],
         TERM2GENE = pathway_gene_set,
         minGSSize = 1
         
)

y1 = GSEA(G1Plasma_ch7_PBSvNP_ranks,
         pvalueCutoff = 0.5,
         pAdjustMethod = "fdr",
         #OrgDb = org.Mm.eg.db, 
         #TERM2GENE = CTD_gene_sets,
         #TERM2GENE = kegg_gene_sets_msig[,c("gs_description", "entrez_gene")],
         #TERM2GENE = biocarta_gene_sets[,c("gs_description", "entrez_gene")],
         #TERM2GENE = Reactome_gene_sets[,c("gs_description", "entrez_gene")],
         TERM2GENE = pathway_gene_set,
         minGSSize = 1
         
)

y2 = GSEA(B_ch7_PBSvNP_PBSvNP_ranks,
          pvalueCutoff = 0.5,
          pAdjustMethod = "fdr",
          #OrgDb = org.Mm.eg.db, 
          #TERM2GENE = CTD_gene_sets,
          #TERM2GENE = kegg_gene_sets_msig[,c("gs_description", "entrez_gene")],
          #TERM2GENE = biocarta_gene_sets[,c("gs_description", "entrez_gene")],
          #TERM2GENE = Reactome_gene_sets[,c("gs_description", "entrez_gene")],
          TERM2GENE = pathway_gene_set,
          minGSSize = 1
          
)

y3 = GSEA(G1B_ch7_PBSvNP_PBSvNP_ranks,
          pvalueCutoff = 0.5,
          pAdjustMethod = "fdr",
          #OrgDb = org.Mm.eg.db, 
          #TERM2GENE = CTD_gene_sets,
          #TERM2GENE = kegg_gene_sets_msig[,c("gs_description", "entrez_gene")],
          #TERM2GENE = biocarta_gene_sets[,c("gs_description", "entrez_gene")],
          #TERM2GENE = Reactome_gene_sets[,c("gs_description", "entrez_gene")],
          TERM2GENE = pathway_gene_set,
          minGSSize = 1
          
)


head(y)

y_filt = y[which(abs(y$qvalue)<0.25),]
#y_filt = as.data.frame(y)

y_filt$Condition = "PBS"
y_filt[which((y_filt$NES)>0),]$Condition = "NP" 

y_filt$Day = "Challenge 7"
y_filt$Subset = "Plasma cells"

y_filt1 = y1[which(abs(y1$qvalue)<0.25),]
#y_filt = as.data.frame(y)

y_filt1$Condition = "PBS"
y_filt1[which((y_filt1$NES)>0),]$Condition = "NP" 

y_filt1$Day = "Challenge 7"
y_filt1$Subset = "IgG1+ Plasma cells"

y_filt2 = y2[which(abs(y2$qvalue)<0.25),]
#y_filt = as.data.frame(y)

y_filt2$Condition = "PBS"
y_filt2[which((y_filt2$NES)>0),]$Condition = "NP" 

y_filt2$Day = "Challenge 7"
y_filt2$Subset = "B cells"

y_filt3 = y3[which(abs(y2$qvalue)<0.25),]
#y_filt = as.data.frame(y)

y_filt3$Condition = "PBS"
y_filt3[which((y_filt3$NES)>0),]$Condition = "NP" 

y_filt3$Day = "Challenge 7"
y_filt3$Subset = "IgG1+ B cells"

y_filt4 = rbind(y_filt, y_filt1, y_filt2, y_filt3)
y_filt4 = y_filt4[which(y_filt4$p.adjust < 0.1),]

y_filt4 = y_filt4[which(abs(y_filt4$NES)>1.9 & y_filt4$qvalue < 0.05),]

ggplot(y_filt4[which(abs(y_filt4$NES)>1.9 & y_filt4$qvalue < 0.05),], aes(x = Subset, y = fct_reorder(Description, abs(NES)))) + 
  geom_point(aes(color = abs(NES), size = qvalue) ) +
  theme_bw(base_size = 14) +
  scale_colour_gradient2( low = "blue",high="red", mid = "white",na.value = "black") +
  ylab(NULL) +
  ggtitle("GSEA at Challenge 7 in T1D_Timepoints") +
  facet_grid(~Condition) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 16),
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x =element_text(size = 16, face = "bold"),
        axis.title.x = element_blank(),
        strip.text.x  = element_text(size = 16, face = "bold")) 
ggsave(file = "BcellPlasmaCellGSEACh7.png", units = "in",width = 25, height = 20, dpi = 400)

data1 = PlasmaNPvPBSCh7
keyvals <- ifelse(
  (data1$avg_log2FC < -1.5 & data1$p_val_adj <0.05), 'royalblue',
  ifelse((data1$avg_log2FC > 1.5 & data1$p_val_adj < 0.05), 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'black'] <- 'Not DEG'
names(keyvals)[keyvals == 'red'] <- 'PBS'
names(keyvals)[keyvals == 'royalblue'] <- 'NP'

#png(file = "AllergyScaf_NPPBSCh4_VolcanoPlot.png",units = "in",width = 10, height = 10, res = 400)
EnhancedVolcano(data1,
                lab = rownames(data1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colCustom = keyvals,
                selectLab = rownames(data1)[which(names(keyvals) %in% c('NP', 'PBS'))],
                title = "DEGs b/w NP and PBS at Challenge 7",
                drawConnectors = F,
                widthConnectors = 0.75,
                boxedLabels = F)
dev.off()


## ----subcluster T cells-------------------------------------------------------------------------------------------------------------------------------
Idents(int.T1D_Timepoints) = int.T1D_Timepoints$Cell.types2

names(int.T1D_Timepoints@graphs)

int.T1D_Timepoints = FindSubCluster(
  object = int.T1D_Timepoints,
  cluster = "T Cells",
  graph.name = "SCT_nn",
  resolution = 0.5,
  subcluster.name = "T Cells Subtypes",
  algorithm = 2
)

unique(int.T1D_Timepoints$`T Cells Subtypes`)

Idents(int.T1D_Timepoints) = int.T1D_Timepoints$`T Cells Subtypes`

DimPlot(int.T1D_Timepoints, reduction = "umap", split.by = c( "Tx"), combine = F, label = T)

table(Idents(Bcells),Bcells@meta.data$sample)



Idents(int.T1D_Timepoints) = int.T1D_Timepoints$Cell.types2
Tcells = subset(x = int.T1D_Timepoints, idents = c("CD4 T Cells", "CD8 T Cells"))
Tcells <- FindNeighbors(Tcells, reduction = "integrated.dr", dims = 1:30)
Tcells <- FindClusters(Tcells, res = 0.3, graph.name = "SCT_nn")
Tcells <- RunUMAP(Tcells, dims = 1:20, reduction = "integrated.dr")

names(Tcells@graphs)

DimPlot(Tcells, reduction = "umap", split.by = c( "Tx"), combine = F, label = T)


get_conserved <- function(cluster){
  Seurat::FindConservedMarkers(Tcells,
                               ident.1 = cluster,
                               grouping.var = "sample",
                               only.pos = TRUE) %>%
    rownames_to_column(var = "gene")  %>%
    #left_join(y = unique(annotations[, c("gene_name", "description")]),
    #          by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
unique(Idents(Tcells))

conserved_markers.Tcells <- map_dfr(0:5, get_conserved)
#conserved_markers.Tcells[is.na(conserved_markers.Tcells)] <- 0

DefaultAssay(Tcells) = "SCT"
head(conserved_markers.Tcells)

# Extract top 10 markers per cluster
conserved_markers.Tcells <- conserved_markers.Tcells %>% 
  mutate(avg_fc = (`PBS PreCh1_avg_log2FC` + 
                     `PBS Ch4_avg_log2FC` + 
                     `NP Ch4_avg_log2FC` +
                     `PBS Ch7_avg_log2FC` +
                     `NP Ch7_avg_log2FC`) /5)

top20.Tcells <- conserved_markers.Tcells %>% 
  mutate(avg_fc = (`PBS PreCh1_avg_log2FC` + 
                     `PBS Ch4_avg_log2FC` + 
                     `NP Ch4_avg_log2FC` +
                     `PBS Ch7_avg_log2FC` +
                     `NP Ch7_avg_log2FC`) /5) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)

Tcells = PrepSCTFindMarkers(Tcells)
Tcell2v4 = FindMarkers(Tcells, ident.1 = "2", ident.2 = "4",assay= "SCT",
                       recorrect_umi = FALSE )

Tcell2v4R = FindMarkers(Tcells, ident.1 = "2", ident.2 = "4",assay= "RNA",
                       recorrect_umi = FALSE )

Tcell0v3R = FindMarkers(Tcells, ident.1 = "0", ident.2 = "3",assay= "SCT",
                        recorrect_umi = FALSE )

Tcell0v4R = FindMarkers(Tcells, ident.1 = "0", ident.2 = "4",assay= "SCT",
                        recorrect_umi = FALSE )


#DefaultAssay(Tcells) = "RNA"

Tcellmarkers = c(
  "Cd3e","Cd3g", "Cd3d" ,"Cd8a","Cd8b1" ,"Cd4",
  #naive T cells:
  "Ccr7", "Cd28",
  #cytotoxic CD8:
  "Prf1","Gzmb","Gzma","Gzmk","Ccl4","Ccl5",
  #Th1:
  "Ccr4","Cxcr3", "Cxcr4","Ccr5", "Fasl","Trp73","Il12rb1","Il12rb2","Il18r1" , "Il18rap","Ifng","Ifngr1", "Stat1","Tbx21", #Tbet 
  #Th2:
  "Il4","Il4ra", "Il5","Il5ra", "Gata3",
  #Th9:
  "Ccr3", "Ccr6", "Spi1","Il9r",#"Il22ra1",
  #Th17: Ccr6, Ccr4 Nk1.1
  "Il17a",  "Il17f" , "Il17re", "Il17rc", "Il17ra" ,"Il17rd", "Il17rb", "Il17d" ,
  "Klrb1c","Il6ra","Tgfb1","Il23a",
  #Treg: CD127 = Il7r
  "Il7r", "Il2ra", "Ctla4", "Il2", "Foxp3", "Il10","Il10ra",
  #Tfh:
  "Cxcr5", "Cd40lg", "Icos", "Bcl6", "Il21r",
  #gammadelta T cell:
  "Il23r", "Rorc","Rora","Tcrg-V6",
  #Memory:
  "Cd44","Sell", #CD62L
  "Cd27", "Itgae", #Cd103
  "Lag3",
  "Izumo1r", #Fr-4
  "Nt5e", #CD73
  "Slamf6", "Xcl1","Cx3cr1", "Tox", "Pdcd1",
  "Trgv2", "Tcrg-C1", "Tcrg-C2","Tcrg-C4","Trdc", "Trdv4", "Trac", 
  "Klrg1", #CD8 memory
  "Klk1", "Nkg7", "Klrc2", "Ncr1", "Klre1", "Itgam", "Itgax", "Klrb1",
  "Ftl1", "Fth1", "Apoe", "Cxcl10", "Cd274",
  "Mki67", "Tuba1b" #dividing
)

DotPlot(Tcells, features = Tcellmarkers, scale = F, assay = "RNA") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  
  theme_prism(base_family = "Arial", base_size = 12) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1))
dev.off()

table(Tcell_data$seurat_clusters, Tcell_data$age_day)

rownames(Tcell_data)
rownames(Tcells)[grepl("^Gzm", rownames(Tcells))]
rownames(Tcells)[grepl("^Fas", rownames(Tcells))]

rownames(Tcell_data)[grepl("^Tcr", rownames(Tcell_data))]
rownames(Tcell_data)[grepl("^Trac", rownames(Tcell_data))]

rownames(Tcell_data)[grepl("^Tcrg", rownames(Tcell_data))]

table(Idents(Tcells), Tcells$sample)

VlnPlot(Tcell_data, features = c("Cd3e","Cd3g", "Cd3d","Cd8a","Cd4" ,"Foxp3", "Itgam", "Itgax",
                                 "Ftl1", "Fth1" ), split.by = "age_day")

VlnPlot(Tcell_data, features = c("Il17a","Il23r", "Rorc","Rora","Tcrg-V6", "Trgv2", "Trdc", "Trdv4"
), split.by = "age_day")


VlnPlot(Tcell_data, features = c("Cd8a","Cd8b1","Cd4" ,"Tcrg-V4", "Tcrg-V6", "Tcrg-C1", "Tcrg-C2", "Tcrg-V1", "Tcrg-C4", 
                                 "Tcrg-V5"), split.by = "age_day")

VlnPlot(Tcells, features = c("Lag3", "Izumo1r","Cd274"), split.by = "sample")

VlnPlot(Tcells, features = c("Cd4", "Cd8a","Gzmk","Nkg7", "Ccl5", "Ccr9", "Itga4"), split.by = "sample") +
  theme(legend.position = 'right')

VlnPlot(Tcells, features = c("Cd4","Il12rb1", "Cxcr3", "Tbx21", "Il18r1", "Ctsw", "Cxcr6", "Il2rb", "Pde7a"), split.by = "sample") +
  theme(legend.position = 'right')



VlnPlot(Tcell_data, features = c("Cd8a","Cd8b1","Cd4", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa"),
        split.by = "age_day")

VlnPlot(Tcell_data, features = c("Cd40lg", "Icos", "Bcl6", "Il21r","Itgam", "Itgax",
                                 "Ccr6", "Cxcr3", "Cxcr5"), split.by = "age_day")

cells.use <- WhichCells(Tcell_data, expression=`H2-Eb1` > 0 |`Cd4` > 0 , idents = c("5"))

table(Tcell_data$integrated_snn_res.0.3, Tcell_data$age_day)

# Cluster 0 is Il18r1 Ctsw CD4
# Cluster 1 is Treg CD4+ Foxp3+ IL10+
# Cluster 2 is Gzma/b+ CD4 T Gzma/b
# Cluster 3 is Naive CD4 T Ccr7+ Sell+
# Cluster 4 is Th2 T 
# Cluster 5 is CD8+ T Ccl5+ Ctsw


# Cluster 0 is CD4 effector memory Cd44+ CD62L- (Sell-)
# Cluster 1 is CD4 effector memory Cd44+ CD62L- (Sell-)
# Cluster 2 is CD4 effector memory Cd44+ CD62L- (Sell-)
# Cluster 3 is Cytotoxic CD8 T Gmza/b+Ccl5+
# Cluster 4 is Naive CD8 T Ccr7 
# Cluster 5 is CD4 CD11b
# Cluster 6 are CD4 CD11b
# Cluster 7 are gd CD4- Cd8- Il17a "Tcrg-C1", "Tcrg-C2","Tcrg-C4", "Trdc", "Trdv4",
# Cluster 8 are Treg CD4+ Foxp3+ Ctla4+

Tcells = RenameIdents(Tcells, 
                          "0" = "Il18r1 CD4 T",
                          "1" = "Treg",
                          "2" = "Gzma/b+ CD4 T",
                          "3" = "Naive CD4 T",
                          "4" = "Th2 T",
                          "5" = "CD8 T"
)

Tcells$TcellIDs = Idents(Tcells)
int.T1D_Timepoints[[]]
head(int.T1D_Timepoints)
# Generate a new column called sub_cluster in the metadata
int.T1D_Timepoints$Cell.types3 <- as.character(Idents(int.T1D_Timepoints))

# Change the information of cells containing sub-cluster information
int.T1D_Timepoints$Cell.types3[Cells(Tcells)] <- as.character(Idents(Tcells))

DimPlot(int.T1D_Timepoints, reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 6, repel = T, 
        group.by = "Cell.types3",
        cols = getPalette(21)) +
  theme_prism(base_size = 12, base_fontface = "bold") + 
  theme(legend.text = element_text(size = 14)) +
  labs( title = "UMAP of Cell Types") 
length(unique(int.T1D_Timepoints$Cell.types3))

saveRDS(int.T1D_Timepoints, file = "./updated_21cluster_int.T1D_Timepoints.rds")


# subcluster B Cells real ----
levels(Idents(int.T1D_Timepoints))
Bcells = subset(x = int.T1D_Timepoints, idents =  c("B Cells", "IgG1+ B Cells"))
length(Bcells$orig.ident)
ElbowPlot(Tcells)

Bcells <- FindNeighbors(Bcells, reduction = "integrated.dr", dims = 1:20)
Bcells <- FindClusters(Bcells, res = 0.3, graph.name = "SCT_nn")
Bcells <- RunUMAP(Bcells, dims = 1:20, reduction = "integrated.dr")

DimPlot(Bcells, reduction = "umap", split.by = c( "sample"), combine = F, label = T)

table(Idents(Bcells), Bcells$sample)

get_conserved <- function(cluster){
  Seurat::FindConservedMarkers(Bcells,
                               ident.1 = cluster,
                               grouping.var = "sample",
                               only.pos = TRUE) %>%
    rownames_to_column(var = "gene")  %>%
    #left_join(y = unique(annotations[, c("gene_name", "description")]),
    #          by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
unique(Idents(Bcells))

conserved_markers.Bcells <- map_dfr(0:5, get_conserved)
conserved_markers.Bcells[is.na(conserved_markers.Bcells)] <- 0

DefaultAssay(Bcells) = "SCT"
head(conserved_markers.Bcells)
conserved_markers$`PBS PreCh1_avg_log2FC`


# Extract top 10 markers per cluster
conserved_markers.Bcells <- conserved_markers.Bcells %>% 
  mutate(avg_fc = (`PBS PreCh1_avg_log2FC` + 
                     `PBS Ch4_avg_log2FC` + 
                     `NP Ch4_avg_log2FC` +
                     `PBS Ch7_avg_log2FC` +
                     `NP Ch7_avg_log2FC`) /5)



top20.Bcells <- conserved_markers.Bcells %>% 
  mutate(avg_fc = (`PBS PreCh1_avg_log2FC` + 
                     `PBS Ch4_avg_log2FC` + 
                     `NP Ch4_avg_log2FC` +
                     `PBS Ch7_avg_log2FC` +
                     `NP Ch7_avg_log2FC`) /5) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)


DotPlot(Bcells, features =c(
  "Ptprc", #CD45
  
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Cd3e", "Trac",
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #Myeloid
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Cd74", "Cd44", "Cxcr4", "Bach2", "Sell", "Ltb", "Prdm1", #"Blimp1",
  "Il4","Il4ra", "Il5","Il5ra", "Il1a", "Il1b", "Il6","Ctla4",
  "Havcr1","Tnfrsf10b",
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3"  #Plasma 
), scale = F, assay = "SCT"
) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")


DotPlot(Bcells, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #Myeloid
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3"  #Plasma 
), scale = F
) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") 


DotPlot(Bcells, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
  "Col4a1","Flt1","Nkx2-3", "Tie1", #Stromal
  "S100a8", "S100a9", #Neut
  "Clec9a",  "Itgax", "Cd209a", # DCs /CD11c+ SiglecH DC-sign(Cd209a/d)
  "Itgae",#CD103
  "Itgam", "Adgre1", "Adgre4","Apoe","Cd14","Ms4a7", "Cd68", #Mac
  "Mrc1", "Msr1", "Fcgr3",  "Ccr3", "Cxcl2",  "Cx3cr1", "Il1a", "Il1b","Ccr2", #Mac
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe"  #Plasma 
), scale = F, assay = "RNA"
) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism(base_family = "Arial", base_size = 12) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1, size = 8, face = "plain"))


colnames(Bcells[[]])

Idents(Bcells)

VlnPlot(Bcells, features = c("Cd19", "Il1a", "Pax5", "Cd80", "Cd86", "Spi1", "Spn", "Sdc1", "Cd1d1", "Itgam", "Ighd"), 
        split.by = "sample", assay = "RNA") +
  theme(legend.position = 'right')

VlnPlot(int.T1D_Timepoints, features = c("Cd19","Ighd", "Ighm", "Ighg1", "Ighg2b", "Cr2", "Fcer2a", "Fcrl1"), 
        split.by = "sample", assay = "RNA") +
  theme(legend.position = 'right')


Bcells = PrepSCTFindMarkers(Bcells)
Bcell0v1 = FindMarkers(Bcells, ident.1 = "0", ident.2 = "1",assay= "SCT",
                       recorrect_umi = FALSE )


Bcells = RenameIdents(Bcells, 

                      "5" = "2",
                      "4" = "3"
)



