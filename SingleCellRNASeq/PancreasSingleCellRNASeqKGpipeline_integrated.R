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
data_dir = "113210_JR_pool/"

sample.names = list.files(path = '113210_JR_pool/')
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
Week18_2 =  remove_soup(Week18_2, data_dir, sample_dir_list[10])
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
Week18_2 =  filter_ddqcr(Week18_2, sample_list[2], sample_dir_list[10])
Week18_3 =  filter_ddqcr(Week18_3, sample_list[3], sample_dir_list[10])
Week18_4 =  filter_ddqcr(Week18_4, sample_list[4], sample_dir_list[11])
Week18_5 =  filter_ddqcr(Week18_5, sample_list[5], sample_dir_list[12])
Week18_6 =  filter_ddqcr(Week18_6, sample_list[6], sample_dir_list[13])

Week12_1 = filter_ddqcr(Week12_1, sample_list[7], sample_dir_list[14])
Week12_2 =  filter_ddqcr(Week12_2, sample_list[8], sample_dir_list[15])
Week12_3 =  filter_ddqcr(Week12_3, sample_list[10], sample_dir_list[16])
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
Week18_2 =  find_doublets(Week18_2, sample_list[2], sample_dir_list[10])
Week18_3 =  find_doublets(Week18_3, sample_list[3], sample_dir_list[10])
Week18_4 =  find_doublets(Week18_4, sample_list[4], sample_dir_list[11])
Week18_5 =  find_doublets(Week18_5, sample_list[5], sample_dir_list[12])
Week18_6 =  find_doublets(Week18_6, sample_list[6], sample_dir_list[13])

Week12_1 = find_doublets(Week12_1, sample_list[7], sample_dir_list[14])
Week12_2 =  find_doublets(Week12_2, sample_list[8], sample_dir_list[15])
Week12_3 =  find_doublets(Week12_3, sample_list[10], sample_dir_list[16])
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

getPalette = colorRampPalette(brewer.pal(10, "Set1"))

#png(file = "UMAP_31clusters.png",    units = "in",width = 11, height = 11, res = 400)
# DimPlot(int.T1D_Timepoints, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 10, repel = T, 
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




## ---- Annotate cluster ID---------------------------------------------------------------------------------------------------------------------------------------



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
# "Cd710a", "Cd110" ,"Ms4a1", "Ighm", 
# "Cd300a", "Ly6c2", "Cd14", "Ccr2", 
# "Smox", "Apoe", "Ms4a7", "Ccr5", "Fcgr1", "Adgre1", "Cd68",
# "Ftl1", "Fth1", 
# "Tmem1110", "Trem2", "Aif1", "Csf1r",
# "S100a10", "Camp", "Retnlg", "Ly6g",
# "Ifitm1", "Siglech", "Klk1", 
# "Zbtb46", "Itgax", "Cd86", "Cd2010a",
# "Bst2", "Cmah", "Ly6a", "Nrp1", "Clec4g", 
# "Cacnb3", "Fscn1", "Syn3", "Tmem150c", 
# "Snca", "Ube2o"
#"Ncr1",  "Klre1",#NK
# Extract top genes for cluster 18
cluster_18_genes <- top20_sct %>% 
  filter(cluster == 18) %>% 
  pull(gene)

# Extract top genes for cluster 110
cluster_110_genes <- top20_sct %>% 
  filter(cluster == 110) %>% 
  pull(gene)

# Print genes for cluster 18
cat("Top 20 genes for Cluster 18:\n")
print(cluster_18_genes)

# Print genes for cluster 110
cat("\nTop 20 genes for Cluster 110:\n")
print(cluster_110_genes)



DotPlot(T1D_Timepoints, features = c( "Ptprc", #CD45
                                      "Trbc1","Trbc2","Cd3g","Cd3e","Cd3d","Lat", #T Cells
                                      "C1qa","C1qb","C1qc","Lyz2","Adgre1", #Macrophages
                                      "Cd710a", "Cd710b","Ms4a1","Cd110",#B Cells
                                      "S100a4","S100a11","Itgax","Zbtb46","Flt3","H2-Aa","Cd74","Mdh2","Gm2a","Pglyrp1","H2-Oa",#cDc
                                      "Siglech","Cd2010a","Cox6a2","Rnase6",#pDC
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
                                  "10" = "T Cell",
                                  "10" = "T Cell",
                                  "11" = "Macrophage",
                                  "12" = "T Cell",
                                  "13" = "T Cell",
                                  "14" = "T Cell",
                                  "15" = "T Cell",
                                  "16" = "T Cell",
                                  "17" = "Acinar Cell",
                                  "18" = "Mesenchymal-like", #TBD
                                  "110" = "Dendritic Cell", #TBD
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

T1D_Timepoints = readRDS("./annotated_T1D_Timepoints_v3.rds")
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

ggsave(file = "CellProportions_10clusters.png",
    units = "in",width = 11, height = 11, dpi = 400)
dev.off()



## ----table of total cells per sample------------------------------------------------------------------------------------------------------------------
table(T1D_Timepoints@meta.data$sample,T1D_Timepoints@meta.data$sample)

#' 
## ----more cluster ids---------------------------------------------------------------------------------------------------------------------------------


## ----Cluster ID dotplot-------------------------------------------------------------------------------------------------------------------------------
#Cluster IDs
DimPlot(T1D_Timepoints)


#' #Subcluster T cells
## ----Subcluster T cell-------------------------------------------------------------------------------------------------------------------------------

names(int.T1D_Timepoints@graphs)

T1D_Timepoints = FindSubCluster(
  object = T1D_Timepoints,
  cluster = "T Cell",
  graph.name = "SCT_nn",
  resolution = 0.5,
  subcluster.name = "T Cell Subtypes",
  algorithm = 2
)

unique(T1D_Timepoints$`T Cell Subtypes`)

#Idents(T1D_Timepoints) = T1D_Timepoints$`T Cell Subtypes`

#DimPlot(T1D_Timepoints, reduction = "umap", group.by = , combine = F, label = T)



Tcells = subset(x = T1D_Timepoints, idents = "T Cell")
ElbowPlot(Tcells, ndims = 50)
Tcells <- RunUMAP(Tcells, dims = 1:30)
Tcells <- FindNeighbors(Tcells, dims = 1:30, verbose = F)
Tcells <- FindClusters(Tcells, res = 0.7)#graph.name = "SCT_nn"


DimPlot(Tcells, reduction = "umap", combine = F, label = T)


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

conserved_markers.Tcells <- map_dfr(0:22, get_conserved)
conserved_markers.Tcells[is.na(conserved_markers.Tcells)] <- 0
head(conserved_markers.Tcells)

DefaultAssay(Tcells) = "SCT"
conserved_markers.Tcells$Week6_3_avg_log2FC
# Extract top 10 markers per cluster
top20_TCells <- conserved_markers.Tcells %>% 
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


Tcells <- PrepSCTFindMarkers(Tcells)
all.markers <- FindAllMarkers(Tcells , only.pos = TRUE,recorrect_umi = FALSE)

top20_sct=all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)

# Save the top 20 markers grouped by cluster as a CSV file
write.csv(top20_sct, file = "top20_sct_Tcells_T1D_Timepoints.csv", row.names = FALSE)


DotPlot(Tcells,assay = "SCT", features = c( "Tnfrsf4","Tbc1d4","Cd4", #CD4
                                      "Cd8a","Blk",#CD8
                                      "Ccr7","Il7r","Lef1","Sell",#Memory
                                      "Gzma","Gzmb","Klrg1","Tbx21","Zeb2",#SLEC
                                      "Ccl5", "Gzmk", "Klrc1", "Nkg7", "Cd38", "Cxcr3", "Cxcr6", "Fasl", "Ifng",#CD8 Effector
                                      "Cd40lg","Il6ra",#Tcon effector
                                      "Nt5e","Izumo1r",#Anergy
                                      "Slamf6","Tcf7",#Progenitor
                                      "Cd200","Il21","Tnfsf8", #Th21
                                      "Rora","Il17f", "Stat3", "Ccr6",#Th17
                                      "Bcl6","Cxcr5","Tox2",#Tfh
                                      "Ctla4","Foxp3","Icos","Ikzf2","Il2ra","Tigit",#Tregs
                                      "Cd610","Egr1","Egr2","Nr4a1",#Activation
                                      "Ifit1","Ifit3","Isg15","Stat1",#Interferon sensing
                                      "Mki67","Stmn1",#Proliferation
                                      "Tox","Pdcd1","Lag3",#Exhaustion
                                      "Klra8","Klrb1a","Eomes","Ncr1","Klre1", #NK
                                      "Tcrg-C2",#Gamm_Delta-,"Trac","Cd3e"
                                      "Hs3st1","Gata3",#ILC2
                                      "Hebp1","Ncoa7","Cd74","Cited4","H2-Ab1","Rorc",#ILC3

                                     
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

# Check average expression of Nkg7 per cluster
AverageExpression(Tcells, features = "Cd8a", group.by = "seurat_clusters")

# Define the T cell subtypes for each cluster
TCellSubtype <- c("CD4 T Cell",  # Tcon central memory
                  "CD8 T Cell",  # CD8 central memory
                  "CD4 T Cell",  # Tcon effector~TBD
                  "CD4 T Cell",  # Tcon progenitor like ~TBD
                  "CD4 T Cell",  # Tregs
                  "CD4 T Cell",  # Tcon central memory /memory~TBD
                  "CD4 T Cell",  # Tcon Recently Activated 
                  "CD4 T Cell",  # Tcon central memory
                  "CD4 T Cell",  # Tcon central memory
                  "CD8 T Cell",  # CD8 central memory
                  "CD8 T Cell",  # CD8 central memory
                  "NK Cell",     # NK Cell
                  "CD4 T Cell",  # Tcon central memory/effector/progenitor~TBD
                  "CD4 T Cell",  # Tcon central memory
                  "CD8 T Cell",  # CD8 central memory
                  "CD8 T Cell",  # CD8 central memory
                  "Gamma Delta T Cell",  # Gamma Delta T Cell
                  "ILC2",        # ILC2
                  "CD8 T Cell",  # CD8 central memory
                  "CD4 T Cell",  # Tcon Interferon Sensing
                  "CD8 T Cell",  # CD8 central memory
                  "CD8 T Cell",  # CD8 effector like SLEC
                  "ILC3")        # ILC3

# Ensure the Idents are in numeric format
ident_integers <- as.integer(Idents(Tcells))

# Check if there are any unexpected levels
unique(ident_integers)

# Add the T cell subtype information to the metadata of the Seurat object
Tcells$TCellType <- TCellSubtype[as.integer(Idents(Tcells))]

# Check the updated metadata
head(Tcells@meta.data)
DimPlot(Tcells,reduction = 'umap',group.by = "TCellType" ,label = T)

### DEGs for Tcon memory cluster ---------------

# Extract top genes for cluster 0
cluster_0_genes <- top20_sct  %>% 
  filter(cluster == 0) %>% 
  pull(gene)
# Print genes for cluster 0
cat("Top 20 genes for Cluster 0:\n")
print(cluster_0_genes)

# Extract top genes for cluster 3
cluster_3_genes <- top20_sct  %>% 
  filter(cluster == 3) %>% 
  pull(gene)
# Print genes for cluster 3
cat("Top 20 genes for Cluster 3:\n")
print(cluster_3_genes)

# Extract top genes for cluster 5
cluster_5_genes <- top20_sct  %>% 
  filter(cluster == 5) %>% 
  pull(gene)
# Print genes for cluster 5
cat("Top 20 genes for Cluster 5:\n")
print(cluster_5_genes)
cluster5<- top20_sct  %>% 
  filter(cluster == 5) 

# Extract top genes for cluster 7
cluster_7_genes <- top20_sct  %>% 
  filter(cluster == 7) %>% 
  pull(gene)
# Print genes for cluster 7
cat("Top 20 genes for Cluster 7:\n")
print(cluster_7_genes)
cluster7<- top20_sct  %>% 
  filter(cluster == 7) 
cluster7

# Extract top genes for cluster 8
cluster_8_genes <- top20_sct  %>% 
  filter(cluster == 8) %>% 
  pull(gene)
# Print genes for cluster 8
cat("Top 20 genes for Cluster 8:\n")
print(cluster_8_genes)
cluster8<- top20_sct  %>% 
  filter(cluster == 8) 
cluster8

# Extract top genes for cluster 12
cluster_12_genes <- top20_sct  %>% 
  filter(cluster == 12) %>% 
  pull(gene)
# Print genes for cluster 12
cat("Top 20 genes for Cluster 12:\n")
print(cluster_12_genes)
cluster12<- top20_sct  %>% 
  filter(cluster == 12) 
cluster12

# Extract top genes for cluster 13
cluster_13_genes <- top20_sct  %>% 
  filter(cluster == 13) %>% 
  pull(gene)
# Print genes for cluster 13
cat("Top 20 genes for Cluster 13:\n")
print(cluster_13_genes)
cluster13<- top20_sct  %>% 
  filter(cluster == 13) 
cluster13
#Dotplot for checking CD4 Tcon memory subtypes
DotPlot(Tcells, features = c( "Tnfrsf4","Tbc1d4","Cd4", #CD4
                              "Cd8a","Blk",#CD8
                              "Ccr7","Il7r","Lef1","Sell",#Memory
                              "Igfbp4","Zbtb7b","Acvrl1","St8sia6","Myo10",#0
                              "Crlf3","Selenop","Tent5c",#3
                              "Trbv2",#5
                              "Trbv4",#7
                              "Trbv19",#8
                              "Trbv1","Pdlim4",#12
                              "Trbv20"#13
                              
                              
                              
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")


### DEGs for CD8 memory cluster ---------------
cluster_1_genes <- top20_sct  %>% 
  filter(cluster == 1) %>% 
  pull(gene)
# Print genes for cluster 1
cat("Top 20 genes for Cluster 1:\n")
print(cluster_1_genes)
cluster1<- top20_sct  %>% 
  filter(cluster == 1) 
cluster1

cluster_10_genes <- top20_sct  %>% 
  filter(cluster == 10) %>% 
  pull(gene)
# Print genes for cluster 10
cat("Top 20 genes for Cluster 10:\n")
print(cluster_10_genes)
cluster10<- top20_sct  %>% 
  filter(cluster == 10) 
cluster10

cluster_10_genes <- top20_sct  %>% 
  filter(cluster == 10) %>% 
  pull(gene)
# Print genes for cluster 10
cat("Top 20 genes for Cluster 10:\n")
print(cluster_10_genes)
cluster10<- top20_sct  %>% 
  filter(cluster == 10) 
cluster10

cluster_14_genes <- top20_sct  %>% 
  filter(cluster == 14) %>% 
  pull(gene)
# Print genes for cluster 14
cat("Top 20 genes for Cluster 14:\n")
print(cluster_14_genes)
cluster14<- top20_sct  %>% 
  filter(cluster == 14) 
cluster14

cluster_15_genes <- top20_sct  %>% 
  filter(cluster == 15) %>% 
  pull(gene)
# Print genes for cluster 15
cat("Top 20 genes for Cluster 15:\n")
print(cluster_15_genes)
cluster15<- top20_sct  %>% 
  filter(cluster == 15) 
cluster15

cluster_18_genes <- top20_sct  %>% 
  filter(cluster == 18) %>% 
  pull(gene)
# Print genes for cluster 18
cat("Top 20 genes for Cluster 18:\n")
print(cluster_18_genes)
cluster18<- top20_sct  %>% 
  filter(cluster == 18) 
cluster18

cluster_20_genes <- top20_sct  %>% 
  filter(cluster == 20) %>% 
  pull(gene)
# Print genes for cluster 20
cat("Top 20 genes for Cluster 20:\n")
print(cluster_20_genes)
cluster20<- top20_sct  %>% 
  filter(cluster == 20) 
cluster20

#Dotplot for checking CD8 memory subtypes
DotPlot(Tcells, features = c( "Tnfrsf4","Tbc1d4","Cd4", #CD4
                                            "Cd8a","Blk",#CD8
                                            "Ccr7","Il7r","Lef1","Sell",#Memory
                                            "Itgae","Ccr10","Adgrg5",#1
                                            "Trbv14",#9
                                            "Trbv29",#10
                                            "Trbv19",#14
                                            "Trbv17", #15 
                                            "Trbv2"#18
                              
                                            
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")


# Define the T cell subtypes for each cluster
TCellDetailedSubtype <- c("Tcon memory",#0
                  "CD8 memory", #1=Tissue Resident 
                  "Tcon exhausted effector-like",  
                  "Tcon memory",#3
                  "Tregs",  # 4
                  "Tcon memory",  #5 Antigen Specific Trbv2 Enriched
                  "Tcon activated ",  #6 
                  "Tcon memory",  # Antigen Specific Trbv4 Enriched
                  "Tcon memory",  # 8-Antigen Specific Trbv19 Enriched
                  "CD8 memory",  # 9- Antigen Specific Trbv14 Enriched
                  "CD8 memory",  # 10- Antigen Specific Trbv29 Enriched
                  "NK Cell",    
                  "Tcon memory",  # Antigen Specific Trbv1 Enriched
                  "Tcon memory",  # Antigen Specific Trbv20 Enriched
                  "CD8 memory",  # 14-Antigen Specific Trbv19 Enriched
                  "CD8 memory",  # 15-Antigen Specific Trbv17 Enriched
                  "Gamma Delta T Cell",  
                  "ILC2",        
                  "CD8 memory",  # 18-Antigen Specific Trbv2 Enriched
                  "Tcon Interferon Sensing",  #19 
                  "CD8 memory",  # 20-Antigen Specific Trbv1 Enriched
                  "CD8 exhausted effector-like",  #21 
                  "ILC3")        

# Ensure the Idents are in numeric format
ident_integers <- as.integer(Idents(Tcells))

# Check if there are any unexpected levels
unique(ident_integers)

# Add the T cell subtype information to the metadata of the Seurat object
Tcells$TCellSubType <- TCellDetailedSubtype[as.integer(Idents(Tcells))]

DimPlot(Tcells,reduction = 'umap',group.by = "TCellSubType" ,label = T)


# Define the T cell subtypes for each cluster
TCellDetailedMemorySubtype <- c("Tcon memory-1",#0
                          "CD8 memory-1", #1=Tissue Resident 
                          "Tcon exhausted effector-like",  
                          "Tcon memory-1",#3
                          "Tregs",  # 4
                          "Tcon memory-2(Trbv2+)",  #5 Antigen Specific Trbv2 Enriched
                          "Tcon activated ",  #6 
                          "Tcon memory-3(Trbv4+)",  # Antigen Specific Trbv4 Enriched
                          "Tcon memory-4(Trbv19+)",  # 8-Antigen Specific Trbv19 Enriched
                          "CD8 memory-2(Trbv14+)",  # 9- Antigen Specific Trbv14 Enriched
                          "CD8 memory-3(Trbv29+)",  # 10- Antigen Specific Trbv29 Enriched
                          "NK Cell",    
                          "Tcon memory-5(Trbv1+)",  # Antigen Specific Trbv1 Enriched
                          "Tcon memory-6(Trbv20+)",  # Antigen Specific Trbv20 Enriched
                          "CD8 memory-4(Trbv19+)",  # 14-Antigen Specific Trbv19 Enriched
                          "CD8 memory-5(Trbv17+)",  # 15-Antigen Specific Trbv17 Enriched
                          "Gamma Delta T Cell",  
                          "ILC2",        
                          "CD8 memory-6(Trbv2+)",  # 18-Antigen Specific Trbv2 Enriched
                          "Tcon Interferon Sensing",  #19 
                          "CD8 memory-7(Trbv1+)",  # 20-Antigen Specific Trbv1 Enriched
                          "CD8 exhausted effector-like",  #21 
                          "ILC3")        

# Ensure the Idents are in numeric format
ident_integers <- as.integer(Idents(Tcells))

# Check if there are any unexpected levels
unique(ident_integers)

# Add the T cell subtype information to the metadata of the Seurat object
Tcells$TCellMemorySubType <- TCellDetailedMemorySubtype[as.integer(Idents(Tcells))]

DimPlot(Tcells,reduction = 'umap',group.by = "TCellMemorySubType" ,label = T)

Tcells$Idents<-Idents(Tcells)

unique(Tcells$Idents)

Tcells = RenameIdents(Tcells, 
                              "0" = "CD4 T Cell",#Tcon central memory
                              "1" = "CD8 T Cell",#CD8 central memory
                              "2" = "CD4 T Cell",#Tcon effector~TBD
                              "3" = "CD4 T Cell",#Tcon progenitor like ~TBD
                              "4" = "CD4 T Cell",#Tregs
                              "5" = "CD4 T Cell",#Tcon central memory /memory~TBD
                              "6" = "CD4 T Cell",#Tcon Recently Activated 
                              "7" = "CD4 T Cell",#Tcon central memory
                              "8" = "CD4 T Cell",#Tcon central memory
                              "10" = "CD8 T Cell",#CD8 central memory
                              "10" = "CD8 T Cell",#CD8 central memory
                              "11" = "NK Cell",
                              "12" = "CD4 T Cell",#Tcon central memory/effector/progenitor~TBD
                              "13" = "CD4 T Cell",#Tcon central memory
                              "14" = "CD8 T Cell",#CD8 central memory
                              "15" = "CD8 T Cell",#CD8 central memory
                              "16" = "Gamma Delta T Cell",
                              "17" = "ILC2",
                              "18" = "CD8 T Cell", #CD8 central memory
                              "19" = "CD4 T Cell", #Tcon Interferon Sensing
                              "20" = "CD8 T Cell",#CD8 central memory
                              "21" = "CD8 T Cell",#CD8 effector like SLEC
                              "22" = "ILC3")

DimPlot(Tcells, reduction = "umap", 
        #split.by = c( "Tx"), 
        combine = F, label = T)
saveRDS(Tcells, file = "./annotated_TCells_T1D_Timepoints_v1.rds")

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

y_filt3 = y_filt3[which(abs(y_filt3$NES)>1.10),]

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
# T1D_Timepoints$Cell.types <-Idents(T1D_Timepoints)
# 
# 
# # pseudobulk cells by stimulation condition AND cell type AND donor
# bulk <- AggregateExpression(T1D_Timepoints, group.by = c("time", "Cell.types", "group"),
#                             assays = "RNA",
#                             features = NULL,
#                             return.seurat = FALSE,
#                             #group.by = "Cell.types",
#                             add.ident = NULL,
#                             slot = "counts",
#                             verbose = TRUE)
# 
# pseudobulk_PBS1 = AggregateExpression(
#   subset(T1D_Timepoints, idents = c("PBS PreCh1")),
#   assays = "RNA",
#   features = NULL,
#   return.seurat = FALSE,
#   group.by = "Cell.types",
#   add.ident = NULL,
#   slot = "counts",
#   verbose = TRUE
# )
# 
# Idents(TNBC) = "Day"
# pseudobulkD10 = AggregateExpression(
#   subset(TNBC,idents = c("D10")),
#   assays = "RNA",
#   features = NULL,
#   return.seurat = FALSE,
#   group.by = "new.cluster.ids",
#   add.ident = NULL,
#   slot = "counts",
#   verbose = TRUE
# )
# Idents(TNBC) = "new.cluster.ids"
# head(pseudobulkD10)
# 
# Idents(TNBC) = "Day"
# pseudobulkD7 = AggregateExpression(
#   subset(TNBC,idents = c("D7")),
#   assays = "RNA",
#   features = NULL,
#   return.seurat = FALSE,
#   group.by = "new.cluster.ids",
#   add.ident = NULL,
#   slot = "counts",
#   verbose = TRUE
# )
# Idents(TNBC) = "new.cluster.ids"
# head(pseudobulkD7)
# 
# 
# signature_genes1 = read.csv("/Users/lailarad/Documents/BI_EAE/aaron_pEAE/Signature_genes1.csv",row.names = NULL)
# signature_genesTx = read.csv("/Users/lailarad/Documents/BI_EAE/Cohort1_RNA_seq/signature_genes_VLA4tx.csv",row.names = NULL)
# signature_genes_day10_pEAE_bulk = read.csv("/Users/lailarad/Documents/BI_EAE/aaron_pEAE/aaronpEAE_signature_genes_day10.csv",row.names = NULL)
# 
# 
# 
# 
# Features(int.T1D_Timepoints)
# Features(int.T1D_Timepoints)[grepl("^H2-", Features(int.T1D_Timepoints))]
# int.T1D_Timepoints[["RNA"]]$counts[grepl("^H2-", Features(int.T1D_Timepoints)),]
# mat[grepl("^H2-", Features(int.T1D_Timepoints)),]
# 
# pseudobulkset = pseudobulk_PBS1
# mat = pseudobulkset$RNA
# mat=log(mat + 1, base = 2)
# 
# head(mat)
# 
# pseudobulkset = bulk
# mat = pseudobulkset$RNA
# mat=log(mat + 1, base = 2)
# mat = mat - rowMeans(mat)
# 
# 
# 
# max(mat)
# 
# 
# 
# sigGenes = (intersect(rownames(mat), signature_genes1$x))
# sigGenes = (intersect(rownames(mat), signature_genes_day10_pEAE_bulk$x))
# sigGenes = (intersect(rownames(mat), signature_genesTx$genes))
# 
# pheatmap::pheatmap(t(mat),
#                    #annotation_col = anno,
#                    show_colnames = T,
#                    fontsize_row = 18,
#                    fontsize_col = 12,
#                    cluster_cols = T,
#                    cluster_rows = T,
#                    #annotation_colors = ann_colors,
#                    #filename = paste(path, "signatureCelltypesD710_vla4tx_log2.png",sep=""),
#                    width = 8,
#                    height = 5
#                    
#                    
# )
# dev.off()
# getwd()
# plots <- VlnPlot(TNBC, features = sigGenes[(25:31)], 
#                  split.by = "sample",
#                  group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
# wrap_plots(plots = plots, ncol = 3)
# 
# DotPlot(TNBC, features = sigGenes) + 
#   RotatedAxis()




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

y_filt4 = y_filt4[which(abs(y_filt4$NES)>1.10 & y_filt4$qvalue < 0.05),]

ggplot(y_filt4[which(abs(y_filt4$NES)>1.10 & y_filt4$qvalue < 0.05),], aes(x = Subset, y = fct_reorder(Description, abs(NES)))) + 
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


## ----Subcluster cells-------------------------------------------------------------------------------------------------------------------------------
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
  #Th10:
  "Ccr3", "Ccr6", "Spi1","Il10r",#"Il22ra1",
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

VlnPlot(Tcells, features = c("Cd4", "Cd8a","Gzmk","Nkg7", "Ccl5", "Ccr10", "Itga4"), split.by = "sample") +
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
  "Cd110",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
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
  "Cd110",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
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
  "S100a8", "S100a10", #Neut
  "Clec10a",  "Itgax", "Cd2010a", # DCs /CD11c+ SiglecH DC-sign(Cd2010a/d)
  "Itgae",#CD103
  "Itgam", "Adgre1", "Adgre4","Apoe","Cd14","Ms4a7", "Cd68", #Mac
  "Mrc1", "Msr1", "Fcgr3",  "Ccr3", "Cxcl2",  "Cx3cr1", "Il1a", "Il1b","Ccr2", #Mac
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
  "Cd110",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe"  #Plasma 
), scale = F, assay = "RNA"
) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism(base_family = "Arial", base_size = 12) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1, size = 8, face = "plain"))


colnames(Bcells[[]])

Idents(Bcells)

VlnPlot(Bcells, features = c("Cd110", "Il1a", "Pax5", "Cd80", "Cd86", "Spi1", "Spn", "Sdc1", "Cd1d1", "Itgam", "Ighd"), 
        split.by = "sample", assay = "RNA") +
  theme(legend.position = 'right')

VlnPlot(int.T1D_Timepoints, features = c("Cd110","Ighd", "Ighm", "Ighg1", "Ighg2b", "Cr2", "Fcer2a", "Fcrl1"), 
        split.by = "sample", assay = "RNA") +
  theme(legend.position = 'right')


Bcells = PrepSCTFindMarkers(Bcells)
Bcell0v1 = FindMarkers(Bcells, ident.1 = "0", ident.2 = "1",assay= "SCT",
                       recorrect_umi = FALSE )


Bcells = RenameIdents(Bcells, 

                      "5" = "2",
                      "4" = "3"
)



