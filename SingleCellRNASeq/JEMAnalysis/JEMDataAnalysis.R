###JEM Data Analysis

##Annotation
library(Seurat)
library(plyr)
library(dplyr)
library(cowplot)

library(ggplot2)
library(openxlsx)
library(stringr)
theme_set(theme_cowplot())

## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
library(Seurat)
packageVersion("Seurat")
library(SoupX)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scRNAseq")
library(scRNAseq)

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
#BiocManager::install("multtest")
library(multtest)
library(metap) 
library(ggprism)
#BiocManager::install('glmGamPoi')
library(glmGamPoi)
 

# Import Raw Files and Prprocess----

setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/RawDataset/")
Week4_expression_matrix1 <- ReadMtx(
  mtx = "GSM4213196_NOD_4w_2734_matrix.mtx.gz", features = "GSM4213196_NOD_4w_2734_genes.tsv.gz",
  cells = "GSM4213196_NOD_4w_2734_barcodes.tsv.gz")

Week4_mouse1<-CreateSeuratObject(counts = Week4_expression_matrix1,min.cells = 3, min.features = 200)


Week4_expression_matrix2 <- ReadMtx(
  mtx = "GSM4213199_NOD_4w_2849_matrix.mtx.gz", features = "GSM4213199_NOD_4w_2849_genes.tsv.gz",
  cells = "GSM4213199_NOD_4w_2849_barcodes.tsv.gz")

Week4_mouse2<-CreateSeuratObject(counts = Week4_expression_matrix2,min.cells = 3, min.features = 200)


Week8_expression_matrix1 <- ReadMtx(
  mtx = "GSM4213200_NOD_8w_2849_matrix.mtx.gz", features = "GSM4213200_NOD_8w_2849_genes.tsv.gz",
  cells = "GSM4213200_NOD_8w_2849_barcodes.tsv.gz")

Week8_mouse1<-CreateSeuratObject(counts = Week8_expression_matrix1,min.cells = 3, min.features = 200)

Week8_expression_matrix2 <- ReadMtx(
  mtx = "GSM4213197_NOD_8w_2734_matrix.mtx.gz", features = "GSM4213197_NOD_8w_2734_genes.tsv.gz",
  cells = "GSM4213197_NOD_8w_2734_barcodes.tsv.gz")

Week8_mouse2<-CreateSeuratObject(counts = Week8_expression_matrix2,min.cells = 3, min.features = 200)


Week15_expression_matrix1 <- ReadMtx(
  mtx = "GSM4213198_NOD_15w_2734_matrix.mtx.gz", features = "GSM4213198_NOD_15w_2734_genes.tsv.gz",
  cells = "GSM4213198_NOD_15w_2734_barcodes.tsv.gz")

Week15_mouse1<-CreateSeuratObject(counts = Week15_expression_matrix1,min.cells = 3, min.features = 200)

Week4_mouse1@meta.data$sample <-"Week4_1"
Week4_mouse2@meta.data$sample <-"Week4_2"
Week8_mouse1@meta.data$sample <-"Week8_1"
Week8_mouse2@meta.data$sample <-"Week8_2"
Week15_mouse1@meta.data$sample <-"Week15_1"

Week4_mouse1@meta.data$time <-"Week4"
Week4_mouse2@meta.data$time <-"Week4"
Week8_mouse1@meta.data$time <-"Week8"
Week8_mouse2@meta.data$time <-"Week8"
Week15_mouse1@meta.data$time <-"Week15"

Week4_mouse1 <- NormalizeData(Week4_mouse1) # normalize the data with log norm
Week4_mouse1 <- ScaleData(Week4_mouse1, verbose = F) # linear transformation prior to

Week4_mouse2 <- NormalizeData(Week4_mouse2) # normalize the data with log norm
Week4_mouse2 <- ScaleData(Week4_mouse2, verbose = F) # linear transformation prior to

Week8_mouse1 <- NormalizeData(Week8_mouse1) # normalize the data with log norm
Week8_mouse1 <- ScaleData(Week8_mouse1, verbose = F) # linear transformation prior to

Week8_mouse2 <- NormalizeData(Week8_mouse2) # normalize the data with log norm
Week8_mouse2 <- ScaleData(Week8_mouse2, verbose = F) # linear transformation prior to

Week15_mouse1 <- NormalizeData(Week15_mouse1) # normalize the data with log norm
Week15_mouse1 <- ScaleData(Week15_mouse1, verbose = F) # linear transformation prior to


JEM_Data <- merge(Week4_mouse1, y = c(Week4_mouse2, Week8_mouse1,Week8_mouse2,Week15_mouse1), project = "T1DProgression")

#Combined_Timepoints <- merge(Week4_mouse1, y = c( Week8_mouse1,Week15_mouse1), add.cell.ids = c("4Weeks_1", "8Weeks_1","15Weeks_1"), project = "T1DProgression")



JEM_Data[["percent.mt"]] <- PercentageFeatureSet(JEM_Data, pattern = "^mt-")

VlnPlot(JEM_Data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)


plot1 <- FeatureScatter(JEM_Data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(JEM_Data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


JEM_Data <- subset(JEM_Data, subset = nFeature_RNA > 200  & percent.mt < 5)
#' scTransform vignette:
#' <https://satijalab.org/seurat/articles/sctransform_vignette.html>
#' install glmGamPoi before using, significantly improves speed

JEM_Data <- SCTransform(JEM_Data, vars.to.regress = "percent.mt", verbose = FALSE)
getwd()
#JEM_Data  <- SCTransform(JEM_Data, vars.to.regress = "percent.mt", verbose=F)
saveRDS(JEM_Data , file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/JEM_Data.rds")

JEM_Data = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData//JEM_Data_sct_filtered.rds")
#JEM_Data <- SCTransform(JEM_Data)
JEM_Data <- RunPCA(JEM_Data)
JEM_Data <- RunUMAP(JEM_Data, dims = 1:30)
JEM_Data[[]]
DimPlot(JEM_Data, reduction = "umap", group.by = c("sample", "seurat_clusters"))
DimPlot(JEM_Data, reduction = "umap", split.by = c("sample"))

ElbowPlot(JEM_Data, ndims = 50)
JEM_Data <- FindNeighbors(JEM_Data, dims = 1:30, verbose = F)
JEM_Data <- FindClusters(JEM_Data, verbose =F, resolution = 0.4)
DimPlot(JEM_Data, group.by = c("time"), label=F)

DimPlot(JEM_Data, reduction = "umap", group.by = c("sample"))

table(Idents(JEM_Data),JEM_Data@meta.data$sample)


## ---- Annotate clustes---------------------------------------------------------------------------------------------------------------------------------------



JEM_Data <- PrepSCTFindMarkers(JEM_Data)
all.markers <- FindAllMarkers(JEM_Data , only.pos = TRUE)

top20_sct=all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)

# Plot the heatmap for the top 20 markers
DoHeatmap(
  object = JEM_Data, 
  features = top20_sct$gene,  # Use the top 20 marker genes
  size = 4,                   # Adjust text size for readability
  group.by = "ident",         # Group cells by cluster identity
  label = TRUE                # Show cluster labels
) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) # Adjust color scale








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



DotPlot(JEM_Data, features = c( "Ptprc", #CD45
                                      "Trbc1","Trbc2","Cd3g","Cd3e","Cd3d","Lat", #T Cells
                                      "C1qa","C1qb","C1qc","Lyz2", #Macrophages
                                      "Cd79a", "Cd79b","Ms4a1","Cd19",#B Cells
                                      "S100a4","S100a11","Mdh2","Gm2a","Pglyrp1","H2-Oa",#cDc
                                      "Siglech","Cd209a","Cd209d","Cox6a2","Rnase6",#pDC
                                      "Igj","Ighg2b","Ighg2c","Igha",#Plasma
                                      "Cd34","Eng","Col15a1","Tie1","Kdr","Pecam1",#Endothelial Cells
                                      "Il34","Pdgfrb","Pdgfra","Col3a1"#Mesenchymal
                                      # "Il34","Pdgfrb","Pdgfa","Col3a1","Foxs1","Mapk10" #Mesenchymal
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

DimPlot(JEM_Data,label = T,label.box = T,label.size = 8)

JEM_Data = RenameIdents(JEM_Data, 
                              "0" = "Endothelial Cells",
                              "1" = "Macrophage",
                              "2" = "T Cell",
                              "3" = "Macrophage",
                              "4" = "Endothelial Cells",
                              "5" = "Mesenchymal",
                              "6" = "T Cell",
                              "7" = "B Cell",
                              "8" = "cDc",
                              "9" = "Endothelial Cells",
                              "10" = "Macrophage",
                              "11" = "Endothelial Cells",
                              "12" = "Macrophage",
                              "13" = "Endothelial Cells",
                              "14" = "Macrophage",
                              "15" = "Plasma",
                              "16" = "pDC",
                              "17" = "Endothelial Cells",
                              "18" = "Endothelial Cells" #TBD

)


DimPlot(JEM_Data,label = T,label.box = T,label.size = 8)
saveRDS(JEM_Data , file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/JEM_Data.rds")




# Subcluster T Cells ----

Tcells_JEM = subset(x = JEM_Data, idents = "T Cell")
ElbowPlot(Tcells_JEM, ndims = 50)
Tcells_JEM <- RunUMAP(Tcells_JEM, dims = 1:30)
Tcells_JEM <- FindNeighbors(Tcells_JEM, dims = 1:30)
Tcells_JEM <- FindClusters(Tcells_JEM, res = 0.7)#graph.name = "SCT_nn"


DimPlot(Tcells_JEM, label = T,
        label.box = TRUE, label.size = 8)



DefaultAssay(Tcells_JEM)
Tcells_JEM <- PrepSCTFindMarkers(Tcells_JEM)
all.markers <- FindAllMarkers(Tcells_JEM , only.pos = TRUE,recorrect_umi = FALSE)

top20_sct=all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)


DotPlot(Tcells_JEM,assay = "RNA", features = c( "Tnfrsf4","Tbc1d4","Cd4", #CD4
                                            "Cd8a","Blk",#CD8
                                            "Ccr7","Il7r","Lef1","Sell",#Memory
                                            "Gzma","Gzmb","Klrg1","Tbx21","Zeb2",#SLEC
                                            "Ccl5", "Gzmk", "Klrc1", "Nkg7", "Cd38", "Cxcr3", "Cxcr6", "Fasl", "Ifng",#CD8 Effector
                                            "Cd40lg","Il6ra",#Tcon effector
                                            "Nt5e",#Anergy
                                            "Slamf6","Tcf7",#Progenitor
                                            "Cd200","Il21","Tnfsf8", #Th21
                                            "Rora","Il17f", "Stat3", "Ccr6",#Th17
                                            "Bcl6","Cxcr5","Tox2",#Tfh
                                            "Ctla4","Foxp3","Icos","Ikzf2","Il2ra","Tigit",#Tregs
                                            "Cd69","Egr1","Egr2","Nr4a1",#Activation
                                            "Ifit1","Ifit3","Isg15","Stat1",#Interferon sensing
                                            "Mki67","Stmn1",#Proliferation
                                            "Tox","Pdcd1","Lag3",#Exhaustion
                                            "Klra8","Klrb1a","Eomes","Ncr1","Klre1", #NK
                                            "Tcrg-C2",#Gamm_Delta-,"Trac","Cd3e"
                                            "Hs3st1","Gata3",#ILC2
                                            "Hebp1","Ncoa7","Cd74","Cited4","H2-Ab1","Rorc"#ILC3
                                            
                                            
)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

# Check average expression of Nkg7 per cluster
AverageExpression(Tcells_JEM, features = "Cd8a", group.by = "seurat_clusters")

# Define the T cell subtypes for each cluster
Tcells_JEMSubtype <- c("CD8 exhausted effector-like",  # 0
                  "Tcon exhausted effector-like",  # 1
                  "CD8 exhausted effector-like",  #2
                  "CD8 proliferating-like",  #3
                  "CD8 exhausted effector-like",  #4
                  "Tregs",  #5 
                  "Tregs",  #6 
                  "CD8 exhausted effector-like",  #7
                  "CD8 proliferating effector-like",  #8
                  "Tcon memory",  #9
                  "NK Cell",  #10
                  "CD8 exhausted effector-like",  #11
                  "Macrophage",  #12
                  "B Cell"  #13
                  
)        

# Ensure the Idents are in numeric format
ident_integers <- as.integer(Idents(Tcells_JEM))

# Check if there are any unexpected levels
unique(ident_integers)

# Add the T cell subtype information to the metadata of the Seurat object
Tcells_JEM$TCellType <- Tcells_JEMSubtype[as.integer(Idents(Tcells_JEM))]

DimPlot(Tcells_JEM,group.by = "TCellType", label = T,
        label.box = TRUE, label.size = 8)




# Change T Cell clustering name in original JEM object ----

cells_to_update <- colnames(Tcells_JEM)
JEM_Data$CellSubType<-Idents(JEM_Data)
DimPlot(JEM_Data,group.by = "CellSubType")

# Ensure Tcells$TCellType is a factor
Tcells_JEM$TCellType <- as.character(Tcells_JEM$TCellType)

# Expand levels in T1D_Timepoints$CellSubType to include new ones
JEM_Data$CellSubType <- factor(JEM_Data$CellSubType, 
                                     levels = unique(c(levels(JEM_Data$CellSubType), unique(Tcells_JEM$TCellType))))


# Assign values without generating invalid factor levels
JEM_Data$CellSubType[cells_to_update] <- Tcells_JEM$TCellType

levels(JEM_Data$CellSubType)
DimPlot(JEM_Data,group.by = "CellSubType",label = T,label.box = T,label.size = 6,repel = T)


saveRDS(JEM_Data , file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/JEM_Data.rds")






# 
# Mac_cells <- subset(Combined_Timepoints, idents = "Mac")
# 
# # new_counts <- as.sparse(t(read.csv("Macrophage_magic_imputed.csv", row.names = 1)))
# # 
# # Mac <- CreateSeuratObject(counts = (new_counts))
# # # Check the structure of the new_counts data frame to ensure it matches the format expected by Seurat
# # str(new_counts)
# # 
# # # Convert the new_counts data frame to a matrix
# # new_counts_matrix <- as.matrix(new_counts)
# # devtools::install_github("KrishnaswamyLab/MAGIC", subdir='Rmagic')
# # # Add the new_counts_matrix to your Seurat object as raw counts
# # Mac_cells@assays$Magic_RNA@data <- new_counts_matrix
# # Mac_cells[["Magic_RNA"]]<- (new_counts_matrix)
# # 
# # Mac_cells@assays$Magic
# 
# mat <- GetAssayData(object = Mac_cells, assay = "RNA", slot = "data")
# mat_matrix<-as.matrix(mat)
# # Convert the matrix to a dataframe
# mat_df <- as.data.frame(mat_matrix)
# mat_df$Genes<-rownames(mat_df)
# getwd()
# # Write the dataframe to an Excel file
# write.xlsx(mat_df, "/Users/jyotirmoyroy/Desktop/Current Projects/MetabolicSensor_Type1Diabetes/SingleCellFlux/Maccell_expression.xlsx",rowNames=TRUE)
# 
# 
# 
# ###  Macrophage SubType Analysis ###
# 
# #Perform linear dimensional reduction
# 
# Mac<-readRDS("Mac_Subtype.rds")
# DimPlot(Mac)
# Mac2<- subset(Mac, idents = "Mac1")
# mat_Mac2 <- GetAssayData(object = Mac2, assay = "RNA", layer  = "data")
# mat_matrix_Mac2<-as.matrix(mat_Mac2 )
# # Convert the matrix to a dataframe
# mat_df_Mac2 <- as.data.frame(mat_matrix_Mac2)
# mat_df_Mac2$Genes<-rownames(mat_df_Mac2)
# getwd()
# # Write the dataframe to an Excel file
# write.xlsx(mat_df_Mac2, "/Users/jyotirmoyroy/Desktop/Current Projects/MetabolicSensor_Type1Diabetes/SingleCellFlux/Mac1_expression.xlsx",rowNames=TRUE)
# 
# 
# #Clustering Analysis
# Mac_cells <- FindNeighbors(Mac_cells, dims = 1:10)
# Mac_cells <- FindClusters(Mac_cells, resolution = 1.2)
# 
# 
# #Run UMAP
# Mac_cells <- RunUMAP(Mac_cells, dims = 1:10)
# 
# # note that you can set `label = TRUE` or use the LabelClusters function to help label
# DimPlot(Mac_cells, reduction = "umap")
# 
# 
# DoHeatmap(object = Mac_cells, features= c("Apoe","Trem2",
#                                           "Birc5","Stmn1","Cdca3","Mki67","Ccna2","Cdk1","Cdca8",
#                                           "Atf3","Ccl3","Egr1","Junb",
#                                           "Cxcl9","Stat1","Ccl5","Tapbp","Tap1","Cd40","Il12b",
#                                           "Prdx1","Il1rn","Lgals1","Lgals3","Cd36","Anxa1","Anxa4"),slot = 'scale.data')
# 
# 
# 
# new.cluster.ids <- c("Mac1","Mac1","Mac1","Mac1","Mac1",
#                      "Mac3","Mac3","Mac2","Mac4","Mac2","Mac3",
#                      "Mac5","Mac1","Mac1")
# names(new.cluster.ids) <- levels(Mac_cells)
# Mac_cells <- RenameIdents(Mac_cells, new.cluster.ids)
# DimPlot(Mac_cells, reduction = "umap", label = TRUE, pt.size = 0.5) 
# 
# 
# 
# DoHeatmap(object = Mac_cells, features= c("Apoe","Trem2",
#                                           "Birc5","Stmn1","Cdca3","Mki67","Ccna2","Cdk1","Cdca8",
#                                           "Atf3","Ccl3","Egr1","Junb",
#                                           "Cxcl9","Stat1","Ccl5","Tapbp","Tap1","Cd40","Il12b",
#                                           "Prdx1","Il1rn","Lgals1","Lgals3","Cd36","Anxa1","Anxa4"),slot = 'scale.data')
# 
# 
# # Retrive Glycolysis Genes
# library("AnnotationDbi")
# library("org.Mm.eg.db")
# 
# columns(org.Mm.eg.db)
# kegg <- org.Mm.egPATH2EG
# mapped <- mappedkeys(kegg)
# kegg2 <- as.list(kegg[mapped])
# # Extract the genes related to glycolysis from the kegg2 list
# glycolysis_genes_entrez <- kegg2[['00010']]
# # Convert Entrez Gene IDs to common gene names
# glycolysis_genes_common <- select(org.Mm.eg.db, 
#                                   keys = unlist(glycolysis_genes_entrez),
#                                   columns = "SYMBOL",
#                                   keytype = "ENTREZID")
# 
# # Print the list of common gene names related to glycolysis
# library(Rmagic)
# library(magic)
# library(devtools)
# #library(phateR)
# #library(dplyr)
# library(cowplot)
# library(plotly)
# 
# 
# 
# reticulate::py_last_error()
# 
# Mac <- magic(Mac_cells)
# 
# Mac_mat <- Mac_cells@assays$Magic_RNA %>% as.matrix %>% t %>% as.data.frame
# 
# reticulate::py_last_error()
# elements_to_remove <- c("Adh7", "G6pc", "Ldhc", "Pck1", "Pdha2", "Pgk2", "Pklr", "Aldob", "Ldhb-ps", "Gapdh-ps15")
# # Remove the elements
# Gly_genes<-glycolysis_genes_common$SYMBOL
# Gly_genes <- Gly_genes[!Gly_genes %in% elements_to_remove]
# Gly_genes
# Mac<-AddModuleScore(Mac,features = list(Gly_genes) ,name="Glycolysis Activation Score")
# Mac<-AddModuleScore(Mac,features = list(c("Cxcl9","Stat1","Ccl5","Tapbp","Tap1","Cd40","Il12b")) ,name="Inflammatory_Macrophage_Score")
# 
# 
# 
# Glycolysis_Activation_Score<-Mac$Glycolysis.Activation.Score1
# Inflammatory_Macrophage_Score<-Mac$Inflammatory.Macrophage.Score1
# Data=data=as.data.frame(Glycolysis_Activation_Score)
# Data$Inflammatory_Macrophage_Score=Inflammatory_Macrophage_Score
# 
# # Calculate the correlation coefficient and R-squared value
# cor_result <- cor.test(Data$Glycolysis_Activation_Score, Data$Inflammatory_Macrophage_Score)
# cor_coef <- cor_result$estimate
# r_squared <- cor_coef^2
# 
# # Create a scatter plot with correlation and R-squared value
# ggplot(Data, aes(x = Glycolysis_Activation_Score, y = Inflammatory_Macrophage_Score)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +
#   labs(
#     title = paste("Correlation Plot (R-squared =", round(r_squared, 2), ")"),
#     x = "Glycolysis Activation Score",
#     y = "Inflammatory Macrophage Score"
#   )
# 
# 
# 
# 
# Mac<-Mac_cells
# # Apply the function to the active.ident column and store results in a new column called "time"
# Cell_name<-rownames(as.data.frame(Mac@active.ident))
# Mac@meta.data$Time <- sapply(Cell_name, extract_time)
# #_T_cells#metadata$time <- sapply(CD4_T_cells$active.ident, extract_time)
# 
# 
# 
# # Subset cells with a value of 4 in the "Time" column
# Mac_4Weeks <- subset(Mac,subset= Time == 4)
# Mac_8Weeks <- subset(Mac,subset= Time == 8)
# Mac_15Weeks <- subset(Mac,subset= Time == 15)
# 
# Mac_15Weeks <- subset(Mac_15Weeks, idents = "Mac3")
# 
# 
# Mac_15Weeks<-AddModuleScore(Mac_15Weeks,features = list(Gly_genes) ,name="Glycolysis Activation Score")
# Mac_15Weeks<-AddModuleScore(Mac_15Weeks,features = list(c("Cxcl9","Stat1","Ccl5","Tapbp","Tap1","Cd40","Il12b")) ,name="Inflammatory Macrophage Score")
# 
# 
# Glycolysis_Activation_Score<-Mac_15Weeks$Glycolysis.Activation.Score1
# Inflammatory_Macrophage_Score<-Mac_15Weeks$Inflammatory.Macrophage.Score1
# Data=data=as.data.frame(Glycolysis_Activation_Score)
# Data$Inflammatory_Macrophage_Score=Inflammatory_Macrophage_Score
# 
# # Calculate the correlation coefficient and R-squared value
# cor_result <- cor.test(Data$Glycolysis_Activation_Score, Data$Inflammatory_Macrophage_Score)
# cor_coef <- cor_result$estimate
# r_squared <- cor_coef^2
# 
# # Create a scatter plot with correlation and R-squared value
# ggplot(Data, aes(x = Glycolysis_Activation_Score, y = Inflammatory_Macrophage_Score)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +
#   labs(
#     title = paste("Correlation Plot (R-squared =", round(r_squared, 2), ")"),
#     x = "Glycolysis Activation Score",
#     y = "Inflammatory Macrophage Score"
#   )

# Cell Chat Analysis ----
# Apply the function to the active.ident column and store results in a new column called "time"
# BiocManager::install("ComplexHeatmap",force = TRUE)
# BiocManager::install("BiocNeighbors")
# 
# remove.packages("CellChat")
# 
# # You need 'remotes' package to install from GitHub
# install.packages("remotes")
# 
# # Now reinstall CellChat from source (latest version, compatible with Seurat v5)
# remotes::install_github("sqjin/CellChat")


library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
install.packages("openxlsx")
library(openxlsx)

devtools::install_github("satijalab/seurat-wrappers")
library(SeuratWrappers)
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))


JEM_Data<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/JEM_Data.rds")
DefaultAssay(JEM_Data)<-"RNA"
# Retain only the 3 essential layers
JEM_Data[["RNA"]]@layers <- JEM_Data[["RNA"]]@layers[c("counts")]
Layers(JEM_Data[["RNA"]])  # Check that it’s now clean

JEM_Data <- NormalizeData(JEM_Data) # normalize the data with log norm
JEM_Data <- ScaleData(JEM_Data) # linear transformation prior to
Layers(JEM_Data[["RNA"]])  # Check that it’s now clean

Idents(JEM_Data)<- JEM_Data$CellSubType
JEM_Data$samples<-JEM_Data$sample
JEM_Data$sample<-NULL
DimPlot(JEM_Data,label = T,label.box = T,label.size = 6,repel = T)


# Subset cells  based on the "Time" column
JEM_4Weeks <- subset(JEM_Data,subset= time == "Week4")
JEM_8Weeks <- subset(JEM_Data,subset= time == "Week8")
JEM_15Weeks <- subset(JEM_Data,subset= time == "Week15")



### Week 4----
JEM_4Weeks_CellChat  <- createCellChat(object = JEM_4Weeks, group.by = "ident",assay = "RNA")
JEM_4Weeks_CellChat<-updateCellChat(JEM_4Weeks_CellChat)
#Set Ligand Receptor Database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
#dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
JEM_4Weeks_CellChat@DB <- CellChatDB.mouse
rm(JEM_Data)
### 2. Preprocessing the expression data 
# subset the expression data of signaling genes for saving computation cost
JEM_4Weeks_CellChat <- subsetData(JEM_4Weeks_CellChat)
dim(JEM_4Weeks_CellChat@data.signaling)
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 30 * 1024^2 * 1024)  # 5 GB

JEM_4Weeks_CellChat <- identifyOverExpressedGenes(JEM_4Weeks_CellChat)
JEM_4Weeks_CellChat <- identifyOverExpressedInteractions(JEM_4Weeks_CellChat)
#The number of highly variable ligand-receptor pairs used for signaling inference is 1214  

gc()  # Free up memory before running the next computation

### 3. Compute the communication probability
ptm = Sys.time()
JEM_4Weeks_CellChat <- computeCommunProb(JEM_4Weeks_CellChat, type = "triMean")
JEM_4Weeks_CellChat <- filterCommunication(JEM_4Weeks_CellChat, min.cells = 10)

### 4. Infer communication at signaling pathway level 

JEM_4Weeks_CellChat <- computeCommunProbPathway(JEM_4Weeks_CellChat)

### 5. Calculate the aggregated cell-cell communication network 
JEM_4Weeks_CellChat <- aggregateNet(JEM_4Weeks_CellChat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
groupSize <- as.numeric(table(JEM_4Weeks_CellChat@idents))
par(mfrow = c(1,1), xpd=TRUE)
#netVisual_circle(cellChat_W6_NonProgressor@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(JEM_4Weeks_CellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength-Week 4")


mat <- JEM_4Weeks_CellChat@net$weight
par(mfrow = c(3,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

### 6. Identify TNF and related contributing signaling

# Compute the network centrality scores
JEM_4Weeks_CellChat <- netAnalysis_computeCentrality(JEM_4Weeks_CellChat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show <- c("TNF")
netAnalysis_signalingRole_network(JEM_4Weeks_CellChat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(JEM_4Weeks_CellChat)
gg1+ xlim(0,8) + ylim(0,9)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(JEM_4Weeks_CellChat, pattern = "outgoing", height = 16)
ht2 <- netAnalysis_signalingRole_heatmap(JEM_4Weeks_CellChat, pattern = "incoming", height = 16)
ht1 + ht2
dev.off()
pathways.show <- c("IFN-II") 
netVisual_aggregate(JEM_4Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

pathways.show <- c("CXCL") 
netVisual_aggregate(JEM_4Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

pathways.show <- c("CCL") 
netVisual_aggregate(JEM_4Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

pathways.show <- c("CD6") 
netVisual_aggregate(JEM_4Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

# library(NMF)
# library(ggalluvial)
# selectK(cellChat_W6_NonProgressor, pattern = "outgoing")
# nPatterns = 2
# cellChat_W6_NonProgressor <- identifyCommunicationPatterns(cellChat_W6_NonProgressor, pattern = "outgoing", k = nPatterns)

saveRDS(JEM_4Weeks_CellChat, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/JEM_4Weeks_CellChat.rds")
sessionInfo()


### Week 8----
JEM_8Weeks_CellChat  <- createCellChat(object = JEM_8Weeks, group.by = "ident",assay = "RNA")
JEM_8Weeks_CellChat<-updateCellChat(JEM_8Weeks_CellChat)
#Set Ligand Receptor Database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
#dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
JEM_8Weeks_CellChat@DB <- CellChatDB.mouse

### 2. Preprocessing the expression data 
# subset the expression data of signaling genes for saving computation cost
JEM_8Weeks_CellChat <- subsetData(JEM_8Weeks_CellChat)
dim(JEM_8Weeks_CellChat@data.signaling)
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 30 * 1024^2 * 1024)  # 5 GB

JEM_8Weeks_CellChat <- identifyOverExpressedGenes(JEM_8Weeks_CellChat)
JEM_8Weeks_CellChat <- identifyOverExpressedInteractions(JEM_8Weeks_CellChat)
#The number of highly variable ligand-receptor pairs used for signaling inference is 1214  

gc()  # Free up memory before running the next computation

### 3. Compute the communication probability
ptm = Sys.time()
JEM_8Weeks_CellChat <- computeCommunProb(JEM_8Weeks_CellChat, type = "triMean")
JEM_8Weeks_CellChat <- filterCommunication(JEM_8Weeks_CellChat, min.cells = 10)

### 4. Infer communication at signaling pathway level 

JEM_8Weeks_CellChat <- computeCommunProbPathway(JEM_8Weeks_CellChat)

### 5. Calculate the aggregated cell-cell communication network 
JEM_8Weeks_CellChat <- aggregateNet(JEM_8Weeks_CellChat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
groupSize <- as.numeric(table(JEM_8Weeks_CellChat@idents))
par(mfrow = c(1,1), xpd=TRUE)
#netVisual_circle(cellChat_W6_NonProgressor@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(JEM_8Weeks_CellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength-Week 8")


mat <- JEM_8Weeks_CellChat@net$weight
par(mfrow = c(3,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

### 6. Identify TNF and related contributing signaling

# Compute the network centrality scores
JEM_8Weeks_CellChat <- netAnalysis_computeCentrality(JEM_8Weeks_CellChat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show <- c("TNF")
netAnalysis_signalingRole_network(JEM_8Weeks_CellChat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(JEM_8Weeks_CellChat)
gg1+ xlim(0,8) + ylim(0,9)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(JEM_8Weeks_CellChat, pattern = "outgoing", height = 16)
ht2 <- netAnalysis_signalingRole_heatmap(JEM_8Weeks_CellChat, pattern = "incoming", height = 16)
ht1 + ht2
dev.off()

pathways.show <- c("IFN-II") 
netVisual_aggregate(JEM_8Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

pathways.show <- c("CXCL") 
netVisual_aggregate(JEM_8Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

pathways.show <- c("CCL") 
netVisual_aggregate(JEM_8Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

pathways.show <- c("IL16") 
netVisual_aggregate(JEM_8Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

# library(NMF)
# library(ggalluvial)
# selectK(cellChat_W6_NonProgressor, pattern = "outgoing")
# nPatterns = 2
# cellChat_W6_NonProgressor <- identifyCommunicationPatterns(cellChat_W6_NonProgressor, pattern = "outgoing", k = nPatterns)

saveRDS(JEM_8Weeks_CellChat, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/JEM_8Weeks_CellChat.rds")
sessionInfo()


### Week 15----
JEM_15Weeks_CellChat  <- createCellChat(object = JEM_15Weeks, group.by = "ident",assay = "RNA")
JEM_15Weeks_CellChat<-updateCellChat(JEM_15Weeks_CellChat)
#Set Ligand Receptor Database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
#dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
JEM_15Weeks_CellChat@DB <- CellChatDB.mouse

### 2. Preprocessing the expression data 
# subset the expression data of signaling genes for saving computation cost
JEM_15Weeks_CellChat <- subsetData(JEM_15Weeks_CellChat)
dim(JEM_15Weeks_CellChat@data.signaling)
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 30 * 1024^2 * 1024)  # 5 GB

JEM_15Weeks_CellChat <- identifyOverExpressedGenes(JEM_15Weeks_CellChat)
JEM_15Weeks_CellChat <- identifyOverExpressedInteractions(JEM_15Weeks_CellChat)
#The number of highly variable ligand-receptor pairs used for signaling inference is 1214  

gc()  # Free up memory before running the next computation

### 3. Compute the communication probability
ptm = Sys.time()
JEM_15Weeks_CellChat <- computeCommunProb(JEM_15Weeks_CellChat, type = "triMean")
JEM_15Weeks_CellChat <- filterCommunication(JEM_15Weeks_CellChat, min.cells = 10)

### 4. Infer communication at signaling pathway level 

JEM_15Weeks_CellChat <- computeCommunProbPathway(JEM_15Weeks_CellChat)

### 5. Calculate the aggregated cell-cell communication network 
JEM_15Weeks_CellChat <- aggregateNet(JEM_15Weeks_CellChat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
groupSize <- as.numeric(table(JEM_15Weeks_CellChat@idents))
par(mfrow = c(1,1), xpd=TRUE)
#netVisual_circle(cellChat_W6_NonProgressor@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(JEM_15Weeks_CellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength-Week 15")


mat <- JEM_15Weeks_CellChat@net$weight
par(mfrow = c(3,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

### 6. Identify TNF and related contributing signaling

# Compute the network centrality scores
JEM_15Weeks_CellChat <- netAnalysis_computeCentrality(JEM_15Weeks_CellChat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show <- c("TNF")
netAnalysis_signalingRole_network(JEM_15Weeks_CellChat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(JEM_15Weeks_CellChat)
gg1+ xlim(0,8) + ylim(0,9)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(JEM_15Weeks_CellChat, pattern = "outgoing", height = 16)
ht2 <- netAnalysis_signalingRole_heatmap(JEM_15Weeks_CellChat, pattern = "incoming", height = 16)
ht1 + ht2
dev.off()

pathways.show <- c("IFN-II") 
netVisual_aggregate(JEM_15Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

pathways.show <- c("CXCL") 
netVisual_aggregate(JEM_15Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

pathways.show <- c("CCL") 
netVisual_aggregate(JEM_15Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

pathways.show <- c("TNF") 
netVisual_aggregate(JEM_15Weeks_CellChat, signaling = pathways.show, layout = "circle",vertex.label.cex = 2.25)

# library(NMF)
# library(ggalluvial)
# selectK(cellChat_W6_NonProgressor, pattern = "outgoing")
# nPatterns = 2
# cellChat_W6_NonProgressor <- identifyCommunicationPatterns(cellChat_W6_NonProgressor, pattern = "outgoing", k = nPatterns)

saveRDS(JEM_15Weeks_CellChat, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JEMData/JEM_15Weeks_CellChat.rds")
sessionInfo()

# Multiple Datasets Comparison ####

## Week 4 vs Week 8 ----

object.list_W4Vs8 <- list(Week4 = JEM_4Weeks_CellChat, Week8 = JEM_8Weeks_CellChat)
cellchat_JEM_W4vs8 <- mergeCellChat(object.list_W4Vs8, add.names = names(object.list_W4Vs8))

cellchat_JEM_W4vs8

### 1. Identify altered interactions and cell populations 

gg1 <- compareInteractions(cellchat_JEM_W4vs8, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_JEM_W4vs8, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat_JEM_W4vs8, weight.scale = T, measure = "weight", vertex.label.cex = 2)

gg1 <- netVisual_heatmap(cellchat_JEM_W4vs8, measure = "weight",font.size = 12, font.size.title = 18)
#> Do heatmap based on a merged object
gg1

num.link <- sapply(object.list_W4Vs8, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list_W4Vs8)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list_W4Vs8[[i]], title = names(object.list_W4Vs8)[i], weight.MinMax = weight.MinMax)+ xlim(0,4) + ylim(0,5)
}
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_JEM_W4vs8, idents.use = "Macrophage",font.size = 18,font.size.title = 18,dot.size = 4.5,label.size = 6)
gg1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_JEM_W4vs8, idents.use = "pDC",font.size = 18,font.size.title = 18,dot.size = 4.5,label.size = 6)
gg2
gg3 <- netAnalysis_signalingChanges_scatter(cellchat_JEM_W4vs8, idents.use = "cDc",font.size = 18,font.size.title = 18,dot.size = 4.5,label.size = 6)
gg3

unique(cellchat_W6_PvsNP@meta$ident)
### 2. Identify altered signaling with distinct interaction strength 
gg1 <- rankNet(cellchat_JEM_W4vs8, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_JEM_W4vs8, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

### 3.Compare outgoing (or incoming) signaling patterns associated with each cell population 
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list_W4Vs8[[i]]@netP$pathways, object.list_W4Vs8[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list_W4Vs8[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list_W4Vs8)[i], width = 10, height = 18,font.size = 11,font.size.title = 14)
ht2 = netAnalysis_signalingRole_heatmap(object.list_W4Vs8[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list_W4Vs8)[i+1], width = 10, height = 18,font.size = 11,font.size.title = 14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list_W4Vs8[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list_W4Vs8)[i], width = 10, height = 18, color.heatmap = "GnBu",font.size = 11,font.size.title = 14)
ht2 = netAnalysis_signalingRole_heatmap(object.list_W4Vs8[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list_W4Vs8)[i+1], width = 10, height = 18, color.heatmap = "GnBu",font.size = 11,font.size.title = 14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list_W4Vs8[[i]], pattern = "all", signaling = pathway.union, title = names(object.list_W4Vs8)[i], width = 10, height = 18, color.heatmap = "OrRd",font.size = 11,font.size.title = 14)
ht2 = netAnalysis_signalingRole_heatmap(object.list_W4Vs8[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list_W4Vs8)[i+1], width = 10, height = 18, color.heatmap = "OrRd",font.size = 11,font.size.title = 14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# ### 4.Identify dysfunctional signaling by comparing the communication probabities
# netVisual_bubble(cellchat_JEM_W4vs8, sources.use = 2, targets.use = c(1,3:17),  comparison = c(1, 2), angle.x = 45)
# netVisual_bubble(cellchat_JEM_W4vs8, sources.use = 5, targets.use = c(1:4,6:17),  comparison = c(1, 2), angle.x = 45)
# 
# 
# gg1 <- netVisual_bubble(cellchat_W6_PvsNP, sources.use = 2, targets.use = c(1,3:17),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Progressor Vs Non-Progressor", angle.x = 45, remove.isolate = T)
# #> Comparing communications on a merged object
# gg2 <- netVisual_bubble(cellchat_W6_PvsNP, sources.use = 2, targets.use = c(1,3:17),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Progressor Vs Non-Progressor", angle.x = 45, remove.isolate = T)
# #> Comparing communications on a merged object
# gg1 + gg2
# 
# gg1 <- netVisual_bubble(cellchat_W6_PvsNP, sources.use = 5, targets.use = c(1:4,6:17),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Progressor Vs Non-Progressor", angle.x = 45, remove.isolate = T)
# #> Comparing communications on a merged object
# gg2 <- netVisual_bubble(cellchat_W6_PvsNP, sources.use = 5, targets.use = c(1:4,6:17),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Progressor Vs Non-Progressor", angle.x = 45, remove.isolate = T)
# #> Comparing communications on a merged object
# gg1 + gg2

# saveRDS(cellchat_W6_PvsNP, file = "cellchat_W6_PvsNP.rds")
# cellchat_W6_PvsNP<-readRDS("cellchat_W6_PvsNP.rds")


## Week 8 vs Week 15 ----

object.list_W8Vs15 <- list(Week8 = JEM_8Weeks_CellChat, Week15 = JEM_15Weeks_CellChat)
cellchat_JEM_W8Vs15 <- mergeCellChat(object.list_W8Vs15, add.names = names(object.list_W8Vs15))

cellchat_JEM_W8Vs15

### 1. Identify altered interactions and cell populations 

gg1 <- compareInteractions(cellchat_JEM_W8Vs15, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_JEM_W8Vs15, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat_JEM_W8Vs15, weight.scale = T, measure = "weight", vertex.label.cex = 2)

gg1 <- netVisual_heatmap(cellchat_JEM_W8Vs15, measure = "weight",font.size = 12, font.size.title = 18)
#> Do heatmap based on a merged object
gg1

num.link <- sapply(object.list_W8Vs15, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list_W8Vs15)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list_W8Vs15[[i]], title = names(object.list_W8Vs15)[i], weight.MinMax = weight.MinMax)+ xlim(0,4) + ylim(0,5)
}
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_JEM_W8Vs15, idents.use = "Macrophage",font.size = 18,font.size.title = 18,dot.size = 4.5,label.size = 6)
gg1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_JEM_W8Vs15, idents.use = "pDC",font.size = 18,font.size.title = 18,dot.size = 4.5,label.size = 6)
gg2
gg3 <- netAnalysis_signalingChanges_scatter(cellchat_JEM_W8Vs15, idents.use = "cDc",font.size = 18,font.size.title = 18,dot.size = 4.5,label.size = 6)
gg3

unique(cellchat_W6_PvsNP@meta$ident)
### 2. Identify altered signaling with distinct interaction strength 
gg1 <- rankNet(cellchat_JEM_W8Vs15, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_JEM_W8Vs15, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

### 3.Compare outgoing (or incoming) signaling patterns associated with each cell population 
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list_W8Vs15[[i]]@netP$pathways, object.list_W8Vs15[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list_W8Vs15[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list_W8Vs15)[i], width = 10, height = 18,font.size = 11,font.size.title = 14)
ht2 = netAnalysis_signalingRole_heatmap(object.list_W8Vs15[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list_W8Vs15)[i+1], width = 10, height = 18,font.size = 11,font.size.title = 14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list_W8Vs15[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list_W8Vs15)[i], width = 10, height = 18, color.heatmap = "GnBu",font.size = 11,font.size.title = 14)
ht2 = netAnalysis_signalingRole_heatmap(object.list_W8Vs15[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list_W8Vs15)[i+1], width = 10, height = 18, color.heatmap = "GnBu",font.size = 11,font.size.title = 14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list_W8Vs15[[i]], pattern = "all", signaling = pathway.union, title = names(object.list_W8Vs15)[i], width = 10, height = 18, color.heatmap = "OrRd",font.size = 11,font.size.title = 14)
ht2 = netAnalysis_signalingRole_heatmap(object.list_W8Vs15[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list_W8Vs15)[i+1], width = 10, height = 18, color.heatmap = "OrRd",font.size = 11,font.size.title = 14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

