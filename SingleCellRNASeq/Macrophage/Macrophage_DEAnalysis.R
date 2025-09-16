## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile")
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


## DESeq and SCPA of Macrophages ----
T1D_Timepoints<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v4.rds")

DimPlot(T1D_Timepoints)
T1D_Timepoints$time
Macrophage <- subset(T1D_Timepoints, idents = "Macrophage")
DimPlot(Macrophage,label=T,label.box = T,label.size = 4)
Macrophage_Week6 <- subset(T1D_Timepoints, idents = "Macrophage", subset = time %in% "Week6")
Macrophage_Week12 <- subset(T1D_Timepoints, idents = "Macrophage", subset = time %in% "Week12")
table(Macrophage$group)
table(Macrophage_Week6$group)
table(Macrophage_Week12$group)

Idents(Macrophage_Week6) <- Macrophage_Week6$group
#Idents(Macrophage_Week12) <- Macrophage_Week12$group

#Macrophage = PrepSCTFindMarkers(Macrophage)
Macrophage_Week6 = PrepSCTFindMarkers(Macrophage_Week6)
#Macrophage_Week12 = PrepSCTFindMarkers(Macrophage_Week12)


Macrophage_Week6_PvNP = FindMarkers(Macrophage_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                            recorrect_umi = FALSE)
Macrophage_Week6_PvNP$p_val_adj<-p.adjust(Macrophage_Week6_PvNP$p_val, method = "BH")



Macrophage_Week6_PvNP$symbol = rownames(Macrophage_Week6_PvNP)
#Macrophage_Week12_PvNP$symbol = rownames(Macrophage_Week12_PvNP)

### Volcano Plot ----


#### Week 6----
Macrophage_Week6_PvNP_data = Macrophage_Week6_PvNP

keyvals <- ifelse(
  (Macrophage_Week6_PvNP_data$avg_log2FC < -1 & Macrophage_Week6_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Macrophage_Week6_PvNP_data$avg_log2FC > 1 & Macrophage_Week6_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'

Macrophage_Week6_PvNP_data$p_val_adj
genes_to_label <- c("Jdp2", "Cebpb", "Fosl2", "Atf3","Hspa1a", "Hspa1b",
                    "Tnf")

EnhancedVolcano(Macrophage_Week6_PvNP_data,
                lab = rownames(Macrophage_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                boxedLabels= T,
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 10,
                colAlpha = 0.75,
                selectLab =genes_to_label,
                title = "Macrophage: Progressor vs Non-Progressor at Week 6",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 7,
                drawConnectors = T,
                widthConnectors = 0.5,
                xlim = c(-6, 6),  # Set x-axis from -8 to 8
                ylim = c(0,12),
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("Volcanoplot_Macrophage_T1DTimepoints_Week6.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(Macrophage_Week6_PvNP_data,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Macrophage_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)




## SCPA Analysis ----

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
msig_h <- msigdbr(species = "Mus musculus", category = "H")
msig_h <- format_pathways(msig_h)


# Load KEGG pathways from C2 category
kegg_gene_sets_msig <- msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG_LEGACY") %>%
  format_pathways()
# Load Gene Ontology Biological Process (GO:BP) pathways
GO_gene_sets <- msigdbr(species = "mouse", subcategory = "GO:BP") %>%
  format_pathways()

# Combine the pathways into a single data frame
# Combine the pathways, ensuring the format remains the same
pathways <- c(msig_h,kegg_gene_sets_msig)

### Week 6

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
  filter(Pval < 0.05,adjPval<1,FC>=1) %>%
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

table(Macrophage_Week6$group)

write.csv(scpa_out_Mac_W6,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Macrophage_T1DTimepoints_SCPA_Week6_Updated.csv",row.names = TRUE)


### Combined

Macrophage_Progressor <- seurat_extract(Macrophage,
                                           meta1 = "group", value_meta1 = "Progressor")
Macrophage_NonProgressor <- seurat_extract(Macrophage,
                                              meta1 = "group", value_meta1 = "Non-Progressor")

pathways <- c(kegg_gene_sets_msig)

scpa_out_Mac <- compare_pathways(samples = list(Macrophage_Progressor, Macrophage_NonProgressor),
                                    pathways = pathways)
head(scpa_out_Mac)




library(ggplot2)

# Filter for significant pathways
filtered_data <- scpa_out_Mac %>%
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

## Macrophage Subclustering ----



Macrophage <- subset(T1D_Timepoints, idents = "Macrophage")
ElbowPlot(Macrophage, ndims = 50)
Macrophage <- RunUMAP(Macrophage, dims = 1:30)
Macrophage <- FindNeighbors(Macrophage, dims = 1:30, verbose = F)
Macrophage <- FindClusters(Macrophage, res = 0.1)#graph.name = "SCT_nn"
table(Macrophage$SCT_snn_res.0.1)

DimPlot(Macrophage, reduction = "umap", combine = F, label = T,repel = TRUE,label.box =TRUE,
        label.size = 8,alpha=1,pt.size = 2)

# get_conserved <- function(cluster){
#   Seurat::FindConservedMarkers(Macrophage,
#                                ident.1 = cluster,
#                                grouping.var = "sample",
#                                only.pos = TRUE) %>%
#     rownames_to_column(var = "gene")  %>%
#     #left_join(y = unique(annotations[, c("gene_name", "description")]),
#     #          by = c("gene" = "gene_name")) %>%
#     cbind(cluster_id = cluster, .)
# }
# unique(Idents(Macrophage))
# 
# conserved_markers.Macrophage <- map_dfr(0:4, get_conserved)
# conserved_markers.Macrophage[is.na(conserved_markers.Macrophage)] <- 0
# head(conserved_markers.Macrophage)

DefaultAssay(Macrophage) = "SCT"
Macrophage <- PrepSCTFindMarkers(Macrophage)
all.markers <- FindAllMarkers(Macrophage , only.pos = TRUE,recorrect_umi = FALSE)
# 
top20_sct <- all.markers %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.5) %>% 
  top_n(n = 20, wt = avg_log2FC)

write.csv(top20_sct, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Annotations/top20_sct_Macrophages_T1DTimepoints.csv", row.names = FALSE)

# top20_sct_pct=all.markers %>% 
#   group_by(cluster) %>% 
#   top_n(n = 20, 
#         wt = pct.1)
# # Save the top 20 markers grouped by cluster as a CSV file
# write.csv(top20_sct_pct, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Annotations/top20_sct_pct_Macrophages_T1DTimepoints.csv", row.names = FALSE)
# 

# DotPlot(Macrophage,assay = "SCT", features = c( "Atf3","Jun","Fos","Ccl3","Egr1","Junb",#NFkB activated
#                                             "Apoe","Trem2",#Not activated
#                                             "Cxcl9","Stat1","Ccl5","Tapbp","Tap1","Tap2","Cd40","Il12b","Psmb8",# IFN-Y and NFkB activated
#                                             "Prdx1","Il1rn","Lgals1","Lgals3","Cd36","Anxa1","Anxa4","Anxa5","Anxa2","Cxcl14",#Anti-inflammatory, dead cell cleared
#                                             "Birc5","Stmn1","Cdca3","Mki7","Ccna2","Cdk1","Cdca8"#Cell-cycle
#                                             
#                                             
# )) + 
#   RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")
# 

DotPlot(Macrophage,assay = "SCT",features = c("Apoe","Selenop","Maf","Mrc1","F13a1","Csf1r",#Cluster 0-Mac1 --Anti-Inflammatory like
                                              "H2-Ab1", "H2-Aa","Cd74","Tyrobp","Itgb2","Cd52","Gm2a", #Cluster 2-Mac 2-Inflammatory, Antigen Presenting
                                                    "Clec9a","Tap1","Psmb8","Igkc","Shisa5","Ets1","Limd2","Lgals3" #Cluster 1-Mac-3- Tissue Resident


)) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")


DotPlot(Macrophage, assay = "SCT", features = c(
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
MacSubtype <- c("Mac-1", #NFKB Activatred
                "Mac-2",
                "Mac-3")        

# Ensure the Idents are in numeric format
ident_integers <- as.integer(Idents(Macrophage))

# Check if there are any unexpected levels
unique(ident_integers)

# Add the Macrophage subtype information to the metadata of the Seurat object
Macrophage$MacType <- MacSubtype[as.integer(Idents(Macrophage))]

all.markers <- FindAllMarkers(Macrophage , only.pos = TRUE,recorrect_umi = FALSE)
top20_sct <- all.markers %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.5) %>% 
  top_n(n = 20, wt = avg_log2FC)

write.csv(top20_sct, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Annotations/top20_sct_Macrophages_T1DTimepoints.csv", row.names = FALSE)

top20_sct <- all.markers %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.2) %>% 
  top_n(n = 20, wt = avg_log2FC)

write.csv(top20_sct, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Annotations/top20_sct_Macrophages_T1DTimepoints.csv", row.names = FALSE)

# Check the updated metadata
DimPlot(Macrophage, reduction = "umap",group.by = "MacType", combine = F, label = T,repel = TRUE,label.box =TRUE,
        label.size = 8,alpha=1,pt.size = 2)

Idents(Macrophage)<-Macrophage$MacType

DotPlot(Macrophage, assay = "SCT", features = c(
  "Apoe", "Selenop", "Maf", "Mrc1", "F13a1", "Csf1r", # Mac1 --Anti-Inflammatory like
  "Fgr","S100a4","Il17ra","Igkc", "Shisa5", "Ets1", "Limd2", "Lgals3",# Mac2 --Anti-Inflammatory like
  "H2-Ab1", "H2-Aa", "Cd74", "Tyrobp", "Itgb2", "Cd52", "Gm2a","Clec9a", "Tlr3", "Tap1", "Psmb8" # Mac 3 and Some Mac2 -Inflammatory, Antigen Presenting
)) + 
  RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")  # Ensure only one color scale is applied



saveRDS(Macrophage,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/Macrophage_v1.rds")

## BarPlot of Macrophage Subtype Composition ----

### Week 6----
Macrophage_w6<- subset(Macrophage, subset= time %in% "Week6")
DimPlot(Macrophage_w6, reduction = "umap", combine = F, label = T,repel = TRUE,label.box =TRUE,
        label.size = 8,alpha=1,pt.size = 2)

# Create frequency tables for Macrophage subtypes of Week 6 using active.ident and group
freq_week6_Mac <- prop.table(x = table(Macrophage_w6@active.ident,Macrophage_w6@meta.data[, "group"]), margin =2)

library(ggsci)
# Extract the colors from the plot
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=3)
#color_list <- pal_npg("nrc", alpha = 0.8)(5)
Macrophage_w6$MacType<-as.factor(Macrophage_w6$MacType)
cell_types<-levels(Macrophage_w6$MacType)
# Create frequency tables for Week 18 using active.ident and group
# Create a vector of colors for each bar
named_colors<-setNames(color_list,cell_types)
bar_colors <- named_colors[rownames(freq_week6_Mac)]
#rep(color_list, length.out = nrow(freq_week6_TCells))
# Set the resolution to 300 DPI
# Plot for Week 6
par(mar = c(5, 6, 4, 6))  # Increasing the right margin to 6

# Plot the bar plot with modified x-axis label rotation
barplot(
  freq_week6_Mac, 
  xlim = c(0, ncol(freq_week6_Mac) + 3),
  col = bar_colors, 
  legend.text = TRUE, 
  args.legend = list(
    x = ncol(freq_week6_Mac) + 2, 
    y = max(colSums(freq_week6_Mac)), 
    bty = "n", 
    cex = 1.5
  ), 
  width = 1.25, 
  ylab = "Fraction of Macrophages",
  cex.lab = 2,   # Increase font size of axis titles
  cex.names = 2,  # Increase font size of x-axis labels
  cex.axis = 1.5     # Increase font size of axis tick labels
)


### Week 12---

Macrophage_w12<- subset(Macrophage, subset= time %in% "Week12")
DimPlot(Macrophage_w12, reduction = "umap", combine = F, label = T,repel = TRUE,label.box =TRUE,
        label.size = 8,alpha=1,pt.size = 2)

# Create frequency tables for Macrophage subtypes of Week 12 using active.ident and group
freq_week12_Mac <- prop.table(x = table(Macrophage_w12@active.ident,Macrophage_w12@meta.data[, "group"]), margin =2)

library(ggsci)
# Extract the colors from the plot
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=3)
Macrophage_w12$MacType<-as.factor(Macrophage_w12$MacType)
cell_types<-levels(Macrophage_w12$MacType)
named_colors<-setNames(color_list,cell_types)
bar_colors <- named_colors[rownames(freq_week12_Mac)]

# Set the resolution to 300 DPI
# Plot for Week 12
par(mar = c(5, 6, 4, 6))  # Increasing the right margin to 6

# Plot the bar plot with modified x-axis label rotation
barplot(
  freq_week12_Mac, 
  xlim = c(0, ncol(freq_week12_Mac) + 3),
  col = bar_colors, 
  legend.text = TRUE, 
  args.legend = list(
    x = ncol(freq_week12_Mac) + 2, 
    y = max(colSums(freq_week12_Mac)), 
    bty = "n", 
    cex = 1.5
  ), 
  width = 1.25, 
  ylab = "Fraction of Macrophages",
  cex.lab = 2,   # Increase font size of axis titles
  cex.names = 2,  # Increase font size of x-axis labels
  cex.axis = 1.5     # Increase font size of axis tick labels
)




## GSEA Analysis ----
# Install packages if not installed
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) install.packages("clusterProfiler")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) install.packages("org.Mm.eg.db")
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
if (!requireNamespace("enrichplot", quietly = TRUE)) install.packages("enrichplot")

### Load libraries
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene annotation
library(msigdbr)
library(enrichplot)

### Mac1 vs Mac2----

deg_results_Mac1Vs2 <- FindMarkers(Macrophage, ident.1 = "Mac-1", ident.2 = "Mac-2",
                           logfc.threshold = 0, assay= "SCT",
                           recorrect_umi = FALSE)


deg_results_Mac1Vs2$p_val_adj<-p.adjust(deg_results_Mac1Vs2$p_val, method = "BH")

deg_results_Mac1Vs2$symbol = rownames(deg_results_Mac1Vs2)

keyvals <- ifelse(
  (deg_results_Mac1Vs2$avg_log2FC < -1 & deg_results_Mac1Vs2$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((deg_results_Mac1Vs2$avg_log2FC > 1 & deg_results_Mac1Vs2$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated: Mac-1'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated: Mac-2'



EnhancedVolcano(deg_results_Mac1Vs2,
                lab = rownames(deg_results_Mac1Vs2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                boxedLabels= T,
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(deg_results_Mac1Vs2)[which(names(keyvals) %in% c('Upregulated: Mac-1', 'Uregulated: Mac-2'))],
                title = "Mac-1 vs Mac-2",
                xlab = expression("log"[2] ~ "Fold Change (Mac-1 / Mac-2)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 15,
                drawConnectors = T,
                widthConnectors = 0.5,
               # xlim = c(-15, 15),  # Set x-axis from -8 to 8
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Volcanoplot_Mac1VsMac2_T1DTimepoints.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(deg_results_Mac1Vs2,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Macrophage_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)


# Order genes by log fold change
deg_results_Mac1Vs2 <- deg_results_Mac1Vs2[order(deg_results_Mac1Vs2$avg_log2FC, decreasing = TRUE), ]

# Prepare ranked gene list
lfc_vector_Mac1VsMac2 <- setNames(deg_results_Mac1Vs2$avg_log2FC, rownames(deg_results_Mac1Vs2))  # Named vector for GSEA

sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names

sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names

kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names

biocarta_gene_sets <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                          biocarta_gene_sets$gs_name) # Pathway names

mm_hallmark_df <- data.frame(
  gs_name = rep(names(mm_hallmark_sets), sapply(mm_hallmark_sets, length)),  # Pathway names
  gene_symbol = unlist(mm_hallmark_sets)  # Flatten the list into a single vector
)

# Combine the pathways into a single TERM2GENE dataset
combined_pathways <- bind_rows(msig_h,kegg_gene_sets_msig,GO_gene_sets)
colnames(combined_pathways) <- c("term", "gene")  # Rename columns to match clusterProfiler format

gsea_results_Mac1VsMac2 <- GSEA(
  geneList = lfc_vector_Mac1VsMac2,  # Ordered ranked gene list
  minGSSize = 5,          # Minimum gene set size
  maxGSSize = 500,        # Maximum gene set size
  pvalueCutoff = 0.05,       # p-value cutoff
  eps = 0,                # Boundary for calculating p-values
  seed = TRUE,            # Reproducibility
  pAdjustMethod = "BH",   # Benjamini-Hochberg correction
  TERM2GENE = combined_pathways  # Combined KEGG & MSigDB pathways
)
# Convert GSEA results to a data frame
gsea_results_Mac1VsMac2_df <- as.data.frame(gsea_results_Mac1VsMac2)


# Color gradient function for NES
color_gradient <- scale_fill_gradient2(
  low = "blue",      # Low NES values (negative)
  mid = "white",    # Midpoint (zero)
  high = "red",    # High NES values (positive)
  midpoint = 0,
  limits = c(min(gsea_results_Mac1VsMac2_df$NES, na.rm = TRUE), max(gsea_results_Mac1VsMac2_df$NES, na.rm = TRUE)),
  name = "NES"
)
top_gsea_Mac1VsMac2 <- gsea_results_Mac1VsMac2_df %>%
  arrange(desc(abs(NES))) %>%  # Sort by absolute NES
  head(20)  # Select top 20

ggplot(top_gsea_Mac1VsMac2, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  coord_flip() +
  color_gradient +
  labs(title = "Top 20 GSEA Results for Mac-1 vs Mac-2",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.7, size = 20, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )


# Ridge plot for top enriched pathways
ridgeplot(gsea_results_Mac1VsMac2) + ggtitle("GSEA Enrichment for Mac-1 vs Mac-2")

# Save results to CSV
write.csv(as.data.frame(gsea_results), "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/GSEA_Mac1VsMac2_results.csv")

### Mac1 vs Mac3----

deg_results_Mac1Vs3 <- FindMarkers(Macrophage, ident.1 = "Mac-1", ident.2 = "Mac-3",
                                   logfc.threshold = 0, assay= "SCT",
                                   recorrect_umi = FALSE)


deg_results_Mac1Vs3$p_val_adj<-p.adjust(deg_results_Mac1Vs3$p_val, method = "BH")

deg_results_Mac1Vs3$symbol = rownames(deg_results_Mac1Vs3)

keyvals <- ifelse(
  (deg_results_Mac1Vs3$avg_log2FC < -1 & deg_results_Mac1Vs3$p_val_adj < 0.05), '#1465AC',  
  ifelse((deg_results_Mac1Vs3$avg_log2FC > 1 & deg_results_Mac1Vs3$p_val_adj < 0.05), '#B31B21',  
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated: Mac-1'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated: Mac-3'



EnhancedVolcano(deg_results_Mac1Vs3,
                lab = rownames(deg_results_Mac1Vs3),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                boxedLabels= T,
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(deg_results_Mac1Vs3)[which(names(keyvals) %in% c('Upregulated: Mac-1', 'Uregulated: Mac-3'))],
                title = "Mac-1 vs Mac-3",
                xlab = expression("log"[2] ~ "Fold Change (Mac-1 / Mac-3)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 15,
                drawConnectors = T,
                widthConnectors = 0.5,
                # xlim = c(-15, 15),  # Set x-axis from -8 to 8
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Volcanoplot_Mac1VsMac3_T1DTimepoints.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(deg_results_Mac1Vs3,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Mac1VsMac3_DEG.csv",row.names = TRUE)


# Order genes by log fold change
deg_results_Mac1Vs3 <- deg_results_Mac1Vs3[order(deg_results_Mac1Vs3$avg_log2FC, decreasing = TRUE), ]

# Prepare ranked gene list
lfc_vector_Mac1VsMac3 <- setNames(deg_results_Mac1Vs3$avg_log2FC, rownames(deg_results_Mac1Vs3))  # Named vector for GSEA

sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names

sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names

kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names

biocarta_gene_sets <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names

mm_hallmark_df <- data.frame(
  gs_name = rep(names(mm_hallmark_sets), sapply(mm_hallmark_sets, length)),  # Pathway names
  gene_symbol = unlist(mm_hallmark_sets)  # Flatten the list into a single vector
)

# Combine the pathways into a single TERM2GENE dataset
combined_pathways <- bind_rows(msig_h,kegg_gene_sets_msig,GO_gene_sets)
colnames(combined_pathways) <- c("term", "gene")  # Rename columns to match clusterProfiler format

gsea_results_Mac1VsMac3 <- GSEA(
  geneList = lfc_vector_Mac1VsMac3,  # Ordered ranked gene list
  minGSSize = 5,          # Minimum gene set size
  maxGSSize = 500,        # Maximum gene set size
  pvalueCutoff = 0.05,       # p-value cutoff
  eps = 0,                # Boundary for calculating p-values
  seed = TRUE,            # Reproducibility
  pAdjustMethod = "BH",   # Benjamini-Hochberg correction
  TERM2GENE = combined_pathways  # Combined KEGG & MSigDB pathways
)
# Convert GSEA results to a data frame
gsea_results_Mac1VsMac3_df <- as.data.frame(gsea_results_Mac1VsMac3)


# Color gradient function for NES
color_gradient <- scale_fill_gradient2(
  low = "blue",      # Low NES values (negative)
  mid = "white",    # Midpoint (zero)
  high = "red",    # High NES values (positive)
  midpoint = 0,
  limits = c(min(gsea_results_Mac1VsMac3_df$NES, na.rm = TRUE), max(gsea_results_Mac1VsMac3_df$NES, na.rm = TRUE)),
  name = "NES"
)
top_gsea_Mac1VsMac3 <- gsea_results_Mac1VsMac3_df %>%
  arrange(desc(abs(NES))) %>%  # Sort by absolute NES
  head(20)  # Select top 20

ggplot(top_gsea_Mac1VsMac3, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  coord_flip() +
  color_gradient +
  labs(title = "Top 20 GSEA Results for Mac-1 vs Mac-3",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.7, size = 20, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )


# Ridge plot for top enriched pathways
ridgeplot(gsea_results_Mac1VsMac3) + ggtitle("GSEA Enrichment for Mac-1 vs Mac-3")

# Save results to CSV
write.csv(as.data.frame(gsea_results), "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/GSEA_Mac1VsMac3_results.csv")


### Mac2 vs Mac3----

deg_results_Mac2Vs3 <- FindMarkers(Macrophage, ident.1 = "Mac-2", ident.2 = "Mac-3",
                                   logfc.threshold = 0, assay= "SCT",
                                   recorrect_umi = FALSE)


deg_results_Mac2Vs3$p_val_adj<-p.adjust(deg_results_Mac2Vs3$p_val, method = "BH")

deg_results_Mac2Vs3$symbol = rownames(deg_results_Mac2Vs3)

keyvals <- ifelse(
  (deg_results_Mac2Vs3$avg_log2FC < -1 & deg_results_Mac2Vs3$p_val_adj < 0.05), '#1465AC',  
  ifelse((deg_results_Mac2Vs3$avg_log2FC > 1 & deg_results_Mac2Vs3$p_val_adj < 0.05), '#B31B21',  
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated: Mac-2'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated: Mac-3'



EnhancedVolcano(deg_results_Mac2Vs3,
                lab = rownames(deg_results_Mac2Vs3),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                boxedLabels= T,
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(deg_results_Mac2Vs3)[which(names(keyvals) %in% c('Upregulated: Mac-2', 'Uregulated: Mac-3'))],
                title = "Mac-2 vs Mac-3",
                xlab = expression("log"[2] ~ "Fold Change (Mac-2 / Mac-3)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 15,
                drawConnectors = T,
                widthConnectors = 0.5,
                # xlim = c(-15, 15),  # Set x-axis from -8 to 8
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Volcanoplot_Mac2VsMac3_T1DTimepoints.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(deg_results_Mac2Vs3,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Mac2VsMac3_DEG.csv",row.names = TRUE)


# Order genes by log fold change
deg_results_Mac2Vs3 <- deg_results_Mac2Vs3[order(deg_results_Mac2Vs3$avg_log2FC, decreasing = TRUE), ]

# Prepare ranked gene list
lfc_vector_Mac2VsMac3 <- setNames(deg_results_Mac2Vs3$avg_log2FC, rownames(deg_results_Mac2Vs3))  # Named vector for GSEA

sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names

sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names

kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names

biocarta_gene_sets <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names

mm_hallmark_df <- data.frame(
  gs_name = rep(names(mm_hallmark_sets), sapply(mm_hallmark_sets, length)),  # Pathway names
  gene_symbol = unlist(mm_hallmark_sets)  # Flatten the list into a single vector
)

# Combine the pathways into a single TERM2GENE dataset
combined_pathways <- bind_rows(msig_h,kegg_gene_sets_msig,GO_gene_sets)
colnames(combined_pathways) <- c("term", "gene")  # Rename columns to match clusterProfiler format

gsea_results_Mac2VsMac3 <- GSEA(
  geneList = lfc_vector_Mac2VsMac3,  # Ordered ranked gene list
  minGSSize = 5,          # Minimum gene set size
  maxGSSize = 500,        # Maximum gene set size
  pvalueCutoff = 0.05,       # p-value cutoff
  eps = 0,                # Boundary for calculating p-values
  seed = TRUE,            # Reproducibility
  pAdjustMethod = "BH",   # Benjamini-Hochberg correction
  TERM2GENE = combined_pathways  # Combined KEGG & MSigDB pathways
)
# Convert GSEA results to a data frame
gsea_results_Mac2VsMac3_df <- as.data.frame(gsea_results_Mac2VsMac3)


# Color gradient function for NES
color_gradient <- scale_fill_gradient2(
  low = "blue",      # Low NES values (negative)
  mid = "white",    # Midpoint (zero)
  high = "red",    # High NES values (positive)
  midpoint = 0,
  limits = c(min(gsea_results_Mac2VsMac3_df$NES, na.rm = TRUE), max(gsea_results_Mac2VsMac3_df$NES, na.rm = TRUE)),
  name = "NES"
)
top_gsea_Mac2VsMac3 <- gsea_results_Mac2VsMac3_df %>%
  arrange(desc(abs(NES))) %>%  # Sort by absolute NES
  head(20)  # Select top 20

ggplot(top_gsea_Mac2VsMac3, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  coord_flip() +
  color_gradient +
  labs(title = "Top 20 GSEA Results for Mac-2 vs Mac-3",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.7, size = 20, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )


# Ridge plot for top enriched pathways
ridgeplot(gsea_results_Mac2VsMac3) + ggtitle("GSEA Enrichment for Mac-2 vs Mac-3")

# Save results to CSV
write.csv(as.data.frame(gsea_results), "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/GSEA_Mac2VsMac3_results.csv")


### Progressor Vs Non-Progressor ----

### Week 6----

#### Mac1----
Macrophage<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/Macrophage_v1.rds")
Mac1 <- subset(Macrophage, idents = "Mac-1",subset = time=='Week6')
Idents(Mac1)<-Mac1$group
unique(Mac1$group)
deg_results_Mac1 <- FindMarkers(Mac1, ident.1 = "Progressor", ident.2 = "Non-Progressor",
                                   logfc.threshold = 0, assay= "SCT",
                                   recorrect_umi = FALSE)


deg_results_Mac1$p_val_adj<-p.adjust(deg_results_Mac1$p_val, method = "BH")

deg_results_Mac1$symbol = rownames(deg_results_Mac1)

keyvals <- ifelse(
  (deg_results_Mac1$avg_log2FC < -1 & deg_results_Mac1$p_val_adj < 0.05), '#1465AC',  
  ifelse((deg_results_Mac1$avg_log2FC > 1 & deg_results_Mac1$p_val_adj < 0.05), '#B31B21',  
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated: Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated: Non-Progressor'

library(EnhancedVolcano)

EnhancedVolcano(deg_results_Mac1,
                lab = rownames(deg_results_Mac1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                boxedLabels= T,
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 6,
                colAlpha = 0.75,
                selectLab = rownames(deg_results_Mac1)[which(names(keyvals) %in% c('Upregulated: Progressor', 'Uregulated: Non-Progressor'))],
                title = "Mac-1: Progressor vs Non-Progressor",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 15,
                drawConnectors = T,
                widthConnectors = 0.5,
                # xlim = c(-15, 15),  # Set x-axis from -8 to 8
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Volcanoplot_Mac1_PvsNP_Week6.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(deg_results_Mac1,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Mac1_PVsNP_DEG_Week6.csv",row.names = TRUE)

# Order genes by log fold change
deg_results_Mac1 <- deg_results_Mac1[order(deg_results_Mac1$avg_log2FC, decreasing = TRUE), ]

# Prepare ranked gene list
lfc_vector_Mac1 <- setNames(deg_results_Mac1$avg_log2FC, rownames(deg_results_Mac1))  # Named vector for GSEA

sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names

sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names

kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names

biocarta_gene_sets <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names

# mm_hallmark_df <- data.frame(
#   gs_name = rep(names(mm_hallmark_sets), sapply(mm_hallmark_sets, length)),  # Pathway names
#   gene_symbol = unlist(mm_hallmark_sets)  # Flatten the list into a single vector
# )
pwl_msigdbr<-c(pwl_hallmark,pwl_kegg,pwl_reactome)
# # Combine the pathways into a single TERM2GENE dataset
# combined_pathways <- bind_rows(msig_h,kegg_gene_sets_msig,GO_gene_sets)
# colnames(combined_pathways) <- c("term", "gene")  # Rename columns to match clusterProfiler format
combined_pathways <- data.frame(
  gs_name = rep(names(pwl_msigdbr), sapply(pwl_msigdbr, length)),  # Pathway names
  gene_symbol = unlist(pwl_msigdbr)  # Flatten the list into a single vector
)

gsea_results_Mac1 <- GSEA(
  geneList = lfc_vector_Mac1,  # Ordered ranked gene list
  minGSSize = 5,          # Minimum gene set size
  maxGSSize = 500,        # Maximum gene set size
  pvalueCutoff = 0.1,       # p-value cutoff
  eps = 0,                # Boundary for calculating p-values
  seed = TRUE,            # Reproducibility
  pAdjustMethod = "BH",   # Benjamini-Hochberg correction
  TERM2GENE = combined_pathways  # Combined KEGG & MSigDB pathways
)
# Convert GSEA results to a data frame
gsea_results_Mac1_df <- as.data.frame(gsea_results_Mac1)


# Color gradient function for NES
color_gradient <- scale_fill_gradient2(
  low = "blue",      # Low NES values (negative)
  mid = "white",    # Midpoint (zero)
  high = "red",    # High NES values (positive)
  midpoint = 0,
  limits = c(min(gsea_results_Mac1_df$NES, na.rm = TRUE), max(gsea_results_Mac1_df$NES, na.rm = TRUE)),
  name = "NES"
)
top_gsea_Mac1 <- gsea_results_Mac1_df %>%
  arrange(desc(abs(NES))) %>%  # Sort by absolute NES
  head(30)  # Select top 20

# ggplot(top_gsea_Mac1, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
#   geom_bar(stat = "identity", show.legend = TRUE) +
#   coord_flip() +
#   color_gradient +
#   labs(title = "Top GSEA Results for Mac-1: Progressor vs Non-Progressor",
#        x = "Pathway",
#        y = "Normalized Enrichment Score (NES)") +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.7, size = 20, face = "bold"),
#     axis.title = element_text(size = 19),
#     axis.text = element_text(size = 12),
#     panel.grid.major = element_line(color = "grey90"),
#     panel.grid.minor = element_blank(),
#     legend.position = "bottom"
#   )

ggplot(top_gsea_Mac1, aes(x = reorder(Description, NES), y = NES, color = NES, size = p.adjust)) +
  geom_point(show.legend = TRUE) +  # Plot points with varying sizes
  scale_size_continuous(range = c(3, 10), name = "Adjusted p-value") +  # Control the range of bubble sizes
  scale_color_gradient2(
    low = "blue",      # Low NES values (negative)
    mid = "white",     # Midpoint (zero)
    high = "red",      # High NES values (positive)
    midpoint = 0,
    limits = c(min(gsea_results_Mac1_df$NES, na.rm = TRUE), max(gsea_results_Mac1_df$NES, na.rm = TRUE)),
    name = "NES"
  ) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1.5, color = "black") +  # Add bold line at y = 0 (NES=0)
  coord_flip() +
  labs(
    title = "GSEA Mac-1: Progressor vs Non-Progressor",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold", family = "Arial"),  # Title with larger size and better font
    axis.title = element_text(size = 20, face = "bold", family = "Arial"),  # Axis titles bold and larger font
    axis.text = element_text(size = 16, family = "Arial"),  # Axis text size increased and font changed
    panel.grid.major.x = element_blank(),  # Removed x-axis gridlines
    panel.grid.major.y = element_line(color = "grey90"),  # Horizontal gridlines retained
    panel.grid.minor = element_blank(),  # Minor gridlines removed
    legend.position = "right",  # Move legend to the right side
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # Bold legend title with improved font
    legend.text = element_text(size = 14, family = "Arial")  # Increased legend text size
  )


# Save results to CSV
write.csv(gsea_results_Mac1_df, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/GSEA_Mac1_ProgressorVsNonProgressor_Week6.csv")

#### Mac2----
Macrophage<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/Macrophage_v1.rds")
Mac2 <- subset(Macrophage, idents = "Mac-2",subset = time=='Week6')
Idents(Mac2)<-Mac2$group
unique(Mac2$group)
deg_results_Mac2 <- FindMarkers(Mac2, ident.1 = "Progressor", ident.2 = "Non-Progressor",
                                logfc.threshold = 0, assay= "SCT",
                                recorrect_umi = FALSE)


deg_results_Mac2$p_val_adj<-p.adjust(deg_results_Mac2$p_val, method = "BH")

deg_results_Mac2$symbol = rownames(deg_results_Mac2)

keyvals <- ifelse(
  (deg_results_Mac2$avg_log2FC < -1 & deg_results_Mac2$p_val_adj < 0.05), '#1465AC',  
  ifelse((deg_results_Mac2$avg_log2FC > 1 & deg_results_Mac2$p_val_adj < 0.05), '#B31B21',  
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated: Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated: Non-Progressor'

library(EnhancedVolcano)

EnhancedVolcano(deg_results_Mac2,
                lab = rownames(deg_results_Mac2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                boxedLabels= T,
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 6,
                colAlpha = 0.75,
                selectLab = rownames(deg_results_Mac2)[which(names(keyvals) %in% c('Upregulated: Progressor', 'Uregulated: Non-Progressor'))],
                title = "Mac-2: Progressor vs Non-Progressor",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 15,
                drawConnectors = T,
                widthConnectors = 0.5,
                # xlim = c(-15, 15),  # Set x-axis from -8 to 8
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Volcanoplot_Mac2_PvsNP_Week6.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(deg_results_Mac2,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Mac2_PVsNP_DEG_Week6.csv",row.names = TRUE)

# Order genes by log fold change
deg_results_Mac2 <- deg_results_Mac2[order(deg_results_Mac2$avg_log2FC, decreasing = TRUE), ]

# Prepare ranked gene list
lfc_vector_Mac2 <- setNames(deg_results_Mac2$avg_log2FC, rownames(deg_results_Mac2))  # Named vector for GSEA

sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names

sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names

kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names

biocarta_gene_sets <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names

# mm_hallmark_df <- data.frame(
#   gs_name = rep(names(mm_hallmark_sets), sapply(mm_hallmark_sets, length)),  # Pathway names
#   gene_symbol = unlist(mm_hallmark_sets)  # Flatten the list into a single vector
# )
pwl_msigdbr<-c(pwl_hallmark,pwl_kegg,pwl_reactome)
# # Combine the pathways into a single TERM2GENE dataset
# combined_pathways <- bind_rows(msig_h,kegg_gene_sets_msig,GO_gene_sets)
# colnames(combined_pathways) <- c("term", "gene")  # Rename columns to match clusterProfiler format
combined_pathways <- data.frame(
  gs_name = rep(names(pwl_msigdbr), sapply(pwl_msigdbr, length)),  # Pathway names
  gene_symbol = unlist(pwl_msigdbr)  # Flatten the list into a single vector
)

gsea_results_Mac2 <- GSEA(
  geneList = lfc_vector_Mac2,  # Ordered ranked gene list
  minGSSize = 5,          # Minimum gene set size
  maxGSSize = 500,        # Maximum gene set size
  pvalueCutoff = 0.1,       # p-value cutoff
  eps = 0,                # Boundary for calculating p-values
  seed = TRUE,            # Reproducibility
  pAdjustMethod = "BH",   # Benjamini-Hochberg correction
  TERM2GENE = combined_pathways  # Combined KEGG & MSigDB pathways
)
# Convert GSEA results to a data frame
gsea_results_Mac2_df <- as.data.frame(gsea_results_Mac2)


# Color gradient function for NES
color_gradient <- scale_fill_gradient2(
  low = "blue",      # Low NES values (negative)
  mid = "white",    # Midpoint (zero)
  high = "red",    # High NES values (positive)
  midpoint = 0,
  limits = c(min(gsea_results_Mac2_df$NES, na.rm = TRUE), max(gsea_results_Mac2_df$NES, na.rm = TRUE)),
  name = "NES"
)
top_gsea_Mac2 <- gsea_results_Mac2_df %>%
  arrange(desc(abs(NES))) %>%  # Sort by absolute NES
  head(30)  # Select top 20

# ggplot(top_gsea_Mac2, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
#   geom_bar(stat = "identity", show.legend = TRUE) +
#   coord_flip() +
#   color_gradient +
#   labs(title = "Top GSEA Results for Mac-1: Progressor vs Non-Progressor",
#        x = "Pathway",
#        y = "Normalized Enrichment Score (NES)") +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.7, size = 20, face = "bold"),
#     axis.title = element_text(size = 19),
#     axis.text = element_text(size = 12),
#     panel.grid.major = element_line(color = "grey90"),
#     panel.grid.minor = element_blank(),
#     legend.position = "bottom"
#   )

ggplot(top_gsea_Mac2, aes(x = reorder(Description, NES), y = NES, color = NES, size = p.adjust)) +
  geom_point(show.legend = TRUE) +  # Plot points with varying sizes
  scale_size_continuous(range = c(3, 10), name = "Adjusted p-value") +  # Control the range of bubble sizes
  scale_color_gradient2(
    low = "blue",      # Low NES values (negative)
    mid = "white",     # Midpoint (zero)
    high = "red",      # High NES values (positive)
    midpoint = 0,
    limits = c(min(gsea_results_Mac2_df$NES, na.rm = TRUE), max(gsea_results_Mac2_df$NES, na.rm = TRUE)),
    name = "NES"
  ) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1.5, color = "black") +  # Add bold line at y = 0 (NES=0)
  coord_flip() +
  labs(
    title = "GSEA Mac-2: Progressor vs Non-Progressor",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold", family = "Arial"),  # Title with larger size and better font
    axis.title = element_text(size = 20, face = "bold", family = "Arial"),  # Axis titles bold and larger font
    axis.text = element_text(size = 16, family = "Arial"),  # Axis text size increased and font changed
    panel.grid.major.x = element_blank(),  # Removed x-axis gridlines
    panel.grid.major.y = element_line(color = "grey90"),  # Horizontal gridlines retained
    panel.grid.minor = element_blank(),  # Minor gridlines removed
    legend.position = "right",  # Move legend to the right side
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # Bold legend title with improved font
    legend.text = element_text(size = 14, family = "Arial")  # Increased legend text size
  )


# Save results to CSV
write.csv(gsea_results_Mac2_df, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/GSEA_Mac2_ProgressorVsNonProgressor_Week6.csv")


#### Mac3----
Macrophage<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/Macrophage_v1.rds")
Mac3 <- subset(Macrophage, idents = "Mac-3",subset = time=='Week6')
Idents(Mac3)<-Mac3$group
unique(Mac3$group)
deg_results_Mac3 <- FindMarkers(Mac3, ident.1 = "Progressor", ident.2 = "Non-Progressor",
                                logfc.threshold = 0, assay= "SCT",
                                recorrect_umi = FALSE)


deg_results_Mac3$p_val_adj<-p.adjust(deg_results_Mac3$p_val, method = "BH")

deg_results_Mac3$symbol = rownames(deg_results_Mac3)

keyvals <- ifelse(
  (deg_results_Mac3$avg_log2FC < -1 & deg_results_Mac3$p_val_adj < 0.05), '#1465AC',  
  ifelse((deg_results_Mac3$avg_log2FC > 1 & deg_results_Mac3$p_val_adj < 0.05), '#B31B21',  
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated: Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated: Non-Progressor'

library(EnhancedVolcano)

EnhancedVolcano(deg_results_Mac3,
                lab = rownames(deg_results_Mac3),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                boxedLabels= T,
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 6,
                colAlpha = 0.75,
                selectLab = rownames(deg_results_Mac3)[which(names(keyvals) %in% c('Upregulated: Progressor', 'Uregulated: Non-Progressor'))],
                title = "Mac-3: Progressor vs Non-Progressor",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 15,
                drawConnectors = T,
                widthConnectors = 0.5,
                # xlim = c(-15, 15),  # Set x-axis from -8 to 8
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Volcanoplot_Mac3_PvsNP_Week6.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(deg_results_Mac3,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Mac3_PVsNP_DEG_Week6.csv",row.names = TRUE)

# Order genes by log fold change
deg_results_Mac3 <- deg_results_Mac3[order(deg_results_Mac3$avg_log2FC, decreasing = TRUE), ]

# Prepare ranked gene list
lfc_vector_Mac3 <- setNames(deg_results_Mac3$avg_log2FC, rownames(deg_results_Mac3))  # Named vector for GSEA

sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names

sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names

kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names

biocarta_gene_sets <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names

# mm_hallmark_df <- data.frame(
#   gs_name = rep(names(mm_hallmark_sets), sapply(mm_hallmark_sets, length)),  # Pathway names
#   gene_symbol = unlist(mm_hallmark_sets)  # Flatten the list into a single vector
# )
pwl_msigdbr<-c(pwl_hallmark,pwl_kegg,pwl_reactome)
# # Combine the pathways into a single TERM2GENE dataset
# combined_pathways <- bind_rows(msig_h,kegg_gene_sets_msig,GO_gene_sets)
# colnames(combined_pathways) <- c("term", "gene")  # Rename columns to match clusterProfiler format
combined_pathways <- data.frame(
  gs_name = rep(names(pwl_msigdbr), sapply(pwl_msigdbr, length)),  # Pathway names
  gene_symbol = unlist(pwl_msigdbr)  # Flatten the list into a single vector
)

gsea_results_Mac3 <- GSEA(
  geneList = lfc_vector_Mac3,  # Ordered ranked gene list
  minGSSize = 5,          # Minimum gene set size
  maxGSSize = 500,        # Maximum gene set size
  pvalueCutoff = 0.1,       # p-value cutoff
  eps = 0,                # Boundary for calculating p-values
  seed = TRUE,            # Reproducibility
  pAdjustMethod = "BH",   # Benjamini-Hochberg correction
  TERM2GENE = combined_pathways  # Combined KEGG & MSigDB pathways
)
# Convert GSEA results to a data frame
gsea_results_Mac3_df <- as.data.frame(gsea_results_Mac3)


# Color gradient function for NES
color_gradient <- scale_fill_gradient2(
  low = "blue",      # Low NES values (negative)
  mid = "white",    # Midpoint (zero)
  high = "red",    # High NES values (positive)
  midpoint = 0,
  limits = c(min(gsea_results_Mac3_df$NES, na.rm = TRUE), max(gsea_results_Mac3_df$NES, na.rm = TRUE)),
  name = "NES"
)
top_gsea_Mac3 <- gsea_results_Mac3_df %>%
  arrange(desc(abs(NES))) %>%  # Sort by absolute NES
  head(30)  # Select top 20

# ggplot(top_gsea_Mac3, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
#   geom_bar(stat = "identity", show.legend = TRUE) +
#   coord_flip() +
#   color_gradient +
#   labs(title = "Top GSEA Results for Mac-1: Progressor vs Non-Progressor",
#        x = "Pathway",
#        y = "Normalized Enrichment Score (NES)") +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.7, size = 20, face = "bold"),
#     axis.title = element_text(size = 19),
#     axis.text = element_text(size = 12),
#     panel.grid.major = element_line(color = "grey90"),
#     panel.grid.minor = element_blank(),
#     legend.position = "bottom"
#   )

ggplot(top_gsea_Mac3, aes(x = reorder(Description, NES), y = NES, color = NES, size = p.adjust)) +
  geom_point(show.legend = TRUE) +  # Plot points with varying sizes
  scale_size_continuous(range = c(3, 10), name = "Adjusted p-value") +  # Control the range of bubble sizes
  scale_color_gradient2(
    low = "blue",      # Low NES values (negative)
    mid = "white",     # Midpoint (zero)
    high = "red",      # High NES values (positive)
    midpoint = 0,
    limits = c(min(gsea_results_Mac3_df$NES, na.rm = TRUE), max(gsea_results_Mac3_df$NES, na.rm = TRUE)),
    name = "NES"
  ) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1.5, color = "black") +  # Add bold line at y = 0 (NES=0)
  coord_flip() +
  labs(
    title = "GSEA Mac-2: Progressor vs Non-Progressor",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold", family = "Arial"),  # Title with larger size and better font
    axis.title = element_text(size = 20, face = "bold", family = "Arial"),  # Axis titles bold and larger font
    axis.text = element_text(size = 16, family = "Arial"),  # Axis text size increased and font changed
    panel.grid.major.x = element_blank(),  # Removed x-axis gridlines
    panel.grid.major.y = element_line(color = "grey90"),  # Horizontal gridlines retained
    panel.grid.minor = element_blank(),  # Minor gridlines removed
    legend.position = "right",  # Move legend to the right side
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # Bold legend title with improved font
    legend.text = element_text(size = 14, family = "Arial")  # Increased legend text size
  )


# Save results to CSV
write.csv(gsea_results_Mac3_df, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/GSEA_Mac3_ProgressorVsNonProgressor_Week6.csv")


### Week 12----

#### Mac1----
Macrophage<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/Macrophage_v1.rds")
Mac1 <- subset(Macrophage, idents = "Mac-1",subset = time=='Week12')
Idents(Mac1)<-Mac1$group
unique(Mac1$group)
deg_results_Mac1 <- FindMarkers(Mac1, ident.1 = "Progressor", ident.2 = "Non-Progressor",
                                logfc.threshold = 0, assay= "SCT",
                                recorrect_umi = FALSE)


deg_results_Mac1$p_val_adj<-p.adjust(deg_results_Mac1$p_val, method = "BH")

deg_results_Mac1$symbol = rownames(deg_results_Mac1)

keyvals <- ifelse(
  (deg_results_Mac1$avg_log2FC < -1 & deg_results_Mac1$p_val_adj < 0.05), '#1465AC',  
  ifelse((deg_results_Mac1$avg_log2FC > 1 & deg_results_Mac1$p_val_adj < 0.05), '#B31B21',  
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated: Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated: Non-Progressor'

library(EnhancedVolcano)

EnhancedVolcano(deg_results_Mac1,
                lab = rownames(deg_results_Mac1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                boxedLabels= T,
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 6,
                colAlpha = 0.75,
                selectLab = rownames(deg_results_Mac1)[which(names(keyvals) %in% c('Upregulated: Progressor', 'Uregulated: Non-Progressor'))],
                title = "Mac-1: Progressor vs Non-Progressor",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 15,
                drawConnectors = T,
                widthConnectors = 0.5,
                # xlim = c(-15, 15),  # Set x-axis from -8 to 8
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Volcanoplot_Mac1_PvsNP_Week12.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(deg_results_Mac1,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Mac1_PVsNP_DEG_Week12.csv",row.names = TRUE)

# Order genes by log fold change
deg_results_Mac1 <- deg_results_Mac1[order(deg_results_Mac1$avg_log2FC, decreasing = TRUE), ]

# Prepare ranked gene list
lfc_vector_Mac1 <- setNames(deg_results_Mac1$avg_log2FC, rownames(deg_results_Mac1))  # Named vector for GSEA

sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names

sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names

kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names

biocarta_gene_sets <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names

# mm_hallmark_df <- data.frame(
#   gs_name = rep(names(mm_hallmark_sets), sapply(mm_hallmark_sets, length)),  # Pathway names
#   gene_symbol = unlist(mm_hallmark_sets)  # Flatten the list into a single vector
# )
pwl_msigdbr<-c(pwl_hallmark,pwl_kegg,pwl_reactome)
# # Combine the pathways into a single TERM2GENE dataset
# combined_pathways <- bind_rows(msig_h,kegg_gene_sets_msig,GO_gene_sets)
# colnames(combined_pathways) <- c("term", "gene")  # Rename columns to match clusterProfiler format
combined_pathways <- data.frame(
  gs_name = rep(names(pwl_msigdbr), sapply(pwl_msigdbr, length)),  # Pathway names
  gene_symbol = unlist(pwl_msigdbr)  # Flatten the list into a single vector
)

gsea_results_Mac1 <- GSEA(
  geneList = lfc_vector_Mac1,  # Ordered ranked gene list
  minGSSize = 5,          # Minimum gene set size
  maxGSSize = 500,        # Maximum gene set size
  pvalueCutoff = 0.1,       # p-value cutoff
  eps = 0,                # Boundary for calculating p-values
  seed = TRUE,            # Reproducibility
  pAdjustMethod = "BH",   # Benjamini-Hochberg correction
  TERM2GENE = combined_pathways  # Combined KEGG & MSigDB pathways
)
# Convert GSEA results to a data frame
gsea_results_Mac1_df <- as.data.frame(gsea_results_Mac1)


# Color gradient function for NES
color_gradient <- scale_fill_gradient2(
  low = "blue",      # Low NES values (negative)
  mid = "white",    # Midpoint (zero)
  high = "red",    # High NES values (positive)
  midpoint = 0,
  limits = c(min(gsea_results_Mac1_df$NES, na.rm = TRUE), max(gsea_results_Mac1_df$NES, na.rm = TRUE)),
  name = "NES"
)
top_gsea_Mac1 <- gsea_results_Mac1_df %>%
  arrange(desc(abs(NES))) %>%  # Sort by absolute NES
  head(30)  # Select top 20

# ggplot(top_gsea_Mac1, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
#   geom_bar(stat = "identity", show.legend = TRUE) +
#   coord_flip() +
#   color_gradient +
#   labs(title = "Top GSEA Results for Mac-1: Progressor vs Non-Progressor",
#        x = "Pathway",
#        y = "Normalized Enrichment Score (NES)") +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.7, size = 20, face = "bold"),
#     axis.title = element_text(size = 19),
#     axis.text = element_text(size = 12),
#     panel.grid.major = element_line(color = "grey90"),
#     panel.grid.minor = element_blank(),
#     legend.position = "bottom"
#   )

ggplot(top_gsea_Mac1, aes(x = reorder(Description, NES), y = NES, color = NES, size = p.adjust)) +
  geom_point(show.legend = TRUE) +  # Plot points with varying sizes
  scale_size_continuous(range = c(3, 10), name = "Adjusted p-value") +  # Control the range of bubble sizes
  scale_color_gradient2(
    low = "blue",      # Low NES values (negative)
    mid = "white",     # Midpoint (zero)
    high = "red",      # High NES values (positive)
    midpoint = 0,
    limits = c(min(gsea_results_Mac1_df$NES, na.rm = TRUE), max(gsea_results_Mac1_df$NES, na.rm = TRUE)),
    name = "NES"
  ) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1.5, color = "black") +  # Add bold line at y = 0 (NES=0)
  coord_flip() +
  labs(
    title = "GSEA Mac-1: Progressor vs Non-Progressor",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold", family = "Arial"),  # Title with larger size and better font
    axis.title = element_text(size = 20, face = "bold", family = "Arial"),  # Axis titles bold and larger font
    axis.text = element_text(size = 16, family = "Arial"),  # Axis text size increased and font changed
    panel.grid.major.x = element_blank(),  # Removed x-axis gridlines
    panel.grid.major.y = element_line(color = "grey90"),  # Horizontal gridlines retained
    panel.grid.minor = element_blank(),  # Minor gridlines removed
    legend.position = "right",  # Move legend to the right side
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # Bold legend title with improved font
    legend.text = element_text(size = 14, family = "Arial")  # Increased legend text size
  )


# Save results to CSV
write.csv(gsea_results_Mac1_df, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/GSEA_Mac1_ProgressorVsNonProgressor_Week12.csv")

#### Mac2----
Macrophage<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/Macrophage_v1.rds")
Mac2 <- subset(Macrophage, idents = "Mac-2",subset = time=='Week12')
Idents(Mac2)<-Mac2$group
unique(Mac2$group)
deg_results_Mac2 <- FindMarkers(Mac2, ident.1 = "Progressor", ident.2 = "Non-Progressor",
                                logfc.threshold = 0, assay= "SCT",
                                recorrect_umi = FALSE)


deg_results_Mac2$p_val_adj<-p.adjust(deg_results_Mac2$p_val, method = "BH")

deg_results_Mac2$symbol = rownames(deg_results_Mac2)

keyvals <- ifelse(
  (deg_results_Mac2$avg_log2FC < -1 & deg_results_Mac2$p_val_adj < 0.05), '#1465AC',  
  ifelse((deg_results_Mac2$avg_log2FC > 1 & deg_results_Mac2$p_val_adj < 0.05), '#B31B21',  
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated: Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated: Non-Progressor'

library(EnhancedVolcano)

EnhancedVolcano(deg_results_Mac2,
                lab = rownames(deg_results_Mac2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                boxedLabels= T,
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 6,
                colAlpha = 0.75,
                selectLab = rownames(deg_results_Mac2)[which(names(keyvals) %in% c('Upregulated: Progressor', 'Uregulated: Non-Progressor'))],
                title = "Mac-2: Progressor vs Non-Progressor",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 15,
                drawConnectors = T,
                widthConnectors = 0.5,
                # xlim = c(-15, 15),  # Set x-axis from -8 to 8
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Volcanoplot_Mac2_PvsNP_Week12.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(deg_results_Mac2,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Mac2_PVsNP_DEG_Week12.csv",row.names = TRUE)

# Order genes by log fold change
deg_results_Mac2 <- deg_results_Mac2[order(deg_results_Mac2$avg_log2FC, decreasing = TRUE), ]

# Prepare ranked gene list
lfc_vector_Mac2 <- setNames(deg_results_Mac2$avg_log2FC, rownames(deg_results_Mac2))  # Named vector for GSEA

sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names

sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names

kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names

biocarta_gene_sets <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names

# mm_hallmark_df <- data.frame(
#   gs_name = rep(names(mm_hallmark_sets), sapply(mm_hallmark_sets, length)),  # Pathway names
#   gene_symbol = unlist(mm_hallmark_sets)  # Flatten the list into a single vector
# )
pwl_msigdbr<-c(pwl_hallmark,pwl_kegg,pwl_reactome)
# # Combine the pathways into a single TERM2GENE dataset
# combined_pathways <- bind_rows(msig_h,kegg_gene_sets_msig,GO_gene_sets)
# colnames(combined_pathways) <- c("term", "gene")  # Rename columns to match clusterProfiler format
combined_pathways <- data.frame(
  gs_name = rep(names(pwl_msigdbr), sapply(pwl_msigdbr, length)),  # Pathway names
  gene_symbol = unlist(pwl_msigdbr)  # Flatten the list into a single vector
)

gsea_results_Mac2 <- GSEA(
  geneList = lfc_vector_Mac2,  # Ordered ranked gene list
  minGSSize = 5,          # Minimum gene set size
  maxGSSize = 500,        # Maximum gene set size
  pvalueCutoff = 0.1,       # p-value cutoff
  eps = 0,                # Boundary for calculating p-values
  seed = TRUE,            # Reproducibility
  pAdjustMethod = "BH",   # Benjamini-Hochberg correction
  TERM2GENE = combined_pathways  # Combined KEGG & MSigDB pathways
)
# Convert GSEA results to a data frame
gsea_results_Mac2_df <- as.data.frame(gsea_results_Mac2)


# Color gradient function for NES
color_gradient <- scale_fill_gradient2(
  low = "blue",      # Low NES values (negative)
  mid = "white",    # Midpoint (zero)
  high = "red",    # High NES values (positive)
  midpoint = 0,
  limits = c(min(gsea_results_Mac2_df$NES, na.rm = TRUE), max(gsea_results_Mac2_df$NES, na.rm = TRUE)),
  name = "NES"
)
top_gsea_Mac2 <- gsea_results_Mac2_df %>%
  arrange(desc(abs(NES))) %>%  # Sort by absolute NES
  head(30)  # Select top 20

# ggplot(top_gsea_Mac2, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
#   geom_bar(stat = "identity", show.legend = TRUE) +
#   coord_flip() +
#   color_gradient +
#   labs(title = "Top GSEA Results for Mac-1: Progressor vs Non-Progressor",
#        x = "Pathway",
#        y = "Normalized Enrichment Score (NES)") +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.7, size = 20, face = "bold"),
#     axis.title = element_text(size = 19),
#     axis.text = element_text(size = 12),
#     panel.grid.major = element_line(color = "grey90"),
#     panel.grid.minor = element_blank(),
#     legend.position = "bottom"
#   )

ggplot(top_gsea_Mac2, aes(x = reorder(Description, NES), y = NES, color = NES, size = p.adjust)) +
  geom_point(show.legend = TRUE) +  # Plot points with varying sizes
  scale_size_continuous(range = c(3, 10), name = "Adjusted p-value") +  # Control the range of bubble sizes
  scale_color_gradient2(
    low = "blue",      # Low NES values (negative)
    mid = "white",     # Midpoint (zero)
    high = "red",      # High NES values (positive)
    midpoint = 0,
    limits = c(min(gsea_results_Mac2_df$NES, na.rm = TRUE), max(gsea_results_Mac2_df$NES, na.rm = TRUE)),
    name = "NES"
  ) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1.5, color = "black") +  # Add bold line at y = 0 (NES=0)
  coord_flip() +
  labs(
    title = "GSEA Mac-2: Progressor vs Non-Progressor",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold", family = "Arial"),  # Title with larger size and better font
    axis.title = element_text(size = 20, face = "bold", family = "Arial"),  # Axis titles bold and larger font
    axis.text = element_text(size = 16, family = "Arial"),  # Axis text size increased and font changed
    panel.grid.major.x = element_blank(),  # Removed x-axis gridlines
    panel.grid.major.y = element_line(color = "grey90"),  # Horizontal gridlines retained
    panel.grid.minor = element_blank(),  # Minor gridlines removed
    legend.position = "right",  # Move legend to the right side
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # Bold legend title with improved font
    legend.text = element_text(size = 14, family = "Arial")  # Increased legend text size
  )


# Save results to CSV
write.csv(gsea_results_Mac2_df, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/GSEA_Mac2_ProgressorVsNonProgressor_Week12.csv")


#### Mac3----
Macrophage<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/Macrophage_v1.rds")
Mac3 <- subset(Macrophage, idents = "Mac-3",subset = time=='Week12')
Idents(Mac3)<-Mac3$group
unique(Mac3$group)
deg_results_Mac3 <- FindMarkers(Mac3, ident.1 = "Progressor", ident.2 = "Non-Progressor",
                                logfc.threshold = 0, assay= "SCT",
                                recorrect_umi = FALSE)


deg_results_Mac3$p_val_adj<-p.adjust(deg_results_Mac3$p_val, method = "BH")

deg_results_Mac3$symbol = rownames(deg_results_Mac3)

keyvals <- ifelse(
  (deg_results_Mac3$avg_log2FC < -1 & deg_results_Mac3$p_val_adj < 0.05), '#1465AC',  
  ifelse((deg_results_Mac3$avg_log2FC > 1 & deg_results_Mac3$p_val_adj < 0.05), '#B31B21',  
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated: Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated: Non-Progressor'

library(EnhancedVolcano)

EnhancedVolcano(deg_results_Mac3,
                lab = rownames(deg_results_Mac3),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                boxedLabels= T,
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 6,
                colAlpha = 0.75,
                selectLab = rownames(deg_results_Mac3)[which(names(keyvals) %in% c('Upregulated: Progressor', 'Uregulated: Non-Progressor'))],
                title = "Mac-3: Progressor vs Non-Progressor",
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p"[adj] * ")"),  # Subscript in y-axis
                max.overlaps = 15,
                drawConnectors = T,
                widthConnectors = 0.5,
                # xlim = c(-15, 15),  # Set x-axis from -8 to 8
                legendPosition = 'right',
                caption = NULL ) # Removes watermark text)



# Save the plot
ggsave("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Volcanoplot_Mac3_PvsNP_Week12.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(deg_results_Mac3,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/Mac3_PVsNP_DEG_Week12.csv",row.names = TRUE)

# Order genes by log fold change
deg_results_Mac3 <- deg_results_Mac3[order(deg_results_Mac3$avg_log2FC, decreasing = TRUE), ]

# Prepare ranked gene list
lfc_vector_Mac3 <- setNames(deg_results_Mac3$avg_log2FC, rownames(deg_results_Mac3))  # Named vector for GSEA

sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names

sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names

kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names

biocarta_gene_sets <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names

# mm_hallmark_df <- data.frame(
#   gs_name = rep(names(mm_hallmark_sets), sapply(mm_hallmark_sets, length)),  # Pathway names
#   gene_symbol = unlist(mm_hallmark_sets)  # Flatten the list into a single vector
# )
pwl_msigdbr<-c(pwl_hallmark,pwl_kegg,pwl_reactome)
# # Combine the pathways into a single TERM2GENE dataset
# combined_pathways <- bind_rows(msig_h,kegg_gene_sets_msig,GO_gene_sets)
# colnames(combined_pathways) <- c("term", "gene")  # Rename columns to match clusterProfiler format
combined_pathways <- data.frame(
  gs_name = rep(names(pwl_msigdbr), sapply(pwl_msigdbr, length)),  # Pathway names
  gene_symbol = unlist(pwl_msigdbr)  # Flatten the list into a single vector
)

gsea_results_Mac3 <- GSEA(
  geneList = lfc_vector_Mac3,  # Ordered ranked gene list
  minGSSize = 5,          # Minimum gene set size
  maxGSSize = 500,        # Maximum gene set size
  pvalueCutoff = 0.1,       # p-value cutoff
  eps = 0,                # Boundary for calculating p-values
  seed = TRUE,            # Reproducibility
  pAdjustMethod = "BH",   # Benjamini-Hochberg correction
  TERM2GENE = combined_pathways  # Combined KEGG & MSigDB pathways
)
# Convert GSEA results to a data frame
gsea_results_Mac3_df <- as.data.frame(gsea_results_Mac3)


# Color gradient function for NES
color_gradient <- scale_fill_gradient2(
  low = "blue",      # Low NES values (negative)
  mid = "white",    # Midpoint (zero)
  high = "red",    # High NES values (positive)
  midpoint = 0,
  limits = c(min(gsea_results_Mac3_df$NES, na.rm = TRUE), max(gsea_results_Mac3_df$NES, na.rm = TRUE)),
  name = "NES"
)
top_gsea_Mac3 <- gsea_results_Mac3_df %>%
  arrange(desc(abs(NES))) %>%  # Sort by absolute NES
  head(30)  # Select top 20

# ggplot(top_gsea_Mac3, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
#   geom_bar(stat = "identity", show.legend = TRUE) +
#   coord_flip() +
#   color_gradient +
#   labs(title = "Top GSEA Results for Mac-1: Progressor vs Non-Progressor",
#        x = "Pathway",
#        y = "Normalized Enrichment Score (NES)") +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.7, size = 20, face = "bold"),
#     axis.title = element_text(size = 19),
#     axis.text = element_text(size = 12),
#     panel.grid.major = element_line(color = "grey90"),
#     panel.grid.minor = element_blank(),
#     legend.position = "bottom"
#   )

ggplot(top_gsea_Mac3, aes(x = reorder(Description, NES), y = NES, color = NES, size = p.adjust)) +
  geom_point(show.legend = TRUE) +  # Plot points with varying sizes
  scale_size_continuous(range = c(3, 10), name = "Adjusted p-value") +  # Control the range of bubble sizes
  scale_color_gradient2(
    low = "blue",      # Low NES values (negative)
    mid = "white",     # Midpoint (zero)
    high = "red",      # High NES values (positive)
    midpoint = 0,
    limits = c(min(gsea_results_Mac3_df$NES, na.rm = TRUE), max(gsea_results_Mac3_df$NES, na.rm = TRUE)),
    name = "NES"
  ) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1.5, color = "black") +  # Add bold line at y = 0 (NES=0)
  coord_flip() +
  labs(
    title = "GSEA Mac-2: Progressor vs Non-Progressor",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold", family = "Arial"),  # Title with larger size and better font
    axis.title = element_text(size = 20, face = "bold", family = "Arial"),  # Axis titles bold and larger font
    axis.text = element_text(size = 16, family = "Arial"),  # Axis text size increased and font changed
    panel.grid.major.x = element_blank(),  # Removed x-axis gridlines
    panel.grid.major.y = element_line(color = "grey90"),  # Horizontal gridlines retained
    panel.grid.minor = element_blank(),  # Minor gridlines removed
    legend.position = "right",  # Move legend to the right side
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # Bold legend title with improved font
    legend.text = element_text(size = 14, family = "Arial")  # Increased legend text size
  )


# Save results to CSV
write.csv(gsea_results_Mac3_df, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/DEG-Pathway/Macrophage/GSEA_Mac3_ProgressorVsNonProgressor_Week12.csv")


# Re-annotate main T1D_Timepoints Object ----
Macrophage<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/Macrophage_v1.rds")
T1D_Timepoints<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v4.rds")

DimPlot(Macrophage,group.by = "MacType")

# Ensure the cell names in Macrophage exist in T1D_Timepoints
cells_to_update <- colnames(Macrophage)

T1D_Timepoints$MacSubType<- T1D_Timepoints$CellSubType
DimPlot(T1D_Timepoints,group.by = "MacSubType")
# Ensure Tcells$TCellSubType is a factor
Macrophage$MacType <- as.character(Macrophage$MacType)

# Expand levels in T1D_Timepoints$MacSubType to include new ones
T1D_Timepoints$MacSubType <- factor(T1D_Timepoints$MacSubType, 
                                     levels = unique(c(levels(T1D_Timepoints$MacSubType), unique(Macrophage$MacType))))


# Assign values without generating invalid factor levels
T1D_Timepoints$MacSubType[cells_to_update] <- Macrophage$MacType

levels(T1D_Timepoints$MacSubType)
DimPlot(T1D_Timepoints,group.by = "MacSubType")
getwd()

saveRDS(T1D_Timepoints, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v5.rds")

