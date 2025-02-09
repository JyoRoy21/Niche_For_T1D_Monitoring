## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq")
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
#' 
# DESeq and GSEA of Dendritic Cells ----
T1D_Timepoints<-readRDS("./annotated_T1D_Timepoints_v4.rds")

DimPlot(T1D_Timepoints)
T1D_Timepoints$time
DC_Week6 <- subset(T1D_Timepoints, idents = "Dendritic Cell", subset = time %in% "Week6")
DC_Week12 <- subset(T1D_Timepoints, idents = "Dendritic Cell", subset = time %in% "Week12")


Idents(DC_Week6) <- DC_Week6$group
Idents(DC_Week12) <- DC_Week12$group


DC_Week6 = PrepSCTFindMarkers(DC_Week6)
DC_Week12 = PrepSCTFindMarkers(DC_Week12)
table(DC_Week6$group)
Idents(DC_Week6)
DC_Week6_PvNP = FindMarkers(DC_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                            recorrect_umi = FALSE,p.adjust.method = "none")
DC_Week12_PvNP = FindMarkers(DC_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                             recorrect_umi = FALSE)
rownames(DC_Week12_PvNP)

DC_Week6_PvNP$symbol = rownames(DC_Week6_PvNP)
DC_Week12_PvNP$symbol = rownames(DC_Week12_PvNP)

### Volcano Plot ----

#### Week 6----
DC_Week6_PvNP_data = DC_Week6_PvNP

keyvals <- ifelse(
  (DC_Week6_PvNP_data$avg_log2FC < -1 & DC_Week6_PvNP_data$p_val < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((DC_Week6_PvNP_data$avg_log2FC > 1 & DC_Week6_PvNP_data$p_val < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated'
names(keyvals)[keyvals == '#1465AC'] <- 'Downregulated'


#png(file = "AllergyScaf_NPPBSCh4_VolcanoPlot.png",units = "in",width = 10, height = 10, res = 400)
EnhancedVolcano(DC_Week6_PvNP_data,
                lab = rownames(DC_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.05,
                FCcutoff = 1,
                colCustom = keyvals,
                selectLab = rownames(DC_Week6_PvNP_data)[which(names(keyvals) %in% c('Upregulated', 'Downregulated'))],
                title = "Dendritic Cells: Progressor vs Non-Progressor at Week 6",
                max.overlaps = 15,
                drawConnectors = T,
                widthConnectors = 0.75,
                ylim = c(0,4),
                boxedLabels = F)
dev.off()
write.csv(DC_Week6_PvNP_data,"DendriticCells_DEG_Week6.csv",row.names = TRUE)

#### Week 12----
DC_Week12_PvNP_data = DC_Week12_PvNP

keyvals <- ifelse(
  (DC_Week12_PvNP_data$avg_log2FC < -1 & DC_Week12_PvNP_data$p_val < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((DC_Week12_PvNP_data$avg_log2FC > 1 & DC_Week12_PvNP_data$p_val < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated'
names(keyvals)[keyvals == '#1465AC'] <- 'Downregulated'


#png(file = "AllergyScaf_NPPBSCh4_VolcanoPlot.png",units = "in",width = 10, height = 10, res = 400)
EnhancedVolcano(DC_Week12_PvNP_data,
                lab = rownames(DC_Week12_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.05,
                FCcutoff = 1,
                colCustom = keyvals,
                selectLab = rownames(DC_Week12_PvNP_data)[which(names(keyvals) %in% c('Upregulated', 'Downregulated'))],
                title = "Dendritic Cells: Progressor vs Non-Progressor at Week 12",
                max.overlaps = 20,
                drawConnectors = T,
                widthConnectors = 0.75,
                ylim = c(0,3),
                boxedLabels = F)
dev.off()

write.csv(DC_Week12_PvNP_data,"DendriticCells_DEG_Week12.csv",row.names = TRUE)

### GSEA Analysis ----

DC_Week6_PvNP_geneList =  DC_Week6_PvNP %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(symbol, avg_log2FC)

DC_Week12_PvNP_geneList =  DC_Week12_PvNP %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(symbol, avg_log2FC)

DC_Week6_PvNP_ranks<- deframe(DC_Week6_PvNP_geneList)
DC_Week12_PvNP_ranks<- deframe(DC_Week12_PvNP_geneList)

# Define Gene Sets for Enrichment Analysis
all_gene_sets = msigdbr(species = "Mus musculus")
unique(all_gene_sets$gs_subcat)
unique(all_gene_sets$gs_cat)
kegg_gene_sets_msig = msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG")
biocarta_gene_sets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:BIOCARTA")
Reactome_gene_sets = msigdbr(species = "mouse", subcategory = "CP:REACTOME")
Wiki_gene_sets = msigdbr(species = "mouse", subcategory = "CP:WIKIPATHWAYS")
HM_gene_sets = msigdbr(species = "mouse", category = "H")
GO_gene_sets = msigdbr(species = "mouse", subcategory = "GO:BP")
combined_data_set = rbind(kegg_gene_sets_msig, 
                          HM_gene_sets)

pathway_gene_set = combined_data_set[,c("gs_name", "gene_symbol")]

DC_Week6_PvNP_GSEA = GSEA(DC_Week6_PvNP_ranks,
         pvalueCutoff = 0.5,
         pAdjustMethod = "fdr",
         TERM2GENE = pathway_gene_set,
         minGSSize = 2
         
)

y1 = GSEA(DC_Week12_PvNP_ranks,
          pvalueCutoff = 0.5,
          pAdjustMethod = "fdr",
          TERM2GENE = pathway_gene_set,
          minGSSize = 2
          
)



# Convert GSEA results into a data frame
DC_Week6_PvNP_results <- as.data.frame(DC_Week6_PvNP_GSEA@result)

# Create dot plot
ggplot(DC_Week6_PvNP_results, aes(x = NES, y = reorder(Description, NES))) +
  geom_point(aes(size = -log10(pvalue), color = NES)) +
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  scale_size(range = c(2, 10)) +  # Adjust dot sizes
  theme_minimal() +
  labs(
    title = "GSEA Dot Plot - DC Week 6 (PBS vs NP)",
    x = "Normalized Enrichment Score (NES)",
    y = "Pathways",
    size = "-log10(p-value)",
    color = "NES"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.position = "right"
  )

head(y)

y_filt = y[which(abs(y$qvalue)<0.25),]
#y_filt = as.data.frame(y)

y_filt$Condition = "Non-Progressor"
y_filt[which((y_filt$NES)>0),]$Condition = "Progressor" 

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

y



