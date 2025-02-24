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
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(SingleCellExperiment))
BiocManager::install("scran")
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))

#' 
## DESeq and GSEA of Dendritic Cells ----
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
DC_Week6_PvNP = FindMarkers(DC_Week6, 
                            ident.1 = "Progressor", 
                            ident.2 = "Non-Progressor", 
                            assay = "SCT",
                            recorrect_umi = FALSE)
# Apply Benjamini-Hochberg correction
DC_Week6_PvNP$p_val_adj <- p.adjust(DC_Week6_PvNP$p_val, method = "BH")

DC_Week12_PvNP = FindMarkers(DC_Week12,
                             ident.1 = "Progressor",
                             ident.2 = "Non-Progressor",
                             assay= "SCT",
                             recorrect_umi = FALSE)
DC_Week12_PvNP$p_val_adj <- p.adjust(DC_Week12_PvNP$p_val, method = "BH")

colnames(DC_Week6_PvNP)

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
pathways <- c(msig_h)



DC_Week6$group
DC_W6_Progressor <- seurat_extract(DC_Week6,
                      meta1 = "group", value_meta1 = "Progressor")
DC_W6_NonProgressor <- seurat_extract(DC_Week6,
                      meta1 = "group", value_meta1 = "Non-Progressor")
scpa_out <- compare_pathways(samples = list(DC_W6_Progressor, DC_W6_NonProgressor), 
                             pathways = pathways)
head(scpa_out)
plot_rank(scpa_out = scpa_out, 
          pathway = "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
          base_point_size = 2, 
          highlight_point_size = 3)

scpa_out <- scpa_out %>%
  filter(FC > 0) %>%  # Remove rows where FC is non-positive
  mutate(
    log10_pval = -log10(Pval),         # Calculate -log10(p-value)
    log2_fc = log2(FC)                 # Calculate log2(FC)
  )

library(ggplot2)

library(ggplot2)

# Highlight additional pathways
highlight_pathways <- c("")

library(ggplot2)
library(RColorBrewer)

# Create the scatter plot with subtle gridlines
ggplot(scpa_out, aes(x = log2_fc, y = log10_pval)) +
  # Plot points, color by log2(FC), with a continuous color scale
  geom_point(aes(color = log2_fc), size = 5) +  # Increased dot size
  
  # Highlight the selected pathways with labels on the left side
  geom_text(aes(label = ifelse(Pathway %in% highlight_pathways, Pathway, ""),
                hjust = 1.1, vjust = 0), size = 4, color = "red") +  # Move label to left
  
  # Add axis lines
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "lightgrey", size = 0.5),  # Light grey gridlines for major ticks
    panel.grid.minor = element_blank(),  # No minor gridlines
    axis.line = element_line(color = "black", size = 0.5),  # Black solid axis lines
    legend.position = "right",  # Show legend for color
    text = element_text(size = 14),  # Increase text size for readability
    axis.title = element_text(size = 16),  # Larger axis title font size
    axis.text = element_text(size = 12)  # Larger axis text font size
  ) +
  
  # Apply color gradient based on log2(FC) using the YlOrRd color palette
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd"), 
                        name = "log2(Fold Change)") +  # Add label for the color legend
  
  # Add labels and title
  labs(
    x = expression("log"[2] * "(Fold Change)"),
    y = expression("-log"[10] * "(p-value)")
  )

