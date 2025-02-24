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
## DESeq and SCPA of Macrophages ----
T1D_Timepoints<-readRDS("./annotated_T1D_Timepoints_v4.rds")

DimPlot(T1D_Timepoints)
T1D_Timepoints$time
Macrophage_Week6 <- subset(T1D_Timepoints, idents = "Macrophage", subset = time %in% "Week6")
Macrophage_Week12 <- subset(T1D_Timepoints, idents = "Macrophage", subset = time %in% "Week12")

table(Macrophage_Week6$group)
Idents(Macrophage_Week6)
Idents(Macrophage_Week6) <- Macrophage_Week6$group
Idents(Macrophage_Week12) <- Macrophage_Week12$group


Macrophage_Week6 = PrepSCTFindMarkers(Macrophage_Week6)
Macrophage_Week12 = PrepSCTFindMarkers(Macrophage_Week12)

Macrophage_Week6_PvNP = FindMarkers(Macrophage_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                            recorrect_umi = FALSE)
Macrophage_Week6_PvNP$p_val_adj<-p.adjust(Macrophage_Week6_PvNP$p_val, method = "BH")


Macrophage_Week12_PvNP = FindMarkers(Macrophage_Week12, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                             recorrect_umi = FALSE)
Macrophage_Week12_PvNP$p_val_adj<-p.adjust(Macrophage_Week12_PvNP$p_val, method = "BH")

Macrophage_Week6_PvNP$symbol = rownames(Macrophage_Week6_PvNP)
Macrophage_Week12_PvNP$symbol = rownames(Macrophage_Week12_PvNP)

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

EnhancedVolcano(Macrophage_Week6_PvNP_data,
                lab = rownames(Macrophage_Week6_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(Macrophage_Week6_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "Macrophage: Progressor vs Non-Progressor at Week 6",
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
ggsave("Volcanoplot_Macrophage_T1DTimepoints_Week6.png", width = 12, height = 8, dpi = 600)

#dev.off()
write.csv(Macrophage_Week6_PvNP_data,"Macrophage_T1DTimepoints_DEG_Week6.csv",row.names = TRUE)

#### Week 12----
Macrophage_Week12_PvNP_data = Macrophage_Week12_PvNP

keyvals <- ifelse(
  (Macrophage_Week12_PvNP_data$avg_log2FC < -1 & Macrophage_Week12_PvNP_data$p_val_adj < 0.05), '#1465AC',  # Non-Progressor (Blue)
  ifelse((Macrophage_Week12_PvNP_data$avg_log2FC > 1 & Macrophage_Week12_PvNP_data$p_val_adj < 0.05), '#B31B21',  # Progressor (Red)
         'grey'))  # Not DEG

keyvals[is.na(keyvals)] <- 'grey'

# Assign names
names(keyvals)[keyvals == 'grey'] <- 'Non Significant'
names(keyvals)[keyvals == '#B31B21'] <- 'Upregulated- Progressor'
names(keyvals)[keyvals == '#1465AC'] <- 'Uregulated- Non-Progressor'
 

EnhancedVolcano(Macrophage_Week12_PvNP_data,
                lab = rownames(Macrophage_Week12_PvNP_data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                pCutoffCol="p_val_adj",
                FCcutoff = 1,
                colCustom = keyvals,
                labSize = 5.5,
                colAlpha = 0.75,
                selectLab = rownames(Macrophage_Week12_PvNP_data)[which(names(keyvals) %in% c('Upregulated- Progressor', 'Uregulated- Non-Progressor'))],
                title = "Macrophage: Progressor vs Non-Progressor at Week 12",
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
ggsave("Volcanoplot_Macrophage_T1DTimepoints_Week12.png", width = 12, height = 8, dpi = 600)


write.csv(Macrophage_Week12_PvNP_data,"Macrophage_T1DTimepoints_DEG_Week12.csv",row.names = TRUE)




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

