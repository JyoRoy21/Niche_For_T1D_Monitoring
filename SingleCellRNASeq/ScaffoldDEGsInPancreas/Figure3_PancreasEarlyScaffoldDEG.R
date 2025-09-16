# 1. Load Libraries ----
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

# 2. Load Single Cell Object ----
getwd()
T1D_Timepoints = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v4.rds")
unique(T1D_Timepoints$CellSubType)
Idents(T1D_Timepoints)<-T1D_Timepoints$CellSubType
week6_data <- subset(T1D_Timepoints,subset =time  %in% c("Week6"))
week6_data$sample_id <- ifelse(week6_data$group == "Non-Progressor",
                               paste0("NonProgressor", gsub("Week6_", "", week6_data$sample)),
                               paste0("Progressor", gsub("Week6_", "", week6_data$sample)))

# Extract raw counts and metadata to create SingleCellExperiment object
counts_w6 <- week6_data@assays$RNA$counts 
metadata_w6 <- week6_data@meta.data
# Set up metadata as desired for aggregation and DE analysis
metadata_w6$cluster_id <- factor(week6_data@active.ident)

# Create single cell experiment object
sce_w6 <- SingleCellExperiment(assays = list(counts = counts_w6), 
                            colData = metadata_w6)

dim(colData(sce_w6))
head(colData(sce_w6))

# 3. Expression of Early Scaffold Genes in CellTypes ----
ScafGenes <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/ScaffoldDESEQ/DEGAnalysis_Volcano_Early.csv", sep=",", header=T) # Metadata file
rownames(ScafGenes)<-ScafGenes$X
ScafGenes$X<-NULL
colnames(ScafGenes)
# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])

# Ensure selected genes exist in the dataset
selected_genes <- selected_genes[selected_genes %in% rownames(week6_data)]
library(Seurat)
# Generate the heatmap
DoHeatmap(
  week6_data,
  features = selected_genes,
  #group.by = "cell_type",  # Replace with your cell type metadata column
  slot = "scale.data"  # Uses scaled expression data
) + scale_fill_gradientn(colors = c("blue", "white", "red"))


dotplot <- DotPlot(
  week6_data,
  features = selected_genes,
  #group.by = "cell_type",
  dot.scale = 8
) + RotatedAxis()  # Rotate x-axis labels for better visibility

# Adjust color gradient manually
dotplot + scale_color_gradientn(colors = c("blue", "white", "red"))+  
  scale_size_continuous(range = c(1, 8))  +  # Min size = 1, Max size = 8
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Increase x-axis text size & rotate
    axis.text.y = element_text(size = 16),  # Increase y-axis text size
    axis.title.x = element_text(size = 18, face = "bold"),  # Increase x-axis title size
    axis.title.y = element_text(size = 18, face = "bold")   # Increase y-axis title size
  )+
  labs(
    x = "Immunological Niche Based Differentially Expressed Genes-Early Stage",  # Custom x-axis label
    y = "Pancreas Cell Types"  # Custom y-axis label
  )

# 5. DE Genes in Cell Types ----


## I. Macrophage ----
Macrophage_Week6 <- subset(week6_data, idents = "Macrophage")
Idents(Macrophage_Week6) <- Macrophage_Week6$group
Macrophage_Week6 = PrepSCTFindMarkers(Macrophage_Week6)

Macrophage_Week6_PvNP = FindMarkers(Macrophage_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                    recorrect_umi = FALSE)

colnames(ScafGenes)
# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])
ScafDEGenes<-ScafGenes[selected_genes,]
# Filter for DE genes in the Macrophage_Week6_PvNP dataset
Macrophage_Week6_PvNP_filtered <- Macrophage_Week6_PvNP[
  (Macrophage_Week6_PvNP$p_val <= 0.05) & 
    (abs(Macrophage_Week6_PvNP$avg_log2FC) >= 1), 
]

# Ensure that the trends (direction) of fold change are the same (both ≥1 or both ≤-1)
same_trend_genes <- intersect(
  rownames(ScafDEGenes[ScafDEGenes$log2FoldChange >= 1,]),
  rownames(Macrophage_Week6_PvNP_filtered[Macrophage_Week6_PvNP_filtered$avg_log2FC >= 1,])
)

# Combine genes with the same trend for both datasets (log2FoldChange <= -1)
same_trend_genes <- union(
  same_trend_genes,
  intersect(
    rownames(ScafDEGenes[ScafDEGenes$log2FoldChange <= -1,]),
    rownames(Macrophage_Week6_PvNP_filtered[Macrophage_Week6_PvNP_filtered$avg_log2FC <= -1,])
  )
)
# Display the common DE genes
same_trend_genes
VlnPlot(Macrophage_Week6, features = c("Blk","Ctla4","Ms4a1"))

options(future.globals.maxSize = 5000 * 1024^2)  # Increase to 5GB

## II. TRegs ----
Treg_Week6 <- subset(week6_data, idents = "Tregs")
Idents(Treg_Week6) <- Treg_Week6$group
Treg_Week6 = PrepSCTFindMarkers(Treg_Week6)

Treg_Week6_PvNP = FindMarkers(Treg_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                    recorrect_umi = FALSE)

# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])
ScafDEGenes<-ScafGenes[selected_genes,]
# Filter for DE genes in the Treg_Week6_PvNP dataset
Treg_Week6_PvNP_filtered <- Treg_Week6_PvNP[
  (Treg_Week6_PvNP$p_val <= 0.05) & 
    (abs(Treg_Week6_PvNP$avg_log2FC) >= 1), 
]

# Ensure that the trends (direction) of fold change are the same (both ≥1 or both ≤-1)
same_trend_genes <- intersect(
  rownames(ScafDEGenes[ScafDEGenes$log2FoldChange>= 1,]),
  rownames(Treg_Week6_PvNP_filtered[Treg_Week6_PvNP_filtered$avg_log2FC >= 1,])
)

# Combine genes with the same trend for both datasets (log2FoldChange <= -1)
same_trend_genes <- union(
  same_trend_genes,
  intersect(
    rownames(ScafDEGenes[ScafDEGenes$log2FoldChange <= -1,]),
    rownames(Treg_Week6_PvNP_filtered[Treg_Week6_PvNP_filtered$avg_log2FC <= -1,])
  )
)
# Display the common DE genes
same_trend_genes
VlnPlot(Treg_Week6, features = c("Ms4a1"))

## III. NKCells ----
NKCell_Week6 <- subset(week6_data, idents = "NK Cell")
Idents(NKCell_Week6) <- NKCell_Week6$group
NKCell_Week6 = PrepSCTFindMarkers(NKCell_Week6)

NKCell_Week6_PvNP = FindMarkers(NKCell_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                              recorrect_umi = FALSE)

# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])
ScafDEGenes<-ScafGenes[selected_genes,]
# Filter for DE genes in the NKCell_Week6_PvNP dataset
NKCell_Week6_PvNP_filtered <- NKCell_Week6_PvNP[
  (NKCell_Week6_PvNP$p_val <= 0.05) & 
    (abs(NKCell_Week6_PvNP$avg_log2FC) >= 1), 
]

# Ensure that the trends (direction) of fold change are the same (both ≥1 or both ≤-1)
same_trend_genes <- intersect(
  rownames(ScafDEGenes[ScafDEGenes$log2FoldChange>= 1,]),
  rownames(NKCell_Week6_PvNP_filtered[NKCell_Week6_PvNP_filtered$avg_log2FC >= 1,])
)

# Combine genes with the same trend for both datasets (log2FoldChange <= -1)
same_trend_genes <- union(
  same_trend_genes,
  intersect(
    rownames(ScafDEGenes[ScafDEGenes$log2FoldChange <= -1,]),
    rownames(NKCell_Week6_PvNP_filtered[NKCell_Week6_PvNP_filtered$avg_log2FC <= -1,])
  )
)
# Display the common DE genes
same_trend_genes
VlnPlot(NKCell_Week6, features = c("Blk"))

## IV. Tcon Activated ----
TconActivate_Week6 <- subset(week6_data, idents = "Tcon activated ")
Idents(TconActivate_Week6) <- TconActivate_Week6$group
TconActivate_Week6 = PrepSCTFindMarkers(TconActivate_Week6)

TconActivate_Week6_PvNP = FindMarkers(TconActivate_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                recorrect_umi = TRUE)

# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])
ScafDEGenes<-ScafGenes[selected_genes,]
# Filter for DE genes in the TconActivate_Week6_PvNP dataset
TconActivate_Week6_PvNP_filtered <- TconActivate_Week6_PvNP[
  (TconActivate_Week6_PvNP$p_val <= 0.05) & 
    (abs(TconActivate_Week6_PvNP$avg_log2FC) >= 1), 
]

# Ensure that the trends (direction) of fold change are the same (both ≥1 or both ≤-1)
same_trend_genes <- intersect(
  rownames(ScafDEGenes[ScafDEGenes$log2FoldChange>= 1,]),
  rownames(TconActivate_Week6_PvNP_filtered[TconActivate_Week6_PvNP_filtered$avg_log2FC >= 1,])
)

# Combine genes with the same trend for both datasets (log2FoldChange <= -1)
same_trend_genes <- union(
  same_trend_genes,
  intersect(
    rownames(ScafDEGenes[ScafDEGenes$log2FoldChange <= -1,]),
    rownames(TconActivate_Week6_PvNP_filtered[TconActivate_Week6_PvNP_filtered$avg_log2FC <= -1,])
  )
)
# Display the common DE genes
same_trend_genes
VlnPlot(TconActivate_Week6, features = c("Ms4a1"))


## V. ILC2 ----
ILC2_Week6 <- subset(week6_data, idents = "ILC2")
Idents(ILC2_Week6) <- ILC2_Week6$group
# Check the class and structure of the object
class(ILC2_Week6)

# Check the dimensions of the object (rows = genes, columns = cells)
dim(ILC2_Week6)

ILC2_Week6 = PrepSCTFindMarkers(ILC2_Week6)

ILC2_Week6_PvNP = FindMarkers(ILC2_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                                      recorrect_umi = FALSE)

# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])
ScafDEGenes<-ScafGenes[selected_genes,]
# Filter for DE genes in the ILC2_Week6_PvNP dataset
ILC2_Week6_PvNP_filtered <- ILC2_Week6_PvNP[
  (ILC2_Week6_PvNP$p_val <= 0.05) & 
    (abs(ILC2_Week6_PvNP$avg_log2FC) >= 1), 
]

# Ensure that the trends (direction) of fold change are the same (both ≥1 or both ≤-1)
same_trend_genes <- intersect(
  rownames(ScafDEGenes[ScafDEGenes$log2FoldChange>= 1,]),
  rownames(ILC2_Week6_PvNP_filtered[ILC2_Week6_PvNP_filtered$avg_log2FC >= 1,])
)

# Combine genes with the same trend for both datasets (log2FoldChange <= -1)
same_trend_genes <- union(
  same_trend_genes,
  intersect(
    rownames(ScafDEGenes[ScafDEGenes$log2FoldChange <= -1,]),
    rownames(ILC2_Week6_PvNP_filtered[ILC2_Week6_PvNP_filtered$avg_log2FC <= -1,])
  )
)
# Display the common DE genes
same_trend_genes
VlnPlot(ILC2_Week6, features = c("Cntn1","Lipe","Cd6","Ctla4","Ms4a1"))

## V. Tcon ----
TconMem_Week6 <- subset(week6_data, idents = "Tcon memory")
Idents(TconMem_Week6) <- TconMem_Week6$group
TconMem_Week6 = PrepSCTFindMarkers(TconMem_Week6)

TconMem_Week6_PvNP = FindMarkers(TconMem_Week6, ident.1 = "Progressor", ident.2 = "Non-Progressor", assay= "SCT",
                              recorrect_umi = FALSE)

# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])
ScafDEGenes<-ScafGenes[selected_genes,]
# Filter for DE genes in the TconMem_Week6_PvNP dataset
TconMem_Week6_PvNP_filtered <- TconMem_Week6_PvNP[
  (TconMem_Week6_PvNP$p_val <= 0.05) & 
    (abs(TconMem_Week6_PvNP$avg_log2FC) >= 1), 
]

# Ensure that the trends (direction) of fold change are the same (both ≥1 or both ≤-1)
same_trend_genes <- intersect(
  rownames(ScafDEGenes[ScafDEGenes$log2FoldChange>= 1,]),
  rownames(TconMem_Week6_PvNP_filtered[TconMem_Week6_PvNP_filtered$avg_log2FC >= 1,])
)

# Combine genes with the same trend for both datasets (log2FoldChange <= -1)
same_trend_genes <- union(
  same_trend_genes,
  intersect(
    rownames(ScafDEGenes[ScafDEGenes$log2FoldChange <= -1,]),
    rownames(TconMem_Week6_PvNP_filtered[TconMem_Week6_PvNP_filtered$avg_log2FC <= -1,])
  )
)
# Display the common DE genes
same_trend_genes
VlnPlot(TconMem_Week6, features = c("Ms4a1"))


# 4.Preparing the single-cell dataset for pseudobulk analysis -----
# Extract unique names of clusters (= levels of cluster_id factor variable)
cluster_names <- levels(colData(sce_w6)$cluster_id)
cluster_names

# Extract unique names of samples (= levels of sample_id factor variable)
sample_names <- unique(colData(sce_w6)$sample_id)
sample_names

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce_w6)[, c("cluster_id", "sample_id")]
head(groups)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce_w6)), 
                                groupings = groups, fun = "sum") 

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]

# Understanding tstrsplit()

# Using which() to look up tstrsplit() output
Mac_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "Macrophage")
Mac_idx

colnames(aggr_counts)[Mac_idx]
aggr_counts[1:10, Mac_idx]

# As a reminder, we stored our cell types in a vector called cluster_names
cluster_names

# Loop over all cell types to extract corresponding counts, and store information in a list

## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)

head(colData(sce_w6))

# Extract sample-level variables
metadata <- colData(sce_w6) %>% 
  as.data.frame() %>% 
  dplyr::select(group,sample_id)
dim(metadata)
head(metadata)
# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]
dim(metadata)
head(metadata)
# Rename rows
rownames(metadata) <- metadata$sample_id
head(metadata)

# Number of cells per sample and cluster
t <- table(colData(sce_w6)$sample_id,
           colData(sce_w6)$cluster_id)

# Creating metadata list

## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}

# Explore the different components of the list
str(metadata_ls)  
# Macrophage ----

## A. DESeq2 Analysis- Sample Level ----

idx <- which(names(counts_ls) == "Macrophage")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Check contents of extracted objects
cluster_counts[1:6, 1:6]
head(cluster_metadata)


# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group)
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
# Perform PCA and store the ggplot object
pca_data <- DESeq2::plotPCA(rld, ntop = 100, intgroup = "group", returnData = TRUE)

# Load ggplot2 for custom labeling
library(ggplot2)

# Calculate proportion of variance explained
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Generate a high-quality PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 6, alpha = 0.8) +  # Increase point size
  geom_text_repel(size = 6, fontface = "bold", box.padding = 0.5) +  # Smarter text labeling
  stat_ellipse(level = 0.68, size = 1.2, linetype = "dashed") +  # 68% confidence interval (One SEM)
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +  # Custom color palette
  theme_classic(base_size = 18) +  # Increase base font size for publication
  labs(
    title = "PCA of Samples",
    x = paste0("PC1 (", percentVar[1], "% Variance)"),
    y = paste0("PC2 (", percentVar[2], "% Variance)"),
    color = "Group"
  ) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5)
  )
# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("group"), drop=F])

## B. DESeq2 Analysis- Main ----

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
# Plot dispersion estimates
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# # Generate results object
# res <- results(dds, 
#                name = "group_Progressor_vs_Non.Progressor",
#                alpha = 0.05)
# 
# # Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
# res <- lfcShrink(dds, 
#                  coef = "group_Progressor_vs_Non.Progressor",
#                  res=res,
#                  type = "apeglm")
# 
# # Turn the DESeq2 results object into a tibble for use with tidyverse functions
# res_tbl <- res %>%
#   data.frame() %>%
#   rownames_to_column(var = "gene") %>%
#   as_tibble() %>%
#   arrange(padj)
# 
# # Check results output
# res_tbl 

results_macrophage <- as.data.frame(results(dds, contrast = c("group", "Progressor", "Non-Progressor")))

# Replace NA p values with 1
results_macrophage$pvalue[is.na(results_macrophage$pvalue)] <- 1
results_macrophage$padj[is.na(results_macrophage$padj)] <- 1

DE_macrophage_genes <- rownames(results_macrophage[
  results_macrophage$pvalue < 0.05 & 
    (abs(results_macrophage$log2FoldChange) >= 1), 
])

# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])
# 3. Find common genes between DE macrophage genes and immunological niche genes
common_genes <- intersect(DE_macrophage_genes, selected_genes)
# 4. Show the common genes
print(common_genes)

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)



# Define color mapping with stronger contrast for visibility
keyvals <- ifelse(results_macrophage$pvalue < 0.05 & results_macrophage$log2FoldChange > 1, "#7E0000",  # Dark red for FDR-significant up
                  ifelse(results_macrophage$pvalue < 0.05 & results_macrophage$log2FoldChange < -1, "#002E66",  # Dark blue for FDR-significant down
                                         "gray"))  # Gray for non-significant genes

# Assign names for legend
names(keyvals) <- ifelse(results_macrophage$pvalue < 0.05 & results_macrophage$log2FoldChange > 1, "Upregulated-Progressor",
                         ifelse(results_macrophage$pvalue < 0.05 & results_macrophage$log2FoldChange < -1, "Upregulated-Non-Progressor",
                                              "Not Significant"))



# Genes to label
# genes_to_label_early <- c("S100a8", "S100a9", "Timp4", "Trpv1", "Ctla4", "Ms4a1",
#                           "Ffar2", "Adipoq", "Gys2", "Lipe", "Pck1", "Scd1", "Tshr")

EnhancedVolcano(results_macrophage,
                lab = rownames(results_macrophage),
                #selectLab = genes_to_label_early,
                x = 'log2FoldChange',
                y = 'pvalue',  
                title = 'Macropage-Early Timepoint',
                pCutoff = 0.05,  
                FCcutoff = 1, 
                boxedLabels = TRUE,
                #pointSize = point_sizes,  # Increase size for FDR-significant points
                labSize = 4.5,
                colAlpha = 0.75,
                legendLabels = c("Not Significant", "Upregulated-Progressor", "Upregulated-Non-Progressor"),
                legendPosition = 'right',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p-value)"),  # Subscript in y-axis
                xlim = c(-8, 8),  # Set x-axis from -8 to 8
                ylim = c(0, 6),  # Set y-axis limit to 6
                colCustom = keyvals,  # Apply custom colors
                caption = NULL  # Removes watermark text
)

# ILC3 ----

## A. DESeq2 Analysis- Sample Level ----

idx <- which(names(counts_ls) == "ILC2")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Check contents of extracted objects
cluster_counts[1:6, 1:6]
head(cluster_metadata)


# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group)
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
# Perform PCA and store the ggplot object
pca_data <- DESeq2::plotPCA(rld, ntop = 100, intgroup = "group", returnData = TRUE)

# Load ggplot2 for custom labeling
library(ggplot2)
#table(week6_data$CellSubType,week6_data$sample)
# Calculate proportion of variance explained
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Generate a high-quality PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 6, alpha = 0.8) +  # Increase point size
  geom_text_repel(size = 6, fontface = "bold", box.padding = 0.5) +  # Smarter text labeling
  stat_ellipse(level = 0.68, size = 1.2, linetype = "dashed") +  # 68% confidence interval (One SEM)
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +  # Custom color palette
  theme_classic(base_size = 18) +  # Increase base font size for publication
  labs(
    title = "PCA of Samples",
    x = paste0("PC1 (", percentVar[1], "% Variance)"),
    y = paste0("PC2 (", percentVar[2], "% Variance)"),
    color = "Group"
  ) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5)
  )

## B. DESeq2 Analysis- Main ----

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
# Plot dispersion estimates
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# # Generate results object
# res <- results(dds, 
#                name = "group_Progressor_vs_Non.Progressor",
#                alpha = 0.05)
# 
# # Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
# res <- lfcShrink(dds, 
#                  coef = "group_Progressor_vs_Non.Progressor",
#                  res=res,
#                  type = "apeglm")
# 
# # Turn the DESeq2 results object into a tibble for use with tidyverse functions
# res_tbl <- res %>%
#   data.frame() %>%
#   rownames_to_column(var = "gene") %>%
#   as_tibble() %>%
#   arrange(padj)
# 
# # Check results output
# res_tbl 

results_ILC2 <- as.data.frame(results(dds, contrast = c("group", "Progressor", "Non-Progressor")))

# Replace NA p values with 1
results_ILC2$pvalue[is.na(results_ILC2$pvalue)] <- 1
results_ILC2$padj[is.na(results_ILC2$padj)] <- 1

DE_ILC2_genes <- rownames(results_ILC2[
  results_ILC2$pvalue < 0.05 & 
    (abs(results_ILC2$log2FoldChange) >= 1), 
])

# Select significant genes
selected_genes <- rownames(ScafGenes[ScafGenes$pvalue < 0.05 & abs(ScafGenes$log2FoldChange) >= 1, ])
# 3. Find common genes between DE macrophage genes and immunological niche genes
common_genes <- intersect(DE_ILC2_genes, selected_genes)
# 4. Show the common genes
print(common_genes)

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)



# Define color mapping with stronger contrast for visibility
keyvals <- ifelse(results_macrophage$pvalue < 0.05 & results_macrophage$log2FoldChange > 1, "#7E0000",  # Dark red for FDR-significant up
                  ifelse(results_macrophage$pvalue < 0.05 & results_macrophage$log2FoldChange < -1, "#002E66",  # Dark blue for FDR-significant down
                         "gray"))  # Gray for non-significant genes

# Assign names for legend
names(keyvals) <- ifelse(results_macrophage$pvalue < 0.05 & results_macrophage$log2FoldChange > 1, "Upregulated-Progressor",
                         ifelse(results_macrophage$pvalue < 0.05 & results_macrophage$log2FoldChange < -1, "Upregulated-Non-Progressor",
                                "Not Significant"))



# Genes to label
# genes_to_label_early <- c("S100a8", "S100a9", "Timp4", "Trpv1", "Ctla4", "Ms4a1",
#                           "Ffar2", "Adipoq", "Gys2", "Lipe", "Pck1", "Scd1", "Tshr")

EnhancedVolcano(results_macrophage,
                lab = rownames(results_macrophage),
                #selectLab = genes_to_label_early,
                x = 'log2FoldChange',
                y = 'pvalue',  
                title = 'Macropage-Early Timepoint',
                pCutoff = 0.05,  
                FCcutoff = 1, 
                boxedLabels = TRUE,
                #pointSize = point_sizes,  # Increase size for FDR-significant points
                labSize = 4.5,
                colAlpha = 0.75,
                legendLabels = c("Not Significant", "Upregulated-Progressor", "Upregulated-Non-Progressor"),
                legendPosition = 'right',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  # Subscript in x-axis
                ylab = expression("-log"[10] ~ "(p-value)"),  # Subscript in y-axis
                xlim = c(-8, 8),  # Set x-axis from -8 to 8
                ylim = c(0, 6),  # Set y-axis limit to 6
                colCustom = keyvals,  # Apply custom colors
                caption = NULL  # Removes watermark text
)

