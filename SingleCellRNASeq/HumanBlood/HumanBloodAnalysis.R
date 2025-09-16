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
T1D_Blood = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/T1D_Seurat_Object_Final.rds")
T1D_Blood$Cluster_Annotation_All <- gsub("_", "-", T1D_Blood$Cluster_Annotation_All)
T1D_Blood <- FindVariableFeatures(T1D_Blood, selection.method = "vst", nfeatures = 2000)
T1D_Blood <- RunPCA(T1D_Blood, features = VariableFeatures(object = T1D_Blood))
ElbowPlot(T1D_Blood)
T1D_Blood <- RunUMAP(T1D_Blood, dims = 1:20)
DimPlot(T1D_Blood)
unique(T1D_Blood$Sample_ID)
unique(T1D_Blood$COND)
unique(T1D_Blood$Cluster_Annotation_All)
table(T1D_Blood$Cluster_Annotation_All,T1D_Blood$COND)
Idents(T1D_Blood)<-T1D_Blood$Cluster_Annotation_All
DimPlot(T1D_Blood)

# Extract raw counts and metadata to create SingleCellExperiment object
counts_blood <- T1D_Blood@assays$RNA$counts 
metadata_blood <- T1D_Blood@meta.data
# Set up metadata as desired for aggregation and DE analysis
metadata_blood$cluster_id <- factor(T1D_Blood@active.ident)

# Create single cell experiment object
sce_blood <- SingleCellExperiment(assays = list(counts = counts_blood), 
                               colData = metadata_blood)

dim(colData(sce_blood))
head(colData(sce_blood))


# 4.Preparing the single-cell dataset for pseudobulk analysis -----
# Extract unique names of clusters (= levels of cluster_id factor variable)
cluster_names <- levels(colData(sce_blood)$cluster_id)
cluster_names

# Extract unique names of samples (= levels of sample_id factor variable)
sample_names <- unique(colData(sce_blood)$Sample_ID)
sample_names

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce_blood)[, c("cluster_id", "Sample_ID")]
head(groups)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce_blood)), 
                                groupings = groups, fun = "sum") 

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]

# Using which() to look up tstrsplit() output
CMono_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "C-Monocyte")
CMono_idx

colnames(aggr_counts)[CMono_idx]
aggr_counts[1:10, CMono_idx]

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

head(colData(sce_blood))

# Extract sample-level variables
metadata <- colData(sce_blood) %>% 
  as.data.frame() %>% 
  dplyr::select(COND,Sample_ID)
dim(metadata)
head(metadata)
# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]
dim(metadata)
head(metadata)
# Rename rows
rownames(metadata) <- metadata$Sample_ID
colnames(metadata)[colnames(metadata) == "Sample_ID"] <- "sample_id"
colnames(metadata)[colnames(metadata) == "COND"] <- "Group"
# Number of cells per sample and cluster
t <- table(colData(sce_blood)$Sample_ID,
           colData(sce_blood)$cluster_id)

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
names(metadata)
names(df)
library(Matrix)
library(dplyr)

# CD8 Naive ----


idx <- which(names(counts_ls) == "CD8-Naive")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Check contents of extracted objects
cluster_counts[1:6, 1:6]
class(cluster_counts)
head(cluster_metadata)




# Convert sparse matrix to dense format
cluster_counts_dense <- as.matrix(cluster_counts)
dim(cluster_counts_dense)
# Remove genes (rows) that have zero counts across all samples (columns)
filtered_counts <- cluster_counts_dense[rowSums(cluster_counts_dense) > 0, ]
# Check dimensions after filtering
dim(filtered_counts)

dds <- DESeqDataSetFromMatrix(filtered_counts, 
                              colData = cluster_metadata, 
                              design = ~ Group)

# Apply variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)  # Normalize counts with VST

# Remove batch effect (if any, assuming 'Batch' variable exists in the metadata)
mat <- assay(vsd)  # Get VST normalized counts
mat <- limma::removeBatchEffect(mat, vsd$Batch)  # Remove batch effect
assay(vsd) <- mat  # Update the assay with batch-corrected counts

# Convert the VST normalized data into a data frame
filtered_normalized_counts <- as.data.frame(assay(vsd))

class(filtered_normalized_counts)

# Extract group information
group1_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[1]]
group2_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[2]]

# Initialize vector to store p-values
p_values <- numeric(nrow(filtered_normalized_counts))

# Loop through each gene and perform t-test
for (i in 1:nrow(filtered_normalized_counts)) {
  gene_expr <- filtered_normalized_counts[i, ]  # Get expression values for the gene
  
  # Perform t-test between the two groups
  test_result <- t.test(gene_expr[group1_samples], gene_expr[group2_samples])
  
  # Store the p-value
  p_values[i] <- test_result$p.value
}

# Create a data frame of results
t_test_results <- data.frame(Gene = rownames(filtered_normalized_counts), P_Value = p_values)

# Filter for significant genes (p-value < 0.05)
significant_genes <- t_test_results[t_test_results$P_Value < 0.05, ]
dim(significant_genes)
# Check results
head(significant_genes)
filtered_significant_counts <- filtered_normalized_counts[significant_genes$Gene, ]
# Check the resulting dataset
dim(filtered_significant_counts)



# Perform PCA on the significant gene dataset
pca_result <- prcomp(t(filtered_significant_counts), center = TRUE, scale. = TRUE)

# Extract PCA components
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

# Add group information to pca_data
pca_data$group <- cluster_metadata$Group
pca_data$name <- colnames(filtered_significant_counts)

# Define custom colors and shapes
group_colors <- c("T1D" = "#E60000",  # Brighter Red
                  "H" = "#3D7A60")  # Brighter Green

group_shapes <- c("T1D" = 16,  # Circle
                  "H" = 17)  # Triangle

# Calculate percentage of variance explained for the first two PCs
percentVar <- 100 * (summary(pca_result)$importance[2, 1:2])

# PCA Plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 5, alpha = 0.7) +  # Points
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group), show.legend = FALSE, level = 0.68) +  # Confidence Ellipses
  scale_color_manual(values = group_colors) +  # Group Colors
  scale_fill_manual(values = group_colors) +  # Ellipse Fill
  scale_shape_manual(values = group_shapes) +  # Group Shapes
  labs(title = "PCA -CD8 Naive",
       x = paste0("PC1 (", round(percentVar[1], 1), "% Variance)"),
       y = paste0("PC2 (", round(percentVar[2], 1), "% Variance)")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),  # Remove Axis Numbers
    axis.ticks = element_blank(),  # Remove Ticks
    axis.line = element_line(color = "black", linewidth = 1),  # Black Axis Lines
    panel.grid = element_blank()  # Remove Grid Lines
  )

write.csv(filtered_significant_counts,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/CD8Naive_filtered_counts.csv" )


head(cluster_metadata)
# Create a new column 'Sample' as row names and drop it from the data
cluster_metadata$Sample <- rownames(cluster_metadata)
cluster_metadata <- cluster_metadata[, c("Sample", "Group")]
# Write the data frame to a CSV file
write.csv(cluster_metadata, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/CD8Naive_metadata_samples.csv", row.names = FALSE)

# C-Monocyte ----


idx <- which(names(counts_ls) == "C-Monocyte")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Check contents of extracted objects
cluster_counts[1:6, 1:6]
class(cluster_counts)
head(cluster_metadata)




# Convert sparse matrix to dense format
cluster_counts_dense <- as.matrix(cluster_counts)
dim(cluster_counts_dense)
# Remove genes (rows) that have zero counts across all samples (columns)
filtered_counts <- cluster_counts_dense[rowSums(cluster_counts_dense) > 0, ]
# Check dimensions after filtering
dim(filtered_counts)

dds <- DESeqDataSetFromMatrix(filtered_counts, 
                              colData = cluster_metadata, 
                              design = ~ Group)

# Apply variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)  # Normalize counts with VST

# Remove batch effect (if any, assuming 'Batch' variable exists in the metadata)
mat <- assay(vsd)  # Get VST normalized counts
mat <- limma::removeBatchEffect(mat, vsd$Batch)  # Remove batch effect
assay(vsd) <- mat  # Update the assay with batch-corrected counts

# Convert the VST normalized data into a data frame
filtered_normalized_counts <- as.data.frame(assay(vsd))

class(filtered_normalized_counts)

# Extract group information
group1_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[1]]
group2_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[2]]

# Initialize vector to store p-values
p_values <- numeric(nrow(filtered_normalized_counts))

# Loop through each gene and perform t-test
for (i in 1:nrow(filtered_normalized_counts)) {
  gene_expr <- filtered_normalized_counts[i, ]  # Get expression values for the gene
  
  # Perform t-test between the two groups
  test_result <- t.test(gene_expr[group1_samples], gene_expr[group2_samples])
  
  # Store the p-value
  p_values[i] <- test_result$p.value
}

# Create a data frame of results
t_test_results <- data.frame(Gene = rownames(filtered_normalized_counts), P_Value = p_values)

# Filter for significant genes (p-value < 0.05)
significant_genes <- t_test_results[t_test_results$P_Value < 0.05, ]
dim(significant_genes)
# Check results
head(significant_genes)
filtered_significant_counts <- filtered_normalized_counts[significant_genes$Gene, ]
# Check the resulting dataset
dim(filtered_significant_counts)



# Perform PCA on the significant gene dataset
pca_result <- prcomp(t(filtered_significant_counts), center = TRUE, scale. = TRUE)

# Extract PCA components
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

# Add group information to pca_data
pca_data$group <- cluster_metadata$Group
pca_data$name <- colnames(filtered_significant_counts)

# Define custom colors and shapes
group_colors <- c("T1D" = "#E60000",  # Brighter Red
                  "H" = "#3D7A60")  # Brighter Green

group_shapes <- c("T1D" = 16,  # Circle
                  "H" = 17)  # Triangle

# Calculate percentage of variance explained for the first two PCs
percentVar <- 100 * (summary(pca_result)$importance[2, 1:2])

# PCA Plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 5, alpha = 0.7) +  # Points
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group), show.legend = FALSE, level = 0.68) +  # Confidence Ellipses
  scale_color_manual(values = group_colors) +  # Group Colors
  scale_fill_manual(values = group_colors) +  # Ellipse Fill
  scale_shape_manual(values = group_shapes) +  # Group Shapes
  labs(title = "PCA -Classical Monocyte",
       x = paste0("PC1 (", round(percentVar[1], 1), "% Variance)"),
       y = paste0("PC2 (", round(percentVar[2], 1), "% Variance)")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),  # Remove Axis Numbers
    axis.ticks = element_blank(),  # Remove Ticks
    axis.line = element_line(color = "black", linewidth = 1),  # Black Axis Lines
    panel.grid = element_blank()  # Remove Grid Lines
  )

write.csv(filtered_significant_counts,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/C-Monocyte_filtered_counts.csv" )


head(cluster_metadata)
# Create a new column 'Sample' as row names and drop it from the data
cluster_metadata$Sample <- rownames(cluster_metadata)
cluster_metadata <- cluster_metadata[, c("Sample", "Group")]
# Write the data frame to a CSV file
write.csv(cluster_metadata, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/C-Monocyte_metadata_samples.csv", row.names = FALSE)

# CD4-EM ----
idx <- which(names(counts_ls) == "CD4-EM")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Check contents of extracted objects
cluster_counts[1:6, 1:6]
class(cluster_counts)
head(cluster_metadata)




# Convert sparse matrix to dense format
cluster_counts_dense <- as.matrix(cluster_counts)
dim(cluster_counts_dense)
# Remove genes (rows) that have zero counts across all samples (columns)
filtered_counts <- cluster_counts_dense[rowSums(cluster_counts_dense) > 0, ]
# Check dimensions after filtering
dim(filtered_counts)

dds <- DESeqDataSetFromMatrix(filtered_counts, 
                              colData = cluster_metadata, 
                              design = ~ Group)

# Apply variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)  # Normalize counts with VST

# Remove batch effect (if any, assuming 'Batch' variable exists in the metadata)
mat <- assay(vsd)  # Get VST normalized counts
mat <- limma::removeBatchEffect(mat, vsd$Batch)  # Remove batch effect
assay(vsd) <- mat  # Update the assay with batch-corrected counts

# Convert the VST normalized data into a data frame
filtered_normalized_counts <- as.data.frame(assay(vsd))

class(filtered_normalized_counts)

# Extract group information
group1_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[1]]
group2_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[2]]

# Initialize vector to store p-values
p_values <- numeric(nrow(filtered_normalized_counts))

# Loop through each gene and perform t-test
for (i in 1:nrow(filtered_normalized_counts)) {
  gene_expr <- filtered_normalized_counts[i, ]  # Get expression values for the gene
  
  # Perform t-test between the two groups
  test_result <- t.test(gene_expr[group1_samples], gene_expr[group2_samples])
  
  # Store the p-value
  p_values[i] <- test_result$p.value
}

# Create a data frame of results
t_test_results <- data.frame(Gene = rownames(filtered_normalized_counts), P_Value = p_values)

# Filter for significant genes (p-value < 0.05)
significant_genes <- t_test_results[t_test_results$P_Value < 0.05, ]
dim(significant_genes)
# Check results
head(significant_genes)
filtered_significant_counts <- filtered_normalized_counts[significant_genes$Gene, ]
# Check the resulting dataset
dim(filtered_significant_counts)



# Perform PCA on the significant gene dataset
pca_result <- prcomp(t(filtered_significant_counts), center = TRUE, scale. = TRUE)

# Extract PCA components
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

# Add group information to pca_data
pca_data$group <- cluster_metadata$Group
pca_data$name <- colnames(filtered_significant_counts)

# Define custom colors and shapes
group_colors <- c("T1D" = "#E60000",  # Brighter Red
                  "H" = "#3D7A60")  # Brighter Green

group_shapes <- c("T1D" = 16,  # Circle
                  "H" = 17)  # Triangle

# Calculate percentage of variance explained for the first two PCs
percentVar <- 100 * (summary(pca_result)$importance[2, 1:2])

# PCA Plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 5, alpha = 0.7) +  # Points
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group), show.legend = FALSE, level = 0.68) +  # Confidence Ellipses
  scale_color_manual(values = group_colors) +  # Group Colors
  scale_fill_manual(values = group_colors) +  # Ellipse Fill
  scale_shape_manual(values = group_shapes) +  # Group Shapes
  labs(title = "PCA -CD4-Effector Memory T Cell",
       x = paste0("PC1 (", round(percentVar[1], 1), "% Variance)"),
       y = paste0("PC2 (", round(percentVar[2], 1), "% Variance)")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),  # Remove Axis Numbers
    axis.ticks = element_blank(),  # Remove Ticks
    axis.line = element_line(color = "black", linewidth = 1),  # Black Axis Lines
    panel.grid = element_blank()  # Remove Grid Lines
  )
dim(filtered_significant_counts)
write.csv(filtered_significant_counts,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/CD4-EM_filtered_counts.csv" )


head(cluster_metadata)
# Create a new column 'Sample' as row names and drop it from the data
cluster_metadata$Sample <- rownames(cluster_metadata)
cluster_metadata <- cluster_metadata[, c("Sample", "Group")]
# Write the data frame to a CSV file
write.csv(cluster_metadata, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/CD4-EM_metadata_samples.csv", row.names = FALSE)

# B-Naive ----
idx <- which(names(counts_ls) == "B-Naive")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Check contents of extracted objects
cluster_counts[1:6, 1:6]
class(cluster_counts)
head(cluster_metadata)




# Convert sparse matrix to dense format
cluster_counts_dense <- as.matrix(cluster_counts)
dim(cluster_counts_dense)
# Remove genes (rows) that have zero counts across all samples (columns)
filtered_counts <- cluster_counts_dense[rowSums(cluster_counts_dense) > 0, ]
# Check dimensions after filtering
dim(filtered_counts)

dds <- DESeqDataSetFromMatrix(filtered_counts, 
                              colData = cluster_metadata, 
                              design = ~ Group)

# Apply variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)  # Normalize counts with VST

# Remove batch effect (if any, assuming 'Batch' variable exists in the metadata)
mat <- assay(vsd)  # Get VST normalized counts
mat <- limma::removeBatchEffect(mat, vsd$Batch)  # Remove batch effect
assay(vsd) <- mat  # Update the assay with batch-corrected counts

# Convert the VST normalized data into a data frame
filtered_normalized_counts <- as.data.frame(assay(vsd))

class(filtered_normalized_counts)

# Extract group information
group1_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[1]]
group2_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[2]]

# Initialize vector to store p-values
p_values <- numeric(nrow(filtered_normalized_counts))

# Loop through each gene and perform t-test
for (i in 1:nrow(filtered_normalized_counts)) {
  gene_expr <- filtered_normalized_counts[i, ]  # Get expression values for the gene
  
  # Perform t-test between the two groups
  test_result <- t.test(gene_expr[group1_samples], gene_expr[group2_samples])
  
  # Store the p-value
  p_values[i] <- test_result$p.value
}

# Create a data frame of results
t_test_results <- data.frame(Gene = rownames(filtered_normalized_counts), P_Value = p_values)

# Filter for significant genes (p-value < 0.05)
significant_genes <- t_test_results[t_test_results$P_Value < 0.05, ]
dim(significant_genes)
# Check results
head(significant_genes)
filtered_significant_counts <- filtered_normalized_counts[significant_genes$Gene, ]
# Check the resulting dataset
dim(filtered_significant_counts)

# Perform PCA on the significant gene dataset
pca_result <- prcomp(t(filtered_significant_counts), center = TRUE, scale. = TRUE)

# Extract PCA components
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

# Add group information to pca_data
pca_data$group <- cluster_metadata$Group
pca_data$name <- colnames(filtered_significant_counts)

# Define custom colors and shapes
group_colors <- c("T1D" = "#E60000",  # Brighter Red
                  "H" = "#3D7A60")  # Brighter Green

group_shapes <- c("T1D" = 16,  # Circle
                  "H" = 17)  # Triangle

# Calculate percentage of variance explained for the first two PCs
percentVar <- 100 * (summary(pca_result)$importance[2, 1:2])

# PCA Plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 5, alpha = 0.7) +  # Points
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group), show.legend = FALSE, level = 0.68) +  # Confidence Ellipses
  scale_color_manual(values = group_colors) +  # Group Colors
  scale_fill_manual(values = group_colors) +  # Ellipse Fill
  scale_shape_manual(values = group_shapes) +  # Group Shapes
  labs(title = "PCA -B-Naive",
       x = paste0("PC1 (", round(percentVar[1], 1), "% Variance)"),
       y = paste0("PC2 (", round(percentVar[2], 1), "% Variance)")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),  # Remove Axis Numbers
    axis.ticks = element_blank(),  # Remove Ticks
    axis.line = element_line(color = "black", linewidth = 1),  # Black Axis Lines
    panel.grid = element_blank()  # Remove Grid Lines
  )
dim(filtered_significant_counts)
write.csv(filtered_significant_counts,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/B-Naive_filtered_counts.csv" )

head(cluster_metadata)
# Create a new column 'Sample' as row names and drop it from the data
cluster_metadata$Sample <- rownames(cluster_metadata)
cluster_metadata <- cluster_metadata[, c("Sample", "Group")]
# Write the data frame to a CSV file
write.csv(cluster_metadata, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/B-Naive_metadata_samples.csv", row.names = FALSE)

# NK ----
idx <- which(names(counts_ls) == "NK")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Check contents of extracted objects
cluster_counts[1:6, 1:6]
class(cluster_counts)
head(cluster_metadata)




# Convert sparse matrix to dense format
cluster_counts_dense <- as.matrix(cluster_counts)
dim(cluster_counts_dense)
# Remove genes (rows) that have zero counts across all samples (columns)
filtered_counts <- cluster_counts_dense[rowSums(cluster_counts_dense) > 0, ]
# Check dimensions after filtering
dim(filtered_counts)

dds <- DESeqDataSetFromMatrix(filtered_counts, 
                              colData = cluster_metadata, 
                              design = ~ Group)

# Apply variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)  # Normalize counts with VST

# Remove batch effect (if any, assuming 'Batch' variable exists in the metadata)
mat <- assay(vsd)  # Get VST normalized counts
mat <- limma::removeBatchEffect(mat, vsd$Batch)  # Remove batch effect
assay(vsd) <- mat  # Update the assay with batch-corrected counts

# Convert the VST normalized data into a data frame
filtered_normalized_counts <- as.data.frame(assay(vsd))

class(filtered_normalized_counts)

# Extract group information
group1_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[1]]
group2_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[2]]

# Initialize vector to store p-values
p_values <- numeric(nrow(filtered_normalized_counts))

# Loop through each gene and perform t-test
for (i in 1:nrow(filtered_normalized_counts)) {
  gene_expr <- filtered_normalized_counts[i, ]  # Get expression values for the gene
  
  # Perform t-test between the two groups
  test_result <- t.test(gene_expr[group1_samples], gene_expr[group2_samples])
  
  # Store the p-value
  p_values[i] <- test_result$p.value
}

# Create a data frame of results
t_test_results <- data.frame(Gene = rownames(filtered_normalized_counts), P_Value = p_values)

# Filter for significant genes (p-value < 0.05)
significant_genes <- t_test_results[t_test_results$P_Value < 0.05, ]
dim(significant_genes)
# Check results
head(significant_genes)
filtered_significant_counts <- filtered_normalized_counts[significant_genes$Gene, ]
# Check the resulting dataset
dim(filtered_significant_counts)

# Perform PCA on the significant gene dataset
pca_result <- prcomp(t(filtered_significant_counts), center = TRUE, scale. = TRUE)

# Extract PCA components
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

# Add group information to pca_data
pca_data$group <- cluster_metadata$Group
pca_data$name <- colnames(filtered_significant_counts)

# Define custom colors and shapes
group_colors <- c("T1D" = "#E60000",  # Brighter Red
                  "H" = "#3D7A60")  # Brighter Green

group_shapes <- c("T1D" = 16,  # Circle
                  "H" = 17)  # Triangle

# Calculate percentage of variance explained for the first two PCs
percentVar <- 100 * (summary(pca_result)$importance[2, 1:2])

# PCA Plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 5, alpha = 0.7) +  # Points
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group), show.legend = FALSE, level = 0.68) +  # Confidence Ellipses
  scale_color_manual(values = group_colors) +  # Group Colors
  scale_fill_manual(values = group_colors) +  # Ellipse Fill
  scale_shape_manual(values = group_shapes) +  # Group Shapes
  labs(title = "PCA -NK",
       x = paste0("PC1 (", round(percentVar[1], 1), "% Variance)"),
       y = paste0("PC2 (", round(percentVar[2], 1), "% Variance)")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),  # Remove Axis Numbers
    axis.ticks = element_blank(),  # Remove Ticks
    axis.line = element_line(color = "black", linewidth = 1),  # Black Axis Lines
    panel.grid = element_blank()  # Remove Grid Lines
  )
dim(filtered_significant_counts)
write.csv(filtered_significant_counts,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/NK_filtered_counts.csv" )

head(cluster_metadata)
# Create a new column 'Sample' as row names and drop it from the data
cluster_metadata$Sample <- rownames(cluster_metadata)
cluster_metadata <- cluster_metadata[, c("Sample", "Group")]
# Write the data frame to a CSV file
write.csv(cluster_metadata, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/NK_metadata_samples.csv", row.names = FALSE)

# T-reg ----
idx <- which(names(counts_ls) == "T-reg")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Check contents of extracted objects
cluster_counts[1:6, 1:6]
class(cluster_counts)
head(cluster_metadata)

# Convert sparse matrix to dense format
cluster_counts_dense <- as.matrix(cluster_counts)
dim(cluster_counts_dense)
# Remove genes (rows) that have zero counts across all samples (columns)
filtered_counts <- cluster_counts_dense[rowSums(cluster_counts_dense) > 0, ]
# Check dimensions after filtering
dim(filtered_counts)

dds <- DESeqDataSetFromMatrix(filtered_counts, 
                              colData = cluster_metadata, 
                              design = ~ Group)

# Apply variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)  # Normalize counts with VST

# Remove batch effect (if any, assuming 'Batch' variable exists in the metadata)
mat <- assay(vsd)  # Get VST normalized counts
mat <- limma::removeBatchEffect(mat, vsd$Batch)  # Remove batch effect
assay(vsd) <- mat  # Update the assay with batch-corrected counts

# Convert the VST normalized data into a data frame
filtered_normalized_counts <- as.data.frame(assay(vsd))

class(filtered_normalized_counts)

# Extract group information
group1_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[1]]
group2_samples <- colnames(filtered_normalized_counts)[cluster_metadata$Group == unique(cluster_metadata$Group)[2]]

# Initialize vector to store p-values
p_values <- numeric(nrow(filtered_normalized_counts))

# Loop through each gene and perform t-test
for (i in 1:nrow(filtered_normalized_counts)) {
  gene_expr <- filtered_normalized_counts[i, ]  # Get expression values for the gene
  
  # Perform t-test between the two groups
  test_result <- t.test(gene_expr[group1_samples], gene_expr[group2_samples])
  
  # Store the p-value
  p_values[i] <- test_result$p.value
}

# Create a data frame of results
t_test_results <- data.frame(Gene = rownames(filtered_normalized_counts), P_Value = p_values)

# Filter for significant genes (p-value < 0.05)
significant_genes <- t_test_results[t_test_results$P_Value < 0.05, ]
dim(significant_genes)
# Check results
head(significant_genes)
filtered_significant_counts <- filtered_normalized_counts[significant_genes$Gene, ]
# Check the resulting dataset
dim(filtered_significant_counts)

# Perform PCA on the significant gene dataset
pca_result <- prcomp(t(filtered_significant_counts), center = TRUE, scale. = TRUE)

# Extract PCA components
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

# Add group information to pca_data
pca_data$group <- cluster_metadata$Group
pca_data$name <- colnames(filtered_significant_counts)

# Define custom colors and shapes
group_colors <- c("T1D" = "#E60000",  # Brighter Red
                  "H" = "#3D7A60")  # Brighter Green

group_shapes <- c("T1D" = 16,  # Circle
                  "H" = 17)  # Triangle

# Calculate percentage of variance explained for the first two PCs
percentVar <- 100 * (summary(pca_result)$importance[2, 1:2])

# PCA Plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 5, alpha = 0.7) +  # Points
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group), show.legend = FALSE, level = 0.68) +  # Confidence Ellipses
  scale_color_manual(values = group_colors) +  # Group Colors
  scale_fill_manual(values = group_colors) +  # Ellipse Fill
  scale_shape_manual(values = group_shapes) +  # Group Shapes
  labs(title = "PCA -Tregs",
       x = paste0("PC1 (", round(percentVar[1], 1), "% Variance)"),
       y = paste0("PC2 (", round(percentVar[2], 1), "% Variance)")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),  # Remove Axis Numbers
    axis.ticks = element_blank(),  # Remove Ticks
    axis.line = element_line(color = "black", linewidth = 1),  # Black Axis Lines
    panel.grid = element_blank()  # Remove Grid Lines
  )
dim(filtered_significant_counts)
write.csv(filtered_significant_counts,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/Tregs_filtered_counts.csv" )

head(cluster_metadata)
# Create a new column 'Sample' as row names and drop it from the data
cluster_metadata$Sample <- rownames(cluster_metadata)
cluster_metadata <- cluster_metadata[, c("Sample", "Group")]
# Write the data frame to a CSV file
write.csv(cluster_metadata, "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Blood/Tregs_metadata_samples.csv", row.names = FALSE)

cluster_names


