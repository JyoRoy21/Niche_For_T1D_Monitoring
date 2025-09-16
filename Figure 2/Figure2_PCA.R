
# Load the Libraries ----

pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, dnapath, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)


# Scaffold Metabolomics ----
# Read the data
data_path <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold/Early Timepoint/Data_normalized_scaf_early.csv"
Scaf_Early_Metab <- read.table(data_path, sep = ",", header = TRUE, check.names = FALSE)

# Extract metadata_scaf_metab (sample group labels from row 1)
metadata_scaf_metab <- Scaf_Early_Metab[1, -1]  # Exclude first column (which contains metabolite names)
metadata_scaf_metab <- gsub("Non-Progressor-Early", "Non-Progressor", metadata_scaf_metab)
metadata_scaf_metab <- gsub("Progressor-Early", "Progressor", metadata_scaf_metab)

# Extract metabolite names from first column, then remove the first row
metabolite_names <- Scaf_Early_Metab[-1, 1]  # Store metabolite names
Scaf_Early_Metab <- Scaf_Early_Metab[-1, -1] # Remove first row (metadata_scaf_metab) and first column (metabolite names)

# Convert to numeric
Scaf_Early_Metab <- as.data.frame(lapply(Scaf_Early_Metab, as.numeric))
rownames(Scaf_Early_Metab) <- metabolite_names  # Assign metabolite names as rownames

# Transpose: Samples as rows, metabolites as columns
Scaf_Early_Metab_t <- as.data.frame(t(Scaf_Early_Metab))
colnames(Scaf_Early_Metab_t) <- metabolite_names  # Restore proper column names
Scaf_Early_Metab_t$Group <- as.factor(unlist(metadata_scaf_metab))  # Assign group labels

# Perform pca_scaf_metab on the transposed data (excluding 'Group' column)
pca_scaf_metab <- prcomp(Scaf_Early_Metab_t[, -ncol(Scaf_Early_Metab_t)], scale. =TRUE)

# Get the variance explained by each principal component
explained_variance_scaf_metab <- summary(pca_scaf_metab)$importance[2, 1:2] * 100  # Variance explained by PC1 and PC2

# Add pca_scaf_metab results to the dataframe
Scaf_Early_Metab_t$pca_scaf_metab_x <- pca_scaf_metab$x[, 1]  # PC1 values
Scaf_Early_Metab_t$pca_scaf_metab_y <- pca_scaf_metab$x[, 2]  # PC2 values


# Define updated color and shape mappings
group_colors <- c("Progressor" = "#E60000",  # Brighter Red
                  "Non-Progressor" = "#3D7A60")  # Brighter Green



group_shapes <- c("Progressor" = 16,  # Circle for Progressor
                  "Non-Progressor" = 17)  # Triangle for Non-Progressor

# Create a pca_scaf_metab plot with custom colors, shapes, and confidence ellipses, without grid lines and axis labels
ggplot(Scaf_Early_Metab_t, aes(x = pca_scaf_metab_x, y = pca_scaf_metab_y, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +  # Scatter plot for each sample
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE,level = 0.68) +  # Confidence interval ellipses
  scale_color_manual(values = group_colors) +  # Set color for groups
  scale_fill_manual(values = group_colors) +  # Set fill color for ellipses
  scale_shape_manual(values = group_shapes) +  # Set shape for groups
  labs(title = "Scaffold Metabolomics",
       x = paste("PC1 (", round(explained_variance_scaf_metab[1], 1), "% Variance)", sep = ""),
       y = paste("PC2 (", round(explained_variance_scaf_metab[2], 1), "% Variance)", sep = "")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),  # Remove axis numbers
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.line = element_line(color = "black", linewidth = 1),  # Use 'linewidth' instead of 'size'
    panel.grid = element_blank()  # Remove grid lines
  )


# Blood Metabolomics ----
# Read the data
data_path2 <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Blood/Early Timepoint/Data_normalized_plasma_early.csv"
Blood_Early_Metab <- read.table(data_path2, sep = ",", header = TRUE, check.names = FALSE)

# Extract metadata_blood_metab (sample group labels from row 1)
metadata_blood_metab <- Blood_Early_Metab[1, -1]  # Exclude first column (which contains metabolite names)
metadata_blood_metab <- gsub("Non-Progressor-Early", "Non-Progressor", metadata_blood_metab)
metadata_blood_metab <- gsub("Progressor-Early", "Progressor", metadata_blood_metab)

# Extract metabolite names from first column, then remove the first row
metabolite_names <- Blood_Early_Metab[-1, 1]  # Store metabolite names
Blood_Early_Metab <- Blood_Early_Metab[-1, -1] # Remove first row (metadata_blood_metab) and first column (metabolite names)

# Convert to numeric
Blood_Early_Metab <- as.data.frame(lapply(Blood_Early_Metab, as.numeric))
rownames(Blood_Early_Metab) <- metabolite_names  # Assign metabolite names as rownames

# Transpose: Samples as rows, metabolites as columns
Blood_Early_Metab_t <- as.data.frame(t(Blood_Early_Metab))
colnames(Blood_Early_Metab_t) <- metabolite_names  # Restore proper column names
Blood_Early_Metab_t$Group <- as.factor(unlist(metadata_blood_metab))  # Assign group labels

# Perform pca_blood_metab on the transposed data (excluding 'Group' column)
pca_blood_metab <- prcomp(Blood_Early_Metab_t[, -ncol(Blood_Early_Metab_t)], scale. =TRUE)

# Get the variance explained by each principal component
explained_variance_blood_metab <- summary(pca_blood_metab)$importance[2, 1:2] * 100  # Variance explained by PC1 and PC2

# Add pca_blood_metab results to the dataframe
Blood_Early_Metab_t$pca_blood_metab_x <- pca_blood_metab$x[, 1]  # PC1 values
Blood_Early_Metab_t$pca_blood_metab_y <- pca_blood_metab$x[, 2]  # PC2 values


# Define updated color and shape mappings
group_colors <- c("Progressor" = "#E60000",  # Brighter Red
                  "Non-Progressor" = "#3D7A60")  # Brighter Green




group_shapes <- c("Progressor" = 16,  # Circle for Progressor
                  "Non-Progressor" = 17)  # Triangle for Non-Progressor

# Create a pca_blood_metab plot with custom colors, shapes, and confidence ellipses, without grid lines and axis labels
ggplot(Blood_Early_Metab_t, aes(x = pca_blood_metab_x, y = pca_blood_metab_y, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +  # Scatter plot for each sample
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE,level = 0.68) +  # Confidence interval ellipses
  scale_color_manual(values = group_colors) +  # Set color for groups
  scale_fill_manual(values = group_colors) +  # Set fill color for ellipses
  scale_shape_manual(values = group_shapes) +  # Set shape for groups
  labs(title = "Blood Metabolomics",
       x = paste("PC1 (", round(explained_variance_blood_metab[1], 1), "% Variance)", sep = ""),
       y = paste("PC2 (", round(explained_variance_blood_metab[2], 1), "% Variance)", sep = "")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),  # Remove axis numbers
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.line = element_line(color = "black", linewidth = 1),  # Use 'linewidth' instead of 'size'
    panel.grid = element_blank()  # Remove grid lines
  )

# Scaffold Transcriptomics ----
library(ggplot2)
library(dplyr)

# Load Data
data_path3 <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_ttestfiltered.csv"
Scaf_Early_Transcriptome <- read.table(data_path3, sep = ",", header = TRUE, check.names = FALSE)
Scaf_Early_Transcriptome <- as.data.frame(Scaf_Early_Transcriptome)

data_path_metadata_scaf_transcriptome <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_metadata.csv"
metadata_scaf_Transcriptome <- read.table(data_path_metadata_scaf_transcriptome, sep = ",", header = TRUE, check.names = FALSE)
metadata_scaf_Transcriptome <-metadata_scaf_Transcriptome[,-1]
# 🛑 Remove Batch 4 Samples from Metadata

# Now, remove Batch 4 samples
metadata_scaf_Transcriptome <- metadata_scaf_Transcriptome %>% 
  filter(Batch != 4)


# 🛑 Remove Batch 4 Samples from Expression Data
samples_to_keep <- metadata_scaf_Transcriptome$Samples  # Get sample names after filtering
# Ensure column names are characters for matching
samples_to_keep <- as.character(samples_to_keep)

# Subset the columns but keep the first column (gene names)
Scaf_Early_Transcriptome <- Scaf_Early_Transcriptome[, c(1, which(colnames(Scaf_Early_Transcriptome) %in% samples_to_keep))]

# Ensure Row Names are Gene Names
rownames(Scaf_Early_Transcriptome) <- Scaf_Early_Transcriptome[, 1]  # Assign gene names
Scaf_Early_Transcriptome <- Scaf_Early_Transcriptome[, -1]  # Remove gene name column

# Transpose: Samples as Rows, Genes as Columns
Scaf_Early_Transcriptome_t <- as.data.frame(t(Scaf_Early_Transcriptome))

# Merge Metadata
Scaf_Early_Transcriptome_t$Group <- metadata_scaf_Transcriptome$Group  # Add group info

# Perform PCA on Expression Data (Excluding Group Column)
pca <- prcomp(Scaf_Early_Transcriptome_t[, -ncol(Scaf_Early_Transcriptome_t)], scale. = TRUE)

# Get Variance Explained by PC1 & PC2
explained_variance <- summary(pca)$importance[2, 1:2] * 100  

# Add PCA Results to Data
Scaf_Early_Transcriptome_t$pca_x <- pca$x[, 1]  # PC1 values
Scaf_Early_Transcriptome_t$pca_y <- pca$x[, 2]  # PC2 values

#  Define Colors & Shapes
group_colors <- c("Progressor" = "#E60000",  # Brighter Red
                  "Non-Progressor" = "#3D7A60")  # Brighter Green

group_shapes <- c("Progressor" = 16,  # Circle
                  "Non-Progressor" = 17)  # Triangle

#PCA Plot
ggplot(Scaf_Early_Transcriptome_t, aes(x = pca_x, y = pca_y, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +  # Points
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE,level = 0.68) +  # Confidence Ellipses
  scale_color_manual(values = group_colors) +  # Group Colors
  scale_fill_manual(values = group_colors) +  # Ellipse Fill
  scale_shape_manual(values = group_shapes) +  # Group Shapes
  labs(title = "Scaffold Transcriptomic",
       x = paste("PC1 (", round(explained_variance[1], 1), "% Variance)", sep = ""),
       y = paste("PC2 (", round(explained_variance[2], 1), "% Variance)", sep = "")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),  # Remove Axis Numbers
    axis.ticks = element_blank(),  # Remove Ticks
    axis.line = element_line(color = "black", linewidth = 1),  # Black Axis Lines
    panel.grid = element_blank()  # Remove Grid Lines
  )


# Scaffold Transcriptomics-DESEQ Based ----

#Metadata Importing
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Pathway Analysis/metadata_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Pathway Analysis/metadata_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Pathway Analysis/metadata_MetabolomicsCohort_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_batch3 <- as.data.frame(meta_batch3)

# Merge metadata by columns (i.e., add samples from Batch 2 to Batch 1)
meta_combined <- rbind(meta_batch1, meta_batch2,meta_batch3)

# Preview the combined metadata
head(meta_combined)


#Counts Data Importing
counts_batch1 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Pathway Analysis/gene_expected_count.annot_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch1 <- na.omit(counts_batch1)

counts_batch2 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Pathway Analysis/gene_expected_count.annot_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch2 <- na.omit(counts_batch2)

counts_batch3 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Pathway Analysis/gene_expected_count.annot_MetabolomicsCohort_Week6.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch3 <- na.omit(counts_batch3)


#Remove duplicate names
counts_batch1 <- counts_batch1[!duplicated(counts_batch1[, 1]), ]
genes <- counts_batch1[, 1]
rownames(counts_batch1) <- genes
counts_batch1 <- counts_batch1[, -1]

counts_batch2 <- counts_batch2[!duplicated(counts_batch2[, 1]), ]
genes <- counts_batch2[, 1]
rownames(counts_batch2) <- genes
counts_batch2 <- counts_batch2[, -1]

counts_batch3 <- counts_batch3[!duplicated(counts_batch3[, 1]), ]
genes <- counts_batch3[, 1]
rownames(counts_batch3) <- genes
counts_batch3 <- counts_batch3[, -1]

#Combine data
# Merge counts_batch1 and counts_batch2
combined_counts <- merge(counts_batch1, counts_batch2, by = "row.names", all = TRUE)
# Merge the result with counts_batch3
combined_counts <- merge(combined_counts, counts_batch3, by.x = "Row.names", by.y = "row.names", all = TRUE)

# Set rownames back to genes
rownames(combined_counts) <- combined_counts$Row.names
combined_counts <- combined_counts[, -1]

# Preview the combined dataset
head(combined_counts)

library(DESeq2)
library(ggplot2)
library(dplyr)



# Assuming case1_f1 is your bulk expression data and meta_combined is the metadata

# Step 1: Run DESeq2 for Early and Intermediate groups
# Subset data for Early and Intermediate

early_data <- combined_counts[, meta_combined$Time == "Early"]
early_data <- na.omit(early_data)
meta_early <- meta_combined[meta_combined$Time == "Early", ]

early_data <- flexiDEG.function1(early_data, meta_early, # Run Function 1
                                 convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                 batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
early_data <- early_data[!grepl("^Gm[0-9]", rownames(early_data)), ]


# DESeq2 analysis for Early Time
dds_early <- DESeqDataSetFromMatrix(countData = early_data, colData = meta_early, design = ~ Batch+Group)
dds_early <- DESeq(dds_early)
results_early <- as.data.frame(results(dds_early, contrast = c("Group", "Progressor", "Non-Progressor")))

# Replace NA p values with 1
results_early$pvalue[is.na(results_early$pvalue)] <- 1
results_early$padj[is.na(results_early$padj)] <- 1


early_data_Batchcorrected <- combined_counts[, meta_combined$Time == "Early"]
early_data_Batchcorrected <- na.omit(early_data_Batchcorrected)
meta_early <- meta_combined[meta_combined$Time == "Early", ]

early_data_Batchcorrected <- flexiDEG.function1(early_data_Batchcorrected, meta_early, # Run Function 1
                                 convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                 batches = T, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
early_data_Batchcorrected <- early_data_Batchcorrected[!grepl("^Gm[0-9]", rownames(early_data_Batchcorrected)), ]

# 1. Subset significant DEGs from DESeq2 results

deg_genes <- rownames(results_early)[
  results_early$pvalue <= 0.05 #& abs(results_early$log2FoldChange) >= 1
]

# 2. Subset these genes from the batch-corrected data
early_data_Batchcorrected_DEGFiltered <- early_data_Batchcorrected[rownames(early_data_Batchcorrected) %in% deg_genes, ]

# Make sure gene names are rownames already
# Transpose the data (samples as rows, genes as columns)
early_DEG_t <- as.data.frame(t(early_data_Batchcorrected_DEGFiltered))

# Add Group info from metadata
early_DEG_t$Group <- meta_early$Group

###PCA----

# Run PCA (exclude Group column)
pca_early <- prcomp(early_DEG_t[, -ncol(early_DEG_t)], scale. = TRUE)

# % Variance explained by PC1 & PC2
explained_var <- summary(pca_early)$importance[2, 1:2] * 100

# Add PC scores
early_DEG_t$pca_x <- pca_early$x[, 1]
early_DEG_t$pca_y <- pca_early$x[, 2]

# Colors and shapes
group_colors <- c("Progressor" = "#E60000", "Non-Progressor" = "#3D7A60")
group_shapes <- c("Progressor" = 16, "Non-Progressor" = 17)

# Plot
ggplot(early_DEG_t, aes(x = pca_x, y = pca_y, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.68) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "PCA: Early Timepoint (DEG-Filtered)",
       x = paste0("PC1 (", round(explained_var[1], 1), "% Variance)"),
       y = paste0("PC2 (", round(explained_var[2], 1), "% Variance)")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )

### LDA ----
# Step 1: PCA
pca_result <- prcomp(lda_input[, -ncol(lda_input)], scale. = TRUE)
top_pcs <- as.data.frame(pca_result$x[, 1:20])  # Select top 20 PCs

# Step 2: Add group info
top_pcs$Group <- lda_input$Group

# Step 3: LDA on top PCs
lda_fit <- lda(Group ~ ., data = top_pcs)
lda_pred <- predict(lda_fit)
library(ggplot2)

# Make plot data
lda_df <- data.frame(LD1 = lda_pred$x[,1], Group = top_pcs$Group)

# Plot LD1 separation
ggplot(lda_df, aes(x = Group, y = LD1, fill = Group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8, aes(color = Group)) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  theme_minimal(base_size = 18) +
  labs(title = "LDA (on Top PCs): LD1 Separates Progressor vs Non-Progressor",
       y = "LD1", x = "") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

ggplot(lda_df, aes(x = LD1, fill = Group)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = group_colors) +
  theme_minimal(base_size = 18) +
  labs(title = "LDA Density Plot: LD1 Discriminates Groups", x = "LD1", y = "Density")

### PLSDA----
devtools::install_github("mixOmicsTeam/mixOmics")
library(mixOmics)
packageVersion("mixOmics")  # Should be 6.24.0 or similar


# Install if needed
# install.packages("mixOmics")

# Restart R session first (Ctrl+Shift+F10 in RStudio)
# Then unload conflicting packages if loaded
detach("package:pls", unload = TRUE)
detach("package:caret", unload = TRUE)
detach("package:plsRglm", unload = TRUE)

# Prepare data
X <- t(early_data_Batchcorrected_DEGFiltered)  # samples x genes
Y <- meta_early$Group                          # class labels

Y <- as.factor(meta_early$Group)  # Convert class labels to factor

# Now run PLS-DA
plsda_model <-  mixOmics::plsda(X, Y, ncomp = 2)
# Extract component scores

scores <- plsda_model$variates$X  # This contains LV1, LV2, ...

# Build data frame for ggplot
plot_df <- data.frame(LV1 = scores[,1],
                      LV2 = scores[,2],
                      Group = Y)

# Define colors & shapes
group_colors <- c("Progressor" = "#E60000", "Non-Progressor" = "#3D7A60")
group_shapes <- c("Progressor" = 16, "Non-Progressor" = 17)

# Plot using ggplot2
library(ggplot2)
ggplot(plot_df, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.5, aes(fill = Group), show.legend = FALSE, level = 0.70) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "PLS-DA: Early Timepoint (DEG-filtered)",
       x = "Latent Variable 1", y = "Latent Variable 2") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )

set.seed(123)  # For reproducibility
tune_results <- tune.plsda(
  X, Y,
  ncomp = 2,
  validation = "Mfold",
  folds = 5,
  dist = "max.dist",  # or "centroids.dist"
  progressBar = TRUE,
  measure = "BER",  # Balanced Error Rate
  test.keepX = seq(5, 100, by = 5)  # Try 5 to 100 genes
)

# Extract VIP scores for all components
vip_scores <- vip(plsda_model)
# Extract VIPs from component 1
vip_lv1 <- vip_scores[, 1]

# Sort descending and take top 20
top_vip_genes <- sort(vip_lv1, decreasing = TRUE)[1:20]

# Build data frame for plotting
vip_df <- data.frame(
  Gene = names(top_vip_genes),
  VIP = as.numeric(top_vip_genes)
)

# Plot
ggplot(vip_df, aes(x = reorder(Gene, VIP), y = VIP)) +
  geom_bar(stat = "identity", fill = "#2C3E50") +
  coord_flip() +
  labs(title = "Top 20 Genes by VIP Score (LV1)",
       x = "Gene", y = "VIP Score") +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )

top_gene_names <- names(sort(vip_scores[,1], decreasing = TRUE)[1:50])
# Subset X to top VIP genes only
X_top <- X[, top_gene_names]

# Rerun PLS-DA
plsda_top_model <- plsda(X_top, Y, ncomp = 2)

# Extract scores
scores_top <- plsda_top_model$variates$X

# Plot
plot_df_top <- data.frame(LV1 = scores_top[,1],
                          LV2 = scores_top[,2],
                          Group = Y)

ggplot(plot_df_top, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.68) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "PLS-DA: Top 20 VIP Genes Only",
       x = "Latent Variable 1", y = "Latent Variable 2") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )

#### Input Dataset for Prediction----
#Metadata Importing
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_MetabolomicsCohort_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch4 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_ScRNASeqCohort.csv", sep=",", header=T) # Metadata file

meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_batch3 <- as.data.frame(meta_batch3)
meta_batch4 <- as.data.frame(meta_batch4)

# Merge metadata by columns (i.e., add samples from Batch 2 to Batch 1)
meta_combined <- rbind(meta_batch1, meta_batch2,meta_batch3,meta_batch4)

# Preview the combined metadata
head(meta_combined)

counts_batch1 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/gene_expected_count.annot_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch1 <- na.omit(counts_batch1)

counts_batch2 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/gene_expected_count.annot_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch2 <- na.omit(counts_batch2)

counts_batch3 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/gene_expected_count.annot_MetabolomicsCohort_Week6.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch3 <- na.omit(counts_batch3)

counts_batch4 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/gene_expected_count.annot_ScRNASeqCohort.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch4 <- na.omit(counts_batch4)

#Remove duplicate names
counts_batch1 <- counts_batch1[!duplicated(counts_batch1[, 1]), ]
genes <- counts_batch1[, 1]
rownames(counts_batch1) <- genes
counts_batch1 <- counts_batch1[, -1]

counts_batch2 <- counts_batch2[!duplicated(counts_batch2[, 1]), ]
genes <- counts_batch2[, 1]
rownames(counts_batch2) <- genes
counts_batch2 <- counts_batch2[, -1]

counts_batch3 <- counts_batch3[!duplicated(counts_batch3[, 1]), ]
genes <- counts_batch3[, 1]
rownames(counts_batch3) <- genes
counts_batch3 <- counts_batch3[, -1]


counts_batch4 <- counts_batch4[!duplicated(counts_batch4[, 1]), ]
genes <- counts_batch4[, 1]
rownames(counts_batch4) <- genes
counts_batch4 <- counts_batch4[, -1]

#Combine data
# Merge counts_batch1 and counts_batch2
combined_counts <- merge(counts_batch1, counts_batch2, by = "row.names", all = TRUE)
# Merge the result with counts_batch3
combined_counts <- merge(combined_counts, counts_batch3, by.x = "Row.names", by.y = "row.names", all = TRUE)
# Merge the result with counts_batch4
combined_counts <- merge(combined_counts, counts_batch4, by.x = "Row.names", by.y = "row.names", all = TRUE)



# Set rownames back to genes
rownames(combined_counts) <- combined_counts$Row.names
combined_counts <- combined_counts[, -1]

# Preview the combined dataset
head(combined_counts)




# Normalizing Genes Input to Machine Learning Model 
combined_counts <- combined_counts[, meta_combined$Samples]  # Ensure Sample_IDs match column names in 

# Filter Data for time
meta_combined_early <- meta_combined[meta_combined$Time == "Early", ]

# Filter combined_counts to keep only the samples in the subset metadata
combined_counts_early <- combined_counts[, colnames(combined_counts) %in% meta_combined_early$Samples]

combined_counts_early <- combined_counts_early[, meta_combined_early$Samples]  # Ensure Sample_IDs match column names in 
# Remove rows with NA values
combined_counts_early <- combined_counts_early[complete.cases(combined_counts_early), ]



case1_f1 <- flexiDEG.function1(combined_counts_early, meta_combined_early, # Run Function 1
                               convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                               batches = T, quality = T, variance = F,use_pseudobulk = F) # Select filters: 2, 0, 15

case1_f1 <- case1_f1[!grepl("^Gm[0-9]", rownames(case1_f1)), ]


case1_f1_allsignificant <- case1_f1[rownames(case1_f1) %in% deg_genes, ]
write.csv(case1_f1_allsignificant  , "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered.csv", row.names = TRUE)




top_25 <- names(sort(vip_scores[,1], decreasing = TRUE)[1:25])
case1_f1_top25 <- case1_f1[rownames(case1_f1) %in% top_25, ]
write.csv(case1_f1_top25  , "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered_PLSDATop25.csv", row.names = TRUE)


top_50 <- names(sort(vip_scores[,1], decreasing = TRUE)[1:50])
case1_f1_top50 <- case1_f1[rownames(case1_f1) %in% top_50, ]
write.csv(case1_f1_top50  , "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered_PLSDATop50.csv", row.names = TRUE)

top_100 <- names(sort(vip_scores[,1], decreasing = TRUE)[1:100])
case1_f1_top100 <- case1_f1[rownames(case1_f1) %in% top_100, ]
write.csv(case1_f1_top100  , "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered_PLSDATop100.csv", row.names = TRUE)


top_200 <- names(sort(vip_scores[,1], decreasing = TRUE)[1:200])
case1_f1_top200 <- case1_f1[rownames(case1_f1) %in% top_200, ]
write.csv(case1_f1_top200  , "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered_PLSDATop200.csv", row.names = TRUE)

top_500 <- names(sort(vip_scores[,1], decreasing = TRUE)[1:500])
case1_f1_top500 <- case1_f1[rownames(case1_f1) %in% top_500, ]
write.csv(case1_f1_top500  , "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered_PLSDATop500.csv", row.names = TRUE)

top_750 <- names(sort(vip_scores[,1], decreasing = TRUE)[1:750])
case1_f1_top750 <- case1_f1[rownames(case1_f1) %in% top_750, ]
write.csv(case1_f1_top750  , "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered_PLSDATop750.csv", row.names = TRUE)
