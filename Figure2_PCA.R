
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



