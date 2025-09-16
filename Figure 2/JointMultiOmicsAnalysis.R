# Preprocessing Data----

#Metadata Importing
meta_batch <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Pathway Analysis/metadata_MetabolomicsCohort_PathwayAnalysis.csv", sep=",", header=T) # Metadata file


#Counts Data Importing
counts_batch <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Pathway Analysis/gene_expected_count.annot_MetabolomicsCohort_Week6.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch <- na.omit(counts_batch)

#Remove duplicate names
counts_batch <- counts_batch[!duplicated(counts_batch[, 1]), ]
genes <- counts_batch[, 1]
rownames(counts_batch) <- genes
counts_batch <- counts_batch[, -1]
head(counts_batch)
dim(counts_batch)
counts_data <- flexiDEG.function1(counts_batch, meta_batch, # Run Function 1
                                 convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                 batches = F, quality = T, variance = F,use_pseudobulk = F) #
head(counts_data)
dim(counts_data)
counts_data <- counts_data[!grepl("^Gm[0-9]", rownames(counts_data)), ]
dim(counts_data)

write.csv(counts_data, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JointMultiOmicsAnalysis/TranscriptomicsData_counts_Early.csv")


# Loading and Preprocessing Data  ----
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JointMultiOmicsAnalysis/")
metabolomics_data <- read.csv("MetabolomicsData_counts_Early.csv",check.names = FALSE)
transcriptomics_data <- read.csv("TranscriptomicsData_counts_Early.csv",check.names = FALSE)
metadata<- read.csv("Samples_Metadata.csv")
rownames(metabolomics_data)<-metabolomics_data$Compound
metabolomics_data$Compound<- NULL
head(metabolomics_data)

# Remove rows with duplicate external_gene_name in transcriptomics_data
transcriptomics_data <- transcriptomics_data[!duplicated(transcriptomics_data$external_gene_name), ]
rownames(transcriptomics_data) <- transcriptomics_data$external_gene_name
transcriptomics_data$`external_gene_name`<- NULL
head(transcriptomics_data)


#Procedure without scaling
X = metabolomics_data
Y = transcriptomics_data

X<-t(X)
Y<-t(Y)

#Ensure the samples are in same order
Y <- Y[rownames(X), , drop = FALSE]

# Scale X (metabolomics) and Y (transcriptomics)
X <- scale(X, center = TRUE, scale = TRUE)
Y <- scale(Y, center = TRUE, scale = TRUE)


# Run O2PLS----
BiocManager::install("ropls")
library(ropls)
library(ggplot2)
library(ggrepel)
ls("package:OmicsPLS")
library(OmicsPLS)



set.seed(123)  # For reproducibility

cv_result <- crossval_o2m(X = X, Y = Y,
                          a = 1:5,       # candidate joint components
                          ax = 0:5,      # specific components in X
                          ay = 0:5,      # specific components in Y
                          nr_folds = 5)  # 5-fold CV


# Get prediction error for best model
which_min <- which(cv_result$press == min(cv_result$press), arr.ind = TRUE)
print(which_min)  # Gives indices for best (n, nX, nY)


X_pca <- prcomp(X, center = TRUE, scale. = TRUE)
X_var <- X_pca$sdev^2
X_var <- X_var / sum(X_var)

qplot(1:length(X_var), X_var, geom = "line") +
  geom_point() +
  ggtitle("Scree Plot: PCA of Metabolomics (X)") +
  xlab("Principal Component") + ylab("Proportion of Variance Explained")

Y_pca <- prcomp(Y, center = TRUE, scale. = TRUE)
Y_var <- Y_pca$sdev^2
Y_var <- Y_var / sum(Y_var)

qplot(1:length(Y_var), Y_var, geom = "line") +
  geom_point() +
  ggtitle("Scree Plot: PCA of Transcriptomics (Y)") +
  xlab("Principal Component") + ylab("Proportion of Variance Explained")

cross_cov <- t(X) %*% Y
svd_vals <- svd(cross_cov)$d
svd_vals <- svd_vals / sum(svd_vals)

qplot(1:length(svd_vals), svd_vals, geom = "line") +
  geom_point() +
  ggtitle("Scree Plot: Cross-covariance XᵀY") +
  xlab("Component") + ylab("Proportion of Cross-Covariance Explained")


library(OmicsPLS)
# Perform O2PLS with 1 predictive and 1 orthogonal component
final_o2pls <- o2m(X, Y, n = 2, nx = 0, ny = 5)

T_joint <- final_o2pls$Tt  # joint scores from X
U_joint <- final_o2pls$U  # joint scores from Y

library(ggplot2)

df_scores <- data.frame(Sample = rownames(X),
                        T1 = T_joint[, 1],
                        T2 = T_joint[, 2],
                        Group = metadata$Group)  # replace with actual metadata column

ggplot(df_scores, aes(x = T1, y = T2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "O2PLS Joint Scores (from X)", x = "T1", y = "T2")


df_u <- data.frame(Sample = rownames(Y),
                   U1 = U_joint[, 1],
                   U2 = U_joint[, 2],
                   Group = metadata$Group)

ggplot(df_u, aes(x = U1, y = U2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "O2PLS Joint Scores (from Y)", x = "U1", y = "U2")

# Loadings for X (metabolites)
head(final_o2pls$W.)

# Loadings for Y (genes)
head(final_o2pls$C.)

# Extract joint and orthogonal components from Y (transcriptomics)
U1 <- final_o2pls$U[, 1]            # Predictive component
Uo1 <- final_o2pls$U_Xosc[, 1]      # Orthogonal component

# Create score data frame
score_df <- data.frame(
  Sample = rownames(Y),
  Predictive = U1,
  Orthogonal = Uo1,
  Group = metadata$Group  # Replace 'Group' if your column name is different
)

# Set custom colors and shapes
group_colors <- c(
  "Progressor" = "#E60000",         # Brighter Red
  "Non-Progressor" = "#3D7A60"      # Brighter Green
)

group_shapes <- c(
  "Progressor" = 16,                # Circle
  "Non-Progressor" = 17             # Triangle
)

# Make sure Group is a factor with desired order
score_df$Group <- factor(score_df$Group, levels = c("Progressor", "Non-Progressor"))


ggplot(score_df, aes(x = Predictive, y = Orthogonal, color = Group, shape = Group)) +
  # Confidence ellipse
  stat_ellipse(geom = "polygon", aes(fill = Group), alpha = 0.2, level = 0.68, show.legend = FALSE) +
  # Sample points
  geom_point(size = 5, alpha = 0.8) +
  # Custom color and shape scales
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  # Titles and axes
  labs(
    #title = "O2PLS Score Plot (Transcriptomics Y)",
    x = "Joint Component",
    y = "Orthogonal Component"
  ) +
  # Style
  theme_minimal(base_size = 20) +
  theme(
    axis.title = element_text(size = 25, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )


## Display Top Coefficient ----
gene_loadings_joint <- final_o2pls$C.[, 1]  # Joint component 1, transcriptomics (Y)
metab_loadings_joint <- final_o2pls$W.[, 1] # Joint component 1, metabolomics (X)
top_genes <- sort(abs(gene_loadings_joint), decreasing = TRUE)[1:30]
top_metabs <- sort(abs(metab_loadings_joint), decreasing = TRUE)[1:30]

names(top_genes)       # Gene symbols
names(top_metabs)      # Metabolite names

# # Extract joint loadings
# loadings_X <- loadings(o2pls_model_scaled, "Xjoint")[, 1]  # Metabolites
# loadings_Y <- loadings(o2pls_model_scaled, "Yjoint")[, 1]  # Genes
# top_n <- 30  # Number of top genes/metabolites to show from each
# 
# # Get top metabolites
# top_metab <- sort(abs(loadings_X), decreasing = TRUE)[1:top_n]
# top_metab_names <- names(top_metab)
# top_metab_vals <- loadings_X[top_metab_names]
# 
# # Get top genes
# top_genes <- sort(abs(loadings_Y), decreasing = TRUE)[1:top_n]
# top_gene_names <- names(top_genes)
# top_gene_vals <- loadings_Y[top_gene_names]
# 
# 
# feature_df <- data.frame(
#   Feature = c(top_metab_names, top_gene_names),
#   Coefficient = c(top_metab_vals, top_gene_vals),
#   Type = c(rep("Metabolite", top_n), rep("Gene", top_n))
# )
# 
# 
# library(ggplot2)
# 
# # Sort features by absolute coefficient (descending)
# feature_df <- feature_df %>%
#   dplyr::arrange(desc((Coefficient))) %>%
#   dplyr::mutate(Feature = factor(Feature, levels = Feature))
# 
# # Nature-compatible color palette
# nature_colors <- c("Gene" = "red", "Metabolite" = "navyblue")
# 
# # Final horizontal bar plot
# ggplot(feature_df, aes(x = Feature, y = Coefficient, fill = Type)) +
#   geom_col(width = 0.7) +
#   scale_fill_manual(values = nature_colors) +
#   labs(
#     title = "Top Contributing Features from O2PLS",
#     x = NULL,
#     y = "O2PLS Loading Coefficient"
#   ) +
#   theme_minimal(base_size = 18) +
#   theme(
#     # Typography
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
#     axis.text.y = element_text(size = 14),
#     
#     # Legend
#     legend.position = "top",
#     legend.title = element_blank(),
#     
#     # Grid and axis styling
#     panel.grid.major.y = element_blank(),  # Remove Y grid lines
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank(),
#     axis.line = element_line(color = "black", linewidth = 1),  # Add black axes
#     axis.ticks = element_line(color = "black")
#   )
# 
# # Correlation Analysis ----
# 
# top_features <- feature_df$Feature
# 
# # Subset from scaled data
# X_sub <- X_scaled[, intersect(top_features, colnames(X_scaled)), drop = FALSE]
# Y_sub <- Y_scaled[, intersect(top_features, colnames(Y_scaled)), drop = FALSE]
# 
# # Combine
# combined_df <- cbind(X_sub, Y_sub)
# # All-by-all Pearson correlation
# full_cor_matrix <- cor(combined_df, method = "pearson")
# library(pheatmap)
# 
# # Define Nature-style color gradient
# heat_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
# 
# # Annotate gene vs metabolite type
# feature_annotation <- data.frame(Type = feature_df$Type)
# rownames(feature_annotation) <- feature_df$Feature
# 
# # Order matrix to match feature order
# full_cor_matrix <- full_cor_matrix[feature_df$Feature, feature_df$Feature]
# 
# # Heatmap
# pheatmap(full_cor_matrix,
#          cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          color = heat_colors,
#          fontsize = 12,
#          border_color = NA,
#          annotation_row = feature_annotation,
#          annotation_col = feature_annotation,
#          main = "Correlation Between Top O2PLS Features")
# 
# 
# #----
# set.seed(12)
# #crossval_o2m_adjR2(X, Y, 1:3, 0:3, 0:3,  nr_cores = 8,nr_folds = nrow(X))
# crossval_o2m(X, Y, 1:3, 0:3, 0:3,  nr_cores = 8, nr_folds = nrow(X))
# 
# fit0 = o2m(X, Y, n=2, nx=1, ny=1)
# fit0
# summary(fit0)
# 
# #Procedure WITH scaling
# X1 = scale(X, scale=F)
# Y1 = scale(Y, scale=F)
# 
# set.seed(1221L)
# crossval_o2m_adjR2(X1, Y1, 1:3, 0:3, 0:3, nr_folds = nrow(X1))
# crossval_o2m(X1, Y1, 1:3, 0:3, 0:3, nr_folds = nrow(X1))