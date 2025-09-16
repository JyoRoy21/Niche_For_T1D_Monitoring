
pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, dnapath, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)

library(dplyr)
# MSIGDBR Pathways ----
BiocManager::install("msigdb")
library(msigdbr)
library(ExperimentHub)
library(GSEABase)
# Needs msigdbr package: https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
msigdbr_collections() # Take a look at all the pathway groups in the msigdbr database

sets_hallmark <- msigdbr(species="Mus musculus", collection="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names
sets_reactome <- msigdbr(species="Mus musculus", subcollection="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names
kegg_gene_sets <- msigdbr(species="Mus musculus", subcollection="CP:KEGG_LEGACY") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names
biocarta_gene_sets <- msigdbr(species="Mus musculus", subcollection="CP:BIOCARTA") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names
pwl_msigdbr <- c(pwl_hallmark, pwl_reactome, pwl_kegg) # Compile them all
length(pwl_msigdbr)

getwd()

# Type 1 Diabetes (Bulk RNA)----

## Import Data----

#Metadata
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_MetabolomicsCohort_PathwayAnalysis.csv", sep=",", header=T) # Metadata file

meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_batch3 <- as.data.frame(meta_batch3)

#Counts Data

counts_batch1 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/gene_expected_count.annot_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch1 <- na.omit(counts_batch1)

counts_batch2 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/gene_expected_count.annot_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch2 <- na.omit(counts_batch2)

counts_batch3 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/gene_expected_count.annot_MetabolomicsCohort_Week6.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
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


# Merge metadata by columns (i.e., add samples from Batch 2 to Batch 1)
meta_combined <- rbind(meta_batch1, meta_batch2,meta_batch3)
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

## Early Timepoint----
early_data <- combined_counts[, meta_combined$Time == "Early"]
early_data <- na.omit(early_data)
meta_early <- meta_combined[meta_combined$Time == "Early", ]

early_data <- flexiDEG.function1(early_data, meta_early, # Run Function 1
                                 convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                 batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
early_data <- early_data[!grepl("^Gm[0-9]", rownames(early_data)), ]


# DESEQ2
dds_early <- DESeqDataSetFromMatrix(countData = early_data, colData = meta_early, design = ~ Batch+Group)
dds_early <- DESeq(dds_early)
results_early <- as.data.frame(results(dds_early, contrast = c("Group", "Progressor", "Non-Progressor")))
# Replace all NA values with 1 in results_early
results_early[is.na(results_early)] <- 1
# Filter for significant genes with log2FC >= 1 or <= -1 and p-value <= 0.05
T1D_DEGs_Early <- results_early[
  (results_early$log2FoldChange >= 1 | results_early$log2FoldChange <= -1) &
    results_early$pvalue <= 0.05,
]

# Cancer 4T1 (Bulk RNA)----

## Import Data----

#Metadata
meta_cancer <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/PanDiseaseInflammation/Cancer4T1_Metadata.csv", sep=",", header=T) # Metadata file
meta_cancer <- as.data.frame(meta_cancer)


#Counts Data
counts_cancer <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/PanDiseaseInflammation/Cancer4T1_RNASeq.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_cancer <- na.omit(counts_cancer)
counts_cancer <- counts_cancer[!duplicated(counts_cancer[, 1]), ]
genes <- counts_cancer[, 1]
rownames(counts_cancer) <- genes
counts_cancer <- counts_cancer[, -1]

counts_cancer <- flexiDEG.function1(counts_cancer, meta_cancer, # Run Function 1
                                 convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                 batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
counts_cancer <- counts_cancer[!grepl("^Gm[0-9]", rownames(counts_cancer)), ]

# DESEQ2
dds_cancer <- DESeqDataSetFromMatrix(countData = counts_cancer, colData = meta_cancer, design = ~ Group)
dds_cancer <- DESeq(dds_cancer)
results_cancer <- as.data.frame(results(dds_cancer, contrast = c("Group", "Diseased", "Healthy")))
# Replace all NA values with 1 in results_early
results_cancer[is.na(results_cancer)] <- 1
# Filter for significant genes with log2FC >= 1 or <= -1 and p-value <= 0.05
Cancer4T1_DEGs <- results_cancer[
  (results_cancer$log2FoldChange >= 1 | results_cancer$log2FoldChange <= -1) &
    results_cancer$pvalue <= 0.05,
]

# EAE (Bulk RNA)----
meta_EAE <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/PanDiseaseInflammation/EAE_RNASeq_Metadata.csv", sep=",", header=T) # Metadata file
meta_EAE <- as.data.frame(meta_EAE)


counts_EAE <- read.csv("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/PanDiseaseInflammation/EAE_RNASeq.csv", 
                           header = TRUE, check.names = FALSE)
counts_EAE <- na.omit(counts_EAE)

counts_EAE <- counts_EAE[!duplicated(counts_EAE[, 1]), ]
genes <- counts_EAE[, 1]
rownames(counts_EAE) <- genes
counts_EAE <- counts_EAE[, -1]

counts_EAE <- flexiDEG.function1(counts_EAE, meta_EAE, # Run Function 1
                                    convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                    batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
counts_EAE <- counts_EAE[!grepl("^Gm[0-9]", rownames(counts_EAE)), ]

# DESEQ2
dds_EAE <- DESeqDataSetFromMatrix(countData = counts_EAE, colData = meta_EAE, design = ~ Group)
dds_EAE <- DESeq(dds_EAE)
results_EAE <- as.data.frame(results(dds_EAE, contrast = c("Group", "Diseased", "Healthy")))
# Replace all NA values with 1 in results_early
results_EAE[is.na(results_EAE)] <- 1
# Filter for significant genes with log2FC >= 1 or <= -1 and p-value <= 0.05
EAE_DEGs <- results_EAE[
  (results_EAE$log2FoldChange >= 1 | results_EAE$log2FoldChange <= -1) &
    results_EAE$pvalue <= 0.05,
]


# # Remove rows with NA 
# counts_EAE_raw <- na.omit(counts_EAE_raw)
# 
# # Step 1: Keep only numeric columns + Gene
# counts_EAE_clean <- counts_EAE_raw[, c("Gene", setdiff(colnames(counts_EAE_raw), c("Gene", "Target")))]
# 
# # Step 2: Collapse duplicate genes using base R aggregate
# counts_EAE <- aggregate(. ~ Gene, data = counts_EAE_clean, FUN = mean)
# 
# # Step 3: Set gene names as rownames
# rownames(counts_EAE) <- counts_EAE$Gene
# counts_EAE <- counts_EAE[, -1]  # drop Gene column
# 
# 
# library(limma)
# 
# # Step 1: Extract housekeeping genes from the original file
# # Re-read to get Target info if needed
# target_info <- read.csv("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/PanDiseaseInflammation/EAE_OpenArray.csv", 
#                         header = TRUE, check.names = FALSE)[, c("Gene", "Target")]
# 
# # Get unique housekeeping gene names
# housekeeping_genes <- unique(target_info$Gene[target_info$Target == "Endogenous Control"])
# 
# # Step 2: Compute ΔCt using mean of housekeeping genes per sample
# hk_means <- colMeans(counts_EAE[rownames(counts_EAE) %in% housekeeping_genes, , drop = FALSE], na.rm = TRUE)
# dct <- sweep(counts_EAE, 2, hk_means, FUN = "-")  # ΔCt = Ct(gene) - mean Ct(HK)
# 
# # Step 3: Convert to -ΔCt so higher = higher expression
# expr_matrix <- -1 * dct
# 
# # Step 4: Match sample order in metadata to columns in expr_matrix
# meta_EAE <- meta_EAE[match(colnames(expr_matrix), meta_EAE$Sample), ]
# 
# # Step 5: Design matrix for limma
# design <- model.matrix(~ 0 + meta_EAE$Group)
# colnames(design) <- levels(factor(meta_EAE$Group))  # Ensure correct names
# 
# # Step 6: Run limma differential expression analysis
# fit <- lmFit(expr_matrix, design)
# contrast_matrix <- makeContrasts(Diseased_vs_Healthy = Diseased - Healthy, levels = design)
# fit2 <- contrasts.fit(fit, contrast_matrix)
# fit2 <- eBayes(fit2)
# 
# # Step 7: Extract DEGs and filter significant ones
# DEGs_EAE <- topTable(fit2, adjust.method = "fdr", number = Inf)
# EAE_DEGs <- DEGs_EAE %>%
#   filter(abs(logFC) >= 1, P.Value <= 0.05)
# 
# # Output head
# head(EAE_DEGs)


# Comparison Across Multiple Disease ----
## DESEQ ----
### Without Directionality ----

# Assuming rownames are gene names
genes_EAE <- rownames(EAE_DEGs)
genes_Cancer <- rownames(Cancer4T1_DEGs)
genes_T1D <- rownames(T1D_DEGs_Early)
#install.packages("VennDiagram")  # Only once





library(VennDiagram)
library(grid)

# Custom colors (calmer, more refined palette)
custom_colors <- c("#3E8EDE", "#D94E5D", "#56B870")  # Blue, red, green

# Create high-quality Venn plot
venn.plot <- venn.diagram(
  x = list(
    `EAE` = genes_EAE,
    `Cancer` = genes_Cancer,
    `Type 1 Diabetes` = genes_T1D
  ),
  filename = NULL,
  fill = custom_colors,
  alpha = 0.7,
  lty = "blank", # No borders
  cex = 2.5, # Inner number font size
  fontface = "bold",
  fontfamily = "Helvetica",
  cat.fontface = "bold",
  cat.fontfamily = "Helvetica",
  cat.cex = 2.2,
  cat.col = custom_colors,
  margin = 0.08,
  main = "Shared and Unique DEGs Across Disease Models",
  main.cex = 2.5,
  main.fontface = "bold"
)

# Draw it in viewer or export to file
grid.newpage()
grid.draw(venn.plot)


# Pairwise intersections
EAE_Cancer <- intersect(genes_EAE, genes_Cancer)
EAE_T1D <- intersect(genes_EAE, genes_T1D)
Cancer_T1D <- intersect(genes_Cancer, genes_T1D)

# Optionally: number of overlaps
length(EAE_Cancer)    # EAE ∩ Cancer
length(EAE_T1D)       # EAE ∩ T1D
length(Cancer_T1D)    # Cancer ∩ T1D

# Show actual gene names
EAE_Cancer
EAE_T1D
Cancer_T1D


# library(UpSetR)
# upset_data <- list(
#   EAE = genes_EAE,
#   Cancer = genes_Cancer,
#   T1D = genes_T1D
# )
# upset(fromList(upset_data), order.by = "freq", sets = c("Cancer", "T1D", "EAE"),
#       sets.bar.color = "#444444", main.bar.color = "#2E86C1",
#       text.scale = 1.8, point.size = 4, line.size = 1.5)


### With Directionality ----

# Step 1: Extract logFCs and signs
signs_EAE <- sign(EAE_DEGs$log2FoldChange)
names(signs_EAE) <- rownames(EAE_DEGs)

signs_Cancer <- sign(Cancer4T1_DEGs$log2FoldChange)
names(signs_Cancer) <- rownames(Cancer4T1_DEGs)

signs_T1D <- sign(T1D_DEGs_Early$log2FoldChange)
names(signs_T1D) <- rownames(T1D_DEGs_Early)

# Step 2: Define a helper function to get shared genes with matching sign
get_overlap_same_sign <- function(a_signs, b_signs) {
  common_genes <- intersect(names(a_signs), names(b_signs))
  same_sign_genes <- common_genes[a_signs[common_genes] == b_signs[common_genes]]
  return(same_sign_genes)
}

# Step 3: Pairwise intersections with same sign
Cancer_T1D_same_sign <- get_overlap_same_sign(signs_Cancer, signs_T1D)
EAE_T1D_same_sign <- get_overlap_same_sign(signs_EAE, signs_T1D)
EAE_Cancer_same_sign <- get_overlap_same_sign(signs_EAE, signs_Cancer)

# Step 4: Optional — intersection across all 3 with same sign
triple_common <- Reduce(intersect, list(Cancer_T1D_same_sign, EAE_T1D_same_sign, EAE_Cancer_same_sign))

# Output
list(
  Cancer_T1D_same_sign = Cancer_T1D_same_sign,
  EAE_T1D_same_sign = EAE_T1D_same_sign,
  EAE_Cancer_same_sign = EAE_Cancer_same_sign,
  All_Three_same_sign = triple_common
)

# Extract gene names by up/down for each dataset
up_EAE <- rownames(EAE_DEGs[EAE_DEGs$log2FoldChange > 0, ])
down_EAE <- rownames(EAE_DEGs[EAE_DEGs$log2FoldChange < 0, ])

up_Cancer <- rownames(Cancer4T1_DEGs[Cancer4T1_DEGs$log2FoldChange > 0, ])
down_Cancer <- rownames(Cancer4T1_DEGs[Cancer4T1_DEGs$log2FoldChange < 0, ])

up_T1D <- rownames(T1D_DEGs_Early[T1D_DEGs_Early$log2FoldChange > 0, ])
down_T1D <- rownames(T1D_DEGs_Early[T1D_DEGs_Early$log2FoldChange < 0, ])

# Upregulated genes in common
common_up <- list(
  Cancer4T1 = up_Cancer,
  T1D = up_T1D,
  EAE = up_EAE
)

# Downregulated genes in common
common_down <- list(
  Cancer4T1 = down_Cancer,
  T1D = down_T1D,
  EAE = down_EAE
)


## SVC Gene Marker ----
### Without Directionality ----

# Assuming rownames are gene names

data_path3 <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_ttestfiltered.csv"
Scaf_Early_Transcriptome <- read.table(data_path3, sep = ",", header = TRUE, check.names = FALSE)
Scaf_Early_Transcriptome <- as.data.frame(Scaf_Early_Transcriptome)

rownames(Scaf_Early_Transcriptome)<-Scaf_Early_Transcriptome[,1]

genes_EAE <- rownames(EAE_DEGs)
genes_Cancer <- rownames(Cancer4T1_DEGs)
genes_T1D <- rownames(Scaf_Early_Transcriptome)
#install.packages("VennDiagram")  # Only once





library(VennDiagram)
library(grid)

# Custom colors (calmer, more refined palette)
custom_colors <- c("#3E8EDE", "#D94E5D", "#56B870")  # Blue, red, green

# Create high-quality Venn plot
venn.plot <- venn.diagram(
  x = list(
    `EAE` = genes_EAE,
    `Cancer` = genes_Cancer,
    `Type 1 Diabetes` = genes_T1D
  ),
  filename = NULL,
  fill = custom_colors,
  alpha = 0.7,
  lty = "blank", # No borders
  cex = 2.5, # Inner number font size
  fontface = "bold",
  fontfamily = "Helvetica",
  cat.fontface = "bold",
  cat.fontfamily = "Helvetica",
  cat.cex = 2.2,
  cat.col = custom_colors,
  margin = 0.08,
  main = "Shared and Unique DEGs Across Disease Models",
  main.cex = 2.5,
  main.fontface = "bold"
)

# Draw it in viewer or export to file
grid.newpage()
grid.draw(venn.plot)


# Pairwise intersections
EAE_Cancer <- intersect(genes_EAE, genes_Cancer)
EAE_T1D <- intersect(genes_EAE, genes_T1D)
Cancer_T1D <- intersect(genes_Cancer, genes_T1D)

# Optionally: number of overlaps
length(EAE_Cancer)    # EAE ∩ Cancer
length(EAE_T1D)       # EAE ∩ T1D
length(Cancer_T1D)    # Cancer ∩ T1D

# Show actual gene names
EAE_Cancer
EAE_T1D
Cancer_T1D


# library(UpSetR)
# upset_data <- list(
#   EAE = genes_EAE,
#   Cancer = genes_Cancer,
#   T1D = genes_T1D
# )
# upset(fromList(upset_data), order.by = "freq", sets = c("Cancer", "T1D", "EAE"),
#       sets.bar.color = "#444444", main.bar.color = "#2E86C1",
#       text.scale = 1.8, point.size = 4, line.size = 1.5)


## GSEA ----
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
# KEGG, REACTOME, HALLMARK
mm_metabolism_df <- data.frame(
  gs_name = rep(names(pwl_msigdbr), sapply(pwl_msigdbr, length)),  # Pathway names
  gene_symbol = unlist(pwl_msigdbr)  # Flatten the list into a single vector
)


### T1D----
results_early$gene <- rownames(results_early)

# Prepare the ranked list for Islet
results_early <- results_early[, c("gene", "log2FoldChange", "padj")]
results_early <- results_early[!is.na(results_early$log2FoldChange), ]
results_early <- results_early[order(results_early$log2FoldChange, decreasing = TRUE), ]
lfc_vector_early <- setNames(results_early$log2FoldChange, results_early$gene)





T1D_GSEA_early <- GSEA(
  geneList = lfc_vector_early, # Your ordered ranked gene list for Islet
  minGSSize = 5, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set size
  pvalueCutoff = 1, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed for reproducibility
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = mm_metabolism_df  # Use the new data frame
)

# Extract results for Islet
T1D_GSEA_early <- as.data.frame(T1D_GSEA_early)

### Cancer----
results_cancer$gene <- rownames(results_cancer)

# Prepare the ranked list for Islet
results_cancer <- results_cancer[, c("gene", "log2FoldChange", "padj")]
results_cancer <- results_cancer[!is.na(results_cancer$log2FoldChange), ]
results_cancer <- results_cancer[order(results_cancer$log2FoldChange, decreasing = TRUE), ]
lfc_vector_cancer <- setNames(results_cancer$log2FoldChange, results_cancer$gene)



Cancer_GSEA <- GSEA(
  geneList = lfc_vector_cancer, # Your ordered ranked gene list for Islet
  minGSSize = 5, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set size
  pvalueCutoff = 1, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed for reproducibility
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = mm_metabolism_df  # Use the new data frame
)

# Extract results for Islet
Cancer_GSEA <- as.data.frame(Cancer_GSEA)

### EAE----
results_EAE$gene <- rownames(results_EAE)

# Prepare the ranked list for EAE
results_EAE <- results_EAE[, c("gene", "log2FoldChange", "padj")]
results_EAE <- results_EAE[!is.na(results_EAE$log2FoldChange), ]
results_EAE <- results_EAE[order(results_EAE$log2FoldChange, decreasing = TRUE), ]
lfc_vector_EAE <- setNames(results_EAE$log2FoldChange, results_EAE$gene)



EAE_GSEA <- GSEA(
  geneList = lfc_vector_EAE, # Your ordered ranked gene list for Islet
  minGSSize = 5, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set size
  pvalueCutoff = 1, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed for reproducibility
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = mm_metabolism_df  # Use the new data frame
)

# Extract results for Islet
EAE_GSEA <- as.data.frame(EAE_GSEA)

### Combined ----
library(dplyr)
library(ggplot2)
library(forcats)

# Add disease labels
T1D_GSEA_early$Disease <- "T1D"
Cancer_GSEA$Disease <- "Cancer"
EAE_GSEA$Disease <- "EAE"

# Select relevant columns
T1D_sel <- T1D_GSEA_early %>% select(ID, NES, p.adjust, Disease)
Cancer_sel <- Cancer_GSEA %>% select(ID, NES, p.adjust, Disease)
EAE_sel <- EAE_GSEA %>% select(ID, NES, p.adjust, Disease)

# Combine into one dataframe
gsea_combined <- bind_rows(T1D_sel, Cancer_sel, EAE_sel)

# # Optional: filter top pathways if needed (e.g., top 15 by p-adjusted value)
# top_pathways <- gsea_combined %>%
#   group_by(Disease) %>%
#   top_n(-10, p.adjust) %>%
#   pull(ID) %>% unique()
# gsea_combined <- gsea_combined %>% filter(ID %in% top_pathways)

filtered_pathways <- c(
  #Cancer
  "KEGG_RIBOSOME",
  "REACTOME_TRANSLATION",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
  "REACTOME_RRNA_PROCESSING",
  "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
  "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_NEUTROPHIL_DEGRANULATION",
  "REACTOME_TNFS_BIND_THEIR_PHYSIOLOGICAL_RECEPTORS",
  "REACTOME_FORMATION_OF_THE_BETA_CATENIN_TCF_TRANSACTIVATING_COMPLEX",
  "REACTOME_PRE_NOTCH_EXPRESSION_AND_PROCESSING",
  
  #T1D
  "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING",
  "REACTOME_SIGNALING_BY_INTERLEUKINS",
  "REACTOME_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  "REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS",
  "REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS",
  "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
  "REACTOME_FCGR_ACTIVATION",
  "REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS",
  "REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "KEGG_INSULIN_SIGNALING_PATHWAY",
  "KEGG_PPAR_SIGNALING_PATHWAY",
  "HALLMARK_ADIPOGENESIS",
  "REACTOME_METABOLISM_OF_CARBOHYDRATES",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY",
  "REACTOME_PYRUVATE_METABOLISM",
  "REACTOME_ION_HOMEOSTASIS",
  "HALLMARK_PEROXISOME",
  "REACTOME_DISEASES_OF_METABOLISM",
  #EAE
  "REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE",
  "REACTOME_KERATINIZATION",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "REACTOME_STRIATED_MUSCLE_CONTRACTION"
  # "REACTOME_MITOCHONDRIAL_TRANSLATION",
  # "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",
  # "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",
  # "REACTOME_SELENOAMINO_ACID_METABOLISM",
  # "REACTOME_INFLUENZA_INFECTION"
)

# Filter the combined GSEA dataframe
gsea_combined <- gsea_combined %>%
  filter(ID %in% filtered_pathways)



# Format for plotting
gsea_combined <- gsea_combined %>%
  mutate(
    neglog10padj = -log10(p.adjust + 1e-10),
    ID = factor(ID, levels = rev(filtered_pathways))  # use rev() if you want top-to-bottom like list
  )

# Plot
ggplot(gsea_combined, aes(x = Disease, y = ID)) +
  geom_point(aes(size = neglog10padj, color = NES)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(2, 10)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Enrichment Across Diseases",
    x = "Disease Model",
    y = "Pathway",
    color = "NES",
    size = "-log10(FDR)"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14, face = "bold", angle = 30, hjust = 1),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )





