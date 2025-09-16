pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, dnapath, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)
library(dplyr)
## MSIGDBR Pathways ----
# Needs msigdbr package: https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
msigdbr_collections() # Take a look at all the pathway groups in the msigdbr database
sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names
sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names
kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG_LEGACY") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names
biocarta_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:BIOCARTA") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names
pwl_msigdbr <- c(pwl_hallmark, pwl_reactome, pwl_kegg, pwl_biocarta) # Compile them all
length(pwl_msigdbr)


getwd()

## NOD Early Stage- Progressor Vs Non-Progressor----

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(dplyr)


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

# Load necessary libraries
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

## Adoptive Transfer - Day 10- BDC2.5 vs OVA----
#Metadata Importing
meta_AT <- read.table("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data/metadata_AdoptiveTransfer.csv", sep=",", header=T) # Metadata file
meta_AT <- as.data.frame(meta_AT)

counts_AT <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data/AdoptiveTransfer_Rawcounts.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_AT <- na.omit(counts_AT)

dim(counts_AT)
dim(meta_AT)
counts_AT <- flexiDEG.function1(counts_AT, meta_AT, # Run Function 1
                                 convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                 batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
counts_AT <- counts_AT[!grepl("^Gm[0-9]", rownames(counts_AT)), ]

# DESeq2 analysis for Adoptive Transfer
# set refs
meta_AT$Group <- relevel(factor(meta_AT$Group), ref = "OVA")
meta_AT$Time  <- relevel(factor(meta_AT$Time),  ref = "D0")

dds_AT <- DESeqDataSetFromMatrix(countData = counts_AT, colData = meta_AT,
                                 design = ~ Time + Group + Time:Group)
dds_AT <- DESeq(dds_AT)
resultsNames(dds_AT)
res_BDC_vs_OVA_d0 <- results(dds_AT, name = "Group_BDC_vs_OVA")

res_BDC_vs_OVA_d10 <- results(
  dds_AT,
  contrast = list(c("Group_BDC_vs_OVA","TimeD10.GroupBDC"))
)

# Import
AT_genes_D10  <- as.data.frame(res_BDC_vs_OVA_d10);  AT_genes_D10$gene  <- rownames(res_BDC_vs_OVA_d10)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_AT_D10  <- AT_genes_D10$stat;  names(lfc_vector_AT_D10)  <- rownames(AT_genes_D10)

# Drop NAs
lfc_vector_AT_D10  <- lfc_vector_AT_D10[!is.na(lfc_vector_AT_D10)]
lfc_vector_AT_D10  <- sort(lfc_vector_AT_D10,  decreasing = TRUE)
# --- Collect each set and convert into 2-column (gs_name, gene_symbol) ---

# C8
CellTypeMSigDB_gene_sets <- msigdbr(species="Mus musculus", category="C8")
mm_c8_sets <- split(CellTypeMSigDB_gene_sets$gene_symbol, CellTypeMSigDB_gene_sets$gs_name)
mm_c8_df <- data.frame(
  gs_name = rep(names(mm_c8_sets), sapply(mm_c8_sets, length)),
  gene_symbol = unlist(mm_c8_sets)
)


# Hallmark
hallmark <- msigdbr(species = "Mus musculus", category  = "H")
mm_h_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
mm_h_df <- data.frame(
  gs_name = rep(names(mm_h_sets), sapply(mm_h_sets, length)),
  gene_symbol = unlist(mm_h_sets)
)


# KEGG
kegg_all <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:KEGG_LEGACY")
mm_kegg_sets <- split(kegg_all$gene_symbol, kegg_all$gs_name)
mm_kegg_df <- data.frame(
  gs_name = rep(names(mm_kegg_sets), sapply(mm_kegg_sets, length)),
  gene_symbol = unlist(mm_kegg_sets)
)

# --- Final combined TERM2GENE data frame ---
mm_all_df <- rbind(mm_c8_df, mm_h_df, mm_kegg_df)

library(clusterProfiler)
library(msigdbr)
# Day 7
gsea_results_AT <- GSEA(
  geneList      = lfc_vector_AT_D10,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  #keyType       = "SYMBOL",       # <- tell it explicitly
  TERM2GENE     = mm_all_df
)
gsea_results_AT_df <- as.data.frame(gsea_results_AT)
