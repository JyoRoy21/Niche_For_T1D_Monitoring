pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, dnapath, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)
library(dplyr)

NOD_BloodMicroArray<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Input/NOD Blood Microarray/NOD_BloodMicroarray_NormalizedGeneExpresssion.csv",check.names = F)
NOD_BloodMicroArray_metadata<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Input/NOD Blood Microarray/NOD_BloodMicroarray_Metadata.csv")

T1DGene_top100<- read.csv("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered_PLSDATop100.csv", sep=",", header=T,check.names = FALSE)
T1DGene_top100<-as.data.frame(T1DGene_top100)

## ---- standardize key column names ----
# Rename the gene column to "Gene" if needed:
if (!"Gene" %in% names(NOD_BloodMicroArray)) {
  # try common alternatives; edit if your file uses something else
  gcol <- intersect(c("Gene_Symbol"), names(NOD_BloodMicroArray))[1]
  stopifnot(length(gcol) == 1)
  names(NOD_BloodMicroArray)[names(NOD_BloodMicroArray) == gcol] <- "Gene"
}

# Extract the 100-gene symbols from the first column of your list
sig_genes <- unique(T1DGene_top100[[1]]) %>% as.character()

## ---- sanity: make sure metadata SampleIDs are present in expression ----
sample_ids <- NOD_BloodMicroArray_metadata$Sample
missing_cols <- setdiff(sample_ids, names(NOD_BloodMicroArray))
if (length(missing_cols)) {
  stop(paste("These SampleIDs are not columns in the expression file:", paste(missing_cols, collapse = ", ")))
}

## ---- subset expression to the 100 genes & order columns by metadata ----
expr100 <- NOD_BloodMicroArray %>%
  filter(Gene %in% sig_genes) %>%
  distinct(Gene, .keep_all = TRUE)

# keep only the Gene column + the 6 samples from metadata, ordered
expr100 <- expr100 %>% select(Gene, all_of(sample_ids))

# set rownames and convert to matrix
rownames(expr100) <- expr100$Gene
expr100$Gene <- NULL
expr100 <- as.matrix(expr100)

## ---- z-score by gene (row-wise) for visualization & PLSDA ----
# scale() operates column-wise, so transpose, scale, transpose back
expr100_z <- t(scale(t(expr100)))  # rows mean=0, sd=1

## =========================
## 1) HEATMAP (100 genes x 6 samples)
## =========================
ann_col <- NOD_BloodMicroArray_metadata %>%
  select(Sample, Group) %>%
  column_to_rownames("Sample")

pheatmap(
  expr100_z,
  annotation_col = ann_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,     # let it try; if you want fixed order, set to FALSE
  show_rownames = TRUE,   # toggle TRUE if you want to see all genes
  fontsize_col = 11,
  main = "Top-100 Scaffold Signature in NOD Blood (z-scored per gene)"
)

## =========================
## 2) PLS-DA (should show no separation if signal is absent)
## =========================
# X must be samples x features
X <- t(expr100_z)   # now 6 x 100
Y <- factor(ann_col$Group, levels = c("Late/Non Progressor","Progressor"))  # set order

set.seed(123)
pls <- plsda(X, Y, ncomp = 2)

# Quick score plot
plotIndiv(
  pls, comp = c(1,2), group = Y, legend = TRUE,
  title = "PLS-DA: 100-gene Scaffold Signature in Blood"
)

## Optional: quick CV/permutation-style feel for separability (with n=3/3)
set.seed(123)
perf_res <- perf(pls, validation = "Mfold", folds = 3, nrepeat = 50, progressBar = FALSE)
print(perf_res$error.rate)  # expect high/near-chance error rates with no separation