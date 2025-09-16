
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
pwl_msigdbr <- c(pwl_hallmark, pwl_reactome, pwl_kegg) # Compile them all
length(pwl_msigdbr)


getwd()

#Metadata Importing
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_MetabolomicsCohort_PathwayAnalysis.csv", sep=",", header=T) # Metadata file

meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_batch3 <- as.data.frame(meta_batch3)

# Merge metadata by columns (i.e., add samples from Batch 2 to Batch 1)
meta_combined <- rbind(meta_batch1, meta_batch2,meta_batch3)

# Preview the combined metadata
head(meta_combined)


# Counts Data Importing----
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




# Normalizing Genes Input to Machine Learning Model ---- 
combined_counts <- combined_counts[, meta_combined$Samples]  # Ensure Sample_IDs match column names in 

# Filter Data for time
meta_combined_early <- meta_combined[meta_combined$Time == "Early", ]

# Filter combined_counts to keep only the samples in the subset metadata
combined_counts_early <- combined_counts[, colnames(combined_counts) %in% meta_combined_early$Samples]

combined_counts_early <- combined_counts_early[, meta_combined_early$Samples]  # Ensure Sample_IDs match column names in 
# Remove rows with NA values
combined_counts_early <- combined_counts_early[complete.cases(combined_counts_early), ]


## Option 1: Only select t-test filtered genes ----

case1_f1 <- flexiDEG.function1(combined_counts_early, meta_combined_early, # Run Function 1
                               convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                               batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 2, 0, 15

case1_f1 <- case1_f1[!grepl("^Gm[0-9]", rownames(case1_f1)), ]



# DESeq2 analysis for Week 6
dds_Early <- DESeqDataSetFromMatrix(countData = case1_f1, colData = meta_combined_early, design = ~ Batch+Group)
dds_Early <- DESeq(dds_Early)
results_Early <- as.data.frame(results(dds_Early, contrast = c("Group", "Progressor", "Non-Progressor")))
results_Early$gene <- rownames(results_Early)

# GSEA ----------------------------------------------------------------- 

# Prepare the ranked list for Islet
results_Early <- results_Early[, c("gene", "log2FoldChange", "padj")]
results_Early <- results_Early[!is.na(results_Early$log2FoldChange), ]
#Order by fold change
results_Early <- results_Early[order(results_Early$log2FoldChange, decreasing = TRUE), ]
results_Early$padj[is.na(results_Early$padj)] <- 1


geneList <- setNames(results_Early$log2FoldChange, results_Early$gene)
de <- names(geneList)[abs(geneList) > 1]
# library(DOSE)
# edo <- enrichDGN(de)

# ## convert gene ID to Symbol
# edox <- setReadable(edo, 'org.Mm.eg.db', 'SYMBOL')
# p1 <- cnetplot(edox, foldChange=geneList)
# ## categorySize can be scaled by 'pvalue' or 'geneNum'
# p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
# p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
# cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))


# Enrichment analysis (for example with GO Biological Processes)
ego <- enrichGO(gene         = de,
                OrgDb        = org.Mm.eg.db,
                keyType      = "SYMBOL",
                ont          = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.01,
                readable     = TRUE)
dev.off()
barplot(ego, showCategory = 10)

#Figure 5D
cnetplot(ego, node_label="category")


# Keep only top 10 pathways
#ego_top <- simplify(ego)  # remove redundant GO terms if GO
ego_top <- ego[1:10, ] 


# Gene-Concept Network Plot
cnetplot(ego_top,
         showCategory = 5,   # Top 10 pathways
         foldChange = NULL,   # Optional: add expression data to color nodes
         circular = FALSE, 
         colorEdge = TRUE)


