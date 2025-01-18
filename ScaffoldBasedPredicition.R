
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
kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names
biocarta_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:BIOCARTA") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names
pwl_msigdbr <- c(pwl_hallmark, pwl_reactome, pwl_kegg, pwl_biocarta) # Compile them all
length(pwl_msigdbr)


getwd()

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


# Counts Data Importing----
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

# # Selecting Genes for Input to Machine Learning Model 
# early_data <- combined_counts[, meta_combined$Time == "Early" & meta_combined$Batch!= '4']
# early_data <- na.omit(early_data)
# meta_early <- meta_combined[meta_combined$Time == "Early" & meta_combined$Batch!= '4', ]
# 
# early_data <- flexiDEG.function1(early_data, meta_early, # Run Function 1
#                                  convert_genes = F, exclude_riken = T, exclude_pseudo = F,
#                                  batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  100000000
# 
# # Remove rows where row names start with "Gm" followed by a digit
# early_data <- early_data[!grepl("^Gm[0-9]", rownames(early_data)), ]
# # DESeq2 analysis for Early Time
# dds_early <- DESeqDataSetFromMatrix(countData = early_data, colData = meta_early, design = ~ Batch+Group)
# dds_early <- DESeq(dds_early)
# results_early <- as.data.frame(results(dds_early, contrast = c("Group", "Progressor", "Non-Progressor")))
# # Replace NA values in the pvalue column with 1
# results_early$pvalue[is.na(results_early$pvalue)] <- 1
# 
# # Verify the changes
# summary(results_early$pvalue)
# 
# 
# # Define the file path and name
# output_file <- "results_early_DESeq2.csv"
# # Save the data frame as a CSV file
# write.csv(results_early, file = output_file, row.names = TRUE)
# 
# 
# # Filter significant genes based on the p-value threshold
# significant_genes <- results_early[results_early$pvalue < 0.05, ]
# # Count the number of significant genes
# num_significant_genes <- nrow(significant_genes)
# # Print the result
# cat("Number of significant genes (p-value < 0.05):", num_significant_genes, "\n")
# 
# # Subset the significant gene rows from 'combined_counts'
# significant_gene_names <- rownames(significant_genes)
# combined_counts <- combined_counts[significant_gene_names, ]


# Normalizing Genes Input to Machine Learning Model ---- 
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

#t-test Based Filtering
# Step 1: Filter `meta_combined_early` to exclude rows with `N/A` in the `Group` column
#meta_combined_early_clean <- meta_combined_early[meta_combined_early$Batch != "4", ]
meta_combined_early_clean <- meta_combined_early[
  !(meta_combined_early$Batch == "4" ), 
]
# Step 2: Get the list of sample names to keep
valid_samples <- meta_combined_early_clean$Samples
# Step 3: Filter `case1_f1` to include only rows corresponding to the valid samples
case1_f1_clean <- case1_f1[,colnames(case1_f1) %in% valid_samples ]

#case1_f3 <- flexiDEG.function2(case1_f1_clean, meta_combined_early_clean) # Run Function 2
case1_f2<- flexiDEG_ttest(case1_f1_clean, meta_combined_early_clean)

# Step 4: Subset rows in case1_f1 that are present in case1_f2
case1_f1_subset <- case1_f1[rownames(case1_f1) %in% rownames(case1_f2), ]
write.csv(case1_f1_subset  , "PredictionInput_ttestfiltered.csv", row.names = TRUE)

# Machine Learning Model P-Alpha Selection ---- 
case1_f1_clean <- case1_f1_clean[rownames(case1_f1_clean) %in% rownames(case1_f2), ]
case1_f4 <- flexiDEG.function4(case1_f1_clean  , meta_combined_early_clean,validation_option = 3) 


EN1 <- na.omit(case1_f1_clean[unique(rownames(case1_f1_clean)[as_vector(case1_f4[[1]])]), ])
EN2 <- na.omit(case1_f1_clean[unique(rownames(case1_f1_clean)[as_vector(case1_f4[[2]])]), ]) 
EN3 <- na.omit(case1_f1_clean[unique(rownames(case1_f1_clean)[as_vector(case1_f4[[3]])]), ]) 

# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors (meta_combined_early_clean)
unique_colSide <- unique(colSide)


#dev.off()
heatmap.2(as.matrix(EN3), scale="row", col=coul_gsva, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN3), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=colSide)

# Machine Learning Model -Prediction ---- 

# Extract Sample names and number
sample_name<- meta_combined_early_clean$Samples #Sample names
sample_num <- length(sample_name) # Number of samples
user_input<-"Samples"
batch <- c(rep(user_input,sample_num))
grouped<- meta_combined_early_clean$Group
n=length(unique(grouped)) 
# Combine lines to convert 'grouped' into a factor with levels 0 to (n-1)
groups <- factor(as.integer(factor(grouped, labels = 1:n)) - 1, levels = 0:(n-1))
# Create a data frame with the factor variable
xfactors <- cbind.data.frame(groups)
# Transpose the data matrix bc_best to create the dataset
bc_best <- case1_f1_clean
dataset <- t(bc_best) # Load the data
dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
set.seed(123) # Split the data into training and test set
training.samples <- dataset$groups %>% createDataPartition(p = 1, list = F)
train.data  <- dataset[training.samples, ]

# Prepare test data
# Step 1: Select all batch 4
meta_combined_testset<- meta_combined_early[meta_combined_early$Batch == "4", ]
# Step 2: Get the list of sample names to keep
test_samples <- meta_combined_testset$Samples
# Step 3: Filter `case1_f1` to include only rows corresponding to the valid samples
case_test <- case1_f1_subset[,colnames(case1_f1_subset) %in% test_samples ]
#case_test <- case_test[rownames(case1_f1) %in% rownames(case1_f2), ]

test.data <- t(case_test)
test.data <- as.data.frame(test.data) # Add factor as the first column
test.data <- test.data[complete.cases(test.data), ] # For cases without missed items

x <- model.matrix(groups~., train.data)[,-1] # Predictor variables
y <- train.data$groups # Outcome variable
num_groups <- length(unique(y)) 
foldid <- sample(1:length(y), size = length(y), replace = TRUE)
# Create an empty vector to store the AUC values for each alpha
alpha_val=0
prediction_results <- data.frame()  # Dataframe to store predictions for all iterations
iter=1000

pb <- txtProgressBar(min = 1, max = iter, style = 3)
for (j in 1:iter) {
  min_samples_per_group <- 4  # Adjust as needed
  criterion_met <- FALSE

  
  while (!criterion_met) {
    # Bootstrap sampling
    boot_indices <- sample(length(y), replace = TRUE)
    samples_per_group <- table(y[boot_indices])
    if (all(samples_per_group >= min_samples_per_group)) {
      criterion_met <- TRUE
    }
  
  }
  
  if (!criterion_met) {
    message("Failed to meet sampling criterion after max attempts. Skipping iteration.")
    next
  }
  
  # Fit Elastic Net model
  cvfit <- tryCatch(
    suppressMessages(cv.glmnet(x[boot_indices, ], y[boot_indices], alpha = alpha_val, foldid = foldid, family = "binomial")),
    error = function(e) NULL
  )
  
  if (is.null(cvfit)) {
    j <- j - 1
    next
  }
  
  best.lambda <- cvfit$lambda.min
  fit <- tryCatch(
    suppressMessages(glmnet(x[boot_indices, ], y[boot_indices], alpha = alpha_val, lambda = best.lambda, family = "binomial")),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    j <- j - 1
    next
  }
  
  # Test data predictions
  x_test <- model.matrix( ~ ., test.data)[, -1]  # Predictor variables for test data
  probabilities <- tryCatch(
    predict(fit, newx = x_test, type = "response", s = best.lambda),
    error = function(e) NULL
  )
  
  if (is.null(probabilities)) {
    j <- j - 1
    next
  }
  
 
  
  # Convert s1 to 1 or 0 based on threshold (s1 > 0.5 is classified as 1, otherwise 0)
  pred_labels <- ifelse(probabilities[, 1] > 0.5, 1, 0)
  
  # Create a dataframe for the current iteration results
  iteration_results <- as.data.frame(t(pred_labels))  # Transpose so each column represents a sample
  colnames(iteration_results) <- rownames(probabilities)  # Set column names to sample names
  iteration_results$iteration <- j  # Add iteration number as a column
  
  # Reorder columns so that the iteration number is the first column
  iteration_results <- iteration_results[, c("iteration", colnames(iteration_results)[1:ncol(iteration_results)-1])]
  
  # Bind the results from this iteration to the overall prediction_results
  prediction_results <- rbind(prediction_results, iteration_results)
  # Update the progress bar
  setTxtProgressBar(pb, j)
}

# Close the progress bar after completion
close(pb)

# Sum the values for each column (excluding the iteration column)
sums_per_sample <- colSums(prediction_results[, -1])

# Display the sums for each sample
print(sums_per_sample)

# Label as 'Progressor' if sum > 500, else 'Non-Progressor'
labels <- ifelse(sums_per_sample > 500, "Progressor", "Non-Progressor")

# Add labels to the dataframe
labels_df <- data.frame(sample = names(sums_per_sample), label = labels)

# Display the labeled results
print(labels_df)

