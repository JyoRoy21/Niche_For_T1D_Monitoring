library(OmicsPLS)

getwd()
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/JointMultiOmicsAnalysis/")
metabolomics_data <- read.csv("MetabolomicsData_counts_Early.csv",check.names = FALSE)
transcriptomics_data <- read.csv("TranscriptomicsData_counts_Early.csv",check.names = FALSE)
rownames(metabolomics_data)<-metabolomics_data$`Compound Method`
metabolomics_data$`Compound Method`<- NULL
#rownames(transcriptomics_data)<-transcriptomics_data$external_gene_name

# Remove rows with duplicate external_gene_name in transcriptomics_data
transcriptomics_data <- transcriptomics_data[!duplicated(transcriptomics_data$external_gene_name), ]
# Set rownames again after removing duplicates
rownames(transcriptomics_data) <- transcriptomics_data$external_gene_name
transcriptomics_data$`external_gene_name`<- NULL



#Procedure without scaling
X = metabolomics_data
Y = transcriptomics_data
# Impute missing values in both X and Y
X_imputed <- impute_matrix(X)
Y_imputed <- impute_matrix(Y)

# Now, you can proceed with the analysis using the imputed data
X <- X_imputed
Y <- Y_imputed

X<-t(X)
Y<-t(Y)
set.seed(12)
#crossval_o2m_adjR2(X, Y, 1:3, 0:3, 0:3,  nr_cores = 8,nr_folds = nrow(X))
crossval_o2m(X, Y, 1:3, 0:3, 0:3,  nr_cores = 8, nr_folds = nrow(X))

fit0 = o2m(X, Y, n=2, nx=1, ny=1)
fit0
summary(fit0)

#Procedure WITH scaling
X1 = scale(X, scale=F)
Y1 = scale(Y, scale=F)

set.seed(1221L)
crossval_o2m_adjR2(X1, Y1, 1:3, 0:3, 0:3, nr_folds = nrow(X1))
crossval_o2m(X1, Y1, 1:3, 0:3, 0:3, nr_folds = nrow(X1))

fit1 = o2m(X1, Y1, ?, ?, ?)
fit1
summary(fit1)