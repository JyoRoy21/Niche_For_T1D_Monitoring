
library(ComplexHeatmap)
library(EnhancedVolcano)
# Install if not already installed
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")

# Load the package
library(circlize)


## Scaffold Metabolomics- Overt Progressor vs Non-Progressor-Late
Metabolites_ScafLateOvertT1D<- read.csv("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Scaffold_Plasma/Results/Scaffold/Late Timepoint/Volcano_Scaf_Late_OvertProgressorVsNonProgressor.csv", sep=",", header=T,check.names = FALSE)
Metabolites_ScafLateOvertT1D<-as.data.frame(Metabolites_ScafLateOvertT1D)
head(Metabolites_ScafLateOvertT1D)
colnames(Metabolites_ScafLateOvertT1D)[1]<-"Metabolites"
colnames(Metabolites_ScafLateOvertT1D)[3]<-"log2FC"
colnames(Metabolites_ScafLateOvertT1D)[5]<-"neglog10p"

# Ensure numeric
Metabolites_ScafLateOvertT1D$log2FC <- as.numeric(Metabolites_ScafLateOvertT1D$log2FC)
Metabolites_ScafLateOvertT1D$raw.pval   <- as.numeric(Metabolites_ScafLateOvertT1D$raw.pval)
rownames(Metabolites_ScafLateOvertT1D)<-Metabolites_ScafLateOvertT1D$Metabolites
# Thresholds
p_thr  <- 0.10
fc_thr <- 1

library(EnhancedVolcano)

df <- Metabolites_ScafLateOvertT1D

# Optional: create adjusted p-values (BH method) for future filtering if needed
df$padj <- p.adjust(df$raw.pval, method = "BH")

# Build keyvals color mapping
keyvals <- ifelse(df$raw.pval < 0.1 & df$log2FC > 1, "#7E0000",     # Dark red high conf up
                  ifelse(df$raw.pval < 0.1 & df$log2FC < -1, "#002E66",    # Dark blue high conf down
                         ifelse(df$log2FC > 1 & df$raw.pval < 0.05, "#FF0000",    # Bright red nominal up
                                ifelse(df$log2FC < -1 & df$raw.pval < 0.05, "#1E90FF",   # Bright blue nominal down
                                       "gray"))))                                               # Non-significant

# Assign names for legend mapping
names(keyvals) <- ifelse(df$raw.pval < 0.1 & df$log2FC > 1, "Upregulated",
                         ifelse(df$raw.pval < 0.1 & df$log2FC < -1, "Downregulated",
                                ifelse(df$log2FC > 1 & df$raw.pval < 0.05, "p < 0.05 & log[2]FC > 1",
                                       ifelse(df$log2FC < -1 & df$raw.pval < 0.05, "p < 0.05 & log[2]FC < -1",
                                              "Not Significant"))))

# Determine x-axis limits
xmax <- max(2, ceiling(max(abs(df$log2FC), na.rm = TRUE)))

# Volcano plot
EnhancedVolcano(
  df,
  lab = df$Metabolites,
  selectLab = df$Metabolites[df$raw.pval < 0.1 & abs(df$log2FC) >= 1],
  x = "log2FC",
  y = "raw.pval",
  pCutoff = 0.1,
  FCcutoff = 1,
  xlab = expression("log"[2]~"Fold Change (Progressor / Non-Progressor)"),
  ylab = expression("-log"[10]~"(p-value)"),
  title = "Late Timepoint: Overt Progressor vs Non-Progressor (Scaffold Metabolomics)",
  xlim = c(-xmax, xmax),
  ylim=c(0,8),
  boxedLabels = TRUE,
  pointSize = 2.5,
  labSize = 6,
  colAlpha = 0.8,
  drawConnectors = TRUE,
  colCustom = keyvals,
  legendPosition = "right"
)

## Plasma Metabolomics- Overt Progressor vs Non-Progressor-Late

Metabolites_PlasmaLateOvertT1D<- read.csv("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Scaffold_Plasma/Results/Blood/Late Timepoint/Volcano_Plasma_OvertProgressorVsNonProgressor.csv", sep=",", header=T,check.names = FALSE)
Metabolites_PlasmaLateOvertT1D<-as.data.frame(Metabolites_PlasmaLateOvertT1D)
head(Metabolites_PlasmaLateOvertT1D)
colnames(Metabolites_PlasmaLateOvertT1D)[1]<-"Metabolites"
colnames(Metabolites_PlasmaLateOvertT1D)[3]<-"log2FC"
colnames(Metabolites_PlasmaLateOvertT1D)[5]<-"neglog10p"

# Ensure numeric
Metabolites_PlasmaLateOvertT1D$log2FC <- as.numeric(Metabolites_PlasmaLateOvertT1D$log2FC)
Metabolites_PlasmaLateOvertT1D$raw.pval   <- as.numeric(Metabolites_PlasmaLateOvertT1D$raw.pval)
rownames(Metabolites_PlasmaLateOvertT1D)<-Metabolites_PlasmaLateOvertT1D$Metabolites
# Thresholds
p_thr  <- 0.10
fc_thr <- 1

library(EnhancedVolcano)

df_Plasma <- Metabolites_PlasmaLateOvertT1D

# Optional: create adjusted p-values (BH method) for future filtering if needed
df_Plasma$padj <- p.adjust(df_Plasma$raw.pval, method = "BH")

# Build keyvals color mapping
keyvals <- ifelse(df_Plasma$raw.pval < 0.1 & df_Plasma$log2FC > 1, "#7E0000",     # Dark red high conf up
                  ifelse(df_Plasma$raw.pval < 0.1 & df_Plasma$log2FC < -1, "#002E66",    # Dark blue high conf down
                         ifelse(df_Plasma$log2FC > 1 & df_Plasma$raw.pval < 0.05, "#FF0000",    # Bright red nominal up
                                ifelse(df_Plasma$log2FC < -1 & df_Plasma$raw.pval < 0.05, "#1E90FF",   # Bright blue nominal down
                                       "gray"))))                                               # Non-significant

# Assign names for legend mapping
names(keyvals) <- ifelse(df_Plasma$raw.pval < 0.1 & df_Plasma$log2FC > 1, "Upregulated",
                         ifelse(df_Plasma$raw.pval < 0.1 & df_Plasma$log2FC < -1, "Downregulated",
                                ifelse(df_Plasma$log2FC > 1 & df_Plasma$raw.pval < 0.05, "p < 0.05 & log[2]FC > 1",
                                       ifelse(df_Plasma$log2FC < -1 & df_Plasma$raw.pval < 0.05, "p < 0.05 & log[2]FC < -1",
                                              "Not Significant"))))

# Determine x-axis limits
xmax <- max(2, ceiling(max(abs(df_Plasma$log2FC), na.rm = TRUE)))

# Volcano plot
EnhancedVolcano(
  df_Plasma,
  lab = df_Plasma$Metabolites,
  selectLab = df_Plasma$Metabolites[df_Plasma$raw.pval < 0.1 & abs(df_Plasma$log2FC) >= 1],
  x = "log2FC",
  y = "raw.pval",
  pCutoff = 0.1,
  FCcutoff = 1,
  xlab = expression("log"[2]~"Fold Change (Progressor / Non-Progressor)"),
  ylab = expression("-log"[10]~"(p-value)"),
  title = "Late Timepoint: Overt Progressor vs Non-Progressor (Plasma Metabolomics)",
  xlim = c(-xmax, xmax),
  ylim=c(0,6),
  boxedLabels = TRUE,
  pointSize = 2.5,
  labSize = 6,
  colAlpha = 0.8,
  drawConnectors = TRUE,
  colCustom = keyvals,
  legendPosition = "right"
)
