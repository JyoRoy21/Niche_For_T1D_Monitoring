
library(ComplexHeatmap)

# Install if not already installed
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")

# Load the package
library(circlize)
T1DGene_top100<- read.csv("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered_PLSDATop100.csv", sep=",", header=T,check.names = FALSE)
T1DGene_top100<-as.data.frame(T1DGene_top100)
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/metadata_MetabolomicsCohort_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
rownames(T1DGene_top100)<-T1DGene_top100[,1]
meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_batch3 <- as.data.frame(meta_batch3)

# Merge metadata by columns (i.e., add samples from Batch 2 to Batch 1)
T1D_metadata_3batch <- rbind(meta_batch1, meta_batch2,meta_batch3)
T1D_metadata_3batch<-T1D_metadata_3batch[T1D_metadata_3batch$Time=="Early",]

T1DGene_top100<-T1DGene_top100[ ,T1D_metadata_3batch$Samples]
colnames(T1DGene_top100)

# Expression matrix (genes x samples)
mat <- as.matrix(T1DGene_top100)

# Metadata: group info for samples
group <- T1D_metadata_3batch$Group
names(group) <- T1D_metadata_3batch$Samples

# Annotation colors
group_col <- c("Progressor" = "#E60000", "Non-Progressor" = "#3D7A60")

# Column annotation
ha_col <- HeatmapAnnotation(
  Group = group[colnames(mat)],
  col = list(Group = group_col),
  annotation_name_side = "left",
  show_annotation_name = TRUE,
  annotation_legend_param = list(title = "Group")
)

# Row-wise Z-score scaling
mat_scaled <- t(scale(t(mat)))
library(colorRamp2)
# Load the package
library(ComplexHeatmap)
# Color function for heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# Heatmap with gene names
ht <- Heatmap(
  mat_scaled,
  name = "Z-score",
  top_annotation = ha_col,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6, fontface = "italic"),  # Italicized gene names
  show_column_names = FALSE,
  heatmap_legend_param = list(title = "Expression", legend_direction = "vertical")
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
## PLSDA----

# Required packages
library(mixOmics)
library(ggplot2)

# Prepare X and Y
X <- t(T1DGene_top100)  # Transpose: samples x genes
Y <- as.factor(T1D_metadata_3batch$Group)

# Run PLS-DA with 2 components
plsda_model <- plsda(X, Y, ncomp = 2)
plsda_model$prop_expl_var
# Extract component scores
scores <- plsda_model$variates$X  # LV1, LV2

# Build ggplot data frame
plot_df <- data.frame(
  LV1 = scores[, 1],
  LV2 = scores[, 2],
  Group = Y
)


# Define group-specific colors & shapes (Nature style)
group_colors <- c("Progressor" = "#E60000", "Non-Progressor" = "#3D7A60")
group_shapes <- c("Progressor" = 16, "Non-Progressor" = 17)
# Compute variance explained manually

xvar <- round(plsda_model$prop_expl_var$X[1] * 100, 1)
yvar <- round(plsda_model$prop_expl_var$X[2] * 100, 1)
xlab <- paste0("PLS Component 1 (", xvar, "%)")
ylab <- paste0("PLS Component 2(", yvar, "%)")


# Plot
ggplot(plot_df, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 1) +
  stat_ellipse(geom = "polygon", aes(fill = Group), alpha = 0.5, level = 0.68, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes)+
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.grid = element_blank(),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +labs(x = xlab, y = ylab)


## Enrichment Analysis
T1DGene_top100<- read.csv("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_DESEqttestfiltered_PLSDATop100.csv", sep=",", header=T,check.names = FALSE)
T1DGene_top100<-as.data.frame(T1DGene_top100)

library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(enrichplot)
library(ggplot2)
library(dplyr)
# 0) Inputs
genes <- T1DGene_top100[[1]] |> as.character() |> unique() |> na.omit()

# Optional but recommended: set a background = all genes you tested (e.g., all expressed genes)
# universe_genes <- all_expressed_mouse_symbols
universe_genes <- NULL

# 1) Build IMMUNE-FOCUSED TERM2GENE from MSigDB

msig_c7 <- msigdbr(species = "Mus musculus", category = "C7")
term2gene_c7 <- msig_c7 %>% dplyr::select(gs_name, gene_symbol)

reactome_all <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
term2gene_reactome <- reactome_all  %>% dplyr::select(gs_name, gene_symbol)

hallmark <- msigdbr(species = "Mus musculus", collection  = "H") 
term2gene_h <- hallmark %>% dplyr::select(gs_name, gene_symbol)

kegg_all <- msigdbr(species="Mus musculus", subcategory="CP:KEGG_LEGACY") 
term2gene_kegg <- kegg_all %>% dplyr::select(gs_name, gene_symbol)


## D) Combine and de-dup
TERM2GENE <- bind_rows(
  mutate(term2gene_c7, source = "C7"),
  mutate(term2gene_reactome, source = "REACTOME"),
  mutate(term2gene_h, source = "HALLMARK"),
  mutate(term2gene_kegg, source = "KEGG"),
) |>
  distinct(gs_name, gene_symbol)

# 2) Over-representation analysis (user-defined gene sets)
EnrichedPathways_Top100 <- enricher(
  gene = genes,
  TERM2GENE = dplyr::select(TERM2GENE, gs_name, gene_symbol),
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff  = 1,
  minGSSize = 5
)
EnrichedPathways_Top100 <- as.data.frame(EnrichedPathways_Top100)
head(EnrichedPathways_Top100, 10)               # first 10 rows in console

library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

plot_pathway_counts <- function(EnrichedPathways_Top100,
                                min_genes = 4,
                                top_n = 25,
                                title = "Biomarkers per Pathway") {
  if (is.null(EnrichedPathways_Top100) || !nrow(EnrichedPathways_Top100)) {
    return(ggplot() + ggtitle(paste0(title, " (no pathways found)")))
  }
  
  df <- EnrichedPathways_Top100 %>%
    filter(Count >= min_genes) %>%
    mutate(Pathway = ifelse(!is.na(Description), Description, ID)) %>%
    arrange(desc(Count)) %>%
    slice_head(n = top_n)
  
  if (!nrow(df)) {
    return(ggplot() + ggtitle(paste0(title, " (no pathways with ≥", min_genes, " genes)")))
  }
  
  df <- df %>%
    mutate(Pathway = str_wrap(Pathway, width = 55),
           Pathway = factor(Pathway, levels = rev(unique(Pathway))))
  
  ggplot(df, aes(x = Pathway, y = Count, fill = Count)) +
    geom_col(width = 0.65) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = c(0, 0.04), breaks = pretty_breaks(4)) +
    # Alternate muted warm palette: deep gold → orange → brick red
    scale_fill_gradientn(
      colours = c("#F4A261", "#E76F51", "#7B2D26"),
      name = "Gene Count"
    ) +
    labs(x = NULL, y = "Number of genes in pathway", title = title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.3),
      axis.line.x = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.text.y = element_text(color = "black", size = 12),     # bigger Y tick labels
      axis.text.x = element_text(color = "black", size = 12),     # bigger X tick labels
      axis.title.y = element_text(size = 14, face = "bold"),      # Y axis title bigger
      axis.title.x = element_text(size = 14, face = "bold"),      # X axis title bigger
      plot.title = element_text(face = "bold", size = 18, hjust = 0),
      legend.title = element_text(face = "bold", size = 14),      # bigger legend title
      legend.text  = element_text(size = 12),                     # bigger legend labels
      legend.key.height = unit(5, "mm"),
      legend.key.width  = unit(3, "mm"),
      plot.margin = margin(8, 12, 8, 12)
    )
}


# Example call
plot_pathway_counts(EnrichedPathways_Top100, min_genes = 4, top_n = 50, title = "Enriched Pathways")
write.csv(EnrichedPathways_Top100,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PredictionInput_Top100_EnrichedPathways.csv")

