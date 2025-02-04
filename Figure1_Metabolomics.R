library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 
library(ggnewscale)
install.packages("viridis")  # Install if not installed
library(viridis)


## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold")
early_metab_pathway <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold/Early Timepoint/pathway_results_scaf_early.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
early_metab_pathway <- na.omit(early_metab_pathway)

intermediate_metab_pathway <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold/Intermediate Timepoint/pathway_results_scaf_intermediate.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
intermediate_metab_pathway <- na.omit(intermediate_metab_pathway)

late_metab_pathway <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold/Late Timepoint/pathway_results_scaf_late.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
late_metab_pathway <- na.omit(late_metab_pathway)


# Add a timepoint column to each dataset
early_metab_pathway$Timepoint <- "Early"
intermediate_metab_pathway$Timepoint <- "Intermediate"
late_metab_pathway$Timepoint <- "Late"

# Bind all data together
all_metab_pathways <- bind_rows(early_metab_pathway, intermediate_metab_pathway, late_metab_pathway)

# Rename first column as "Pathway" if it's unnamed
colnames(all_metab_pathways)[1] <- "Pathway"

# Select pathways that have at least one Raw p < 0.05
filtered_pathways <- all_metab_pathways %>%
  group_by(Pathway) %>%
  filter(any(`FDR` < 0.05)) %>%
  ungroup()

# Create the dotplot
# Plot
# PLOT
library(ggplot2)
library(ggnewscale)
library(viridis)
library(RColorBrewer)  # For Nature-compatible color palettes

ggplot(filtered_pathways, aes(x = -log10(FDR), y = Pathway, size = Impact)) +
  
  # Early timepoint (Yellow to Red)
  geom_point(data = subset(filtered_pathways, Timepoint == "Early"), aes(color = -log10(`Raw p`)), alpha = 0.8) +
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "Early") +
  new_scale_color() +
  
  # Intermediate timepoint (Reversed Viridis: Light to Deep Green)
  geom_point(data = subset(filtered_pathways, Timepoint == "Intermediate"), aes(color = -log10(`Raw p`)), alpha = 0.8) +
  scale_color_gradient(low = "#D5E21A", high = "#008000", name = "Intermediate") +  # Light to deep green
  new_scale_color() +
  
  # Late timepoint (Blues, correctly applied)
  geom_point(data = subset(filtered_pathways, Timepoint == "Late"), aes(color = -log10(`Raw p`)), alpha = 0.8) +
  scale_color_gradient(low = "#BFD3E6", high = "#08306B", name = "Late") +  # Light to deep blue
  
  # Dot size scale
  scale_size(range = c(2, 8)) +
  
  # Improve plot theme
  theme_minimal(base_size = 14) +
  labs(title = "Metabolomics Pathway Enrichment",
       x = expression(-log[10](FDR)),
       y = "Pathways",
       size = "Impact") +
  
  theme(
    axis.text.y = element_text(size = 12, family = "Arial"),
    axis.text.x = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_blank(),  # Remove only x-axis grid
    panel.grid.major.y = element_line(color = "grey80"),  # Keep y-axis grid
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )

