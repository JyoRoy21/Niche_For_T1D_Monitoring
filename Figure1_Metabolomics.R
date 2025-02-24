## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 
library(ggnewscale)
install.packages("viridis")  # Install if not installed
library(viridis)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(ggplot2)
library(ggnewscale)
library(viridis)
library(RColorBrewer)  # For Nature-compatible color palettes
# Function to Expand Metabolites (Convert from "X;Y;Z" format to long format)
library(tidyr)
## ----DotPlot Pathways-----------------------------------------------------------------------------------------------------------------------------------

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
library(dplyr)

# Calculate -log10(FDR) for Early timepoint
early_fdr <- filtered_pathways %>%
  filter(Timepoint == "Early") %>%
  mutate(log10_FDR = -log10(FDR))

# Reorder Pathway factor levels based on descending log10_FDR
filtered_pathways <- filtered_pathways %>%
  mutate(Pathway = factor(Pathway, levels = early_fdr$Pathway[order(early_fdr$log10_FDR, decreasing = FALSE)]))

# Assuming 'filtered_pathways' is your dataframe
filtered_pathways <- filtered_pathways %>%
  group_by(Pathway) %>%
  filter(any(Impact > 0.0)) %>%
  ungroup()

# Create the dotplot
# Plot
# PLOT


filtered_pathways$`-log10(p)`
ggplot(filtered_pathways, aes(x = `-log10(p)`, y = Pathway, size = Impact)) +
  
  # Early timepoint (Yellow to Red)
  geom_point(data = subset(filtered_pathways, Timepoint == "Early"), aes(color = -log10(`Raw p`)), alpha = 0.8) +
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "Early") +
  new_scale_color() +
  
  # Intermediate timepoint (Reversed Viridis: Light to Deep Green)
  geom_point(data = subset(filtered_pathways, Timepoint == "Intermediate"), aes(color = -log10(`Raw p`)), alpha = 0.8) +
  scale_color_gradient(low = "#D5E21A", high = "#008000", name = "Intermediate") +  # Light to deep green
  new_scale_color() +
  
  # Late timepoint (Blues)
  geom_point(data = subset(filtered_pathways, Timepoint == "Late"), aes(color = -log10(`Raw p`)), alpha = 0.8) +
  scale_color_gradient(low = "#BFD3E6", high = "#08306B", name = "Late") +  # Light to deep blue
  
  # Dot size scale with custom breaks and labels
  scale_size_continuous(
    name = "Impact",
    breaks = c(0.2, 0.3, 0.4, 0.5),
    labels = c("≤0.2", "0.3", "0.4", "0.5"),
    range = c(2, 8)
  ) +
  
  # Improve plot theme
  theme_minimal(base_size = 14) +
  labs(title = "Metabolomics Pathway Enrichment",
       x = expression(-log[10](P[adj])),
       y = "Pathways") +
  
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

## ----Raindrop plot-----------------------------------------------------------------------------------------------------------------------------------

setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold")

# Load Metabolite-Pathway Mapping
pathway_metabolite_list <- read_excel("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold/Scaffold_SigMetabPathways_Metabolites.xls", sheet = 1)

pathway_metabolite_long <- pathway_metabolite_list %>%
  separate_rows(Metabolite, sep = ";") %>%
  mutate(Metabolite = trimws(Metabolite))  # Remove extra spaces


intermediate_metabs <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold/Intermediate Timepoint/Volcano_scaf_intermediate_unfiltered_FDR.csv", 
                                  sep=",", header=TRUE, check.names=FALSE)
colnames(intermediate_metabs)[1] <- "Metabolite"
intermediate_metabs <- intermediate_metabs %>% 
  na.omit() %>% 
  mutate(Timepoint = "Intermediate")


# Load Metabolite Data for Each Timepoint
early_metabs <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold/Early Timepoint/Volcano_scaf_early_unfiltered_FDR.csv",
                           sep=",", header=T, check.names = FALSE) 
colnames(early_metabs)[1] <- "Metabolite"
early_metabs <- early_metabs %>% 
  na.omit() %>% 
  mutate(Timepoint = "Early")


intermediate_metabs <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold/Intermediate Timepoint/Volcano_scaf_intermediate_unfiltered_FDR.csv", 
                                  sep=",", header=TRUE, check.names=FALSE)
colnames(intermediate_metabs)[1] <- "Metabolite"
intermediate_metabs <- intermediate_metabs %>% 
  na.omit() %>% 
  mutate(Timepoint = "Intermediate")


late_metabs <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold/Late Timepoint/Volcano_scaf_Late_unfiltered_FDR.csv", 
                          sep=",", header=T, check.names = FALSE)
colnames(late_metabs)[1] <- "Metabolite"
late_metabs <- late_metabs %>% 
  na.omit() %>% 
  mutate(Timepoint = "Late")


# Combine all datasets
all_metabs <- bind_rows(early_metabs, intermediate_metabs, late_metabs)
colnames(pathway_metabolite_long)
colnames(all_metabs)

pathway_metabolite_long$Metabolite <- as.character(trimws(pathway_metabolite_long$Metabolite))
all_metabs$Metabolite <- as.character(trimws(all_metabs$Metabolite))

# Merge Metabolites with Pathways
all_metabs_with_pathways <- all_metabs %>%
  inner_join(pathway_metabolite_long, by = "Metabolite")


# Adjust x_pos calculation to increase space between timepoints within each pathway
all_metabs_with_pathways <- all_metabs_with_pathways %>%
  mutate(Pathway = factor(Pathway, levels = unique(Pathway)),
         Timepoint = factor(Timepoint, levels = c("Early", "Intermediate", "Late")),
         x_pos = as.numeric(Pathway) * 10 + as.numeric(Timepoint) * 2 - 4)  # Centering timepoints within each pathway

# Define breaks and labels
breaks <- seq(10, max(all_metabs_with_pathways$x_pos), by = 10)
labels <- levels(all_metabs_with_pathways$Pathway)


# Plot
ggplot(all_metabs_with_pathways, aes(x = x_pos, y = `log2(FC)`, size = -log10(p.ajusted))) +
  
  # Early Timepoint (Yellow to Red)
  geom_point(data = subset(all_metabs_with_pathways, Timepoint == "Early"), aes(color = abs(`log2(FC)`)), alpha = 0.8) +
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "Early") +
  new_scale_color() +
  
  # Intermediate Timepoint (Reversed Green Gradient)
  geom_point(data = subset(all_metabs_with_pathways, Timepoint == "Intermediate"), aes(color = abs(`log2(FC)`)), alpha = 0.8) +
  scale_color_gradient(low = "#D5E21A", high = "#008000", name = "Intermediate") +
  new_scale_color() +
  
  # Late Timepoint (Blue Shades)
  geom_point(data = subset(all_metabs_with_pathways, Timepoint == "Late"), aes(color = abs(`log2(FC)`)), alpha = 0.8) +
  scale_color_gradient(low = "#BFD3E6", high = "#08306B", name = "Late") +
  
  # Add vertical separator lines between pathways (solid lines)
  geom_vline(xintercept = seq(15, max(all_metabs_with_pathways$x_pos), by = 10), linetype = "solid", color = "grey50") +
  
  # Add a denser dotted horizontal line at y = 0 with custom linetype
  geom_hline(yintercept = 0, linetype = "1212", color = "grey50", size = 1.5) +
  
  scale_size(range = c(2, 8), name = bquote(bold("-log10(" * p[adj] * ")")))+

  # Custom x-axis labels (Only showing pathway names)
  scale_x_continuous(breaks = breaks, labels = labels) +
  
  # Improve plot theme
  theme_minimal(base_size = 14) +
  labs(title="",
    x = "Pathways",
       y = expression(log[2]("FC Progressor Vs Non-Progressor"))) +
  
  theme(
    axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "Arial"),  # Adjusted hjust for right shift
    axis.text.y = element_text(size = 16, family = "Arial"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank(),  # Remove vertical grid
    panel.grid.major.y = element_blank(),  # Remove all horizontal gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border around the plot
  )


