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
Mac_T1D_CombinedTimepoints_Cellwise_Factors_pvalue_log2FC_
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/Spectra/T1D_Timepoints/Macrophage/")

Macrophage_spectra_pathway <- as.data.frame(read.table("Mac_T1D_CombinedTimepoints_Cellwise_Factors_pvalue_log2FC_.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
Macrophage_spectra_pathway <- na.omit(Macrophage_spectra_pathway)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Ensure the dataset is properly formatted
Macrophage_spectra_pathway <- Macrophage_spectra_pathway %>%
  arrange(log2FC)  # Sorting by log2FC for better visualization


# Filter pathways with FDR < 0.05 and remove specific factors
filtered_pathways <- Macrophage_spectra_pathway %>%
  filter(FDR < 0.05 & 
           !Factor %in% c("52-X-global-X-52", "4-X-global-X-all_circadian-rhythm", 
                          "137-X-global-X-137", "48-X-global-X-48", "135-X-global-X-135",
                          "70-X-global-X-all_MET_metabolism")) %>%
  arrange(log2FC)  # Sort by log2FC

# Remove "ABC-X-global-X-" from Factor names using gsub
# Remove "ABC-X-global-X-" and "all_" from Factor names, and capitalize first letter

# Define colors for up/downregulated pathways
filtered_pathways$color <- ifelse(filtered_pathways$log2FC > 0, "#D73027", "#4575B4")


# Plot the barplot with enhanced aesthetics
ggplot(filtered_pathways, aes(x = reorder(Factor, log2FC), y = log2FC, fill = color)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_identity() +  # Uses the manually defined colors
  coord_flip() +  # Flip axis for better readability
  labs(x = "Pathways", 
       y = expression(paste("log" [2], "(FC - Progressor vs Non-Progressor)")), 
       title = "Significant Pathways (FDR < 0.05)") +
  theme_minimal(base_size = 16) +  # Set a larger base font size
  theme(axis.text.y = element_text(size = 14, face = "bold"),  # Improve y-axis labels
        axis.title = element_text(size = 16, face = "bold"),  # Bold the axis titles
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Center and bold the title
        axis.text.x = element_text(size = 14),  # Improve x-axis labels
        panel.grid = element_blank(),  # Remove gridlines
        axis.ticks = element_line(size = 0.5),  # Adjust axis tick size
        plot.margin = margin(20, 20, 20, 20))  # Add some padding around the plot

