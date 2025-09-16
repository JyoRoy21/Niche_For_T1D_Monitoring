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

## Progresssor Vs Non-Progressor ----
### ----DotPlot Pathways-----------------------------------------------------------------------------------------------------------------------------------

setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results")
Macrophage_metab_pathway <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/pathway_results_Macrophage_ProgressorVsNonProgressor.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
Macrophage_metab_pathway <- na.omit(Macrophage_metab_pathway)

Macrophage_metab_pathway
colnames(Macrophage_metab_pathway)

library(ggplot2)
library(RColorBrewer)


# Rename the second column to "filtered_Pathway"
colnames(Macrophage_metab_pathway)[1] <- "filtered_Pathway"

# Filter pathways with Raw p < 0.05 and Impact > 0
filtered_pathways <- Macrophage_metab_pathway[Macrophage_metab_pathway$`Raw p` < 0.05 & Macrophage_metab_pathway$Impact > 0, ]


# Define custom color limits to start at a slightly higher shade
min_logP <- min(-log10(filtered_pathways$`Raw p`))
max_logP <- max(-log10(filtered_pathways$`Raw p`))

ggplot(filtered_pathways, aes(x = -log10(FDR), y = reorder(filtered_Pathway, -log10(FDR)), size = Impact, color = -log10(`Raw p`))) +
  geom_point(alpha = 0.8) +
  
  # Adjust color scale: shift the lower bound to a higher shade
  scale_color_gradientn(
    colors = brewer.pal(9, "YlOrRd"),
    name = expression(bold(-log[10](P))),  # Make legend title bold
    limits = c(min_logP * 0.5, max_logP),  
    oob = scales::squish  
  ) +
  
  # Adjust dot size scale
  scale_size_continuous(
    name = "Impact",
    range = c(6, 12)
  ) +
  
  # Improve plot aesthetics
  theme_minimal(base_size = 22) +
  labs(title = "Macrophage Metabolomics",
       x = expression(-log[10](P[adj])),
       y = "Pathways") +
  
  theme(
    axis.text.y = element_text(size = 22, family = "Arial"),
    axis.text.x = element_text(size = 18, family = "Arial"),
    axis.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )

### ----Raindrop plot-----------------------------------------------------------------------------------------------------------------------------------


# Load Metabolite-Pathway Mapping
pathway_metabolite_list <- read_excel("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/Macrophage_SignificantPathways_Metabolites.xlsx", sheet = 1)

pathway_metabolite_long <- pathway_metabolite_list %>%
  separate_rows(Metabolite, sep = ";") %>%
  mutate(Metabolite = trimws(Metabolite))  # Remove extra spaces


Macrophage_metabs <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/volcano_Macrophage_ProgressorVsNonProgressor.csv", 
                                  sep=",", header=TRUE, check.names=FALSE)
colnames(Macrophage_metabs)[1] <- "Metabolite"
colnames(Macrophage_metabs)[5] <- "-log10(padj)"


# Merge with pathway annotations
all_metabs_with_pathways <- pathway_metabolite_long %>%
  inner_join(Macrophage_metabs, by = "Metabolite") 

# Adjust x-axis positions for pathways
all_metabs_with_pathways <- all_metabs_with_pathways %>%
  mutate(Pathway = factor(Pathway, levels = unique(Pathway)),
         x_pos = as.numeric(Pathway) * 10)  # Spread out pathways

# Define breaks and labels for x-axis
breaks <- seq(10, max(all_metabs_with_pathways$x_pos), by = 10)
labels <- levels(all_metabs_with_pathways$Pathway)

# Load RColorBrewer
library(RColorBrewer)

# Custom color gradient without white (based on reversed RdBu palette)
custom_colors <- rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#92C5DE", "#4393C3", "#2166AC", "#053061"))

# Raindrop Plot with Custom Red-Blue Color Transition
ggplot(all_metabs_with_pathways, aes(x = x_pos, y = `log2(FC)`, size = `-log10(padj)`, color = abs(`log2(FC)`))) +
  geom_point(alpha = 0.8, position = position_jitter(width = 3, height = 0)) +  # Increased jitter width
  
  # Reverse gradient with no white
  scale_color_gradientn(colors = custom_colors, name = bquote(bold("|log"[2]*"(FC)|"))) +
  scale_size(range = c(4, 8), name = bquote(bold("-log10(" * p[adj] * ")"))) +
  
  # Add vertical separators between pathways
  geom_vline(xintercept = seq(15, max(all_metabs_with_pathways$x_pos), by = 10), linetype = "solid", color = "grey50") +
  
  # Horizontal reference line at log2(FC) = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 1.2) +
  
  # Custom x-axis labels
  scale_x_continuous(breaks = breaks, labels = labels) +
  
  # Improve plot theme
  theme_minimal(base_size = 14) +
  labs(
    x = "Pathways",
    y = expression(log[2]("FC Progressor Vs Non-Progressor"))
  ) +
  
  theme(
    axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "Arial"),
    axis.text.y = element_text(size = 16, family = "Arial"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

### ---Heatmap ----------------------------------------------------------------------------------------------------------------------------------
# Load necessary libraries
library(pheatmap)
library(RColorBrewer)

# Read the dataset, skipping the first row (metadata row)
Macrophage_Normlaized_metabolities <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/data_normalized_Macrophage_ProgressorVsNonProgressor.csv", sep=",", header=T, check.names = FALSE))

# Remove the first row which contains metadata, and store the metabolite names in a separate object
metadata <- Macrophage_Normlaized_metabolities[1, ]  # Store the metadata
Macrophage_Normlaized_metabolities <- Macrophage_Normlaized_metabolities[-1, ]  # Remove the metadata row

# Set the first column (metabolite names) as rownames
rownames(Macrophage_Normlaized_metabolities) <- Macrophage_Normlaized_metabolities[, 1]
Macrophage_Normlaized_metabolities <- Macrophage_Normlaized_metabolities[, -1]  # Remove the first column which is metabolite names

# Convert the data to numeric (force numeric conversion for the entire dataset)
Macrophage_Normlaized_metabolities[] <- lapply(Macrophage_Normlaized_metabolities, function(x) as.numeric(as.character(x)))

# Check for any NA values created by non-numeric values and replace them (optional)
Macrophage_Normlaized_metabolities[is.na(Macrophage_Normlaized_metabolities)] <- 0  # Replace NA values with 0

# Define the groups for metabolites (matching with your previous categorization)
metabolite_groups <- c(
  "Melibiose" = "Carbohydrates/Sugars",
  "D-Mannose" = "Carbohydrates/Sugars",
  "myo-Inositol" = "Carbohydrates/Sugars",
  "GlcNAc" = "Carbohydrates/Sugars",
  "D-Fructose 1,6-biphosphate" = "Carbohydrates/Sugars",
  
  "L-Aspartic Acid" = "Amino Acids and Derivatives",
  "N-acetylaspartate" = "Amino Acids and Derivatives",
  "L-asparagine" = "Amino Acids and Derivatives",
  "L-Glutamic acid" = "Amino Acids and Derivatives",
  "L-Glutamine" = "Amino Acids and Derivatives",
  "L-Tyrosine" = "Amino Acids and Derivatives",
  "L-Histidine" = "Amino Acids and Derivatives",
  
  "Uridine 5-diphosphogalactose" = "Nucleotides and Derivatives",
  "beta-Nicotinamide adenine dinucleotide" = "Nucleotides and Derivatives",
  "Nicotinic acid" = "Nucleotides and Derivatives",
  "Quinolinate" = "Nucleotides and Derivatives",
  
  "Phosphoenolpyruvic acid" = "Organic Acids and TCA Cycle Intermediates",
  "Succinic acid" = "Organic Acids and TCA Cycle Intermediates",
  "L-Malic acid" = "Organic Acids and TCA Cycle Intermediates",
  "Citric acid" = "Organic Acids and TCA Cycle Intermediates",
  "Pyruvic acid" = "Organic Acids and TCA Cycle Intermediates"
)

# Get the metabolites in the dataset
metabolites_in_data <- rownames(Macrophage_Normlaized_metabolities)

# Ensure the metabolite names match
matching_metabolites <- metabolites_in_data %in% names(metabolite_groups)
matched_metabolites <- metabolites_in_data[matching_metabolites]

# Subset data with only the matching metabolites
Macrophage_Normlaized_metabolities <- Macrophage_Normlaized_metabolities[matching_metabolites, ]

# Assign group labels for each metabolite
group_labels <- metabolite_groups[matched_metabolites]

# Sort the metabolites by their group
ordered_metabolites <- matched_metabolites[order(factor(group_labels, levels = unique(group_labels)))]

# Reorder the data based on the metabolite group
Macrophage_Normlaized_metabolities <- Macrophage_Normlaized_metabolities[ordered_metabolites, ]
# Row-wise Z-normalization
Macrophage_Normlaized_metabolities <- t(scale(t(Macrophage_Normlaized_metabolities)))

# Extract metadata for coloring the samples (Progressor vs Non-Progressor)
sample_metadata <- colnames(Macrophage_Normlaized_metabolities)
sample_group <- ifelse(grepl("NP", sample_metadata), "Non-Progressor", "Progressor")

# Create a color palette for the sample groups
sample_group_colors <- c("Progressor" = "blue", "Non-Progressor" = "red")

# Define a new color palette for the heatmap (avoiding yellow)
heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

# Plot heatmap using pheatmap
pheatmap(Macrophage_Normlaized_metabolities,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = data.frame(Group = sample_group, row.names = colnames(Macrophage_Normlaized_metabolities)),
         annotation_colors = list(Group = sample_group_colors),
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = heatmap_colors,
         fontsize = 12,
         fontsize_row = 14,
         fontsize_col = 14,  # Increased font size for x-axis labels
         main = "Macrophage Metabolomics Heatmap")


## M1-M2-M0 Vs Progresssor-Non-Progressor ----
### ----DotPlot Pathways-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results")
Macrophage_PVsNP <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/ProgressorVsNonProgressor/pathway_results_Macrophage_ProgressorVsNonProgressor.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
Macrophage_PVsNP  <- na.omit(Macrophage_PVsNP)

Macrophage_M1VsM0 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/M1VsM0/pathway_results_M1VsM0.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
Macrophage_M1VsM0 <- na.omit(Macrophage_M1VsM0)

Macrophage_M2VsM0 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/M2VsM0/pathway_results_M2VsM0.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
Macrophage_M2VsM0 <- na.omit(Macrophage_M2VsM0)


# Add a timepoint column to each dataset
Macrophage_PVsNP$Type <- "PVsNP"
Macrophage_M1VsM0$Type <- "M1VsM0"
Macrophage_M2VsM0$Type <- "M2VsM0"

# Bind all data together
all_metab_pathways <- bind_rows(Macrophage_PVsNP, Macrophage_M1VsM0, Macrophage_M2VsM0)

# Rename first column as "Pathway" if it's unnamed
colnames(all_metab_pathways)[1] <- "Pathway"

filtered_pathways <- all_metab_pathways %>%
  filter(Pathway %in% c(
    "Galactose metabolism",
    "Fructose and mannose metabolism",
    "Amino sugar and nucleotide sugar metabolism",
    "Nicotinate and nicotinamide metabolism",
    "Phenylalanine, tyrosine and tryptophan biosynthesis",
    "Glycolysis or Gluconeogenesis",
    "Citrate cycle (TCA cycle)",
    "Alanine, aspartate and glutamate metabolism",
    "Histidine metabolism"
  ))

unique(filtered_pathways$Pathway)

# Calculate -log10(FDR) for Early timepoint
PVsNP_fdr <- filtered_pathways %>%
  filter(Type == "PVsNP") %>%
  mutate(log10_FDR = -log10(FDR))

# Reorder Pathway factor levels based on descending log10_FDR
filtered_pathways <- filtered_pathways %>%
  mutate(Pathway = factor(Pathway, levels = PVsNP_fdr$Pathway[order(PVsNP_fdr$log10_FDR, decreasing = FALSE)]))

# Assuming 'filtered_pathways' is your dataframe
filtered_pathways <- filtered_pathways %>%
  group_by(Pathway) %>%
  ungroup()
filtered_pathways <- filtered_pathways %>%
  mutate(log10_FDR = -log10(FDR))
# Create the dotplot
# Plot
# PLOT
filtered_pathways$log10_FDR

filtered_pathways
ggplot(filtered_pathways, aes(x = log10_FDR, y = Pathway, size = Impact)) +
  
  #PVsNP timepoint (Yellow to Red)
  geom_point(data = subset(filtered_pathways, Type == "PVsNP"), aes(color = -log10(`Raw p`)), alpha = 0.8) +
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "Progressor Vs Non-Progressor") +
  new_scale_color() +
  
  # Intermediate timepoint (Reversed Viridis: Light to Deep Green)
  geom_point(data = subset(filtered_pathways, Type == "M1VsM0"), aes(color = -log10(`Raw p`)), alpha = 0.8) +
  scale_color_gradient(low = "#D5E21A", high = "#008000", name = "M1 Vs M0") +  # Light to deep green
  new_scale_color() +
  
  # Late timepoint (Blues)
  geom_point(data = subset(filtered_pathways, Type == "M2VsM0"), aes(color = -log10(`Raw p`)), alpha = 0.8) +
  scale_color_gradient(low = "#BFD3E6", high = "#08306B", name = "M2 Vs M0") +  # Light to deep blue
  
  # Dot size scale with custom breaks and labels
  scale_size_continuous(
    name = "Impact",
    breaks = c(0.2, 0.3, 0.4, 0.5),
    labels = c("≤0.2", "0.3", "0.4", "0.5"),
    range = c(2, 8)
  ) +
  
  # Improve plot theme
  theme_minimal(base_size = 14) +
  labs(title = "Macrophage  Metabolomics Pathway Enrichment",
       x = expression(-log[10](P[adj])),
       y = "Pathways") +
  
  theme(
    axis.text.y = element_text(size = 16, family = "Arial"),
    axis.text.x = element_text(size = 16, family = "Arial"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_blank(),  # Remove only x-axis grid
    panel.grid.major.y = element_line(color = "grey80"),  # Keep y-axis grid
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )

### ----Raindrop plot-----------------------------------------------------------------------------------------------------------------------------------

#setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Results/Scaffold")

# Load Metabolite-Pathway Mapping
pathway_metabolite_list <- read_excel("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/ProgressorVsNonProgressor/Macrophage_SignificantPathways_Metabolites.xlsx", sheet = 1)

pathway_metabolite_long <- pathway_metabolite_list %>%
  separate_rows(Metabolite, sep = ";") %>%
  mutate(Metabolite = trimws(Metabolite))  # Remove extra spaces


# Load Metabolite Data for Each Type
Macrophage_metabs_PVsNP <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/ProgressorVsNonProgressor/volcano_Macrophage_ProgressorVsNonProgressor.csv", 
                                sep=",", header=TRUE, check.names=FALSE)
colnames(Macrophage_metabs_PVsNP)[1] <- "Metabolite"
colnames(Macrophage_metabs_PVsNP)[5] <- "-log10(padj)"
Macrophage_metabs_PVsNP <- Macrophage_metabs_PVsNP %>% 
  na.omit() %>% 
  mutate(Type = "PVsNP")


Macrophage_metabs_M1VsM0 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/M1VsM0/volcano_M1VsM0.csv", 
                                      sep=",", header=TRUE, check.names=FALSE)
colnames(Macrophage_metabs_M1VsM0)[1] <- "Metabolite"
colnames(Macrophage_metabs_M1VsM0)[5] <- "-log10(padj)"
Macrophage_metabs_M1VsM0 <- Macrophage_metabs_M1VsM0 %>% 
  na.omit() %>% 
  mutate(Type = "M1VsM0")


Macrophage_metabs_M2VsM0 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Metabolomics/Macrophage/Results/M2VsM0/volcano_M2VsM0.csv", 
                                       sep=",", header=TRUE, check.names=FALSE)
colnames(Macrophage_metabs_M2VsM0)[1] <- "Metabolite"
colnames(Macrophage_metabs_M2VsM0)[5] <- "-log10(padj)"
Macrophage_metabs_M2VsM0 <- Macrophage_metabs_M2VsM0 %>% 
  na.omit() %>% 
  mutate(Type = "M2VsM0")


# Identify common metabolites across all three comparisons
common_metabolites <- Reduce(intersect, list(
  Macrophage_metabs_PVsNP$Metabolite, 
  Macrophage_metabs_M1VsM0$Metabolite, 
  Macrophage_metabs_M2VsM0$Metabolite
))
all_metabolites <- unique(c(
  Macrophage_metabs_PVsNP$Metabolite, 
  Macrophage_metabs_M1VsM0$Metabolite, 
  Macrophage_metabs_M2VsM0$Metabolite
))
# Identify metabolites that are NOT common to all three
non_common_metabolites <- setdiff(all_metabolites, common_metabolites)

# Print non-common metabolites
print(non_common_metabolites)

# Filter each dataset to retain only the common metabolites
Macrophage_metabs_PVsNP <- Macrophage_metabs_PVsNP %>% filter(Metabolite %in% common_metabolites)
Macrophage_metabs_M1VsM0 <- Macrophage_metabs_M1VsM0 %>% filter(Metabolite %in% common_metabolites)
Macrophage_metabs_M2VsM0 <- Macrophage_metabs_M2VsM0 %>% filter(Metabolite %in% common_metabolites)

# Combine all datasets
all_metabs <- bind_rows(Macrophage_metabs_PVsNP, Macrophage_metabs_M1VsM0, Macrophage_metabs_M2VsM0)

# Check column names
colnames(all_metabs)

# Combine all datasets
all_metabs <- bind_rows(Macrophage_metabs_PVsNP, Macrophage_metabs_M1VsM0, Macrophage_metabs_M2VsM0)
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
         Type= factor(Type, levels = c("PVsNP", "M1VsM0", "M2VsM0")),
         x_pos = as.numeric(Pathway) * 10 + as.numeric(Type) * 2 - 4)  # Centering timepoints within each pathway

# Define breaks and labels
breaks <- seq(10, max(all_metabs_with_pathways$x_pos), by = 10)
labels <- levels(all_metabs_with_pathways$Pathway)


# Plot
ggplot(all_metabs_with_pathways, aes(x = x_pos, y = `log2(FC)`, size = -log10(p.ajusted))) +
  
  # Early Type (Yellow to Red)
  geom_point(data = subset(all_metabs_with_pathways, Type == "PVsNP"), aes(color = abs(`log2(FC)`)), alpha = 0.8) +
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "Progressor Vs Non-Progressor") +
  new_scale_color() +
  
  # Intermediate Type (Reversed Green Gradient)
  geom_point(data = subset(all_metabs_with_pathways, Type == "M1VsM0"), aes(color = abs(`log2(FC)`)), alpha = 0.8) +
  scale_color_gradient(low = "#D5E21A", high = "#008000", name = "M1 Vs M0") +
  new_scale_color() +
  
  # Late Type (Blue Shades)
  geom_point(data = subset(all_metabs_with_pathways, Type == "M2VsM0"), aes(color = abs(`log2(FC)`)), alpha = 0.8) +
  scale_color_gradient(low = "#BFD3E6", high = "#08306B", name = "M2 Vs M0") +
  
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
       y = expression(log[2]("Fold Change"))) +
  
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



