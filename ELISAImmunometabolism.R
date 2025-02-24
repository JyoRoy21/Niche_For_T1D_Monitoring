# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/ELISA/")


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Import the data
data <- read.csv("ELISANODMousePlottingData.csv")

# Convert relevant columns to factors
data$Group <- factor(data$Group, levels = c("Non-Progressor", "Progressor"))
data$Time <- factor(data$Time, levels = c("Early", "Intermediate", "Late"))

# Keep '>' values as is, but track them for separate plotting
data$AverageConcentration <- as.character(data$AverageConcentration)

# Identify censored values and create a new column for censored
data$censored <- grepl(">", data$AverageConcentration)

# For statistical analysis, replace '>' values with NA (they won't be considered)
data$AverageConcentration[data$censored] <- NA

# Convert to numeric values (non-'>' values will convert, NAs will stay)
data$AverageConcentration <- as.numeric(data$AverageConcentration)

# Remove 'Late' timepoint for plotting (Optional)
data <- data %>% filter(Time != "Late")

# Perform t-tests for each timepoint (excluding censored values)
ttests <- data %>%
  group_by(Time) %>%
  summarise(p_value = t.test(AverageConcentration ~ Group, data = .)$p.value)

# Plotting with improved aesthetics
# Plotting with boxplots, individual data points, significance line, and y-axis range until 600
library(ggplot2)
library(ggsignif)

ggplot(data, aes(x = Group, y = AverageConcentration, fill = Group)) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.4, size = 1) + # Bold black border and no transparency for edges
  # Show individual data points (including censored ones)
  geom_point(aes(color = censored), shape = 16, size = 2, position = position_jitter(width = 0.1, height = 0)) + 
  scale_color_manual(values = c("black", "red")) + # Black for non-censored, red for censored
  facet_wrap(~ Time, scales = "free_y") + # Facet for each timepoint
  theme_classic() +  # Clean theme for publication quality
  labs(title = "Comparison of Concentrations Between Progressors and Non-Progressors", 
       y = "Average Concentration (Ng/ml)",
       x = "Group") +
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5), # Axis label formatting
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18, face = "bold"), # Facet label styling
        legend.position = "none", # Remove legend for cleaner look
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), # Title styling
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Border for the entire plot area
        ylim = c(0, 650)) + # Set y-axis range up to 600
  geom_signif(comparisons = list(c("Progressor", "Non-Progressor")), 
              map_signif_level = TRUE, 
              textsize = 5, 
              vjust = -0.5, # Adjust p-value label position
              tip_length = 0.03, # Length of the line above the p-value
              color = "black") # Color for the significance line
