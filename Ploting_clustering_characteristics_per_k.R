##Ploting_clustering_characteristics_per_k####
###author:haim krupkin#
###date: 04/14/2025
###description####


###packages###
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(scales)   # For number formatting

###analysis###
per_clustering_size<-read.csv("C:/Users/haim_krupkin/Downloads/pdb100_cluster_tuning(1).csv")
###
long_data <- per_clustering_size %>%
  pivot_longer(cols = -k, names_to = "Label", values_to = "Value") %>%
  mutate(Label = gsub("\\.", " ", Label))  # Replace '.' with ' '

# Create the plot
plot_of_intrest<-long_data%>%
  mutate(Label=ifelse(Label=="Within cluster sum of squares",
                      "Within-cluster sum of squares",
                      Label))%>%
  dplyr::filter(Label %in% c("Within-cluster sum of squares",
                             "Folds per cluster",
                             "SCOP entropy",
                             "Samples per cluster"))%>%
  mutate(Label = fct_relevel(Label, 
                             "Within-cluster sum of squares", 
                             "Folds per cluster", 
                             "SCOP entropy", 
                             "Samples per cluster"))%>%
  ggplot( aes(x = k, y = Value)) +
  geom_line(size = 1,color="cyan") +
  geom_point(size = 2,color="cyan") +
  facet_wrap(~ Label, scales = "free")+#, nrow = 6) +  # Adjust number of rows as needed
  labs(title = "", x = "k", y = "") +
  #scale_x_log10() +  # Log scale for x-axis
  scale_color_brewer(palette = "Set1") +  # Use a color palette similar to Seaborn
  theme_cowplot() +  # Minimal theme with base font size
  scale_y_continuous(labels = label_number()) +  # Format y-axis labels as numbers
  
  theme(
    strip.background = element_blank(),  # Hide strip background
    
    strip.text = element_text(size = 10),  # Title size for facets
    axis.text = element_text(size = 10),   # Axis text size
    axis.title = element_text(size = 12),  # Axis title size
    
    #panel.grid.major = element_line(color = "grey80"),  # Major grid lines
    #panel.grid.minor = element_blank()  # No minor grid lines
  )

plot_of_intrest
# Save the plot
ggsave(plot=plot_of_intrest,
       file="C:/Users/haim_krupkin/Desktop/phd/AltmanRotation/Plots_for_Collapse_Clustering_revision/Figure_S1.png",
       dpi = 900, width =  6, height = 4)

###modified plot figure S1
plot_of_intrest_v2<-long_data%>%
  mutate(Label=ifelse(Label=="Within cluster sum of squares",
                      "Within-cluster sum of squares",
                      Label))%>%
  mutate(Label=ifelse(Label=="Folds per cluster",
                      "Average Folds Per Cluster",
                      Label))%>%
  dplyr::filter(Label %in% c("Within-cluster sum of squares",
                             "Average Folds Per Cluster",
                             "SCOP entropy",
                             "Samples per cluster"))%>%
  mutate(Label = fct_relevel(Label, 
                             "Within-cluster sum of squares", 
                             "Average Folds Per Cluster", 
                             "SCOP entropy", 
                             "Samples per cluster"))%>%
  ggplot( aes(x = k, y = Value)) +
  geom_line(size = 1,color="navyblue") +
  geom_point(size = 2,color="navyblue") +
  facet_wrap(~ Label, scales = "free")+#, nrow = 6) +  # Adjust number of rows as needed
  labs(title = "", x = "k", y = "") +
  #scale_x_log10() +  # Log scale for x-axis
  scale_color_brewer(palette = "Set1") +  # Use a color palette similar to Seaborn
  theme_cowplot() +  # Minimal theme with base font size
  scale_y_continuous(trans = 'log10',labels = label_number()) +  # Format y-axis labels as numbers
  
  theme(
    strip.background = element_blank(),  # Hide strip background
    
    strip.text = element_text(size = 11),  # Title size for facets
    axis.text = element_text(size = 8),   # Axis text size
    axis.title.x = element_text(size = 16),   # Axis text size
    
    axis.title = element_text(size = 12),  # Axis title size
    
    #panel.grid.major = element_line(color = "grey80"),  # Major grid lines
    #panel.grid.minor = element_blank()  # No minor grid lines
  )

plot_of_intrest_v2
ggsave(plot=plot_of_intrest_v2,
       file="C:/Users/haim_krupkin/Desktop/phd/AltmanRotation/Plots_for_Collapse_Clustering_revision/Figure_S1_v3.png",
       dpi = 900, width =  6, height = 4)
