###Analysing_AfM_enrichment_Analysis###
#author: Haim Krupkin
#date: 05/05/2025

###description####
#this code is intended to naalyse the results of the AfM pathogenicity enrichemtn analysis



###packages#####
#install.packages("reticulate",lib ="/home/users/hkrupkin/R_packages")
library(reticulate,lib.loc="/home/users/hkrupkin/R_packages")
reticulate::py_install('pandas', pip = TRUE)
#install.packages("tidyverse",lib ="/home/users/hkrupkin/R_packages")
library(dplyr,lib.loc="/home/users/hkrupkin/R_packages")
library(ggplot2,lib.loc="/home/users/hkrupkin/R_packages")
library(ggrepel,lib.loc="/home/users/hkrupkin/R_packages") 
library(labeling,lib.loc="/home/users/hkrupkin/R_packages")
library(farver,lib.loc="/home/users/hkrupkin/R_packages")
library(scales,lib.loc="/home/users/hkrupkin/R_packages")  # Load the scales package for formatting numbers
#install.packages("ggridges",lib="/home/users/hkrupkin/R_packages")
library(ggridges,lib.loc="/home/users/hkrupkin/R_packages")
#install.packages("patchwork",lib="/home/users/hkrupkin/R_packages")
library(patchwork,lib.loc="/home/users/hkrupkin/R_packages")
#install.packages("tidyr",lib="/home/users/hkrupkin/R_packages")
library(tidyr,lib.loc="/home/users/hkrupkin/R_packages")
#install.packages("data.table",lib="/home/users/hkrupkin/R_packages")
library(data.table,lib.loc="/home/users/hkrupkin/R_packages")
library(parallel)


#AfM<-read.delim("/oak/stanford/groups/rbaltman/alptartici/substrateSpec/sota_vep_data/AlphaMissense_aa_substitutions.tsv",
#                skip=3)
#table_of_pathogenicity<-table(AfM$am_class)

#proportion_pathogenic_alpha_missense_wide<-table_of_pathogenicity[[3]]/sum(table_of_pathogenicity)
  
  
#write.csv(proportion_pathogenic_alpha_missense_wide,
#            "/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/proportion_pathogenic_alpha_missense_wide.csv")
proportion_pathogenic_alpha_missense_wide<-read.csv("/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/proportion_pathogenic_alpha_missense_wide.csv")
proportion_pathogenic_alpha_missense_wide<-proportion_pathogenic_alpha_missense_wide$x
#0.4337479
###analysis###
#loading the data#
summarised_per_cluster_pathogeniciry_AfM_with_p.value<-read.csv(
          "/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/summarised_per_cluster_pathogeniciry_AfM_with_p_value.csv")
####
###plotting the results###
qqplot_clustering_enrichemtn_resutls<-ggplot(summarised_per_cluster_pathogeniciry_AfM_with_p.value,
                                             aes(x = -log10(p_expected), y = -log10(p_value))) +
  geom_point(color = "black", size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(
    title = "QQ Plot of Cluster Enriched or depleted p-values",
    x = "Expected -log10(p)",
    y = "Observed -log10(p)"
  ) +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 12),
    
    axis.ticks = element_line(color = "black")
  )+
  scale_y_continuous(labels = comma)
qqplot_clustering_enrichemtn_resutls
ggsave(plot=qqplot_clustering_enrichemtn_resutls,
       file="/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/qqplot_clustering_enrichemtn_resutls.png",
       dpi=900,
       width=5,
       height=5)

Odds_Ratio_plot<-ggplot(summarised_per_cluster_pathogeniciry_AfM_with_p.value,
                        aes(x = log2_odds_ratio,
                            y = -log10(p_value))) +
  geom_point(color = "black", size = 2, alpha = 0.7) +
  labs(
    title = "QQ Plot of Cluster Enrichment p-values",
    y = "-log10(p.value)",
    x = "log2(Odds Ratio)"
  ) +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 12),
    
    axis.ticks = element_line(color = "black")
  )+
  scale_y_continuous(labels = comma)
Odds_Ratio_plot
ggsave(plot=Odds_Ratio_plot,
       file="/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/Odds_Ratio_plot.png",
       dpi=900,
       width=5,
       height=5)
###finding log2odds ratio thresholds####
p0 <- 0.43
odds0 <- p0 / (1 - p0)

# Define the two new probabilities
p1_low <- 0.025
p1_high <- 0.975

# Compute odds for each
odds1_low <- p1_low / (1 - p1_low)
odds1_high <- p1_high / (1 - p1_high)

# Compute odds ratios
or_low <- odds1_low / odds0
or_high <- odds1_high / odds0

# Convert to log2 odds ratios
log2or_low <- log2(or_low)
log2or_high <- log2(or_high)

# Print
log2or_low
log2or_high

summarised_per_cluster_pathogeniciry_AfM_with_p.value<-summarised_per_cluster_pathogeniciry_AfM_with_p.value%>%
  dplyr::mutate(odds_ratio_thresholded=ifelse(log2_odds_ratio>log2or_high & p_value_adj<0.05,"High Enrichment",
                       ifelse(log2_odds_ratio<log2or_low &p_value_adj<0.05 ,
                              "High Depletion","Altered")))

color_palette <- c(
  "High Enrichment" = "red",   # Color for High Enrichment
  "Altered" = "grey",            # Color for Altered
  "High Depletion" = "blue"       # Color for High Depletion
)

Odds_Ratio_plot<-summarised_per_cluster_pathogeniciry_AfM_with_p.value%>%
  dplyr::filter(!is.na(odds_ratio_thresholded))%>%
   ggplot(aes(x = log2_odds_ratio,
                            y = -log10(p_value_adj))) +
  geom_point(aes(color = odds_ratio_thresholded), size = 2, alpha = 0.6) +
  labs(
    y = "-log10(P.Value Adjusted)",
    x = expression(log[2]~frac((frac("Pathogenic Mutation in Cluster","Non-Pathogenic Mutation in cluster")),
                               (frac("Pathogenic Mutation AlphaMissense Wide","Non-Pathogenic Mutation AlphaMissense Wide")))),
    color="AlphaMissense\nPathogenicity\nCategory"
  ) +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.ticks = element_line(color = "black"),
    legend.position = "none",
    legend.justification = "center",
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 12)   # Increase legend title size
  ) +
  scale_color_manual(values = color_palette) +  # Use the custom color palette
  scale_y_continuous(labels = comma)+
  geom_hline(yintercept =-log10(0.05) ,color="magenta",linetype="dashed")+
  geom_vline(xintercept =log2or_high ,color="yellow",linetype="dashed")+
  geom_vline(xintercept =log2or_low ,color="yellow",linetype="dashed")

Odds_Ratio_plot
ggsave(plot=Odds_Ratio_plot,
       file="/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/Volcano_cluster_Odds_Ratio_plot_alpha_missense.png",
       dpi=900,
       width=6,
       height=5)

table(summarised_per_cluster_pathogeniciry_AfM_with_p.value$odds_ratio_thresholded)
###Get number of residue that we might be able to predict using this approach ####
highly_predictive_clusters<-summarised_per_cluster_pathogeniciry_AfM_with_p.value%>%
  dplyr::filter(odds_ratio_thresholded %in% c("High Depletion","High Enrichment"))
highly_predictive_clusters<-highly_predictive_clusters%>%
  mutate(number_of_resdieus_for_real_in_cluster=n_cluster/19)
total_number_of_residues_we_can_predict_well<-sum(highly_predictive_clusters$n_cluster/19)
total_number_of_residues_in_human_proteome<-sum(summarised_per_cluster_pathogeniciry_AfM_with_p.value$n_cluster/19)
print(total_number_of_residues_we_can_predict_well/total_number_of_residues_in_human_proteome)

###comparing to humsavar and clinvar####
humsavar_clinvar<-read.csv("/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/merged_humsavar_and_clinvar_clusters_pathogenicity.csv")

most_pathogenic_clusters<-humsavar_clinvar%>%
  dplyr::filter((P.value_corrected_Clinvar <0.05 & Number_of_residues_in_cluster_pathogenic_Clinvar>4 )| (P.value_corrected_humsavar<0.05 & Number_of_residues_in_cluster_pathogenic_humsavar>4 ))


summarised_per_cluster_pathogeniciry_AfM_with_p.value<-summarised_per_cluster_pathogeniciry_AfM_with_p.value%>%
  mutate(enrichment_status=ifelse(proportion_pathogenic>0.43 & p_value_adj<0.05,
                "Enriched_significant","Not"))



humsavar_clinvar_data_AfM_data<-summarised_per_cluster_pathogeniciry_AfM_with_p.value%>%
  dplyr::filter(cluster_id %in% most_pathogenic_clusters$Cluster_ID)

###
label_data <- data.frame(
  previosuly_discovered = c("Identified Pathogenic Cluster by ClinVar/Humsavar", 
                            "Not Identified As Pathogenic Cluster by ClinVar/Humsavar"),
  label="AlphaMissense Wide",
  x = rep(proportion_pathogenic_alpha_missense_wide, 2),  # Same x position for both labels
  y = c(25, 1400)  # Adjust y positions as needed for each facet
)

plot_histograms<-summarised_per_cluster_pathogeniciry_AfM_with_p.value %>%
  mutate(previosuly_discovered = ifelse(cluster_id %in% most_pathogenic_clusters$Cluster_ID,
                                        "Identified Pathogenic Cluster by ClinVar/Humsavar",
                                        "Not Identified As Pathogenic Cluster by ClinVar/Humsavar")) %>%
  ggplot(aes(x = proportion_pathogenic, fill = previosuly_discovered)) +
  geom_histogram(binwidth = 0.02, alpha = 0.7, color = "black", position = "identity") +
  facet_wrap(~ previosuly_discovered, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(x = "Proportion Of Pathogenic Mutations in Cluster",
       y = "Number Of Clusters",
       fill = "Cluster Type",
       title = "") +
  theme(
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 20),
    axis.ticks = element_line(color = "black"),
    strip.text = element_text(size = 14, color = "black"),
    legend.position = "none"
  )+
  geom_vline(xintercept=proportion_pathogenic_alpha_missense_wide,
             linetype="dashed",
             color="black")+
  geom_text(data = label_data, 
            aes(label = label, x = x, y = y),
            size = 4,           # Adjust size as needed
            color = "black",     # Change color for better visibility
            fontface = "bold",  # Make the text bold
            hjust = 1,        # Center the text horizontally
            vjust = 0)
plot_histograms
ggsave(plot=plot_histograms,
       file="/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/histograms_humsavar_clinvar_alpha_missense_cluster_predictions.png",
       dpi=900,
       width=6,
       height=4)

###plotting things###
number_of_clusters_enriched_for_alpha_missense_variants<-table(summarised_per_cluster_pathogeniciry_AfM_with_p.value$enrichment_status)[[1]]
number_of_clusters_sampled<-length(unique(summarised_per_cluster_pathogeniciry_AfM_with_p.value$cluster_id))
proportion_pathogenic_AfM_clusters=number_of_clusters_enriched_for_alpha_missense_variants/number_of_clusters_sampled
print("high confidnece clusters: ")
print(table(summarised_per_cluster_pathogeniciry_AfM_with_p.value$odds_ratio_thresholded))

print("number of highly predictive pathogenic clsuters also found in humsavar/clinvar")
length(intersect(most_pathogenic_clusters$Cluster_ID,
                 highly_predictive_clusters[highly_predictive_clusters$odds_ratio_thresholded=="High Enrichment",]$cluster_id))
print("number of overlap of clinvar/humsavar in direction with alpha missense")
length(intersect(most_pathogenic_clusters$Cluster_ID,
                 summarised_per_cluster_pathogeniciry_AfM_with_p.value[summarised_per_cluster_pathogeniciry_AfM_with_p.value$proportion_pathogenic>0.43 &summarised_per_cluster_pathogeniciry_AfM_with_p.value$p_value_adj<0.05 ,]$cluster_id))

print("number_of_clusters_enriched_for_alpha_missense_variants")
print(number_of_clusters_enriched_for_alpha_missense_variants)
print("proportion_pathogenic_AfM_clusters")
print(proportion_pathogenic_AfM_clusters)

###Merged_ClinVar_Humsavar_AlphaFold####
###creating a merged dataframe of clinvar,humsavar, and alphafold
#head(humsavar_clinvar)
for_merge_humsavar_clinvar<-humsavar_clinvar%>%
  dplyr::select(-X,-Mutations_humsavar,
                -Number_of_residues_total_humsavar,-Number_of_pathogenic_residues_overall_humsavar,
                -Number_of_residues_in_cluster_humsavar,-Number_of_residues_in_cluster_pathogenic_humsavar,
                -Mutations_Clinvar,
                -Number_of_residues_total_Clinvar,-Number_of_pathogenic_residues_overall_Clinvar,
                -Number_of_residues_in_cluster_Clinvar,-Number_of_residues_in_cluster_pathogenic_Clinvar)
for_merge_summarised_per_cluster_pathogeniciry_AfM_with_p.value<-summarised_per_cluster_pathogeniciry_AfM_with_p.value%>%
  mutate(number_of_residues_in_cluster=cluster_size/19)%>%
  dplyr::select(2,number_of_residues_in_cluster,10,16,19,26,27,28,33,34,36)

merged_humsavar_clinvar_AfM<-merge(humsavar_clinvar,
      summarised_per_cluster_pathogeniciry_AfM_with_p.value,
      all=T,
      by.x="Cluster_ID",
      by.y="cluster_id")

###writing the merged dataframe###
write.csv(merged_humsavar_clinvar_AfM,
          "/home/users/hkrupkin/merged_humsavar_clinvar_AfM.csv")

###getting the clusters that are humsavar and clinvar associated
merged_humsavar_clinvar_AfM<-read.csv("/home/users/hkrupkin/merged_humsavar_clinvar_AfM.csv")
#############################################################################################
clusters_of_high_agreement_ClinVar_Humsavar_df<-merged_humsavar_clinvar_AfM%>%
  dplyr::filter((P.value_corrected_humsavar<0.05 & Number_of_residues_in_cluster_pathogenic_humsavar>4) & (P.value_corrected_Clinvar<0.05& Number_of_residues_in_cluster_pathogenic_Clinvar>4))
length(unique(clusters_of_high_agreement_ClinVar_Humsavar_df$Cluster_ID))

clusters_of_high_agreement_ClinVar_Humsavar<-unique(clusters_of_high_agreement_ClinVar_Humsavar_df$Cluster_ID)
###checking what proportion of proteins have resiudes rfom these clusters###
#loading cluster defintions###
long_df<-read.csv("/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/human_collapse_residues_and_clustering.csv")
head(long_df)
all_proteins_in_uniprot<-unique(long_df$uniprot)
length(all_proteins_in_uniprot)
all_uniprot_residues_from_pathogenic_clusters<-long_df%>%
  dplyr::filter(cluster_id %in% clusters_of_high_agreement_ClinVar_Humsavar)
nrow(all_uniprot_residues_from_pathogenic_clusters)

number_of_proteins_with_agreed_upon_pathogenic_clinvar_humsavar<-length(unique(all_uniprot_residues_from_pathogenic_clusters$uniprot))


number_of_proteins_in_AF_database<-21834
proportion_of_proteins_with_clusters_of_agreement_pathogenic_clinvar_humsavar<-number_of_proteins_with_agreed_upon_pathogenic_clinvar_humsavar*100/number_of_proteins_in_AF_database
print("percentage of proteins with agreed upon pathogenic clinvar_humsavar is : ")
print(proportion_of_proteins_with_clusters_of_agreement_pathogenic_clinvar_humsavar)
