###Tanimotos_all_by_all###
#author:haim krupkin
#date:04/24/2025
##description
###this script is insteede to get the tanimotot claclatuoin of all vs all###


###
#BiocManager::install("ChemmineR")
#BiocManager::install("ChemmineOB")

library(jsonlite)
library(dplyr)
library(tidyr)
library(ChemmineOB)
library(ChemmineR)
library(ggplot2)
###loading_ligand_defintions####
ligand_binders<-read.csv("C:/Users/haim_krupkin/Downloads/only_significant_biolip.csv")

###we are going to check for ligands for whom co occurance occured###

over10_occruances_ligands<-read.table("C:/Users/haim_krupkin/Desktop/phd/AltmanRotation/ligands_with_over_10_cluster_associated_with_them.tsv")

#####
ligand_binders<-ligand_binders%>%
  dplyr::filter(k>4)%>%
  dplyr::mutate(percent_of_residues=k/N)
ligand_binders<-ligand_binders%>%
  dplyr::filter(percent_of_residues<1)
ions <- c(
  "CA",  # Calcium
  "ZN",  # Zinc
  "MG",  # Magnesium
  "MN",  # Manganese
  "FE",  # Iron
  "CU",  # Copper
  "CO",  # Cobalt
  "CL",  # Chloride
  "NI",  # Nickel
  "AG",  # Silver
  "AU",  # Gold
  "HG",  # Mercury
  "PB",  # Lead
  "BA",  # Barium
  "SR",  # Strontium
  "TL",  # Thallium
  "LI",  # Lithium
  "K",   # Potassium
  "NA",  # Sodium
  "H",   # Hydrogen (in some contexts)
  "OH",  # Hydroxide (not a metal ion, but often included)
  "BR",  # Bromide
  "I",   # Iodide
  "F",    # Fluoride
  "XE", #Xenon
  "CD", #"Cadmium"
  "OS" #osmium
)
ligand_binders<-ligand_binders%>%
  mutate(Molecule_Class = case_when(
    ligand == "peptide" ~ "peptide",
    ligand == "dna" ~ "DNA",
    ligand == "rna" ~ "RNA",
    (ligand) %in% ions~ "Ion",
    TRUE ~ "Small_molecule"  # Default case if none of the conditions are met
  ))
#
q_biolip <- fromJSON("C:/Users/haim_krupkin/Downloads/ligand_defintions.json")
####merging###
ligand_binders_with_ids<-ligand_binders%>%
  merge(q_biolip,
        by.x="ligand",
        by.y="ligid",
        all.x=T)

###
smiles_and_definitions<-ligand_binders_with_ids%>%
  dplyr::filter(Molecule_Class=="Small_molecule")%>%
  dplyr::mutate(smiles=as.character(smiles))%>%
  dplyr::select(ligand,smiles)%>%
  distinct()
smiles_and_definitions<-smiles_and_definitions%>%
  dplyr::filter(smiles!="NA" & !is.na(smiles))%>%
  dplyr::filter(ligand %in% over10_occruances_ligands$x)
head(smiles_and_definitions)
nrow(smiles_and_definitions)
###creating small test df###
smiles_and_definitions<-smiles_and_definitions[,]




# Loop through each SMILES entry and convert to SDF
#smiles_and_definitions$sdf_compound <- 
smiles_and_definitions$sdf_compound<-lapply(smiles_and_definitions$smiles, function(smiles_value) {
  tryCatch(
    {
      ChemmineR::smiles2sdf(smiles_value)  # Attempt conversion
    },
    error = function(e) {
      return("invalid_smiles")  # Return "invalid_smiles" if conversion fails
    }
  )
})

smiles_and_definitions$ap_compound<-lapply(smiles_and_definitions$sdf_compound, function(sdf_compound) {
  tryCatch(
    {
      ChemmineR::sdf2ap(sdf_compound)  # Attempt conversion
    },
    error = function(e) {
      return("invalid_sdf")  # Return "invalid_smiles" if conversion fails
    }
  )
})

smiles_and_definitions_backup<-smiles_and_definitions

smiles_and_definitions<-smiles_and_definitions%>%
  dplyr::filter(ap_compound!="invalid_sdf")

nrow(smiles_and_definitions)

df_all_tanimoto_results<-data.frame()
running_big_index=0
for (row_index_ligand1 in 1:nrow(smiles_and_definitions)){
  print(row_index_ligand1)
  print(running_big_index)
  smiles_and_definitions_row_ligand1<-smiles_and_definitions[row_index_ligand1,]
  print(smiles_and_definitions_row_ligand1)
  for (row_index_ligand2 in 1:nrow(smiles_and_definitions)){
    running_big_index=running_big_index+1
    #print(row_index_ligand2)
    smiles_and_definitions_row_ligand2<-smiles_and_definitions[row_index_ligand2,]
    #print(smiles_and_definitions_row_ligand1)
    tanimoto_similaritiy<-cmp.similarity(smiles_and_definitions_row_ligand1$ap_compound[[1]],
                   smiles_and_definitions_row_ligand2$ap_compound[[1]])
    #print(tanimoto_similaritiy)
    df_tanimoto_result<-data.frame(ligand1=smiles_and_definitions_row_ligand1$ligand,
                                   ligand2=smiles_and_definitions_row_ligand2$ligand,
                                   tanimoto_similaritiy=tanimoto_similaritiy)
    df_all_tanimoto_results<-rbind(df_all_tanimoto_results,df_tanimoto_result)
    
  }
  }

#df_all_tanimoto_results
#View(df_all_tanimoto_results)
df_all_tanimoto_results%>%
  ggplot(aes(x=ligand1,y=ligand2))+
  geom_tile(aes(fill=tanimoto_similaritiy))+
  labs(x = "Ligand", y = "Ligand", fill = "Tanimoto Similarity") +  # Set labels
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.key.size = unit(1.5, "cm"),
    legend.position = "top"  # Place legend at the top
  )
###
write.csv(df_all_tanimoto_results,"C:/Users/haim_krupkin/Desktop/phd/AltmanRotation/df_all_tanimoto_results.csv")
###
similarity_matrix <- df_all_tanimoto_results %>%
  pivot_wider(names_from = ligand2, values_from = tanimoto_similaritiy)

# Convert to matrix and set row names
ligand_names_row<-similarity_matrix$ligand1
rownames(similarity_matrix) <- similarity_matrix$ligand1
similarity_matrix <- as.matrix(similarity_matrix[, -1])  # Remove the ligand1 column
rownames(similarity_matrix)<-ligand_names_row
# Perform hierarchical clustering
dist_matrix <- dist(1 - similarity_matrix, method = "euclidean")  # Use 1 - similarity for distance
hc <- hclust(dist_matrix)

# Reorder the matrix based on clustering
ordered_labels <- hc$labels[hc$order]

# Create the heatmap with ggplot
p <- ggplot(df_all_tanimoto_results, aes(x = factor(ligand1, levels = ordered_labels), 
                            y = factor(ligand2, levels = ordered_labels))) +
  geom_tile(aes(fill = tanimoto_similaritiy)) +
  labs(x = "Ligand", y = "Ligand", fill = "Tanimoto Similarity") +  # Set labels
  scale_fill_gradient(low = "white", high = "red") +  # White to red color scale
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12)#,  # Increase legend text size
    #legend.key.size = unit(0.4, "cm"),
    #legend.justification = "center",
    #legend.position = "top"  # Place legend at the top
  )
p
ggsave(plot=p,
       filename = "C:/Users/haim_krupkin/Desktop/phd/AltmanRotation/tanimoto_similaritiy_similar_ligands_ligands_with_over_10_cluster_associated.png",
       dpi=900,
       width=7,
       height=5)
library(ggpubr)
no_identical_molecule_with_paris_no_repetitive<-read.csv("C:/Users/haim_krupkin/Desktop/phd/AltmanRotation/no_identical_molecule_with_paris_no_repetitive.csv")
data_for_ploting<-no_identical_molecule_with_paris_no_repetitive%>%
  mutate(similarity_category = case_when(
    tanimoto_similaritiy < 0.4 ~ "Low Similarity",
    tanimoto_similaritiy >= 0.4 & tanimoto_similaritiy < 0.6 ~ "Moderate Similarity",
    tanimoto_similaritiy >= 0.6 & tanimoto_similaritiy < 0.85 ~ "High Similarity",
    tanimoto_similaritiy >= 0.85 ~ "Very High Similarity",
    TRUE ~ "Undefined"  # This handles any NA or unexpected values
  ))%>%
  mutate(similarity_category = factor(similarity_category, 
                                      levels = c("Very Low Similarity", 
                                                 "Low Similarity", 
                                                 "Moderate Similarity", 
                                                 "High Similarity", 
                                                 "Very High Similarity")))

install.packages("DescTools")
library(DescTools)

# Assuming your data frame is named 'data_for_ploting'
# Perform the Jonckheere-Terpstra test
jt_test <- JonckheereTerpstraTest(normalized_co_occurrence ~ similarity_category,
                                  data = data_for_ploting)

# Extract the p-value
p_value <- jt_test$p.value

# Create the plot
boxplotingv2<-ggplot(data_for_ploting, aes(x = similarity_category,
                                           y = normalized_co_occurrence,
                                     fill = similarity_category)) +
  geom_boxplot(outlier.shape = 1, alpha = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("Very Low Similarity" = "#FFCCCC",  # Light red
                               "Low Similarity" = "#FF9999",   # Medium light red
                               "Moderate Similarity" = "#FF6666", # Medium red
                               "High Similarity" = "#FF3333",  # Darker red
                               "Very High Similarity" = "#FF0000", # Bright red
                               "Undefined" = "#FFCCCC")) +
  labs(x = "Tanimoto Similarity Category Ligand 1 and Ligand 2", 
       y = "Co Occurrence ligand 1 and Ligand 2") +
  theme(axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.title.y = element_text(size = 16),
        legend.position = "none", 
        legend.justification = "none",
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black")) +
  # Add the p-value to the plot
  annotate("text", x = 1, y = max(data_for_ploting$normalized_co_occurrence, na.rm = TRUE) * 0.9, 
           label = paste0("Jonckheere-Terpstra test P.Value<",
                          ifelse(p_value==0,
                          "2.2e-16",
                          p_value)),
                         size = 6,
           hjust = 0)
boxplotingv2
ggsave(plot=boxplotingv2,
       filename = "C:/Users/haim_krupkin/Desktop/phd/AltmanRotation/tanimoto_vs_cooocuance_groupsV2.png",
       dpi=900,
       width=8,
       height=6)
