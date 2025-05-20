###Ligand co occuance analysis####
#description: this code is meant to find co occurances of ligands
#author: haim krupkin
#date:04/24/2025



library(dplyr,lib.loc="/home/users/hkrupkin/R_packages")
library(parallel)
library(ggplot2,lib.loc="/home/users/hkrupkin/R_packages")
library(ggrepel,lib.loc="/home/users/hkrupkin/R_packages") 
library(labeling,lib.loc="/home/users/hkrupkin/R_packages")
library(farver,lib.loc="/home/users/hkrupkin/R_packages")
library(scales,lib.loc="/home/users/hkrupkin/R_packages")


ligand_binders<-read.csv("/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/only_significant_biolip.csv")

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



###UpsetR####
count_co_occurrences <- function(data, ligand1, ligand2) {
  ligand1_cluster_labels <- data %>%
    filter((ligand1==ligand)) %>% pull(cluster_label)
  ligand2_cluster_labels <- data %>%
    filter((ligand2==ligand)) %>% pull(cluster_label)
  return(length(intersect(ligand1_cluster_labels, ligand2_cluster_labels)))
}

# Print the results
list_of_ligands <- unique(ligand_binders$ligand)
data_frame_saving_combinations <- data.frame()
length(list_of_ligands)
running_index_ligand1 = 0

# Check for ligands with at least 2 appearances
ligand_counts <- table(ligand_binders$ligand)
valid_ligands <- names(ligand_counts[ligand_counts >= 10])
length(valid_ligands)
num_valid_ligands <- length(valid_ligands)

total_pairs <- (num_valid_ligands * (num_valid_ligands + 1)) / 2
cat("Total number of ligand pairs to be processed (including self-pairs):", total_pairs, "\n")

start_time <- Sys.time()
# Estimate time
estimated_time_per_pair <- 0.1  # Adjust this based on your previous experience
total_estimated_time <- total_pairs * estimated_time_per_pair
cat("Estimated time to complete (in seconds):", total_estimated_time, "\n")

# Start timing
start_time <- Sys.time()

# Set up parallel processing
num_cores <- detectCores() - 1  # Leave one core free
results_list <- mclapply(seq_along(valid_ligands), function(i) {
  ligand1 <- valid_ligands[i]
  results <- list()
  
  # Count ligand1-ligand1
  co_occurrence_ligand1_and_ligand1 <- count_co_occurrences(ligand_binders, ligand1, ligand1)
  results[[1]] <- data.frame(ligand1 = ligand1,
                             ligand2 = ligand1,
                             co_occurrence_ligand1_and_ligand2 = co_occurrence_ligand1_and_ligand1)
  
  for (j in (i + 1):length(valid_ligands)) {
    ligand2 <- valid_ligands[j]
    co_occurrence_ligand1_and_ligand2 <- count_co_occurrences(ligand_binders, ligand1, ligand2)
    results[[j - i]] <- data.frame(ligand1 = ligand1,
                                   ligand2 = ligand2,
                                   co_occurrence_ligand1_and_ligand2 = co_occurrence_ligand1_and_ligand2)
  }
  
  # Print progress
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  estimated_time_per_ligand <- elapsed_time / (i + 1)  # Average time per ligand
  remaining_ligands <- num_valid_ligands - (i + 1)
  estimated_time_left <- remaining_ligands * estimated_time_per_ligand
  
  cat(sprintf("Processed ligand %d of %d: %s | Estimated time left: %.2f seconds\n", 
              i + 1, num_valid_ligands, ligand1, estimated_time_left))
  
  return(do.call(rbind, results))
}, mc.cores = num_cores)

# Combine results
data_frame_saving_combinations <- do.call(rbind, results_list)


data_frame_saving_combinations_before_naing <- do.call(rbind, results_list)


# End timing
end_time <- Sys.time()
cat("Total time taken (in seconds):", as.numeric(difftime(end_time, start_time, units = "secs")), "\n")

write.csv(data_frame_saving_combinations,
          "/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/ligands_cooccuance_in_same_cluster.csv")

data_frame_saving_combinations<-read.csv("/oak/stanford/groups/rbaltman/hkrupkin/COLLAPSE_WORK/ligands_cooccuance_in_same_cluster.csv")

data_frame_saving_combinations$co_occurrence_ligand1_and_ligand2<-as.numeric(data_frame_saving_combinations$co_occurrence_ligand1_and_ligand2)

#data_frame_saving_combinations%>%
#  dplyr::filter(is.na(co_occurrence_ligand1_and_ligand2))%>%
#  View()

data_frame_saving_combinations%>%
  dplyr::filter(ligand1!=ligand2)%>%
  View()

data_frame_saving_combinations<-data_frame_saving_combinations%>%
  mutate(co_occurrence_ligand1_and_ligand2=as.numeric(co_occurrence_ligand1_and_ligand2))
####checkign specific examples

most_common_ligand1<-data_frame_saving_combinations%>%
  dplyr::arrange(co_occurrence_ligand1_and_ligand2)%>%
  top_n(70)%>%pull(ligand1)
most_common_ligand2<-data_frame_saving_combinations%>%
  dplyr::arrange(co_occurrence_ligand1_and_ligand2)%>%
  top_n(70)%>%pull(ligand2)

ligands_common<-c(most_common_ligand1,most_common_ligand2)
ligands_common
data_frame_saving_combinations%>%
  dplyr::filter(ligand1 %in%ligands_common & ligand2 %in% ligands_common)%>%
  View()

data_frame_saving_combinations%>%
  dplyr::filter(ligand1 %in%ligands_common & ligand2 %in% ligands_common)%>%
  ggplot(aes(x=ligand1 ,y=ligand2 ))+
  geom_tile(aes(fill=co_occurrence_ligand1_and_ligand2))+
  theme_bw()+
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "top",
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="Ligand",y="Ligand",fill="Number Of clusters Binding Same Ligand")
