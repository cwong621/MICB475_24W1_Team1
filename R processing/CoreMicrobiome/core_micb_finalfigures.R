# Load libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library("sf")

# Load data
load("HIV Phyloseq/HIV_final.RData") #entire dataset
load("HIV Phyloseq/leftphyloseq.RData") #left cohort
load("HIV Phyloseq/rightphyloseq.RData") #right cohort

# See metadata
SAMP <- sample_data(HIV_final) #entire dataset
SAMP_left <- sample_data(ps_left) #left cohort
SAMP_right <- sample_data(ps_right) #right cohort

#convert phyloseq object to relative abundance
HIV_RA <- transform_sample_counts(HIV_final, fun=function(x) x/sum(x))
HIV_RA_left <- transform_sample_counts(ps_left, fun=function(x) x/sum(x))
HIV_RA_right <- transform_sample_counts(ps_right, fun=function(x) x/sum(x))

# Filter dataset (responsive, nonresponsive and healthy before and after (regardless of viral load))
### Note: before and after mean at visit 2 and 3, respectively

#responsive visit 2 (regardless of viral load) both cohorts
res_2_both <- subset_samples(HIV_RA, `response_patient_by_visit`=="responsive 2")
#responsive visit 3 (regardless of viral load) both cohorts
res_3_both <- subset_samples(HIV_RA, `response_patient_by_visit`=="responsive 3")
#responsive visit 2 (regardless of viral load) left cohort
res_2_left <- subset_samples(HIV_RA_left, `response_patient_by_visit`=="responsive 2")
#responsive visit 3 (regardless of viral load) left cohort
res_3_left <- subset_samples(HIV_RA_left, `response_patient_by_visit`=="responsive 3")
#responsive visit 2 (regardless of viral load) right cohort
res_2_right <- subset_samples(HIV_RA_right, `response_patient_by_visit`=="responsive 2")
#responsive visit 3 (regardless of viral load) right cohort
res_3_right <- subset_samples(HIV_RA_right, `response_patient_by_visit`=="responsive 3")

#nonresponsive visit 2 (regardless of viral load) both cohorts
nonres_2_both <- subset_samples(HIV_RA, `response_patient_by_visit`=="nonresponsive 2")
#nonresponsive visit 3 (regardless of viral load) both cohorts
nonres_3_both <- subset_samples(HIV_RA, `response_patient_by_visit`=="nonresponsive 3")
#nonresponsive visit 2 (regardless of viral load) left cohort
nonres_2_left <- subset_samples(HIV_RA_left, `response_patient_by_visit`=="nonresponsive 2")
#nonresponsive visit 3 (regardless of viral load) left cohort
nonres_3_left <- subset_samples(HIV_RA_left, `response_patient_by_visit`=="nonresponsive 3")
#nonresponsive visit 2 (regardless of viral load) right cohort
nonres_2_right <- subset_samples(HIV_RA_right, `response_patient_by_visit`=="nonresponsive 2")
#nonresponsive visit 3 (regardless of viral load) right cohort
nonres_3_right <- subset_samples(HIV_RA_right, `response_patient_by_visit`=="nonresponsive 3")


#Healthy at visit 2  both cohorts
healthy_2_both <- subset_samples(HIV_RA, `response_patient_by_visit`=="Healthy 2" & `start_viral_load_patient_by_visit`=="Healthy 2")
#Healthy at visit 3  both cohorts
healthy_3_both <- subset_samples(HIV_RA, `response_patient_by_visit`=="Healthy 3" & `start_viral_load_patient_by_visit`=="Healthy 3")
#Healthy at visit 2 left cohort
healthy_2_left <- subset_samples(HIV_RA_left, `response_patient_by_visit`=="Healthy 2" & `start_viral_load_patient_by_visit`=="Healthy 2")
#Healthy at visit 3 left cohort
healthy_3_left <- subset_samples(HIV_RA_left, `response_patient_by_visit`=="Healthy 3" & `start_viral_load_patient_by_visit`=="Healthy 3")
#Healthy at visit 2 right cohort
healthy_2_right <- subset_samples(HIV_RA_right, `response_patient_by_visit`=="Healthy 2" & `start_viral_load_patient_by_visit`=="Healthy 2")
#Healthy at visit 3 right cohort
healthy_3_right <- subset_samples(HIV_RA_right, `response_patient_by_visit`=="Healthy 3" & `start_viral_load_patient_by_visit`=="Healthy 3")



### Abundance and prevalence - Find ASVs for the different categories ###

#responsive visit 2 (regardless of viral load) both cohorts
res_2_ASVs_both <- core_members(res_2_both, detection=0.001, prevalence=0.5)
#responsive visit 3 (regardless of viral load) both cohorts
res_3_ASVs_both <- core_members(res_3_both, detection=0.001, prevalence=0.5)
#responsive visit 2 (regardless of viral load) left cohort
res_2_ASVs_left <- core_members(res_2_left, detection=0.001, prevalence=0.5)
#responsive visit 3 (regardless of viral load) left cohort
res_3_ASVs_left <- core_members(res_3_left, detection=0.001, prevalence=0.5)
#responsive visit 2 (regardless of viral load) right cohort
res_2_ASVs_right <- core_members(res_2_right, detection=0.001, prevalence=0.5)
#responsive visit 3 (regardless of viral load) right cohort
res_3_ASVs_right <- core_members(res_3_right, detection=0.001, prevalence=0.5)


#nonresponsive visit 2 (regardless of viral load) both cohorts
nonres_2_ASVs_both <- core_members(nonres_2_both, detection=0.001, prevalence=0.5)
#nonresponsive visit 3 (regardless of viral load) both cohorts
nonres_3_ASVs_both <- core_members(nonres_3_both, detection=0.001, prevalence=0.5)
#nonresponsive visit 2 (regardless of viral load) left cohort
nonres_2_ASVs_left <- core_members(nonres_2_left, detection=0.001, prevalence=0.5)
#nonresponsive visit 3 (regardless of viral load) left cohort
nonres_3_ASVs_left <- core_members(nonres_3_left, detection=0.001, prevalence=0.5)
#nonresponsive visit 2 (regardless of viral load) right cohort
nonres_2_ASVs_right <- core_members(nonres_2_right, detection=0.001, prevalence=0.5)
#nonresponsive visit 3 (regardless of viral load) right cohort
nonres_3_ASVs_right <- core_members(nonres_3_right, detection=0.001, prevalence=0.5)


#Healthy at visit 2 both cohorts
healthy_2_ASVs_both <- core_members(healthy_2_both, detection=0.001, prevalence=0.5)
#Healthy at visit 3 both cohorts
healthy_3_ASVs_both <- core_members(healthy_3_both, detection=0.001, prevalence=0.5)
#Healthy at visit 2 left cohort
healthy_2_ASVs_left <- core_members(healthy_2_left, detection=0.001, prevalence=0.5)
#Healthy at visit 3 left cohort
healthy_3_ASVs_left <- core_members(healthy_3_left, detection=0.001, prevalence=0.5)
#Healthy at visit 2 right cohort
healthy_2_ASVs_right <- core_members(healthy_2_right, detection=0.001, prevalence=0.5)
#Healthy at visit 3 right cohort
healthy_3_ASVs_right <- core_members(healthy_3_right, detection=0.001, prevalence=0.5)



### Make lists to use later for Venn diagrams ###


#Responsive, nonresponsive, and healthy comparison at visit 2 both cohorts
list_before_both <- list(Nonresponsive = nonres_2_ASVs_both, Responsive = res_2_ASVs_both, Healthy = healthy_2_ASVs_both)
#Responsive, nonresponsive, and healthy comparison at visit 3 both cohorts
list_after_both <- list(Nonresponsive = nonres_3_ASVs_both, Responsive = res_3_ASVs_both, Healthy = healthy_3_ASVs_both)
#Responsive, nonresponsive, and healthy comparison at visit 2 left cohort
list_before_left <- list(Nonresponsive = nonres_2_ASVs_left, Responsive = res_2_ASVs_left, Healthy = healthy_2_ASVs_left)
#Responsive, nonresponsive, and healthy comparison at visit 3 left cohort
list_after_left <- list(Nonresponsive = nonres_3_ASVs_left, Responsive = res_3_ASVs_left, Healthy = healthy_3_ASVs_left)
#Responsive, nonresponsive, and healthy comparison at visit 2 right cohort
list_before_right <- list(Nonresponsive = nonres_2_ASVs_right, Responsive = res_2_ASVs_right, Healthy = healthy_2_ASVs_right)
#Responsive, nonresponsive, and healthy comparison at visit 3 right cohort
list_after_right <- list(Nonresponsive = nonres_3_ASVs_right, Responsive = res_3_ASVs_right, Healthy = healthy_3_ASVs_right)


dir.create("venn_diagrams_updated")



### Make a Venn-diagram ###

#Responsive, nonresponsive, and healthy comparison at visit 2 both cohorts
venn_before_both <- ggVennDiagram(x = list_before_both)
ggsave("venn_diagrams_updated/Both Cohorts Before Responsive, Nonresponsive, and Healthy at visit 2.png", venn_before_both)
#Responsive, nonresponsive, and healthy comparison at visit 3 both cohorts
venn_after_both <- ggVennDiagram(x = list_after_both)
ggsave("venn_diagrams_updated/Both Cohorts After Responsive, Nonresponsive, and Healthy at visit 3.png", venn_after_both)


#Responsive, nonresponsive, and healthy comparison at visit 2 left cohort
venn_before_left <- ggVennDiagram(x = list_before_left)
ggsave("venn_diagrams_updated/Left Cohort Before Responsive, Nonresponsive, and Healthy at visit 2.png", venn_before_left)
#Responsive, nonresponsive, and healthy comparison at visit 3 left cohort
venn_after_left <- ggVennDiagram(x = list_after_left)
ggsave("venn_diagrams_updated/Left Cohort After Responsive, Nonresponsive, and Healthy at visit 2.png", venn_after_left)

#Responsive, nonresponsive, and healthy comparison at visit 2 right cohort
venn_before_right <- ggVennDiagram(x = list_before_right)
ggsave("venn_diagrams_updated/Right Cohort Before Responsive, Nonresponsive, and Healthy at visit 2.png", venn_before_right)
#Responsive, nonresponsive, and healthy comparison at visit 3 right cohort
venn_after_right <- ggVennDiagram(x = list_after_right)
ggsave("venn_diagrams_updated/Right Cohort After Responsive, Nonresponsive, and Healthy at visit 3.png", venn_after_right)


### Find Unique and Overlapping ASVs for each category ###

# Extract taxonomy table
taxonomy_df <- as.data.frame(tax_table(HIV_final))
taxonomy_df$ASV_ID <- rownames(taxonomy_df)

# Create a function to find overlaps and uniques with taxonomy
find_overlaps_and_uniques_with_taxonomy <- function(asv_lists, taxonomy) {
  overlaps <- list(
    Responsive_Nonresponsive = intersect(asv_lists$Responsive, asv_lists$Nonresponsive),
    Responsive_Healthy = intersect(asv_lists$Responsive, asv_lists$Healthy),
    Nonresponsive_Healthy = intersect(asv_lists$Nonresponsive, asv_lists$Healthy),
    All = Reduce(intersect, asv_lists) # Intersection of all three groups
  )
  uniques <- list(
    Unique_Responsive = setdiff(asv_lists$Responsive, union(asv_lists$Nonresponsive, asv_lists$Healthy)),
    Unique_Nonresponsive = setdiff(asv_lists$Nonresponsive, union(asv_lists$Responsive, asv_lists$Healthy)),
    Unique_Healthy = setdiff(asv_lists$Healthy, union(asv_lists$Responsive, asv_lists$Nonresponsive))
  )
  
  # Map taxonomy to overlaps and uniques
  map_to_taxonomy <- function(asv_list) {
    taxonomy_subset <- merge(data.frame(ASV_ID = asv_list), taxonomy, by = "ASV_ID", all.x = TRUE)
    return(taxonomy_subset)
  }
  
  overlaps_taxonomy <- lapply(overlaps, map_to_taxonomy)
  uniques_taxonomy <- lapply(uniques, map_to_taxonomy)
  
  return(list(Overlaps = overlaps_taxonomy, Uniques = uniques_taxonomy))
}

# Define ASV lists
all_datasets <- list(
  Before_Both = list(Nonresponsive = nonres_2_ASVs_both, Responsive = res_2_ASVs_both, Healthy = healthy_2_ASVs_both),
  After_Both = list(Nonresponsive = nonres_3_ASVs_both, Responsive = res_3_ASVs_both, Healthy = healthy_3_ASVs_both),
  Before_Left = list(Nonresponsive = nonres_2_ASVs_left, Responsive = res_2_ASVs_left, Healthy = healthy_2_ASVs_left),
  After_Left = list(Nonresponsive = nonres_3_ASVs_left, Responsive = res_3_ASVs_left, Healthy = healthy_3_ASVs_left),
  Before_Right = list(Nonresponsive = nonres_2_ASVs_right, Responsive = res_2_ASVs_right, Healthy = healthy_2_ASVs_right),
  After_Right = list(Nonresponsive = nonres_3_ASVs_right, Responsive = res_3_ASVs_right, Healthy = healthy_3_ASVs_right)
)

# Initialize results storage
all_results <- list()

# Process each dataset
for (name in names(all_datasets)) {
  # Calculate overlaps and uniques with taxonomy
  dataset_results <- find_overlaps_and_uniques_with_taxonomy(all_datasets[[name]], taxonomy_df)
  all_results[[name]] <- dataset_results
}

# Save taxonomy results into CSV
output_dir <- "taxonomy_results_with_info"
dir.create(output_dir, showWarnings = FALSE)

for (name in names(all_results)) {
  for (region in names(all_results[[name]])) {
    for (group in names(all_results[[name]][[region]])) {
      file_name <- paste0(output_dir, "/", name, "_", region, "_", group, "_taxonomy.csv")
      write.csv(all_results[[name]][[region]][[group]], file_name, row.names = FALSE)
    }
  }
}



### Make a Venn-diagram with normalized colors and visible headings ###

# Create directory for normalized Venn diagrams
dir.create("venn_diagrams_normalized")

# Function to adjust text placement
adjust_set_labels <- function(plot) {
  plot + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),  # Center and resize title
      text = element_text(size = 12),                    # Adjust overall text size
      legend.position = "right",                         # Keep legend consistent
      plot.margin = margin(20, 20, 20, 40)               # Add margin to prevent cutoff
    ) +
    scale_fill_gradient(low = "white", high = "blue", limits = c(0, 31))    # Normalize color gradient
}

# Before ART - Both Cohorts
venn_before_both <- ggVennDiagram(list_before_both) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 31)) +
  labs(title = "Before ART - Both Cohorts") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    text = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 40)
  )
ggsave("venn_diagrams_normalized/Before_ART_Both_Cohorts.png", venn_before_both)

# After ART - Both Cohorts
venn_after_both <- ggVennDiagram(list_after_both) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 31)) +
  labs(title = "After ART - Both Cohorts") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    text = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 40)
  )
ggsave("venn_diagrams_normalized/After_ART_Both_Cohorts.png", venn_after_both)

# Before ART - Left Cohort
venn_before_left <- ggVennDiagram(list_before_left) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 31)) +
  labs(title = "0 Weeks - Cohort A") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    text = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 40)
  )
ggsave("venn_diagrams_normalized/Before_ART_Left_Cohort.png", venn_before_left)

# After ART - Left Cohort
venn_after_left <- ggVennDiagram(list_after_left) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 31)) +
  labs(title = "24 Weeks - Cohort A") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    text = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 40)
  )
ggsave("venn_diagrams_normalized/After_ART_Left_Cohort.png", venn_after_left)

# Before ART - Right Cohort
venn_before_right <- ggVennDiagram(list_before_right) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 31)) +
  labs(title = "0 Weeks - Cohort B") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    text = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 40)
  )
ggsave("venn_diagrams_normalized/Before_ART_Right_Cohort.png", venn_before_right)

# After ART - Right Cohort
venn_after_right <- ggVennDiagram(list_after_right) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 31)) +
  labs(title = "24 Weeks - Cohort B") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    text = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 40)
  )
ggsave("venn_diagrams_normalized/After_ART_Right_Cohort.png", venn_after_right)


# Load libraries for making combined figure
library(ggVennDiagram)
library(patchwork)

# Combine the diagrams into a figure
combined_plot <- (venn_before_left | venn_after_left) / (venn_before_right | venn_after_right) +
  plot_layout(guides = "collect") # Share legend across plots

# Display the combined plot
print(combined_plot)

# Save the combined plot as a single image
ggsave("venn_diagrams_normalized/venn_diagrams_panel.png", combined_plot, width = 16, height = 16)
