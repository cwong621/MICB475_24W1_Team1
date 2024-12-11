###### Running a pairwise PERMANOVA on Jaccard left cohort ######

# Load libraries
library(pairwiseAdonis)
library(dplyr)
library(ggplot2)

# Run pairwise PERMANOVA
pairwise_results <- pairwise.adonis2(
  jc_dm_left ~ response_patient_by_visit,
  data = sample_data_left_wPD,
  permutations = 999
)
pairwise_results

# Select for statistically significant values
pairwise_df <- as.data.frame(pairwise_results)
                 
pr_columns <- grep("\\.Pr\\.\\.F\\.$", colnames(pairwise_df), value = TRUE)

significant_results <- list()

for (col in pr_columns) {
  significant_indices <- which(pairwise_df[[col]] < 0.05)
  if (length(significant_indices) > 0) {
    significant_results[[col]] <- significant_indices
  }
}

# Create a data frame with significant pairings
significant_results <- as.data.frame(significant_results)

significant_results <- significant_results %>%
  pivot_longer(
    cols = everything(),
    names_to = "Comparison"
  )

significant_results <- significant_results %>%
  mutate(c(0.014,0.012,0.001,0.002,0.025))

colnames(significant_results) <- c("Comparison", "value", "Pr(<F)")

significant_pairs <- data.frame(
  Comparison = c("Healthy.3_vs_responsive.3", "Healthy.3_vs_nonresponsive.3", 
                 "Healthy.2_vs_responsive.3", "Healthy.2_vs_nonresponsive.3", 
                 "responsive.3_vs_nonresponsive.3"),
  Pr_F = c(0.014, 0.012, 0.001, 0.002, 0.025)
)

significant_pairs <- significant_pairs %>%
  mutate(Group1 = str_extract(Comparison, "^[^_]+"),  # Extract text before '_vs_'
         Group2 = str_extract(Comparison, "(?<=_vs_).*"))  # Extract text after '_vs_'

significant_pairs <- significant_pairs %>%
  mutate(
    Group1 = gsub("\\.", " ", Group1),  # Replace '.' with a space
    Group2 = gsub("\\.", " ", Group2)   # Replace '.' with a space
)

significant_pairs <- significant_pairs %>%
  mutate(Pair_label = paste0("Pair_", row_number()))

significant_groups <- unique(c(significant_pairs$Group1, significant_pairs$Group2))

# Subset for significant groups
pcoa_filtered <- pcoa_coords_w_metadata %>%
  filter(response_patient_by_visit %in% significant_groups)

pairs_long <- significant_pairs %>%
  select(Group1, Group2, Pair_label) %>%
  pivot_longer(cols = c(Group1, Group2), names_to = "Group_Type", values_to = "response_patient_by_visit") %>%
  select(-Group_Type)
  
# Join pairs_long with pcoa_filtered
pcoa_filtered_with_pairs <- pairs_long %>%
  left_join(pcoa_filtered, by = "response_patient_by_visit")

# Plot the filtered PCA plot

ggplot(pcoa_healthy_3_responsive_3, aes(x = PC1, y= PC2, color = response_patient_by_visit)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill=response_patient_by_visit), type = "t", alpha = 0.2) +
  labs(
    title = "PCoA Plot for Healthy (after) and Responsive (after)",
    x = "PC1 [6.2%]",
    y = "PC2 [4.4%]"
  ) 
