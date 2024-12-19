# Load libraries

library(phyloseq)
library(vegan)
library(tidyverse)

# Load data
load("./leftphyloseq_rare.RData")
load("./rightphyloseq_rare.RData")

# Adding phylogenetic distance to phyloseq objects

sample_data_left_wPD <- data.frame(sample_data(leftphyloseq_rare), estimate_richness(leftphyloseq_rare))
sample_data_right_wPD <- data.frame(sample_data(rightphyloseq_rare), estimate_richness(rightphyloseq_rare))

###### PERMANOVA on Weighted Unifrac ######

# Calculate weighted unifrac distance matrix (left)
wu_dm_left <- UniFrac(leftphyloseq_rare, weighted = TRUE)
pcoa_wu_left <- ordinate(leftphyloseq_rare, method="PCoA", distance=wu_dm_left)
plot_ordination(leftphyloseq_rare, pcoa_wu_left, color = "response_patient_by_visit")

adonis2(wu_dm_left ~ `response_patient_by_visit`*`start_viral_load_patient_by_visit`, data = sample_data_left_wPD)

adonis2(wu_dm_left ~ `response_patient_by_visit`, data = sample_data_left_wPD)
# Pr(>F) 0.927

# Calculate weighted unifrac distance matrix (right)
wu_dm_right <- UniFrac(rightphyloseq_rare, weighted = TRUE)
pcoa_wu_right <- ordinate(rightphyloseq_rare, method="PCoA", distance=wu_dm_right)
plot_ordination(rightphyloseq_rare, pcoa_wu_right, color = "response_patient_by_visit")

adonis2(wu_dm_right ~ `response_patient_by_visit`*`start_viral_load_patient_by_visit`, data = sample_data_right_wPD)

adonis2(wu_dm_right ~ `response_patient_by_visit`, data = sample_data_right_wPD)
# Pr(>F) 0.827

###### PERMANOVA on Bray-Curtis #####

# Calculate Bray-Curtis distance matrix (left)
bc_dm_left <- distance(leftphyloseq_rare, method="bray")
pcoa_bc_left <- ordinate(leftphyloseq_rare, method="PCoA", distance=bc_dm_left)
plot_ordination(leftphyloseq_rare, pcoa_bc_left, color = "response_patient_by_visit")

adonis2(bc_dm_left ~ `response_patient_by_visit`*`start_viral_load_patient_by_visit`, data = sample_data_left_wPD)

adonis2(bc_dm_left ~ `response_patient_by_visit`, data = sample_data_left_wPD)
# Pr(>F) 0.324

# Calculate Bray-Curtis distance matrix (right)
bc_dm_right <- distance(rightphyloseq_rare, method="bray")
pcoa_bc_right <- ordinate(rightphyloseq_rare, method="PCoA", distance=bc_dm_right)
plot_ordination(rightphyloseq_rare, pcoa_bc_right, color = "response_patient_by_visit")

adonis2(bc_dm_right ~ `response_patient_by_visit`*`start_viral_load_patient_by_visit`, data = sample_data_right_wPD)

adonis2(bc_dm_right ~ `response_patient_by_visit`, data = sample_data_right_wPD)
# Pr(>F) 0.118

###### PERMANOVA on Jaccard #####

# Calculate Jaccard distance matrix (left)
jc_dm_left <- distance(leftphyloseq_rare, method="jaccard", binary = TRUE)
pcoa_jc_left <- ordinate(leftphyloseq_rare, method="PCoA", distance=jc_dm_left)
plot_ordination(leftphyloseq_rare, pcoa_jc_left, color = "response_patient_by_visit")

jaccard_left_stat <- adonis2(jc_dm_left ~ `response_patient_by_visit`*`start_viral_load_patient_by_visit`, data = sample_data_left_wPD)

adonis2(jc_dm_left ~ `response_patient_by_visit`, data = sample_data_left_wPD)
# Pr(>F) 0.047 *

# Calculate Jaccard distance matrix (right)
jc_dm_right <- distance(rightphyloseq_rare, method="jaccard", binary = TRUE)
pcoa_jc_right <- ordinate(rightphyloseq_rare, method="PCoA", distance=jc_dm_right)
plot_ordination(rightphyloseq_rare, pcoa_jc_right, color = "response_patient_by_visit")

adonis2(jc_dm_right ~ `response_patient_by_visit`*`start_viral_load_patient_by_visit`, data = sample_data_right_wPD)

adonis2(jc_dm_right ~ `response_patient_by_visit`, data = sample_data_right_wPD) 
# Pr(>F) 0.2

######### PERMANOVA ON UNWEIGHTED UNIFRAC ########

# calculate distance matrix right
uu_dm_right <- UniFrac(rightphyloseq_rare, weighted = FALSE)
pcoa_uu_right <- ordinate(rightphyloseq_rare, method="PCoA", distance=uu_dm_right)
plot_ordination(rightphyloseq_rare, pcoa_uu_right, color = "response_patient_by_visit")

adonis2(uu_dm_right ~ `response_patient_by_visit`*`start_viral_load_patient_by_visit`, data = sample_data_right_wPD)

adonis2(uu_dm_right ~ `response_patient_by_visit`, data = sample_data_right_wPD) 
# Pr(>F) 0.602

# calculate distance matrix left
uu_dm_left <- UniFrac(leftphyloseq_rare, weighted = FALSE)
pcoa_uu_left <- ordinate(leftphyloseq_rare, method="PCoA", distance=uu_dm_left)
plot_ordination(leftphyloseq_rare, pcoa_uu_left, color = "response_patient_by_visit")

adonis2(uu_dm_left ~ `response_patient_by_visit`*`start_viral_load_patient_by_visit`, data = sample_data_left_wPD)

adonis2(uu_dm_left ~ `response_patient_by_visit`, data = sample_data_left_wPD) 
# Pr(>F) 0.132


