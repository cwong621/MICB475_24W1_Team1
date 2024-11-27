# Load libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library("sf")

# Load data
load("HIV Phyloseq/rightphyloseq.RData")

# See metadata
SAMP <- sample_data(ps_right)

#convert phyloseq object to relative abundance
HIV_final_RA <- transform_sample_counts(ps_right, fun=function(x) x/sum(x))


# Filter dataset (high load responsive before and after, low load responsive before and after, 
#high load nonresponsive before and after, low load nonresponsive before and after, 
#responsive and nonresponsive before and after (regardless of viral load), 
#healthy controls before and after, ART experienced nonresponsive before and after)
### Note: before and after mean at visit 2 and 3, respectively

#high load responsive visit 2 (before ART)
high_res_2 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="responsive 2" & `start_viral_load_patient_by_visit`=="HIV high pretreatment 2")
#high load responsive visit 3 (after ART)
high_res_3 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="responsive 3" & `start_viral_load_patient_by_visit`=="HIV high pretreatment 3")
#low load responsive visit 2 (before ART)
low_res_2 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="responsive 2" & `start_viral_load_patient_by_visit`=="HIV low pretreatment 2")
#low load responsive visit 3 (after ART)
low_res_3 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="responsive 3" & `start_viral_load_patient_by_visit`=="HIV low pretreatment 3")
#high load nonresponsive visit 2 (before ART)
high_nonres_2 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="nonresponsive 2" & `start_viral_load_patient_by_visit`=="HIV high pretreatment 2")
#high load nonresponsive visit 3 (after ART)
high_nonres_3 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="nonresponsive 3" & `start_viral_load_patient_by_visit`=="HIV high pretreatment 3")
#low load nonresponsive visit 2 (before ART)
low_nonres_2 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="nonresponsive 2" & `start_viral_load_patient_by_visit`=="HIV low pretreatment 2")
#low load nonresponsive visit 3 (after ART)
low_nonres_3 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="nonresponsive 3" & `start_viral_load_patient_by_visit`=="HIV low pretreatment 3")
#responsive visit 3 (regardless of viral load)
res_3 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="responsive 3")
#nonresponsive visit 3 (regardless of viral load)
nonres_3 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="nonresponsive 3")

#responsive visit 2 (regardless of viral load)
res_2 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="responsive 2")
#nonresponsive visit 2 (regardless of viral load)
nonres_2 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="nonresponsive 2")

#Healthy at visit 2
healthy_2 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="Healthy 2" & `start_viral_load_patient_by_visit`=="Healthy 2")
#Healthy at visit 3
healthy_3 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="Healthy 3" & `start_viral_load_patient_by_visit`=="Healthy 3")

#ART experienced nonresponsive at visit 2
art_exp_nonres_2 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="ART-experienced nonresponsive V2 2" & `start_viral_load_patient_by_visit`=="NA 2")
#ART experienced nonresponsive at visit 3
art_exp_nonres_3 <- subset_samples(HIV_final_RA, `response_patient_by_visit`=="ART-experienced nonresponsive V3 3" & `start_viral_load_patient_by_visit`=="NA 3")


### Abundance and prevalence - Find ASVs for the different categories ###
#high load responsive visit 2 (before ART)
high_res_2_ASVs <- core_members(high_res_2, detection=0.001, prevalence=0.5)
#high load responsive visit 3 (after ART)
high_res_3_ASVs <- core_members(high_res_3, detection=0.001, prevalence=0.5)
#low load responsive visit 2 (before ART)
low_res_2_ASVs <- core_members(low_res_2, detection=0.001, prevalence=0.5)
#low load responsive visit 3 (after ART)
low_res_3_ASVs <- core_members(low_res_3, detection=0.001, prevalence=0.5)
#high load nonresponsive visit 2 (before ART)
high_nonres_2_ASVs <- core_members(high_nonres_2, detection=0.001, prevalence=0.5)
#high load nonresponsive visit 3 (after ART)
high_nonres_3_ASVs <- core_members(high_nonres_3, detection=0.001, prevalence=0.5)
#low load nonresponsive visit 2 (before ART)
low_nonres_2_ASVs <- core_members(low_nonres_2, detection=0.001, prevalence=0.5)
#high load nonresponsive visit 3 (after ART)
low_nonres_3_ASVs <- core_members(low_nonres_3, detection=0.001, prevalence=0.5)
#responsive visit 3 (regardless of viral load)
res_3_ASVs <- core_members(res_3, detection=0.001, prevalence=0.5)
#nonresponsive visit 3 (regardless of viral load)
nonres_3_ASVs <- core_members(nonres_3, detection=0.001, prevalence=0.5)

#responsive visit 2 (regardless of viral load)
res_2_ASVs <- core_members(res_3, detection=0.001, prevalence=0.5)
#nonresponsive visit 2 (regardless of viral load)
nonres_2_ASVs <- core_members(nonres_3, detection=0.001, prevalence=0.5)

#Healthy at visit 2
healthy_2_ASVs <- core_members(healthy_2, detection=0.001, prevalence=0.5)
#Healthy at visit 3
healthy_3_ASVs <- core_members(healthy_3, detection=0.001, prevalence=0.5)

#ART experienced nonresponsive at visit 2
art_exp_nonres_2_ASVs <- core_members(art_exp_nonres_2, detection=0.001, prevalence=0.5)
#ART experienced nonresponsive at visit 3
art_exp_nonres_3_ASVs <- core_members(art_exp_nonres_3, detection=0.001, prevalence=0.5)



### Make lists to use later for Venn diagrams ###

#HIV high load responsive comparison visit 2 and 3
list_high_resp_2_and_3 <- list(HIV_high_responsive_before_ART = high_res_2_ASVs, HIV_high_responsive_after_ART = high_res_3_ASVs)

#HIV low load responsive comparison visit 2 and 3
list_low_resp_2_and_3 <- list(HIV_low_responsive_before_ART = low_res_2_ASVs, HIV_low_responsive_after_ART = low_res_3_ASVs)

#Responsive comparison visit 2 and 3
list_resp_2_and_3 <- list(responsive_before_ART = res_2_ASVs, responsive_after_ART = res_3_ASVs)


#HIV high load nonresponsive comparison visit 2 and 3
list_high_nonresp_2_and_3 <- list(HIV_high_nonresponsive_before_ART = high_nonres_2_ASVs, HIV_high_nonresponsive_after_ART = high_nonres_3_ASVs)

#HIV low load nonresponsive comparison visit 2 and 3
list_low_nonresp_2_and_3 <- list(HIV_low_nonresponsive_before_ART = low_nonres_2_ASVs, HIV_low_nonresponsive_after_ART = low_nonres_3_ASVs)

#Nonresponsive comparison visit 2 and 3
list_nonresp_2_and_3 <- list(nonresponsive_before_ART = nonres_2_ASVs, nonresponsive_after_ART = nonres_3_ASVs)

#HIV high load responsive, nonresponsive, and healthy comparison at visit 3
list_HIV_high_healthy_resp_nonresp <- list(HIV_high_nonresponsive_after_ART = high_nonres_3_ASVs, HIV_high_responsive_after_ART = high_res_3_ASVs, Healthy = healthy_3_ASVs)

#HIV low load responsive, nonresponsive, and healthy comparison at visit 3
list_HIV_low_healthy_resp_nonresp <- list(HIV_low_nonresponsive_after_ART = low_nonres_3_ASVs, HIV_low_responsive_after_ART = low_res_3_ASVs, Healthy = healthy_3_ASVs)

#Responsive, nonresponsive, and healthy comparison at visit 3
list_healthy_resp_nonresp <- list(Nonresponsive = nonres_3_ASVs, Responsive = res_3_ASVs, Healthy = healthy_3_ASVs)

#Healthy comparison at visit 2 and 3
list_healthy_2and3 <- list(healthy_visit_2 = healthy_2_ASVs, healthy_visit_3 = healthy_3_ASVs)

#ART experienced nonresponsive comparison at visit 2 and 3
list_artexp_2and3 <- list(ART_experienced_nonresponsive_2 = art_exp_nonres_2_ASVs, ART_experienced_nonresponsive_3 = art_exp_nonres_2_ASVs)

#HIV high responsive, nonresponsive, and ART experienced at visit 3
list_high_resp_nonresp_artexp <- list(HIV_high_responsive_after_ART = high_res_3_ASVs, HIV_high_nonresponsive_after_ART = high_nonres_3_ASVs, ART_experienced_nonresponsive = art_exp_nonres_3_ASVs)

#HIV low responsive, nonresponsive, and ART experienced at visit 3
list_low_resp_nonresp_artexp <- list(HIV_low_responsive_after_ART = low_res_3_ASVs, HIV_low_nonresponsive_after_ART = low_nonres_3_ASVs, ART_experienced_nonresponsive = art_exp_nonres_3_ASVs)

#Responsive, nonresponsive, and ART experienced at visit 3 (regardless of viral load)
list_resp_nonresp_artexp <- list(responsive_after_ART = res_3_ASVs, nonresponsive_after_ART = nonres_3_ASVs, ART_experienced_nonresponsive = art_exp_nonres_3_ASVs)

#Healthy, responsive, nonresponsive, and ART experienced at visit 3 (regardless of viral load)
list_healthy_resp_nonresp_artexp <- list(Healthy = healthy_3_ASVs, responsive_after_ART = res_3_ASVs, nonresponsive_after_ART = nonres_3_ASVs, ART_experienced_nonresponsive = art_exp_nonres_3_ASVs)

#Responsive, nonresponsive, and ART experienced at visit 2 (regardless of viral load)
list_resp_nonresp_artexp_2 <- list(responsive_after_ART = res_2_ASVs, nonresponsive_after_ART = nonres_2_ASVs, ART_experienced_nonresponsive = art_exp_nonres_2_ASVs)

#Healthy, responsive, nonresponsive, and ART experienced at visit 2 (regardless of viral load)
list_healthy_resp_nonresp_artexp_2 <- list(Healthy = healthy_2_ASVs, responsive_after_ART = res_2_ASVs, nonresponsive_after_ART = nonres_2_ASVs, ART_experienced_nonresponsive = art_exp_nonres_2_ASVs)

#Responsive, nonresponsive, and healthy comparison at visit 2
list_healthy_resp_nonresp_2 <- list(Nonresponsive = nonres_2_ASVs, Responsive = res_2_ASVs, Healthy = healthy_2_ASVs)


dir.create("venn_diagrams_right_cohort")




### Make a Venn-diagram ###

#HIV high load responsive comparison visit 2 and 3
right_venn_high_resp_2_and_3 <- ggVennDiagram(x = list_high_resp_2_and_3)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn HIV high responsive at visit 2 and 3.png", right_venn_high_resp_2_and_3)

#HIV low load responsive comparison visit 2 and 3
right_venn_low_resp_2_and_3 <- ggVennDiagram(x = list_low_resp_2_and_3)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn HIV low responsive at visit 2 and 3.png", right_venn_low_resp_2_and_3)

#HIV high load nonresponsive comparison visit 2 and 3
right_venn_high_nonresp_2_and_3 <- ggVennDiagram(x = list_high_nonresp_2_and_3)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn HIV high nonresponsive at visit 2 and 3.png", right_venn_high_nonresp_2_and_3)


#HIV low load nonresponsive comparison visit 2 and 3
right_venn_low_nonresp_2_and_3 <- ggVennDiagram(x = list_low_nonresp_2_and_3)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn HIV low nonresponsive at visit 2 and 3.png", right_venn_low_nonresp_2_and_3)


#HIV high load responsive, nonresponsive, and healthy comparison at visit 3
right_venn_HIV_high_healthy_resp_nonresp <- ggVennDiagram(x = list_HIV_high_healthy_resp_nonresp)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn HIV high responsive, nonresponsive and healthy at visit 3.png", right_venn_HIV_high_healthy_resp_nonresp)


#HIV low load responsive, nonresponsive, and healthy comparison at visit 3
right_venn_HIV_low_healthy_resp_nonresp <- ggVennDiagram(x = list_HIV_low_healthy_resp_nonresp)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn HIV low responsive, nonresponsive and healthy at visit 3.png", right_venn_HIV_low_healthy_resp_nonresp)


#Responsive, nonresponsive, and healthy comparison at visit 3 (regardless of viral load)
right_venn_healthy_resp_nonresp <- ggVennDiagram(x = list_healthy_resp_nonresp)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn Responsive, nonresponsive and healthy at visit 3.png", right_venn_healthy_resp_nonresp)


#Healthy comparison at visit 2 and 3
right_venn_healthy_2and3 <- ggVennDiagram(x = list_healthy_2and3)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn Healthy at visits 2 and 3.png", right_venn_healthy_2and3)


#ART experienced nonresponsive comparison at visit 2 and 3
right_venn_artexp_2and3 <- ggVennDiagram(x = list_artexp_2and3)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn ART experienced nonresponsive at visits 2 and 3.png", right_venn_artexp_2and3)


#HIV high responsive, nonresponsive, and ART experienced at visit 3
right_venn_high_resp_nonresp_artexp <- ggVennDiagram(x = list_high_resp_nonresp_artexp)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn HIV high responsive, nonresponsive, and ART experienced at visit 3.png", right_venn_high_resp_nonresp_artexp)


#HIV low responsive, nonresponsive, and ART experienced at visit 3
right_venn_low_resp_nonresp_artexp <- ggVennDiagram(x = list_low_resp_nonresp_artexp)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn HIV low responsive, nonresponsive, and ART experienced at visit 3.png", right_venn_low_resp_nonresp_artexp)


#Responsive, nonresponsive, and ART experienced at visit 3 (regardless of viral load)
right_venn_resp_nonresp_artexp <- ggVennDiagram(x = list_resp_nonresp_artexp)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn Responsive, nonresponsive and ART experienced nonresponsive at visit 3.png", right_venn_resp_nonresp_artexp)


#Healthy, responsive, nonresponsive, and ART experienced at visit 3 (regardless of viral load)
right_venn_healthy_resp_nonresp_artexp <- ggVennDiagram(x = list_healthy_resp_nonresp_artexp)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn Responsive, nonresponsive, ART experienced nonresponsive, and Healthy at visit 3.png", right_venn_healthy_resp_nonresp_artexp)

#Responsive, nonresponsive, and ART experienced at visit 2 (regardless of viral load)
right_venn_resp_nonresp_artexp_2 <- ggVennDiagram(x = list_resp_nonresp_artexp_2)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn Responsive, nonresponsive and ART experienced nonresponsive at visit 2.png", right_venn_resp_nonresp_artexp_2)

#Healthy, responsive, nonresponsive, and ART experienced at visit 2 (regardless of viral load)
right_venn_healthy_resp_nonresp_artexp_2 <- ggVennDiagram(x = list_healthy_resp_nonresp_artexp_2)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn Responsive, nonresponsive, ART experienced nonresponsive, and Healthy at visit 2.png", right_venn_healthy_resp_nonresp_artexp_2)

#Responsive, nonresponsive, and healthy comparison at visit 2
right_venn_healthy_resp_nonresp_2 <- ggVennDiagram(x = list_healthy_resp_nonresp_2)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn Responsive, nonresponsive, and Healthy at visit 2.png", right_venn_healthy_resp_nonresp_2)


#Nonresponsive comparison visit 2 and 3
right_venn_nonresp_2_and_3 <- ggVennDiagram(x = list_nonresp_2_and_3)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn Nonresponsive at visit 2 and 3.png", right_venn_nonresp_2_and_3)


#Responsive comparison visit 2 and 3
right_venn_resp_2_and_3 <- ggVennDiagram(x = list_resp_2_and_3)
ggsave("venn_diagrams_right_cohort/Right Cohort Venn Responsive at visit 2 and 3.png", right_venn_resp_2_and_3)







