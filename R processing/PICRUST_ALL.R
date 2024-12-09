
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)

# Import PICRUST2 pathway file
abundance_data <- read_delim("pathabun_exported/pathway_abundance.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE) %>%
  as.data.frame()

# Import metadata and filter groups for pairwise comparison
left <- read_delim("left cohort.csv")
right <- read_delim("right cohort.csv")

left$response_patient_by_visit <- factor(left$response_patient_by_visit, 
                                         levels = c("Healthy 2", "nonresponsive 2", "responsive 2","Healthy 3", "nonresponsive 3", "responsive 3" ))
right$response_patient_by_visit <- factor(right$response_patient_by_visit, 
                                          levels = c("Healthy 2", "nonresponsive 2", "responsive 2","Healthy 3", "nonresponsive 3", "responsive 3" ))

# Load PICRUST analysis function
source("PICRUST_Function.R")

# Week 0
left_v2_NvH <- filter(left, response_patient_by_visit=="nonresponsive 2" | response_patient_by_visit=="Healthy 2")
left_v2_NvHgg <- PICRUST_function(left_v2_NvH, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Left Cohort, Nonresponsive vs Healthy, Week 0")
ggsave(left_v2_NvHgg, file="PICRUST/left_v2_NvH.jpg", width=11)

left_v2_RvH <- filter(left, response_patient_by_visit=="responsive 2" | response_patient_by_visit=="Healthy 2")
left_v2_RvHgg <- PICRUST_function(left_v2_RvH, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Left Cohort, Responsive vs Healthy, Week 0")
ggsave(left_v2_RvHgg, file="PICRUST/left_v2_RvH.jpg", width=11)

left_v2_RvN <- filter(left, response_patient_by_visit=="responsive 2" | response_patient_by_visit=="nonresponsive 2")
left_v2_RvNgg <- PICRUST_function(left_v2_RvN, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Left Cohort, Responsive vs Nonresponsive, Week 0")
ggsave(left_v2_RvNgg, file="PICRUST/left_v2_RvN.jpg", width=11)

right_v2_NvH <- filter(right, response_patient_by_visit=="nonresponsive 2" | response_patient_by_visit=="Healthy 2")
right_v2_NvHgg <- PICRUST_function(right_v2_NvH, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Right Cohort, Nonresponsive vs Healthy, Week 0")
ggsave(right_v2_NvHgg, file="PICRUST/right_v2_NvH.jpg", width=11)

right_v2_RvH <- filter(right, response_patient_by_visit=="responsive 2" | response_patient_by_visit=="Healthy 2")
right_v2_RvHgg <- PICRUST_function(right_v2_RvH, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Right Cohort, Responsive vs Healthy, Week 0")
ggsave(right_v2_RvHgg, file="PICRUST/right_v2_RvH.jpg", width=11)

right_v2_RvN <- filter(right, response_patient_by_visit=="responsive 2" | response_patient_by_visit=="nonresponsive 2")
right_v2_RvNgg <- PICRUST_function(right_v2_RvN, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Right Cohort, Responsive vs Nonresponsive, Week 0")
ggsave(right_v2_RvNgg, file="PICRUST/right_v2_RvN.jpg", width=11)


## Week 24

left_v3_NvH <- filter(left, response_patient_by_visit=="nonresponsive 3" | response_patient_by_visit=="Healthy 3")
left_v3_NvHgg <- PICRUST_function(left_v3_NvH, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Left Cohort, Nonresponsive vs Healthy, Week 24")
ggsave(left_v3_NvHgg, file="PICRUST/left_v3_NvH.jpg", width=11)

left_v3_RvH <- filter(left, response_patient_by_visit=="responsive 3" | response_patient_by_visit=="Healthy 3")
left_v3_RvHgg <- PICRUST_function(left_v3_RvH, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Left Cohort, Responsive vs Healthy, Week 24")
ggsave(left_v3_RvHgg, file="PICRUST/left_v3_RvH.jpg", width=11)

left_v3_RvN <- filter(left, response_patient_by_visit=="responsive 3" | response_patient_by_visit=="nonresponsive 3")
left_v3_RvNgg <- PICRUST_function(left_v3_RvN, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Left Cohort, Responsive vs Nonresponsive, Week 24")
ggsave(left_v3_RvNgg, file="PICRUST/left_v3_RvN.jpg", width=11)

right_v3_NvH <- filter(right, response_patient_by_visit=="nonresponsive 3" | response_patient_by_visit=="Healthy 3")
right_v3_NvHgg <- PICRUST_function(right_v3_NvH, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Right Cohort, Nonresponsive vs Healthy, Week 24")
ggsave(right_v3_NvHgg, file="PICRUST/right_v3_NvH.jpg", width=11)

right_v3_RvH <- filter(right, response_patient_by_visit=="responsive 3" | response_patient_by_visit=="Healthy 3")
right_v3_RvHgg <- PICRUST_function(right_v3_RvH, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Right Cohort, Responsive vs Healthy, Week 24")
ggsave(right_v3_RvHgg, file="PICRUST/right_v3_RvH.jpg", width=11)

right_v3_RvN <- filter(right, response_patient_by_visit=="responsive 3" | response_patient_by_visit=="nonresponsive 3")
right_v3_RvNgg <- PICRUST_function(right_v3_RvN, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Right Cohort, Responsive vs Nonresponsive, Week 24")
ggsave(right_v3_RvNgg, file="PICRUST/right_v3_RvN.jpg", width=11)

# Week 24 vs 0
right_v2v3_R <- filter(right, response_patient_by_visit=="responsive 2" | response_patient_by_visit=="responsive 3")
right_v2v3_Rgg <- PICRUST_function(right_v2v3_R, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Right Cohort, Responsive, Week 24 vs Week 0")
ggsave(right_v2v3_Rgg, file="PICRUST/right_v2v3_R.jpg", width=11)

right_v2v3_NR <- filter(right, response_patient_by_visit=="nonresponsive 2" | response_patient_by_visit=="nonresponsive 3")
right_v2v3_NRgg <- PICRUST_function(right_v2v3_NR, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Right Cohort, Nonresponsive, Week 24 vs Week 0")
ggsave(right_v2v3_NRgg, file="PICRUST/right_v2v3_NR.jpg", width=11)

right_v2v3_H <- filter(right, response_patient_by_visit=="Healthy 2" | response_patient_by_visit=="Healthy 3")
right_v2v3_Hgg <- PICRUST_function(right_v2v3_H, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Right Cohort, Healthy, Week 24 vs Week 0")
ggsave(right_v2v3_Hgg, file="PICRUST/right_v2v3_H.jpg", width=11)

left_v2v3_R <- filter(left, response_patient_by_visit=="responsive 2" | response_patient_by_visit=="responsive 3")
left_v2v3_Rgg <- PICRUST_function(left_v2v3_R, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Left Cohort, Responsive, Week 24 vs Week 0")
ggsave(left_v2v3_Rgg, file="PICRUST/left_v2v3_R.jpg", width=11)

left_v2v3_NR <- filter(left, response_patient_by_visit=="nonresponsive 2" | response_patient_by_visit=="nonresponsive 3")
left_v2v3_NRgg <- PICRUST_function(left_v2v3_NR, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Left Cohort, Nonresponsive, Week 24 vs Week 0")
ggsave(left_v2v3_NRgg, file="PICRUST/left_v2v3_NR.jpg", width=11)

left_v2v3_H <- filter(left, response_patient_by_visit=="Healthy 2" | response_patient_by_visit=="Healthy 3")
left_v2v3_Hgg <- PICRUST_function(left_v2v3_H, abundance_data) +
  labs(x = "Log2FoldChange", y="Pathways", title="Left Cohort, Healthy, Week 24 vs Week 0")
ggsave(left_v2v3_Hgg, file="PICRUST/left_v2v3_H.jpg", width=11)

