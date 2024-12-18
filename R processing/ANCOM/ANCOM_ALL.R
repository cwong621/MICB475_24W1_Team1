# Load packages
library(ANCOMBC)
library(tidyverse)

# Load phyloseq objects
load("./leftphyloseq.RData")
load("./rightphyloseq.RData")

# Set levels in response_patient_by_visit metadata column
ps_left@sam_data$response_patient_by_visit <- factor(ps_left@sam_data$response_patient_by_visit, 
            levels = c("Healthy 2", "nonresponsive 2", "responsive 2", "ART-experienced nonresponsive 2",  "Healthy 3", "nonresponsive 3", "responsive 3", "ART-experienced nonresponsive 3" ))
ps_right@sam_data$response_patient_by_visit <- factor(ps_right@sam_data$response_patient_by_visit, 
                                                     levels = c("Healthy 2", "nonresponsive 2", "responsive 2", "ART-experienced nonresponsive 2",  "Healthy 3", "nonresponsive 3", "responsive 3", "ART-experienced nonresponsive 3" ))
# Load ANCOM functions
source("ANCOMBC_Function.R")

####PAIRWISE COMPARISONS ####

## Nonresponsive/Responsive vs Healthy
ps_leftv2 <- subset_samples(ps_left, Visit==2 & Cohort_Short!="B")
leftv2 <- ANCOMBC_HealthyRef(ps_leftv2) +
  labs(x = NULL, y = NULL, title = "Cohort A, Week 0") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/leftv2.jpg", leftv2)

ps_rightv2 <- subset_samples(ps_right, Visit==2 & Cohort_Short!="B")
rightv2 <- ANCOMBC_HealthyRef(ps_rightv2) +
  labs(x = NULL, y = NULL, title = "Cohort B, Week 0") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/rightv2.jpg", rightv2)

ps_rightv3 <- subset_samples(ps_right, Visit==3 & Cohort_Short!="B")
rightv3 <- ANCOMBC_HealthyRef(ps_rightv3) +
  labs(x = NULL, y = NULL, title = "Cohort B, Week 24") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/rightv3.jpg", rightv3)

ps_leftv3 <- subset_samples(ps_left, Visit==3 & Cohort_Short!="B")
leftv3 <- ANCOMBC_HealthyRef(ps_leftv3) +
  labs(x = NULL, y = NULL, title = "Cohort A, Week 24") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/leftv3.jpg", leftv3)

## Responsive vs Nonresponsive

ps_leftv2NRRef <- subset_samples(ps_left, Visit==2 & Cohort_Short=="A")
leftv2NRRef <- ANCOMBC_NRRef(ps_leftv2NRRef) +
  labs(x = NULL, y = NULL, title = "Cohort A, Week 0") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/leftv2NRRef.jpg", leftv2NRRef)

ps_rightv2NRRef <- subset_samples(ps_right, Visit==2 & Cohort_Short=="A")
rightv2NRRef <- ANCOMBC_NRRef(ps_rightv2NRRef) +
  labs(x = NULL, y = NULL, title = "Cohort B, Week 0") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/rightv2NRRef.jpg", rightv2NRRef)

ps_rightv3NRRef <- subset_samples(ps_right, Visit==3 & Cohort_Short=="A")
rightv3NRRef <- ANCOMBC_NRRef(ps_rightv3NRRef) +
  labs(x = NULL, y = NULL, title = "Cohort B, Week 24") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/rightv3NRRef.jpg", rightv3NRRef)

ps_leftv3NRRef <- subset_samples(ps_left, Visit==3 & Cohort_Short=="A")
leftv3NRRef <- ANCOMBC_NRRef(ps_leftv3NRRef) +
  labs(x = NULL, y = NULL, title = "Cohort A, Week 24") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/leftv3NRRef.jpg", leftv3NRRef)

## Week 24 vs Week 0 comparisons
ps_leftNR <- subset_samples(ps_left, response_patient=="nonresponsive" & Cohort_Short!="B")
leftNR <- ANCOMBC_V3V2(ps_leftNR) +
  labs(x = NULL, y = NULL, title = "Cohort A, Nonresponsive") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/leftNR.jpg", leftNR)

ps_leftR <- subset_samples(ps_left, response_patient=="responsive" & Cohort_Short!="B")
leftR <- ANCOMBC_V3V2(ps_leftR) +
  labs(x = NULL, y = NULL, title = "Cohort A, Responsive") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/leftR.jpg", leftR)

ps_leftH <- subset_samples(ps_left, response_patient=="Healthy" & Cohort_Short!="B")
leftH <- ANCOMBC_V3V2(ps_leftH) +
  labs(x = NULL, y = NULL, title = "Cohort A, Healthy") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/leftH.jpg", leftH)

ps_rightNR <- subset_samples(ps_right, response_patient=="nonresponsive" & Cohort_Short!="B")
rightNR <- ANCOMBC_V3V2(ps_rightNR) +
  labs(x = NULL, y = NULL, title = "Cohort B, Nonresponsive") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/rightNR.jpg", rightNR)

ps_rightR <- subset_samples(ps_right, response_patient=="responsive" & Cohort_Short!="B")
rightR <- ANCOMBC_V3V2(ps_rightR) +
  labs(x = NULL, y = NULL, title = "Cohort B, Responsive") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/rightR.jpg", rightR)

ps_rightH <- subset_samples(ps_right, response_patient=="Healthy" & Cohort_Short!="B")
rightH <- ANCOMBC_V3V2(ps_rightH) +
  labs(x = NULL, y = NULL, title = "Cohort B, Healthy") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("ANCOM/rightH.jpg", rightH)