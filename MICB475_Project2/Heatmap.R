library(tidyverse)
library(reshape2)

#### Load data ####
load("./meta_redef.RData")

#make metadata table for heatmap (only ART naive samples from pretreatment timepoint)
meta_heat <- filter(meta_redef, `Cohort_Short`=="A") %>%
  filter(`Visit`==2)

#viral load vs response (trend but not significant)
response_load <- table(meta_heat$response_patient, meta_heat$start_viral_load)
response_load <- response_load[,c(2,1)]
heatmap_load <- ggplot(melt(response_load), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="Starting HIV viral load", y="ART response", fill="Count")

#cd4 vs response **significant
response_cd4 <- table(meta_heat$response_patient, meta_heat$`CD4+_Cat`)
response_cd4 <- response_cd4[,c(1,3,2)]
heatmap_cd4 <- ggplot(melt(response_cd4), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="CD4 count", y="ART response", fill="Count")

#sex vs response (no difference)
response_sex <- table(meta_heat$response_patient, meta_heat$`Gender`)
heatmap_sex <- ggplot(melt(response_sex), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="Sex", y="ART response", fill="Count")

#bmi vs response (no clear difference)
response_bmi <- table(meta_heat$response_patient, meta_heat$`BMI_Cat`)
response_bmi <- response_bmi[,c(5,1,2,4,3)]
heatmap_bmi <- ggplot(melt(response_bmi), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="BMI category", y="ART response", fill="Count")

#il6 vs response (no difference)
response_il6 <- table(meta_heat$response_patient, meta_heat$`IL6_Cat`)
heatmap_il6 <- ggplot(melt(response_il6), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="IL-6 category", y="ART response", fill="Count")

#crp vs response (no difference)
response_crp <- table(meta_heat$response_patient, meta_heat$`CRP_Cat`)
heatmap_crp <- ggplot(melt(response_crp), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="CRP", y="ART response", fill="Count")

#location vs response
response_location <- table(meta_heat$response_patient, meta_heat$`Arm_Short`)
heatmap_location <- ggplot(melt(response_location), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="Location", y="ART response", fill="Count")

#cotrimoxazole use vs response (trend not significant)
response_cotri <- table(meta_heat$response_patient, meta_heat$`Cotrimoxazole_HIV`)
heatmap_cotri <- ggplot(melt(response_cotri), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="Cotrimoxazole status", y="ART response", fill="Count")

#save heatmaps
ggsave(filename = "./Heatmaps/viralload_response_heatmap.png"
       , heatmap_load
       , height=4, width=5)
ggsave(filename = "./Heatmaps/cd4_response_heatmap.png"
       , heatmap_cd4
       , height=4, width=5)
ggsave(filename = "./Heatmaps/cotri_response_heatmap.png"
       , heatmap_cotri
       , height=4, width=5)
ggsave(filename = "./Heatmaps/location_response_heatmap.png"
       , heatmap_location
       , height=4, width=5)
ggsave(filename = "./Heatmaps/crp_response_heatmap.png"
       , heatmap_crp
       , height=4, width=5)
ggsave(filename = "./Heatmaps/il6_response_heatmap.png"
       , heatmap_il6
       , height=4, width=5)
ggsave(filename = "./Heatmaps/sex_response_heatmap.png"
       , heatmap_sex
       , height=4, width=5)
ggsave(filename = "./Heatmaps/bmi_response_heatmap.png"
       , heatmap_bmi
       , height=4, width=5)