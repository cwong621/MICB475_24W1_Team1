library(tidyverse)
library(reshape2)

# load metadata
meta <- read_delim(file="hiv_metadata.tsv", delim="\t")

#### METADATA WRANGLING ####

# define thresholds
v2threshold <- 25000 #some sources say below 10000 is low if not on ART (may adjust later)
v3threshold <- 200 #CDC defines <200 as viral suppression under ART

# clean HIV-1 viral load column and create response column
meta_redef <- mutate(meta, `HIV-1_viral_load` = ifelse(grepl("ND", `HIV-1_viral_load`), "0", `HIV-1_viral_load`)) %>% #ND (not detectable) viral loads converted to 0
  mutate(`HIV-1_viral_load` = ifelse(grepl("Negative", `HIV-1_viral_load`), NA, `HIV-1_viral_load`)) %>% #values for HIV negative samples coverted to NA
  mutate(`HIV-1_viral_load` = as.numeric(`HIV-1_viral_load`)) %>%
  filter(`Cohort_Short`!="B" | `HIV-1_viral_load`>=v3threshold) %>% #remove ART-experienced individuals responsive to therapy
  group_by(`PID`) %>%
  filter(n()>1) %>% #remove individuals with only one visit
  ungroup() %>%
  mutate(response = case_when( #create response column 
    `Visit`==2 & `HIV-1_viral_load`<v2threshold & `Cohort_Short`=="A" ~ "HIV low pretreatment", 
    `Visit`==2 & `HIV-1_viral_load`>=v2threshold & `Cohort_Short`=="A" ~ "HIV high pretreatment",
    `Visit`==3 & `HIV-1_viral_load`<v3threshold & `Cohort_Short`=="A" ~ "responsive",
    `Visit`==3 & `HIV-1_viral_load`>=v3threshold & `Cohort_Short`=="A" ~ "nonresponsive",
    `Visit`==3 & `Cohort_Short`=="B" ~ "ART-experienced nonresponsive V3",
    `Visit`==2 & `Cohort_Short`=="B" ~ "ART-experienced nonresponsive V2",
    is.na(`HIV-1_viral_load`) ~ "Healthy"
  ))%>%
  mutate(IL6_Cat = case_when( #create column categorizing IL6 level
    `IL-6_pg_mL`<=5 ~ "normal",
    `IL-6_pg_mL`>5 ~ "high"
  ))%>%
  mutate(CRP_Cat = case_when( #create column categorizing CRP level
    `CRP_mg_L`<=3 ~ "normal",
    `CRP_mg_L`>3 ~ "high"
  ))%>%
  mutate(response_patient = case_when( #make response column that assigns response value to both timepoints
    `Visit`==3 & `HIV-1_viral_load`<v3threshold & `Cohort_Short`=="A" ~ "responsive",
    `Visit`==3 & `HIV-1_viral_load`>=v3threshold & `Cohort_Short`=="A" ~ "nonresponsive",
    `Visit`==3 & `Cohort_Short`=="B" ~ "ART-experienced nonresponsive V3",
    `Visit`==2 & `Cohort_Short`=="B" ~ "ART-experienced nonresponsive V2",
    is.na(`HIV-1_viral_load`) ~ "Healthy"
  ))  %>%
  mutate(start_viral_load = case_when( #make column for starting viral load high/low
    `Visit`==2 & `HIV-1_viral_load`<v2threshold & `Cohort_Short`=="A" ~ "HIV low pretreatment", 
    `Visit`==2 & `HIV-1_viral_load`>=v2threshold & `Cohort_Short`=="A" ~ "HIV high pretreatment",
    `Visit`==2 & `Cohort_Short`=="B" ~ NA,
    is.na(`HIV-1_viral_load`) ~ "Healthy"
  )) %>%
  group_by(`PID`) %>% 
  mutate(response_patient = ifelse(is.na(response_patient), na.locf(response_patient), response_patient)) %>% # assign response value to both timepoints
  mutate(start_viral_load = ifelse(is.na(start_viral_load), na.locf(start_viral_load), start_viral_load)) %>% # assign start viral load value to both timepoints
  ungroup() %>%
  unite(response_patient_by_visit, response_patient, Visit, sep=" ", remove= F )%>% 
  unite(start_viral_load_patient_by_visit, start_viral_load, Visit, sep=" ", remove = F)

# save modified metadata file to uplaod to server
write.table(meta_redef, file="hiv_metadata_filt.tsv", quote=FALSE, sep="\t", row.names=FALSE)


# load manifest and filter samples to match redefined metadata
manifest_filt <- read_delim(file="hiv_manifest.tsv", delim="\t") %>%
  subset(`sample-id` %in% meta_redef$`sample-id`)

#save filtered manifest to upload to server
write.table(manifest_filt, file="hiv_manifest_filt.tsv", quote=FALSE, sep="\t", row.names=FALSE)

#### HEATMAP ####
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
heatmap_crp <- ggplot(melt(response_location), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="Location", y="ART response", fill="Count")

#cotrimoxazole use vs response (trend not significant)
response_cotri <- table(meta_heat$response_patient, meta_heat$`Cotrimoxazole_HIV`)
heatmap_cotri <- ggplot(melt(response_cotri), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="Cotrimoxazole status", y="ART response", fill="Count")

#save heatmap
ggsave(filename = "viralload_response_heatmap.png"
       , heatmap
       , height=4, width=5)
