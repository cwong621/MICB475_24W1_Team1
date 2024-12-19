library(tidyverse)
library(zoo)

# load metadata
meta <- read_delim(file="./hiv_metadata.tsv", delim="\t")

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
save(meta_redef, file="meta_redef.RData")

# load manifest and filter samples to match redefined metadata
manifest_filt <- read_delim(file="hiv_manifest.tsv", delim="\t") %>%
  subset(`sample-id` %in% meta_redef$`sample-id`)

#save filtered manifest to upload to server
write.table(manifest_filt, file="hiv_manifest_filt.tsv", quote=FALSE, sep="\t", row.names=FALSE)
