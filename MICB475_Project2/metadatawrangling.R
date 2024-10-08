library(tidyverse)

# load metadata
meta <- read_delim(file="hiv_metadata.tsv", delim="\t")

#filter out ART-experienced samples and patients with only 1 visit
meta_filt <- filter(meta, `Cohort_Short`!="B") %>%
  group_by(`PID`) %>%
  filter(n()>1)

# make response column

v2threshold <- 10000
v3threshold <- 200

meta_redef <- select(meta_filt,`sample-id`, `PID`,`Cohort_Short`, `Visit`, `HIV-1_viral_load`)%>%
  mutate(`HIV-1_viral_load` = ifelse(grepl("ND", `HIV-1_viral_load`), "0", `HIV-1_viral_load`)) %>% #ND viral loads converted to 0
  mutate(`HIV-1_viral_load` = ifelse(grepl("Negative", `HIV-1_viral_load`), NA, `HIV-1_viral_load`)) %>%
  mutate(`HIV-1_viral_load` = as.numeric(`HIV-1_viral_load`)) %>%
  mutate(`Visit` = as.factor(`Visit`)) %>%
  mutate(response = case_when(
    `Visit`=="2" & `HIV-1_viral_load`<v2threshold ~ "HIVlow",
    `Visit`=="2" & `HIV-1_viral_load`>v2threshold ~ "HIVhigh",
    `Visit`=="3" & `HIV-1_viral_load`<v3threshold ~ "responsive",
    `Visit`=="3" & `HIV-1_viral_load`>v3threshold ~ "nonresponsive",
    is.na(`HIV-1_viral_load`) ~ "negative"
  ))

meta_heat <- filter(meta_redef, `Cohort_Short`=="A") %>%
  select(`PID`,`response`)
heatmap(meta_heat)

