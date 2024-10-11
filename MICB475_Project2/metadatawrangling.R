library(tidyverse)
library(reshape2)

# load metadata
meta <- read_delim(file="hiv_metadata.tsv", delim="\t")

#### METADATA WRANGLING ####
#filter out ART-experienced samples and patients with only 1 visit
meta_filt <- filter(meta, `Cohort_Short`!="B") %>%
  group_by(`PID`) %>%
  filter(n()>1) %>%
  ungroup()

# make response column
v2threshold <- 10000 #some sources say below 10000 is low if not on ART (may adjust later)
v3threshold <- 200 #CDC defines <200 as viral suppression under ART

meta_redef <- select(meta_filt,`sample-id`,,`Gender`, `PID`,`Cohort_Short`, `Visit`, `HIV-1_viral_load`, `CRP_mg_L`, `IL-6_pg_mL`,`CD4+_count`,`CD4+_Cat`, `CD4_CD103_Percent`, `CD4_PD1_Percent`, `CD4_HLADRpos_CD38pos_Percent`)%>%
  mutate(`HIV-1_viral_load` = ifelse(grepl("ND", `HIV-1_viral_load`), "0", `HIV-1_viral_load`)) %>% #ND (not detectable) viral loads converted to 0
  mutate(`HIV-1_viral_load` = ifelse(grepl("Negative", `HIV-1_viral_load`), NA, `HIV-1_viral_load`)) %>% #values for HIV negative samples coverted to NA
  mutate(`HIV-1_viral_load` = as.numeric(`HIV-1_viral_load`)) %>%
  mutate(`Visit` = as.factor(`Visit`)) %>%
  mutate(response = case_when(
    `Visit`=="2" & `HIV-1_viral_load`<v2threshold ~ "HIV low", 
    `Visit`=="2" & `HIV-1_viral_load`>=v2threshold ~ "HIV high",
    `Visit`=="3" & `HIV-1_viral_load`<v3threshold ~ "responsive",
    `Visit`=="3" & `HIV-1_viral_load`>=v3threshold ~ "nonresponsive",
    is.na(`HIV-1_viral_load`) ~ "negative"
  ))


# load manifest and filter samples to match redefined metadata
samples_filt <- select(meta_redef, `sample-id`)
manifest_filt <- read_delim(file="hiv_manifest.tsv", delim="\t") %>%
  right_join(samples_filt, join_by(`sample-id`))

#save filtered manifest
write.csv(manifest_filt, file="hiv_manifest_filt.tsv", row.names=FALSE)

#### HEATMAP ####
#make metadata table for heatmap
meta_heat <- group_by(meta_redef, `PID`) %>%
  arrange(`PID`, `response`) %>%
  summarise(response_comb = paste(`response`, collapse = ","), .groups = "drop") %>%
  separate(response_comb, into=c("start_viral_load", "response"), sep=",") %>%
  filter(`start_viral_load`!="negative")

# make heatmap
response_table <- table(meta_heat$response, meta_heat$start_viral_load) 
response_table <- response_table[,c(2,1)]
heatmap <- ggplot(melt(response_table), aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "lightskyblue1", high = "lightskyblue4") +
  labs(x="Starting HIV viral load", y="ART response", fill="Count")

#save heatmap
ggsave(filename = "art_response_heatmap.png"
       , heatmap
       , height=4, width=4)
