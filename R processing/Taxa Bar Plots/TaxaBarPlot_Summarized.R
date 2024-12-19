#### Taxa Bar Plot Recreation

library(phyloseq)
library(tidyverse)

## Load in .RData ##
load("leftphyloseq.RData")
load("rightphyloseq.RData")

# Convert to relative abundance
left_RA <- transform_sample_counts(ps_left, function(x) x/sum(x))
right_RA <- transform_sample_counts(ps_right, function(x) x/sum(x))

# "glom" by phylum to remove black bars
left_phylum <- tax_glom(left_RA, taxrank = "Phylum", NArm=FALSE)
right_phylum <- tax_glom(right_RA, taxrank = "Phylum", NArm=FALSE)

left_melt <- psmelt(left_phylum) 
left_summ <- left_melt %>%
  group_by(response_patient_by_visit, Phylum) %>%
  summarise_at("Abundance", mean)

## Creating the taxa bar plot for the response of each patient

gg_taxa_rpbv_left <- ggplot(left_summ, aes(x=response_patient_by_visit, y=Abundance, fill=Phylum)) +
  geom_bar(color="black", stat="identity", position="fill") +
  theme(axis.text.x = element_text(angle=45, hjust = 1))
gg_taxa_rpbv_left

# Saving the plot as a .png
ggsave("plot_taxonomy_rpbv_left_summary.png"
       , gg_taxa_rpbv_left
       , height=8, width =15)


## Right
right_melt <- psmelt(right_phylum) 
right_summ <- right_melt %>%
  group_by(response_patient_by_visit, Phylum) %>%
  summarise_at("Abundance", mean)

## Creating the taxa bar plot for the response of each patient

gg_taxa_rpbv_right <- ggplot(right_summ, aes(x=response_patient_by_visit, y=Abundance, fill=Phylum)) +
  geom_bar(color="black", stat="identity", position="fill")+
  theme(axis.text.x = element_text(angle=45, hjust = 1))
gg_taxa_rpbv_right

# Saving the plot as a .png
ggsave("plot_taxonomy_rpbv_right_summary.png"
       , gg_taxa_rpbv_right
       , height=8, width =15)
