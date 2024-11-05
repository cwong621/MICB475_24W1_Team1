#### Taxa Bar Plot Recreation

library(phyloseq)
library(tidyverse)
# library(picante)
# library(ape)


## Load in .RData ##
load("HIV_rare.RData")

# Plot bar plot of taxonomy
plot_bar(HIV_rare, fill="Phylum")

# Convert to relative abundance
HIV_RA <- transform_sample_counts(HIV_rare, function(x) x/sum(x))

# "glom" by phylum to remove black bars
HIV_phylum <- tax_glom(HIV_RA, taxrank = "Phylum", NArm=FALSE)

## Creating the taxa bar plot for the response of each patient
gg_taxa_rpbv <- plot_bar(HIV_phylum, fill="Phylum") + 
  facet_wrap(.~response_patient_by_visit, scales = "free_x")
gg_taxa_rpbv

# Saving the plot as a .png
ggsave("plot_taxonomy_rpbv.png"
       , gg_taxa_rpbv
       , height=8, width =15)

## Create the taxa bar plot for the starting viral load
## of each patient by visit
gg_taxa_svlpbv <- plot_bar(HIV_phylum, fill="Phylum") + 
  facet_wrap(.~start_viral_load_patient_by_visit, scales = "free_x")
gg_taxa_svlpbv

# Saving the plot as a .png
ggsave("plot_taxonomy_svlpbv.png"
       , gg_taxa_svlpbv
       , height=8, width =12)