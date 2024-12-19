#### Taxa Bar Plot Recreation

library(phyloseq)
library(tidyverse)

## Load in .RData ##
load("leftphyloseq.RData")
load("rightphyloseq.RData")


# Plot bar plot of taxonomy
plot_bar(ps_left, fill="Phylum")

plot_bar(ps_right, fill="Phylum")

# Convert to relative abundance
left_RA <- transform_sample_counts(ps_left, function(x) x/sum(x))
right_RA <- transform_sample_counts(ps_right, function(x) x/sum(x))

# "glom" by phylum to remove black bars
left_phylum <- tax_glom(left_RA, taxrank = "Phylum", NArm=FALSE)
right_phylum <- tax_glom(right_RA, taxrank = "Phylum", NArm=FALSE)

## Creating the taxa bar plot for the response of each patient
gg_taxa_rpbv_left <- plot_bar(left_phylum, fill="Phylum") + 
  facet_wrap(.~response_patient_by_visit, scales = "free_x")
gg_taxa_rpbv_left

# Saving the plot as a .png
ggsave("plot_taxonomy_rpbv_left.png"
       , gg_taxa_rpbv_left
       , height=8, width =15)

gg_taxa_rpbv_right <- plot_bar(right_phylum, fill="Phylum") + 
  facet_wrap(.~response_patient_by_visit, scales = "free_x")
gg_taxa_rpbv_right

# Saving the plot as a .png
ggsave("plot_taxonomy_rpbv_right.png"
       , gg_taxa_rpbv_right
       , height=8, width =15)

## Create the taxa bar plot for the starting viral load
## of each patient by visit
gg_taxa_svlpbv_left <- plot_bar(left_phylum, fill="Phylum") + 
  facet_wrap(.~start_viral_load_patient_by_visit, scales = "free_x")
gg_taxa_svlpbv_left

# Saving the plot as a .png
ggsave("plot_taxonomy_svlpbv_left.png"
       , gg_taxa_svlpbv_left
       , height=8, width =15)

gg_taxa_svlpbv_right <- plot_bar(right_phylum, fill="Phylum") + 
  facet_wrap(.~start_viral_load_patient_by_visit, scales = "free_x")
gg_taxa_svlpbv_right

# Saving the plot as a .png
ggsave("plot_taxonomy_svlpbv_right.png"
       , gg_taxa_svlpbv_right
       , height=8, width =15)
