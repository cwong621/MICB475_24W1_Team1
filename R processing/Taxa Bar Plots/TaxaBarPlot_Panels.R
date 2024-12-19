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

# Renaming the response categories
left_summ_renamed <- mutate(left_summ, 
                        response_patient_by_visit = recode(response_patient_by_visit,
                                             "ART-experienced nonresponsive V2 2" = "ART-Experienced Nonresponsive (Week 0)",
                                             "ART-experienced nonresponsive V3 3" = "ART-Experienced Nonresponsive (Week 24)",
                                             "Healthy 2" = "Healthy (Week 0)",
                                             "Healthy 3" = "Healthy (Week 24)",
                                             "nonresponsive 2" = "Nonresponsive (Week 0)",
                                             "nonresponsive 3" = "Nonresponsive (Week 24)",
                                             "responsive 2" = "Responsive (Week 0)",
                                             "responsive 3" = "Responsive (Week 24)") )
## Filtering out ART-experienced samples
left_summ_renamed_noARTexpV2 <- filter(left_summ_renamed, response_patient_by_visit != "ART-Experienced Nonresponsive (Week 0)")
left_summ_renamed_noARTexp <- filter(left_summ_renamed_noARTexpV2, response_patient_by_visit != "ART-Experienced Nonresponsive (Week 24)")


## Creating the taxa bar plot for the response of each patient

gg_taxa_rpbv_left <- ggplot(left_summ_renamed_noARTexp, aes(x=response_patient_by_visit, y=Abundance, fill=Phylum)) +
  geom_bar(color="black", stat="identity", position="fill") +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  xlab(NULL)
gg_taxa_rpbv_left

# Saving the plot as a .png
ggsave("plot_taxonomy_rpbv_left_summary.png"
       , gg_taxa_rpbv_left
       , height=8, width =8)


## Right
right_melt <- psmelt(right_phylum) 
right_summ <- right_melt %>%
  group_by(response_patient_by_visit, Phylum) %>%
  summarise_at("Abundance", mean)


# Renaming the response categories again
right_summ_renamed <- mutate(right_summ, 
                            response_patient_by_visit = recode(response_patient_by_visit,
                                                               "ART-experienced nonresponsive V2 2" = "ART-Experienced Nonresponsive (Week 0)",
                                                               "ART-experienced nonresponsive V3 3" = "ART-Experienced Nonresponsive (Week 24)",
                                                               "Healthy 2" = "Healthy (Week 0)",
                                                               "Healthy 3" = "Healthy (Week 24)",
                                                               "nonresponsive 2" = "Nonresponsive (Week 0)",
                                                               "nonresponsive 3" = "Nonresponsive (Week 24)",
                                                               "responsive 2" = "Responsive (Week 0)",
                                                               "responsive 3" = "Responsive (Week 24)") )
# Filtering out ART-experienced samples
right_summ_renamed_noARTexpV2 <- filter(right_summ_renamed, response_patient_by_visit != "ART-Experienced Nonresponsive (Week 0)")
right_summ_renamed_noARTexp <- filter(right_summ_renamed_noARTexpV2, response_patient_by_visit != "ART-Experienced Nonresponsive (Week 24)")


## Creating the taxa bar plot for the response of each patient

gg_taxa_rpbv_right <- ggplot(right_summ_renamed_noARTexp, aes(x=response_patient_by_visit, y=Abundance, fill=Phylum)) +
  geom_bar(color="black", stat="identity", position="fill")+
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  xlab(NULL)
gg_taxa_rpbv_right

# Saving the plot as a .png
ggsave("plot_taxonomy_rpbv_right_summary.png"
       , gg_taxa_rpbv_right
       , height=8, width =8)


## Making the plots into one panel -> the right panel will be squished

taxabar_paneled <- grid.arrange(
  arrangeGrob(
    gg_taxa_rpbv_left + theme(legend.position = "none", plot.title = element_text(hjust = 0)) + labs(title = bquote(bold("A"))), 
    gg_taxa_rpbv_right + theme(legend.position = "right", plot.title element_text(hjust = 0)) + labs(title = bquote(bold("B"))),
    nrow = 2, ncol = 2), nrow = 2,
  heights = c(15, 1),
  widths = c(10, 1)# Adjust heights to give more space to plots
)

taxabar_paneled

ggsave(filename = "taxabar_paneled.png"
       , taxabar_paneled
       , height=15, width=10
       , )
