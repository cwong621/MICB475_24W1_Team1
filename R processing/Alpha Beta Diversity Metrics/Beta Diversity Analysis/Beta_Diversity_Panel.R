# Load libraries
library(gridExtra)
library(cowplot)
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggplot2)
library(ggsignif)

# Load in data
load("leftphyloseq_rare.RData")
load("rightphyloseq_rare.RData")

#create list of the samples to remove from the rarified phyloseq object
Samples_toRemove <- c("ART-experienced nonresponsive V2 2", 
                      "ART-experienced nonresponsive V3 3")

#remove those from your phyloseq object
phyloseq_B_noARTnonresponse <- subset_samples(rightphyloseq_rare, !(response_patient_by_visit %in% Samples_toRemove))
phyloseq_A_noARTnonresponse <- subset_samples(leftphyloseq_rare, !(response_patient_by_visit %in% Samples_toRemove))

#save Rdata
save(phyloseq_B_noARTnonresponse, file="phyloseq_B_noARTnonresponse.RData")
save(phyloseq_B_noARTnonresponse, file="phyloseq_B_noARTnonresponsee.RData")

######################## BETA DIVERSITY ###################################

############ bray curtis method ##############

#BC Cohort B
bc_dm_B <- distance(phyloseq_B_noARTnonresponse, method="bray")

pcoa_bc_B <- ordinate(phyloseq_B_noARTnonresponse, method="PCoA", distance=bc_dm_B)

plot_ordination(phyloseq_B_noARTnonresponse, pcoa_bc_B, color = "response_patient_by_visit")

gg_pcoa_B <- plot_ordination(phyloseq_B_noARTnonresponse, pcoa_bc_B, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit", title = "Bray Curtis Cohort B") +
  scale_color_discrete(name = "Response Status by Visit", 
                       labels = c( "Healthy (Week 0)",
                                   "Healthy (Week 24)",
                                   "Nonresponsive (Week 0)",
                                   "Nonresponsive (Week 24)",
                                   "Responsive (Week 0)",
                                   "Responsive (Week 24)")) +
  stat_ellipse()

gg_pcoa_B

ggsave("plot_beta_bray_response_cohort_B.png"
       , gg_pcoa_B
       , height=5, width=7)


#BC Cohort A
bc_dm_A <- distance(phyloseq_A_noARTnonresponse, method="bray")

pcoa_bc_A <- ordinate(phyloseq_A_noARTnonresponse, method="PCoA", distance=bc_dm_A)

plot_ordination(phyloseq_A_noARTnonresponse, pcoa_bc_A, color = "response_patient_by_visit")

gg_pcoa_A <- plot_ordination(phyloseq_A_noARTnonresponse, pcoa_bc_A, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit", title = "Bray Curtis Cohort A") +
  scale_color_discrete(name = "Response Status by Visit", 
                       labels = c( "Healthy (Week 0)",
                                   "Healthy (Week 24)",
                                   "Nonresponsive (Week 0)",
                                   "Nonresponsive (Week 24)",
                                   "Responsive (Week 0)",
                                   "Responsive (Week 24)")) +
  stat_ellipse()

gg_pcoa_A

ggsave("plot_beta_bray_response_cohort_A.png"
       , gg_pcoa_A
       , height=5, width=7)



############ Jaccard #############

#JC Cohort B
jc_dm_B <- distance(phyloseq_B_noARTnonresponse, method="jaccard", binary = TRUE)

pcoa_jc_B <- ordinate(phyloseq_B_noARTnonresponse, method="PCoA", distance=jc_dm_B)

plot_ordination(phyloseq_B_noARTnonresponse, pcoa_jc_B, color = "response_patient_by_visit")

gg_jaccard_pcoa_B <- plot_ordination(phyloseq_B_noARTnonresponse, pcoa_jc_B, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit", title = "Jaccard Cohort B" ) +
  scale_color_discrete(name = "Response Status by Visit", 
                       labels = c( "Healthy (Week 0)",
                                  "Healthy (Week 24)",
                                  "Nonresponsive (Week 0)",
                                  "Nonresponsive (Week 24)",
                                  "Responsive (Week 0)",
                                  "Responsive (Week 24)")) +
  stat_ellipse()

gg_jaccard_pcoa_B

ggsave("plot_beta_jaccard_response_cohort_B.png"
       , gg_jaccard_pcoa_B
       , height=5, width=7)


#JC Cohort A
jc_dm_A <- distance(phyloseq_A_noARTnonresponse, method="jaccard", binary = TRUE)

pcoa_jc_A <- ordinate(phyloseq_A_noARTnonresponse, method="PCoA", distance=jc_dm_A)

plot_ordination(phyloseq_A_noARTnonresponse, pcoa_jc_A, color = "response_patient_by_visit")

gg_jaccard_pcoa_A <- plot_ordination(phyloseq_A_noARTnonresponse, pcoa_jc_A, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit", title = "Jaccard Cohort A" ) +
  scale_color_discrete(name = "Response Status by Visit", 
                       labels = c( "Healthy (Week 0)",
                                   "Healthy (Week 24)",
                                   "Nonresponsive (Week 0)",
                                   "Nonresponsive (Week 24)",
                                   "Responsive (Week 0)",
                                   "Responsive (Week 24)")) +
  stat_ellipse()

gg_jaccard_pcoa_A

ggsave("plot_beta_jaccard_response_cohort_A.png"
       , gg_jaccard_pcoa_A
       , height=5, width=7)


############### Unweighted Unifrac Beta Diversity ##############################

#UU Cohort B
uu_dm_B <- UniFrac(phyloseq_B_noARTnonresponse, weighted = FALSE)

pcoa_uu_B <- ordinate(phyloseq_B_noARTnonresponse, method="PCoA", distance=uu_dm_B)

plot_ordination(phyloseq_B_noARTnonresponse, pcoa_uu_B, color = "response_patient_by_visit")

gg_unweighted_unifrac_pcoa_B <- plot_ordination(phyloseq_B_noARTnonresponse, pcoa_uu_B, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit", title = "Unweighted Unifrac Cohort B") +
  scale_color_discrete(name = "Response Status by Visit", 
                       labels = c( "Healthy (Week 0)",
                                   "Healthy (Week 24)",
                                   "Nonresponsive (Week 0)",
                                   "Nonresponsive (Week 24)",
                                   "Responsive (Week 0)",
                                   "Responsive (Week 24)")) +
  stat_ellipse()

gg_unweighted_unifrac_pcoa_B

ggsave("plot_beta_unweighted_unifrac_response_cohort_B.png"
       , gg_unweighted_unifrac_pcoa_B
       , height=5, width=7)


#UU Cohort A

uu_dm_A <- UniFrac(phyloseq_A_noARTnonresponse, weighted = FALSE)

pcoa_uu_A <- ordinate(phyloseq_A_noARTnonresponse, method="PCoA", distance=uu_dm_A)

plot_ordination(phyloseq_A_noARTnonresponse, pcoa_uu_A, color = "response_patient_by_visit")

gg_unweighted_unifrac_pcoa_A <- plot_ordination(phyloseq_A_noARTnonresponse, pcoa_uu_A, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit", title = " Unweighted Unifrac Cohort A") +
  scale_color_discrete(name = "Response Status by Visit", 
                      labels = c("Healthy (Week 0)",
                                 "Healthy (Week 24)",
                                 "Nonresponsive (Week 0)",
                                "Nonresponsive (Week 24)",
                                "Responsive (Week 0)",
                                "Responsive (Week 24)")) +
  stat_ellipse()

gg_unweighted_unifrac_pcoa_A

ggsave("plot_beta_unweighted_unifrac_response_cohort_A.png"
       , gg_unweighted_unifrac_pcoa_A
       , height=5, width=7)


#################### weighted unifrac beta diversity ###########################

#WU Cohort B
wu_dm_B <- UniFrac(phyloseq_B_noARTnonresponse, weighted = TRUE)

pcoa_wu_B <- ordinate(phyloseq_B_noARTnonresponse, method="PCoA", distance=wu_dm_B)

plot_ordination(phyloseq_B_noARTnonresponse, pcoa_wu_B, color = "response_patient_by_visit")

gg_weighted_unifrac_pcoa_B <- plot_ordination(phyloseq_B_noARTnonresponse, pcoa_wu_B, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit", title = "Weighted Unifrac Cohort B") +
  scale_color_discrete(name = "Response Status by Visit", 
                       labels = c("Healthy (Week 0)",
                                  "Healthy (Week 24)",
                                  "Nonresponsive (Week 0)",
                                  "Nonresponsive (Week 24)",
                                  "Responsive (Week 0)",
                                  "Responsive (Week 24)")) +
  stat_ellipse()

gg_weighted_unifrac_pcoa_B

ggsave("plot_beta_weighted_unifrac_response_cohort_B.png"
       , gg_weighted_unifrac_pcoa_B
       , height=5, width=7)


#WU Cohort A
wu_dm_A <- UniFrac(phyloseq_A_noARTnonresponse, weighted = TRUE)

pcoa_wu_A <- ordinate(phyloseq_A_noARTnonresponse, method="PCoA", distance=wu_dm_A)

plot_ordination(phyloseq_A_noARTnonresponse, pcoa_wu_A, color = "response_patient_by_visit")

gg_weighted_unifrac_pcoa_A <- plot_ordination(phyloseq_A_noARTnonresponse, pcoa_wu_A, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit", title = "Weighted Unifrac Cohort A") +
  scale_color_discrete(name = "Response Status by Visit", 
                       labels = c("Healthy (Week 0)",
                                  "Healthy (Week 24)",
                                  "Nonresponsive (Week 0)",
                                  "Nonresponsive (Week 24)",
                                  "Responsive (Week 0)",
                                  "Responsive (Week 24)")) +
  stat_ellipse()

gg_weighted_unifrac_pcoa_A

ggsave("plot_beta_weighted_unifrac_response_cohort_A.png"
       , gg_weighted_unifrac_pcoa_A
       , height=5, width=7)


############ PANEL FIGURE ############

# extract legend
legend_2 <- cowplot::get_plot_component(gg_unweighted_unifrac_pcoa_A, 'guide-box-right', return_all = TRUE) 
                                         
cowplot::ggdraw(legend_2)

ggsave("legend.png", legend_2, width = 3, height = 2)


# make panel figure
Figure_1_Beta_Diversity <- grid.arrange(
  arrangeGrob(
    gg_unweighted_unifrac_pcoa_A + theme(legend.position = "none", plot.title = element_text(hjust = -0.15)) + labs(title = bquote(bold("A") ~ "  Unweighted Unifrac Cohort A")), 
    gg_unweighted_unifrac_pcoa_B + theme(legend.position = "none", plot.title = element_text(hjust = -0.15)) + labs(title = bquote(bold("B") ~ "  Unweighted Unifrac Cohort B")), 
    gg_weighted_unifrac_pcoa_A + theme(legend.position = "none", plot.title = element_text(hjust = -0.16)) + labs(title = bquote(bold("C") ~ "  Weighted Unifrac Cohort A")), 
    gg_weighted_unifrac_pcoa_B + theme(legend.position = "none", plot.title = element_text(hjust = -0.16)) + labs(title = bquote(bold("D") ~ "  Weighted Unifrac Cohort B")), 
    gg_jaccard_pcoa_A + theme(legend.position = "none", plot.title = element_text(hjust = -0.09)) + labs(title = bquote(bold("E") ~ "  Jaccard Cohort A")),
    gg_jaccard_pcoa_B + theme(legend.position = "none", plot.title = element_text(hjust = -0.09)) + labs(title = bquote(bold("F") ~ "  Jaccard Cohort B")),
    cowplot::ggdraw(legend_2),
    nrow = 4, ncol = 2), nrow = 1,
  heights = unit(c(40), "cm"))

# save as image
ggsave("Figure_1_Beta_Diversity.png", plot = Figure_1_Beta_Diversity, width = 10, height = 16)

