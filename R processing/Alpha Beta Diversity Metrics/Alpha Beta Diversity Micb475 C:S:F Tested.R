# Load libraries

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggplot2)
library(rbiom)
library(ggsignif)

# Load in data

load("HIV_rare.RData")
load("leftphyloseq_rare.RData")
load("rightphyloseq_rare.RData")

# response_patient_by_visit

# Alpha Diversity

plot_right <- plot_richness(rightphyloseq_rare)

plot_richness(rightphyloseq_rare, measures = c("Shannon","Chao1")) 


plot_left <- plot_richness(leftphyloseq_rare)



alpha_response_richness_right <- plot_richness(rightphyloseq_rare, x = "response_patient_by_visit", measures = c("Shannon","Chao1")) +
  xlab("Response Status by Visit") +
  geom_boxplot()
alpha_response_richness_right

alpha_response_richness_left <- plot_richness(leftphyloseq_rare, x = "response_patient_by_visit", measures = c("Shannon","Chao1")) +
  xlab("Response Status by Visit") +
  geom_boxplot()
alpha_response_richness_left


ggsave(filename = "chao1_shannon_alpha_response_richness_right.png"
       , alpha_response_richness_right
       , height=5, width=6)

ggsave(filename = "chao1_shannon_alpha_response_richness_left.png"
       , alpha_response_richness_left
       , height=5, width=6)


# alpha_viral_load_richness <- plot_richness(HIV_rare, x = "start_viral_load_patient_by_visit", measures = c("Shannon","Chao1")) +
  # xlab("Initial Viral Load by Visit") +
  # geom_boxplot()
#alpha_viral_load_richness


#ggsave(filename = "plot_alpha_viral_load_richness.png"
#       , alpha_viral_load_richness
#       , height=5, width=6)


alpha_right <- estimate_richness(rightphyloseq_rare)
alpha_left <- estimate_richness(leftphyloseq_rare)

# Getting the phyloseq objects into dataframe formats for future use
leftphylo_dat <- sample_data(leftphyloseq_rare)
leftphylo_datfrm <- data.frame(leftphylo_dat)

rightphylo_dat <- sample_data(rightphyloseq_rare)
rightphylo_datfrm <- data.frame(rightphylo_dat)

# Making dataframes including richness measures for Shannon/Chao1 testing
leftphylo_datfrm_wdiv <- data.frame(leftphylo_datfrm, alpha_left)
rightphylo_datfrm_wdiv <- data.frame(rightphylo_datfrm, alpha_right)


### Statistical testing of Shannon/Chao1 - Kruskal-Wallis ###

# Chao1 testing
kruskal_chao1_left <- kruskal.test(Chao1 ~ `response_patient_by_visit`, data = leftphylo_datfrm_wdiv)
kruskal_chao1_left  # something significant here!

kruskal_chao1_right <- kruskal.test(Chao1 ~ `response_patient_by_visit`, data = rightphylo_datfrm_wdiv)
kruskal_chao1_right  # nothing significant here

# Shannon testing
kruskal_shannon_left <- kruskal.test(Shannon ~ `response_patient_by_visit`, data = leftphylo_datfrm_wdiv)
kruskal_shannon_left  # something significant here!

kruskal_shannon_right <- kruskal.test(Shannon ~ `response_patient_by_visit`, data = rightphylo_datfrm_wdiv)
kruskal_shannon_right  # nothing significant here


### Log-transforming and performing ANOVAs on left cohort measures - as was described in course material

# Chao1
left_chao1_vs_rpbv_logtrsfm <- lm(log(Chao1) ~ `response_patient_by_visit`, data=leftphylo_datfrm_wdiv)
anova_left_chao1_vs_rpbv_logtrsfm <- aov(left_chao1_vs_rpbv_logtrsfm)
summary(anova_left_chao1_vs_rpbv_logtrsfm)

# Seeing which comparison is significant
TukeyHSD(anova_left_chao1_vs_rpbv_logtrsfm) # nonresponsive 3 vs both Healthy!


# Shannon
left_shannon_vs_rpbv_logtrsfm <- lm(log(Shannon) ~ `response_patient_by_visit`, data=leftphylo_datfrm_wdiv)
anova_left_shannon_vs_rpbv_logtrsfm <- aov(left_shannon_vs_rpbv_logtrsfm)
summary(anova_left_shannon_vs_rpbv_logtrsfm)

# Seeing which comparison is significant
TukeyHSD(anova_left_shannon_vs_rpbv_logtrsfm) # again, nonresponsive 3 vs both Healthy!



## Re-plotting left cohort Chao1 and Shannon plots to show significant comparisons

# Filtering the data to not include the ART-experienced nonresponsive cohort
leftphylo_datfrm_wdiv_noARTexpV2 <- filter(leftphylo_datfrm_wdiv, response_patient_by_visit != "ART-experienced nonresponsive V2 2")
leftphylo_datfrm_wdiv_noARTexp <- filter(leftphylo_datfrm_wdiv_noARTexpV2, response_patient_by_visit != "ART-experienced nonresponsive V3 3")

# Chao1
chao1_reponse_left_plot_sig <- ggplot(leftphylo_datfrm_wdiv_noARTexp, aes(response_patient_by_visit, Chao1, ylim =40)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("nonresponsive 3","Healthy 2"), c("nonresponsive 3", "Healthy 3")),
              y_position = c(500,550),
              annotations = c("*", "*")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Left Cohort Response Status") +
  ylim(100, 600) + # tweak axis to include significance symbols
  ylab("Alpha Diversity Measure")

chao1_reponse_left_plot_sig

# Shannon
shannon_reponse_left_plot_sig <- ggplot(leftphylo_datfrm_wdiv_noARTexp, aes(response_patient_by_visit, Shannon, ylim =40)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("nonresponsive 3","Healthy 2"), c("nonresponsive 3", "Healthy 3")),
              y_position = c(5.5,6),
              annotations = c("*", "*")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Left Cohort Response Status") +
  ylim(3,6.2) + # tweak axis to include significance symbols
  ylab("Alpha Diversity Measure")

shannon_reponse_left_plot_sig

### Regenerating right cohort plots to have the same axes/aesthetics

# Filtering the data to not include the ART-experienced nonresponsive cohort
rightphylo_datfrm_wdiv_noARTexpV2 <- filter(rightphylo_datfrm_wdiv, response_patient_by_visit != "ART-experienced nonresponsive V2 2")
rightphylo_datfrm_wdiv_noARTexp <- filter(rightphylo_datfrm_wdiv_noARTexpV2, response_patient_by_visit != "ART-experienced nonresponsive V3 3")

# Chao1
chao1_reponse_right_plot_sig <- ggplot(rightphylo_datfrm_wdiv_noARTexp, aes(response_patient_by_visit, Chao1, ylim =40)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Right Cohort Response Status") +
  ylab("Alpha Diversity Measure")

chao1_reponse_right_plot_sig

# Shannon
shannon_reponse_right_plot_sig <- ggplot(rightphylo_datfrm_wdiv_noARTexp, aes(response_patient_by_visit, Shannon, ylim =40)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Right Cohort Response Status") +
  #ylim(3,6.2) + # tweak axis to include significance symbols
  ylab("Alpha Diversity Measure")

shannon_reponse_right_plot_sig

## Save all the Chao1/Shannon graphs as images

# Left cohorts
ggsave(filename = "chao1_reponse_left_plot_sig.png"
       , chao1_reponse_left_plot_sig
       , height=5, width=6)

ggsave(filename = "shannon_reponse_left_plot_sig.png"
       , shannon_reponse_left_plot_sig
       , height=5, width=6)

# Right cohorts
ggsave(filename = "chao1_reponse_right_plot_sig.png"
       , chao1_reponse_right_plot_sig
       , height=5, width=6)

ggsave(filename = "shannon_reponse_right_plot_sig.png"
       , shannon_reponse_right_plot_sig
       , height=5, width=6)






##### Phylogenetic diversity  #####

# calculate Faith's phylogenetic diversity as PD
phylo_dist_right <- pd(t(otu_table(rightphyloseq_rare)), phy_tree(rightphyloseq_rare),
                 include.root=F) 

phylo_dist_left <- pd(t(otu_table(leftphyloseq_rare)), phy_tree(leftphyloseq_rare),
                  include.root=F) 

# add PD to metadata table
sample_data(rightphyloseq_rare)$PD <- phylo_dist_right$PD

sample_data(leftphyloseq_rare)$PD <- phylo_dist_left$PD


# plot any metadata category against the PD
faith_reponse_right_plot.pd <- ggplot(sample_data(rightphyloseq_rare), aes(response_patient_by_visit, PD)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Right Cohort Response Status") +
  ylab("Phylogenetic Diversity")

faith_reponse_left_plot.pd <- ggplot(sample_data(leftphyloseq_rare), aes(response_patient_by_visit, PD)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Left Cohort Response Status") +
  ylab("Phylogenetic Diversity")

# view plots
faith_reponse_right_plot.pd
faith_reponse_left_plot.pd

######### Run Kruskal-Wallis statistical testing #########


# Testing for significance in at least one comparison for each cohort
kruskal_phylodiv_left <- kruskal.test( PD ~ response_patient_by_visit, data = leftphylo_datfrm)
kruskal_phylodiv_left   ## there is something significant here!

kruskal_phylodiv_right <- kruskal.test( PD ~ response_patient_by_visit, data = rightphylo_datfrm)
kruskal_phylodiv_right  ## nothing significant here

### Log-transforming and performing an ANOVA on the left cohort - as was described in course material

left_pd_vs_rpbv_logtrsfm <- lm(log(PD) ~ `response_patient_by_visit`, data=leftphylo_datfrm)
anova_left_pd_vs_rpbv_logtrsfm <- aov(left_pd_vs_rpbv_logtrsfm)
summary(anova_left_pd_vs_rpbv_logtrsfm)

# Seeing which comparison is significant
TukeyHSD(anova_left_pd_vs_rpbv_logtrsfm)

### Generating polished plots with statistical test annotations

# Filtering the data to not include the ART-experienced nonresponsive cohort
leftphylo_datfrm_noARTexpV2 <- filter(leftphylo_datfrm, response_patient_by_visit != "ART-experienced nonresponsive V2 2")
leftphylo_datfrm_noARTexp <- filter(leftphylo_datfrm_noARTexpV2, response_patient_by_visit != "ART-experienced nonresponsive V3 3")

faith_reponse_left_plot.pd <- ggplot(leftphylo_datfrm_noARTexp, aes(response_patient_by_visit, PD, ylim =40)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("nonresponsive 3","Healthy 2")),
              y_position = c(30),
              annotations = c("*")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Left Cohort Response Status") +
  ylim(c(10,32)) + #tweaked the axis a bit since the asterisk was getting cut off
  ylab("Phylogenetic Diversity")

faith_reponse_left_plot.pd

# Filtering the right plot to not include the ART-experienced group as well
rightphylo_datfrm_noARTexpV2 <- filter(rightphylo_datfrm, response_patient_by_visit != "ART-experienced nonresponsive V2 2")
rightphylo_datfrm_noARTexp <- filter(rightphylo_datfrm_noARTexpV2, response_patient_by_visit != "ART-experienced nonresponsive V3 3")

# Regenerating right plot without that ART-experienced group
faith_reponse_right_plot.pd <- ggplot(rightphylo_datfrm_noARTexp, aes(response_patient_by_visit, PD)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Right Cohort Response Status") +
  ylab("Phylogenetic Diversity")

faith_reponse_right_plot.pd

############

# save graphs as images
ggsave(filename = "plot_faith_reponse_richness_right.png"
       , faith_reponse_right_plot.pd
       , height=5, width=6)

ggsave(filename = "plot_faith_reponse_richness_left.png"
       , faith_reponse_left_plot.pd
       , height=5, width=6
       , )


###### BETA DIVERSITY #########

############ bray curtis method ##############
bc_dm_right <- distance(rightphyloseq_rare, method="bray")

pcoa_bc <- ordinate(rightphyloseq_rare, method="PCoA", distance=bc_dm_right)

plot_ordination(rightphyloseq_rare, pcoa_bc, color = "response_patient_by_visit")

gg_pcoa_right <- plot_ordination(rightphyloseq_rare, pcoa_bc, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit") +
  stat_ellipse()
gg_pcoa_right

ggsave("plot_beta_bray_response_right_pcoa_2.png"
       , gg_pcoa_right
       , height=5, width=7)


# left
bc_dm_left <- distance(leftphyloseq_rare, method="bray")

pcoa_bc_left <- ordinate(leftphyloseq_rare, method="PCoA", distance=bc_dm_left)

plot_ordination(leftphyloseq_rare, pcoa_bc_left, color = "response_patient_by_visit", shape = "start_viral_load_patient_by_visit")

gg_pcoa_left <- plot_ordination(leftphyloseq_rare, pcoa_bc_left, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit") +
  stat_ellipse()
gg_pcoa_left

ggsave("plot_beta_bray_response_left_pcoa_2.png"
       , gg_pcoa_left
       , height=5, width=7)


############ Jaccard #############
jc_dm_right <- distance(rightphyloseq_rare, method="jaccard", binary = TRUE)

pcoa_jc_right <- ordinate(rightphyloseq_rare, method="PCoA", distance=jc_dm_right)

plot_ordination(rightphyloseq_rare, pcoa_jc_right, color = "response_patient_by_visit")

gg_jaccard_pcoa_right <- plot_ordination(rightphyloseq_rare, pcoa_jc_right, color = "response_patient_by_visit", shape="start_viral_load_patient_by_visit") +
  labs(col = "Response Status by Visit", shape = "start_viral_load_patient_by_visit" )
gg_jaccard_pcoa_right

ggsave("plot_beta_jaccard_response_right_pcoa.png"
       , gg_jaccard_pcoa_right
       , height=5, width=7)

# left cohort
jc_dm_left <- distance(leftphyloseq_rare, method="jaccard", binary = TRUE)

pcoa_jc_left <- ordinate(leftphyloseq_rare, method="PCoA", distance=jc_dm_left)

plot_ordination(leftphyloseq_rare, pcoa_jc_left, color = "response_patient_by_visit")

gg_jaccard_pcoa_left <- plot_ordination(leftphyloseq_rare, pcoa_jc_left, color = "response_patient_by_visit", shape="start_viral_load_patient_by_visit") +
  labs(col = "Response Status by Visit", shape = "start_viral_load_patient_by_visit" )
gg_jaccard_pcoa_left

ggsave("plot_beta_jaccard_response_left_pcoa.png"
       , gg_jaccard_pcoa_left
       , height=5, width=7)



############### Unweighted Unifrac Beta Diversity ##############################
?distance

# whole cohort

uu_dm_whole <- distance(HIV_rare, method = "uunifrac")

uu_distance <- UniFrac(HIV_rare, weighted = FALSE, normalized = TRUE, parallel = FALSE)
# both of the codes abovegive the same output and error

pcoa_uu_whole <- ordinate(HIV_rare, method="PCoA", distance=uu_dm_whole)

plot_ordination(HIV_rare, pcoa_uu_whole, color = "response_patient_by_visit")

gg_unweighted_unifrac_pcoa_whole <- plot_ordination(HIV_rare, pcoa_uu_whole, color = "response_patient_by_visit", axes = c(1,2)) +
  labs(col = "Response Status by Visit") +
  stat_ellipse()

gg_unweighted_unifrac_pcoa_whole

ggsave("plot_beta_unweighted_unifrac_response_whole_pcoa.png"
       , gg_unweighted_unifrac_pcoa_whole
       , height=5, width=7)


#right cohort
uu_dm_right <- distance(rightphyloseq_rare, "uunifrac")

pcoa_uu_right <- ordinate(rightphyloseq_rare, method="PCoA", distance=uu_dm_right)

plot_ordination(rightphyloseq_rare, pcoa_uu_right, color = "response_patient_by_visit")

gg_unweighted_unifrac_pcoa_right <- plot_ordination(rightphyloseq_rare, pcoa_uu_right, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit") +
  stat_ellipse()

gg_unweighted_unifrac_pcoa_right

ggsave("plot_beta_unweighted_unifrac_response_right_pcoa.png"
       , gg_unweighted_unifrac_pcoa_right
       , height=5, width=7)


#left cohort
uu_dm_left <- distance(leftphyloseq_rare, method="uunifrac")

pcoa_uu_left <- ordinate(leftphyloseq_rare, method="PCoA", distance=uu_dm_left)

plot_ordination(leftphyloseq_rare, pcoa_uu_left, color = "response_patient_by_visit")

gg_unweighted_unifrac_pcoa_left <- plot_ordination(leftphyloseq_rare, pcoa_uu_left, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit") +
  stat_ellipse()

gg_unweighted_unifrac_pcoa_left

ggsave("plot_beta_unweighted_unifrac_response_left_pcoa.png"
       , gg_unweighted_unifrac_pcoa_left
       , height=5, width=7)





#################### weighted unifrac beta diversity ###########################
# whole
wu_dm_whole <- distance(HIV_rare, method="wunifrac")

pcoa_wu_whole <- ordinate(HIV_rare, method="PCoA", distance=wu_dm_whole)

plot_ordination(HIV_rare, pcoa_wu_whole, color = "response_patient_by_visit")

gg_weighted_unifrac_pcoa_whole <- plot_ordination(HIV_rare, pcoa_wu_whole, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit") +
  stat_ellipse()

gg_weighted_unifrac_pcoa_whole

ggsave("plot_beta_weighted_unifrac_response_whole_pcoa.png"
       , gg_weighted_unifrac_pcoa_whole
       , height=5, width=7)



# right cohort
wu_dm_right <- distance(rightphyloseq_rare, method="wunifrac")

pcoa_wu_right <- ordinate(rightphyloseq_rare, method="PCoA", distance=wu_dm_right)

plot_ordination(rightphyloseq_rare, pcoa_wu_right, color = "response_patient_by_visit")

gg_weighted_unifrac_pcoa_right <- plot_ordination(rightphyloseq_rare, pcoa_wu_right, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit") +
  stat_ellipse()

gg_weighted_unifrac_pcoa_right

ggsave("plot_beta_weighted_unifrac_response_right_pcoa.png"
       , gg_weighted_unifrac_pcoa_right
       , height=5, width=7)



# left
wu_dm_left <- distance(leftphyloseq_rare, method="wunifrac")

pcoa_wu_left <- ordinate(leftphyloseq_rare, method="PCoA", distance=wu_dm_left)

plot_ordination(leftphyloseq_rare, pcoa_wu_left, color = "response_patient_by_visit")

gg_weighted_unifrac_pcoa_left <- plot_ordination(leftphyloseq_rare, pcoa_wu_left, color = "response_patient_by_visit") +
  labs(col = "Response Status by Visit") +
  stat_ellipse()

gg_weighted_unifrac_pcoa_left

ggsave("plot_beta_weighted_unifrac_response_left_pcoa.png"
       , gg_weighted_unifrac_pcoa_left
       , height=5, width=7)


