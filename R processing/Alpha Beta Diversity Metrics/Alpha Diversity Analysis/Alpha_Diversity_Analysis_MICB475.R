### Alpha Diversity Analysis ###
# MICB 475 Final Project - Team 1 


# Load libraries

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggplot2)
library(rbiom)
library(ggsignif)
library(gridExtra)


# Load in data

load("HIV_rare.RData")
load("leftphyloseq_rare.RData")
load("rightphyloseq_rare.RData")


#### Preliminary Analysis ####


# Obtain alpha diversity measures

plot_A <- plot_richness(rightphyloseq_rare)
plot_richness(rightphyloseq_rare, measures = c("Shannon","Chao1")) 
     # preliminary visualization to see what we are working with

plot_B <- plot_richness(leftphyloseq_rare)
      
      # From here on, the "left" and "right" cohorts will be labeled as cohorts
      # "A" and "B" instead.
  
# Comparing Shannon's (abundance and evenness) and Chao1 (total richness)
# between each of cohort A and B
alpha_response_richness_B <- plot_richness(rightphyloseq_rare, 
  x = "response_patient_by_visit", measures = c("Shannon","Chao1")) +
  xlab("Response Status by Visit") +
  geom_boxplot()
alpha_response_richness_B


alpha_response_richness_A <- plot_richness(leftphyloseq_rare, 
  x = "response_patient_by_visit", measures = c("Shannon","Chao1")) +
  xlab("Response Status by Visit") +
  geom_boxplot()
alpha_response_richness_A

  # Seems to be some significant measures in cohort A.


# Save these images for discussion, if we want to continue analysing these 
# metrics.
ggsave(filename = "chao1_shannon_alpha_response_richness_B.png"
       , alpha_response_richness_B
       , height=5, width=6)

ggsave(filename = "chao1_shannon_alpha_response_richness_A.png"
       , alpha_response_richness_A
       , height=5, width=6)



## Preparation for further analysis

# Obtaining richness estimates from the phyloseq objects
alpha_B <- estimate_richness(rightphyloseq_rare)
alpha_A <- estimate_richness(leftphyloseq_rare)

# Getting the phyloseq objects into dataframe formats for ease of manipulation
A_phylo_dat <- sample_data(leftphyloseq_rare)
A_phylo_datfrm <- data.frame(A_phylo_dat)

B_phylo_dat <- sample_data(rightphyloseq_rare)
B_phylo_datfrm <- data.frame(B_phylo_dat)

# Adding richness measures onto the new dataframes for Shannon/Chao1 testing
A_phylo_datfrm_wdiv <- data.frame(A_phylo_datfrm, alpha_A)
B_phylo_datfrm_wdiv <- data.frame(B_phylo_datfrm, alpha_B)




#### Shannon/Chao1 Analysis ####


# Chao1 significance testing

kruskal_chao1_A <- kruskal.test(Chao1 ~ `response_patient_by_visit`, 
                                data = A_phylo_datfrm_wdiv)
kruskal_chao1_A  # something significant here!


kruskal_chao1_B <- kruskal.test(Chao1 ~ `response_patient_by_visit`, 
                                data = B_phylo_datfrm_wdiv)
kruskal_chao1_B  # nothing significant here...


# Shannon significance testing

kruskal_shannon_A <- kruskal.test(Shannon ~ `response_patient_by_visit`, 
                                  data = A_phylo_datfrm_wdiv)
kruskal_shannon_A  # something significant here!


kruskal_shannon_B <- kruskal.test(Shannon ~ `response_patient_by_visit`, 
                                  data = B_phylo_datfrm_wdiv)
kruskal_shannon_B  # nothing significant here


## Log-transforming and performing ANOVAs on significant cohort measures
## to see with TukeyHSD() which comparisons are the significant ones


# Chao1
A_chao1_vs_rpbv_logtrsfm <- lm(log(Chao1) ~ `response_patient_by_visit`, 
                               data=A_phylo_datfrm_wdiv)
anova_A_chao1_vs_rpbv_logtrsfm <- aov(A_chao1_vs_rpbv_logtrsfm)
summary(anova_A_chao1_vs_rpbv_logtrsfm)

# Seeing which comparison is significant
TukeyHSD(anova_A_chao1_vs_rpbv_logtrsfm) 
      # nonresponsive 3 vs both Healthy!


# Shannon
A_shannon_vs_rpbv_logtrsfm <- lm(log(Shannon) ~ `response_patient_by_visit`, 
                                 data=A_phylo_datfrm_wdiv)
anova_A_shannon_vs_rpbv_logtrsfm <- aov(A_shannon_vs_rpbv_logtrsfm)
summary(anova_A_shannon_vs_rpbv_logtrsfm)

# Seeing which comparison is significant
TukeyHSD(anova_A_shannon_vs_rpbv_logtrsfm) 
      # again, nonresponsive 3 vs both Healthy!






#### Publication Plot Generation - Shannon (and Chao1) ####


## Data polishing - Cohort A

    # Filtering the data to not include the ART-experienced nonresponsive cohort,
    # as their sample number is low and show no significant variations

A_phylo_datfrm_wdiv_noARTexpV2 <- filter(A_phylo_datfrm_wdiv, 
                                      response_patient_by_visit != "ART-experienced nonresponsive V2 2")
A_phylo_datfrm_wdiv_noARTexp <- filter(A_phylo_datfrm_wdiv_noARTexpV2, 
                                      response_patient_by_visit != "ART-experienced nonresponsive V3 3")

    # Renaming response group names
A_phylo_datfrm_wdiv_noARTexp_renamed <- mutate(A_phylo_datfrm_wdiv_noARTexp, 
                                                 response_patient_by_visit = recode(response_patient_by_visit,
                                                "Healthy 2" = "Healthy (Week 0)",
                                                "Healthy 3" = "Healthy (Week 24)",
                                                "nonresponsive 2" = "Nonresponsive (Week 0)",
                                                "nonresponsive 3" = "Nonresponsive (Week 24)",
                                                "responsive 2" = "Responsive (Week 0)",
                                                "responsive 3" = "Responsive (Week 24)") )
A_phylo_datfrm_wdiv_noARTexp_renamed


# Plotting Chao1 measures
chao1_reponse_A_plot_sig <- ggplot(A_phylo_datfrm_wdiv_noARTexp_renamed, aes(response_patient_by_visit, Chao1, ylim =40)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Nonresponsive (Week 24)","Healthy (Week 0)"), c("Nonresponsive (Week 24)", "Healthy (Week 24)")),
              y_position = c(500,550),
              annotations = c("*", "*")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Cohort A Response Status") +
  ylim(100, 600) + # tweak axis to include significance symbols
  ylab("Alpha Diversity Measure")

chao1_reponse_A_plot_sig

# Plotting Shannon Measures
shannon_reponse_A_plot_sig <- ggplot(A_phylo_datfrm_wdiv_noARTexp_renamed, aes(response_patient_by_visit, Shannon, ylim =40)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Nonresponsive (Week 24)","Healthy (Week 24)"), c("Nonresponsive (Week 24)", "Healthy (Week 0)")),
              y_position = c(5.5,6),
              annotations = c("*", "*")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab(NULL) + # Removed X-axis labels since this panel will go in the manuscript and have its axis described in the figure caption
  ylim(3,6.2) + # tweak axis to include significance symbols
  ylab("Shannon Index")
shannon_reponse_A_plot_sig


### Data Polishing again - Cohort B

# Filtering the data to not include the ART-experienced nonresponsive cohort
B_phylo_datfrm_wdiv_noARTexpV2 <- filter(B_phylo_datfrm_wdiv, 
                                         response_patient_by_visit != "ART-experienced nonresponsive V2 2")
B_phylo_datfrm_wdiv_noARTexp <- filter(B_phylo_datfrm_wdiv_noARTexpV2, 
                                       response_patient_by_visit != "ART-experienced nonresponsive V3 3")

# Again, renaming the response groups for publication
B_phylo_datfrm_wdiv_noARTexp_renamed <- mutate(B_phylo_datfrm_wdiv_noARTexp, 
                                                 response_patient_by_visit = recode(response_patient_by_visit,
                                                 "Healthy 2" = "Healthy (Week 0)",
                                                 "Healthy 3" = "Healthy (Week 24)",
                                                 "nonresponsive 2" = "Nonresponsive (Week 0)",
                                                 "nonresponsive 3" = "Nonresponsive (Week 24)",
                                                 "responsive 2" = "Responsive (Week 0)",
                                                 "responsive 3" = "Responsive (Week 24)") )


# Plotting Chao1 measures
chao1_reponse_B_plot_sig <- ggplot(B_phylo_datfrm_wdiv_noARTexp_renamed, aes(response_patient_by_visit, Chao1, ylim =40)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Cohort B Response Status") +
  ylab("Alpha Diversity Measure")

chao1_reponse_B_plot_sig

# Plotting Shannon measures
shannon_reponse_B_plot_sig <- ggplot(B_phylo_datfrm_wdiv_noARTexp_renamed, aes(response_patient_by_visit, Shannon, ylim =40)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab(NULL) +  # Removed X axis label, since this plot is going in the manuscript and will have its axis described in the figure caption
  ylab("Shannon Index")

shannon_reponse_B_plot_sig



## Save all the Chao1/Shannon graphs as images

# Cohort A
ggsave(filename = "chao1_reponse_A_plot_sig.png"
       , chao1_reponse_A_plot_sig
       , height=5, width=6)

ggsave(filename = "shannon_reponse_A_plot_sig.png"
       , shannon_reponse_A_plot_sig
       , height=5, width=6)

# Cohort B
ggsave(filename = "chao1_reponse_B_plot_sig.png"
       , chao1_reponse_B_plot_sig
       , height=5, width=6)

ggsave(filename = "shannon_reponse_B_plot_sig.png"
       , shannon_reponse_B_plot_sig
       , height=5, width=6)






##### Phylogenetic diversity - Faith's #####

# calculate Faith's phylogenetic diversity as PD 
phylo_dist_B <- pd(t(otu_table(rightphyloseq_rare)), phy_tree(rightphyloseq_rare),
                 include.root=F) 

phylo_dist_A <- pd(t(otu_table(leftphyloseq_rare)), phy_tree(leftphyloseq_rare),
                  include.root=F) 

# add PD to metadata table
sample_data(rightphyloseq_rare)$PD <- phylo_dist_B$PD

sample_data(leftphyloseq_rare)$PD <- phylo_dist_A$PD


#### Faith's Preliminary analysis ####

# plot any metadata category against the PD
faith_reponse_B_plot_pd <- ggplot(sample_data(rightphyloseq_rare), aes(response_patient_by_visit, PD)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Right Cohort Response Status") +
  ylab("Phylogenetic Diversity")

faith_reponse_A_plot_pd <- ggplot(sample_data(leftphyloseq_rare), aes(response_patient_by_visit, PD)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Left Cohort Response Status") +
  ylab("Phylogenetic Diversity")

# view plots
faith_reponse_B_plot_pd
faith_reponse_A_plot_pd


# Add phylogenetic distance to the dataframes with the other alpha diversity measures

A_phylo_datfrm_wdiv$PD <- phylo_dist_A$PD
B_phylo_datfrm_wdiv$PD <- phylo_dist_B$PD

### Run Kruskal-Wallis statistical testing ###


# Testing for significance in at least one comparison for each cohort
kruskal_phylodiv_A <- kruskal.test( PD ~ response_patient_by_visit, data = A_phylo_datfrm_wdiv)
kruskal_phylodiv_A   ## there is something significant here!

kruskal_phylodiv_right <- kruskal.test( PD ~ response_patient_by_visit, data = B_phylo_datfrm_wdiv)
kruskal_phylodiv_right  ## nothing significant here

### Log-transforming and performing an ANOVA on the left cohort - as was described in course material

A_pd_vs_rpbv_logtrsfm <- lm(log(PD) ~ `response_patient_by_visit`, data=A_phylo_datfrm_wdiv)
anova_A_pd_vs_rpbv_logtrsfm <- aov(A_pd_vs_rpbv_logtrsfm)
summary(anova_A_pd_vs_rpbv_logtrsfm)

# Seeing which comparison is significant
TukeyHSD(anova_A_pd_vs_rpbv_logtrsfm)
      # nonresponsive after ART and healthy at week 0 only

### Generating polished plots with statistical test annotations

# Filtering the data to not include the ART-experienced nonresponsive cohort
A_phylo_datfrm_wdiv_noARTexpV2 <- filter(A_phylo_datfrm_wdiv, response_patient_by_visit != "ART-experienced nonresponsive V2 2")
A_phylo_datfrm_wdiv_noARTexp <- filter(A_phylo_datfrm_wdiv_noARTexpV2, response_patient_by_visit != "ART-experienced nonresponsive V3 3")

# Again, filter this so that the response names are better for publication

A_phylo_datfrm_wdiv_noARTexp_renamed <- mutate(A_phylo_datfrm_wdiv_noARTexp, 
                                               response_patient_by_visit = recode(response_patient_by_visit,
                                                                                  "Healthy 2" = "Healthy (Week 0)",
                                                                                  "Healthy 3" = "Healthy (Week 24)",
                                                                                  "nonresponsive 2" = "Nonresponsive (Week 0)",
                                                                                  "nonresponsive 3" = "Nonresponsive (Week 24)",
                                                                                  "responsive 2" = "Responsive (Week 0)",
                                                                                  "responsive 3" = "Responsive (Week 24)") )


faith_reponse_A_plot_pd <- ggplot(A_phylo_datfrm_wdiv_noARTexp_renamed, aes(response_patient_by_visit, PD, ylim =40)) + 
  geom_boxplot() +                 # used the dataframe with diversity metrics to avoid renaming again 
  geom_signif(comparisons = list(c("Nonresponsive (Week 24)","Healthy (Week 0)")),
              y_position = c(30),
              annotations = c("*")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab(NULL) +
  ylim(c(10,32)) + #tweaked the axis a bit since the asterisk was getting cut off
  ylab("Phylogenetic Diversity")

faith_reponse_A_plot_pd

# Filtering the right plot to not include the ART-experienced group as well
B_phylo_datfrm_wdiv_noARTexpV2 <- filter(B_phylo_datfrm_wdiv, response_patient_by_visit != "ART-experienced nonresponsive V2 2")
B_phylo_datfrm_wdiv_noARTexp <- filter(B_phylo_datfrm_wdiv_noARTexpV2, response_patient_by_visit != "ART-experienced nonresponsive V3 3")

# once again, rename response categories

B_phylo_datfrm_wdiv_noARTexp_renamed <- mutate(B_phylo_datfrm_wdiv_noARTexp, 
                                               response_patient_by_visit = recode(response_patient_by_visit,
                                                                                  "Healthy 2" = "Healthy (Week 0)",
                                                                                  "Healthy 3" = "Healthy (Week 24)",
                                                                                  "nonresponsive 2" = "Nonresponsive (Week 0)",
                                                                                  "nonresponsive 3" = "Nonresponsive (Week 24)",
                                                                                  "responsive 2" = "Responsive (Week 0)",
                                                                                  "responsive 3" = "Responsive (Week 24)") )


# Regenerating right plot without that ART-experienced group
faith_reponse_B_plot_pd <- ggplot(B_phylo_datfrm_wdiv_noARTexp_renamed, aes(response_patient_by_visit, PD)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab(NULL) +
  ylab("Phylogenetic Diversity")

faith_reponse_B_plot_pd

############

# save graphs as images
ggsave(filename = "plot_faith_reponse_richness_B.png"
       , faith_reponse_B_plot_pd
       , height=5, width=6)

ggsave(filename = "plot_faith_reponse_richness_A.png"
       , faith_reponse_A_plot_pd
       , height=5, width=6
       , )

### Generating Paneled Figure for Publication

alpha_diversity_paneled <- grid.arrange(
  arrangeGrob(
    shannon_reponse_A_plot_sig + theme(legend.position = "none", plot.title = element_text(hjust = 0)) + labs(title = bquote(bold("A"))), 
    shannon_reponse_B_plot_sig + theme(legend.position = "none", plot.title = element_text(hjust = 0)) + labs(title = bquote(bold("B"))), 
    faith_reponse_A_plot_pd + theme(legend.position = "none", plot.title = element_text(hjust = 0)) + labs(title = bquote(bold("C"))), 
    faith_reponse_B_plot_pd + theme(legend.position = "none", plot.title = element_text(hjust = 0)) + labs(title = bquote(bold("D"))),
    nrow = 2, ncol = 2), nrow = 2,
  heights = c(20, 1)  # Adjust heights to give more space to plots
)
alpha_diversity_paneled

ggsave(filename = "alpha_diversity_paneled.png"
       , alpha_diversity_paneled
       , height=8, width=10
       , )
