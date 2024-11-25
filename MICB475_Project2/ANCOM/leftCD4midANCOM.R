#200-499 left 

library(DT)
library(ANCOMBC)
library(nlme)
library(tidyverse)
library(compositions)
library(readr)
library(ALDEx2)
library(GUniFrac)
library(phyloseq)
library(phyloseqCompanion)



####ANCOM####
#https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
# Load data 
load("./HIV_final.RData")
HIV_final@sam_data$response_patient_by_visit <- factor(HIV_final@sam_data$response_patient_by_visit, levels = c("Healthy 2", "Healthy 3", "responsive 2", "nonresponsive 2", "responsive 3", "nonresponsive 3", "ART-experienced nonresponsive 2", "ART-experienced nonresponsive 3" ))
left <- read_delim("left cohort.csv", delim=",")

HIV_final@sam_data$response_patient_by_visit <- factor(HIV_final@sam_data$response_patient_by_visit, 
                                                       levels = c("nonresponsive 2", "responsive 2", "Healthy 2", "ART-experienced nonresponsive 2", "nonresponsive 3", "responsive 3", "Healthy 3", "ART-experienced nonresponsive 3" ))


ps_left <- subset_samples(HIV_final, `original.sample.id` %in% left$original.sample.id)


aids <- subset_samples(ps_left, Visit==2 & CD4._Cat=="200-499" & Cohort_Short=="A")
set.seed(1)
out = ancombc(data = aids, tax_level = "Genus", 
              formula = "response_patient_by_visit", 
              lib_cut = 1000, 
              group = "response_patient_by_visit", struc_zero = TRUE, neg_lb = TRUE,
              conserve = TRUE, alpha = 0.05, global = TRUE)
res = out$res
res_global = out$res_global

tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "responsive - nonresponsive")

colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_lfc) = col_name
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())

colnames(df_se) = col_name
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

df_fig_v2 = df_lfc %>% 
  filter(abs(`responsive - nonresponsive`) > 1.5) %>%
  transmute(Taxon, 
            `responsive - nonresponsive` = round(`responsive - nonresponsive`, 2)) %>%
  pivot_longer(cols = `responsive - nonresponsive`, 
               names_to = "group", values_to = "value") %>%
  arrange(Taxon)

lo = floor(min(df_fig_v2$value))
up = ceiling(max(df_fig_v2$value))
mid = 0
p_v2 = df_fig_v2 %>%
  ggplot(aes(x = group, y = Taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "Log Fold Change") +
  geom_text(aes(group, Taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Left, CD4 count 200-499, Before Treatment") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_v2

ggsave("ANCOM/leftCD4midnaive.jpg", p_v2)
