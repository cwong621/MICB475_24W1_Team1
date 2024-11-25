## left V3 resp vs healthy
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

load("./HIV_final.RData")
right <- read_delim("left cohort.csv", delim=",")

HIV_final@sam_data$response_patient_by_visit <- factor(HIV_final@sam_data$response_patient_by_visit, 
                                                       levels = c("nonresponsive 2", "responsive 2", "Healthy 2", "ART-experienced nonresponsive 2", "nonresponsive 3", "responsive 3", "Healthy 3", "ART-experienced nonresponsive 3" ))


ps_left <- subset_samples(HIV_final, `original.sample.id` %in% right$original.sample.id)

v2 <- subset_samples(ps_left, Visit==3 & Cohort_Short!="B" & response_patient!="nonresponsive")
set.seed(1)
out = ancombc(data = v2, tax_level = "Genus", 
              formula = "response_patient_by_visit", 
              lib_cut = 1000, 
              group = "response_patient_by_visit", struc_zero = TRUE, neg_lb = TRUE,
              conserve = TRUE, alpha = 0.05, global = TRUE)
res = out$res
res_global = out$res_global

tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "healthy - responsive")

colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log 2 Fold Changes from the Primary Result") %>%
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

df_fig_v2 <- df_lfc %>% 
  filter(abs(`healthy - responsive`) > 1.5) %>%
  transmute(Taxon, 
            `healthy - responsive` = round(`healthy - responsive`, 2)) %>%
  pivot_longer(cols = `healthy - responsive`,
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
                       name = "Log 2 Fold Change") +
  geom_text(aes(group, Taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Left cohort, Visit 2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_v2

ggsave("ANCOM/leftv2healthyvsresp.jpg", p_v2)
