ANCOMBC_HealthyRef = function(phylo) {
  set.seed(1)
  
  out = ancombc(data = phylo, tax_level = "Genus", 
                formula = "response_patient_by_visit", 
                lib_cut = 1000, 
                group = "response_patient_by_visit", struc_zero = TRUE, neg_lb = TRUE,
                conserve = TRUE, alpha = 0.05)
  res = out$res
  
  col_name = c("Taxon", "Intercept", "nonresponsive - healthy", "responsive - healthy")
  
  df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
    mutate(taxon_id = res$diff_abn$taxon) %>%
    dplyr::select(taxon_id, everything())
  colnames(df_lfc) = col_name
  
  df_fig_v2 = df_lfc %>% 
    filter(abs(`nonresponsive - healthy`) > 1.5 | abs(`responsive - healthy`) > 1.5) %>%
    transmute(Taxon, 
              `nonresponsive - healthy` = round(`nonresponsive - healthy`, 2),
              `responsive - healthy` = round(`responsive - healthy`, 2)) %>%
    pivot_longer(cols = `nonresponsive - healthy`:`responsive - healthy`,
                 names_to = "group", values_to = "value") %>%
    arrange(Taxon)
  
  lo = floor(min(df_fig_v2$value))
  up = ceiling(max(df_fig_v2$value))
  mid = 0
  heatmap = df_fig_v2 %>%
    ggplot(aes(x = group, y = Taxon, fill = value)) + 
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         na.value = "white", midpoint = 0, limit = c(-3, 3),
                         name = "Log 2 Fold Change") +
    geom_text(aes(group, Taxon, label = value), color = "black", size = 4) +
    theme_minimal()
  
}

ANCOMBC_NRRef = function(phylo) {
  set.seed(1)
  
  out = ancombc(data = phylo, tax_level = "Genus", 
                formula = "response_patient_by_visit", 
                lib_cut = 1000, 
                group = "response_patient_by_visit", struc_zero = TRUE, neg_lb = TRUE,
                conserve = TRUE, alpha = 0.05)
  res = out$res
  
  col_name = c("Taxon", "Intercept", "responsive - nonresponsive")
  
  df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
    mutate(taxon_id = res$diff_abn$taxon) %>%
    dplyr::select(taxon_id, everything())
  colnames(df_lfc) = col_name
  
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
  heatmap = df_fig_v2 %>%
    ggplot(aes(x = group, y = Taxon, fill = value)) + 
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         na.value = "white", midpoint = 0, limit = c(-3, 3),
                         name = "Log 2 Fold Change") +
    geom_text(aes(group, Taxon, label = value), color = "black", size = 4) +
    theme_minimal()
  
}

ANCOMBC_V3V2 = function(phylo) {
  set.seed(1)
  
  out = ancombc(data = phylo, tax_level = "Genus", 
                formula = "response_patient_by_visit", 
                lib_cut = 1000, 
                group = "response_patient_by_visit", struc_zero = TRUE, neg_lb = TRUE,
                conserve = TRUE, alpha = 0.05)
  res = out$res
  
  col_name = c("Taxon", "Intercept", "Week 24 - Week 0")
  
  df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
    mutate(taxon_id = res$diff_abn$taxon) %>%
    dplyr::select(taxon_id, everything())
  colnames(df_lfc) = col_name
  
  df_fig_v2 = df_lfc %>% 
    filter(abs(`Week 24 - Week 0`) > 1.5) %>%
    transmute(Taxon, 
              `Week 24 - Week 0` = round(`Week 24 - Week 0`, 2)) %>%
    pivot_longer(cols = `Week 24 - Week 0`,
                 names_to = "group", values_to = "value") %>%
    arrange(Taxon)
  
  heatmap = df_fig_v2 %>%
    ggplot(aes(x = group, y = Taxon, fill = value)) + 
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         na.value = "white", midpoint = 0, limit = c(-3, 3),
                         name = "Log 2 Fold Change") +
    geom_text(aes(group, Taxon, label = value), color = "black", size = 4) +
    theme_minimal()
  
}