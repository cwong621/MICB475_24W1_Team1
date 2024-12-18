PICRUST_function = function(metadata, abundance_data) {

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = metadata$'...1'
sample_names = append(sample_names, "#OTU ID")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] 

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata = metadata[metadata$`...1` %in% abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq ####
#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                        metadata = metadata, group = "response_patient_by_visit", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(`#OTU ID` %in% feature_with_p_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
abundance_desc = abundance_desc[,-c(nrow(abundance_desc):ncol(abundance_desc))] 

# Lead the function in
source("DESeq2_function.R")

# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, metadata, "response_patient_by_visit")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05 & abs(log2FoldChange)>2)
# You can also filter by Log2fold change

sig_res <- sig_res[order(sig_res$log2FoldChange),]

rightgg <- ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity") + 
  theme_bw()
}

