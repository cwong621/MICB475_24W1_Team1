## Generate phyloseq object

# Load packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

# Load files

meta <- read_delim(file="./hiv_export/hiv_metadata_filt.tsv", delim = "\t")

otu <- read_delim(file="./hiv_export/table_export/feature-table.txt", delim = "\t", skip = 1)

tax <- read_delim(file="./hiv_export/taxonomy_export/taxonomy.tsv", delim = "\t")

phylotree <- read_delim("./hiv_export/rooted_tree_export/tree.nwk")

# Format OTU table

# Save everything except #OTU ID as matrix
otu_mat <- as.matrix(otu[,-1])
# Make #OTU ID the row names
rownames(otu_mat) <- otu$`#OTU ID`
# Save as table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

# Format metadata

# Save everything except sample_id as dataframe
meta_df <- as.data.frame(meta[,-1])
# Make sample_ids into row names
rownames(meta_df) <- meta$`sample-id`
# Make phyloseq sample data
META <- sample_data(meta_df)

# Format taxonomy table

# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(col=Taxon, sep="; ",
           into = c("Domain", "Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
# Save everything except feature IDs
tax_mat <- tax_mat[,-1]
# Make feature IDs the row nanmes
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)

# Merge everything to create phyloseq object HIV

HIV <- phyloseq(OTU, META, TAX, phylotree)
