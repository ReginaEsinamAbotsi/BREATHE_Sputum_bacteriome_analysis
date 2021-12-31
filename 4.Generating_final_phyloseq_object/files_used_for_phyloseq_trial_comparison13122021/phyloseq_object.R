
#13th December 2021
#Creating phyloseq object used for all downstream
library(phyloseq)


#Creating phyloseq object with the phylogenetic tree. 
original_abundances <- read.csv("ASVcount948V2.csv", header = TRUE, row.names = 1)
original_taxonomy <- read.csv("taxdat948V2.csv", row.names = 1)
original_metadata <- read.csv("metadata_longer_namesfinaldiv.csv", header = TRUE, row.names = 1)
mytree<- readRDS("fitGTR_regina_generated.RDS")


ps_original <- phyloseq(otu_table(original_abundances, taxa_are_rows = T),
                        tax_table(as.matrix(original_taxonomy)), # use as.matrix to convert to matrix, otherwise ranks or taxon IDs/ASVs are not well preserved
                        sample_data(original_metadata))

ps_original

saveRDS(ps_original, "breathe_sputum_complete_phyloseq13122021.RDS") 




