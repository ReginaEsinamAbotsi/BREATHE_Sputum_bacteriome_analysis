met<-read.csv("meta_edited13122021.csv", row.names = 1, stringsAsFactors = T)

phy <- readRDS("phyloseq_object_trial_only28thNov2021_CLR.RDS")
sample_data(phy) <- met
colnames(sample_data(phy))
str(sample_data(phy2)$visit)
phy2 <- readRDS("phyloseq_object_trial_only28thNov2021.RDS")
sample_data(phy2) <- met
sample_data(phy) <- met


saveRDS(phy, "breathe_sputum_phyloseq_trial_13Dec2021CLR.RDS") 
saveRDS(phy2, "breathe_sputum_phyloseq_trial_13Dec2021.RDS") 
