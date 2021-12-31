#http://www2.decipher.codes/AlignSequences.html
# load the DECIPHER library in R
library(DECIPHER)

# specify the path to the FASTA file (in quotes)
fas <- "/Users/reginaabotsi/Documents/Files/REGINA/MY_PROJECTS/FOLLOW_UP_MICROBIOME/16S_analysis_my_attempt11thMay2021/dna_fasta.fasta"

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet(fas)

# look at some of the sequences (optional)
seqs

library(ape)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)


# view the alignment in a browser (optional)
BrowseSeqs(alignment, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(alignment,
                file="/Users/reginaabotsi/Documents/Files/REGINA/MY_PROJECTS/FOLLOW_UP_MICROBIOME/16S_analysis_my_attempt11thMay2021/alignment_regina2.fasta")

library(phangorn)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

tree<-phy_tree(fitGTR$tree)
saveRDS(fitGTR, "fitGTR_regina_generatedCORRECT.RDS")

saveRDS(phy, "breathe_sputum_complete_phyloseq.RDS") 
readRDS("breathe_sputum.RDS")
