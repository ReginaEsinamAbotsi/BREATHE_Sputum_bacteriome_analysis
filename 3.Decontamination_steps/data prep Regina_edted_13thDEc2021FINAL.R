### Produced by Sugnet Lubbe, Department of Statistics and ActSci, 
### Stellenbosch University, May 2020
### --------------------------------------------------
### Revised by Shantelle Claassen-Weitz, Department of Pathology, 
### University of Cape Town, October 2020
### --------------------------------------------------
### Revised by Regina Abotsi, Department of Molecular and Cell Biology
### University of Cape Town, 26th January 2021


pathname <- NULL
pathname <- "/Users/reginaabotsi/Documents/Files/REGINA/MY_PROJECTS/MY_R_projects/Microbiome/Decontamination_Regina_edit25012020/UpdatedfromShantelle/Decontamination_13thDEc2021/"
#=============================================================================================================
# --- Read data 
#=============================================================================================================


datmat <- read.csv(paste(pathname,"data/ASVs_counts.csv",sep=""), header=T, row.names=1)
head(datmat[,1:3])
datmat <- as.matrix(datmat)
dim(datmat) #3665 ASVs  1151 SampleIDs

taxdat <- read.csv(paste(pathname,"data/ASVs_taxonomy.csv",sep=""), header=T, row.names=1)
head(taxdat)
dim(taxdat) #3665     7

ASV.order <- order(taxdat[,1], taxdat[,2], taxdat[,3], taxdat[,4], taxdat[,5], taxdat[,6], taxdat[,7])
taxdat <- taxdat[ASV.order,]
datmat <- datmat[match(rownames(taxdat),rownames(datmat)),]  

#select ASVs=Bacteria at Kingdom level only
datmat0 <- datmat[!is.na(taxdat$Kingdom) & taxdat$Kingdom=="Bacteria",]
taxdat0 <- taxdat[!is.na(taxdat$Kingdom) & taxdat$Kingdom=="Bacteria",]
dim(datmat0) #3481 ASVs  1151 SampleIDs
dim(taxdat0) #3481  7

#=============================================================================================================
# --- Identify biological samples, spikes, mocks and negative controls
#=============================================================================================================

labdata <- read.csv (paste(pathname,"data/Runs01to03_labdataUpdated.csv",sep=""))
dim(labdata)
head(labdata)

labdata$sample_id_run <- gsub("-",".",labdata$sample_id_run)

# === Biological samples
# ----------------------
samples <- labdata$sample_id_run[(labdata$sample_type=="SP" | labdata$sample_type!="sequencing_control")]
length(samples) #960 SP specimens
table(is.na(labdata$templateDNA_16SqPCR_copies), labdata$sample_type, exclude=NULL)
samples <- samples[!is.na(labdata$templateDNA_16SqPCR_copies[match(samples, labdata$sample_id_run)])]
length(samples) #960 SP specimens with 16S copy number data

# === Negative controls
# ---------------------
table(labdata$nontemplate, exclude=NULL)
negcontrols <- labdata$sample_id_run[!is.na(labdata$nontemplate) & labdata$nontemplate=="Primestore"]
length(negcontrols) #43 Primestore specimens
table(is.na(labdata$templateDNA_16SqPCR_copies), labdata$nontemplate, exclude=NULL)
negcontrols <- negcontrols[!is.na(labdata$templateDNA_16SqPCR_copies[match(negcontrols, labdata$sample_id_run)])]
length(negcontrols) #43 Primestore specimens with 16S copy number data

# === Mocks
# ---------------------
table(labdata$mock, exclude=NULL)
mocks <- labdata$sample_id_run[!is.na(labdata$mock) & (labdata$mock=="Zymobiomics_cells" | labdata$mock=="Zymobiomics_DNA")]
length(mocks) #24 Zymo controls specimens
table(is.na(labdata$templateDNA_16SqPCR_copies), labdata$mock, exclude=NULL)
mocks <- mocks[!is.na(labdata$templateDNA_16SqPCR_copies[match(mocks, labdata$sample_id_run)])]
length(mocks) #24 mocks with 16S copy number data

# === check for missing samples
# --- Biological samples
samples[is.na(match(samples,colnames(datmat0)))]
#----- no missing samples

# --- Negative controls
negcontrols[is.na(match(negcontrols,colnames(datmat0)))]
# --- no missing samples

# --- Mocks
mocks[is.na(match(mocks,colnames(datmat0)))]
# --- no missing samples

# === List of SampleIDs prior to QC steps:
#write.csv(samples, paste(pathname, "SampleIDs_Samples_All.csv", sep=""))
#write.csv(negcontrols, paste(pathname, "SampleIDs_Negatives.csv", sep=""))
#write.csv(mocks, paste(pathname, "SampleIDs_Mocks.csv", sep=""))
#write.csv(datmat0, paste(pathname, "ASV_bacteria.csv", sep=""))

#=============================================================================================================
# --- Mocks
#=============================================================================================================

col.legend <- function (col.vec, names, per.col=50, per.page=12)
{
  dev.new()
  old.mai <- par("mai")
  on.exit(par(mai=old.mai))
  par(mai=rep(0,4))
  n <- length(col.vec)
  num.cols <- ceiling(n/per.col)
  j <- 0
  while(num.cols > 0)
  {  cc <- per.page
  plot (c(0.95,cc+0.95),c(0,per.col),type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  i <- 1+j
  x <- 1
  y <- per.col
  while (i <= n)
  {
    points (x,y,pch=15,col=col.vec[i])
    text (x,y,names[i],pos=4,cex=0.55)
    i <- i + 1
    y <- y - 1
    if (y ==0) { y <- per.col
    x <- x + 1
    }
  }
  num.cols <- num.cols-per.page
  j <- j+per.col*per.page
  if (num.cols>0) { dev.new()
    par(mai=rep(0,4))
  }
  }
}

seqmat <- datmat0
seq.mocks <- seqmat[,match(mocks, colnames(seqmat))]
seq.mocks <- seq.mocks[apply(seq.mocks,1,sum)>0,]
dim(seq.mocks) #13 ASVs 24 mocks
mock.type <- labdata$mock[match(mocks,labdata$sample_id_run)]
mock.type  
seq.mocks <- seq.mocks[,order(mock.type)]
ASV.totals <- t(apply(seq.mocks, 1, function(x) tapply(x, mock.type, sum)))
ASV.order <- rev(order(ASV.totals[,1]+ASV.totals[,2]))
output <- data.frame(seq.mocks[ASV.order,],ASV.totals[ASV.order,])
output <- data.frame (output, taxdat0[match(rownames(output),rownames(taxdat0)),])
#write.csv (output, paste(pathname,"Output_mocks.csv",sep=""))

ZymoCells.prop <- c(0.157,0.133,0.1,0.0,0.113,0.046,0.159,0.188,0.104)
names(ZymoCells.prop) <- c("ASV_56","ASV_46","ASV_51","ASV_153","ASV_54","ASV_61","ASV_92","ASV_82","ASV_85")

ZymoDNA.prop <- c(0.157,0.133,0.1,0.0,0.113,0.046,0.159,0.188,0.104)
names(ZymoDNA.prop) <- c("ASV_56","ASV_46","ASV_51","ASV_153","ASV_54","ASV_61","ASV_92","ASV_82","ASV_85")

seq.mocks.prop <- apply(seq.mocks, 2, function(x)x/sum(x))
#seq.mocks.prop <- rbind (seq.mocks.prop, acnes=0, faecalis=0, odontolyticus=0, radiodurans=0)
seq.mocks.prop <- as.data.frame(cbind (ZymoCells=0, dummy1=NA, seq.mocks.prop, dummy2=NA, ZymoDNA=0))
seq.mocks.prop$ZymoCells[match(names(ZymoCells.prop),rownames(seq.mocks.prop))] <- ZymoCells.prop#?
seq.mocks.prop$ZymoDNA[match(names(ZymoDNA.prop),rownames(seq.mocks.prop))] <- ZymoDNA.prop#?

mock.cols <- rep(rgb(red=166,green=166,blue=166,maxColorValue=255),nrow(seq.mocks))
mock.cols[match(c("ASV_51","ASV_153","ASV_54","ASV_92",
                  "ASV_56","ASV_61","ASV_85",
                  "ASV_82","ASV_46",
                  "ASV_1","ASV_1827","ASV_2",
                  "ASV_18"),rownames(seq.mocks.prop))] <- 
  rgb(rbind(c(0,166,205),c(0,199,204),c(0,199,204),c(150,0,0),
            c(192,50,0),c(2,76,240),c(200,18,35),
            c(214,0,0),c(240,115,0),
            c(166,166,166),c(166,166,166),c(166,166,166),
            c(166,166,166)),maxColorValue=255)

names.vec <- colnames(seq.mocks.prop)
names.vec[match(c("dummy1","dummy2"),names.vec)] <- ""
barplot (as.matrix(seq.mocks.prop), col=mock.cols, names.arg=names.vec, las=2, cex.names=0.5)
col.legend (rev(mock.cols), rev(rownames(seq.mocks.prop)))

###Comment (SC): 
# DNA ZBmocks are more comparable to theoretical compositions than extracted ZBmocks
# Compare these profiles to ZBmock profiles as per manufacturers reports: 
# Zymo.prop  <- ASV_56        :0.157
#               ASV_46        :0.133
#               ASV_51        :0.100
#               ASV_54/ASV_153:0.113
#               ASV_61        :0.046
#               ASV_92        :0.159
#               ASV_82        :0.188
#               ASV_85        :0.104

#=============================================================================================================
# --- Reproducibility
#=============================================================================================================

##----WRrepeats (73x2;1x3) and BRrepeats (28x2)

table(labdata$WRrepeats, exclude=NULL)
WRrep.IDs <- labdata$sample_id_run[!is.na(labdata$WRrepeats)]
WRrep.labels <- labdata$WRrepeats[!is.na(labdata$WRrepeats)]

get.Rsq <- function(rep.IDs, rep.labels)
{
  Rsq <- NULL
  label.vec <- levels(factor(rep.labels))
  for (i in 1:length(label.vec))
  {
    this.rep <- label.vec[i]
    these.IDs <- rep.IDs[rep.labels==this.rep]
    while (length(these.IDs)>1)
    {
      count1 <- seqmat[,match(these.IDs[1],colnames(seqmat))]
      count2 <- seqmat[,match(these.IDs[2],colnames(seqmat))]
      if (sd(count1)==0) print (c(this.rep, these.IDs[1],sum(count1)))
      if (sd(count2)==0) print (c(this.rep, these.IDs[2],sum(count2)))
      
      Rsq <- rbind (Rsq, c(cor(as.vector(count1),as.vector(count2))^2,
                           labdata$age_years_numeric[match(these.IDs[1],labdata$sample_id_run)],
                           labdata$templateDNA_16SqPCR_copies[match(these.IDs[1],labdata$sample_id_run)],
                           labdata$templateDNA_16SqPCR_copies[match(these.IDs[2],labdata$sample_id_run)],
                           labdata$FinalReads_ASVpipeline[match(these.IDs[1],labdata$sample_id_run)],
                           labdata$FinalReads_ASVpipeline[match(these.IDs[2],labdata$sample_id_run)]))
      rownames(Rsq)[nrow(Rsq)] <- this.rep
      these.IDs <- these.IDs[-1]
    }
  }
  colnames(Rsq) <- c("R^2","age", "copies1", "copies2", "reads1", "reads2")
  Rsq
}

WR.Rsq <- as.data.frame(get.Rsq (WRrep.IDs, WRrep.labels))

# === List of SampleIDs WRrepeats:
#write.csv(WRrep.IDs, paste(pathname, "SampleIDs_WRrepeats.csv", sep=""))
#write.csv(WR.Rsq, paste(pathname, "Output_WRrepeats.csv", sep=""))

  table(labdata$BRrepeats, labdata$run, exclude=NULL)
  table(labdata$BRrepeats, exclude=NULL)
  BRrep.IDs <- labdata$sample_id_run[!is.na(labdata$BRrepeats)]
  BRrep.labels <- labdata$BRrepeats[!is.na(labdata$BRrepeats)]
  BR.Rsq <- as.data.frame(get.Rsq (BRrep.IDs, BRrep.labels))

# === List of SampleIDs BRrepeats:
#write.csv(BRrep.IDs, paste(pathname, "HSampleIDs_BRrepeats.csv", sep=""))
#write.csv(BR.Rsq, paste(pathname, "Output_BRrepeats.csv", sep=""))

Rsq.plots <- function (WR.mat, BR.mat, age.lim=range(c(WR.mat$age,BR.mat$age)), 
                       copy.lim=range(c(WR.mat$copies1,WR.mat$copies2,BR.mat$copies1,BR.mat$copies2)),
                       read.lim=range(c(WR.mat$reads1,WR.mat$reads2,BR.mat$reads1,BR.mat$reads2)))
{
  # R.sq vs age
  dev.new()
  plot (rbind(WR.mat,BR.mat)[,c(2,1)], type="n", xlab="Age days numeric", ylab="Reproducibility r-square", main="Repeats", xlim=age.lim)
  points (WR.mat$age, WR.mat[,1], col="deepskyblue", pch=1)
  points (BR.mat$age, BR.mat[,1], col="violetred1", pch=1)
  legend ("bottomright", c("WR repeats", "BR repeats"), pch=1, col=c("deepskyblue","violetred1"))
  
  # R.sq vs copies
  dev.new()
  plot (range(rbind(WR.mat,BR.mat)[,3:4]), range(rbind(WR.mat,BR.mat)[,1], na.rm=T), 
        type="n", xlab="Copy number", ylab="Reproducibility r-square", main="Repeats", xlim=copy.lim)
  points (WR.mat$copies1, WR.mat[,1], col="deepskyblue", pch=1)
  points (WR.mat$copies2, WR.mat[,1], col="deepskyblue", pch=1)
  points (BR.mat$copies1, BR.mat[,1], col="violetred1", pch=1)
  points (BR.mat$copies2, BR.mat[,1], col="violetred1", pch=1)
  for (i in 1:nrow(WR.mat))
    lines (WR.mat[i,3:4],rep(WR.mat[i,1],2),col="deepskyblue")
  for (i in 1:nrow(BR.mat))
    lines (BR.mat[i,3:4],rep(BR.mat[i,1],2),col="violetred1")
  legend ("bottomright", c("WR repeats", "BR repeats"), pch=1, col=c("deepskyblue","violetred1"))
  
  # R.sq vs final reads
  dev.new()
  plot (range(rbind(WR.mat,BR.mat)[,5:6]), range(rbind(WR.mat,BR.mat)[,1], na.rm=T), 
        type="n", xlab="Final reads", ylab="Reproducibility r-square", main="Repeats", xlim=read.lim)
  points (WR.mat$reads1, WR.mat[,1], col="deepskyblue", pch=1)
  points (WR.mat$reads2, WR.mat[,1], col="deepskyblue", pch=1)
  points (BR.mat$reads1, BR.mat[,1], col="violetred1", pch=1)
  points (BR.mat$reads2, BR.mat[,1], col="violetred1", pch=1)
  for (i in 1:nrow(WR.mat))
    lines (WR.mat[i,5:6],rep(WR.mat[i,1],2),col="deepskyblue")
  for (i in 1:nrow(BR.mat))
    lines (BR.mat[i,5:6],rep(BR.mat[i,1],2),col="violetred1")
  legend ("bottomright", c("WR repeats", "BR repeats"), pch=1, col=c("deepskyblue","violetred1"))
}

Rsq.plots (WR.Rsq, BR.Rsq)

#which.WR <- labdata[match(WRrep.IDs,labdata$sample_id_run),]
#selected.WR <- which.WR$sample_id_run[which.WR$age_years_numeric > 10 & 
#                                        which.WR$templateDNA_16SqPCR_copies>300 & 
#                                        which.WR$FinalReads_ASVpipeline>1000]
#WRrep.labels.2 <- WRrep.labels[match(selected.WR,WRrep.IDs)]
#WRrep.IDs.2 <- WRrep.IDs[match(selected.WR,WRrep.IDs)]
#WR.Rsq.2 <- as.data.frame(get.Rsq (WRrep.IDs.2, WRrep.labels.2))

#which.BR <- labdata[match(BRrep.IDs,labdata$sample_id_run),]
#selected.BR <- which.BR$sample_id_run[which.BR$age_years_numeric > 10 & 
#                                        which.BR$templateDNA_16SqPCR_copies>300 & 
#                                        which.BR$FinalReads_ASVpipeline>1000]
#BRrep.labels.2 <- BRrep.labels[match(selected.BR,BRrep.IDs)]
#BRrep.IDs.2 <- BRrep.IDs[match(selected.BR,BRrep.IDs)]
#BR.Rsq.2 <- as.data.frame(get.Rsq (BRrep.IDs.2, BRrep.labels.2))

#Rsq.plots (WR.Rsq.2, BR.Rsq.2)

# === List of SampleIDs WRrepeats and BRrepeats after removing "problem samples":
#write.csv(WRrep.IDs.2, paste(pathname, "SampleIDs_WRrepeats_noproblems.csv", sep=""))
#write.csv(WR.Rsq.2, paste(pathname, "Output_WRrepeats_noproblems.csv", sep=""))
#write.csv(BRrep.IDs.2, paste(pathname, "SampleIDs_BRrepeats_noproblems.csv", sep=""))
#write.csv(BR.Rsq.2, paste(pathname, "Output_BRrepeats_noproblems.csv", sep=""))



##----WRrepeats: WRrepeats_WP (54X2) and WRrepeats_BP (20X2)

table(labdata$WRrepeats_WP, exclude=NULL)
WPrep.IDs <- labdata$sample_id_run[!is.na(labdata$WRrepeats_WP)]
WPrep.labels <- labdata$WRrepeats_WP[!is.na(labdata$WRrepeats_WP)]

get.Rsq <- function(rep.IDs, rep.labels)
{
  Rsq <- NULL
  label.vec <- levels(factor(rep.labels))
  for (i in 1:length(label.vec))
  {
    this.rep <- label.vec[i]
    these.IDs <- rep.IDs[rep.labels==this.rep]
    while (length(these.IDs)>1)
    {
      count1 <- seqmat[,match(these.IDs[1],colnames(seqmat))]
      count2 <- seqmat[,match(these.IDs[2],colnames(seqmat))]
      if (sd(count1)==0) print (c(this.rep, these.IDs[1],sum(count1)))
      if (sd(count2)==0) print (c(this.rep, these.IDs[2],sum(count2)))
      
      Rsq <- rbind (Rsq, c(cor(as.vector(count1),as.vector(count2))^2,
                           labdata$age_years_numeric[match(these.IDs[1],labdata$sample_id_run)],
                           labdata$templateDNA_16SqPCR_copies[match(these.IDs[1],labdata$sample_id_run)],
                           labdata$templateDNA_16SqPCR_copies[match(these.IDs[2],labdata$sample_id_run)],
                           labdata$FinalReads_ASVpipeline[match(these.IDs[1],labdata$sample_id_run)],
                           labdata$FinalReads_ASVpipeline[match(these.IDs[2],labdata$sample_id_run)]))
      rownames(Rsq)[nrow(Rsq)] <- this.rep
      these.IDs <- these.IDs[-1]
    }
  }
  colnames(Rsq) <- c("R^2","age", "copies1", "copies2", "reads1", "reads2")
  Rsq
}

WP.Rsq <- as.data.frame(get.Rsq (WPrep.IDs, WPrep.labels))

# === List of SampleIDs WRrepeats_WP:
#write.csv(WPrep.IDs, paste(pathname, "SampleIDs_WPrepeats.csv", sep=""))
#write.csv(WP.Rsq, paste(pathname, "Output_WPrepeats.csv", sep=""))

table(labdata$WRrepeats_BP, labdata$run, exclude=NULL)
table(labdata$WRrepeats_BP, exclude=NULL)
BPrep.IDs <- labdata$sample_id_run[!is.na(labdata$WRrepeats_BP)]
BPrep.labels <- labdata$WRrepeats_BP[!is.na(labdata$WRrepeats_BP)]
BP.Rsq <- as.data.frame(get.Rsq (BPrep.IDs, BPrep.labels))

# === List of SampleIDs WRrepeats_BP:
#write.csv(BPrep.IDs, paste(pathname, "SampleIDs_BPrepeats.csv", sep=""))
#write.csv(BP.Rsq, paste(pathname, "Output_BPrepeats.csv", sep=""))

Rsq.plots <- function (WP.mat, BP.mat, age.lim=range(c(WP.mat$age,BP.mat$age)), 
                       copy.lim=range(c(WP.mat$copies1,WP.mat$copies2,BP.mat$copies1,BP.mat$copies2)),
                       read.lim=range(c(WP.mat$reads1,WP.mat$reads2,BP.mat$reads1,BP.mat$reads2)))
{
  # R.sq vs age
  dev.new()
  plot (rbind(WP.mat,BP.mat)[,c(2,1)], type="n", xlab="Age days numeric", ylab="Reproducibility r-square", main="Repeats", xlim=age.lim)
  points (WP.mat$age, WP.mat[,1], col="mediumorchid4", pch=1)
  points (BP.mat$age, BP.mat[,1], col="violetred1", pch=1)
  legend ("bottomright", c("WP repeats", "BP repeats"), pch=1, col=c("mediumorchid4","violetred1"))
  
  # R.sq vs copies
  dev.new()
  plot (range(rbind(WP.mat,BP.mat)[,3:4]), range(rbind(WP.mat,BP.mat)[,1], na.rm=T), 
        type="n", xlab="Copy number", ylab="Reproducibility r-square", main="Repeats", xlim=copy.lim)
  points (WP.mat$copies1, WP.mat[,1], col="mediumorchid4", pch=1)
  points (WP.mat$copies2, WP.mat[,1], col="mediumorchid4", pch=1)
  points (BP.mat$copies1, BP.mat[,1], col="violetred1", pch=1)
  points (BP.mat$copies2, BP.mat[,1], col="violetred1", pch=1)
  for (i in 1:nrow(WP.mat))
    lines (WP.mat[i,3:4],rep(WP.mat[i,1],2),col="mediumorchid4")
  for (i in 1:nrow(BP.mat))
    lines (BP.mat[i,3:4],rep(BP.mat[i,1],2),col="violetred1")
  legend ("bottomright", c("WP repeats", "BP repeats"), pch=1, col=c("mediumorchid4","violetred1"))
  
  # R.sq vs final reads
  dev.new()
  plot (range(rbind(WP.mat,BP.mat)[,5:6]), range(rbind(WP.mat,BP.mat)[,1], na.rm=T), 
        type="n", xlab="Final reads", ylab="Reproducibility r-square", main="Repeats", xlim=read.lim)
  points (WP.mat$reads1, WP.mat[,1], col="mediumorchid4", pch=1)
  points (WP.mat$reads2, WP.mat[,1], col="mediumorchid4", pch=1)
  points (BP.mat$reads1, BP.mat[,1], col="violetred1", pch=1)
  points (BP.mat$reads2, BP.mat[,1], col="violetred1", pch=1)
  for (i in 1:nrow(WP.mat))
    lines (WP.mat[i,5:6],rep(WP.mat[i,1],2),col="mediumorchid4")
  for (i in 1:nrow(BP.mat))
    lines (BP.mat[i,5:6],rep(BP.mat[i,1],2),col="violetred1")
  legend ("bottomright", c("WP repeats", "BP repeats"), pch=1, col=c("mediumorchid4","violetred1"))
}

Rsq.plots (WP.Rsq, BP.Rsq)

#which.WP <- labdata[match(WPrep.IDs,labdata$sample_id_run),]
#selected.WP <- which.WP$sample_id_run[which.WP$age_years_numeric > 10 & 
#                                        which.WP$templateDNA_16SqPCR_copies>300 & 
#                                        which.WP$FinalReads_ASVpipeline>1000]
#WPrep.labels.2 <- WPrep.labels[match(selected.WP,WPrep.IDs)]
#WPrep.IDs.2 <- WPrep.IDs[match(selected.WP,WPrep.IDs)]
#WP.Rsq.2 <- as.data.frame(get.Rsq (WPrep.IDs.2, WPrep.labels.2))

#which.BP <- labdata[match(BPrep.IDs,labdata$sample_id_run),]
#selected.BP <- which.BP$sample_id_run[which.BP$age_years_numeric > 10 & 
#                                        which.BP$templateDNA_16SqPCR_copies>300 & 
#                                        which.BP$FinalReads_ASVpipeline>1000]
#BPrep.labels.2 <- BPrep.labels[match(selected.BP,BPrep.IDs)]
#BPrep.IDs.2 <- BPrep.IDs[match(selected.BP,BPrep.IDs)]
#BP.Rsq.2 <- as.data.frame(get.Rsq (BPrep.IDs.2, BPrep.labels.2))

#Rsq.plots (WP.Rsq.2, BP.Rsq.2)

# === List of SampleIDs WRrepeats_WP and WRrepeats_BP after removing "problem samples":
#write.csv(WPrep.IDs.2, paste(pathname, "SampleIDs_WPrepeats_noproblems.csv", sep=""))
#write.csv(WP.Rsq.2, paste(pathname, "Output_WPrepeats_noproblems.csv", sep=""))
#write.csv(BPrep.IDs.2, paste(pathname, "SampleIDs_BPrepeats_noproblems.csv", sep=""))
#write.csv(BP.Rsq.2, paste(pathname, "Output_BPrepeats_noproblems.csv", sep=""))


###Comment (SC): 
# I would not suggest to exclude any samples based on WRrepeats and BRrepeats reproducibility
# All, except 1 repeat had R2 values >0.92

#=============================================================================================================
# --- labdata analysis for biological samples and negative controls
#=============================================================================================================

# --- identify samples and negative controls
seq.samples <- seqmat[,match(samples,colnames(seqmat))] # 3481 ASVs x 960 specimens
dim(seq.samples)
seq.negcontrols <- seqmat[,match(negcontrols,colnames(seqmat))] # 3481 ASVs x 43 specimens
dim(seq.negcontrols)
seq.comb <- cbind (seq.samples, seq.negcontrols) # 3481 ASVs x 1003 specimens
dim(seq.comb)
seq.comb <- seq.comb[apply(seq.comb,1,sum)>0,] # 3219 ASVs>0 x 1003 specimens
dim(seq.comb)
seq.samples <- seq.samples[match(rownames(seq.comb),rownames(seq.samples)),] # 3219 x 960
dim(seq.samples)
seq.negcontrols <- seq.negcontrols[match(rownames(seq.comb),rownames(seq.negcontrols)),] # 3219 x 43
dim(seq.negcontrols)
labdata2 <- labdata[match(colnames(seq.comb),labdata$sample_id_run),] # 1003 x 60
dim(labdata2)
rm("seq.comb")

n <- nrow(labdata2)
copy.cutoff <- 300

# --- copies vs final reads
plot (c(1,n),range(labdata2$FinalReads_ASVpipeline,na.rm=T), type="n", xlab="copies", ylab="Final reads")
points((1:n)[labdata2$templateDNA_16SqPCR_copies>=copy.cutoff], labdata2$FinalReads_ASVpipeline[labdata2$templateDNA_16SqPCR_copies>=copy.cutoff], 
       col="deepskyblue", pch=c(17,1)[as.numeric(is.na(labdata2$nontemplate[labdata2$templateDNA_16SqPCR_copies>=copy.cutoff]))+1])
points((1:n)[labdata2$templateDNA_16SqPCR_copies<copy.cutoff], labdata2$FinalReads_ASVpipeline[labdata2$templateDNA_16SqPCR_copies<copy.cutoff], 
       col="red4", pch=c(17,1)[as.numeric(is.na(labdata2$nontemplate[labdata2$templateDNA_16SqPCR_copies<copy.cutoff]))+1])
legend("topleft",c(paste("copies <",copy.cutoff,sep=""),paste("copies >=",copy.cutoff,sep="")), col=c("red4","deepskyblue"),pch=18)
legend("topright",c("biological samples","negative controls"),col="black",pch=c(1,17))

plot (range(labdata2$templateDNA_16SqPCR_copies,na.rm=T), range(labdata2$FinalReads_ASVpipeline,na.rm=T), type="n", log="x", xlab="copies", ylab="Final reads")
points (labdata2$templateDNA_16SqPCR_copies,labdata2$FinalReads_ASVpipeline,col=c("red4","deepskyblue")[as.numeric(is.na(labdata2$nontemplate))+1])
panel.smooth(labdata2$templateDNA_16SqPCR_copies[is.na(labdata2$nontemplate)],labdata2$FinalReads_ASVpipeline[is.na(labdata2$nontemplate)],col="black")
cor(labdata2$templateDNA_16SqPCR_copies[is.na(labdata2$nontemplate)],labdata2$FinalReads_ASVpipeline[is.na(labdata2$nontemplate)])

# --- copies vs diversity
require(vegan)
div.vals <- diversity(t(cbind(seq.samples,seq.negcontrols)))
temp.mat <- cbind(seq.samples,seq.negcontrols)
temp.mat <- apply(temp.mat, 2, function(x)x/sum(x))
dim(temp.mat)
div2 <- diversity(t(temp.mat))
plot (div.vals, div2)
rm("temp.mat")
plot (c(1,n),range(div.vals,na.rm=T), type="n", xlab="copies", ylab="Shannon diversity")
points((1:n)[labdata2$templateDNA_16SqPCR_copies>=copy.cutoff], div.vals[labdata2$templateDNA_16SqPCR_copies>=copy.cutoff], 
       col="deepskyblue", pch=c(17,1)[as.numeric(is.na(labdata2$nontemplate[labdata2$templateDNA_16SqPCR_copies>=copy.cutoff]))+1])
points((1:n)[labdata2$templateDNA_16SqPCR_copies<copy.cutoff], div.vals[labdata2$templateDNA_16SqPCR_copies<copy.cutoff], 
       col="red4", pch=c(17,1)[as.numeric(is.na(labdata2$nontemplate[labdata2$templateDNA_16SqPCR_copies<copy.cutoff]))+1])
legend("topleft",c(paste("copies <",copy.cutoff,sep=""),paste("copies >=",copy.cutoff,sep="")), col=c("red4","deepskyblue"),pch=18)
legend("topright",c("biological samples","negative controls"),col="black",pch=c(1,17))

plot (range(labdata2$templateDNA_16SqPCR_copies,na.rm=T), range(div.vals,na.rm=T), type="n", log="x", xlab="copies", ylab="Shannon diversity")
points (labdata2$templateDNA_16SqPCR_copies,div.vals,col=c("red4","deepskyblue")[as.numeric(is.na(labdata2$nontemplate))+1])
panel.smooth(labdata2$templateDNA_16SqPCR_copies[is.na(labdata2$nontemplate)],div.vals[is.na(labdata2$nontemplate)],col="black")
cor(labdata2$templateDNA_16SqPCR_copies[is.na(labdata2$nontemplate)],div.vals[is.na(labdata2$nontemplate)])

# --- copies vs age
plot (c(1,n),range(labdata2$age_years_numeric,na.rm=T), type="n", xlab="copies", ylab="age")
points((1:n)[labdata2$templateDNA_16SqPCR_copies>=copy.cutoff], labdata2$age_years_numeric[labdata2$templateDNA_16SqPCR_copies>=copy.cutoff], 
       col="deepskyblue", pch=c(17,1)[as.numeric(is.na(labdata2$nontemplate[labdata2$templateDNA_16SqPCR_copies>=copy.cutoff]))+1])
points((1:n)[labdata2$templateDNA_16SqPCR_copies<copy.cutoff], labdata2$age_years_numeric[labdata2$templateDNA_16SqPCR_copies<copy.cutoff], 
       col="red4", pch=c(17,1)[as.numeric(is.na(labdata2$nontemplate[labdata2$templateDNA_16SqPCR_copies<copy.cutoff]))+1])
legend("topleft",c(paste("copies <",copy.cutoff,sep=""),paste("copies >=",copy.cutoff,sep="")), col=c("red4","deepskyblue"),pch=18)
legend("topright",c("biological samples","negative controls"),col="black",pch=c(1,17))

plot (range(labdata2$templateDNA_16SqPCR_copies,na.rm=T), range(labdata2$age_years_numeric,na.rm=T), type="n", log="x", xlab="copies", ylab="age")
points (labdata2$templateDNA_16SqPCR_copies,labdata2$age_years_numeric,col=c("red4","deepskyblue")[as.numeric(is.na(labdata2$nontemplate))+1])
panel.smooth(labdata2$templateDNA_16SqPCR_copies[is.na(labdata2$nontemplate)],labdata2$age_years_numeric[is.na(labdata2$nontemplate)],col="black")
cor(labdata2$templateDNA_16SqPCR_copies[is.na(labdata2$nontemplate)],labdata2$age_years_numeric[is.na(labdata2$nontemplate)])

# ===================================================================================================================
# log ratio biplots

LRbiplot.calc <- function (datmat)
{
  log.zero.adjust <- function(c.vec)
  { # uniform replacement
    k <- sum(c.vec==0)
    c.vec[c.vec==0] <- runif(k, 0.1, 1)
    c.vec
  }
  
  N    <- t(apply (datmat, 1, log.zero.adjust))
  P    <- N / sum(N)
  rm   <- apply(P, 1, sum)
  cm   <- apply(P, 2, sum)
  Y    <- as.matrix(log(P))
  mc   <- t(Y) %*% as.vector(rm)
  Y    <- Y - rep(1, nrow(P)) %*% t(mc)
  mr   <- Y %*% as.vector(cm)
  Y    <- Y - mr %*% t(rep(1,ncol(P)))
  Z    <- diag(sqrt(rm)) %*% Y %*% diag(sqrt(cm))
  svdZ <- svd(Z)
  Fmat <- diag(1/sqrt(rm)) %*% svdZ$u[,1:2] %*% diag(svdZ$d[1:2])
  Gmat <- diag(1/sqrt(cm)) %*% svdZ$v[,1:2]
  # --- implement lambda-scaling: Gower, Lubbe & le Roux, section 2.3.1
  lambda <- ((nrow(Fmat)/nrow(Gmat)) * (sum(Gmat^2)/sum(Fmat^2)))^0.25
  Fmat <- lambda*Fmat
  Gmat <- Gmat/lambda
  list (rows=Fmat, cols=Gmat, qual=svdZ$d^2/sum(svdZ$d^2))
}

PCArarefy.calc <- function (countmat)
{
  require(vegan)
  mat <- countmat[apply(countmat,1,sum)>1000,]
  rmat <- rrarefy (mat, 1000)
  
  X <- scale(rmat,scale=F)
  svd.out <- svd(X)
  PCscores <- svd.out$u[,1:2] %*% diag(svd.out$d[1:2])
  Gmat <- svd.out$v[,1:2]
  Fmat <- matrix (NA, nrow=nrow(countmat),ncol=2)
  rownames(Fmat) <- rownames(countmat)
  Fmat[match(rownames(rmat),rownames(Fmat)),] <- PCscores
  list (rows=Gmat, cols=Fmat)
}

PCoArarefy.calc <- function (countmat)
{
  require(vegan)
  mat <- countmat[apply(countmat,1,sum)>1000,]
  rmat <- rrarefy (mat, 1000)
  
  DD <- vegdist(rmat)
  MDSscores <- cmdscale(DD)
  Gmat <- NULL
  Fmat <- matrix (NA, nrow=nrow(countmat),ncol=2)
  rownames(Fmat) <- rownames(countmat)
  Fmat[match(rownames(rmat),rownames(Fmat)),] <- MDSscores
  list (rows=Gmat, cols=Fmat)
}

PCoAbeta.calc <- function (countmat, do.rarefy=F)
{
  require(vegan)
  if (!do.rarefy) mat <- countmat
  else
  {
    mat <- countmat[apply(countmat,1,sum)>1000,]
    mat <- rrarefy (mat, 1000)
  }
  
  DD <- betadiver(mat,method="w")
  MDSscores <- cmdscale(DD)
  Gmat <- NULL
  Fmat <- matrix (NA, nrow=nrow(countmat),ncol=2)
  rownames(Fmat) <- rownames(countmat)
  Fmat[match(rownames(mat),rownames(Fmat)),] <- MDSscores
  list (rows=Gmat, cols=Fmat)
}

LR.group.plot <- function(mat, cat.vec, col.vec, sample.type, ...)
{
  plot (range(mat[,1],na.rm=T),range(mat[,2],na.rm=T), type="n", xaxt="n", yaxt="n", xlab="", ylab="", ...)
  for (i in nlevels(cat.vec):1)
    points (mat[as.numeric(cat.vec)==i,1], mat[as.numeric(cat.vec)==i,2], 
            col=col.vec[i], pch=ifelse(is.na(sample.type[as.numeric(cat.vec)==i]),1,17))
  if (any(is.na(cat.vec))) points (mat[is.na(cat.vec),1], mat[is.na(cat.vec),2], col="black", pch=17)
   legend ("bottomright", names(col.vec), col=col.vec, pch=15, cex=0.7)
   legend ("bottomleft", c("biological samples","negative controls"), pch=c(1,17), col="black")
}

calc.rarecurve <- function (x, step = 1, sample) 
{
  ### --- based on rarecurve() function from package vegan
  ### --- to return computations without plotting
  
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  out <- lapply(seq_len(nr), function(i) { n <- seq(1, tot[i], by = step)
  if (n[length(n)] != tot[i]) 
    n <- c(n, tot[i])
  drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  out2 <- sapply(out, function(unit) { mat <- cbind(attr(unit,"Subsample"),unit)
  rownames(mat) <- NULL
  colnames(mat) <- c("reads","ASVs.OTUs")
  mat
  })
  names(out2) <- rownames(x)
  out2 
}

plot.rarecurve <- function(rare.out, col.vec=NULL, ...)
{
  n <- length(rare.out)
  if (is.null(col.vec)) col.vec <- rep("black",n)
  while (length(col.vec)<n) col.vec <- c(col.vec, col.vec)
  x.range <- range(sapply(rare.out, function(mat) range(mat[,1],na.rm=T)))
  y.range <- range(sapply(rare.out, function(mat) range(mat[,2],na.rm=T)))
  plot (x.range, y.range, type="n", ...)
  for (i in 1:n)
    lines (rare.out[[i]], col=col.vec[i])
}

my.rarefaction <- function(count.table, cat.vec, col.vec, step=100, ...)
{
  rare.out <- calc.rarecurve(t(count.table), step=step)
  levels(cat.vec) <- c(levels(cat.vec),"missing")
  cat.vec[is.na(cat.vec)] <- "missing"
  temp.table <- table(cat.vec)
  col.vec <- col.vec[temp.table>0]
  cat.labels <- names(temp.table[temp.table>0])
  n <- length(cat.labels)
  dev.new()
  #  row.num <- floor(sqrt(n))
  #  col.num <- ceiling(n/row.num)
  #  par(mfrow=c(row.num, col.num))
  m <- length(rare.out)
  
  for (k in 1:n)
  {
    j <- 1
    sub.list <- NULL
    for (i in 1:m)
      if (cat.vec[i]==cat.labels[k]) { sub.list[[j]] <- rare.out[[i]]
      j <- j+1
      }
    plot.rarecurve(sub.list, col=col.vec[k], ...)
    lines (rep(1000,2),par("usr")[3:4])
  }
  
}

LR <- LRbiplot.calc(cbind(seq.samples, seq.negcontrols))$cols 
PCA <- PCArarefy.calc(t(cbind(seq.samples, seq.negcontrols)))$cols 
PCoA <- PCoArarefy.calc(t(cbind(seq.samples, seq.negcontrols)))$cols

# --- LR biplot for copies
#extra cols options : "red4","red1","orange","gold","yellow","greenyellow","lawngreen","seagreen3","lightseagreen",
#                     "deepskyblue","royalblue","slateblue4","navyblue","midnightblue","black")
#extra cols options: "#920023","#F6815C","#E7F3F9","#BCDDEE","#85B8D7","#3D669D","#242872","#12143A","#05061c"

range(labdata2$templateDNA_16SqPCR_copies, na.rm=T)
copy.breaks <- c(0,300,500,1000,2000,5000,10000,50000,100000,10000000)
copy.cols <- c("red4","red1","orange","yellow","lawngreen","lightseagreen",
               "deepskyblue","royalblue","midnightblue")
copy.cat <- cut(labdata2$templateDNA_16SqPCR_copies,copy.breaks)
table(copy.cat)
names(copy.cols) <- levels(copy.cat)
names(copy.cols) <- c("<300",">=300-500",">=500-1,000",">=1,000-2,000",">=2,000-5,000",
                      ">=5,000-10,000",">=10,000-50,000",">=50,000-100,000",">=100,000")

require(vegan)

# --- LR biplot for Final Reads
range(labdata2$FinalReads_ASVpipeline, na.rm=T)
FinalRead.breaks <- c(0,1000,2000,3000,4000,5000,10000,20000,50000,100000)
FinalRead.cols <- c("red4","red1","orange","yellow","lawngreen","lightseagreen",
                    "deepskyblue","royalblue","midnightblue")
reads.cat <- cut(labdata2$FinalReads_ASVpipeline,FinalRead.breaks, include.lowest=T)
table(reads.cat)
names(FinalRead.cols) <- c("<1,000",">=1,000-2,000",">=2,000-3,000",">=3,000-4,000",">=4,000-5,000",">=5,000-10,000",
                           ">=10,000-20,000",">=20,000-50,000",">=50,000")

# --- LR biplot for age
range(labdata2$age_years_numeric, na.rm=T)
age.breaks <- c(6,7,8,9,10,11,12,13,14,15,16,17,18,22)
age.cols <- c("red4","red1","orange","yellow","greenyellow","lawngreen","seagreen3","lightseagreen",
              "deepskyblue","royalblue","slateblue4","navyblue","midnightblue")
age.cat <- cut(labdata2$age_years_numeric,age.breaks)
table(age.cat)
names(age.cols) <- levels(age.cat)
names(age.cols) <- c("6-7",">7-8",">8-9",">9-10",">10-11",">11-12",
                     "12-13",">13-14",">14-15",">15-16",">16-17",">17-18",
                     "18-22")

# --- LR biplot for Primestore vs SP
table(labdata2$sample_type, exclude=NULL) # 960 + 43
sample.type.cols <- c("deepskyblue","violetred")
names(sample.type.cols) <- levels(labdata2$sample_type)

# --- LR biplot for Runs
table(labdata2$run, exclude=NULL)
run.cols <- c("red","deepskyblue","goldenrod1")
names(run.cols) <- levels(labdata2$run)

# --- LR biplot for country
table(labdata2$country, exclude=NULL)
country.cols <- c("orange","purple","black")
names(country.cols) <- levels(labdata2$country)

# --- LR biplot for timepoint
table(labdata2$timepoint, exclude=NULL)
timepoint.cols <- c("orchid1","orchid4","plum1","grey","black")
names(timepoint.cols) <- levels(labdata2$timepoint)

# --- LR biplot for sampling_induced_expectorated
table(labdata2$sampling_induced_expectorated, exclude=NULL)
sampling_induced_expectorated.cols <- c("mediumseagreen","red","grey","black")
names(sampling_induced_expectorated.cols) <- levels(labdata2$sampling_induced_expectorated)


LR.group.plot (mat=LR, cat.vec=copy.cat, col.vec=copy.cols, sample.type=labdata2$nontemplate, main="Copies")
LR.group.plot (mat=LR, cat.vec=reads.cat, col.vec=FinalRead.cols, sample.type=labdata2$nontemplate, main="Final Reads")
LR.group.plot (mat=LR, cat.vec=age.cat, col.vec=age.cols, sample.type=labdata2$nontemplate, main="Age")
LR.group.plot (mat=LR, cat.vec=labdata2$sample_type, col.vec=sample.type.cols, sample.type=labdata2$nontemplate, main="Primestore vs SP")
LR.group.plot (mat=LR, cat.vec=labdata2$run, col.vec=run.cols, sample.type=labdata2$nontemplate, main="run")
LR.group.plot (mat=LR, cat.vec=labdata2$country, col.vec=country.cols, sample.type=labdata2$nontemplate, main="country")
LR.group.plot (mat=LR, cat.vec=labdata2$timepoint, col.vec=timepoint.cols, sample.type=labdata2$nontemplate, main="timepoint")
LR.group.plot (mat=LR, cat.vec=labdata2$sampling_induced_expectorated, col.vec=sampling_induced_expectorated.cols, sample.type=labdata2$nontemplate, main="sampling_induced_expectorated")

LR.group.plot (mat=PCA, cat.vec=copy.cat, col.vec=copy.cols, sample.type=labdata2$nontemplate, main="Copies rarefied PCA")
LR.group.plot (mat=PCA, cat.vec=reads.cat, col.vec=FinalRead.cols, sample.type=labdata2$nontemplate, main="Final Reads rarefied PCA")
LR.group.plot (mat=PCA, cat.vec=age.cat, col.vec=age.cols, sample.type=labdata2$nontemplate, main="Age rarefied PCA")
LR.group.plot (mat=PCA, cat.vec=labdata2$sample_type, col.vec=sample.type.cols, sample.type=labdata2$nontemplate, main="Primestore vs SP rarefied PCA")
LR.group.plot (mat=PCA, cat.vec=labdata2$run, col.vec=run.cols, sample.type=labdata2$nontemplate, main="run rarefied PCA")
LR.group.plot (mat=PCA, cat.vec=labdata2$country, col.vec=country.cols, sample.type=labdata2$nontemplate, main="country rarefied PCA")
LR.group.plot (mat=PCA, cat.vec=labdata2$timepoint, col.vec=timepoint.cols, sample.type=labdata2$nontemplate, main="timepoint rarefied PCA")
LR.group.plot (mat=PCA, cat.vec=labdata2$sampling_induced_expectorated, col.vec=sampling_induced_expectorated.cols, sample.type=labdata2$nontemplate, main="sampling_induced_expectorated rarefied PCA")

LR.group.plot (mat=PCoA, cat.vec=copy.cat, col.vec=copy.cols, sample.type=labdata2$nontemplate, main="Copies rarefied PCoA")
LR.group.plot (mat=PCoA, cat.vec=reads.cat, col.vec=FinalRead.cols, sample.type=labdata2$nontemplate, main="Final Reads rarefied PCoA")
LR.group.plot (mat=PCoA, cat.vec=age.cat, col.vec=age.cols, sample.type=labdata2$nontemplate, main="Age rarefied PCoA")
LR.group.plot (mat=PCoA, cat.vec=labdata2$sample_type, col.vec=sample.type.cols, sample.type=labdata2$nontemplate, main="Primestore vs SP rarefied PCoA")
LR.group.plot (mat=PCoA, cat.vec=labdata2$run, col.vec=run.cols, sample.type=labdata2$nontemplate, main="run rarefied PCoA")
LR.group.plot (mat=PCoA, cat.vec=labdata2$country, col.vec=country.cols, sample.type=labdata2$nontemplate, main="country rarefied PCoA")
LR.group.plot (mat=PCoA, cat.vec=labdata2$timepoint, col.vec=timepoint.cols, sample.type=labdata2$nontemplate, main="timepoint rarefied PCoA")
LR.group.plot (mat=PCoA, cat.vec=labdata2$sampling_induced_expectorated, col.vec=sampling_induced_expectorated.cols, sample.type=labdata2$nontemplate, main="sampling_induced_expectorated rarefied PCoA")

my.rarefaction (cbind(seq.samples,seq.negcontrols), cat.vec=copy.cat, col.vec=copy.cols, 
                xlab="Final reads", ylab="ASVs", main="Copies", ylim=c(0,200))
my.rarefaction (cbind(seq.samples,seq.negcontrols), cat.vec=reads.cat, col.vec=FinalRead.cols, 
                xlab="Final reads", ylab="ASVs", main="Final reads", ylim=c(0,200))
my.rarefaction (cbind(seq.samples,seq.negcontrols), cat.vec=age.cat, col.vec=age.cols, 
                xlab="Final reads", ylab="ASVs", main="Age", ylim=c(0,200))
my.rarefaction (cbind(seq.samples,seq.negcontrols), cat.vec=labdata2$run, col.vec=run.cols, 
                xlab="Final reads", ylab="ASVs", main="Runs", ylim=c(0,200))

###Comment (SC): 
# I don't think we are seeing an effect of read numbers, 16S copies or age here
# Perhaps we might want to consider excluding samples with <1000 reads and <300 16S copies??


# ===================================== Low biomass additional investigation

dim(cbind(seq.samples, seq.negcontrols))    # 3219 ASVs x 1003 specimens
seq.comb <- cbind(seq.samples, seq.negcontrols)[,!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$templateDNA_16SqPCR_copies<10000)]
dim(seq.comb)      # 3219 ASVs x 157 low biomass specimens 
seq.comb <- seq.comb[apply(seq.comb,1,sum)>0,] 
dim(seq.comb) # 1359 ASVs>0 x 157 low biomass specimens
LR.low.bio <- LRbiplot.calc(seq.comb)$cols
PCA.low.bio <- PCArarefy.calc(t(seq.comb))$cols  
PCoA.low.bio <- PCoArarefy.calc(t(seq.comb))$cols  
labdata.low.bio <- labdata2[!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$templateDNA_16SqPCR_copies<10000),]
dim(labdata.low.bio)    # 157 low biomass specimens x 60

# --- LR biplot for copies

copy.breaks <- c(0,100,300,500,1000,2000,5000,10000)
copy.cols <- c("red4","orange","yellow","lawngreen","seagreen","lightseagreen","deepskyblue")
copy.cat <- cut(labdata.low.bio$templateDNA_16SqPCR_copies,copy.breaks)
table(copy.cat) 
sum(table(copy.cat)) # 157
names(copy.cols) <- levels(copy.cat)
names(copy.cols) <- c("<100",">=100-300",">=300-500",">=500-1000",">=1000-2000",">=2000-5000",">=5000-10000")

LR.group.plot (mat=LR.low.bio, cat.vec=copy.cat, col.vec=copy.cols, sample.type=labdata.low.bio$nontemplate, main="Copies log-ratio")
LR.group.plot (mat=PCA.low.bio, cat.vec=copy.cat, col.vec=copy.cols, sample.type=labdata.low.bio$nontemplate, main="Copies rarefied PCA")
LR.group.plot (mat=PCoA.low.bio, cat.vec=copy.cat, col.vec=copy.cols, sample.type=labdata.low.bio$nontemplate, main="Copies rarefied PCoA")

# ===================================== Low read count additional investigation

dim(cbind(seq.samples, seq.negcontrols))    # 3219 ASVs x 1003 specimens
seq.comb <- cbind(seq.samples, seq.negcontrols)[,!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$FinalReads_ASVpipeline<15000)]
dim(seq.comb)      # 3219 ASVs x 416 low read count specimens 
seq.comb <- seq.comb[apply(seq.comb,1,sum)>0,] 
dim(seq.comb) # 1692 ASVs>0 x 416 low read count specimens
LR.low.bio <- LRbiplot.calc(seq.comb)$cols
PCA.low.bio <- PCArarefy.calc(t(seq.comb))$cols  
PCoA.low.bio <- PCoArarefy.calc(t(seq.comb))$cols  
labdata.low.bio <- labdata2[!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$FinalReads_ASVpipeline<15000),]
dim(labdata.low.bio)    # 416 low biomass specimens x 60

# --- LR biplot for low read counts

read.breaks <- c(0,1000,2000,5000,10000,15000)
read.cols <- c("red4","orange","yellow","lawngreen","seagreen")
read.cat <- cut(labdata.low.bio$FinalReads_ASVpipeline,read.breaks)
table(read.cat) 
sum(table(read.cat)) # 394
names(read.cols) <- levels(read.cat)
names(read.cols) <- c("<1000",">=1000-2000",">=2000-5000",">=5000-10000",">=10000-15000")

LR.group.plot (mat=LR.low.bio, cat.vec=read.cat, col.vec=read.cols, sample.type=labdata.low.bio$nontemplate, main="Reads log-ratio")
LR.group.plot (mat=PCA.low.bio, cat.vec=read.cat, col.vec=read.cols, sample.type=labdata.low.bio$nontemplate, main="Reads rarefied PCA")
LR.group.plot (mat=PCoA.low.bio, cat.vec=read.cat, col.vec=read.cols, sample.type=labdata.low.bio$nontemplate, main="Reads rarefied PCoA")

# ===================================== Remove low reads and low 16S copy numbers

# === Remove samples based on age years (no samples to remove)
seq.comb <- cbind(seq.samples, seq.negcontrols)[,!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$age_years_numeric>1)]
dim(seq.comb)
# 3219 x 1003 (remove no samples based on age_years)
seq.comb <- seq.comb[apply(seq.comb,1,sum)>0,] # 3219 x 1003
labdata3 <- labdata2[!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$age_years_numeric>1),]
dim(labdata3)
# 1003 x 60

# === Remove copies<=100
seq.comb <- seq.comb[,!is.na(labdata3$nontemplate) | labdata3$templateDNA_16SqPCR_copies>100]
dim(seq.comb)
# 3219 x 1001 # remove 2 samples
seq.comb <- seq.comb[apply(seq.comb,1,sum)>0,] # 3216 x 1001
labdata3 <- labdata3[!is.na(labdata3$nontemplate) | labdata3$templateDNA_16SqPCR_copies>100,]
dim(labdata3)
# 1001 x 60

# === Remove reads<=1000
seq.comb <- seq.comb[,!is.na(labdata3$nontemplate) | labdata3$FinalReads_ASVpipeline>1000]
dim(seq.comb)
# 3216 x 996 # remove 5 samples
seq.comb <- seq.comb[apply(seq.comb,1,sum)>0,] # 3216 x 996
labdata3 <- labdata3[!is.na(labdata3$nontemplate) | labdata3$FinalReads_ASVpipeline>1000,]
dim(labdata3)
# 996 x 58 (60)


# ==============================================================================================================
# Barplots

top25.col <- colors()[c(552,131,90,97,501,142,653,83,87,657,497,614,610,102,8,72,45,128,26,30,84,544,450,524,40)]

# --- low copy number vs negative controls

dim(seq.samples)
dim(seq.negcontrols)
sample.copies <- labdata2$templateDNA_16SqPCR_copies[is.na(labdata2$nontemplate)]
dim(sample.copies)
seq.low.copy <- seq.samples[,sample.copies>100 & sample.copies<1000]
seq.low.copy <- seq.low.copy[,rev(order(sample.copies[sample.copies>100 & sample.copies<1000]))]
dim(seq.low.copy)
seq.comb <- cbind (seq.low.copy, seq.negcontrols)
seq.comb <- seq.comb[apply(seq.comb,1,sum)>0,]
dim(seq.comb)
rev(sort(apply(seq.comb,2,sum)))
seq.comb <- seq.comb[,apply(seq.comb,2,sum)>0]
prop.comb <- apply(seq.comb, 2, function(x)x/sum(x))
dim(prop.comb)

max.ra <- apply(prop.comb,1,max)
prop.other <- prop.comb[max.ra>=0.01,]
other <- apply(prop.comb[max.ra<0.01,],2,sum)
prop.other <- rbind (prop.other, other)
dim(prop.other)

rm("seq.low.copy")

tax.low.copy <- taxdat0[match(rownames(prop.other),rownames(taxdat0)),]
head(tax.low.copy)
table(tax.low.copy$Phylum, exclude=NULL)

phylum.col <- NULL
#            Actinobacteriota: 35               
phylum.col <- c(phylum.col, colorRampPalette(c(colors()[652],colors()[502]))(35))
#                Bacteroidota: 24
phylum.col <- c(phylum.col, colorRampPalette(c(colors()[86],colors()[85]))(24))
#             Campilobacterota: 2
phylum.col <- c(phylum.col, colors()[c(620,624)])
#              Cyanobacteria: 3
phylum.col <- c(phylum.col, colorRampPalette(c("white",colors()[457]))(4)[-1])
#                 Deinococcota: 1
phylum.col <- c(phylum.col, colors()[427])
#                 Desulfobacterota: 1
phylum.col <- c(phylum.col, colors()[413])
#                   Firmicutes: 34
phylum.col <- c(phylum.col, colorRampPalette(c(colors()[60],colors()[552],colors()[137]))(34))
#                Fusobacteria: 12
phylum.col <- c(phylum.col, colorRampPalette(c("white",colors()[638]))(7)[-1])
#              Patescibacteria: 2
phylum.col <- c(phylum.col, colorRampPalette(c("white",colors()[457]))(3)[-1])
#               Proteobacteria: 107
phylum.col <- c(phylum.col, colorRampPalette(c(colors()[601],colors()[26],colors()[73]))(107))
#            no Phylum: 1
phylum.col <- c(phylum.col, grey(0.5))
# other = "palegoldenrod"
phylum.col[length(phylum.col)] <- "palegoldenrod"

dim(prop.other)
length(phylum.col)

barplot (prop.other, col=phylum.col, las=2, cex.names=0.5, border=NA)
col.legend (phylum.col, rownames(prop.other))

mean.ra <- apply(prop.other[-nrow(prop.other),],1,mean)
top.ASVs <- names(rev(sort(mean.ra)))[1:25]
col.vec2 <- colorRampPalette(c(grey(0.9),grey(0.1)))(nrow(prop.other)-1)
col.vec2 <- c(col.vec2, "palegoldenrod")
col.vec2[match(top.ASVs,rownames(prop.other))] <- top25.col

barplot (prop.other, col=col.vec2, las=2, cex.names=0.5, border=NA)
col.legend (col.vec2, rownames(prop.other))
taxdat0[match(top.ASVs,rownames(taxdat0)),]

# --- very low copy number vs negative controls

dim(seq.samples)
dim(seq.negcontrols)
sample.copies <- labdata2$templateDNA_16SqPCR_copies[is.na(labdata2$nontemplate)]
seq.low.copy <- seq.samples[,sample.copies<=100]
seq.low.copy <- seq.low.copy[,rev(order(sample.copies[sample.copies<=100]))]
dim(seq.low.copy)
seq.comb <- cbind (seq.low.copy, seq.negcontrols)
seq.comb <- seq.comb[apply(seq.comb,1,sum)>0,]
dim(seq.comb)
rev(sort(apply(seq.comb,2,sum)))
seq.comb <- seq.comb[,apply(seq.comb,2,sum)>0]
prop.comb <- apply(seq.comb, 2, function(x)x/sum(x))
dim(prop.comb)

max.ra <- apply(prop.comb,1,max)
prop.other <- prop.comb[max.ra>=0.01,]
other <- apply(prop.comb[max.ra<0.01,],2,sum)
prop.other <- rbind (prop.other, other)
dim(prop.other)

rm("seq.low.copy")

tax.low.copy <- taxdat0[match(rownames(prop.other),rownames(taxdat0)),]
head(tax.low.copy)
table(tax.low.copy$Phylum, exclude=NULL)

phylum.col <- NULL
#            Actinobacteriota: 32               
phylum.col <- c(phylum.col, colorRampPalette(c(colors()[652],colors()[502]))(32))
#                Bacteroidota: 15
phylum.col <- c(phylum.col, colorRampPalette(c(colors()[86],colors()[85]))(15))
#             Campilobacterota: 1                 
phylum.col <- c(phylum.col, colorRampPalette(c("white",colors()[614]))(3)[-1])
#                Cyanobacteria: 3
phylum.col <- c(phylum.col, colorRampPalette(c("white",colors()[638]))(4)[-1])
#                 Deinococcota: 1
phylum.col <- c(phylum.col, colors()[427])
#                 Desulfobacterota: 1
phylum.col <- c(phylum.col, colors()[413])
#              Firmicutes: 25
phylum.col <- c(phylum.col, colorRampPalette(c("white",colors()[457]))(26)[-1])
#             Fusobacteria: 3
phylum.col <- c(phylum.col, colors()[c(620,622,624)])
#               Proteobacteria: 102
phylum.col <- c(phylum.col, colorRampPalette(c(colors()[601],colors()[26],colors()[73]))(102))
#            no Phylum: 1
phylum.col <- c(phylum.col, grey(0.5))
# other = "palegoldenrod"
phylum.col[length(phylum.col)] <- "palegoldenrod"

dim(prop.other)
length(phylum.col)

barplot (prop.other, col=phylum.col, las=2, cex.names=0.5, border=NA)
col.legend (phylum.col, rownames(prop.other))

mean.ra <- apply(prop.other[-nrow(prop.other),],1,mean)
top.ASVs <- names(rev(sort(mean.ra)))[1:25]
col.vec2 <- colorRampPalette(c(grey(0.9),grey(0.1)))(nrow(prop.other)-1)
col.vec2 <- c(col.vec2, "palegoldenrod")
col.vec2[match(top.ASVs,rownames(prop.other))] <- top25.col

barplot (prop.other, col=col.vec2, las=2, cex.names=0.5, border=NA)
col.legend (col.vec2, rownames(prop.other))
taxdat0[match(top.ASVs,rownames(taxdat0)),]

# --- very low  read counts vs negative controls

dim(seq.samples)
dim(seq.negcontrols)
sample.reads <- labdata2$FinalReads_ASVpipeline[is.na(labdata2$nontemplate)]
seq.low.reads <- seq.samples[,sample.reads<=100]
seq.low.reads <- seq.low.reads[,rev(order(sample.reads[sample.reads<=1000]))]
dim(seq.low.reads)
seq.comb <- cbind (seq.low.reads, seq.negcontrols)
seq.comb <- seq.comb[apply(seq.comb,1,sum)>0,]
dim(seq.comb)
rev(sort(apply(seq.comb,2,sum)))
seq.comb <- seq.comb[,apply(seq.comb,2,sum)>0]
prop.comb <- apply(seq.comb, 2, function(x)x/sum(x))
dim(prop.comb)

max.ra <- apply(prop.comb,1,max)
prop.other <- prop.comb[max.ra>=0.01,]
other <- apply(prop.comb[max.ra<0.01,],2,sum)
prop.other <- rbind (prop.other, other)
dim(prop.other)

rm("seq.low.reads")

tax.low.reads <- taxdat0[match(rownames(prop.other),rownames(taxdat0)),]
head(tax.low.reads)
table(tax.low.reads$Phylum, exclude=NULL)

phylum.col <- NULL
#            Actinobacteriota: 32               
phylum.col <- c(phylum.col, colorRampPalette(c(colors()[652],colors()[502]))(32))
#                Bacteroidota: 13
phylum.col <- c(phylum.col, colorRampPalette(c(colors()[86],colors()[85]))(13))
#             Campilobacterota: 1                 
phylum.col <- c(phylum.col, colors()[96])
#             Cyanobacteria: 2                 
phylum.col <- c(phylum.col, colorRampPalette(c("white",colors()[614]))(3)[-1])
#                 Deinococcota: 1
phylum.col <- c(phylum.col, colors()[427])
#                 Desulfobacterota: 1
phylum.col <- c(phylum.col, colors()[413])
#              Firmicutes: 22
phylum.col <- c(phylum.col, colorRampPalette(c("white",colors()[457]))(23)[-1])
#             Fusobacteria: 3
phylum.col <- c(phylum.col, colors()[c(620,622,624)])
#               Proteobacteria: 97
phylum.col <- c(phylum.col, colorRampPalette(c(colors()[601],colors()[26],colors()[73]))(97))
#            no Phylum: 1
phylum.col <- c(phylum.col, grey(0.5))
# other = "palegoldenrod"
phylum.col[length(phylum.col)] <- "palegoldenrod"

dim(prop.other)
length(phylum.col)

barplot (prop.other, col=phylum.col, las=2, cex.names=0.5, border=NA)
col.legend (phylum.col, rownames(prop.other))

mean.ra <- apply(prop.other[-nrow(prop.other),],1,mean)
top.ASVs <- names(rev(sort(mean.ra)))[1:25]
col.vec2 <- colorRampPalette(c(grey(0.9),grey(0.1)))(nrow(prop.other)-1)
col.vec2 <- c(col.vec2, "palegoldenrod")
col.vec2[match(top.ASVs,rownames(prop.other))] <- top25.col

barplot (prop.other, col=col.vec2, las=2, cex.names=0.5, border=NA)
col.legend (col.vec2, rownames(prop.other))
taxdat0[match(top.ASVs,rownames(taxdat0)),]


#=============================================================================================================
# --- Decontam
#=============================================================================================================

dim(seq.samples) # 3219 x 960
dim(seq.negcontrols) # 3219 x 43
dim(labdata2) # 1003   60

#  We use 43 negative controls + 960 biological samples 
#      (remove samples with <=100 copies and <=1000 reads)

seq.comb <- cbind (seq.samples, seq.negcontrols)
dim(seq.comb) # 3219 x 1003
seq.comb <- seq.comb[,!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$age_years_numeric>1)]
labdata2 <- labdata2[!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$age_years_numeric>1),]
dim(seq.comb) # 3219 x 1003 (no samples removed based on age)
seq.comb <- seq.comb[,!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$templateDNA_16SqPCR_copies>100)]
labdata2 <- labdata2[!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$templateDNA_16SqPCR_copies>100),]
dim(seq.comb) # 3219 x 1001 (removed 2 low copy samples)
seq.comb <- seq.comb[,!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$FinalReads_ASVpipeline>1000)]
labdata2 <- labdata2[!is.na(labdata2$nontemplate) | (is.na(labdata2$nontemplate) & labdata2$FinalReads_ASVpipeline>1000),]
dim(labdata2) # 996 x 60 (removed 5 low read number samples)

seq.comb <- seq.comb[apply(seq.comb,1,sum)>0,]
dim(seq.comb) # 3216 x 996 (remove 3 ASVs with zero reads for these 996 samples)
seq.samples <- seq.comb[,is.na(labdata2$nontemplate)]
dim(seq.samples) # 3216 x 953 biological samples
# --- decontam
require(decontam)

out <- isContaminant (t(seq.comb), conc=labdata2$templateDNA_16SqPCR_copies,
                      neg=!is.na(labdata2$nontemplate), threshold=0.40)
length(out$contaminant)  #3216
sum(out$contaminant)     #70 potential contaminant ASVs
#write.csv(out,"ContaminantsR.csv")
#write.csv(out, paste(pathname, "Contaminants.csv", sep=""))

prop.table <- apply(seq.comb, 2, function(x) x/sum(x))
neg <- !is.na(labdata2$nontemplate)
mat <- apply (prop.table, 1, function(x)tapply(x,neg,function(y)sum(y>0)/length(y)))

plot (c(0,max(mat)),c(0,max(mat)),type="n",asp=1,xlab="prevalence in neg controls",ylab="prevalence in biological")
text (mat[2,!out$contaminant],mat[1,!out$contaminant],colnames(mat)[!out$contaminant],col="olivedrab",cex=0.5)
text (mat[2,out$contaminant],mat[1,out$contaminant],colnames(mat)[out$contaminant],col="red",cex=0.5)

contaminants <- rownames(out)[out$contaminant]
for (i in 1:length(contaminants))
  seq.samples <- seq.samples[-match(contaminants[i],rownames(seq.samples)),]
dim(seq.samples)         # 3146 953
seq.samples <- seq.samples[apply(seq.samples,1,sum)>0,]
dim(seq.samples)         # 2829 953

# --- remove spurious ASVs (<=10 reads across the dataset without problem samples)
num.counts <- apply(seq.samples,1,sum)
prop.samples <- apply(seq.samples, 2, function(x)x/sum(x))
max.per.sample <- apply(prop.samples,1,max)
table(num.counts)[1:20]

for (i in 2:20)
{
  cat("========== total count = ", i, "==========\n") 
  sub.vec <- max.per.sample[num.counts==i]
  dev.new()
  hist(sub.vec*100, xlab="%", main=paste("Total count =",i))
  num.reads <- num.counts[num.counts<=i]
  ASVs <- rownames(seq.samples)[num.counts<=i]
  max.props <- max.per.sample[num.counts<=i]
  detected <- apply(seq.samples[num.counts<=i,,drop=F],1,function(x)sum(x>0))
  dat <- data.frame(num.reads,ASVs,max.props,detected)
  write.csv (dat, file=paste(pathname,"spurious",i,".csv",sep=""))
}

seq.samples2 <- seq.samples[apply(seq.samples,1,sum)>10,]
dim(seq.samples2)         # 1668 953

taxdat0 <- taxdat0[match(rownames(seq.samples2),rownames(taxdat0)),]
dim(taxdat0)       # 1668 ASVs with 7 columns (7 taxon levels)
labdata2 <- labdata2[is.na(labdata2$nontemplate),]
dim(labdata2)     # only 953 samples, remove negative controls

table(taxdat0$Phylum, exclude=NULL)   # 15 ASVs with no Phylum
prop.samples <- apply(seq.samples2, 2, function(x)x/sum(x))

taxdat0[,2] <- as.character(taxdat0[,2])
taxdat0[is.na(taxdat0[,2]),2] <- rownames(taxdat0)[is.na(taxdat0[,2])]
taxdat0[,3] <- as.character(taxdat0[,3])
taxdat0[is.na(taxdat0[,3]),3] <- rownames(taxdat0)[is.na(taxdat0[,3])]
taxdat0[,4] <- as.character(taxdat0[,4])
taxdat0[is.na(taxdat0[,4]),4] <- rownames(taxdat0)[is.na(taxdat0[,4])]
taxdat0[,5] <- as.character(taxdat0[,5])
taxdat0[is.na(taxdat0[,5]),5] <- rownames(taxdat0)[is.na(taxdat0[,5])]
taxdat0[,6] <- as.character(taxdat0[,6])
taxdat0[is.na(taxdat0[,6]),6] <- rownames(taxdat0)[is.na(taxdat0[,6])]

seq.imputed <- seq.samples2
set.seed (730317)
seq.imputed[seq.imputed==0] <- runif(sum(seq.samples2==0),0.1,1)

#write.csv (taxdat0, paste(pathname,"taxdat953.csv",sep=""))
#write.csv (seq.samples2, paste(pathname,"ASVcount953.csv",sep=""))
#write.csv (prop.samples, paste(pathname,"ASVprop953.csv",sep=""))
#write.csv (labdata2, paste(pathname,"labdata953.csv",sep=""))
#write.csv (seq.imputed, paste(pathname,"ASVcount953_imputed.csv",sep=""))

#5 samples were excluded  from ASVcount953.csv to produce ASVcount948.csv
#1. SP.BSP0066.BREATHe.R01.P4.E01_S337_L001_R
#2. SP.B1001.trial.BSP0227.18m.BREATHe.R01.P4.H11_S383_L001_R
#3. SP.C305.comparison.group.BYF104T3.baseline.Malawi.BREATHe.R03.P4.H02_S374_L001_R
#4. SP.C309.comparison.group.BYF108T3.baseline.Malawi.BREATHe.R03.P4.H03_S375_L001_R
#5. SP.C316.comparison.group.BYF10GT3.baseline.Malawi.BREATHe.R03.P4.H01_S373_L001_R

#The removal of these five samples resulted in the exclusion of some ASVs that occur only in these 5
#Hence these were removed from  taxdat953.csv to produce taxdat948.csv used for downstream analyses


