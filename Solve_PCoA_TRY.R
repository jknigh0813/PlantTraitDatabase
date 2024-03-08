#####################################################################
# Solve PCoAs from TRY Data
# keeps the 50 highest by variance explained
# Note: this step is computationally expensive. This code is best adapted for parallel execution on HPC
# Created: James Knighton (7/25/2023)
#####################################################################

library(ape)
library(phytools)
library(MPSEM)

inpath = "C:/PhyloTraitEst/"
outpath = "C:/PhyloTraitEst/RF_Model/"
TraitName = "LeafN"

#1.Read in Newick Tree and TRY data
phylotree <- read.newick(paste(inpath,"TRY_Values_NewickTrees/TRY_",TraitName,"_Newick.txt",sep=""))
Traits = read.csv(paste(inpath,"TRY_Values_NewickTrees/TRY_",TraitName,"_Values.csv",sep=""))

#2. Solve distance matrix
distmat = cophenetic.phylo(phylotree)
  
#3. Solve PCoAs
PCOAS <- pcoa(distmat)
PCA_Vals = as.data.frame(PCOAS$vectors)
PCA_Vals$spec.name = gsub("_", " ", row.names(PCA_Vals))

#4. Export eigenvectors for all species for RF trait prediction
Trait_PCOAS <- merge(Traits,PCA_Vals,by="spec.name",all=T)
Train = data.frame(Trait_PCOAS$spec.name, Trait_PCOAS$Value, Trait_PCOAS[,7:57])
write.table(Train,paste(outpath,"RF_Train_",TraitName,"_PCOA.csv",sep=""),row.names = FALSE,sep=",")
