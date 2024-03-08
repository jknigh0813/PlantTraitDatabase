#####################################################################
# Solve PEMs from TRY Data
# keeps the most rank-correlated columns for RF training
# Note: this step is computationally expensive. This code is best adapted for parallel execution on HPC
# Created: James Knighton (7/25/2023)
#####################################################################

library(phytools)
library(MPSEM)

inpath = "C:/PhyloTraitEst/"
outpath = "C:/PhyloTraitEst/RF_Model/"
TraitName = "gsmax"

#1.Read in Newick Tree and TRY data
phylotree <- read.newick(paste(inpath,"TRY_Values_NewickTrees/TRY_",TraitName,"_Newick.txt",sep=""))
Traits = read.csv(paste(inpath,"TRY_Values_NewickTrees/TRY_",TraitName,"_Values.csv",sep=""))

#2. Build PEM graph
Traits.pgraph <- Phylo2DirectedGraph(phylotree)
#Traits.pgraph$edge$distance[Traits.pgraph$edge$distance > 200] = 200
Eigens = as.data.frame(PEM.build(Traits.pgraph))
Eigens$spec.name = gsub("_", " ", row.names(Eigens))

#3. Merge PEMs with TRY trait data
Trait_Eigen <- merge(Traits,Eigens,by="spec.name",all=T)
Trait_Eigen = Trait_Eigen[!is.na(Trait_Eigen$V_1),]

#4. Find n-most rank-correlated PEMs, discard the rest
result = matrix(nrow=ncol(Trait_Eigen)-6,ncol=2)
for (j in 7:(ncol(Trait_Eigen)))
{
  result[j-6,1] = j
  result[j-6,2] = cor(Trait_Eigen$Value, Trait_Eigen[,j], method = "spearman",use="pairwise.complete.obs")
  print(j)
}
result[,2] = abs(result[,2])
result <- result[order(-result[,2]),]

#5. Export PEMs and observations for RF trait prediction
Train = data.frame(Trait_Eigen$spec.name, Trait_Eigen$Value, Trait_Eigen[,result[1:50,1]])
write.table(Train,paste(outpath,"RF_Train_",TraitName,".csv",sep=""),row.names = FALSE,sep=",")
