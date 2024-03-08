#####################################################################
# Pagel's lambda test for TRY observations
# Created: James Knighton (7/25/2023)
#####################################################################

library(phytools)
library(Taxonstand)
library(diversitree)
library(ape)
library(caper)

inpath = "C:/PhyloTraitEst/"
outpath = "C:/PhyloTraitEst/RF_Model/"
TraitName = "LeafN"

#1.Read in Newick Tree and TRY data
phylotree <- read.newick(paste(inpath,"TRY_Values_NewickTrees/TRY_",TraitName,"_Newick.txt",sep=""))
Traits = read.csv(paste(inpath,"TRY_Values_NewickTrees/TRY_",TraitName,"_Values.csv",sep=""))
Traits$MatchString = gsub(" ", "_", Traits$spec.name)
                      
#2. Compute median species trait value from TRY database
for (i in 1:length(phylotree$tip.label))
{
  speciesname = phylotree$tip.label[i]
  Records = Traits[Traits$MatchString == speciesname,]
  phylotree$trait[i] = median(Records$Value)
}

#3. Pagel's lambda test for phylogenetic signal
L1 = phylosig(phylotree, phylotree$trait, method="lambda",test=TRUE,niter=100, se=NULL, start=NULL,
              control=list())

#K1 = phylosig(phylotree, phylotree$trait, method="K",test=TRUE,nsim=1000, se=NULL, start=NULL,
#              control=list())


