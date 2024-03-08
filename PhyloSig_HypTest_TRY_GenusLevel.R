#####################################################################
# Pagel's lambda test | Genus level analysis
# Created: James Knighton (7/25/2023)
#####################################################################

library(phytools)
library(Taxonstand)
library(diversitree)
library(ape)
library(caper)

inpath = "C:/PhyloTraitEst/TRY_Values_NewickTrees/"
outpath = "C:/PhyloTraitEst/TRY_Values_NewickTrees/"
TraitName = "P50"

#1.Read in Newick Tree and TRY data
phylotree <- read.newick(paste(inpath,"TRY_",TraitName,"_Newick_GenusLevel.txt",sep=""))
Traits = read.csv(paste(inpath,"TRY_",TraitName,"_Values_GenusLevel.csv",sep=""))
Traits$MatchString = gsub(" ", "_", Traits$spec.name)
                      
#2. Median species trait value from TRY database
for (i in 1:length(phylotree$tip.label))
{
  speciesname = phylotree$tip.label[i]
  Records = Traits[Traits$MatchString == speciesname,]
  phylotree$trait[i] = median(Records$Value)
}

#3. Pagel's lambda tests for phylogenetic signal

L1 = phylosig(phylotree, phylotree$trait, method="lambda",test=TRUE,niter=10, se=NULL, start=NULL,
              control=list())

#K1 = phylosig(phylotree, phylotree$trait, method="K",test=TRUE,nsim=1000, se=NULL, start=NULL,
#              control=list())
