#####################################################################
# The following code computes phylogenetic eigenvectors from TRY data
#Created: James Knighton (7/25/2023)
#####################################################################

library(rtry)
library(WorldFlora)
library(V.PhyloMaker2)
library(ape)
library(dplyr)
traitname = "LeafN"

inpath = "C:/PhyloTraitEst/TRY_Values_NewickTrees/"

#Import TRY data (downloaded 7/25/2023)
TRYdata1 <- read.csv(paste(inpath,"TRY_",traitname,"_Values.csv",sep=""))
genera = unique(TRYdata1$Genus)
for (g in 1:length(genera))
{
    subset = TRYdata1[TRYdata1$Genus == genera[g],]
    specname = subset$Species[1]
    TRYdata1$Species[TRYdata1$Genus == genera[g]] = specname
    TRYdata1$spec.name[TRYdata1$Genus == genera[g]] = paste(genera[g]," ",specname,sep="")
    
}

#Run V.Phylomaker2 scenario S3
WFO_Names = data.frame(species = TRYdata1$spec.name, genus = TRYdata1$Genus, family = TRYdata1$Family)
phylotree <- phylo.maker(WFO_Names, scenarios=c("S3"))

#Write Newick Tree
write.tree(phylotree$scenario.3,paste(inpath,"/TRY_",traitname,"_Newick_GenusLevel.txt",sep=""))

#Write Raw Data
TraitExport = unique(data.frame(TRYdata1$DataID,TRYdata1$spec.name,TRYdata1$Family,TRYdata1$Genus,TRYdata1$Species,TRYdata1$Value))
colnames(TraitExport) <- c('DataID','spec.name','Family','Genus','Species','Value')
write.table(TraitExport,paste(inpath,"/TRY_",traitname,"_Values_GenusLevel.csv",sep=""),row.names = FALSE,sep=",",quote = FALSE)
