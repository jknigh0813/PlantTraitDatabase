#####################################################################
# Reformat TRY observations for GraPhlan plotting
#Created: James Knighton (7/25/2023)
#####################################################################

#library(phytools)
library(MPSEM)
library(WorldFlora)
library(V.PhyloMaker2)
library(dplyr)
library(ape)

inpath = "C:/PhyloTraitEst/"
outpath = "C:/PhyloTraitEst/Graphlan/"
#TraitNames = c("gsmax","P12","P50","P88","rdmax","WUE","height","SLA","LeafN")
TraitNames = c("height")
GraphLan = data.frame("Species","Genus","Family")
colnames(GraphLan) <- c("Species","Genus","Family")

for (j in 1:length(TraitNames))
{
  TraitName = TraitNames[j]
  
  #1.Read in Newick Tree and TRY data
  Traits = read.csv(paste(inpath,"TRY_Values_NewickTrees/TRY_",TraitName,"_Values.csv",sep=""))
  Traits$MatchString = paste(Traits$Genus,Traits$Species)
  specieslist = unique(Traits$MatchString)
    
  Out = data.frame(matrix(ncol=4,nrow=length(unique(specieslist))))
  
  #2. Compute median species trait value from TRY database
  for (i in 1:length(specieslist))
  {
    speciesname = specieslist[i]
    Records = Traits[Traits$MatchString == speciesname,]
    Out[i,1] = unique(Traits$Species[Traits$MatchString == speciesname])
    Out[i,2] = unique(Traits$Genus[Traits$MatchString == speciesname])
    Out[i,3] = unique(Traits$Family[which(Traits$MatchString == speciesname)])
    Out[i,4] = median(Records$Value)
    
  }
  colnames(Out) <- c("Species","Genus","Family",'Value')
  
  GraphLan <- merge(GraphLan, Out, by = c("Species","Genus","Family"), all.x = TRUE, all.y = TRUE)
  
}

#colnames(GraphLan) <- c("species","genus","family","gsmax","P12","P50","P88","rdmax","WUE")
colnames(GraphLan) <- c("species","genus","family","height")
GraphLan$species = paste(GraphLan$genus,GraphLan$species)

#3. Run V.Phylomaker2 scenario S3
phylotree <- phylo.maker(GraphLan, scenarios=c("S3"))
GraphLan$MatchString = gsub(" ", "_", GraphLan$species)
GraphLanOut = GraphLan[GraphLan$MatchString %in% phylotree$scenario.3$tip.label,]

#4. Write GraPhlan Newick Tree
write.tree(phylotree$scenario.3,paste(outpath,"TRY_height_Newick.txt",sep=""))

#5. Write GraPhlan Annotation File
write.table(GraphLanOut,paste(outpath,"TRY_height_Annotation.csv",sep=""),row.names = FALSE,sep=",")

