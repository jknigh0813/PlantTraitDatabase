#####################################################################
# Reformat RF residuals for GraPhlan plotting
# Created: James Knighton (7/25/2023)
#####################################################################

#library(phytools)
library(MPSEM)
library(WorldFlora)
library(V.PhyloMaker2)
library(dplyr)
library(ape)

inpath = "C:/PhyloTraitEst/RF_Model/"
outpath = "C:/PhyloTraitEst/Graphlan/"
TraitNames = c("gsmax","P12","P50","P88","rdmax","WUE","height","SLA","LeafN")

for (j in 1:7)
{
  GraphLan = data.frame("Species")
  colnames(GraphLan) <- c("Species")
  
  TraitName = TraitNames[j]
  
  #1.Read in Newick Tree and TRY data
  Traits <- read.table(paste(inpath,"RF_Test_Predict_",TraitName,".csv",sep=""),sep=",",header=FALSE)
  colnames(Traits) <- c("Species","Sim","Pred")
  Traits$Res = Traits$Sim - Traits$Pred
  specieslist = unique(Traits$Species)
  
  Out = data.frame(matrix(ncol=2,nrow=length(unique(specieslist))))
  
  #2. Median species residual values
  for (i in 1:length(specieslist))
  {
    speciesname = specieslist[i]
    Records = Traits[Traits$Species == speciesname,]
    Out[i,1] = speciesname
    Out[i,2] = median(Records$Res)
    
  }
  colnames(Out) <- c("Species",'Value')
  
  GraphLan <- merge(GraphLan, Out, by = c("Species"), all.x = TRUE, all.y = TRUE)
  GraphLan = GraphLan[which(GraphLan$Species != "Species"),]
  colnames(GraphLan) <- c("spec.name","Residual")  
  
  #3. Match Names
  NameMatch = WFO.match(spec.data=GraphLan, WFO.file = "C:/PhyloTraitEst/WFO/classification.csv", counter=1, verbose=TRUE, Fuzzy=0)
  NameMatch = NameMatch[which(NameMatch$Subseq == 1),]
  NameMatch = NameMatch[which(NameMatch$taxonRank == "species"),]
  NameMatch = NameMatch[which(NameMatch$scientificName != "Quercus × hispanica"),]
  
  GraphLanOut = data.frame(NameMatch$scientificName,NameMatch$genus,NameMatch$family,NameMatch$Residual)
  colnames(GraphLanOut) <- c("species","genus","family","Residual")
  
  #4. Run V.Phylomaker2 scenario S3
  phylotree <- phylo.maker(GraphLanOut, scenarios=c("S3"))
  
  GraphLanOut$MatchString = gsub(" ", "_", GraphLanOut$species)
  GraphLanOut = GraphLanOut[GraphLanOut$MatchString %in% phylotree$scenario.3$tip.label,]
  GraphLanOut$Annotate = 0.1 + 0.9*abs(GraphLanOut$Residual)/max(abs(GraphLanOut$Residual))
  
  #5. Write GraPhlan Newick Tree
  write.tree(phylotree$scenario.3,paste(outpath,"RF_",TraitName,"_Newick.txt",sep=""))
  
  #5. Write GraPhlan Annotation File
  write.table(GraphLanOut,paste(outpath,"RF_",TraitName,"_Annotation.csv",sep=""),row.names = FALSE,sep=",")
  
  #Clean up if memory is getting stressed
  rm(Traits)
  rm(NameMatch)
  rm(phylotree)
  rm(GraphLan)
  rm(GraphLanOut)
  gc()
}





