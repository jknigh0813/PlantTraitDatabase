#####################################################################
# Validate BCGI tree names against WFO
#Created: James Knighton (7/25/2023)
#####################################################################

library(rtry)
library(WorldFlora)
library(V.PhyloMaker2)
library(ape)
library(dplyr)

inpath = "C:/PhyloTraitEst/"

#1. Import BCGI tree names (downloaded 7/25/2023)
Trees <- read.csv(paste(inpath,"TRY_RawData/global_tree_search_trees_1_7.csv",sep=""))
keep = grepl('^\\w+\\s\\w+$', Trees$TaxonName)
Trees = as.data.frame(Trees[keep,])
colnames(Trees) <- c('TaxonName')
Trees[,2:3] <- data.frame(do.call('rbind', strsplit(as.character(Trees$TaxonName),' ',fixed=TRUE)))
colnames(Trees) <- c('spec.name','Genus','Species')
Trees = Trees[which(Trees$Species != "sp"),]

#2. Match names against WFO
NameMatch = WFO.match(spec.data=Trees, WFO.file = "C:/PhyloTraitEst/WFO/classification.csv", counter=1, verbose=TRUE,Fuzzy=0)
NameMatch = NameMatch[which(NameMatch$taxonomicStatus=="Accepted"),]

#3. Write corrected names
TraitExport = data.frame(NameMatch$scientificName,NameMatch$family,NameMatch$Genus,NameMatch$specificEpithet)
colnames(TraitExport) <- c('spec.name','Family','Genus','Species')
write.table(TraitExport,paste(inpath,"TRY_Values_NewickTrees/GlobalTrees_WFONames.csv",sep=""),row.names = FALSE,sep=",",quote = FALSE)
