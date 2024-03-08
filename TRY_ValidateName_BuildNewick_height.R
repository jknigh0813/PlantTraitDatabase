#####################################################################
# Reformat height TRY data for further analysis
# Created: James Knighton (7/25/2023)
#####################################################################

library(rtry)
library(WorldFlora)
library(V.PhyloMaker2)
library(ape)
library(dplyr)

inpath = "C:/PhyloTraitEst/"

#0. Import TRY data (downloaded 7/25/2023)
TRYdata1 <- rtry_import(paste(inpath,"TRY_RawData/28262.txt",sep=""))

#1. Extract TRY data and clean
Trait = data.frame(TRYdata1$ObsDataID,TRYdata1$AccSpeciesName, TRYdata1$StdValue, TRYdata1$DataName, TRYdata1$TraitID)
colnames(Trait) <- c('DataID','FullName','Value','DataName','TraitID')
Trait = Trait[which(Trait$TraitID == 3106),]
datatypes = unique(Trait$DataName)
traittypes = unique(Trait$TraitID)
keeps = c("Plant height maximum","Maximum height","Maximum height max","Plant height observed max","Plant height vegetative max.")
Trait = Trait[which(Trait$DataName %in% keeps),]
Trait = Trait[which(!is.na(Trait$Value)),]
Trait = Trait[which(Trait$Value > 1),] #keep only trees > 1 m
Trait = Trait[which(Trait$FullName != "unknown"),]
keep = grepl('^\\w+\\s\\w+$', Trait$FullName)
Trait = Trait[keep,]
Trait[,6:7] <- data.frame(do.call('rbind', strsplit(as.character(Trait$FullName),' ',fixed=TRUE)))
colnames(Trait) <- c('DataID','spec.name','Value','DataName','TraitID','Genus','Species')
Trait = Trait[which(Trait$Species != "sp"),]

#remove plants with very high variance observations (i.e., potentially mislabeled observations?)
temp = Trait %>% group_by(spec.name) %>% summarise(SD = sd(Value,na.rm=TRUE), Mean = mean(Value,na.rm=TRUE), N = length(Value))
temp$cv = temp$SD/temp$Mean
highvar = unique(temp$spec.name[which(temp$cv > 0.5)])
Trait = Trait[(which(!Trait$spec.name %in% highvar)),]
trees = read.csv("C:/PhyloTraitEst/TRY_Values_NewickTrees/GlobalTrees_WFONames.csv")
Trait = Trait[which(Trait$spec.name %in% trees$spec.name),]

#2. Match names against WFO
NameMatch = WFO.match(spec.data=Trait, WFO.file = "C:/PhyloTraitEst/WFO/classification.csv", counter=1, verbose=TRUE, Fuzzy=0)
NameMatch = NameMatch[which(NameMatch$taxonomicStatus=="Accepted"),]
NameMatch = NameMatch[!duplicated(NameMatch$DataID),]

#3. Run V.Phylomaker2 scenario S3
WFO_Names = data.frame(species = NameMatch$scientificName, genus = NameMatch$genus, family = NameMatch$family)
phylotree <- phylo.maker(WFO_Names, scenarios=c("S3"))

#4. Write Newick Tree
write.tree(phylotree$scenario.3,paste(inpath,"TRY_Values_NewickTrees/TRY_height_Newick.txt",sep=""))

#5. Write Raw Data
TraitExport = unique(data.frame(NameMatch$DataID,NameMatch$scientificName,NameMatch$family,NameMatch$Genus,NameMatch$specificEpithet,NameMatch$Value))
colnames(TraitExport) <- c('DataID','spec.name','Family','Genus','Species','Value')
write.table(TraitExport,paste(inpath,"TRY_Values_NewickTrees/TRY_height_Values.csv",sep=""),row.names = FALSE,sep=",",quote = FALSE)
