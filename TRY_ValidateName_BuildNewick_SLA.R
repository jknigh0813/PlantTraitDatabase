#####################################################################
# Reformat SLA TRY data for further analysis
# Created: James Knighton (7/25/2023)
#####################################################################

library(rtry)
library(WorldFlora)
library(V.PhyloMaker2)
library(ape)

inpath = "D:/PhyloTraitEst/"

#0. Import TRY data (downloaded 7/25/2023)
TRYdata1 <- rtry_import(paste(inpath,"TRY_RawData/32298.txt",sep=""))

#1. Extract TRY data and clean
Trait = data.frame(TRYdata1$ObsDataID,TRYdata1$AccSpeciesName, TRYdata1$StdValue, TRYdata1$DataName, TRYdata1$TraitID, TRYdata1$ErrorRisk)
colnames(Trait) <- c('DataID','FullName','Value','DataName','TraitID','ErrorRisk')
Trait = Trait[which(Trait$TraitID == 3117),]
Trait = Trait[which(Trait$ErrorRisk <= 5),]
datatypes = unique(Trait$DataName)
#keeps = c("SLA: undefined if petiole in- or excluded","SLA: undefined if petiole in- or excluded (1)")
keeps = c("SLA: undefined if petiole in- or excluded")
Trait = Trait[which(Trait$DataName %in% keeps),]
Trait = Trait[which(!is.na(Trait$Value)),]
Trait = Trait[which(Trait$Value > 0),]
Trait = Trait[which(Trait$FullName != "unknown"),]
Trait = Trait[which(Trait$FullName != "Coffea arabica" & Trait$FullName != "Glycine max"),]
keep = grepl('^\\w+\\s\\w+$', Trait$FullName)
Trait = Trait[keep,]
Trait[,6:7] <- data.frame(do.call('rbind', strsplit(as.character(Trait$FullName),' ',fixed=TRUE)))
colnames(Trait) <- c('DataID','spec.name','Value','DataName','TraitID','Genus','Species')
Trait = Trait[which(Trait$Species != "sp"),]

#2. Match names against WFO
NameMatch = WFO.match(spec.data=Trait, WFO.file = "C:/PhyloTraitEst/WFO/classification.csv", no.dates = TRUE, counter=1, verbose=TRUE,Fuzzy=0)
NameMatch = NameMatch[which(NameMatch$taxonomicStatus=="Accepted"),]
NameMatch = NameMatch[!duplicated(NameMatch$DataID),]

#3. Run V.Phylomaker2 scenario S3
WFO_Names = data.frame(species = NameMatch$scientificName, genus = NameMatch$Genus, family = NameMatch$family)
phylotree <- phylo.maker(WFO_Names, scenarios=c("S3"))

#4. Write Newick Tree
write.tree(phylotree$scenario.3,paste(inpath,"TRY_Values_NewickTrees/TRY_SLA_Newick.txt",sep=""))

#5. Write Raw Data
TraitExport = unique(data.frame(NameMatch$DataID,NameMatch$scientificName,NameMatch$family,NameMatch$genus,NameMatch$specificEpithet,NameMatch$Value))
colnames(TraitExport) <- c('DataID','spec.name','Family','Genus','Species','Value')
write.table(TraitExport,paste(inpath,"TRY_Values_NewickTrees/TRY_SLA_Values.csv",sep=""),row.names = FALSE,sep=",",quote = FALSE)
