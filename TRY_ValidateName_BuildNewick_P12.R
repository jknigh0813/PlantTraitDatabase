#####################################################################
# Reformat P12 TRY data for further analysis
#Created: James Knighton (7/25/2023)
#####################################################################

library(rtry)
library(WorldFlora)
library(V.PhyloMaker2)
library(ape)

inpath = "C:/PhyloTraitEst/"

#0. Import TRY data (downloaded 7/25/2023)
TRYdata1 <- rtry_import(paste(inpath,"TRY_RawData/27968.txt",sep=""))

#1. Extract TRY data and clean
Trait = data.frame(TRYdata1$ObsDataID,TRYdata1$AccSpeciesName, TRYdata1$StdValue, TRYdata1$DataName, TRYdata1$TraitID, TRYdata1$ErrorRisk)
colnames(Trait) <- c('DataID','FullName','Value','DataName','TraitID','ErrorRisk')
trait_vals = unique(Trait$TraitID)
Trait = Trait[which(Trait$ErroRisk < 7),]
Trait = Trait[which(Trait$TraitID == 719),]
Trait = Trait[which(Trait$Value <= -0.25),]
datatypes = unique(Trait$DataName)
keeps = c("Xylem water potential at which 12% of conductivity is lost (P12)")
Trait = Trait[which(Trait$DataName %in% keeps),]

Trait = Trait[which(!is.na(Trait$Value)),]
Trait = Trait[which(Trait$FullName != "unknown"),]
Trait = Trait[which(Trait$FullName != "Coffea arabica" & Trait$FullName != "Glycine max"),]
Trait[,6:7] <- data.frame(do.call('rbind', strsplit(as.character(Trait$FullName),' ',fixed=TRUE)))
colnames(Trait) <- c('DataID','spec.name','Value','DataName','TraitID','Genus','Species')
keep = grepl('^\\w+\\s\\w+$', Trait$spec.name)
Trait = Trait[keep,]

#2. Match names against WFO
NameMatch = WFO.match(spec.data=Trait, WFO.file = "C:/PhyloTraitEst/WFO/classification.csv", counter=1, verbose=TRUE,Fuzzy=0)
NameMatch = NameMatch[which(NameMatch$taxonomicStatus=="Accepted"),]
NameMatch = NameMatch[!duplicated(NameMatch$DataID),]

#3. Run V.Phylomaker2 scenario S3
WFO_Names = data.frame(species = NameMatch$scientificName, genus = NameMatch$Genus, family = NameMatch$family)
phylotree <- phylo.maker(WFO_Names, scenarios=c("S3"))

#4. Write Newick Tree
write.tree(phylotree$scenario.3,paste(inpath,"TRY_Values_NewickTrees/TRY_P12_Newick.txt",sep=""))

#5. Write Raw Data
TraitExport = unique(data.frame(NameMatch$DataID,NameMatch$scientificName,NameMatch$family,NameMatch$genus,NameMatch$specificEpithet,NameMatch$Value))
colnames(TraitExport) <- c('DataID','spec.name','Family','Genus','Species','Value')
write.table(TraitExport,paste(inpath,"TRY_Values_NewickTrees/TRY_P12_Values.csv",sep=""),row.names = FALSE,sep=",",quote = FALSE)
