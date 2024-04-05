#####################################################################
# Reformat gsMAX TRY data for further analysis
#Created: James Knighton (7/25/2023)
#####################################################################

library(rtry)
library(WorldFlora)
library(V.PhyloMaker2)
library(ape)
library(dplyr)

inpath = "D:/PhyloTraitEst/"

#0. Import TRY data (downloaded 7/25/2023)
TRYdata1 <- rtry_import(paste(inpath,"TRY_RawData/27967.txt",sep=""))

#1a. Extract TRY data and clean
Trait = data.frame(TRYdata1$ObsDataID,TRYdata1$AccSpeciesName, TRYdata1$StdValue, TRYdata1$DataName, TRYdata1$TraitID,TRYdata1$DatasetID,TRYdata1$Reference, TRYdata1$ErrorRisk)
colnames(Trait) <- c('DataID','FullName','Value','DataName','TraitID','DatasetID','Ref','ErrorRisk')
Trait = Trait[which(Trait$ErrorRisk <= 5),]
Trait = Trait[which(Trait$TraitID == 45),]
datatypes = unique(Trait$DataName)
keeps = c("Stomata conductance to water vapour per leaf area","Stomata conductance per leaf area at Asat measurement","stomata conductance to water vapour under conditions of light saturation, high relative humidity, and growth CO2 concentration.","stomata conductance per leaf area","Stomata conductance per leaf area max","Max of max stomata conductance","Maximum stomata conductance per leaf area","stomata conductance in the limit of saturating light, zero VPD, maximum soil water content and at the reference value of Cs = 350 ppm","Stomata conductance per leaf area at Amax measurement (saturating light and saturating CO2)","Stomat conductance at Amax")
Trait = Trait[which(Trait$DataName %in% keeps),]
Trait = Trait[which(Trait$Value>0),]
Trait = Trait[which(!is.na(Trait$Value)),]
Trait = Trait[which(Trait$FullName != "unknown"),]
keep = grepl('^\\w+\\s\\w+$', Trait$FullName)
Trait = Trait[keep,]
Trait[,7:8] <- data.frame(do.call('rbind', strsplit(as.character(Trait$FullName),' ',fixed=TRUE)))
colnames(Trait) <- c('DataID','spec.name','Value','DataName','TraitID','DatasetID','Genus','Species')
Trait = Trait[which(Trait$Species != "sp"),]

#1b. remove species with very high variance gsmax observations (i.e., potentially mislabeled observations?)
temp = Trait %>% group_by(spec.name) %>% summarise(SD = sd(Value,na.rm=TRUE), Mean = mean(Value,na.rm=TRUE), N = length(Value))
temp$cv = temp$SD/temp$Mean
highvar = unique(temp$spec.name[which(temp$cv > 0.6)])
Trait = Trait[(which(!Trait$spec.name %in% highvar)),]

#2. Match names against WFO
NameMatch = WFO.match(spec.data=Trait, WFO.file = "C:/PhyloTraitEst/WFO/classification.csv", counter=1, verbose=TRUE, Fuzzy=0)
NameMatch = NameMatch[which(NameMatch$taxonomicStatus=="Accepted"),]
NameMatch = NameMatch[!duplicated(NameMatch$DataID),]

#3. Run V.Phylomaker2 scenario S3
WFO_Names = data.frame(species = NameMatch$scientificName, genus = NameMatch$Genus, family = NameMatch$family)
phylotree <- phylo.maker(WFO_Names, scenarios=c("S3"))

#4. Write Newick Tree
write.tree(phylotree$scenario.3,paste(inpath,"TRY_Values_NewickTrees/TRY_gsmax_Newick.txt",sep=""))

#5. Write Raw Data
TraitExport = unique(data.frame(NameMatch$DataID,NameMatch$scientificName,NameMatch$family,NameMatch$genus,NameMatch$specificEpithet,NameMatch$Value))
colnames(TraitExport) <- c('DataID','spec.name','Family','Genus','Species','Value')
write.table(TraitExport,paste(inpath,"TRY_Values_NewickTrees/TRY_gsmax_Values.csv",sep=""),row.names = FALSE,sep=",",quote = FALSE)
