#####################################################################
# Reformat P50 TRY data for further analysis
#Created: James Knighton (7/25/2023)
#####################################################################

library(rtry)
library(WorldFlora)
library(V.PhyloMaker2)
library(ape)

inpath = "C:/PhyloTraitEst/"

#0. Import TRY data (downloaded 7/25/2023)
XFTdata <- read.csv(paste(inpath,"TRY_RawData/XFT_full_database_download_20240324-005823.csv",sep=""))

#1. Extract TRY data and clean
Trait = data.frame(XFTdata$rowID,XFTdata$Cleaned.binomial, XFTdata$Developmental.stage, XFTdata$Plant.organ, XFTdata$Curve, XFTdata$P50..MPa.,XFTdata$P12..MPa.,XFTdata$P88..MPa.)
Trait$XFTdata.P50..MPa. = as.numeric(Trait$XFTdata.P50..MPa.)
Trait = Trait[which(Trait$XFTdata.Developmental.stage == 'A'),]
Trait = Trait[which(Trait$XFTdata.Plant.organ == 'S'),]
Trait = Trait[which(Trait$XFTdata.Curve == 'S'),]
Trait = Trait[which(!is.nan(Trait$XFTdata.P50..MPa.)),]
Trait = Trait[which(!is.nan(Trait$XFTdata.Cleaned.binomial)),]
Trait = Trait[which(Trait$XFTdata.P50..MPa. <= -0.5),]
Out = data.frame(Trait$XFTdata.rowID,Trait$XFTdata.Cleaned.binomial,Trait$XFTdata.P50..MPa.,Trait$XFTdata.rowID,Trait$XFTdata.rowID)
Out[,6:7] <- data.frame(do.call('rbind', strsplit(as.character(Out$Trait.XFTdata.Cleaned.binomial),' ',fixed=TRUE)))
colnames(Out) <- c('DataID','spec.name','Value','DataName','TraitID','Genus','Species')
keep = grepl('^\\w+\\s\\w+$', Out$spec.name)
Trait = Trait[keep,]

#2. Match names against WFO
NameMatch = WFO.match(spec.data=Out, WFO.file = "C:/PhyloTraitEst/WFO/classification.csv", counter=1, verbose=TRUE,Fuzzy=0)
NameMatch = NameMatch[which(NameMatch$taxonomicStatus=="Accepted"),]
NameMatch = NameMatch[!duplicated(NameMatch$DataID),]

#3.Run V.Phylomaker2 scenario S3
WFO_Names = data.frame(species = NameMatch$scientificName, genus = NameMatch$Genus, family = NameMatch$family)
phylotree <- phylo.maker(WFO_Names, scenarios=c("S3"))

#4. Write Newick Tree
write.tree(phylotree$scenario.3,paste(inpath,"TRY_Values_NewickTrees/XFT_P50_Newick.txt",sep=""))

#5. Write Raw Data
TraitExport = unique(data.frame(NameMatch$DataID,NameMatch$scientificName,NameMatch$family,NameMatch$genus,NameMatch$specificEpithet,NameMatch$Value))
colnames(TraitExport) <- c('DataID','spec.name','Family','Genus','Species','Value')
write.table(TraitExport,paste(inpath,"TRY_Values_NewickTrees/XFT_P50_Values.csv",sep=""),row.names = FALSE,sep=",",quote = FALSE)
