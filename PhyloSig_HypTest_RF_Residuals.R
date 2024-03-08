#####################################################################
# Pagel's lambda test for RF model trait residuals
# Created: James Knighton (7/25/2023)
#####################################################################

library(phytools)
library(Taxonstand)
library(diversitree)
library(ape)
library(caper)
library(V.PhyloMaker2)
library(WorldFlora)

inpath = "C:/PhyloTraitEst/RF_Model/"
outpath = "C:/PhyloTraitEst/RF_Model/"
TraitName = "gsmax"

#1.Read in Newick Tree and TRY data
Trait <- read.table(paste(inpath,"RF_Test_Predict_",TraitName,".csv",sep=""),sep=",",header=FALSE)
Trait$MatchString = gsub(" ", "_", Trait$V1)
colnames(Trait) <- c('spec.name','Obs','Sim',"MatchString")

#2. Run V.Phylomaker2 scenario S3
NameMatch = WFO.match(spec.data=Trait, WFO.file = "C:/PhyloTraitEst/WFO/classification.csv", counter=1, verbose=TRUE, Fuzzy=0)

#3. Validate against WFO (should already be handled in earlier scripts, but just a double check on H2O output)
WFO_Names = data.frame(species = NameMatch$scientificName, genus = NameMatch$genus, family = NameMatch$family)
phylotree <- phylo.maker(WFO_Names, scenarios=c("S3"))
phylotree = phylotree$scenario.3

#4. Compute median species trait value from TRY database
for (i in 1:length(phylotree$tip.label))
{
  speciesname = phylotree$tip.label[i]
  Records = Trait[Trait$MatchString == speciesname,]
  phylotree$trait[i] = median(Records$Obs - Records$Sim)
}

#5. Pagel's lambda tests for phylogenetic signal
L1 = phylosig(phylotree, phylotree$trait, method="lambda",test=TRUE,nsim=1000, se=NULL, start=NULL,
              control=list())

#K1 = phylosig(phylotree, phylotree$trait, method="K",test=TRUE,nsim=1000, se=NULL, start=NULL,
#              control=list())
