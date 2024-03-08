#######################################################################################################
# Impute trait values based on PEMs for BCGI tree names
# Solving PEMs for BCGI is memory limited due to the large number of species (>55K)
# This script circumvents this limitation by batching the database
# Note: this step is very very very computationally expensive. This code is best adapted for parallel execution on HPC
# Created: James Knighton (7/25/2023)
#######################################################################################################

library(h2o)

h2o.init(nthreads = 6, #Number of threads -1 means use all cores on your machine
         max_mem_size = "4G")  #max mem size is the maximum memory to allocate to H2O

inpath = "C:/PhyloTraitEst/RF_Model/GlobalPredictSet/"
TraitName = "WUE"

for (chunk in 0:58)
{
  print(paste("Solving ",chunk," of 58",sep=""))

  #1. Import training and prediction files
  if (chunk == 0)
  {
    trainfile <- paste(inpath,"RF_",TraitName,"_1_train.csv",sep="")
    predfile <- paste(inpath,"RF_",TraitName,"_1_train.csv",sep="")
  }
  else
  {
    trainfile <- paste(inpath,"RF_",TraitName,"_",chunk,"_train.csv",sep="")
    predfile <- paste(inpath,"RF_",TraitName,"_",chunk,"_predict.csv",sep="")
  }  
  
  traindata <- h2o.importFile(trainfile)
  preddata <- h2o.importFile(predfile)    
  
  splits <- h2o.splitFrame(data = traindata, 
                           ratios = c(0.7, 0.15),  #partition data into 70%, 15%, 15% chunks
                           seed = 5)  #setting a seed will guarantee reproducibility
  train <- splits[[1]]
  valid <- splits[[2]]
  test <- splits[[3]]
  
  #2. Identify response and predictor variables
  y <- "Trait_Eigen.Value"
  x <- setdiff(names(traindata),c(y,"Trait_Eigen.spec.name"))
  
  #3. Random Forest
  RF <- h2o.randomForest(x = x,
                              y = y,
                              training_frame = train,
                              model_id = "RF_Trait",
                              nfolds = 10,
                              validation_frame = valid,
                              ntrees = 300,
                              max_depth = 30,
                              #stopping_rounds = 10,
                              stopping_tolerance = 1e-3,
                              stopping_metric = "MSE",
                              seed = 200)
  
  predtraits = h2o.predict(object=RF, newdata=preddata)
  
  #4. Random Forest Predictions
  Output = as.data.frame(preddata$Trait_Eigen.spec.name)
  Output[,2] = as.data.frame(predtraits)
  
  if (chunk == 0) {All = Output}
  if (chunk > 0) {All = rbind(All,Output)}
}

#5. Write trait predictions
write.table(All,paste(inpath,"RF_Global_Trees_",TraitName,".csv",sep=""),row.names = FALSE,sep=",")

