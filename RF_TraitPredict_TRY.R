###########################################################################
# Train RF models to impute trait values based on PEMs and TRY observations
# Created: James Knighton (7/25/2023)
###########################################################################

library(h2o)

inpath = "C:/PhyloTraitEst/RF_Model/"
TraitName = "gsVPD"

h2o.init(nthreads = 6, #Number of threads -1 means use all cores on your machine
         max_mem_size = "2G")  #max mem size is the maximum memory to allocate to H2O

#1. Import TRY training dataset
TraitData <- paste(inpath,"RF_Train_",TraitName,".csv",sep="")
data <- h2o.importFile(TraitData)
dim(data)

splits <- h2o.splitFrame(data = data, 
                         ratios = c(0.7, 0.15),  #partition data into 70%, 15%, 15% chunks
                         seed = 6)  #setting a seed will guarantee reproducibility
train <- splits[[1]]
valid <- splits[[2]]
test <- splits[[3]]

#2. Identify response and predictor variables
y <- "Trait_Eigen.Value"
x <- setdiff(names(data),c(y,"Trait_Eigen.spec.name"))
print(x)

#3. Random Forest
RF <- h2o.randomForest(x = x,
                            y = y,
                            training_frame = train,
                            model_id = "RF_Trait",
                            nfolds = 8,
                            validation_frame = valid,  #only used if stopping_rounds > 0
                            ntrees = 300,
                            max_depth = 50,
                            stopping_rounds = 10,
                            stopping_tolerance = 1e-3,
                            stopping_metric = "MSE",
                            seed = 200)

RF_perf <- h2o.performance(model = RF, newdata = test)
RF_perf

summary(RF)

#4. Random Forest Predictions
pred_train = h2o.predict(object=RF, newdata=train)
pred_valid = h2o.predict(object=RF, newdata=valid)
pred_test = h2o.predict(object=RF, newdata=test)

Output = as.data.frame(train$Trait_Eigen.Value)
Output[,2] = as.data.frame(pred_train)
Output_V = as.data.frame(valid$Trait_Eigen.Value)
Output_V[,2] = as.data.frame(pred_valid)
Output_T = as.data.frame(test$Trait_Eigen.spec.name)
Output_T[,2] = as.data.frame(test$Trait_Eigen.Value)
Output_T[,3] = as.data.frame(pred_test)

#5a. Scatter plots of observed versus simulated
plotmin = 0
plotmax = 1e-2
par(mfrow=c(1,3))
plot(c(plotmin,plotmax),c(plotmin,plotmax),type="l")
points(Output$Trait_Eigen.Value,Output$predict,xlim=c(plotmin,plotmax),ylim=c(plotmin,plotmax))
plot(c(plotmin,plotmax),c(plotmin,plotmax),type="l")
points(Output_V$Trait_Eigen.Value,Output_V$predict,xlim=c(plotmin,plotmax),ylim=c(plotmin,plotmax))
plot(c(plotmin,plotmax),c(plotmin,plotmax),type="l")
points(Output_T$Trait_Eigen.Value,Output_T$predict,xlim=c(plotmin,plotmax),ylim=c(plotmin,plotmax))

#5b. Write test dataset performance
write.table(Output_T,paste(inpath,"RF_Test_Predict_",TraitName,".csv",sep=""),row.names = FALSE,col.names = FALSE,sep=",")
