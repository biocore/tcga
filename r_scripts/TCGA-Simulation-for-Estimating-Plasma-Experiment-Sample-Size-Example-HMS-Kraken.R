# simulationLOOCVPred1VsAllFullDataBDN.R
# Author: Greg Poore
# Date: Sept 28, 2019
# Purpose: Simulate model performance 

#-------------------------------#
# Load dependencies
require(devtools)
require(doMC)
require(tibble)
require(gbm)
require(splitstackshape)
require(reshape2)
require(ggpubr)
require(caret) # for model building
require(pROC) # for AUC calculations
require(DMwR) # for SMOTE class imbalance correction
require(ROSE) # for ROSE class imbalance correction
require(purrr) # for functional programming using map()
require(plyr) # for data manipulation
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(gbm) # for machine learning
require(tibble) # for df operations
require(cowplot) # for plotting
require(PRROC) # for precision-recall curves
require(MLmetrics) # for multi-class learning

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
# THE IMPORTED DATA SHOULD CONTAIN 3 R OBJECTS:
# 1. METADATA FOR A PARTICULAR SEQUENCING CENTER AND ASSIGNMENT ALGORITHM (E.G. HMS-KRAKEN METADATA, AS BELOW)
# 2. VOOM-SNM NORMALIZED DATA FOR A PARTICULAR SEQUENCING CENTER AND ASSIGNMENT ALGORITHM (E.G. HMS-KRAKEN DATA, AS BELOW)
# 3. A LIST OF SAMPLE-IDS FOR EACH PERMUTATION
# NB: #3 COULD BE GENERATED HERE AS WELL, BUT NOTE THAT R SOMETIMES HAS PROBLEMS MAKING RANDOM STRATIFIED
# SAMPLES WHEN RUNNING PARALLEL COMPUTING, WHICH IS NECESSARY TO SPEED UP MODEL TRAINING/TESTING
# NB: AN .RDATA FILE COULD ALSO BE LOADED, AS BELOW:
# load("simulationMaterials_Final_Nov21_KRAKEN.RData")

#------------------------------------------------------


## To avoid the error below when training the GBM learner, a custom GBM grid can be entered
# The dataset size is too small or subsampling rate is too large: nTrain*bag.fraction <= n.minobsinnode
defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                               n.trees = floor((1:3) * 50),
                               shrinkage = 0.1,
                               n.minobsinnode = 5)

customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                              n.trees = floor((1:3) * 50),
                              shrinkage = 0.1,
                              n.minobsinnode = 1)

# Model setup
sampleSizes <- seq(5,40,5)
metadataSim <- simulationMetadataQCHMS_Kraken ## MODIFY AS NEEDED
snmDataSim <- snmDataSimulationHMS_Kraken ## MODIFY AS NEEDED
caretTuneGrid <- customGBMGrid ## MODIFY AS NEEDED
numKFold <- 4
numResampleIter <- 1
iterSize <- 10

# SETTING UP LIST VARIABLES FOR LATER INPUT
multiClassSummaryStats <- list()
multiClassSummaryStatsDist <- list()
varImpBestModelDF2OrderedNonzeroListStackMedianLOOCV <- list()
varImpBestModelDF2OrderedNonzeroListStackMedianDistMedian <- list()

for(kk in 1:8){ # 5-40 by steps of 5

  sampleSizeSelected <- sampleSizes[kk]
  
  for(jj in 1:iterSize){

    sampleIDs <- krakenHMS_LIST[[kk]][[jj]] ## MODIFY AS NEEDED
    
    mlDataY <- droplevels(metadataSim[sampleIDs,])
    mlDataX <- snmDataSim[rownames(mlDataY),]
    dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
    
    indexSuper <- 1:dim(mlDataY)[1]
    predProbs <- list()
    obsClass <- vector()
    predClass <- vector()
    varImpBestModelDF2OrderedNonzeroList <- list()
    
    # Do LOOCV model building and testing
    for(ii in 1:length(indexSuper)){
      index <- indexSuper[ii]
      print(index)
      trainX <- mlDataX[-index,]
      trainY <- mlDataY[-index,]$disease_type_consol
      testX <- mlDataX[index,,drop=FALSE]
      testY <- mlDataY[index,,drop=FALSE]$disease_type_consol
      print(testY)
      
      refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
      refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
      
      obsClass[ii] <- as.character(refactoredTestY)
      
      set.seed(42)
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = multiClassSummary,
                           classProbs = TRUE,
                           verboseIter = TRUE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      mlModel <- train(x = trainX,
                       y = refactoredTrainY,
                       method = "gbm",
                       preProcess = c("scale","center"),
                       trControl = ctrl,
                       verbose = TRUE,
                       # metric = "ROC",
                       tuneGrid = caretTuneGrid)
      
      predProbs[ii] <- list(predict(mlModel, newdata = testX, type = "prob"))
      predClass[ii] <- as.character(predict(mlModel, newdata = testX, type = "raw"))
      
      varImpBestModelDF <- as.data.frame(varImp( mlModel$finalModel, scale = FALSE ))
      varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
      varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
      colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
      varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
      varImpBestModelDF2OrderedNonzeroList[[ii]] <- varImpBestModelDF2OrderedNonzero

      rm(mlModel)
    }
    
    # CLASSES AS FOLLOWS: NSCLC = GROUPED LUNG CANCER | PRAD = PROSTATE CANCER | SKCM = MELANOMA
    classes <- c("NSCLC", "PRAD", "SKCM")
    loocvPreds <- cbind(obs = factor(obsClass,
                                     levels = classes),
                        pred = factor(predClass,
                                      levels = classes),
                        do.call(rbind,predProbs))
    multiClassSummaryStats[[jj]] <- multiClassSummary(loocvPreds, lev = classes)
    print(multiClassSummaryStats[[jj]])
  }
  
  multiClassSummaryStatsDist[[kk]] <- data.frame(do.call(rbind, multiClassSummaryStats),
                                                 sampleSize = sampleSizeSelected)
  print(multiClassSummaryStatsDist[[kk]])
  csvFilename <- paste0("multiClassSummaryStatsDist_SampleSize_",sampleSizeSelected,
                        ".csv")
  write.csv(multiClassSummaryStatsDist[[kk]],file = csvFilename)
  
  multiClassSummaryStatsDistInterim <- data.frame(do.call(rbind, multiClassSummaryStatsDist))
  multiClassSummaryStatsDistInterim.melted <- melt(multiClassSummaryStatsDistInterim,
                                               id.vars = "sampleSize")
  
  multiClassSummaryStatsDistInterimPerfPlot <- 
    ggline(multiClassSummaryStatsDistInterim.melted,
         x = "sampleSize",
         y = "value",
         xlab = "Sample Size",
         add = c("mean_se"),
         facet.by = "variable")
  
  ggsave(filename = "BBB_multiClassSummaryStatsDistInterimPerfPlot.png",
         plot = multiClassSummaryStatsDistInterimPerfPlot,
         dpi = "retina",
         width = 12, 
         units = "in")
  
}

multiClassSummaryStatsDistAll <- data.frame(do.call(rbind, multiClassSummaryStatsDist))
multiClassSummaryStatsDistAll.melted <- melt(multiClassSummaryStatsDistAll,
                                             id.vars = "sampleSize")
csvFilenameFinal <- "AAA_multiClassSummaryStatsDistAll.csv"
write.csv(multiClassSummaryStatsDistAll,file = csvFilenameFinal)

multiClassSummaryStatsDistAllPerfPlot <- 
  ggline(multiClassSummaryStatsDistAll.melted,
            x = "sampleSize",
            y = "value",
            xlab = "Sample Size",
            add = c("mean_se"),
            facet.by = "variable")

ggsave(filename = "AAA_multiClassSummaryStatsDistAllPerfPlot.png",
       plot = multiClassSummaryStatsDistAllPerfPlot,
       dpi = "retina",
       width = 12, 
       units = "in")




#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------