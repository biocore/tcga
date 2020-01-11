# age_regression_cfdna_shogun.R

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
require(purrr) # for functional programming using map()
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(gbm) # for machine learning
require(tibble) # for df operations
require(cowplot) # for plotting
require(PRROC) # for precision-recall curves
require(MLmetrics) # for multi-class learning
require(gmodels) # for CI estimates


numCores <- detectCores()
registerDoMC(cores=numCores)

defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                               n.trees = floor((1:3) * 50),
                               shrinkage = 0.1,
                               n.minobsinnode = 5)
customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                              n.trees = floor((1:3) * 50),
                              shrinkage = 0.1,
                              n.minobsinnode = 1)

## Load data
# load("snmKrakenAndMetadataFiltered_Dec2_Final.RData")
# load("snmCfdnaShogunAndMetadata_Dec2_Final.RData")

load("shogun_vbAndMetadataFiltered_Dec2_Final.RData")


loocvAge <- function(metaData, readData, caretTuneGrid = defaultGBMGrid){
  
  # Do LOOCV model building and testing
  
  multiClassSummaryStats <- list()
  multiClassSummaryStatsDist <- list()
  numKFold <- 4
  numResampleIter <- 1
  iterSize <- 1
  
  for(jj in 1:iterSize){
    
    mlDataY <- droplevels(metaData[metaData$HvsC == "Control",])
    mlDataX <- readData[rownames(mlDataY),]
    dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
    
    # Create data partitions
    # set.seed(42)
    indexSuper <- 1:dim(mlDataY)[1]
    predProbs <- list()
    predAge <- list()
    obsClass <- list()
    predClass <- list()
    varImpBestModelDF2OrderedNonzeroList <- list()
    
    for(ii in 1:length(indexSuper)){
      print(sprintf("Iteration: %d/%d", ii, length(indexSuper)))
      index <- indexSuper[ii]
      # print(index)
      trainX <- mlDataX[-index,]
      trainY <- mlDataY[-index,]$host_age
      testX <- mlDataX[index,,drop=FALSE]
      testY <- mlDataY[index,,drop=FALSE]$host_age
      # print(testY)
      
      refactoredTrainY <- trainY
      refactoredTestY <- testY
      
      obsClass[ii] <- refactoredTestY
      
      set.seed(42)
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           # sampling = "up",
                           # summaryFunction = multiClassSummary,
                           # classProbs = TRUE,
                           verboseIter = FALSE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      mlModel <- train(x = trainX,
                       y = refactoredTrainY,
                       method = "gbm",
                       preProcess = c("scale","center"),
                       trControl = ctrl,
                       verbose = TRUE,
                       metric = "MAE",
                       tuneGrid = caretTuneGrid)
      
      predAge[ii] <- list(predict(mlModel, newdata = testX))
      predClass[ii] <- as.character(predict(mlModel, newdata = testX, type = "raw"))
      
      rm(mlModel)
    }
    
    predAgePerf <<- data.frame(obs = do.call(c, obsClass),
                              pred = do.call(c, predAge))
    print(postResample(pred = predAgePerf$pred, obs = predAgePerf$obs))
    ggplot(predAgePerf, aes(x = obs, y = pred)) + geom_point() + coord_equal() +
      geom_smooth(method='lm', formula= y~x) + 
      labs(x = "Observed age (yr)", y = "Predicted age (yr)", 
           title = "Predicting age using raw data and GBMs (LOOCV)") -> p
    
    print(predAgePerf)
    print(p)
    

  }
  
  return(predAgePerf)
}

## Kraken age predictions
loocvAge(metaData = metadataPSMatchedDPQCFiltered,
         readData = vbMergedDataKrakenDPCFDecontam)
krakenAgePer <- predAgePerf

save(krakenAgePer, file = "krakenAgePer.RData")

postResample(pred = krakenAgePer$pred, obs = krakenAgePer$obs)

ggscatter(krakenAgePer, 
          x = "obs", 
          y = "pred",
          xlab = "Observed age (yr)",
          ylab = "Predicted age (yr)") +
  coord_equal() + geom_smooth(method='lm', formula= y~x)

## Shogun age predictions
loocvAge(metaData = metadataPSMatchedDPQCFiltered,
         readData = shogundataPSUniqueDPDecontamQC)
shogunAgePer <- predAgePerf

save(shogunAgePer, file = "shogunAgePer.RData")

postResample(pred = shogunAgePer$pred, obs = shogunAgePer$obs)

ggscatter(shogunAgePer, 
          x = "obs", 
          y = "pred",
          xlab = "Observed age (yr)",
          ylab = "Predicted age (yr)") +
  coord_equal() + geom_smooth(method='lm', formula= y~x)

## Combined age predictions
allAgePerf <- rbind(data.frame(krakenAgePer, dataset = "Kraken"),
                    data.frame(shogunAgePer, dataset = "Shogun"))

write.csv(allAgePerf, file = "../AARF/edf8g.csv")

ggscatter(allAgePerf, 
          x = "obs", 
          y = "pred",
          palette = "nejm",
          color = "dataset",
          legend.title = "Dataset",
          xlab = "Observed age (yr)",
          ylab = "Predicted age (yr)") +
  annotate(geom = "text", x = 50, y = 32, 
           label = "Kraken MAE: 11.9650 yr\nShogun MAE: 10.6174 yr") -> allAgeRegressionPlot

ggsave(plot = allAgeRegressionPlot,
       filename = "allAgeRegressionPlot_Dec13.svg",
       width = 9,
       units = "in",
       dpi = "retina")



