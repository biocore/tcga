# loocv_subsampling100.R
# Gregory Poore
# Dec 9 2019
# Purpose: See how subsampling affects the perf of LOO machine learning models

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
require(caret) # for machine learning
require(gmodels) # for CI interval calc

numCores <- detectCores()
registerDoMC(cores=numCores)

## Load data
# LOAD KRAKEN- AND SHOGUN-DERIVED PLASMA DATA R OBJECTS
# CAN LOAD MATRICES / DATA FRAMES OR TOGETHER USING A .RDATA FILE, AS SHOWN BELOW:
# load("snmCfdnaShogunAndMetadata_Dec2_Final.RData")
# load("snmKrakenAndMetadataFiltered_Dec2_Final.RData")

## DEFINE HYPERPARAMETER SEARCH GRIDS FOR MACHINE LEARNING
# THE defaultGBMGrid WAS USED BELOW BUT CAN BE ALTERED, AS SHOWN
# BY THE customGBMGrid EXAMPLE
defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                               n.trees = floor((1:3) * 50),
                               shrinkage = 0.1,
                               n.minobsinnode = 5)
customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                              n.trees = floor((1:3) * 50),
                              shrinkage = 0.1,
                              n.minobsinnode = 1)
numKFold <- 4
numResampleIter <- 1
iterSize <- 100 # I.E. SUBSAMPLING ITERATIONS

#-------------------------------------#
## Main function 
#-------------------------------------#

# NB: SMALLEST COHORT (SKCM) CONTAINED 16 TOTAL SAMPLES --> DEFINES THE SUBSAMPLING SIZE
loocvDTs_subsample <- function(snmData, samplingSize = 16, DTs, iterSize = 100, caretTuneGrid = defaultGBMGrid,
                               filenameString = paste(DTs,collapse = "__"), HvsCFlag = FALSE){
  
  if(HvsCFlag){ # HvsCFlag IS TO COMPARE HEALTHY CONTROLS VS GROUPED CANCER
    # metadataPSMatchedDPQCFiltered IS A LOADED METADATA DATA FRAME CONTAINING A COLUMN CALLED disease_type_consol
    # WITH THE SAMPLE TYPE DESIGNATION (E.G. PRAD). HvsC IS ANOTHER COLUMN CONTAINING A "HEALTHY" OR "CANCER" LABEL
    # FOR ALL SAMPLES
    metaTmpX <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in% DTs),])
    metaTmpX$disease_type_consol <- metaTmpX$HvsC
    classes <- gsub(" ","",levels(metaTmpX$disease_type_consol))
  } else{
    metaTmpX <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in% DTs),])
    classes <- gsub(" ","",DTs)
  }
  
  rownameSubsampleList <- t(lapply(1:iterSize, function(x) as.data.frame(stratified(metaTmpX,
                                                                   group = "disease_type_consol",
                                                                   size = samplingSize,
                                                                   keep.rownames = TRUE,
                                                                   replace = FALSE,
                                                                   bothSets = FALSE))$rn))
  
  # Do LOOCV model building and testing
  
  multiClassSummaryStats <- list()
  multiClassSummaryStatsDist <- list()
  numKFold <- 4
  numResampleIter <- 1
  metaData <- metaTmpX
  snmData <- snmData 
  iterSize <- iterSize
  for(jj in 1:iterSize){
    print(sprintf("Sub-Sampling Iteration: %d/%d", jj, iterSize))

    mlDataY <- droplevels(metaData[rownameSubsampleList[[jj]],])
    mlDataX <- snmData[rownames(mlDataY),]
    dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
    
    print(rownames(mlDataY))
    
    indexSuper <- 1:dim(mlDataY)[1]
    predProbs <- list()
    obsClass <- vector()
    predClass <- vector()
    varImpBestModelDF2OrderedNonzeroList <- list()
    
    for(ii in 1:length(indexSuper)){
      print(sprintf("LOOCV Iteration: %d/%d", ii, length(indexSuper)))
      index <- indexSuper[ii]
      # print(index)
      trainX <- mlDataX[-index,]
      trainY <- mlDataY[-index,]$disease_type_consol
      testX <- mlDataX[index,,drop=FALSE]
      testY <- mlDataY[index,,drop=FALSE]$disease_type_consol
      # print(testY)

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
                           verboseIter = FALSE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)

      mlModel <- train(x = trainX,
                       y = refactoredTrainY,
                       method = "gbm",
                       preProcess = c("scale","center"),
                       trControl = ctrl,
                       verbose = FALSE,
                       metric = "ROC",
                       tuneGrid = caretTuneGrid)

      predProbs[ii] <- list(predict(mlModel, newdata = testX, type = "prob"))
      predClass[ii] <- as.character(predict(mlModel, newdata = testX, type = "raw"))


      rm(mlModel)
    }

    # CALCULATE PERFORMANCE SUMMARY STATISTICS
    loocvPreds <- cbind(obs = factor(obsClass,
                                     levels = classes),
                        pred = factor(predClass,
                                      levels = classes),
                        do.call(rbind,predProbs))

    multiClassSummaryStats[[jj]] <- multiClassSummary(loocvPreds, lev = classes)

    print(confusionMatrix(loocvPreds$obs, loocvPreds$pred))
    print(multiClassSummaryStats[[jj]])
    
  }
  
  # MERGE PERFORMANCE SUMMARY STATISTICS
  multiClassSummaryStatsDist <- data.frame(do.call(rbind, multiClassSummaryStats))
  print(multiClassSummaryStatsDist)
  filenamePerfDist <- paste0(filenameString,"__iter",iterSize,"__PerfDist.csv")
  write.csv(multiClassSummaryStatsDist, file = filenamePerfDist)
  
  # CALCULATE 95% CONFIDENCE INTERVALS OF PERFORMANCE
  aucrocCI <- ci(multiClassSummaryStatsDist$AUC)
  auprCI <- ci(multiClassSummaryStatsDist$prAUC)
  
  print(paste0("Mean AUCROC = ",round(aucrocCI[1],4)," | 95% CI: [",
               round(aucrocCI[2],4),",",round(aucrocCI[3],4),"]"))
  
  print(paste0("Mean AUPR = ",round(auprCI[1],4)," | 95% CI: [",
               round(auprCI[2],4),",",round(auprCI[3],4),"]"))

  return(multiClassSummaryStatsDist)
  
}

#-------------------------------------#
## RUN FUNCTION
#-------------------------------------#

# RUN ABOVE FUNCTION FOR SUBSAMPLING EFFECTS ON PROSTATE CANCER VS LUNG CANCER DISCRIMINATION
prad_nsclc_ss16 <- loocvDTs_subsample(snmData = snmDataKrakenCFDecontamDPQC,
                   samplingSize = 16, 
                   iterSize = iterSize,
                   DTs = c("PRAD","NSCLC"),
                   caretTuneGrid = defaultGBMGrid)
save(prad_nsclc_ss16, file = "prad_nsclc_ss16_100.RData")

# RUN ABOVE FUNCTION FOR SUBSAMPLING EFFECTS ON PROSTATE CANCER VS HEALTHY CONTROLS DISCRIMINATION
prad_control_ss16 <- loocvDTs_subsample(snmData = snmDataKrakenCFDecontamDPQC,
                                      samplingSize = 16, 
                                      iterSize = iterSize,
                                      DTs = c("PRAD","Control"),
                                      caretTuneGrid = defaultGBMGrid)
save(prad_control_ss16, file = "prad_control_ss16_100.RData")

# RUN ABOVE FUNCTION FOR SUBSAMPLING EFFECTS ON LUNG CANCER VS HEALTHY CONTROLS DISCRIMINATION
nsclc_control_ss16 <- loocvDTs_subsample(snmData = snmDataKrakenCFDecontamDPQC,
                                        samplingSize = 16, 
                                        iterSize = iterSize,
                                        DTs = c("NSCLC","Control"),
                                        caretTuneGrid = defaultGBMGrid)
save(nsclc_control_ss16, file = "nsclc_control_ss16_100.RData")

# RUN ABOVE FUNCTION FOR SUBSAMPLING EFFECTS ON MELANOMA VS HEALTHY CONTROLS DISCRIMINATION
skcm_control_ss16 <- loocvDTs_subsample(snmData = snmDataKrakenCFDecontamDPQC,
                                         samplingSize = 16, 
                                         iterSize = iterSize,
                                         DTs = c("SKCM","Control"),
                                         caretTuneGrid = defaultGBMGrid)
save(skcm_control_ss16, file = "skcm_control_ss16_100.RData")

# MERGE PERFORMANCE DATA FOR PLOTTING
hVsCSubsample <- rbind(data.frame(prad_control_ss16, CT = "PC"),
                       data.frame(nsclc_control_ss16, CT = "LC"),
                       data.frame(skcm_control_ss16, CT = "SKCM"))
hVsCSubsample.melted <- melt(hVsCSubsample, id.vars = "CT")
levels(hVsCSubsample.melted$variable) <- c("Log Loss",
                                     "Mean AUROC",
                                     "Mean AUPR",
                                     "Accuracy",
                                     "Kappa",
                                     "Mean F1",
                                     "Mean Sensitivity",
                                     "Mean Specificity",
                                     "Mean Pos Pred Value",
                                     "Mean Neg Pred Value",
                                     "Mean Precision",
                                     "Mean Recall",
                                     "Mean Detection Rate",
                                     "Mean Balanced Accuracy"
)

# PLOT PERFORMANCE DISTRIBUTIONS FOR ALL 100 ITERATIONS AS BOXPLOTS
ggboxplot(hVsCSubsample.melted,
          x = "CT",
          y = "value",
          fill = "CT",
          legend.title = "Cancer Type",
          ylim = c(0.2,1.35),
          ncol = 7,
          panel.labs.font = list(size = 10),
          ylab = "Performance Value",
          facet.by = "variable") +
  stat_compare_means(comparisons = list(c("PC","LC"), c("PC","SKCM"), c("LC","SKCM")),
                     label = "p.signif") -> subsamplePlot

# SAVE PLOT DATA AS AN .RDATA FILE
save(hVsCSubsample.melted, subsamplePlot, file = "subsamplePlotAndData.RData")

# SAVE PLOT IN SVG FORMAT
ggsave(filename = "subsamplePlot_Dec9.svg",
       plot = subsamplePlot,
       dpi = "retina",
       width = 14, 
       height = 6,
       units = "in")