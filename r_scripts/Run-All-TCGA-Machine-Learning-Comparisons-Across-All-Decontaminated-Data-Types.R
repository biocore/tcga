# mlFinalRe-Run_Barnacle_v3_Nov7_2019.R
# Author: Greg Poore
# Date: Nov 7, 2019
# Purpose: Final re-run

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
require(purrr) # for functional programming using map()
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

# LOAD AS MUCH OR AS LITTLE DATA AS DESIRED.
# ALL VOOM-SNM NORMALIZED DATA SHOULD BE MATCHED TO ONE METADATA DATA FRAME
# ANOTHER METADATA DATA FRAME CAN BE LOADED WITH GROUPED PATHOLOGY STAGE INFORMATION
# NB: MANY TCGA SAMPLES DID NOT CONTAIN PATHOLOGY INFORMATION AND THUS WERE DISCARDED
# WHEN EXAMINING STAGE I VS STAGE IV EFFECTS
# IT IS ALSO POSSIBLE TO LOAD .RDATA FILES, AS DONE BELOW:
# load("tcgaVbDataAndMetadataAndSNM.RData")
# load("snmDataSampleTypeWithExpStrategyFINAL.RData")
# load("snmDataDecontam4Datasets_Nov112019.RData")
# load("snmDataAndMetadataDecontamPlateCenter.RData")
# load("tcgaPathBinnedData.RData")

# GBM HYPERPARAMETER SEARCH GRID BELOW (DEFAULT PER CARET PACKAGE)
defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                               n.trees = floor((1:3) * 50),
                               shrinkage = 0.1,
                               n.minobsinnode = 5)

# Model setup -- MODIFY AS NEEDED; ALL THESE COMPARISONS WILL BE TESTED
sampleTypeList <- c("Stage I vs IV",
                    "Primary Tumor vs Solid Tissue Normal",
                    "Blood Derived Normal",
                    "Primary Tumor")
# MODIFY AS NEEDED
datasetList <- list(snmDataSampleTypeWithExpStrategy,
  snmDataLikelyContamRemovedFA,
  snmDataAllPutativeContamRemovedFA,
  snmDataPlateCenterContamRemovedFA,
  snmDataMostStringentContamRemovedFA)

# MODIFY AS NEEDED
datasetListNames <- c("snmDataSampleTypeWithExpStrategy",
  "snmDataLikelyContamRemovedFA",
  "snmDataAllPutativeContamRemovedFA",
  "snmDataPlateCenterContamRemovedFA",
  "snmDataMostStringentContamRemovedFA")

# MATCHED METADATA DATA FRAME FOR ALL COUNT DATA:
metaTmpQC <- metadataSamplesAllQC # MODIFY AS NEEDED
# PATHOLOGY METADATA DATA FRAME:
metaTmpPath <- metadataSamplesAllQCPath # MODIFY AS NEEDED
caretTuneGrid <- defaultGBMGrid
numKFold <- 4
numResampleIter <- 1
prroc_roc <- list()
prroc_pr <- list()
perf <- list()
perfTmp <- list()
perfTmp2 <- list()

for(jj in seq_along(datasetList)){
  dataTmp <- datasetList[[jj]]

  datasetName <- datasetListNames[[jj]]

  for(kk in seq_along(sampleTypeList)){
  st <- sampleTypeList[kk]
  
  if(st == "Primary Tumor vs Solid Tissue Normal"){
    metaTmp <- metaTmpQC
    metaTmp2 <- droplevels(metaTmp[metaTmp$sample_type %in% c("Primary Tumor",
                                                              "Solid Tissue Normal"),])
  } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
    metaTmp <- metaTmpQC
    metaTmp2 <- droplevels(metaTmp[metaTmp$sample_type %in% st,])
  } else if(st == "Stage I vs IV"){
    metaTmp <- metaTmpPath
    metaTmp2 <- droplevels(metaTmp[(metaTmp$sample_type %in% "Primary Tumor") &
                                     (metaTmp$pathologic_stage_label_binned %in% c("Stage1","Stage4")),])
  }

    for(ii in seq_along(levels(metaTmp2$disease_type))){
      
      dt <- levels(metaTmp2$disease_type)[ii]
      print(dt)
      
      if(st == "Primary Tumor vs Solid Tissue Normal"){
        metaTmp3 <- droplevels(metaTmp2[metaTmp2$disease_type == dt,])
        metaTmp3$predY <- factor(gsub('([[:punct:]])|\\s+','',metaTmp3$sample_type))
        positiveClass <- "PrimaryTumor"
        negativeClass <- "SolidTissueNormal"
      } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
        metaTmp3 <- metaTmp2
        metaTmp3$predY <- factor(ifelse(metaTmp2$disease_type == dt, 
                                        yes = dt, 
                                        no = "OtherCancerType"),
                                 levels = c(dt, "OtherCancerType"))
        positiveClass <- gsub('([[:punct:]])|\\s+','',dt)
        negativeClass <- "OtherCancerType"
      } else if(st == "Stage I vs IV"){
        metaTmp3 <- droplevels(metaTmp2[metaTmp2$disease_type == dt,])
        metaTmp3$predY <- metaTmp3$pathologic_stage_label_binned
        positiveClass <- "Stage4"
        negativeClass <- "Stage1"
      }
      
      print(table(metaTmp3$predY))
      
      # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
      if(length(table(metaTmp3$predY)) < 2){next}
      
      # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 20 SAMPLES IN EITHER CLASS
      if(any(table(metaTmp3$predY) < 20)){next}
      
      minorityClassSize <- min(table((metaTmp3$predY)))
      majorityClassSize <- max(table((metaTmp3$predY)))

      minorityClassName <- names(table(metaTmp3$predY)[which(table(metaTmp3$predY) == min(table(metaTmp3$predY)))])
      majorityClassName <- names(table(metaTmp3$predY)[which(table(metaTmp3$predY) == max(table(metaTmp3$predY)))])
      
      mlDataY <- metaTmp3
      mlDataX <- dataTmp[rownames(mlDataY),]
      
      # USE 70% OF DATA FOR TRAINING AND 30% FOR TESTING
      set.seed(42)
      index <- createDataPartition(mlDataY$predY, p = 0.7, list = FALSE)
      trainX <- mlDataX[index,]
      trainY <- mlDataY[index,]$predY
      testX <- mlDataX[-index,]
      testY <- mlDataY[-index,]$predY
      # print(testY)
      
      refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
      refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
      
      set.seed(42)
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = twoClassSummary,
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
                       metric = "ROC",
                       tuneGrid = caretTuneGrid)
      
      predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
      fg <- predProbs[refactoredTestY == positiveClass]
      bg <- predProbs[refactoredTestY == negativeClass]
      
      # CALCULATE ROC AND PR CURVE INFORMATION
      prroc_roc[[ii]] <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
      prroc_pr[[ii]] <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)

      rocCurveData <- cbind(as.data.frame(prroc_roc[[ii]]$curve), disease_type = dt, sample_type = st)
      prCurveData <- cbind(as.data.frame(prroc_pr[[ii]]$curve), disease_type = dt, sample_type = st)
      
      # SUMMARIZE MODEL PERFORMANCES
      perf[[ii]] <- data.frame(diseaseType = dt,
                               sampleType = st,
                               datasetName = datasetName,
                               minorityClassSize = minorityClassSize,
                               majorityClassSize = majorityClassSize,
                               minorityClassName = minorityClassName,
                               majorityClassName = majorityClassName,
                               aucroc = prroc_roc[[ii]]$auc,
                               aupr = prroc_pr[[ii]]$auc.integral)
      
      print(perf[[ii]])
      
      confusionMatrix <- confusionMatrix(predict(mlModel, newdata = testX, type = "raw"), 
                                         refactoredTestY, 
                                         positive = positiveClass)
      
      print(confusionMatrix)
      
      #--------------------------------------#
      # Save performance into relevant files #
      #--------------------------------------#
      
      filepathPerfPlots <- paste0("./perfPlots__",datasetName)
      filepathPerfPlotsDataPR <- paste0("./dataPR__",datasetName)
      filepathPerfPlotsDataROC <- paste0("./dataROC__",datasetName)
      filepathFeatures <- paste0("./features__",datasetName)
      filepathPerfStats <- paste0("./stats__",datasetName)
      filepathConfusionMatrix <- paste0("./confusionMatrix__",datasetName)
      
      if(!( dir.exists( file.path(filepathPerfPlots)))){
        dir.create(file.path(filepathPerfPlots))
      }
      if(!( dir.exists( file.path(filepathPerfPlotsDataPR)))){
        dir.create(file.path(filepathPerfPlotsDataPR))
      }
      if(!( dir.exists( file.path(filepathPerfPlotsDataROC)))){
        dir.create(file.path(filepathPerfPlotsDataROC))
      }
      if(!( dir.exists( file.path(filepathFeatures)))){
        dir.create(file.path(filepathFeatures))
      }
      if(!( dir.exists( file.path(filepathPerfStats)))){
        dir.create(file.path(filepathPerfStats))
      }
      if(!( dir.exists( file.path(filepathConfusionMatrix)))){
        dir.create(file.path(filepathConfusionMatrix))
      }
      
      filenameROC <- paste0(filepathPerfPlots,"/",
                            dt,
                            " -- ",
                            st,
                            " -- ROC.png")
      
      filenamePR <- paste0(filepathPerfPlots,"/",
                           dt,
                           " -- ",
                           st,
                            " -- PR.png")

      filenameROCData <- paste0(filepathPerfPlotsDataROC,"/",
                            dt,
                            " -- ",
                            st,
                            " -- ROC.csv")
      
      filenamePRData <- paste0(filepathPerfPlotsDataPR,"/",
                           dt,
                           " -- ",
                           st,
                            " -- PR.csv")
      
      filenameCM <- paste0(filepathConfusionMatrix,"/",
                           dt,
                           " -- ",
                           st,
                           " -- CM.txt")
      
      filenameFeatures <- paste0(filepathFeatures,"/",
                            dt,
                            " -- ",
                            st,
                            " -- Features.csv")
      
      filenameCSV <- paste0(filepathPerfStats,"/",
                            dt,
                            " -- ",
                            st,
                            " -- Perf.csv")

      write.csv(perf[[ii]], file = filenameCSV)
      
      png(filename=filenameROC, width = 6, height = 4, units = 'in', res = 300)
      plot(prroc_roc[[ii]])
      dev.off()
      
      png(filename=filenamePR, width = 6, height = 4, units = 'in', res = 300)
      plot(prroc_pr[[ii]])
      dev.off()

      write.table(prCurveData, sep=",", file = filenamePRData, col.names = FALSE)
      write.table(rocCurveData, sep=",", file = filenameROCData, col.names = FALSE)
      
      # EXTRACT AND SAVE RANKED FEATURE IMPORTANCE INFORMATION
      varImpBestModelDF <- as.data.frame(varImp( mlModel$finalModel, scale = FALSE ))
      varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
      varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
      colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
      varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
      write.csv(varImpBestModelDF2OrderedNonzero, file = filenameFeatures, row.names = FALSE)
      
      sink(file = filenameCM)
      print(dt)
      print(confusionMatrix)
      sink()
      
      rm(mlModel)
      
    }

  # SUMMARIZE RESULTS
  perfTmp[[kk]] <- do.call(rbind, perf)

  }

  # SUMMARIZE RESULTS
  perfTmp2[[jj]] <- do.call(rbind, perfTmp)

  write.csv(perfTmp2[[jj]], file = paste0("perfFinalRe-Run_Nov11_",datasetName,".csv"))


}

# SUMMARIZE RESULTS
perf1VsAll <- do.call(rbind, perfTmp2)

write.csv(perf1VsAll, file = "perfFinalRe-Run_Nov11.csv")



#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------