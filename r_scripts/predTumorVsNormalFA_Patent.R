## Goal: Predict tumor vs normal
predTumorVsNormalFA_Patent <- function(qcMLMetadata, 
                                qcMLDataSNM,
                                cancerTypeString = NULL,
                                nlTissueComparison = "Solid Tissue Normal",
                                tumorTissueComparison = "Primary Tumor", 
                                flagForAlternativeNormal = FALSE,
                                cancerTypeStringAlternativeNormal = NULL,
                                caretModel = "gbm",
                                samplingStrategy = "up",
                                numResampleIter = 2,
                                numKFold = 3,
                                trainSetProp = 0.7,
                                caretTuneGrid = NULL,
                                plotPath = "./provisional-patent-roc-predicting-cancer-tumor-vs-normal",
                                modelFeaturesPath = "./provisional-patent-csv-model-features-cancer-tumor-vs-normal",
                                functionOutputPath = "./provisional-patent-function-output-cancer-tumor-vs-normal",
                                confusionMatrixPath = "./provisional-patent-confusion-matrices-cancer-tumor-vs-normal",
                                confusionMatrixStatsPath = "./provisional-patent-stats-confusion-matrices-cancer-tumor-vs-normal",
                                ...){
  
  require(caret) # for model building
  require(pROC) # for AUC calculations
  require(DMwR) # for SMOTE class imbalance correction
  require(ROSE) # for ROSE class imbalance correction
  require(purrr) # for functional programming using map()
  require(dplyr) # for data manipulation
  require(doMC) # for parallel computing
  require(gbm) # for machine learning
  require(tibble) # for df operations
  require(cowplot) # for plotting
  
  numCores <- detectCores()
  registerDoMC(cores=numCores)

  if(!( dir.exists( file.path(plotPath)))){
    dir.create(file.path(plotPath))
  }

  if(!( dir.exists( file.path(modelFeaturesPath)))){
    dir.create(file.path(modelFeaturesPath))
  }
  if(!( dir.exists( file.path(functionOutputPath)))){
    dir.create(file.path(functionOutputPath))
  }

  if(!( dir.exists( file.path(confusionMatrixPath)))){
    dir.create(file.path(confusionMatrixPath))
  }
  if(!( dir.exists( file.path(confusionMatrixStatsPath)))){
    dir.create(file.path(confusionMatrixStatsPath))
  }

  sinkFilename <- sprintf("%s/%s - %s vs. %s (CV k-fold of %d|Train proportion of %d) Sink.txt", 
    functionOutputPath, cancerTypeString, tumorTissueComparison, nlTissueComparison, numKFold, trainSetProp*100)
  modelFeaturesFilename <- sprintf("%s/%s - %s vs. %s (CV k-fold of %d|Train proportion of %d) Model Features.csv", 
    modelFeaturesPath, cancerTypeString, tumorTissueComparison, nlTissueComparison, numKFold, trainSetProp*100)
  confusionMatrixFilename <- sprintf("%s/%s - %s vs. %s (CV k-fold of %d|Train proportion of %d) Confusion Matrix.csv", 
    confusionMatrixPath, cancerTypeString, tumorTissueComparison, nlTissueComparison, numKFold, trainSetProp*100)
  confusionMatrixStatsFilename <- sprintf("%s/%s - %s vs. %s (CV k-fold of %d|Train proportion of %d) Confusion Matrix Stats.txt", 
    confusionMatrixStatsPath, cancerTypeString, tumorTissueComparison, nlTissueComparison, numKFold, trainSetProp*100)
  
  sink(file = sinkFilename)
  
  cancerTypeRows <- (qcMLMetadata$disease_type == cancerTypeString)
  if(flagForAlternativeNormal == TRUE){
    cancerTypeAlternativeNormalRows <- ((qcMLMetadata$disease_type == cancerTypeStringAlternativeNormal) & (qcMLMetadata$sample_type == nlTissueComparison))
    cancerTypeRowsTotal <- (cancerTypeRows | cancerTypeAlternativeNormalRows)
    # rocTitle <- paste(paste("ROC curves",
    #                         paste0("(Learner: ",caretModel,"/","rep=",as.character(numResampleIter),"/","kFold=",as.character(numKFold),")\n")),
    #                   cancerTypeString,"\n",
    #                   paste0('(',tumorTissueComparison,' vs. ',nlTissueComparison,' from ',cancerTypeStringAlternativeNormal,')') )
  }
  else{
    cancerTypeRowsTotal <- cancerTypeRows
    # rocTitle <- paste(paste("ROC curves",
    #                         paste0("(Learner: ",caretModel,"/","rep=",as.character(numResampleIter),"/","kFold=",as.character(numKFold),")\n")),
    #                   cancerTypeString,"\n",
    #                   paste0('(',tumorTissueComparison,' vs. ',nlTissueComparison,')') )
  }
  
  # Set normal comparison
  if(nlTissueComparison == "Blood Derived Normal"){
    nlTissueRows <- (qcMLMetadata$sample_type == "Blood Derived Normal")
  }
  else if(nlTissueComparison == "Solid Tissue Normal"){
    nlTissueRows <- (qcMLMetadata$sample_type == "Solid Tissue Normal")
  }
  else if(nlTissueComparison == "Primary Tumor"){ # Permits comparing metastatic and recurrent tumors to primary tumors
    nlTissueRows <- (qcMLMetadata$sample_type == "Primary Tumor")
  }
  
  # Set tumor tissue comparison
  if(tumorTissueComparison == "Primary Tumor"){
    tumorTissueRows <- (qcMLMetadata$sample_type == "Primary Tumor")
  }
  else if(tumorTissueComparison == "Metastatic"){
    tumorTissueRows <- (qcMLMetadata$sample_type == "Metastatic")
  }
  else if(tumorTissueComparison == "Recurrent Tumor"){
    tumorTissueRows <- (qcMLMetadata$sample_type == "Recurrent Tumor")
  }
  else if(tumorTissueComparison == "Primary Blood Derived Cancer - Peripheral Blood"){
    tumorTissueRows <- (qcMLMetadata$sample_type == "Primary Blood Derived Cancer - Peripheral Blood")
  }
  else if(tumorTissueComparison == "Additional - New Primary"){
    tumorTissueRows <- (qcMLMetadata$sample_type == "Additional - New Primary")
  }

  
  extractedTumorVsNormalRows <- (nlTissueRows | tumorTissueRows)
  
  cancerMetadata <- droplevels(qcMLMetadata[(extractedTumorVsNormalRows & cancerTypeRowsTotal),])
  
  ###########
  
  mlDataY <- cancerMetadata
  mlDataX <- qcMLDataSNM[rownames(mlDataY),]
  dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
  
  # Examine imbalances
  print("The class imbalances are given below...")
  print(table(mlDataY$disease_type))
  print(table(mlDataY$sample_type))
  
  numNl <- as.character(table(mlDataY$sample_type)[1])
  numTumor <- as.character(table(mlDataY$sample_type)[2])
  
  # Create data partitions
  set.seed(42)
  index <- createDataPartition(mlDataY$sample_type, p = trainSetProp, list = FALSE)
  trainX <- mlDataX[index,]
  trainY <- mlDataY[index,]$sample_type
  testX <- mlDataX[-index,]
  testY <- mlDataY[-index,]$sample_type
  
  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
  refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
  
  refactoredTrainY <- relevel(refactoredTrainY, ref = gsub('([[:punct:]])|\\s+','',nlTissueComparison))
  refactoredTestY <- relevel(refactoredTestY, ref = gsub('([[:punct:]])|\\s+','',nlTissueComparison))
  
  test_roc <- function(model, data, classes) {
    
    roc(classes,
        predict(model, data, type = "prob")[, gsub('([[:punct:]])|\\s+','',tumorTissueComparison)])
    
  }

  if(caretModel == "glmnet"){
    set.seed(42)
    ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
    sprintf("Now training model with %s sampling...", samplingStrategy)
    
    # Build up-sampled model
    ctrl$sampling <- samplingStrategy
    
    print("Now training model with up sampling...")
    mlModel <- train(x = trainX,
                   y = refactoredTrainY,
                   method = caretModel,
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   metric = "ROC",
                   tuneGrid = caretTuneGrid)
    
    mlModel %>%
      test_roc(data = testX, classes = refactoredTestY) %>%
      auc() %>% print()
    }
  else{
    set.seed(42)
    ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       verboseIter = TRUE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
    sprintf("Now training model with %s sampling...", samplingStrategy)
  
    # Build up-sampled model
    ctrl$sampling <- samplingStrategy
    
    print("Now training model with up sampling...")
    mlModel <- train(x = trainX,
                   y = refactoredTrainY,
                   method = caretModel,
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   verbose = TRUE,
                   metric = "ROC",
                   tuneGrid = caretTuneGrid)
    
    mlModel %>%
      test_roc(data = testX, classes = refactoredTestY) %>%
      auc() %>% print()
  }
  
  
  
  # # Build SMOTE model
  # smotest <- list(name = "SMOTE with more neighbors!",
  #               func = function (x, y) {
  #                 library(DMwR)
  #                 dat <- if (is.data.frame(x)) x else as.data.frame(x)
  #                 dat$.y <- y
  #                 dat <- SMOTE(.y ~ ., data = dat, perc.over = 200, k = 5, perc.under = 200)
  #                 list(x = dat[, !grepl(".y", colnames(dat), fixed = TRUE)], 
  #                      y = dat$.y)
  #                 },
  #               first = TRUE)

  
  # Examine results for test set
  positiveClass <- gsub('([[:punct:]])|\\s+','',tumorTissueComparison)
  confusionMatrix <- confusionMatrix(predict(mlModel, newdata = testX, type = "raw"), 
                                    refactoredTestY, 
                                    positive = positiveClass)

  write.csv(as.matrix(confusionMatrix), file = confusionMatrixFilename)

  print("Here is the confusion matrix using the test set...")
  print(confusionMatrix)
  
  model_list <- list(mlModel = mlModel)
  
  model_list_roc <- model_list %>%
    map(test_roc, data = testX, classes = refactoredTestY)
  
  model_list_roc %>%
    map(auc) %>% print()
  
  model_list_roc %>%
    map(auc) %>% array() %>% unlist() %>% max() -> maxROC
  
  model_list_roc %>%
    map(auc) %>% array() %>% unlist() %>% which.max() -> nameIndex
  
  model_list_roc %>%
    map(auc) %>% list() -> rocList
  
  bestSamplingStrategy <- names(rocList[[1]])[nameIndex]
  
  rocInsetText <- paste("AUC:\n",sprintf("%1.4f",maxROC),"\n",paste0('(',samplingStrategy,')'))
  
  results_list_roc <- list()
  num_mod <- 1
  
  for(the_roc in model_list_roc){
    
    results_list_roc[[num_mod]] <- 
      data.frame(tpr = the_roc$sensitivities,
                 fpr = 1 - the_roc$specificities,
                 model = names(model_list)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  results_df_roc <- do.call("rbind",results_list_roc)
  
  custom_col <- c("#000000", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#65ABBC")
  if(flagForAlternativeNormal == TRUE){
    rocTitle <- paste(paste("ROC curve",
                            paste0("(Learner: ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")\n")),
                      cancerTypeString,"|",sprintf("Trained on %d%% of data | Tested on %d%% of data", trainSetProp*100, (100-trainSetProp*100)),"\n",
                      paste0('(',numTumor,' ',tumorTissueComparison,' vs. ',numNl,' ',nlTissueComparison,' from ',cancerTypeStringAlternativeNormal,')') )
    rocPlotFileTitle <- paste(paste("ROC curve",
                                    paste0("(Learner ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")")),
                              cancerTypeString,
                              paste0('(',numTumor,' ',tumorTissueComparison,' vs. ',numNl,' ',nlTissueComparison,' from ',cancerTypeStringAlternativeNormal,').png') )
  }
  else{
    rocTitle <- paste(paste("ROC curve",
                            paste0("(Learner: ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")\n")),
                      cancerTypeString,"|",sprintf("Trained on %d%% of data | Tested on %d%% of data", trainSetProp*100, (100-trainSetProp*100)),"\n",
                      paste0('(',numTumor,' ',tumorTissueComparison,' vs. ',numNl,' ',nlTissueComparison,')') )
    rocPlotFileTitle <- paste(paste("ROC curve",
                                    paste0("(Learner ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")")),
                              cancerTypeString,
                              paste0('(',numTumor,' ',tumorTissueComparison,' vs. ',numNl,' ',nlTissueComparison,').png') )
  }
  
  g <- ggplot(aes(x = fpr,  y = tpr, group = model), data = results_df_roc) +
    geom_path(aes(color = model), size = 1) +
    scale_color_manual(values = custom_col) +
    geom_abline(intercept = 0, slope = 1, color = "gray", size = 1) +
    # theme_bw(base_size = 18) +
    coord_equal(ratio=1) + xlim(0, 1) + ylim(0,1) +
    labs(x = "False Positive Rate", y = "True Positive Rate", title = rocTitle) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), legend.position="none") +
    annotate("text", x = 0.75, y = 0.25, label = rocInsetText, size = 4)
  
  ggsave(plot = g, 
         filename = rocPlotFileTitle,
         path = plotPath,
         device = "png",
         width = 16.2,
         height = 5.29,
         units = "in",
         dpi = "retina")
  

  varImpBestModelDF <- as.data.frame(varImp( model_list[[nameIndex]]$finalModel, scale = FALSE ))
  varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
  varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
  colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
  varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]

  write.csv(varImpBestModelDF2OrderedNonzero, file = modelFeaturesFilename, row.names = FALSE)
  
  print(sprintf("Number of non-zero features used by the tuned model: %d", dim(varImpBestModelDF2OrderedNonzero)[1]))

  sink()
  
  print(g)
  print("Here is the confusion matrix using the test set...")
  print(confusionMatrix)

  sink(file = confusionMatrixStatsFilename)
  print("Here is the confusion matrix using the test set...")
  print(confusionMatrix)
  sink()
  
  results <- list(
    # trainX = trainX,
    # trainY = refactoredTrainY,
    testX = testX,
    testY = refactoredTestY,
    rocPlot = g,
    rocPlotData = results_df_roc,
    rocModelData = model_list_roc,
    maxROC = maxROC,
    modelList = model_list,
    modelNonzeroVariableImp = varImpBestModelDF2OrderedNonzero,
    modelAllVariableImp = varImpBestModelDF2Ordered,
    confusionMatrix = confusionMatrix)
  
  return(results)
}