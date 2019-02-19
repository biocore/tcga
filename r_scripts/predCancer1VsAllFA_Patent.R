# predCancer1VsAll.R
predCancer1VsAllFA_Patent <- function(qcMLMetadata, 
                                 qcMLDataSNM,
                                 cancerTypeString = NULL,
                                 sampleTypeComparison = "Primary Tumor",
                                 caretModel = "gbm",
                                 samplingStrategy = "up",
                                 numResampleIter = 2,
                                 numKFold = 3,
                                 trainSetProp = 0.7,
                                 caretTuneGrid = NULL,
                                 plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all",
                                 modelFeaturesPath = "./provisional-patent-csv-model-features-cancer-type-1-vs-all",
                                 functionOutputPath = "./provisional-patent-function-output-cancer-type-1-vs-all",
                                 confusionMatrixPath = "./provisional-patent-confusion-matrices-cancer-type-1-vs-all",
                                 confusionMatrixStatsPath = "./provisional-patent-stats-confusion-matrices-cancer-type-1-vs-all",
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

  sinkFilename <- sprintf("%s/%s - %s (CV k-fold of %d|Train proportion of %d) Sink.txt", 
    functionOutputPath, cancerTypeString, sampleTypeComparison, numKFold, trainSetProp*100)
  modelFeaturesFilename <- sprintf("%s/%s - %s (CV k-fold of %d|Train proportion of %d) Model Features.csv", 
    modelFeaturesPath, cancerTypeString, sampleTypeComparison, numKFold, trainSetProp*100)
  confusionMatrixFilename <- sprintf("%s/%s - %s (CV k-fold of %d|Train proportion of %d) Confusion Matrix.csv", 
    confusionMatrixPath, cancerTypeString, sampleTypeComparison, numKFold, trainSetProp*100)
  confusionMatrixStatsFilename <- sprintf("%s/%s - %s (CV k-fold of %d|Train proportion of %d) Confusion Matrix Stats.txt", 
    confusionMatrixStatsPath, cancerTypeString, sampleTypeComparison, numKFold, trainSetProp*100)
  
  sink(file = sinkFilename)

  
  extractRowsSampleTypeRows <- qcMLMetadata$sample_type == sampleTypeComparison
  
  cancerTypeComparison <- qcMLMetadata$disease_type
  cancerTypeComparisonFactor <- factor(ifelse(cancerTypeComparison == cancerTypeString, yes = cancerTypeString, no = "OtherCancerType"),
                                       levels = c(cancerTypeString, "OtherCancerType"))
  qcMLMetadata$cancerTypeComparison <- cancerTypeComparisonFactor
  
  cancerMetadata <- droplevels(qcMLMetadata[extractRowsSampleTypeRows,])
  
  mlDataY <- cancerMetadata
  mlDataX <- qcMLDataSNM[rownames(mlDataY),]
  dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
  
  # Examine available samples
  print(sprintf("Chosen cancer type of interest: %s", cancerTypeString))
  print(sprintf("Chosen sample type of comparison: %s", sampleTypeComparison))
  print("Breakdown of samples:")
  print(table(mlDataY$cancerTypeComparison))
  
  numCancerType <- as.character(table(mlDataY$cancerTypeComparison)[1])
  numOther <- as.character(table(mlDataY$cancerTypeComparison)[2])
  
  # Create data partitions
  set.seed(42)
  index <- createDataPartition(mlDataY$cancerTypeComparison, p = trainSetProp, list = FALSE)
  trainX <- mlDataX[index,]
  trainY <- mlDataY[index,]$cancerTypeComparison
  testX <- mlDataX[-index,]
  testY <- mlDataY[-index,]$cancerTypeComparison

  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
  refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))

  refactoredTrainY <- relevel(refactoredTrainY, ref = gsub('([[:punct:]])|\\s+','',cancerTypeString))
  refactoredTestY <- relevel(refactoredTestY, ref = gsub('([[:punct:]])|\\s+','',cancerTypeString))

  test_roc <- function(model, data, classes) {

    roc(classes,
        predict(model, data, type = "prob")[, gsub('([[:punct:]])|\\s+','',cancerTypeString)])

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

  # Examine results for test set
  positiveClass <- gsub('([[:punct:]])|\\s+','',cancerTypeString)
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
  rocTitle <- paste(paste("ROC curves",
                          paste0("(Learner: ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")\n")),
                    cancerTypeString,"vs all other cancer types |",sprintf("Trained on %d%% of data | Tested on %d%% of data", trainSetProp*100, (100-trainSetProp*100)),"\n",
                    paste0('(',numCancerType,' ',cancerTypeString,' vs. ',numOther,' other cancer types (All ',sampleTypeComparison,')',')') )
  rocPlotFileTitle <- paste(cancerTypeString,paste("vs all other cancer types - ROC curves",
                          paste0("(Learner ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")")),
                          paste0('(',numCancerType,' vs. ',numOther,' other cancer types (All ',sampleTypeComparison,')',').png') )

  
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
  # file.show(sinkFilename)

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
    modelList = model_list,
    maxROC = maxROC,
    modelNonzeroVariableImp = varImpBestModelDF2OrderedNonzero,
    modelAllVariableImp = varImpBestModelDF2Ordered,
    confusionMatrix = confusionMatrix)

  return(results)
}