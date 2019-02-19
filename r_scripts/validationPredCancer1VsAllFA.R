# validationPredCancer1VsAll.R
validationPredCancer1VsAllFA <- function(qcMLMetadataS1,
                                qcMLMetadataS2,
                                qcMLDataSNMS1,
                                qcMLDataSNMS2,
                                cancerTypeString = NULL,
                                 sampleTypeComparison = "Primary Tumor",
                                 caretModel = "gbm",
                                 samplingStrategy = "up",
                                 numResampleIter = 2,
                                 numKFold = 3,
                                 trainSetProp = 0.7,
                                 caretTuneGrid = NULL,
                                 plotPath = "./roc-ggplots-predicting-cancer-type-1-vs-all",
                                 originalModel = NULL,
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
  require(ggsci) # for plot colors
  
  numCores <- detectCores()
  registerDoMC(cores=numCores)

  #---------- Split 1 -------------
  
  extractRowsSampleTypeRowsS1 <- qcMLMetadataS1$sample_type == sampleTypeComparison
  
  cancerTypeComparisonS1 <- qcMLMetadataS1$disease_type
  cancerTypeComparisonFactorS1 <- factor(ifelse(cancerTypeComparisonS1 == cancerTypeString, yes = cancerTypeString, no = "OtherCancerType"),
                                       levels = c(cancerTypeString, "OtherCancerType"))
  qcMLMetadataS1$cancerTypeComparison <- cancerTypeComparisonFactorS1
  
  cancerMetadataS1 <- droplevels(qcMLMetadataS1[extractRowsSampleTypeRowsS1,])

  #---------- Split 2 -------------

  extractRowsSampleTypeRowsS2 <- qcMLMetadataS2$sample_type == sampleTypeComparison
  
  cancerTypeComparisonS2 <- qcMLMetadataS2$disease_type
  cancerTypeComparisonFactorS2 <- factor(ifelse(cancerTypeComparisonS2 == cancerTypeString, yes = cancerTypeString, no = "OtherCancerType"),
                                       levels = c(cancerTypeString, "OtherCancerType"))
  qcMLMetadataS2$cancerTypeComparison <- cancerTypeComparisonFactorS2
  
  cancerMetadataS2 <- droplevels(qcMLMetadataS2[extractRowsSampleTypeRowsS2,])

  ###########
  
  mlDataYS1 <- cancerMetadataS1
  mlDataXS1 <- qcMLDataSNMS1[rownames(mlDataYS1),]
  dim(mlDataYS1)[1] == dim(mlDataXS1)[1] # Sanity check

  mlDataYS2 <- cancerMetadataS2
  mlDataXS2 <- qcMLDataSNMS2[rownames(mlDataYS2),]
  dim(mlDataYS2)[1] == dim(mlDataXS2)[1] # Sanity check
  
  # Examine available samples
  print(sprintf("Chosen cancer type of interest: %s", cancerTypeString))
  print(sprintf("Chosen sample type of comparison: %s", sampleTypeComparison))
  print("Breakdown of samples:")
  print("For split 1...")
  print(table(mlDataYS1$cancerTypeComparison))
  print("For split 2...")
  print(table(mlDataYS2$cancerTypeComparison))
  
  numCancerTypeS1 <- as.character(table(mlDataYS1$cancerTypeComparison)[1])
  numCancerTypeS2 <- as.character(table(mlDataYS2$cancerTypeComparison)[1])
  numOtherS1 <- as.character(table(mlDataYS1$cancerTypeComparison)[2])
  numOtherS2 <- as.character(table(mlDataYS2$cancerTypeComparison)[2])
  
  # Create data partitions
  set.seed(42)
  # index <- createDataPartition(mlDataY$cancerTypeComparison, p = trainSetProp, list = FALSE)
  trainXS1 <- mlDataXS1
  trainYS1 <- mlDataYS1$cancerTypeComparison
  testXS1 <- mlDataXS2
  testYS1 <- mlDataYS2$cancerTypeComparison

  trainXS2 <- mlDataXS2
  trainYS2 <- mlDataYS2$cancerTypeComparison
  testXS2 <- mlDataXS1
  testYS2 <- mlDataYS1$cancerTypeComparison
  
  refactoredTrainYS1 <- factor(gsub('([[:punct:]])|\\s+','',trainYS1))
  refactoredTestYS1 <- factor(gsub('([[:punct:]])|\\s+','',testYS1))

  refactoredTrainYS2 <- factor(gsub('([[:punct:]])|\\s+','',trainYS2))
  refactoredTestYS2 <- factor(gsub('([[:punct:]])|\\s+','',testYS2))
  
  refactoredTrainYS1 <- relevel(refactoredTrainYS1, ref = gsub('([[:punct:]])|\\s+','',cancerTypeString))
  refactoredTestYS1 <- relevel(refactoredTestYS1, ref = gsub('([[:punct:]])|\\s+','',cancerTypeString))

  refactoredTrainYS2 <- relevel(refactoredTrainYS2, ref = gsub('([[:punct:]])|\\s+','',cancerTypeString))
  refactoredTestYS2 <- relevel(refactoredTestYS2, ref = gsub('([[:punct:]])|\\s+','',cancerTypeString))

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
    
    # Build up-sampled model
    ctrl$sampling <- samplingStrategy
    
    print(sprintf("Now training model 1 with %s sampling...", samplingStrategy))
    mlModelS1 <- train(x = trainXS1,
                   y = refactoredTrainYS1,
                   method = caretModel,
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   metric = "ROC",
                   tuneGrid = caretTuneGrid)

    mlModelS1 %>%
      test_roc(data = testXS1, classes = refactoredTestYS1) %>%
      auc() %>% print()

    print(sprintf("Now training model 2 with %s sampling...", samplingStrategy))
    mlModelS2 <- train(x = trainXS2,
                   y = refactoredTrainYS2,
                   method = caretModel,
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   metric = "ROC",
                   tuneGrid = caretTuneGrid)
    
    mlModelS2 %>%
      test_roc(data = testXS2, classes = refactoredTestYS2) %>%
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
  
    # Build up-sampled model
    ctrl$sampling <- samplingStrategy
    
    print(sprintf("Now training model 1 with %s sampling...", samplingStrategy))
    mlModelS1 <- train(x = trainXS1,
                   y = refactoredTrainYS1,
                   method = caretModel,
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   verbose = TRUE,
                   metric = "ROC",
                   tuneGrid = caretTuneGrid)
    
    mlModelS1 %>%
      test_roc(data = testXS1, classes = refactoredTestYS1) %>%
      auc() %>% print()

    print(sprintf("Now training model 2 with %s sampling...", samplingStrategy))
    mlModelS2 <- train(x = trainXS2,
                   y = refactoredTrainYS2,
                   method = caretModel,
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   verbose = TRUE,
                   metric = "ROC",
                   tuneGrid = caretTuneGrid)
    
    mlModelS2 %>%
      test_roc(data = testXS2, classes = refactoredTestYS2) %>%
      auc() %>% print()
  }

  
  # Examine results for test set
  positiveClass <- gsub('([[:punct:]])|\\s+','',cancerTypeString)
  confusionMatrixS1 <- confusionMatrix(predict(mlModelS1, newdata = testXS1, type = "raw"), 
                                    refactoredTestYS1, 
                                    positive = positiveClass)
  confusionMatrixS2 <- confusionMatrix(predict(mlModelS2, newdata = testXS2, type = "raw"), 
                                    refactoredTestYS2, 
                                    positive = positiveClass)

  confusionMatrices <- list(cmModelS1 = confusionMatrixS1,
                            cmModelS2 = confusionMatrixS2)

  print("Here is the confusion matrix for the model built on data from split 1 and tested on data from split 2...")
  print(confusionMatrixS1)

  print("Here is the confusion matrix for the model built on data from split 2 and tested on data from split 1...")
  print(confusionMatrixS2)
  
  model_list <- list(mlModelS1 = mlModelS1,
                     mlModelS2 = mlModelS2)

  mlModelS1 %>%
      test_roc(data = testXS1, classes = refactoredTestYS1) %>%
      auc() -> mlModelS1_AUC

  mlModelS2 %>%
  test_roc(data = testXS2, classes = refactoredTestYS2) %>%
  auc() -> mlModelS2_AUC
  
  rocInsetText <- paste0("AUC Model Split 1: ",sprintf("%1.3f",mlModelS1_AUC),"\n","AUC Model Split 2: ",sprintf("%1.3f",mlModelS2_AUC))

  mlModelS1 %>%
      test_roc(data = testXS1, classes = refactoredTestYS1) %>%
      list() -> rocMlModelS1

  mlModelS2 %>%
    test_roc(data = testXS2, classes = refactoredTestYS2) %>%
    list() -> rocMlModelS2

  model_list_roc <- c(rocMlModelS1, rocMlModelS2)
  
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

  if(!is.null(originalModel)){
    results_df_roc <- rbind(originalModel$rocPlotData, results_df_roc)
    rocInsetText <- paste0("AUC Original: ",sprintf("%1.3f",originalModel$maxROC),"\n",
                          "AUC Model Split 1: ",sprintf("%1.3f",mlModelS1_AUC),"\n",
                          "AUC Model Split 2: ",sprintf("%1.3f",mlModelS2_AUC))
  }
  
  # custom_col <- c("#000000", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#65ABBC")
  rocTitle <- "Test"
  rocPlotFileTitle <- "Test.png"

  # rocTitle <- paste(paste("ROC curves",
  #                         paste0("(Learner: ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")\n")),
  #                   cancerTypeString,"vs all other cancer types |",sprintf("Trained on %d%% of data | Tested on %d%% of data", trainSetProp*100, (100-trainSetProp*100)),"\n",
  #                   paste0('(',numCancerType,' ',cancerTypeString,' vs. ',numOther,' other cancer types (All ',sampleTypeComparison,')',')') )
  # rocPlotFileTitle <- paste(cancerTypeString,paste("vs all other cancer types - ROC curves",
  #                         paste0("(Learner ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")")),
  #                         paste0('(',numCancerType,' vs. ',numOther,' other cancer types (All ',sampleTypeComparison,')',').png') )

  
  g <- ggplot(aes(x = fpr,  y = tpr, color = model), data = results_df_roc) +
    geom_path(size = 1) +
    # scale_color_manual(values = custom_col) +
    geom_abline(intercept = 0, slope = 1, color = "gray", size = 1) +
    # theme_bw(base_size = 18) +
    coord_equal(ratio=1) + xlim(0, 1) + ylim(0,1) +
    labs(x = "False Positive Rate", y = "True Positive Rate", title = rocTitle) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + #, legend.position="none"
    annotate("text", x = 0.70, y = 0.25, label = rocInsetText, size = 4) +
    scale_color_nejm(name = "Model Type", labels = c("mlModelS1" = "Trained on split 1 | Tested on split 2", 
                                                      "mlModel" = "Trained and tested on original Voom-SNM Data",
                                                        "mlModelS2" = "Trained on split 2 | Tested on split 1"))

  print(g)

  if(!( dir.exists( file.path(plotPath)))){
    dir.create(file.path(plotPath))
  }
  
  ggsave(plot = g, 
         filename = rocPlotFileTitle,
         path = plotPath,
         device = "png",
         width = 16.2,
         height = 5.29,
         units = "in",
         dpi = "retina")
  

  varImpBestModelDFS1 <- as.data.frame(varImp( model_list[[1]]$finalModel, scale = FALSE ))
  varImpBestModelDF2S1 <- rownames_to_column(varImpBestModelDFS1, "Taxa")
  varImpBestModelDF2OrderedS1 <- varImpBestModelDF2S1[order(-varImpBestModelDF2S1$Overall),]
  colnames(varImpBestModelDF2OrderedS1)[2] <- "varImp"
  varImpBestModelDF2OrderedNonzeroS1 <- varImpBestModelDF2OrderedS1[varImpBestModelDF2OrderedS1$varImp != 0,]
  
  print(sprintf("Number of non-zero features used by the tuned model 1: %d", dim(varImpBestModelDF2OrderedNonzeroS1)[1]))

  varImpBestModelDFS2 <- as.data.frame(varImp( model_list[[1]]$finalModel, scale = FALSE ))
  varImpBestModelDF2S2 <- rownames_to_column(varImpBestModelDFS2, "Taxa")
  varImpBestModelDF2OrderedS2 <- varImpBestModelDF2S2[order(-varImpBestModelDF2S2$Overall),]
  colnames(varImpBestModelDF2OrderedS2)[2] <- "varImp"
  varImpBestModelDF2OrderedNonzeroS2 <- varImpBestModelDF2OrderedS2[varImpBestModelDF2OrderedS2$varImp != 0,]
  
  print(sprintf("Number of non-zero features used by the tuned model 2: %d", dim(varImpBestModelDF2OrderedNonzeroS2)[1]))

  intersectingModelFeatures <- intersect(varImpBestModelDF2OrderedNonzeroS1$Taxa, varImpBestModelDF2OrderedNonzeroS2$Taxa)
  print(sprintf("Number of intersecting features between tuned model 1 and tuned model 2: %d", length(intersectingModelFeatures)))

  results <- list(
    # trainX = trainX,
    # trainY = refactoredTrainY,
    testXS1 = testXS1,
    testYS1 = refactoredTestYS1,
    testXS2 = testXS2,
    testYS2 = refactoredTestYS2,
    rocPlot = g,
    rocPlotData = results_df_roc,
    rocModelData = model_list_roc,
    modelList = model_list,
    modelNonzeroVariableImpS1 = varImpBestModelDF2OrderedNonzeroS1,
    modelAllVariableImpS1 = varImpBestModelDF2OrderedS1,
    modelNonzeroVariableImpS2 = varImpBestModelDF2OrderedNonzeroS2,
    modelAllVariableImpS2 = varImpBestModelDF2OrderedS2,
    intersectingModelFeatures = intersectingModelFeatures,
    confusionMatrices = confusionMatrices)

  return(results)
}