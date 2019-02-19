## Goal: Predict tumor vs normal
validationPredTumorVsNormalFA <- function(qcMLMetadataS1,
                                qcMLMetadataS2,
                                qcMLDataSNMS1,
                                qcMLDataSNMS2,
                                cancerTypeString = NULL,
                                nlTissueComparison = "Solid Tissue Normal",
                                tumorTissueComparison = "Primary Tumor", 
                                flagForAlternativeNormal = FALSE,
                                cancerTypeStringAlternativeNormal = NULL,
                                caretModel = "gbm",
                                samplingStrategy = "up",
                                numResampleIter = 2,
                                numKFold = 3,
                                caretTuneGrid = NULL,
                                plotPath = "./roc-ggplots-validation-predicting-tumor-vs-normal",
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
  cancerTypeRowsS1 <- (qcMLMetadataS1$disease_type == cancerTypeString)
  if(flagForAlternativeNormal == TRUE){
    cancerTypeAlternativeNormalRowsS1 <- ((qcMLMetadataS1$disease_type == cancerTypeStringAlternativeNormal) & (qcMLMetadataS1$sample_type == nlTissueComparison))
    cancerTypeRowsTotalS1 <- (cancerTypeRowsS1 | cancerTypeAlternativeNormalRowsS1)
    # rocTitle <- paste(paste("ROC curves",
    #                         paste0("(Learner: ",caretModel,"/","rep=",as.character(numResampleIter),"/","kFold=",as.character(numKFold),")\n")),
    #                   cancerTypeString,"\n",
    #                   paste0('(',tumorTissueComparison,' vs. ',nlTissueComparison,' from ',cancerTypeStringAlternativeNormal,')') )
  }
  else{
    cancerTypeRowsTotalS1 <- cancerTypeRowsS1
    # rocTitle <- paste(paste("ROC curves",
    #                         paste0("(Learner: ",caretModel,"/","rep=",as.character(numResampleIter),"/","kFold=",as.character(numKFold),")\n")),
    #                   cancerTypeString,"\n",
    #                   paste0('(',tumorTissueComparison,' vs. ',nlTissueComparison,')') )
  }
  
  # Set normal comparison
  if(nlTissueComparison == "Blood Derived Normal"){
    nlTissueRowsS1 <- (qcMLMetadataS1$sample_type == "Blood Derived Normal")
  }
  else if(nlTissueComparison == "Solid Tissue Normal"){
    nlTissueRowsS1 <- (qcMLMetadataS1$sample_type == "Solid Tissue Normal")
  }
  else if(nlTissueComparison == "Primary Tumor"){ # Permits comparing metastatic and recurrent tumors to primary tumors
    nlTissueRowsS1 <- (qcMLMetadataS1$sample_type == "Primary Tumor")
  }
  
  # Set tumor tissue comparison
  if(tumorTissueComparison == "Primary Tumor"){
    tumorTissueRowsS1 <- (qcMLMetadataS1$sample_type == "Primary Tumor")
  }
  else if(tumorTissueComparison == "Metastatic"){
    tumorTissueRowsS1 <- (qcMLMetadataS1$sample_type == "Metastatic")
  }
  else if(tumorTissueComparison == "Recurrent Tumor"){
    tumorTissueRowsS1 <- (qcMLMetadataS1$sample_type == "Recurrent Tumor")
  }
  else if(tumorTissueComparison == "Primary Blood Derived Cancer - Peripheral Blood"){
    tumorTissueRowsS1 <- (qcMLMetadataS1$sample_type == "Primary Blood Derived Cancer - Peripheral Blood")
  }
  else if(tumorTissueComparison == "Additional - New Primary"){
    tumorTissueRowsS1 <- (qcMLMetadataS1$sample_type == "Additional - New Primary")
  }

  
  extractedTumorVsNormalRowsS1 <- (nlTissueRowsS1 | tumorTissueRowsS1)
  
  cancerMetadataS1 <- droplevels(qcMLMetadataS1[(extractedTumorVsNormalRowsS1 & cancerTypeRowsTotalS1),])

  #---------- Split 2 -------------

  cancerTypeRowsS2 <- (qcMLMetadataS2$disease_type == cancerTypeString)
  if(flagForAlternativeNormal == TRUE){
    cancerTypeAlternativeNormalRowsS2 <- ((qcMLMetadataS2$disease_type == cancerTypeStringAlternativeNormal) & (qcMLMetadataS2$sample_type == nlTissueComparison))
    cancerTypeRowsTotalS2 <- (cancerTypeRowsS2 | cancerTypeAlternativeNormalRowsS2)
    # rocTitle <- paste(paste("ROC curves",
    #                         paste0("(Learner: ",caretModel,"/","rep=",as.character(numResampleIter),"/","kFold=",as.character(numKFold),")\n")),
    #                   cancerTypeString,"\n",
    #                   paste0('(',tumorTissueComparison,' vs. ',nlTissueComparison,' from ',cancerTypeStringAlternativeNormal,')') )
  }
  else{
    cancerTypeRowsTotalS2 <- cancerTypeRowsS2
    # rocTitle <- paste(paste("ROC curves",
    #                         paste0("(Learner: ",caretModel,"/","rep=",as.character(numResampleIter),"/","kFold=",as.character(numKFold),")\n")),
    #                   cancerTypeString,"\n",
    #                   paste0('(',tumorTissueComparison,' vs. ',nlTissueComparison,')') )
  }
  
  # Set normal comparison
  if(nlTissueComparison == "Blood Derived Normal"){
    nlTissueRowsS2 <- (qcMLMetadataS2$sample_type == "Blood Derived Normal")
  }
  else if(nlTissueComparison == "Solid Tissue Normal"){
    nlTissueRowsS2 <- (qcMLMetadataS2$sample_type == "Solid Tissue Normal")
  }
  else if(nlTissueComparison == "Primary Tumor"){ # Permits comparing metastatic and recurrent tumors to primary tumors
    nlTissueRowsS2 <- (qcMLMetadataS2$sample_type == "Primary Tumor")
  }
  
  # Set tumor tissue comparison
  if(tumorTissueComparison == "Primary Tumor"){
    tumorTissueRowsS2 <- (qcMLMetadataS2$sample_type == "Primary Tumor")
  }
  else if(tumorTissueComparison == "Metastatic"){
    tumorTissueRowsS2 <- (qcMLMetadataS2$sample_type == "Metastatic")
  }
  else if(tumorTissueComparison == "Recurrent Tumor"){
    tumorTissueRowsS2 <- (qcMLMetadataS2$sample_type == "Recurrent Tumor")
  }
  else if(tumorTissueComparison == "Primary Blood Derived Cancer - Peripheral Blood"){
    tumorTissueRowsS2 <- (qcMLMetadataS2$sample_type == "Primary Blood Derived Cancer - Peripheral Blood")
  }
  else if(tumorTissueComparison == "Additional - New Primary"){
    tumorTissueRowsS2 <- (qcMLMetadataS2$sample_type == "Additional - New Primary")
  }

  
  extractedTumorVsNormalRowsS2 <- (nlTissueRowsS2 | tumorTissueRowsS2)
  
  cancerMetadataS2 <- droplevels(qcMLMetadataS2[(extractedTumorVsNormalRowsS2 & cancerTypeRowsTotalS2),])
  
  ###########
  
  mlDataYS1 <- cancerMetadataS1
  mlDataXS1 <- qcMLDataSNMS1[rownames(mlDataYS1),]
  dim(mlDataYS1)[1] == dim(mlDataXS1)[1] # Sanity check

  mlDataYS2 <- cancerMetadataS2
  mlDataXS2 <- qcMLDataSNMS2[rownames(mlDataYS2),]
  dim(mlDataYS2)[1] == dim(mlDataXS2)[1] # Sanity check
  
  # Examine imbalances
  print("The class imbalances are given below...")
  print("For split 1...")
  print(table(mlDataYS1$disease_type))
  print(table(mlDataYS1$sample_type))

  print("For split 2...")
  print(table(mlDataYS2$disease_type))
  print(table(mlDataYS2$sample_type))
  
  numNlS1 <- as.character(table(mlDataYS1$sample_type)[1])
  numTumorS1 <- as.character(table(mlDataYS1$sample_type)[2])

  numNlS2 <- as.character(table(mlDataYS2$sample_type)[1])
  numTumorS2 <- as.character(table(mlDataYS2$sample_type)[2])
  
  # Create data partitions
  set.seed(42)
  # index <- createDataPartition(mlDataY$sample_type, p = trainSetProp, list = FALSE)

  trainXS1 <- mlDataXS1
  trainYS1 <- mlDataYS1$sample_type
  testXS1 <- mlDataXS2
  testYS1 <- mlDataYS2$sample_type

  trainXS2 <- mlDataXS2
  trainYS2 <- mlDataYS2$sample_type
  testXS2 <- mlDataXS1
  testYS2 <- mlDataYS1$sample_type
  
  refactoredTrainYS1 <- factor(gsub('([[:punct:]])|\\s+','',trainYS1))
  refactoredTestYS1 <- factor(gsub('([[:punct:]])|\\s+','',testYS1))

  refactoredTrainYS2 <- factor(gsub('([[:punct:]])|\\s+','',trainYS2))
  refactoredTestYS2 <- factor(gsub('([[:punct:]])|\\s+','',testYS2))
  
  refactoredTrainYS1 <- relevel(refactoredTrainYS1, ref = gsub('([[:punct:]])|\\s+','',nlTissueComparison))
  refactoredTestYS1 <- relevel(refactoredTestYS1, ref = gsub('([[:punct:]])|\\s+','',nlTissueComparison))

  refactoredTrainYS2 <- relevel(refactoredTrainYS2, ref = gsub('([[:punct:]])|\\s+','',nlTissueComparison))
  refactoredTestYS2 <- relevel(refactoredTestYS2, ref = gsub('([[:punct:]])|\\s+','',nlTissueComparison))
  
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
  positiveClass <- gsub('([[:punct:]])|\\s+','',tumorTissueComparison)
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
  if(flagForAlternativeNormal == TRUE){
    rocTitle <- "Test"
    rocPlotFileTitle <- "Test"
    # rocTitle <- paste(paste("ROC curve",
    #                         paste0("(Learner: ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")\n")),
    #                   cancerTypeString,"|",sprintf("Trained on %d%% of data | Tested on %d%% of data", trainSetProp*100, (100-trainSetProp*100)),"\n",
    #                   paste0('(',numTumor,' ',tumorTissueComparison,' vs. ',numNl,' ',nlTissueComparison,' from ',cancerTypeStringAlternativeNormal,')') )
    # rocPlotFileTitle <- paste(paste("ROC curve",
    #                                 paste0("(Learner ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")")),
    #                           cancerTypeString,
    #                           paste0('(',numTumor,' ',tumorTissueComparison,' vs. ',numNl,' ',nlTissueComparison,' from ',cancerTypeStringAlternativeNormal,').png') )
  }
  else{
    rocTitle <- paste(paste("ROC curve",
                            paste0("(Learner: ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")\n")),
                      cancerTypeString,"\n","Comparing models built on a split pipeline for validation (50% stratified splits)","\n",
                      paste0('(',numTumorS1,'|',numTumorS2,' ',tumorTissueComparison,' vs. ',numNlS1,'|',numNlS2,' ',nlTissueComparison,')') )
    rocPlotFileTitle <- paste(paste(cancerTypeString,
                                    "- ROC curve",
                                    paste0("(Learner ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")")),
                              paste0('(',numTumorS1,'|',numTumorS2,' ',tumorTissueComparison,' vs. ',numNlS1,'|',numNlS2,' ',nlTissueComparison,').png') )
  }
  
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
  # save(model_list, file  = paste0(caretModel,"_model_list.RData"))
  
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