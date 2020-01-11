# permutation_testing.R
# Greg Poore
# Nov 18 2019
# Purpose: Permute ranks to identify differences for rank-based normalization

## Load dependencies
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
require(gtools) # for permute function
require(limma) # for voom normalization
require(edgeR) # for voom normalization
require(snm) # for snm normalization
require(gmodels) # for confidence interval calc


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

caretTuneGrid <- defaultGBMGrid
numKFold <- 4
numResampleIter <- 1

# metadataPSMatchedQC, vbMergedDataKrakenCFDecontamQC,
# file = "kraken_vbAndMetadata.RData"

load("kraken_vbAndMetadataFiltered_Dec2.RData")

#-------------------------------------------------

kraken_age_permute <- permPerfHvsC_Age(iterPerm = 100,
                                   permuteRankOrNominalAge = "nominal",
                                   metaData = metadataPSMatchedDPQC,
                                   readData = vbMergedDataKrakenDPCFDecontamQC,
                                   permuteFlag = TRUE)

kraken_age_static <- permPerfHvsC_Age(iterPerm = 100,
                                  permuteRankOrNominalAge = "nominal",
                                  metaData = metadataPSMatchedDPQC,
                                  readData = vbMergedDataKrakenDPCFDecontamQC,
                                  permuteFlag = FALSE)

kraken_age_ps_data <- rbind(data.frame(kraken_age_permute, permuted = "TRUE"),
                            data.frame(kraken_age_static, permuted = "FALSE"))
kraken_age_ps_data.melted <- melt(kraken_age_ps_data)

save(kraken_age_ps_data, file = "kraken_age_ps_data_100.RData")

write.csv(kraken_age_ps_data.melted, file = "../AARF/edf8h.csv")

ggboxplot(kraken_age_ps_data.melted,
          x = "permuted",
          y = "value",
          palette = "nejm",
          ylim  = c(0.5,1.02),
          xlab = "Permuted Age",
          ylab = "Performance Value",
          panel.labs = list(variable = c("AUROC", "AUPR")),
          legend = "none",
          fill = "permuted",
          facet.by = "variable") +
  theme(strip.text.x = element_text(size = 14)) +
  stat_compare_means(comparisons = list(c("TRUE","FALSE")), 
                     # label.y = 0.95,
                     label = "p.format") #-> kraken_age_permute_plot

ggsave(plot = kraken_age_permute_plot,
       filename = "kraken_age_permute_plot_100.svg",
       width = 3.5,
       height = 3.5,
       units = "in",
       dpi = "retina")


#-------------------------------------------------

kraken_sex_permute <- permPerfHvsC_Sex(iterPerm = 100,
                                   metaData = metadataPSMatchedDPQC,
                                   readData = vbMergedDataKrakenDPCFDecontamQC,
                                   permuteFlag = TRUE)

kraken_sex_static <- permPerfHvsC_Sex(iterPerm = 100,
                                       metaData = metadataPSMatchedDPQC,
                                       readData = vbMergedDataKrakenDPCFDecontamQC,
                                       permuteFlag = FALSE)

kraken_sex_ps_data <- rbind(data.frame(kraken_sex_permute, permuted = "TRUE"),
                            data.frame(kraken_sex_static, permuted = "FALSE"))
kraken_sex_ps_data.melted <- melt(kraken_sex_ps_data)

save(kraken_sex_ps_data, file = "kraken_sex_ps_data_100.RData")

write.csv(kraken_sex_ps_data.melted, file = "../AARF/edf8i.csv")

ggboxplot(kraken_sex_ps_data.melted,
          x = "permuted",
          y = "value",
          palette = "nejm",
          ylim  = c(0.5,1.02),
          xlab = "Permuted Sex",
          ylab = "Performance Value",
          panel.labs = list(variable = c("AUROC", "AUPR")),
          legend = "none",
          fill = "permuted",
          facet.by = "variable") +
  theme(strip.text.x = element_text(size = 14)) +
  stat_compare_means(comparisons = list(c("TRUE","FALSE")), 
                     # label.y = 0.95,
                     label = "p.format") #-> kraken_sex_permute_plot

ggsave(plot = kraken_sex_permute_plot,
       filename = "kraken_sex_permute_plot_100.svg",
       width = 3.5,
       height = 3.5,
       units = "in",
       dpi = "retina")

#-------------------------------------------------


kraken_ageAndSex_permute <- permPerfHvsC_AgeAndSex(iterPerm = 100,
                                             permuteRankOrNominalAge = "nominal",
                                             metaData = metadataPSMatchedDPQC,
                                             readData = vbMergedDataKrakenDPCFDecontamQC,
                                             permuteFlagAge = TRUE,
                                             permuteFlagSex = TRUE)

kraken_ageAndSex_static <- permPerfHvsC_AgeAndSex(iterPerm = 100,
                                                   permuteRankOrNominalAge = "nominal",
                                                   metaData = metadataPSMatchedDPQC,
                                                   readData = vbMergedDataKrakenDPCFDecontamQC,
                                                   permuteFlagAge = FALSE,
                                                   permuteFlagSex = FALSE)

kraken_ageAndSex_ps_data <- rbind(data.frame(kraken_ageAndSex_permute, permuted = "TRUE"),
                            data.frame(kraken_ageAndSex_static, permuted = "FALSE"))
kraken_ageAndSex_ps_data.melted <- melt(kraken_ageAndSex_ps_data)

save(kraken_ageAndSex_ps_data, file = "kraken_ageAndSex_ps_data_100.RData")

write.csv(kraken_ageAndSex_ps_data.melted, file = "../AARF/edf8j.csv")

ggboxplot(kraken_ageAndSex_ps_data.melted,
          x = "permuted",
          y = "value",
          palette = "nejm",
          ylim  = c(0.5,1.05),
          xlab = "Permuted Age & Sex",
          ylab = "Performance Value",
          panel.labs = list(variable = c("AUROC", "AUPR")),
          legend = "none",
          fill = "permuted",
          facet.by = "variable") +
  theme(strip.text.x = element_text(size = 14)) +
  stat_compare_means(comparisons = list(c("TRUE","FALSE")), 
                     # label.y = 0.95,
                     label = "p.format") #-> kraken_ageAndSex_permute_plot

ggsave(plot = kraken_ageAndSex_permute_plot,
       filename = "kraken_ageAndSex_permute_plot_100.svg",
       width = 3.5,
       height = 3.5,
       units = "in",
       dpi = "retina")


#-------------------------------------------------


tmp <- metadataPSMatchedQC$host_age
tmp2 <- metadataPSMatchedQC$host_age_ordered
ageRank <- match(tmp, sort(unique(tmp)))

tmp3 <- permute(ageRank)

metadataPSMatchedQC$host_age_perm <- ordered(tmp3)

perm100Results <- permPerf(iterPerm = 100)

permPerfHvsC(iterPerm = 5, permuteRankOrNominalAge = "nominal", permuteFlag = FALSE)

load("kraken_vbAndMetadataFiltered_Nov19.RData")
permPerfHvsC_sex(iterPerm = 5, permuteFlag = TRUE)

metadataPSMatchedQC$host_age_perm <- ordered(match(metadataPSMatchedQC$host_age, sort(unique(metadataPSMatchedQC$host_age))))

vsnm()

mlHvsC(snmData = snmDataKrakenCFDecontamQC)

loocvDTs(snmData = snmDataKrakenCFDecontamQC,
         samplingSize = 15, 
         DTs = c("Prostate Cancer","Skin Cutaneous Melanoma","Lung Adenocarcinoma"))

permPerfLOOCV(iterPerm = 3, permuteRankOrNominalAge = "nominal", permuteFlag = TRUE,
              samplingSize = 5, 
              DTs = c("Prostate Cancer","Skin Cutaneous Melanoma","Lung Adenocarcinoma"),
              caretTuneGrid = customGBMGrid)

#------------------------------------------------


## Main function ##
permPerfHvsC_Age <- function(iterPerm = 10, 
                         permuteRankOrNominalAge = "rank", 
                         metaData,
                         readData,
                         permuteFlag = TRUE){
  
  predProbs <- list()
  fg <- list()
  bg <- list()
  perf <- list()
  prroc_roc <- list()
  prroc_pr <- list()
  
  if(tolower(permuteRankOrNominalAge) == "rank"){
    agePrePermuteRank <- match(metaData$host_age, sort(unique(metaData$host_age)))
    agePermList <- t(lapply(1:iterPerm, 
                            function(x) sample(agePrePermuteRank, 
                                               length(agePrePermuteRank), 
                                               replace = FALSE)))
  } else if (tolower(permuteRankOrNominalAge) == "nominal"){
    agePermList <- t(lapply(1:iterPerm, 
                            function(x) sample(metaData$host_age, 
                                               length(metaData$host_age), 
                                               replace = FALSE)))
  } 
  
  # if(permuteFlag == FALSE){
  #   iterPerm <- 1
  # }
  
  qcMetadata <- metaData # metadataPSPLMatchedQC
  
  # set.seed(floor(runif(1)*500))
  
  # if(tolower(permuteRankOrNominalAge) == "rank"){
  #   agePrePermute <- match(metaData$host_age, sort(unique(metaData$host_age)))
  #   qcMetadata$host_age_perm <-  ordered(sample(agePrePermute, length(agePrePermute)))
  # } else if (tolower(permuteRankOrNominalAge) == "nominal"){
  #   agePrePermute <- metaData$host_age
  #   qcMetadata$host_age_perm <- sample(agePrePermute, length(agePrePermute))
  # }
  
  for(ii in 1:iterPerm){
    
    print(sprintf("Iteration: %d/%d", ii, iterPerm))
    
    # qcMetadata <- metaData # metadataPSPLMatchedQC
    
    # set.seed(floor(runif(1)*500))
    
    if(tolower(permuteRankOrNominalAge) == "rank"){
      # agePrePermute <- match(metaData$host_age, sort(unique(metaData$host_age)))
      qcMetadata$host_age_perm <-  ordered(agePermList[[ii]])
    } else if (tolower(permuteRankOrNominalAge) == "nominal"){
      # agePrePermute <- metaData$host_age
      qcMetadata$host_age_perm <- agePermList[[ii]]
    }
    
    if(permuteFlag == FALSE & tolower(permuteRankOrNominalAge) == "nominal"){
      qcMetadata$host_age_perm <- metaData$host_age
    } else if (permuteFlag == FALSE & tolower(permuteRankOrNominalAge) == "rank"){
      qcMetadata$host_age_perm <- ordered(metaData$host_age)
    }
    
    print(qcMetadata$host_age_perm)
    
    qcData <- readData # dataPSPLUniqueQC
    
    # Set up design matrix
    covDesignNorm <- model.matrix(~0 + disease_type +
                                    host_age_perm,
                                  data = qcMetadata)
    
    # Check row dimensions
    dim(covDesignNorm)[1] == dim(qcData)[1]
    
    # print(colnames(covDesignNorm))
    colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
    # print(colnames(covDesignNorm))
    
    # Set up counts matrix
    counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
    
    # Normalize using edgeR and then plug into voom
    dge <- DGEList(counts = counts)
    vdge__perm <- voom(dge, design = covDesignNorm, 
                       plot = TRUE, save.plot = TRUE, 
                       normalize.method="quantile")
    
    # Apply
    bio.var <- model.matrix(~disease_type,
                            data=qcMetadata)
    
    adj.var <- model.matrix(~host_age_perm,
                            data=qcMetadata)
    
    colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
    colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
    
    snmDataObjOnly__perm <- snm(raw.dat = vdge__perm$E, 
                                                             bio.var = bio.var, 
                                                             adj.var = adj.var, 
                                                             rm.adj=TRUE,
                                                             verbose = TRUE,
                                                             diagnose = TRUE)
    snmData__perm <- t(snmDataObjOnly__perm$norm.dat)
    
    metaTmp1 <- metaData
    metaTmp1$disease_type <- metaData$HvsC
    
    mlDataY <- metaTmp1
    mlDataX <- snmData__perm[rownames(mlDataY),]
    
    if(permuteFlag){
      seedVal <- 42
      print(sprintf("Random Number Seed Value: %d",seedVal))
    } else{
      seedVal <- sample(1:1e3, 1)
      print(sprintf("Random Number Seed Value: %d",seedVal))
    }
    
    set.seed(seedVal)
    index <- createDataPartition(mlDataY$disease_type, p = 0.7, list = FALSE)
    trainX <- mlDataX[index,]
    trainY <- mlDataY[index,]$disease_type
    testX <- mlDataX[-index,]
    testY <- mlDataY[-index,]$disease_type
    
    refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
    refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
    
    set.seed(seedVal)
    ctrl <- trainControl(method = "repeatedcv",
                         number = numKFold,
                         repeats = numResampleIter,
                         sampling = "up",
                         summaryFunction = twoClassSummary,
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
                     tuneGrid = defaultGBMGrid)
    
    positiveClass <- "Cancer"
    negativeClass <- "Control"
    predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
    fg[[ii]] <- predProbs[refactoredTestY == positiveClass]
    bg[[ii]] <- predProbs[refactoredTestY == negativeClass]
    
    prroc_roc[[ii]] <- roc.curve(scores.class0 = fg[[ii]], scores.class1 = bg[[ii]], curve = T)
    plot(prroc_roc[[ii]])
    prroc_pr[[ii]] <- pr.curve(scores.class0 = fg[[ii]], scores.class1 = bg[[ii]], curve = T, rand.compute=T)
    plot(prroc_pr[[ii]])
    
    perf[[ii]] <- data.frame(aucroc = prroc_roc[[ii]]$auc,
                             aupr = prroc_pr[[ii]]$auc.integral)
    
    rm(mlModel)
    
  }
  
  perfAll <- do.call(rbind, perf)
  print(perfAll)
  
  require(gmodels)
  print(ci(perfAll$aucroc))
  print(ci(perfAll$aupr))
  
  return(perfAll)

}

#------------------------------------------------


## Main function ##
permPerfHvsC_AgeAndSex <- function(iterPerm = 10, 
                         permuteRankOrNominalAge = "rank", 
                         metaData,
                         readData,
                         permuteFlagAge = TRUE,
                         permuteFlagSex = FALSE){
  
  predProbs <- list()
  fg <- list()
  bg <- list()
  perf <- list()
  prroc_roc <- list()
  prroc_pr <- list()
  
  sexPermList <- t(lapply(1:iterPerm, 
                          function(x) sample(metaData$sex, 
                                             length(metaData$sex), 
                                             replace = FALSE)))
  
  if(tolower(permuteRankOrNominalAge) == "rank"){
    agePrePermuteRank <- match(metaData$host_age, sort(unique(metaData$host_age)))
    agePermList <- t(lapply(1:iterPerm, 
                            function(x) sample(agePrePermuteRank, 
                                               length(agePrePermuteRank), 
                                               replace = FALSE)))
  } else if (tolower(permuteRankOrNominalAge) == "nominal"){
    agePermList <- t(lapply(1:iterPerm, 
                            function(x) sample(metaData$host_age, 
                                               length(metaData$host_age), 
                                               replace = FALSE)))
  } 
  
  # if(permuteFlag == FALSE){
  #   iterPerm <- 1
  # }
  
  qcMetadata <- metaData # metadataPSPLMatchedQC
  
  # set.seed(floor(runif(1)*500))
  
  # if(tolower(permuteRankOrNominalAge) == "rank"){
  #   agePrePermute <- match(metaData$host_age, sort(unique(metaData$host_age)))
  #   qcMetadata$host_age_perm <-  ordered(sample(agePrePermute, length(agePrePermute)))
  # } else if (tolower(permuteRankOrNominalAge) == "nominal"){
  #   agePrePermute <- metaData$host_age
  #   qcMetadata$host_age_perm <- sample(agePrePermute, length(agePrePermute))
  # }
  
  for(ii in 1:iterPerm){
    
    print(sprintf("Iteration: %d/%d", ii, iterPerm))
    
    # qcMetadata <- metaData # metadataPSPLMatchedQC
    
    # set.seed(floor(runif(1)*500))
    
    if(tolower(permuteRankOrNominalAge) == "rank"){
      # agePrePermute <- match(metaData$host_age, sort(unique(metaData$host_age)))
      qcMetadata$host_age_perm <-  ordered(agePermList[[ii]])
    } else if (tolower(permuteRankOrNominalAge) == "nominal"){
      # agePrePermute <- metaData$host_age
      qcMetadata$host_age_perm <- agePermList[[ii]]
    }
    
    if(permuteFlagAge == FALSE & tolower(permuteRankOrNominalAge) == "nominal"){
      qcMetadata$host_age_perm <- metaData$host_age
    } else if (permuteFlagAge == FALSE & tolower(permuteRankOrNominalAge) == "rank"){
      qcMetadata$host_age_perm <- ordered(metaData$host_age)
    }
    
    print(qcMetadata$host_age_perm)
    
    if(permuteFlagSex){
      qcMetadata$sex_perm <- factor(sexPermList[[ii]])
    } else{
      qcMetadata$sex_perm <- factor(metaData$sex)
    }
    
    print(qcMetadata$sex_perm)
    
    qcData <- readData # dataPSPLUniqueQC
    
    # Set up design matrix
    covDesignNorm <- model.matrix(~0 + disease_type +
                                    host_age_perm + 
                                    sex_perm,
                                  data = qcMetadata)
    
    # Check row dimensions
    dim(covDesignNorm)[1] == dim(qcData)[1]
    
    # print(colnames(covDesignNorm))
    colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
    # print(colnames(covDesignNorm))
    
    # Set up counts matrix
    counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
    
    # Normalize using edgeR and then plug into voom
    dge <- DGEList(counts = counts)
    vdge__perm <- voom(dge, design = covDesignNorm, 
                       plot = TRUE, save.plot = TRUE, 
                       normalize.method="quantile")
    
    # Apply
    bio.var <- model.matrix(~disease_type,
                            data=qcMetadata)
    
    adj.var <- model.matrix(~host_age_perm + 
                              sex_perm,
                            data=qcMetadata)
    
    colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
    colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
    
    snmDataObjOnly__perm <- snm(raw.dat = vdge__perm$E, 
                                bio.var = bio.var, 
                                adj.var = adj.var, 
                                rm.adj=TRUE,
                                verbose = TRUE,
                                diagnose = TRUE)
    snmData__perm <- t(snmDataObjOnly__perm$norm.dat)
    
    metaTmp1 <- metaData
    metaTmp1$disease_type <- metaData$HvsC
    
    mlDataY <- metaTmp1
    mlDataX <- snmData__perm[rownames(mlDataY),]
    
    if(permuteFlagAge | permuteFlagSex){
      seedVal <- 42
      print(sprintf("Random Number Seed Value: %d",seedVal))
    } else{
      seedVal <- sample(1:1e3, 1)
      print(sprintf("Random Number Seed Value: %d",seedVal))
    }
    
    set.seed(seedVal)
    index <- createDataPartition(mlDataY$disease_type, p = 0.7, list = FALSE)
    trainX <- mlDataX[index,]
    trainY <- mlDataY[index,]$disease_type
    testX <- mlDataX[-index,]
    testY <- mlDataY[-index,]$disease_type
    
    refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
    refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
    
    set.seed(seedVal)
    ctrl <- trainControl(method = "repeatedcv",
                         number = numKFold,
                         repeats = numResampleIter,
                         sampling = "up",
                         summaryFunction = twoClassSummary,
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
                     tuneGrid = defaultGBMGrid)
    
    positiveClass <- "Cancer"
    negativeClass <- "Control"
    predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
    fg[[ii]] <- predProbs[refactoredTestY == positiveClass]
    bg[[ii]] <- predProbs[refactoredTestY == negativeClass]
    
    prroc_roc[[ii]] <- roc.curve(scores.class0 = fg[[ii]], scores.class1 = bg[[ii]], curve = T)
    plot(prroc_roc[[ii]])
    prroc_pr[[ii]] <- pr.curve(scores.class0 = fg[[ii]], scores.class1 = bg[[ii]], curve = T, rand.compute=T)
    plot(prroc_pr[[ii]])
    
    perf[[ii]] <- data.frame(aucroc = prroc_roc[[ii]]$auc,
                             aupr = prroc_pr[[ii]]$auc.integral)
    
    rm(mlModel)
    
  }
  
  perfAll <- do.call(rbind, perf)
  print(perfAll)
  
  require(gmodels)
  print(ci(perfAll$aucroc))
  print(ci(perfAll$aupr))
  
  return(perfAll)
  
}

#------------------------------------------------


## Main function ##
permPerfHvsC_Sex <- function(iterPerm = 10,  
                             metaData,
                             readData,
                             permuteFlag = TRUE){
  
  predProbs <- list()
  fg <- list()
  bg <- list()
  perf <- list()
  prroc_roc <- list()
  prroc_pr <- list()
  
  sexPermList <- t(lapply(1:iterPerm, 
                            function(x) sample(metaData$sex, 
                                               length(metaData$sex), 
                                               replace = FALSE)))
  
  # if(permuteFlag == FALSE){
  #   iterPerm <- 1
  # }
  
  qcMetadata <- metaData # metadataPSPLMatchedQC
  
  # set.seed(floor(runif(1)*500))
  
  # if(tolower(permuteRankOrNominalAge) == "rank"){
  #   agePrePermute <- match(metadataPSMatchedQC$host_age, sort(unique(metadataPSMatchedQC$host_age)))
  #   qcMetadata$host_age_perm <-  ordered(sample(agePrePermute, length(agePrePermute)))
  # } else if (tolower(permuteRankOrNominalAge) == "nominal"){
  #   agePrePermute <- metadataPSMatchedQC$host_age
  #   qcMetadata$host_age_perm <- sample(agePrePermute, length(agePrePermute))
  # }
  
  for(ii in 1:iterPerm){
    
    print(sprintf("Iteration: %d/%d", ii, iterPerm))
    
    # qcMetadata <- metadataPSMatchedQC # metadataPSPLMatchedQC
    
    # set.seed(floor(runif(1)*500))
    
    qcMetadata$sex_perm <- factor(sexPermList[[ii]])
    
    if(permuteFlag == FALSE){
      qcMetadata$sex_perm <- factor(metaData$sex)
    }
    
    print(qcMetadata$sex_perm)
    
    qcData <- readData # dataPSPLUniqueQC
    
    # Set up design matrix
    covDesignNorm <- model.matrix(~0 + disease_type +
                                    sex_perm,
                                  data = qcMetadata)
    
    # Check row dimensions
    dim(covDesignNorm)[1] == dim(qcData)[1]
    
    # print(colnames(covDesignNorm))
    colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
    # print(colnames(covDesignNorm))
    
    # Set up counts matrix
    counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
    
    # Normalize using edgeR and then plug into voom
    dge <- DGEList(counts = counts)
    vdge__perm <- voom(dge, design = covDesignNorm, 
                       plot = TRUE, save.plot = TRUE, 
                       normalize.method="quantile")
    
    # Apply
    bio.var <- model.matrix(~disease_type,
                            data=qcMetadata)
    
    adj.var <- model.matrix(~sex_perm,
                            data=qcMetadata)
    
    colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
    colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
    
    snmDataObjOnly__perm <- snm(raw.dat = vdge__perm$E, 
                                 bio.var = bio.var, 
                                 adj.var = adj.var, 
                                 rm.adj=TRUE,
                                 verbose = TRUE,
                                 diagnose = TRUE)
    snmData__perm <- t(snmDataObjOnly__perm$norm.dat)
    
    metaTmp1 <- metaData
    metaTmp1$disease_type <- metaData$HvsC
    
    mlDataY <- metaTmp1
    mlDataX <- snmData__perm[rownames(mlDataY),]
    
    if(permuteFlag){
      seedVal <- 42
      print(sprintf("Random Number Seed Value: %d",seedVal))
    } else{
      seedVal <- sample(1:1e3, 1)
      print(sprintf("Random Number Seed Value: %d",seedVal))
    }
    
    set.seed(seedVal)
    index <- createDataPartition(mlDataY$disease_type, p = 0.7, list = FALSE)
    trainX <- mlDataX[index,]
    trainY <- mlDataY[index,]$disease_type
    testX <- mlDataX[-index,]
    testY <- mlDataY[-index,]$disease_type
    
    refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
    refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
    
    set.seed(seedVal)
    ctrl <- trainControl(method = "repeatedcv",
                         number = numKFold,
                         repeats = numResampleIter,
                         sampling = "up",
                         summaryFunction = twoClassSummary,
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
                     tuneGrid = defaultGBMGrid)
    
    positiveClass <- "Cancer"
    negativeClass <- "Control"
    predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
    fg[[ii]] <- predProbs[refactoredTestY == positiveClass]
    bg[[ii]] <- predProbs[refactoredTestY == negativeClass]
    
    prroc_roc[[ii]] <- roc.curve(scores.class0 = fg[[ii]], scores.class1 = bg[[ii]], curve = T)
    plot(prroc_roc[[ii]])
    prroc_pr[[ii]] <- pr.curve(scores.class0 = fg[[ii]], scores.class1 = bg[[ii]], curve = T, rand.compute=T)
    plot(prroc_pr[[ii]])
    
    perf[[ii]] <- data.frame(aucroc = prroc_roc[[ii]]$auc,
                             aupr = prroc_pr[[ii]]$auc.integral)
    
    rm(mlModel)
    
  }
  
  perfAll <- do.call(rbind, perf)
  print(perfAll)
  
  require(gmodels)
  print(ci(perfAll$aucroc))
  print(ci(perfAll$aupr))
  
  return(perfAll)
  
}




vsnm <- function(orderedBoolean){
  ## Load packages ##
  require(limma)
  require(edgeR)
  require(dplyr)
  require(snm)
  require(doMC)
  require(tibble)
  require(gbm)
  
  numCores <- detectCores()
  registerDoMC(cores=numCores)
  
  qcMetadata <- metadataPSMatchedQC # metadataPSPLMatchedQC
  qcData <- vbMergedDataKrakenCFDecontamQC # dataPSPLUniqueQC
  
  # Set up design matrix
  covDesignNorm <- model.matrix(~0 + disease_type +
                                  host_age_perm,
                                  data = qcMetadata)
  
  # Check row dimensions
  dim(covDesignNorm)[1] == dim(qcData)[1]
  
  print(colnames(covDesignNorm))
  colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
  print(colnames(covDesignNorm))
  
  # Set up counts matrix
  counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
  
  # Normalize using edgeR and then plug into voom
  dge <- DGEList(counts = counts)
  vdge_vbMergedDataKrakenCFDecontamQC <<- voom(dge, design = covDesignNorm, 
                                               plot = TRUE, save.plot = TRUE, 
                                               normalize.method="quantile")
  
  # Apply
  bio.var <- model.matrix(~disease_type,
                          data=qcMetadata)
  
  adj.var <- model.matrix(~host_age_perm,
                            data=qcMetadata)
  
  colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
  colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
  print(dim(adj.var))
  print(dim(bio.var))
  print(dim(t(vdge_vbMergedDataKrakenCFDecontamQC$E)))
  print(dim(covDesignNorm))
  
  snmDataObjOnlyVBMergedDataKrakenCFDecontamQC <- snm(raw.dat = vdge_vbMergedDataKrakenCFDecontamQC$E, 
                                                      bio.var = bio.var, 
                                                      adj.var = adj.var, 
                                                      rm.adj=TRUE,
                                                      verbose = TRUE,
                                                      diagnose = TRUE)
  snmDataKrakenCFDecontamQC <<- t(snmDataObjOnlyVBMergedDataKrakenCFDecontamQC$norm.dat)
  
}

mlHvsC <- function(snmData){
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
  
  table(metadataPSMatchedQC$disease_type)
  
  caretTuneGrid <- defaultGBMGrid
  numKFold <- 4
  numResampleIter <- 1
  metaTmp1 <- droplevels(metadataPSMatchedQC[(metadataPSMatchedQC$disease_type %in% c("Prostate Cancer",
                                                                                      "Skin Cutaneous Melanoma",
                                                                                      "Lung Adenocarcinoma")) |
                                               (metadataPSMatchedQC$disease_type %in% c("healthy control") & 
                                                  metadataPSMatchedQC$hiv_status_clean == "HIV-"),])
  
  tmp <- metaTmp1
  tmp$disease_type <- factor(ifelse(tmp$disease_type == "healthy control", yes = "healthy", no = "cancer"))
  
  mlDataY <- tmp
  mlDataX <- snmData[rownames(mlDataY),]
  
  set.seed(42)
  index <- createDataPartition(mlDataY$disease_type, p = 0.7, list = FALSE)
  trainX <- mlDataX[index,]
  trainY <- mlDataY[index,]$disease_type
  testX <- mlDataX[-index,]
  testY <- mlDataY[-index,]$disease_type
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
                   tuneGrid = defaultGBMGrid)
  
  positiveClass <- "cancer"
  negativeClass <- "healthy"
  predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
  fg <- predProbs[refactoredTestY == positiveClass]
  bg <- predProbs[refactoredTestY == negativeClass]
  
  prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  
  plot(prroc_roc)
  plot(prroc_pr)
  
  predClass <- predict(mlModel, newdata = testX)
  print(confusionMatrix(data = predClass, reference = refactoredTestY, positive = positiveClass))
}

loocvDTs <- function(snmData, samplingSize = 15, DTs, caretTuneGrid = defaultGBMGrid){
  
  defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                                 n.trees = floor((1:3) * 50),
                                 shrinkage = 0.1,
                                 n.minobsinnode = 5)
  customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                                n.trees = floor((1:3) * 50),
                                shrinkage = 0.1,
                                n.minobsinnode = 1)
  
  metaTmpX <- droplevels(metadataPSMatchedQC[(metadataPSMatchedQC$disease_type %in% DTs),])
  
  # Do LOOCV model building and testing
  classes <- gsub(" ","",DTs)
  multiClassSummaryStats <- list()
  multiClassSummaryStatsDist <- list()
  numKFold <- 4
  numResampleIter <- 1
  # caretTuneGrid <- defaultGBMGrid #defaultGBMGrid # customGBMGrid
  metaData <- metaTmpX
  snmData <- snmData # dataPSUniqueDecontamQC # 
  iterSize <- 1
  for(jj in 1:iterSize){
    metadataSimSampled <- as.data.frame(stratified(metaData,
                                                   group = "disease_type",
                                                   size = samplingSize,
                                                   keep.rownames = TRUE,
                                                   replace = FALSE,
                                                   bothSets = FALSE))
    rownames(metadataSimSampled) <- metadataSimSampled$rn
    mlDataY <- metadataSimSampled
    mlDataX <- snmData[rownames(mlDataY),]
    dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
    
    # Create data partitions
    # set.seed(42)
    indexSuper <- 1:dim(mlDataY)[1]
    predProbs <- list()
    obsClass <- vector()
    predClass <- vector()
    varImpBestModelDF2OrderedNonzeroList <- list()
    
    for(ii in 1:length(indexSuper)){
      index <- indexSuper[ii]
      print(index)
      trainX <- mlDataX[-index,]
      trainY <- mlDataY[-index,]$disease_type
      testX <- mlDataX[index,,drop=FALSE]
      testY <- mlDataY[index,,drop=FALSE]$disease_type
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
    
    loocvPreds <- cbind(obs = factor(obsClass,
                                     levels = classes),
                        pred = factor(predClass,
                                      levels = classes),
                        do.call(rbind,predProbs))
    # multiClassSummaryStats <- multiClassSummary(loocvPreds, lev = classes)
    # print(multiClassSummaryStats)
    
    multiClassSummaryStats[[jj]] <- multiClassSummary(loocvPreds, lev = classes)
    print(multiClassSummaryStats[[jj]])
  }
  
  print(confusionMatrix(loocvPreds$obs, loocvPreds$pred))
  multiClassSummaryStatsDist <- data.frame(do.call(rbind, multiClassSummaryStats))
}

#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------

permPerfLOOCV <- function(iterPerm = 10, permuteRankOrNominalAge = "rank", permuteFlag = TRUE,
                          samplingSize = 15, DTs, caretTuneGrid = defaultGBMGrid){
  
  predProbs <- list()
  fg <- list()
  bg <- list()
  perf <- list()
  prroc_roc <- list()
  prroc_pr <- list()
  multiClassSummaryStats <- list()
  multiClassSummaryStatsDist <- list()
  
  if(tolower(permuteRankOrNominalAge) == "rank"){
    agePrePermuteRank <- match(metadataPSMatchedQC$host_age, sort(unique(metadataPSMatchedQC$host_age)))
    agePermList <- t(lapply(1:iterPerm, 
                            function(x) sample(agePrePermuteRank, 
                                               length(agePrePermuteRank), 
                                               replace = FALSE)))
  } else if (tolower(permuteRankOrNominalAge) == "nominal"){
    agePermList <- t(lapply(1:iterPerm, 
                            function(x) sample(metadataPSMatchedQC$host_age, 
                                               length(metadataPSMatchedQC$host_age), 
                                               replace = FALSE)))
  } 
  
  if(permuteFlag == FALSE){
    iterPerm <- 1
  }
  
  qcMetadata <- metadataPSMatchedQC # metadataPSPLMatchedQC
  
  for(kk in 1:iterPerm){
    
    print(sprintf("Iteration: %d/%d", kk, iterPerm))
    
    if(tolower(permuteRankOrNominalAge) == "rank"){
      qcMetadata$host_age_perm <-  ordered(agePermList[[kk]])
    } else if (tolower(permuteRankOrNominalAge) == "nominal"){
      qcMetadata$host_age_perm <- agePermList[[kk]]
    }
    
    if(permuteFlag == FALSE & tolower(permuteRankOrNominalAge) == "nominal"){
      qcMetadata$host_age_perm <- metadataPSMatchedQC$host_age
    } else if (permuteFlag == FALSE & tolower(permuteRankOrNominalAge) == "rank"){
      qcMetadata$host_age_perm <- ordered(metadataPSMatchedQC$host_age)
    }
    
    print(qcMetadata$host_age_perm)
    
    qcData <- vbMergedDataKrakenCFDecontamQC # dataPSPLUniqueQC
    
    # Set up design matrix
    covDesignNorm <- model.matrix(~0 + disease_type +
                                    host_age_perm,
                                  data = qcMetadata)
    
    colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
    
    # Set up counts matrix
    counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
    
    # Normalize using edgeR and then plug into voom
    dge <- DGEList(counts = counts)
    vdge_vbMergedDataKrakenCFDecontamQC_perm <- voom(dge, design = covDesignNorm, 
                                                     plot = TRUE, save.plot = TRUE, 
                                                     normalize.method="quantile")
    
    # Apply
    bio.var <- model.matrix(~disease_type,
                            data=qcMetadata)
    
    adj.var <- model.matrix(~host_age_perm,
                            data=qcMetadata)
    
    colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
    colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
    
    snmDataObjOnlyVBMergedDataKrakenCFDecontamQC_perm <- snm(raw.dat = vdge_vbMergedDataKrakenCFDecontamQC_perm$E, 
                                                             bio.var = bio.var, 
                                                             adj.var = adj.var, 
                                                             rm.adj=TRUE,
                                                             verbose = TRUE,
                                                             diagnose = TRUE)
    snmDataKrakenCFDecontamQC_perm <- t(snmDataObjOnlyVBMergedDataKrakenCFDecontamQC_perm$norm.dat)
    
    metaTmpX <- droplevels(metadataPSMatchedQC[(metadataPSMatchedQC$disease_type %in% DTs),])
    
    # Do LOOCV model building and testing
    classes <- gsub(" ","",DTs)
    numKFold <- 4
    numResampleIter <- 1
    # caretTuneGrid <- defaultGBMGrid #defaultGBMGrid # customGBMGrid
    metaData <- metaTmpX
    snmData <- snmDataKrakenCFDecontamQC_perm # dataPSUniqueDecontamQC # 
    iterSize <- 1
    for(jj in 1:iterSize){
      metadataSimSampled <- as.data.frame(stratified(metaData,
                                                     group = "disease_type",
                                                     size = samplingSize,
                                                     keep.rownames = TRUE,
                                                     replace = FALSE,
                                                     bothSets = FALSE))
      rownames(metadataSimSampled) <- metadataSimSampled$rn
      mlDataY <- metadataSimSampled
      mlDataX <- snmData[rownames(mlDataY),]
      dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
      
      # Create data partitions
      # set.seed(42)
      indexSuper <- 1:dim(mlDataY)[1]
      predProbs <- list()
      obsClass <- vector()
      predClass <- vector()
      varImpBestModelDF2OrderedNonzeroList <- list()
      
      for(ii in 1:length(indexSuper)){
        print(sprintf("Sub-iteration: %d/%d",ii,length(indexSuper)))
        index <- indexSuper[ii]
        # print(index)
        trainX <- mlDataX[-index,]
        trainY <- mlDataY[-index,]$disease_type
        testX <- mlDataX[index,,drop=FALSE]
        testY <- mlDataY[index,,drop=FALSE]$disease_type
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
      
      loocvPreds <- cbind(obs = factor(obsClass,
                                       levels = classes),
                          pred = factor(predClass,
                                        levels = classes),
                          do.call(rbind,predProbs))
      # multiClassSummaryStats <- multiClassSummary(loocvPreds, lev = classes)
      # print(multiClassSummaryStats)
      
      multiClassSummaryStats[[jj]] <- multiClassSummary(loocvPreds, lev = classes)
      print(multiClassSummaryStats[[jj]])
    }
    
    print(confusionMatrix(loocvPreds$obs, loocvPreds$pred))
    multiClassSummaryStatsDist[[kk]] <- data.frame(do.call(rbind, multiClassSummaryStats))
    
  }
  
  perfAll <- do.call(rbind, multiClassSummaryStatsDist)
  print(perfAll)
  
  # require(gmodels)
  # print(ci(perfAll$aucroc))
  # print(ci(perfAll$aupr))
  
  return(perfAll)
  
}

