# Shogun_Bootstrap500_cfDNA_Nov21.R
# Author: Greg Poore
# Date: Dec 2, 2019
# Purpose: Bootstrap of healthy vs cancer

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
# YOU SHOULD IMPORT A METADATA DATA FRAME AND A VOOM-SNM NORMALIZED COUNT MATRIX (OR DATA FRAME)
# YOU COULD ALSO LOAD A .RDATA FILE HERE, AS SHOWN BELOW:
# load("snmKrakenAndMetadataFiltered_Dec2_Final.RData") 

# A DEFAULT GBM GRIDSEARCH IS GIVEN BELOW (DEFAULT PER THE CARET PACKAGE). IT CAN BE CUSTOMIZED AS DESIRED.
defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                               n.trees = floor((1:3) * 50),
                               shrinkage = 0.1,
                               n.minobsinnode = 5)
# MAIN FUNCTION
bootstrapPerf <- function(metaData, snmData, iterSizeBootstrap=500, fileName){
	tmp <- metaData
	tmp$disease_type <- metaData$HvsC
	numKFold <- 4
	numResampleIter <- 1
	caretTuneGrid <- defaultGBMGrid

	mlDataY <- tmp
	mlDataX <- snmData[rownames(mlDataY),]
	dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check

	predProbs <- list()
	obsClass <- list()
	predClass <- list()
	varImpBestModelDF2OrderedNonzeroList <- list()

	fg <- list()
	bg <- list()
	perf <- list()
	prroc_roc <- list()
	prroc_pr <- list()

	for(ii in 1:iterSizeBootstrap){
	  print(sprintf("Iteration: %d", ii))
	  index <- createDataPartition(mlDataY$disease_type, p = 0.7, list = FALSE)
	  trainX <- mlDataX[index,]
	  trainY <- mlDataY[index,]$disease_type
	  testX <- mlDataX[-index,]
	  testY <- mlDataY[-index,]$disease_type
	  
	  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
	  refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
	  
	  obsClass[[ii]] <- as.character(refactoredTestY)
	  
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

	require(gmodels)
	require(ggplot2)
	require(gridExtra)

	perfAll <- do.call(rbind, perf)

	write.csv(perfAll, file = paste0("perfAll__iter",iterSizeBootstrap,"__",fileName,".csv"))
	print(perfAll)
	
	aucrocCI <- ci(perfAll$aucroc)
	auprCI <- ci(perfAll$aupr)

	tmpROC <- data.frame(do.call(rbind,lapply(prroc_roc, '[[', 3))[,-3])
	colnames(tmpROC) <- c("FPR","Sensitivity")
	ggplot(tmpROC,aes(x=FPR,y=Sensitivity)) + coord_equal() + 
	  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
	  scale_fill_distiller(palette=4, direction=-1) +
	  scale_x_continuous(expand = c(0, 0), limits = c(0,1)) +
	  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
	  theme(legend.position='none') +
	  annotate(geom="text", 
	label = paste0("Mean AUCROC = ",round(aucrocCI[1],3),"\n95% CI: [",
	               round(aucrocCI[2],3),",",round(aucrocCI[3],3),"]"), 
	           color = "white", x = 0.6, y = 0.3, size = 6) +
	  ggtitle(paste0("ROC Density Plot (",iterSizeBootstrap," iterations)")) -> ggROC

	tmpPR <- data.frame(do.call(rbind,lapply(prroc_pr, '[[', 4))[,-3])
	colnames(tmpPR) <- c("Recall","Precision")
	ggplot(tmpPR,aes(x=Recall,y=Precision)) + coord_equal() + 
	  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
	  scale_fill_distiller(palette=4, direction=-1) +
	  scale_x_continuous(expand = c(0, 0), limits = c(0,1)) +
	  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
	  theme(legend.position='none') +
	  annotate(geom="text", 
	           label = paste0("Mean AUPR = ",round(auprCI[1],3),"\n95% CI: [",
	                          round(auprCI[2],3),",",round(auprCI[3],3),"]"), 
	           color = "white", x = 0.4, y = 0.3, size = 6) +
	  ggtitle(paste0("PR Density Plot (",iterSizeBootstrap," iterations)")) -> ggPR

	# Save density plot data
	write.csv(tmpROC, file = paste0("rocDensityData__iter",iterSizeBootstrap,"__",fileName,".csv"))
	write.csv(tmpPR, file = paste0("prDensityData__iter",iterSizeBootstrap,"__",fileName,".csv"))

	ggROCPR <- ggarrange(ggROC, ggPR, labels = c("a","b"), 
	          ncol = 2, nrow = 1)

	ggROCPRannotated <- annotate_figure(ggROCPR,
	                fig.lab.face = "bold",
	                top = paste0("Discriminating healthy vs cancer samples using decontaminated\nplasma cell-free microbial DNA (",
	                             table(tmp$disease_type)[1]," cancer samples | ",table(tmp$disease_type)[2]," healthy samples)"))

	save(ggROCPRannotated, file = paste0("ggROCPRannotated_iter",iterSizeBootstrap,fileName,".RData"))

	ggsave(plot = ggROCPRannotated, 
		filename = paste0("ggROCPRannotated_iter",iterSizeBootstrap,fileName,".svg"), 
		width = 14, units = "in", dpi = "retina")


}

# THE METADATA DATA FRAME LOADED (ABOVE) WAS CALLED: "metadataPSMatchedDPQCFiltered"
# IT HAD A COLUMN CALLED "disease_type_consol" THAT CONTAINED SAMPLE TYPE NAMES FOR EACH OF THE SAMPLES
# KEY: PRAD = PROSTATE CANCER | NSCLC = LUNG CANCER (BROADLY) | Control = Control
metaTmpPRADvsHealthy <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in%
 																							c("PRAD", "Control")),])

metaTmpNSCLCvsHealthy <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in%
 																							c("NSCLC", "Control")),])

# CALL FUNCTION AND RUN
bootstrapPerf(metaData = metadataPSMatchedDPQCFiltered, snmData = snmDataKrakenCFDecontamDPQC, fileName = "AllDTs_Dec2")

bootstrapPerf(metaData = metaTmpPRADvsHealthy, snmData = snmDataKrakenCFDecontamDPQC, fileName = "PRADvsHealthy_Dec2")

bootstrapPerf(metaData = metaTmpNSCLCvsHealthy, snmData = snmDataKrakenCFDecontamDPQC, fileName = "NSCLCvsHealthy_Dec2")

