# All_Tumor_sequence_center_effects.R
# Author: Greg Poore
# Date: Aug 21, 2018
# Purpose: To explore the contribution of sequence center to cancer microbiome profiles

#-------------------------------#
# Load dependencies
require(devtools)
require(ggbiplot)
require(doMC)
require(plotly)

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("tcgaVbDataAndMetadataAndSNM.RData")
load("tcgaPathBinnedData.RData") # Loads snmDataSampleTypeWithExpStrategyPath, metadataSamplesAllQCPath
load("snmDataSampleTypeWithExpStrategyFINAL.RData")
snmDataSampleTypeWithExpStrategyPath <- snmDataSampleTypeWithExpStrategy[rownames(metadataSamplesAllQCPath),]

metadataSamplesAllQCPath %>%
  filter(sample_type == "Primary Tumor") %>%
  # filter(pathologic_stage_label_binned == "Stage1" | pathologic_stage_label_binned == "Stage4") %>%
  group_by(disease_type,pathologic_stage_label_binned) %>%
  summarise(n=n()) %>%
  data.frame() -> tmp

#---------------------------------------------------------------------#
defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                               n.trees = floor((1:3) * 50),
                               shrinkage = 0.1,
                               n.minobsinnode = 5)
customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                              n.trees = floor((1:3) * 50),
                              shrinkage = 0.1,
                              n.minobsinnode = 1)

source("predcaretLvsHPathFA.R")
caretLvsHPathFABIC <- caretLvsHPathFANoAge(cancerTypeString = "Breast Invasive Carcinoma", 
                                      allMetadataPath = metadataSamplesAllQCPath, 
                                     snmDataML = snmDataSampleTypeWithExpStrategyPath,
                                     caretModel = "gbm",
                                     lowGradePath = "Stage1",
                                     highGradePath = "Stage4",
                                     caretTuneGrid = customGBMGrid,
                                     numResampleIter = 1,
                                     numKFold = 4,
                                     testSetProp = 0.7)

caretLvsHPathFACOAD <- caretLvsHPathFANoAge(cancerTypeString = "Colon Adenocarcinoma", 
                                 allMetadataPath = metadataSamplesAllQCPath, 
                                 snmDataML = snmDataSampleTypeWithExpStrategyPath,
                                 caretModel = "gbm",
                                 lowGradePath = "Stage1",
                                 highGradePath = "Stage4",
                                 caretTuneGrid = defaultGBMGrid,
                                 numResampleIter = 1,
                                 numKFold = 4,
                                 testSetProp = 0.7)

caretLvsHPathFASTAD <- caretLvsHPathFANoAge(cancerTypeString = "Stomach Adenocarcinoma", 
                                   allMetadataPath = metadataSamplesAllQCPath, 
                                   snmDataML = snmDataSampleTypeWithExpStrategyPath,
                                   caretModel = "gbm",
                                   lowGradePath = "Stage1",
                                   highGradePath = "Stage4",
                                   caretTuneGrid = customGBMGrid,
                                   numResampleIter = 1,
                                   numKFold = 4,
                                   testSetProp = 0.7)

source("predcaretLvsHPathFA.R")
source("predcaretLvsHPathFANoAge.R")
caretLvsHPathFATHCA <- caretLvsHPathFANoAge(cancerTypeString = "Thyroid Carcinoma", 
                                   allMetadataPath = metadataSamplesAllQCPath, 
                                   snmDataML = snmDataSampleTypeWithExpStrategyPath,
                                   caretModel = "gbm",
                                   lowGradePath = "Stage1",
                                   highGradePath = "Stage4",
                                   samplingStrategy = "up",
                                   caretTuneGrid = customGBMGrid,
                                   numResampleIter = 1,
                                   numKFold = 4,
                                   testSetProp = 0.7)

customGBMGridPancancer <-  expand.grid(interaction.depth = seq(1,10),
                                  n.trees = floor((1:3) * 50),
                                  shrinkage = 0.05,
                                  n.minobsinnode = 1)
caretLvsHPathFAX <- caretLvsHPathFA(cancerTypeString = "Pan Cancer", 
                                   allMetadataPath = metadataSamplesAllQCPath, 
                                   snmDataML = snmDataSampleTypeWithExpStrategyPath,
                                   tissueType = "Primary Tumor",
                                   caretModel = "gbm",
                                   lowGradePath = "Stage1",
                                   highGradePath = "Stage4",
                                   caretTuneGrid = customGBMGridPancancer,
                                   numResampleIter = 1,
                                   numKFold = 1,
                                   testSetProp = 0.7)



