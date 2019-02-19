# All_Tumor_sequence_center_effects.R
# Author: Greg Poore
# Date: Aug 21, 2018
# Purpose: To explore the contribution of sequence center to cancer microbiome profiles

#-------------------------------#
# Load dependencies
require(devtools)
require(ggbiplot)
require(doMC)
# require(caretEnsemble)
# require(plotly)

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("tcgaVbDataAndMetadataAndSNM.RData")
load("snmDataSampleTypeWithExpStrategyFINAL.RData")

#----------------------------------------------------

## To avoid the error below when training the GBM learner, a custom GBM grid can be entered
# The dataset size is too small or subsampling rate is too large: nTrain*bag.fraction <= n.minobsinnode
defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                        n.trees = floor((1:3) * 50),
                        shrinkage = 0.1,
                        n.minobsinnode = 5)

customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                        n.trees = floor((1:3) * 50),
                        shrinkage = 0.1,
                        n.minobsinnode = 1)

glmnetGrid <- expand.grid(.alpha = seq(.05, 1, length = 15),
                          .lambda = c((1:5)/10))

### Testing ###
# 1252 T / 129 N

# source("predTumorVsNormalEnsembleFA.R")

source("predTumorVsNormalFA_Patent.R")
bicPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                  qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                  cancerTypeString = "Breast Invasive Carcinoma",
                                                  nlTissueComparison = "Solid Tissue Normal",
                                                  tumorTissueComparison = "Primary Tumor", 
                                                  flagForAlternativeNormal = FALSE,
                                                  cancerTypeStringAlternativeNormal = NULL,
                                                  caretModel = "gbm",
                                                  samplingStrategy = "up",
                                                  numResampleIter = 1,
                                                  numKFold = 2,
                                                  trainSetProp = 0.7,
                                                  plotPath = "./roc-ggplots-predicting-tumor-vs-normal",
                                                  caretTuneGrid = customGBMGrid)

# tmpPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
#                                                   qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
#                                                   cancerTypeString = "Uterine Corpus Endometrial Carcinoma",
#                                                   nlTissueComparison = "Solid Tissue Normal",
#                                                   tumorTissueComparison = "Primary Tumor", 
#                                                   flagForAlternativeNormal = FALSE,
#                                                   cancerTypeStringAlternativeNormal = NULL,
#                                                   caretModel = "gbm",
#                                                   samplingStrategy = "up",
#                                                   numResampleIter = 1,
#                                                   numKFold = 2,
#                                                   trainSetProp = 0.7,
#                                                   plotPath = "./roc-ggplots-predicting-tumor-vs-normal",
#                                                   caretTuneGrid = customGBMGrid)


#---------- Predict tumor vs normal using cancer microbiome in Acute Myeloid Leukemia ----------#

## NB: AML is not in the QC metadata

#---------- Predict tumor vs normal using cancer microbiome in Adrenocortical Carcinoma ----------#

# 79 T / 110 N (RCC)
adrenocorticalcarcinomaPredTumorVsRCCNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                                        qcMLDataSNM = snmDataSampleType,
                                                                        cancerTypeString = "Adrenocortical Carcinoma",
                                                                        nlTissueComparison = "Solid Tissue Normal",
                                                                        tumorTissueComparison = "Primary Tumor", 
                                                                        flagForAlternativeNormal = TRUE,
                                                                        cancerTypeStringAlternativeNormal = "Kidney Renal Clear Cell Carcinoma",
                                                                        caretModel = "gbm",
                                                                        numResampleIter = 1,
                                                                        numKFold = 2,
                                                                        testSetProp = 0.7,
                                                                        caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Bladder Urothelial Carcinoma ----------#

# 557 T / 48 N
bladderPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                                    qcMLDataSNM = snmDataSampleType,
                                                                    cancerTypeString = "Bladder Urothelial Carcinoma",
                                                                    nlTissueComparison = "Solid Tissue Normal",
                                                                    tumorTissueComparison = "Primary Tumor", 
                                                                    flagForAlternativeNormal = FALSE,
                                                                    cancerTypeStringAlternativeNormal = NULL,
                                                                    caretModel = "gbm",
                                                                    numResampleIter = 1,
                                                                    numKFold = 2,
                                                                    testSetProp = 0.7,
                                                                    caretTuneGrid = defaultGBMGrid)

# 557 T / 124 NB
bladderPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                    qcMLDataSNM = snmDataSampleType,
                                                    cancerTypeString = "Bladder Urothelial Carcinoma",
                                                    nlTissueComparison = "Blood Derived Normal",
                                                    tumorTissueComparison = "Primary Tumor", 
                                                    flagForAlternativeNormal = FALSE,
                                                    cancerTypeStringAlternativeNormal = NULL,
                                                    caretModel = "gbm",
                                                    numResampleIter = 1,
                                                    numKFold = 2,
                                                    testSetProp = 0.7,
                                                    caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Brain Lower Grade Glioma ----------#

# 607 T / 92 NB
blggPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                    qcMLDataSNM = snmDataSampleType,
                                                    cancerTypeString = "Brain Lower Grade Glioma",
                                                    nlTissueComparison = "Blood Derived Normal",
                                                    tumorTissueComparison = "Primary Tumor", 
                                                    flagForAlternativeNormal = FALSE,
                                                    cancerTypeStringAlternativeNormal = NULL,
                                                    caretModel = "gbm",
                                                    numResampleIter = 1,
                                                    numKFold = 2,
                                                    testSetProp = 0.7,
                                                    caretTuneGrid = defaultGBMGrid)

# 34 RT / 607 T 
blggPredPrimaryTumorVsRecurrentTumor <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                   qcMLDataSNM = snmDataSampleType,
                                                   cancerTypeString = "Brain Lower Grade Glioma",
                                                   nlTissueComparison = "Primary Tumor",
                                                   tumorTissueComparison = "Recurrent Tumor", 
                                                   flagForAlternativeNormal = FALSE,
                                                   cancerTypeStringAlternativeNormal = NULL,
                                                   caretModel = "gbm",
                                                   numResampleIter = 1,
                                                   numKFold = 2,
                                                   testSetProp = 0.7,
                                                   caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Breast Invasive Carcinoma ----------#

# 1252 T / 129 N
bicPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                               qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                               cancerTypeString = "Breast Invasive Carcinoma",
                                               nlTissueComparison = "Solid Tissue Normal",
                                               tumorTissueComparison = "Primary Tumor", 
                                               flagForAlternativeNormal = FALSE,
                                               cancerTypeStringAlternativeNormal = NULL,
                                               caretModel = "gbm",
                                               numResampleIter = 1,
                                               numKFold = 2,
                                               testSetProp = 0.7,
                                               caretTuneGrid = defaultGBMGrid)

# 1252 T / 107 NB
bicPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Breast Invasive Carcinoma",
                                                nlTissueComparison = "Blood Derived Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.7,
                                                caretTuneGrid = defaultGBMGrid)

# 9 M / 1252 T ************
bicGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                            n.trees = floor((1:3) * 50),
                            shrinkage = 0.1,
                            n.minobsinnode = 1)
bicPredPrimaryTumorVsMetastaticTumor <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                             qcMLDataSNM = snmDataSampleType,
                                                             cancerTypeString = "Breast Invasive Carcinoma",
                                                             nlTissueComparison = "Primary Tumor",
                                                             tumorTissueComparison = "Metastatic", 
                                                             flagForAlternativeNormal = FALSE,
                                                             cancerTypeStringAlternativeNormal = NULL,
                                                             caretModel = "gbm",
                                                             numResampleIter = 1,
                                                             numKFold = 2,
                                                             testSetProp = 0.7,
                                                             caretTuneGrid = bicGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma ----------#

# 374 T / 5 N ***************
cervGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                              n.trees = floor((1:3) * 50),
                              shrinkage = 0.1,
                              n.minobsinnode = 1)
cervicalPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",
                                                nlTissueComparison = "Solid Tissue Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.7,
                                                caretTuneGrid = cervGBMGrid)

# 374 T / 70 NB
cervicalPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                               qcMLDataSNM = snmDataSampleType,
                                               cancerTypeString = "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",
                                               nlTissueComparison = "Blood Derived Normal",
                                               tumorTissueComparison = "Primary Tumor", 
                                               flagForAlternativeNormal = FALSE,
                                               cancerTypeStringAlternativeNormal = NULL,
                                               caretModel = "gbm",
                                               numResampleIter = 1,
                                               numKFold = 2,
                                               testSetProp = 0.7,
                                               caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Cholangiocarcinoma ----------#

# 36 T / 9 N *************
cholGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                            n.trees = floor((1:3) * 50),
                            shrinkage = 0.1,
                            n.minobsinnode = 1)
cholPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Cholangiocarcinoma",
                                                nlTissueComparison = "Solid Tissue Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.7,
                                                caretTuneGrid = cholGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Colon Adenocarcinoma ----------#

# 837 T / 70 N
colonPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                 qcMLDataSNM = snmDataSampleType,
                                                 cancerTypeString = "Colon Adenocarcinoma",
                                                 nlTissueComparison = "Solid Tissue Normal",
                                                 tumorTissueComparison = "Primary Tumor", 
                                                 flagForAlternativeNormal = FALSE,
                                                 cancerTypeStringAlternativeNormal = NULL,
                                                 caretModel = "gbm",
                                                 numResampleIter = 1,
                                                 numKFold = 2,
                                                 testSetProp = 0.7,
                                                 caretTuneGrid = defaultGBMGrid)

# 837 T / 109 NB
colonPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                  qcMLDataSNM = snmDataSampleType,
                                                  cancerTypeString = "Colon Adenocarcinoma",
                                                  nlTissueComparison = "Blood Derived Normal",
                                                  tumorTissueComparison = "Primary Tumor", 
                                                  flagForAlternativeNormal = FALSE,
                                                  cancerTypeStringAlternativeNormal = NULL,
                                                  caretModel = "gbm",
                                                  numResampleIter = 1,
                                                  numKFold = 2,
                                                  testSetProp = 0.7,
                                                  caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Esophageal Carcinoma ----------#

# 253 T / 37 N
esophagealcarcinomaPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                  qcMLDataSNM = snmDataSampleType,
                                                  cancerTypeString = "Esophageal Carcinoma",
                                                  nlTissueComparison = "Solid Tissue Normal",
                                                  tumorTissueComparison = "Primary Tumor", 
                                                  flagForAlternativeNormal = FALSE,
                                                  cancerTypeStringAlternativeNormal = NULL,
                                                  caretModel = "gbm",
                                                  numResampleIter = 1,
                                                  numKFold = 2,
                                                  testSetProp = 0.7,
                                                  caretTuneGrid = defaultGBMGrid)

# 253 T / 48 NB
esophagealcarcinomaPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                 qcMLDataSNM = snmDataSampleType,
                                                 cancerTypeString = "Esophageal Carcinoma",
                                                 nlTissueComparison = "Blood Derived Normal",
                                                 tumorTissueComparison = "Primary Tumor", 
                                                 flagForAlternativeNormal = FALSE,
                                                 cancerTypeStringAlternativeNormal = NULL,
                                                 caretModel = "gbm",
                                                 numResampleIter = 1,
                                                 numKFold = 2,
                                                 testSetProp = 0.7,
                                                 caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Glioblastoma Multiforme ----------#

# 405 T / 5 N *********************
# gbmGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
#                             n.trees = floor((1:3) * 50),
#                             shrinkage = 0.1,
#                             n.minobsinnode = 1)
# gbmPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC,
#                                                                 qcMLDataSNM = snmDataSampleType,
#                                                                 cancerTypeString = "Glioblastoma Multiforme",
#                                                                 nlTissueComparison = "Solid Tissue Normal",
#                                                                 tumorTissueComparison = "Primary Tumor",
#                                                                 flagForAlternativeNormal = FALSE,
#                                                                 cancerTypeStringAlternativeNormal = NULL,
#                                                                 caretModel = "gbm",
#                                                                 numResampleIter = 1,
#                                                                 numKFold = 2,
#                                                                 testSetProp = 0.7,
#                                                                 caretTuneGrid = gbmGBMGrid)

# 405 T / 74 NB
gbmPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                               qcMLDataSNM = snmDataSampleType,
                                                               cancerTypeString = "Glioblastoma Multiforme",
                                                               nlTissueComparison = "Blood Derived Normal",
                                                               tumorTissueComparison = "Primary Tumor", 
                                                               flagForAlternativeNormal = FALSE,
                                                               cancerTypeStringAlternativeNormal = NULL,
                                                               caretModel = "gbm",
                                                               numResampleIter = 1,
                                                               numKFold = 2,
                                                               testSetProp = 0.7,
                                                               caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Head and Neck Squamous Cell Carcinoma ----------#

# 692 T / 70 N
hnsccPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Head and Neck Squamous Cell Carcinoma",
                                                nlTissueComparison = "Solid Tissue Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.7,
                                                caretTuneGrid = defaultGBMGrid)

# 692 T / 143 NB
hnsccPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                               qcMLDataSNM = snmDataSampleType,
                                               cancerTypeString = "Head and Neck Squamous Cell Carcinoma",
                                               nlTissueComparison = "Blood Derived Normal",
                                               tumorTissueComparison = "Primary Tumor", 
                                               flagForAlternativeNormal = FALSE,
                                               cancerTypeStringAlternativeNormal = NULL,
                                               caretModel = "gbm",
                                               numResampleIter = 1,
                                               numKFold = 2,
                                               testSetProp = 0.7,
                                               caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Kidney Chromophobe ----------#

# 116 T / 66 N
kidneychromophobePredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Kidney Chromophobe",
                                                nlTissueComparison = "Solid Tissue Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.7,
                                                caretTuneGrid = defaultGBMGrid)

# 116 T / 9 NB ************
# kidneychromophobePredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
#                                                qcMLDataSNM = snmDataSampleType,
#                                                cancerTypeString = "Kidney Chromophobe",
#                                                nlTissueComparison = "Blood Derived Normal",
#                                                tumorTissueComparison = "Primary Tumor", 
#                                                flagForAlternativeNormal = FALSE,
#                                                cancerTypeStringAlternativeNormal = NULL,
#                                                caretModel = "gbm",
#                                                numResampleIter = 1,
#                                                numKFold = 2,
#                                                testSetProp = 0.7,
#                                                caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Kidney Renal Clear Cell Carcinoma ----------#

# 1026 T / 110 N
rccPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                  qcMLDataSNM = snmDataSampleType,
                                                  cancerTypeString = "Kidney Renal Clear Cell Carcinoma",
                                                  nlTissueComparison = "Solid Tissue Normal",
                                                  tumorTissueComparison = "Primary Tumor", 
                                                  flagForAlternativeNormal = FALSE,
                                                  cancerTypeStringAlternativeNormal = NULL,
                                                  caretModel = "gbm",
                                                  numResampleIter = 1,
                                                  numKFold = 2,
                                                  testSetProp = 0.7,
                                                  caretTuneGrid = defaultGBMGrid)

# 1026 T / 5 NB **************
# rccPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
#                                                  qcMLDataSNM = snmDataSampleType,
#                                                  cancerTypeString = "Kidney Renal Clear Cell Carcinoma",
#                                                  nlTissueComparison = "Blood Derived Normal",
#                                                  tumorTissueComparison = "Primary Tumor", 
#                                                  flagForAlternativeNormal = FALSE,
#                                                  cancerTypeStringAlternativeNormal = NULL,
#                                                  caretModel = "gbm",
#                                                  numResampleIter = 1,
#                                                  numKFold = 2,
#                                                  testSetProp = 0.7,
#                                                  caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Kidney Renal Papillary Cell Carcinoma ----------#

# 328 T / 36 N
rpccPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                  qcMLDataSNM = snmDataSampleType,
                                                  cancerTypeString = "Kidney Renal Papillary Cell Carcinoma",
                                                  nlTissueComparison = "Solid Tissue Normal",
                                                  tumorTissueComparison = "Primary Tumor", 
                                                  flagForAlternativeNormal = FALSE,
                                                  cancerTypeStringAlternativeNormal = NULL,
                                                  caretModel = "gbm",
                                                  numResampleIter = 1,
                                                  numKFold = 2,
                                                  testSetProp = 0.7,
                                                  caretTuneGrid = defaultGBMGrid)

# 328 T / 35 NB
rpccPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                 qcMLDataSNM = snmDataSampleType,
                                                 cancerTypeString = "Kidney Renal Papillary Cell Carcinoma",
                                                 nlTissueComparison = "Blood Derived Normal",
                                                 tumorTissueComparison = "Primary Tumor", 
                                                 flagForAlternativeNormal = FALSE,
                                                 cancerTypeStringAlternativeNormal = NULL,
                                                 caretModel = "gbm",
                                                 numResampleIter = 1,
                                                 numKFold = 2,
                                                 testSetProp = 0.7,
                                                 caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Liver Hepatocellular Carcinoma ----------#

# 425 T / 72 N
hccPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                 qcMLDataSNM = snmDataSampleType,
                                                 cancerTypeString = "Liver Hepatocellular Carcinoma",
                                                 nlTissueComparison = "Solid Tissue Normal",
                                                 tumorTissueComparison = "Primary Tumor", 
                                                 flagForAlternativeNormal = FALSE,
                                                 cancerTypeStringAlternativeNormal = NULL,
                                                 caretModel = "gbm",
                                                 numResampleIter = 1,
                                                 numKFold = 2,
                                                 testSetProp = 0.7,
                                                 caretTuneGrid = defaultGBMGrid)

# 425 T / 32 NB
hccPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Liver Hepatocellular Carcinoma",
                                                nlTissueComparison = "Blood Derived Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.7,
                                                caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Lung Adenocarcinoma ----------#

# 711 T / 134 N
ladPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                    qcMLDataSNM = snmDataSampleType,
                                                    cancerTypeString = "Lung Adenocarcinoma",
                                                    nlTissueComparison = "Solid Tissue Normal",
                                                    tumorTissueComparison = "Primary Tumor", 
                                                    flagForAlternativeNormal = FALSE,
                                                    cancerTypeStringAlternativeNormal = NULL,
                                                    caretModel = "gbm",
                                                    numResampleIter = 1,
                                                    numKFold = 2,
                                                    testSetProp = 0.7,
                                                    caretTuneGrid = defaultGBMGrid)

# 711 T / 104 NB
ladPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                    qcMLDataSNM = snmDataSampleType,
                                                    cancerTypeString = "Lung Adenocarcinoma",
                                                    nlTissueComparison = "Blood Derived Normal",
                                                    tumorTissueComparison = "Primary Tumor", 
                                                    flagForAlternativeNormal = FALSE,
                                                    cancerTypeStringAlternativeNormal = NULL,
                                                    caretModel = "gbm",
                                                    numResampleIter = 1,
                                                    numKFold = 2,
                                                    testSetProp = 0.7,
                                                    caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Lung Squamous Cell Carcinoma ----------#

# 552 T / 83 N
lsccPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Lung Squamous Cell Carcinoma",
                                                nlTissueComparison = "Solid Tissue Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.7,
                                                caretTuneGrid = defaultGBMGrid)

# 552 T / 18 NB **************
lsccGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                             n.trees = floor((1:3) * 50),
                             shrinkage = 0.1,
                             n.minobsinnode = 1)
lsccPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                               qcMLDataSNM = snmDataSampleType,
                                               cancerTypeString = "Lung Squamous Cell Carcinoma",
                                               nlTissueComparison = "Blood Derived Normal",
                                               tumorTissueComparison = "Primary Tumor", 
                                               flagForAlternativeNormal = FALSE,
                                               cancerTypeStringAlternativeNormal = NULL,
                                               caretModel = "gbm",
                                               numResampleIter = 1,
                                               numKFold = 2,
                                               testSetProp = 0.7,
                                               caretTuneGrid = lsccGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Lymphoid Neoplasm Diffuse Large B-cell Lymphoma ----------#

# 55 T / 7 NB **************
dlbclGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                          n.trees = floor((1:3) * 50),
                          shrinkage = 0.1,
                          n.minobsinnode = 1)
dlbclPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                 qcMLDataSNM = snmDataSampleType,
                                                 cancerTypeString = "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
                                                 nlTissueComparison = "Blood Derived Normal",
                                                 tumorTissueComparison = "Primary Tumor", 
                                                 flagForAlternativeNormal = FALSE,
                                                 cancerTypeStringAlternativeNormal = NULL,
                                                 caretModel = "gbm",
                                                 numResampleIter = 1,
                                                 numKFold = 2,
                                                 testSetProp = 0.7,
                                                 caretTuneGrid = dlbclGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Mesothelioma ----------#

# 87 T / 134 N (LAD)
mesotheliomaPredTumorVsLADNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                             qcMLDataSNM = snmDataSampleType,
                                                             cancerTypeString = "Mesothelioma",
                                                             nlTissueComparison = "Solid Tissue Normal",
                                                             tumorTissueComparison = "Primary Tumor", 
                                                             flagForAlternativeNormal = TRUE,
                                                             cancerTypeStringAlternativeNormal = "Lung Adenocarcinoma",
                                                             caretModel = "gbm",
                                                             numResampleIter = 1,
                                                             numKFold = 2,
                                                             testSetProp = 0.7,
                                                             caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Ovarian Serous Cystadenocarcinoma ----------#

# 919 T / 16 N ************
ovGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                           n.trees = floor((1:3) * 50),
                           shrinkage = 0.1,
                           n.minobsinnode = 1)
ovarianPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Ovarian Serous Cystadenocarcinoma",
                                                nlTissueComparison = "Solid Tissue Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.7,
                                                caretTuneGrid = ovGBMGrid)

# 919 T / 68 NB
ovarianPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                               qcMLDataSNM = snmDataSampleType,
                                               cancerTypeString = "Ovarian Serous Cystadenocarcinoma",
                                               nlTissueComparison = "Blood Derived Normal",
                                               tumorTissueComparison = "Primary Tumor", 
                                               flagForAlternativeNormal = FALSE,
                                               cancerTypeStringAlternativeNormal = NULL,
                                               caretModel = "gbm",
                                               numResampleIter = 1,
                                               numKFold = 2,
                                               testSetProp = 0.7,
                                               caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Pancreatic Adenocarcinoma ----------#

# 178 T / 4 N ****************
padGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                          n.trees = floor((1:3) * 50),
                          shrinkage = 0.1,
                          n.minobsinnode = 1)
padPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                    qcMLDataSNM = snmDataSampleType,
                                                    cancerTypeString = "Pancreatic Adenocarcinoma",
                                                    nlTissueComparison = "Solid Tissue Normal",
                                                    tumorTissueComparison = "Primary Tumor", 
                                                    flagForAlternativeNormal = FALSE,
                                                    cancerTypeStringAlternativeNormal = NULL,
                                                    caretModel = "gbm",
                                                    numResampleIter = 1,
                                                    numKFold = 2,
                                                    testSetProp = 0.7,
                                                    caretTuneGrid = padGBMGrid,
                                                    bag.fraction = 0.75)

#---------- Predict tumor vs normal using cancer microbiome in Pheochromocytoma and Paraganglioma ----------#

# 179 T / 3 N ****************
pheoGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                           n.trees = floor((1:3) * 50),
                           shrinkage = 0.1,
                           n.minobsinnode = 1)
pheoPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Pheochromocytoma and Paraganglioma",
                                                nlTissueComparison = "Solid Tissue Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.5,
                                                caretTuneGrid = pheoGBMGrid,
                                                bag.fraction = 0.9)

#---------- Predict tumor vs normal using cancer microbiome in Prostate Adenocarcinoma ----------#

# 641 T / 67 N
prostatePredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                 qcMLDataSNM = snmDataSampleType,
                                                 cancerTypeString = "Prostate Adenocarcinoma",
                                                 nlTissueComparison = "Solid Tissue Normal",
                                                 tumorTissueComparison = "Primary Tumor", 
                                                 flagForAlternativeNormal = FALSE,
                                                 cancerTypeStringAlternativeNormal = NULL,
                                                 caretModel = "gbm",
                                                 numResampleIter = 1,
                                                 numKFold = 2,
                                                 testSetProp = 0.7,
                                                 caretTuneGrid = defaultGBMGrid)

# 641 T / 121 NB
prostatePredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Prostate Adenocarcinoma",
                                                nlTissueComparison = "Blood Derived Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.7,
                                                caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Rectum Adenocarcinoma ----------#

# 641 T / 67 N
rectumGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                               n.trees = floor((1:3) * 50),
                               shrinkage = 0.1,
                               n.minobsinnode = 1)
rectumPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                     qcMLDataSNM = snmDataSampleType,
                                                     cancerTypeString = "Rectum Adenocarcinoma",
                                                     nlTissueComparison = "Solid Tissue Normal",
                                                     tumorTissueComparison = "Primary Tumor", 
                                                     flagForAlternativeNormal = FALSE,
                                                     cancerTypeStringAlternativeNormal = NULL,
                                                     caretModel = "gbm",
                                                     numResampleIter = 1,
                                                     numKFold = 2,
                                                     testSetProp = 0.7,
                                                     caretTuneGrid = rectumGBMGrid)

# 641 T / 121 NB
rectumPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                    qcMLDataSNM = snmDataSampleType,
                                                    cancerTypeString = "Rectum Adenocarcinoma",
                                                    nlTissueComparison = "Blood Derived Normal",
                                                    tumorTissueComparison = "Primary Tumor", 
                                                    flagForAlternativeNormal = FALSE,
                                                    cancerTypeStringAlternativeNormal = NULL,
                                                    caretModel = "gbm",
                                                    numResampleIter = 1,
                                                    numKFold = 2,
                                                    testSetProp = 0.7,
                                                    caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Sarcoma ----------#

# 299 T / 7 N *************
sarcomaGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                            n.trees = floor((1:3) * 50),
                            shrinkage = 0.1,
                            n.minobsinnode = 1)
sarcomaPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                   qcMLDataSNM = snmDataSampleType,
                                                   cancerTypeString = "Sarcoma",
                                                   nlTissueComparison = "Solid Tissue Normal",
                                                   tumorTissueComparison = "Primary Tumor", 
                                                   flagForAlternativeNormal = FALSE,
                                                   cancerTypeStringAlternativeNormal = NULL,
                                                   caretModel = "gbm",
                                                   numResampleIter = 1,
                                                   numKFold = 2,
                                                   testSetProp = 0.7,
                                                   caretTuneGrid = sarcomaGBMGrid)

# 299 T / 36 NB
sarcomaPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                  qcMLDataSNM = snmDataSampleType,
                                                  cancerTypeString = "Sarcoma",
                                                  nlTissueComparison = "Blood Derived Normal",
                                                  tumorTissueComparison = "Primary Tumor", 
                                                  flagForAlternativeNormal = FALSE,
                                                  cancerTypeStringAlternativeNormal = NULL,
                                                  caretModel = "gbm",
                                                  numResampleIter = 1,
                                                  numKFold = 2,
                                                  testSetProp = 0.7,
                                                  caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Skin Cutaneous Melanoma ----------#

# 122 T / 158 N
skincutaneousmelanomaPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                    qcMLDataSNM = snmDataSampleType,
                                                    cancerTypeString = "Skin Cutaneous Melanoma",
                                                    nlTissueComparison = "Blood Derived Normal",
                                                    tumorTissueComparison = "Primary Tumor", 
                                                    flagForAlternativeNormal = FALSE,
                                                    cancerTypeStringAlternativeNormal = NULL,
                                                    caretModel = "gbm",
                                                    numResampleIter = 1,
                                                    numKFold = 2,
                                                    testSetProp = 0.7,
                                                    caretTuneGrid = defaultGBMGrid)

# 511 M / 122 T
skincutaneousmelanomaPredMetastaticTumorVsPrimaryTumor <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                                               qcMLDataSNM = snmDataSampleType,
                                                                               cancerTypeString = "Skin Cutaneous Melanoma",
                                                                               nlTissueComparison = "Primary Tumor",
                                                                               tumorTissueComparison = "Metastatic", 
                                                                               flagForAlternativeNormal = FALSE,
                                                                               cancerTypeStringAlternativeNormal = NULL,
                                                                               caretModel = "gbm",
                                                                               numResampleIter = 1,
                                                                               numKFold = 2,
                                                                               testSetProp = 0.7,
                                                                               caretTuneGrid = defaultGBMGrid)


#---------- Predict tumor vs normal using cancer microbiome in Stomach Adenocarcinoma ----------#

# 862 T / 114 N
stomachPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                qcMLDataSNM = snmDataSampleType,
                                                cancerTypeString = "Stomach Adenocarcinoma",
                                                nlTissueComparison = "Solid Tissue Normal",
                                                tumorTissueComparison = "Primary Tumor", 
                                                flagForAlternativeNormal = FALSE,
                                                cancerTypeStringAlternativeNormal = NULL,
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                testSetProp = 0.7,
                                                caretTuneGrid = defaultGBMGrid)

# 862 T / 116 NB
stomachPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                               qcMLDataSNM = snmDataSampleType,
                                               cancerTypeString = "Stomach Adenocarcinoma",
                                               nlTissueComparison = "Blood Derived Normal",
                                               tumorTissueComparison = "Primary Tumor", 
                                               flagForAlternativeNormal = FALSE,
                                               cancerTypeStringAlternativeNormal = NULL,
                                               caretModel = "gbm",
                                               numResampleIter = 1,
                                               numKFold = 2,
                                               testSetProp = 0.7,
                                               caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Testicular Germ Cell Tumors ----------#
          
# 150 T / 121 NB (Prostate Adeno)
testicularPredTumorVsProstateAdenoNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                   qcMLDataSNM = snmDataSampleType,
                                                   cancerTypeString = "Testicular Germ Cell Tumors",
                                                   nlTissueComparison = "Blood Derived Normal",
                                                   tumorTissueComparison = "Primary Tumor", 
                                                   flagForAlternativeNormal = TRUE,
                                                   cancerTypeStringAlternativeNormal = "Prostate Adenocarcinoma",
                                                   caretModel = "gbm",
                                                   numResampleIter = 1,
                                                   numKFold = 2,
                                                   testSetProp = 0.7,
                                                   caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Thymoma ----------#

# 120 T / 79 N (Thyroid Carcinoma)
thymomaPredTumorVsThyroidCarcinomaNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                                   qcMLDataSNM = snmDataSampleType,
                                                                   cancerTypeString = "Thymoma",
                                                                   nlTissueComparison = "Solid Tissue Normal",
                                                                   tumorTissueComparison = "Primary Tumor", 
                                                                   flagForAlternativeNormal = TRUE,
                                                                   cancerTypeStringAlternativeNormal = "Thyroid Carcinoma",
                                                                   caretModel = "gbm",
                                                                   numResampleIter = 1,
                                                                   numKFold = 2,
                                                                   testSetProp = 0.7,
                                                                   caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Thyroid Carcinoma ----------#

# 655 T / 79 N
thyroidcarcinomaPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                    qcMLDataSNM = snmDataSampleType,
                                                    cancerTypeString = "Thyroid Carcinoma",
                                                    nlTissueComparison = "Solid Tissue Normal",
                                                    tumorTissueComparison = "Primary Tumor", 
                                                    flagForAlternativeNormal = FALSE,
                                                    cancerTypeStringAlternativeNormal = NULL,
                                                    caretModel = "gbm",
                                                    numResampleIter = 1,
                                                    numKFold = 2,
                                                    testSetProp = 0.7,
                                                    caretTuneGrid = defaultGBMGrid)

# 655 T / 137 NB
thyroidcarcinomaPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                   qcMLDataSNM = snmDataSampleType,
                                                   cancerTypeString = "Thyroid Carcinoma",
                                                   nlTissueComparison = "Blood Derived Normal",
                                                   tumorTissueComparison = "Primary Tumor", 
                                                   flagForAlternativeNormal = FALSE,
                                                   cancerTypeStringAlternativeNormal = NULL,
                                                   caretModel = "gbm",
                                                   numResampleIter = 1,
                                                   numKFold = 2,
                                                   testSetProp = 0.7,
                                                   caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Uterine Carcinosarcoma ----------#

# 57 T / 45 N (Uterine Corpus Endometrial Carcinoma)
uterinecarcinosarcomaPredTumorVsEndometrialCarcinomaNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                                    qcMLDataSNM = snmDataSampleType,
                                                                    cancerTypeString = "Uterine Carcinosarcoma",
                                                                    nlTissueComparison = "Solid Tissue Normal",
                                                                    tumorTissueComparison = "Primary Tumor", 
                                                                    flagForAlternativeNormal = TRUE,
                                                                    cancerTypeStringAlternativeNormal = "Uterine Corpus Endometrial Carcinoma",
                                                                    caretModel = "gbm",
                                                                    numResampleIter = 1,
                                                                    numKFold = 2,
                                                                    testSetProp = 0.7,
                                                                    caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Uterine Corpus Endometrial Carcinoma ----------#

# 1032 T / 45 N
endometrialcarcinomaPredTumorVsNormalTissue <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                             qcMLDataSNM = snmDataSampleType,
                                                             cancerTypeString = "Uterine Corpus Endometrial Carcinoma",
                                                             nlTissueComparison = "Solid Tissue Normal",
                                                             tumorTissueComparison = "Primary Tumor", 
                                                             flagForAlternativeNormal = FALSE,
                                                             cancerTypeStringAlternativeNormal = NULL,
                                                             caretModel = "gbm",
                                                             numResampleIter = 1,
                                                             numKFold = 2,
                                                             testSetProp = 0.7,
                                                             caretTuneGrid = defaultGBMGrid)

# 1032 T / 161 NB
endometrialcarcinomaPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                            qcMLDataSNM = snmDataSampleType,
                                                            cancerTypeString = "Uterine Corpus Endometrial Carcinoma",
                                                            nlTissueComparison = "Blood Derived Normal",
                                                            tumorTissueComparison = "Primary Tumor", 
                                                            flagForAlternativeNormal = FALSE,
                                                            cancerTypeStringAlternativeNormal = NULL,
                                                            caretModel = "gbm",
                                                            numResampleIter = 1,
                                                            numKFold = 2,
                                                            testSetProp = 0.7,
                                                            caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Uveal Melanoma ----------#

# 131 T / 51 NB
uvealmelanomaPredTumorVsNormalBlood <- predTumorVsNormalFA_Patent(qcMLMetadata = metadataSamplesAllQC, 
                                                                qcMLDataSNM = snmDataSampleType,
                                                                cancerTypeString = "Uveal Melanoma",
                                                                nlTissueComparison = "Blood Derived Normal",
                                                                tumorTissueComparison = "Primary Tumor", 
                                                                flagForAlternativeNormal = FALSE,
                                                                cancerTypeStringAlternativeNormal = NULL,
                                                                caretModel = "gbm",
                                                                numResampleIter = 1,
                                                                numKFold = 2,
                                                                testSetProp = 0.7,
                                                                caretTuneGrid = defaultGBMGrid)

# save.image(file = "mlPredTumorVsNormalAllTumors.RData")

mlModelListTumorVsNormalMicrobes <- list(
  adrenocorticalcarcinomaPredTumorVsRCCNormalTissue = adrenocorticalcarcinomaPredTumorVsRCCNormalTissue,
  bladderPredTumorVsNormalTissue = bladderPredTumorVsNormalTissue,
  bladderPredTumorVsNormalBlood = bladderPredTumorVsNormalBlood,
  blggPredTumorVsNormalBlood = blggPredTumorVsNormalBlood,
  blggPredPrimaryTumorVsRecurrentTumor = blggPredPrimaryTumorVsRecurrentTumor,
  bicPredTumorVsNormalTissue = bicPredTumorVsNormalTissue,
  bicPredTumorVsNormalBlood = bicPredTumorVsNormalBlood,
  bicPredPrimaryTumorVsMetastaticTumor = bicPredPrimaryTumorVsMetastaticTumor,
  cervicalPredTumorVsNormalTissue = cervicalPredTumorVsNormalTissue,
  cervicalPredTumorVsNormalBlood = cervicalPredTumorVsNormalBlood,
  cholPredTumorVsNormalTissue = cholPredTumorVsNormalTissue,
  colonPredTumorVsNormalTissue = colonPredTumorVsNormalTissue,
  colonPredTumorVsNormalBlood = colonPredTumorVsNormalBlood,
  esophagealcarcinomaPredTumorVsNormalTissue = esophagealcarcinomaPredTumorVsNormalTissue,
  esophagealcarcinomaPredTumorVsNormalBlood = esophagealcarcinomaPredTumorVsNormalBlood,
  # gbmPredTumorVsNormalTissue = gbmPredTumorVsNormalTissue, # too few normal tissue samples on which to use for ML
  gbmPredTumorVsNormalBlood = gbmPredTumorVsNormalBlood,
  hnsccPredTumorVsNormalTissue = hnsccPredTumorVsNormalTissue,
  hnsccPredTumorVsNormalBlood = hnsccPredTumorVsNormalBlood,
  kidneychromophobePredTumorVsNormalTissue = kidneychromophobePredTumorVsNormalTissue,
  rccPredTumorVsNormalTissue = rccPredTumorVsNormalTissue,
  rpccPredTumorVsNormalTissue = rpccPredTumorVsNormalTissue,
  rpccPredTumorVsNormalBlood = rpccPredTumorVsNormalBlood,
  hccPredTumorVsNormalTissue = hccPredTumorVsNormalTissue,
  hccPredTumorVsNormalBlood = hccPredTumorVsNormalBlood,
  ladPredTumorVsNormalTissue = ladPredTumorVsNormalTissue,
  ladPredTumorVsNormalBlood = ladPredTumorVsNormalBlood,
  lsccPredTumorVsNormalTissue = lsccPredTumorVsNormalTissue,
  lsccPredTumorVsNormalBlood = lsccPredTumorVsNormalBlood,
  dlbclPredTumorVsNormalBlood = dlbclPredTumorVsNormalBlood,
  mesotheliomaPredTumorVsLADNormalTissue = mesotheliomaPredTumorVsLADNormalTissue,
  ovarianPredTumorVsNormalTissue = ovarianPredTumorVsNormalTissue,
  ovarianPredTumorVsNormalBlood = ovarianPredTumorVsNormalBlood,
  padPredTumorVsNormalTissue = padPredTumorVsNormalTissue,
  pheoPredTumorVsNormalTissue = pheoPredTumorVsNormalTissue,
  prostatePredTumorVsNormalTissue = prostatePredTumorVsNormalTissue,
  prostatePredTumorVsNormalBlood = prostatePredTumorVsNormalBlood,
  rectumPredTumorVsNormalTissue = rectumPredTumorVsNormalTissue,
  rectumPredTumorVsNormalBlood = rectumPredTumorVsNormalBlood,
  sarcomaPredTumorVsNormalTissue = sarcomaPredTumorVsNormalTissue,
  sarcomaPredTumorVsNormalBlood = sarcomaPredTumorVsNormalBlood,
  skincutaneousmelanomaPredTumorVsNormalBlood = skincutaneousmelanomaPredTumorVsNormalBlood,
  skincutaneousmelanomaPredMetastaticTumorVsPrimaryTumor = skincutaneousmelanomaPredMetastaticTumorVsPrimaryTumor,
  stomachPredTumorVsNormalTissue = stomachPredTumorVsNormalTissue,
  stomachPredTumorVsNormalBlood = stomachPredTumorVsNormalBlood,
  testicularPredTumorVsProstateAdenoNormalBlood = testicularPredTumorVsProstateAdenoNormalBlood,
  thymomaPredTumorVsThyroidCarcinomaNormalTissue = thymomaPredTumorVsThyroidCarcinomaNormalTissue,
  thyroidcarcinomaPredTumorVsNormalTissue = thyroidcarcinomaPredTumorVsNormalTissue,
  thyroidcarcinomaPredTumorVsNormalBlood = thyroidcarcinomaPredTumorVsNormalBlood,
  uterinecarcinosarcomaPredTumorVsEndometrialCarcinomaNormalTissue = uterinecarcinosarcomaPredTumorVsEndometrialCarcinomaNormalTissue,
  endometrialcarcinomaPredTumorVsNormalTissue = endometrialcarcinomaPredTumorVsNormalTissue,
  endometrialcarcinomaPredTumorVsNormalBlood = endometrialcarcinomaPredTumorVsNormalBlood,
  uvealmelanomaPredTumorVsNormalBlood = uvealmelanomaPredTumorVsNormalBlood
)

mlGridListTumorVsNormalMicrobes <- list(
  defaultGBMGrid = defaultGBMGrid,
  bicGBMGrid = bicGBMGrid,
  cervGBMGrid = cervGBMGrid,
  cholGBMGrid = cholGBMGrid,
  # gbmGBMGrid = gbmGBMGrid,
  dlbclGBMGrid = dlbclGBMGrid,
  ovGBMGrid = ovGBMGrid,
  padGBMGrid = padGBMGrid,
  pheoGBMGrid = pheoGBMGrid,
  sarcomaGBMGrid = sarcomaGBMGrid
)