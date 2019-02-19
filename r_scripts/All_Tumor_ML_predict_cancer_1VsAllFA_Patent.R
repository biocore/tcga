# All_Tumor_sequence_center_effects.R
# Author: Greg Poore
# Date: Aug 21, 2018
# Purpose: To explore the contribution of sequence center to cancer microbiome profiles

#-------------------------------#
# Load dependencies
require(devtools)
require(doMC)
require(tibble)
require(gbm)

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
# load("amlVbDataAndMetadataAndSNMByGender.RData")
load("tcgaVbDataAndMetadataAndSNM.RData")
load("snmDataSampleTypeWithExpStrategyFINAL.RData")

#------------------------------------------------------
# Predict covariates using ELN / RF while accounting for class imbalances

## Goal: Predict tumor vs normal


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

# require(pheatmap)
# modelTest <- brcaPredCancer1VsAllBDN
# stringTest <- "Breast Invasive Carcinoma"
# tmpMeta <- metadataSamplesAllQC[metadataSamplesAllQC$sample_type == "Blood Derived Normal",]
# tmpData <- snmDataSampleTypeWithExpStrategy[rownames(tmpMeta),modelTest$modelNonzeroVariableImp$Taxa]
# tmpLabel <- data.frame(factor(ifelse(tmpMeta$disease_type == stringTest,
#                                      yes = stringTest, no = "Other")))
# rownames(tmpLabel) <- rownames(tmpMeta)
# 
# tmpLabel <- data.frame(factor(tmpMeta$disease_type))
# colnames(tmpLabel)[1] <- "diseaseType"
# rownames(tmpLabel) <- rownames(tmpMeta)
# pheatmap(tmpData,
#          annotation_row = tmpLabel,
#          show_rownames = FALSE,
#          show_colnames = FALSE)
# 
# tmpDF <- data.frame(Vals = tmpData[,modelTest$modelNonzeroVariableImp$Taxa[1]], tmpLabel = tmpLabel)
# colnames(tmpDF)[2] <- "Label"
# ggplot(tmpDF, aes(x = Vals, fill = Label)) + geom_density()


#---------- Predict tumor vs normal using cancer microbiome in Acute Myeloid Leukemia ----------#

## NB: AML is not in the QC metadata
# # 258 T / 1866 total BDN
# source("predTumorVsNormalAML.R")
# amlPredCancerVsBDN <- predTumorVsNormalAML(qcMLMetadata = metadataSamplesAllQCAML,
#                                                      qcMLDataSNM = snmDataGenderWithAML,
#                                                      cancerTypeString = "Acute Myeloid Leukemia",
#                                                      nlTissueComparison = "Blood Derived Normal",
#                                                      tumorTissueComparison = "Primary Blood Derived Cancer - Peripheral Blood",
#                                                      flagForAlternativeNormal = FALSE,
#                                                      cancerTypeStringAlternativeNormal = NULL,
#                                                      caretModel = "gbm",
#                                                      numResampleIter = 2,
#                                                      numKFold = 3,
#                                                      trainSetProp = 0.7,
#                                                      caretTuneGrid = defaultGBMGrid,
#                                                      plotPath = "./roc-ggplots-predicting-aml-vs-all-blood")

# save(amlPredCancerVsBDN, file = "amlPredCancerVsBDNModelK3R2.RData")



source("predCancer1VsAllFA_Patent.R")
#---------- Calculate DE microbes in Adrenocortical Carcinoma ----------#
accPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Adrenocortical Carcinoma",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Bladder Urothelial Carcinoma ----------#
bucPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                          qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                          cancerTypeString = "Bladder Urothelial Carcinoma",
                                          sampleTypeComparison = "Primary Tumor",
                                          caretModel = "gbm",
                                          numResampleIter = 1,
                                          numKFold = 2,
                                          trainSetProp = 0.7,
                                          caretTuneGrid = defaultGBMGrid,
                                          plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

# glmnetTuneGrid = expand.grid(.alpha = seq(.5, 1, length = 15), .lambda = seq(0.001,0.1,10))
bucPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Bladder Urothelial Carcinoma",
                                            sampleTypeComparison = "Blood Derived Normal",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")


#---------- Calculate DE microbes in Brain Lower Grade Glioma ----------#
blggPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                          qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                          cancerTypeString = "Brain Lower Grade Glioma",
                                          sampleTypeComparison = "Primary Tumor",
                                          caretModel = "gbm",
                                          numResampleIter = 1,
                                          numKFold = 2,
                                          trainSetProp = 0.7,
                                          caretTuneGrid = defaultGBMGrid,
                                          plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")
# source("predCancer1VsAllFA_Patent.R")
blggPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Brain Lower Grade Glioma",
                                           sampleTypeComparison = "Blood Derived Normal",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = customGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Breast Invasive Carcinoma ----------#
brcaPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Breast Invasive Carcinoma",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

brcaPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Breast Invasive Carcinoma",
                                            sampleTypeComparison = "Blood Derived Normal",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma ----------#
cervPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

cervPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",
                                            sampleTypeComparison = "Blood Derived Normal",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Cholangiocarcinoma ----------#
cholPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Cholangiocarcinoma",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Colon Adenocarcinoma ----------#
coadPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Colon Adenocarcinoma",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

coadPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Colon Adenocarcinoma",
                                            sampleTypeComparison = "Blood Derived Normal",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Esophageal Carcinoma ----------#
esophPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Esophageal Carcinoma",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

esophPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Esophageal Carcinoma",
                                            sampleTypeComparison = "Blood Derived Normal",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.5,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Glioblastoma Multiforme ----------#
gbmPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Glioblastoma Multiforme",
                                            sampleTypeComparison = "Primary Tumor",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")
source("predCancer1VsAllFA_Patent.R")
gbmPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                              qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                              cancerTypeString = "Glioblastoma Multiforme",
                                              sampleTypeComparison = "Blood Derived Normal",
                                              caretModel = "gbm",
                                              numResampleIter = 1,
                                              numKFold = 2,
                                              trainSetProp = 0.7,
                                              caretTuneGrid = defaultGBMGrid,
                                              plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Head and Neck Squamous Cell Carcinoma ----------#
hnsccPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                          qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                          cancerTypeString = "Head and Neck Squamous Cell Carcinoma",
                                          sampleTypeComparison = "Primary Tumor",
                                          caretModel = "gbm",
                                          numResampleIter = 1,
                                          numKFold = 2,
                                          trainSetProp = 0.7,
                                          caretTuneGrid = defaultGBMGrid,
                                          plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

hnsccPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Head and Neck Squamous Cell Carcinoma",
                                           sampleTypeComparison = "Blood Derived Normal",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Kidney Chromophobe ----------#
kidneychromophobePredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Kidney Chromophobe",
                                            sampleTypeComparison = "Primary Tumor",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

# Only 9 blood derived normals available for Kidney Chromophobe
kidneychromophobePredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                             qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                             cancerTypeString = "Kidney Chromophobe",
                                             sampleTypeComparison = "Blood Derived Normal",
                                             caretModel = "gbm",
                                             numResampleIter = 1,
                                             numKFold = 2,
                                             trainSetProp = 0.7,
                                             caretTuneGrid = customGBMGrid,
                                             plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")


#---------- Calculate DE microbes in Kidney Renal Clear Cell Carcinoma ----------#
rccPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                        qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                        cancerTypeString = "Kidney Renal Clear Cell Carcinoma",
                                                        sampleTypeComparison = "Primary Tumor",
                                                        caretModel = "gbm",
                                                        numResampleIter = 1,
                                                        numKFold = 2,
                                                        trainSetProp = 0.7,
                                                        caretTuneGrid = defaultGBMGrid,
                                                        plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

# Only 5 blood derived normals available for RCC
rccPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                         qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                         cancerTypeString = "Kidney Renal Clear Cell Carcinoma",
                                                         sampleTypeComparison = "Blood Derived Normal",
                                                         caretModel = "gbm",
                                                         numResampleIter = 1,
                                                         numKFold = 2,
                                                         trainSetProp = 0.7,
                                                         caretTuneGrid = customGBMGrid,
                                                         plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Kidney Renal Papillary Cell Carcinoma ----------#
rpccPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                          qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                          cancerTypeString = "Kidney Renal Papillary Cell Carcinoma",
                                          sampleTypeComparison = "Primary Tumor",
                                          caretModel = "gbm",
                                          numResampleIter = 1,
                                          numKFold = 2,
                                          trainSetProp = 0.7,
                                          caretTuneGrid = defaultGBMGrid,
                                          plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

# NB: Only 35 blood derived normals available for Kidney Renal Papillary Cell Carcinoma
rpccPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Kidney Renal Papillary Cell Carcinoma",
                                           sampleTypeComparison = "Blood Derived Normal",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = customGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Liver Hepatocellular Carcinoma ----------#
hccPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Liver Hepatocellular Carcinoma",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

# NB: Only 32 blood derived normals available for HCC
hccPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Liver Hepatocellular Carcinoma",
                                            sampleTypeComparison = "Blood Derived Normal",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = customGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Lung Adenocarcinoma ----------#
luadPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                          qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                          cancerTypeString = "Lung Adenocarcinoma",
                                          sampleTypeComparison = "Primary Tumor",
                                          caretModel = "gbm",
                                          numResampleIter = 1,
                                          numKFold = 2,
                                          trainSetProp = 0.7,
                                          caretTuneGrid = defaultGBMGrid,
                                          plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

luadPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Lung Adenocarcinoma",
                                           sampleTypeComparison = "Blood Derived Normal",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Lung Squamous Cell Carcinoma ----------#
lsccPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                          qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                          cancerTypeString = "Lung Squamous Cell Carcinoma",
                                          sampleTypeComparison = "Primary Tumor",
                                          caretModel = "gbm",
                                          numResampleIter = 1,
                                          numKFold = 2,
                                          trainSetProp = 0.7,
                                          caretTuneGrid = defaultGBMGrid,
                                          plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

# NB: Only 18 blood derived normals available for LSCC
lsccPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Lung Squamous Cell Carcinoma",
                                           sampleTypeComparison = "Blood Derived Normal",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = customGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Lymphoid Neoplasm Diffuse Large B-cell Lymphoma ----------#
dlbclPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

# NB: Only 7 blood derived normals available for DLBCL
dlbclPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
                                            sampleTypeComparison = "Blood Derived Normal",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = customGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Mesothelioma ----------#
mesoPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Mesothelioma",
                                            sampleTypeComparison = "Primary Tumor",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Ovarian Serous Cystadenocarcinoma ----------#
ovPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Ovarian Serous Cystadenocarcinoma",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = customGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

ovPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Ovarian Serous Cystadenocarcinoma",
                                            sampleTypeComparison = "Blood Derived Normal",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Pancreatic Adenocarcinoma ----------#
padPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                         qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                         cancerTypeString = "Pancreatic Adenocarcinoma",
                                         sampleTypeComparison = "Primary Tumor",
                                         caretModel = "gbm",
                                         numResampleIter = 1,
                                         numKFold = 2,
                                         trainSetProp = 0.7,
                                         caretTuneGrid = defaultGBMGrid,
                                         plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Pheochromocytoma and Paraganglioma ----------#
pheoPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                          qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                          cancerTypeString = "Pheochromocytoma and Paraganglioma",
                                          sampleTypeComparison = "Primary Tumor",
                                          caretModel = "gbm",
                                          numResampleIter = 1,
                                          numKFold = 2,
                                          trainSetProp = 0.7,
                                          caretTuneGrid = defaultGBMGrid,
                                          plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Prostate Adenocarcinoma ----------#
prostatePredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                         qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                         cancerTypeString = "Prostate Adenocarcinoma",
                                         sampleTypeComparison = "Primary Tumor",
                                         caretModel = "gbm",
                                         numResampleIter = 1,
                                         numKFold = 2,
                                         trainSetProp = 0.7,
                                         caretTuneGrid = defaultGBMGrid,
                                         plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

prostatePredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                          qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                          cancerTypeString = "Prostate Adenocarcinoma",
                                          sampleTypeComparison = "Blood Derived Normal",
                                          caretModel = "gbm",
                                          numResampleIter = 1,
                                          numKFold = 2,
                                          trainSetProp = 0.7,
                                          caretTuneGrid = defaultGBMGrid,
                                          plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Rectum Adenocarcinoma ----------#
rectumPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                               qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                               cancerTypeString = "Rectum Adenocarcinoma",
                                               sampleTypeComparison = "Primary Tumor",
                                               caretModel = "gbm",
                                               numResampleIter = 1,
                                               numKFold = 2,
                                               trainSetProp = 0.7,
                                               caretTuneGrid = defaultGBMGrid,
                                               plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

rectumPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                cancerTypeString = "Rectum Adenocarcinoma",
                                                sampleTypeComparison = "Blood Derived Normal",
                                                caretModel = "gbm",
                                                numResampleIter = 1,
                                                numKFold = 2,
                                                trainSetProp = 0.7,
                                                caretTuneGrid = defaultGBMGrid,
                                                plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Sarcoma ----------#
sarcomaPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                             qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                             cancerTypeString = "Sarcoma",
                                             sampleTypeComparison = "Primary Tumor",
                                             caretModel = "gbm",
                                             numResampleIter = 1,
                                             numKFold = 2,
                                             trainSetProp = 0.7,
                                             caretTuneGrid = defaultGBMGrid,
                                             plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

# NB: Only 36 blood derived normals exist for Sarcoma
sarcomaPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                              qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                              cancerTypeString = "Sarcoma",
                                              sampleTypeComparison = "Blood Derived Normal",
                                              caretModel = "gbm",
                                              numResampleIter = 1,
                                              numKFold = 2,
                                              trainSetProp = 0.7,
                                              caretTuneGrid = customGBMGrid,
                                              plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Skin Cutaneous Melanoma ----------#
cutaneousmelanomaPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                             qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                             cancerTypeString = "Skin Cutaneous Melanoma",
                                             sampleTypeComparison = "Primary Tumor",
                                             caretModel = "gbm",
                                             numResampleIter = 1,
                                             numKFold = 2,
                                             trainSetProp = 0.7,
                                             caretTuneGrid = defaultGBMGrid,
                                             plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

cutaneousmelanomaPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                              qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                              cancerTypeString = "Skin Cutaneous Melanoma",
                                              sampleTypeComparison = "Blood Derived Normal",
                                              caretModel = "gbm",
                                              numResampleIter = 1,
                                              numKFold = 2,
                                              trainSetProp = 0.7,
                                              caretTuneGrid = defaultGBMGrid,
                                              plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

cutaneousmelanomaPredCancer1VsAllMet <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                         qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                         cancerTypeString = "Skin Cutaneous Melanoma",
                                                         sampleTypeComparison = "Metastatic",
                                                         caretModel = "gbm",
                                                         numResampleIter = 1,
                                                         numKFold = 2,
                                                         trainSetProp = 0.7,
                                                         caretTuneGrid = defaultGBMGrid,
                                                         plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Stomach Adenocarcinoma ----------#
stadPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                        qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                        cancerTypeString = "Stomach Adenocarcinoma",
                                                        sampleTypeComparison = "Primary Tumor",
                                                        caretModel = "gbm",
                                                        numResampleIter = 1,
                                                        numKFold = 2,
                                                        trainSetProp = 0.7,
                                                        caretTuneGrid = defaultGBMGrid,
                                                        plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

stadPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                         qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                         cancerTypeString = "Stomach Adenocarcinoma",
                                                         sampleTypeComparison = "Blood Derived Normal",
                                                         caretModel = "gbm",
                                                         numResampleIter = 1,
                                                         numKFold = 2,
                                                         trainSetProp = 0.7,
                                                         caretTuneGrid = defaultGBMGrid,
                                                         plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Testicular Germ Cell Tumors ----------#
testiculargermcellPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Testicular Germ Cell Tumors",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Thymoma ----------#
thymomaPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                         qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                         cancerTypeString = "Thymoma",
                                                         sampleTypeComparison = "Primary Tumor",
                                                         caretModel = "gbm",
                                                         numResampleIter = 1,
                                                         numKFold = 2,
                                                         trainSetProp = 0.7,
                                                         caretTuneGrid = defaultGBMGrid,
                                                         plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Thyroid Carcinoma ----------#
thyroidcarcinomaPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                           cancerTypeString = "Thyroid Carcinoma",
                                           sampleTypeComparison = "Primary Tumor",
                                           caretModel = "gbm",
                                           numResampleIter = 1,
                                           numKFold = 2,
                                           trainSetProp = 0.7,
                                           caretTuneGrid = defaultGBMGrid,
                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

thyroidcarcinomaPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Thyroid Carcinoma",
                                            sampleTypeComparison = "Blood Derived Normal",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 2,
                                            trainSetProp = 0.7,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Uterine Carcinosarcoma ----------#
uterinecarcinosarcomaPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                       qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                       cancerTypeString = "Uterine Carcinosarcoma",
                                                       sampleTypeComparison = "Primary Tumor",
                                                       caretModel = "gbm",
                                                       numResampleIter = 1,
                                                       numKFold = 2,
                                                       trainSetProp = 0.7,
                                                       caretTuneGrid = defaultGBMGrid,
                                                       plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Uterine Corpus Endometrial Carcinoma ----------#
endometrialcarcinomaPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                       qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                       cancerTypeString = "Uterine Corpus Endometrial Carcinoma",
                                                       sampleTypeComparison = "Primary Tumor",
                                                       caretModel = "gbm",
                                                       numResampleIter = 1,
                                                       numKFold = 2,
                                                       trainSetProp = 0.7,
                                                       caretTuneGrid = defaultGBMGrid,
                                                       plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

endometrialcarcinomaPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                        qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                        cancerTypeString = "Uterine Corpus Endometrial Carcinoma",
                                                        sampleTypeComparison = "Blood Derived Normal",
                                                        caretModel = "gbm",
                                                        numResampleIter = 1,
                                                        numKFold = 2,
                                                        trainSetProp = 0.7,
                                                        caretTuneGrid = defaultGBMGrid,
                                                        plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

#---------- Calculate DE microbes in Uveal Melanoma ----------#
uvealmelanomaPredCancer1VsAllPT <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                           qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                           cancerTypeString = "Uveal Melanoma",
                                                           sampleTypeComparison = "Primary Tumor",
                                                           caretModel = "gbm",
                                                           numResampleIter = 1,
                                                           numKFold = 2,
                                                           trainSetProp = 0.7,
                                                           caretTuneGrid = defaultGBMGrid,
                                                           plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")

uvealmelanomaPredCancer1VsAllBDN <- predCancer1VsAllFA_Patent(qcMLMetadata = metadataSamplesAllQC,
                                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                            cancerTypeString = "Uveal Melanoma",
                                                            sampleTypeComparison = "Blood Derived Normal",
                                                            caretModel = "gbm",
                                                            numResampleIter = 1,
                                                            numKFold = 2,
                                                            trainSetProp = 0.7,
                                                            caretTuneGrid = defaultGBMGrid,
                                                            plotPath = "./provisional-patent-roc-predicting-cancer-type-1-vs-all")


