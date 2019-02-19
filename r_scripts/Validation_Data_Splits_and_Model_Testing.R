# All_Tumor_sequence_center_effects.R
# Author: Greg Poore
# Date: Aug 21, 2018
# Purpose: To explore the contribution of sequence center to cancer microbiome profiles

#-------------------------------#
# Load dependencies
require(devtools)
require(ggbiplot)
require(doMC)
# require(plotly)

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
vbDataBarnDFReconciled <- read.csv(file = 'vbDataBarnDFReconciled_18116_files.csv',
                                   row.names = 1,
                                   stringsAsFactors = FALSE)
metadataSamplesAll <- read.csv(file = 'tcgaMetadataPickledReconciled_18116_files.csv',
                               header = TRUE,
                               row.names = 2,
                               stringsAsFactors = FALSE,
                               strip.white = TRUE,
                               na.strings=c(""," ","NA"))

cols2Factor <- c("disease_type",
                 "ethnicity",
                 "experimental_strategy",
                 "gender",
                 "investigation",
                 "platform",
                 "primary_site",
                 "race",
                 "reference_genome",
                 "sample_type",
                 "vital_status",
                 "tissue_source_site_label",
                 "data_submitting_center_label",
                 "country_of_sample_procurement",
                 "histological_diagnosis_label",
                 "pathologic_t_label",
                 "pathologic_n_label",
                 "pathologic_stage_label",
                 "icd03_histology_label",
                 "icd03_histology_site",
                 "portion_is_ffpe",
                 "new_tumor_event_after_initial_trtmt",
                 "primary_therapy_outcome_success_label",
                 "analyte_type_label")

# Convert categorical labels to factors
metadataSamplesAll[cols2Factor] <- lapply(metadataSamplesAll[cols2Factor], factor)
metadataSamplesAll$sample_type <- relevel(metadataSamplesAll$sample_type, ref = "Solid Tissue Normal")

## Examine metadata for samples with very poor quality entries (i.e. lots of NA entries)
sort(sapply(metadataSamplesAll, function(x) sum(is.na(x)))) # <---- shows 79 very poor entries

# Remove samples with very poor metadata entry
metadataSamplesAllQC <- droplevels(metadataSamplesAll[! (is.na(metadataSamplesAll$race) | 
                                                is.na(metadataSamplesAll$portion_is_ffpe) | 
                                                is.na(metadataSamplesAll$icd10) |
                                                is.na(metadataSamplesAll$analyte_amount) |
                                                is.na(metadataSamplesAll$age_at_diagnosis)),])
# metadataSamplesAllQCAML <- droplevels(metadataSamplesAll[! (is.na(metadataSamplesAll$race) | 
#                                                            is.na(metadataSamplesAll$portion_is_ffpe) |
#                                                            is.na(metadataSamplesAll$age_at_diagnosis)),])
# metadataSamplesAllQCSurvival <- droplevels(metadataSamplesAll[! (is.na(metadataSamplesAll$race) | 
#                                                         is.na(metadataSamplesAll$portion_is_ffpe) | 
#                                                         is.na(metadataSamplesAll$icd10) |
#                                                         is.na(metadataSamplesAll$analyte_amount) |
#                                                         is.na(metadataSamplesAll$age_at_diagnosis) | 
#                                                         is.na(metadataSamplesAll$days_to_death)),])


sort(sapply(metadataSamplesAllQC, function(x) sum(is.na(x)))) # <---- recheck
# sort(sapply(metadataSamplesAllQCSurvival, function(x) sum(is.na(x)))) # <---- recheck

vbDataBarnDFReconciledQC <- vbDataBarnDFReconciled[rownames(metadataSamplesAllQC),]
# vbDataBarnDFReconciledQCSurvival <- vbDataBarnDFReconciled[rownames(metadataSamplesAllQCSurvival),]

# Look at the metadata for counts
# sample_type, pathologic_stage_label, new_tumor_event_after_initial_trtmt
table(metadataSamplesAll$data_submitting_center_label)
table(metadataSamplesAllQC$data_submitting_center_label)
table(metadataSamplesAll$disease_type)
table(metadataSamplesAllQC$disease_type)

# Plot distributions
ggplot(as.data.frame(table(metadataSamplesAllQC$data_submitting_center_label)),
       aes(x = Var1, y = Freq)) + geom_bar(stat = "identity")

#-----------------------------------------------------------------#
## Splitting of data                     ##
#-----------------------------------------------------------------#
require(data.table)
require(splitstackshape)
require(ggplot2)
require(cowplot)
require(ggsci)
stratSamplingMetadataQC <- stratified(metadataSamplesAllQC, 
                                     group = c("data_submitting_center_label","sample_type","disease_type"), 
                                     size = 0.5, 
                                     keep.rownames = TRUE, 
                                     bothSets = TRUE, 
                                     replace = FALSE)
split1 <- as.data.frame(stratSamplingMetadataQC[[1]])
rownames(split1) <- split1$rn
split2 <- as.data.frame(stratSamplingMetadataQC[[2]])
rownames(split2) <- split2$rn

seqCenterCompare <- data.frame(Split = c(rep("Split 1", dim(split1)[1]),
                                         rep("Split 2", dim(split2)[1])),
                               SeqCenter = c(as.character(split1$data_submitting_center_label),
                                             as.character(split2$data_submitting_center_label)))
ggplot(seqCenterCompare, aes(x = SeqCenter, fill = Split)) + geom_bar(position = "dodge",stat = "count") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) + scale_fill_nejm() +
  geom_text(stat='count', aes(label=..count..), vjust=-1, position = position_dodge(0.9)) + 
  labs(x = "Sequencing Center", y = "Count", title = "Validation split - Sequence center distribution")

sampleTypeCompare <- data.frame(Split = c(rep("Split 1", dim(split1)[1]),
                                         rep("Split 2", dim(split2)[1])),
                               SampleType = c(as.character(split1$sample_type),
                                             as.character(split2$sample_type)))
ggplot(sampleTypeCompare, aes(x = SampleType, fill = Split)) + geom_bar(position = "dodge",stat = "count") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) + scale_fill_nejm() +
  geom_text(stat='count', aes(label=..count..), vjust=-1, position = position_dodge(0.9)) + 
  labs(x = "Sample Type", y = "Count", title = "Validation split - Sample type distribution")

diseaseTypeCompare <- data.frame(Split = c(rep("Split 1", dim(split1)[1]),
                                          rep("Split 2", dim(split2)[1])),
                                DiseaseType = c(as.character(split1$investigation),
                                               as.character(split2$investigation)))
ggplot(diseaseTypeCompare, aes(x = DiseaseType, fill = Split)) + geom_bar(position = "dodge",stat = "count") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_fill_nejm() +
  geom_text(stat='count', aes(label=..count..), hjust=-0.25, position = position_dodge(0.9), angle = 90) + 
  labs(x = "Disease Type", y = "Count", title = "Validation split - Disease type distribution")

# Subset vb data using splits

split1MetadataQC <- split1
split2MetadataQC <- split2

split1VbDataQC <- vbDataBarnDFReconciledQC[rownames(split1MetadataQC),]
split2VbDataQC <- vbDataBarnDFReconciledQC[rownames(split2MetadataQC),]
save(split1MetadataQC, split2MetadataQC, split1VbDataQC, split2VbDataQC, 
     file = "validationSplitVbDataAndMetadata.RData")

## *** See Jupyter Notebook titled "TCGA Validation Spits Batch Correction -- Final Analysis"
# for batch effect correction ***

# Load batch effect corrected data
load("validationSplits_snmDataSampleTypeWithExpStrategyFAS1S2.RData")
load("tcgaVbDataAndMetadataAndSNM.RData")
load("snmDataSampleTypeWithExpStrategyFINAL.RData")



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

source("predTumorVsNormalFA.R")
lihcPredTumorVsNormalTissue <- predTumorVsNormalFA(qcMLMetadata = metadataSamplesAllQC, 
                                                  qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                                  cancerTypeString = "Liver Hepatocellular Carcinoma",
                                                  nlTissueComparison = "Solid Tissue Normal",
                                                  tumorTissueComparison = "Primary Tumor", 
                                                  flagForAlternativeNormal = FALSE,
                                                  cancerTypeStringAlternativeNormal = NULL,
                                                  caretModel = "gbm",
                                                  samplingStrategy = "up",
                                                  numResampleIter = 1,
                                                  numKFold = 4,
                                                  trainSetProp = 0.5,
                                                  plotPath = "./roc-ggplots-predicting-tumor-vs-normal",
                                                  caretTuneGrid = customGBMGrid)

source("validationPredTumorVsNormalFA.R")
valLihcPredTumorVsNormalTissue <- validationPredTumorVsNormalFA(qcMLMetadataS1 = split1MetadataQC,
                                                   qcMLMetadataS2 = split2MetadataQC,
                                                   qcMLDataSNMS1 = snmDataSampleTypeWithExpStrategyFAS1,
                                                   qcMLDataSNMS2 = snmDataSampleTypeWithExpStrategyFAS2,
                                                   cancerTypeString = "Liver Hepatocellular Carcinoma",
                                                   nlTissueComparison = "Solid Tissue Normal",
                                                   tumorTissueComparison = "Primary Tumor", 
                                                   flagForAlternativeNormal = FALSE,
                                                   cancerTypeStringAlternativeNormal = NULL,
                                                   caretModel = "gbm",
                                                   samplingStrategy = "up",
                                                   numResampleIter = 1,
                                                   numKFold = 4,
                                                   plotPath = "./roc-ggplots-predicting-tumor-vs-normal",
                                                   originalModel = lihcPredTumorVsNormalTissue,
                                                   caretTuneGrid = customGBMGrid)

source("predCancer1VsAllFA.R")
tmpPredCancer1VsAllBDN <- predCancer1VsAllFA(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Lung Adenocarcinoma",
                                            sampleTypeComparison = "Blood Derived Normal",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 4,
                                            trainSetProp = 0.5,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./roc-ggplots-predicting-cancer-type-1-vs-all")

source("validationPredCancer1VsAllFA.R")
valTmpPredTumorVsNormalTissue <- validationPredCancer1VsAllFA(qcMLMetadataS1 = split1MetadataQC,
                                                                qcMLMetadataS2 = split2MetadataQC,
                                                                qcMLDataSNMS1 = snmDataSampleTypeWithExpStrategyFAS1,
                                                                qcMLDataSNMS2 = snmDataSampleTypeWithExpStrategyFAS2,
                                                                cancerTypeString = "Lung Adenocarcinoma",
                                                                sampleTypeComparison = "Blood Derived Normal",
                                                                caretModel = "gbm",
                                                                samplingStrategy = "up",
                                                                numResampleIter = 1,
                                                                numKFold = 4,
                                                                plotPath = "./roc-ggplots-predicting-cancer-type-1-vs-all",
                                                                originalModel = tmpPredCancer1VsAllBDN,
                                                                caretTuneGrid = defaultGBMGrid)

source("predCancer1VsAllFA.R")
lggPredCancer1VsAllPT <- predCancer1VsAllFA(qcMLMetadata = metadataSamplesAllQC,
                                             qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                             cancerTypeString = "Brain Lower Grade Glioma",
                                             sampleTypeComparison = "Primary Tumor",
                                             caretModel = "gbm",
                                             numResampleIter = 1,
                                             numKFold = 4,
                                             trainSetProp = 0.5,
                                             caretTuneGrid = defaultGBMGrid,
                                             plotPath = "./roc-ggplots-predicting-cancer-type-1-vs-all")

source("validationPredCancer1VsAllFA.R")
valLGGPredCancer1VsAllPT <- validationPredCancer1VsAllFA(qcMLMetadataS1 = split1MetadataQC,
                                                              qcMLMetadataS2 = split2MetadataQC,
                                                              qcMLDataSNMS1 = snmDataSampleTypeWithExpStrategyFAS1,
                                                              qcMLDataSNMS2 = snmDataSampleTypeWithExpStrategyFAS2,
                                                              cancerTypeString = "Brain Lower Grade Glioma",
                                                              sampleTypeComparison = "Primary Tumor",
                                                              caretModel = "gbm",
                                                              samplingStrategy = "up",
                                                              numResampleIter = 1,
                                                              numKFold = 4,
                                                              plotPath = "./roc-ggplots-predicting-cancer-type-1-vs-all",
                                                              originalModel = lggPredCancer1VsAllPT,
                                                              caretTuneGrid = defaultGBMGrid)

source("predCancer1VsAllFA.R")
stadPredCancer1VsAllPT <- predCancer1VsAllFA(qcMLMetadata = metadataSamplesAllQC,
                                            qcMLDataSNM = snmDataSampleTypeWithExpStrategy,
                                            cancerTypeString = "Stomach Adenocarcinoma",
                                            sampleTypeComparison = "Primary Tumor",
                                            caretModel = "gbm",
                                            numResampleIter = 1,
                                            numKFold = 4,
                                            trainSetProp = 0.5,
                                            caretTuneGrid = defaultGBMGrid,
                                            plotPath = "./roc-ggplots-predicting-cancer-type-1-vs-all")

# save(stadPredCancer1VsAllPT, file = "stadPredCancer1VsAllPT.RData")
source("validationPredCancer1VsAllFA.R")
valSTADPredCancer1VsAllPT <- validationPredCancer1VsAllFA(qcMLMetadataS1 = split1MetadataQC,
                                                         qcMLMetadataS2 = split2MetadataQC,
                                                         qcMLDataSNMS1 = snmDataSampleTypeWithExpStrategyFAS1,
                                                         qcMLDataSNMS2 = snmDataSampleTypeWithExpStrategyFAS2,
                                                         cancerTypeString = "Stomach Adenocarcinoma",
                                                         sampleTypeComparison = "Primary Tumor",
                                                         caretModel = "gbm",
                                                         samplingStrategy = "up",
                                                         numResampleIter = 1,
                                                         numKFold = 4,
                                                         plotPath = "./roc-ggplots-predicting-cancer-type-1-vs-all",
                                                         originalModel = stadPredCancer1VsAllPT,
                                                         caretTuneGrid = defaultGBMGrid)



#---------- Predict tumor vs normal using cancer microbiome in Acute Myeloid Leukemia ----------#

## NB: AML is not in the QC metadata
# # 258 T / 77 N
# amlPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
#                                                      qcMLDataSNM = snmDataSampleType,
#                                                      cancerTypeString = "Acute Myeloid Leukemia",
#                                                      nlTissueComparison = "Solid Tissue Normal",
#                                                      tumorTissueComparison = "Primary Blood Derived Cancer - Peripheral Blood", 
#                                                      flagForAlternativeNormal = FALSE,
#                                                      cancerTypeStringAlternativeNormal = NULL,
#                                                      caretModel = "gbm",
#                                                      numResampleIter = 1,
#                                                      numKFold = 2,
#                                                      testSetProp = 0.7,
#                                                      caretTuneGrid = defaultGBMGrid)

#---------- Predict tumor vs normal using cancer microbiome in Adrenocortical Carcinoma ----------#

# 79 T / 110 N (RCC)
adrenocorticalcarcinomaPredTumorVsRCCNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
bladderPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
bladderPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
blggPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
blggPredPrimaryTumorVsRecurrentTumor <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
bicPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
bicPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
bicPredPrimaryTumorVsMetastaticTumor <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
cervicalPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
cervicalPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
cholPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
colonPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
colonPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
esophagealcarcinomaPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
esophagealcarcinomaPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
# gbmPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
gbmPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
hnsccPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
hnsccPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
kidneychromophobePredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
# kidneychromophobePredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
rccPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
# rccPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
rpccPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
rpccPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
hccPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
hccPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
ladPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
ladPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
lsccPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
lsccPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
dlbclPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
mesotheliomaPredTumorVsLADNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
ovarianPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
ovarianPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
padPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
pheoPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
prostatePredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
prostatePredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
rectumPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
rectumPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
sarcomaPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
sarcomaPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
skincutaneousmelanomaPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
skincutaneousmelanomaPredMetastaticTumorVsPrimaryTumor <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
stomachPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
stomachPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
testicularPredTumorVsProstateAdenoNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
thymomaPredTumorVsThyroidCarcinomaNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
thyroidcarcinomaPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
thyroidcarcinomaPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
uterinecarcinosarcomaPredTumorVsEndometrialCarcinomaNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
endometrialcarcinomaPredTumorVsNormalTissue <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
endometrialcarcinomaPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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
uvealmelanomaPredTumorVsNormalBlood <- predTumorVsNormal(qcMLMetadata = metadataSamplesAllQC, 
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

# save(mlModelListTumorVsNormalMicrobes, file = "mlModelListTumorVsNormalMicrobes.RData")
save(mlGridListTumorVsNormalMicrobes, file = "mlGridListTumorVsNormalMicrobes.RData")

# final <- data.frame(actual = testY,
#                     predict(model_rf, newdata = testX, type = "prob"))

## Goal: Predict tumor type using tumor-only data
mlDataY <- subsetMetadata4CancersQC[subsetMetadata4CancersQC$sample_type == "Primary Tumor",]
mlDataX <- snmData[rownames(mlDataY),]
dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check

# Examine imbalances
table(mlDataY$disease_type) # Cervical: 374 // Lung: 710 // Ovarian: 917 // Stomach: 862

# Create data partitions
set.seed(42)
index <- createDataPartition(mlDataY$disease_type, p = 0.7, list = FALSE)
trainX <- mlDataX[index,]
trainY <- mlDataY[index,]$disease_type
testX <- mlDataX[-index,]
testY <- mlDataY[-index,]$disease_type

set.seed(42)
modelRF <- train(x = trainX,
                 y = trainY,
                 method = "rf",
                 preProcess = c("scale","center"),
                 trControl = trainControl(method = "repeatedcv",
                                          number = 10,
                                          repeats = 10,
                                          verboseIter = TRUE))
final <- data.frame(actual = testY,
                    predict(model_rf, newdata = testX, type = "prob"))



#------------------------------------------------------


# calcDEMicrobesQuantile <- function(countData, allMetadata, adenoRows, outputFolderName, nlTissueComparison){
#   require(limma)
#   require(Glimma)
#   
#   cancerMetadata <- droplevels(allMetadata[(nlTissueComparison & adenoRows),])
#   
#   covDesignShort <- model.matrix(~0 + sample_type + data_submitting_center_label,
#                                  data = cancerMetadata)
#   
#   colnames(covDesignShort) <- gsub('([[:punct:]])|\\s+','_',colnames(covDesignShort))
#   colnames(covDesignShort)[1:2] <- c("Normal_Tissue","Primary_Tumor")
#   
#   counts <- t(countData[rownames(cancerMetadata),])
#   
#   v <- voom(counts = counts, design = covDesignShort, plot = TRUE, normalize.method = "quantile")
#   vfit <- lmFit(v, covDesignShort)
#   vfit <- eBayes(vfit)
#   contrast.matrix <- makeContrasts(Primary_Tumor - Normal_Tissue, levels = covDesignShort)
#   vfit2 <- contrasts.fit(vfit, contrasts = contrast.matrix)
#   vfit2 <- eBayes(vfit2)
#   dt <- decideTests(vfit2)
# 
#   require(Glimma)
#   glMDPlot(vfit2, coef = 1, status = dt,
#            groups = cancerMetadata$sample_type,
#            counts = v, side.main="Symbols",
#            folder=outputFolderName)
# }
# 
# 
# calcDEMicrobesQuantile(countData = vbBarnFilt4,
#                     allMetadata = subsetMetadata4Cancers,
#                     adenoRows = cervicalAdenoRows,
#                     outputFolderName = "Glimma_Cervical_Quantile",
#                     nlTissueComparison = nlVsTumorRows)

# Separate lung cancer metadata
lungCancerMetadata <- droplevels(subsetMetadata4Cancers[(nlVsTumorRows & lungAdenoRows),])

# Form design matrix: lung primary tumor vs. normal tissue
# covDesign <- model.matrix(~0 + sample_type + 
#                             data_submitting_center_label + 
#                             experimental_strategy +
#                             platform,
#                             data = lungCancerMetadata)
# colnames(covDesign) <- c("Solid_Tissue_Normal",
#                          "Primary_Tumor",
#                          "data_sub_Harvard",
#                          "data_sub_UNC",
#                          "exp_strategy_WGS",
#                          "platform_Illumina_HiSeq")
# head(covDesign)

covDesignShort <- model.matrix(~0 + sample_type + data_submitting_center_label,
                                 data = lungCancerMetadata)
head(covDesignShort)
colnames(covDesignShort) <- c("Solid_Tissue_Normal",
                              "Primary_Tumor",
                              "data_sub_Harvard",
                              "data_sub_UNC")
head(covDesignShort)

# lungSampleType <- factor(lungCancerMetadata$sample_type)
# design <- model.matrix(~0+lungSampleType)
# colnames(design) <- levels(lungSampleType)
# head(design)

# Extract lung cancer data with matching metadata
counts <- t(vbBarnFilt4[rownames(lungCancerMetadata),])

# Normalize using edgeR and then plug into voom
dge <- DGEList(counts = counts)
keep <- filterByExpr(dge, covDesignShort)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
vdge <- voom(dge, design = covDesignShort, plot = TRUE, normalize.method="none")
vdgeFit <- lmFit(vdge, covDesignShort)
vdgeFit <- eBayes(vdgeFit)
contrast.matrix <- makeContrasts(Primary_Tumor - Solid_Tissue_Normal, levels = covDesignShort)
vdgeFit2 <- contrasts.fit(vdgeFit, contrasts = contrast.matrix)
vdgeFit2 <- eBayes(vdgeFit2)
vdgeFitDT <- decideTests(vdgeFit2)
sum(abs(vdgeFitDT))

require(Glimma)
glMDPlot(vdgeFit2, coef = 1, status = vdgeFitDT,
         groups = lungCancerMetadata$sample_type,
         counts = vdge, side.main="Symbols",
         folder="GlimmaTestDGE")

# #  Alternatively, use voom directly while accounting for noisy data
# v <- voom(counts = counts, design = covDesign, plot = TRUE, normalize.method = "quantile")
# vfit <- lmFit(v, covDesign)
# vfit <- eBayes(vfit)
# topTable(vfit, coef=1)
# 
# contrast.matrix <- makeContrasts(Primary_Tumor - Solid_Tissue_Normal, levels = covDesign)
# vfit2 <- contrasts.fit(vfit, contrasts = contrast.matrix)
# vfit2 <- eBayes(vfit2)
# topTable(vfit2, coef=1)


### Testing ###
v <- voom(counts = counts, design = covDesignShort, plot = TRUE, normalize.method = "quantile")
vfit <- lmFit(v, covDesignShort)
vfit <- eBayes(vfit)
tmp <- topTable(vfit, coef=2, p.value = 0.05, number = Inf)
dim(tmp)

contrast.matrix <- makeContrasts(Primary_Tumor - Solid_Tissue_Normal, levels = covDesignShort)
vfit2 <- contrasts.fit(vfit, contrasts = contrast.matrix)
vfit2 <- eBayes(vfit2)
dt <- decideTests(vfit2)
tmp2 <- topTreat(vfit2, coef=1, p.value = 0.05, number = Inf)
dim(tmp2)

require(Glimma)
glMDPlot(vfit2, coef = 1, status = dt,
         groups = lungCancerMetadata$sample_type,
         counts = v, side.main="Symbols",
         folder="GlimmaTest")
# library(gplots)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale="row",
#           labRow=v$genes$SYMBOL[i], labCol=group,
#           col=mycol, trace="none", density.info="none", 
#           margin=c(8,6), lhei=c(2,10), dendrogram="column")

###########################################
###########################################
###########################################

lungCancerMetadataFiltered <- lungCancerMetadata[!is.na(lungCancerMetadata$pathologic_n_label),]

## SNM correction using all 4 batches; supervised
counts <- t(vbBarnFilt4[rownames(lungCancerMetadataFiltered),])
dim(counts)

covDesignShort <- model.matrix(~0 + sample_type +
                                 data_submitting_center_label +
                                 platform +
                                 tissue_source_site_label +
                                 reference_genome +
                                 country_of_sample_procurement +
                                 portion_is_ffpe +
                                 aliquot_concentration +
                                 analyte_amount +
                                 analyte_type_label,
                                 data=lungCancerMetadataFiltered)
dim(covDesignShort)
colnames(covDesignShort)
colnames(covDesignShort) <- gsub('([[:punct:]])|\\s+','_',colnames(covDesignShort))
colnames(covDesignShort)[1:2] <- c("Normal_Tissue","Primary_Tumor")
colnames(covDesignShort)


v <- voom(counts = counts, design = covDesignShort, plot = TRUE, normalize.method = "quantile")


require(snm)
sort(sapply(lungCancerMetadataFiltered, function(x) sum(is.na(x))))
bio.var <- model.matrix(~sample_type +
                          pathologic_stage_label +
                          histological_diagnosis_label,
                          data=lungCancerMetadataFiltered)
adj.var <- model.matrix(~data_submitting_center_label +
                          platform +
                          tissue_source_site_label +
                          portion_is_ffpe,
                          data=lungCancerMetadataFiltered)
# adj.var <- model.matrix(~data_submitting_center_label +
#                           platform +
#                           tissue_source_site_label +
#                           reference_genome +
#                           country_of_sample_procurement +
#                           portion_is_ffpe +
#                           #aliquot_concentration +
#                           #analyte_amount +
#                           analyte_type_label,
#                         data=lungCancerMetadata)
colnames(bio.var) <- gsub('([[:punct:]])|\\s+','_',colnames(bio.var))
colnames(adj.var) <- gsub('([[:punct:]])|\\s+','_',colnames(adj.var))
dim(adj.var)
dim(bio.var)
dim(t(v$E))
dim(covDesignShort)
snmDataObj <- snm(raw.dat = v$E, bio.var = bio.var, adj.var = adj.var, rm.adj=TRUE,
                  verbose = TRUE)
snmData <- t(snmDataObj$norm.dat)

pcaPlotting <- function(pcaObject,pcChoices, dataLabels, factorString, titleString){
  require(ggbiplot)
  theme_update(plot.title = element_text(hjust = 0.5))
  g <- ggbiplot(pcaObject,pcChoices, obs.scale = 1, var.scale = 1,
                groups = dataLabels, ellipse = TRUE,
                circle = TRUE,var.axes=FALSE) + scale_color_discrete(name = factorString) +
        theme_bw() + 
        #theme(legend.direction = "horizontal", legend.position = "top") +
        ggtitle(titleString) + theme(plot.title = element_text(hjust = 0.5))
  
  print(g)
}

unnormalizedPCAPlot <- pcaPlotting(pcaObject = prcomp(t(v$E)),
                                   pcChoices = 1:2,
                                   dataLabels = lungCancerMetadata$data_submitting_center_label,
                                   factorString = "Batch",
                                   titleString = "PCA w/o Batch Correction")

snmPCAPlot <- pcaPlotting(pcaObject = prcomp(snmData),
                          pcChoices = 1:2,
                          dataLabels = lungCancerMetadata$data_submitting_center_label,
                          factorString = "Sequencing Center",
                          titleString = "PCA w/ SNM Correction")

grid.arrange(unnormalizedPCAPlot,snmPCAPlot,
             nrow=1, ncol=2,
             top = textGrob("Sequencing Center Effects on log-normalized data",gp=gpar(fontsize=20,font=3)))

unnormalizedPCAPlotTT <- pcaPlotting(pcaObject = prcomp(t(v$E)),
                                   pcChoices = 1:2,
                                   dataLabels = lungCancerMetadataFiltered$pathologic_stage_label,
                                   factorString = "Disease type",
                                   titleString = "TT: PCA w/o Batch Correction")

snmPCAPlotTT <- pcaPlotting(pcaObject = prcomp(snmData),
                          pcChoices = 1:2,
                          dataLabels = lungCancerMetadataFiltered$pathologic_stage_label,
                          factorString = "Disease type",
                          titleString = "TT: PCA w/ SNM Correction")

grid.arrange(unnormalizedPCAPlotTT,snmPCAPlotTT,
             nrow=1, ncol=2,
             top = textGrob("Effects of batch correction on disease type",gp=gpar(fontsize=20,font=3)))


#-------------------------------------------------------------#
require(glmnet)
# 
# nlVsTumorRows <- (subsetMetadata4Cancers$sample_type == "Solid Tissue Normal" | subsetMetadata4Cancers$sample_type == "Primary Tumor")
# bloodVsTumorRows <- (subsetMetadata4Cancers$sample_type == "Blood Derived Normal" | subsetMetadata4Cancers$sample_type == "Primary Tumor")
# # Extract samples from specific cancers
# lungAdenoRows <- (subsetMetadata4Cancers$disease_type == "Lung Adenocarcinoma")
# cervicalAdenoRows <- (subsetMetadata4Cancers$disease_type == "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma")
# ovarianAdenoRows <- (subsetMetadata4Cancers$disease_type == "Ovarian Serous Cystadenocarcinoma")
# stomachAdenoRows <- (subsetMetadata4Cancers$disease_type == "Stomach Adenocarcinoma")
# 


nlVsTumorRowsLungFiltered <- (lungCancerMetadataFiltered$sample_type == "Solid Tissue Normal" | lungCancerMetadataFiltered$sample_type == "Primary Tumor")
# bloodVsTumorRowsLungFiltered <- (lungCancerMetadataFiltered$sample_type == "Blood Derived Normal" | lungCancerMetadataFiltered$sample_type == "Primary Tumor")

cancerMetadataLungFiltered <- droplevels(lungCancerMetadataFiltered[(nlVsTumorRowsLungFiltered),])
cancerExprdataLungFiltered <- snmData[(nlVsTumorRowsLungFiltered),]

mlDataX <- snmData # dimension (nobs x nvars)
mlDataY <- lungCancerMetadataFiltered$sample_type

tumorVsNormalModel <- cv.glmnet(x= mlDataX, y = mlDataY,
                             family = "binomial",
                             nfolds = 10, # 10-fold CV
                             alpha = 1, #alpha = 0.5 for elastic net
                             type.measure="class",
                             parallel = T)
plot(tumorVsNormalModel)

# Contributing non-zero coefs
sigBactVsBaseProbes <- coef(bactVsBaseModel, s="lambda.1se")[which(coef(bactVsBaseModel)!=0),]
bactVsBaseIDs <- names(sigBactVsBaseProbes[-1])
bactVsBaseIDs