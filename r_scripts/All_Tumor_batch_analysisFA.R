# All_Tumor_sequence_center_effects.R
# Author: Greg Poore
# Date: Aug 21, 2018
# Purpose: To explore the contribution of sequence center to cancer microbiome profiles

#-------------------------------#
# Load dependencies
require(devtools)
require(ggbiplot)
require(doMC)

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
# load("snmObjectsAllCancers.RData")
# load("tcgaVbDataAndMetadataAndSNM.RData")

#-----------------------------------------------------
# Oct 10, 2018 --> PVCA
require(pvca)
require(ggplot2)
require(ggpubr)
require(tibble)
require(reshape2)
require(cowplot)
require(ggsci)

# load("pvcaSampleTypePostSNM.RData")
# load("pvcaSampleTypeVoomNoSNM.RData")
# load("pvcaSampleSNMwithDzType.RData")
# load("pvcaVbRawNoVoomNoSNMwithDzType.RData")
# load("pvcaVoomNoSNMwithDzType.RData")
load("pvcaVbRawNoVoomNoSNM_ExtendedFiltered_FA.RData")
load("pvcaVoomNoSNM_ExtendedFiltered_FA.RData")
load("pvcaSampleWithExpStrategySNM_ExtendedFiltered_FA.RData")

# tmp2 <- as.data.frame(rbind(tmpVoomNoSNM,tmpSampleSNM))
# tmp3 <- rownames_to_column(tmp2, "group")
# data.m <- melt(tmp3, id.vars = "group")
# data.m$group <- factor(data.m$group, levels = c("tmpVoomNoSNM", "tmpSampleSNM"))

tmp2 <- as.data.frame(rbind(pvcaVbRawNoVoomNoSNM_ExtendedFiltered_FA,
                            pvcaVoomNoSNM_ExtendedFiltered_FA,
                            pvcaSampleWithExpStrategySNM_ExtendedFiltered_FA))
tmp3 <- rownames_to_column(tmp2, "group")
tmp3$group <- c("Raw Counts", "Voom", "Voom-SNM")
data.m <- melt(tmp3, id.vars = "group")
data.m$group <- factor(data.m$group, levels = c("Raw Counts",
                                                "Voom",
                                                "Voom-SNM"))

g <- ggplot(data.m, aes(x = variable, y = value, fill = group)) + 
  geom_bar(aes(fill = group), position = "dodge", stat = "identity") + 
  geom_text(aes(label=round(value,3)), position=position_dodge(width=0.9), vjust=-0.25) + 
  labs(x = "Technical & Biological Effects",
       y = "Weighted average proportion variance",
       title = "Principal variance component analysis of batch effect correction procedures") +
  # theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.3, 0.75)) + 
  scale_x_discrete(labels=c("sample_type" = "Sample Type", 
                            "disease_type" = "Disease Type", 
                            "data_submitting_center_label" = "Sequencing Center",
                            "platform" = "Sequencing Platform",
                            "experimental_strategy" = "Experimental strategy",
                            "tissue_source_site_label" = "Tissue Source Site",
                            "portion_is_ffpe" = "FFPE Fixation",
                            "resid" = "Residual \n(not explained by technical variation)")) +
  scale_fill_nejm(name = "Data types", labels = c("Raw Count Data", 
                                                      "Voom Normalized Data",
                                                      "Voom Normalized & SNM Corrected Data"))
g
save_plot(g, filename = "PVCA batch effect correction v2.png",
          base_width = 18,
          units = "in",
          dpi = "retina")



#-----------------------------------------------------


snmDataObjSampleType <- snmObjectsAllCancers$snmDataObjSampleTypeAllCancers

hist(log10(snmDataObjSampleType$pval))

snm.plot(snmDataObjSampleType, col.by=snmDataObjSampleType$bio.var)

tmp.fit <- fitted(snmDataObjSampleType)

# tmp <- lm.fit(x = t(snmDataObjSampleType$raw.dat), y = snmDataSampleType)
# snmDataSampleTypeLMFit <- tmp
# save(snmDataSampleTypeLMFit, file = "snmDataSampleTypeLMFit.RData")

# require(pheatmap)
# pheatmap(snmDataSampleTypeLMFit$coefficients,
#          show_rownames = FALSE,
#          show_colnames = FALSE)
require(matrixStats)
divSNMDataSampleType <- snmDataSampleType / t(snmDataObjSampleType$raw.dat)
taxaMedians <- data.frame(Medians = colMedians(divSNMDataSampleType), 
                          Taxa = colnames(divSNMDataSampleType),
                          pval = factor(ifelse(snmDataObjSampleType$pval <=0.05, yes = "P-value <= 0.05", no = "P-value > 0.05")))
sampleMedians <- data.frame(Medians = rowMedians(divSNMDataSampleType), 
                            Samples = rownames(divSNMDataSampleType),
                            SeqCenter = metadataSamplesAllQC$data_submitting_center_label,
                            SampleType = metadataSamplesAllQC$sample_type,
                            CancerType = metadataSamplesAllQC$disease_type)
g <- ggplot(taxaMedians, aes(x = reorder(Taxa, -Medians), y = Medians, fill = pval)) + 
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
g
gs <- ggplot(sampleMedians, aes(x = reorder(Samples, -Medians), y = Medians, fill = CancerType)) + 
  geom_bar(stat = "identity") + coord_flip() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  scale_y_log10() + labs(y = "Median of Normalizing Ratios Per Sample", x = "Samples", fill='Cancer Type')  
gs
#--------------------------------------------#

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

# Subset metadata and count data for 4 filtered tumor types (lung, cervical, ovarian, stomach)
lungAdenoRows <- (metadataSamplesAll$disease_type == "Lung Adenocarcinoma")
cervicalAdenoRows <- (metadataSamplesAll$disease_type == "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma")
ovarianAdenoRows <- (metadataSamplesAll$disease_type == "Ovarian Serous Cystadenocarcinoma")
stomachAdenoRows <- (metadataSamplesAll$disease_type == "Stomach Adenocarcinoma")
filtered4TumorRows <- (lungAdenoRows | cervicalAdenoRows | ovarianAdenoRows | stomachAdenoRows)

subsetMetadata4Cancers <- metadataSamplesAll[filtered4TumorRows,]
subsetVBData4Cancers <- vbDataBarnDFReconciled[filtered4TumorRows,]

## Examine metadata for samples with very poor quality entries (i.e. lots of NA entries)
sort(sapply(metadataSamplesAll, function(x) sum(is.na(x)))) # <---- shows 79 very poor entries
sort(sapply(subsetMetadata4Cancers, function(x) sum(is.na(x)))) # <---- shows one very poor entry

# Remove samples with very poor metadata entry
metadataSamplesAllQC <- droplevels(metadataSamplesAll[! (is.na(metadataSamplesAll$race) | 
                                                is.na(metadataSamplesAll$portion_is_ffpe) | 
                                                is.na(metadataSamplesAll$icd10) |
                                                is.na(metadataSamplesAll$analyte_amount) |
                                                is.na(metadataSamplesAll$age_at_diagnosis)),])
subsetMetadata4CancersQC <- droplevels(subsetMetadata4Cancers[! (is.na(subsetMetadata4Cancers$race) | 
                                                        is.na(subsetMetadata4Cancers$icd10) |
                                                        is.na(subsetMetadata4Cancers$age_at_diagnosis)),])
metadataSamplesAllQCSurvival <- droplevels(metadataSamplesAll[! (is.na(metadataSamplesAll$race) | 
                                                        is.na(metadataSamplesAll$portion_is_ffpe) | 
                                                        is.na(metadataSamplesAll$icd10) |
                                                        is.na(metadataSamplesAll$analyte_amount) |
                                                        is.na(metadataSamplesAll$age_at_diagnosis) | 
                                                        is.na(metadataSamplesAll$days_to_death)),])


sort(sapply(metadataSamplesAllQC, function(x) sum(is.na(x)))) # <---- recheck
sort(sapply(subsetMetadata4CancersQC, function(x) sum(is.na(x)))) # <---- recheck
sort(sapply(metadataSamplesAllQCSurvival, function(x) sum(is.na(x)))) # <---- recheck

vbDataBarnDFReconciledQC <- vbDataBarnDFReconciled[rownames(metadataSamplesAllQC),]
vbDataBarnDFReconciledQCSurvival <- vbDataBarnDFReconciled[rownames(metadataSamplesAllQCSurvival),]
subsetVBData4CancersQC <- subsetVBData4Cancers[rownames(subsetMetadata4CancersQC),]

# Look at the metadata for counts
# sample_type, pathologic_stage_label, new_tumor_event_after_initial_trtmt
table(metadataSamplesAll$data_submitting_center_label)
table(metadataSamplesAllQC$data_submitting_center_label)
table(metadataSamplesAll$disease_type)
table(metadataSamplesAllQC$disease_type)

## Plot distributions
# ggplot(as.data.frame(table(metadataSamplesAllQC$data_submitting_center_label)),
#        aes(x = Var1, y = Freq)) + geom_bar(stat = "identity")

#-----------------------------------------------------------------#
## Batch correction of voom-normalized data                      ##
#-----------------------------------------------------------------#

require(limma)
require(edgeR)
require(dplyr)
require(snm)

## NB: Need to use QC metadata and associated samples

qcMetadata <- metadataSamplesAllQC
qcData <- vbDataBarnDFReconciledQC

# Set up design matrix
covDesignNorm <- model.matrix(~0 + sample_type +
                                data_submitting_center_label,
                              data = qcMetadata)
colnames(covDesignNorm)
colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
colnames(covDesignNorm)

# Set up counts matrix
counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)

# Normalize using edgeR and then plug into voom
dge <- DGEList(counts = counts)
keep <- filterByExpr(dge, covDesignNorm)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "TMM")
vdge <- voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE, normalize.method="none")

# Apply
bio.var.sample.type <- model.matrix(~sample_type, # histological_diagnosis_label and disease_type tried but cause function to fail
                        data=qcMetadata)
bio.var.disease.type <- model.matrix(~disease_type,
                                    data=qcMetadata)
bio.var.path.stage <- model.matrix(~pathologic_stage_label,
                                     data=qcMetadata)
dim(bio.var.sample.type)
dim(bio.var.disease.type)
dim(bio.var.path.stage)
adj.var <- model.matrix(~data_submitting_center_label +
                          platform +
                          tissue_source_site_label +
                          portion_is_ffpe,
                        data=qcMetadata)
adj.var.disease.type <- model.matrix(~data_submitting_center_label,
                                    data=qcMetadata)
adj.var.path.stage <- model.matrix(~data_submitting_center_label,
                                  data=qcMetadata)
dim(adj.var)
dim(adj.var.disease.type)
dim(adj.var.path.stage)
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
colnames(bio.var.sample.type) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var.sample.type))
colnames(bio.var.disease.type) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var.disease.type))
colnames(bio.var.path.stage) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var.path.stage))
colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
colnames(adj.var.disease.type) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var.disease.type))
colnames(adj.var.path.stage) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var.path.stage))
dim(adj.var)
dim(bio.var.sample.type)
dim(bio.var.disease.type)
dim(bio.var.path.stage)
dim(t(vdge$E))
dim(covDesignNorm)
snmDataObjSampleType <- snm(raw.dat = vdge$E, 
                        bio.var = bio.var.sample.type, 
                        adj.var = adj.var, 
                        rm.adj=TRUE,
                        verbose = TRUE)
snmDataSampleType <- t(snmDataObjSampleType$norm.dat)
# snmDataObjDiseaseType <- snm(raw.dat = vdge$E, 
#                             bio.var = bio.var.disease.type, 
#                             adj.var = adj.var.disease.type, 
#                             rm.adj=TRUE,
#                             num.iter = 2,
#                             lmer.max.iter=50,
#                             nbins = 10,
#                             verbose = TRUE)
# snmDataDiseaseType <- t(snmDataObjDiseaseType$norm.dat)
# snmDataObjPathStage <- snm(raw.dat = vdge$E,
#                              bio.var = bio.var.path.stage,
#                              adj.var = adj.var.path.stage,
#                              rm.adj=TRUE,
#                              lmer.max.iter=100,
#                              nbins = 10,
#                              verbose = TRUE)
# snmDataPathStage <- t(snmDataObjPathStage$norm.dat)

##########

pcaPlotting <- function(pcaObject,pcChoices, dataLabels, factorString, titleString){
  require(ggbiplot)
  require(ggsci)
  theme_update(plot.title = element_text(hjust = 0.5))
  g <- ggbiplot(pcaObject,pcChoices, obs.scale = 1, var.scale = 1,
                groups = dataLabels, ellipse = TRUE,
                circle = TRUE,var.axes=FALSE) + scale_color_discrete(name = factorString) +
    theme_bw() + 
    #theme(legend.direction = "horizontal", legend.position = "top") +
    ggtitle(titleString) + theme(plot.title = element_text(hjust = 0.5))
  
  print(g)
}

tmpMetadata <- droplevels(metadataSamplesAllQC[metadataSamplesAllQC$sample_type == "Primary Tumor",])
tmpData <- snmDataSampleTypeWithExpStrategy[rownames(tmpMetadata),]

tmpPCAPlot <- pcaPlotting(pcaObject = prcomp(tmpData),
                                   pcChoices = 1:2,
                                   dataLabels = tmpMetadata$disease_type,
                                   factorString = "Disease Type",
                                   titleString = "PTs by Disease Type")

unnormalizedPCAPlot <- pcaPlotting(pcaObject = prcomp(t(vdge$E)),
                                   pcChoices = 1:2,
                                   dataLabels = qcMetadata$data_submitting_center_label,
                                   factorString = "Batch",
                                   titleString = "PCA w/o Batch Correction")

snmPCAPlotSampleType <- pcaPlotting(pcaObject = prcomp(snmDataSampleType),
                          pcChoices = 1:2,
                          dataLabels = qcMetadata$data_submitting_center_label,
                          factorString = "Sequencing Center",
                          titleString = "PCA w/ SNM Correction (Target: Sample Type)")

# snmPCAPlotDiseaseType <- pcaPlotting(pcaObject = prcomp(snmDataDiseaseType),
#                           pcChoices = 1:2,
#                           dataLabels = qcMetadata$data_submitting_center_label,
#                           factorString = "Sequencing Center",
#                           titleString = "PCA w/ SNM Correction (Target: Disease Type)")
# 
# snmPCAPlotPathStage <- pcaPlotting(pcaObject = prcomp(snmDataPathStage),
#                                      pcChoices = 1:2,
#                                      dataLabels = qcMetadata$data_submitting_center_label,
#                                      factorString = "Sequencing Center",
#                                      titleString = "PCA w/ SNM Correction (Target: Pathological Stage)")

require(gridExtra)
# grid.arrange(unnormalizedPCAPlot,
#              snmPCAPlotSampleType,
#              snmPCAPlotDiseaseType,
#              nrow=3, ncol=1,
#              top = textGrob("Sequencing Center Effects on voom-normalized data (4 Cancers: LAD, OV, CERV, STAD)",gp=gpar(fontsize=20,font=3)))

grid.arrange(unnormalizedPCAPlot,
             snmPCAPlotSampleType,
             snmPCAPlotDiseaseType,
             snmPCAPlotPathStage,
             nrow=2, ncol=2,
             top = textGrob("Sequencing Center Effects on voom-normalized data (All Cancers)",gp=gpar(fontsize=20,font=3)))

# snmObjects4Cancers <- list(
#   snmDataObjSampleType4Cancers = snmDataObjSampleType,
#   snmDataObjDiseaseType4Cancers = snmDataObjDiseaseType
# )

snmObjectsAllCancers <- list(
  snmDataObjSampleTypeAllCancers = snmDataObjSampleType,
  # snmDataObjDiseaseTypeAllCancers = snmDataObjDiseaseType,
  # snmDataObjPathStageAllCancers = snmDataObjPathStage,
  # snmPCAPlotDiseaseType = snmPCAPlotDiseaseType,
  # snmPCAPlotPathStage = snmPCAPlotPathStage,
  unnormalizedPCAPlot = unnormalizedPCAPlot,
  snmPCAPlotSampleType = snmPCAPlotSampleType
)

# save(snmObjects4Cancers, file = "snmObjects4Cancers.RData")
save(snmObjectsAllCancers, file = "snmObjectsAllCancers_v2.RData")
load("snmObjectsAllCancers.RData")
snmDataSampleType <- t(snmObjectsAllCancers$snmDataObjSampleTypeAllCancers$norm.dat)

# #------------------------------------------------------
# # Account for negative values before exporting to BIOM table
# 
# hist(snmData) # approx. normally distributed
# hist(snmData-min(snmData)) # adding constant merely shifts distribution
# snmDataPositiveOnly <- snmData-min(snmData)
# snmPCAPlotPositiveOnly <- pcaPlotting(pcaObject = prcomp(snmDataPositiveOnly),
#                                       pcChoices = 1:2,
#                                       dataLabels = subsetMetadata4CancersQC$data_submitting_center_label,
#                                       factorString = "Sequencing Center",
#                                       titleString = "PCA w/ SNM Correction")
# snmPCAPlotPositiveOnly # <-- no effect on PCA
# 
# #------------------------------------------------------
# 
# # Save data as BIOM table
# require(biomformat)
# 
# qiimeBIOM_SNM_PositiveOnly <- make_biom(t(snmDataPositiveOnly), sample_metadata = NULL, id = "QIIME 4 tumors SNM", matrix_element_type = "float")
# write_biom(qiimeBIOM_SNM_PositiveOnly, biom_file = "qiime4tumorDataSNM_PositiveOnly.biom")
# write.csv(qiime4tumorMetadata, file = "qiime4tumorMetadata.csv")

#------------------------------------------------------
# Predict covariates using ELN / RF while accounting for class imbalances

## Goal: Predict tumor vs normal
predTumorVsNormal <- function(qcMLMetadata, 
                              qcMLDataSNM,
                              cancerTypeString = NULL,
                              nlTissueComparison = "Solid Tissue Normal",
                              tumorTissueComparison = "Primary Tumor", 
                              flagForAlternativeNormal = FALSE,
                              cancerTypeStringAlternativeNormal = NULL,
                              caretModel = "gbm",
                              numResampleIter = 5,
                              numKFold = 3,
                              testSetProp = 0.7,
                              caretTuneGrid = NULL,
                              ...){
  
  require(caret) # for model building
  require(pROC) # for AUC calculations
  require(DMwR) # for SMOTE class imbalance correction
  require(ROSE) # for ROSE class imbalance correction
  require(purrr) # for functional programming using map()
  require(dplyr) # for data manipulation
  require(doMC) # for parallel computing
  
  numCores <- detectCores()
  registerDoMC(cores=numCores)
  
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
  index <- createDataPartition(mlDataY$sample_type, p = testSetProp, list = FALSE)
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
  
  set.seed(42)
  ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       verboseIter = TRUE,
                       allowParallel=TRUE)
  
  print("Now training model with equally weighted samples...")
  modelVanilla <- train(x = trainX,
                        y = refactoredTrainY,
                        method = caretModel,
                        preProcess = c("scale","center"),
                        metric = "ROC",
                        verbose = TRUE,
                        trControl = ctrl,
                        tuneGrid = caretTuneGrid)
  
  modelVanilla %>%
    test_roc(data = testX, classes = refactoredTestY) %>%
    auc() %>% print()
  
  
  # Fix imbalance by sample weighting
  
  model_weights <- ifelse(refactoredTrainY == gsub('([[:punct:]])|\\s+','',nlTissueComparison),
                          
                          (1/table(refactoredTrainY)[1]) * 0.5,
                          (1/table(refactoredTrainY)[2]) * 0.5
                          
                          )
  sum(model_weights) == 1 # Sanity check
  
  # Use the same seed to ensure same cross-validation splits
  ctrl$seeds <- modelVanilla$control$seeds
  
  # Build weighted model
  print("Now training model with weighted samples...")
  weightedFit <- train(x = trainX,
                       y = refactoredTrainY,
                       method = caretModel,
                       preProcess = c("scale","center"),
                       trControl = ctrl,
                       metric = "ROC",
                       verbose = TRUE,
                       weights = model_weights,
                       tuneGrid = caretTuneGrid)
  
  weightedFit %>%
    test_roc(data = testX, classes = refactoredTestY) %>%
    auc() %>% print()
  
  # Build down-sampled model
  ctrl$sampling <- "down"
  
  print("Now training model with down sampling...")
  downFit <- train(x = trainX,
                   y = refactoredTrainY,
                   method = caretModel,
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   verbose = TRUE,
                   tuneGrid = caretTuneGrid,
                   metric = "ROC",
                   ...)
  
  downFit %>%
    test_roc(data = testX, classes = refactoredTestY) %>%
    auc() %>% print()
  
  # Build up-sampled model
  ctrl$sampling <- "up"
  
  print("Now training model with up sampling...")
  upFit <- train(x = trainX,
                 y = refactoredTrainY,
                 method = caretModel,
                 preProcess = c("scale","center"),
                 trControl = ctrl,
                 verbose = TRUE,
                 metric = "ROC",
                 tuneGrid = caretTuneGrid)
  
  upFit %>%
    test_roc(data = testX, classes = refactoredTestY) %>%
    auc() %>% print()
  
  # Build SMOTE model
  ctrl$sampling <- "smote"
  
  print("Now training model with SMOTE sampling...")
  smoteFit <- train(x = trainX,
                    y = refactoredTrainY,
                    method = caretModel,
                    preProcess = c("scale","center"),
                    trControl = ctrl,
                    verbose = TRUE,
                    metric = "ROC",
                    tuneGrid = caretTuneGrid)
  
  smoteFit %>%
    test_roc(data = testX, classes = refactoredTestY) %>%
    auc() %>% print()
  
  # Build ROSE model
  ctrl$sampling <- "rose"
  
  print("Now training model with ROSE sampling...")
  roseFit <- train(x = trainX,
                   y = refactoredTrainY,
                   method = caretModel,
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   verbose = TRUE,
                   metric = "ROC",
                   tuneGrid = caretTuneGrid)
  
  roseFit %>%
    test_roc(data = testX, classes = refactoredTestY) %>%
    auc() %>% print()
  
  # Examine results for test set
  
  modelVanillaCM <- confusionMatrix(predict(modelVanilla, newdata = testX, type = "raw"), 
                                    refactoredTestY, 
                                    positive = gsub('([[:punct:]])|\\s+','',tumorTissueComparison))
  weightedFitCM <- confusionMatrix(predict(weightedFit, newdata = testX, type = "raw"), 
                                    refactoredTestY, 
                                    positive = gsub('([[:punct:]])|\\s+','',tumorTissueComparison))
  downFitCM <- confusionMatrix(predict(downFit, newdata = testX, type = "raw"), 
                                    refactoredTestY, 
                                    positive = gsub('([[:punct:]])|\\s+','',tumorTissueComparison))
  upFitCM <- confusionMatrix(predict(upFit, newdata = testX, type = "raw"), 
                                    refactoredTestY, 
                                    positive = gsub('([[:punct:]])|\\s+','',tumorTissueComparison))
  smoteFitCM <- confusionMatrix(predict(smoteFit, newdata = testX, type = "raw"), 
                                    refactoredTestY, 
                                    positive = gsub('([[:punct:]])|\\s+','',tumorTissueComparison))
  roseFitCM <- confusionMatrix(predict(roseFit, newdata = testX, type = "raw"), 
                                refactoredTestY, 
                                positive = gsub('([[:punct:]])|\\s+','',tumorTissueComparison))
  
  confusionMatrices <- list(original = modelVanillaCM,
                            weighted = weightedFitCM,
                            down = downFitCM,
                            up = upFitCM,
                            SMOTE = smoteFitCM,
                            ROSE = roseFitCM)
  
  model_list <- list(original = modelVanilla,
                     weighted = weightedFit,
                     down = downFit,
                     up = upFit,
                     SMOTE = smoteFit,
                     ROSE = roseFit)
  
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
  
  rocInsetText <- paste("Max test set AUC:\n",sprintf("%1.4f",maxROC),"\n",paste0('(',bestSamplingStrategy,')'))
  
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
    rocTitle <- paste(paste("ROC curves",
                            paste0("(Learner: ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")\n")),
                      cancerTypeString,"\n",
                      paste0('(',numTumor,' ',tumorTissueComparison,' vs. ',numNl,' ',nlTissueComparison,' from ',cancerTypeStringAlternativeNormal,')') )
    rocPlotFileTitle <- paste(paste("ROC curves",
                            paste0("(Learner ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")")),
                      cancerTypeString,
                      paste0('(',numTumor,' ',tumorTissueComparison,' vs. ',numNl,' ',nlTissueComparison,' from ',cancerTypeStringAlternativeNormal,').png') )
    }
  else{
    rocTitle <- paste(paste("ROC curves",
                            paste0("(Learner: ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")\n")),
                      cancerTypeString,"\n",
                      paste0('(',numTumor,' ',tumorTissueComparison,' vs. ',numNl,' ',nlTissueComparison,')') )
    rocPlotFileTitle <- paste(paste("ROC curves",
                            paste0("(Learner ",caretModel,"|","rep=",as.character(numResampleIter),"|","kFold=",as.character(numKFold),")")),
                      cancerTypeString,
                      paste0('(',numTumor,' ',tumorTissueComparison,' vs. ',numNl,' ',nlTissueComparison,').png') )
  }
  
  g <- ggplot(aes(x = fpr,  y = tpr, group = model), data = results_df_roc) +
    geom_path(aes(color = model), size = 1) +
    scale_color_manual(values = custom_col) +
    geom_abline(intercept = 0, slope = 1, color = "gray", size = 1) +
    theme_bw(base_size = 18) + coord_equal(ratio=1) + xlim(0, 1) + ylim(0,1) +
    labs(color="Sampling strategy\nduring model training", x = "False Positive Rate", y = "True Positive Rate", title = rocTitle) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    annotate("text", x = 0.75, y = 0.25, label = rocInsetText)
  
  ggsave(plot = g, 
         filename = rocPlotFileTitle,
         path = "./roc-ggplots-predicting-tumor-vs-normal",
         device = "png",
         dpi = "retina")
  
  
  resampling <- resamples(model_list)
  trellis.par.set(caretTheme())
  p <- bwplot(resampling)
  
  
  save(model_list, file = paste0(caretModel,"_model_list.RData"))
  
  results <- list(
    # trainX = trainX,
    # trainY = refactoredTrainY,
    testX = testX,
    testY = refactoredTestY,
    resamplingPlot = p,
    resampling = resampling,
    rocPlot = g,
    rocPlotData = results_df_roc,
    rocModelData = model_list_roc,
    modelList = model_list,
    modelWeights = model_weights,
    confusionMatrices = confusionMatrices)
  
  return(results)
}

## To avoid the error below when training the GBM learner, a custom GBM grid can be entered
# The dataset size is too small or subsampling rate is too large: nTrain*bag.fraction <= n.minobsinnode
defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                        n.trees = floor((1:3) * 50),
                        shrinkage = 0.1,
                        n.minobsinnode = 5)

# customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
#                         n.trees = floor((1:3) * 50),
#                         shrinkage = 0.1,
#                         n.minobsinnode = 1)

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
                                               qcMLDataSNM = snmDataSampleType,
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