# All_Tumor_clinical_analyses.R
# Author: Greg Poore
# Date: Oct 8, 2018
# Purpose: To explore clinically-oriented differential expression and machine learning analyses

# Load dependencies
require(ggplot2)
require(ggsci)
require(limma)
require(Glimma)
require(edgeR)
require(maftools)
require(dplyr)
require(TCGAmutations)
require(doMC)

numCores <- detectCores()
registerDoMC(cores=numCores)

#------------------------------------------------------

# Load data
load("tcgaVbDataAndMetadataAndSNM.RData")
load("alphaDiversityMetrics.RData")
load("snmDataSampleTypeWithExpStrategyFINAL.RData")

#------------------------------------------------------

# require(sevenbridges)
# a <- Auth(token = "e5664ac8582f46b5b68c08a88381fbea",
#           platform = "cgc")
# a$api(path = "projects", method = "GET")
# a$project()
# 
# p <- a$project(id = "jkanbar/tcga-kraken")
# p$file()
# tmp <- p$file(fields="_all", complete = TRUE)
# 
# 
# cols2Keep <- c(
#   "reference_genome",
#   "case_id",
#   "experimental_strategy",
#   "disease_type",
#   "aliquot_uuid",
#   "gender",
#   "aliquot_id",
#   "data_subtype",
#   "sample_uuid",
#   "platform",
#   "investigation",
#   "case_uuid",
#   "data_format",
#   "data_type",
#   "sample_type",
#   "primary_site",
#   "sample_id"
# )
# tmp3 <- list()
# for(ii in 1:length(tmp)){
#   tmp3[[ii]] <- data.frame(filename = tmp[[ii]]$name, as.data.frame(tmp[[ii]]$`.->metadata`)[,cols2Keep])
# }
# 
# tmp4 <- do.call("rbind", tmp3)
# cgcMetadataKrakenProj <- tmp4
# save(cgcMetadataKrakenProj, file = "cgcMetadataKrakenProj.RData")

#------------------------------------------------------
load("cgcMetadataKrakenProj.RData")
load("cgcAPIMetadataJoined.RData")
# metadataSamplesAllCGC <- left_join(metadataSamplesAll, 
#                                    cgcMetadataKrakenProj[, -which(names(cgcMetadataKrakenProj) %in% 
#                                                                     c("reference_genome",
#                                                                       "experimental_strategy",
#                                                                       "disease_type",
#                                                                       "aliquot_uuid",
#                                                                       "gender",
#                                                                       "sample_uuid",
#                                                                       "platform",
#                                                                       "investigation",
#                                                                       "case_uuid",
#                                                                       "sample_type",
#                                                                       "primary_site"))],
#                                    by = "filename")
# rownames(metadataSamplesAllCGC) <- rownames(metadataSamplesAll)
# head(metadataSamplesAllCGC,2)
# 
# metadataSamplesAllQCCGC <- left_join(metadataSamplesAllQC, 
#                                    cgcMetadataKrakenProj[, -which(names(cgcMetadataKrakenProj) %in% 
#                                                                     c("reference_genome",
#                                                                       "experimental_strategy",
#                                                                       "disease_type",
#                                                                       "aliquot_uuid",
#                                                                       "gender",
#                                                                       "sample_uuid",
#                                                                       "platform",
#                                                                       "investigation",
#                                                                       "case_uuid",
#                                                                       "sample_type",
#                                                                       "primary_site"))],
#                                    by = "filename")
# rownames(metadataSamplesAllQCCGC) <- rownames(metadataSamplesAllQC)
# head(metadataSamplesAllQCCGC,2)
# 
# metadataSamplesAllQCSurvivalCGC <- left_join(metadataSamplesAllQCSurvival, 
#                                      cgcMetadataKrakenProj[, -which(names(cgcMetadataKrakenProj) %in% 
#                                                                       c("reference_genome",
#                                                                         "experimental_strategy",
#                                                                         "disease_type",
#                                                                         "aliquot_uuid",
#                                                                         "gender",
#                                                                         "sample_uuid",
#                                                                         "platform",
#                                                                         "investigation",
#                                                                         "case_uuid",
#                                                                         "sample_type",
#                                                                         "primary_site"))],
#                                      by = "filename")
# rownames(metadataSamplesAllQCSurvivalCGC) <- rownames(metadataSamplesAllQCSurvival)
# head(metadataSamplesAllQCSurvivalCGC,2)
# 
# save(metadataSamplesAllCGC, metadataSamplesAllQCCGC, metadataSamplesAllQCSurvivalCGC,
#      file = "cgcAPIMetadataJoined.RData")
#------------------------------------------------------

# Load clinical metadata from TCGAmutations
## Load available TCGA datasets and extract radiation information
availableStudies <- TCGAmutations::tcga_available()
for(tcgaStudy in head(availableStudies$Study_Abbreviation,-1)){
  tcga_load(study = tcgaStudy)
}
accClinicalMetadata <- as.data.frame(getClinicalData(tcga_acc_mc3))
blcaClinicalMetadata <- as.data.frame(getClinicalData(tcga_blca_mc3)) #
brcaClinicalMetadata <- as.data.frame(getClinicalData(tcga_brca_mc3)) #
cescClinicalMetadata <- as.data.frame(getClinicalData(tcga_cesc_mc3))
cholClinicalMetadata <- as.data.frame(getClinicalData(tcga_chol_mc3))

coadClinicalMetadata <- as.data.frame(getClinicalData(tcga_coad_mc3)) #
dlbcClinicalMetadata <- as.data.frame(getClinicalData(tcga_dlbc_mc3))
escaClinicalMetadata <- as.data.frame(getClinicalData(tcga_esca_mc3))
gbmClinicalMetadata <- as.data.frame(getClinicalData(tcga_gbm_mc3))

hnscClinicalMetadata <- as.data.frame(getClinicalData(tcga_hnsc_mc3)) #
kichClinicalMetadata <- as.data.frame(getClinicalData(tcga_kich_mc3))
kircClinicalMetadata <- as.data.frame(getClinicalData(tcga_kirc_mc3)) #
kirpClinicalMetadata <- as.data.frame(getClinicalData(tcga_kirp_mc3))
lamlClinicalMetadata <- as.data.frame(getClinicalData(tcga_laml_mc3))

lggClinicalMetadata <- as.data.frame(getClinicalData(tcga_lgg_mc3)) #
lihcClinicalMetadata <- as.data.frame(getClinicalData(tcga_lihc_mc3)) #
luadClinicalMetadata <- as.data.frame(getClinicalData(tcga_luad_mc3)) #
luscClinicalMetadata <- as.data.frame(getClinicalData(tcga_lusc_mc3))

mesoClinicalMetadata <- as.data.frame(getClinicalData(tcga_meso_mc3))
ovClinicalMetadata <- as.data.frame(getClinicalData(tcga_ov_mc3)) #
paadClinicalMetadata <- as.data.frame(getClinicalData(tcga_paad_mc3))
pcpgClinicalMetadata <- as.data.frame(getClinicalData(tcga_pcpg_mc3))
pradClinicalMetadata <- as.data.frame(getClinicalData(tcga_prad_mc3)) #

readClinicalMetadata <- as.data.frame(getClinicalData(tcga_read_mc3))
sarcClinicalMetadata <- as.data.frame(getClinicalData(tcga_sarc_mc3))
skcmClinicalMetadata <- as.data.frame(getClinicalData(tcga_skcm_mc3)) #
stadClinicalMetadata <- as.data.frame(getClinicalData(tcga_stad_mc3)) #
tgctClinicalMetadata <- as.data.frame(getClinicalData(tcga_tgct_mc3))

thcaClinicalMetadata <- as.data.frame(getClinicalData(tcga_thca_mc3)) #
thymClinicalMetadata <- as.data.frame(getClinicalData(tcga_thym_mc3))
ucecClinicalMetadata <- as.data.frame(getClinicalData(tcga_ucec_mc3)) #
ucsClinicalMetadata <- as.data.frame(getClinicalData(tcga_ucs_mc3))
uvmClinicalMetadata <- as.data.frame(getClinicalData(tcga_uvm_mc3))

# cols2Keep <- c(
#   "Tumor_Sample_Barcode",
#   "patient_id",
#   "bcr_patient_uuid", ## NB: This aligns with the case_uuid in the QIIME mapping file
#   "postoperative_rx_tx",
#   "radiation_therapy"
# )
# ind <- 1
# listClinData <- list()
# for(studyClinData in grep("ClinicalMetadata$", ls(), value = TRUE)){
#   dfClinData <- get(studyClinData)
#   dfClinData$study <- factor(studyClinData)
#   dfClinData$bcr_patient_uuid <- toupper(dfClinData$bcr_patient_uuid)
#   if(any(grepl("^radiation_therapy$",colnames(dfClinData), ignore.case = TRUE))){
#     print(studyClinData)
#     listClinData[[ind]] <- dfClinData[,cols2Keep]
#     ind <- ind + 1
#   }
#   else{next}
# }
# 
# mergedClinicalDataDF <- do.call("rbind", listClinData)
# dim(mergedClinicalDataDF)
# sum(duplicated(mergedClinicalDataDF)) # Shows lots of duplicated rows, so make unique
# uniqueMergedClinicalDataDF <- unique(mergedClinicalDataDF) 
# table(uniqueMergedClinicalDataDF$radiation_therapy)
# 
# dim(uniqueMergedClinicalDataDF)

#------------- COAD: MSI vs. MSS -------------#

# Data alignment
coadClinicalCols2Keep <- c(
  "Tumor_Sample_Barcode",
  "bcr_patient_uuid", ## NB: This aligns with the case_uuid in the QIIME mapping file
  "microsatellite_instability",
  "kras_mutation_found",
  "colon_polyps_present",
  "loss_expression_of_mismatch_repair_proteins_by_ihc"
)

coadMetadataCGC <- metadataSamplesAllCGC[metadataSamplesAllCGC$disease_type == "Colon Adenocarcinoma",]
coadMetadataCGC$case_uuid <- toupper(coadMetadataCGC$case_uuid)
coadClinicalMetadata$bcr_patient_uuid <- toupper(coadClinicalMetadata$bcr_patient_uuid)

coadMetadataCGCClinical <- left_join(coadMetadataCGC, 
                                     coadClinicalMetadata[,coadClinicalCols2Keep], 
                                     by = c("case_uuid" = "bcr_patient_uuid"))
rownames(coadMetadataCGCClinical) <- rownames(coadMetadataCGC)

# Subset data
coadMetadataCGCClinical_microsatellite <- coadMetadataCGCClinical[!is.na(coadMetadataCGCClinical$microsatellite_instability),]
coadMetadataCGCClinical_microsatellitePT <- droplevels(coadMetadataCGCClinical_microsatellite[coadMetadataCGCClinical_microsatellite$sample_type == "Primary Tumor",])
#---------------------------
voomMetadata <- coadMetadataCGCClinical_microsatellitePT
voomCountData <- t(vbDataBarnDFReconciled[rownames(voomMetadata),])
# Differential abundance analysis
limmaFormula <- formula(~0 + microsatellite_instability + data_submitting_center_label + platform)

covDesignShort <- model.matrix(limmaFormula, data = voomMetadata)
colnames(covDesignShort) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignShort))
colnames(covDesignShort)

dge <- DGEList(counts = voomCountData)
keep <- filterByExpr(dge, covDesignShort)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "TMM")
vdge <- voom(dge, design = covDesignShort, plot = TRUE, save.plot = TRUE, normalize.method="none")
vdgeFit <- lmFit(vdge, covDesignShort)
vdgeFit <- eBayes(vdgeFit)
contrast.matrix <- makeContrasts(microsatelliteinstabilityYES - microsatelliteinstabilityNO, levels = covDesignShort)
vdgeFit2 <- contrasts.fit(vdgeFit, contrasts = contrast.matrix)
vdgeFit2 <- eBayes(vdgeFit2)
vdgeFitDT <- decideTests(vdgeFit2)

# results <- list(dge = dge, vdge = vdge, vdgeFit2 = vdgeFit2, vdgeFitDT = vdgeFitDT, limmaFormula = limmaFormula)

titleXY <- "Volcano Plot - COAD Primary Tumors - Microsatellite Instability (Y|N = 18|100)"
glXYPlot(x = vdgeFit2$coef, y = vdgeFit2$lod,
         # coef = 1,
         xlab="logFC",
         ylab="log-odds",
         status = vdgeFitDT,
         groups = coadMetadataCGCClinical_microsatellitePT$microsatellite_instability,
         main = titleXY,
         side.ylab = "Voom Normalized Abundance",
         html = titleXY,
         folder = "Glimma-plots-clinical-Volcano",
         counts = vdge)

# predMSSvsMSI.R

coadMetadataQCCGC <- metadataSamplesAllQCCGC[metadataSamplesAllQCCGC$disease_type == "Colon Adenocarcinoma",]
coadMetadataQCCGC$case_uuid <- toupper(coadMetadataQCCGC$case_uuid)
coadClinicalMetadata$bcr_patient_uuid <- toupper(coadClinicalMetadata$bcr_patient_uuid)

coadMetadataQCCGCClinical <- left_join(coadMetadataQCCGC, 
                                     coadClinicalMetadata[,coadClinicalCols2Keep], 
                                     by = c("case_uuid" = "bcr_patient_uuid"))
rownames(coadMetadataQCCGCClinical) <- rownames(coadMetadataQCCGC)

# Subset data
coadMetadataQCCGCClinical_microsatellite <- coadMetadataQCCGCClinical[!is.na(coadMetadataQCCGCClinical$microsatellite_instability),]
coadMetadataQCCGCClinical_microsatellitePT <- droplevels(coadMetadataQCCGCClinical_microsatellite[coadMetadataQCCGCClinical_microsatellite$sample_type == "Primary Tumor",])
#---------------------------

source("predMSSvsMSI.R")
customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                              n.trees = floor((1:3) * 50),
                              shrinkage = 0.1,
                              n.minobsinnode = 1)
tmp <- predMSSvsMSI(qcMLMetadata = coadMetadataQCCGCClinical_microsatellitePT, 
                         qcMLDataSNM = snmDataSampleType[rownames(coadMetadataQCCGCClinical_microsatellitePT),],
                         cancerTypeString = "Colon Adenocarcinoma",
                         sampleTypeComparison = "Primary Tumor",
                         caretModel = "gbm",
                         numResampleIter = 3,
                         numKFold = 4,
                         trainSetProp = 0.5,
                         caretTuneGrid = customGBMGrid,
                         ggPath = "./roc-ggplots-clinical")

#-----------------------------------------#
#------------- Sanity checks -------------#
#-----------------------------------------#
require(dplyr)
require(ggplot2)
require(ggpubr)
require(bigrquery)

#------------- HPV status across all cancers -------------#

clinical_table = "[isb-cgc:tcga_201607_beta.Clinical_data]"
cloud_project_workshop = "hybrid-coyote-219120"
sqlQuery = paste("SELECT ParticipantBarcode, Study, hpv_calls, hpv_status ", 
                 "FROM ", clinical_table,sep="")
sqlQuery
hpv_table = query_exec(sqlQuery,project = cloud_project_workshop)

isbcgcHPV <- hpv_table
hpvPancancerMeta <- left_join(metadataSamplesAllQCCGC, isbcgcHPV, by = c("case_id" = "ParticipantBarcode") )
hpvPancancerData <- data.frame(HPV = snmDataSampleTypeWithExpStrategy[rownames(metadataSamplesAllQCCGC),"k__Viruses.f__Papillomaviridae.g__Alphapapillomavirus"])
hpvPancancerCombined <- droplevels(cbind(hpvPancancerMeta, hpvPancancerData))

interactVec <- as.character(interaction(hpvPancancerCombined$investigation, hpvPancancerCombined$hpv_status, sep = " "))
interactVec[which(is.na(interactVec))] <- as.character(hpvPancancerCombined$investigation[which(is.na(interactVec))])
interactVec <- factor(interactVec)
hpvPancancerCombined$hpvInteract <- interactVec

# hpvPancancerCombined %>%
#   filter((sample_type %in% c("Primary Tumor")) &
#            !(hpv_status %in% c("Indeterminate", "Positive"))) %>%
#   summarise(mean = mean(HPV)) -> ptGrandMeans
# 
# hpvPancancerCombined %>%
#   filter((sample_type %in% c("Solid Tissue Normal")) &
#            !(hpv_status %in% c("Indeterminate", "Positive"))) %>%
#   summarise(mean = mean(HPV)) -> stnGrandMeans
# 
# hpvPancancerCombined %>%
#   filter((sample_type %in% c("Blood Derived Normal")) &
#            !(hpv_status %in% c("Indeterminate", "Positive"))) %>%
#   summarise(mean = mean(HPV)) -> bdnGrandMeans

hpvPancancerCombined %>%
  filter((sample_type %in% c("Primary Tumor")) &
           (hpv_status %in% c("Negative"))) %>%
  summarise(mean = mean(HPV)) -> ptGrandMeans

hpvPancancerCombined %>%
  filter((sample_type %in% c("Solid Tissue Normal")) &
           (hpv_status %in% c("Negative"))) %>%
  summarise(mean = mean(HPV)) -> stnGrandMeans

hpvPancancerCombined %>%
  filter((sample_type %in% c("Blood Derived Normal")) &
           (hpv_status %in% c("Negative"))) %>%
  summarise(mean = mean(HPV)) -> bdnGrandMeans

allGrandMeans <- data.frame(sample_type = c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal"),
                            means = rbind(ptGrandMeans, stnGrandMeans, bdnGrandMeans))
hpvPancancerCombined %>%
  filter((sample_type %in% c("Blood Derived Normal", "Solid Tissue Normal", "Primary Tumor")) &
           !(hpv_status %in% c("Indeterminate"))) %>%
  ggboxplot(x = "hpvInteract", y = "HPV", 
            color = "sample_type",
            # add = "median_iqr",
            # facet.by = "sample_type",
            # palette = "Blues",
            xlab = "Disease Type", ylab = "SNM Normalized Abundance", 
            title = "Pancancer Comparison of Alphapapillomavirus Genus Abundance",
            # legend = "right",
            # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  geom_hline(data = allGrandMeans, aes(yintercept = mean),
             linetype = 2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(angle = 30) -> p

facet(p, facet.by = "sample_type", nrow = 3, ncol = 1)

hpvCervicalCancerComparisons <- list( c("TCGA-CESC Positive", "TCGA-CESC Negative"))
hpvPancancerCombined %>%
  filter((sample_type %in% c("Blood Derived Normal", "Primary Tumor")) &
           !(hpv_status %in% c("Indeterminate")) &
           (disease_type == "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma")) %>%
  ggboxplot(x = "hpvInteract", y = "HPV", 
            # color = "sample_type",
            add = "jitter",
            facet.by = "sample_type",
            # palette = pal_nejm(),
            xlab = "Clinical HPV Status", ylab = "SNM Normalized Abundance",
            ylim = c(-11, 20),
            title = "Pancancer Comparison of Alphapapillomavirus Genus Abundance in Cervical Cancer",
            # legend = "right",
            # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(labels=c("TCGA-CESC Negative" = "Negative",
                            "TCGA-CESC Positive" = "Positive")) +
  scale_color_nejm() +
  rotate_x_text(angle = 30) +
  stat_compare_means(comparisons = hpvCervicalCancerComparisons, label = "p.signif", method.args = list(alternative = "greater")) -> p# Add pairwise comparisons p-value
  
ggsave(p, filename = "HPV in CESC.png", 
       path = "./Clinical Validation Plots",
       dpi = "retina",
       units = "in",
       height = 5, width = 4)

#------------- LCV status in stomach cancers -------------#

stadMasterPatientTable <- read.csv(file = "STAD_Master_Patient_Table_20140207.csv", stringsAsFactors = FALSE)

stadPancancerMeta <- left_join(metadataSamplesAllQCCGC,
                               stadMasterPatientTable, by = c("case_id" = "TCGA.barcode") )
rownames(stadPancancerMeta) <- rownames(metadataSamplesAllQCCGC)
lcvPancancerData <- data.frame(LCV = snmDataSampleTypeWithExpStrategy[rownames(stadPancancerMeta),"k__Viruses.o__Herpesvirales.f__Herpesviridae.g__Lymphocryptovirus"],
                               HPylori = snmDataSampleTypeWithExpStrategy[rownames(stadPancancerMeta),"k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales.f__Helicobacteraceae.g__Helicobacter"])
lcvPancancerCombined <- droplevels(cbind(stadPancancerMeta, lcvPancancerData))

lcvComparisons <- list( c("EBV", "CIN"), c("EBV", "GS"), c("EBV","MSI") )
lcvPancancerCombined %>%
  filter((sample_type %in% c("Blood Derived Normal", "Solid Tissue Normal", "Primary Tumor")) &
           !(is.na(Molecular.Subtype)) &
           (disease_type == "Stomach Adenocarcinoma")) %>%
  ggboxplot(x = "Molecular.Subtype", y = "LCV", 
            # color = "sample_type",
            add = "jitter",
            facet.by = "sample_type",
            # palette = "lancet",
            ylim = c(-5, 22),
            xlab = "STAD Molecular subtype (The Cancer Genome Atlas Research Network, 2014. Nature)", ylab = "SNM Normalized Abundance", 
            title = "Pancancer Comparison of Lymphocryptovirus Genus Abundance in Stomach Adenocarcinoma",
            # legend = "right",
            # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  scale_color_nejm() +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = lcvComparisons, label = "p.signif", method = "wilcox.test") -> p # Add pairwise comparisons p-value

ggsave(p, filename = "EBV in STAD.png", path = "./Clinical Validation Plots",dpi = "retina", units = "in",
       height = 4, width = 7)

hpyloriComparisons <- list( c("Primary Tumor", "Solid Tissue Normal"))
lcvPancancerCombined %>%
  filter((sample_type %in% c("Solid Tissue Normal", "Primary Tumor")) &
           (disease_type == "Stomach Adenocarcinoma")) %>%
  # filter(experimental_strategy == "WGS") %>%
  ggboxplot(x = "sample_type", y = "HPylori", 
            # color = "sample_type",
            add = "jitter",
            line.color = "gray",
            # facet.by = "sample_type",
            palette = "lancet",
            # ylim = c(-5, 25),
            xlab = "Sample Type", ylab = "SNM Normalized Abundance", 
            title = "Comparison of Helicobacter Genus Abundance in Stomach Adenocarcinoma",
            # legend = "right",
            # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(angle = 15) + 
  stat_compare_means(comparisons = hpyloriComparisons,
                     method = "wilcox.test",
                     # ref.group = "Solid Tissue Normal",
                     # comparisons = hpyloriComparisons,
                     method.args = list(alternative = "less")) -> p# Add pairwise comparisons p-value

ggsave(p, filename = "H Pylori in STAD.png", path = "./Clinical Validation Plots",dpi = "retina", units = "in",
       height = 4, width = 4)

# hpyloriComparisons <- list( c("Primary Tumor", "Solid Tissue Normal"))
# lcvPancancerCombined %>%
#   filter((sample_type %in% c("Solid Tissue Normal", "Primary Tumor")) &
#            (disease_type == "Stomach Adenocarcinoma")) %>%
#   group_by(sample_id) %>%
#   filter(n() >= 2) %>% # & experimental_strategy == "RNA-Seq"
#   ggpaired(x = "sample_type", y = "HPylori", 
#             color = "sample_type",
#             id = "case_uuid",
#             add = "jitter",
#             # facet.by = "sample_type",
#             palette = "lancet",
#             # ylim = c(-5, 25),
#             xlab = "Sample Type", ylab = "SNM Normalized Abundance", 
#             title = "Comparison of Helicobacter Genus Abundance in Stomach Adenocarcinoma",
#             # legend = "right",
#             # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
#             legend.title = "Sample Type") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   stat_compare_means(comparisons = hpyloriComparisons, method = "t.test") # Add pairwise comparisons p-value


#------------- HPV status in cervical cancer -------------#

cervicalHPVMeta <- droplevels(metadataSamplesAllQC[metadataSamplesAllQC$disease_type == "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",])
cervicalHPVData <- data.frame(HPV = snmDataSampleTypeWithExpStrategy[rownames(cervicalHPVMeta),"k__Viruses.f__Papillomaviridae.g__Alphapapillomavirus"])
cervicalHPVCombined <- cbind(cervicalHPVMeta, cervicalHPVData)

cervicalHPVcomparisons <- list( c("Primary Tumor", "Solid Tissue Normal"), c("Primary Tumor", "Blood Derived Normal"), c("Solid Tissue Normal", "Blood Derived Normal"))
cervicalHPVCombined %>%
  filter(sample_type %in% c("Blood Derived Normal", "Solid Tissue Normal", "Primary Tumor")) %>%
  ggboxplot(x = "sample_type", y = "HPV", 
            color = "sample_type", 
            add = "jitter", 
            palette = "lancet",
            xlab = "Sample Type", ylab = "SNM Normalized Abundance", 
            title = "Comparison of Alphapapillomavirus Genus Abundance in Cervical Cancer",
            legend = "right",
            order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = cervicalHPVcomparisons, label = "p.signif") + # Add pairwise comparisons p-value
  stat_compare_means(label.y = -10)                   # Add global p-value
  
#------------- HPV status in HNSC -------------#

# Data alignment
hnscClinicalCols2Keep <- c(
  "Tumor_Sample_Barcode",
  "bcr_patient_uuid", ## NB: This aligns with the case_uuid in the QIIME mapping file
  "hpv_status_by_p16_testing",
  "hpv_status_by_ish_testing"
)

hnscMetadataQCCGC <- metadataSamplesAllQCCGC[metadataSamplesAllQCCGC$disease_type == "Head and Neck Squamous Cell Carcinoma",]
hnscMetadataQCCGC$case_uuid <- toupper(hnscMetadataQCCGC$case_uuid)
hnscClinicalMetadata$bcr_patient_uuid <- toupper(hnscClinicalMetadata$bcr_patient_uuid)

hnscMetadataQCCGCClinical <- left_join(hnscMetadataQCCGC, 
                                       hnscClinicalMetadata[,hnscClinicalCols2Keep], 
                                       by = c("case_uuid" = "bcr_patient_uuid"))
rownames(hnscMetadataQCCGCClinical) <- rownames(hnscMetadataQCCGC)

# Subset data
hnscMetadataQCCGCClinical_HPVp16 <- droplevels(hnscMetadataQCCGCClinical[!is.na(hnscMetadataQCCGCClinical$hpv_status_by_p16_testing),])
hnscMetadataQCCGCClinical_HPVish <- droplevels(hnscMetadataQCCGCClinical[!is.na(hnscMetadataQCCGCClinical$hpv_status_by_ish_testing),])

hnscHPVp16Data <- data.frame(HPV = snmDataSampleTypeWithExpStrategy[rownames(hnscMetadataQCCGCClinical_HPVp16),"k__Viruses.f__Papillomaviridae.g__Alphapapillomavirus"])
hnscHPVishData <- data.frame(HPV = snmDataSampleTypeWithExpStrategy[rownames(hnscMetadataQCCGCClinical_HPVish),"k__Viruses.f__Papillomaviridae.g__Alphapapillomavirus"])
hnscHPVp16Combined <- cbind(hnscMetadataQCCGCClinical_HPVp16, hnscHPVp16Data)
hnscHPVishCombined <- cbind(hnscMetadataQCCGCClinical_HPVish, hnscHPVishData)

testType <- factor(c(rep("p16 Testing",dim(hnscHPVp16Combined)[1]), rep("ISH Testing", dim(hnscHPVishCombined)[1])))
testValue <- factor(c(as.character(hnscHPVp16Combined$hpv_status_by_p16_testing), as.character(hnscHPVishCombined$hpv_status_by_ish_testing)))
hnscHPVbothCombined <- cbind(rbind(hnscHPVp16Combined,hnscHPVishCombined),testType, testValue)
hnscHPVcomparisons <- list( c("Negative", "Positive"))
hnscHPVbothCombined %>%
  filter(sample_type == "Primary Tumor") %>%
  filter(!is.na(hpv_status_by_p16_testing)) %>%
  filter(!is.na(hpv_status_by_ish_testing)) %>%
  ggboxplot(x = "testValue", y = "HPV", 
            # color = "testValue",
            facet.by = "testType",
            add = "jitter",
            # palette = "lancet",
            xlab = "Clinical Testing for HPV", ylab = "SNM Normalized Abundance", 
            ylim = c(-3, 18),
            title = "Comparison of Alphapapillomavirus Genus Abundance in\nHead and Neck Squamous Cell Carcinoma Primary Tumors",
            legend = "right",
            legend.title = "Clinical HPV Status",
            font.label = list(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_nejm() +
  rotate_x_text(angle = 30) +
  stat_compare_means(comparisons = hnscHPVcomparisons, 
                     label = "p.signif",
                     method.args = list(alternative = "less"),
                     method = "t.test") -> p # Add pairwise comparisons p-value

ggsave(p, filename = "HPV in HNSCC.png", 
       path = "./Clinical Validation Plots",
       dpi = "retina",
       units = "in",
       height = 4, width = 3)

# my_comparisons <- list( c("Primary Tumor", "Solid Tissue Normal"), c("Primary Tumor", "Blood Derived Normal"), c("Solid Tissue Normal", "Blood Derived Normal"))

hnscHPVcomparisons <- list( c("Negative", "Positive"))
hnscHPVp16Combined %>%
  filter(sample_type == "Primary Tumor") %>%
ggboxplot(x = "hpv_status_by_p16_testing", y = "HPV", 
          color = "hpv_status_by_p16_testing",
          add = "jitter",
          palette = "lancet",
          xlab = "Clinical Testing for HPV", ylab = "SNM Normalized Abundance", 
          title = "Comparison of Alphapapillomavirus Genus Abundance in\nHead and Neck Squamous Cell Carcinoma Primary Tumors",
          legend = "right",
          legend.title = "HPV Status by\nP16 Testing",
          font.label = list(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = hnscHPVcomparisons, 
                     method.args = list(alternative = "less"),
                     method = "t.test") # Add pairwise comparisons p-value

hnscHPVcomparisons <- list( c("Negative", "Positive"))
hnscHPVishCombined %>%
  filter(sample_type == "Primary Tumor") %>%
ggboxplot(x = "hpv_status_by_ish_testing", y = "HPV", 
          color = "hpv_status_by_ish_testing",
          add = "jitter",
          palette = "lancet",
          xlab = "Clinical Testing for HPV", ylab = "SNM Normalized Abundance", 
          title = "Comparison of Alphapapillomavirus Genus Abundance in\nHead and Neck Squamous Cell Carcinoma Primary Tumors",
          legend = "right",
          legend.title = "HPV Status by\nISH Testing",
          font.label = list(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = hnscHPVcomparisons, 
                     method.args = list(alternative = "less"),
                     method = "t.test") # Add pairwise comparisons p-value

#------------- HBV/HCV status in LIHC -------------#

# Data alignment
lihcClinicalCols2Keep <- c(
  "Tumor_Sample_Barcode",
  "bcr_patient_uuid", ## NB: This aligns with the case_uuid in the QIIME mapping file
  "history_hepato_carcinoma_risk_factors"
)

lihcMetadataQCCGC <- metadataSamplesAllQCCGC[metadataSamplesAllQCCGC$disease_type == "Liver Hepatocellular Carcinoma",]
lihcMetadataQCCGC$case_uuid <- toupper(lihcMetadataQCCGC$case_uuid)
lihcClinicalMetadata$bcr_patient_uuid <- toupper(lihcClinicalMetadata$bcr_patient_uuid)

lihcMetadataQCCGCClinical <- left_join(lihcMetadataQCCGC, 
                                       lihcClinicalMetadata[,lihcClinicalCols2Keep], 
                                       by = c("case_uuid" = "bcr_patient_uuid"))
rownames(lihcMetadataQCCGCClinical) <- rownames(lihcMetadataQCCGC)

# Subset data
lihcMetadataQCCGCClinical_Riskfactors <- droplevels(lihcMetadataQCCGCClinical[!is.na(lihcMetadataQCCGCClinical$history_hepato_carcinoma_risk_factors),])

hbvHcvGenera <- c("k__Viruses.f__Flaviviridae.g__Hepacivirus", "k__Viruses.f__Hepadnaviridae.g__Orthohepadnavirus")
lihcHepData <- snmDataSampleTypeWithExpStrategy[rownames(lihcMetadataQCCGCClinical_Riskfactors),hbvHcvGenera]
lihcHepDataCombined <- droplevels(cbind(lihcMetadataQCCGCClinical_Riskfactors, lihcHepData))

# my_comparisons <- list( c("Primary Tumor", "Solid Tissue Normal"), c("Primary Tumor", "Blood Derived Normal"), c("Solid Tissue Normal", "Blood Derived Normal"))
lihcHepComparisons <- list( c("Hepatitis_B", "Hepatitis_C"), c("Hepatitis_B","Alcohol_consumption"), c("Hepatitis_C","Alcohol_consumption"))
lihcHepDataCombined %>%
  filter((history_hepato_carcinoma_risk_factors %in% c("Hepatitis_C", "Hepatitis_B","Alcohol_consumption")) &
           (sample_type %in% c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal"))) %>%
  ggboxplot(x = "history_hepato_carcinoma_risk_factors", y = "k__Viruses.f__Hepadnaviridae.g__Orthohepadnavirus", 
          # color = "sample_type",
          facet.by = "sample_type",
          palette = "lancet",
          add = c("jitter"), # add = "jitter", 
          # shape = "sample_type",
          xlab = "Clinically Assessed Patient History Risk Factors for Hepatocellular Carcinoma", ylab = "SNM Normalized Abundance", 
          title = "Comparison of Orthohepadnavirus Genus Abundance in Liver Hepatocellular Carcinoma",
          legend = "right",
          legend.title = "Sample Type",
          font.label = list(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("Alcohol_consumption" = "EtOH", 
                            "Hepatitis_B" = "Hep B",
                            "Hepatitis_C" = "Hep C")) +
  # rotate_x_text(angle = 30) +
  stat_compare_means(comparisons = lihcHepComparisons, label = "p.signif") -> p # Add pairwise comparisons p-value

ggsave(p, filename = "HBV in LIHC.png", path = "./Clinical Validation Plots",dpi = "retina", units = "in",
       height = 4, width = 7)

# # my_comparisons <- list( c("Primary Tumor", "Solid Tissue Normal"), c("Primary Tumor", "Blood Derived Normal"), c("Solid Tissue Normal", "Blood Derived Normal"))
# lihcHepComparisons <- list( c("Hepatitis_B", "Hepatitis_C"), c("Hepatitis_B","Alcohol_consumption"), c("Hepatitis_C","Alcohol_consumption"))
# lihcHepDataCombined %>%
#   filter((history_hepato_carcinoma_risk_factors %in% c("Hepatitis_C", "Hepatitis_B","Alcohol_consumption")) &
#            (sample_type %in% c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal"))) %>%
#   filter(experimental_strategy == "WGS") %>%
#   ggboxplot(x = "history_hepato_carcinoma_risk_factors", y = "k__Viruses.f__Flaviviridae.g__Hepacivirus", 
#             color = "sample_type",
#             facet.by = "sample_type",
#             palette = "lancet",
#             add = c("jitter"), # add = "jitter", 
#             # shape = "sample_type",
#             xlab = "Clinically Assessed Patient History Risk Factors for Hepatocellular Carcinoma", ylab = "SNM Normalized Abundance", 
#             title = "Comparison of Hepacivirus Genus Abundance in Liver Hepatocellular Carcinoma",
#             legend = "right",
#             legend.title = "Sample Type",
#             font.label = list(size = 14, face = "bold")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_x_discrete(labels=c("Alcohol_consumption" = "EtOH", 
#                             "Hepatitis_B" = "Hep B",
#                             "Hepatitis_C" = "Hep C")) +
#   # rotate_x_text(angle = 30) +
#   stat_compare_means(comparisons = lihcHepComparisons, label = "p.signif") # Add pairwise comparisons p-value

#------------- H pylori in STAD -------------#

# Data alignment
stadClinicalCols2Keep <- c(
  "Tumor_Sample_Barcode",
  "bcr_patient_uuid", ## NB: This aligns with the case_uuid in the QIIME mapping file
  "h_pylori_infection"
)

stadMetadataQCCGC <- metadataSamplesAllQCCGC[metadataSamplesAllQCCGC$disease_type == "Stomach Adenocarcinoma",]
stadMetadataQCCGC$case_uuid <- toupper(stadMetadataQCCGC$case_uuid)
stadClinicalMetadata$bcr_patient_uuid <- toupper(stadClinicalMetadata$bcr_patient_uuid)

stadMetadataQCCGCClinical <- left_join(stadMetadataQCCGC, 
                                       stadClinicalMetadata[,stadClinicalCols2Keep], 
                                       by = c("case_uuid" = "bcr_patient_uuid"))
rownames(stadMetadataQCCGCClinical) <- rownames(stadMetadataQCCGC)

# Subset data
stadMetadataQCCGCClinical_HPylori <- droplevels(stadMetadataQCCGCClinical[!is.na(stadMetadataQCCGCClinical$h_pylori_infection),])

hpyloriGenus <- "k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales.f__Helicobacteraceae.g__Helicobacter"
stadHPyloriData <- data.frame(HPylori = snmDataSampleTypeWithExpStrategy[rownames(stadMetadataQCCGCClinical_HPylori),
                                     "k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales.f__Helicobacteraceae.g__Helicobacter"])
stadHPyloriDataCombined <- droplevels(cbind(stadMetadataQCCGCClinical_HPylori, stadHPyloriData))

# my_comparisons <- list( c("Primary Tumor", "Solid Tissue Normal"), c("Primary Tumor", "Blood Derived Normal"), c("Solid Tissue Normal", "Blood Derived Normal"))
# lihcHepComparisons <- list( c("Hepatitis_B", "Hepatitis_C"), c("Hepatitis_B","Alcohol_consumption"), c("Hepatitis_C","Alcohol_consumption"))

stadComparisons <- list( c( "Yes","No"))
stadHPyloriDataCombined %>%
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal")) %>%
  ggboxplot(x = "h_pylori_infection", y = "HPylori", 
            color = "sample_type",
            facet.by = "sample_type",
            palette = "lancet",
            add = c("jitter"), # add = "jitter", 
            # shape = "sample_type",
            xlab = "Clinical H Pylori Testing Result", ylab = "SNM Normalized Abundance", 
            title = "Comparison of Helicobacter Genus Abundance in Stomach Adenocarcinoma",
            legend = "right",
            legend.title = "Sample Type",
            font.label = list(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # scale_x_discrete(labels=c("Alcohol_consumption" = "EtOH", 
  #                           "Hepatitis_B" = "Hep B",
  #                           "Hepatitis_C" = "Hep C")) +
  # rotate_x_text(angle = 30) +
  stat_compare_means(comparisons = stadComparisons, label = "p.signif") # Add pairwise comparisons p-value

#------------- Gender differences -------------#

genderBug <- "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Micrococcales.f__Intrasporangiaceae.g__Tetrasphaera"
genderBugData <- data.frame(Tetrasphaera = snmDataSampleTypeWithExpStrategy[,genderBug])
genderDataCombined <- droplevels(cbind(metadataSamplesAllQC, genderBugData))

genderComparisons <- list( c( "MALE","FEMALE"))
genderCancers <- c("Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",
                   "Breast Invasive Carcinoma",
                   "Ovarian Serous Cystadenocarcinoma",
                   "Prostate Adenocarcinoma",
                   "Testicular Germ Cell Tumors",
                   "Uterine Carcinosarcoma",
                   "Uterine Corpus Endometrial Carcinoma")
genderDataCombined %>%
  filter((sample_type %in% c("Solid Tissue Normal")) &
           !(disease_type %in% genderCancers)) %>%
  ggboxplot(x = "gender", y = "Tetrasphaera", 
            color = "gender",
            facet.by = "disease_type",
            palette = "lancet",
            add = c("jitter"), # add = "jitter", 
            # shape = "sample_type",
            xlab = "Gender", ylab = "SNM Normalized Abundance", 
            ylim = c(10,18),
            title = "Comparison of Tetrasphaera Genus Abundance Across Genders in Solid Tissue Normals",
            # legend = "right",
            legend.title = "Gender") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("FEMALE" = "Female",
                            "MALE" = "Male")) +
  # rotate_x_text(angle = 30) +
  stat_compare_means(comparisons = genderComparisons, label.y = 17, label = "p.signif") # Add pairwise comparisons p-value


#------------- COAD: KRAS Mutation -------------#

# Subset data
coadMetadataCGCClinical_kras <- coadMetadataCGCClinical[!is.na(coadMetadataCGCClinical$kras_mutation_found),]
coadMetadataCGCClinical_krasPT <- droplevels(coadMetadataCGCClinical_kras[coadMetadataCGCClinical_kras$sample_type == "Primary Tumor",])

voomMetadata <- coadMetadataCGCClinical_krasPT
voomCountData <- counts <- t(vbDataBarnDFReconciled[rownames(coadMetadataCGCClinical_krasPT),])
# Differential abundance analysis
limmaFormula <- formula(~0 + kras_mutation_found + data_submitting_center_label + platform)

covDesignShort <- model.matrix(limmaFormula, data = voomMetadata)
colnames(covDesignShort) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignShort))
colnames(covDesignShort)

dge <- DGEList(counts = voomCountData)
keep <- filterByExpr(dge, covDesignShort)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "TMM")
vdge <- voom(dge, design = covDesignShort, plot = TRUE, save.plot = TRUE, normalize.method="none")
vdgeFit <- lmFit(vdge, covDesignShort)
vdgeFit <- eBayes(vdgeFit)
contrast.matrix <- makeContrasts(krasmutationfoundYES - krasmutationfoundNO, levels = covDesignShort)
vdgeFit2 <- contrasts.fit(vdgeFit, contrasts = contrast.matrix)
vdgeFit2 <- eBayes(vdgeFit2)
vdgeFitDT <- decideTests(vdgeFit2)

# results <- list(dge = dge, vdge = vdge, vdgeFit2 = vdgeFit2, vdgeFitDT = vdgeFitDT, limmaFormula = limmaFormula)

titleXY <- "Volcano Plot - COAD Primary Tumors - KRAS Mutation (Y|N = 28|30)"
glXYPlot(x = vdgeFit2$coef, y = vdgeFit2$lod,
         # coef = 1,
         xlab="logFC",
         ylab="log-odds",
         status = vdgeFitDT,
         groups = voomMetadata$kras_mutation_found,
         main = titleXY,
         side.ylab = "Voom Normalized Abundance",
         html = titleXY,
         folder = "Glimma-plots-clinical-Volcano",
         counts = vdge)

#------------- LUAD: Smoking vs nonsmoking -------------#

# Checking which columns to keep
table(luadClinicalMetadata$kras_mutation_found) # Y|N 23|39
table(luadClinicalMetadata$tobacco_smoking_history)

# Data alignment
luadClinicalCols2Keep <- c(
  "Tumor_Sample_Barcode",
  "bcr_patient_uuid", ## NB: This aligns with the case_uuid in the QIIME mapping file
  "tobacco_smoking_history", # See meaning of values here: https://groups.google.com/forum/#!topic/cbioportal/irEXZRj9Who
  "number_pack_years_smoked",
  "kras_mutation_found"
)

luadMetadataCGC <- metadataSamplesAllCGC[metadataSamplesAllCGC$disease_type == "Lung Adenocarcinoma",]
luadMetadataCGC$case_uuid <- toupper(luadMetadataCGC$case_uuid)
luadClinicalMetadata$bcr_patient_uuid <- toupper(luadClinicalMetadata$bcr_patient_uuid)

luadMetadataCGCClinical <- left_join(luadMetadataCGC, 
                                     luadClinicalMetadata[,luadClinicalCols2Keep], 
                                     by = c("case_uuid" = "bcr_patient_uuid"))
rownames(luadMetadataCGCClinical) <- rownames(luadMetadataCGC)

# Subset data
luadMetadataCGCClinical_smoking <- luadMetadataCGCClinical[!is.na(luadMetadataCGCClinical$tobacco_smoking_history),]
luadMetadataCGCClinical_kras <- luadMetadataCGCClinical[!is.na(luadMetadataCGCClinical$kras_mutation_found),]
table(luadMetadataCGCClinical_smoking$sample_type)
# luadMetadataCGCClinical_smoking <- luadMetadataCGCClinical_smoking[luadMetadataCGCClinical_smoking$tobacco_smoking_history %in% c("1","3"),]

luadMetadataCGCClinical_smoking$smokingHistory <- ordered(ifelse(luadMetadataCGCClinical_smoking$tobacco_smoking_history == 1, yes = "Nonsmoker", no = "Smoker"),
                                                          levels = c("Nonsmoker","Smoker"))
luadMetadataCGCClinical_smokingPT <- droplevels(luadMetadataCGCClinical_smoking[luadMetadataCGCClinical_smoking$sample_type == "Primary Tumor",])
luadMetadataCGCClinical_smokingSTN <- droplevels(luadMetadataCGCClinical_smoking[luadMetadataCGCClinical_smoking$sample_type == "Solid Tissue Normal",])
luadMetadataCGCClinical_smokingBDN <- droplevels(luadMetadataCGCClinical_smoking[luadMetadataCGCClinical_smoking$sample_type == "Blood Derived Normal",])

luadMetadataCGCClinical_krasPT <- droplevels(luadMetadataCGCClinical_kras[luadMetadataCGCClinical_kras$sample_type == "Primary Tumor",])
luadMetadataCGCClinical_krasSTN <- droplevels(luadMetadataCGCClinical_kras[luadMetadataCGCClinical_kras$sample_type == "Solid Tissue Normal",])
luadMetadataCGCClinical_krasBDN <- droplevels(luadMetadataCGCClinical_kras[luadMetadataCGCClinical_kras$sample_type == "Blood Derived Normal",])

voomMetadata <- luadMetadataCGCClinical_krasPT
voomCountData <- t(vbDataBarnDFReconciled[rownames(voomMetadata),])
table(voomMetadata$smokingHistory)
table(voomMetadata$kras_mutation_found)
# Differential abundance analysis
limmaFormula <- formula(~0 + kras_mutation_found + data_submitting_center_label)

covDesignShort <- model.matrix(limmaFormula, data = voomMetadata)
colnames(covDesignShort) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignShort))
colnames(covDesignShort)

dge <- DGEList(counts = voomCountData)
keep <- filterByExpr(dge, covDesignShort)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "TMM")
vdge <- voom(dge, design = covDesignShort, plot = TRUE, save.plot = TRUE, normalize.method="none")
vdgeFit <- lmFit(vdge, covDesignShort)
vdgeFit <- eBayes(vdgeFit)
contrast.matrix <- makeContrasts(krasmutationfoundYES - krasmutationfoundNO, levels = covDesignShort)
vdgeFit2 <- contrasts.fit(vdgeFit, contrasts = contrast.matrix)
vdgeFit2 <- eBayes(vdgeFit2)
vdgeFitDT <- decideTests(vdgeFit2)

# results <- list(dge = dge, vdge = vdge, vdgeFit2 = vdgeFit2, vdgeFitDT = vdgeFitDT, limmaFormula = limmaFormula)

titleXY <- "Volcano Plot - LUAD Primary Tumor - KRAS Mutation (Yes|No = 26|55)"
glXYPlot(x = vdgeFit2$coef, y = vdgeFit2$lod,
         # coef = 1,
         xlab="logFC",
         ylab="log-odds",
         status = vdgeFitDT,
         groups = voomMetadata$kras_mutation_found,
         main = titleXY,
         side.ylab = "Voom Normalized Abundance",
         html = titleXY,
         folder = "Glimma-plots-clinical-Volcano",
         counts = vdge)
