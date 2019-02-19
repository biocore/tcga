# Tumor_Microbe_Subtyping.R
# Author: Greg Poore
# Date: Sept 27, 2018
# Purpose: To identify subtype of cancer microbiome (pancancer and per cancer)

#-------------------------------#
# Load dependencies
require(doMC)
require(plyr)
require(dplyr)
require(tibble)
require(CancerSubtypes)
require(SNFtool)
# require(iClusterPlus)
numCores <- detectCores()
registerDoMC(cores=numCores)

## Import cancer microbiome data
load("tcgaMetadataRad.RData") # This has the tumor barcode of interest necessary to merge with RM's data
load("tcgaVbDataAndMetadataAndSNM.RData")
load("immunoOncData.RData")
load("snmDataSampleTypeWithExpStrategyFINAL.RData")
load("cgcMetadataKrakenProj.RData")
load("cgcAPIMetadataJoined.RData")

immunityPaperCibersort <- read.table("TCGA.Kallisto.fullIDs.cibersort.relative.tsv",
                                     sep = "\t",
                                     header = TRUE,
                                     strip.white = TRUE,
                                     stringsAsFactors = FALSE)
immunityPaperCibersortFilt <- immunityPaperCibersort[,-c(2,25:27)]

####

aliquotID <- metadataSamplesAllQCSurvivalCGC$aliquot_id
aliquotIDSub <- gsub("-",".",aliquotID)
metadataSamplesAllQCSurvivalCGC$aliquot_id <- aliquotIDSub
aliquotIDBoolean <- aliquotIDSub %in% immunityPaperCibersort$SampleID

metadataSamplesAllQCSurvivalCGCCibersort <- metadataSamplesAllQCSurvivalCGC[aliquotIDBoolean,]

# #--------------------------- For Gibs ---------------------------#
# 
# aliquotID <- metadataSamplesAllQCCGC$aliquot_id
# aliquotIDSub <- gsub("-",".",aliquotID)
# metadataSamplesAllQCCGC$aliquot_id <- aliquotIDSub
# aliquotIDBoolean <- aliquotIDSub %in% immunityPaperCibersort$SampleID
# 
# metadataSamplesAllQCCGCCibersort <- metadataSamplesAllQCSurvivalCGC[aliquotIDBoolean,]
# 
# aliquotID <- metadataSamplesAllCGC$aliquot_id
# aliquotIDSub <- gsub("-",".",aliquotID)
# metadataSamplesAllCGC$aliquot_id <- aliquotIDSub
# aliquotIDBoolean <- aliquotIDSub %in% immunityPaperCibersort$SampleID
# 
# metadataSamplesAllCGCCibersort <- metadataSamplesAllQCSurvivalCGC[aliquotIDBoolean,]

#------------------------------------------------------#

# coadMetadata <- metadataSamplesAllQC[metadataSamplesAllQC$disease_type == "Colon Adenocarcinoma",]
dzMetadata <- droplevels(metadataSamplesAllQCSurvivalCGCCibersort[(metadataSamplesAllQCSurvivalCGCCibersort$disease_type == "Stomach Adenocarcinoma") &
                                                           (metadataSamplesAllQCSurvivalCGCCibersort$sample_type == "Primary Tumor"),])
dzDataPT <- snmDataSampleTypeWithExpStrategy[rownames(dzMetadata),]

tmp <- data.frame(sampleID = rownames(dzMetadata), aliquotID = dzMetadata$aliquot_id)

tmp2 <- left_join(tmp, immunityPaperCibersortFilt, by = c("aliquotID" = "SampleID"))
tmp3 <- tmp2[!duplicated(tmp2[,c("sampleID")]),]
rownames(tmp3) <- tmp3$sampleID
cibersortDzSampleID <- tmp3[,-c(1:2)]

identical(rownames(dzDataPT),rownames(cibersortDzSampleID)) # Sanity check
dzDataIntersect <- dzDataPT
dzMetadataIntersect <- dzMetadata
dzCibersortIntersect <- cibersortDzSampleID

#-------------------------- SNF+CC --------------------------#
# STAD and RCC stand out
result <- ExecuteSNF.CC(datasets = list(t(dzDataIntersect), 
                                        t(dzCibersortIntersect)#,
                                        # t(coadMutationsIntersect)
                                        ), 
                        clusterNum = 3, K = 25, alpha = 0.5, t = 20,
                        maxK = 10, pItem = 0.8, reps = 1000)
sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
plot(sil)
group=result$group
distanceMatrix=result$distanceMatrix
p_value=survAnalysis(mainTitle="Testing",time = dzMetadataIntersect$days_to_death,
                     status = rep(1, times = length(dzMetadataIntersect$days_to_death)),
                     group,
                     distanceMatrix,similarity=TRUE)

save(result, sil, group, distanceMatrix, p_value, file = "Tumor_Microbe_SubtypingFA_Results_012419.RData")
load("Tumor_Microbe_SubtypingFA_Results_012419.RData")

# tmpgroup <- group

# Reformat figures for NEJM
require(ggsci)
require(pheatmap)
require(factoextra)
require(survminer)
require(survival)

sData <- data.frame(days_to_death = dzMetadataIntersect$days_to_death,
                    status = rep(1, times = length(dzMetadataIntersect$days_to_death)),
                    group = group)

sFit <- survfit(Surv(days_to_death, status) ~ group,
                data = sData )
ggsurvplot(sFit, pval = TRUE, 
           # conf.int = TRUE, 
           palette = c("#0072B5FF","#BC3C29FF", "#E18727FF"),
           break.time.by = 100,
           surv.median.line = "hv",
           xlim = c(0, 2200),
           xlab = "Time (Days)",
           risk.table = TRUE,
           legend.labs = c("Plasma cell high", "APC high / T-cell low", "Plasma cell low"),
           # cumevents = TRUE,
           pval.coord = c(1200,0.75),
           risk.table.col = "strata",
           risk.table.height = 0.5#Useful when you have multiple groups
)

modDistanceMatrix <- distanceMatrix
rownames(modDistanceMatrix) <- names(group)
colnames(modDistanceMatrix) <- names(group)
anno <- data.frame(group=group)
rownames(anno) <- names(group)

pheatmap(distanceMatrix, treeheight_row = 0, treeheight_col = 0)

ggSil <- fviz_silhouette(sil, 
                         legend = "none",
                         ticks = FALSE,
                         tickslab = FALSE,
                         palette = c("#0072B5FF","#BC3C29FF", "#E18727FF"),
                         orientation = "horizontal",
                         ggtheme = theme_pubr())
ggSil + theme(plot.title = element_text(hjust = 0.5))

ggSil + coord_flip() + scale_color_manual(labels = c("T999", "T888", "XXX"), values = c("#0072B5FF","#BC3C29FF", "#E18727FF"))
# scale_color_nejm() +


### Looking for inter-group differences
require(ggpubr)
require(ggsci)
require(reshape2)
require(dplyr)

dzCibersortIntersectConcat <- cbind(dzCibersortIntersect, group)
tmpCB <- melt(dzCibersortIntersectConcat, id.vars = "group")
tmpCB$group <- as.character(paste("Group",tmpCB$group,sep=""))
tmpCB$variable <- as.factor(gsub("\\.", " ", tmpCB$variable))
head(tmpCB)
dim(tmpCB)

# cibersortComparisons <- list( c("Group1", "Group2"), c("Group1", "Group3"), c("Group2", "Group3"))
# tmpCB %>%
#   filter(variable %in% c("Macrophages M2", "T cells gamma delta")) %>%
#   ggboxplot(x = "group", y = "value", 
#             color = "group",
#             facet.by = "variable",
#             add = "jitter",
#             # palette = "lancet",
#             xlab = "Immune Cell Types", ylab = "CIBERSORT Normalized Abundance", 
#             title = "Comparison of Immune Cell Abundances Among STAD Immuno-Oncology-Microbiome (IOM) Subtypes",
#             legend = "right",
#             legend.title = "IOM Group",
#             font.label = list(size = 14, face = "bold")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   rotate_x_text(angle = 45) +
#   scale_color_nejm() +
#   stat_compare_means(mapping = aes(label = "color"), 
#                      # method = "anova",
#                      comparisons = cibersortComparisons)#,
#                      # label = "p.signif",
#                      # ref.group = "Group2",
#                      # label.y = -0.01)

cibersortComparisons <- list( c("Group1", "Group3"), c("Group1", "Group2"), c("Group2", "Group3"))
immuneCellListofInterest <- c("Neutrophils",
                              "Dendritic cells activated",
                              "Eosinophils",
                              "Macrophages M0",
                              "Macrophages M2",
                              "Plasma cells",
                              "T cells CD8",
                              "T cells follicular helper",
                              "T cells regulatory  Tregs "
                              )
tmpCBPseudo <- tmpCB
# tmpCBPseudo$value <- log((tmpCB$value+0.000001))
tmpCB %>%
  filter(!(variable %in% c("T cells CD4 naive", "T cells gamma delta"))) %>%
  # filter(variable %in% immuneCellListofInterest) %>%
  ggboxplot(x = "group", y = "value", 
            color = "group",
            facet.by = "variable",
            add = "jitter",
            # ylim = c(-20, 5),
            ylim = c(-0.01, 0.7),
            # palette = "lancet",
            xlab = "IOM STAD Subtype Group", ylab = "CIBERSORT Relative Abundance", 
            title = "Comparison of Immune Cell Abundances Among STAD Immuno-Oncology-Microbiome (IOM) Subtypes",
            legend = "none",
            legend.title = "IOM Group",
            font.label = list(size = 20, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(angle = 45) +
  # scale_color_nejm(labels = c("T999", "T888", "XXX", values = c("#0072B5FF","#BC3C29FF", "#E18727FF"))) +
  scale_color_manual(labels = c("T999", "T888", "XXX"), values = c("#0072B5FF","#BC3C29FF", "#E18727FF")) +
  scale_x_discrete(labels=c("Group1" = "Plasma cell high",
                            "Group2" = "APC high / T-cell low",
                            "Group3" = "Plasma cell low")) +
  rotate_x_text(angle = 30) +
  # yscale("log10", .format = FALSE) +
  stat_compare_means(comparisons = cibersortComparisons,
                     label = "p.signif",
                     label.y = c(0.5, 0.60, 0.65))

M <- 1795
pv <- sapply(1:M, function(i){
  mydataframe <- data.frame(y=dzDataIntersect[,i], ig=group)
  fit <- aov(y ~ ig, data=mydataframe)
  summary(fit)[[1]][["Pr(>F)"]][1]
})
names(pv) <- colnames(dzDataIntersect)
pVal <- data.frame(rawPVal = pv, pAdj = p.adjust(pv, method = "BH"))
pValOrdered <- pVal[order(pVal$pAdj),]

topXMicrobeIOMPlotSTAD <- function(numToPlot, groupVar = group){
  topXMicrobeNames <- rownames(pValOrdered)[1:numToPlot]
  topXMicrobeData <- data.frame(cbind(dzDataIntersect[,topXMicrobeNames], group= group))
  topXMicrobeData.melted <- melt(topXMicrobeData, id.vars = c("group"))
  topXMicrobeData.melted$variable <- factor(ldply(strsplit(as.character(topXMicrobeData.melted$variable),"g__*"))[[2]])
  topXMicrobeData.melted$group <- as.character(paste("Group",topXMicrobeData.melted$group,sep=""))
  head(topXMicrobeData.melted)
  
  microbeComparisons <- list( c("Group1", "Group3"), c("Group1", "Group2"), c("Group2", "Group3"))
  topXMicrobeData.melted %>%
    ggboxplot(x = "group", y = "value", 
              color = "group",
              facet.by = "variable",
              add = "jitter",
              ylim = c(-6, 15),
              # ylim = c(-0.01, 0.7),
              # palette = "lancet",
              xlab = "IOM STAD Subtype Group", ylab = "SNM Normalized Abundance", 
              title = "Comparison of Microbial Abundances Among STAD Immuno-Oncology-Microbiome (IOM) Subtypes",
              legend = "none",
              legend.title = "IOM Group",
              font.label = list(size = 14, face = "bold")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(labels = c("T999", "T888", "XXX"), values = c("#0072B5FF","#BC3C29FF", "#E18727FF")) +
    rotate_x_text(angle = 30) +
    scale_x_discrete(labels=c("Group1" = "Plasma cell high",
                              "Group2" = "APC high / T-cell low",
                              "Group3" = "Plasma cell low")) +
    stat_compare_means(comparisons = microbeComparisons,
                       label = "p.signif")
}

topXMicrobeIOMPlotSTAD(20)




# ###
# 
# result2 <- ExecuteSNF.CC(datasets = list(t(dzDataIntersect), 
#                                         t(dzCibersortIntersect)), 
#                          clusterNum = 4, K = 25, alpha = 0.5, t = 20,
#                          maxK = 10, pItem = 0.8, reps = 50)
# sil2=silhouette_SimilarityMatrix(result2$group, result2$distanceMatrix)
# plot(sil2)
# group2=result2$group
# distanceMatrix2=result2$distanceMatrix
# p_value2=survAnalysis(mainTitle="Testing",time = dzMetadataIntersect$days_to_death,
#                      status = rep(1, times = length(dzMetadataIntersect$days_to_death)),
#                      group2,
#                      distanceMatrix2,similarity=TRUE)

# save(result, sil, group, distanceMatrix, p_value, file = "snf.cc.4clusterCibersortCorrected.RData")
# load("snf.cc.4cluster.RData")

# sigClustRes <- sigclustTest(Data = t(coadCibersortIntersect),
#                             group = result$group,
#                             nsim = 100,
#                             icovest = 1)
# sigClustRes

#-------------------------- Apply clusters to rest of data --------------------------#
# 
# coadMetadataQC <- metadataSamplesAllQC[metadataSamplesAllQC$disease_type == "Colon Adenocarcinoma",]
# coadDataQC <- snmDataSampleTypeWithExpStrategy[rownames(coadMetadataQC),]
# coadDataQCPT <- coadDataQC[rownames(coadMetadataQC[coadMetadataQC$sample_type == "Primary Tumor",]),]
# 
# sampleIntersectQC <- Reduce(intersect,
#                           list(rownames(coadDataQCPT),
#                                # rownames(mhc2ConservativeWithDummiesMergedCCSample),
#                                rownames(cibersortTCGAMergedCCSample)#,
#                                # rownames(mutationsTCGAMergedCCSample)
#                           ))
# length(sampleIntersectQC)
# 
# coadDataQCIntersect <- coadDataQCPT[sampleIntersectQC,]
# coadCibersortQCIntersect <- cibersortTCGAMergedCCSample[sampleIntersectQC,]
# coadMutationsQCIntersect <- mutationsTCGAMergedCCSample[sampleIntersectQC,]
# coadMetadataQCIntersect <- coadMetadataQC[sampleIntersectQC,]
# 
# trainData <- list(coadDataIntersect, 
#                   as.matrix(coadCibersortIntersect))
# testData <- list(coadDataQCIntersect[!(rownames(coadDataQCIntersect) %in% rownames(coadDataIntersect)),],
#                  as.matrix( coadCibersortQCIntersect[!(rownames(coadDataQCIntersect) %in% rownames(coadDataIntersect)),] ) )
# trainGroups = as.vector(group)
# 
# str(trainData)
# str(testData)
# 
# newLabel = groupPredict(train = trainData,
#                         test = testData,
#                         groups = trainGroups,
#                         K = 25,
#                         alpha = 0.5,
#                         t = 20,
#                         method = 1) # method = 0 means to use local and global consistency; method = 1 means to use label propagation
# newLabel
# 
# #-------------------------- Rank features by NMI --------------------------#
# 
# K = 25
# alpha = 0.5
# t = 20
# Dist1 = dist2(as.matrix(coadDataIntersect),as.matrix(coadDataIntersect))
# Dist2 = dist2(as.matrix(coadCibersortIntersect),as.matrix(coadCibersortIntersect))
# W1 = affinityMatrix(Dist1, K, alpha)
# W2 = affinityMatrix(Dist2, K, alpha)
# W = SNF(list(W1,W2), K, t)
# estimateNumberOfClustersGivenGraph(W, NUMC = 2:10)
# 
# source("rankFeaturesByNMI_GP.R")
# featureNMI <- rankFeaturesByNMI_GP(data = list(coadDataIntersect,coadCibersortIntersect), 
#                                    W = W,
#                                    pickClusterEst = 4) 
# taxaFeatureNMI <- data.frame(Taxa = colnames(coadDataIntersect), 
#                                  NMI = featureNMI[[1]][[1]],
#                                  Rank = featureNMI[[2]][[1]])
# taxaFeatureNMIOrdered <- taxaFeatureNMI[order(taxaFeatureNMI[,"Rank"]),]
# cibersortFeatureNMI <- data.frame(CellType = colnames(coadCibersortIntersect), 
#                                       NMI = featureNMI[[1]][[2]],
#                                       Rank = featureNMI[[2]][[2]])
# cibersortFeatureNMIOrdered <- cibersortFeatureNMI[order(cibersortFeatureNMI[,"Rank"]),]
# 
# write.csv(taxaFeatureNMIOrdered, file = "taxaFeatureNMIOrdered.csv")
# write.csv(cibersortFeatureNMIOrdered, file = "cibersortFeatureNMIOrdered.csv")
