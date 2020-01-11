#-------------------------------#
# Finding contaminants in model features
#-------------------------------#

require(readtext)
require(stringr)

# NOTES:
# - NON-ZERO FEATURES FOR EACH DATA TYPE AND COMPARISON WERE SAVED AS CSV FILES DURING MACHINE LEARNING
# - THESE CSV FILES ARE THEN LOADED PER DATA TYPE (SHOWN BELOW)
# - AFTER LOADING, SPIKED PSEUDO-CONTAMINANTS ARE IDENTIFIED PER CSV FILE AND THEN THEIR FEATURE IMPORTANCE SCORES SUMMED PER CSV FILE
# - THE SUMMED FEATURE IMPORTANCE SCORE FROM THE SPIKED PSEUDO-CONTAMINANTS IS THEN DIVDED BY THE TOTAL SUMMED FEAT. IMP. SCORE
# - THIS PROVIDES AN ESTIMATE OF THE CONTRIBUTION OF SPIKED PSEUDO-CONTAMINANTS TO THE MODEL'S PREDICTIONS
# - A HIGH PROPORTION CONTRIBUTION OF SPIKED PSEUDO-CONTAMINANTS WOULD SUGGEST AN UNRELIABLE MODEL
# - ANOTHER APPROACH COULD USE THE RANKINGS TO REMOVE ANY TAXA RANKED LOWER THAN SPIKED PSEUDO-CONTAMINANTS

## Likely Contam
filePathLikelyContam <- "~/Downloads/features__snmDataLikelyContamRemovedFA/*.csv"
modelFeatLikelyContam <- readtext(filePathLikelyContam)
modelFeatLikelyContam$group <- sapply(strsplit(modelFeatLikelyContam$doc_id, split = "--\\sF"),'[',1)
data.frame(tail(modelFeatLikelyContam))

modelFeatLikelyContam %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(group_sum = sum(varImp), count = dplyr::n()) %>%
  data.frame() -> modelFeatLikelyContamSUM
modelFeatLikelyContamSUM %>% write.csv(file = "pseudoSUMFEATURESmodelFeatLikelyContam_Dec3.csv")

contamModelFeatLikelyContam <- data.frame(modelFeatLikelyContam[grep(pattern = ".*contaminant.*", modelFeatLikelyContam$text),])
contamModelFeatLikelyContam$group <- sapply(strsplit(contamModelFeatLikelyContam$doc_id, split = "--\\sF"),'[',1)
write.csv(contamModelFeatLikelyContam,
          file = "pseudoContamModelFeatLikelyContam_Dec3.csv")

contamModelFeatLikelyContam %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(group_sum = sum(varImp), count = dplyr::n()) %>%
  data.frame() -> contamModelFeatLikelyContamSUM

left_join(modelFeatLikelyContamSUM, contamModelFeatLikelyContamSUM, by = "group") %>%
  mutate(percentContam = group_sum.y*100/group_sum.x) -> percentContamLikelyContam
percentContamLikelyContam %>% write.csv("percentContamLikelyContam_Dec3.csv")

## All Putative
filePathAllPutative <- "~/Downloads/features__snmDataAllPutativeContamRemovedFA/*.csv"
modelFeatAllPutative <- readtext(filePathAllPutative)
modelFeatAllPutative$group <- sapply(strsplit(modelFeatAllPutative$doc_id, split = "--\\sF"),'[',1)

modelFeatAllPutative %>%
  group_by(group) %>%
  dplyr::summarise(group_sum = sum(varImp), count = dplyr::n()) %>%
  data.frame() -> modelFeatAllPutativeSUM
modelFeatAllPutativeSUM %>% write.csv(file = "pseudoSUMFEATURESmodelFeatAllPutative_Dec3.csv")

contamModelFeatAllPutative <- data.frame(modelFeatAllPutative[grep(pattern = ".*contaminant.*", modelFeatAllPutative$text),])
contamModelFeatAllPutative$group <- sapply(strsplit(contamModelFeatAllPutative$doc_id, split = "--\\sF"),'[',1)
write.csv(contamModelFeatAllPutative,
          file = "pseudoContamModelFeatAllPutative_Dec3.csv")

contamModelFeatAllPutative %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(group_sum = sum(varImp), count = dplyr::n()) %>%
  data.frame() -> contamModelFeatAllPutativeSUM

left_join(modelFeatAllPutativeSUM, contamModelFeatAllPutativeSUM, by = "group") %>%
  mutate(percentContam = group_sum.y*100/group_sum.x) -> percentContamAllPutative
percentContamAllPutative %>% write.csv("percentContamAllPutative_Dec3.csv")

## PlateCenter
filePathPlateCenter <- "~/Downloads/features__snmDataPlateCenterContamRemovedFA/*.csv"
modelFeatPlateCenter <- readtext(filePathPlateCenter)
modelFeatPlateCenter$group <- sapply(strsplit(modelFeatPlateCenter$doc_id, split = "--\\sF"),'[',1)

modelFeatPlateCenter %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(group_sum = sum(varImp), count = dplyr::n()) %>%
  data.frame() -> modelFeatPlateCenterSUM
modelFeatPlateCenterSUM %>% write.csv(file = "pseudoSUMFEATURESmodelFeatPlateCenter_Dec3.csv")

contamModelFeatPlateCenter <- data.frame(modelFeatPlateCenter[grep(pattern = ".*contaminant.*", modelFeatPlateCenter$text),])
contamModelFeatPlateCenter$group <- sapply(strsplit(contamModelFeatPlateCenter$doc_id, split = "--\\sF"),'[',1)
write.csv(contamModelFeatPlateCenter,
          file = "pseudoContamModelFeatPlateCenter_Dec3.csv")

contamModelFeatPlateCenter %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(group_sum = sum(varImp), count = dplyr::n()) %>%
  data.frame() -> contamModelFeatPlateCenterSUM

left_join(modelFeatPlateCenterSUM, contamModelFeatPlateCenterSUM, by = "group") %>%
  mutate(percentContam = group_sum.y*100/group_sum.x) -> percentContamPlateCenter
percentContamPlateCenter %>% write.csv("percentContamPlateCenter_Dec3.csv")

## MostStringent
filePathMostStringent <- "~/Downloads/features__snmDataMostStringentContamRemovedFA/*.csv"
modelFeatMostStringent <- readtext(filePathMostStringent)
modelFeatMostStringent$group <- sapply(strsplit(modelFeatMostStringent$doc_id, split = "--\\sF"),'[',1)

modelFeatMostStringent %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(group_sum = sum(varImp), count = dplyr::n()) %>%
  data.frame() -> modelFeatMostStringentSUM
modelFeatMostStringentSUM %>% write.csv(file = "pseudoSUMFEATURESmodelFeatMostStringent_Dec3.csv")

contamModelFeatMostStringent <- data.frame(modelFeatMostStringent[grep(pattern = ".*contaminant.*", modelFeatMostStringent$text),])
contamModelFeatMostStringent$group <- sapply(strsplit(contamModelFeatMostStringent$doc_id, split = "--\\sF"),'[',1)
write.csv(contamModelFeatMostStringent,
          file = "pseudoContamModelFeatMostStringent_Dec3.csv")

contamModelFeatMostStringent %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(group_sum = sum(varImp), count = dplyr::n()) %>%
  data.frame() -> contamModelFeatMostStringentSUM

left_join(modelFeatMostStringentSUM, contamModelFeatMostStringentSUM, by = "group") %>%
  mutate(percentContam = group_sum.y*100/group_sum.x) -> percentContamMostStringent
percentContamMostStringent %>% write.csv("percentContamMostStringent_Dec3.csv")

percentContamLikelyContam_AllPutative <- left_join(percentContamLikelyContam, percentContamAllPutative, 
                  by = "group", suffix = c(".likelyContam",".allPutative"))
percentContamAll3 <- left_join(percentContamLikelyContam_AllPutative, percentContamMostStringent, by = "group")

filter(percentContamAll3, grepl("*Primary Tumor[[:space:]]+$",percentContamAll3$group)) %>%
  write.csv(file = "contamContribution_PT_Dec3.csv")

filter(percentContamAll3, grepl("*Primary Tumor vs Solid Tissue Normal[[:space:]]+$",percentContamAll3$group)) %>%
  write.csv(file = "contamContribution_TvsN_Dec3.csv")

filter(percentContamAll3, grepl("*Blood Derived Normal[[:space:]]+$",percentContamAll3$group)) %>%
  write.csv(file = "contamContribution_BDN_Dec3.csv")

filter(percentContamAll3, grepl("*Stage I vs IV[[:space:]]+$",percentContamAll3$group)) %>%
  write.csv(file = "contamContribution_Stage_Dec3.csv")

#-------------------------------#
# Extended figs -- bubble plots
#-------------------------------#
require(reshape2)
library(dplyr)
require(ggpubr)
require(ggsci)
# contamContribData <- read.csv("extFigBubblePlotPseudoContaminantContribution.csv", stringsAsFactors = FALSE)
# contamContribDataPercent <- contamContribData[,c(1:4)]
# contamContribDataPercent.melt <- melt(contamContribDataPercent, id.vars = "DT")
# 
# ggplot(contamContribDataPercent.melt,
#        aes(x = DT, y = variable, 
#            size = value, label = value, fill = value)) +
#   geom_point(aes(fill = value), shape = 21) +
#   geom_text(size = 5, nudge_y = -0.35) +
#   scale_size(range = c(10, 30), guide = F) +
#   scale_y_discrete(labels = c("Likely\nContaminants\nRemoved",
#                               "All Putative\nContaminants\nRemoved",
#                               "Most Stringent\nFiltering")) +
#   theme_pubr() + scale_fill_gradient(low="blue",high="darkorange", name = "Percentage Contribution of Pseudo-Contaminants to Model Predictions") +
#   theme(axis.title.x=element_blank(), axis.title.y=element_blank()) -> contamContribDataPercentPlot
# # theme(panel.grid.major = element_line(linetype = 2, color = "black"))
# ggsave(contamContribDataPercentPlot, width = 18, units = "in",
#        filename = "contamContribDataPercentPlot.svg",
#        dpi = "retina")

## PT Only
contamContribDataPTOnly <- read.csv("contamContribution_PT_Dec3_Formatted.csv", stringsAsFactors = FALSE)
# contamContribDataPercentPTOnly <- contamContribDataPTOnly[,c(2:5)]
contamContribDataPercentPTOnly.melt <- melt(contamContribDataPTOnly, id.vars = "DT")
contamContribDataPercentPTOnly.melt$value <- round(contamContribDataPercentPTOnly.melt$value, digits = 2)

write.csv(contamContribDataPercentPTOnly.melt, file = "./AARF/edf7a.csv")

ggplot(contamContribDataPercentPTOnly.melt,
       aes(x = DT, y = variable, 
           size = value, label = value, fill = value)) +
  geom_point(aes(fill = value), shape = 21) +
  geom_text(size = 4, nudge_y = -0.3) +
  scale_size(range = c(5, 25), guide = F) +
  scale_y_discrete(labels = c("Likely\nContaminants\nRemoved",
                              "Decontamination\nBy Plate-Center\nCombination",
                              "All Putative\nContaminants\nRemoved",
                              "Most Stringent\nFiltering")) +
  ggpubr::rotate_x_text() +
  # coord_flip() +
  theme_pubr() + scale_fill_gradient(low="blue",high="darkorange", name = "Percentage Contribution of Pseudo-Contaminants to Model Predictions") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) -> contamContribDataPercentPlotPTOnly
# theme(panel.grid.major = element_line(linetype = 2, color = "black"))

ggsave(contamContribDataPercentPlotPTOnly, 
       width = 18, 
       units = "in",
       filename = "contamContribDataPercentPlotPTOnly_Dec3.svg",
       dpi = "retina")

## BDN Only
contamContribDataBDNOnly <- read.csv("contamContribution_BDN_Dec3_Formatted.csv", stringsAsFactors = FALSE)
# contamContribDataPercentPTOnly <- contamContribDataPTOnly[,c(2:5)]
contamContribDataBDNOnly.melt <- melt(contamContribDataBDNOnly, id.vars = "DT")
contamContribDataBDNOnly.melt$value <- round(contamContribDataBDNOnly.melt$value, digits = 2)

write.csv(contamContribDataBDNOnly.melt, file = "./AARF/edf7b.csv")

ggplot(contamContribDataBDNOnly.melt,
       aes(x = DT, y = variable, 
           size = value, label = value, fill = value)) +
  geom_point(aes(fill = value), shape = 21) +
  geom_text(size = 4, nudge_y = -0.3) +
  scale_size(range = c(5, 25), guide = F) +
  scale_y_discrete(labels = c("Likely\nContaminants\nRemoved",
                              "Decontamination\nBy Plate-Center\nCombination",
                              "All Putative\nContaminants\nRemoved",
                              "Most Stringent\nFiltering")) +
  ggpubr::rotate_x_text() +
  # coord_flip() +
  theme_pubr() + scale_fill_gradient(low="blue",high="darkorange", name = "Percentage Contribution of Pseudo-Contaminants to Model Predictions") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) -> contamContribDataPercentPlotBDNOnly
# theme(panel.grid.major = element_line(linetype = 2, color = "black"))

ggsave(contamContribDataPercentPlotBDNOnly, 
       width = 18, 
       units = "in",
       filename = "contamContribDataPercentPlotBDNOnly_Dec3.svg",
       dpi = "retina")
