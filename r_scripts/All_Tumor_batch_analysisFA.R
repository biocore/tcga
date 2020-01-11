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

