# Contamination_simulation_FA.R
# Author: Greg Poore
# Date: Aug 16 2019
# Purpose: To decontaminate TCGA data and see if/how it affects the machine learning analyses

# Load dependencies
require(ggplot2)
require(ggsci)
require(limma)
require(Glimma)
require(edgeR)
require(dplyr)
require(doMC)

numCores <- detectCores()
registerDoMC(cores=numCores)

#---------------------------------------------
# Using plate number for batch to decontaminate
#---------------------------------------------

# Function for extracting last n characters from R string
# URL: https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

# NB: metadataSamplesAllQCCGC IS A METADATA FILE OF TCGA SAMPLES CONTAINING ALIQUOT IDS FOR ALL SAMPLES
# IN A COLUMN CALLED "aliquot_id"
tmp <- as.character(metadataSamplesAllQCCGC$aliquot_id)
metadataSamplesAllQCCGC$PlateCenter <- factor(substrRight(tmp, 7))
# NB: SINCE DECONTAM ESSENTIALLY PERFORMS LINEAR REGRESSION BETWEEN READ FRACTIONS AND 
# ANALYTE CONCENTRATIONS, AT LEAST 10 SAMPLES ARE REQUIRED PER PLATE-CENTER COMBINATION
# TO BE PROCESSED FOR IDENTIFYING PUTATIVE CONTAMINANTS. NOTE ALSO THAT ANY CONTAMINANT
# IDENTIFIED IN ANY ONE PLATE-CENTER BATCH WILL BE REMOVED FROM THE WHOLE DATASET
booleanPlateCenter <- as.logical(table(metadataSamplesAllQCCGC$PlateCenter)>10)
sufficientPlateCenter <- names(table(metadataSamplesAllQCCGC$PlateCenter))[booleanPlateCenter]
metadataSamplesAllQCCGC$PlateCenterFlag <- (metadataSamplesAllQCCGC$PlateCenter %in% sufficientPlateCenter)
metadataSamplesAllQCCGC_PlateCenterSubset <- droplevels(metadataSamplesAllQCCGC[metadataSamplesAllQCCGC$PlateCenterFlag,])

# NB: vbContaminatedDataQC CONTAINS RAW TCGA DATA FROM KRAKEN AND THE SPIKED PSEUDOCONTAMINANTS
vbContaminatedDataQC_PlateCenterSubset <- vbContaminatedDataQC[rownames(metadataSamplesAllQCCGC_PlateCenterSubset),]

# Decontam
require(decontam)
countData <- vbContaminatedDataQC_PlateCenterSubset
countMetadata <- metadataSamplesAllQCCGC_PlateCenterSubset

contamdf.freq.plateCenter <- isContaminant(seqtab = as.matrix(countData), 
                                             conc = countMetadata$aliquot_concentration, 
                                             method = "frequency", 
                                             batch = countMetadata$PlateCenter,
                                             threshold = 0.1) # DEFAULT VALUE IS 0.1
save(contamdf.freq.plateCenter, file = "contamdf.freq.plateCenter.RData")
table(contamdf.freq.plateCenter$contaminant)


contamSum <- colSums(as.matrix(vbContaminatedDataQC)[,contamdf.freq.plateCenter$contaminant])
sum(contamSum)/sum(colSums(as.matrix(vbContaminatedDataQC))) #--> 0.02140984
plateCenterSplitGenera <- strsplit(rownames(contamdf.freq.plateCenter)[contamdf.freq.plateCenter$contaminant],split = ".g__")
generaNamesPlateCenter <- sapply(plateCenterSplitGenera, "[",2)
tail(generaNamesPlateCenter)
# NB: This ^ cut off the spiked contaminants, which should be included for downstream processing. FIX WITH FOLLOWING CODE:
generaNamesPlateCenter[(length(generaNamesPlateCenter)-1):(length(generaNamesPlateCenter))] <- 
  c("contaminant4RandomSpikesHarvard", "contaminant5RandomSpikes1000")

# RESULTANT PLATE-CENTER DECONTAMINATED DATA:
vbContaminatedDataQCPlateCenterContamRemoved <- vbContaminatedDataQC[,!(generaNamesAllWithSpikedContaminants %in% 
                                                                            c(generaNamesPlateCenter, 
                                                                              contamReagents$Genera,
                                                                              generaCrossContamLikely))]

# Evaluate how much data was lost
1-sum(colSums(vbContaminatedDataQCPlateCenterContamRemoved))/sum(colSums(vbContaminatedDataQC))
# = 0.6714914