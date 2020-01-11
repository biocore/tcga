#-----------------------------------------------
#           VSNM Normalization                 #
#-----------------------------------------------

vsnm <- function(){
  ## Load packages ##
  require(limma)
  require(edgeR)
  require(dplyr)
  require(snm)
  require(doMC)
  require(tibble)
  require(gbm)
  
  numCores <- detectCores()
  registerDoMC(cores=numCores)
  
  qcMetadata <- INPUT_METADATA # ADAPT THIS AS NEEDED
  qcData <- INPUT_COUNT_DATA # ADAPT THIS AS NEEDED
  
  # Set up design matrix
  covDesignNorm <- model.matrix(~0 + disease_type_consol +
                                    host_age + # host_age should be numeric
                                    sex, # sex should be a factor
                                  data = qcMetadata)
  
  # Check row dimensions
  dim(covDesignNorm)[1] == dim(qcData)[1]
  
  print(colnames(covDesignNorm))
  # The following corrects for column names that are incompatible with downstream processing
  colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
  print(colnames(covDesignNorm))
  
  # Set up counts matrix
  counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
  
  # Quantile normalize and plug into voom
  dge <- DGEList(counts = counts)
  vdge <<- voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE, 
                                              normalize.method="quantile")
  
  # List biological and normalization variables in model matrices
  bio.var <- model.matrix(~disease_type_consol,
                          data=qcMetadata)
  
  adj.var <- model.matrix(~host_age +
                              sex,
                            data=qcMetadata)
  
  colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
  colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
  print(dim(adj.var))
  print(dim(bio.var))
  print(dim(t(vdge$E)))
  print(dim(covDesignNorm))
  
  snmDataObjOnly <- snm(raw.dat = vdge$E, 
                        bio.var = bio.var, 
                        adj.var = adj.var, 
                        rm.adj=TRUE,
                        verbose = TRUE,
                        diagnose = TRUE)
  snmData <<- t(snmDataObjOnly$norm.dat)
  
}