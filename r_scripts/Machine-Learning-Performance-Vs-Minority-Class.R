#-------------------------------#
# ML specificity vs minority class size
#-------------------------------#

# IMPORT DATA (DONE HERE FROM A CSV FILE)
allPerfNov11 <- read.csv(" ADD_YOUR_FILE_PATH_HERE ",
                         stringsAsFactors = FALSE,
                         strip.white = TRUE)

ptPerfNov11 <- allPerfNov11[allPerfNov11$sampleType == "Primary Tumor",]

# ESTIMATE LINEAR REGRESSION AND FIT STATISTICS
summary(lm(auroc_fullData ~ log10(minorityClassSize), ptPerfNov11))
summary(lm(aupr_fullData ~ log10(minorityClassSize), ptPerfNov11))

# MAIN REGRESSION FIT AND PLOTTING FUNCTION
plotMinorityClassSizeVsPerf <- function(valData, 
                                        fileName = "tmp",
                                        auprFlag = TRUE,
                                        # xlowerlim = 0.5, 
                                        ylowerlim = 0.5, 
                                        names = NULL,
                                        title = "Classifier performance vs minority class size"){
  require(ggplot2)
  if(auprFlag){
    df <- valData[,-grep("AUROC",colnames(valData),ignore.case = TRUE)]
    colnames(df)[grep("AUPR",colnames(df),ignore.case = TRUE)] <- "perf"
    typeString <- "AUPR"
  } else{
    df <- valData[,-grep("AUPR",colnames(valData),ignore.case = TRUE)]
    colnames(df)[grep("AUROC",colnames(df),ignore.case = TRUE)] <- "perf"
    typeString <- "AUROC"
  }
  print(df)
  fit <- lm(perf ~ log10(minorityClassSize), df)
  print(summary(fit))
  write.csv(df, file = paste0("./AARF/",fileName))
  p <- ggplot(df, aes(x = log10(minorityClassSize), y = perf)) +
    # guides(colour = guide_legend(override.aes = list(size = 8)))+
    # theme(legend.key=element_rect(fill=NA)) +
    geom_point(alpha = 0.8) + 
    theme_pubr() +
    geom_smooth(method = 'lm', formula = y ~ x) +
    ylim(c(ylowerlim,1)) +
    # geom_abline(linetype = 2) +
    labs(x = expression(Log[10](Minority~Class~Size)),
         y = paste(typeString,"of classifiers trained/tested\n on full VSNM data"),
         title = title) +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5),
          text = element_text(size=20), axis.text = element_text(color = "black"))
  if(is.null(names)){
    p <- p
  }
  else{
    require(ggrepel)
    p <- p + geom_label_repel(aes(label=DT),
                              force = 1,
                              segment.alpha = 0.8,
                              label.padding = 0.1,
                              hjust=-0.15, vjust=0.15, size = 5, 
                              segment.color = 'grey50', show.legend = FALSE)
  }
  print(p)
  return(p)
}

# PLOT AUROC VS MINORITY CLASS
ptPerfPlotAUROC <- plotMinorityClassSizeVsPerf(ptPerfNov11[,c("DT",
                                                         "minorityClassSize",
                                                         "auroc_fullData",
                                                         "aupr_fullData")],
                                               fileName = "minorityClassAUROC.csv",
                                          auprFlag = FALSE, ylowerlim = 0, names = TRUE, title = "")

# SAVE PLOT
ggsave(plot = ptPerfPlotAUROC,
       filename = "./Clinical Validation Plots/minorityclassPt1vsAllPerfPlot_Dec2_AUROC.svg",
       width = 10,
       height = 5,
       units = "in")

# PLOT AUPR VS MINORITY CLASS
ptPerfPlotAUPR <- plotMinorityClassSizeVsPerf(ptPerfNov11[,c("DT",
                                            "minorityClassSize",
                                            "auroc_fullData",
                                            "aupr_fullData")],
                                            fileName = "minorityClassAUPR.csv",
                            auprFlag = TRUE, ylowerlim = 0, names = TRUE, title = "")

# SAVE PLOT
ggsave(plot = ptPerfPlotAUPR,
       filename = "./Clinical Validation Plots/minorityclassPt1vsAllPerfPlot_Dec2_AUPR.svg",
       width = 10,
       height = 5,
       units = "in")

