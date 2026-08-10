################################################################################
################################################################################
# Scripts for plotting assocaited with predictive classification
################################################################################
################################################################################

makePredStatPlt <- function(inpPredDfLst, pltCols=NULL, histCols = NULL, remHist = NULL, fileAdd = NULL, binBool = NULL, multParamBool = NULL, fullDfBool=NULL){
  if(is.null(pltCols)){
    pltCols <- c("SID","Histotype", 
                 "Sensitivity", "Specificity", 
                 "Accuracy","F1Score",  
                 "Low", "High",
                 "AUC", "PRAUC")
  }
  if(is.null(histCols)){
    histCols <- c("CCC", "EC", "HGSC", "MC")
  }
  if(!is.null(remHist)){
    histCols <- histCols[!histCols %in% remHist]
  }
  pltDf <- data.frame(matrix(nrow=0, ncol=length(pltCols)))
  colnames(pltDf) <- pltCols
  if(is.null(binBool)){
    colnames(pltDf)[grep("Low", colnames(pltDf))] <- "Acc_Lower"
    colnames(pltDf)[grep("High", colnames(pltDf))] <- "Acc_Higher"
  }else{
    colnames(pltDf)[grep("Low", colnames(pltDf))] <- "CI_Lower"
    colnames(pltDf)[grep("High", colnames(pltDf))] <- "CI_Higher"
  }
  grepVec <- paste(histCols, collapse="|")
  
  for(i in 1:length(inpPredDfLst)){
    inpPredDf <- inpPredDfLst[[i]]
    if(!"SID" %in% colnames(inpPredDf)){
      inpPredDf$SID <- names(inpPredDfLst)[i]
    }
    if(length(which(colnames(inpPredDf) %in% "AUPRC"))>0){
      colnames(inpPredDf)[which(colnames(inpPredDf) %in% "AUPRC")] <- "PRAUC"
    }
    tmpSid <- names(inpPredDfLst)[i]
    colIndVec <- c()
    for(z in 1:length(pltCols)){
      tmpCol <- pltCols[z]
      tmpInds <- grep(tmpCol, colnames(inpPredDf))
      colIndVec <- append(colIndVec, tmpInds)
    }
    colIndVec <- colIndVec[!duplicated(colIndVec)]
    for(j in 1:nrow(inpPredDf)){
      tmpRank <-inpPredDf[j,]
      tmpVals <- tmpRank[, colIndVec]
      hVals <- tmpVals[grep(grepVec, colnames(tmpVals))]
      nhVals <- tmpVals[!colnames(tmpVals) %in% colnames(hVals)]
      for(k in 1:length(histCols)){
        tmpH <- histCols[k]
        hCols <- hVals[grep(tmpH, colnames(hVals))]
        hCols <- unlist(hCols)
        nhVals <- unlist(nhVals)
        hCols <- append(hCols, nhVals)
        hCols <- append(hCols, tmpH)
        names(hCols)[length(hCols)] <- "Histotype"
        names(hCols) <- gsub(paste(tmpH, "_", sep=""), "", names(hCols))
        nRow <- hCols[match(pltCols, names(hCols))]
        pltDf[nrow(pltDf)+1,] <- nRow
      }
    }
  }
  if(!is.null(multParamBool)){
    if(is.null(binBool)){
      pltDf$Histotype <- "Model"
    }
  }
  skipCols <- c("Histotype", "SID")
  pltDfLong <-  pltDf %>% 
    pivot_longer(
      cols = pltCols[!pltCols %in% skipCols], 
      names_to = "Measurement",
      values_to = "Value"
    )
  if(length(grep("Lower", pltDfLong$Measurement))>0){
    if(is.null(binBool)){
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% "Acc_Lower")] <- "Acc_Range"
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% c("Acc_Higher", "Acc_Upper"))] <- "Acc_Range"
    }else{
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% "CI_Lower")] <- "CI_Range"
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% c("CI_Higher", "CI_Upper"))] <- "CI_Range"
    }
  }
  pltDfLong <- data.frame(pltDfLong)
  pltDfLong$Value <- gsub(",", ".", pltDfLong$Value)
  pltDfLong$Value <- as.numeric(pltDfLong$Value)
  
  pltDfLong <- pltDfLong[!pltDfLong$Value %in% c(0, "0"),]
  
  if(!is.null(multParamBool)){
    tmpPlt <- ggplot(pltDfLong, aes(fill=Histotype, y=Value, x=Histotype)) + 
      geom_boxplot() +
      #geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
      xlab("Histotype") +
      ylab("Percent") + 
      ggtitle(paste("Mean values for 5-fold cross-validated iterations (n=100)\nEC excluded due to lack of HSPs" , sep="")) + 
      scale_fill_viridis(discrete = TRUE) +
      theme(text = element_text(size=24), 
            axis.text.x = element_text(size=16, face="bold"),
            legend.text=element_text(size=16),
            plot.title = element_text(hjust = 0.5)) +
      facet_grid(rows = vars(SID), cols = vars(Measurement), 
                 scales="free_y", switch = 'y')
  }else{
    pltDfLong$Histotype <- factor(pltDfLong$Histotype, levels=c("CCC", "EC", "HGSC", "MC"))
    tmpPlt <- ggplot(pltDfLong, aes(fill=Histotype, y=Value, x=Histotype)) + 
      geom_boxplot() +
      #geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
      xlab("Histotype") +
      ylab("Percent") + 
      ggtitle(paste("Mean values for 5-fold cross-validated iterations (n=100)\nEC excluded due to lack of HSPs" , sep="")) + 
      scale_fill_viridis(discrete = TRUE) +
      theme(text = element_text(size=24), 
            axis.text.x = element_text(size=16, face="bold"),
            legend.text=element_text(size=16),
            plot.title = element_text(hjust = 0.5)) +
      facet_grid(rows = vars(SID), cols = vars(Measurement), 
                 scales="free_y", switch = 'y')
  }
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(is.null(fileAdd)){
    outFile <- paste(outDir, "Pred_Stats_Plt.pdf", sep="")
  }else{
    outFile <- paste(outDir, "Pred_Stats_Plt_", fileAdd, ".pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt, width=40, height=20, units = "cm")
}

makePredEpochPlt <- function(inpPredDfLst, pltCols=NULL, histCols = NULL, remHist = NULL, fileAdd = NULL, binBool = NULL, multParamBool = NULL, fullDfBool=NULL){
  if(is.null(pltCols)){
    pltCols <- c("SID","Histotype", 
                 "Sensitivity", "Specificity", 
                 "Accuracy","F1Score",  
                 "Low", "High",
                 "AUC", "PRAUC")
  }
  if(is.null(histCols)){
    histCols <- c("CCC", "EC", "HGSC", "MC")
  }
  if(!is.null(remHist)){
    histCols <- histCols[!histCols %in% remHist]
  }
  if(!"Epoch" %in% pltCols){
    pltCols <- append(pltCols, "Epoch")
  }
  pltDf <- data.frame(matrix(nrow=0, ncol=length(pltCols)))
  colnames(pltDf) <- pltCols
  if(is.null(binBool)){
    colnames(pltDf)[grep("Low", colnames(pltDf))] <- "Acc_Lower"
    colnames(pltDf)[grep("High", colnames(pltDf))] <- "Acc_Higher"
  }else{
    colnames(pltDf)[grep("Low", colnames(pltDf))] <- "CI_Lower"
    colnames(pltDf)[grep("High", colnames(pltDf))] <- "CI_Higher"
  }
  grepVec <- paste(histCols, collapse="|")
  
  for(i in 1:length(inpPredDfLst)){
    inpPredDf <- inpPredDfLst[[i]]
    if(length(which(colnames(inpPredDf) %in% "AUPRC"))>0){
      colnames(inpPredDf)[which(colnames(inpPredDf) %in% "AUPRC")] <- "PRAUC"
    }
    tmpSid <- names(inpPredDfLst)[i]
    colIndVec <- c()
    for(z in 1:length(pltCols)){
      tmpCol <- pltCols[z]
      tmpInds <- grep(tmpCol, colnames(inpPredDf))
      colIndVec <- append(colIndVec, tmpInds)
    }
    colIndVec <- colIndVec[!duplicated(colIndVec)]
    for(j in 1:nrow(inpPredDf)){
      tmpRank <-inpPredDf[j,]
      tmpVals <- tmpRank[, colIndVec]
      hVals <- tmpVals[grep(grepVec, colnames(tmpVals))]
      nhVals <- tmpVals[!colnames(tmpVals) %in% colnames(hVals)]
      #tmpCats <- colnames(tmpVals)
      for(k in 1:length(histCols)){
        tmpH <- histCols[k]
        hCols <- hVals[grep(tmpH, colnames(hVals))]
        hCols <- unlist(hCols)
        nhVals <- unlist(nhVals)
        hCols <- append(hCols, nhVals)
        hCols <- append(hCols, tmpH)
        names(hCols)[length(hCols)] <- "Histotype"
        names(hCols) <- gsub(paste(tmpH, "_", sep=""), "", names(hCols))
        nRow <- hCols[match(pltCols, names(hCols))]
        pltDf[nrow(pltDf)+1,] <- nRow
      }
    }
  }
  if(!is.null(multParamBool)){
    pltDf$Histotype <- "Model"
  }
  skipCols <- c("Histotype", "SID", "Epoch")
  pltDfLong <-  pltDf %>% 
    pivot_longer(
      cols = pltCols[!pltCols %in% skipCols], 
      names_to = "Measurement",
      values_to = "Value"
    )
  if(!is.null(multParamBool)){
    pltDfLong <- pltDfLong[!pltDfLong$Value == 0, ]
    if(is.null(binBool)){
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% "Acc_Lower")] <- "Acc_Range"
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% c("Acc_Higher", "Acc_Upper"))] <- "Acc_Range"
    }else{
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% "CI_Lower")] <- "CI_Range"
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% c("CI_Higher", "CI_Upper"))] <- "CI_Range"
    }
  }
  pltDfLong <- data.frame(pltDfLong)
  pltDfLong$Value <- gsub(",", ".", pltDfLong$Value)
  pltDfLong$Value <- as.numeric(pltDfLong$Value)
  
  pltDfLong$Epoch <- as.numeric(pltDfLong$Epoch)
  pltDfLong$Epoch <- factor(pltDfLong$Epoch, levels = min(pltDfLong$Epoch):max(pltDfLong$Epoch))
  pltDfLong <- pltDfLong[!pltDfLong$Value %in% c(0, "0"), ]
  
  for(l in 1:length(table(pltDfLong$SID))){
    tmpSid <- names(table(pltDfLong$SID))[l]
    outDir <- paste(plotPath, "/", sep="")
    ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
    if(is.null(fileAdd)){
      outFile <- paste(outDir, "Epoch_Predstats_Plt_", tmpSid, ".pdf", sep="")
    }else{
      outFile <- paste(outDir, "Epoch_Predstats_Plt_", fileAdd,"_", tmpSid, ".pdf", sep="")
    }
    
    if(!is.null(multParamBool)){
      tmpPlt <- ggplot(pltDfLong, aes(fill=Measurement, y=Value, x=Epoch)) + 
        geom_boxplot() +
        #geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
        xlab("Histotype") +
        ylab("Percent") + 
        ggtitle(paste("Model performance for cohort",tmpSid,"\nEC excluded due to lack of HSPs" , sep="")) + 
        scale_fill_viridis(discrete = TRUE) +
        theme(text = element_text(size=24), 
              axis.text.x = element_blank(),
              legend.text=element_text(size=16),
              plot.title = element_text(hjust = 0.5)) +
        facet_grid(rows = vars(SID), cols = vars(Measurement), 
                   scales="free_y", switch = 'y')
    }else{
      tmpPlt <- ggplot(pltDfLong, aes(x=Epoch, y=Value, fill=Histotype)) + 
        geom_boxplot() +
        xlab("Histotype") +
        ylab("Percent") + 
        ggtitle(paste("Model performance for cohort",tmpSid, "\nEC excluded due to lack of HSPs" , sep="")) + 
        scale_fill_viridis(discrete = TRUE) +
        theme(text = element_text(size=24), 
              axis.text.x = element_blank(),
              legend.text=element_text(size=16),
              plot.title = element_text(hjust = 0.5)) +
        facet_grid(rows = vars(Histotype), cols = vars(Measurement), 
                   scales="free_y", switch = 'y')
    }
    ggsave(outFile, plot=tmpPlt, width=60, height=40, units = "cm")
  }
}

makePredStatPltFlipped <- function(inpPredDfLst, pltCols=NULL, histCols = NULL, remHist = NULL, fileAdd = NULL, binBool = NULL, multParamBool = NULL, fullDfBool=NULL){
  if(is.null(pltCols)){
    pltCols <- c("SID","Histotype", 
                 "Sensitivity", "Specificity", 
                 "Accuracy","F1Score",  
                 "Low", "High",
                 "AUC", "PRAUC")
  }
  if(is.null(histCols)){
    histCols <- c("CCC", "EC", "HGSC", "MC")
  }
  if(!is.null(remHist)){
    histCols <- histCols[!histCols %in% remHist]
  }
  pltDf <- data.frame(matrix(nrow=0, ncol=length(pltCols)))
  colnames(pltDf) <- pltCols
  if(is.null(binBool)){
    colnames(pltDf)[grep("Low", colnames(pltDf))] <- "Acc_Lower"
    colnames(pltDf)[grep("High", colnames(pltDf))] <- "Acc_Higher"
  }else{
    colnames(pltDf)[grep("Low", colnames(pltDf))] <- "CI_Lower"
    colnames(pltDf)[grep("High", colnames(pltDf))] <- "CI_Higher"
  }
  grepVec <- paste(histCols, collapse="|")
  
  for(i in 1:length(inpPredDfLst)){
    inpPredDf <- inpPredDfLst[[i]]
    if(!"SID" %in% colnames(inpPredDf)){
      inpPredDf$SID <- names(inpPredDfLst)[i]
    }
    if(length(which(colnames(inpPredDf) %in% "AUPRC"))>0){
      colnames(inpPredDf)[which(colnames(inpPredDf) %in% "AUPRC")] <- "PRAUC"
    }
    tmpSid <- names(inpPredDfLst)[i]
    colIndVec <- c()
    for(z in 1:length(pltCols)){
      tmpCol <- pltCols[z]
      tmpInds <- grep(tmpCol, colnames(inpPredDf))
      colIndVec <- append(colIndVec, tmpInds)
    }
    colIndVec <- colIndVec[!duplicated(colIndVec)]
    for(j in 1:nrow(inpPredDf)){
      tmpRank <-inpPredDf[j,]
      tmpVals <- tmpRank[, colIndVec]
      hVals <- tmpVals[grep(grepVec, colnames(tmpVals))]
      nhVals <- tmpVals[!colnames(tmpVals) %in% colnames(hVals)]
      for(k in 1:length(histCols)){
        tmpH <- histCols[k]
        hCols <- hVals[grep(tmpH, colnames(hVals))]
        hCols <- unlist(hCols)
        nhVals <- unlist(nhVals)
        hCols <- append(hCols, nhVals)
        hCols <- append(hCols, tmpH)
        names(hCols)[length(hCols)] <- "Histotype"
        names(hCols) <- gsub(paste(tmpH, "_", sep=""), "", names(hCols))
        nRow <- hCols[match(pltCols, names(hCols))]
        pltDf[nrow(pltDf)+1,] <- nRow
      }
    }
  }
  if(!is.null(multParamBool)){
    if(is.null(binBool)){
      pltDf$Histotype <- "Model"
    }
  }
  skipCols <- c("Histotype", "SID")
  pltDfLong <-  pltDf %>% 
    pivot_longer(
      cols = pltCols[!pltCols %in% skipCols], 
      names_to = "Measurement",
      values_to = "Value"
    )
  if(length(grep("Lower", pltDfLong$Measurement))>0){
    if(is.null(binBool)){
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% "Acc_Lower")] <- "Acc_Range"
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% c("Acc_Higher", "Acc_Upper"))] <- "Acc_Range"
    }else{
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% "CI_Lower")] <- "CI_Range"
      pltDfLong$Measurement[which(pltDfLong$Measurement %in% c("CI_Higher", "CI_Upper"))] <- "CI_Range"
    }
  }
  pltDfLong <- data.frame(pltDfLong)
  pltDfLong$Value <- gsub(",", ".", pltDfLong$Value)
  pltDfLong$Value <- as.numeric(pltDfLong$Value)
  
  pltDfLong <- pltDfLong[!pltDfLong$Value %in% c(0, "0"),]
  
  if(!is.null(multParamBool)){
    tmpPlt <- ggplot(pltDfLong, aes(fill=Histotype, y=Histotype, x=Value)) + 
      geom_boxplot() +
      #geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
      xlab("Percent") +
      ylab("Histotype") + 
      ggtitle(paste("Mean values for 5-fold cross-validated iterations (n=100)\nEC excluded due to lack of HSPs" , sep="")) + 
      scale_fill_viridis(discrete = TRUE) +
      theme(text = element_text(size=24), 
            axis.text.x = element_text(size=16, face="bold"),
            legend.text=element_text(size=16),
            plot.title = element_text(hjust = 0.5)) +
      facet_grid(rows = vars(Measurement), cols = vars(SID),  
                 scales="free_x", switch = 'x')
  }else{
    pltDfLong$Histotype <- factor(pltDfLong$Histotype, levels=c("CCC", "EC", "HGSC", "MC"))
    tmpPlt <- ggplot(pltDfLong, aes(fill=Histotype, y=Histotype, x=Value)) + 
      geom_boxplot() +
      #geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
      xlab("Percent") +
      ylab("Histotype") +  
      ggtitle(paste("Mean values for 5-fold cross-validated iterations (n=100)\nEC excluded due to lack of HSPs" , sep="")) + 
      scale_fill_viridis(discrete = TRUE) +
      theme(text = element_text(size=24), 
            axis.text.x = element_text(size=16, face="bold"),
            legend.text=element_text(size=16),
            plot.title = element_text(hjust = 0.5)) +
      facet_grid(rows = vars(Measurement), cols = vars(SID), 
                 scales="free_x", switch = 'x')
  }
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(is.null(fileAdd)){
    outFile <- paste(outDir, "Pred_Stats_Plt.pdf", sep="")
  }else{
    outFile <- paste(outDir, "Pred_Stats_Plt_", fileAdd, ".pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt, width=40, height=20, units = "cm")
}

################################################################################
# Deprecated
################################################################################

# makePredStatPlt_OLD <- function(inpPredDfLst, fileAdd = NULL, binBool = NULL){
#   # Script for plotting predictive classification statistics
#   pltDf <- data.frame(matrix(nrow=0, ncol=5))
#   colnames(pltDf) <- c("Histotype", "Sensitivity", "Specificity", "AUC","SID")
#   for(i in 1:length(inpPredDfLst)){
#     inpPredDf <- inpPredDfLst[[i]]
#     tmpSid <- names(inpPredDfLst)[i]
#     sensCols <- grep("Sensitivity",colnames(inpPredDf))
#     specCols <- grep("Specificity",colnames(inpPredDf))
#     aucCols <- grep("AUC",colnames(inpPredDf))
#     for(j in 1:nrow(inpPredDf)){
#       tmpRank <-inpPredDf[j,]
#       sensVals <- tmpRank[, sensCols]
#       specVals <- tmpRank[, specCols]
#       colnames(sensVals) <- gsub("_Sensitivity", "", colnames(sensVals) )
#       colnames(specVals) <- gsub("_Specificity", "", colnames(specVals) )
#       for(k in 1:length(sensVals)){
#         tmpH <- names(sensVals)[k]
#         tmpSens <- sensVals[[tmpH]]
#         tmpSpec <- specVals[[tmpH]]
#         if(!is.null(aucCols)){
#           tmpAUC <- tmpRank[, aucCols]
#           if(length(tmpAUC) > 1){
#             colnames(tmpAUC) <- gsub("_AUC", "", colnames(tmpAUC))
#             tmpAUC <- tmpAUC[[tmpH]]
#           }
#         }else{
#           tmpRank[, aucCols] <- gsub(",",".", tmpRank[, aucCols])
#           tmpRank[, aucCols] <- as.numeric(tmpRank[, aucCols])
#           tmpAUC <- tmpRank[, aucCols]
#         }
#         pltDf[nrow(pltDf)+1,] <- c(tmpH, tmpSens, tmpSpec, tmpAUC, tmpSid)
#       }
#     }
#   }
#   pltDf <- pltDf[!pltDf$Histotype %in% "EC", ]
#   pltDfLong <-  pltDf %>% 
#     pivot_longer(
#       cols = c("Sensitivity", "Specificity", "AUC"), 
#       names_to = "Measurement",
#       values_to = "Value"
#     )
#   pltDfLong <- data.frame(pltDfLong)
#   pltDfLong$Value <- gsub(",", ".", pltDfLong$Value)
#   pltDfLong$Value <- as.numeric(pltDfLong$Value)
#   pltDfLong$Histotype <- factor(pltDfLong$Histotype, levels=c("CCC", "EC", "HGSC", "MC"))
#   
#   tmpPlt <- ggplot(pltDfLong, aes(fill=Histotype, y=Value, x=Histotype)) + 
#     geom_boxplot() +
#     #geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
#     xlab("Histotype") +
#     ylab("Percent") + 
#     ggtitle(paste("Sensitivity and Specificity for multiclass models\nEC excluded due to lack of HSPs" , sep="")) + 
#     scale_fill_viridis(discrete = TRUE) +
#     theme(text = element_text(size=24), 
#           axis.text.x = element_text(size=16, face="bold"),
#           legend.text=element_text(size=16),
#           plot.title = element_text(hjust = 0.5)) +
#     facet_grid(rows = vars(SID), cols = vars(Measurement), 
#                scales="free_y", switch = 'y')
#   outDir <- paste(plotPath, "/", sep="")
#   ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
#   if(is.null(fileAdd)){
#     outFile <- paste(outDir, "StepWisePred_Stats_Plt.pdf", sep="")
#   }else{
#     outFile <- paste(outDir, "StepWisePred_Stats_Plt_", fileAdd, ".pdf", sep="")
#   }
#   ggsave(outFile, plot=tmpPlt, width=60, height=40, units = "cm")
# }