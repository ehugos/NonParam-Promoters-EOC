################################################################################
# DMP scripts
################################################################################

makeChampDmpRes <- function(inpBeta, inpPheno, focusCol=NULL, inpArr=NULL, pCut=NULL, bCut = NULL, focH=NULL){
  # Function for generating DMP-results using ChAMP
  message("Initiating ChAMP DMP-analysis")
  if(is.null(inpArr)){
    inpArr <- "EPIC"
  }
  if(is.null(focusCol)){
    focusCol <- "Histotype"
  }
  if(is.null(pCut)){
    pCut <- 0.05
  }
  if(is.null(bCut)){
    bCut <- 0.2
  }
  dmpChamp <- list() 
  focusLoc <- which(colnames(inpPheno)==focusCol)
  phenoCol <- inpPheno[,focusLoc]
  contVec <- makeCatCombs(phenoCol)
  if(!is.null(focH)){
    contNames <- sapply(contVec, function(x) strsplit(x, " - ")[[1]][[1]], USE.NAMES=FALSE)
    contVec <- contVec[which(contNames %in% focH)]
  }
  # Perform dmp-analysis for contrasts
  for(i in 1:length(contVec)){
    # We use try, as some of the comparisons will fail due to insufficient data for analysis
    # Convert DMRcate contrasts into single entries for champ
    tmpCont <- contVec[[i]]
    splitCont <- strsplit(contVec[i],split=' - ')
    splitCont <- strsplit(splitCont[[1]]," ")
    cont_1 <- splitCont[[1]]
    cont_2 <- splitCont[[2]]
    #compTitle <- paste("DMP_",  cont_1, "_", cont_2, sep="")
    # tryCatch to avoid script stopping due to error
    skip_to_next <- NULL
    tmpCont <- gsub(" - ", "_", tmpCont)
    message(paste("Running DMP analaysis for contrast: ", tmpCont, sep=""))
    tryCatch(
      dmpRes <-  ChAMP::champ.DMP(beta = inpBeta, 
                                  pheno =  phenoCol,
                                  adjPVal = pCut, 
                                  adjust.method = "BH",
                                  compare.group = c(cont_1, cont_2), 
                                  arraytype = inpArr)
      , error=function(e){
        skip_to_next <<- TRUE
        print(paste("Error: ", e, sep=""))
        message(paste("Error encountered for contrast: ", tmpCont, sep=""))
      })
    if(!is.null(skip_to_next)){
      dfCols <- c("logFC", "AveExpr", "t", "P.Value",   "adj.P.Val",  "B",   
                  paste(cont_1, "AVG", sep="_"),  paste(cont_2, "AVG", sep="_"), "deltaBeta", "CHR", "MAPINFO",
                  "Strand", "Type", "gene feature",  "cgi","feat.cgi", "UCSC_Islands_Name", "SNP_ID SNP_DISTANCE")
      falseDf <- data.frame(matrix(nrow=0, ncol=length(dfCols)))
      colnames(falseDf) <- dfCols
      dmpChamp[[i]] <- falseDf
      names(dmpChamp)[i] <- tmpCont
      next()
    }else{
      dmpRes <- dmpRes[[1]]
      dmpRes <- dmpRes[which(abs(dmpRes[,9]) > bCut),]
      dmpRes <- dmpRes[which(dmpRes[,5] < pCut),]
      dmpChamp[[i]] <- dmpRes
      names(dmpChamp)[i] <- tmpCont
    } 
  }
  return(dmpChamp)
}

makeMinfiDmpRes <- function(inpBeta, inpPheno, focusCol=NULL, inpArr=NULL, pCut=NULL, focH=NULL){
  # Script for findings minfi DMPs using DMPfinder
  message("Initiating Minfi DMP-analysis")
  if(is.null(inpArr)){
    inpArr <- "EPIC"
  }
  if(is.null(focusCol)){
    focusCol <- "Histotype"
  }
  if(is.null(pCut)){
    pCut <- 0.05
  }
  dmpChamp <- list() 
  focusLoc <- which(colnames(inpPheno)==focusCol)
  phenoCol <- inpPheno[,focusLoc]
  contVec <- makeCatCombs(phenoCol)
  if(!is.null(focH)){
    contNames <- sapply(contVec, function(x) strsplit(x, " - ")[[1]][[2]], USE.NAMES=FALSE)
    contVec <- contVec[which(contNames %in% focH)]
  }
  # Perform dmp-analysis for contrasts
  for(i in 1:length(contVec)){
    # We use try, as some of the comparisons will fail due to insufficient data for analysis
    # Convert DMRcate contrasts into single entries for champ
    tmpCont <- contVec[[i]]
    splitCont <- strsplit(contVec[i],split=' - ')
    splitCont <- strsplit(splitCont[[1]]," ")
    cont_1 <- splitCont[[1]]
    cont_2 <- splitCont[[2]]
    
    minfiP <- inpPheno[which(inpPheno$Histotype %in% c(cont_1, cont_2)),]
    minfiB <- inpBeta[, match(minfiP$barcode, colnames(inpBeta))]
    
    #compTitle <- paste("DMP_",  cont_1, "_", cont_2, sep="")
    # tryCatch to avoid script stopping due to error
    skip_to_next <- NULL
    tmpCont <- gsub(" - ", "_", tmpCont)
    message(paste("Running DMP analaysis for contrast: ", tmpCont, sep=""))
    tryCatch(
      dmpRes <-  minfi::dmpFinder(dat= as.matrix(minfiB), 
                                  pheno =  factor(minfiP$Histotype, levels = c(cont_1, cont_2)),
                                  type = "categorical",
                                  qCutoff = pCut)
      , error=function(e){
        skip_to_next <<- TRUE
        print(paste("Error: ", e, sep=""))
        message(paste("Error encountered for contrast: ", tmpCont, sep=""))
      })
    if(!is.null(skip_to_next)){
      dmpChamp[[i]] <- dmpRes
      names(dmpChamp)[i] <- tmpCont
      next()
    }else{
      # Keep only DMP-results with delta-beta above or below abs(0.2)
      # dmpRes <- dmpRes[which(abs(dmpRes$intercept) > 0.2),]
      dmpRes <- dmpRes[which(dmpRes$qval < pCut),]
      dmpChamp[[i]] <- dmpRes
      names(dmpChamp)[i] <- tmpCont
    } 
  }
  return(dmpChamp)
}

makeDmpFreqDf <- function(inpDmp, inpPheno, excLst=NULL, freqBool=NULL){
  # Updated script for finding HSGs, creates a summary dataframe complete with 
  # frequency of the DEG and its l2fc if available for a contrast
  tmpHisto <- names(table(inpPheno$Histotype))
  if(!is.null(excLst)){
    tmpHisto <- tmpHisto[!tmpHisto %in% excLst]
  }
  # inpDmp <- inpDmp[!is.na(names(inpDmp))]
  outDmpLst <- list()
  for(i in 1:length(tmpHisto)){
    tmpH <- tmpHisto[i]
    comps <-  getFocusedComps(inpDmp, tmpH, revBool = TRUE)
    # Retrieve all entries that contain/match the histotype comparison
    tmpDmpLst <- inpDmp[names(inpDmp) %in% comps]
    for(j in 1:length(tmpDmpLst)){
      tmpDmpRes <- tmpDmpLst[[j]]
      tmpDmpRes <- tibble::rownames_to_column(tmpDmpRes, "cpgID")
      tmpDmpLst[[j]] <- tmpDmpRes
    }
    # Combine matrix into one combined matrix
    tmpDmpMat <- do.call(rbind, tmpDmpLst)
    unqDmp <- unique(tmpDmpMat$cpgID)
    hDmpDf <- data.frame(matrix(nrow=0, 
                                ncol=3+length(tmpDmpLst)))
    # Create DMP frequency dataframe
    colnames(hDmpDf) <- c("cpg", "noDmp", "dmpFreq", comps)
    for(j in 1:length(unqDmp)){
      tmpCpg <- unqDmp[[j]]
      dmpRows <- tmpDmpMat[tmpDmpMat$cpgID %in% tmpCpg, ]
      dmpRows <- tibble::rownames_to_column(dmpRows, "Contrast")
      dmpRows$Contrast <- gsub("\\.[0-9]{1,4}","",dmpRows$Contrast)
      dmpFreq <- round(nrow(dmpRows)/length(comps),3)
      rowInd <- nrow(hDmpDf)+1
      hDmpDf[rowInd, ] <- c(tmpCpg, nrow(dmpRows), dmpFreq, rep(NA, length(comps)))
      for(k in 1:nrow(dmpRows)){
        if(dmpRows$Contrast[k] %in% colnames(hDmpDf)){
          tmpColInd <- which(colnames(hDmpDf) %in% dmpRows$Contrast[k])
          hDmpDf[rowInd , tmpColInd] <- dmpRows$deltaBeta[k]
        }
      }
    }
    hDmpDf <- hDmpDf[order(hDmpDf$dmpFreq, decreasing = TRUE), ]
    if(!is.null(freqBool)){
      hDmpDf <- hDmpDf[hDmpDf$dmpFreq >= freqBool,]
    }
    outDmpLst[[i]] <- hDmpDf
    names(outDmpLst)[i] <- tmpH
  }
  return(outDmpLst)
}

# Function for creating contrast vectors for use with DMRcate
makeCatCombs <- function(phenoInp){
  catLst <- unique(phenoInp)
  cat_combs <- list() 
  catCross <- crossing(phenoInp, phenoInp)
  for(i in 1:nrow(catCross)){
    if(catCross[i,1] == catCross[i,2]){
      next
    }else{
      cat_combs[i] <- paste(catCross[i,1], catCross[i,2], sep= " - ")  
    }
  }
  cat_combs <- unlist(cat_combs)
  cat_combs[cat_combs == "NULL"] = NA
  cat_combs <- na.omit(cat_combs)
  cat_combs <- unique(cat_combs)
  return(cat_combs)
}

makeDMPMatch <- function(inpMeth, cpgLocs, inpGenes){
  inpGenes <- na.omit(unique(inpGenes))
  outLst <- list()
  # Need to adjust script to factor in conserved DMRs (i.e. to extract DMRs matching DEgs)
  for(i in 1:length(histotypes)){
    tmpH <- histotypes[i]
    tmpConts <- sapply(names(inpExp), function(x) strsplit(x, "_")[[1]][1], USE.NAMES=FALSE)
    contInds <- which(tmpConts %in% tmpH)
    expLst <- inpExp[contInds]
    geneDf <- data.frame(matrix(nrow=length(inpGenes), ncol = 4))
    rownames(geneDf) <- inpGenes
    colnames(geneDf) <- c("Gene", names(expLst))
    geneDf$Gene <- rownames(geneDf)
    for(j in 1:length(expLst)){
      tmpExp <- expLst[[j]]
      tmpName <- names(expLst)[j]
      matchInds <- which(geneDf$Gene %in% tmpExp$external_gene_name)
      geneDf[matchInds, tmpName] <- 1
    }
    geneDf <- geneDf[rowSums(is.na(geneDf[, names(expLst)])) < 3,]
    geneDf <- geneDf[order(rowSums(geneDf[, names(expLst)])),]
    outLst[[i]] <- geneDf
    names(outLst)[i] <-tmpH
  }
  return(outLst)
}

################################################################################
# DMR analysis
################################################################################

# Function for creating DMR-annotation for use with DMRcate
makeDmrAnno <- function(inpBeta, phenoInp, focusCol, pCut = NULL){
  if(is.null(pCut)){
    pCut <- 0.05
  }
  # Get index position of focus column in phenotype input matrix
  focusLoc <- which(colnames(phenoInp) %in% focusCol)
  # Make contrast vector for DMRcate to use in comparison
  conts <- makeCatCombs(phenoInp[,focusLoc])
  # Import diagnosis & histotype as factor-vectors for comparative analysis
  hisFactor <- factor(phenoInp[,focusLoc])
  hisFactor <- relevel(hisFactor, ref="HGSC")
  # Design matrix for Focus group vs. Control without intercept
  tmpDesMat <- model.matrix(~0+hisFactor, data=phenoInp)
  colnames(tmpDesMat) <- gsub("hisFactor","",colnames(tmpDesMat))
  # Create contrast matrixes for DMRcate
  tmpContMatrix <- limma::makeContrasts(contrasts = conts,levels=tmpDesMat)
  # Iterate over all given contrasts in the matrix and extract their individual results
  # CPG annotate allows for a contrast matrix but only for one contrast at a given time
  annoLst <- list()
  # Loop through annotation object, store results in resLst
  for (i in 1:length(conts)){
    contName <- conts[i]
    message(paste("Now annotating: ", contName, sep=""))
    resName <- DMRcate::cpg.annotate("array",
                                     inpBeta,
                                     what="Beta",
                                     arraytype = "EPICv1",
                                     analysis.type="differential", 
                                     contrast=TRUE, 
                                     cont.matrix=tmpContMatrix, 
                                     design=tmpDesMat,
                                     coef=contName,
                                     fdr=pCut,
                                     annotation= c(array = "IlluminaHumanMethylationEPIC", 
                                                   annotation = "ilm10b5.hg38"))
    annoLst[[i]] <- resName
    names(annoLst)[i] <- contName
  }
  return(annoLst)
}

makeDmrResult <- function(annoLst, minCpg = NULL, betaCut = NULL, lambdaInp = NULL, cInp = NULL){ 
  if(is.null(minCpg)){
    minCpg <- 3
  }
  if(is.null(betaCut)){
    betaCut <- 0.2
  }
  if(is.null(lambdaInp)){
    lambdaInp <- 1000
  }
  if(is.null(cInp)){
     cInp <- 2
  }
  # Makes DMR-results from annotated meth-object by DMRcate, saves into named list
  dmrResLst <- list()
  for (i in 1:length(annoLst)){
    skip_to_next <- FALSE
    tmpAnn <- annoLst[[i]]
    tmpName <- names(annoLst)[i]
    message(paste("Generating dmr-results for contrast: ", tmpName, sep=""))
    # Trycatch clause for entries with no results
    tryCatch(
      # Absolute delta-beta > 0.2 (20%) for significant regions
      # Minimum number of consecutive cpg's for DMR
      tmpDMR <- DMRcate::dmrcate(tmpAnn, 
                                 lambda=lambdaInp, 
                                 C=cInp, 
                                 betacutoff=betaCut,
                                 min.cpg=minCpg)
      , error=function(e){
        print(paste("Error: ", e, sep=""))
        skip_to_next <<- TRUE
      })
    if(skip_to_next){
      message(paste("Contrast: ", tmpName, " skipped due to error", sep=""))
      next
    }else{
      message(paste("Contrast : ", tmpName, " demarcated", sep=""))
      dmrResLst[[i]] <- tmpDMR
      names(dmrResLst)[i] <- tmpName
    }
  }
  return(dmrResLst)
}

# Extract ranges from DMRcate results
makeDmrRanges <- function(dmrResLst){
  # Requires DMR-results objects in the form of named list as input
  rangeLst <- list() 
  for(i in 1:length(dmrResLst)){
    skip_to_next <- FALSE
    tmpDmr <- dmrResLst[[i]]
    tmpName <- names(dmrResLst)[i]
    message(paste("Extracting ranges for contrast: ", tmpName, sep=""))
    tryCatch(
      tmpRanges <- DMRcate::extractRanges(tmpDmr, genome = "hg19")
      , error=function(e){
        print(paste("Error: ", e, sep=""))
        skip_to_next <- TRUE
      })
    if(skip_to_next){
      message(paste("Ranges could not be extracted due to error for: ", tmpName, sep=""))
      next
    }else{
      tmpRanges <- tmpRanges[order(-abs(tmpRanges$maxdiff)),]
      rangeLst[[i]] <- tmpRanges
      names(rangeLst)[i] <- tmpName 
      # Save ranges to local directory
      rangePath <- paste(outPath,"/" ,tmpName, '_Ranges.csv', sep="")
      write.csv(tmpRanges, file=rangePath)
    } 
  }
  return(rangeLst)
}

# Extract ranges from DMRcate results object(s) in the form of named list, saves into named list
makeFunAnno <- function(annoLst, cpgInp, collectionInp, annoInp){
  goResLst <- list()
  for(i in 1:length(annoLst)){
    skip_to_next <- FALSE
    # Requires named list of cpg.annotate result objects
    tmpRegion <- annoLst[[i]]
    tmpName <- names(annoLst)[i]
    message(paste("Functional annotation initiated for: ", tmpName, sep=""))
    tryCatch(
      # GO-annotation using KEGG or GO
      tmpGo <- missMethyl::goregion(
        tmpRegion,
        all.cpg = cpgInp,
        collection = collectionInp,
        array.type = "EPIC",
        prior.prob = TRUE,
        anno = annoInp,
        fract.counts = TRUE,
        genomic.features = c("TSS200","TSS1500"),
        sig.genes = TRUE)
      , error=function(e){
        print(paste("Error: ", e, sep=""))
        skip_to_next <- TRUE
      })
    if(skip_to_next){
      message(paste("Functional annotation skipped for: ", tmpName, " due to error!", sep=""))
      next
    }else{
      goResLst[[i]] <- tmpGo
      names(goResLst)[i] <- tmpName
    }
    # Get top 20 results
    # topHisGO[[i]] <- topGSA(tempGO, number = 20, sort = TRUE)
  }
  return(goResLst)
}

# Get CpG IDs for each DMR
makeDmrCpgId <- function(inpRanges, cpgRanges){
  outLst <- list()
  for(i in 1:length(inpRanges)){
    contName <- names(inpRanges)[i]
    tmpDmr <- inpRanges[[i]]
    tmpDmr$cpgIds <- NA
    message(paste("Retrieving CpG's for DMR's for contrast: ", contName, sep=""))
    # Iterate over each DMR, retrieve overlapping CpG's, save ID's as list in dataframe
    for(j in 1:length(tmpDmr)){
      skip_to_next <- FALSE
      tryCatch(
        tmpDmrCpg <- names(subsetByOverlaps(cpgRanges, tmpDmr[j]))
        , error=function(e){
          print(paste("Error: ", e, sep=""))
          skip_to_next <- TRUE
        })
      if(skip_to_next){
        next
      }else{
        tmpCpgLst <- list(tmpDmrCpg)
        tmpDmr$cpgIds[j] <- tmpCpgLst 
      }
    }
    # Must be in character format if we wish to write to .csv format
    tmpDmr$cpgIds = as.character(tmpDmr$cpgIds)
    outLst[[i]] <- tmpDmr
    names(outLst)[i] <- contName
  }
  return(outLst)
}

makeDmpRegCovDf <- function(inpCpgPos, inpDmpLst){
  hspDmpDf <- data.frame(matrix(nrow=0, ncol=8))
  colnames(hspDmpDf) <- c("Histotype", "Ref_Sig_CpG", 
                          "Cov_1", "Freq_1",
                          "Cov_2", "Freq_2",
                          "Cov_3", "Freq_3")
  # Simple summary script of frequency of DMP results for a phenotype when compared to a reference
  for(i in 1:length(inpCpgPos)){
    tmpCpgs <- inpCpgPos[[i]]
    tmpH <- names(inpCpgPos)[i]
    if(is.null(inpCpgPos[[i]])){
      next()
    }
    contNames <- sapply(names(inpDmpLst), function(x) strsplit(x, "_")[[1]][[2]], USE.NAMES=FALSE)
    hConts <- inpDmpLst[which(contNames %in% tmpH)]
    matchVec <- c()
    for(j in 1:length(hConts)){
      tmpCont <- hConts[[j]]
      matchDmp <- rownames(tmpCont)[which(rownames(tmpCont) %in% tmpCpgs)]
      matchVec <- append(matchVec, matchDmp)
    }
    tabCov <- table(matchVec)
    cov_1 <- length(which(tabCov >= 1))
    cov_2 <- length(which(tabCov >= 2))
    cov_3 <- length(which(tabCov >= 3))
    hspDmpDf[nrow(hspDmpDf)+1, ] <- c(tmpH, length(tmpCpgs), 
                          cov_1, round(cov_1/length(tmpCpgs),2),
                          cov_2, round(cov_2/length(tmpCpgs),2),
                          cov_3, round(cov_3/length(tmpCpgs),2))
  }
  return(hspDmpDf)
}

makeMissDmpType <- function(inpBeta, inpPheno, inpCont, missCpg, revBool = NULL){
  # input beta dataframe, 
  # input-pheno with histotype column and barcode column
  # input contrast (groups present in histotype columns)
  # missing CpG sites as vector
  refHist <- inpCont[[1]]
  contHist <- inpCont[[2]]
  refSamps <- inpPheno$barcode[which(inpPheno$Histotype %in% refHist)]
  contSamps <- inpPheno$barcode[which(inpPheno$Histotype %in% contHist)]
  missB <- inpBeta[missCpg,]
  compType <- data.frame(matrix(nrow=1, ncol=6))
  compType[1,] <- c(0,0,0,0,0,0)
  colnames(compType) <- c("Ref","Comp", "nMissing", 
                          "Normal_Normal", "Normal_NotNormal", "NotNormal_NotNormal")
  if(is.null(missCpg) | length(missCpg) == 0){
    message(message("ERROR: No missing CpGs supplemented for contrast: ", refHist, "_", contHist, sep=""))
    return(compType)
  }
  compType[1, "Ref"] <- refHist
  compType[1, "Comp"] <- contHist
  # Revbool allows us to also check no. missing or found DMPs 
  nIt <- length(missCpg)
  compType[1, "nMissing"] <- nIt
  message(paste("Investigating distribution type of missing CpGs for reference group: ", refHist, 
                " and contrast group: ", contHist, sep=""))
  for(k in 1:nIt){
    if(abs(k)%%round(nIt/10) == 0 || is.na(abs(k)%%round(nIt/10))){
      message(paste("Processing CpG: [", k, "/", nIt,"]"))
    }
    tmpC <- missCpg[k]
    refVec <- missB[tmpC, which(colnames(missB) %in% refSamps)]
    contVec <- missB[tmpC, which(colnames(missB) %in% contSamps)]
    # Normality test
    refNorm <-  shapiro.test(as.numeric(refVec))
    contNorm <- shapiro.test(as.numeric(contVec))
    if(refNorm$p.value >= 0.05 & contNorm$p.value >= 0.05){
        compType[1, "Normal_Normal"] <- compType[1, "Normal_Normal"] + 1
    }else if(refNorm$p.value >= 0.05 & contNorm$p.value <= 0.05 || contNorm$p.value >= 0.05 & refNorm$p.value <= 0.05){
        compType[1, "Normal_NotNormal"] <- compType[1, "Normal_NotNormal"] + 1
    }else{
        compType[1, "NotNormal_NotNormal"] <- compType[1, "NotNormal_NotNormal"] + 1
    }
  }
  return(compType)
}

makeDMPRegCovCompPlt <- function(inpCovLst, noHits = NULL, fileExt=NULL){
  # Function for plotting DMP coverage
  
  covDf <- do.call(rbind.data.frame, inpCovLst)
  covDf <- rownames_to_column(covDf, var="Program")
  covDf$Program <- gsub("\\.[1-9]", "", covDf$Program)
  percCols <- grep("Freq", colnames(covDf))
  pltDf <-covDf %>% 
    pivot_longer(
      cols = percCols, 
      names_to = "Coverage_perc",
      values_to = "Percent"
    )
  pltDf <- data.frame(pltDf)
  pltDf$Coverage_perc <- gsub("Freq_", "", pltDf$Coverage_perc)
  pltDf$Percent <- as.numeric(pltDf$Percent)
  pltDf$Percent <- 100*pltDf$Percent
  pltDf$Histotype <- factor(pltDf$Histotype, levels = names(table(pltDf$Histotype)))
  
  colCats <- unique(pltDf$Program)
  colvals <- viridis(n=length(colCats), option="B")
  colVec <- colvals
  names(colVec) <- colCats
  
  tmpPlt <- ggplot(pltDf, aes(fill=Program, y=Percent, x=Coverage_perc)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    geom_text(aes(label = Percent), 
              position = position_dodge(width = 1), 
              size = 6, 
              vjust = 1.5, 
              colour = "grey") +
    xlab("EOC comparisons (n=3)") +
    ylab("Percent of HSP CpG's") + 
    ggtitle(paste("Percentage of significant CpGs in \nnp-promoter regions found as DMPs" , sep="")) + 
    scale_fill_manual(values = colVec) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    guides(fill=guide_legend(title="Program")) +
    facet_wrap(~factor(Histotype), nrow=3, strip.position = "bottom")
  
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(!is.null(fileExt)){
    outFile <- paste(outDir, "extDMPCov_all_", fileExt,".pdf", sep="")
  }else{
    outFile <- paste(outDir, "extDMPCov_all.pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
}

makeMissDMPPlt <- function(missLst, refCpgLst, noHits = NULL, fileAdd = NULL, titleCat = NULL){
  if(is.null(titleCat)){
    titleCat <- "Promoter"
  }
  # Combine missing promoters in contrasts as one dataframe
  missLst <- missLst[!names(missLst) %in% ""]
  mergeMissDf <- Reduce(function(...) merge(..., all=T), missLst)
  
  missCCC <- as.numeric(mergeMissDf$nMiss[mergeMissDf$Ref %in% "CCC"])
  missHGSC <- as.numeric(mergeMissDf$nMiss[mergeMissDf$Ref %in% "HGSC"])
  missMC <- as.numeric(mergeMissDf$nMiss[mergeMissDf$Ref %in% "MC"])
  
  totCCC <- length(refCpgLst[["CCC"]])
  totHGSC <- length(refCpgLst[["HGSC"]])
  totMC <- length(refCpgLst[["MC"]])
  
  missCCCR <- paste(min(missCCC), max(missCCC), sep="-")
  missHGSCR <- paste(min(missHGSC), max(missHGSC), sep="-")
  missMCR <- paste(min(missMC), max(missMC), sep="-")
  totVec <- c(rep(totCCC[1], length(missCCC)), 
               rep(totHGSC[1], length(missHGSC)),
               rep(totMC[1], length(missMC)))
  mergeMissDf$nRef <- totVec
  missVec <- c(rep(missCCCR[1], length(missCCC)), 
               rep(missHGSCR[1], length(missHGSC)),
               rep(missMCR[1], length(missMC)))
  mergeMissDf$nMiss <- missVec
  
  pltDf <- mergeMissDf %>% pivot_longer(
    cols = c(4:6), 
    names_to = "DistComp",
    values_to = "Hits")
  pltDf <- data.frame(pltDf)
  pltDf$Comp <- factor(pltDf$Comp)
  
  colCats <- names(table(pltDf$DistComp))
  colvals <- viridis(n=length(colCats), option="C")
  colVec <- colvals
  names(colVec) <- colCats
  
  tmpPlt <- ggplot(pltDf, aes(fill=DistComp, y=Hits, x=Comp)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    # geom_bar(stat="identity", width=.5, position = "dodge") 
    #geom_text(aes(label = Percent), size = 8, vjust = 1.5, colour = "white") +
    xlab("EOC group compared with reference") +
    ylab("Number of CpGs") + 
    ggtitle(paste("Distribution type of significant", titleCat, " CpG's not\nfound as DMPs in any EOC comparison" , sep="")) + 
    scale_fill_manual(values = colVec) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    guides(fill=guide_legend(title="Distribution\nTypes")) +
    facet_wrap(~factor(Ref), 
               nrow=4, 
               strip.position = "bottom", scales = "free_y")
  
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(!is.null(fileAdd)){
    outFile <- paste(outDir,fileAdd, "_missingDMPPlt.pdf", sep="")
  }else{
    outFile <- paste(outDir, "missingDMPPlt.pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt, width=30, height=30, units = "cm")
}

makeHitTypeDf <- function(inpCov, noHits = NULL, fileExt = NULL){
  # Create summary dataframe for the percentage of coverage in no. hits for each dmr-caller
  wDf <- data.frame(matrix(nrow=0, ncol=3))
  wDf <- inpCov %>% pivot_longer(
    cols = c(6:ncol(inpCov)), 
    names_to = "Program",
    values_to = "Hits"
  )
  # pltDf$Hits <- ifelse(pltDf$Hits > 0, 1, pltDf$Hits)
  wDf$totCpgNew <- ifelse(wDf$totCpgNew >= 5, "5+", 
                          ifelse(wDf$totCpgNew == 4, 4, "3-"))
  wDf$totCpgNew <- as.factor(wDf$totCpgNew)
  wDf$Hits <- as.factor(wDf$Hits)
  
  hspCats <- nrow(wDf)/4
  wDf$Program <- ifelse(wDf$Program %in% "DMRCATE", "DMRcate", wDf$Program)
  wDf$Program <- ifelse(wDf$Program %in% "BUMPHUNTER", "Bumphunter",wDf$Program)
  wDf$Program <- ifelse(wDf$Program %in% "PROBELASSO", "ProbeLasso",wDf$Program)
  wDf$Program <- ifelse(wDf$Program %in% "SESAME", "Sesame",wDf$Program)
  wDf$Program <- factor(wDf$Program,  levels=c("Bumphunter","DMRcate", "ProbeLasso", "Sesame"))
  # Create dataframe of the total number of regions, and in how many comparisons these regions were found
  hitTypePltDf <- data.frame(matrix(nrow=0, ncol=9))
  colnames(hitTypePltDf) <- c("Program", "no_Regions", "no_Cov", 
                              "1_hits", "1_perc", 
                              "2_hits", "2_perc", 
                              "3_hits", "3_perc")
  for(i in 1:length(names(table(wDf$Program)))){
    tmpP <- names(table(wDf$Program))[i]
    tmpCats <- table(wDf$Hits[wDf$Program %in% tmpP])
    tmpCats <- tmpCats[!names(tmpCats) %in% 0]
    p1 <- sum(tmpCats)/hspCats
    p2 <- sum(tmpCats[which(names(tmpCats) > 1)])/hspCats
    p3 <- sum(tmpCats[which(names(tmpCats) > 2)])/hspCats
    hitTypePltDf[nrow(hitTypePltDf)+1, ] <- c(tmpP, hspCats, sum(tmpCats), 
                                              sum(tmpCats), p1, 
                                              sum(tmpCats[which(names(tmpCats) > 1)]), p2, 
                                              sum(tmpCats[which(names(tmpCats) > 2)]), p3)
  }
  
  histCats <- table(wDf$Histotype)/4
  histNames <- names(histCats)
  # Create dataframe showing how many of the significant regions for each histotype that were found
  hitTypePltDf2 <- data.frame(matrix(nrow=0, ncol=3+(2*length(histNames))))
  colnames(hitTypePltDf2) <- c("Program", "no_Regions",
                               paste(histNames, "Regions", sep="_"), 
                               "no_Cov", 
                               histNames)
  for(j in 1:length(names(table(wDf$Program)))){
    tmpP <- names(table(wDf$Program))[j]
    tmpCats <- table(wDf$Histotype[wDf$Program %in% tmpP & !wDf$Hits %in% 0])
    tmpCats <- tmpCats[!names(tmpCats) %in% 0]
    if(length(which(!histNames %in% names(tmpCats))) > 0){
      tmpCats <- append(tmpCats, rep(0, length(which(!histNames %in% names(tmpCats)))))
      names(tmpCats)[(length(which(histNames %in% names(tmpCats)))+1):length(tmpCats)] <-  histNames[which(!histNames %in% names(tmpCats))]
      tmpCats <- tmpCats[match(histNames, names(tmpCats))]
    }
    hitTypePltDf2[nrow(hitTypePltDf2)+1, ] <- c(tmpP, hspCats, 
                                                histCats, 
                                                sum(tmpCats), 
                                                tmpCats)
  }
  return(list("HitType"=hitTypePltDf, "HistType"= hitTypePltDf2))
}

makeDmpStatPlt <- function(inpDmpLst, fileExt = NULL){
  pltDf <- data.frame(matrix(ncol=3, nrow=0))
  colnames(pltDf) <- c("Program", "Contrast","nDMP")
  for(i in 1:length(inpDmpLst)){
    tmpProg <- names(inpDmpLst)[i]
    tmpDmp <- inpDmpLst[[i]]
    for(j in 1:length(tmpDmp)){
      tmpCont <- names(tmpDmp)[j]
      contDmp <- tmpDmp[[j]]
      nDmp <- nrow(contDmp)
      pltDf[nrow(pltDf)+1, ] <- c(tmpProg, tmpCont, nDmp)
    }
  }
  colCats <- unique(pltDf$Program)
  colvals <- viridis(n=length(colCats), option="B")
  colVec <- colvals
  names(colVec) <- colCats
  pltDf$nDMP <- as.numeric(pltDf$nDMP)
  
  tmpPlt <- ggplot(pltDf, aes(fill=Program, y=nDMP, x=Contrast)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    xlab("Contrast") +
    ylab("Number of DMPs in contrast") + 
    ggtitle("Number of DMP's found in contrasts in\nEOC histotype comparisons") + 
    scale_fill_manual(values = colVec) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(angle = 90, size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5))
  if(!is.null(fileExt)){
    outFile <- paste(plotPath,"/", fileExt, "_dmpContCounts.pdf", sep="")
  }else{
    outFile <- paste(plotPath, "dmpContCounts.pdf", sep="/")
  }
  ggsave(filename = outFile, plot=tmpPlt, width=30, height=20, units = "cm")
}

makeDmpDiffDf <- function(inpDmpLst, fileExt = NULL, pltBool = NULL){
  pltDf <- data.frame(matrix(ncol=9, nrow=0))
  colnames(pltDf) <- c("Program", "Contrast", "nDMP", "Comparison", 
                       "nCont", "OL", "Missing", 
                       "percOl", "percContOl")
  missDf <- data.frame(matrix(ncol=4, nrow=0))
  colnames(missDf) <- c("Program", "Contrast", "Comparison", "Missing")
  for(i in 1:length(inpDmpLst)){
    tmpProg <- names(inpDmpLst)[i]
    tmpDmp <- inpDmpLst[[i]]
    for(j in 1:length(tmpDmp)){
      tmpCont <- names(tmpDmp)[j]
      tmpRef <- tmpDmp[[tmpCont]]
      nDmp <- nrow(tmpRef)
      compProgs <- names(inpDmpLst)[-which(names(inpDmpLst) %in% tmpProg)]
      for(k in 1:length(compProgs)){
        compP <- compProgs[[k]]
        dmpComp <- inpDmpLst[[compP]][[tmpCont]]
        nCDmp <- length(rownames(dmpComp))
        nOl <- length(which(rownames(tmpRef) %in% rownames(dmpComp)))
        olPerc <- round((nOl/nDmp),3)
        contPerc <- round((nOl/nCDmp),3)
        nMiss <- length(which(!rownames(tmpRef) %in% rownames(dmpComp)))
        missVec <- rownames(tmpRef)[which(!rownames(tmpRef) %in% rownames(dmpComp))]
        missVec <- paste(missVec, collapse = ",")
        pltDf[nrow(pltDf)+1, ] <- c(tmpProg, tmpCont, nDmp, 
                                    compP, nCDmp,
                                    nOl, nMiss, 
                                    olPerc, contPerc)
        missDf[nrow(missDf)+1, ] <- c(tmpProg, tmpCont, compP, missVec)
      }
    }
  }
  pltDf$percOl <- 100*as.numeric(pltDf$percOl)
  if(!is.null(pltBool)){
    colCats <- unique(pltDf$Program)
    colvals <- viridis(n=length(colCats), option="B")
    colVec <- colvals
    names(colVec) <- colCats
    tmpPlt <- ggplot(pltDf, aes(fill=Comparison, y=percOl, x=Contrast)) + 
      geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
      xlab("Contrast") +
      ylab("Percent of DMPs") + 
      ggtitle("Percent of DMPs found in other DMP callers\nin EOC histotype comparisons") + 
      scale_fill_manual(values = colVec) +
      ylim(0, 100) +
      theme(text = element_text(size=24), 
            axis.text.x = element_text(angle = 90, size=16, face="bold"),
            legend.text=element_text(size=16),
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~factor(Program), nrow=4, strip.position = "bottom")
    if(!is.null(fileExt)){
      outFile <- paste(plotPath,"/", fileExt, "_dmpDiff.pdf", sep="")
    }else{
      outFile <- paste(plotPath, "dmpDiff.pdf", sep="/")
    }
    ggsave(filename = outFile, plot=tmpPlt, width=30, height=20, units = "cm")
  }
  return(list("diffDf" = pltDf, "missDf" = missDf))
}

makeDMROls <- function(inpDmrLst){
  # Compare DMR's based on genomic coordinates
  extDMROlPercDf <- data.frame(matrix(nrow=0, ncol=9))
  colnames(extDMROlPercDf) <- c("Reference","Comparison", "Contrast", 
                                "refDmrCount","compDmrCount", 
                                "overlaps", "unique", 
                                "percOl", "percUnq")
  for(i in 1:length(inpDmrLst)){
    refProg <- names(inpDmrLst)[i]
    refProgLst <- inpDmrLst[[refProg]]
    nProg <- names(inpDmrLst)[!names(inpDmrLst) %in% refProg]
    nProgLst <- inpDmrLst[nProg]
    for(j in 1:length(names(refProgLst))){
      tmpCont <- names(refProgLst)[j]
      for(k in 1:length(nProgLst)){
        compProg <- names(nProgLst)[k]
        refDmr <- refProgLst[[tmpCont]]
        contDmr <- nProgLst[[k]]
        if(!tmpCont %in% names(contDmr)){
          extDMROlPercDf[nrow(extDMROlPercDf)+1,] <- c(refProg, compProg, tmpCont, 
                                                       nrow(refDmr),0,
                                                       0, 0, 0, 0)
          next()
        }
        contDmr <- contDmr[[tmpCont]]
        if(is.null(refProg)){
          refParamVec <- c("chr"="chr","start"="start", "end"="end")
        }else if(toupper(refProg) %in% "DMRCATE"){
          refParamVec <- c("chr"="chr","start"="start", "end"="end")
        }else if(toupper(refProg) %in% "BUMPHUNTER"){
          refParamVec <- c("chr"="seqnames","start"="start", "end"="end")
        }else if(toupper(refProg) %in% "PROBELASSO"){
          refParamVec <- c("chr"="seqnames","start"="start", "end"="end")
        }else if(toupper(refProg) %in% "SESAME"){
          refParamVec <- c("chr"="Seg_Chrm","start"="Seg_Start","end"="Seg_End")
        }
        
        if(is.null(compProg)){
          compParamVec <- c("chr"="chr","start"="start", "end"="end")
        }else if(toupper(compProg) %in% "DMRCATE"){
          compParamVec <- c("chr"="chr","start"="start", "end"="end")
        }else if(toupper(compProg) %in% "BUMPHUNTER"){
          compParamVec <- c("chr"="seqnames","start"="start", "end"="end")
        }else if(toupper(compProg) %in% "PROBELASSO"){
          compParamVec <- c("chr"="seqnames","start"="start", "end"="end")
        }else if(toupper(compProg) %in% "SESAME"){
          compParamVec <- c("chr"="Seg_Chrm","start"="Seg_Start","end"="Seg_End")
        }
        
        refRanges <- makeGRangesFromDataFrame(refDmr,
                                              ignore.strand=TRUE,
                                              seqinfo=NULL,
                                              seqnames.field=refParamVec[["chr"]],
                                              start.field=refParamVec[["start"]],
                                              end.field=refParamVec[["end"]],
                                              starts.in.df.are.0based=FALSE)
        contRanges <- makeGRangesFromDataFrame(contDmr,
                                               ignore.strand=TRUE,
                                               seqinfo=NULL,
                                               seqnames.field=compParamVec[["chr"]],
                                               start.field=compParamVec[["start"]],
                                               end.field=compParamVec[["end"]],
                                               starts.in.df.are.0based=FALSE)
        dmrOLS <- findOverlaps(refRanges, 
                               contRanges)
        olsCount <- length(dmrOLS@from)
        unqCount <- length(unique(dmrOLS@from))
        percOLS <- 100*length(dmrOLS@from)/nrow(refDmr)
        unqPercOLS <- 100*length(unique(dmrOLS@from))/nrow(refDmr)
        extDMROlPercDf[nrow(extDMROlPercDf)+1,] <- c(refProg, compProg, tmpCont, 
                                                     nrow(refDmr),nrow(contDmr),
                                                     olsCount, unqCount, percOLS, unqPercOLS)
        
      }
    }
  }
  extDMROlPercDf$percUnq <- as.numeric(extDMROlPercDf$percUnq)
  return(extDMROlPercDf)
}

makeSesSegDmr_MULT <- function(inpSesDmrLst, inpCpgAnn, minCpg = NULL, filtDist = NULL){
  if(is.null(minCpg)){
    minCpg <- 1
  }
  # Sesame DMR's are currently only 1 CpG in length at its lowest
  # Adjust so that the smallest region contains 2 or more significant CpG's (i.e. smallest detected)
  # Group region into "DMR"s rather then showcase them as individual CpGs
  adjSesDMR <- list()
  for(i in 1:length(inpSesDmrLst)){
    tmpSes <- inpSesDmrLst[[i]]
    contName <- names(inpSesDmrLst)[i]
    if(nrow(tmpSes) == 0){
      message(paste("No DMRs found for: ", contName, sep=""))
      next()
    }
    contName <- names(inpSesDmrLst)[i]
    message(paste("Adjusting DMR results from sesame for contrast: ", contName,   " to comparable format", sep=""))
    if(!is.null(filtDist)){
      message(paste("Removing all CpGs with delta-beta < ", filtDist, sep=""))
      tmpSes <- tmpSes[which(abs(tmpSes$Seg_Est) >= filtDist), ]
    }
    allSegs <- table(tmpSes$Seg_ID)
    keepSegs <- names(allSegs[which(allSegs >= minCpg)])
    tmpSes <- tmpSes[which(tmpSes$Seg_ID %in% keepSegs), ]
    nIt <- length(keepSegs)
    sesDMRLst <- list()
    sesDMRLst <- future_lapply(1:nIt, function(j){
      if(abs(j)%%round(nIt/10) == 0){
        message(paste("Processing segment: [", j, "/", nIt, "] with n >=", minCpg, " CpGs", sep =""))
      }
      tmpSeg <- keepSegs[j]
      tmpSegs <- tmpSes[which(tmpSes$Seg_ID %in% tmpSeg), ]
      # Filter on Beta-value separation
      tmpSegs$start <- NA
      tmpSegs$end <- NA
      tmpSegs$strand <- NA
      tmpSegs$chr <- NA
      for(m in 1:nrow(tmpSegs)){
          tmpAnn <- inpCpgAnn[tmpSegs$Probe_ID[m],]
          tmpSegs$start[m] <- tmpAnn$start
          tmpSegs$end[m] <- tmpAnn$end
          tmpSegs$strand[m] <- tmpAnn$strand
          tmpSegs$chr[m] <- tmpAnn$chr
      }
      cpgStart <- tmpSegs$start[which.min(tmpSegs$start)]
      cpgEnd <- tmpSegs$end[which.max(tmpSegs$end)]
      cpgStrand <- tmpSegs$strand[1]
      tmpCpg <- tmpSegs$Probe_ID
      nCpg <- length(tmpCpg)
      estMean <- mean(tmpSegs$Seg_Est)
      sesDMRLst[[j]] <- c(tmpSegs$Seg_Chrm[1], tmpSegs$Seg_Start[1],tmpSegs$Seg_End[1], 
                            abs(as.numeric(tmpSegs$Seg_End[1]) - as.numeric(tmpSegs$Seg_Start[1])), estMean,
                            cpgStart, cpgEnd, cpgStrand, nCpg, paste(tmpCpg, collapse=","))
        # sesDMRDf[nrow(sesDMRDf)+1, ] <- c(tmpSegs$Seg_Chrm[1], tmpSegs$Seg_Start[1],tmpSegs$Seg_End[1], 
        #                                   abs(as.numeric(tmpSegs$Seg_End[1]) - as.numeric(tmpSegs$Seg_Start[1])),tmpSegs$Seg_Est[1],
        #                                   cpgStart, cpgEnd, cpgStrand, paste(tmpCpg, collapse=","))
    })
    sesDMRDf <- do.call(rbind.data.frame, sesDMRLst)
    colnames(sesDMRDf) <- c("Seg_Chrm","Seg_Start","Seg_End", "Width", "Mean_DeltaBeta",
                                                  "CpG_Start", "CpG_End","CpG_Strand", "nCpG","CpGs")
    adjSesDMR[[i]] <- sesDMRDf
    names(adjSesDMR)[i] <- contName
  }
  return(adjSesDMR)
}

makeDmrStatPlt <- function(inpDmrLst, fileExt = NULL){
  pltDf <- data.frame(matrix(ncol=length(inpDmrLst[[1]]), nrow=length(inpDmrLst)))
  rownames(pltDf) <- names(inpDmrLst)
  colnames(pltDf) <- names(inpDmrLst[[1]])
  for(i in 1:length(inpDmrLst)){
    tmpProg <- names(inpDmrLst)[i]
    tmpDmr <- inpDmrLst[[i]]
    for(j in 1:length(tmpDmr)){
      tmpCont <- names(tmpDmr)[j]
      contDmr <- tmpDmr[[j]]
      nDmr <- nrow(contDmr)
      pltDf[tmpProg, tmpCont] <- nDmr
    }
  }
  pltDf <- rownames_to_column(pltDf, "Program")
  pltDfLong <- pltDf %>% pivot_longer(
    cols = c(2:ncol(pltDf)), 
    names_to = "Contrast",
    values_to = "nDMR"
  )
  colCats <- unique(pltDfLong$Program)
  colvals <- viridis(n=length(colCats), option="B")
  colVec <- colvals
  names(colVec) <- colCats
  
  tmpPlt <- ggplot(pltDfLong, aes(fill=Program, y=nDMR, x=Contrast)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    xlab("Contrast") +
    ylab("Number of DMRs in contrast") + 
    ggtitle("Number of DMR's found in contrasts in\nEOC histotype comparisons") + 
    scale_fill_manual(values = colVec) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(angle = 90, size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5))
  if(!is.null(fileExt)){
    outFile <- paste(plotPath,"/", fileExt, "_dmrContCounts.pdf", sep="")
  }else{
    outFile <- paste(plotPath, "dmrContCounts.pdf", sep="/")
  }
  ggsave(filename = outFile, plot=tmpPlt, width=30, height=20, units = "cm")
}

# Get stats over HSP coverage in external DMR callers
makeDMRHSPstatDf <- function(inpDMRLst, inpHsp, inpHspLocs, inpBeta, inpSigs){
  outDf <- data.frame(matrix(nrow=0, ncol=8))
  colnames(outDf) <- c("Gene","totCpg", "sigCpg", "Histotype", "DMRCATE", "BUMPHUNTER", "PROBELASSO", "SESAME")
  for(i in 1:length(names(inpHsp))){
    refHist <- names(inpHsp)[i]
    tmpHsp <- inpHsp[[refHist]]
    tmpLocs <- inpHspLocs[[refHist]]
    tmpSigs <- inpSigs[[refHist]]
    for(j in 1:nrow(tmpHsp)){
      tmpGene <- tmpHsp[["external_gene_name"]][j]
      tmpEns <- tmpHsp[["ensembl_gene_id"]][j]
      geneLocs <- tmpLocs[[tmpEns]]
      # nCpG <- nrow(geneLocs)
      nCpG <- length(which(rownames(geneLocs) %in% rownames(inpBeta)))
      nSigCpg <- length(which(rownames(geneLocs) %in% tmpSigs))
      sigRegRanges <- makeGRangesFromDataFrame(geneLocs)
      outDf[nrow(outDf)+1, ] <- c(tmpGene, nCpG,nSigCpg,  refHist, 0,0,0,0)
      rownames(outDf)[nrow(outDf)] <- tmpGene
      for(k in 1:length(inpDMRLst)){
        tmpExt <- inpDMRLst[[k]]
        inpProg <- toupper(names(inpDMRLst)[k])
        contNames <- names(tmpExt)
        compNames <- sapply(contNames, function(x) strsplit(x, "_")[[1]][[2]], USE.NAMES=FALSE)
        contInds <- which(compNames %in% refHist)
        dmrConts <- tmpExt[contInds]
        
        matchInd <- 0
        for(l in 1:length(dmrConts)){
          tmpDmrs <- dmrConts[[l]]
          # Setup different columns for creating ranges based on input program
          if(toupper(inpProg) %in% "DMRCATE"){
            dmrParamVec <- c("chr"="seqnames","start"="start", "end"="end")
          }else if(toupper(inpProg) %in% "BUMPHUNTER"){
            dmrParamVec <- c("chr"="seqnames","start"="start", "end"="end")
          }else if(toupper(inpProg) %in% "PROBELASSO"){
            dmrParamVec <- c("chr"="seqnames","start"="start", "end"="end")
          }else if(toupper(inpProg) %in% "SESAME"){
            dmrParamVec <- c("chr"="Seg_Chrm","start"="Seg_Start","end"="Seg_End")
          }
          dmrRanges <- GenomicRanges::makeGRangesFromDataFrame(tmpDmrs,
                                                               ignore.strand=FALSE,
                                                               seqinfo=NULL,
                                                               seqnames.field=dmrParamVec["chr"],
                                                               start.field=dmrParamVec["start"],
                                                               end.field=dmrParamVec["end"],
                                                               starts.in.df.are.0based=FALSE)
          hspDmrOls <- GenomicRanges::findOverlaps(sigRegRanges, dmrRanges, ignore.strand=TRUE)
          if(length(hspDmrOls) >0){
            matchInd <- matchInd + 1
          }
        }
        outDf[tmpGene,inpProg] <- matchInd 
      }
    }
  }
  return(outDf)
}

makeDmrRegionstatDf <- function(inpBeta, inpDMRLst, inpRegions, inpRegionsPosLst, cpgCoords){
  # Get stats over coverage of defined regions in results from external DMR callers
  outDf <- data.frame(matrix(nrow=0, ncol=5+length(inpDMRLst)))
  colnames(outDf) <- c("Gene","totCpgOld", "sigCpgOld", "totCpgNew", "Histotype", 
                       toupper(names(inpDMRLst)))
  # Iterate through list of dataframes with CpGs
  for(i in 1:length(names(inpRegionsPosLst))){
    refHist <- names(inpRegions)[i]
    tmpRegs <- inpRegions[[refHist]]
    tmpLocs <- inpRegionsPosLst[[refHist]]
    for(j in 1:nrow(tmpRegs)){
      tmpReg <- tmpRegs[j, ]
      tmpEns <- tmpReg[, "ensembl_gene_id"]
      tmpGene <- tmpReg[, "external_gene_name"]
      if(tmpGene %in% ""){
        tmpGene <- tmpEns
      }
      # Get dataframe with CpG locations and id
      nSigCpg <- tmpReg[,"sigCpg"]
      if(!tmpEns %in% names(tmpLocs)){
        message(paste(tmpEns, " not found in inpRegionsPosLst", sep=""))
        next()
      }
      geneLocs <- tmpLocs[[tmpEns]]
      geneLocs <- geneLocs[which(rownames(geneLocs) %in% rownames(inpBeta)), ]
      nCpgOld <- nrow(geneLocs)
      # Match overlapping CpG names to the coordinate vector to get ranges
      newCoordLocs <- cpgCoords[rownames(geneLocs), ]
      # Check the number of available CpGs in the other given CpG annotation
      nCpGNew <- nrow(newCoordLocs)
      # Create genomic ranges from new coordinates
      sigRegRanges <- makeGRangesFromDataFrame(newCoordLocs)
      outDf[nrow(outDf)+1, ] <- c(tmpGene, nCpgOld, nSigCpg, nCpGNew,  refHist, 0,0,0,0)
      rownames(outDf)[nrow(outDf)] <- tmpGene
      # Loop through DMRs, compare their coordinates to that of the region of interest 
      for(k in 1:length(inpDMRLst)){
        tmpExt <- inpDMRLst[[k]]
        inpProg <- toupper(names(inpDMRLst)[k])
        contNames <- names(tmpExt)
        compNames <- sapply(contNames, function(x) strsplit(x, "_")[[1]][[2]], USE.NAMES=FALSE)
        contInds <- which(compNames %in% refHist)
        dmrConts <- tmpExt[contInds]
        
        matchInd <- 0
        for(l in 1:length(dmrConts)){
          tmpDmrs <- dmrConts[[l]]
          # Setup different columns for creating ranges based on input program
          if(toupper(inpProg) %in% "DMRCATE"){
            dmrParamVec <- c("chr"="chr","start"="start", "end"="end")
          }else if(toupper(inpProg) %in% "BUMPHUNTER"){
            dmrParamVec <- c("chr"="seqnames","start"="start", "end"="end")
          }else if(toupper(inpProg) %in% "PROBELASSO"){
            dmrParamVec <- c("chr"="seqnames","start"="start", "end"="end")
          }else if(toupper(inpProg) %in% "SESAME"){
            dmrParamVec <- c("chr"="Seg_Chrm","start"="Seg_Start","end"="Seg_End")
          }
          dmrRanges <- GenomicRanges::makeGRangesFromDataFrame(tmpDmrs,
                                                               ignore.strand=TRUE,
                                                               seqinfo=NULL,
                                                               seqnames.field=dmrParamVec["chr"],
                                                               start.field=dmrParamVec["start"],
                                                               end.field=dmrParamVec["end"],
                                                               starts.in.df.are.0based=FALSE)
          # Get overlaps between input region and external DMRs
          hspDmrOls <- GenomicRanges::findOverlaps(sigRegRanges, 
                                                   dmrRanges, 
                                                   ignore.strand=TRUE)
          if(length(hspDmrOls) >0){
            matchInd <- matchInd + 1
          }else{
            next()
          }
        }
        outDf[tmpGene,inpProg] <- matchInd 
      }
    }
  }
  return(outDf)
}

makeDMRRegion_CovCompPlt <- function(inpCov, noHits = NULL, fileExt = NULL){
  # Plot the coverage of HSPs in DMRs seen to no. comparisons they were found in
  wDf <- data.frame(matrix(nrow=0, ncol=3))
  wDf <- inpCov %>% pivot_longer(
    cols = c(6:ncol(inpCov)), 
    names_to = "Program",
    values_to = "Hits"
  )
  colnames(wDf)[which(colnames(wDf) %in% "totCpgNew")] <- "totCpg"
  wDf$totCpg <- ifelse(wDf$totCpg >= 5, "5+", ifelse(wDf$totCpg == 4, 4, "3-"))
  wDf$totCpg <- as.factor(wDf$totCpg)
  wDf$Hits <- as.factor(wDf$Hits)
  
  hspCats <- nrow(wDf)/4
  wDf$Program <- ifelse(wDf$Program %in% "DMRCATE", "DMRcate",wDf$Program)
  wDf$Program <- ifelse(wDf$Program %in% "BUMPHUNTER", "Bumphunter",wDf$Program)
  wDf$Program <- ifelse(wDf$Program %in% "PROBELASSO", "ProbeLasso",wDf$Program)
  wDf$Program <- ifelse(wDf$Program %in% "SESAME", "Sesame",wDf$Program)
  wDf$Program <- factor(wDf$Program,  levels=c("Bumphunter","DMRcate", "ProbeLasso", "Sesame"))
  
  hitTypePltDf <- data.frame(matrix(nrow=0, ncol=9))
  colnames(hitTypePltDf) <- c("Program", "no_HSP", "no_Cov", "1_hits", "1_perc", "2_hits", "2_perc", "3_hits", "3_perc")
  for(i in 1:length(names(table(wDf$Program)))){
    tmpP <- names(table(wDf$Program))[i]
    tmpCats <- table(wDf$Hits[wDf$Program %in% tmpP])
    tmpCats <- tmpCats[!names(tmpCats) %in% 0]
    p1 <- sum(tmpCats)/sum(hspCats)
    p2 <- sum(tmpCats[which(names(tmpCats) > 1)])/sum(hspCats)
    p3 <- sum(tmpCats[which(names(tmpCats) > 2)])/sum(hspCats)
    hitTypePltDf[nrow(hitTypePltDf)+1, ] <- c(tmpP, hspCats, sum(tmpCats), 
                                              sum(tmpCats), p1, 
                                              sum(tmpCats[which(names(tmpCats) > 1)]), p2, 
                                              sum(tmpCats[which(names(tmpCats) > 2)]), p3)
  }
  
  pltDf <- data.frame(matrix(nrow=0, ncol=3))
  colnames(pltDf) <- c("Program", "Percent", "Hits")
  
  colnames(hitTypePltDf) <- gsub("_perc", "",colnames(hitTypePltDf))
  inpCovLong <-  hitTypePltDf %>% 
    pivot_longer(
      cols = c("1", "2", "3"), 
      names_to = "Hits",
      values_to = "Percent"
    )
  inpCovLong$Percent <- as.numeric(inpCovLong$Percent)
  inpCovLong$Percent <- 100*inpCovLong$Percent
  for(k in 1:nrow(inpCovLong)){
    perc <- round(inpCovLong[k,"Percent"],2)
    pltDf[nrow(pltDf)+1,] <- c(inpCovLong[k,"Program"], perc, inpCovLong[k,"Hits"])
  }
  inpCovLong$Program <- factor(inpCovLong$Program,  levels=c("Bumphunter","DMRcate", "ProbeLasso", "Sesame"))
  
  pltDf$Hits <- ifelse(pltDf$Hits %in% 1, "1/3", ifelse(pltDf$Hits %in% 2, "2/3","3/3"))
  tmpPlt <- ggplot(pltDf, aes(fill=Program, y=Percent, x=Program)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    geom_text(aes(label = Percent), size = 8, vjust = -1, colour = "black") +
    xlab("DMR-caller") +
    ylab("Percent (of total HSPs)") + 
    ylim(0,100) + 
    ggtitle(paste("Percentage of regions n=(",sum(table(wDf$totCpg)/4), ") in DMR comparisons" , sep="")) + 
    scale_fill_manual(values = c("Bumphunter" = viridis(n=6, option="F")[2],
                                 "DMRcate"=viridis(n=6, option="F")[3],
                                 "ProbeLasso"=viridis(n=6, option="F")[4],
                                 "Sesame"=viridis(n=6, option="F")[5]),
    ) +
    #scale_fill_viridis(discrete = TRUE, option = "F") +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    guides(fill=guide_legend(title="DMR-caller")) +
    facet_wrap(~factor(Hits), nrow=3, strip.position = "bottom")
  
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(is.null(fileExt)){
    outFile <- paste(outDir, "extDMR_CpG_Cov.pdf", sep="")
  }else{
    outFile <- paste(outDir, fileExt, "_extDMR_CpG_Cov.pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
}

makeDMRRegion_CpGLenCompPlt <- function(inpCov, noHits = NULL, fileExt = NULL){
  wDf <- data.frame(matrix(nrow=0, ncol=3))
  wDf <- inpCov %>% pivot_longer(
    cols = c(6:ncol(inpCov)), 
    names_to = "Program",
    values_to = "Hits"
  )
  colnames(wDf)[which(colnames(wDf) %in% "totCpgNew")] <- "totCpg"
  # pltDf$Hits <- ifelse(pltDf$Hits > 0, 1, pltDf$Hits)
  wDf$totCpg <- ifelse(wDf$totCpg >= 5, "5+", ifelse(wDf$totCpg == 4, 4, "3-"))
  wDf$totCpg <- as.factor(wDf$totCpg)
  wDf$Hits <- as.factor(wDf$Hits)
  
  hspCats <- table(wDf$totCpg)/4
  wDf$Program <- ifelse(wDf$Program %in% "DMRCATE", "DMRcate",wDf$Program)
  wDf$Program <- ifelse(wDf$Program %in% "BUMPHUNTER", "Bumphunter",wDf$Program)
  wDf$Program <- ifelse(wDf$Program %in% "PROBELASSO", "ProbeLasso",wDf$Program)
  wDf$Program <- ifelse(wDf$Program %in% "SESAME", "Sesame",wDf$Program)
  wDf$Program <- factor(wDf$Program,  levels=c("Bumphunter","DMRcate", "ProbeLasso", "Sesame"))
  
  covTypePltDf <- data.frame(matrix(nrow=0, ncol=9))
  colnames(covTypePltDf) <- c("Program", "no_HSP", "no_Cov", "3-", "3_perc", "4=", "4_perc", "5+", "5_perc")
  for(i in 1:length(names(table(wDf$Program)))){
    tmpP <- names(table(wDf$Program))[i]
    tmpCats <- table(wDf$totCpg[wDf$Program %in% tmpP & !wDf$Hits %in% 0])
    tmpPerc <- round(100*(tmpCats/hspCats),2)
    covTypePltDf[nrow(covTypePltDf)+1, ] <- c(tmpP, length(wDf$Histotype)/4, sum(tmpCats), 
                                              tmpCats[1], tmpPerc[1], tmpCats[2], tmpPerc[2], tmpCats[3], tmpPerc[3])
  }
  
  pltDf <- data.frame(matrix(nrow=0, ncol=3))
  colnames(pltDf) <- c("Program", "Percent", "CpGs")
  
  colnames(covTypePltDf) <- gsub("_perc", "",colnames(covTypePltDf))
  inpCovLong <-  covTypePltDf %>% 
    pivot_longer(
      cols = c("3", "4", "5"), 
      names_to = "CpGs",
      values_to = "Percent"
    )
  inpCovLong$Percent <- as.numeric(inpCovLong$Percent)
  for(k in 1:nrow(inpCovLong)){
    perc <- round(inpCovLong[k,"Percent"],2)
    pltDf[nrow(pltDf)+1,] <- c(inpCovLong[k,"Program"], perc, inpCovLong[k,"CpGs"])
  }
  inpCovLong$Program <- factor(inpCovLong$Program,  levels=c("Bumphunter","DMRcate", "ProbeLasso", "Sesame"))
  
  pltDf$CpGs <- ifelse(pltDf$CpGs %in% 3, "3-", ifelse(pltDf$CpGs %in% 4, "4","5+"))
  tmpPlt <- ggplot(pltDf, aes(fill=Program, y=Percent, x=Program)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    geom_text(aes(label = Percent), size = 8, vjust = -1, colour = "black") +
    xlab("DMR-caller") +
    ylab("Percent (of total HSPs)") + 
    ylim(0,100) + 
    ggtitle(paste("Percentage of HSPs with  n=(",sum(table(wDf$totCpg)/4), ") in DMR comparisons" , sep="")) + 
    scale_fill_manual(values = c("Bumphunter" = viridis(n=6, option="F")[2],
                                 "DMRcate"=viridis(n=6, option="F")[3],
                                 "ProbeLasso"=viridis(n=6, option="F")[4],
                                 "Sesame"=viridis(n=6, option="F")[5]),
    ) +
    #scale_fill_viridis(discrete = TRUE, option = "F") +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    guides(fill=guide_legend(title="DMR-caller")) +
    facet_wrap(~factor(CpGs), nrow=3, strip.position = "bottom")
  
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(is.null(fileExt)){
    outFile <- paste(outDir, "extDMR_CpG_Len_Cov.pdf", sep="")
  }else{
    outFile <- paste(outDir, fileExt, "_extDMR_CpG_Len_Cov.pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
}

makeNonParamDMP_MULT <- function(inpBeta, inpPheno, inpCont, sigBool=NULL, inpM = NULL, pCut = NULL, bCut = NULL, adjMeth = NULL){
  refHist <- inpCont[1]
  contHist <- inpCont[2]
  message(paste("Performing DMP analysis between reference group: ", refHist, 
                " and contrast group: ", contHist, sep=""))
  if(is.null(pCut)){
    pCut <- 0.05
  }
  if(is.null(bCut)){
    bCut <- 0.2
  }
  if(is.null(adjMeth)){
    adjMeth <- "BH"
  }
  
  refSamps <- inpPheno$barcode[which(inpPheno$barcode %in% inpPheno$barcode[which(inpPheno$Histotype %in% refHist)])]
  nRefSamps <- inpPheno$barcode[which(inpPheno$barcode %in% inpPheno$barcode[which(inpPheno$Histotype %in% contHist)])]
  testBeta <- inpBeta[, testSamps] 
  message(paste("Calculating difference in beta value between groups: [", refHist, "_", contHist, "]", sep=""))
  phenoBeta <- testBeta[, refSamps]
  nPhenoBeta <- testBeta[, nRefSamps]
  # Apply trimean and get difference in beta values
  tmRef <- future.apply::future_apply(phenoBeta, 1, FUN = makeTm)
  tmNRef <- future.apply::future_apply(nPhenoBeta, 1, FUN = makeTm)
  tmDiffs <- tmRef - tmNRef
  message(paste("Performing Mann-Whitney U test for significance between groups: [", refHist, "_", contHist, "]", sep=""))
  # Perform the mann-whitney U-test o
  mDf <- log2(testBeta/(1-testBeta))
  refM <- mDf[, refSamps]
  nRefM <- mDf[, nRefSamps]
  mwM <- cbind(refM, nRefM)
  colnames(mwM) <- inpPheno$Histotype[match(colnames(mwM),inpPheno$barcode)]
  mwUTest <- future.apply::future_apply(mwM, 1, FUN = makeMannWhitP) 
  # Combine dataframes into result dataframe
  resDf <- cbind(tmDiffs, mwUTest)
  resDf <- data.frame(resDf)
  colnames(resDf) <- c("DeltaBeta", "adjP")
  resDf <- rownames_to_column(resDf, var="CpG") 
  resDf$Ref <- refHist
  resDf$Contrast <- contHist
  
  if(!is.null(sigBool)){
    resDf <- resDf[which(resDf$adjP < pCut),] 
  }
  if(!is.null(bCut)){
    resDf <- resDf[which(abs(resDf$DeltaBeta) >= bCut),]
  }
  return(resDf)
}

################################################################################
################################################################################
################################################################################
# Experimental
# Attempt to speed up DMP function 
################################################################################
################################################################################
################################################################################

nParamDmpFun <- function(bRow){
  # Get order of histotypes
  hVec <- as.character(bRow[which(names(bRow) %in% "SortCol")])
  hVec <- strsplit(hVec, ",")[[1]]
  ref <- hVec[1]
  nRef <- names(table(hVec))[!names(table(hVec)) %in% ref]
  refInds <- which(hVec %in% ref)
  nRefInds <- which(hVec %in% nRef)
  # Calculate delta-beta
  bValRow <- as.numeric(bRow[-which(names(bRow) %in% "SortCol")])
  tmRef <- makeTm(bValRow[refInds])
  nTmRef <- makeTm(bValRow[nRefInds])
  tDiff <- tmRef - nTmRef
  # Calculate brunner-munzel's p for m-values
  # Non-parametric, does not assume equal variance or size between populations
  mRow <- log2(bValRow/(1-bValRow))
  brunmunzP<- brunnermunzel::brunnermunzel.test(x=mRow[refInds], 
                                                 y=mRow[nRefInds])
  #valV <- as.numeric(mRow)
  #grpV <- factor(hVec, 
  #               levels = c( ref, nRef))
  # Retrieve unadjusted p-value 
  #brunnDf <- data.frame(valV, grpV)
  #colnames(brunnDf) <- c("M", "Group")
  #brunmunzP <- brunnermunzel::brunnermunzel.test(M ~ Group, 
  #                                           data = brunnDf)
  #
  pVal <- brunmunzP$p.value
  return(c(tDiff, pVal))
}

makeNonParamDMP_MULT_V2 <- function(inpBeta, inpPheno, inpCont, sigBool=NULL, pCut = NULL, bCut = NULL, adjMeth = NULL){
  refHist <- inpCont[1]
  contHist <- inpCont[2]
  message(paste("Performing DMP analysis between reference group: ", refHist, 
                " and contrast group: ", contHist, sep=""))
  if(is.null(pCut)){
    pCut <- 0.05
  }
  if(is.null(adjMeth)){
    adjMeth <- "BH"
  }
  # Select reference and contrast samples, bind them together
  refBeta <- inpBeta[, which(colnames(inpBeta) %in% inpPheno$barcode[which(inpPheno$Histotype %in% refHist)])]
  nRefBeta <- inpBeta[, which(colnames(inpBeta) %in% inpPheno$barcode[which(inpPheno$Histotype %in% contHist)])]
  dmpBeta <- cbind(refBeta, nRefBeta)
  # Create annotation vector to be used in downstream function, add to dataframe
  hVec <- paste(inpPheno$Histotype[match(colnames(dmpBeta),inpPheno$barcode)], collapse = ",")
  dmpBeta$SortCol <- hVec
  # Apply function for non-parametric DMP calling using future for multiprocessing
  # Warning: Long runtimes despite paralell processing
  dmpVals <- future.apply::future_apply(dmpBeta, 1, 
                                        FUN = nParamDmpFun)
  # Combine dataframes into result dataframe
  resDf <- data.frame(t(dmpVals))
  colnames(resDf) <- c("DeltaBeta", "p")
  # Correct p-values based on number of performed p-value estimations 
  resDf$adjP <- p.adjust(resDf$p, 
                         method = adjMeth,
                         n = nrow(resDf))
  # Sort based on smallest adjusted p-values
  resDf <- resDf[order(resDf$adjP, 
                       decreasing = FALSE), ]
  resDf <- rownames_to_column(resDf, var="CpG") 
  resDf$Ref <- refHist
  resDf$Contrast <- contHist
  
  # Filter away sites below adjusted p of certain value
  if(!is.null(sigBool)){
    resDf <- resDf[which(resDf$adjP < pCut),] 
  }
  # Filter away sites below absolute value of delta B between groups (trimean)
  if(!is.null(bCut)){
    resDf <- resDf[which(abs(resDf$DeltaBeta) >= bCut),]
  }
  rownames(resDf) <- resDf$CpG
  resDf <- resDf[order(resDf$adjP, decreasing = FALSE), ]
  return(resDf)
}

makeMissTypeDfPlt <- function(inpMissDmpDf){
  colCats <- unique(allMissDfLong$Program)
  colvals <- viridis(n=length(colCats), option="B")
  colVec <- colvals
  names(colVec) <- colCats
  
  inpMissDmpDf$Perc_Normal_Normal <- inpMissDmpDf$Normal_Normal/inpMissDmpDf$nMissing
  inpMissDmpDf$Perc_Normal_NotNormal <- inpMissDmpDf$Normal_NotNormal/inpMissDmpDf$nMissing
  inpMissDmpDf$Perc_NotNormal_NotNormal <- inpMissDmpDf$NotNormal_NotNormal/inpMissDmpDf$nMissing
  
  is.nan.data.frame <- function(x){
    do.call(cbind, lapply(x, is.nan))
  }
  
  inpMissDmpDf[is.nan.data.frame(inpMissDmpDf)] <- 0
  percCols <- grep("Perc", colnames(inpMissDmpDf))
  inpMissDmpDfLong <- inpMissDmpDf %>% 
    pivot_longer(
      cols = percCols, 
      names_to = "Missing_Type",
      values_to = "Percent"
    )
  inpMissDmpDfLong  <- data.frame(inpMissDmpDfLong )
  inpMissDmpDfLong$Percent <- as.numeric(inpMissDmpDfLong$Percent)
  inpMissDmpDfLong$Percent <- round(100*inpMissDmpDfLong$Percent, 2)
  inpMissDmpDfLong$Cont <- factor(inpMissDmpDfLong$Cont)
  inpMissDmpDfLong$Missing_Type <- gsub("Perc_", "",inpMissDmpDfLong$Missing_Type)
  
  tmpPlt <- ggplot(inpMissDmpDfLong, aes(fill=Program, y=Percent, x=Cont)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    geom_text(aes(label = Percent),
              position = position_dodge(width = 1),
              size = 7, 
              vjust = -1, 
              colour = "black") +
    xlab("Contrast") +
    ylab("Percent of CpG's") + 
    ylim(0,100) + 
    ggtitle(paste("Percentage of DMPs discovered by the non-parametric approach\nnot deemed significant DMPs" , sep="")) + 
    scale_fill_manual(values = colVec) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    guides(fill=guide_legend(title="Program")) +
    facet_wrap(~factor(Missing_Type), nrow=4, strip.position = "bottom")
  outFile <- paste(plotPath, "allMissingDmp.pdf", sep="")
  ggsave(outFile, 
         plot=tmpPlt, 
         width=60, height=40, units = "cm")
}