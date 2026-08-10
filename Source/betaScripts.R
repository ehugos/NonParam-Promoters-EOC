################################################################################
################################################################################
################################################################################
# Scripts associated with CpG or Beta/M-value analysis
################################################################################
################################################################################
################################################################################

makePhenoBetaFilter <- function(inpB, inpP){
  # Simple function for removing any missing samples in pheno/beta
  keepSamps <- intersect(colnames(inpB), inpP$barcode)
  inpB <- inpB[,keepSamps]
  inpP <- inpP[inpP$barcode %in% keepSamps, ]
  inpB <- inpB[,match(inpP$barcode, colnames(inpB))]
  return(list(inpB, inpP))
}

makeTm <- function(x){
  # Simple function for performing Tukeys trimean on a given vector
  # Split into three quartiles
  q <- quantile(as.numeric(x), 
                prob=c(.25,.5,.75))
  tuk <- (q[1]+(2*q[2])+q[3])/4
  tukRes <- tuk[[1]]
  return(tukRes)
}

makeMannWhitP <- function(inpRow, inpCont=NULL){
  # Simple function for performing the Mann-Whitney u-test given an input dataframe row
  # Requires column names to be in Reference and non-reference group before being submitted
  # Unless input-phenotype is supplemented (to be implemented)
  if(is.null(inpCont)){
    ref <- names(inpRow)[1]
    nRef <- names(table(names(inpRow)))[!names(table(names(inpRow))) %in% ref]
  }else{
    ref <- inpCont[[1]]
    nRef <- inpCont[[2]]
  }
  grpV <- names(inpRow)
  # Reformat the colnames of the dataframe into reference and contrast group
  grpV <- ifelse(grpV %in% ref, 
                 "ref", "cont")
  grpV <- factor(grpV, 
                 levels = c( "ref", "cont"))
  valV <- as.numeric(inpRow)
  names(valV) <- grpV
  mannWhitP <- pairwise.wilcox.test(x = valV,
                                    g = grpV,
                                    p.adjust.method = "BH",
                                    paired = FALSE)
  pVal <- mannWhitP$p.value
  return(pVal)
}

makeKW <- function(inpM, inpPheno, inpH){
  # Function for performing dunns test
  tmpPheno <- inpPheno[which(inpPheno$barcode %in% colnames(inpM)),]
  hVec <- tmpPheno$Histotype[match(colnames(inpM), tmpPheno$barcode)]
  conts <- names(table(hVec))[!names(table(hVec)) %in% inpH]
  kwPDf <- data.frame(matrix(nrow=0, 
                             ncol=2))
  colnames(kwPDf) <- c("CpG", "pVal")
  
  for(m in 1:nrow(inpM)){
    tmpCpg <- rownames(inpM)[m]
    tmpC <- inpM[m,]
    kwDf <- rbind(tmpC, hVec)
    kwDf <- data.frame(t(kwDf))
    colnames(kwDf) <- c("M", "Histotype")
    # Set levels (first level is reference)
    kwDf$Histotype <- factor(kwDf$Histotype, 
                             levels=c(inpH, conts))
    kwTest <- kruskal.test(x=as.numeric(kwDf$M), 
                           g=kwDf$Histotype)
    kwP <- kwTest$p.value
    # If there is a significance, we perform Dunn's test to check differences to the reference, between groups
    kwPDf[nrow(kwPDf)+1, ] <- c(tmpCpg, 
                                kwP) 
  }
  return(kwPDf)
}

makeDunn <- function(inpM, inpPheno, inpH, inpMethod = NULL){
  if(is.null(inpMethod)){
    inpMethod <- "none"
  }
  # Function for performing dunns test
  tmpPheno <- inpPheno[which(inpPheno$barcode %in% colnames(inpM)),]
  hVec <- tmpPheno$Histotype[match(colnames(inpM), tmpPheno$barcode)]
  conts <- names(table(hVec))[!names(table(hVec)) %in% inpH]
  dunnProbeDf <- data.frame(matrix(nrow=0, 
                                   ncol=3))
  colnames(dunnProbeDf) <- c("CpG", "Contrast", "pVal")
  
  for(m in 1:nrow(inpM)){
    tmpCpg <- rownames(inpM)[m]
    tmpC <- inpM[m,]
    dDf <- rbind(tmpC, hVec)
    dDf <- data.frame(t(dDf))
    colnames(dDf) <- c("M", "Histotype")
    # Set levels (first level is reference)
    dDf$Histotype <- factor(dDf$Histotype, 
                            levels=c(inpH, conts))
    dunnCW <- DunnTest(x=as.numeric(dDf$M), 
                       g=dDf$Histotype, 
                       method = inpMethod)
    dDf <- data.frame(dunnCW[[1]])
    for(n in 1:nrow(dDf)){
      dunnProbeDf[nrow(dunnProbeDf)+1, ] <- c(tmpCpg, rownames(dDf)[n],dDf$pval[n]) 
    }
    #tmpP <- keepDunn$pval
    #dunnProbeDf[m, ] <-  tmpP
  }
  if(is.null(inpMethod)){
    # Adjust p-value based on the total number of tests performed in the promoter region
    dunnProbeDf$pAdj <- p.adjust(dunnProbeDf$pVal, method="BH")
    # Select tests involving the reference histotype
  }else{
    dunnProbeDf$pAdj <- dunnProbeDf$pVal
  }
  dunnRows <- dunnProbeDf[grep(inpH, dunnProbeDf$Contrast), ]
  return(dunnRows)
}

makeBetaGrpDists <- function(phenoBetas, nPhenoBetas, tmpPheno, tmpHistotypes, tmpName, histBool=NULL, distType=NULL){
  if(is.null(distType)){
    distType <- "TRIMEAN"
  }
  if(!is.null(histBool)){
    distConts <- tmpHistotypes[!tmpHistotypes %in% tmpName]
    distDf <- data.frame(matrix(nrow=nrow(phenoBetas), ncol=(length(distConts))))
    colnames(distDf) <- distConts
  }else{
    distDf <- data.frame(matrix(nrow=nrow(phenoBetas), ncol=1))
    distConts <- c("all")
    colnames(distDf) <- distConts
  }
  rownames(distDf) <- rownames(phenoBetas)
  for(l in 1:length(distConts)){
    tmpCont <- distConts[l]
    if(tmpCont %in% "all"){
      tmpNPB <- nPhenoBetas
    }else{
      tmpNPB <- nPhenoBetas[,match(tmpPheno$Sample_ID[tmpPheno$Histotype %in% tmpCont], colnames(nPhenoBetas))]
    }
    for(m in 1:nrow(phenoBetas)){
      tmpPBr <- as.data.frame(phenoBetas[m,])
      tmpPBr  <- as.numeric(tmpPBr[ , colSums(is.na(tmpPBr)) == 0])
      tmpNPBr <- as.data.frame(tmpNPB[m,])
      tmpNPBr <- as.numeric(tmpNPBr[ , colSums(is.na(tmpNPBr)) == 0])
      # colnames(tmpPBr) <- colnames(tmpNPBr)
      if(toupper(distType) %in% "MEDIAN"){
        # Get distance between the median methylated probes 
        medDist <- stats::dist(rbind(median(tmpPBr), median(tmpNPBr)), method = "euclidean")
        distDf[m,tmpCont] <- medDist[1]
      }else if(toupper(distType) %in% "GEOMEAN"){
        # Take the geometric mean from beta-distribution (less sensitive to outliers)
        pBmean <- exp(mean(log(as.numeric(tmpPBr[as.numeric(tmpPBr)>0]))))
        nPBmean <- exp(mean(log(as.numeric(tmpNPBr[as.numeric(tmpNPBr)>0]))))
        meanDist <- stats::dist(rbind(as.numeric(pBmean), as.numeric(nPBmean)), method = "euclidean")
        distDf[m,tmpCont] <- meanDist[1]
      }else if(toupper(distType) %in% "MEAN"){
        pBmean <- mean(tmpPBr)
        nPBmean <- mean(tmpNPBr)
        meanDist <- stats::dist(rbind(as.numeric(pBmean), as.numeric(nPBmean)), method = "euclidean")
        distDf[m,tmpCont] <- meanDist[1]
      }else if(toupper(distType) %in% "TRIMEAN"){
        # Get the trimean (q1+2q2+q3)/4
        # q1 = lower 25%, q2 = 25-75%, q3 = upper 75%  
        pBmean <- litteR::trimean(tmpPBr)
        nPBmean <- litteR::trimean(tmpNPBr)
        meanDist <- stats::dist(rbind(as.numeric(pBmean), 
                                      as.numeric(nPBmean)), 
                                method = "euclidean")
        distDf[m,tmpCont] <- meanDist[1]
      }
    }
  }
  return(distDf)
}

addPromoToGene_Mult <- function(inpGenes, nCores = NULL){
  if(is.null(nCores)){
    nCores <- 1
  }
  plan(multisession, workers = nCores)
  # Add promoter information to genes
  inpGenes$promo_start <- NA
  inpGenes$promo_end <- NA
  nIt <- nrow(inpGenes)
  nItDec <- nIt/10
  inpGenes <- future_lapply(1:nIt, function(i){
    if(abs(i)%%round(nItDec) == 0 || is.na(abs(i)%%round(nItDec))){
      message(paste("Processing gene: [", i, "/", nIt,"]"))
    }
    geneRow <- inpGenes[i,]
    # Add promoter (2000bp Upstream of TSS, 200bp downstream)
    if(geneRow$strand == "-"){
        geneRow$promo_start <- geneRow$end - 200 
        geneRow$promo_end <- geneRow$end + 2000
      }else{
        geneRow$promo_start <- geneRow$start - 2000
        geneRow$promo_end <- geneRow$start + 200
      }
    inpGenes[i,] <- geneRow
  })
  inpGenes <- do.call(rbind.data.frame, inpGenes) 
  return(inpGenes)
}

makeHg19GeneCoords <- function(inpGenes){
  # Function for mapping hg38 coordinates to their hg19 positions
  # Requires a dataframe with ensembl-gene id's as input
  tmpMart = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                   host="https://grch37.ensembl.org", 
                   path="/biomart/martservice", 
                   dataset="hsapiens_gene_ensembl")
  # Get coordinate data for all genes in the input list
  hg19Coords <- getBM(attributes = c("ensembl_gene_id",
                                     "chromosome_name",
                                     "start_position",
                                     "end_position",
                                     "strand"),
  filters = "ensembl_gene_id", 
  values = inpGenes$ensembl_gene_id,
  mart = tmpMart)
  keepGenes <- intersect(hg19Coords$ensembl_gene_id, inpGenes$ensembl_gene_id)
  hg19Coords <- hg19Coords[hg19Coords$ensembl_gene_id %in% keepGenes, ]
  inpGenes <- inpGenes[inpGenes$ensembl_gene_id %in% keepGenes, ]
  hg19Coords$strand[hg19Coords$strand == "1"] <- "+"
  hg19Coords$strand[hg19Coords$strand == "-1"] <- "-"
  # Transer new data to input dataframe
  inpGenes <- inpGenes[!is.na(inpGenes$ensembl_gene_id),]
  hg19Coords <- hg19Coords[match(inpGenes$ensembl_gene_id, hg19Coords$ensembl_gene_id),]
  inpGenes$start <- hg19Coords$start_position
  inpGenes$end <- hg19Coords$end_position
  inpGenes$strand <- hg19Coords$strand
  inpGenes$chr <- paste("chr", hg19Coords$chromosome_name,sep="")
  return(inpGenes)
}

# Create list of promoters from dataframe compatible with scripts
makePromoLst <- function(inpDeg, inpPromoDf, cpgInp){
  tmpGenes <- inpPromoDf[inpPromoDf$Gene %in% inpDeg$ensembl_gene_id,]
  tmpLst <- list()
  for(i in 1:nrow(tmpGenes)){
    if(!nrow(tmpGenes)<10){
      if(abs(i%%round(nrow(tmpGenes))/10) == 0){
        message(paste("Processing gene: [", i, "/", nrow(tmpGenes),"]"))
      }
    }
    tmpRow <- tmpGenes[i,]
    tmpName <- tmpRow$Gene
    cpgInd <- grep("cpg", colnames(tmpRow))
    cpgCols <- tmpRow[cpgInd]
    cpgCols <- cpgCols[!is.na(cpgCols)]
    cpgDat <- cpgInp[cpgCols,]
    tmpLst[[i]] <- cpgDat
    names(tmpLst)[i] <- tmpName
  }
  return(tmpLst)
}

makeBetaGrpVar <- function(cpgDf, tmpHistotypes, tmpPheno, distType=NULL){
  # Return the standard deviation between groups in a dataframe with respect to the chosen method for averaging
  if(is.null(distType)){
    distType <- "TRIMEAN"
  }
  varDf <- data.frame(matrix(nrow=nrow(cpgDf), ncol=(length(tmpHistotypes))))
  colnames(varDf) <- tmpHistotypes
  rownames(varDf) <- rownames(cpgDf)
  for(l in 1:length(tmpHistotypes)){
    tmpH <- tmpHistotypes[l]
    phenoBetas <- cpgDf[,match(tmpPheno$barcode[tmpPheno$Histotype %in% tmpH], colnames(cpgDf))]
    for(m in 1:nrow(phenoBetas)){
      tmpPB <- as.data.frame(phenoBetas[m,])
      tmpPB  <- as.numeric(tmpPB[ , colSums(is.na(tmpPB)) == 0])
      # Get median of histotype CpG site methylation
      if(toupper(distType) %in% "MEDIAN"){
        pBm <- median(tmpPB)
      }else if(toupper(distType) %in% "GEOMEAN"){
        pBm  <- exp(mean(log(as.numeric(tmpPB[as.numeric(tmpPB)>0]))))
      }else if(toupper(distType) %in% "TRIMEAN"){
        pBm <- makeTm(tmpPB)
      }else if(toupper(distType) %in% "MEAN"){
        pBm <- mean(tmpPB)
      }
      nMedPB <- tmpPB[!tmpPB %in% pBm]
      medPBVar <- sum(abs(nMedPB-pBm)/length(nMedPB))
      varDf[m,tmpH] <- medPBVar
    }
  }
  return(varDf)
}

makeTriMeanDf <- function(inpB, inpPheno){
  # Simple script for creating a dataframe of trimean values
  outDf <- matrix(nrow=nrow(inpB), ncol=length(table(inpPheno$Histotype)))
  colnames(outDf) <- names(table(inpPheno$Histotype))
  rownames(outDf) <- rownames(inpB)
  for(i in 1:length(table(inpPheno$Histotype))){
    tmpH <- names(table(inpPheno$Histotype))[i]
    tmpColNames <- inpPheno$barcode[which(inpPheno$Histotype %in% tmpH)]
    tmpB <- inpB[, tmpColNames]
    for(j in 1:nrow(tmpB)){
      tmpName <- rownames(tmpB)[j]
      betaCol <- tmpB[j,]
      betaCol <- na.omit(betaCol)
      tmpP <- makeTm(as.numeric(betaCol))
      outDf[j, tmpH] <- tmpP
    }
  }
  return(outDf)
}

makeCpgDf <- function(inpCpg, tmpBeta, tmpPheno, inpHist){
  tmpPheno <- tmpPheno[tmpPheno$barcode %in% colnames(tmpBeta), ]
  tmpBeta <- tmpBeta[,colnames(tmpBeta) %in% tmpPheno$barcode]
  # cpgCols <- which(colnames(inpCpg) %like% "cpg_")
  cpgCols <- which(grepl("cpg_",colnames(inpCpg)))
  cpgVec <- c()
  for(i in 1:nrow(inpCpg)){
    tmpGRow <- inpCpg[i,]
    cpgOls <- tmpGRow[,cpgCols]
    cpgOls <- cpgOls[!is.na(cpgOls)]
    cpgVec <- append(cpgVec, cpgOls)
  }
  cpgVec <- cpgVec[!duplicated(cpgVec)]
  cpgVec <- cpgVec[cpgVec %in% rownames(tmpBeta)]
  outBeta <- tmpBeta[cpgVec,]
  return(outBeta)
}

makeGeneMatch <- function(inpExp, inpGenes){
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
    geneDf <- geneDf[rowSums(is.na(geneDf)) < 3,]
    geneDf <- geneDf[order(rowSums(geneDf[, names(expLst)])),]
    outLst[[i]] <- geneDf
    names(outLst)[i] <-tmpH
  }
  return(outLst)
}

getGeneCpgBeta <- function(geneCpgInp, inpBeta, inpPheno, cpgInf, 
                           inpFoc = NULL, revBool=NULL, printBool = NULL, allBool = NULL){
  # Requires inpPheno to have a grouping column names "Histotype"
  # Iterates through promoter CpG dataframe, returns list of dataframes with matched beta values for a histotype of interest
  
  # First check if we are selecting all available samples
  if(!is.null(allBool)){
    revBool <- NULL
    inpFoc <- names(table(inpPheno$Histotype))
  }else{
    # If we are not choosing all samples, check that we have a phenotype group of interest
    if(is.null(inpFoc)){
      message("No phenotype of interest chosen, selecting the most prevalent group")
      inpFoc <- names(sort(table(inpPheno$Histotype), decreasing = TRUE))[1]
    }
    # Check if inpFoc is of single histotype or contrast type
    if(grepl("_", inpFoc, fixed=TRUE)){
      inpFoc <- strsplit(inpFoc, "_")[[1]]
    }
  }
  betaLst <- list()
  # Only select samples of the chosen histotype(s)
  # RevBool selects all histotypes NOT in the chosen histotype 
  # To retrieve ALL CpG-Betas, we thus name our input object (inpFoc) something not in histotypes
  if(!is.null(revBool)){
    infSamples <- inpPheno[!inpPheno$Histotype %in% inpFoc,]
  }else{
    infSamples <- inpPheno[inpPheno$Histotype %in% inpFoc,]
  }
  for(j in 1:nrow(geneCpgInp)){
    if(is.null(printBool)){
      if(is.na(abs(j)%%round(nrow(geneCpgInp)/10)) || abs(j)%%round(nrow(geneCpgInp)/10) == 0){
        message(paste("Processing gene: [", j, "/", nrow(geneCpgInp),"]"))
      }
    }
    # Get CpG-Gene dataframe for gene j
    geneRow <- geneCpgInp[j, ]
    # Extract CpG names
    geneName <- geneRow$Gene
    geneCpgs <- geneRow[,!colnames(geneRow) %in% c("Gene")]
    geneCpgs <- geneCpgs[!is.na(geneCpgs)]
    geneCpgs <- geneCpgs[which(geneCpgs %in% rownames(inpBeta))]
    if(length(geneCpgs) == 0){
      #message(paste("No CpGs found in dataframe for gene: ", geneName, sep=""))
      next()
    }
    # Extract cpg values from the beta matrix using rownames (CpGs) for gene
    betaTmp <-  inpBeta[geneCpgs,]
    # Remove NA-rows (i.e. no such probe in inpBeta)
    betaTmp <- betaTmp[!is.na(row.names(betaTmp)), ]
    # Grab positional data for the CpGs
    infDf <- cpgInf[rownames(betaTmp), ]
    # Keep only columns in beta-matrix which are of the correct histotype(s)
    betaTmp <- betaTmp[,which(colnames(betaTmp) %in% infSamples$barcode)]
    # Remove samples with NA
    betaTmp <- betaTmp[, colSums(is.na(betaTmp)) == 0]
    #  Add histotype to sample (column) name
    if(!is.null(revBool)){
      # If we have all samples except the histotype of interest, we must add them with respect to each samples histotype
      for(k in 1:ncol(betaTmp)){
        cHisto <- infSamples$Histotype[which(infSamples$barcode == colnames(betaTmp)[k])]
        colnames(betaTmp)[k] <- paste(cHisto, colnames(betaTmp)[k], sep="_")
      }
    }else{
      if(!is.null(allBool)){
        hSampVec <- inpPheno$Histotype[match(colnames(betaTmp), inpPheno$barcode)] 
        colnames(betaTmp) <- paste(hSampVec, colnames(betaTmp), sep="_")
      }else{
        # If not, we simply add the focus group of interest
        colnames(betaTmp) <- paste(inpFoc, colnames(betaTmp), sep="_")
      }
    }
    # Merge dataframes by CpG id
    geneCpgs  <- merge(infDf, betaTmp, by = 'row.names', all = TRUE)
    # Order seen to position
    geneCpgs <- geneCpgs[order(geneCpgs$pos),]
    betaLst[[j]] <- geneCpgs   
    names(betaLst)[j] <- geneName
  }
  # Remove empty entries from list
  betaLst <- betaLst[!names(betaLst) %in% ""]
  betaLst <- betaLst[!is.na(names(betaLst))]
  return(betaLst)
}

getGeneCpgBeta_MULT <- function(geneCpgInp, inpBeta, inpPheno, cpgInf, 
                                inpFoc = NULL, revBool=NULL, allBool = NULL, 
                                printBool = NULL, nCores=NULL){
  # Uses gene CpG dataframe as input
  # Requires inpPheno to have a grouping column names "Histotype"
  # Iterates through promoter CpG dataframe, returns list of dataframes with matched beta values for a histotype of interest
  # First check if we are selecting all available samples
  if(!is.null(allBool)){
    revBool <- NULL
    inpFoc <- names(table(inpPheno$Histotype))
    if(is.null(printBool)){
      message(paste("Getting CpG Beta values for all histotypes"))
    }
  }else{
    # If we are not choosing all samples, check that we have a phenotype group of interest
    if(is.null(inpFoc)){
      message("No phenotype of interest chosen, selecting the most prevalent group")
      inpFoc <- names(sort(table(inpPheno$Histotype), decreasing = TRUE))[1]
    }
    # Check if inpFoc is of single histotype or contrast type
    if(grepl("_", inpFoc, fixed=TRUE)){
      inpFoc <- strsplit(inpFoc, "_")[[1]]
    }
    if(is.null(printBool)){
      message(paste("Getting CpG Beta values for histotype:", inpFoc))
    }
  }
  # Only select samples of the chosen histotype(s)
  # RevBool selects all histotypes NOT in the chosen histotype 
  # To retrieve ALL CpG-Betas, we thus name our input object (inpFoc) something not in histotypes
  if(!is.null(revBool)){
    infSamples <- inpPheno[!inpPheno$Histotype %in% inpFoc,]
  }else{
    infSamples <- inpPheno[inpPheno$Histotype %in% inpFoc,]
  }
  nIt <- nrow(geneCpgInp)
  refCpg <- rownames(inpBeta)
  
  betaLst <- list()
  betaLst <- future_lapply(1:nIt, function(j){
    if(is.na(abs(j)%%round(nrow(geneCpgInp)/10)) || abs(j)%%round(nrow(geneCpgInp)/10) == 0){
      if(is.null(printBool)){
        message(paste("Processing gene: [", j, "/", nrow(geneCpgInp),"]"))
      }
    }
    # Get CpG-Gene dataframe for gene j
    geneRow <- geneCpgInp[j, ]
    # Extract CpG names
    geneName <- geneRow$Gene
    geneCpgs <- geneRow[,!colnames(geneRow) %in% c("Gene")]
    geneCpgs <- geneCpgs[!is.na(geneCpgs)]
    geneCpgs <- geneCpgs[which(geneCpgs %in% refCpg)]
    if(length(geneCpgs) == 0){
      #if(is.null(printBool)){
      #    message(paste("No CpGs found in dataframe for gene: ", geneName, sep=""))
      #}
      betaLst[[j]] <- NULL
    }else{
      # Extract cpg values from the beta matrix using rownames (CpGs) for gene
      betaTmp <-  inpBeta[geneCpgs,]
      # Keep only columns in beta-matrix which are of the correct histotype(s)
      betaTmp <- betaTmp[,which(colnames(betaTmp) %in% infSamples$barcode)]
      # Remove NA-rows (i.e. no such probe in inpBeta)
      betaTmp <- na.omit(betaTmp)
      # Grab positional data for the CpGs
      infDf <- cpgInf[which(rownames(cpgInf) %in% rownames(betaTmp)), ]
      infDf <- infDf[match(rownames(betaTmp), rownames(infDf)), ]
      #  Add histotype to sample (column) name
      if(!is.null(revBool)){
        # If we have all samples except the histotype of interest, we must add them with respect to each samples histotype
        for(k in 1:ncol(betaTmp)){
          cHisto <- infSamples$Histotype[which(infSamples$barcode == colnames(betaTmp)[k])]
          colnames(betaTmp)[k] <- paste(cHisto, colnames(betaTmp)[k], sep="_")
        }
      }else{
        if(!is.null(allBool)){
          hSampVec <- inpPheno$Histotype[match(colnames(betaTmp), inpPheno$barcode)] 
          colnames(betaTmp) <- paste(hSampVec, colnames(betaTmp), sep="_")
        }else{
          # If not, we simply add the focus group of interest
          colnames(betaTmp) <- paste(inpFoc, colnames(betaTmp), sep="_")
        }
      }
      # Merge dataframes by CpG id
      geneCpgs  <- merge(infDf, 
                         betaTmp, by = 'row.names', all = TRUE)
      # Remove first column, set as rownames
      cpgNames <- geneCpgs[,1]
      geneCpgs <- geneCpgs[,-1]
      rownames(geneCpgs) <- cpgNames
      # Order seen to position
      keepCpgs <- geneCpgs[order(geneCpgs$pos),]
      betaLst[[j]] <- keepCpgs  
    }
  })
  names(betaLst) <- geneCpgInp$Gene
  # Remove empty entries from list
  remInds <- which(lengths(betaLst) == 0)
  betaLst <- betaLst[-remInds]
  return(betaLst)
}

getGeneInf <- function(geneName, geneInf, inpTx, filterType="HGNC"){
  # Retrieve gene data from ensembl (grch38) 
  # Match input type to filter from biomart query
  if(toupper(filterType) =="ENSEMBL"){
    filterVar <- "ensembl_gene_id"
  }else if(toupper(filterType) =="HGNC"){
    filterVar <- "external_gene_name"
  }else if(toupper(filterType) =="UCSC"){
    filterVar <- "txGene"
  }else{
    stop("Incorrect filter type chosen, please choose HGNC, UCSC or ENSEMBL")
  }
  # Get gene-data from UCSC
  txRowInd <- which(names(inpTx) == geneName) 
  # Get gene-data from ensembl
  ensInd <- which(colnames(geneInf) %in% filterVar)
  # Get chromosome, gene coordinates for DEG genes
  ensRowInd <- which(geneInf[ensInd] == geneName)
  if(!length(txRowInd) == 0 ){
    txRow <- inpTx[txRowInd]
    txChr <- txRow@seqnames
    txStart <- txRow@ranges@start
    txWidth <- txRow@ranges@width
    txEnd <- txStart + txWidth -1
    txRow <- list("chr"= txChr, "start"= txStart, "end"=txEnd, "width"=txWidth, "txGene"=geneName)
    txDf <- data.frame(txRow)
    if(!length(ensRowInd) == 0 ){
      ensRow <- geneInf[ensRowInd,]
      ensRow$width <- abs(ensRow$start - ensRow$end)
      txLoc <- which(names(inpTx) == geneName)
      if(length(txLoc) == 0){
        tmpStrand <- "*"
      }else{
        tmpStrand <- as.character(inpTx[txLoc]@strand@values) 
      }
      ensRow$strand <- tmpStrand
    }
    if(length(ensRowInd) == 0 && length(txRowInd) ==0){
      message(paste("Gene-data missing for gene", geneName))
      return(NA)
    }
    
    if(filterVar == "txGene"){
      colnames(ensRow) <- paste("ENS.", colnames(ensRow), sep="")
      geneDat <- cbind(txDf, ensRow)
    }else{
      colnames(txDf) <- paste("TX.", colnames(txDf), sep="")
      geneDat <- cbind(ensRow, txDf)
    }
    return(geneDat)
  }
}

makeGeneRegOverlap <- function(degDf, dmrDf){
  # Function for matching genes to regions
  degDmrLst <- list()
  ind <- 1
  # Remove chr part of chromosome name to enable matching
  dmrDf$seqnames <- gsub("chr", "", dmrDf$seqnames)
  # Adjust strand to be in the correct format
  degDf$strand_orientation <- gsub("-1", "-", degDf$strand_orientation)
  degDf$strand_orientation <- gsub("1", "+", degDf$strand_orientation)
  degDf$chromosome_name <- as.character(degDf$chromosome_name)
  
  # Keep only data with a shared chromosome between DEG/DMR to avoid analsysis with different levels
  sharedSeqs <- sort(intersect(unique(dmrDf$seqnames), unique(degDf$chromosome_name)))
  degDf <- degDf[degDf$chromosome_name %in% sharedSeqs,]
  dmrDf <- dmrDf[dmrDf$seqnames %in% sharedSeqs,]
  
  # Create G-ranges object for DEG
  degRanges <- makeGRangesFromDataFrame(degDf,
                                        keep.extra.columns=FALSE,
                                        ignore.strand=FALSE,
                                        seqinfo=NULL,
                                        seqnames.field=c("chromosome_name"),
                                        start.field="start_position",
                                        end.field=c("end_position"),
                                        strand.field="strand_orientation",
                                        starts.in.df.are.0based=FALSE)
  names(degRanges) <- degDf$ensembl_gene_id
  # Create G-ranges object for DMR
  dmrRanges <- makeGRangesFromDataFrame(dmrDf,
                                        keep.extra.columns=FALSE,
                                        ignore.strand=FALSE,
                                        seqinfo=NULL,
                                        seqnames.field=c("seqnames"),
                                        start.field="start",
                                        end.field=c("end"),
                                        strand.field="strand",
                                        starts.in.df.are.0based=FALSE)
  names(dmrRanges) <- dmrDf$overlapping.genes
  indOverlaps <- findOverlaps(degRanges, dmrRanges)
  # We sometimes hit multiple highly similar ranges, thus we filter these out
  indOverlaps <- indOverlaps[!duplicated(indOverlaps@from), ]
  # Get Rows
  degRows <- degDf[indOverlaps@from,]
  dmrRows <- dmrDf[indOverlaps@to,]
  # As row-matrix is of identical dimensionality, we simply cbind it
  degDmrDf <- as.data.frame(cbind(degRows, dmrRows))
  return(degDmrDf)
}

makeCpgGeneOverlap_MULT <- function(inpGenes = NULL, cpgInp, inpGeneInf, 
                                    type=NULL, strandBool = NULL){
  # Function for extracting the cpg-positions from DMRcate genomicranges output
  # To utilize the filter function of dplyr, we need to have it as a dataframe and not DFrame
  if(is.null(inpGenes)){
    inpGenes <- inpGeneInf$ensembl_gene_id
    type <- NULL
  }
  if(is.null(type)){
    type="GENE"
  }
  
  degDf <- inpGeneInf[inpGeneInf$ensembl_gene_id %in% inpGenes, ]
  if(nrow(degDf) == 0){
    message("No supplied genes found in geneInf")
    return()
  }
  cpgInp <- data.frame(cpgInp)
  # Requires named list of dataframe converted Granges results 
  if(!"strand" %in% colnames(degDf)){
    stop("Error, no strand information supplied in input-matrix")
  } 
  # Extract name, split into separate contrast names for comparative purposes
  message("Getting CpG matches for chosen genes")
  degDf <- degDf[!is.na(degDf$ensembl_gene_id),]
  nIt <- nrow(degDf)
  if(is.null(strandBool)){
    strandIgn <- TRUE
  }else{
    strandIgn = FALSE
  }
  # Create ranges object based on CpG input before the loop to save time
  cpgRanges <- makeGRangesFromDataFrame(cpgInp,
                                        ignore.strand=strandIgn,
                                        seqinfo=NULL,
                                        start.field="start",
                                        end.field="end",
                                        strand.field = "strand",
                                        seqnames.field="chr") 
  geneCpgLst <- list()
  geneCpgLst <- future_lapply(1:nIt, function(i){
    if(abs(i)%%round(nIt/10) == 0 || is.na(abs(i)%%round(nIt/10))){
      message(paste("Processing gene: [", i, "/", nIt,"]"))
    }
    # Extract gene row
    geneRow <- degDf[i,]
    if(toupper(type)=="PROMOTER"){
      # Get start position, end position of gene
      # For PROMOTER, we choose the region that is within 5200 bp upstream of TSS
      geneStart <- geneRow$promo_start
      geneEnd <- geneRow$promo_end
    }else if(toupper(type)=="GENE"){
      if(geneRow$strand == "-"){
        # If gene, we want to include gene body and TSS[-200,2000]
        geneStart <- geneRow$start
        geneEnd <- geneRow$promo_end
      }else{
        # If gene, we want to include gene body and TSS[-200,2000]
        geneStart <- geneRow$promo_start
        geneEnd <- geneRow$end
      }
    }else if(toupper(type)=="BODY"){
      # If gene, we want to include gene body and TSS[-200,2000]
      geneStart <- geneRow$start
      geneEnd <- geneRow$end
    }else{
      stop("Incorrect overlap type, please choose GENE or PROMOTER as input")
    }
    # Call getGeneCpgs to retrieve overlapping CPGs
    posInp <- c(geneStart, geneEnd)
    cpgMatches <- getGeneCpgs(geneRow=geneRow, 
                              cpgInp=cpgInp,
                              posInp=posInp,
                              rangesInp = cpgRanges,
                              strandBool = strandBool)
    if(!is.null(cpgMatches)){
      cpgMatches <- cpgMatches[order(cpgMatches$pos, decreasing = FALSE),]
    }
    geneCpgLst[[i]] <- cpgMatches
    ## Check if we have results, if we do save result under gene-name 
    #if(!is.null(cpgMatches)){
    #  # Sort CpG matches seen to their position (i.e. from smallest to largest)
    #  cpgMatches <- cpgMatches[order(cpgMatches$pos, decreasing = FALSE),]
    #  geneCpgLst[j] <- cpgMatches
    #  # names(geneCpgLst)[j] <- geneRow$ensembl_gene_id
    #}
  })
  names(geneCpgLst) <- degDf$ensembl_gene_id
  # Filter away null entries 
  geneCpgLst <- geneCpgLst[lengths(geneCpgLst) != 0]
  return(geneCpgLst)
}

makeGeneCpGLocDataframe <- function(inpLst){
  # Script for turning list of gene dataframes into one comprehensive dataframe listing all promoter Cpg's
  promoDf <- data.frame(matrix(ncol=1, nrow=length(inpLst)))
  colnames(promoDf) <- "Gene"
  for(i in 1:length(inpLst)){
    tmpGene <- rownames(inpLst[[i]])
    if(length(tmpGene) > ncol(promoDf)-1){
      pDiff <- abs(ncol(promoDf)-1 - length(tmpGene))
      for(j in 1:pDiff){
        colN <- paste("cpg", ncol(promoDf), sep="_")
        promoDf <- promoDf %>% tibble::add_column(!! paste0("cpg_", ncol(promoDf)) := NA, .after = ncol(promoDf))
      }
    }
    if(length(tmpGene) < ncol(promoDf)-1){
      diff <- ncol(promoDf)-1-length(tmpGene)
      tmpGene <- c(tmpGene, rep(NA,diff))
    }
    tmpGene <- c(names(inpLst)[i],tmpGene)
    promoDf[i,] <- tmpGene
  }
  return(promoDf)
}

makeGeneCpGLocDataframe_MULT <- function(inpLst){
  # Script for turning list of gene dataframes into one comprehensive dataframe listing all promoter Cpg's
  promoDf <- data.frame(matrix(ncol=1, nrow=length(inpLst)))
  colnames(promoDf) <- "Gene"
  nIt <- length(inpLst)
  cpgLst <- list()
  cpgLst <- future_lapply(1:nIt, function(i){
    if(abs(i)%%round(nIt/10) == 0){
      message(paste("Processing gene: [", i, "/", nIt,"]"))
    }
    tmpGene <- rownames(inpLst[[i]])
    #if(length(tmpGene) > ncol(promoDf)-1){
    #  pDiff <- abs(ncol(promoDf)-1 - length(tmpGene))
    #  for(j in 1:pDiff){
    #    colN <- paste("cpg", ncol(promoDf), sep="_")
    #    promoDf <- promoDf %>% tibble::add_column(!! paste0("cpg_", ncol(promoDf)) := NA, .after = ncol(promoDf))
    #  }
    #}
    #if(length(tmpGene) < ncol(promoDf)-1){
    #  diff <- ncol(promoDf)-1-length(tmpGene)
    #  tmpGene <- c(tmpGene, rep(NA,diff))
    #}
    #tmpGene <- c(names(inpLst)[i],tmpGene)
    cpgLst[[i]] <- tmpGene
  })
  maxL <- max(sapply(cpgLst, length))
  ## Add NA values to list elements so vectors are of the same length
  cpgLst <- lapply(cpgLst, function(v) { c(v, rep(NA, maxL-length(v)))})
  promoDf <- do.call(rbind, cpgLst)
  # Add CpG index as column names
  colnames(promoDf) <- paste("cpg_", 1:ncol(promoDf), sep="")
  promoDf <- as.data.frame(promoDf)
  # Add gene-name as first column
  promoDf <- cbind("Gene" = names(inpLst), promoDf)
  return(promoDf)
}

rankGeneMeanVar <- function(inpGeneLst, inpGeneInf){
  geneDf <- data.frame(matrix(nrow=length(inpGeneLst), ncol=7))
  colnames(geneDf) <- c("ensembl_gene_id","external_gene_name", "maxStdGene", "maxDiff", "stdMean", "meanValVec", "sidName")
  for(i in 1:length(inpGeneLst)){
    tmpEns <- names(inpGeneLst)[i]
    tmpSym <- inpGeneInf$external_gene_name[which(inpGeneInf$ensembl_gene_id %in% tmpEns)]
    tmpBetaDf <- inpGeneLst[[tmpEns]]
    meanCols <- grep("Mean", colnames(tmpBetaDf))
    if(length(meanCols) <= 1){
      next()
    }
    stdCols <- grep("STD", colnames(tmpBetaDf))
    meanDf <- tmpBetaDf[, meanCols]
    rownames(meanDf) <- tmpBetaDf$CpG
    # stdDf <-  tmpBetaDf[, stdCols]
    # Get train CpGs and not-tran CPGs
    meanTrain <- data.frame(meanDf[,which(colnames(meanDf) %in% "TRAIN_Mean")]) 
    colnames(meanTrain) <- "TRAIN_Mean"
    rownames(meanTrain) <- rownames(meanDf) 
    meanNTrain <- meanDf[,which(!colnames(meanDf) %in% "TRAIN_Mean")]
    stdDf <- data.frame(matrix(ncol=4, nrow=length(meanTrain)))
    colnames(stdDf) <- c("CpG", "TrainMean", "DevExt", "extMeanVec")
    # Calculate standard deviance based on the geometric mean between train and ext
    # I.e. we calculate the standard deviance using the train-data rather then the actual STD 
    # That is to say the deviance from the training cohort, not the deviance from the population mean
    
    if(is.null(ncol(meanNTrain))){
      meanNTrain <- data.frame(meanNTrain)
      rownames(meanNTrain) <- rownames(meanDf)
      colnames(meanNTrain) <- colnames(meanDf)[which(!colnames(meanDf) %in% "TRAIN_Mean")]
    }
    for(k in 1:nrow(meanTrain)){
      distVec <- c()
      tmpTrainMean <- meanTrain[k,]
      tmpCpgName <- rownames(meanNTrain)[k]
      tmpExtMeanVec <- meanNTrain[k,]
      for(l in 1:length(tmpExtMeanVec)){
        tmpExtMean <- tmpExtMeanVec[l]
        if(is.na(tmpExtMean) | is.na(tmpTrainMean)){
          tmpDev <- NA
        }else{
          tmpDev <- sqrt((as.numeric(tmpTrainMean) - as.numeric(tmpExtMean))^2)
        }
        distVec <- append(distVec, tmpDev)
      }
      distVec <- na.omit(distVec)
      extStd <- sum(distVec/length(distVec))
      stdDf[k,] <- c(tmpCpgName, tmpTrainMean, extStd, paste(tmpExtMeanVec, collapse=","))
    }
    # Remove entries without match in Train (as we use Train as reference)
    stdDf <- stdDf[!is.na(stdDf$TrainMean),]
    # We do NOT use geom-m mean here, as we want it to be sensitive for outliers
    meanStd <- litteR::trimean(as.numeric(stdDf$DevExt))
    maxStd <- max(stdDf$DevExt)
    # Calculate geometric mean for each dataset
    triMeanVec <- c()
    for(j in 1:length(meanDf)){
      tmpCol <- meanDf[,j]
      tmpCol <- as.numeric(tmpCol[!is.na(tmpCol)])
      triMCol <- litteR::trimean(tmpCol)
      triMeanVec  <- append(triMeanVec, triMCol)
    }
    mDiff <- max(abs(triMeanVec[1] - triMeanVec[-1]))
    geneDf[i,] <- c(tmpEns, 
                    tmpSym, 
                    maxStd,
                    mDiff,
                    meanStd,
                    paste(triMeanVec,collapse=","),
                    paste(colnames(meanDf),collapse=","))
  }
  geneDf <- geneDf[order(geneDf$maxStdGene, decreasing = FALSE),]
  geneDf <- geneDf[!is.na(geneDf$ensembl_gene_id), ]
  return(geneDf)
}

makeExtBetaMean <- function(hGeneBetas, inpPheno, inpH, inpMeth=NULL){
  # Function for getting central tendency for external data cohorts CpG sites 
  if(is.null(inpMeth)){
    inpMeth <- "TRIMEAN"
  }
  tmpB <- hGeneBetas[ ,which(colnames(hGeneBetas) %in% inpPheno$barcode[which(inpPheno$Histotype %in% inpH)])]
  bDf <- data.frame(matrix(nrow=0, ncol=3))
  colnames(bDf) <- c("CpG","Mean", "STD")
  for(j in 1:nrow(tmpB)){
    tmpBr <- tmpB[j,]
    if(toupper(inpMeth) %in% "TRIMEAN"){
      # Get trimean
      meanVal <- makeTm(as.matrix(tmpBr))
    }else if(inpMeth %in% "ARIMEAN"){
      meanVal <- mean(tmpBr)
    }else if(inpMeth %in% "GEOMEAN"){
      meanVal <-exp(mean(log(as.numeric(tmpBr))))
    }
    # Get standard deviation of distribution
    tmpStd <- sd(as.numeric(tmpBr))
    bDf[nrow(bDf)+1,] <- c(rownames(tmpB)[j], meanVal, tmpStd)
  }
  return(bDf)
} 

makeGeneMeanLst <- function(inpBetaLst, inpPhenoLst, inpCpgLocsEPIC, inpCpgLocs450K, inpH, arrTypeLst, methType = NULL){
  # Function for getting trimean beta values for CpG sites in genes for multiple cohorts
  if(is.null(methType)){
    methType <- "TRIMEAN"
  }
  outLst <- list()
  tmpLocs <- inpCpgLocsEPIC[[which(names(inpCpgLocsEPIC) %in% inpH)]]
  for(i in 1:length(tmpLocs)){
    tmpGene <- names(tmpLocs)[i]
    geneBetaDf <- data.frame(matrix(ncol=0, nrow=0))
    for(j in 1:length(inpBetaLst)){
      tmpSid <- names(inpBetaLst)[j]
      tmpB <- inpBetaLst[[tmpSid]]
      tmpPheno <- inpPhenoLst[[tmpSid]]
      arrType <- arrTypeLst[[tmpSid]]
      if(arrType %in% "EPIC"){
        tmpCpgOls <- inpCpgLocsEPIC[which(names(inpCpgLocsEPIC) %in% inpH)][[1]]
      }else{
        tmpCpgOls <- inpCpgLocs450K[which(names(inpCpgLocs450K) %in% inpH)][[1]]
      }
      if(!tmpGene %in% names(tmpCpgOls)){
        next()
      }
      geneOls <- tmpCpgOls[[which(names(tmpCpgOls) %in% tmpGene)]]
      if(nrow(geneBetaDf)==0){
        geneBetaDf <- data.frame(matrix(ncol=0, nrow=nrow(geneOls)))
        geneBetaDf[,ncol(geneBetaDf) + 1] <- rownames(geneOls)
        geneBetaDf[,ncol(geneBetaDf) + 1] <- geneOls$chr
        geneBetaDf[,ncol(geneBetaDf) + 1] <- geneOls$pos
        geneBetaDf[,ncol(geneBetaDf) + 1] <- geneOls$strand
        colnames(geneBetaDf) <- c("CpG", "chr", "pos", "strand")
      }
      geneBetas <- tmpB[which(rownames(tmpB) %in% rownames(geneOls)),]
      geneBetas <- geneBetas[, tmpPheno$barcode[which(tmpPheno$Histotype %in% inpH)]]
      if(nrow(geneBetas) == 0 || ncol(geneBetas) == 0){
        next()
      }
      gBetaDf <- makeExtBetaMean(hGeneBetas = geneBetas, 
                                 inpPheno = tmpPheno, 
                                 inpH = inpH, 
                                 inpMeth = methType)
      colnames(gBetaDf)[2:3] <- c(paste(tmpSid, 
                                        colnames(gBetaDf)[2:3], sep="_"))
      nColInd <- which(!colnames(gBetaDf) %in% colnames(geneBetaDf))
      for(k in 1:length(nColInd)){
        tmpInd <- nColInd[k]
        geneBetaDf[,ncol(geneBetaDf)+1] <- NA
        colnames(geneBetaDf)[ncol(geneBetaDf)] <- colnames(gBetaDf)[tmpInd]
        for(l in 1:nrow(gBetaDf)){
          tmpRow <- gBetaDf[l,]
          tmpRInd <- which(geneBetaDf$CpG %in% tmpRow$CpG)
          geneBetaDf[tmpRInd,ncol(geneBetaDf)] <- gBetaDf[l,tmpInd]
        }
      }
    }
    geneBetaDf <- geneBetaDf[,c(1:4, grep("Mean", colnames(geneBetaDf)))]
    outLst[[i]] <- geneBetaDf
    names(outLst)[i] <- tmpGene
  }
  return(outLst)
}

makePromoBetaRank_MULT <- function(inpRegLocLst, inpPheno, inpGeneInf, allBetas, allM = NULL, focGrp = NULL,
                              minCpg=NULL, sigCpg = NULL, histBool=NULL,
                              totCutThresh=NULL, histCutThresh=NULL,
                              distType=NULL,
                              cpgCutThresh=NULL, cpgCutFreq=NULL,
                              bCut=NULL, pCut = NULL, noCat=NULL,
                              noPBool = NULL, noBBool = NULL,
                              brownBool=NULL, verbBool = NULL,
                              nCores = NULL){
  if(is.null(distType)){
    distType <- "TRIMEAN"
  }
  if(is.null(sigCpg)){
    sigCpg <- 2
  }
  if(is.null(minCpg)){
    minCpg <- 3
  }
  if(is.null(cpgCutFreq)){
    cpgCutFreq <- 0.33
  }
  if(is.null(pCut)){
    pCut <- 0.05
  }
  if(is.null(bCut)){
    bCut <- 0.2
  }
  if(is.null(allM)){
    allM <- log2(allBetas/(1-allBetas))
  }
  outLst <- list()
  tmpPheno <- inpPheno
  tmpPheno$Sample_ID <- tmpPheno$barcode
  phenoGrp <- names(table(tmpPheno$Histotype))
  # Set up the number of categories required for significance
  if(is.null(noCat)){
    noCat <- length(phenoGrp)-1
  }
  if(is.null(focGrp)){
    focGrp <- names(sort(table(tmpPheno$Histotype), decreasing = TRUE))[1]
    message(paste("No focusgroup for analysis submitted, selecting largest group in given cohort: ", focGrp, sep=""))
  }
  # Get sample id's for pheno and non-pheno samp
  phenoIdx <- tmpPheno[which(tmpPheno$Histotype %in% focGrp),]
  nPhenoIdx <- tmpPheno[which(!tmpPheno$Histotype %in% focGrp),]
  nPhenoTab <- names(table(nPhenoIdx$Histotype))
  message(paste("Retrieving regions showing the largest separation seen to beta for histotype: ", focGrp, sep=""))
  refCpgs <- rownames(allBetas)
  nIt <- length(inpRegLocLst)
  resRankLst <- list()
  resRankLst <- future_lapply(1:nIt, function(j){
    # Print out every 10% of genes covered by function
    if(abs(j)%%round(nIt/10) == 0){
      message(paste("Processing gene: [", j, "/",nIt,"]"))
    }
    # Get CpG-Gene dataframe for gene j
    geneBDf <- inpRegLocLst[[j]]
    tmpEns <- names(inpRegLocLst)[j]
    tmpSym <- inpGeneInf$external_gene_name[match(tmpEns, inpGeneInf$ensembl_gene_id)]
    geneBDf <- geneBDf[which(rownames(geneBDf) %in% refCpgs),]
    if(nrow(geneBDf) == 0 | is.null(nrow(geneBDf)) | nrow(geneBDf) < minCpg){
      # if(!is.null(verbBool)){
      #  message(paste("CpGs found to overlap promoter for gene: ", tmpEns, " did not fulfifll requirements", sep=""))
      # }
      resRankLst[[j]] <- NULL
    }else{
      tmpCpgs <- rownames(geneBDf)
      # Perform kruskal-wallis test for significance
      tmpM <- allM[tmpCpgs, ]
      kwDf <- makeKW(inpM = tmpM,
                     inpPheno = tmpPheno,
                     inpH = focGrp)
      # Perform Dunns test for significance in pairwise comparisons (do not correct p-values)
      dunnDf <- makeDunn(inpM = tmpM,
                     inpPheno = tmpPheno,
                     inpH = focGrp)
      dunnDf$Contrast <- gsub(paste("-", focGrp, sep=""), "", dunnDf$Contrast)
      # Perform trimean calculation for effect size
      tmpB <- allBetas[tmpCpgs, ]
      phenoBetas <- tmpB[,phenoIdx$barcode]
      nPhenoBetas <- tmpB[,nPhenoIdx$barcode]
      distDf <- makeBetaGrpDists(phenoBetas = phenoBetas,
                                 nPhenoBetas = nPhenoBetas,
                                 tmpPheno = tmpPheno,
                                 tmpHistotypes = phenoGrp,
                                 tmpName = focGrp,
                                 histBool=TRUE,
                                 distType=distType)
      ebVec <- c()
      for(k in 1:length(nPhenoTab)){
        # Summarize raw p-values using empirical browns for each histotype
        # I.e. get the significance for each comparison over the promoter region
        tmpNH <- nPhenoTab[[k]]
        ebM <- tmpM[,which(colnames(tmpM) %in% tmpPheno$barcode[which(tmpPheno$Histotype %in% c(focGrp, tmpNH))])]
        tmpEb <- EmpiricalBrownsMethod::empiricalBrownsMethod(ebM,
                                                              as.numeric(dunnDf[which(dunnDf$Contrast %in% tmpNH),
                                                                                "pVal"]))
        ebVec <- append(ebVec, 
                        tmpEb)
      }
      names(ebVec) <- nPhenoTab
      # Convert result dataframe to wide format compatible with tm distance format
      dunnDfW <- pivot_wider(dunnDf,
                             id_cols = CpG,
                             names_from = Contrast,
                             values_from = pAdj)
      dunnDfW <- data.frame(dunnDfW)
      # Check to see how many CpGs pass significance thresholds
      passDunn <- dunnDfW[which(rowSums(dunnDfW < pCut) >= noCat),]
      sigPass <- passDunn$CpG
      passTM <- distDf[which(rowSums(distDf > bCut) >= noCat),]
      if(!is.null(bCut)){
        sigPass <- intersect(sigPass, rownames(passTM))
      }
      # Get percentage of passes of CpG's promoter
      passPerc <- length(sigPass)/nrow(dunnDfW)
      passEb <- ebVec[which(ebVec < pCut)]

      if(length(sigPass) == 0 || length(sigPass) < sigCpg || length(passEb) < sigCpg || passPerc < cpgCutThresh){
        resRankLst[[j]] <- NULL
      }else{
        # Summarize the trimean distance of the promoter separation for each histotype
        # For genes that pass the filtering criterion of promoter separation, we amplify their distance
        # With respect to the variation within the CpG probes (i.e. the "consistency")
        histDist <- apply(distDf, 2, makeTm)
        histSum <- round(sum(histDist), 2)
        pSum <- apply(dunnDfW[,2:ncol(dunnDfW)], 2, makeTm)
        nCpg <- length(tmpCpgs)
        brownSum <- sum(ebVec)
        resRankLst[[j]] <- c(tmpEns, tmpSym, nCpg, length(sigPass),
                             histSum, histDist, paste(distDf, sep=", "),
                             pSum, paste(dunnDfW[,2:ncol(dunnDfW)], sep=", "),
                             ebVec, brownSum)
      }
    }
  })
  resRankLst <- resRankLst[lengths(resRankLst) != 0]
  if(length(resRankLst)== 0){
    message(paste("No genes passed filtering criterion"))
    return()
  }else{
    outDf <- do.call(rbind.data.frame, resRankLst)
    # Create result dataframe
    colnames(outDf) <- c("ensembl_gene_id", "external_gene_name", "nCpg","sigCpg","sumDist",
                         nPhenoTab,
                         paste("all_dist",nPhenoTab, sep="_"),
                         paste("sum_adjP",nPhenoTab, sep="_"),
                         paste("all_adjP",nPhenoTab, sep="_"),
                         paste("EmpiricalBrown", nPhenoTab, sep="_"),
                         "brownSum")
    outDf <- outDf[order(outDf[,"brownSum"], decreasing = FALSE),]
  }
  return(outDf)
}

makePromoBetaRank_V3 <- function(inpRegLocLst, inpPheno, inpGeneInf, allBetas, allM = NULL, focGrp = NULL, 
                                 minCpg=NULL, sigCpg = NULL, noCat=NULL,  
                                 bCut=NULL,  pCut = NULL, brownBool=NULL, cpgCutFreq=NULL, 
                                 distType=NULL, verbBool = NULL){
  if(is.null(distType)){
    distType <- "TRIMEAN"
  }
  if(is.null(sigCpg)){
    sigCpg <- 2
  }
  if(is.null(minCpg)){
    minCpg <- 3
  }
  if(is.null(cpgCutFreq)){
    cpgCutFreq <- 0.33
  }
  if(is.null(pCut)){
    pCut <- 0.05
  }
  if(is.null(bCut)){
    bCut <- 0.2
  }
  
  if(is.null(allM)){
    allM <- log2(allBetas/(1-allBetas))
  }
  
  outLst <- list()
  tmpPheno <- inpPheno
  tmpPheno$Sample_ID <- tmpPheno$barcode
  phenoGrps <- names(table(tmpPheno$Histotype))
  # Set up the number of categories required for significance
  if(is.null(noCat)){
    noCat <- length(phenoGrps)-1
  }
  if(is.null(focGrp)){
    focGrp <- names(sort(table(tmpPheno$Histotype), decreasing = TRUE))[1]
    message(paste("No focusgroup for analysis submitted, selecting largest group in given cohort: ", focGrp, sep=""))
  }
  # Get sample id's for pheno and non-pheno samp
  phenoIdx <- tmpPheno[which(tmpPheno$Histotype %in% focGrp),]
  nPhenoIdx <- tmpPheno[which(!tmpPheno$Histotype %in% focGrp),]
  nPhenoTab <- names(table(nPhenoIdx$Histotype))
  # Create result dataframe
  geneDf <- data.frame(matrix(nrow=0, ncol=6+(5*length(nPhenoTab))))
  colnames(geneDf) <- c("ensembl_gene_id", "external_gene_name", "nCpg","sigCpg","sumDist",
                        nPhenoTab,
                        paste("all_dist",nPhenoTab, sep="_"), 
                        paste("sum_adjP",nPhenoTab, sep="_"), 
                        paste("all_adjP",nPhenoTab, sep="_"),
                        paste("EmpiricalBrown", nPhenoTab, sep="_"),
                        "brownSum")
  message(paste("Retrieving regions showing the largest separation seen to beta for histotype: ", focGrp, sep=""))
  # resRankLst <- list()
  refCpgs <- rownames(allBetas)
  nIt <- length(inpRegLocLst)
  for(j in 1:nIt){
    # Print out every 10% of genes covered by function
    if(abs(j)%%round(nIt/10) == 0){
      message(paste("Processing gene: [", j, "/",nIt,"]"))
    }
    # Get CpG-Gene dataframe for gene j
    tmpEns <- names(inpRegLocLst)[j]
    geneBDf <- inpRegLocLst[[j]]
    tmpSym <- inpGeneInf$external_gene_name[which(inpGeneInf$ensembl_gene_id %in% tmpEns)]
    if(is.null(tmpSym)){
      tmpSym <- tmpEns
    }
    tmpCpgs <- rownames(geneBDf)[which(rownames(geneBDf) %in% refCpgs)]
    # geneBDf <- geneBDf[which(rownames(geneBDf) %in% refCpgs),]
    if(length(tmpCpgs) == 0 | length(tmpCpgs) < minCpg){
      next()
    }else{
      # Perform trimean calculation for effect size 
      tmpB <- allBetas[which(rownames(allBetas) %in% tmpCpgs), ]
      # Perform kruskal-wallis test for significance
      tmpM <- allM[which(rownames(allM) %in% tmpCpgs), ]
      kwDf <- makeKW(inpM = tmpM, 
                     inpPheno = tmpPheno, 
                     inpH = focGrp)
      # If we lack enough probes showing significant differences in p, move to the next
      # See which histotypes pass the empirical browns method significance
      passKw <- kwDf[which(as.numeric(kwDf$pVal) < pCut), ]
      if(nrow(passKw) < sigCpg){
          next()
      }else{
        # Perform Dunns test for significance in pairwise comparisons (Obtain raw-values)
        dunnDf <- makeDunn(inpM = tmpM, 
                           inpPheno = tmpPheno, 
                           inpH = focGrp)
        # Obtain adjusted values
        #dunnDfAdj <- makeDunn(inpM = tmpM, 
        #                   inpPheno = tmpPheno, 
        #                   inpH = focGrp,
        #                   inpMethod = "BH")
        phenoBetas <- tmpB[ ,which(colnames(tmpB) %in% phenoIdx$barcode)]
        nPhenoBetas <- tmpB[ ,which(colnames(tmpB) %in% nPhenoIdx$barcode)]
        distDf <- makeBetaGrpDists(phenoBetas = phenoBetas, 
                                   nPhenoBetas = nPhenoBetas, 
                                   tmpPheno = tmpPheno, 
                                   tmpHistotypes = phenoGrps,
                                   tmpName = focGrp, 
                                   histBool=TRUE, 
                                   distType=distType)
        dunnDf$Contrast <- gsub(paste("-", focGrp, sep=""), "", dunnDf$Contrast)
        dunnDfAdj$Contrast <- gsub(paste("-", focGrp, sep=""), "", dunnDfAdj$Contrast)
        # Summarize raw p-values using empirical browns for each histotype
        # I.e. get the significance for each comparison over the promoter region
        ebVec <- c()
        for(k in 1:length(nPhenoTab)){
          tmpNH <- nPhenoTab[[k]]
          ebM <- tmpM[,which(colnames(tmpM) %in% tmpPheno$barcode[which(tmpPheno$Histotype %in% c(focGrp, tmpNH))])]
          tmpEb <- EmpiricalBrownsMethod::empiricalBrownsMethod(ebM, as.numeric(dunnDf[which(dunnDf$Contrast %in% tmpNH), "pVal"]))
          ebVec <- append(ebVec, tmpEb)
        }
        names(ebVec) <- nPhenoTab
        # Convert result dataframe to wide format compatible with tm distance format
        dunnDfW <- pivot_wider(dunnDf, 
                               id_cols = CpG, 
                               names_from = Contrast, 
                               values_from = pAdj)
        dunnDfW <- data.frame(dunnDfW)
        dunnDfW <- column_to_rownames(dunnDfW, "CpG")
        # Check to see how many CpGs pass significance thresholds
        passDunn <- dunnDfW[which(rowSums(dunnDfW < pCut) >= noCat),]
        sigPassCpg <- rownames(passDunn)
        passTM <- distDf[which(rowSums(distDf > bCut) >= noCat),]
        if(!is.null(bCut)){
          sigPassCpg <- intersect(sigPassCpg, rownames(passTM))
        }
        # Significance thresholds
        if(length(sigPassCpg) == 0 || length(sigPassCpg) < sigCpg){
          next()
        }
        # Get percentage of passes of CpG's promoter
        passPerc <- length(sigPassCpg)/nrow(dunnDfW)
        if(!is.null(cpgCutFreq)){
          if(passPerc < cpgCutFreq){
            next()
          }
        }
        # See which histotypes pass the empirical browns method significance
        passEb <- ebVec[which(ebVec < pCut)]
        if(!is.null(brownBool)){
          if(length(passEb) < sigCpg){
            next()
          }
        }
        # Summarize the trimean distance of the promoter separation for each histotype
        histDist <- apply(distDf, 2, makeTm)
        histSum <- round(sum(histDist), 2)
        pSum <- apply(dunnDfW, 2, makeTm)
        nCpg <- length(tmpCpgs)
        brownSum <- sum(ebVec)
        geneDf[nrow(geneDf)+1, ] <- c(tmpEns, tmpSym, nCpg, length(sigPassCpg),
                                      histSum, histDist, paste(distDf, sep=", "),
                                      pSum, paste(dunnDfW, sep=", "),
                                      ebVec, brownSum)
      }
    }
  }
  if(nrow(geneDf) == 0){
    message(paste("No genes passed filtering criterion"))
  }else{
    if(nrow(geneDf) == 0){
      message(paste("No genes passed filtering criterion"))
    }else{
    #   ebCols <- grep("EmpiricalBrown", 
    #                  colnames(geneDf))
    #   for(l in 1:lenght(ebCols)){
    #     geneDf[,ebCols[l]] <- p.adjust(as.numeric(ebCols[l]), n = nrow(geneDf), method = "BH")
    #   }
    #   # Correct EB p-values and filter based on  
    #   if(!is.null(brownBool)){
    #     # See which histotypes pass the empirical browns method significance
    #     passEb <- ebVec[which(ebVec < pCut)]
    #     if(length(passEb) < sigCpg){
    #       next()
    #     }
    #   }
    #   geneDf <- geneDf[order(geneDf[,"brownSum"], decreasing = FALSE),]
    # }
    geneDf <- geneDf[order(geneDf[,"brownSum"], decreasing = FALSE),]
  }
  return(geneDf)
  }
  
makePromoBetaRank_V4 <- function(inpRegLocLst, inpPheno, inpGeneInf, allBetas, allM = NULL, focGrp = NULL, 
                                   minCpg=NULL, sigCpg = NULL, noCat=NULL,  
                                   bCut=NULL,  pCut = NULL, brownBool=NULL, cpgCutFreq=NULL, 
                                   distType=NULL, verbBool = NULL){
    if(is.null(distType)){
      distType <- "TRIMEAN"
    }
    if(is.null(sigCpg)){
      sigCpg <- 2
    }
    if(is.null(minCpg)){
      minCpg <- 3
    }
    if(is.null(cpgCutFreq)){
      cpgCutFreq <- 0.33
    }
    if(is.null(pCut)){
      pCut <- 0.05
    }
    if(is.null(bCut)){
      bCut <- 0.2
    }
    
    if(is.null(allM)){
      allM <- log2(allBetas/(1-allBetas))
    }
    
    outLst <- list()
    tmpPheno <- inpPheno
    tmpPheno$Sample_ID <- tmpPheno$barcode
    phenoGrps <- names(table(tmpPheno$Histotype))
    # Set up the number of categories required for significance
    if(is.null(noCat)){
      noCat <- length(phenoGrps)-1
    }
    if(is.null(focGrp)){
      focGrp <- names(sort(table(tmpPheno$Histotype), decreasing = TRUE))[1]
      message(paste("No focusgroup for analysis submitted, selecting largest group in given cohort: ", focGrp, sep=""))
    }
    # Get sample id's for pheno and non-pheno samp
    phenoIdx <- tmpPheno[which(tmpPheno$Histotype %in% focGrp),]
    nPhenoIdx <- tmpPheno[which(!tmpPheno$Histotype %in% focGrp),]
    nPhenoTab <- names(table(nPhenoIdx$Histotype))
    # Create result dataframe
    geneDf <- data.frame(matrix(nrow=0, ncol=6+(5*length(nPhenoTab))))
    colnames(geneDf) <- c("ensembl_gene_id", "external_gene_name", "nCpg","sigCpg","sumDist",
                          nPhenoTab,
                          paste("all_dist",nPhenoTab, sep="_"), 
                          paste("sum_adjP",nPhenoTab, sep="_"), 
                          paste("all_adjP",nPhenoTab, sep="_"),
                          paste("EmpiricalBrown", nPhenoTab, sep="_"),
                          "brownSum")
    message(paste("Retrieving regions showing the largest separation seen to beta for histotype: ", focGrp, sep=""))
    # resRankLst <- list()
    refCpgs <- rownames(allBetas)
    nIt <- length(inpRegLocLst)
    for(j in 1:nIt){
      # Print out every 10% of genes covered by function
      if(abs(j)%%round(nIt/10) == 0){
        message(paste("Processing gene: [", j, "/",nIt,"]"))
      }
      # Get CpG-Gene dataframe for gene j
      tmpEns <- names(inpRegLocLst)[j]
      geneBDf <- inpRegLocLst[[j]]
      tmpSym <- inpGeneInf$external_gene_name[which(inpGeneInf$ensembl_gene_id %in% tmpEns)]
      if(is.null(tmpSym)){
        tmpSym <- tmpEns
      }
      tmpCpgs <- rownames(geneBDf)[which(rownames(geneBDf) %in% refCpgs)]
      # geneBDf <- geneBDf[which(rownames(geneBDf) %in% refCpgs),]
      if(length(tmpCpgs) == 0 | length(tmpCpgs) < minCpg){
        next()
      }else{
        # Perform trimean calculation for effect size 
        tmpB <- allBetas[which(rownames(allBetas) %in% tmpCpgs), ]
        # Perform kruskal-wallis test for significance
        tmpM <- allM[which(rownames(allM) %in% tmpCpgs), ]
        kwDf <- makeKW(inpM = tmpM, 
                       inpPheno = tmpPheno, 
                       inpH = focGrp)
        # If we lack enough probes showing significant differences in p, move to the next
        # See which histotypes pass the empirical browns method significance
        passKw <- kwDf[which(as.numeric(kwDf$pVal) < pCut), ]
        if(nrow(passKw) < sigCpg){
          next()
        }else{
          # Perform Dunns test for significance in pairwise comparisons (Obtain raw-values)
          dunnDf <- makeDunn(inpM = tmpM, 
                             inpPheno = tmpPheno, 
                             inpH = focGrp)
          # Obtain adjusted values
          dunnDfAdj <- makeDunn(inpM = tmpM, 
                             inpPheno = tmpPheno, 
                             inpH = focGrp,
                             inpMethod = "BH")
          phenoBetas <- tmpB[ ,which(colnames(tmpB) %in% phenoIdx$barcode)]
          nPhenoBetas <- tmpB[ ,which(colnames(tmpB) %in% nPhenoIdx$barcode)]
          distDf <- makeBetaGrpDists(phenoBetas = phenoBetas, 
                                     nPhenoBetas = nPhenoBetas, 
                                     tmpPheno = tmpPheno, 
                                     tmpHistotypes = phenoGrps,
                                     tmpName = focGrp, 
                                     histBool=TRUE, 
                                     distType=distType)
          dunnDf$Contrast <- gsub(paste("-", focGrp, sep=""), "", dunnDf$Contrast)
          dunnDfAdj$Contrast <- gsub(paste("-", focGrp, sep=""), "", dunnDfAdj$Contrast)
          # Summarize raw p-values using empirical browns for each histotype
          # I.e. get the significance for each comparison over the promoter region
          ebVec <- c()
          for(k in 1:length(nPhenoTab)){
            tmpNH <- nPhenoTab[[k]]
            ebM <- tmpM[,which(colnames(tmpM) %in% tmpPheno$barcode[which(tmpPheno$Histotype %in% c(focGrp, tmpNH))])]
            tmpEb <- EmpiricalBrownsMethod::empiricalBrownsMethod(ebM, as.numeric(dunnDf[which(dunnDf$Contrast %in% tmpNH), "pVal"]))
            ebVec <- append(ebVec, tmpEb)
          }
          names(ebVec) <- nPhenoTab
          # Convert result dataframe to wide format compatible with tm distance format
          dunnDfW <- pivot_wider(dunnDfAdj, 
                                 id_cols = CpG, 
                                 names_from = Contrast, 
                                 values_from = pAdj)
          dunnDfW <- data.frame(dunnDfW)
          dunnDfW <- column_to_rownames(dunnDfW, "CpG")
          # Check to see how many CpGs pass significance thresholds
          passDunn <- dunnDfW[which(rowSums(dunnDfW < pCut) >= noCat),]
          sigPassCpg <- rownames(passDunn)
          passTM <- distDf[which(rowSums(distDf > bCut) >= noCat),]
          if(!is.null(bCut)){
            sigPassCpg <- intersect(sigPassCpg, rownames(passTM))
          }
          # Significance thresholds
          if(length(sigPassCpg) == 0 || length(sigPassCpg) < sigCpg){
            next()
          }
          # Get percentage of passes of CpG's promoter
          passPerc <- length(sigPassCpg)/nrow(dunnDfW)
          if(!is.null(cpgCutFreq)){
            if(passPerc < cpgCutFreq){
              next()
            }
          }
          # See which histotypes pass the empirical browns method significance
          passEb <- ebVec[which(ebVec < pCut)]
          if(!is.null(brownBool)){
            if(length(passEb) < sigCpg){
              next()
            }
          }
          # Summarize the trimean distance of the promoter separation for each histotype
          histDist <- apply(distDf, 2, makeTm)
          histSum <- round(sum(histDist), 2)
          pSum <- apply(dunnDfW, 2, makeTm)
          nCpg <- length(tmpCpgs)
          brownSum <- sum(ebVec)
          geneDf[nrow(geneDf)+1, ] <- c(tmpEns, tmpSym, nCpg, length(sigPassCpg),
                                        histSum, histDist, paste(distDf, sep=", "),
                                        pSum, paste(dunnDfW, sep=", "),
                                        ebVec, brownSum)
        }
      }
    }
    if(nrow(geneDf) == 0){
      message(paste("No genes passed filtering criterion"))
    }else{
      if(nrow(geneDf) == 0){
        message(paste("No genes passed filtering criterion"))
      }else{
        ebCols <- grep("EmpiricalBrown", 
                          colnames(geneDf))
        for(l in 1:lenght(ebCols)){
             geneDf[,ebCols[l]] <- p.adjust(as.numeric(ebCols[l]), n = nrow(geneDf), method = "BH")
        }
        # Correct EB p-values and filter based on  
        if(!is.null(brownBool)){
            # See which histotypes pass the empirical browns method significance
          passEb <- ebVec[which(ebVec < pCut)]
          if(length(passEb) < sigCpg){
            next()
          }
        }
        geneDf <- geneDf[order(geneDf[,"brownSum"], decreasing = FALSE),]
      }
      geneDf <- geneDf[order(geneDf[,"brownSum"], decreasing = FALSE),]
    }
  }
  return(geneDf)
}

makePromoBetaRank <- function(inpBetas, inpPheno, inpGeneInf, allBetas, focGrp = NULL, 
                                allM = NULL, minCpg=NULL, sigCpg = NULL, histBool=NULL, 
                                totCutThresh=NULL, histCutThresh=NULL, 
                                distType=NULL, sumType=NULL, cpgCutThresh=NULL,
                                bCut=NULL, pCut = NULL, noCat=NULL,
                                noPBool = NULL, noBBool = NULL, 
                                brownBool=NULL, verbBool = NULL){
  if(is.null(distType)){
    distType <- "TRIMEAN"
  }
  if(is.null(sigCpg)){
    sigCpg <- 2
  }
  if(is.null(minCpg)){
    minCpg <- 3
  }
  if(is.null(pCut)){
    pCut <- 0.05
  }
  if(is.null(bCut)){
    bCut <- 0.2
  }
  outLst <- list()
  tmpPheno <- inpPheno
  phenoGrp <- names(table(tmpPheno$Histotype))
  if(is.null(focGrp)){
    focGrp <- names(sort(table(tmpPheno$Histotype), decreasing = TRUE))[1]
    message(paste("No focusgroup for analysis submitted, selecting largest group in given cohort: ", focGrp, sep=""))
  }
  if(!focGrp %in% phenoGrp){
    message(paste("ERROR! The given focusgroup: ", focGrp, "does not exist within cohort data, please select another primary phenotype"))
    return(NULL)
  }
  tmpGenes <- names(inpBetas)
  # Get sample id's for pheno and non-pheno samp
  phenoIdx <- tmpPheno[which(tmpPheno$Histotype %in% focGrp),]
  nPhenoIdx <- tmpPheno[which(!tmpPheno$Histotype %in% focGrp),]
  nPhenoTab <- table(nPhenoIdx$Histotype)
  nPhenoHists <- names(nPhenoTab)
  # Create result dataframe
  geneDf <- data.frame(matrix(nrow=0, ncol=8+(length(nPhenoHists)*2)+(length(phenoGrp)*2)+(length(phenoGrp)-1)))
  colnames(geneDf) <- c("ensembl_gene_id", "external_gene_name", "nCpg","sigCpg","sumDist","brownSum",
                        nPhenoHists, 
                        paste("all_dist", nPhenoHists,sep="_"), 
                        paste("var",phenoGrp,sep="_"), 
                        "varSum", "varSumDist", "norm", 
                        paste("sum_adjP",nPhenoHists,sep="_"), 
                        paste("all_adjP",nPhenoHists,sep="_"))
  message(paste("Retrieving regions showing the largest separation seen to beta for histotype: ", focGrp, sep=""))
  for(j in 1:length(inpBetas)){
    # Print out every 10% of genes covered by function
    if(abs(j)%%round(length(inpBetas)/10) == 0){
      message(paste("Processing gene: [", j, "/",length(inpBetas),"]"))
    }
    tmpEns <- names(inpBetas)[j]
    tmpSym <- inpGeneInf$external_gene_name[match(tmpEns, inpGeneInf$ensembl_gene_id)]
    geneBDf <- inpBetas[[j]]
    if(nrow(geneBDf) == 0 | is.null(nrow(geneBDf)) | nrow(geneBDf) < minCpg){
      if(!is.null(verbBool)){
        message(paste("CpGs found to overlap promoter for gene: ", tmpEns, " did not fulfifll requirements", sep=""))
      }
      next()
    }
    tmpGeneRow <- inpGeneInf[inpGeneInf$ensembl_gene_id %in% tmpEns,]
    tmpSym <- tmpGeneRow$external_gene_name
    if(is.na(tmpSym) || tmpSym == "" || is.null(tmpSym)){
      tmpSym <- tmpEns
    }
    # Retrieve Beta values for given promoter region
    cpgDf <- allBetas[which(rownames(allBetas) %in% rownames(geneBDf)),]
    # Remove potential NA rows
    cpgDf <- na.omit(cpgDf)
    if(nrow(cpgDf) <= 1){
      if(!is.null(verbBool)){
        message("Gene: ", tmpEns, "has n<=1 CpGs and will be skipped")
      }
      next()
    }
    # If m-values are supplied, we use them instead of beta-values when testing for significance
    if(!is.null(allM)){
      mDf <- allM[rownames(cpgDf),]
    }else{
      mDf <- cpgDf
    }
    # Get beta-values for phenotype of interest, and all others
    phenoBetas <- cpgDf[,colnames(cpgDf) %in% phenoIdx$barcode]
    colnames(phenoBetas) <- phenoIdx$Sample_ID[match(colnames(phenoBetas), phenoIdx$barcode)]
    nPhenoBetas <- cpgDf[,colnames(cpgDf) %in% nPhenoIdx$barcode]
    colnames(nPhenoBetas) <- nPhenoIdx$Sample_ID[match(colnames(nPhenoBetas), nPhenoIdx$barcode)]
    phenoBetas <- phenoBetas[,order(colnames(phenoBetas))]
    nPhenoBetas <- nPhenoBetas[,order(colnames(nPhenoBetas))]
    # Get two dataframes consisting of the distance between the beta values for histotypes, and one with the internal variance of probes for histotypes
    # distDf <- makeBetaGrpDists(phenoBetas = phenoBetas, nPhenoBetas = nPhenoBetas, 
    #                            tmpPheno=tmpPheno, inpNFoc = phenoGrp,inpFoc = focGrp, 
    #                            histBool=TRUE, distType=distType)
    distDf <- makeBetaGrpDists(phenoBetas = phenoBetas, 
                               nPhenoBetas = nPhenoBetas, 
                               tmpPheno=tmpPheno,
                               tmpHistotypes = phenoGrp,
                               tmpName = focGrp, 
                               histBool=TRUE, 
                               distType=distType)
    varDf <- makeBetaGrpVar(cpgDf= cpgDf, 
                            tmpHistotypes = phenoGrp, tmpPheno = tmpPheno,
                            distType=distType)
    # Create vector to store Kruskal-Wallis results, and dataframe for Dunn's test for significance
    kwVec <- c()
    dunnProbeDf <- data.frame(matrix(nrow=nrow(cpgDf), ncol=length(nPhenoTab)))
    rownames(dunnProbeDf) <- rownames(cpgDf)
    colnames(dunnProbeDf) <- names(nPhenoTab)
    # Create dataframe to store uncorrected p-values
    dunnRawDf <- dunnProbeDf
    tmpH <- inpPheno$Histotype[match(colnames(mDf), inpPheno$barcode)]
    for(m in 1:nrow(mDf)){
      # Perform Kruskal-Wallis Multiple Comparisons test for significance between groups
      # Utilize M-values to avoid the issue of heteroscedasticity found in beta-values
      tmpC <- mDf[m,]
      cwDf <- rbind(tmpC, tmpH)
      cwDf <- data.frame(t(cwDf))
      colnames(cwDf) <- c("M", "Histotype")
      # Set levels (first level is reference)
      cwDf$Histotype <- factor(cwDf$Histotype, levels = c(focGrp, names(nPhenoTab)))
      cwTest <- kruskal.test(x=as.numeric(cwDf$M) , g=cwDf$Histotype)
      kwVec <- append(kwVec, cwTest$p.value)
      # Perform Dunn's test to check significance in pairwise comparisons, adjusted for multiple comparisons (BH)
      dunnRes <- DunnTest(x=as.numeric(cwDf$M), 
                          g=cwDf$Histotype,
                          alternative = "two.sided",
                          method = "BH") 
      tmpDunn <- data.frame(dunnRes[[1]])
      dunnRows <- sapply(rownames(tmpDunn), function(x) strsplit(x, "-")[[1]][[2]], USE.NAMES=FALSE)
      keepDunn <- tmpDunn[which(dunnRows %in% focGrp),]
      dunnProbeDf[m, ] <-  keepDunn$pval
      # Repeat for uncorrected p-values
      dunnRaw <- DunnTest(x=as.numeric(cwDf$M),
                          g=cwDf$Histotype,
                          alternative = "two.sided",
                          method = "none")
      tmpRaw <- data.frame(dunnRaw[[1]])
      dunnRawRows <- sapply(rownames(tmpRaw), function(x) strsplit(x, "-")[[1]][[2]], USE.NAMES=FALSE)
      keepRaw <- tmpRaw[which(dunnRawRows %in% focGrp),]
      dunnRawDf[m, ] <-  keepRaw$pval
    }
    # Use Empirical Browns method to aggregate p-values, requires uncorrected p-values
    # We use EB's method as CpG-sites in promoter regions experience comethylation, and thus are dependent
    # Perform for Kruskal Wallis (unadjusted)
    kwBrown <- EmpiricalBrownsMethod::empiricalBrownsMethod(cpgDf, kwVec)
    # Perform for unadjusted Dunn
    brownVec <- c() 
    for(o in 1:ncol(dunnRawDf)){
      tmpC <- dunnRawDf[,o]
      tmpH <- colnames(dunnRawDf)[o]
      # Get samples used to generate p-values (i.e. the focus-phenotype and the other group in the pairwise comparison)
      tmpSamps <- cpgDf[,inpPheno$barcode[which(inpPheno$Histotype %in% c(focGrp, tmpH))]]
      brownVal <- EmpiricalBrownsMethod::empiricalBrownsMethod(tmpSamps, tmpC)
      brownVec <- append(brownVec, brownVal)
      names(brownVec)[length(brownVec)] <- tmpH
    }
    # Set up the number of categories required for significance
    if(is.null(noCat)){
      noCat <- ncol(distDf)-1
    }
    # Filter for significance using browns method (i.e. for the entire promoter region)
    if(!is.null(brownBool)){
      if(!length(which(brownVec < pCut)) >= noCat){
        next()
      }
    }
    
    # Filter for significance using adjusted p-values from Dunn's method
    passKw <- dunnProbeDf[which(rowSums(dunnProbeDf < pCut) >= noCat),]
    sigDistDf <- distDf[rownames(distDf) %in% rownames(passKw),]
    # Filter for effect size based on CpGs passing significance
    sigDf <- sigDistDf[which(rowSums(sigDistDf > bCut) >= noCat),]
    if(nrow(sigDf) < sigCpg){
      next()
    }
    # Get percentage of passes of CpG's promoter
    passPerc <- nrow(sigDf)/nrow(dunnProbeDf)
    if(passPerc < cpgCutThresh){
      next()
    }
    
    # If criteria passed, get the norm of each gene-row
    normVec <- c()
    for(l in 1:nrow(distDf)){
      tmpRow <- distDf[l,]
      normVec <- append(normVec, norm(tmpRow, type="2") )
    }
    # Summarize the trimean distance of the promoter separation for each histotype
    # For genes that pass the filtering criterion of promoter separation, we amplify their distance
    # With respect to the variation within the CpG probes (i.e. the "consistency")
    histDist <- c()
    for(c in 1:(ncol(distDf))){
      tmpCol <- distDf[,c]
      # Get the  geometric mean of the promoter methylation
      tM <- litteR::trimean(tmpCol)
      histDist <- append(histDist, tM)
      names(histDist)[c] <- colnames(distDf)[c]
    }
    
    varHist <- c()
    for(d in 1:ncol(varDf)){
      tmpVCol <- varDf[,d]
      geomVM <- litteR::trimean(tmpVCol)
      varHist <- append(varHist, geomVM)
      names(varHist)[d] <- colnames(varDf)[d]
    }
    norm <- litteR::trimean(normVec)
    # Multiply (i.e. divide by) by the internal distance between points
    # This gives us a score for distance from pheno to n-pheno, with respect to separation between histotypes
    histSum <- sum(histDist)
    varSum <-  sum(varHist)
    filtDist <- round(histSum/varSum,2)
    histDist <- round(histDist,2)
    histSum <- round(histSum,2)
    varSum <-  round(varSum,2)
    brownSum <- sum(brownVec)
    tmpDists <- distDf
    #geneDf[nrow(geneDf)+1,] <- rep(0, ncol(geneDf))
    #geneDf[nrow(geneDf), "ensembl_gene_id"] <- tmpEns
    #geneDf[nrow(geneDf), "external_gene_name"] <- tmpSym
    #geneDf[nrow(geneDf), "nCpg"] <- nrow(distDf)
    #geneDf[nrow(geneDf), "sigCpg"] <- nrow(sigDf)
    #geneDf[nrow(geneDf), "sumDist"] <- histSum
    #geneDf[nrow(geneDf), "brownSum"] <- brownSum
    # geneDf[nrow(geneDf), phenoGrp[!phenoGrp %in% focGrp],sep="_") <- histDist
    geneDf[nrow(geneDf)+1,] <- c(tmpEns, tmpSym, nrow(distDf), nrow(sigDf),
                                 histSum, brownSum, 
                                 histDist, 
                                 paste(tmpDists, sep=", "),
                                 varHist, varSum, filtDist, 
                                 norm,  
                                 brownVec, paste(dunnProbeDf, sep=", "))
  }
  
  if(nrow(geneDf) == 0){
    message(paste("No genes passed filtering criterion"))
  }else{
    geneDf <- geneDf[order(geneDf[,"brownSum"], decreasing = FALSE),]
  }
  return(geneDf)
}

makePromoCovDf <- function(inpDf, inpCpgs){
  allGenes <- inpCpgs$Gene
  promoCov_cpg_cols <- inpCpgs[,grep("cpg",colnames(inpCpgs))]
  covVec <- c()
  nIt <- nrow(promoCov_cpg_cols)
  for(i in 1:nIt){
    if(abs(i)%%round(nIt/10) == 0 || is.na(abs(i)%%round(nIt/10))){
      message(paste("Processing row: [", i, "/", nIt,"]"))
    }
    tmpCov <- length(which(promoCov_cpg_cols[i,] %in% rownames(inpDf)))
    covVec <- append(covVec, tmpCov)
    names(covVec)[length(covVec)] <- allGenes[i]
  }
  cov_0 <- length(which(covVec == 0))
  cov_1 <- length(which(covVec == 1))
  cov_2 <- length(which(covVec == 2))
  cov_3 <- length(which(covVec == 3))
  cov_4 <- length(which(covVec == 4))
  cov_5 <- length(which(covVec == 5))
  cov_1m <- length(which(covVec >= 1 & covVec < 3))
  cov_3m <- length(which(covVec >= 3 & covVec < 5))
  cov_5m <- length(which(covVec >= 5))
  return(c("ALL" = nrow(promoCov_cpg_cols),"0"=cov_0,"1"=cov_1, "2"=cov_2,"3"=cov_3, "4"=cov_4, "5"=cov_5,"1+"=cov_1m, "3+"=cov_3m, "5+"=cov_5m))
}

makePromoCovDf_MULT <- function(inpCpgs, inpInf, inpDf = NULL){
  allGenes <- inpCpgs$Gene
  promoCov_cpg_cols <- inpCpgs[,grep("cpg",colnames(inpCpgs))]
  nIt <- nrow(promoCov_cpg_cols)
  covLst <- list()
  if(!is.null(inpDf)){
    allRef <- rownames(inpDf)
  }
  covLst <- future_lapply(1:nIt, function(i){
    if(abs(i)%%round(nIt/10) == 0 || is.na(abs(i)%%round(nIt/10))){
      message(paste("Processing row: [", i, "/", nIt,"]"))
    }
    pCovRow <- promoCov_cpg_cols[i,]
    if(!is.null(inpDf)){
      pCovRow <- pCovRow[which(pCovRow %in% allRef)]
    }
    nCov <- length(pCovRow)
    covLst[[i]] <- nCov
  })
  covDf <- do.call(rbind.data.frame, covLst)
  colnames(covDf) <- "Cov"
  cov_0 <- length(which(covDf$Cov == 0))
  cov_1 <- length(which(covDf$Cov == 1))
  cov_2 <- length(which(covDf$Cov == 2))
  cov_3 <- length(which(covDf$Cov == 3))
  cov_4 <- length(which(covDf$Cov == 4))
  cov_5 <- length(which(covDf$Cov >= 5))
  cov_1_3m <- length(which(covDf$Cov >= 1 & covDf$Cov <= 3))
  cov_4_6m <- length(which(covDf$Cov > 3 & covDf$Cov <= 6))
  cov_7_10m <- length(which(covDf$Cov > 6 & covDf$Cov <= 10))
  cov_10plus <- length(which(covDf$Cov > 10))
  resVec <- c("ALL" = nrow(promoCov_cpg_cols),
              "0"=cov_0,"1"=cov_1, "2"=cov_2,"3"=cov_3, "4"=cov_4, "5"=cov_5,
              "1-3"=cov_1_3m, "4-6"=cov_4_6m, "7-10"=cov_7_10m, "10+"=cov_10plus)
  return(resVec)
}

makeCpGDistTypeDf <-function(inpBeta, inpPheno){
  # Test cpg sites in input beta matrix, for each histotype in the input phenotype data for normality, skewness and distribution type
  cpgTypeDf <- data.frame(matrix(nrow=0, ncol=7))
  colnames(cpgTypeDf) <- c("ShapNormal","ShapNotNormal", "LightSkew", "HeavySkew", "LeptoKurtic",  "PlatyKurtic", "MesoKurtic")
  for(i in 1:length(names(table(inpPheno$Histotype)))){
    tmpH <- names(table(inpPheno$Histotype))[i]
    message(paste("Retrieving distribution type data for phenotype:", tmpH, sep=" "))
    cpgTypeDf[nrow(cpgTypeDf)+1, ] <- c(0,0,0,0,0,0,0)
    rownames(cpgTypeDf)[nrow(cpgTypeDf)] <- tmpH
    tmpB <- inpBeta[,which(colnames(inpBeta) %in% inpPheno$barcode[which(inpPheno$Histotype %in% tmpH)])]
    nIt <- nrow(tmpB)
    for(j in 1:nIt){
      if(abs(j)%%round(nIt/10) == 0){
        message(paste("Processing CpG site: [", j, "/", nIt,"]"))
      }
      tmpCrow <- tmpB[j,]
      # If p-value is lower then 0.05, data significantly deviates from normal dist
      normRes <- shapiro.test(as.numeric(tmpCrow))
      if(normRes$p.value > 0.05){
        cpgTypeDf[tmpH,"ShapNormal"] <- cpgTypeDf[tmpH,"ShapNormal"] + 1
      }else{
        cpgTypeDf[tmpH,"ShapNotNormal"] <- cpgTypeDf[tmpH,"ShapNotNormal"] + 1
      }
      # Skewness between -1, 1 is seen as normal
      # Skewness between -2, 2 is seen as lightly skewed
      # Skewness beyond this is substantially different
      skewRes <- e1071::skewness(as.numeric(tmpCrow))
      if(abs(skewRes) >= 1 &  abs(skewRes) < 2){
        cpgTypeDf[tmpH,"LightSkew"] <- cpgTypeDf[tmpH,"LightSkew"] + 1
      }else if(abs(skewRes) >= 2){
        cpgTypeDf[tmpH,"HeavySkew"] <- cpgTypeDf[tmpH,"HeavySkew"] + 1
      }
      # e1071 uses type 3 = excess kurtosis where we change boundary values by a factor of -3
      # That is to say, instead of 3 representing normal distribution, 0 is normal distribution
      kurtRes <- e1071::kurtosis(as.numeric(tmpCrow), type = 3)
      if(kurtRes > 0){
        cpgTypeDf[tmpH,"LeptoKurtic"] <- cpgTypeDf[tmpH,"LeptoKurtic"] + 1
      }else if(kurtRes < 0){
        cpgTypeDf[tmpH,"PlatyKurtic"] <- cpgTypeDf[tmpH,"PlatyKurtic"] + 1
      }else{
        cpgTypeDf[tmpH,"MesoKurtic"] <- cpgTypeDf[tmpH,"MesoKurtic"] + 1
      }
    }
  }
  return(cpgTypeDf)
}

makeCpGDistTypeDf_MULT <-function(inpBeta, inpPheno, nCores = NULL){
  # Start multicore session
  if(is.null(nCores)){
    nCores <- 1
  }
  # plan(multisession, workers = nCores+2)
  
  # Test cpg sites in input beta matrix, for each histotype in the input phenotype data for normality, skewness and distribution type
  cpgTypeDf <- data.frame(matrix(nrow=0, ncol=7))
  colnames(cpgTypeDf) <- c("ShapNormal","ShapNotNormal", "LightSkew", "HeavySkew", "LeptoKurtic",  "PlatyKurtic", "MesoKurtic")
  for(i in 1:length(names(table(inpPheno$Histotype)))){
    tmpH <- names(table(inpPheno$Histotype))[i]
    message(paste("Retrieving distribution type data for phenotype:", tmpH, sep=" "))
    cpgTypeDf[nrow(cpgTypeDf)+1, ] <- c(0,0,0,0,0,0,0)
    rownames(cpgTypeDf)[nrow(cpgTypeDf)] <- tmpH
    tmpB <- inpBeta[,which(colnames(inpBeta) %in% inpPheno$barcode[which(inpPheno$Histotype %in% tmpH)])]
    nIt <- nrow(tmpB)
    cpgDistLst <- list()
    cpgDistLst <- future_lapply(1:nIt, function(j){
      
      if(abs(j)%%round(nIt/10) == 0 || is.na(abs(j)%%round(nIt/10))){
        message(paste("Processing gene: [", j, "/", nIt,"]"))
      }
      tmpCrow <- tmpB[j,]
      # If p-value is lower then 0.05, data significantly deviates from normal dist
      normRes <- shapiro.test(as.numeric(tmpCrow))
      if(normRes$p.value > 0.05){
        shapTest <- "ShapNormal"
      }else{
        shapTest <- "ShapNotNormal"
      }
      # Skewness between -1, 1 is seen as normal
      # Skewness between -2, 2 is seen as lightly skewed
      # Skewness beyond this is substantially different
      skewRes <- e1071::skewness(as.numeric(tmpCrow))
      if(abs(skewRes) >= 1 &  abs(skewRes) < 2){
        skewRes <- "LightSkew"
      }else if(abs(skewRes) >= 2){
        skewRes <- "HeavySkew"
      }else{
        skewRes <- "NoSkew"
      }
      # e1071 uses type 3 = excess kurtosis where we change boundary values by a factor of -3
      # That is to say, instead of 3 representing normal distribution, 0 is normal distribution
      kurtRes <- e1071::kurtosis(as.numeric(tmpCrow), type = 3)
      if(kurtRes > 0){
        kurtRes <- "LeptoKurtic"
      }else if(kurtRes < 0){
        kurtRes <- "PlatyKurtic"
      }else{
        kurtRes <- "MesoKurtic"
      }
      cpgDistLst[[j]] <- c(shapTest, skewRes, kurtRes) 
    })
    cpgDistDf <- do.call(rbind.data.frame, cpgDistLst)
    rownames(cpgDistDf) <- rownames(tmpB)
    colnames(cpgDistDf) <- c("ShapTest", "SkewTest", "KurtTest")
    cpgTypeDf[tmpH, "ShapNormal"] <- sum(cpgDistDf$ShapTest %in%"ShapNormal") 
    cpgTypeDf[tmpH, "ShapNotNormal"] <- sum(cpgDistDf$ShapTest %in%"ShapNotNormal")
    cpgTypeDf[tmpH, "LightSkew"] <- sum(cpgDistDf$SkewTest %in% "LightSkew")
    cpgTypeDf[tmpH, "HeavySkew"] <- sum(cpgDistDf$SkewTest %in% "HeavySkew")
    cpgTypeDf[tmpH, "LeptoKurtic"] <- sum(cpgDistDf$KurtTest %in% "LeptoKurtic")
    cpgTypeDf[tmpH, "PlatyKurtic"] <- sum(cpgDistDf$KurtTest %in% "PlatyKurtic")
    cpgTypeDf[tmpH, "MesoKurtic"] <- sum(cpgDistDf$KurtTest %in% "MesoKurtic")
  }
  return(cpgTypeDf)
}

makeBetaDfMethType <- function(inpBeta, inpPheno){
  # Test the trimean distribution for each histotype
  cpgBetaDf <- data.frame(matrix(nrow=0, ncol=6))
  colnames(cpgBetaDf) <- c("Hypo","Hemi", "Hyper", "Hypo_perc","Hemi_perc", "Hyper_perc")
  for(i in 1:length(names(table(inpPheno$Histotype)))){
    tmpH <- names(table(inpPheno$Histotype))[i]
    cpgBetaDf[nrow(cpgBetaDf)+1, ] <- c(0,0,0,0,0,0)
    rownames(cpgBetaDf)[nrow(cpgBetaDf)] <- tmpH
    tmpB <- inpBeta[,which(colnames(inpBeta) %in% inpPheno$barcode[which(inpPheno$Histotype %in% tmpH)])]
    message(paste("Calculating CpG methylation type for phenotype: ", tmpH, sep=""))
    for(j in 1:nrow(tmpB)){
      if(abs(j)%%round(nrow(tmpB)/10) == 0 || is.na(abs(j)%%round(nrow(tmpB)/10))){
        message(paste("Processing CpG: [", j, "/", nrow(tmpB),"]"))
      }
      tmpCrow <- tmpB[j,]
      for(k in 1:ncol(tmpCrow)){
        if(tmpCrow[,k] < 0.3){
          cpgBetaDf[tmpH,"Hypo"] <- cpgBetaDf[tmpH,"Hypo"] + 1
        }else if(tmpCrow[,k] > 0.3 & tmpCrow[,k] < 0.7){
          cpgBetaDf[tmpH,"Hemi"] <- cpgBetaDf[tmpH,"Hemi"] + 1
        }else if(tmpCrow[,k] > 0.7){
          cpgBetaDf[tmpH,"Hyper"] <- cpgBetaDf[tmpH,"Hyper"] + 1
        }
      }
    }
  }
  for(l in 1:length(table(inpPheno$Histotype))){
    tmpH2 <- names(table(inpPheno$Histotype))[l]
    cpgBetaDf[tmpH2, "Hypo_perc"] <- cpgBetaDf[tmpH2, "Hypo"]/(nrow(tmpB)*length(which(inpPheno$Histotype %in% tmpH2)))
    cpgBetaDf[tmpH2, "Hemi_perc"] <- cpgBetaDf[tmpH2, "Hemi"]/(nrow(tmpB)*length(which(inpPheno$Histotype %in% tmpH2)))
    cpgBetaDf[tmpH2, "Hyper_perc"] <- cpgBetaDf[tmpH2, "Hyper"]/(nrow(tmpB)*length(which(inpPheno$Histotype %in% tmpH2)))
  }
  return(cpgBetaDf)
}

makeBetaDfMethType_MULT <- function(inpBeta, inpPheno, nCores = NULL){
  # Start multicore session
  if(is.null(nCores)){
    nCores <- 1
  }
  # Test the trimean distribution for each histotype
  cpgBetaDf <- data.frame(matrix(nrow=0, ncol=6))
  colnames(cpgBetaDf) <- c("Hypo","Hemi", "Hyper", "Hypo_perc","Hemi_perc", "Hyper_perc")
  for(i in 1:length(names(table(inpPheno$Histotype)))){
    tmpH <- names(table(inpPheno$Histotype))[i]
    cpgBetaDf[nrow(cpgBetaDf)+1, ] <- c(0,0,0,0,0,0)
    rownames(cpgBetaDf)[nrow(cpgBetaDf)] <- tmpH
    tmpB <- inpBeta[,which(colnames(inpBeta) %in% inpPheno$barcode[which(inpPheno$Histotype %in% tmpH)])]
    message(paste("Calculating CpG methylation type for phenotype: ", tmpH, sep=""))
    nIt <- nrow(tmpB)
    countRes <- list()
    countRes <- future_lapply(1:nIt, function(j){
      if(abs(j)%%round(nIt/10) == 0 || is.na(abs(j)%%round(nIt/10))){
        message(paste("Processing CpG: [", j, "/", nIt,"]"))
      }
      tmpCrow <- tmpB[j,]
      nHemi <- length(which(tmpCrow > 0.3 & tmpCrow < 0.7))
      nHypo <- length(which(tmpCrow < 0.3))
      nHyper <- length(which(tmpCrow > 0.7))
      # countRes[nrow(countRes)+1,]
      countRes[[j]] <- c(nHypo, nHemi, nHyper)
    })
    countDf <- do.call(rbind.data.frame, countRes)
    # Add stats to DF for given histotype
    cpgBetaDf[tmpH,"Hypo"] <- sum(countDf[,1])
    cpgBetaDf[tmpH,"Hemi"] <- sum(countDf[,2])
    cpgBetaDf[tmpH,"Hyper"] <- sum(countDf[,3])
    totCounts <- nIt*length(which(inpPheno$Histotype %in% tmpH))
    cpgBetaDf[tmpH, "Hypo_perc"] <- cpgBetaDf[tmpH, "Hypo"]/totCounts
    cpgBetaDf[tmpH, "Hemi_perc"] <- cpgBetaDf[tmpH, "Hemi"]/totCounts
    cpgBetaDf[tmpH, "Hyper_perc"] <- cpgBetaDf[tmpH, "Hyper"]/totCounts
  }
  return(cpgBetaDf)
}

makeDiffBetaDf <- function(inpBeta, inpPheno, sigBool = NULL, meanFunc = NULL, pltBool = NULL, fileExt = NULL){
  # Get difference between mean/median and trimean summaries of CpG site beta values for different phenotypic groups
  # I.e. for any given CpG site, how much does the trimean summary value of the site for group A, differ compared to the mean/median summary?
  resDf <- data.frame(matrix(nrow=0, ncol=7))
  
  if(is.null(meanFunc)){
    meanFunc <- "Mean"
  }else{
    meanFunc <- "Median" 
  }
  
  colnames(resDf) <- c("Ref","Cont",
                       "sigMean","sigTM", 
                       "sigMean_sigTM","sigMean_nSigTM","nSigMean_sigTM")
  for(i in 1:length(table(inpPheno$Histotype))){
    tmpRef <- names(table(inpPheno$Histotype))[i]
    nHists <- names(table(inpPheno$Histotype))[!names(table(inpPheno$Histotype)) %in% tmpRef]
    
    if(!is.null(sigBool)){
      tmpTitle <- paste("Difference between ", meanFunc, " mean and trimean beta for different\ndistribution type comparisons in significant CpG sites" , sep="")
      if(!is.null(fileExt)){
        outFile <- paste(plotPath, "/", fileExt, "_", tmpRef, "_", meanFunc, "_distVarBoxPlot_SigCpg.pdf", sep= "")
      }else{
        outFile <- paste(plotPath, "/", tmpRef, "_", meanFunc, "_distVarBoxPlot_SigCpg.pdf", sep= "")
      }
    }else{
      tmpTitle <- paste("Difference between ", meanFunc," and trimean beta for different\ndistribution type comparisons in CpG sites" , sep="")
      if(!is.null(fileExt)){
        outFile <- paste(plotPath, "/", fileExt, "_", tmpRef, "_", meanFunc, "_distVarBoxPlot_Cpg.pdf", sep= "")
      }else{
        outFile <- paste(plotPath, "/", tmpRef, "_", meanFunc, "_distVarBoxPlot_Cpg.pdf", sep= "")
      }
    }
    histLst <- list()
    for(j in 1:length(nHists)){
      contHist <- nHists[j]
      refBeta <- inpBeta[,inpPheno$barcode[which(inpPheno$Histotype %in% tmpRef)]]
      contBeta <- inpBeta[,inpPheno$barcode[which(inpPheno$Histotype %in% contHist)]]
      # Get mean/median difference between groups
      if(meanFunc %in% "Mean"){
        meanDiffs <- abs(apply(refBeta, 1, function(x) mean(x)) - apply(contBeta, 1, function(x) mean(x)))
      }else if(meanFunc %in% "Median"){
        meanDiffs <- abs(apply(refBeta, 1, function(x) median(x)) - apply(contBeta, 1, function(x) median(x)))
      }else{
        message("Error in meanFunc, this should not be possible, troubleshoot makeDiffBetaDf function")
        return()
      }
      # Get trimean difference between groups
      tmDiffs <- abs(apply(refBeta, 1, function(x) litteR::trimean(x)) - apply(contBeta, 1, function(x) litteR::trimean(x)))
      # Test for normality
      distTypeRef <- apply(refBeta, 1, function(x) shapiro.test(as.numeric(x))$p.value)
      distTypeCont <- apply(contBeta, 1, function(x) shapiro.test(as.numeric(x))$p.value)
      # Test for skewness
      skewRef <- apply(refBeta, 1, function(x) e1071::skewness(as.numeric(x)))
      skewCont <- apply(contBeta, 1, function(x) e1071::skewness(as.numeric(x)))
      # Save results into dataframe
      diffDf <- data.frame(matrix(nrow=length(meanDiffs), ncol=7))
      colnames(diffDf) <- c("CpG", meanFunc, "TM","NormRef", "NormCont", "SkewRef", "SkewCont")
      diffDf[,1] <- names(meanDiffs)
      diffDf[,2] <- meanDiffs
      diffDf[,3] <- tmDiffs
      diffDf[,4] <- distTypeRef
      diffDf[,5] <- distTypeCont
      diffDf[,6] <- skewRef
      diffDf[,7] <- skewCont
      # Keep only rows that are significant for either mean or for trimean
      if(!is.null(sigBool)){
        sigDf <- diffDf[diffDf[,2] >= 0.2 | diffDf[,3] >= 0.2, ]
      }else{
        sigDf <- diffDf
      }
      sigDf$NormRef <- ifelse(sigDf$NormRef > 0.05, "Normal", "NotNormal")
      sigDf$NormCont <- ifelse(sigDf$NormCont > 0.05, "Normal", "NotNormal")
      sigDf$SkewRef <- ifelse(abs(sigDf$SkewRef) >= 1, "Skewed", "NotSkewed")
      sigDf$SkewCont <- ifelse(abs(sigDf$SkewCont) >= 1, "Skewed", "NotSkewed")
      sigDf$normType <- paste(sigDf$NormRef, sigDf$NormCont, sep="_")
      sigDf$normType <- ifelse(sigDf$normType %in% "NotNormal_Normal", "Normal_NotNormal", sigDf$normType)
      sigDf$Mean_TM_Diff <- abs(sigDf[,which(colnames(sigDf) %in% meanFunc)] - sigDf$TM)
      sigDf$Contrast <- paste(tmpRef, contHist, sep="_")
      write.csv(sigDf, paste(outPath, "/TriMean_vs_", meanFunc, "_sig_sites_", paste(tmpRef, contHist, sep="_"), ".csv", sep=""))
      histLst[[length(histLst)+1]] <- sigDf
      names(histLst)[length(histLst)] <- contHist
    }
    # Merge dataframes in list into long format for plotting
    allSigDf = Reduce(function(...) merge(..., all=T), histLst)
    if(!is.null(pltBool)){
      allSigDfLong <- allSigDf %>% 
        pivot_longer(
          cols = c(meanFunc, "TM"), 
          names_to = "Type",
          values_to = "Difference"
        )
      allSigDfLong$Type <- ifelse(allSigDfLong$Type %in% "TM", "TriMean", meanFunc)
      tmpPlt <- ggplot(allSigDfLong, aes(fill=Contrast, y=Mean_TM_Diff, x=Contrast)) + 
        geom_boxplot() + 
        xlab("Contrast") +
        ylab("Delta-beta") + 
        ggtitle(tmpTitle) + 
        scale_fill_viridis(discrete = TRUE) +
        theme(text = element_text(size=24), 
              axis.text.x = element_text(size=16, face="bold"),
              legend.text=element_text(size=16),
              plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~normType, scales = "free_x")
      ggsave(outFile, plot=tmpPlt, width=40, height=30, units = "cm")
    }
    
    sigCats <- names(table(allSigDf$Contrast))
    for(k in 1:length(sigCats)){
      tmpCat <- sigCats[k]
      tmpCats <- strsplit(tmpCat, "_")
      refCat <- tmpCats[[1]][1]
      contCat <- tmpCats[[1]][2]
      
      diffDf <- allSigDf[which(allSigDf$Contrast %in% tmpCat), ]
      
      meanSig <- diffDf[which(diffDf[,2] >= 0.2), ]
      tmSig <- diffDf[which(diffDf[,3] >= 0.2), ]
      
      #normNormMean <- litteR::trimean(as.numeric(diffDf[which(diffDf$normType %in% "Normal_Normal"),"Mean_TM_Diff"]))
      #normNNormMean <- litteR::trimean(as.numeric(diffDf[which(diffDf$normType %in% "Normal_NotNormal"),"Mean_TM_Diff"]))
      #nNormNNormMean <- litteR::trimean(as.numeric(diffDf[which(diffDf$normType %in% "NotNormal_NotNormal"),"Mean_TM_Diff"]))
      
      # Get CpG's which are significant for both trimean, and mean
      mean_TM_allSig <- diffDf[which(diffDf[,2] >= 0.2 & diffDf[,3] >= 0.2), ]
      # Get CpG's which are deemed significant by mean, but not trimean
      meanSig_noTm <- diffDf[which(diffDf[,3] < 0.2 & diffDf[,2] >= 0.2), ]
      # Get CpG's which are deemed significant by trimean, but not mean
      noMean_TmSig <- diffDf[which(diffDf[,2] < 0.2 & diffDf[,3] >= 0.2), ]
      # Add statistics to dataframe
      resDf[nrow(resDf)+1,] <- c(refCat, contCat, 
                                 nrow(meanSig), nrow(tmSig),
                                 nrow(mean_TM_allSig),
                                 nrow(meanSig_noTm),
                                 nrow(noMean_TmSig))
    }
  }
  return(resDf)
}

makeTMVarDiff <- function(inpBeta, inpPheno, fileExt = NULL, pltBool = NULL){
  # Calculate variance between Mean, Median and Trimean of CpG site beta values for different phenotypic groups
  meanValsDf <- data.frame(matrix(nrow=nrow(inpBeta), ncol=(3*length(table(inpPheno$Histotype)))))
  colnames(meanValsDf) <- c(paste(names(table(inpPheno$Histotype)), "Mean", sep="_"),
                            paste(names(table(inpPheno$Histotype)), "TriMean", sep="_"),
                            paste(names(table(inpPheno$Histotype)), "Median", sep="_"))
  rownames(meanValsDf) <- rownames(inpBeta)
  for(i in 1:nrow(inpBeta)){
    if(abs(i)%%round(nrow(inpBeta)/10) == 0){
      message(paste("Processing row: [", i, "/", nrow(inpBeta),"]"))
    }
    for(j in 1:length(names(table(inpPheno$Histotype)))){
      tmpH <- names(table(inpPheno$Histotype))[j]
      hBeta <- inpBeta[i,inpPheno$barcode[which(inpPheno$Histotype %in% tmpH)]]
      hBMat <- as.matrix(hBeta)
      cMean <- mean(hBMat)
      cMed <- median(hBMat)
      cTM <- litteR::trimean(hBMat)
      meanValsDf[i, 
                 c(paste(tmpH, "Mean", sep="_"),
                   paste(tmpH, "TriMean", sep="_"),
                   paste(tmpH, "Median", sep="_"))] <- c(cMean, cTM, cMed)
    }
  }
  triMeanVar <- abs(meanValsDf[,grep("_Mean", colnames(meanValsDf))] - meanValsDf[,grep("_TriMean", colnames(meanValsDf))])
  triMedVar <- abs(meanValsDf[,grep("_Median", colnames(meanValsDf))] - meanValsDf[,grep("_TriMean", colnames(meanValsDf))])
  
  if(!is.null(pltBool)){
    pltDfLongMean <-  triMeanVar %>% 
      pivot_longer(
        cols = colnames(triMeanVar), 
        names_to = "Histotype",
        values_to = "TM_Mean")
    pltDfLongMean <- data.frame(pltDfLongMean) 
    
    pltDfLongMed <-  triMedVar %>% 
      pivot_longer(
        cols = colnames(triMedVar), 
        names_to = "Histotype",
        values_to = "TM_Med")
    pltDfLongMed <- data.frame(pltDfLongMed)  
    
    pltDfLongMean$TM_Med <- pltDfLongMed$TM_Med
    pltDfLongDf <-  pltDfLongMean %>% 
      pivot_longer(
        cols = c("TM_Mean", "TM_Med"), 
        names_to = "VarType",
        values_to = "Var")
    pltDfLongDf <- data.frame(pltDfLongDf)  
    pltDfLongDf$Histotype <- gsub("_Mean", "", pltDfLongDf$Histotype)
    tmpPlt <- ggplot(pltDfLongDf, aes(fill=Histotype, y=Var, x=Histotype)) + 
      geom_boxplot() +
      #geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
      xlab("Histotype") +
      ylab("Variance (Beta)") + 
      ggtitle(paste("Variance (Beta) between Trimean and Mean, Median" , sep="")) + 
      scale_fill_viridis(discrete = TRUE) +
      theme(text = element_text(size=24), 
            axis.text.x = element_text(size=16, face="bold"),
            legend.text=element_text(size=16),
            plot.title = element_text(hjust = 0.5)) +
      facet_grid(rows = vars(VarType), 
                 scales="free_y", switch = 'y')
    outDir <- paste(plotPath, "/", sep="")
    ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
    if(is.null(fileExt)){
      outFile <- paste(outDir, "TM_Mean_Var.pdf", sep="/")
    }else{
      outFile <- paste(outDir, "/TM_Mean_Var_", fileExt, ".pdf", sep="")
    }
    ggsave(outFile, plot=tmpPlt, width=60, height=40, units = "cm")
  }
  
  colnames(triMeanVar) <- gsub("_Mean", "", colnames(triMeanVar))
  colnames(triMedVar) <- gsub("_Median", "", colnames(triMedVar))
  TriMeanVarLong <- 
    TriMedVarLong <- 
    meanMean <- colMeans(triMeanVar)
  meanMed <- colMeans(triMedVar)
  #meanRange <- c(colMins(triMeanVar), colMaxs(triMeanVar))
  #medRange <- c(colMins(as.matrix(triMedVar)), colMaxs(as.matrix(triMedVar)))
  
  outDf <- data.frame(matrix(nrow=6, ncol=4))
  rownames(outDf) <- c("Mean_Mean_vs_TM",
                       "Mean_Med_vs_TM",
                       "Min_Mean_vs_TM",
                       "Max_Mean_vs_TM",
                       "Min_Med_vs_TM",
                       "Max_Med_vs_TM")
  colnames(outDf) <- names(table(inpPheno$Histotype))
  outDf[1, ] <- c(meanMean) 
  outDf[2, ] <- c(meanMed) 
  outDf[3, ] <- colMins(as.matrix(triMeanVar)) 
  outDf[4, ] <- colMaxs(as.matrix(triMeanVar))
  outDf[5, ] <- colMins(as.matrix(triMedVar))
  outDf[6, ] <- colMaxs(as.matrix(triMedVar))
  
  if(is.null(fileExt)){
    write.csv(outDf, paste(outPath, "TM_Mean_Med_Var_Sum.csv"))
  }else{
    write.csv(outDf, paste(outPath, fileExt, "_TM_Mean_Med_Var_Sum.csv"))
  }
  return(outDf)
}

# makeRegSigCpg <- function(inpBeta, inpM, inpPheno, inpGenes, inpPromoCpgs){
#   # Get all CpGs which pass the CW and Dunn requirements
#   allPromoSigCpgs <- list()
#   # Get CpG's identified as significant by the HSP-method
#   for(i in 1:length(inpGenes)){
#     tmpHsp <- inpGenes[[i]]
#     if(nrow(tmpHsp) == 0){
#       next()
#     }
#     tmpH <- names(inpGenes)[i]
#     tmpHistotypes <- names(table(inpPheno$Histotype))[!names(table(inpPheno$Histotype)) %in% tmpH]
#     tmpName <- tmpH
#     
#     # Get genes and associated sample columns
#     tmpBetas <- inpPromoCpgs[which(names(inpPromoCpgs) %in% tmpHsp$ensembl_gene_id)]
#     hBarcodes <- inpPheno$barcode[which(inpPheno$Histotype %in% tmpH)]
#     
#     tmpPheno <- inpPheno
#     allM <- inpM
#     
#     histoFoc <- unique(tmpPheno$Histotype[!tmpPheno$Histotype %in% tmpName])
#     histoFoc <- histoFoc[order(histoFoc)]
#     
#     sigCpgs <- c()
#     for(j in 1:length(tmpBetas)){
#       tmpEns <- names(tmpBetas)[j]
#       geneBDf <- tmpBetas[[j]]
#       if(nrow(geneBDf) == 0 | is.null(nrow(geneBDf)) | nrow(geneBDf) < 3){
#         next()
#       }
#       tmpGeneRow <- tmpHsp[tmpHsp$ensembl_gene_id %in% tmpEns,]
#       tmpSym <- tmpGeneRow$external_gene_name
#       if(is.na(tmpSym) || tmpSym == "" || is.null(tmpSym)){
#         tmpSym <- tmpEns
#       }
#       # Get sample id's for pheno and non-pheno samp
#       phenoIdx <- tmpPheno[which(tmpPheno$Histotype %in% tmpName),]
#       nPhenoIdx <- tmpPheno[which(!tmpPheno$Histotype %in% tmpName),]
#       cpgDf <- inpBeta[which(rownames(inpBeta) %in% rownames(geneBDf)),]
#       # If any rows have NA, remove these as we do not know hwat theirmethylation would have been
#       cpgDf <- na.omit(cpgDf)
#       if(nrow(cpgDf) <= 1){
#         next()
#       }
#       tmpM <- allM[which(rownames(allM) %in% rownames(geneBDf)),]
#       phenoBetas <- cpgDf[,colnames(cpgDf) %in% phenoIdx$barcode]
#       colnames(phenoBetas) <- phenoIdx$Sample_ID[match(colnames(phenoBetas), phenoIdx$barcode)]
#       nPhenoBetas <- cpgDf[,colnames(cpgDf) %in% nPhenoIdx$barcode]
#       colnames(nPhenoBetas) <- nPhenoIdx$Sample_ID[match(colnames(nPhenoBetas), nPhenoIdx$barcode)]
#       phenoBetas <- phenoBetas[,order(colnames(phenoBetas))]
#       nPhenoBetas <- nPhenoBetas[,order(colnames(nPhenoBetas))]
#       # Get two dataframes consisting of the distance between histotypes, and one with the internal variance of probes for histotypes
#       distDf <- makeBetaGrpDists(phenoBetas, nPhenoBetas, tmpPheno,  tmpHistotypes, tmpName, histBool=TRUE, distType="TRIMEAN")
#       # Perform wilcoxon signed rank test for both vector as a whole, and for individual CpG sites
#       # For whole promoter
#       nPhenoTab <- table(tmpPheno[which(!tmpPheno$Histotype %in% tmpName),"Histotype"])
#       # Perform wilcoxon rank sum test for individual promoters
#       wcProbeDf <- data.frame(matrix(nrow=nrow(cpgDf), ncol=length(nPhenoTab)))
#       rownames(wcProbeDf) <- rownames(cpgDf)
#       colnames(wcProbeDf) <- names(nPhenoTab)
#       cwVec <- c()
#       # Cruscal Wallis requires an even "spread" so we use the m-values instead of the beta-values
#       # In statistics, Dunnett's test is a multiple comparison procedure to compare each of a number of treatments with a single control
#       # We use M-values here to avoid the issue of heteroscedasticity found in beta-values
#       hMatch <- inpPheno$Histotype[match(colnames(tmpM), inpPheno$barcode)]
#       for(m in 1:nrow(tmpM)){
#         tmpC <- tmpM[m,]
#         cwDf <- rbind(tmpC, hMatch)
#         cwDf <- data.frame(t(cwDf))
#         colnames(cwDf) <- c("M", "Histotype")
#         # Set levels (first level is reference)
#         cwDf$Histotype <- factor(cwDf$Histotype, levels = c(tmpName, histoFoc))
#         cwTest <- kruskal.test(x=as.numeric(cwDf$M) , g=cwDf$Histotype)
#         cwVec <- append(cwVec, cwTest$p.value)
#         # If there is a significance, we perform Dunn's test to check differences to the reference, between groups
#         dunnCW <- DunnTest(x=as.numeric(cwDf$M) , g=cwDf$Histotype, method = "BH") 
#         tmpCW <- data.frame(dunnCW[[1]])
#         dunnRows <- sapply(rownames(tmpCW), function(x) strsplit(x, "-")[[1]][[2]], USE.NAMES=FALSE)
#         keepDunn <- tmpCW[which(dunnRows %in% tmpName),]
#         tmpP <- keepDunn$pval
#         wcProbeDf[m, ] <-  tmpP
#       }
#       passWc <- wcProbeDf[which(rowSums(wcProbeDf < 0.05) >=3),]
#       distDf <- distDf[rownames(passWc), ]
#       passTM <- distDf[which(rowSums(distDf >0.2) >=3),]
#       passCpG <- rownames(passWc)[which(rownames(passWc) %in% rownames(passTM))]
#       sigCpgs <- append(sigCpgs, passCpG)
#     }
#     allPromoSigCpgs[[i]] <- sigCpgs
#     names(allPromoSigCpgs)[i] <- tmpH
#   }
#   return(allPromoSigCpgs)
# }

makeRegSigCpg <- function(inpReg,
                          inpBeta,
                          inpPheno,
                          inpH,
                          noCats = NULL, 
                          minCpg = NULL, 
                          pCut = NULL, 
                          bCut = NULL,
                          sigThresh = NULL){
  # Retrieves the significant promoters within a pre-defined region dataframe
  if(is.null(noCats)){
    noCats <- length(names(table(inpPheno$Histotype)))-1
  }
  if(is.null(minCpg)){
    minCpg <- 3
  }
  if(is.null(pCut)){
    pCut <- 0.05
  }
  if(is.null(bCut)){
    bCut <- 0.2
  }
  passCpg <- c()
  inpPheno$Sample_ID <- inpPheno$barcode
  # Get all CpGs which pass the CW and Dunn requirements
  tmpHistotypes <- names(table(inpPheno$Histotype))
  histoFoc <- unique(inpPheno$Histotype[!inpPheno$Histotype %in% inpH])
  histoFoc <- histoFoc[order(histoFoc)]
  hBarcodes <- inpPheno$barcode[which(inpPheno$Histotype %in% inpH)]
  geneBDf <- inpBeta[rownames(inpReg), ]
  geneBDf <- na.omit(geneBDf)
  if(nrow(geneBDf) == 0 | is.null(nrow(geneBDf)) | nrow(geneBDf) < minCpg){
    next()
  }
  # Get sample id's for pheno and non-pheno samp
  phenoIdx <- inpPheno[which(inpPheno$Histotype %in% inpH),]
  nPhenoIdx <- inpPheno[which(!inpPheno$Histotype %in% inpH),]
  nPhenoTab <- table(inpPheno[which(!inpPheno$Histotype %in% inpH),"Histotype"])
  # Get histotype specific beta-values
  phenoBetas <- geneBDf[,colnames(geneBDf) %in% phenoIdx$barcode]
  colnames(phenoBetas) <- phenoIdx$Sample_ID[match(colnames(phenoBetas), phenoIdx$barcode)]
  nPhenoBetas <- geneBDf[,colnames(geneBDf) %in% nPhenoIdx$barcode]
  colnames(nPhenoBetas) <- nPhenoIdx$Sample_ID[match(colnames(nPhenoBetas), nPhenoIdx$barcode)]
  phenoBetas <- phenoBetas[,order(colnames(phenoBetas))]
  nPhenoBetas <- nPhenoBetas[,order(colnames(nPhenoBetas))]
  # Get two dataframes consisting of the distance between histotypes, and one with the internal variance of probes for histotypes
  distDf <- makeBetaGrpDists(phenoBetas = phenoBetas, 
                             nPhenoBetas = nPhenoBetas, 
                             tmpPheno = inpPheno,
                             tmpHistotypes = tmpHistotypes,
                             tmpName = inpH,  
                             histBool=TRUE, 
                             distType="TRIMEAN")
  # Perform wilcoxon signed rank test for both vector as a whole, and for individual CpG sites
  # For whole promoter
  # Perform wilcoxon rank sum test for individual promoters
  # Kruscal Wallis requires an even "spread" so we use the m-values instead of the beta-values
  # In statistics, Dunnett's test is a multiple comparison procedure to compare each of a number of treatments with a single control
  # We use M-values here to avoid the issue of heteroscedasticity found in beta-values
  tmpM <- log2(geneBDf/(1-geneBDf))
  kwDf <- makeKW(inpM = tmpM, 
                 inpPheno = inpPheno, 
                 inpH= inpH)
  dunnDf <- makeDunn(inpM = tmpM, 
                     inpPheno = inpPheno,
                     inpH = inpH)
  # Convert result dataframe to wide format compatible with tm distance format
  dunnDf$Contrast <- gsub(paste("-",inpH, sep=""), "", dunnDf$Contrast)
  dunnDfW <- pivot_wider(dunnDf, 
                         id_cols = CpG, 
                         names_from = Contrast, 
                         values_from = pAdj)
  dunnDfW <- data.frame(dunnDfW)
  dunnDfW <- column_to_rownames(dunnDfW, "CpG")
  passDunn <- dunnDfW[which(rowSums(dunnDfW < pCut) >= noCats),]
  passTM <- distDf[which(rowSums(distDf > bCut) >= noCats),]
  if(is.null(sigThresh)){
    passCpg <- intersect(rownames(passDunn), rownames(passTM))
  }else if(toupper(sigThresh) %in% "DIST"){
    passCpg <- rownames(passTM)
  }else if(toupper(sigThresh) %in% "PVAL"){
    passCpg <- rownames(passDunn)
  }
  return(passCpg)
}

# makeRegSigCpg <- function(inpReg,
#                           inpBeta,
#                           inpPheno,
#                           inpH,
#                           noCats = NULL, 
#                           minCpg = NULL, 
#                           pCut = NULL, 
#                           distCut = NULL,
#                           sigThresh = NULL){
#   # Retrieves the significant promoters within a pre-defined region dataframe
#   if(is.null(noCats)){
#     noCats <- length(names(table(inpPheno$Histotype)))-1
#   }
#   if(is.null(minCpg)){
#     minCpg <- 3
#   }
#   if(is.null(pCut)){
#     pCut <- 0.05
#   }
#   if(is.null(distCut)){
#     distCut <- 0.2
#   }
#   
#   passCpg <- c()
#   # Get all CpGs which pass the CW and Dunn requirements
#   tmpHistotypes <- names(table(inpPheno$Histotype))
#   histoFoc <- unique(inpPheno$Histotype[!inpPheno$Histotype %in% inpH])
#   histoFoc <- histoFoc[order(histoFoc)]
#   hBarcodes <- inpPheno$barcode[which(inpPheno$Histotype %in% inpH)]
#   geneBDf <- inpBeta[rownames(inpReg), ]
#   geneBDf <- na.omit(geneBDf)
#   if(nrow(geneBDf) == 0 | is.null(nrow(geneBDf)) | nrow(geneBDf) < minCpg){
#     next()
#   }
#   tmpM <- log2(geneBDf/(1-geneBDf))
#   # Get sample id's for pheno and non-pheno samp
#   phenoIdx <- inpPheno[which(inpPheno$Histotype %in% inpH),]
#   nPhenoIdx <- inpPheno[which(!inpPheno$Histotype %in% inpH),]
#   nPhenoTab <- table(inpPheno[which(!inpPheno$Histotype %in% inpH),"Histotype"])
#   # Get histotype specific beta-values
#   phenoBetas <- geneBDf[,colnames(geneBDf) %in% phenoIdx$barcode]
#   colnames(phenoBetas) <- phenoIdx$Sample_ID[match(colnames(phenoBetas), phenoIdx$barcode)]
#   nPhenoBetas <- geneBDf[,colnames(geneBDf) %in% nPhenoIdx$barcode]
#   colnames(nPhenoBetas) <- nPhenoIdx$Sample_ID[match(colnames(nPhenoBetas), nPhenoIdx$barcode)]
#   phenoBetas <- phenoBetas[,order(colnames(phenoBetas))]
#   nPhenoBetas <- nPhenoBetas[,order(colnames(nPhenoBetas))]
#   # Get two dataframes consisting of the distance between histotypes, and one with the internal variance of probes for histotypes
#   distDf <- makeBetaGrpDists(phenoBetas = phenoBetas, 
#                              nPhenoBetas = nPhenoBetas, 
#                              tmpPheno = inpPheno,
#                              tmpHistotypes = tmpHistotypes,
#                              tmpName = inpH,  
#                              histBool=TRUE, 
#                              distType="TRIMEAN")
#   # Perform wilcoxon signed rank test for both vector as a whole, and for individual CpG sites
#   # For whole promoter
#   # Perform wilcoxon rank sum test for individual promoters
#   dunnProbeDf <- data.frame(matrix(nrow=nrow(geneBDf), 
#                                        ncol=length(nPhenoTab)))
#   rownames(dunnProbeDf) <- rownames(geneBDf)
#   colnames(dunnProbeDf) <- names(nPhenoTab)
#   cwVec <- c()
#   # Cruscal Wallis requires an even "spread" so we use the m-values instead of the beta-values
#   # In statistics, Dunnett's test is a multiple comparison procedure to compare each of a number of treatments with a single control
#   # We use M-values here to avoid the issue of heteroscedasticity found in beta-values
#   hMatch <- inpPheno$Histotype[match(colnames(tmpM), 
#                                      inpPheno$barcode)]
#   for(m in 1:nrow(tmpM)){
#     tmpC <- tmpM[m,]
#     cwDf <- rbind(tmpC, hMatch)
#     cwDf <- data.frame(t(cwDf))
#     colnames(cwDf) <- c("M", "Histotype")
#     # Set levels (first level is reference)
#     cwDf$Histotype <- factor(cwDf$Histotype, 
#                              levels=c(inpH, histoFoc))
#     cwTest <- kruskal.test(x=as.numeric(cwDf$M), 
#                            g=cwDf$Histotype)
#     cwVec <- append(cwVec, cwTest$p.value)
#     # If there is a significance, we perform Dunn's test to check differences to the reference, between groups
#     #dunnCW <- DunnTest(x=as.numeric(cwDf$M), 
#     #                       g=cwDf$Histotype, 
#     #                       method = "BH") 
#     dunnCW <- DunnTest(x=as.numeric(cwDf$M), 
#                        g=cwDf$Histotype, 
#                        method = "none")
#     tmpCW <- data.frame(dunnCW[[1]])
#     dunnRows <- sapply(rownames(tmpCW), function(x) strsplit(x, "-")[[1]][[2]], USE.NAMES=FALSE)
#     keepDunn <- tmpCW[which(dunnRows %in% inpH),]
#     tmpP <- keepDunn$pval
#     dunnProbeDf[m, ] <-  tmpP
#   }
#   passDunn <- dunnProbeDf[which(rowSums(dunnProbeDf < pCut) >= noCats),]
#   passTM <- distDf[which(rowSums(distDf > distCut) >= noCats),]
#   if(is.null(sigThresh)){
#     passCpg <- intersect(rownames(passDunn), rownames(passTM))
#   }else if(toupper(sigThresh) %in% "DIST"){
#     passCpg <- rownames(passTM)
#   }else if(toupper(sigThresh) %in% "PVAL"){
#     passCpg <- rownames(passDunn)
#   }
#   return(passCpg)
# }

getGeneCpgs <- function(geneRow, cpgInp, posInp=NULL, 
                        grBool=TRUE, 
                        verbBool=FALSE, 
                        rangesInp = NULL,
                        strandBool = NULL){
  # Using a row representing a range (gene, dmr etc.) get CpG's overlapping seen to position
  if(!sum(str_detect(geneRow$chr, 'chr')) > 0){
    # Get chromosome index, add chr for comparative purposes
    # Change chromosome name to character for comparative purposes
    geneRow$chr <- as.character(paste("chr", geneRow$chr, sep=""))
  }
  if(geneRow$strand == 1){
    geneRow$strand <- gsub("1", "+", geneRow$strand)
  }else if(geneRow$strand == -1){
    geneRow$strand <- gsub("-1", "-", geneRow$strand)
  }
  # Promoters gets the 2000bp upstream, 200bp downstream from the gene-start
  # While this works in theory, in practice it does not take into account directionality
  # If clause for manual start/end positions
  if(!is.null(posInp)){
    geneStart <- posInp[1]
    geneEnd <- posInp[2]
  }else if(length(which(colnames(geneRow) %like% "promo")) == 0){
    geneStart <- geneRow$start
    geneEnd <- geneRow$end
  }else{
    # If negative strand, TSS is located at "end" of gene
    if(geneRow$strand == "-"){
      geneStart <- geneRow$start
      geneEnd <- geneRow$promo_end
    }else{
      geneStart <- geneRow$promo_start
      geneEnd <- geneRow$end
    }
  }
  geneRow$start <- geneStart
  geneRow$end <- geneEnd
  if(is.null(strandBool)){
    args1 <- c("start", "end", "chr",  "")
    args2 <- args1
  }else if(geneRow$strand == "*"){
    args1 <- c("start", "end", "chr",  NA)
  }else{
    args1 <- c("start", "end", "chr",  "strand")
    args2 <- c("start", "end", "chr",  "strand")
  }
  if(is.null(rangesInp)){
    # Use granges to match geneRow coordinates with CpG reference dataframe
    cpgOls <- makeGrangesOverlaps(df1=geneRow, df2=cpgInp, 
                                  args1=args1, args2=args2)
  }else{
    # Use granges to match geneRow coordinates with input ranges dataframe
    cpgOls <- makeGrangesOverlaps(df1=geneRow, df2=rangesInp, 
                                  args1=args1, args2 = NULL)
  }  
  cpgMatches <- cpgInp[cpgOls@to,]
  if(!nrow(cpgMatches)==0){
    cpgMatches <- as.data.frame(cpgMatches)
    cpgMatches <- cpgMatches[order(cpgMatches$pos, decreasing = FALSE),]
    return(cpgMatches)
  }else{
    if(verbBool){
      message(paste("Error: No CPG matches found for gene: ", geneRow$ensembl_gene_id))
    }
    return(NULL)
  }
}

makeCpgMethTypePerc <- function(inpGenes, inpCpgs, inpPheno, inpBeta, inpFoc = NULL){
  # Make percentage of methylation type of CpGs associated with genes, with respect to phenotype of interest
  if(is.null(inpFoc)){
    inpFoc <- names(order(table(inpPheno$Histotype), decreasing = TRUE))[1]
  }
  resLst <- list()
  tmpGenes <- inpGenes
  tmpH <- inpFoc
  tmpCpgs <- inpCpgs
  geneDf <- data.frame(matrix(nrow=0, ncol=1+(3*(length(names(table(inpPheno$Histotype)))))))
  colnames(geneDf) <- c("Gene", 
                          paste(names(table(inpPheno$Histotype)), "PercHyper", sep="_"),
                          paste(names(table(inpPheno$Histotype)), "PercHemi", sep="_"),
                          paste(names(table(inpPheno$Histotype)), "PercHypo", sep="_"))
  for(j in 1:nrow(tmpGenes)){
    tmpEns <- tmpGenes$ensembl_gene_id[j]
    geneDf[nrow(geneDf)+1, ] <- c(tmpEns, rep(0, 3*(length(names(table(inpPheno$Histotype))))))
    ensCpgs <- tmpCpgs[which(tmpCpgs$Gene %in% tmpEns), ]
    ensCpgs <- ensCpgs[-1]
    ensCpgs <- ensCpgs[!is.na(ensCpgs)]
      
    for(k in 1:length(names(table(inpPheno$Histotype)))){
        tmpHist <- names(table(inpPheno$Histotype))[k]
        tmpBeta <- inpBeta[ensCpgs, inpPheno$barcode[which(inpPheno$Histotype %in% tmpHist)]]
        
        tmpBeta <- na.omit(tmpBeta)
        tmCpgs <- apply(tmpBeta, 1, function(x) litteR::trimean(as.numeric(x)))
        
        percHypo <- length(which(tmCpgs <= 0.3))/length(tmCpgs)
        geneDf[which(geneDf$Gene %in% tmpEns), paste(tmpHist, "PercHypo", sep="_")] <- as.numeric(percHypo)
        
        percHemi <- length(which(tmCpgs < 0.7 & tmCpgs > 0.3))/length(tmCpgs)
        geneDf[which(geneDf$Gene %in% tmpEns), paste(tmpHist, "PercHemi", sep="_")] <- as.numeric(percHemi)
        
        percHyper <- length(which(tmCpgs >= 0.7))/length(tmCpgs) 
        geneDf[which(geneDf$Gene %in% tmpEns), paste(tmpHist, "PercHyper", sep="_")] <- as.numeric(percHyper)
      }
  }
  geneDf <- geneDf[,c(1, (1+order(colnames(geneDf)[-1])))]
  geneDf[nrow(geneDf)+1,] <- c("Sum",  colSums(mutate_all(geneDf[,-1], function(x) as.numeric(as.character(x)))))
  return(geneDf)
}

makePromoBetaMeanTypeStats <- function(inpPromoPos, inpBeta, inpPheno, nCpg = NULL){
  if(is.null(nCpg)){
    nCpg <- 1
  }
  # Summarize HSPs seen to metylation type 
  outDf <- data.frame(matrix(nrow=length(inpPromoPos), ncol=1+length(names(table(inpPheno$Histotype)))))
  colnames(outDf) <- c("ensembl_gene_id", names(table(inpPheno$Histotype)))
  
  countDf <- data.frame(matrix(nrow=length(names(table(inpPheno$Histotype))), ncol=4))
  colnames(countDf) <- c("n=Hypo", "n=Hemi", "n=Hyper", "n=Total")
  rownames(countDf) <- names(table(inpPheno$Histotype))
  countDf <- countDf %>% replace(is.na(.), 0)
  for(i in 1:length(inpPromoPos)){
    if(abs(i)%%round(length(inpPromoPos)/10) == 0){
      message(paste("Processing row: [", i, "/", length(inpPromoPos),"]"))
    }
    tmpG <- inpPromoPos[[i]]
    tmpEns <- names(inpPromoPos)[i]
    outDf[i,"ensembl_gene_id"] <- tmpEns
    if(length(which(rownames(tmpG) %in% rownames(inpBeta))) < nCpg){
      next()
    }
    tmpB <- inpBeta[rownames(tmpG),]
    tmpB <- na.omit(tmpB)
    for(j in 1:length(names(table(inpPheno$Histotype)))){
      bCopy <- tmpB
      tmpH <- names(table(inpPheno$Histotype))[j]
      bCopy  <- bCopy[, which(colnames(bCopy) %in% inpPheno$barcode[which(inpPheno$Histotype %in% tmpH)])]
      hBVec <- c()
      for(k in 1:nrow(bCopy)){
        cB <- litteR::trimean(as.matrix(bCopy[k,]))
        hBVec <- append(hBVec, cB)
      }
      # Add to counters
      countDf[tmpH,"n=Hypo"] <- countDf[tmpH,"n=Hypo"]  + length(which(hBVec < 0.3))
      countDf[tmpH,"n=Hemi"] <- countDf[tmpH,"n=Hemi"]  + length(which(hBVec > 0.3 & hBVec < 0.7))
      countDf[tmpH,"n=Hyper"] <- countDf[tmpH,"n=Hyper"]  + length(which(hBVec > 0.7))
      countDf[tmpH,"n=Total"] <- countDf[tmpH,"n=Total"]  + length(hBVec)
      # Type definition
      typeStat <- ifelse(hBVec < 0.3, "Hypo", ifelse(hBVec > 0.7, "Hyper", "Hemi"))
      mainStat <- names(table(typeStat))[which.max(table(typeStat))]
      outDf[i,tmpH] <- mainStat
    }
  }
  outDf <- na.omit(outDf)
  resDf <- data.frame(matrix(nrow=length(names(table(inpPheno$Histotype))), ncol=7))
  colnames(resDf) <- c("n=", c("Hemi", "Hyper", "Hypo"), paste("Percent", c("Hemi", "Hyper", "Hypo"), sep="_"))
  rownames(resDf) <- names(table(inpPheno$Histotype))
  for(l in 1:length(names(table(inpPheno$Histotype)))){
    tmpHH <- names(table(inpPheno$Histotype))[l]
    tabStats <- table(outDf[,tmpHH])
    if(length(which(!c("Hemi", "Hyper", "Hypo") %in% names(tabStats)))){
      missVals <- c("Hemi", "Hyper", "Hypo")[which(!c("Hemi", "Hyper", "Hypo") %in% names(tabStats))]
      zeroVec <- rep(0, length(missVals))
      names(zeroVec) <- missVals
      tabStats <- append(tabStats, zeroVec)
    }
    tabPerc <- round(tabStats/nrow(outDf),2)
    resDf[tmpHH, ] <- c(nrow(outDf), tabStats, tabPerc)
  }
  resDf <- tibble::rownames_to_column(resDf, "Histotype")
  countDf <- tibble::rownames_to_column(countDf, "Histotype")
  resDf <- dplyr::left_join(resDf, countDf, by = "Histotype")
  merge(resDf, countDf, by = "Histotype", all = TRUE)
  return(resDf)
}

makePromoBetaMeanTypeStats_MULT <- function(inpPromoPos,
                                        inpBeta, 
                                        inpPheno, 
                                        inpCpgs = NULL,
                                        nCores = NULL,
                                        nCpg = NULL) {
  if(is.null(nCpg)){
    nCpg <- 1
  }
  # Summarize HSPs seen to metylation type 
  outDf <- data.frame(matrix(nrow=length(inpPromoPos), ncol=1+length(names(table(inpPheno$Histotype)))))
  colnames(outDf) <- c("ensembl_gene_id", names(table(inpPheno$Histotype)))
  outDf$ensembl_gene_id <- names(inpPromoPos)
  
  countDf <- data.frame(matrix(nrow=length(names(table(inpPheno$Histotype))), ncol=4))
  colnames(countDf) <- c("n=Hypo", "n=Hemi", "n=Hyper", "n=Total")
  rownames(countDf) <- names(table(inpPheno$Histotype))
  countDf <- countDf %>% replace(is.na(.), 0)
  
  # InpCpgs reduces the object size of the dataframe for the future package
  # If not, the objects stored may be too large for the multicore to handle
  if(!is.null(inpCpgs)){
    inpBeta <- inpBeta[allPromoCpgVec, ]
    inpBeta <- na.omit(inpBeta)
  }
  # Make count of promoter stats
  for(i in 1:length(names(table(inpPheno$Histotype)))){
    tmpH <- names(table(inpPheno$Histotype))[i]
    nIt <- length(inpPromoPos)
    hBet <- inpBeta[inpPheno$barcode[which(inpPheno$Histotype %in% tmpH)]]
    # Call future to use multicore
    cpgDistLst <- list()
    cpgDistLst <- future_lapply(1:nIt, function(j){
      if(is.na(abs(j)%%round(nIt/10)) || abs(j)%%round(nIt/10) == 0){
        message(paste("Processing CpG: [", j, "/", nIt,"]"))
      }
      tmpG <- inpPromoPos[[j]]
      tmpEns <- names(inpPromoPos)[j]
      if(length(which(rownames(tmpG) %in% rownames(hBet))) < nCpg){
        cpgDistLst[[j]] <- NULL
      }else{
        gBet <- gBet[rownames(tmpG), ]
        gBet <- na.omit(gBet)
        hBVec <- c()
        for(k in 1:nrow(gBet)){
          cB <- litteR::trimean(as.matrix(gBet[k,]))
          hBVec <- append(hBVec, cB)
        }
        # Add to counters
        hypoC <- length(which(hBVec < 0.3))
        hemiC <- length(which(hBVec > 0.3 & hBVec < 0.7))
        hyperC <- length(which(hBVec > 0.7))
        totC <- length(hBVec)
        # Type definition
        typeStat <- ifelse(hBVec < 0.3, "Hypo", ifelse(hBVec > 0.7, "Hyper", "Hemi"))
        mainStat <- names(table(typeStat))[which.max(table(typeStat))]
        cpgDistLst[[j]] <- c(hypoC, hemiC, hyperC, totC, mainStat) 
      }
    })
    distDf <- do.call(rbind.data.frame, cpgDistLst)
    colnames(distDf) <- c("Hypo", "Hemi", "Hyper", "Total","Main")
    rownames(distDf) <- names(inpPromoPos)
    outDf[, tmpH] <- distDf$Main
    # Add stats to DF for given histotype
    countDf[tmpH,"n=Hypo"] <- sum(as.numeric(distDf[,"Hypo"]))
    countDf[tmpH,"n=Hemi"] <- sum(as.numeric(distDf[,"Hemi"]))
    countDf[tmpH,"n=Hyper"] <- sum(as.numeric(distDf[,"Hyper"]))
    countDf[tmpH,"n=Total"] <- sum(as.numeric(distDf[,"Total"]))
  }
  # Create summary DF
  resDf <- data.frame(matrix(nrow=length(names(table(inpPheno$Histotype))), ncol=7))
  colnames(resDf) <- c("n=", c("Hemi", "Hyper", "Hypo"), paste("Percent", c("Hemi", "Hyper", "Hypo"), sep="_"))
  rownames(resDf) <- names(table(inpPheno$Histotype))
  for(l in 1:length(names(table(inpPheno$Histotype)))){
    tmpHH <- names(table(inpPheno$Histotype))[l]
    tabStats <- table(outDf[,tmpHH])
    if(length(which(!c("Hemi", "Hyper", "Hypo") %in% names(tabStats)))){
      missVals <- c("Hemi", "Hyper", "Hypo")[which(!c("Hemi", "Hyper", "Hypo") %in% names(tabStats))]
      zeroVec <- rep(0, length(missVals))
      names(zeroVec) <- missVals
      tabStats <- append(tabStats, zeroVec)
    }
    tabPerc <- round(tabStats/nrow(outDf),2)
    resDf[tmpHH, ] <- c(nrow(outDf), tabStats, tabPerc)
  }
  resDf <- tibble::rownames_to_column(resDf, "Histotype")
  countDf <- tibble::rownames_to_column(countDf, "Histotype")
  resDf <- dplyr::left_join(resDf, countDf, by = "Histotype")
  resDf <- merge(resDf, countDf, by = "Histotype", all = TRUE)
  return(resDf)
}

makeBetaTM_Mean_Var <- function(inpCpg, inpPheno, inpBeta, inpGeneInf, sigDiffVal = NULL){
  # Function for calculating the difference between Trimean and Mean
  if(is.null(sigDiffVal)){
    sigDiffVal <- 0.1
  }
  diffLst <- list()
  for(i in 1:length(inpCpg)){
    tmpH <- names(inpCpg)[i]
    if(!tmpH %in% names(table(inpPheno$Histotype))){
      next()
    }
    tmpHistotypes <- names(table(inpPheno$Histotype))[!names(table(inpPheno$Histotype)) %in% tmpH]
    tmpCDf <- inpCpg[[tmpH]]
    if(nrow(tmpCDf) == 0){
      next()
    }
    hDf <- data.frame(matrix(nrow=nrow(tmpCDf), ncol=3+(4*length(table(inpPheno$Histotype)[!names(table(inpPheno$Histotype)) %in% tmpH]))))
    colnames(hDf) <- c("ensembl_gene_id", "external_gene_name", "Histotype",
                       paste(names(table(inpPheno$Histotype))[!names(table(inpPheno$Histotype)) %in% tmpH], "mean_NormNorm", sep="_"),
                       paste(names(table(inpPheno$Histotype))[!names(table(inpPheno$Histotype)) %in% tmpH], "mean_NormNotNorm", sep="_"),
                       paste(names(table(inpPheno$Histotype))[!names(table(inpPheno$Histotype)) %in% tmpH], "range_NormNorm", sep="_"),
                       paste(names(table(inpPheno$Histotype))[!names(table(inpPheno$Histotype)) %in% tmpH], "range_NormNotNorm", sep="_"))
    cpgCols <- grep("cpg", colnames(tmpCDf), fixed = TRUE)
    # Prepare reference dataframe
    refDf <- makeCpgDf(tmpCDf, inpBeta, inpPheno, tmpH)
    refDf <- data.frame(t(refDf))
    refDf$Histotype <- inpPheno$Histotype[match(rownames(refDf), inpPheno$barcode)]
    phenoBetas <- refDf[which(refDf$Histotype %in% tmpH),]
    phenoBetas <- phenoBetas[,-which(colnames(phenoBetas) %in% "Histotype")]
    phenoBetas <- t(phenoBetas)
    # Prepare dataframe of other histotypes
    colnames(phenoBetas) <- inpPheno$Sample_ID[match(colnames(phenoBetas), inpPheno$barcode)]
    nPhenoBetas <- refDf[-which(refDf$Histotype %in% tmpH),]
    nPhenoBetas <- nPhenoBetas[,-which(colnames(nPhenoBetas) %in% "Histotype")]
    nPhenoBetas <- t(nPhenoBetas)
    colnames(nPhenoBetas) <- inpPheno$Sample_ID[match(colnames(nPhenoBetas), inpPheno$barcode)]
    tmDists <- makeBetaGrpDists(phenoBetas, nPhenoBetas, inpPheno,  tmpHistotypes, tmpH, histBool=TRUE, distType="TRIMEAN")
    meanDists <- makeBetaGrpDists(phenoBetas, nPhenoBetas, inpPheno,  tmpHistotypes, tmpH, histBool=TRUE, distType="MEAN")
    for(j in 1:nrow(tmpCDf)){
      tmpGr <- tmpCDf[j,]
      tmpEns <- tmpGr$Gene
      
      tmpSym <- inpGeneInf$external_gene_name[which(inpGeneInf$ensembl_gene_id %in% tmpEns)]
      if(tmpSym == "" | is.na(tmpSym)){
        tmpSym <- tmpEns
      }
      
      hDf[j,c("ensembl_gene_id", "external_gene_name", "Histotype")] <- c(tmpEns, tmpSym, tmpH) 
      
      tmpGr <- tmpGr[,cpgCols]
      tmpGr <- tmpGr[!is.na(tmpGr)]
      tmpGr <- tmpGr[which(tmpGr %in% rownames(inpBeta))]
      if(length(tmpGr) <= 1){
        next()
      }
      
      distTypeDf <- data.frame(matrix(nrow=length(tmpGr), ncol=length(table(inpPheno$Histotype))))
      skewTypeDf <- data.frame(matrix(nrow=length(tmpGr), ncol=length(table(inpPheno$Histotype))))
      rownames(distTypeDf) <- tmpGr
      rownames(skewTypeDf) <- tmpGr
      colnames(distTypeDf) <- names(table(inpPheno$Histotype))
      colnames(skewTypeDf)<- names(table(inpPheno$Histotype))
      for(l in 1:length(table(inpPheno$Histotype))){
        colsH <- names(table(inpPheno$Histotype))[l]
        tmpCols <- inpPheno$barcode[which(inpPheno$Histotype %in% colsH)]
        colsDf <- refDf[tmpCols,tmpGr]
        for(m in 1:ncol(colsDf)){
          tmpC <- colsDf[,m]
          tmpCName <- colnames(colsDf)[m]
          normRes <- shapiro.test(as.numeric(tmpC))
          if(normRes$p.value > 0.05){
            distTypeDf[tmpCName,colsH] <- "Normal"
          }else{
            distTypeDf[tmpCName,colsH] <- "NotNormal"
          }
          skewRes <- e1071::skewness(as.numeric(tmpC))
          if(abs(skewRes) >= 1){
            skewTypeDf[tmpCName,colsH]  <- "Skew"
          }else{
            skewTypeDf[tmpCName,colsH] <- "NotSkew"
          }
        }
      }
      
      rTM <- tmDists[tmpGr,]
      rM <- meanDists[tmpGr,]
      diffDf <- abs(rTM -  rM)
      
      distTypeComp <- distTypeDf
      
      for(s in 1:nrow(distTypeDf)){
        for(t in 1:ncol(distTypeDf)){
          distTypeComp[s,t] <- paste(distTypeDf[s,1], distTypeDf[s,t],sep="_")
        }
      }
      distTypeComp <- distTypeComp[,-1]
      distTypeCompLong <- distTypeComp %>% 
        pivot_longer(
          cols = c(1:3), 
          names_to = "Histotype",
          values_to = "Comp"
        )
      distTypeCompLong <- data.frame(distTypeCompLong)
      
      if(!length(table(distTypeCompLong$Comp)) >1){
        next()
      }
      distTypeCompLong$Comp <-  gsub("NotNormal_Normal", "Normal_NotNormal", distTypeCompLong$Comp)
      
      mCopy <- rM 
      tmCopy <- rTM
      
      mLong <-mCopy %>% 
        pivot_longer(
          cols = c(1:3), 
          names_to = "Histotype",
          values_to = "Beta"
        )
      mLong <- data.frame(mLong)
      
      tmLong <-tmCopy %>% 
        pivot_longer(
          cols = c(1:3), 
          names_to = "Histotype",
          values_to = "Beta"
        )
      tmLong <- data.frame(tmLong)
      
      mLong$Comp <- distTypeCompLong$Comp
      tmLong$Comp <- distTypeCompLong$Comp
      mLong$Beta <- as.numeric(mLong$Beta)
      tmLong$Beta <- as.numeric(tmLong$Beta)
      
      tmLong$Mean <- mLong$Beta
      tmLong$Diff <- abs(mLong$Beta - tmLong$Beta)
      
      diffRange <- c(min(tmLong$Diff), max(tmLong$Diff))
      for(t in 1:length(table(tmLong$Histotype))){
        compH <- names(table(tmLong$Histotype))[t]
        longH <- tmLong[which(tmLong$Histotype %in% compH), ]  
        
        if(!"Normal_Normal" %in% longH$Comp){
          normRange <- c(0,0)
          notNormRange <- c(min(longH$Diff[which(longH$Comp %in% "Normal_NotNormal")]),
                            max(longH$Diff[which(longH$Comp %in% "Normal_NotNormal")]))
          compNotNorm <- mean(longH$Diff)
          compNorm <- 0
        }else if(!"Normal_NotNormal" %in% longH$Comp){
          normRange <- c(min(longH$Diff[which(longH$Comp %in% "Normal_Normal")]),
                         max(longH$Diff[which(longH$Comp %in% "Normal_Normal")]))
          notNormRange <- c(0,0)
          compNotNorm <- 0
          compNorm <- mean(longH$Diff)
        }else{
          normRange <- c(min(longH$Diff[which(longH$Comp %in% "Normal_Normal")]),
                         max(longH$Diff[which(longH$Comp %in% "Normal_Normal")]))
          notNormRange <- c(min(longH$Diff[which(longH$Comp %in% "Normal_NotNormal")]),
                            max(longH$Diff[which(longH$Comp %in% "Normal_NotNormal")]))
          diffCompVals <- as.data.frame(longH[,c("Comp", "Diff")] %>% 
                                          group_by(Comp) %>%  
                                          summarise_all(.funs = mean))
          compNorm <- diffCompVals$Diff[which(diffCompVals$Comp %in% "Normal_Normal")]
          compNotNorm <- diffCompVals$Diff[which(diffCompVals$Comp %in% "Normal_NotNormal")]
        }
        hDf[j,c(paste(compH, "mean_NormNorm", sep="_"),
                paste(compH, "mean_NormNotNorm", sep="_"),
                paste(compH, "range_NormNorm", sep="_"),
                paste(compH, "range_NormNotNorm", sep="_"))] <- c(compNorm,
                                                                  compNotNorm,
                                                                  paste(normRange, collapse=":"),
                                                                  paste(notNormRange, collapse=":"))
      }
    }
    diffLst[[i]] <- hDf
    names(diffLst)[i] <- tmpH
  }
  return(diffLst)
}

makeMeanTMDiffBetaDf <- function(inpMat, inpPheno, inpH, limRows = NULL, pltBool = NULL, fileExt = NULL){
  if(!is.null(limRows)){
    inpMat <- inpMat[1:limRows, ]
  }
  inpMat <- inpMat[,which(colnames(inpMat) %in% inpPheno$barcode)]
  inpPheno <- inpPheno[match(colnames(inpMat),inpPheno$barcode), ]
  nHists <- names(table(inpPheno$Histotype))[!names(table(inpPheno$Histotype)) %in% inpH]
  refCols <- colnames(inpMat)[which(colnames(inpMat) %in% inpPheno$barcode[which(inpPheno$Histotype %in% inpH)])]
  contLst <- list()
  for(j in 1:length(nHists)){
    nRef <- nHists[j]
    message(paste("Running pairwise central tendency calculation between: ", inpH, " vs. ", nRef, sep=""))
    # First, we get distribution types for each CpG in the dataframe (with respect to contrast and reference)
    # Then, we perform the students t-test for each CpG
    contCols <- colnames(inpMat)[which(colnames(inpMat) %in% inpPheno$barcode[which(inpPheno$Histotype %in% nRef)])]
    # Finally, we must loop through all rows to get Dunn's p-value
    pDf <- data.frame(matrix(nrow=0,ncol=3))
    colnames(pDf) <- c("CpG", "Mean", "TM")
    nIt <- nrow(inpMat)
    for(k in 1:nIt){
      if(is.na(abs(k)%%round(nIt/10)) || abs(k)%%round(nIt/10) == 0){
        message(paste("Processing CpG: [", k, "/", nIt,"]"))
      }
      tmpCpg <- rownames(inpMat)[k]
      mRow <- inpMat[tmpCpg,]
      refVals <- as.numeric(mRow[refCols])
      nRefVals <- as.numeric(mRow[contCols])
      # Calculate mean difference
      meanRef <- mean(refVals)
      meanNRef <- mean(nRefVals)
      meanDiff <- meanRef - meanNRef
      # Calculate tm difference
      tmRef <- makeTm(refVals)
      tmNRef <- makeTm(nRefVals)
      tmDiff <- tmRef - tmNRef
      pDf[nrow(pDf)+1,] <- c(tmpCpg, meanDiff,  tmDiff)
    }
    # After we have done the t-test, dunn comparison, continue to determine distribution types
    distTypeRef <- apply(inpMat[, refCols], 1, function(x) shapiro.test(as.numeric(x))$p.value)
    distTypeCont <- apply(inpMat[, contCols], 1, function(x) shapiro.test(as.numeric(x))$p.value)
    # Create dataframe of values
    diffDf <- data.frame(matrix(nrow=nrow(pDf), ncol=5))
    colnames(diffDf) <- c("CpG", "Mean", "TM", "NormRef", "NormCont")
    diffDf[,1] <- pDf$CpG
    diffDf[,2] <- pDf$Mean
    diffDf[,3] <- pDf$TM
    diffDf[,4] <- distTypeRef
    diffDf[,5] <- distTypeCont
    # Translate comparison type
    diffDf$NormRef <- ifelse(diffDf$NormRef > 0.05, "Normal", 
                             "NotNormal")
    diffDf$NormCont <- ifelse(diffDf$NormCont > 0.05, "Normal",
                              "NotNormal")
    diffDf$normType <- paste(diffDf$NormRef, diffDf$NormCont, sep="_")
    diffDf$normType <- ifelse(diffDf$normType %in% "NotNormal_Normal", "Normal_NotNormal", diffDf$normType)
    # Get difference in p-values
    diffDf$diffVal <- abs(as.numeric(diffDf$Mean) - as.numeric(diffDf$TM))
    cont <- paste(nRef, inpH, sep="_")
    diffDf$Contrast <- cont
    contLst[[length(contLst)+1]] <- diffDf
    names(contLst)[length(contLst)] <- cont
  }
  # Merge all comparisons into one dataframe
  allSigDf = Reduce(function(...) merge(..., all=T), contLst)
  allSigDfLong <- allSigDf %>% 
    pivot_longer(
      cols = c("Mean", "TM"), 
      names_to = "Type",
      values_to = "Val"
    )
  allSigDfLong$diffVal <- as.numeric(allSigDfLong$diffVal)
  allSigDfLong$Type <- ifelse(allSigDfLong$Type %in% "Mean", "Mean", "TriMean")
  # Remove duplicate rows (we only care about the difference between the tests here)
  allSigDfLong <- allSigDfLong[!allSigDfLong$Type %in% "TriMean",]
  # Get the average difference for each comparisons type
  meanDiff <- data.frame(allSigDfLong %>% 
                            group_by(normType, Contrast) %>% 
                            summarise(meanValType = mean(diffVal, na.rm = TRUE)))
  meanDiff$meanValType <- round(meanDiff$meanValType, 4) 
  for(l in 1:nrow(meanDiff)){
    pCat <- meanDiff$normType[l]
    contType <- meanDiff$Contrast[l]
    nRow <- length(which(allSigDfLong$normType %in% pCat & allSigDfLong$Contrast %in% contType))
    allSigDfLong$Contrast <- ifelse(allSigDfLong$normType %in% pCat & allSigDfLong$Contrast %in% contType, 
                                    paste(allSigDfLong$Contrast, "\nmean diff=", meanDiff$meanValType[l], "\nno cpg=", nRow, sep=""), allSigDfLong$Contrast)
  }
  
  if(!is.null(pltBool)){
    if(!is.null(fileExt)){
      outFile <- paste(plotPath,  inpH, "_", fileExt, "_Tm_Mean_Diff_SigCpg.pdf", sep= "")
    }else{
      outFile <- paste(plotPath,  inpH, "_Tm_Mean_Diff_.pdf", sep= "")
    }
    
    # Plot differences
    tmpTitle <- paste("Difference between Mean and TriMean in different\n distribution type comparisons in CpG sites", sep="")
    # Manually color categories
    colCats <- names(table(allSigDfLong$Contrast))
    colvals <- viridis(n=length(colCats)/length(nHists), option="B")
    colVec <- rep(colvals, length(nHists))
    colVec <- colVec[order(colVec)]
    names(colVec) <- colCats
    tmpPlt <- ggplot(allSigDfLong, aes(fill=Contrast, y=diffVal, x=Contrast)) + 
      geom_boxplot() + 
      xlab("Contrast") +
      ylab("Difference in adjusted p-value") + 
      ggtitle(tmpTitle) + 
      scale_fill_manual(values = colVec) + 
      theme(text = element_text(size=24), 
            axis.text.x = element_text(size=16, face="bold"),
            legend.text=element_text(size=16),
            plot.title = element_text(hjust = 0.5)) +
      theme(legend.position="none") +
      facet_wrap(~normType, scales = "free_x")
    ggsave(outFile, plot=tmpPlt, width=50, height=30, units = "cm")
  }
  # Get difference in p-value in pairwise comparisons in a contrast against all other possible contrasts
  resDf <- data.frame(matrix(nrow=0, ncol=7))
  colnames(resDf) <- c("Ref","Cont",
                       "sigTT","sigBrunnM", 
                       "sigTT_sigBrunnM","sigTT_nSigBrunnM","nSigTT_sigBrunnM")
  sigCats <- names(table(allSigDf$Contrast))
  for(m in 1:length(sigCats)){
    tmpCat <- sigCats[m]
    tmpCats <- strsplit(tmpCat, "_")
    refCat <- tmpCats[[1]][1]
    contCat <- tmpCats[[1]][2]
    diffDf <- allSigDf[which(allSigDf$Contrast %in% tmpCat), ]
    meanSig <- diffDf[which(diffDf[,2] > 0.2), ]
    tmSig <- diffDf[which(diffDf[,3] > 0.2), ]
    # Get CpG's which are significant for both trimean, and mean
    Mean_TM_allSig <- diffDf[which(diffDf[,2] > 0.2 & diffDf[,3] > 0.2), ]
    # Get CpG's which are deemed significant by mean, but not trimean
    MeanSig_noTM <- diffDf[which(diffDf[,2] > 0.2 & diffDf[,3] < 0.2), ]
    # Get CpG's which are deemed significant by trimean, but not mean
    noMean_TMSig <- diffDf[which(diffDf[,2] < 0.2 & diffDf[,3] > 0.2), ]
    # Add statistics to dataframe
    resDf[nrow(resDf)+1,] <- c(refCat, contCat, 
                               nrow(meanSig), nrow(tmSig),
                               nrow(Mean_TM_allSig),
                               nrow(MeanSig_noTM),
                               nrow(noMean_TMSig))
  }
  return(list("ResDf" = resDf, "allValDf" = allSigDf))
}

makePairPvalDiffBetaDf <- function(inpMat, inpPheno, inpH, limRows = NULL, pltBool = NULL, fileExt = NULL){
  if(!is.null(limRows)){
    inpMat <- inpMat[1:limRows, ]
  }
  inpMat <- inpMat[,which(colnames(inpMat) %in% inpPheno$barcode)]
  inpPheno <- inpPheno[match(colnames(inpMat),inpPheno$barcode), ]
  
  nHists <- names(table(inpPheno$Histotype))[!names(table(inpPheno$Histotype)) %in% inpH]
  refCols <- colnames(inpMat)[which(colnames(inpMat) %in% inpPheno$barcode[which(inpPheno$Histotype %in% inpH)])]
  contLst <- list()
  for(j in 1:length(nHists)){
    nRef <- nHists[j]
    message(paste("Running pairwise comparison between: ", inpH, " vs. ", nRef, sep=""))
    # First, we get distribution types for each CpG in the dataframe (with respect to contrast and reference)
    # Then, we perform the students t-test for each CpG
    contCols <- colnames(inpMat)[which(colnames(inpMat) %in% inpPheno$barcode[which(inpPheno$Histotype %in% nRef)])]
    # Finally, we must loop through all rows to get Dunn's p-value
    pDf <- data.frame(matrix(nrow=0,ncol=3))
    colnames(pDf) <- c("CpG", "TTest", "BM")
    nIt <- nrow(inpMat)
    for(k in 1:nIt){
      if(is.na(abs(k)%%round(nIt/10)) || abs(k)%%round(nIt/10) == 0){
          message(paste("Processing CpG: [", k, "/", nIt,"]"))
      }
      tmpCpg <- rownames(inpMat)[k]
      mRow <- inpMat[tmpCpg,]
      ttRow <- mRow[c(refCols, contCols)]
      mLab <- inpPheno$Histotype[match(names(ttRow) , 
                                       inpPheno$barcode)]
      mLab <- factor(mLab, 
                     levels=c(inpH, nRef))
      # Perform pairwise t-test
      pTTest <- pairwise.t.test(as.numeric(ttRow), 
                                mLab, 
                                p.adjust.method = "none")
      ttP <- pTTest$p.value
      # Perform Brunner-Munzels test
      refVals <- as.numeric(mRow[refCols])
      nRefVals <- as.numeric(mRow[contCols])
      brunnM <- brunnermunzel::brunnermunzel.test(x = refVals,
                                                  y = nRefVals)
      brunnMP <- brunnM$p.value
      pDf[nrow(pDf)+1,] <- c(tmpCpg, ttP,  brunnMP)
    }
    # Adjust p-values based on total number of pairwise comparisons made (i.e. number of CpGs in Df) 
    pDf$ttAdjP <- p.adjust(as.numeric(pDf$TTest), 
                           method = "BH", 
                           n = nrow(pDf))
    pDf$bmAdjP <- p.adjust(as.numeric(pDf$BM), 
                           method = "BH", 
                           n = nrow(pDf))
    # After we have done the t-test, dunn comparison, continue to determine distribution types
    distTypeRef <- apply(inpMat[, refCols], 1, function(x) shapiro.test(as.numeric(x))$p.value)
    distTypeCont <- apply(inpMat[, contCols], 1, function(x) shapiro.test(as.numeric(x))$p.value)
    # Create dataframe of values
    diffDf <- data.frame(matrix(nrow=nrow(pDf), ncol=7))
    colnames(diffDf) <- c("CpG", "TTp", "BMp", "adjustedTT", "adjustedBrunn", "NormRef", "NormCont")
    diffDf[,1] <- pDf$CpG
    diffDf[,2] <- pDf$TTest
    diffDf[,3] <- pDf$BM
    diffDf[,4] <- pDf$ttAdjP
    diffDf[,5] <- pDf$bmAdjP
    diffDf[,6] <- distTypeRef
    diffDf[,7] <- distTypeCont
    # Translate comparison type
    diffDf$NormRef <- ifelse(diffDf$NormRef > 0.05, "Normal", 
                            "NotNormal")
    diffDf$NormCont <- ifelse(diffDf$NormCont > 0.05, "Normal",
                             "NotNormal")
    diffDf$normType <- paste(diffDf$NormRef, diffDf$NormCont, sep="_")
    diffDf$normType <- ifelse(diffDf$normType %in% "NotNormal_Normal", "Normal_NotNormal", diffDf$normType)
    # Get difference in p-values
    diffDf$diffP <- abs(as.numeric(diffDf$TTp) - as.numeric(diffDf$BMp))
    diffDf$diffAdjP <- abs(as.numeric(diffDf$adjustedTT) - as.numeric(diffDf$adjustedBrunn))
    cont <- paste(nRef,inpH, sep="_")
    diffDf$Contrast <- cont
    contLst[[length(contLst)+1]] <- diffDf
    names(contLst)[length(contLst)] <- cont
  }
  # Merge all comparisons into one dataframe
  allSigDf = Reduce(function(...) merge(..., all=T), contLst)
  allSigDfLong <- allSigDf %>% 
      pivot_longer(
        cols = c("BMp", "TTp"), 
        names_to = "Type",
        values_to = "pVal"
      )
  allSigDfLong$Type <- ifelse(allSigDfLong$Type %in% "TTp", "t-test", "Brunner-Munzel")
  # Remove duplicate rows (we only care about the difference between the tests here)
  allSigDfLong <- allSigDfLong[!allSigDfLong$Type %in% "t-test",]
  # Get the average difference for each comparisons type
  meanPDiff <- data.frame(allSigDfLong %>% 
                            group_by(normType, Contrast) %>% 
                            summarise(meanPType = mean(diffAdjP, na.rm = TRUE)))
  meanPDiff$meanPType <- round(meanPDiff$meanPType, 4) 
  for(l in 1:nrow(meanPDiff)){
    pCat <- meanPDiff$normType[l]
    contType <- meanPDiff$Contrast[l]
    nRow <- length(which(allSigDfLong$normType %in% pCat & allSigDfLong$Contrast %in% contType))
    allSigDfLong$Contrast <- ifelse(allSigDfLong$normType %in% pCat & allSigDfLong$Contrast %in% contType, 
                                    paste(allSigDfLong$Contrast, "\nmean diff=", meanPDiff$meanPType[l], "\nno cpg=", nRow, sep=""), allSigDfLong$Contrast)
  }
  
  if(!is.null(pltBool)){
    if(!is.null(fileExt)){
      outFile <- paste(plotPath, "/", inpH, "_", fileExt, "_pvalVarBoxPlot_SigCpg.pdf", sep= "")
    }else{
      outFile <- paste(plotPath, "/", inpH, "_pvalVarBoxPlot.pdf", sep= "")
    }
    
    # Plot differences
    tmpTitle <- paste("Difference between adjusted p-values (BH) for students t-test and Brunner-Munzels\n test in different distribution type comparisons in CpG sites", sep="")
    # Manually color categories
    colCats <- names(table(allSigDfLong$Contrast))
    colvals <- viridis(n=length(colCats)/length(nHists), option="B")
    colVec <- rep(colvals, length(nHists))
    colVec <- colVec[order(colVec)]
    names(colVec) <- colCats
    tmpPlt <- ggplot(allSigDfLong, aes(fill=Contrast, y=diffAdjP, x=Contrast)) + 
      geom_boxplot() + 
      xlab("Contrast") +
      ylab("Difference in adjusted p-value") + 
      ggtitle(tmpTitle) + 
      scale_fill_manual(values = colVec) + 
      theme(text = element_text(size=24), 
            axis.text.x = element_text(size=16, face="bold"),
            legend.text=element_text(size=16),
            plot.title = element_text(hjust = 0.5)) +
      theme(legend.position="none") +
      facet_wrap(~normType, scales = "free_x")
    ggsave(outFile, plot=tmpPlt, width=50, height=30, units = "cm")
  }
  # Get difference in p-value in pairwise comparisons in a contrast against all other possible contrasts
  resDf <- data.frame(matrix(nrow=0, ncol=7))
  colnames(resDf) <- c("Ref","Cont",
                       "sigTT","sigBrunnM", 
                       "sigTT_sigBrunnM","sigTT_nSigBrunnM","nSigTT_sigBrunnM")
  sigCats <- names(table(allSigDf$Contrast))
  for(m in 1:length(sigCats)){
      tmpCat <- sigCats[m]
      tmpCats <- strsplit(tmpCat, "_")
      refCat <- tmpCats[[1]][1]
      contCat <- tmpCats[[1]][2]
      diffDf <- allSigDf[which(allSigDf$Contrast %in% tmpCat), ]
      ttSig <- diffDf[which(diffDf[,2] < 0.05), ]
      bmSig <- diffDf[which(diffDf[,3] < 0.05), ]
      # Get CpG's which are significant for both trimean, and mean
      TT_BM_allSig <- diffDf[which(diffDf[,2] < 0.05 & diffDf[,3] < 0.05), ]
      # Get CpG's which are deemed significant by mean, but not trimean
      TTSig_noBM <- diffDf[which(diffDf[,2] < 0.05 & diffDf[,3] >= 0.05), ]
      # Get CpG's which are deemed significant by trimean, but not mean
      noTT_BMSig <- diffDf[which(diffDf[,2] > 0.05 & diffDf[,3] < 0.05), ]
      # Add statistics to dataframe
      resDf[nrow(resDf)+1,] <- c(refCat, contCat, 
                                 nrow(ttSig), nrow(bmSig),
                                 nrow(TT_BM_allSig),
                                 nrow(TTSig_noBM),
                                 nrow(noTT_BMSig))
  }
  return(list("ResDf" = resDf, "allValDf" = allSigDf))
}

makeMultPvalDiffBetaDf <- function(inpMat, inpPheno, inpH,  pltBool = NULL, limRows = NULL, fileExt = NULL){
  if(!is.null(limRows)){
    inpMat <- inpMat[1:limRows, ]
  }
  
  inpMat <- inpMat[,which(colnames(inpMat) %in% inpPheno$barcode)]
  inpPheno <- inpPheno[match(colnames(inpMat),inpPheno$barcode), ]
  
  nHists <- names(table(inpPheno$Histotype))[!names(table(inpPheno$Histotype)) %in% inpH]
  refCols <- colnames(inpMat)[which(colnames(inpMat) %in% inpPheno$barcode[which(inpPheno$Histotype %in% inpH)])]
  contLst <- list()
  message(paste("Running pairwise comparison between: ", inpH, " and other histotypes: ", paste(nHists, collapse=", "), sep=""))
  allSigDf <- data.frame(matrix(nrow=0,ncol=9))
  colnames(allSigDf) <- c("CpG", "AOV", "KW", "Tukey", "Dunn", "NormRef", "NormCont", "diffAdjP", "Contrast")
  nIt <- nrow(inpMat)
  for(k in 1:nIt){
    if(is.na(abs(k)%%round(nIt/10)) || abs(k)%%round(nIt/10) == 0){
       message(paste("Processing CpG: [", k, "/", nIt,"]"))
    }
    tmpCpg <- rownames(inpMat)[k]
    mRow <- inpMat[tmpCpg,]
    mLab <- inpPheno$Histotype[match(names(mRow) , 
                                       inpPheno$barcode)]
    mLab <- factor(mLab, 
                     levels=c(inpH, nHists))
    
    mDf <- data.frame(as.numeric(mRow), mLab)
    colnames(mDf) <- c("Val", "Histotype")
    
    # Perform anova
    aovRes <- aov(Val ~ Histotype, data = mDf)
    aovP <- summary(aovRes)[[1]][["Pr(>F)"]][1]
    # Get pairwise p.values using tukeys test
    tukeyRes <- TukeyHSD(aovRes)
    tukeyRes <- tukeyRes[[1]]
    # Perform kruskal-wallis
    kwRes <- kruskal.test(x=as.numeric(mDf$Val), 
                 g=mDf$Histotype)
    kwP <- kwRes$p.value
    dunnRes <- DunnTest(x=as.numeric(mDf$Val), 
                        g=mDf$Histotype, 
                        method = "BH")
    dunnRes <- dunnRes[[1]]
    for(j in 1:length(nHists)){
      nRef <- nHists[j]
      cRow <- paste(nRef, inpH, sep="-")
      tmpCont <- paste(nRef, inpH, sep="_")
      tkP <- tukeyRes[cRow, ]
      tkP <- tkP[[4]]
      dP <- dunnRes[cRow, ]
      dP <- dP[[2]]
      # Get distribution type of comparison
      refDist <- shapiro.test(as.numeric(mRow[which(names(mRow) %in% inpPheno$barcode[which(inpPheno$Histotype %in% inpH)])]))
      nRefDist <- shapiro.test(as.numeric(mRow[which(names(mRow) %in% inpPheno$barcode[which(inpPheno$Histotype %in% nRef)])]))
      pDiff <-abs(as.numeric(tkP)-as.numeric(dP))
      allSigDf[nrow(allSigDf)+1, ] <- c(tmpCpg, aovP, kwP, tkP, dP, refDist$p.value, nRefDist$p.value, pDiff, tmpCont)
    }
  }
  # Translate comparison type
  allSigDf$NormRef <- ifelse(allSigDf$NormRef > 0.05, "Normal", 
                             "NotNormal")
  allSigDf$NormCont <- ifelse(allSigDf$NormCont > 0.05, "Normal",
                              "NotNormal")
  allSigDf$normType <- paste(allSigDf$NormRef, allSigDf$NormCont, sep="_")
  allSigDf$normType <- ifelse(allSigDf$normType %in% "NotNormal_Normal", "Normal_NotNormal", allSigDf$normType)
  # Merge all comparisons into one dataframe
  allSigDfLong <- allSigDf %>% 
    pivot_longer(
      cols = c("Tukey", "Dunn"), 
      names_to = "Type",
      values_to = "pVal"
    )
  allSigDfLong$Type <- ifelse(allSigDfLong$Type %in% "Tukey", "Tukeys", "DunnTest")
  # Remove duplicate rows (we only care about the difference between the tests here)
  allSigDfLong <- allSigDfLong[!allSigDfLong$Type %in% "DunnTest",]
  allSigDfLong$diffAdjP <- as.numeric(allSigDfLong$diffAdjP)
  # Get the average difference for each comparisons type
  meanPDiff <- data.frame(allSigDfLong %>% 
                            group_by(normType, Contrast) %>% 
                            summarise(meanPType = mean(diffAdjP, 
                                                       na.rm = TRUE)))
  meanPDiff$meanPType <- round(meanPDiff$meanPType, 4) 
  for(l in 1:nrow(meanPDiff)){
    pCat <- meanPDiff$normType[l]
    contType <- meanPDiff$Contrast[l]
    nRow <- length(which(allSigDfLong$normType %in% pCat & allSigDfLong$Contrast %in% contType))
    allSigDfLong$Contrast <- ifelse(allSigDfLong$normType %in% pCat & allSigDfLong$Contrast %in% contType, 
                                    paste(allSigDfLong$Contrast, "\nmean diff=", meanPDiff$meanPType[l], "\nno cpg=", nRow, sep=""), allSigDfLong$Contrast)
  }
  
  if(!is.null(pltBool)){
    if(!is.null(fileExt)){
      outFile <- paste(plotPath, "/", inpH, "_", fileExt, "_pvalVar_Multiple_Comp_BoxPlot_SigCpg.pdf", sep= "")
    }else{
      outFile <- paste(plotPath, "/", inpH, "_pvalVar_Multiple_Comp_BoxPlot.pdf", sep= "")
    }
    
    # Plot differences
    tmpTitle <- paste("Difference between adjusted p-values (BH) for Tukeys and Dunns test for multiple comparisons\n in different distribution type comparisons in CpG sites", sep="")
    # Manually color categories
    colCats <- names(table(allSigDfLong$Contrast))
    colvals <- viridis(n=length(colCats)/length(nHists), option="B")
    colVec <- rep(colvals, length(nHists))
    colVec <- colVec[order(colVec)]
    names(colVec) <- colCats
    tmpPlt <- ggplot(allSigDfLong, aes(fill=Contrast, y=diffAdjP, x=Contrast)) + 
      geom_boxplot() + 
      xlab("Contrast") +
      ylab("Difference in adjusted p-value") + 
      ggtitle(tmpTitle) + 
      scale_fill_manual(values = colVec) + 
      theme(text = element_text(size=24), 
            axis.text.x = element_text(size=16, face="bold"),
            legend.text=element_text(size=16),
            plot.title = element_text(hjust = 0.5)) +
      theme(legend.position="none") +
      facet_wrap(~normType, scales = "free_x")
    ggsave(outFile, plot=tmpPlt, width=50, height=30, units = "cm")
  }
  # Get difference in p-value in pairwise comparisons in a contrast against all other possible contrasts
  resDf <- data.frame(matrix(nrow=0, ncol=7))
  colnames(resDf) <- c("Ref","Cont",
                       "sigTuk","sigDunn", 
                       "sigTuk_sigDunn","sigTuk_nSigDunn","nSigTuk_sigDunn")
  sigCats <- names(table(allSigDf$Contrast))
  for(m in 1:length(sigCats)){
    tmpCat <- sigCats[m]
    tmpCats <- strsplit(tmpCat, "_")
    refCat <- tmpCats[[1]][1]
    contCat <- tmpCats[[1]][2]
    diffDf <- allSigDf[which(allSigDf$Contrast %in% tmpCat), ]
    tukSig <- diffDf[which(diffDf[,4] < 0.05), ]
    dunnSig <- diffDf[which(diffDf[,5] < 0.05), ]
    # Get CpG's which are significant for both trimean, and mean
    Tuk_Dunn_allSig <- diffDf[which(diffDf[,4] < 0.05 & diffDf[,5] < 0.05), ]
    # Get CpG's which are deemed significant by mean, but not trimean
    TukSig_noDunn <- diffDf[which(diffDf[,4] < 0.05 & diffDf[,5] >= 0.05), ]
    # Get CpG's which are deemed significant by trimean, but not mean
    noTuk_DunnSig <- diffDf[which(diffDf[,4] > 0.05 & diffDf[,5] < 0.05), ]
    # Add statistics to dataframe
    resDf[nrow(resDf)+1,] <- c(refCat, contCat, 
                               nrow(tukSig), nrow(dunnSig),
                               nrow(Tuk_Dunn_allSig),
                               nrow(TukSig_noDunn),
                               nrow(noTuk_DunnSig))
  }
  return(list("ResDf" = resDf, "allValDf" = allSigDf))
}

################################################################################
# External validation scripts
################################################################################

makeExtRankFilt <- function(extRankLst, inpCpgCov, inpPromoLst){
  # Filter out uneccessary columns from promoter ranking output
  extRankFilt <- list()
  for(i in 1:length(extRankLst)){
    tmpRanks <- extRankLst[[i]]
    tmpSid <- names(extRankLst)[i]
    for(j in 1:length(tmpRanks)){
      skipCols <- c("ensembl_gene_id","external_gene_name", "sumDist", "brownSum",
                    "histVar", "nHistVarSum", "varSumDist")
      tmpExtRank <- tmpRanks[[j]]
      varCols <- colnames(tmpExtRank)[which(grepl("var",colnames(tmpExtRank)))]
      # varCols <- colnames(tmpExtRank)[which(colnames(tmpExtRank) %like% "var")]
      skipCols <- append(skipCols, varCols)
      incCols <- which(!colnames(tmpExtRank) %in% skipCols)
      # Filter based on differential DNA-methylation
    }
    extRankFilt[[i]] <- tmpRanks
    names(extRankFilt)[i] <- tmpSid
  }
  
  freqExtSumLst <- list()
  for(i in 1:length(histotypes)){
    filtRankLst <- extRankLst[c("GSE51820", "GSE226823")]
    #filtRankLst <- extRankLst[c("GSE51820")]
    tmpHLst <- list()
    tmpH <- histotypes[i]
    if(!tmpH %in% names(inpPromoLst)){
      next()
    }
    hColNames <- c()
    histoFilt <- histotypes[!histotypes %in% tmpH]
    for(k in 1:length(names(filtRankLst))){
      tmpName <- names(filtRankLst)[k]
      hColNames <- append(hColNames, paste(tmpName, histoFilt , sep="_"))
    }
    sumExtDf <- data.frame(matrix(nrow=nrow(inpPromoLst[[tmpH]]) , ncol=5+length(filtRankLst)*3+length(hColNames)))
    colnames(sumExtDf) <- c("ensembl_gene_id", "external_gene_name","freq", "cpgEPIC", "cpg450K", paste("CpG", names(filtRankLst), sep="_") ,names(filtRankLst), paste(names(filtRankLst), "sumDist", sep="_"), hColNames)
    sumExtDf$ensembl_gene_id <- inpPromoLst[[tmpH]]$ensembl_gene_id
    sumExtDf$external_gene_name <- inpPromoLst[[tmpH]]$external_gene_name
    sumExtDf$cpgEPIC <- inpCpGCov[[tmpH]]$cpg_EPIC[match(sumExtDf$ensembl_gene_id, inpCpGCov[[tmpH]]$ensembl_gene_id)]
    sumExtDf$cpg450K <- inpCpGCov[[tmpH]]$cpg_450K[match(sumExtDf$ensembl_gene_id, inpCpGCov[[tmpH]]$ensembl_gene_id)]
    for(j in 1:length(filtRankLst)){
      tmpRank <- filtRankLst[[j]][[tmpH]]
      if(is.null(tmpRank)){
        next()
      }
      filtCols <- which(colnames(tmpRank) %in% histoFilt)
      tmpSid <- names(filtRankLst)[j]
      colInd <- which(colnames(sumExtDf) %in% names(filtRankLst)[j])
      # Set frequency to 0, then update to 1 if tmpRank contains the gene afer filtering
      sumExtDf[, colInd] <- 0
      sumExtDf[match(tmpRank$ensembl_gene_id, sumExtDf$ensembl_gene_id), colInd] <- 1
      tmpColName <- paste(tmpSid, "sumDist", sep="_")
      # Add information to rows of interest
      sumExtDf[match(tmpRank$ensembl_gene_id, sumExtDf$ensembl_gene_id), tmpColName] <- as.numeric(tmpRank$sumDist)
      
      tmpColName_2 <- paste("CpG", tmpSid, sep="_")
      sumExtDf[match(tmpRank$ensembl_gene_id, sumExtDf$ensembl_gene_id), tmpColName_2] <- as.numeric(tmpRank$nCpg)
      hCols <- tmpRank[, which(colnames(tmpRank) %in% histoFilt)]
      rownames(hCols) <- tmpRank$ensembl_gene_id
      colnames(hCols) <- paste(tmpSid, colnames(hCols), sep="_")
      sumExtDf[match(rownames(hCols), sumExtDf$ensembl_gene_id), match(colnames(hCols), colnames(sumExtDf))] <- hCols
      tmpHNames <- paste(tmpSid, histoFilt, sep="_")
    }
    if(tmpH %in% "MC"){
      sumExtDf <- sumExtDf[,-which(grepl("GSE226823",colnames(sumExtDf)))]
      cLengths <- length(filtRankLst)-1
      sumExtDf[is.na(sumExtDf)] <- 0
      sumExtDf$freq <- sumExtDf[,(6+cLengths)] 
      #sumExtDf <- sumExtDf[!sumExtDf$freq == 0,]
      sumDistDf <- sumExtDf[,which(grepl("sumDist",colnames(sumExtDf)))]
      sumExtDf$totSumDist <- sumDistDf
    }else{
      cLengths <- length(filtRankLst)
      sumExtDf[is.na(sumExtDf)] <- 0
      sumExtDf$freq <- rowSums(sumExtDf[,(5+(2*cLengths)-1):(5+(2*cLengths))])/cLengths 
      # sumExtDf <- sumExtDf[!sumExtDf$freq == 0,]
      sumDistDf <- sumExtDf[,which(grepl("sumDist",colnames(sumExtDf)))]
      sumDistDf[is.na(sumDistDf)] <- 0
      sumExtDf$totSumDist <- rowSums(sumDistDf)
    }
    sumExtDf$freq <- as.numeric(sumExtDf$freq)
    sumExtDf$freqDist <- sumExtDf$totSumDist*sumExtDf$freq
    # We order by totSumDist, since it will order similar to freq, but based on separation
    sumExtDf <- sumExtDf[order(sumExtDf$freqDist, decreasing = TRUE),]
    sumExtDf  <- sumExtDf %>% relocate(totSumDist, .after = freq)
    sumExtDf  <- sumExtDf %>% relocate(freqDist, .after = totSumDist)
    freqExtSumLst[[length(freqExtSumLst)+1]] <- sumExtDf
    names(freqExtSumLst)[length(freqExtSumLst)] <- tmpH
  }
  freqExtSumLst <- freqExtSumLst[sapply(freqExtSumLst, nrow) > 0]
  return(freqExtSumLst)
}


