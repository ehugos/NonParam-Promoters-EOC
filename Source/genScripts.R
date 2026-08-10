################################################################################
################################################################################
################################################################################
# Collection of scripts not associated with a specific point of analysis
################################################################################
################################################################################
################################################################################

unregister_dopar <- function() {
  # Fixes the error; Error in summary.connection(connection) : invalid connection
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

detachAllPackages <- function(){
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if(length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

# Function for creating named list of dataframes
makeDfLst <- function(fileLst, sepC=NULL){
  if(is.null(sepC)){
    sepC=","
  }
  outLst <- list()
  for(i in 1:length(fileLst)){
    fPath <- fileLst[i]
    tmpDf <- read.csv(fPath, row.names=1, as.is = TRUE, sep = sepC)
    fSplt <- strsplit(fPath, "/")[[1]]
    tmpName <- fSplt[(length(fSplt)-1)]
    if("X" %in% colnames(tmpDf)){
      rownames(tmpDf) <- tmpDf[,"X"]
      tmpDf <- tmpDf[,-(which(colnames(tmpDf) %in% "X"))]
    }
    outLst[[i]] <- tmpDf
    names(outLst)[i] <- toupper(tmpName)
  }
  return(outLst)
}

getFocusedComps <- function(inpDeg, histo, revBool=NULL){
  # Remove duplicate entries of DEG results (as direction does not matter for Deseq2)
  # We thus prefer to only keep the entries where the histo is used as reference
  contLst <- names(inpDeg)
  contTmp <- contLst
  for(i in 1:(length(contLst))){
    spltcName <- strsplit(contLst[i], "_")
    if(!is.null(revBool)){
      if(spltcName[[1]][1] != histo){
        tmpInd <- which(contTmp == contLst[i])
        contTmp <- contTmp[-tmpInd]
      }
    }else{
      if(spltcName[[1]][2] != histo){
        tmpInd <- which(contTmp == contLst[i])
        contTmp <- contTmp[-tmpInd]
      }
    }
  }
  return(contTmp)
}

makeListNames <- function(inpLst){
  # New function for quickly "fixing" a list name seen to histotype vs meth results 
  outLst <- inpLst
  degCont <- NA
  dmpCont <- NA
  for(i in 1:length(outLst)){
    tmpCont <- names(outLst)[i]
    # Check so that we have not added lists in the wrong order
    if(is.na(stringr::str_extract(tmpCont, regex("[A-Z]{2,4}_VS_[A-Z]{2,4}"))[[1]])){
      dmpCont <- tmpCont
    }else if(is.na(stringr::str_extract(tmpCont, regex("[A-Z]{2,4} - [A-Z]{2,4}"))[[1]])){
      degCont <- tmpCont
    }else{
      message(paste("ERROR: Element", i, "in list 1", tmpCont ,"does not confer to associated results patterns, exiting!", sep=" "))
      break
    } 
    if(is.na(degCont)){
      # Change name of lists to be equal for downstream analysis
      dmpSplt <- strsplit(dmpCont, " - ")
      cont1 <- dmpSplt[[1]][1]
      cont2 <- dmpSplt[[1]][2]
    }else{
      degSplt <- strsplit(degCont, "_")
      cont1 <- degSplt[[1]][2]
      cont2 <- degSplt[[1]][4]
    }
    newName <- paste(cont1, cont2, sep="_")
    names(outLst)[i] <- newName
  }
  return(outLst)
}

makeGRangesCompatible <- function(inpDf){
  type <- detLstType(inpDf)
  newCols <- c("start", "end", "chr", "strand")
  if(all(newCols %in% colnames(inpDf))){
    startCol <- "start"
    endCol <- "end"
    chrCol <- "chr"
    strandCol <- "strand"
    message("Dataframe already made compatible for Granges comparison")
  }else{
    if(type == "DMP"){
      startCol <- "MAPINFO"
      inpDf$MAPEND <- inpDf[,startCol] + 1
      endCol <- "MAPEND"
      chrCol <- "CHR"
      strandCol <- "Strand"
    }else if(type == "DMR"){
      startCol <- "start"
      endCol <- "end"
      chrCol <- "seqnames"
      strandCol <- "strand"
    }else if(type == "DEG"){
      startCol <- "start_position"
      endCol <- "end_position"
      chrCol <- "chromosome_name"
      strandCol <- "strand"
    }else{
      message("Error: Dataframe of unknown input-type")
      return()
    }
  }
  # Check if "chr" is missing from chr column, if so add it
  if(sum(str_detect(inpDf[,chrCol], 'chr')) == 0 | is.na(sum(str_detect(inpDf[,chrCol], 'chr')))){
    # Get chromosome index, add chr for comparative purposes
    inpDf[,chrCol] <- as.character(paste("chr", inpDf[,chrCol], sep=""))
  }
  
  # Change strand parameter if needed to accepted granges format
  if(!strandCol %in% colnames(inpDf)){
    inpDf$strand <- rep("*", nrow(inpDf))
  }else if(sum(str_detect(replace_na(as.character(inpDf[,strandCol]), ""), 'F')) > 0 | sum(str_detect(replace_na(as.character(inpDf[,strandCol]), ""), 'R')) > 0){
    inpDf[strandCol][inpDf[strandCol] == "F"] <- "+"
    inpDf[strandCol][inpDf[strandCol] == "R"] <- "-"
  }else if(sum(str_detect(replace_na(as.character(inpDf[,strandCol]), ""), '1')) > 0 | sum(str_detect(replace_na(as.character(inpDf[,strandCol]), ""), '-1')) > 0){
    inpDf[strandCol][inpDf[strandCol] == "1"] <- "+"
    inpDf[strandCol][inpDf[strandCol] == "-1"] <- "-"
  }
  
  # Check function for FALSE entries within strand parameter,
  # if present we need to ignor strand
  if(sum(str_detect(replace_na(as.character(inpDf[,strandCol]), ""), 'FALSE')) > 0 | sum(is.na(inpDf[,strandCol])) > 0){
    inpDf[strandCol][inpDf[strandCol] == FALSE] <- "*"
    inpDf[strandCol][is.na(inpDf[strandCol])] <- "*"
  }
  names(inpDf)[names(inpDf) == startCol] <-  newCols[1]
  names(inpDf)[names(inpDf) == endCol] <-  newCols[2]
  names(inpDf)[names(inpDf) == chrCol] <-  newCols[3]
  names(inpDf)[names(inpDf) == strandCol] <-  newCols[4]
  # Finally, remove any rows that contain NA's fr start/end
  inpDf <- inpDf[!is.na(inpDf$start),]
  inpDf <- inpDf[!is.na(inpDf$end),]
  return(inpDf)
}

detLstType <- function(inpDf){
  type <- NULL
  # List-type determination, checks for column names unique to list-type to be present
  if(length(which(colnames(inpDf) %in% c("external_gene_name", "ensembl_gene_id"))) == 2){
    type <- "DEG"
  }else if(length(which(colnames(inpDf) %in% c("overlapping.genes", "no.cpgs"))) == 2){
    type <- "DMR"
  }else if(length(which(colnames(inpDf) %in% c("gene", "B"))) == 2){
    type <- "DMP"
  }else if(length(which(colnames(inpDf) %in% c("start", "end", "chr", "strand"))) == 4){
    type <- "PROCESSED"
  }
  if(is.null(type)){
    stop("Error: Unknown list type, not of type: DMP, DMR or DEG")
  }else{
    return(type)
  }
}

makeListFolders <- function(inpLst, subDir = "", plotDir = TRUE, datDir = TRUE){
  # Create directories for output if they do not already exist
  for(i in 1:length(names(inpLst))){
    if(!subDir == ""){
      newOut <- paste(outPath, names(inpLst)[i], subDir, sep="/")
      newPlot <- paste(plotPath, names(inpLst)[i], subDir, sep="/")
    }else{
      newOut <- paste(outPath, names(inpLst)[i], sep="/")
      newPlot <- paste(plotPath, names(inpLst)[i], sep="/") 
    }
    if(plotDir){
      ifelse(!dir.exists(file.path(newPlot)), dir.create(file.path(newPlot)), FALSE)
    }
    if(datDir){
      ifelse(!dir.exists(file.path(newOut)), dir.create(file.path(newOut)), FALSE)
    }
  }
}

getGeneCounts <- function(inpGenes, inpCounts){
  # Script for getting sample specific counts for genes
  outLst <- list()
  for(i in 1:length(inpGenes)){
    tmpName <- names(inpGenes)[i]
    tmpGenes <- inpGenes[[i]]
    # If this is a contrast gene-list, take reference as "focus"
    if(grepl("_", tmpName, fixed = TRUE)){
      tmpCont1 <- strsplit(tmpName, "_")[[1]][1]
      tmpCont2 <- strsplit(tmpName, "_")[[1]][2]
      pCols <- trainPheno$Sample_ID[which(trainPheno$Histotype %in% c(tmpCont1, tmpCont2))]
    }else{
      pCols <- trainPheno$Sample_ID[which(trainPheno$Histotype %in% tmpName)]
    }
    pCols <- gsub("OV_", "", pCols)
    pCounts <- inpCounts[,which(colnames(inpCounts) %in% pCols)] 
    tmpSigCounts <- pCounts[rownames(pCounts) %in% tmpGenes$ensembl_gene_id,]
    outLst[[i]] <- tmpSigCounts
    names(outLst)[i] <- tmpName
  }
  return(outLst)
}

makeCsvSave <- function(resLst, sType, rBool = NULL, dirExtBool = NULL){
  if(is.null(rBool)){
    rBool <- TRUE
  }
  # Function for saving output dataframe into .csv format
  for(i in 1:length(resLst)){
    outRes <- resLst[[i]]
    if(is.null(nrow(outRes))){
      next()
    }else{
      chaRes <- apply(outRes,2,as.character)
      # In the offchance that conversion goes wrong, save as the original format
      if(is.null(dim(chaRes))){
        chaRes <- outRes
      }
      row.names(chaRes) <- row.names(outRes)
      sName <- names(resLst)[i]
      
      if(!is.null(dirExtBool)){
        tmpDir <- paste(outPath, dirExtBool, "/", sep="")
        ifelse(!dir.exists(tmpDir), dir.create(tmpDir), FALSE)
        outDir <- paste(tmpDir, "/", sName, sep="")
      }else{
        outDir <- paste(outPath, sName, sep="")
      }
      ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
      write.csv(chaRes, paste(outDir, "/", sName, "_", sType, ".csv", sep=""),  row.names=rBool)
    }
  }
}

makeGrangesOverlaps <- function(df1, df2, args1=NULL, args2=NULL, df1IsRanges = NULL){
  # Simple function for doing a g-ranges comparison to identify coordinate overlaps
  # Args must be in the format c(start, end, chr, strand)
  # If strand = NULL, ignore strand will be set to TRUE for ranges construction 
  if(is.null(df1IsRanges)){
    if(is.null(args1)){
      message("ERROR in makeGrangesOverlaps: Need args1 to create ranges1")
      return()
    }
  }
  strand1Bool <- FALSE
  if(is.na(args1[[4]]) | args1[[4]] == ""){
    strand1Bool <- TRUE
  }
  # Create granges object
  ranges1 <- makeGRangesFromDataFrame(df1,
                                      ignore.strand=strand1Bool,
                                      seqinfo=NULL,
                                      start.field=args1[[1]],
                                      end.field=args1[[2]],
                                      seqnames.field=args1[[3]],
                                      strand.field=args1[[4]])
  # Create G-ranges object for DMR
  strand2Bool <- FALSE
  if(is.null(args2)){
    # If args2 is not submitted, we will assume that the other DF is a Granges object
    ranges2 <- df2
  }else{
    # First check if strand is a submitted argument
    if(is.na(args2[[4]]) | args2[[4]] == ""){
      strand2Bool <- TRUE
    }
    # If args2 has "cpg", we update the dataframe to allow for CpG matching
    if(toupper(args2[[1]]) == "CPG" | toupper(args2[[2]]) == "CPG"){
      # As CpGs are at a specific location "pos", we make the granges to be between [pos,pos+2]
      # This is since each probe is a CpG site (i.e. CG pairing)
      df2$start <- df2$pos
      df2$end <- df2$pos+2
      df2 <- within(df2, rm("pos"))
      args2 <- c("start","end","chr","strand")
    }
    ranges2 <- makeGRangesFromDataFrame(df2,
                                          ignore.strand=strand2Bool,
                                          seqinfo=NULL,
                                          start.field=args2[[1]],
                                          end.field=args2[[2]],
                                          seqnames.field=args2[[3]],
                                          strand.field=args2[[4]]) 
  }
  dfOverlaps <- GenomicRanges::findOverlaps(ranges1, ranges2)
  return(dfOverlaps)
}

makeDegFreqDf <- function(inpDeg, excLst=NULL){
  # Updated script for finding HSGs, creates a summary dataframe complete with 
  # frequency of the DEG and its l2fc if available for a contrast
  tmpHisto <- histotypes
  if(!is.null(excLst)){
    tmpHisto <- tmpHisto[!tmpHisto %in% excLst]
  }
  outDegLst <- list()
  for(i in 1:length(tmpHisto)){
    tmpH <- tmpHisto[i]
    comps <-  getFocusedComps(inpDeg, tmpH)
    if(!is.null(excLst)){
      for(j in 1:length(excLst)){
        comps <- comps[!comps %like% excLst[j]]
      }
    }
    # Retrieve all entries that contain/match the histotype comparison
    tmpDegLst <- inpDeg[names(inpDeg) %in% comps]
    # Combine matrix into one combined matrix
    tmpDegMat <- do.call(rbind,tmpDegLst)
    unqEns <- unique(tmpDegMat$ensembl_gene_id)
    hDegDf <- data.frame(matrix(nrow=0, ncol=4+length(tmpDegLst)))
    # Create DEG frequency dataframe
    colnames(hDegDf) <- c("ensembl_gene_id", "external_gene_name","noDeg", "degFreq", comps)
    for(j in 1:length(unqEns)){
      tmpEns <- unqEns[[j]]
      ensRows <- tmpDegMat[tmpDegMat$ensembl_gene_id %in% tmpEns, ]
      ensRows <- tibble::rownames_to_column(ensRows, "Contrast")
      ensRows$Contrast <- gsub("\\.[0-9]{1,4}","",ensRows$Contrast)
      degFreq <- round(nrow(ensRows)/length(comps),3)
      rowInd <- nrow(hDegDf)+1
      hDegDf[rowInd, ] <- c(tmpEns, ensRows$external_gene_name[1], nrow(ensRows), degFreq, rep(NA, length(comps)))
      for(k in 1:nrow(ensRows)){
        if(ensRows$Contrast[k] %in% colnames(hDegDf)){
          tmpColInd <- which(colnames(hDegDf) %in% ensRows$Contrast[k])
          hDegDf[rowInd , tmpColInd] <- ensRows$log2FoldChange[k]
        }
      }
    }
    hDegDf <- hDegDf[order(hDegDf$degFreq, decreasing = TRUE), ]
    outDegLst[[i]] <- hDegDf
    names(outDegLst)[i] <- tmpH
  }
  return(outDegLst)
}
