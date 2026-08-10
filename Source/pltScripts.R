################################################################################
################################################################################
################################################################################
# Scripts for plotting not assocaited with predictive classification
################################################################################
################################################################################
################################################################################

makeFullGeneCpgBoxPlot <- function(inpGenes, inpBeta, inpAnno, cpgInp, geneInfInp, focusGrp="Histotype", fileExt=NULL, dirExt=NULL, allBool=FALSE, pltBool=NULL){
  # Function for creating CpG boxplots for gene in provided dataframe
  # Get overlapping CpGs for each gene in input
  if(is.null(pltBool)){
    pltBool <- TRUE
  }
  
  for(i in 1:length(inpGenes)){
    tmpGenes <- inpGenes[[i]]
    contName <- names(inpGenes)[i]
    if(nrow(tmpGenes) == 0){
      message(paste("No genes found for contrast:", contName, "Proceeding to next contrast", sep=" "))
      next
    }
    spltCont <- strsplit(contName, "_")
    if(length(spltCont) == 1){
      if(allBool){
        tmpCont <- paste(histotypes[order(histotypes)], collapse='_')
      }else{
        tmpCont <- contName
      }
    }else{
      tmpCont <- contName
    }
    
    tmpGLst <- list(tmpGenes)
    names(tmpGLst) <- tmpCont
    cpgGeneOverlap <- makeCpgGeneOverlap(tmpGLst, cpgInp, geneInfInp)
    cpgBetas <- getGeneCpgBeta(cpgGeneOverlap, inpBeta, inpAnno)
    if(length((cpgBetas)[[1]]) > 0){
      tmpCpgs <- cpgBetas[[tmpCont]]
    }else{
      message(paste("No matching CpG sites for genes in contrast:", contName, "Proceeding to next contrast", sep=" "))
      next
    }
    colInd <- which(names(inpAnno) %in% focusGrp)
    catColInd <- which(names(catColProf) %in% focusGrp)
    # Extract DDS object based on the reference in contrast name
    # Function requires input to be in list format, convert
    # tmpGenes <- tmpGenes[order(tmpGenes$padj), ]
    noGenes <- nrow(tmpGenes)
    message(paste("Creating boxplots for: [", noGenes, "] Entries for contrast: ", tmpCont, sep=""))
    pltLst <- list()
    for(j in 1:nrow(tmpGenes)){
      message(paste("Plotting gene [", j, "/", noGenes, "]", sep=""))
      tmpLst <- list()
      tmpRow <- tmpGenes[j,]
      ensName <- tmpRow$ensembl_gene_id
      geneSym <- tmpRow$external_gene_name
      cpgPos <- tmpCpgs[[ensName]]
      # Make sure that ALL overlapping CpGs for gene are on the correct strand
      # cpgPos <- cpgPos[which(cpgPos$strand == tmpRow$strand),]
      if(is.null(cpgPos)){
        message(paste("No overlapping CpGs found for gene:", geneSym, "proceeding to next gene",sep=" "))
        next
      }else if(nrow(cpgPos) < 2){
        message(paste("Insufficient number of overlapping CpGs:", nrow(cpgPos), "for gene", geneSym, "proceeding to next gene",sep=" "))
        next
      }else{
        if(is.null(dirExt)){
          outDir <- paste(plotPath,"/",contName,"/cpgBoxPlot/DEG/", geneSym,"/", sep="")
        }else{
          outDir <- paste(plotPath,"/",contName,"/cpgBoxPlot/DEG/", dirExt,"/",geneSym,"/", sep="")
        }
        ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
        if(!is.null(fileExt)){
          outFile <- paste(outDir, geneSym, "_", fileExt, "_","Gene_Cpg_BoxPlot_", contName, ".pdf", sep="")
        }else{
          outFile <- paste(outDir, geneSym, "_","Gene_Cpg_BoxPlot_", contName, ".pdf", sep="")
        }
        # Filter our plot-vector
        skipVec <- c("Row.names", "chr", "pos", "strand")
        pltDf <- as.data.frame(t(cpgPos[,!colnames(cpgPos) %in% skipVec]))
        colnames(pltDf) <- cpgPos$Row.names
        # Create matching column for samples in input-annotation
        inpAnno$matchCol <- paste(inpAnno$Histotype, inpAnno$barcode , sep="_")
        pltDf$Histotype <- inpAnno$Histotype[match(row.names(pltDf), inpAnno$matchCol)]
        mergeDf <- tidyr::gather(pltDf, key="CpG", value="Beta",  colnames(pltDf)[1]:colnames(pltDf)[ncol(pltDf)-1], factor_key=TRUE)
        # Add CpG position as factor
        mergeDf$pos <- cpgPos$pos[match(mergeDf$CpG , cpgPos$Row.names)]
        mergeDf$cpg_class <- ifelse(between(mergeDf$pos, 
                                            tmpRow$promo_start, 
                                            tmpRow$promo_end), "promo", "body")
        # Sort dataframe, then add an ordered position vector of factors to represent the gene coordinates
        mergeDf <- mergeDf[order(factor(mergeDf$pos, levels = levels(factor(mergeDf$pos)), ordered = TRUE)),]
        mergeDf$pos <- as.factor(mergeDf$pos)
        # Convert to factor, sort separately
        # mergeDf$pos_ordered <- factor(mergeDf$pos, levels = levels(factor(mergeDf$pos)), ordered = TRUE)
        # mergeDf$pos <- as.factor(mergeDf$pos)
        
        if(length(which(mergeDf$cpg_class == "promo")) > 0 & length(which(mergeDf$cpg_class == "body")) > 0){
          # If we have a promoter entry, we create two separate plots to scale them
          promoMat <- mergeDf[which(mergeDf$cpg_class == "promo"),]
          bodyMat <- mergeDf[which(!mergeDf$cpg_class == "promo"),]
          # Change the point size, and shape
          promoPlot <- ggplot(promoMat, aes(x=pos, y=Beta, fill=Histotype)) + 
            geom_boxplot() +
            labs(# x=paste(geneRow$chr, geneRow$start, geneRow$end, geneRow$strand, sep="::"), 
              y="Beta",
              title="") +
            scale_fill_manual(values = catColProf[[catColInd]]) +
            #facet_wrap(~Histotype) +
            # stat_compare_means(comparisons=list(c(cont1, cont2))) +
            scale_y_continuous(expand = expansion(mult = .15)) +
            scale_x_discrete() +
            geom_hline(yintercept=0.2) +
            geom_hline(yintercept=0.7) +
            theme(text = element_text(size=18), 
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
                  legend.position = "none",
                  plot.title = element_text(hjust = 0.5)) + 
            facet_wrap(~factor(cpg_class, c("promo", "body")), nrow =1, strip.position = "bottom")
          
          # Change the point size, and shape
          bodyPlot <- ggplot(bodyMat, aes(x=pos, y=Beta, fill=Histotype)) + 
            geom_boxplot() +
            labs(x=paste(tmpRow$chr, tmpRow$start, tmpRow$end, tmpRow$strand, sep="::"), 
                 y="Beta",
                 title=paste(contName, ": HSP-CPG: ", geneSym, sep="")) +
            scale_fill_manual(values = catColProf[[catColInd]]) +
            #facet_wrap(~Histotype) +
            # stat_compare_means(comparisons=list(c(cont1, cont2))) +
            scale_y_continuous(expand = expansion(mult = .15)) +
            scale_x_discrete() +
            geom_hline(yintercept=0.2) +
            geom_hline(yintercept=0.7) +
            theme(text = element_text(size=18), 
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
                  legend.text=element_text(size=12),
                  plot.title = element_text(hjust = 0.5)) + 
            facet_wrap(~factor(cpg_class, c("promo", "body")), nrow =1, strip.position = "bottom")
          if(tmpGeneRow$strand == "-"){
            tmpPlt <- ggarrange(bodyPlot, promoPlot,
                                ncol = 2, nrow = 1,
                                widths = c(1, 0.5))
          }else{
            tmpPlt <- ggarrange(promoPlot, bodyPlot,
                                ncol = 2, nrow = 1,
                                widths = c(0.5, 1))
          }
        }else{
          tmpPlt <- ggplot(mergeDf, aes(x=pos, y=Beta, fill=Histotype)) + 
            geom_boxplot() +
            labs(x=paste(tmpRow$chr, tmpRow$start, tmpRow$end, tmpRow$strand, sep="::"), 
                 y="Beta",
                 title=paste(contName, ": DEG-CPG: ", geneSym, sep="")) +
            scale_fill_manual(values = catColProf[[catColInd]]) +
            #facet_wrap(~Histotype) +
            # stat_compare_means(comparisons=list(c(cont1, cont2))) +
            scale_y_continuous(expand = expansion(mult = .15)) +
            scale_x_discrete() +
            geom_hline(yintercept=0.2) +
            geom_hline(yintercept=0.7) +
            theme(text = element_text(size=18), 
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
                  legend.text=element_text(size=12),
                  plot.title = element_text(hjust = 0.5)) + 
            facet_wrap(~factor(cpg_class, c("promo", "body")), nrow =1, strip.position = "bottom")
        }
        if(pltBool){
          ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
        }
      }
    }
  }
}

makeCpgBoxPlot <- function(inpCpgs, inpBeta, inpAnno, focusGrp="Histotype", grpName=NULL, fileExt=NULL, dirExt=NULL, pltBool=NULL){
  # Function for plotting indivudual CpG sites as boxplots with respect to histotype
  if(is.null(pltBool)){
    pltBool = TRUE
  }
  if(is.null(grpName)){
    grpName <- "CpGs"
  }
  # Function for creating CpG boxplots for gene in provided dataframe
  # Get overlapping CpGs for each gene in input
  tmpCpgs <- inpCpgs
  cpgBetas <- inpBeta[tmpCpgs,]
  colInd <- which(names(inpAnno) %in% focusGrp)
  catColInd <- which(names(catColProf) %in% focusGrp)
  noCpgs <- nrow(cpgBetas)
  message(paste("Creating boxplots for: ", grpName, ": [", noCpgs, "] ", "CpGs", sep=""))
  pltLst <- list()
  for(j in 1:nrow(cpgBetas)){
    message(paste("Plotting gene [", j, "/", noCpgs, "]", sep=""))
    tmpLst <- list()
    cpgRow <- cpgBetas[j,]
    cpgName <- rownames(cpgRow)
    if(is.null(dirExt)){
      outDir <- paste(plotPath,"/cpgSingleBoxPlot/", grpName, "/", cpgName,"/", sep="")
    }else{
      outDir <- paste(plotPath,"/cpgSingleBoxPlot/", grpName, "/", dirExt,"/", cpgName,"/", sep="")
    }
    ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
    if(!is.null(fileExt)){
      outFile <- paste(outDir, cpgName, "_", fileExt, "_","Single_Cpg_BoxPlot.pdf", sep="")
    }else{
      outFile <- paste(outDir, cpgName, "_","Single_Cpg_BoxPlot.pdf", sep="")
    }
    cpgRow <- t(cpgRow)
    cpgRow <- as.data.frame(cpgRow)
    # Create matching column for samples in input-annotation
    cpgRow$Histotype <- inpAnno$Histotype[match(row.names(cpgRow), inpAnno$barcode)]
    colnames(cpgRow)[1] <- "Beta"
    tmpPlt <- ggplot(cpgRow, aes(x=Histotype, y=Beta, fill=Histotype)) + 
      geom_boxplot() +
      labs(x=cpgName, 
           y="Beta",
           title=paste("Beta values for CPG: ", cpgName, sep="")) +
      scale_fill_manual(values = catColProf[[catColInd]]) +
      #facet_wrap(~Histotype) +
      # stat_compare_means(comparisons=list(c(cont1, cont2))) +
      scale_y_continuous(expand = expansion(mult = .15)) +
      scale_x_discrete() +
      geom_hline(yintercept=0.2) +
      geom_hline(yintercept=0.7) +
      theme(text = element_text(size=18), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
            legend.text=element_text(size=12),
            plot.title = element_text(hjust = 0.5)) 
    if(pltBool){
      ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
    }
    pltLst[[j]] <- tmpPlt
    names(pltLst)[j] <- cpgName
  }
  return(pltLst)
}

makeCpgDensityPlot <- function(inpCpgs, inpBeta, inpAnno, focusGrp="Histotype", grpName=NULL, fileExt=NULL, dirExt=NULL, pltBool=NULL){
  # Function for plotting indivudual CpG sites as boxplots with respect to histotype
  if(is.null(pltBool)){
    pltBool = TRUE
  }
  if(is.null(grpName)){
    grpName <- "CpGs"
  }
  # Function for creating CpG boxplots for gene in provided dataframe
  # Get overlapping CpGs for each gene in input
  tmpCpgs <- inpCpgs
  cpgBetas <- inpBeta[tmpCpgs,]
  colInd <- which(colnames(inpAnno) %in% focusGrp)
  catColInd <- which(names(catColProf) %in% focusGrp)
  noCpgs <- nrow(cpgBetas)
  # message(paste("Creating boxplots for: ", grpName, ": [", noCpgs, "] ", "CpGs", sep=""))
  pltLst <- list()
  for(j in 1:nrow(cpgBetas)){
    message(paste("Plotting gene [", j, "/", noCpgs, "]", sep=""))
    tmpLst <- list()
    cpgRow <- cpgBetas[j,]
    cpgName <- rownames(cpgRow)
    if(is.null(dirExt)){
      outDir <- paste(plotPath,"/cpgDensityPlot/", grpName, "/", cpgName,"/", sep="")
    }else{
      outDir <- paste(plotPath,"/cpgDensityPlot/", grpName, "/", dirExt,"/", cpgName,"/", sep="")
    }
    ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
    if(!is.null(fileExt)){
      outFile <- paste(outDir, cpgName, "_", fileExt, "_","cpgDensityPlot.pdf", sep="")
    }else{
      outFile <- paste(outDir, cpgName, "_","cpgDensityPlot.pdf", sep="")
    }
    cpgRow <- t(cpgRow)
    cpgRow <- as.data.frame(cpgRow)
    colnames(cpgRow)[1] <- "Beta"
    # Create matching column for samples in input-annotation
    cpgRow$Histotype <- inpAnno$Histotype[match(row.names(cpgRow), inpAnno$barcode)]
    statDf <- data.frame(matrix(nrow=4, ncol=3))
    rownames(statDf) <- names(table(cpgRow$Histotype))
    colnames(statDf) <- c("Mean", "Median", "TriMean")
    for(k in 1:length(table(cpgRow$Histotype))){
      tmpHist <- names(table(cpgRow$Histotype))[k]
      tmpBetas <- cpgRow[cpgRow$Histotype %in% tmpHist,]
      ariMean <- mean(tmpBetas$Beta)
      triMean <- litteR::trimean(tmpBetas$Beta)
      med <- median(tmpBetas$Beta)
      statDf$Mean[k] <- ariMean
      statDf$Median[k] <- med
      statDf$TriMean[k] <- triMean
    }
    statDf <- rownames_to_column(statDf, "Histotype")
    statDfLong <- data.frame(statDf %>% pivot_longer(cols=c("Mean", "Median", "TriMean"),
                                                     names_to='StatType',
                                                     values_to='Vals'))
    pltDf <- full_join(cpgRow, statDf, by = c("Histotype"))
    
    tmpPlt <- ggplot(pltDf, aes(x=Beta, group=Histotype, fill=Histotype, after_stat(scaled))) + 
      geom_density() +
      labs(x=cpgName, 
           y="Density Estimate",
           title=paste("Beta values for CPG: ", cpgName, sep="")) +
      scale_fill_manual(values = catColProf[[catColInd]]) +
      scale_y_continuous(expand = expansion(mult = .15)) +
      scale_x_continuous() +
      theme(text = element_text(size=18), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
            legend.text=element_text(size=12),
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~factor(Histotype, c(names(table(cpgRow$Histotype)))), nrow =2, strip.position = "bottom") +
      #geom_vline(data=statDf, aes(xintercept=Mean, colour="Mean"), linewidth = 1) +
      #geom_vline(data=statDf, aes(xintercept=Median, colour="Median"), linewidth =1) +
      geom_vline(data=statDfLong, aes(xintercept=Vals, colour=StatType), linewidth = 1) + 
      scale_color_manual(values = c(Mean = catColProf$Sur_Grp[[1]],
                                    Median = catColProf$Sur_Grp[[2]], 
                                    TriMean = catColProf$Sur_Grp[[3]]),
                         labels = unique(statDfLong$StatType))
    
    
    if(pltBool){
      ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
    }
    pltLst[[j]] <- tmpPlt
    names(pltLst)[j] <- cpgName
  }
  return(pltLst)
}

makeGeneLstCpgBoxPlot <- function(inpGenes, inpBeta, inpPheno, cpgInp, geneInfInp, focusGrp="Histotype", grpName=NULL, fileExt=NULL, dirExt=NULL, pltBool=NULL, promoBool = NULL){
  if(is.null(pltBool)){
    pltBool = TRUE
  }
  if(is.null(grpName)){
    grpName <- "All_Phenotypes"
  }
  if(!is.null(promoBool)){
    promoBool <- "PROMOTER"
  }
  
  # Function for creating CpG boxplots for gene in provided dataframe
  # Get overlapping CpGs for each gene in input
  tmpGenes <- inpGenes
  noGenes <- nrow(tmpGenes)
  # Get CpG's overlapping genes 
  cpgGeneOverlap <- makeCpgGeneOverlap_MULT(inpGenes = tmpGenes$ensembl_gene_id, 
                                            cpgInp = cpgInp, 
                                            inpGeneInf = geneInfInp, 
                                            type = promoBool)
  # Convert list of overlapping CpGs into dataframe
  cpgGeneDf <- makeGeneCpGLocDataframe(cpgGeneOverlap)
  # Get list of dataframes containing beta-values
  cpgBetas <- getGeneCpgBeta(geneCpgInp = cpgGeneDf, 
                             inpBeta = inpBeta, 
                             inpPheno = inpPheno,
                             cpgInf = cpgInp, 
                             allBool = TRUE)
  colInd <- which(names(inpPheno) %in% focusGrp)
  catColInd <- which(names(catColProf) %in% focusGrp)
  message(paste("Creating boxplots for: ", grpName, ": [", noGenes, "] ", "Genes", sep=""))
  pltLst <- list()
  for(j in 1:noGenes){
    message(paste("Plotting gene [", j, "/", noGenes, "]", sep=""))
    tmpLst <- list()
    geneRow <- tmpGenes[j,]
    ensName <- geneRow$ensembl_gene_id
    geneSym <- geneRow$external_gene_name
    cpgPos <- cpgBetas[[ensName]]
    # Make sure that ALL overlapping CpGs for gene are on the correct strand
    # cpgPos <- cpgPos[which(cpgPos$strand == geneRow$strand),]
    if(is.null(cpgPos)){
      message(paste("No overlapping CpGs found for gene:", geneSym, "proceeding to next gene",sep=" "))
      next
    }else if(nrow(cpgPos) < 2){
      message(paste("Insufficient number of overlapping CpGs:", nrow(cpgPos), "for gene", geneSym, "proceeding to next gene",sep=" "))
      next
    }else{
      if(is.null(dirExt)){
        outDir <- paste(plotPath, "cpgBoxPlot/", grpName, "/", geneSym,"/", sep="")
      }else{
        outDir <- paste(plotPath, "cpgBoxPlot/", grpName, "/", dirExt,"/",geneSym,"/", sep="")
      }
      ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
      if(!is.null(fileExt)){
        outFile <- paste(outDir, geneSym, "_", fileExt, "_","Gene_Cpg_BoxPlot.pdf", sep="")
      }else{
        outFile <- paste(outDir, geneSym, "_","Gene_Cpg_BoxPlot.pdf", sep="")
      }
      # Filter our plot-vector
      skipVec <- c("Row.names", "chr", "pos", "strand", "start", "end")
      pltDf <- as.data.frame(t(cpgPos[,!colnames(cpgPos) %in% skipVec]))
      colnames(pltDf) <- cpgPos$Row.names
      # Create matching column for samples in input-annotation
      inpPheno$matchCol <- paste(inpPheno$Histotype, inpPheno$barcode , sep="_")
      pltDf$Histotype <- inpPheno$Histotype[match(row.names(pltDf), inpPheno$matchCol)]
      mergeDf <- tidyr::gather(pltDf, key="CpG", value="Beta",  colnames(pltDf)[!colnames(pltDf) %in% "Histotype"], factor_key=TRUE)
      # Add CpG position as factor
      mergeDf$pos <- cpgPos$pos[match(mergeDf$CpG , cpgPos$Row.names)]
      mergeDf$cpg_class <- ifelse(between(mergeDf$pos, 
                                          geneRow$promo_start, 
                                          geneRow$promo_end), "promo", "body")
      # Sort dataframe, then add an ordered position vector of factors to represent the gene coordinates
      mergeDf <- mergeDf[order(factor(mergeDf$pos, levels = levels(factor(mergeDf$pos)), ordered = TRUE)),]
      mergeDf$pos <- as.factor(mergeDf$pos)
      # Convert to factor, sort separately
      # mergeDf$pos_ordered <- factor(mergeDf$pos, levels = levels(factor(mergeDf$pos)), ordered = TRUE)
      # mergeDf$pos <- as.factor(mergeDf$pos)
      
      if(length(which(mergeDf$cpg_class == "promo")) > 0 & length(which(mergeDf$cpg_class == "body")) > 0){
        # If we have a promoter entry, we create two separate plots to scale them
        promoMat <- mergeDf[which(mergeDf$cpg_class == "promo"),]
        bodyMat <- mergeDf[which(!mergeDf$cpg_class == "promo"),]
        # Change the point size, and shape
        promoPlot <- ggplot(promoMat, aes(x=pos, y=Beta, fill=Histotype)) + 
          geom_boxplot() +
          labs(x=paste(geneRow$chr, geneRow$start, geneRow$end, geneRow$strand, sep="::"), 
               y="Beta",
               title="") +
          scale_fill_manual(values = catColProf[[catColInd]]) +
          #facet_wrap(~Histotype) +
          # stat_compare_means(comparisons=list(c(cont1, cont2))) +
          scale_y_continuous(limits = c(0, 1), expand = expansion(mult = .15)) +
          scale_x_discrete() +
          geom_hline(yintercept=0.2) +
          geom_hline(yintercept=0.7) +
          theme(text = element_text(size=18), 
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
                legend.position = "none",
                plot.title = element_text(hjust = 0.5)) + 
          facet_wrap(~factor(cpg_class, c("promo", "body")), nrow =1, strip.position = "bottom")
        
        # Change the point size, and shape
        bodyPlot <- ggplot(bodyMat, aes(x=pos, y=Beta, fill=Histotype)) + 
          geom_boxplot() +
          labs(x=paste(geneRow$chr, geneRow$start, geneRow$end, geneRow$strand, sep="::"), 
               y="Beta",
               title=paste("Gene-CPG: ", geneSym, sep="")) +
          scale_fill_manual(values = catColProf[[catColInd]]) +
          #facet_wrap(~Histotype) +
          # stat_compare_means(comparisons=list(c(cont1, cont2))) +
          scale_y_continuous(limits = c(0, 1), expand = expansion(mult = .15)) +
          scale_x_discrete() +
          geom_hline(yintercept=0.2) +
          geom_hline(yintercept=0.7) +
          theme(text = element_text(size=18), 
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
                legend.text=element_text(size=12),
                plot.title = element_text(hjust = 0.5)) + 
          facet_wrap(~factor(cpg_class, c("promo", "body")), nrow =1, strip.position = "bottom")
        if(geneRow$strand == "-"){
          tmpPlt <- ggarrange(bodyPlot, promoPlot,
                              ncol = 2, nrow = 1,
                              widths = c(1, 0.5))
        }else{
          tmpPlt <- ggarrange(promoPlot, bodyPlot,
                              ncol = 2, nrow = 1,
                              widths = c(0.5, 1))
        }
      }else{
        tmpPlt <- ggplot(mergeDf, aes(x=pos, y=Beta, fill=Histotype)) + 
          geom_boxplot() +
          labs(x=paste(geneRow$chr, geneRow$start, geneRow$end, geneRow$strand, sep="::"), 
               y="Beta",
               title=paste("DEG-CPG: ", geneSym, sep="")) +
          scale_fill_manual(values = catColProf[[catColInd]]) +
          #facet_wrap(~Histotype) +
          # stat_compare_means(comparisons=list(c(cont1, cont2))) +
          scale_y_continuous(limits = c(0, 1), expand = expansion(mult = .15)) +
          scale_x_discrete() +
          geom_hline(yintercept=0.2) +
          geom_hline(yintercept=0.7) +
          theme(text = element_text(size=18), 
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
                legend.text=element_text(size=12),
                plot.title = element_text(hjust = 0.5)) + 
          facet_wrap(~factor(cpg_class, c("promo", "body")), nrow =1, strip.position = "bottom")
      }
      if(pltBool){
        ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
      }
      pltLst[[j]] <- tmpPlt
      names(pltLst)[j] <- geneSym
    }
  }
  return(pltLst)
}

makeGeneLstCpgDotPlot <- function(inpGenes, inpBeta, inpPheno, cpgInp, geneInfInp, focusGrp="Histotype", 
                                  grpName=NULL, fileExt=NULL, dirExt=NULL, pltBool=NULL, promoBool = NULL){
  if(is.null(pltBool)){
    pltBool = TRUE
  }
  if(is.null(grpName)){
    grpName <- "All_Phenotypes"
  }
  if(!is.null(promoBool)){
    promoBool <- "PROMOTER"
  }
  # Function for creating CpG boxplots for gene in provided dataframe
  # Get overlapping CpGs for each gene in input
  tmpGenes <- inpGenes
  noGenes <- nrow(tmpGenes)
  # Get CpG's overlapping genes 
  cpgGeneOverlap <- makeCpgGeneOverlap_MULT(inpGenes = tmpGenes$ensembl_gene_id, 
                                            cpgInp = cpgInp, 
                                            inpGeneInf = geneInfInp, 
                                            type = promoBool)
  # Convert list of overlapping CpGs into dataframe
  cpgGeneDf <- makeGeneCpGLocDataframe(cpgGeneOverlap)
  # Get list of dataframes containing beta-values
  cpgBetas <- getGeneCpgBeta(geneCpgInp = cpgGeneDf, 
                             inpBeta = inpBeta, 
                             inpPheno = inpPheno,
                             cpgInf = cpgInp, 
                             allBool = TRUE)
  colInd <- which(names(inpPheno) %in% focusGrp)
  catColInd <- which(names(catColProf) %in% focusGrp)
  message(paste("Creating dotplots for: ", grpName, ": [", noGenes, "] ", "Genes", sep=""))
  pltLst <- list()
  for(j in 1:noGenes){
    # print(j)
    message(paste("Plotting gene [", j, "/", noGenes, "]", sep=""))
    tmpLst <- list()
    geneRow <- tmpGenes[j,]
    ensName <- geneRow$ensembl_gene_id
    geneSym <- geneRow$external_gene_name
    if(geneSym == ""){
      geneSym <- ensName
    }
    cpgPos <- cpgBetas[[ensName]]
    # Make sure that ALL overlapping CpGs for gene are on the correct strand
    # cpgPos <- cpgPos[which(cpgPos$strand == geneRow$strand),]
    if(is.null(cpgPos)){
      message(paste("No overlapping CpGs found for gene:", geneSym, "proceeding to next gene",sep=" "))
      next
    }else if(nrow(cpgPos) < 2){
      message(paste("Insufficient number of overlapping CpGs:", nrow(cpgPos), "for gene", geneSym, "proceeding to next gene",sep=" "))
      next
    }else{
      if(is.null(dirExt)){
        outDir <- paste(plotPath, "cpgDotPlot/", grpName, "/", geneSym,"/", sep="")
      }else{
        outDir <- paste(plotPath, "cpgDotPlot/", grpName, "/", dirExt,"/",geneSym,"/", sep="")
      }
      ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
      if(!is.null(fileExt)){
        outFile <- paste(outDir, geneSym, "_", fileExt, "_","Gene_Cpg_DotPlot.pdf", sep="")
      }else{
        outFile <- paste(outDir, geneSym, "_","Gene_Cpg_DotPlot.pdf", sep="")
      }
      # Filter our plot-vector
      skipVec <- c("Row.names", "chr", "pos", "strand", "start", "end")
      pltDf <- as.data.frame(t(cpgPos[,!colnames(cpgPos) %in% skipVec]))
      colnames(pltDf) <- cpgPos$Row.names
      # Create matching column for samples in input-annotation
      inpPheno$matchCol <- paste(inpPheno$Histotype, inpPheno$barcode , sep="_")
      pltDf$Histotype <- inpPheno$Histotype[match(row.names(pltDf), inpPheno$matchCol)]
      mergeDf <- tidyr::gather(pltDf, 
                               key="CpG", 
                               value="Beta",  
                               colnames(pltDf)[1]:colnames(pltDf)[ncol(pltDf)-1], 
                               factor_key=TRUE)
      # Add CpG position as factor
      mergeDf$pos <- cpgPos$pos[match(mergeDf$CpG , cpgPos$Row.names)]
      mergeDf$cpg_class <- ifelse(between(mergeDf$pos, 
                                          geneRow$promo_start, 
                                          geneRow$promo_end), "promo", "body")
      # Sort dataframe, then add an ordered position vector of factors to represent the gene coordinates
      mergeDf <- mergeDf[order(factor(mergeDf$pos, levels = levels(factor(mergeDf$pos)), ordered = TRUE)),]
      mergeDf$pos <- as.factor(mergeDf$pos)
      # Convert to factor, sort separately
      # mergeDf$pos_ordered <- factor(mergeDf$pos, levels = levels(factor(mergeDf$pos)), ordered = TRUE)
      # mergeDf$pos <- as.factor(mergeDf$pos)
      
      if(length(which(mergeDf$cpg_class == "promo")) > 0 & length(which(mergeDf$cpg_class == "body")) > 0){
        # If we have a promoter entry, we create two separate plots to scale them
        promoMat <- mergeDf[which(mergeDf$cpg_class == "promo"),]
        bodyMat <- mergeDf[which(!mergeDf$cpg_class == "promo"),]
        # Create median matrix
        promoMed <- data.frame(matrix(nrow=0, ncol=ncol(promoMat)))
        bodyMed <- data.frame(matrix(nrow=0, ncol=ncol(promoMat)))
        colnames(promoMed) <- colnames(promoMat)
        colnames(bodyMed) <- colnames(promoMed)
        for(k in 1:length(unique(promoMat$Histotype))){
          medH <- unique(promoMat$Histotype)[k]
          promoRows <- promoMat[which(promoMat$Histotype %in% medH), ]  
          # bodyRows <- promoMat[which(bodyMat$Histotype %in% medH), ] 
          bodyRows <- bodyMat[which(bodyMat$Histotype %in% medH), ]
          
          pBetaMed <- as.data.frame(promoRows[,c("CpG", "Beta")] %>% 
                                      group_by(CpG) %>%  
                                      summarise_all(.funs = median))
          bBetaMed <- as.data.frame(bodyRows[,c("CpG", "Beta")] %>% 
                                      group_by(CpG) %>%  
                                      summarise_all(.funs = median))
          for(l in 1:nrow(pBetaMed)){
            promoMed[nrow(promoMed)+1,] <- c(medH, as.character(pBetaMed$CpG)[l], as.numeric(as.character(pBetaMed$Beta[l])),
                                             as.numeric(as.character(promoRows$pos[which(promoRows$CpG %in% as.character(pBetaMed$CpG)[l])][1])), "promo" )
          }
          for(m in 1:nrow(bBetaMed)){
            bodyMed[nrow(bodyMed)+1,] <- c(medH, as.character(bBetaMed$CpG)[m],as.numeric(as.character(bBetaMed$Beta[m])), 
                                           as.numeric(as.character(bodyRows$pos[which(bodyRows$CpG %in% as.character(bBetaMed$CpG)[m])][1])), "body" )
          }
        }
        
        # Change the point size, and shape
        promoPlot <- ggplot(promoMed, aes(x=as.factor(pos), y=as.numeric(Beta), col=Histotype)) + 
          geom_point(size=4) +
          labs(x=paste(geneRow$chr, geneRow$start, geneRow$end, geneRow$strand, sep="::"), 
               y="Beta (trimean)",
               title="") +
          scale_color_manual(values = catColProf[[catColInd]]) +
          #facet_wrap(~Histotype) +
          # stat_compare_means(comparisons=list(c(cont1, cont2))) +
          scale_y_continuous(limits = c(0, 1), expand = expansion(mult = .15)) +
          scale_x_discrete() +
          geom_hline(yintercept=0.2) +
          geom_hline(yintercept=0.7) +
          theme(text = element_text(size=22), 
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
                axis.text.y = element_text(size=12),
                legend.position = "none",
                plot.title = element_text(hjust = 0.5)) + 
          facet_wrap(~factor(cpg_class, c("promo", "body")), nrow =1, strip.position = "bottom")
        
        # Change the point size, and shape
        bodyPlot <- ggplot(bodyMed, aes(x=as.factor(pos), y=as.numeric(Beta), col=Histotype)) + 
          geom_point(size=4) +
          labs(x=paste(geneRow$chr, geneRow$start, geneRow$end, geneRow$strand, sep="::"), 
               y="Beta (trimean)",
               title=paste("Gene-CPG: ", geneSym, sep="")) +
          scale_color_manual(values = catColProf[[catColInd]]) +
          #facet_wrap(~Histotype) +
          # stat_compare_means(comparisons=list(c(cont1, cont2))) +
          scale_y_continuous(limits = c(0, 1), expand = expansion(mult = .15)) +
          scale_x_discrete() +
          geom_hline(yintercept=0.2) +
          geom_hline(yintercept=0.7) +
          theme(text = element_text(size=22), 
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
                axis.text.y = element_text(size=12),
                legend.text=element_text(size=18),
                plot.title = element_text(hjust = 0.5)) + 
          facet_wrap(~factor(cpg_class, c("promo", "body")), nrow =1, strip.position = "bottom")
        if(geneRow$strand == "-"){
          tmpPlt <- ggarrange(bodyPlot, promoPlot,
                              ncol = 2, nrow = 1,
                              widths = c(1, 0.5))
        }else{
          tmpPlt <- ggarrange(promoPlot, bodyPlot,
                              ncol = 2, nrow = 1,
                              widths = c(0.5, 1))
        }
      }else{
        medDf <- data.frame(matrix(nrow=0, ncol=ncol(mergeDf)))
        colnames(medDf) <- colnames(mergeDf)
        for(k in 1:length(unique(mergeDf$Histotype))){
          medH <- unique(mergeDf$Histotype)[k]
          tmpRows <- mergeDf[which(mergeDf$Histotype %in% medH), ]  
          betaMed <- as.data.frame(tmpRows[,c("CpG", "Beta")] %>% 
                                     group_by(CpG) %>%  
                                     summarise_all(.funs = litteR::trimean))
          for(l in 1:nrow(betaMed)){
            medDf[nrow(medDf)+1,] <- c(medH, as.character(betaMed$CpG)[l], as.numeric(as.character(betaMed$Beta[l])),
                                       as.numeric(as.character(tmpRows$pos[which(tmpRows$CpG %in% as.character(betaMed$CpG)[l])][1])), 
                                       tmpRows$cpg_class[which(tmpRows$CpG %in% as.character(betaMed$CpG)[l])][1])
          }
        }
        
        tmpPlt <- ggplot(medDf, aes(x=as.factor(pos), y=as.numeric(Beta), col=Histotype)) + 
          geom_point(size=4) +
          labs(x=paste(geneRow$chr, geneRow$start, geneRow$end, geneRow$strand, sep="::"), 
               y="Beta (triMean)",
               title=paste("DEG-CPG: ", geneSym, sep="")) +
          scale_color_manual(values = catColProf[[catColInd]]) +
          scale_y_continuous(limits = c(0, 1), expand = expansion(mult = .15)) +
          scale_x_discrete() +
          geom_hline(yintercept=0.2) +
          geom_hline(yintercept=0.7) +
          theme(text = element_text(size=22), 
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
                axis.text.y = element_text(size=12),
                legend.text=element_text(size=18),
                plot.title = element_text(hjust = 0.5)) + 
          facet_wrap(~factor(cpg_class, c("promo", "body")), nrow =1, strip.position = "bottom")
      }
      if(pltBool){
        ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
      }
      pltLst[[j]] <- tmpPlt
      names(pltLst)[j] <- geneSym
    }
  }
  return(pltLst)
}

makeExtCompCpgBoxPlot <- function(inpCpgGeneLst, inpBetaLst, inpPhenoLst, inp450k, inpEPIC, inpH, inpGeneInf, fileExt=NULL, dirExt=NULL, pltBool=NULL){
  pltLst <- list()
  if(is.null(pltBool)){
    pltBool = TRUE
  }
  # Function for creating CpG boxplots for gene in provided dataframe
  # Get overlapping CpGs for each gene in input
  for(i in 1:length(inpCpgGeneLst)){
    tmpCpg <- inpCpgGeneLst[[i]] 
    tmpEns <- names(inpCpgGeneLst)[i]
    tmpSym <- inpGeneInf$external_gene_name[which(inpGeneInf$ensembl_gene_id %in% tmpEns)]
    if(length(tmpSym) == 0){
      tmpSym <- tmpEns
    }
    # Get beta values for each entry in beta-list
    gBLst <- list()
    for(j in 1:length(inpBetaLst)){
      tmpB <- inpBetaLst[[j]]
      tmpSid <- names(inpBetaLst)[j]
      tmpP <- inpPhenoLst[[tmpSid]]
      filtDf <- makePhenoBetaFilter(tmpB, tmpP)
      tmpB <- filtDf[[1]]
      tmpP <- filtDf[[2]]
      tmpB <- tmpB[, match(tmpP$barcode[which(tmpP$Histotype %in% inpH)], colnames(tmpB))]
      pCpg <- rownames(tmpCpg)[which(rownames(tmpCpg) %in% rownames(tmpB))]
      if(length(pCpg) == 0){
        next()
      }else if(length(pCpg) == 1){
        tmpCpgB <- tmpB[pCpg,]
        tmpCpgB <- data.frame(tmpCpgB)
        if(nrow(tmpCpgB) > length(pCpg)){
          tmpCpgB <- t(tmpCpgB)
        }
        rownames(tmpCpgB) <- pCpg
      }else{
        tmpCpgB <- tmpB[pCpg,]
      }
      tmpCpgB <- na.omit(tmpCpgB)
      gBLst[[j]] <- tmpCpgB
      names(gBLst)[j] <- tmpSid
    }
    # Remove empty entries from list (i.e. empty dataframes)
    gBLst <-  gBLst[lengths(gBLst) != 0]
    if(is.null(dirExt)){
      outDir <- paste(plotPath,"extCompCpgBoxPlot/", inpH, "/", tmpSym,"/", sep="")
    }else{
      outDir <- paste(plotPath,"extCompCpgBoxPlot/", inpH, "/", dirExt,"/",tmpSym,"/", sep="")
    }
    ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
    if(!is.null(fileExt)){
      outFile <- paste(outDir, tmpSym, "_", fileExt, "_","extCompCpgBoxPlot.pdf", sep="")
    }else{
      outFile <- paste(outDir, tmpSym, "_","extCompCpgBoxPlot.pdf", sep="")
    }
    # Extend so that vectors are of equal length by padding to the max length
    maxL <- 0
    for(k in 1:length(gBLst)){
      tmpB <- gBLst[[k]]
      if(ncol(tmpB) > maxL){
        maxL <- ncol(tmpB)
      }
    }
    for(l in 1:length(gBLst)){
      tmpB <- gBLst[[l]]
      tmpSid <- names(gBLst)[l]
      if(ncol(tmpB) <= maxL){
        tmpDiff <- maxL - ncol(tmpB)
        tmpDf <- data.frame(matrix(nrow=nrow(tmpB), ncol=tmpDiff))
        rownames(tmpDf) <- rownames(tmpB)
        # tmpB <- tibble::rownames_to_column(tmpB, "cpg")
        tmpB <- cbind(tmpB, tmpDf)
        tmpB <- data.frame(t(tmpB))
        tmpB$SID <- tmpSid
        gBLst[[l]] <- tmpB
      }else{
        tmpB$SID <- tmpSid
      }
    }
    # Merge list of dataframes into one dataframe
    mergeDf <- bind_rows(gBLst, .id = "SID")
    mergeDf <- mergeDf[,!colnames(mergeDf) %like% "NA"]
    # mergeDf <- na.omit(mergeDf)
    # mergeDf$cpg_class <- "promo"
    mergeDf <- tidyr::gather(mergeDf, 
                             key="CpG", 
                             value="Beta",  
                             colnames(mergeDf)[!colnames(mergeDf) %in% "SID"], 
                             factor_key=TRUE)
    cpgOrder <-  rownames(tmpCpg)
    mergeDf$cpg_class <- "promo"
    mergeDf$pos <- tmpCpg$pos[match(mergeDf$CpG , rownames(tmpCpg))]
    # Sort dataframe, then add an ordered position vector of factors to represent the gene coordinates
    mergeDf <- mergeDf[order(factor(mergeDf$pos, levels = levels(factor(mergeDf$pos)), ordered = TRUE)),]
    mergeDf$pos <- as.factor(mergeDf$pos)
    # Create plot information
    cpgDatVec <- c(paste(tmpCpg$chr[1], min(tmpCpg$start), max(tmpCpg$end), tmpCpg$strand[1], sep="::"))
    colVec <- viridis::plasma(length(gBLst))
    names(colVec) <- names(gBLst)
    # Change the point size, and shape
    compPlt <- ggplot(mergeDf, aes(x=pos, y=Beta, fill=SID)) + 
      geom_boxplot() +
      labs(x= cpgDatVec, 
           y="Beta",
           title=paste("Comparative boxplot for promoter of gene:", tmpSym, sep=" ")) +
      scale_fill_manual(values = colVec) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = .15)) +
      scale_x_discrete() +
      geom_hline(yintercept=0.2) +
      geom_hline(yintercept=0.7) +
      theme(text = element_text(size=18), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~factor(cpg_class, c("promo", "body")), nrow = 1, strip.position = "bottom")
    if(pltBool){
      ggsave(outFile, plot=compPlt, width=30, height=20, units = "cm")
    }
    pltLst[[i]] <- compPlt 
    names(pltLst)[i] <- tmpSym
  }
  return(pltLst)
}

makeExtCompCpgDotPlot <- function(inpCpgGeneLst, inpBetaLst, inpPhenoLst, inp450k, inpEPIC, inpH, inpGeneInf, fileExt=NULL, dirExt=NULL, pltBool=NULL, promoBool= NULL){
  pltLst <- list()
  if(is.null(pltBool)){
    pltBool = TRUE
  }
  # Function for creating CpG boxplots for gene in provided dataframe
  # Get overlapping CpGs for each gene in input
  for(i in 1:length(inpCpgGeneLst)){
    tmpCpg <- inpCpgGeneLst[[i]] 
    tmpCpg <- tmpCpg[which(rownames(tmpCpg) %in% rownames(inpBetaLst$TRAIN)),]
    tmpEns <- names(inpCpgGeneLst)[i]
    tmpSym <- inpGeneInf$external_gene_name[which(inpGeneInf$ensembl_gene_id %in% tmpEns)]
    if(length(tmpSym) == 0){
      tmpSym <- tmpEns
    }
    # Get beta values for each entry in beta-list
    gBLst <- list()
    for(j in 1:length(inpBetaLst)){
      tmpB <- inpBetaLst[[j]]
      tmpSid <- names(inpBetaLst)[j]
      tmpP <- inpPhenoLst[[tmpSid]]
      filtDf <- makePhenoBetaFilter(tmpB, tmpP)
      tmpB <- filtDf[[1]]
      tmpP <- filtDf[[2]]
      tmpB <- tmpB[, match(tmpP$barcode[which(tmpP$Histotype %in% inpH)], colnames(tmpB))]
      pCpg <- rownames(tmpCpg)[which(rownames(tmpCpg) %in% rownames(tmpB))]
      if(length(pCpg) == 0){
        next()
      }else if(length(pCpg) == 1){
        tmpCpgB <- tmpB[pCpg,]
        tmpCpgB <- data.frame(tmpCpgB)
        if(nrow(tmpCpgB) > length(pCpg)){
          tmpCpgB <- t(tmpCpgB)
        }
        rownames(tmpCpgB) <- pCpg
      }else{
        tmpCpgB <- tmpB[pCpg,]
      }
      tmpCpgB <- na.omit(tmpCpgB)
      gBLst[[j]] <- tmpCpgB
      names(gBLst)[j] <- tmpSid
    }
    # Remove empty entries from list (i.e. empty dataframes)
    gBLst <-  gBLst[lengths(gBLst) != 0]
    if(is.null(dirExt)){
      outDir <- paste(plotPath,"/extCompCpgDotPlot/", inpH, "/", tmpSym,"/", sep="")
    }else{
      outDir <- paste(plotPath,"/extCompCpgDotPlot/", inpH, "/", dirExt,"/",tmpSym,"/", sep="")
    }
    ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
    if(!is.null(fileExt)){
      outFile <- paste(outDir, tmpSym, "_", fileExt, "_","extCompCpgDotPlot.pdf", sep="")
    }else{
      outFile <- paste(outDir, tmpSym, "_","extCompCpgDotPlot.pdf", sep="")
    }
    # Extend so that vectors are of equal length by padding to the max length
    maxL <- 0
    for(k in 1:length(gBLst)){
      tmpB <- gBLst[[k]]
      if(ncol(tmpB) > maxL){
        maxL <- ncol(tmpB)
      }
    }
    for(k in 1:length(gBLst)){
      tmpB <- gBLst[[k]]
      tmpSid <- names(gBLst)[k]
      if(ncol(tmpB) <= maxL){
        tmpDiff <- maxL - ncol(tmpB)
        tmpDf <- data.frame(matrix(nrow=nrow(tmpB), ncol=tmpDiff))
        rownames(tmpDf) <- rownames(tmpB)
        # tmpB <- tibble::rownames_to_column(tmpB, "cpg")
        tmpB <- cbind(tmpB, tmpDf)
        tmpB <- data.frame(t(tmpB))
        tmpB$SID <- tmpSid
        gBLst[[k]] <- tmpB
      }else{
        tmpB$SID <- tmpSid
      }
    }
    # Merge list of dataframes into one dataframe
    mergeDf <- bind_rows(gBLst, .id = "SID")
    mergeDf <- mergeDf[,!colnames(mergeDf) %like% "NA"]
    # Remove rows with NA (i.e. padding, all other rows should have a value for beta)
    mergeDf <- mergeDf[rowSums(is.na(mergeDf[,!colnames(mergeDf) %in% "SID"])) != ncol(mergeDf)-1, ]
    # Create a "long" format dataframe
    mergeDf <- tidyr::gather(mergeDf, key="CpG", value="Beta",  colnames(mergeDf)[!colnames(mergeDf) %in% "SID"], factor_key=TRUE)
    # Add additional information to vector
    cpgOrder <-  rownames(tmpCpg)
    mergeDf$cpg_class <- "promo"
    mergeDf$pos <- tmpCpg$pos[match(mergeDf$CpG , rownames(tmpCpg))]
    # Sort dataframe, then add an ordered position vector of factors to represent the gene coordinates
    mergeDf <- mergeDf[order(factor(mergeDf$pos, levels = levels(factor(mergeDf$pos)), ordered = TRUE)),]
    mergeDf$pos <- as.factor(mergeDf$pos)
    # Create plot information
    cpgDatVec <- c(paste(tmpCpg$chr[1], min(tmpCpg$start), max(tmpCpg$end), tmpCpg$strand[1], sep="::"))
    colVec <- viridis(length(gBLst), option = "B")
    names(colVec) <- names(gBLst)
    # Create new trimean dataframe from data
    promoTriMean <- data.frame(matrix(nrow=0, ncol=ncol(mergeDf)))
    colnames(promoTriMean) <- colnames(mergeDf)
    for(l in 1:length(unique(mergeDf$SID))){
      medS <- unique(mergeDf$SID)[l]
      tmpRows <- mergeDf[which(mergeDf$SID %in% medS), ]  
      tmpRows <- tmpRows[!is.na(tmpRows$Beta),]
      pBetaMean <- as.data.frame(tmpRows[,c("CpG", "Beta")] %>% 
                                   group_by(CpG) %>%  
                                   summarise_all(.funs = litteR::trimean))
      for(m in 1:nrow(pBetaMean)){
        promoTriMean[nrow(promoTriMean)+1,] <- c(medS, 
                                                 as.character(pBetaMean$CpG)[m], 
                                                 as.numeric(as.character(pBetaMean$Beta[m])),
                                                 "promo",
                                                 as.numeric(as.character(tmpRows$pos[which(tmpRows$CpG %in% as.character(pBetaMean$CpG)[m])][1])))
      }
    }
    # Change the point size, and shape
    compPlt <- ggplot(promoTriMean, aes(x=pos, y=as.numeric(Beta), col=SID)) + 
      geom_point(size=4) +
      labs(x= cpgDatVec, 
           y="Beta",
           title=paste("Comparative dotplot of trimean for promoter of gene:", tmpSym, sep=" ")) +
      scale_color_manual(values = colVec) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = .15)) +
      scale_x_discrete() +
      geom_hline(yintercept=0.2) +
      geom_hline(yintercept=0.7) +
      theme(text = element_text(size=22), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
            axis.text.y = element_text(size=12),
            plot.title = element_text(hjust = 0.5),
            legend.text=element_text(size=18)) +
      facet_wrap(~factor(cpg_class, c("promo", "body")), nrow = 1, strip.position = "bottom")
    if(pltBool){
      ggsave(outFile, plot=compPlt, width=30, height=20, units = "cm")
    }
    pltLst[[i]] <- compPlt 
    names(pltLst)[i] <- tmpSym
  }
  return(pltLst)
}

makeMethBetaFreqPlt <- function(inpBeta, inpPheno, inpSid, fileExt=NULL){
  inpBeta <- inpBeta[,colnames(inpBeta) %in% inpPheno$barcode]
  inpPheno <- inpPheno[inpPheno$barcode %in% colnames(inpBeta), ]
  # Function for creating frequency-plot of beta-values for groups in cohort
  #inpBeta <- t(inpBeta)
  # Convert beta-values into long-format for plotting
  betaLong <- inpBeta %>% 
    pivot_longer(
      cols = colnames(inpBeta), 
      names_to = "Sample",
      values_to = "Beta"
    )
  betaLong <- na.omit(betaLong)
  # Add histotype based on barcodes
  betaLong$Histotype <- inpPheno$Histotype[match(betaLong$Sample, inpPheno$barcode)]
  # Create color-vector 
  colVec <- viridis(length(table(inpPheno$Histotype)))
  names(colVec) <- names(table(inpPheno$Histotype))
  # Plot beta-density
  message(paste("Plotting beta densities for:", inpSid, sep=" "))
  tmpPlt <- ggplot(data = betaLong,
                   aes(x=Beta, y=after_stat(scaled), colour=Histotype)) +  
    #geom_line(stat = 'Beta') + 
    geom_density(alpha = 0.1, linewidth = 1.25) +
    # scale_fill_manual(values = colVec) +
    stat_density(geom="line", position="identity") +
    scale_colour_manual(values = colVec) + 
    xlim(0, 1) + 
    labs(title=paste("Beta density for: ", inpSid, sep="")) +
    #ggtitle(paste("Expression PCA-plot: ", focusCol, sep="")) + 
    xlab("Beta") + 
    ylab("Density") + 
    theme(text = element_text(size=18), 
          legend.text=element_text(size=12),
          axis.text.x = element_text(size=16),
          plot.title = element_text(hjust = 0.5))
  if(!is.null(fileExt)){
    outFile <- paste(plotPath,"/", fileExt, "_", inpSid, "_Beta_Density.pdf",sep="")
  }else{
    outFile <- paste(plotPath,"/", inpSid, "_Beta_Density.pdf",sep="")
  }
  ggsave(filename =  outFile, plot= tmpPlt, width=30, height=20, units = "cm")
}

# Function for creating PCA plots  from methylation data (M-values)
makeMethPcaPlot <- function(inpBeta, inpPheno, focusCol = "Histotype", nameBool=NULL, noProbes=NULL){
  if(is.null(noProbes)){
    noProbes <- 500
  }
  inpPheno <- inpPheno[inpPheno$barcode %in% colnames(inpBeta),]
  inpBeta <- inpBeta[,colnames(inpBeta) %in% inpPheno$barcode]
  inpBeta <- as.matrix(inpBeta)
  bestInd <- head(order(-rowVars(inpBeta)), noProbes)
  bestBeta <- inpBeta[bestInd,]
  # Match phenotypic data to expression data
  inpPheno <- inpPheno[match(colnames(bestBeta),inpPheno$barcode),]
  inpPheno <- inpPheno[!is.na(inpPheno$barcode), ]
  # Perform the principal component analysis
  pca <- stats::prcomp(t(bestBeta[, matrixStats::colSds(as.matrix(bestBeta), na.rm = T) > 0]))
  varRatio <- round(100*pca$sdev^2/sum(pca$sdev^2),1) 
  plot_df <- as.data.frame(pca$x)
  # Centres the plot on the data points while keeping the plot square
  limits_x <- c(min(plot_df$PC1), max(plot_df$PC1))
  limits_y <- c(min(plot_df$PC2), max(plot_df$PC2))
  width_x <- limits_x[2] - limits_x[1]
  width_y <- limits_y[2] - limits_y[1]
  padding <- abs(width_y - width_x) / 2
  if (width_y > width_x) {
    limits_x <- limits_x + c(-padding, padding)
  } else {
    limits_y <- limits_y + c(-padding, padding)
  }
  # Set up phenotypic representation of data
  plot_df$Groups <- inpPheno$Histotype
  aesthetic_map <- ggplot2::aes(plot_df$PC1, plot_df$PC2, colour = plot_df$Groups)
  colVec <- viridis(length(table(plot_df$Groups)))
  names(colVec) <- names(table(plot_df$Groups))
  # If names are in catColProf, we replace these
  colVec <- replace(colVec, names(catColProf$Histotype), catColProf$Histotype) 
  
  if(!is.null(nameBool)){
    outFile <- paste(plotPath, "/",nameBool, "_PCA_Plot_Meth_",focusCol, ".pdf", sep="")
    tmpTitle <- paste(nameBool, ": Methylation PCA-plot", sep="")
  }else{
    outFile <- paste(plotPath, "/",focusCol, "_PCA_Plot_Meth.pdf", sep="")
    tmpTitle <- paste(focusCol, ": Methylation PCA-plot", sep="")
  }
  tmpPlt <- ggplot2::ggplot(plot_df, x = plot_df$PC1, y = plot_df$PC2, aesthetic_map) +
    ggplot2::geom_point(size=4) + 
    ggplot2::coord_fixed(ratio = 1, xlim = limits_x, ylim = limits_y) +
    labs(title=tmpTitle,
         subtitle=paste("Top:", noProbes, "variable probes", sep=" ")) +
    #ggtitle(paste("Expression PCA-plot: ", focusCol, sep="")) + 
    xlab(paste0("PC1, VarExp: ", varRatio[1], "%")) + 
    ylab(paste0("PC2, VarExp: ", varRatio[2], "%")) + 
    labs(colour=focusCol) +
    theme(text = element_text(size=18), 
          legend.text=element_text(size=12),
          axis.text.x = element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colVec) 
  ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
}

makeTopHMap <- function(inpExp, inpAnno, noGenes = NULL, fileExt = NULL, pltBool = NULL, titleVal = NULL, nameVal=NULL, rowAnnInpDf=NULL, heatBool = NULL){
  # Produces a heatmap showing the top n most variable genes and histotype dendogram 
  if(is.null(noGenes)){
    noGenes <- 500
  }
  if(is.null(titleVal)){
    titleVal <- "Variable Probes"
  }
  if(is.null(nameVal)){
    nameVal <- "Z-score"
  }
  if(is.null(pltBool)){
    pltBool = TRUE
  }
  if(is.null(fileExt)){
    fileExt <- ""
  }
  inpExp <- inpExp[,which(colnames(inpExp) %in% inpAnno$barcode)]
  inpExp <- inpExp %>% mutate_if(is.character, as.numeric)
  inpExp[is.na(inpExp)] <- 0
  inpExp <- as.matrix(inpExp)
  bestInd <- head(order(-rowVars(inpExp)), noGenes)
  bestExp <- inpExp[bestInd,]
  if(!is.null(heatBool)){
    heat <- bestExp
    # Remove zero-variance rows
    heat <- heat[which(!is.nan(heat[,1])),]
    zInts <- c(max(heat, na.rm = TRUE), min(heat, na.rm=TRUE))
  }else{
    # Scale VST for plotting
    heat <- t(scale(t(bestExp)))
    # Remove zero-variance rows
    heat <- heat[which(!is.nan(heat[,1])),]
    if(abs(round(max(heat, na.rm = TRUE))) > abs(round(min(heat, na.rm = TRUE)))){
      zInts <- c(abs(round(min(heat, na.rm = TRUE))), round(min(heat, na.rm=TRUE)))
    }else{
      zInts <- c(round(max(heat, na.rm = TRUE)), -round(max(heat, na.rm=TRUE)))
    }
  }
  # Set breaks for colormap, should be equal to number of datapoints in color-profile
  myBreaks <- seq(zInts[1], zInts[2], length.out = 15)
  # Grab annotation for samples in DF
  annDf <- data.frame(inpAnno$Histotype[match(colnames(bestExp), inpAnno$barcode)])
  colnames(annDf) <- c("Histotype")
  # Adjust script legends after presence of histotype in data 
  tmpHisto <- names(table(annDf$Histotype))
  # annDf$Sample_ID <- inpAnno$Sample_ID[match(colnames(bestExp), inpAnno$Sample_ID)]
  # colnames(tmpVst) <- annDf$Histotype
  # Create annotation object for heatmap
  colAnn <- HeatmapAnnotation(
    df = annDf,
    # 'col' (samples) or 'row' (gene) annotation?
    which = 'col',
    # default colour for any NA values in the annotation data-frame, 'ann'
    na_col = 'white',
    col = catColProf,
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    show_annotation_name = FALSE,
    gap = unit(1, 'mm'),
    annotation_legend_param = list( 
      Histotype = list(
        # number of rows across which the legend will be arranged
        nrow = length(tmpHisto), 
        title = 'Samples',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 10, fontface = 'bold'),
        labels_gp = gpar(fontsize = 10, fontface = 'bold'))
    )
  )
  if(!is.null(rowAnnInpDf)){
    # Adds row-annotation showing which phenotype the rows belong to
    # Grab annotation for samples in DF
    rowAnnDf <- data.frame(rowAnnInpDf$Histotype[match(rownames(heat), rowAnnInpDf$rowID)])
    colnames(rowAnnDf) <- c("Histotype")
    rowAnn <- HeatmapAnnotation(
      df = rowAnnDf,
      # 'col' (samples) or 'row' (gene) annotation?
      which = 'row',
      # default colour for any NA values in the annotation data-frame, 'ann'
      na_col = 'white',
      col = catColProf,
      annotation_height = 0.6,
      annotation_width = unit(1, 'cm'),
      show_annotation_name = FALSE,
      gap = unit(1, 'mm'),
      annotation_legend_param = list( 
        Histotype = list(
          # number of rows across which the legend will be arranged
          nrow = length(table(rowAnnDf)), 
          title = 'Probes',
          title_position = 'topcenter',
          legend_direction = 'vertical',
          title_gp = gpar(fontsize = 10, fontface = 'bold'),
          labels_gp = gpar(fontsize = 10, fontface = 'bold'))
      )
    )
    # message(paste("Creating heatmap for top:", noGenes, "DEGs for contrast:", tmpName, sep=" "))
    # Create complex heatmap for histotype-histotype comparison
    chMap <- ComplexHeatmap::Heatmap(heat,
                                     name = nameVal, 
                                     col = colorRamp2(myBreaks, magma(15, direction = -1)),
                                     na_col = 'white',
                                     # Set settings for legend in plot 
                                     heatmap_legend_param = list(
                                       color_bar = 'continuous',
                                       legend_direction = 'vertical',
                                       legend_width = unit(5, 'cm'),
                                       legend_height = unit(8.0, 'cm'),
                                       legend_gp = gpar(fontsize = 11),
                                       title_gp = gpar(fontsize = 14)),
                                     cluster_columns = TRUE,
                                     cluster_rows = TRUE,
                                     show_column_dend = TRUE,
                                     show_column_names = FALSE,
                                     show_row_names = FALSE,
                                     show_row_dend = FALSE,
                                     # Set title for column space
                                     column_title = paste("Top: ", noGenes, " ", titleVal, sep=""),
                                     column_title_gp = gpar(fontsize = 16, fontface = 'bold'),
                                     column_title_rot = 0,
                                     column_names_gp = gpar(fontsize = 6),
                                     row_names_gp = gpar(fontsize=8),
                                     column_names_max_height = unit(10, 'cm'),
                                     column_dend_height = unit(25,'mm'),
                                     # cluster methods for rows and columns
                                     clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                                     # clustering_distance_columns = "canberra", 
                                     clustering_method_columns = 'ward.D2',
                                     clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                                     clustering_method_rows = 'ward.D2',
                                     top_annotation = colAnn,
                                     left_annotation = rowAnn)
  }else{
    # message(paste("Creating heatmap for top:", noGenes, "DEGs for contrast:", tmpName, sep=" "))
    # Create complex heatmap for histotype-histotype comparison
    chMap <- ComplexHeatmap::Heatmap(heat,
                                     name = nameVal, 
                                     col = colorRamp2(myBreaks,magma(15, direction = -1)),
                                     na_col = 'white',
                                     # Set settings for legend in plot 
                                     heatmap_legend_param = list(
                                       color_bar = 'continuous',
                                       legend_direction = 'vertical',
                                       legend_width = unit(5, 'cm'),
                                       legend_height = unit(8.0, 'cm')),
                                     cluster_columns = TRUE,
                                     show_column_dend = TRUE,
                                     show_column_names = FALSE,
                                     show_row_names = FALSE,
                                     show_row_dend = FALSE,
                                     # Set title for column space
                                     column_title = paste("Top: ", nrow(noGenes), " ", titleVal, sep=""),
                                     column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                                     column_title_rot = 0,
                                     column_names_gp = gpar(fontsize = 6),
                                     row_names_gp = gpar(fontsize=8),
                                     column_names_max_height = unit(10, 'cm'),
                                     column_dend_height = unit(25,'mm'),
                                     # cluster methods for rows and columns
                                     clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                                     # clustering_distance_columns = "canberra", 
                                     clustering_method_columns = 'ward.D2',
                                     clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                                     clustering_method_rows = 'ward.D2',
                                     top_annotation = colAnn)
  }
  if(!fileExt ==""){
    outFile <- paste(plotPath,"/", fileExt, "_Top_", noGenes, "_HeatMap.pdf", sep="")
  }else{
    outFile <- paste(plotPath,"/Top_", noGenes, "_HeatMap.pdf", sep="")
  }
  if(pltBool){
    # Create heatmaps for top genes
    pdf(outFile)
    print(chMap)
    dev.off()
  }
  return(chMap)
}

makeDistBar <- function(inpDist, inpBeta, fileExt=NULL){
  contNames <- rownames(inpDist)
  pltDf <- data.frame(matrix(nrow=0, ncol=3))
  colnames(pltDf) <- c("Histotype","Type", "Percent")
  for(j in 1:nrow(inpDist)){
    tmpRow <- inpDist[j,]
    percNorm <- round(tmpRow$ShapNormal/nrow(inpBeta),2)
    percNNorm <- round(tmpRow$ShapNotNormal/nrow(inpBeta),2)
    percLS <- round(tmpRow$LightSkew/nrow(inpBeta),2)
    percHS <- round(tmpRow$HeavySkew/nrow(inpBeta),2)
    pltDf[nrow(pltDf)+1,] <- c(rownames(tmpRow), "Normal", percNorm)
    pltDf[nrow(pltDf)+1,] <- c(rownames(tmpRow), "Non-Normal", percNNorm)
    pltDf[nrow(pltDf)+1,] <- c(rownames(tmpRow), "LightSkew", percLS)
    pltDf[nrow(pltDf)+1,] <- c(rownames(tmpRow), "HeavySkew", percHS)
  }
  pltDf$Type <- factor(pltDf$Type,  levels=c("Normal", "Non-Normal", "LightSkew", "HeavySkew"))
  pltDf$Histotype <- factor(pltDf$Histotype, levels=c("CCC", "EC", "HGSC", "MC"))
  pltDf$Percent <- as.numeric(pltDf$Percent)
  
  tmpPlt <- ggplot(pltDf, aes(fill=Histotype, y=Percent, x=Histotype)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    xlab("Histotype") +
    ylab("Percent") + 
    ggtitle(paste("Percentage of CpG sites in histotype, (n=", nrow(inpBeta),")" , sep="")) + 
    scale_fill_viridis(discrete = TRUE) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Type, scales = "free_x")
  
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(!is.null(fileExt)){
    outFile <- paste(outDir, fileExt, "_Distribution_Types_Consensus_Barplot.pdf", sep="")
  }else{
    outFile <- paste(outDir, "_Distribution_Types_Consensus_Barplot.pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
}

makeBetaBar <- function(inpDist, inpBeta, fileExt = NULL){
  # Script for plotting the methylation types seen to beta values of genes
  contNames <- rownames(inpDist)
  pltDf <- data.frame(matrix(nrow=0, ncol=3))
  colnames(pltDf) <- c("Histotype","Type", "Percent")
  for(j in 1:nrow(inpDist)){
    tmpRow <- inpDist[j,]
    percHypo <- round(tmpRow$Hypo_perc,2)*100
    percHemi <- round(tmpRow$Hemi_perc,2)*100
    percHyper <- round(tmpRow$Hyper_perc,2)*100
    pltDf[nrow(pltDf)+1,] <- c(rownames(tmpRow), "Hypo", percHypo)
    pltDf[nrow(pltDf)+1,] <- c(rownames(tmpRow), "Hemi", percHemi)
    pltDf[nrow(pltDf)+1,] <- c(rownames(tmpRow), "Hyper", percHyper)
  }
  pltDf$Type <- factor(pltDf$Type,  levels=c("Hypo", "Hemi", "Hyper"))
  pltDf$Histotype <- factor(pltDf$Histotype, levels=c("CCC", "EC", "HGSC", "MC"))
  pltDf$Percent <- as.numeric(pltDf$Percent)
  
  tmpPlt <- ggplot(pltDf, aes(fill=Histotype, y=Percent, x=Histotype)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") +
    geom_text(aes(label = Percent), size = 8, vjust = -0.5, colour = "black") +
    xlab("Histotype") +
    ylab("Percent") + 
    ggtitle(paste("Methylation type of CpG sites (n=685650)\n for all samples in histotype" , sep="")) + 
    scale_fill_viridis(discrete = TRUE) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Type, scales = "free_x")
  
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(!is.null(fileExt)){
    outFile <- paste(outDir, fileExt, "_Methylation_Types_Consensus_Barplot.pdf", sep="")
  }else{
    outFile <- paste(outDir, "Methylation_Types_Consensus_Barplot.pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt, width=40, height=30, units = "cm")
}

makeDistBarLst <- function(inpDistLst, inpMLst, inpBetaLst){
  # Create list of distribution plots
  pltDf <- data.frame(matrix(nrow=0, ncol=5))
  colnames(pltDf) <- c("SID", "Value", "Histotype","Distribution", "Percent")
  for(i in 1:length(inpDistLst)){
    inpDist <- inpDistLst[[i]]
    tmpSid <- names(inpDistLst)[i]
    inpM <- inpMLst[[tmpSid]]
    tmpBeta <- inpBetaLst[[tmpSid]]
    rowCount <- nrow(tmpBeta)
    for(j in 1:nrow(inpDist)){
      tmpRow <- inpDist[j,]
      percNorm <- round(tmpRow$ShapNormal/rowCount,2)
      percnNorm <- round(tmpRow$ShapNotNormal/rowCount,2)
      percLS <- round(tmpRow$LightSkew/rowCount,2)
      percHS <- round(tmpRow$HeavySkew/rowCount,2)
      pltDf[nrow(pltDf)+1,] <- c(tmpSid, "Beta", rownames(tmpRow), "Normal", percNorm)
      pltDf[nrow(pltDf)+1,] <- c(tmpSid, "Beta", rownames(tmpRow), "NotNormal", percnNorm)
      pltDf[nrow(pltDf)+1,] <- c(tmpSid, "Beta", rownames(tmpRow), "LightSkew", percLS)
      pltDf[nrow(pltDf)+1,] <- c(tmpSid, "Beta", rownames(tmpRow), "HeavySkew", percHS)
    }
    for(k in 1:nrow(inpM)){
      tmpRow <- inpM[k,]
      percNorm <- round(tmpRow$ShapNormal/rowCount,2)
      percnNorm <- round(tmpRow$ShapNotNormal/rowCount,2)
      percLS <- round(tmpRow$LightSkew/rowCount,2)
      percHS <- round(tmpRow$HeavySkew/rowCount,2)
      pltDf[nrow(pltDf)+1,] <- c(tmpSid, "M", rownames(tmpRow), "Normal", percNorm)
      pltDf[nrow(pltDf)+1,] <- c(tmpSid, "M", rownames(tmpRow), "NotNormal", percnNorm)
      pltDf[nrow(pltDf)+1,] <- c(tmpSid, "M", rownames(tmpRow), "LightSkew", percLS)
      pltDf[nrow(pltDf)+1,] <- c(tmpSid, "M", rownames(tmpRow), "HeavySkew", percHS)
    }
  }
  pltDf$Distribution <- factor(pltDf$Distribution,  levels=c("Normal","NotNormal", "LightSkew", "HeavySkew"))
  pltDf$Histotype <- factor(pltDf$Histotype, levels=c("CCC", "EC", "HGSC", "MC"))
  pltDf$Percent <- as.numeric(pltDf$Percent)
  
  bDf <- pltDf[which(pltDf$Value %in% "Beta"), ]
  mDf <- pltDf[which(pltDf$Value %in% "M"), ]
  
  bPlt <- ggplot(bDf, aes(fill=SID, y=Percent, x=Histotype)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    xlab("Histotype") +
    ylab("Percent") + 
    ggtitle(paste("Percentage of CpG sites in histotype (n=685650)" , sep="")) + 
    scale_fill_viridis(discrete = TRUE) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Distribution, scales = "free_x")
  
  mPlt <- ggplot(mDf, aes(fill=SID, y=Percent, x=Histotype)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    xlab("Histotype") +
    ylab("Percent") + 
    ggtitle(paste("Percentage of CpG sites in histotype (n=685650)" , sep="")) + 
    scale_fill_viridis(discrete = TRUE) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Distribution, scales = "free_x")
  
  mergePlt <- bPlt / mPlt
  
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  outFile <- paste(outDir, "Distribution_Types_Consensus_Barplot_Lst.pdf", sep="")
  ggsave(outFile, plot=mergePlt, width=30, height=40, units = "cm")
}

makeHSPModelPlot <- function(inpBeta, inpPheno, inpCpgAnnDf, fileExt = NULL, dirExt = NULL){
  # Plot any given genomic region given in a model input format
  inpCpgAnnDf <- inpCpgAnnDf[which(inpCpgAnnDf$rowID %in% rownames(inpBeta)), ]
  cpgBeta <- inpBeta[inpCpgAnnDf$rowID, ]
  if(is.null(dirExt)){
    outDir <- paste(plotPath)
  }else{
    outDir <- paste(plotPath, dirExt, sep="/")
  }
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(!is.null(fileExt)){
    outFile <- paste(outDir,"/HSP_CpG_Model_", fileExt, ".pdf", sep="")
  }else{
    outFile <- paste(outDir,"/HSP_CpG_Model.pdf", sep="")
  }
  # Create matching column for samples in input-annotation
  inpPheno$matchCol <- paste(inpPheno$Histotype, inpPheno$barcode , sep="_")
  pltDf <- data.frame(t(cpgBeta))
  pltDf$Histotype <- inpPheno$Histotype[match(row.names(pltDf), inpPheno$barcode)]
  mergeDf <- tidyr::gather(pltDf, key="CpG", value="Beta",  colnames(pltDf)[1]:colnames(pltDf)[ncol(pltDf)-1], factor_key=TRUE)
  mergeDf$Trimean <- NA
  for(k in 1:length(unique(mergeDf$Histotype))){
    medH <- names(table(mergeDf$Histotype))[k]
    tmpRows <- mergeDf[which(mergeDf$Histotype %in% medH), ]  
    betaTM <- as.data.frame(tmpRows[,c("CpG", "Beta")] %>% 
                              group_by(CpG) %>%  
                              summarise_all(.funs = litteR::trimean))
    for(l in 1:nrow(betaTM)){
      mergeDf$Trimean[which(mergeDf$CpG %in% betaTM$CpG[l] & mergeDf$Histotype %in% medH)] <- betaTM$Beta[l]
    }
  }
  mergeDf$Histotype_CpG <- inpCpgAnnDf$Histotype[match(mergeDf$CpG, inpCpgAnnDf$rowID)]
  mergeDf$combCol <- paste(mergeDf$Histotype, mergeDf$CpG, sep="_")
  mergeDf <-  mergeDf[!duplicated(mergeDf$combCol),]
  # Change the point size, and shape
  tmpPlt <- ggplot(mergeDf, aes(x=CpG, y=as.numeric(Trimean), col=Histotype)) + 
    geom_point(size=4) +
    labs(x="CpG's in HSP model", 
         y="Beta (trimean)",
         title="CpG Trimean in HSP model") +
    scale_color_manual(values = catColProf$Histotype) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = .15)) +
    scale_x_discrete() +
    geom_hline(yintercept=0.2) +
    geom_hline(yintercept=0.7) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) + 
    facet_wrap(~factor(Histotype_CpG), nrow =1, strip.position = "bottom", scales = "free")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  ggsave(outFile, plot=tmpPlt, width=60, height=40, units = "cm")
}

makeHSPModelPlotLstBased <- function(inpBetaLst, inpPhenoLst, inpCpgAnnDf, fileExt = NULL, dirExt = NULL){
  if(is.null(dirExt)){
    outDir <- paste(plotPath)
  }else{
    outDir <- paste(plotPath, dirExt, sep="/")
  }
  inpCpgAnnDf <- inpCpgAnnDf[!duplicated(inpCpgAnnDf$rowID),]
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(!is.null(fileExt)){
    outFile <- paste(outDir,"/HSP_CpG_Model_Lst_", fileExt, ".pdf", sep="")
  }else{
    outFile <- paste(outDir,"/HSP_CpG_Model_Lst.pdf", sep="")
  }
  
  mergeDfLst <- list()
  for(i in 1:length(inpBetaLst)){
    inpBeta <- inpBetaLst[[i]]
    tmpSid <- names(inpBetaLst)[i]
    inpPheno <- inpPhenoLst[[tmpSid]]
    
    inpCpgAnnDf <- inpCpgAnnDf[which(inpCpgAnnDf$rowID %in% rownames(inpBeta)), ]
    cpgBeta <- inpBeta[inpCpgAnnDf$rowID, ]
    
    # Create matching column for samples in input-annotation
    inpPheno$matchCol <- paste(inpPheno$Histotype, inpPheno$barcode , sep="_")
    pltDf <- data.frame(t(cpgBeta))
    pltDf$Histotype <- inpPheno$Histotype[match(row.names(pltDf), inpPheno$barcode)]
    mergeDf <- tidyr::gather(pltDf, key="CpG", value="Beta",  colnames(pltDf)[1]:colnames(pltDf)[ncol(pltDf)-1], factor_key=TRUE)
    mergeDf$Trimean <- NA
    mergeDf$Sid <- tmpSid
    for(k in 1:length(unique(mergeDf$Histotype))){
      medH <- names(table(mergeDf$Histotype))[k]
      tmpRows <- mergeDf[which(mergeDf$Histotype %in% medH), ]  
      betaTM <- as.data.frame(tmpRows[,c("CpG", "Beta")] %>% 
                                group_by(CpG) %>%  
                                summarise_all(.funs = litteR::trimean))
      for(l in 1:nrow(betaTM)){
        mergeDf$Trimean[which(mergeDf$CpG %in% betaTM$CpG[l] & mergeDf$Histotype %in% medH)] <- betaTM$Beta[l]
      }
    }
    mergeDf$Histotype_CpG <- inpCpgAnnDf$Histotype[match(mergeDf$CpG, inpCpgAnnDf$rowID)]
    mergeDf$combCol <- paste(mergeDf$Histotype, mergeDf$CpG, sep="_")
    mergeDf <-  mergeDf[!duplicated(mergeDf$combCol),]
    mergeDfLst[[i]] <- mergeDf
    names(mergeDfLst)[i] <- tmpSid
  }
  allMerge = Reduce(function(...) merge(..., all=T), mergeDfLst)
  
  allMerge <- allMerge[which(allMerge$Histotype %in% names(table(inpPhenoLst[[1]]$Histotype))), ]
  allMerge$Sid <- gsub("_[0-9]{1}", "", allMerge$Sid)
  allMerge$Sid <- factor(allMerge$Sid)
  # Change the point size, and shape
  tmpPlt <- ggplot(allMerge, aes(x=CpG, y=as.numeric(Trimean), col=Histotype, shape=Sid)) + 
    geom_point(size = 4, stroke = 2) +
    labs(x="CpG's in HSP model", 
         y="Beta (trimean)",
         title="HSP model CpG site Trimean values") +
    scale_color_manual(values = catColProf$Histotype) +
    scale_shape_manual(values=seq(0,length(inpBetaLst))) + 
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = .15)) +
    scale_x_discrete() +
    geom_hline(yintercept=0.2) +
    geom_hline(yintercept=0.7) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) + 
    facet_wrap(~factor(Histotype_CpG), nrow =1, strip.position = "bottom", scales = "free")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  ggsave(outFile, plot=tmpPlt, width=60, height=40, units = "cm")
}

makeDevBoxPlt <- function(inpRanked, fileExt=NULL){
  pltDf <- data.frame(matrix(nrow=0, ncol=3))
  colnames(pltDf) <- c("Histotype", "Max", "Mean")
  for(i in 1:length(inpRanked)){
    tmpRank <-inpRanked[[i]]
    tmpH <- names(inpRanked)[i]
    hStr <- paste(tmpH," (n=", nrow(tmpRank), " HSPs)", sep="")
    for(j in 1:nrow(tmpRank)){
      tmpR <- tmpRank[j,]
      pltDf[nrow(pltDf)+1,] <- c(hStr, tmpR$maxStdGene, tmpR$stdMean)
    }
  }
  pltDfLong <-  pltDf %>% 
    pivot_longer(
      cols = c("Max", "Mean"), 
      names_to = "SD_Type",
      values_to = "Value"
    )
  pltDfLong <- data.frame(pltDfLong)
  pltDfLong$Value <- as.numeric(pltDfLong$Value)
  tmpPlt <- ggplot(pltDfLong, aes(x=SD_Type, y=Value, fill=SD_Type)) + 
    geom_boxplot() +
    labs(title=paste("SD (Beta) for HSP CpG's in external cohorts\nMax=Highest SD in promoter\nMean=Mean SD of promoter", sep=""),
         x ="SD-type", 
         y = "SD (Beta)") +
    scale_fill_manual(values = c("Max" = viridis(n=4, option="B")[2], "Mean"=viridis(n=4, option="B")[3])) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Histotype, scales = "free_x") + 
    guides(fill=guide_legend(title="SD-type"))
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(!is.null(fileExt)){
    outFile <- paste(outDir,fileExt, "_Deviance_Ranked_Plt.pdf", sep="")
  }else{
    outFile <- paste(outDir, "Deviance_Ranked_Plt.pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
}

makeTTestDunnCompPlt <- function(inpBeta, inpPheno, fileAdd = NULL, mBool = NULL){
  # Calculate TM differences across cohort, also include students t-test
  if(!is.null(mBool)){
    inpM <- log2(inpBeta/(1-inpBeta))
  }else{
    inpM <- inpBeta
  }
  histLst <- list()
  for(j in 1:length(names(table(inpPheno$Histotype)))){
    tmpH <- names(table(inpPheno$Histotype))[j]
    pValSumDf <- data.frame(matrix(nrow=0, ncol=8))
    colnames(pValSumDf) <- c("Histotype", "Dun", "TTest", "WC", "CW", "Method", "p", "CpG")
    for(i in 1:nrow(inpM)){
      if(abs(i)%%round(nrow(inpBeta)/10) == 0){
        message(paste("Processing row: [", i, "/", nrow(inpM),"]"))
      }
      tmpCpg <- rownames(inpM)[i]
      tmpRow <- data.frame(inpM[i,], check.names = FALSE)
      hM <- inpM[i,inpPheno$barcode[which(inpPheno$Histotype %in% tmpH)]]
      nHm <- inpM[i,inpPheno$barcode[which(!inpPheno$Histotype %in% tmpH)]]
      # Perform pairwise comparisons
      nHists <- names(table(inpPheno$Histotype))[which(!names(table(inpPheno$Histotype)) %in% tmpH)]
      
      dunn2Pheno <- inpPheno
      dunn2Pheno$Histotype <- factor(dunn2Pheno$Histotype, levels=c(tmpH, nHists))
      
      # Perform Kruskal Wallis test
      cwTest <- kruskal.test(x=as.numeric(tmpRow) , g=dunn2Pheno$Histotype)
      # Perform dunn's test
      dunn2 <- DunnTest(x=as.numeric(tmpRow) , g=dunn2Pheno$Histotype, method = "BH") 
      dunn2 <- data.frame(dunn2[[1]])
      dunnRows <- sapply(rownames(dunn2), function(x) strsplit(x, "-")[[1]][[2]], USE.NAMES=FALSE)
      keepDunn <- dunn2[which(dunnRows %in% tmpH),]
      
      pValDf <- data.frame(matrix(nrow=length(nHists), ncol=5))
      colnames(pValDf) <- c("Dun", "TTest", "WC", "TTDiff", "WCDiff")
      rownames(pValDf) <- nHists
      for(k in 1:length(nHists)){
        tmpHist2 <- nHists[k]
        # Dunn
        # In statistics, Dunnett's test is a multiple comparison procedure to compare each of a number of treatments with a single control
        # We use M-values here to avoid the issue of heteroscedasticity found in beta-values
        nHistM <- nHm[, inpPheno$barcode[which(inpPheno$Histotype %in% tmpHist2)]]
        dunnRow <- keepDunn[grep(tmpHist2, rownames(keepDunn)),]
        dunnP <- dunnRow$pval
        # Students t-test
        tmpTtest <- t.test(as.numeric(hM), as.numeric(nHistM), alternative = "two.sided")
        ttP <- tmpTtest$p.value 
        # Wilcoxon rank sum test
        tmpWc <- wilcox.test(as.numeric(hM), as.numeric(nHistM), alternative = "two.sided")
        wcP <- tmpWc$p.value 
        pValDf[tmpHist2,] <- c(dunnP, ttP, wcP, abs(ttP-dunnP), abs(wcP-dunnP))
      }
      pValDf$CW <- cwTest$p.value
      pValDf$CWDiff <- abs(pValDf$CW - pValDf$Dun)
      pValDf <- tibble::rownames_to_column(pValDf, "Histotype")
      diffCols <- grep("Diff", colnames(pValDf))
      pValLong <- data.frame(pValDf %>% 
                               pivot_longer(
                                 cols = diffCols, 
                                 names_to = "Method",
                                 values_to = "p"))
      pValLong$CpG <- tmpCpg
      pValSumDf <- rbind(pValSumDf, pValLong)
    }
    histLst[[j]] <- pValSumDf
    names(histLst)[j] <- tmpH
  }
  
  # Merge lists, use reference histotype to group plots
  allMerge <- data.frame(rbindlist(histLst, idcol = TRUE ))
  colnames(allMerge)[which(colnames(allMerge) %in% ".id")] <- "Reference"
  allMerge$p <- -log10(allMerge$p)
  
  tmpPlt <- ggplot(allMerge, aes(fill=Histotype, y=p, x=Method)) + 
    geom_boxplot() +
    #geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    xlab("Method") +
    ylab("-log10(p)") + 
    ggtitle(paste("p-value difference vs. Dunn's test" , sep="")) + 
    scale_fill_viridis(discrete = TRUE) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    facet_grid(rows = vars(Reference), 
               scales="free_y", switch = 'y')
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(is.null(fileAdd)){
    outFile <- paste(outDir, "p_val_var.pdf", sep="/")
  }else{
    outFile <- paste(outDir, "/p_val_var_", fileAdd, ".pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt, width=60, height=40, units = "cm")
  
  allMerge_2 <-  data.frame(allMerge %>% 
                              pivot_longer(
                                cols = c(3:5), 
                                names_to = "pTest",
                                values_to = "p2"))
  allMerge_2$p2 <- -log10(allMerge_2$p2)
  tmpPlt2 <- ggplot(allMerge_2, aes(fill=Histotype, y=p2, x=pTest)) + 
    geom_boxplot() +
    #geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    xlab("Method") +
    ylab("-log10(p)") + 
    ggtitle(paste("p-value difference vs. Dunn's test" , sep="")) + 
    scale_fill_viridis(discrete = TRUE) +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    facet_grid(rows = vars(Reference), 
               scales="free_y", switch = 'y')
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  if(is.null(fileAdd)){
    outFile <- paste(outDir, "p_val_var_all.pdf", sep="/")
  }else{
    outFile <- paste(outDir, "/p_val_var_all_", fileAdd, ".pdf", sep="")
  }
  ggsave(outFile, plot=tmpPlt2, width=60, height=40, units = "cm")
}

makePromoCovBarPlt <- function(inpCov){
  # Create dataframe showing the number of promoter regions with a set number of CpG sites represented in array
  pltDf <- data.frame(matrix(nrow=0, ncol=3))
  rownames(inpCov)[1] <- paste(rownames(inpCov)[1], ", n=865859 CpGs", sep="")
  rownames(inpCov)[2] <- paste(rownames(inpCov)[2], ", n=485512 CpGs", sep="")
  colnames(pltDf) <- c("Platform", "Coverage", "Percent")
  for(j in 1:nrow(inpCov)){
    perc1 <- round(inpCov[j,"1-3"]/inpCov[j,"ALL"],2)
    perc2 <- round(inpCov[j,"4-6"]/inpCov[j,"ALL"],2)
    perc3 <- round(inpCov[j,"7-10"]/inpCov[j,"ALL"],2)
    perc4 <- round(inpCov[j,"10+"]/inpCov[j,"ALL"],2)
    pltDf[nrow(pltDf)+1,] <- c(rownames(inpCov)[j], "1-3", perc1)
    pltDf[nrow(pltDf)+1,] <- c(rownames(inpCov)[j], "4-6", perc2)
    pltDf[nrow(pltDf)+1,] <- c(rownames(inpCov)[j], "7-10", perc3)
    pltDf[nrow(pltDf)+1,] <- c(rownames(inpCov)[j], "10+", perc4)
  }
  pltDf$Platform<- factor(pltDf$Platform,  levels=c("EPIC, n=865859 CpGs", "450K, n=485512 CpGs"))
  pltDf$Coverage <- factor(pltDf$Coverage, levels=c("1-3", "4-6", "7-10", "10+"))
  pltDf$Percent <- 100*as.numeric(pltDf$Percent)
  
  tmpPlt <- ggplot(pltDf, aes(fill=Coverage, y=Percent, x=Coverage)) + 
    geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
    xlab("Histotype") +
    ylab("Percent") + 
    ggtitle(paste("Percentage genes with X CpG sites\n in promoter region (n=24979 genes)" , sep="")) + 
    scale_fill_viridis(discrete = TRUE, option = "F") +
    theme(text = element_text(size=24), 
          axis.text.x = element_text(size=16, face="bold"),
          legend.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Platform, scales = "free_x") + 
    guides(fill=guide_legend(title="No. CpG's \nin promoter"))
  
  outDir <- paste(plotPath, "/", sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  outFile <- paste(outDir, "CpG_Promoter_Coverage_Barplot.pdf", sep="")
  ggsave(outFile, plot=tmpPlt, width=30, height=20, units = "cm")
  return(pltDf)
}


#################################################################
# Deprecated
#################################################################

# makeMethMDS <- function(inpBeta, inpPheno, focusCol = "Histotype", nameBool=NULL){
#   # Function for creating PCA plots  from methylation data (M-values)
#   inpPheno <- inpPheno[inpPheno$barcode %in% colnames(inpBeta),]
#   inpBeta <- inpBeta[,colnames(inpBeta) %in% inpPheno$barcode]
#   bestBeta <- as.matrix(inpBeta)
#   # Match phenotypic data to expression data
#   inpPheno <- inpPheno[match(colnames(bestBeta),inpPheno$barcode),]
#   inpPheno <- inpPheno[!is.na(inpPheno$barcode), ]
#   if(!is.null(nameBool)){
#     outFile <- paste(plotPath, "/", "PCA_Plot_Meth_",focusCol,"_", nameBool, ".png", sep="")
#     tmpTitle <- paste("Methylation PCA-plot: ", nameBool, sep="")
#   }else{
#     outFile <- paste(plotPath, "/", "PCA_Plot_Meth_",focusCol,".png", sep="")
#     tmpTitle <- paste("Methylation PCA-plot: ", focusCol, sep="")
#   }
#   htmltools::capturePlot(
#     minfi::mdsPlot(dat = bestBeta, numPositions = 1000, #sampNames = inpPheno$barcode,
#                    sampGroups = inpPheno$Histotype, pch = 16,
#                    pal = catColProf$Histotype[order(catColProf$Histotype)], legendPos = "bottomleft"),
#     filename = outFile,
#     device =defaultPngDevice(),
#     width = 1200,
#     height = 1200,
#     res = 144
#   )
# }
