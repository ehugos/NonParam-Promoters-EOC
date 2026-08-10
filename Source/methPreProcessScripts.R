################################################################################
################################################################################
################################################################################
# Scripts involved in the preprocessing and QC of beta matrices from .idat files
################################################################################
################################################################################
################################################################################

makeProbeConversionRate <- function(rgInp, phenoInp, inpDir, focusCol, pltBool = NULL){
  # Function for generating probe conversion rate and plot it for rgSet
  focusLoc <- which(colnames(phenoInp)==focusCol)
  # Percentage probe conversion rate 
  bsConRes <- wateRmelon::bscon(rgInp)
  bsMat <- cbind(bsConRes, phenoInp[,focusLoc])
  colnames(bsMat) <- c("PercMeth", focusCol)
  bsDf <- as.data.frame(bsMat)
  bsDf <- tibble::rownames_to_column(bsDf, "SampleName")
  bsDf$PercMeth <- as.numeric(bsDf$PercMeth)
  bsDf$PercMeth <- round(bsDf$PercMeth, digits = 0)
  # Plot probe conversion rate
  plt <- ggplot(bsDf, 
                aes(x = SampleName ,
                    y = PercMeth,
                    fill = get(focusCol))) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = catColVec) +
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank()) + 
    xlab(paste0("Sample")) + 
    ylab(paste0("Percent CpG Methylated")) 
  if(!is.null(pltBool)){
    ggsave(paste(inpDir, focusCol, "_Probe_Conversion_plot.pdf",sep=""), plot=plt, width=60, height=40, units = "cm")
  }
  # Return probe conversion rate dataframe
  return(bsDf)
}

makeMethBetaFreqPlt <- function(inpBeta, inpPheno, inpSid, inpDir){
  # Function for creating frequency-plot of beta-values for groups in cohort
  inpBeta <- inpBeta[,colnames(inpBeta) %in% inpPheno$barcode]
  inpPheno <- inpPheno[inpPheno$barcode %in% colnames(inpBeta), ]
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
    ylab("Frequency (%)") + 
    theme(text = element_text(size=18), 
          legend.text=element_text(size=12),
          axis.text.x = element_text(size=16),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste(inpDir, inpSid, "_Beta_Density.pdf",sep=""), plot= tmpPlt, width=30, height=20, units = "cm")
}

makeBetaMatrix <- function(inpPheno, inpArrType, inpName, inpDir, writeBool = NULL, inpCores = NULL){
  if(is.null(inpCores)){
    inpCores <- 6
  }
  # Read in the raw data from the IDAT files; warnings can be ignored.
  rgSet <- read.metharray.exp(targets=inpPheno, verbose = TRUE, force=TRUE, extended = TRUE)
  # Get manifest for rgSet using package function
  manfst <- getManifest(rgSet)
  # Get probe-information associated with the manifest-object
  probeDat <- getProbeInfo(manfst)
  message("Data import completed!")
  # ---------------------------------------
  # Outlier detection & probe filtering
  # ---------------------------------------
  message("Preprocessing initiated")
  # Filter out any failed probes from the rgSet
  detP <- detectionP(rgSet)
  # Remove probes with failed conversion rates
  remove <- apply(detP, 1, function (x) any(x > 0.01))
  rgSet <- rgSet[!rownames(rgSet) %in% names(which(remove)),]
  # Remove samples that fail detection p thresholds (0.01) for more then 5% of samples
  keep_samples <- names(colMeans(detP < 0.01) > 0.95)
  rgSet <- rgSet[,which(colnames(rgSet) %in% keep_samples)]
  # Remove cross-reactive probes
  xloci <- maxprobes::xreactive_probes(array_type = inpArrType)
  rgSet <- rgSet[!rownames(rgSet) %in% xloci,]
  # Identify samples with a low probe conversion rate, remove these from DF
  bsDf <- makeProbeConversionRate(rgInp = rgSet, phenoInp = inpPheno, inpDir = inpDir, focusCol="Histotype")
  # Remove samples with probe conversion rate < 80%
  lowMetSamples <- bsDf[which(bsDf$PercMeth < 80),]
  # If any are found, we remove these
  if(nrow(lowMetSamples)>=1){
    message(paste(nrow(lowMetSamples), " Samples were found to have a probe conversion rate below 80% and were removed from analysis", sep=""))
    rgSet <- rgSet[,!colnames(rgSet) %in% lowMetSamples$SampleName]
    detP <- detP[,!colnames(detP) %in% lowMetSamples$SampleName]
  }
  # Get beta-values from rgSet
  rawBeta <- minfi::getBeta(rgSet)
  # Transform red/green values into a processed methyl set for beta-values
  MSet <- preprocessRaw(rgSet)
  rawPheno <- pData(MSet)
  
  message("Preprocessing completed!")
  # ---------------------------------------
  # Additional Filtering steps (ChAMP)
  # ---------------------------------------
  message("ChAMP probe filtering initiated!")
  champP <- detP
  champBeta <- rawBeta
  colnames(champP) <- rawPheno$Sample_ID
  colnames(champBeta) <- rawPheno$Sample_ID
  champFilter <- champ.filter(beta = champBeta,
                              pd=rawPheno,
                              detP=champP,
                              autoimpute=TRUE,
                              filterDetP=TRUE,
                              SampleCutoff=0.1,
                              detPcut=0.05,
                              filterBeads=TRUE,
                              beadCutoff=0.05,
                              filterXY = FALSE,
                              filterNoCG = TRUE,
                              filterSNPs = TRUE,
                              filterMultiHit = TRUE,
                              fixOutlier = TRUE,
                              arraytype = inpArrType)
  message("ChAMP probe filtering completed!")
  #------------------------
  # Noob + BMIQ
  #------------------------
  message("Noob background correction")
  # Noob background correction
  noobSet <- minfi::preprocessNoob(rgSet,
                                   dyeCorr = TRUE,
                                   verbose = TRUE)
  noobSet <- noobSet[rownames(noobSet) %in% rownames(champFilter$beta),]
  message("Noob completed")
  # Probe liftover? 
  message("BMIQ normalization")
  noobAnno <- minfi::getAnnotation(noobSet)
  #noobMean <- apply(noobSet,1,mean)
  #filter_gene <- rownames(subset(noobMean ,noobMean ==0))
  betaNoob <- getBeta(noobSet)
  nBmiqSet <- champ.norm(beta=betaNoob,
                         resultsDir= paste(inpDir, inpName, sep="/"),
                         method="BMIQ",
                         plotBMIQ=FALSE,
                         arraytype=inpArrType,
                         cores=noCores)
  message("BMIQ completed")
  if(!is.null(writeBool)){
    # Save Beta, M value matrix into separate .csv files for later use in study
    write.csv(nBmiqSet, paste(inpDir, inpName, "/", inpName, "_Beta.csv", sep=""))
  }
  return(nBmiqSet)
}