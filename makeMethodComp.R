################################################################################
################################################################################
################################################################################
# Script for comparing source method to external packages
################################################################################
################################################################################
################################################################################

################################################################################
# Install packages used in scripts
################################################################################
bcPacks <- c("biomaRt", "GenomicRanges","GenomicFeatures", "TxDb.Hsapiens.UCSC.hg38.knownGene", 
             "minfi", "ChAMP", "sesame", "GEOquery", "goseq", "conumee", "EmpiricalBrownsMethod", 
             "ComplexHeatmap", "BiocParallel")
miss_bc_packs <- bcPacks[!(bcPacks %in% installed.packages()[,"Package"])]
if(length(miss_bc_packs)){
  BiocManager::install(miss_bc_packs)
} 
rPacks <- c("devtools","tidymodels", "rtracklayer", "ggplot2", "RColorBrewer", "viridis",
            "ggpubr", "gridExtra", "dplyr", "data.table", "stringr", "tidyr", "tibble",
            "forcats", "cowplot", "colorRamp2", "tibble", "rsample", "yardstick", 
            "DescTools", "patchwork", "htmltools", "future","future.apply", "brunnermunzel")
miss_r_packs <- rPacks[!(rPacks %in% installed.packages()[,"Package"])]
if(length(miss_r_packs)){
  install.packages(miss_r_packs)
} 

################################################################################
# Load packages
################################################################################
invisible(lapply(rPacks, library, character.only = TRUE))
invisible(lapply(bcPacks, library, character.only = TRUE))
# Install hg38 annotation from github
#install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)

################################################################################
################################################################################
# Load global parameters
message("Loading Parameters")
################################################################################
################################################################################
# Disable scientific parameters to show exact values
options(scipen = 999)
histotypes = c("CCC","EC", "HGSC", "MC")
focusGrp <- "Histotype"
alphaV <- 0.05

################################################################################
# Global paths
################################################################################
wd <- getwd()
dataPath <- paste(wd, "Data/", sep="/")
outPath <- paste(wd, "Export/", sep="/")
plotPath <- paste(outPath, "Plots/", sep="")
target_path <- paste(dataPath, "Train/trainPheno.csv", sep="")

################################################################################
# Plot color parameters
################################################################################
colProf <- viridis(50)
catColProf <-  list(Histotype = c("CCC" = viridis(20)[1],
                                  "EC" = viridis(20)[7],
                                  "HGSC" = viridis(20)[13],
                                  "MC" =  viridis(20)[19]),
                    Sur_Grp = c("x<3" = plasma(15)[2],
                                "6>x>3" = plasma(15)[7],
                                "x>6" = plasma(15)[12]))

# If Paths do not exist, create them
ifelse(!dir.exists(file.path(outPath)), dir.create(file.path(outPath)), FALSE)
ifelse(!dir.exists(file.path(plotPath)), dir.create(file.path(plotPath)), FALSE)

# Set up parallel processing (Win - Based on AMD processor)
totCores = length(availableWorkers(logical = FALSE))
# Use 70% of available cores for parallel processing
nCores = round(totCores/2)
# Set the memory allocation in gigabytes for each worker 
options(future.globals.maxSize = 3000 * 1024^2)
plan(multisession, workers = nCores) 

################################################################################
################################################################################
# Source functions
################################################################################
################################################################################

source(paste(wd, "/Source/betaScripts.R", sep=""))
source(paste(wd, "/Source/genScripts.R", sep=""))
source(paste(wd, "/Source/pltScripts.R", sep=""))
source(paste(wd, "/Source/extPkgScripts.R", sep=""))

################################################################################
################################################################################
# Load data
################################################################################
################################################################################

hg38GenePromoterCpgs <- read.csv(paste(outPath, "hg38GenePromoterCpgs.csv", sep=""), row.names = 1)
hg38CpgLocsEPIC <- read.csv(paste(outPath, "hg38CpgLocsEPIC.csv", sep=""), row.names = 1)
hg19CpgLocsEPIC <- read.csv(paste(outPath, "hg19CpgLocsEPIC.csv", sep=""), row.names = 1)
transTableEPIC <- read.csv(paste(outPath, "transTableEPIC.csv", sep=""), row.names = 1)
geneInf <- read.csv(paste(outPath, "hg38GeneInf.csv", sep=""), row.names = 1)
hg19GeneInf <- read.csv(paste(outPath, "hg19GeneInf.csv", sep=""), row.names = 1)

################################################################################
# Read in phenotype data
################################################################################
trainPheno <- read.csv(paste(dataPath, "Train/Trainpheno.csv", sep=""), row.names = 1)

################################################################################
# Read in beta-matrices
################################################################################
trainBeta <- read.csv(paste(dataPath, "Train/Train_Beta.csv",sep=""), row.names = 1)
trainM <- log2(trainBeta/(1-trainBeta))

################################################################################
# Get significant promoter regions
################################################################################
promoRankLst_Files <- list.files(path = outPath, 
                                 pattern = "promoters_ranked.csv",
                                 recursive = TRUE,
                                 full.names = TRUE)
allSigPromo <- makeDfLst(promoRankLst_Files)
allSigPromo <- allSigPromo[!names(allSigPromo) %in% "EC"]

################################################################################
# Get significant CpGs in signficant promoter regions
################################################################################
allSigCpgPromo <- list()
for(i in 1:length(allSigPromo)){
  tmpSigs <- allSigPromo[[i]]
  tmpH <- names(allSigPromo)[i]
  sigPromoVec <- c()
  for(j in 1:nrow(tmpSigs)){
    tmpEns <- tmpSigs$ensembl_gene_id[j]
    tmpReg <- hg38GenePromoterCpgs[hg38GenePromoterCpgs$Gene %in% tmpEns,]
    tmpDf <-  makeCpgDf(data.frame(tmpReg),
                        trainBeta, 
                        trainPheno, 
                        tmpH)
    tmpCpg <- makeRegSigCpg(inpReg = tmpDf,
                            inpBeta = trainBeta,
                            inpPheno = trainPheno,
                            inpH = tmpH,
                            minCpg = 2)
    if(length(tmpCpg) > 0){
      sigPromoVec <- append(sigPromoVec, tmpCpg)
    }
  }
  allSigCpgPromo[[i]] <- sigPromoVec
  names(allSigCpgPromo)[i] <- tmpH
}
makeCsvSave(allSigCpgPromo, "allSigCpgPromo")

allSigCpgPromo_Files <- list.files(path = outPath, 
                                 pattern = "allSigCpgPromo.csv",
                                 recursive = TRUE,
                                 full.names = TRUE)
allSigCpgPromo <- makeDfLst(allSigCpgPromo_Files )

################################################################################
################################################################################
# Run DMP analysis on training cohort
message("Running ChAMP and Minfi DMP analysis on training cohort")
################################################################################
################################################################################

################################################################################
# Run ChAMP DMP analysis (limma)
################################################################################
champDmpLst_M <- list()
for(i in 1:length(names(table(trainPheno$Histotype)))){
  tmpH <- names(table(trainPheno$Histotype))[i]
  message(paste("Now runing DMP analysis for histotype:", tmpH))
  tmpDmp <- makeChampDmpRes(trainM, 
                            trainPheno, 
                            focH = tmpH,
                            bCut = 0)
  for(j in 1:length(tmpDmp)){
    champDmpLst_M[[length(champDmpLst_M)+1]] <- tmpDmp[[j]]
    names(champDmpLst_M)[length(champDmpLst_M)] <- names(tmpDmp)[j]
  }
}

makeCsvSave(resLst = champDmpLst_M, 
            sType = "champDMP_M", 
            dirExtBool = "DMP")

################################################################################
# Run minfi DMP analysis
################################################################################
minfiDmpLst_M  <- list()
for(i in 1:length(names(table(trainPheno$Histotype)))){
  tmpH <- names(table(trainPheno$Histotype))[i]
  message(paste("Now runing DMP analysis for histotype:", tmpH))
  tmpDmp <- makeMinfiDmpRes(inpBeta = trainM, 
                            inpPheno = trainPheno, 
                            focH = tmpH,
                            pCut = 0.05)
  for(j in 1:length(tmpDmp)){
    minfiDmpLst_M[[length(minfiDmpLst_M)+1]] <- tmpDmp[[j]]
    names(minfiDmpLst_M)[length(minfiDmpLst_M)] <- names(tmpDmp)[j]
  }
}

makeCsvSave(resLst = minfiDmpLst_M, 
            sType = "minfiDMP_M", 
            dirExtBool = "DMP")

################################################################################
################################################################################
# Run DMP analysis using non-parametric methods 
# As we do not consider covariates in this study, we use the mann-whitney u-test to test for significance
# And tukeys trimean to estimate effect size

# For complex models, we should consider stratified multivariate mann-whitney (sanon)
# Alternatively the van-elteren test
# Currently very slow, needs to be optimized
message("Running non-parametric DMP analysis")
################################################################################
################################################################################

nParamDmpLst <- list()
for(i in 1:length(names(table(trainPheno$Histotype)))){
  tmpH <- names(table(trainPheno$Histotype))[i]
  nHs <- names(table(trainPheno$Histotype))[!names(table(trainPheno$Histotype)) %in% tmpH]
  for(j in 1:length(nHs)){
    tmpC <- nHs[j]
    # inpBeta, inpPheno, sigBool=NULL, inpM = NULL, pCut = NULL, bCut = NULL
    hspDMPTest <- makeNonParamDMP_MULT_V2(inpBeta = trainBeta, 
                                       inpPheno = trainPheno, 
                                       inpCont = c(tmpH, tmpC),
                                       sigBool = TRUE,
                                       bCut = NULL,
                                       pCut = 0.05)
    nParamDmpLst[[length(nParamDmpLst)+1]] <- hspDMPTest
    names(nParamDmpLst)[length(nParamDmpLst)] <- paste(tmpH, tmpC, sep="_")
  }
}
nParamDmpLst <- nParamDmpLst[order(names(nParamDmpLst))]

makeCsvSave(nParamDmpLst,
            "nonParamDmp",
            dirExtBool = "DMP")

################################################################################
################################################################################
# Compare CpGs seen as significant in promoter regions by source
# To that of the generated DMP results to characterize overlaps
message("Running HSP-DMP comparisons analysis")
################################################################################
################################################################################

champDmpFiles <- list.files(path = outPath, 
                        pattern = "champDMP_M",
                        recursive = TRUE,
                        full.names = TRUE)
champDmpLst <- makeDfLst(champDmpFiles )

minfiDmpFiles <- list.files(path = outPath, 
                             pattern = "minfiDMP_M",
                             recursive = TRUE,
                             full.names = TRUE)
minfiDmpLst <- makeDfLst(minfiDmpFiles)

nParamDmpFiles <- list.files(path = outPath, 
                            pattern = "nonParamDmp",
                            recursive = TRUE,
                            full.names = TRUE)
nParamDmpLst <- makeDfLst(nParamDmpFiles)

################################################################################
# Plot number of DMPs identified in each contrast
################################################################################
dmpAllLst <- list("ChAMP" = champDmpLst,
                  "minfi" = minfiDmpLst)
makeDmpStatPlt(dmpAllLst)

################################################################################
# Get the overlaps and differences between each program and each comparison
################################################################################
dmpDiff <- makeDmpDiffDf(dmpAllLst, pltBool = TRUE)
dmpDiffDf <- dmpDiff[[1]]
dmpMissDf <- dmpDiff[[2]]

write.csv(dmpDiffDf, paste(outPath, "DMP_Overlaps_DMP_Callers.csv",sep=""))

################################################################################
# Investigate the missing DMPs in each comparison
################################################################################

################################################################################
# Rerun without significance cutoff
################################################################################
champDmpAll <- list()
for(i in 1:length(names(table(trainPheno$Histotype)))){
  tmpH <- names(table(trainPheno$Histotype))[i]
  message(paste("Now runing DMP analysis for histotype:", tmpH))
  tmpDmp <- makeChampDmpRes(trainM, 
                            trainPheno, 
                            focH = tmpH,
                            pCut = 1,
                            bCut = 0)
  for(j in 1:length(tmpDmp)){
    champDmpAll[[length(champDmpAll)+1]] <- tmpDmp[[j]]
    names(champDmpAll)[length(champDmpAll)] <- names(tmpDmp)[j]
  }
}

minfiDmpAll <- list()
for(i in 1:length(names(table(trainPheno$Histotype)))){
  tmpH <- names(table(trainPheno$Histotype))[i]
  message(paste("Now runing DMP analysis for histotype:", tmpH))
  tmpDmp <- makeMinfiDmpRes(inpBeta =  trainM,
                            inpPheno =  trainPheno,
                            focH = tmpH,
                            pCut = 1)
  for(j in 1:length(tmpDmp)){
    minfiDmpAll[[length(minfiDmpAll)+1]] <- tmpDmp[[j]]
    names(minfiDmpAll)[length(minfiDmpAll)] <- names(tmpDmp)[j]
  }
}

missDmpChampLst <- list()
missDmpMinfiLst <- list()
for(i in 1:length(nParamDmpLst)){
  tmpnParam <- nParamDmpLst[[i]]
  tmpCont <- names(nParamDmpLst)[i]
  inpCont <- strsplit(tmpCont, "_")[[1]]
  inpCont <- c(inpCont[1], inpCont[2])
  tmpChamp <- champDmpLst[[tmpCont]]
  tmpMinfi <- minfiDmpLst[[tmpCont]]
  missChamp <- rownames(tmpnParam)[which(!rownames(tmpnParam) %in% rownames(tmpChamp))]
  missMinfi <- rownames(tmpnParam)[which(!rownames(tmpnParam) %in% rownames(tmpMinfi))]
  
  missDmpChampLst[[i]] <- makeMissDmpType(inpBeta = trainBeta,
                                  inpPheno = trainPheno, 
                                  inpCont = inpCont,
                                  missCpg = missChamp)
  names(missDmpChampLst)[i] <- tmpCont
  
  missDmpMinfiLst[[i]] <- makeMissDmpType(inpBeta = trainBeta,
                                  inpPheno = trainPheno, 
                                  inpCont = inpCont,
                                  missCpg = missMinfi)
  names(missDmpMinfiLst)[i] <- tmpCont
}

missDmpChampDf <- do.call(rbind.data.frame, missDmpChampLst)
missDmpChampDf <- rownames_to_column(missDmpChampDf, "Cont")
missDmpMinfiDf <- do.call(rbind.data.frame, missDmpMinfiLst)
missDmpMinfiDf <- rownames_to_column(missDmpMinfiDf, "Cont")

missDmpLst <- list("ChAMP" = missDmpChampDf,
                   "minfi" = missDmpMinfiDf)
allMissDf <- do.call(rbind.data.frame, missDmpLst)
allMissDf <- rownames_to_column(allMissDf, "Program")
allMissDf$Program <- gsub("\\.[0-9]{1,2}", "", allMissDf$Program)
makeMissTypeDfPlt(allMissDf)

write.csv(missDmpChampDf, 
          paste(outPath, "DMP/missDmpChampDf.csv",sep=""))
write.csv(missDmpMinfiDf, 
          paste(outPath, "DMP/missDmpMinfiDf.csv",sep=""))

################################################################################
# Generate dataframes with the missing DMPs alongside the champ and limma p-values
################################################################################
compMissDmpLst <- list()
for(i in 1:length(nParamDmpLst)){
  tmpCont <- names(nParamDmpLst)[i]
  tmpNP <- nParamDmpLst[[tmpCont]]
  tmpChamp <- champDmpAll[[tmpCont]]
  tmpMinfi <- minfiDmpAll[[tmpCont]]
  tmpChamp <- tmpChamp[match(rownames(tmpNP), rownames(tmpChamp)),]
  tmpMinfi <- tmpMinfi[match(rownames(tmpNP), rownames(tmpMinfi)),]
  tmpNP$Champ_p <- tmpChamp$P.Value
  tmpNP$Champ_adjP <- tmpChamp$adj.P.Val
  tmpNP$Minfi_adjP <- tmpMinfi$qval
  tmpNP <- tmpNP[which(tmpNP$Champ_adjP > 0.05 | tmpNP$Minfi_adjP > 0.05), ]
  compMissDmpLst[[i]] <- tmpNP
  names(compMissDmpLst)[i] <- tmpCont
}

################################################################################
# See in how many comparisons signifcant CpGs in the promoter regions were found
# As DMPs by ChAMP and minfi
################################################################################
champDmpCov <- makeDmpRegCovDf(allSigCpgPromo, 
                            champDmpLst)
minfiDmpCov <- makeDmpRegCovDf(allSigCpgPromo, 
                            minfiDmpLst)
nParamDmpCov <- makeDmpRegCovDf(inpCpgPos = allSigCpgPromo, 
                             inpDmpLst = nParamDmpLst)

write.csv(champDmpCov, paste(outPath, "DMP/champ_SigPromo_Cov.csv",sep="/"))
write.csv(minfiDmpCov, paste(outPath, "DMP/minfi_SigPromo_Cov.csv",sep="/"))
write.csv(nParamDmpCov, paste(outPath, "DMP/nParam_SigPromo_Cov.csv",sep="/"))

################################################################################
# Plot significant CpGs in significant promoters DMP coverage in ChAMP/minfi
################################################################################
inpCovLst <- list("ChAMP" = champDmpCov, 
                  "minfi" = minfiDmpCov, 
                  "nParam"=nParamDmpCov)
makeDMPRegCovCompPlt(inpCovLst, 
                   fileExt = "sigRegDmpCov")

makeDMPRegCovCompPlt(list("ChAMP" = champDmpCov, 
                          "minfi" = minfiDmpCov), 
                     fileExt = "sigRegDmpCov_NoNP")

################################################################################
# Look at the comparisons distribution type for the significant CpGs that are
# not identified as DMPs in all 3 comparisons
# See missing CpG per comparison, and their associatied distribution type comparison
message("Running SigPromo CpG - DMP missingness analysis")
################################################################################
missSigPromoTypeCpgChampLst <- list()
missSigPromoTypeCpgMinfiLst <- list()
for(i in 1:length(champDmpLst)){
  tmpCont <- names(champDmpLst)[i]
  tmpChampDmp <- champDmpLst[[tmpCont]]
  tmpMinfiDmp <- minfiDmpLst[[tmpCont]]
  inpCont <- strsplit(tmpCont, "_")[[1]]
  inpCont <- c(inpCont[1], inpCont[2])
  sigCpg <- allSigCpgPromo[[inpCont[[1]]]]
  missChampCpg <- sigCpg [which(!sigCpg %in% 
                                 rownames(tmpChampDmp))]
  missMinfiCpg <- sigCpg [which(!sigCpg %in% 
                                 rownames(tmpMinfiDmp))]
  if(length(missChampCpg) > 0){
    message(inpCont)
    tmpChampMiss <- makeMissDmpType(inpBeta = trainBeta,
                                    inpPheno = trainPheno, 
                                    inpCont = inpCont,
                                    missCpg = missChampCpg)
    missSigPromoTypeCpgChampLst[[length(missSigPromoTypeCpgChampLst)+1]] <- tmpChampMiss
    names(missSigPromoTypeCpgChampLst)[length(missSigPromoTypeCpgChampLst)] <- tmpCont
  }
    
  if(length(missMinfiCpg) > 0){
    tmpMinfiMiss <- makeMissDmpType(inpBeta = trainBeta,
                                    inpPheno = trainPheno, 
                                    inpCont = inpCont,
                                    missCpg = missMinfiCpg)
    missSigPromoTypeCpgMinfiLst[[length(missSigPromoTypeCpgMinfiLst)+1]] <- tmpMinfiMiss
    names(missSigPromoTypeCpgMinfiLst)[length(missSigPromoTypeCpgMinfiLst)] <- tmpCont
  }
}

################################################################################
# Create barplot of CpGs in significant promoters not characterized as DMPs 
################################################################################
makeMissDMPPlt(missLst = missSigPromoTypeCpgChampLst,
               refCpgLst = allSigCpgPromo,
               fileAdd = "ChAMP")
makeMissDMPPlt(missLst = missSigPromoTypeCpgMinfiLst, 
               refCpgLst = allSigCpgPromo,
               fileAdd = "minfi")

################################################################################
################################################################################
# Run DMR analysis
message("Running DMR analysis")
################################################################################
################################################################################
conts <- makeCatCombs(trainPheno$Histotype)
conts <- gsub(" - ", "_", conts)

################################################################################
# DMRcate
message("Running DMR analysis using DMRcate")
################################################################################

# Convert DF into Matrix for subsequent use with dmrCate 
betaMat <- as.matrix(trainBeta)

# Generate contrast vectors
hisContVec <- makeCatCombs(trainPheno$Histotype)
# Convert DF into Matrix for subsequent use with dmrCate 
betaMat <- as.matrix(trainBeta)

dmrPheno <- trainPheno
rownames(dmrPheno) <- trainPheno$barcode
betaGrSet <- minfi::makeGenomicRatioSetFromMatrix(
  mat = betaMat,
  pData = dmrPheno,
  array = "IlluminaHumanMethylationEPIC", 
  annotation = "ilm10b5.hg38",
  what="Beta"
)

# Annotate cpg's for contrasts

hisAnnoLst <- makeDmrAnno(inpBeta = betaGrSet, 
                          phenoInp = dmrPheno, 
                          focusCol = "Histotype")
message("Annotate ranges: Pass")

# Create DMR-result objects
hisDmrLst <- makeDmrResult(hisAnnoLst)
message("DMR results: Pass")

# Extract genomic ranges from annotated DMR's
hisRanges <- makeDmrRanges(hisDmrLst)
names(hisRanges) <- gsub(" - ", "_", names(hisRanges))
message("Extract ranges: Pass")

dmrCateDmrLst <- list()
for(i in 1:length(hisRanges)){
  tmpRange <- hisRanges[[i]]
  tmpStart <- tmpRange@ranges@start
  tmpEnd <- tmpRange@ranges@start + tmpRange@ranges@width
  tmpDf <- data.frame(tmpRange)
  tmpDf$start <- tmpStart 
  tmpDf$end <- tmpEnd
  colnames(tmpDf)[which(colnames(tmpDf) %in% "seqnames")] <- "chr"
  dmrCateDmrLst[[i]] <- tmpDf
  names(dmrCateDmrLst)[i] <- names(hisRanges)[i]
}
makeCsvSave(dmrCateDmrLst, 
            "dmrCateDMR", 
            dirExtBool = "DMR")

################################################################################
# Sesame - Warning, long runtime
message("Running DMR analysis using Sesame")
################################################################################

sesameDataCache("EPIC.address")

sesamePheno <- trainPheno 
rownames(sesamePheno) <-  sesamePheno$barcode
# availConts <- colnames(attr(tmpDml, "model.matrix")) 
# Sesames standard setting is hg38 if we choose epic, we must map it to hg19 to enable matching vs other callers
sesCoords <- makeGRangesFromDataFrame(hg19CpgLocsEPIC)

sesameLst <- list()
for(i in 1:length(names(table(sesamePheno$Histotype)))){  
  tmpCont <- names(table(sesamePheno$Histotype))[i]
  # Make contrast vector for DMRcate to use in comparison
  message(paste("Running DML for: ", tmpCont, sep=""))
  # Import diagnosis & histotype as factor-vectors for comparative analysis
  sesamePheno$Histotype <- factor(sesamePheno$Histotype)
  sesamePheno$Histotype <- relevel(sesamePheno$Histotype, 
                                   ref=tmpCont)
  otherConts <- names(table(sesamePheno$Histotype))[!names(table(sesamePheno$Histotype)) %in% tmpCont]
  tmpDML <- sesame::DML(betas = as.matrix(trainBeta),
                        fm = ~Histotype,
                        meta = sesamePheno,
                        BPPARAM = SnowParam(workers = nCores, 
                                            type = "SOCK"))
  for(j in 1:length(otherConts)){
    c2 <- otherConts[j]
    message(paste("Running DML for the contrast: ", c2, "_", tmpCont, sep=""))
    newForm <- paste("Histotype", c2, sep="")
    sesDMR <- sesame::DMR(betas = as.matrix(trainBeta),
                          tmpDML, 
                          contrast = newForm,
                          platform = c("EPIC"),
                          seg.per.locus = 0.5,
                          probe.coords = sesCoords)
    sesDMR <- sesDMR[sesDMR$Seg_Pval_adj < 0.05,]
    sesameLst[[length(sesameLst)+1]] <- sesDMR
    names(sesameLst)[length(sesameLst)] <- paste(c2, tmpCont, sep="_")
  }
}

makeCsvSave(sesameLst,
            dirExtBool = "DMR",
            "sesameDMR")

################################################################################
# Adjust sesame DMRs to be compatible with other packages for comparative purposes 
################################################################################
sesRawFiles <- list.files(path = outPath, 
                          pattern = "sesameDMR",
                          recursive = TRUE,
                          full.names = TRUE)
sesameLst <- makeDfLst(sesRawFiles)

adjSesDmrLst <-  makeSesSegDmr_MULT(inpSesDmrLst = sesameLst, 
                               inpCpgAnn = hg19CpgLocsEPIC,
                               minCpg = 3,
                               filtDist = 0.2)
makeCsvSave(adjSesDmrLst, 
            "AdjSesDMR",
            dirExtBool = "DMR")

################################################################################
# Bumphunter - ChAMP
# Bumphunter is exceedingly buggy, using ChAMP instead
message("Running DMR analysis using ChAMP - Bumphunter")
################################################################################
# Bumphunter works poorly, the champ extension of the module works better
conts <- conts <- makeCatCombs(trainPheno$Histotype)
conts <- gsub(" - ", "_", conts)
bumpHunterLst <- list()
# conts <- conts[!conts %in% c("CCC_EC","CCC_HGSC","CCC_MC","EC_CCC","EC_HGSC","EC_MC")]
for(i in 1:length(conts)){
  tmpCont <- conts[i]
  cs <- strsplit(tmpCont, "_")[[1]]
  contVec <- c(cs[1], cs[2])
  # Set up DMR analysis
  tmpPheno <- trainPheno
  tmpPheno <- tmpPheno[tmpPheno$Histotype %in% contVec,]
  tmpPheno$Histotype <- factor(tmpPheno$Histotype, levels=c(contVec[1], contVec[2]))
  tmpBeta <- trainBeta[, colnames(trainBeta) %in% tmpPheno$barcode]
  message(paste("NOW RUNNING DMR ANALYSIS FOR CONTRAST: ", contVec[1], "_", contVec[2], sep=""))
  tmpBP <- ChAMP::champ.DMR(beta = as.matrix(tmpBeta),
                            pheno = tmpPheno$Histotype,
                            arraytype = "EPIC",
                            method = "Bumphunter",
                            minProbes = 3,
                            adjPvalDmr=0.05,
                            B=100,
                            pickCutoff = 0.95,
                            compare.group=contVec,
                            cores = 10)
  bumpHunterLst[[i]] <- tmpBP
  names(bumpHunterLst)[i] <- tmpCont
}
for(i in 1:length(bumpHunterLst)){
  bumpHunterLst[[i]] <- bumpHunterLst[[i]][[1]] 
}

makeCsvSave(bumpHunterLst, 
            "BumpHunterChamp", 
            dirExtBool = "DMR")

################################################################################
# ProbeLasso
message("Running DMR analysis using probelasso")
################################################################################
pbLassLst <- list()
for(i in 1:length(conts)){
  tmpCont <- conts[i]
  cs <- strsplit(tmpCont, "_")[[1]]
  contVec <- c(cs[1], cs[2])
  tmpPheno <- trainPheno[trainPheno$Histotype %in% contVec,]
  tmpBeta <- trainBeta[, colnames(trainBeta) %in% tmpPheno$barcode]
  message(paste("NOW RUNNING DMR ANALYSIS FOR CONTRAST: ", contVec[1], "_", contVec[2], sep=""))
  pbLass <- ChAMP::champ.DMR(beta = tmpBeta, 
                             pheno = tmpPheno$Histotype,
                             arraytype="EPIC",
                             compare.group=contVec,
                             method = "ProbeLasso",
                             minProbes=3,
                             adjPvalDmr=0.05,
                             cores = 6)
  pbLassLst[[i]] <- pbLass
  names(pbLassLst)[i] <- tmpCont
}

for(i in 1:length(pbLassLst)){
  pbLassLst[[i]] <- pbLassLst[[i]][[1]] 
}
makeCsvSave(pbLassLst, "ProbeLassoDMR")

################################################################################
################################################################################
# DMR overlap analysis
message("investigating DMR overlaps")
################################################################################
################################################################################
dmrCateFiles <- list.files(path = outPath, 
                           pattern = "dmrCateDMR",
                           recursive = TRUE,
                           full.names = TRUE)
dmrCateDmrLst <-  makeDfLst(dmrCateFiles)

adjSesFiles <- list.files(path = outPath, 
                          pattern = "AdjSesDMR",
                          recursive = TRUE,
                          full.names = TRUE)
adjSesDmrLst <- makeDfLst(adjSesFiles)

probeLassoFiles <- list.files(path = outPath, 
                              pattern = "ProbeLassoDMR.csv",
                              recursive = TRUE,
                              full.names = TRUE)
probeLassoDmrLst <- makeDfLst(probeLassoFiles)

bumpHunterFiles <- list.files(path = outPath, 
                         pattern = "BumpHunterChamp.csv",
                         recursive = TRUE,
                         full.names = TRUE)
bumpHunterDmrLst <- makeDfLst(bumpHunterFiles)

extDmrLst <- list("DMRcate"= dmrCateDmrLst, 
                  "bumphunter"= bumpHunterDmrLst,
                  "probelasso"= probeLassoDmrLst,
                  "sesame"= adjSesDmrLst)

########################################################################################
# Plot DMR results next to one another to illustrate the number of hits in each contrast
########################################################################################
makeDmrStatPlt(inpDmrLst = extDmrLst, 
               fileExt = "DmrOverlaps")

#############################################################################################
# Compare DMR caller overlaps (i.e. how many regions in one caller are found found by others)
#############################################################################################
extDMROlPercDf <- makeDMROls(inpDmrLst = extDmrLst)
write.csv(extDMROlPercDf, 
          paste(outPath, "dmrCompOverlapDf.csv", sep="/"))

#############################################################################################
# Plot DMR caller DMR region overlaps
#############################################################################################
colCats <- unique(extDMROlPercDf$Comparison)
colvals <- viridis(n=length(colCats), option="B")
colVec <- colvals
names(colVec) <- colCats
extDMROlPlt <- ggplot(extDMROlPercDf, aes(fill=Comparison, y=percUnq, x=Contrast)) + 
  geom_bar(position = position_dodge2(preserve = "single", padding = 0.05), stat="identity") + 
  xlab("Contrast") +
  ylab("Percent DMR overlaps") + 
  ggtitle("DMR caller overlaps in EOC contrasts\nPercentage of overlapping DMRs (Overlaps/Total DMRs in reference)") + 
  ylim(0, 100) +
  scale_fill_manual(values = colVec) +
  theme(text = element_text(size=24), 
        axis.text.x = element_text(angle = 90, size=16, face="bold"),
        legend.text=element_text(size=16),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Reference)

ggsave(paste(plotPath,"dmrCompOverlapResPlot.pdf", sep="/"),
       plot=extDMROlPlt, width=60, height=40, units = "cm")

################################################################################
################################################################################
# Compare DMR callers with significant promoters, how many HSPs are found as DMRs
message("Comparing HSP coordinates to DMR results")

# Promoter overlaps differ based on if we use hg19 or hg38
# However the probe id and its location is the same irregardless of the mapping
# As such we opt to translate hg38 coordinates to hg19 and then compare
################################################################################
################################################################################

allSigPromoCpgPosHg38 <- list()
for(i in 1:length(allSigPromo)){
  tmpSig <- allSigPromo[[i]]
  tmpH <- names(allSigPromo)[i]
  allPromoCpgPosHg38 <- makeCpgGeneOverlap_MULT(inpGenes = tmpSig$ensembl_gene_id, 
                                            cpgInp = hg38CpgLocsEPIC, 
                                            inpGeneInf = geneInf, 
                                            type="PROMOTER")
  allSigPromoCpgPosHg38[[i]] <- allPromoCpgPosHg38
  names(allSigPromoCpgPosHg38)[i] <- tmpH
}

sigPromoDmrStats <- makeDmrRegionstatDf(inpBeta = trainBeta,
                                        inpDMRLst = extDmrLst,
                                        inpRegions = allSigPromo,
                                        inpRegionsPosLst = allSigPromoCpgPosHg38,
                                        cpgCoords = hg19CpgLocsEPIC)

# Subset so that we only look at regions with 3 or more significant CpGs
# As this is the minimum setting for DMR callers
sigPromoDmrStats <- sigPromoDmrStats[sigPromoDmrStats$sigCpgOld >=3, ]

write.csv(sigPromoDmrStats, 
          paste(outPath, "DMR/sigPromoSigCpgDmrOverlapStats.csv",sep=""))

################################################################################
# Plot coverage of significant regions in X/Y comparisons in DMR caller results
################################################################################
makeDMRRegion_CovCompPlt(inpCov = sigPromoDmrStats, 
                   fileExt = "SigPromo")

################################################################################
# Investigate how many overlaps are found in DMR callers output for significant regions
################################################################################
hitTypeDf <- makeHitTypeDf(inpCov = sigPromoDmrStats)
write.csv(hitTypeDf[[1]], paste(outPath, "sigPromo_NoHitsInExt.csv",sep="/"))
write.csv(hitTypeDf[[2]], paste(outPath, "sigPromo_PhenoStatsInExt.csv",sep="/"))

################################################################################
# Create region dataframe from DMR results, run sigpromo to see how many are classed as sig
################################################################################


