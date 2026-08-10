################################################################################
################################################################################
# Script for comparison against external datasets (long runtimes due to size of datasets)
################################################################################
################################################################################

################################################################################
################################################################################
message("Loading packages")
################################################################################
################################################################################

# Check r-packages, install if missing
r_packs <- c("ggplot2", "RColorBrewer", "viridis", "colorRamp2","ggpubr",
             "gridExtra", "forcats", "cowplot", "tibble", "dplyr", "data.table",
             "stringr", "tidyr", "tidymodels", "e1071", "rsample", "yardstick", "DescTools", "future", "future.apply")
miss_r_packs <- r_packs[!(r_packs %in% installed.packages()[,"Package"])]
if(length(miss_r_packs)){
  install.packages(miss_r_packs)
} 
invisible(lapply(r_packs, library, character.only = TRUE))

# Check bioconductor packages, install if missing
bioc_packs <- c("ChAMP", "GEOquery", "minfi", "minfiData",
                "BiocManager", "GenomicRanges", "GenomicFeatures",
                "minfi", "GEOquery")
miss_bc_packs <- bioc_packs[!(bioc_packs %in% installed.packages()[,"Package"])]
if(length(miss_bc_packs)){
  BiocManager::install(miss_bc_packs)
} 
invisible(lapply(bioc_packs, library, character.only = TRUE))

################################################################################
################################################################################
# Define global parameters
################################################################################
################################################################################

# Disable scientific parameters to get exact values
options(scipen = 999)

histotypes = c("CCC","EC", "HGSC", "MC")
focusGrp <- "Histotype"
alphaV <- 0.05

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

# Get half of the available cores for multiprocessing
# To reset, use "plan(sequential)"
nCores = length(availableWorkers(logical = FALSE))/3
options(future.globals.maxSize = 5000 * 1024^2)
plan(multisession, workers = nCores) 

################################################################################
# Define global paths
################################################################################
wd <- getwd()
dataPath <- paste(wd, "/Data/", sep="/")
outPath <- paste(wd, "/Export/", sep="/")
plotPath <- paste(outPath, "Plots/", sep="")
# If Paths do not exist, create them
ifelse(!dir.exists(file.path(outPath)), dir.create(file.path(outPath)), FALSE)
ifelse(!dir.exists(file.path(plotPath)), dir.create(file.path(plotPath)), FALSE)

################################################################################
################################################################################
# Source functions
################################################################################
################################################################################

source(paste(wd, "/Source/genScripts.R", sep=""))
source(paste(wd, "/Source/betaScripts.R", sep=""))
source(paste(wd, "/Source/pltScripts.R", sep=""))

################################################################################
# Load Cpg-site location mapping data
################################################################################
geneInf <- read.csv(paste(outPath, "hg38GeneInf.csv", sep=""), row.names = 1)
hg38CpgLocsEPIC <- read.csv(paste(outPath, "/hg38CpgLocsEPIC.csv",sep=""), row.names = 1)
hg38CpgLocs450K <- read.csv(paste(outPath, "hg38CpgLocs450K.csv", sep=""), row.names = 1)
promoter_cpgs <- read.csv(paste(outPath, "hg38GenePromoterCpgs.csv", sep=""))
cpgAnno <- read.csv(paste(outPath, "cpgAnno.csv", sep=""), row.names = 1)

################################################################################
# Get gene-CpG overlaps for 450k/EPIC coordinates for significant promoter regions
################################################################################
keepPromoInfLstFiles <- list.files(path = outPath, 
                                   pattern = "sigPromoInfDf",
                                   recursive = TRUE,
                                   full.names = TRUE)
keepPromoInfLst <- makeDfLst(keepPromoInfLstFiles)

cpgOlsHg38450K  <- list()
cpgOlsHg38EPIC <- list()
for(i in 1:length(keepPromoInfLst)){
  tmpInf <- keepPromoInfLst[[i]]
  tmpH <- names(tmpInf)
  ol450K <- makeCpgGeneOverlap_MULT(tmpInf$ensembl_gene_id, 
                                                 hg38CpgLocs450K, 
                                                 geneInf, 
                                                 type="PROMOTER")
  cpgOlsHg38450K[[i]] <- ol450K
  olEPIC <-  makeCpgGeneOverlap_MULT(tmpInf$ensembl_gene_id, 
                                                  hg38CpgLocsEPIC, 
                                                  geneInf, 
                                                  type="PROMOTER", )
  cpgOlsHg38EPIC[[i]] <- olEPIC
}
names(cpgOlsHg38450K) <- names(keepPromoInfLst)
names(cpgOlsHg38EPIC) <- names(keepPromoInfLst)

################################################################################
# Load SigPromo files
################################################################################

################################################################################
################################################################################
# Load datasets
message("Loading datasets used in study")
################################################################################
################################################################################
# Define array-type of dataframes
arrTypeLst <- list("TRAIN" = "EPIC", "GSE133556" = "EPIC", 
                   "GSE211686_1" = "450K", "GSE211686_2" = "EPIC",
                   "GSE226823" = "450K", "GSE263434" = "EPIC", 
                   "GSE267068" = "EPIC", "GSE155760" = "EPIC",
                   "GSE185008_1" = "450K", "GSE185008_2" = "EPIC", "GSE51820" = "450K")

################################################################################
# Load Pheno
################################################################################

################################################################################
# All
################################################################################
gse51820_pheno <- read.csv(paste(dataPath, "GSE51820/gse51820_pheno.csv",sep=""))
gse226823_pheno <- read.csv(paste(dataPath, "GSE226823/gse226823_pheno.csv", sep=""))
trainPheno <- read.csv(paste(dataPath, "Train/trainPheno.csv", sep=""), row.names = 1)

################################################################################
# CCC
################################################################################
gse185008_1_pheno <- read.csv(paste(dataPath, "GSE185008/GSE185008_1_pheno.csv", sep=""), row.names=1)
gse185008_2_pheno <- read.csv(paste(dataPath, "GSE185008/GSE185008_2_pheno.csv", sep=""), row.names=1)

################################################################################
# HGSC
################################################################################
gse133556_pheno <- read.csv(paste(dataPath, "GSE133556/GSE133556_pheno.csv",sep=""), row.names=1)
gse211686_1_pheno <- read.csv(paste(dataPath, "GSE211686/GSE211686_1_pheno.csv", sep=""), row.names=1)
gse211686_2_pheno <- read.csv(paste(dataPath, "GSE211686/GSE211686_2_pheno.csv", sep=""), row.names=1)
gse263434_pheno <- read.csv(paste(dataPath, "GSE263434/GSE263434_pheno.csv", sep=""), row.names=1)
gse267068_pheno <- read.csv(paste(dataPath, "GSE267068/GSE267068_pheno.csv", sep=""), row.names=1)
gse155760_pheno <- read.csv(paste(dataPath, "GSE155760/GSE155760_pheno.csv", sep=""), row.names=1)

################################################################################
# Load Beta 
# Warning: Unless you have signficant RAM available, load for one histotype at a time
################################################################################

################################################################################
# All
################################################################################
gse51820_beta <- read.csv(paste(dataPath, "GSE51820/GSE51820_Beta.csv", sep=""), row.names=1)
gse226823_beta <- read.csv(paste(dataPath, "GSE226823/GSE226823_Beta.csv", sep=""), row.names = 1)
trainBeta <- read.csv(paste(dataPath, "Train/Train_Beta.csv", sep=""), row.names = 1)

################################################################################
# CCC
################################################################################
gse185008_1_beta <- read.csv(paste(dataPath, "GSE185008_1/GSE185008_1_Beta.csv", sep=""), row.names=1)
gse185008_2_beta <- read.csv(paste(dataPath, "GSE185008_2/GSE185008_2_Beta.csv", sep=""), row.names=1)

################################################################################
# HGSC
################################################################################
gse133556_beta <- read.csv(paste(dataPath, "GSE133556/GSE133556_Beta.csv",sep=""), row.names=1)
gse211686_1_beta <- read.csv(paste(dataPath, "GSE211686_1/GSE211686_1_Beta.csv", sep=""), row.names=1)
gse211686_2_beta <- read.csv(paste(dataPath, "GSE211686_2/GSE211686_2_Beta.csv", sep=""), row.names=1)
gse263434_beta <- read.csv(paste(dataPath, "GSE263434/GSE263434_Beta.csv", sep=""), row.names=1)
gse267068_beta <- read.csv(paste(dataPath, "GSE267068/GSE267068_Beta.csv", sep=""), row.names=1)
gse155760_beta <- read.csv(paste(dataPath, "GSE155760/GSE155760_Beta.csv", sep=""), row.names=1)

allChrtPheno <- list("TRAIN"= trainPheno, 
                     "GSE226823" = gse226823_pheno, 
                     "GSE51820" = gse51820_pheno,
                     "GSE211686_1" = gse211686_1_pheno, 
                     "GSE211686_2" = gse211686_2_pheno, 
                     "GSE133556" = gse133556_pheno, 
                     "GSE263434" =  gse263434_pheno, 
                     "GSE267068" = gse267068_pheno, 
                     "GSE155760" = gse155760_pheno, 
                     "GSE185008_1" = gse185008_1_pheno, 
                     "GSE185008_2" =gse185008_2_pheno)

allChrtBeta <- list("TRAIN"= trainBeta, 
                    "GSE226823" =  gse226823_beta, 
                    "GSE51820" = gse51820_beta, 
                    "GSE211686_1" = gse211686_1_beta,
                    "GSE211686_2" = gse211686_2_beta, 
                    "GSE133556" = gse133556_beta, 
                    "GSE263434" =  gse263434_beta, 
                    "GSE267068" = gse267068_beta, 
                    "GSE155760" = gse155760_beta, 
                    "GSE185008_1" = gse185008_1_beta, 
                    "GSE185008_2" = gse185008_2_beta)

for(i in 1:length(allChrtBeta)){
  tmpB <- allChrtBeta[[i]]
  tmpP <- allChrtPheno[[i]]
  chrtFilt <- makePhenoBetaFilter(inpB = tmpB,
                                  inpP = tmpP)
  allChrtBeta[[i]] <- chrtFilt[[1]]
  allChrtPheno[[i]] <- chrtFilt[[2]]
}

# save.image(paste(wd, "extComp.RData", sep="/"))
load(paste(wd, "extComp.RData", sep="/"))

################################################################################
# Make CCC vectors
################################################################################
cccNames <- c("GSE185008_1", "GSE185008_2")
cccBetaLst <- allChrtBeta[cccNames]
cccPhenoLst <- allChrtPheno[cccNames]

# # Make PCA-plots for CCC sets
for(i in 1:length(cccBetaLst)){
  tmpB <- cccBetaLst[[i]]
  tmpSid <- names(cccBetaLst)[i]
  tmpPheno <- cccPhenoLst[[tmpSid]]
  makeMethPcaPlot(tmpB, tmpPheno, nameBool = tmpSid, noProbes = 1000)
  makeMethBetaFreqPlt(tmpB, tmpPheno, inpSid = tmpSid)
}

cccBetaCompLst <- allChrtBeta[c(cccNames, "GSE61820", "GSE226823", "TRAIN")]
cccPhenoCompLst <- allChrtPheno[c(cccNames, "GSE61820", "GSE226823", "TRAIN")]
# Get mean of genes
cccMeanLst <- makeGeneMeanLst(inpBetaLst = cccBetaCompLst, 
                              inpPhenoLst= cccPhenoCompLst, 
                              inpCpgLocsEPIC = cpgOlsHg38EPIC, 
                              inpCpgLocs450K = cpgOlsHg38450K,
                              arrTypeLst = arrTypeLst,
                              inpH = "CCC")
cccStdRank <- rankGeneMeanVar(inpGeneLst = cccMeanLst, 
                              inpGeneInf = geneInf)
write.csv(cccStdRank, 
          paste(outPath, "CCC/cccAllExtStdRank.csv", sep=""))


################################################################################
# Make HGSC vectors
################################################################################
hgscNames <- c("GSE133556", "GSE211686_1", "GSE211686_2", "GSE263434", "GSE267068", "GSE155760")
hgscBetaLst <- allChrtBeta[hgscNames]
hgscPhenoLst <- allChrtPheno[hgscNames]

# Make PCA-plots for HGSC sets
for(i in 1:length(hgscBetaLst)){
  tmpB <- hgscBetaLst[[i]]
  tmpSid <- names(hgscBetaLst)[i]
  tmpPheno <- hgscPhenoLst[[tmpSid]]
  makeMethPcaPlot(tmpB, tmpPheno, nameBool = tmpSid, noProbes = 1000)
  makeMethBetaFreqPlt(tmpB, tmpPheno, inpSid = tmpSid)
}

hgscBetaCompLst <- allChrtBeta[c(hgscNames, "GSE61820", "GSE226823", "TRAIN")]
hgscPhenoCompLst <- allChrtPheno[c(hgscNames, "GSE61820", "GSE226823", "TRAIN")]
# Get mean of genes
hgscMeanLst <- makeGeneMeanLst(inpBetaLst = hgscBetaCompLst, 
                               inpPhenoLst= hgscPhenoCompLst, 
                               inpCpgLocsEPIC = cpgOlsHg38EPIC, 
                               inpCpgLocs450K = cpgOlsHg38450K,
                               arrTypeLst = arrTypeLst,
                               inpH = "HGSC")
hgscStdRank <- rankGeneMeanVar(inpGeneLst = hgscMeanLst, 
                              inpGeneInf = geneInf)
write.csv(hgscStdRank, 
          paste(outPath, "HGSC/hgscAllExtStdRank.csv", sep=""))

# ################################################################################
# # Make EC vectors
# ################################################################################
# ecNames <- c("TRAIN", "GSE226823", "GSE51820", "GSE263434")
# ecBetaLst <- list(trainBeta, gse226823_beta, gse51820_beta, gse263434_beta)
# ecPhenoLst <- list(trainPheno, gse226823_pheno, gse51820_pheno, gse263434_pheno)
# names(ecBetaLst) <- ecNames
# names(ecPhenoLst) <- ecNames
# 
# ################################################################################
# # Make MC vectors
# ################################################################################
# mcNames <- c("TRAIN", "GSE51820")
# mcBetaLst <- list(trainBeta,  gse51820_beta)
# mcPhenoLst <- list(trainPheno, gse51820_pheno)
# names(mcBetaLst) <- mcNames
# names(mcPhenoLst) <- mcNames

mcMeanLst <- makeGeneMeanLst(inpBetaLst = allChrtBeta[c("GSE51820", "TRAIN")], 
                             inpPhenoLst = allChrtPheno[c("GSE51820", "TRAIN")], 
                             inpCpgLocsEPIC = cpgOlsHg38EPIC, 
                             inpCpgLocs450K = cpgOlsHg38450K,
                             arrTypeLst = arrTypeLst,
                             inpH = "MC")
mcStdRank <- rankGeneMeanVar(inpGeneLst = mcMeanLst, 
                               inpGeneInf = geneInf)
write.csv(mcStdRank, 
          paste(outPath, "MC/mcAllExtStdRank.csv", sep=""))

################################################################################
################################################################################
# Plot model in external datasets
message("Plotting individual genes in all cohorts")
################################################################################
################################################################################

################################################################################
# Create annotation for predictive classification model 
################################################################################
stepWiseSigPromoModel_450K <- read.csv(paste(outPath, "cpgModelStepWise_450K_sigPromo.csv", sep=""), row.names = 1)
stepWiseSigPromoModelHists_450K <- rownames(stepWiseSigPromoModel_450K)
stepWiseSigPromoModelCpGs_450K <- stepWiseSigPromoModel_450K[,1]
stepWiseSigPromoModel_450K <- stepWiseSigPromoModelCpGs_450K
names(stepWiseSigPromoModel_450K) <- stepWiseSigPromoModelHists_450K

cpgAnnDf <- data.frame(matrix(nrow=length(stepWiseSigPromoModel_450K), ncol=2))
colnames(cpgAnnDf) <- c("rowID", "Histotype")
cpgAnnDf$rowID <- stepWiseSigPromoModel_450K
cpgAnnDf$Histotype <- names(stepWiseSigPromoModel_450K)
cpgAnnDf$Histotype <- gsub('[0-9]+', '', cpgAnnDf$Histotype)

################################################################################
# Plot model
################################################################################
makeHSPModelPlotLstBased(inpBetaLst = allChrtBeta, 
                         inpPhenoLst = allChrtPheno, 
                         inpCpgAnnDf = cpgAnnDf, 
                         fileExt = "final_multclass_model_allExt")

################################################################################
# Plot individual genes of interest
################################################################################
compBetaLst <- allChrtBeta[c("GSE51820", 
                              "GSE226823",
                              "TRAIN")]
compPhenoLst <- allChrtPheno[c("GSE51820", 
                              "GSE226823",
                              "TRAIN")]

# LRRC41
lrrc41Coords <- rownames(cpgOlsHg38EPIC$CCC$ENSG00000132128) 
cpgAnnDf_lrrc41 <- data.frame(matrix(nrow=length(lrrc41Coords), ncol=2))
colnames(cpgAnnDf_lrrc41) <- c("rowID", "Histotype")
cpgAnnDf_lrrc41$rowID <- lrrc41Coords 
cpgAnnDf_lrrc41$Histotype <- "CCC"
makeHSPModelPlotLstBased(inpBetaLst = compBetaLst, 
                         inpPhenoLst = compPhenoLst, 
                         inpCpgAnnDf = cpgAnnDf_lrrc41, 
                         fileExt = "LRRC41")

# PYY
pyyCoords <- rownames(cpgOlsHg38EPIC$CCC$ENSG00000131096) 
cpgAnnDf_pyy <- data.frame(matrix(nrow=length(pyyCoords), ncol=2))
colnames(cpgAnnDf_pyy) <- c("rowID", "Histotype")
cpgAnnDf_pyy$rowID <- pyyCoords 
cpgAnnDf_pyy$Histotype <- "CCC"
makeHSPModelPlotLstBased(inpBetaLst = compBetaLst, 
                         inpPhenoLst = compPhenoLst, 
                         inpCpgAnnDf = cpgAnnDf_pyy, 
                         fileExt = "PYY")

# CMTM2
cmtm2Coords <- rownames(cpgOlsHg38EPIC$HGSC$ENSG00000140932) 
cpgAnnDf_cmtm2 <- data.frame(matrix(nrow=length(cmtm2Coords ), ncol=2))
colnames(cpgAnnDf_cmtm2) <- c("rowID", "Histotype")
cpgAnnDf_cmtm2$rowID <- cmtm2Coords 
cpgAnnDf_cmtm2$Histotype <- "HGSC"
makeHSPModelPlotLstBased(inpBetaLst = compBetaLst, 
                         inpPhenoLst = compPhenoLst, 
                         inpCpgAnnDf = cpgAnnDf_cmtm2, 
                         fileExt = "CMTM2")

# CMKLR1
cmklr1Coords <- rownames(cpgOlsHg38EPIC$CCC$ENSG00000174600) 
cpgAnnDf_cmklr1 <- data.frame(matrix(nrow=length(cmklr1Coords), ncol=2))
colnames(cpgAnnDf_cmklr1) <- c("rowID", "Histotype")
cpgAnnDf_cmklr1$rowID <- cmklr1Coords 
cpgAnnDf_cmklr1$Histotype <- "CCC"
makeHSPModelPlotLstBased(inpBetaLst = compBetaLst, 
                         inpPhenoLst = compPhenoLst, 
                         inpCpgAnnDf = cpgAnnDf_cmklr1, 
                         fileExt = "CMKLR1")

################################################################################
################################################################################
# External validation of HSP method (runs analysis on external datasets)
message("Running SigPromo method on external data")
################################################################################
################################################################################
allSigPromoFiles <- list.files(path = outPath, 
                               pattern = "promoters_ranked.csv",
                               recursive = TRUE,
                               full.names = TRUE)
allSigPromo <- makeDfLst(allSigPromoFiles)
allSigPromo  <- allSigPromo[!names(allSigPromo) %in% "EC"]

################################################################################
# Perform promoter ranking of signifcant regions in external cohorts
# We use lower thresholds to see what their values are, and not be filtered away
# This is necessary as the CpG coverage is lower in 450K arrays
################################################################################
extBetaLst <- allChrtBeta[c("GSE51820", "GSE226823")]
extPhenoLst <- allChrtPheno[c("GSE51820", "GSE226823")]
extPromoRankLst <- list() 
for(i in 1:length(extBetaLst)){
  tmpBetaMatches <- list()
  tmpExtBeta <- extBetaLst[[i]]
  tmpExtM <- log2(tmpExtBeta /(1-tmpExtBeta ))
  tmpSid <- names(extBetaLst)[i]
  tmpExtPheno <- extPhenoLst[[tmpSid]]
  extRankLst <- list()
  message(paste("Ranking promoter regions for: ", tmpSid, sep=""))
  for(j in 1:length(cpgOlsHg38450K)){
    tmpH <- names(cpgOlsHg38450K)[j]
    tmpGB <- cpgOlsHg38450K[[j]]
    tmpRank <- makePromoBetaRank_V3(inpRegLocLst = tmpGB,
                         inpPheno = tmpExtPheno, 
                         inpGeneInf = geneInf,
                         allBetas = tmpExtBeta,
                         focGrp = tmpH,
                         minCpg = 2, 
                         sigCpg = 2, 
                         noCat = length(table(tmpExtPheno$Histotype))-1,
                         bCut = 0.2, 
                         pCut = 0.05,
                         cpgCutFreq = 0.25, 
                         distType = "TRIMEAN", 
                         brownBool = TRUE)
    if(is.null(tmpRank)){
      next()
    }
    extRankLst[[j]] <- tmpRank
    names(extRankLst)[j] <- tmpH
  }
  extPromoRankLst[[i]] <- extRankLst
  names(extPromoRankLst)[i] <- tmpSid
}
  
makeCsvSave(extPromoRankLst$GSE226823, 
            "SigPromo_Ext_GSE226823")
makeCsvSave(extPromoRankLst$GSE51820, 
            "SigPromo_Ext_GSE51820")

################################################################################
################################################################################
# External validation plots
message("Plotting results from external validation")
################################################################################
################################################################################
deviance_files <- list.files(path = outPath, pattern = "AllExtStdRank",
                             recursive = TRUE,
                             full.names = TRUE)
devianceRankLst <- makeDfLst(deviance_files)

################################################################################
# Plot standard deviance from training cohort
################################################################################
makeDevBoxPlt(devianceRankLst, 
              fileExt = "allExtDeviance")

################################################################################
################################################################################
# Plot CpG-sites from cohorts in list side by side 
message("Creating comparative boxplots for genes of interest")
################################################################################
################################################################################

hg38GenePromoterCpgs <- read.csv(paste(outPath, "hg38GenePromoterCpgs.csv", sep=""), row.names = 1)

allSigPromoDfLst <- list() 
allSigPromoCpGSumDfLst <- list()
for(i in 1:length(allSigPromo)){
  tmpPromoGene <- allSigPromo[[i]]
  tmpName <- names(allSigPromo)[i]
  tmpPromo <- hg38GenePromoterCpgs[which(hg38GenePromoterCpgs$Gene %in% tmpPromoGene$ensembl_gene_id),]
  tmpInf <- geneInf[which(geneInf$ensembl_gene_id %in% tmpPromoGene$ensembl_gene_id),]
  tmpCpgPos <- makePromoLst(tmpInf, 
                            tmpPromo,
                            hg38CpgLocsEPIC)
  allSigPromoCpGSumDfLst[[i]] <- tmpPromo
  names(allSigPromoCpGSumDfLst)[i] <- tmpName
  allSigPromoDfLst[[i]] <- tmpCpgPos
  names(allSigPromoDfLst)[i] <- tmpName
}

cccSigGeneInf <- geneInf[match(allSigPromo$CCC$ensembl_gene_id, geneInf$ensembl_gene_id),]
hgscSigGeneInf <- geneInf[match(allSigPromo$HGSC$ensembl_gene_id, geneInf$ensembl_gene_id),]
mcSigGeneInf <- geneInf[match(allSigPromo$MC$ensembl_gene_id, geneInf$ensembl_gene_id),]

################################################################################
# Plot significant promoter regions in CCC cohorts
################################################################################
cccCompPlotLstBox <- makeExtCompCpgBoxPlot(inpCpgGeneLst = allSigPromoDfLst$CCC, 
                                           inpBetaLst = cccBetaCompLst, 
                                           inpPhenoLst = cccPhenoCompLst, 
                                           inpEPIC = cpgOlsHg38EPIC, 
                                           inp450k = cpgOlsHg38450K, 
                                           inpH = "CCC", 
                                           inpGeneInf = geneInf)

################################################################################
# Plot significant promoter regions in HGSC cohorts
################################################################################
hgscCompPlotLstBox <- makeExtCompCpgBoxPlot(inpCpgGeneLst = allSigPromoDfLst$HGSC, 
                                            inpBetaLst = hgscBetaCompLst, 
                                            inpPhenoLst = hgscPhenoCompLst, 
                                            inpEPIC = cpgOlsHg38EPIC, 
                                            inp450k = cpgOlsHg38450K, 
                                            inpH = "HGSC", 
                                            inpGeneInf = geneInf)

################################################################################
# Plot significant promoter regions in MC cohorts
################################################################################
mcCompPlotLstBox <- makeExtCompCpgBoxPlot(inpCpgGeneLst = allSigPromoDfLst$MC, 
                                          inpBetaLst = list("TRAIN" = trainBeta, 
                                                            "GSE51820" = gse51820_beta), 
                                          inpPhenoLst = list("TRAIN" = trainPheno, 
                                                             "GSE51820" = gse51820_pheno), 
                                          inpEPIC = cpgOlsHg38EPIC, 
                                          inp450k = cpgOlsHg38450K, 
                                          inpH = "MC", 
                                          inpGeneInf = geneInf)

################################################################################
# Create list of plots for in training cohort
################################################################################
cccBoxLst <- makeGeneLstCpgBoxPlot(inpGenes = geneInf[which(geneInf$ensembl_gene_id %in% allSigPromo$CCC$ensembl_gene_id),],
                                   inpBeta = trainBeta, 
                                   inpPheno = trainPheno, 
                                   cpgInp = hg38CpgLocsEPIC, 
                                   geneInfInp = geneInf, 
                                   pltBool = FALSE, 
                                   promoBool = TRUE)

hgscBoxLst <- makeGeneLstCpgBoxPlot(inpGenes = geneInf[which(geneInf$ensembl_gene_id %in% allSigPromo$HGSC$ensembl_gene_id),],
                                    trainBeta, 
                                    trainPheno, 
                                    hg38CpgLocsEPIC, 
                                    geneInf, 
                                    pltBool = FALSE, 
                                    promoBool = TRUE)

mcBoxLst <- makeGeneLstCpgBoxPlot(geneInf[geneInf$ensembl_gene_id %in% allSigPromo$MC$ensembl_gene_id,], 
                                  trainBeta, 
                                  trainPheno, 
                                  hg38CpgLocsEPIC, 
                                  geneInf, 
                                  pltBool = FALSE, 
                                  promoBool = TRUE)

################################################################################
# Plot training cohort and external cohorts next to one another
################################################################################
boxCompLst <- list("CCC" = cccCompPlotLstBox,
                   "HGSC"= hgscCompPlotLstBox,
                   "MC"= mcCompPlotLstBox)

boxIntLst <- list("CCC" = cccBoxLst,
                  "HGSC" = hgscBoxLst, 
                  "MC" = mcBoxLst)

mergeBoxLst <- list()
for(i in 1:length(boxIntLst)){
  tmpBoxLst <- boxIntLst[[i]]
  tmpH <- names(boxCompLst)[i]
  tmpCompLst <- boxCompLst[[tmpH]]
  tmpHLst <- list()
  outDir <- paste(plotPath, "cpgBoxCompPlot/", tmpH, sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  for(j in 1:length(tmpBoxLst)){
    tmpPlt <-tmpCompLst[[j]]
    tmpSym <- names(tmpCompLst)[j]
    if(!tmpSym %in% names(tmpBoxLst)){
      next()
    }else{
      tmpTrainPlt <- tmpBoxLst[[tmpSym]]
      combPlt <- cowplot::plot_grid(tmpTrainPlt, tmpPlt, ncol=2)
      outFile <- paste(outDir, "/", tmpSym, "_compCpgBoxPlot_allExt.pdf",sep="")
      cowplot::save_plot(outFile, combPlt, base_width=24, base_height=16)
      tmpHLst[[j]] <- combPlt
      names(tmpHLst)[j] <- tmpSym
    }
  }
  mergeBoxLst[[i]] <- tmpHLst
  names(mergeBoxLst)[i] <- tmpH
}

################################################################################
################################################################################
# Plot CpG-sites from cohorts in list side by side 
message("Creating comparative dotplots for genes of interest")
################################################################################
################################################################################

# Create comparative dot-plots for best separated promoter genes
cccCompPlotLstDot <- makeExtCompCpgDotPlot(inpCpgGeneLst = allSigPromoDfLst$CCC, 
                                           inpBetaLst = cccBetaCompLst, 
                                           inpPhenoLst = cccPhenoCompLst, 
                                           inpEPIC = cpgOlsHg38EPIC, 
                                           inp450k = cpgOlsHg38450K, 
                                           inpH = "CCC", 
                                           inpGeneInf = geneInf)

hgscCompPlotLstDot <- makeExtCompCpgDotPlot(inpCpgGeneLst = allSigPromoDfLst$HGSC, 
                                            inpBetaLst = hgscBetaCompLst, 
                                            inpPhenoLst = hgscPhenoCompLst, 
                                            inpEPIC = cpgOlsHg38EPIC, 
                                            inp450k = cpgOlsHg38450K, 
                                            inpH = "HGSC", 
                                            inpGeneInf = geneInf)

mcCompPlotLstDot <- makeExtCompCpgDotPlot(inpCpgGeneLst = allSigPromoDfLst$MC, 
                                          inpBetaLst = list("TRAIN" = trainBeta, 
                                                            "GSE51820" = gse51820_beta), 
                                          inpPhenoLst = list("TRAIN" = trainPheno, 
                                                             "GSE51820" = gse51820_pheno), 
                                          inpEPIC = cpgOlsHg38EPIC, 
                                          inp450k = cpgOlsHg38450K, 
                                          inpH = "MC", 
                                          inpGeneInf = geneInf)

################################################################################
# Create list of plots for in training cohort
################################################################################

cccDotLst <- makeGeneLstCpgDotPlot(inpGenes = geneInf[geneInf$ensembl_gene_id %in% allSigPromo$CCC$ensembl_gene_id,],
                                   trainBeta, 
                                   trainPheno, 
                                   hg38CpgLocsEPIC, 
                                   geneInf, 
                                   pltBool = FALSE, 
                                   promoBool = TRUE)

hgscDotLst <- makeGeneLstCpgDotPlot(inpGenes = geneInf[geneInf$ensembl_gene_id %in% allSigPromo$HGSC$ensembl_gene_id,],
                                    trainBeta, 
                                    trainPheno, 
                                    hg38CpgLocsEPIC, 
                                    geneInf, 
                                    pltBool = FALSE, 
                                    promoBool = TRUE)

mcDotLst <- makeGeneLstCpgDotPlot(inpGenes = geneInf[geneInf$ensembl_gene_id %in% allSigPromo$MC$ensembl_gene_id,],
                                  trainBeta, 
                                  trainPheno, 
                                  hg38CpgLocsEPIC, 
                                  geneInf, 
                                  pltBool = FALSE, 
                                  promoBool = TRUE)

################################################################################
# Plot training cohort and external cohorts next to one another
################################################################################
dotCompLst <- list("CCC" = cccCompPlotLstDot, 
                   "HGSC" = hgscCompPlotLstDot, 
                   "MC" = mcCompPlotLstDot)

dotIntLst <- list("CCC" = cccDotLst, 
                  "HGSC" = hgscDotLst, 
                  "MC" = mcDotLst)

mergeDotLst <- list()
# Merge plots based on shared names (trainplot on left, compplot on right)
for(i in 1:length(dotIntLst)){
  tmpDotLst <- dotIntLst[[i]]
  tmpH <- names(dotCompLst)[i]
  tmpCompLst <- dotCompLst[[tmpH]]
  tmpHLst <- list()
  outDir <- paste(plotPath, "cpgDotCompPlot/", tmpH, sep="")
  ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir), recursive = TRUE), FALSE)
  for(j in 1:length(tmpDotLst)){
    tmpPlt <-tmpCompLst[[j]]
    tmpSym <- names(tmpCompLst)[j]
    if(!tmpSym %in% names(tmpDotLst)){
      next()
    }else{
      tmpTrainPlt <- tmpDotLst[[tmpSym]]
      combPlt <- cowplot::plot_grid(tmpTrainPlt, tmpPlt, ncol=2)
      outFile <- paste(outDir, "/", tmpSym, "_compCpgDotPlot_allExt.pdf.pdf",sep="")
      cowplot::save_plot(outFile, combPlt, base_width=24, base_height=16)
      tmpHLst [[j]] <- combPlt
      names(tmpHLst)[j] <- tmpSym
    }
  }
  mergeDotLst[[i]] <- tmpHLst
  names(mergeDotLst)[i] <- tmpH
}

################################################################################
# Distribution type analysis
# Warning: Long runtime
# message("Running distribution tpye analysis for full validation cohort")
################################################################################
allChrtDistTypes_FULL <- list()
for(i in 1:length(allChrtBeta)){
  tmpB <- allChrtBeta[[i]]
  tmpSid <- names(allChrtBeta)[i]
  message(paste("Running promoter methylation analysis for cohort:", tmpSid), sep=" ")
  tmpPheno <- allChrtPheno[[tmpSid]]
  distStats <- makeCpGDistTypeDf_MULT(inpBeta = tmpB,
                                      inpPheno = tmpPheno)
  allChrtDistTypes_FULL[[i]] <- distStats
  names(allChrtDistTypes_FULL)[i] <- tmpSid
}
makeCsvSave(allChrtDistTypes_FULL, "BetaDistributionTypes")

################################################################################
################################################################################
message("Running predictive classification using model")
################################################################################
################################################################################

################################################################################
# Run binary classification on external cohorts with few histotypes
################################################################################

source(paste(wd, "/Source/predScripts.R", sep=""))
source(paste(wd, "/Source/predPltScripts.R", sep=""))

stepWiseSigPromoModel_450K <- read.csv(paste(outPath, "cpgModelStepWise_450K_sigPromo.csv", sep=""), row.names = 1)
stepWiseSigPromoModelHists_450K <- rownames(stepWiseSigPromoModel_450K)
stepWiseSigPromoModelCpGs_450K <- stepWiseSigPromoModel_450K[,1]
stepWiseSigPromoModel_450K <- stepWiseSigPromoModelCpGs_450K
names(stepWiseSigPromoModel_450K) <- stepWiseSigPromoModelHists_450K
ensModel_450K <- read.csv(paste(outPath, "ensModelStepWise_450K_sigPromo.csv", sep=""), row.names = 1)

bestModel450KAnn <- data.frame(stepWiseSigPromoModel_450K,stepWiseSigPromoModelHists_450K, ensModel_450K)
colnames(bestModel450KAnn) <- c("CpG", "Histotype", "Gene")
write.csv(bestModel450KAnn, paste(outPath, "bestModelStepWise_450K_Ann.csv", sep=""))

allSigPromoFiles <- list.files(path = outPath, 
                               pattern = "sigPromoInfDf.csv",
                               recursive = TRUE,
                               full.names = TRUE)
allSigPromo <- makeDfLst(allSigPromoFiles)
allSigPromo  <- allSigPromo[!names(allSigPromo) %in% "EC"]


allSigPromoCpgPos <- list() 
allSigPromoCpgLst <- list()
for(i in 1:length(allSigPromo)){
  tmpPromoGene <- allSigPromo[[i]]
  tmpName <- names(allSigPromo)[i]
  tmpPromo <- promoter_cpgs[which(promoter_cpgs$Gene %in% tmpPromoGene$ensembl_gene_id),]
  tmpInf <- geneInf[which(geneInf$ensembl_gene_id %in% tmpPromoGene$ensembl_gene_id),]
  tmpCpgPos <- makePromoLst(tmpInf, 
                            tmpPromo, 
                            hg38CpgLocsEPIC)
  allSigPromoCpgLst[[i]] <- tmpPromo
  names(allSigPromoCpgLst)[i] <- tmpName
  allSigPromoCpgPos[[i]] <- tmpCpgPos
  names(allSigPromoCpgPos)[i] <- tmpName
}

extPredLst <- allChrtBeta[c("GSE133556", "GSE263434", "GSE155760")]
extPredPheno <- allChrtPheno[c("GSE133556", "GSE263434", "GSE155760")]
# Perform multicass classification using combined model
# Train on train-dataset, apply model directly on external dataset
bestModelExtPredRes_450K_MULT <- makeDirExtBinPred(inpExtBetaLst = extPredLst,
                                                    inpExtPhenoLst = extPredPheno,
                                                    inpCModel = stepWiseSigPromoModel_450K, 
                                                    refPheno = trainPheno,
                                                    refBeta = trainBeta,
                                                    sigCpgLocs = allSigPromoCpgLst[!names(allSigPromoCpgLst) %in% "EC"],
                                                    noEpochs = 100,
                                                    nFold = 5,
                                                    paramBool = FALSE)

write.csv(bestModelExtPredRes_450K_MULT$Merge, paste(outPath, "bestModelExtPredRes_450K_Mult_sumDf.csv", sep=""))
write.csv(data.frame(bestModelExtPredRes_450K_MULT$CpGModel), paste(outPath, "bestModelExtPredRes_450K_Mult_CpGModel.csv", sep=""))
write.csv(data.frame(bestModelExtPredRes_450K_MULT$ensModel), paste(outPath, "bestModelExtPredRes_450K_Mult_EnsModel.csv", sep=""))
makeCsvSave(bestModelExtPredRes_450K_MULT$cvAllLst, "bestModelExtPredRes_450K_Mult_allValsDf")

################################################################################
# Binary classificiation (One vs. All)
################################################################################
# Perform multicass classification using combined model
# Train on train-dataset, apply model directly on external dataset
bestModelExtPredRes_450K_Bin <- makeDirExtBinPred(inpExtBetaLst = extBetaLst,
                                                  inpExtPhenoLst = extPhenoLst,
                                                  inpCModel = stepWiseSigPromoModel_450K,
                                                  refPheno = trainPheno,
                                                  refBeta = trainBeta, 
                                                  sigCpgLocs = allSigPromoCpgLst[!names(allSigPromoCpgLst) %in% "EC"],
                                                  noEpochs = 100,
                                                  nFold = 5,
                                                  paramBool = FALSE)

