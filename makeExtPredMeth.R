################################################################################
################################################################################
################################################################################
# Script for predictive classification using external datasets
################################################################################
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
             "stringr", "tidyr", "tidymodels", "e1071", "glmnet", "MASS", "caret", 
             "kernlab", "ranger", "xgboost", "parsnip", "recipes", "rsample",
             "tune", "yardstick", "DescTools", "Boruta", "litteR", "fmsb", 
             "future", "future.apply", "MLmetrics")
miss_r_packs <- r_packs[!(r_packs %in% installed.packages()[,"Package"])]
if(length(miss_r_packs)){
  install.packages(miss_r_packs)
} 
invisible(lapply(r_packs, library, character.only = TRUE))

# Check bioconductor packages, install if missing
bioc_packs <- c("ChAMP", "GEOquery", "minfi", "minfiData",
                "BiocManager", "GenomicRanges", "GenomicFeatures", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                "minfi", "GEOquery")
miss_bc_packs <- bioc_packs[!(bioc_packs %in% installed.packages()[,"Package"])]
if(length(miss_bc_packs)){
  BiocManager::install(miss_bc_packs)
} 
invisible(lapply(bioc_packs, library, character.only = TRUE))


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
dataPath <- paste(wd, "/Data/", sep="/")
outPath <- paste(wd, "/Export/", sep="/")
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
nCores = round(0.5 * totCores)
options(future.globals.maxSize = 3000 * 1024^2)
plan(multisession, workers = nCores) 

################################################################################
################################################################################
# Source functions
################################################################################
################################################################################

source(paste(wd, "/Source/genScripts.R", sep=""))
source(paste(wd, "/Source/betaScripts.R", sep=""))
source(paste(wd, "/Source/pltScripts.R", sep=""))
source(paste(wd, "/Source/predScripts.R", sep=""))
source(paste(wd, "/Source/predPltScripts.R", sep=""))

################################################################################
# Load rdata
################################################################################
# load(paste(wd, "ExtPredData.RData",sep=""))

################################################################################
# Read in phenotype data
################################################################################
gse51820_pheno <- read.csv(paste(dataPath, "GSE51820/gse51820_pheno.csv",sep=""))
gse226823_pheno <- read.csv(paste(dataPath, "GSE226823/gse226823_pheno.csv", sep=""))
trainPheno <- read.csv(paste(dataPath, "Train/Trainpheno.csv", sep=""))

################################################################################
# Read in beta-matrices
################################################################################
gse51820_beta <- read.csv(paste(dataPath, "GSE51820/gse51820_Beta.csv",sep=""), row.names = 1)
gse226823_beta <- read.csv(paste(dataPath, "GSE226823/gse226823_Beta.csv",sep=""), row.names = 1)
trainBeta <- read.csv(paste(dataPath, "Train/Train_Beta.csv",sep=""), row.names = 1)
trainM <- log2(trainBeta/(1-trainBeta))

################################################################################
# Create data lists
################################################################################
keepSamps51820 <- intersect(gse51820_pheno$barcode, colnames(gse51820_beta))
gse51820_beta <- gse51820_beta[,keepSamps51820]
gse51820_pheno <- gse51820_pheno[gse51820_pheno$barcode %in% keepSamps51820, ]
gse51820_beta <- gse51820_beta[,match(gse51820_pheno$barcode, colnames(gse51820_beta))]
keepSamps226823 <- intersect(gse226823_pheno$barcode, colnames(gse226823_beta))
gse226823_beta <- gse226823_beta[,keepSamps226823]
gse226823_pheno <- gse226823_pheno[gse226823_pheno$barcode %in% keepSamps226823, ]
gse226823_beta <- gse226823_beta[,match(gse226823_pheno$barcode, colnames(gse226823_beta))]

extPhenoLst <- list("GSE51820" = gse51820_pheno,
                    "GSE226823" = gse226823_pheno)
extBetaLst <- list("GSE51820" = gse51820_beta,
                   "GSE226823" = gse226823_beta)

################################################################################
# Load Cpg-site location mapping data
################################################################################
hg38CpgLocsEPIC <- read.csv(paste(outPath, "/hg38CpgLocsEPIC.csv",sep=""), row.names = 1)
hg38CpgLocs450K <- read.csv(paste(outPath, "hg38CpgLocs450K.csv", sep=""), row.names = 1)
promoter_cpgs <- read.csv(paste(outPath, "hg38GenePromoterCpgs.csv", sep=""), row.names = 1)

cpgAnno <- read.csv(paste(outPath, "cpgAnno.csv", sep=""), row.names = 1)
geneInf <- read.csv(paste(outPath, "hg38GeneInf.csv", sep=""), row.names = 1)
# Create list of dataframes 
allPromoCpgDfLst <- makePromoLst(geneInf, 
                                 promoter_cpgs,
                                 hg38CpgLocsEPIC)
################################################################################
################################################################################
# Load data associated with significant promoters in main
message("Load data associated with significant promoters generated in main")
################################################################################
################################################################################

allSigPromoFiles <- list.files(path = outPath, 
                                   pattern = "promoters_ranked.csv",
                                   recursive = TRUE,
                                   full.names = TRUE)
allSigPromo <- makeDfLst(allSigPromoFiles)
allSigPromo  <- allSigPromo[!names(allSigPromo) %in% "EC"]

################################################################################
# Gene information dataframes (reduces loadtimes)
################################################################################
sigPromoInfFiles <- list.files(path = outPath, 
                               pattern = "sigPromoInfDf",
                               recursive = TRUE,
                               full.names = TRUE)
keepPromoInfLst <- makeDfLst(sigPromoInfFiles)

################################################################################
# Get cpg locations for best-probes
################################################################################
allSigPromoCpgPos <- list() 
allSigPromoCpgLst <- list()
for(i in 1:length(allSigPromo)){
  tmpPromoGene <- allSigPromo[[i]]
  tmpName <- names(allSigPromo)[i]
  tmpPromo <- promoter_cpgs[which(promoter_cpgs$Gene %in% tmpPromoGene$ensembl_gene_id),]
  tmpInf <- geneInf[which(geneInf$ensembl_gene_id %in% tmpPromoGene$ensembl_gene_id),]
  tmpCpgPos <- makePromoLst(tmpInf, tmpPromo, hg38CpgLocsEPIC)
  allSigPromoCpgLst[[i]] <- tmpPromo
  names(allSigPromoCpgLst)[i] <- tmpName
  allSigPromoCpgPos[[i]] <- tmpCpgPos
  names(allSigPromoCpgPos)[i] <- tmpName
}

sigPromoCpgFiles <- list.files(path = outPath, 
                               pattern = "sigPromoCpgVecs",
                               recursive = TRUE,
                               full.names = TRUE)
sigPromoCpgLst <- makeDfLst(sigPromoCpgFiles)
sigPromoCpgVec <- read.csv(paste(outPath, "allSigPromoCpgVec.csv", sep=""))

sigPromoCpGLocDfFiles <- list.files(path = outPath, 
                                    pattern = "sigPromoCpGLocDf",
                                    recursive = TRUE,
                                    full.names = TRUE)
keepPromoCpgPos <- makeDfLst(sigPromoCpGLocDfFiles)

################################################################################
# Get cpg location data for CpGs in promoter regions
################################################################################
sigPromoCpgLocsLst <- list()
for(i in 1:length(sigPromoCpgLst)){
  tmpH <- names(sigPromoCpgLst)[i]
  tmpCpgs <- sigPromoCpgLst[[i]]
  tmpCpgs <- as.vector(tmpCpgs[[1]])
  sigPromoCpgLst[[i]] <- tmpCpgs
  tmpLocs <- hg38CpgLocsEPIC[tmpCpgs,]
  sigPromoCpgLocsLst[[i]] <- tmpLocs
  names(sigPromoCpgLocsLst)[i] <- tmpH
}

################################################################################  
################################################################################
# Promoter genes coordinate mapping
message("Mapping coordinate to HSPs")
################################################################################
################################################################################

################################################################################
# Get gene-CpG overlaps for 450k/EPIC coordinates (useful if running script at later occurence)
################################################################################
cpgOlsHg38450K  <- list()
cpgOlsHg38EPIC <- list()
for(i in 1:length(keepPromoInfLst)){
  tmpH <- names(keepPromoInfLst)[i]
  cpgOlsHg38450K[[i]] <- makeCpgGeneOverlap_MULT(keepPromoInfLst[[i]]$ensembl_gene_id, 
                                                 hg38CpgLocs450K, 
                                                 geneInf, 
                                                 type="PROMOTER")
  names(cpgOlsHg38450K)[i] <- tmpH
  cpgOlsHg38EPIC[[i]] <-  makeCpgGeneOverlap_MULT(keepPromoInfLst[[i]]$ensembl_gene_id, 
                                                  hg38CpgLocsEPIC, 
                                                  geneInf, 
                                                  type="PROMOTER")
  names(cpgOlsHg38EPIC)[i] <- tmpH
}

################################################################################
################################################################################
# Stepwise Classification + multi and bin for using HSPs
message("Beginning predictive classification using HSPs")
################################################################################
################################################################################

################################################################################
# Filter training dataset to include CpG sites in both external datasets
# i.e. we keep genes present in both EPIC and 450K array and with represenation in our data
################################################################################
# Get rownames present in both external cohorts
extRowKeep <- intersect(rownames(extBetaLst[[1]]), rownames(extBetaLst[[2]]))
modRowKeep <- intersect(rownames(trainBeta), extRowKeep)
predBetaDf <- trainBeta[modRowKeep,]

cpgOls450KKeep <- list()
for(i in 1:length(cpgOlsHg38450K)){
  tmpLst <- cpgOlsHg38450K[[i]]
  nLst <- list()
  for(j in 1:length(tmpLst)){
    tmpDf <- tmpLst[[j]]
    # Only keep entries which have their rownames in our reference dataframe
    tmpDf <- tmpDf[which(rownames(tmpDf) %in% rownames(predBetaDf)),]
    # Remove genes with less then 3 CpG sites
    if(nrow(tmpDf) < 3){
      next()
    }else{
      nLst[[length(nLst)+1]] <- tmpDf
      names(nLst)[length(nLst)] <- names(tmpLst)[j]
    }
  }
  cpgOls450KKeep[[i]] <- nLst
  names(cpgOls450KKeep)[i] <- names(cpgOlsHg38450K)[i]
}

sigPromo450KFiltKeepLst <- list()
for(i in 1:length(cpgOls450KKeep)){
  tmpH <- names(cpgOls450KKeep)[i]
  tmpPromo <- cpgOls450KKeep[[i]]
  tmpInf <- allSigPromo[[tmpH]]
  tmpFilt <- tmpInf[which(tmpInf$ensembl_gene_id %in% names(tmpPromo)), ]
  sigPromo450KFiltKeepLst[[i]] <- tmpFilt
  names(sigPromo450KFiltKeepLst)[i] <- tmpH
}

makeCsvSave(sigPromo450KFiltKeepLst, 
            "sigPromo450KFiltKeep")

sigPromo450KFiltKeepLstFiles <- list.files(path = outPath, 
                                    pattern = "sigPromo450KFiltKeep",
                                    recursive = TRUE,
                                    full.names = TRUE)
sigPromo450KFiltKeepLst <- makeDfLst(sigPromo450KFiltKeepLstFiles)

################################################################################
################################################################################
# Separate script for selecting a model through stepwise iteration
message("Running stepwise model selection for multiclass model")
################################################################################
################################################################################

################################################################################
# Retrieve the best model combined through stepwise addition for each histotype individually
# Warning, extremely long runtime, can be adressed by lowering noEpochs
################################################################################
# stepWiseModel <- makeStepwiseMultModel(inpCpgDfLst = keepPromoCpgPos, 
#                                                refBeta = predBetaDf, 
#                                                refPheno = trainPheno,
#                                                inpCpgPos = cpgOls450KKeep,
#                                                paramBool = NULL,
#                                                noEpochs = 100,
#                                                minCpg = 2,
#                                                nFold = 5)

stepWiseModel <- makeStepwiseMultModel_MULT(inpCpgDfLst = keepPromoCpgPos, 
                           refBeta = predBetaDf, 
                           refPheno = trainPheno,
                           inpCpgPos = cpgOls450KKeep,
                           paramBool = NULL,
                           noEpochs = 5,
                           minCpg = 2,
                           nFold = 5)
#cpgModelStepwise <- stepWiseModel[[1]]
#ensModelStepwise <- stepWiseModel[[2]]
#statDfStepwise <- stepWiseModel[[3]]
# write.csv(statDfStepwise, paste(outPath, "statDfStepWise_450K_sigPromo.csv", sep=""))

# Retrieve the best combination of the created top n stepwise models using grid search
bestModel <- makeGridBestModel(stepWiseModel$modelLst, 
                               inpBeta = trainBeta, 
                               inpPheno = trainPheno, 
                               sigCpgLocs = allSigPromoCpgLst[!names(allSigPromoCpgLst) %in% "EC"])
annBestModel <- getAnnotatedCpgModel(bestModel, allSigPromoCpgLst[!names(allSigPromoCpgLst) %in% "EC"])

write.csv(data.frame(annBestModel$CpG), 
          paste(outPath, "cpgModelStepWise_450K_sigPromo.csv", sep=""))
write.csv(data.frame(annBestModel$Ens), 
          paste(outPath, "ensModelStepWise_450K_sigPromo.csv", sep=""))

################################################################################
# Load best model
################################################################################
stepWiseSigPromoModel_EPIC <- read.csv(paste(outPath, "cpgModelStepWise_EPIC_sigPromo.csv", sep=""), row.names = 1)
stepWiseSigPromoModelHists_EPIC <- rownames(stepWiseSigPromoModel_EPIC)
stepWiseSigPromoModelCpGs_EPIC <- stepWiseSigPromoModel_EPIC[,1]
stepWiseSigPromoModel_EPIC <- stepWiseSigPromoModelCpGs_EPIC
names(stepWiseSigPromoModel_EPIC) <- stepWiseSigPromoModelHists_EPIC

stepWiseSigPromoModel_450K <- read.csv(paste(outPath, "cpgModelStepWise_450K_sigPromo.csv", sep=""), row.names = 1)
stepWiseSigPromoModelHists_450K <- rownames(stepWiseSigPromoModel_450K)
stepWiseSigPromoModelCpGs_450K <- stepWiseSigPromoModel_450K[,1]
stepWiseSigPromoModel_450K <- stepWiseSigPromoModelCpGs_450K
names(stepWiseSigPromoModel_450K) <- stepWiseSigPromoModelHists_450K
ensModel_450K <- read.csv(paste(outPath, "ensModelStepWise_450K_sigPromo.csv", sep=""), row.names = 1)

bestModel450KAnn <- data.frame(stepWiseSigPromoModel_450K,stepWiseSigPromoModelHists_450K, ensModel_450K)
colnames(bestModel450KAnn) <- c("CpG", "Histotype", "Gene")
write.csv(bestModel450KAnn, paste(outPath, "bestModelStepWise_450K_Ann.csv", sep=""))
################################################################################
################################################################################
message("Running predictive classification using model")
################################################################################
################################################################################

################################################################################
# Multiclass classification
################################################################################
# Perform multicass classification using combined model
# Train on train-dataset, apply model directly on external dataset
bestModelExtPredRes_450K_MULT <- makeDirExtMultPred(inpExtBetaLst = extBetaLst,
                                                      inpExtPhenoLst = extPhenoLst,
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

write.csv(bestModelExtPredRes_450K_Bin$Merge, paste(outPath, "bestModelExtPredRes_450K_Bin_sumDf.csv",sep=""))
write.csv(data.frame(bestModelExtPredRes_450K_Bin$CpGModel), paste(outPath, "bestModelExtPredRes_450K_Bin_CpGModel.csv", sep=""))
write.csv(data.frame(bestModelExtPredRes_450K_Bin$ensModel), paste(outPath, "bestModelExtPredRes_450K_Bin_EnsModel.csv", sep=""))
makeCsvSave(bestModelExtPredRes_450K_Bin$cvAllLst, "bestModelExtPredRes_450K_Bin_allValsDf")

################################################################################
################################################################################
message("Plotting classification results")
################################################################################
################################################################################

################################################################################
# Plot multiclass predictive classification results
################################################################################
makePredStatPlt(inpPredDfLst = bestModelExtPredRes_450K_MULT$MeanCvLst, 
                remHist = "EC",
                pltCols = c("SID", "Histotype", "Sensitivity", "Specificity", "Accuracy", "F1Score"),
                fileAdd = "bestModelSigPromo450K_Mult",
                multParamBool = NULL)

# Plot the full epochs for the predictive classification
makePredEpochPlt(inpPredDfLst = bestModelExtPredRes_450K_MULT$cvAllLst, 
                 remHist = "EC",
                 pltCols = c("SID", "Histotype", "Sensitivity", "Specificity", "Accuracy", "F1Score"),
                 fileAdd = "bestModelSigPromo450K_Mult",
                 multParamBool = NULL)

################################################################################
# Horizontal oritentation (easier to interpret)
################################################################################
makePredStatPltFlipped(inpPredDfLst = bestModelExtPredRes_450K_MULT$cvAllLst, 
                       remHist = "EC",
                       pltCols = c("SID", "Histotype", "Sensitivity", "Specificity", "Accuracy", "F1Score"),
                       fileAdd = "bestModelSigPromo450K_Mult_Horisontal",
                       multParamBool = NULL,
                       binBool = NULL)

makePredStatPltFlipped(inpPredDfLst = bestModelExtPredRes_450K_MULT$cvAllLst, 
                       remHist = "EC",
                       pltCols = c("SID", "Histotype", "Acc_Overall", "AUC", "PRAUC"),
                       fileAdd = "bestModelSigPromo450K_Mult_SumStats_Horisontal",
                       multParamBool = TRUE)

################################################################################
# Plot binary predictive classification results
################################################################################

makePredStatPltFlipped(inpPredDfLst = bestModelExtPredRes_450K_Bin$cvAllLst, 
                       remHist = "EC",
                       pltCols = c("SID", "Histotype", "Sensitivity", "Specificity", "Accuracy", "F1Score"),
                       fileAdd = "bestModelSigPromo450K_Bin_Horisontal",
                       multParamBool = NULL,
                       binBool = NULL)

################################################################################
################################################################################
# Plot model from classification
message("Plotting classification model")
################################################################################
################################################################################

################################################################################
# Create DF required by function 
################################################################################
cpgAnnDf <- data.frame(matrix(nrow=length(stepWiseSigPromoModel_450K), ncol=2))
colnames(cpgAnnDf) <- c("rowID", "Histotype")
cpgAnnDf$rowID <- stepWiseSigPromoModel_450K
cpgAnnDf$Histotype <- names(stepWiseSigPromoModel_450K)
cpgAnnDf$Histotype <- gsub('[0-9]+', '', cpgAnnDf$Histotype)

################################################################################
# Plot model in Training cohort and external cohorts
################################################################################
makeHSPModelPlot(trainBeta, 
                 trainPheno, 
                 cpgAnnDf, 
                 fileExt = "Train_stepWiseSigPromoModel_450K")

for(i in 1:length(extBetaLst)){
  tmpSid <- names(extBetaLst)[[i]]
  tmpB <- extBetaLst[[tmpSid]]
  tmpP <- extPhenoLst[[tmpSid ]]
  makeHSPModelPlot(tmpB, 
                   tmpP, 
                   cpgAnnDf, 
                   fileExt = paste(tmpSid, "stepWiseSigPromoModel_450K", sep="_"))
}

save.image(paste(wd, "/ExtPredData.RData",sep=""))