################################################################################
################################################################################
################################################################################
# Main script for study
# Generates significant differentially methylated promoter regions for genes
################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
# Load packages
################################################################################
################################################################################

################################################################################
# Install packages used in scripts
################################################################################
bcPacks <- c("biomaRt", "GenomicRanges","GenomicFeatures", "TxDb.Hsapiens.UCSC.hg38.knownGene", "minfi", "GEOquery", "goseq", "conumee", "EmpiricalBrownsMethod", "ComplexHeatmap")
miss_bc_packs <- bcPacks[!(bcPacks %in% installed.packages()[,"Package"])]
if(length(miss_bc_packs)){
  BiocManager::install(miss_bc_packs)
} 
rPacks <- c("devtools","tidymodels", "rtracklayer", "ggplot2", "RColorBrewer", "ggpubr", "gridExtra", "dplyr", "data.table", "stringr", "tidyr", "viridis", "tibble",
            "forcats", "cowplot", "colorRamp2", "tibble", "glmnet", "ROCR", "MASS", "caret", "kernlab", "ranger", "xgboost", "adabag", "parsnip", "workflows", "recipes", "rsample", "tidyr", "tune",
            "yardstick", "tuneRanger", "DescTools", "Boruta", "patchwork", "htmltools", "future.apply")
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
# install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)

################################################################################
################################################################################
# Load global parameters
################################################################################
################################################################################

histotypes = c("CCC","EC", "HGSC", "MC")
focusGrp <- "Histotype"
alphaV <- 0.05

################################################################################
# Define plot color parameters
################################################################################
colProf <- viridis(50)
catColProf <-  list(Histotype = c("CCC" = viridis(20)[1],
                                  "EC" = viridis(20)[7],
                                  "HGSC" = viridis(20)[13],
                                  "MC" =  viridis(20)[19]),
                    Sur_Grp = c("x<3" = plasma(15)[2],
                                "6>x>3" = plasma(15)[7],
                                "x>6" = plasma(15)[12]))

################################################################################
# Define global paths
################################################################################
wd <- getwd()
dataPath <- paste(wd, "/Data/", sep="")
outPath <- paste(wd, "/Export/", sep="")
plotPath <- paste(outPath, "Plots/", sep="")
target_path <- paste(dataPath, "Train/trainPheno.csv", sep="")
dmrNewDir <- "C:/Users/xswehu/OneDrive - University of Gothenburg/xswehu_project/R_out/degDmrOut/"
degDir <- "C:/Users/xswehu/OneDrive - University of Gothenburg/xswehu_project/R_out/rnaOut"
countPath <- "C:/Users/xswehu/OneDrive - University of Gothenburg/xswehu_project/Data/RNAseq/FeatureCounts/sum_ftCounts.tab"
target_path <- "C:/Users/xswehu/OneDrive - University of Gothenburg/xswehu_project/Data/DNAmet/idat_dir/sampleData.csv"

# If Paths do not exist, create them
ifelse(!dir.exists(file.path(outPath)), dir.create(file.path(outPath)), FALSE)
ifelse(!dir.exists(file.path(plotPath)), dir.create(file.path(plotPath)), FALSE)

# Load saved workspace
#load(paste(wd,  'loadedData.RData',sep=""))

################################################################################
# Set up parallel processing (furure.apply for Win)
################################################################################
# Get half of the available cores for multiprocessing
nCores = length(availableWorkers(logical = FALSE))*0.5
options(future.globals.maxSize = 4000 * 1024^2)
plan(multisession, workers = nCores) 

################################################################################
################################################################################
# Source functions
################################################################################
################################################################################

source(paste(wd, "/Source/genScripts.R", sep=""))
source(paste(wd, "/Source/betaScripts.R", sep=""))
source(paste(wd, "/Source/pltScripts.R", sep=""))

options(scipen = 999)

################################################################################
################################################################################
# Setup
message("Build mapping data")
# Only needs to be ran once
################################################################################
################################################################################

################################################################################
# Get canonical hg38 ensembl transcripts of protein coding genes on chromosomes 1:22 and X
################################################################################
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl", 
                host = 'www.ensembl.org')
geneInf <- getBM(attributes = c("ensembl_gene_id",
                     "external_gene_name",
                     "entrezgene_id",
                     "chromosome_name",
                     "start_position",
                     "end_position",
                     "strand",
                     "transcript_is_canonical", 
                     "gene_biotype"), 
                 mart = mart)
geneInf <- subset(geneInf,
                  transcript_is_canonical == 1)
# pseudoNames <- names(table(geneInf$gene_biotype[grep("pseudo", geneInf$gene_biotype)]))
geneInf <- geneInf[geneInf$gene_biotype %in% c("protein_coding"), ]
geneInf <- geneInf[geneInf$chromosome_name %in% c(1:22, "X"), ]
geneInf <- geneInf[!duplicated(geneInf$ensembl_gene_id),]
geneInf <- makeGRangesCompatible(geneInf)
# Add promoter coordinates to gene dataframe
geneInf <- addPromoToGene_Mult(geneInf, nCores = 6)
write.csv(geneInf, paste(outPath, "hg38GeneInf.csv", sep=""))

################################################################################
# Get the hg37 gene coordinates for hg38 genes
################################################################################
#hg19GeneInf <- makeHg19GeneCoords(inpGenes = geneInf)
#hg19GeneInf <- makeGRangesCompatible(hg19GeneInf)
#hg19GeneInf <- addPromoToGene_Mult(hg19GeneInf, nCores=6)
#write.csv(hg19GeneInf, paste(outPath, "hg19GeneInf.csv", sep=""))
# hg19GeneInf <- read.csv(paste(outPath, "hg19GeneInf.csv"))

################################################################################
# Get full CpG annotation for the EPIC array
################################################################################
# The standard locations here are given in hg19, with the hg38 columns given as optional inputs not used by most packages
cpgAnno <- cpgAnno <- IlluminaHumanMethylationEPICanno.ilm10b5.hg38 %>% 
  getAnnotation %>% 
  as.data.frame
write.csv(cpgAnno, paste(outPath, "cpgAnno.csv", sep=""))

################################################################################
# Build positional dataframe for the EPIC array
################################################################################
hg19CpgLocsEPIC <- as.data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Locations)
hg19CpgLocsEPIC$start <- hg19CpgLocsEPIC$pos
hg19CpgLocsEPIC$end <- hg19CpgLocsEPIC$pos+2

################################################################################
# Build positional dataframe for the 450K array
################################################################################
hg19CpgLocs450K <- as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations)
hg19CpgLocs450K$start <- hg19CpgLocs450K$pos
hg19CpgLocs450K$end <- hg19CpgLocs450K$start + 2
write.csv(hg19CpgLocs450K, paste(outPath, "hg19CpgLocs450K.csv", sep=""))
# hg19CpgLocs450K <- read.csv(paste(outPath, "hg19CpgLocs450K.csv"))

################################################################################
# Create hg38 annotation (i.e. we remap values from hg19 to hg38)
# 450k and EPIC CpGs map to the same coordinates, as such the hg38 coordinates for EPIC are the hg38 coordinates for 450K 
################################################################################
hg38CpgLocsEPIC <- hg19CpgLocsEPIC
cpgAnno <- cpgAnno[match(rownames(hg38CpgLocsEPIC), rownames(cpgAnno)),]
hg38CpgLocsEPIC$pos <- cpgAnno$Start_hg38
hg38CpgLocsEPIC$start <- cpgAnno$Start_hg38
hg38CpgLocsEPIC$chr <- cpgAnno$CHR_hg38
hg38CpgLocsEPIC$end <- cpgAnno$End_hg38
hg38CpgLocsEPIC$strand <- cpgAnno$Strand_hg38

# Some probes lack hg38 mapping, we filter these away
hg38CpgLocsEPIC <- hg38CpgLocsEPIC[!is.na(hg38CpgLocsEPIC$pos),]

# Create hg39 based dataframe for translating CpGs
sharedCpgHg19Hg38 <- intersect(rownames(hg38CpgLocsEPIC), rownames(hg19CpgLocsEPIC))
transTableEpic <- hg38CpgLocsEPIC[sharedCpgHg19Hg38,]
transTableHg19 <- hg19CpgLocsEPIC[sharedCpgHg19Hg38,]
transTableHg19 <- transTableHg19[match(rownames(transTableEpic), rownames(transTableHg19)),]
transTableEpic$hg19_pos <- transTableHg19$pos
transTableEpic$hg19_start <- transTableHg19$start
transTableEpic$hg19_end <- transTableHg19$end
transTableEpic$hg19_chr <- transTableHg19$chr
transTableEpic$hg19_strand <- transTableHg19$strand

write.csv(hg38CpgLocsEPIC, paste(outPath, "hg38CpgLocsEPIC.csv", sep=""))
write.csv(hg19CpgLocsEPIC, paste(outPath, "hg19CpgLocsEPIC.csv", sep=""))
write.csv(transTableEpic, paste(outPath, "transTableEpic.csv", sep=""))
# hg19CpgLocsEPIC <- read.csv(paste(outPath, "hg19CpgLocsEPIC.csv", sep=""))
# hg38CpgLocsEPIC <- read.csv(paste(outPath, "hg19CpgLocsEPIC.csv", sep=""))
# transTableEpic <- read.csv(paste(outPath, "transTableEpic.csv", sep=""))

################################################################################
# Create hg38 annotation for 450K
################################################################################
sharedCpG <- intersect(rownames(hg19CpgLocs450K), rownames(hg38CpgLocsEPIC))
hg38CpgLocs450K <- hg38CpgLocsEPIC[sharedCpG, ]
hg38CpgLocs450K <- na.omit(hg38CpgLocs450K)
write.csv(hg38CpgLocs450K, paste(outPath, "hg38CpgLocs450K.csv", sep=""))

################################################################################
# Load phenotype data for train
################################################################################
trainPheno <- read.csv(paste(dataPath, "Train/TrainPheno.csv", sep=""),  row.names = 1)
rownames(trainPheno) <- trainPheno$Sample_ID
# Add the barcode parameter for matching in script (must be equal to that of column names in beta-matrix)
trainPheno$barcode <- paste(trainPheno$Sample_ID, trainPheno$barcode.ch1, sep="_")
write.csv(trainPheno ,paste(dataPath, "Train/TrainPheno.csv", sep=""))

################################################################################
# Load Beta-matrix, and M-matrix for train
################################################################################
trainBeta <- read.csv(paste(dataPath, "Train/Train_Beta.csv", sep=""),  row.names = 1)
# Make sure that all samples are accounted for
trainBeta <- trainBeta[,colnames(trainBeta) %in% paste(trainPheno$Sample_ID, trainPheno$barcode.ch1, sep="_")]
# Reorder the dataframe to match the order of the phenotype data if needed
trainBeta <- trainBeta[,match(colnames(trainBeta), paste(trainPheno$Sample_ID, trainPheno$barcode.ch1, sep="_"))]
# Convert beta values into m-values
trainM <- log2(trainBeta/(1-trainBeta))
# Create vector of all genes 
allGenes <- geneInf$ensembl_gene_id

################################################################################
################################################################################
# Setup
message("Begin loading of initial data")
################################################################################
################################################################################

################################################################################
# Load Cpg-site location mapping data
################################################################################
geneInf <- read.csv(paste(outPath, "hg38GeneInf.csv", sep=""), row.names = 1)
hg38CpgLocsEPIC <- read.csv(paste(outPath, "/hg38CpgLocsEPIC.csv",sep=""), row.names = 1)
hg38CpgLocs450K <- read.csv(paste(outPath, "hg38CpgLocs450K.csv", sep=""), row.names = 1)
promoter_cpgs <- read.csv(paste(outPath, "hg38GenePromoterCpgs.csv", sep=""))
cpgAnno <- read.csv(paste(outPath, "cpgAnno.csv", sep=""), row.names = 1)

################################################################################
# Load Training cohort data
################################################################################
trainPheno <- read.csv(paste(dataPath, "Train/Trainpheno.csv", sep=""), row.names = 1)
trainBeta <- read.csv(paste(dataPath, "Train/Train_Beta.csv",sep=""), row.names = 1)
trainM <- log2(trainBeta/(1-trainBeta))

################################################################################
################################################################################
################################################################################
# Quality control plots
message("Plotting PCA, MDS and Densities for training/test cohorts")
################################################################################
################################################################################
################################################################################

################################################################################
# Make PCA-plots
################################################################################
makeMethPcaPlot(trainBeta, trainPheno, nameBool = "Train", noProbes = 1000)

################################################################################
# Plot beta-densities for histotypes in cohorts
################################################################################
makeMethBetaFreqPlt(trainBeta, trainPheno, inpSid = "Train")

################################################################################
# Plot heatmaps for external beta lists and training cohort
################################################################################
makeTopHMap(trainBeta, trainPheno, noGenes = 1000, fileExt = "Train", pltBool = TRUE, 
            titleVal = "Top 1000 variable probes", rowAnnInpDf = NULL, nameVal = "Beta", heatBool = FALSE)

################################################################################
################################################################################
# Get the methylation type of CpG sites in dataframe (hypo, hemi, hyper)
# Warning!!! extremely long runtime, only run this if the information is required
message("Running methylation Beta, M, type analysis")
################################################################################
################################################################################

trainBetaMethType <- makeBetaDfMethType(trainBeta, trainPheno)
# write.csv(trainBetaMethType, paste(outPath, "trainBetaMethType.csv", sep="/"))
# Plot methylation types as barplot
makeBetaBar(trainBetaMethType, trainBeta)

################################################################################
################################################################################
# Get CpG sites overlapping the promoter region
################################################################################
################################################################################

################################################################################
# Create CpG genome coordinate overlaps for all genes of interest 
# Warning: Long runtime, preferably should only be ran once!
################################################################################
# For hg38 genes
hg38GenePromoterCpgDfLst <- makeCpgGeneOverlap_MULT(inpGenes = geneInf$ensembl_gene_id, 
                                     cpgInp = hg38CpgLocsEPIC, 
                                     inpGeneInf = geneInf, 
                                     type="PROMOTER")
#makeCsvSave(hg38GenePromoterCpgDfLst,
#            sType = "hg38GenePromoterCpgDf", 
#            dirExtBool = "GeneDf")

# For hg19 genes
#hg19AllPromoCpgPos <- makeCpgGeneOverlap_MULT(inpGenes = hg19GeneInf$ensembl_gene_id, 
#                                      cpgInp = hg19CpgLocsEPIC, 
#                                      inpGeneInf = hg19GeneInf, 
#                                      type="PROMOTER",
#                                      nCores = 6)

################################################################################
# Create CpG site promoter location dataframe for future reference in script
# Warning: Long runtime
################################################################################
hg38GenePromoterCpgs <- makeGeneCpGLocDataframe_MULT(inpLst = hg38GenePromoterCpgDfLst)
write.csv(hg38GenePromoterCpgs ,paste(outPath, "hg38GenePromoterCpgs.csv", sep=""))
hg38GenePromoterCpgs <- read.csv(paste(outPath, "hg38GenePromoterCpgs.csv", sep=""), row.names = 1)

#hg19Promoter_cpgs <- makeGeneCpGLocDataframe(hg19AllPromoCpgPos)
#write.csv(hg19Promoter_cpgs ,paste(outPath, "hg19GenePromoterCpgs.csv", sep=""))
# hg19Promoter_cpgs <- read.csv(paste(outPath, "hg19GenePromoterCpgs.csv", sep=""), row.names = 1)

################################################################################
# Subset beta dataframe using promoter CpGs (can be used to lower load-times)
################################################################################
allPromoCpgVec <- c()
for(i in 1:length(hg38GenePromoterCpgDfLst)){
  tmpCpgs <- hg38GenePromoterCpgDfLst[[i]]
  allPromoCpgVec <- append(allPromoCpgVec, rownames(tmpCpgs))
}

allPromoBetaVals <- trainBeta[which(rownames(trainBeta) %in% allPromoCpgVec),]
allPromoBetaVals  <- na.omit(allPromoBetaVals)
write.csv(allPromoBetaVals, paste(outPath, "allPromoBetaVals.csv",sep=""))
# allPromoBetaVals <- read.csv(paste(outPath, "allPromoBetaVals.csv",sep=""), row.names = 1)

allPromoMVals <- trainM[which(rownames(trainM) %in% allPromoCpgVec),]
allPromoMVals <- na.omit(allPromoMVals)
write.csv(allPromoMVals, paste(outPath, "allPromoMVals.csv",sep=""))

################################################################################
# Get mean methylation type of promoter region CpG's with 3 or more CpG's
################################################################################
promoBetaStats <- makePromoBetaMeanTypeStats(hg38GenePromoterCpgDfLst, 
                                         trainBeta, 
                                         trainPheno,
                                         nCpg = 3)
write.csv(promoBetaStats, paste(outPath, 
                                "allPromoBetaMeanTypeStats.csv",sep=""))

################################################################################
# Get the statistics for the methylation tpye of all samples, individually
################################################################################
trainPromoBetaMethType <- makeBetaDfMethType_MULT(allPromoBetaVals,
                                                  trainPheno, 
                                                  nCores = 6)

write.csv(trainPromoBetaMethType, paste(outPath, "trainPromoBetaMethType.csv"))
# trainPromoBetaMethType <- read.csv(paste(outPath, "trainPromoBetaMethType.csv", sep="/"))

################################################################################
# Plot methylation types of promoter CpGs as barplot
################################################################################
makeBetaBar(trainPromoBetaMethType, 
            trainBeta, 
            fileExt = "Promoter_Cpgs")

################################################################################
# Get beta values for promoter CpG positions
################################################################################
allPromoBetas <- getGeneCpgBeta_MULT(geneCpgInp = hg38GenePromoterCpgs, 
                                  inpBeta = allPromoBetaVals, 
                                  inpPheno = trainPheno,
                                  cpgInf = hg38CpgLocsEPIC,
                                  allBool = TRUE)

# promoBetaStats <- read.csv(paste(outPath, "promoBetaStats.csv", sep="/"))
################################################################################
################################################################################
# Get distribution type of CpG sites in dataframe
message("Running distribution type analysis for training cohort")
################################################################################
################################################################################

# Warning, long runtime as we are checking ALL relevant CpG sites!
trainCDist <- makeCpGDistTypeDf_MULT(trainBeta, 
                                     trainPheno, 
                                     nCores = 6)
trainMCDist <- makeCpGDistTypeDf_MULT(trainM, 
                                      trainPheno, 
                                      nCores = 6)

# Run the same script for the promoter CpG sites
promoDistTypes <-  makeCpGDistTypeDf_MULT(allPromoBetaVals, 
                                          trainPheno, 
                                          nCores = 6)

# Write results to wd
write.csv(trainCDist, 
          paste(outPath, "trainCpgDistribution.csv", sep="/"))
write.csv(trainMCDist, 
          paste(outPath, "trainCpgDistribution_M_Vals.csv", sep="/"))
write.csv(promoDistTypes, 
          paste(outPath, "trainPromoDistType.csv", sep="/"))

trainCDist <- read.csv(paste(outPath, "trainCpgDistribution.csv", sep="/"), sep=";")
trainMDist <- read.csv(paste(outPath, "trainCpgDistribution_M_Vals.csv", sep="/"))

################################################################################
# Plot distribution types
################################################################################
makeDistBar(trainCDist, 
            trainBeta, 
            fileExt = "CpG_Coverage_Beta")
makeDistBar(trainMCDist, 
            trainM, 
            fileExt = "CpG_Coverage_M")
makeDistBar(promoDistTypes, 
            allPromoBetaVals, 
            fileExt = "CpG_Coverage_Promo")

################################################################################
################################################################################
# Promoter coverage statistics
message("Running promoter coverage analysis")
################################################################################
################################################################################

################################################################################
# Get the no of promoter regions of Ensembl genes in EPIC and 450k arrays containting x CpGs, respectively
################################################################################
covEpic <- makePromoCovDf_MULT(inpCpgs = hg38GenePromoterCpgs,
                               inpDf = hg38CpgLocsEPIC)
cov450k <- makePromoCovDf_MULT(inpCpgs = hg38GenePromoterCpgs,
                               inpDf = hg38CpgLocs450K)

# Add both result vectors to dataframe
promoCovDf <- data.frame(matrix(nrow=2, ncol=11))
rownames(promoCovDf) <- c("EPIC", "450K")
colnames(promoCovDf) <- c("ALL", 0,1,2,3,4,5,"1-3","4-6", "7-10", "10+")
promoCovDf[1,] <- covEpic
promoCovDf[2,] <- cov450k
write.csv(promoCovDf, paste(outPath, "promoCovDf.csv"))

################################################################################
# Save coverage as barplot
################################################################################
promoCovPercDf <- makePromoCovBarPlt(promoCovDf)
write.csv(promoCovPercDf, paste(outPath, "promoCovDfPerc.csv"))

################################################################################
################################################################################
# Comparative analysis of arithmethric mean, median and trimean
message("Comparing mean, median and trimean in training cohort")
################################################################################
################################################################################

################################################################################
# Compare mean and trimean for all CpG sites for phenotypic groups
################################################################################
diffMeanTmDf <- makeDiffBetaDf(inpBeta = trainBeta, 
                               inpPheno = trainPheno, 
                               meanFunc = NULL, 
                               pltBool = TRUE)
write.csv(diffMeanTmDf, paste(outPath, "meanTMDifference.csv", sep="/"))
# Compare median and trimean for all CpG sites for phenotypic groups
diffMedianTmDf <- makeDiffBetaDf(inpBeta = trainBeta, 
                                 inpPheno = trainPheno, 
                                 meanFunc = "Median", 
                                 pltBool = TRUE)
write.csv(diffMedianTmDf, paste(outPath, "medianTMDifference.csv", sep="/"))
# Calculate variance between trimean and median/mean over all CpGs
tmVarDf <- makeTMVarDiff(trainBeta, 
                         trainPheno)
write.csv(tmVarDf, paste(outPath, "TriMean_VarDf_All_Betas.csv", sep=""))

################################################################################
# Repeat for promoter CpGs in EPIC array
################################################################################
# Compare mean and trimean for all CpG sites for phenotypic groups
diffMeanTmDf_Promo <- makeDiffBetaDf(inpBeta = allPromoBetaVals, 
                               inpPheno = trainPheno, 
                               meanFunc = NULL, 
                               pltBool = TRUE)
write.csv(diffMeanTmDf, paste(outPath, "meanTMDifference.csv", sep="/"))
# Compare median and trimean for all CpG sites for phenotypic groups
diffMedianTmDf_Promo <- makeDiffBetaDf(inpBeta = allPromoBetaVals, 
                                 inpPheno = trainPheno, 
                                 meanFunc = "Median", 
                                 pltBool = TRUE)

tmVarDf_Promo <- makeTMVarDiff(allPromoBetaVals, 
                               trainPheno)
write.csv(tmVarDf_Promo, paste(outPath, "TriMean_VarDf_Promo_Betas.csv", sep=""))

################################################################################
################################################################################
# Perform HSP analysis on Train
message("Running HSP analysis")
################################################################################
################################################################################

################################################################################
# Rank genes seen to dna-methylation separation for overlapping CpG-sites with at least 3 overlapping CpG sites within 2000bp upstream, 200bp downstreat of start
# Filter away genes where either the total distance is under 0.6 (i.e. x+z+y < 0.6) or where each histotypes median distance is under 0.2  
################################################################################

# promoB <- trainBeta[which(rownames(trainBeta) %in% allPromoCpgVec),]
# promoM <- trainM[which(rownames(trainM) %in% allPromoCpgVec),]
promoRankLst <- list() 
for(i in 1:length(names(table(trainPheno$Histotype)))){
  tmpH <- names(table(trainPheno$Histotype))[i]
  promoBetaRankDf <- makePromoBetaRank_V3(inpRegLocLst = hg38GenePromoterCpgDfLst,
                                          inpPheno = trainPheno, 
                                          inpGeneInf = geneInf,
                                          allBetas = trainBeta, 
                                          focGrp = tmpH,
                                          minCpg = 3, 
                                          sigCpg = 2, 
                                          noCat = 3,
                                          bCut = 0.2, 
                                          pCut = 0.05,
                                          brownBool = TRUE,
                                          cpgCutFreq = 0.25, 
                                          distType = "TRIMEAN")
  promoRankLst[[i]] <- promoBetaRankDf
  names(promoRankLst)[i] <- tmpH
}
names(promoRankLst) <- names(table(trainPheno$Histotype))

# Remove empty entries from list
promoRankLst <- promoRankLst[map(promoRankLst, function(x) dim(x)[1]) > 0]
makeCsvSave(promoRankLst, "promoters_ranked")
allSigPromo <- promoRankLst

################################################################################  
################################################################################
# Map map coordinates to HSPs
message("Mapping coordinate to HSPs")
################################################################################
################################################################################

# Summarize hsp categories
sigPromoDf <- data.frame(matrix(ncol=2, nrow=4))
colnames(sigPromoDf) = c("Histotype", "nSigPromo")
rownames(sigPromoDf) = names(table(trainPheno$Histotype)) 
for(i in 1:length(allSigPromo)){
  kfSum <- nrow(allSigPromo[[i]])
  sigPromoDf[which(rownames(sigPromoDf) %in% names(allSigPromo)[i]),] <- c(names(allSigPromo)[i], kfSum)
}
sigPromoDf[is.na(sigPromoDf)] <- 0
write.csv(sigPromoDf, paste(outPath, "sigPromoDf.csv",sep="/"))

################################################################################
# Make lists of cpg-sites for involved promoters/genes
################################################################################
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

# Create list of vectors of all CpGs of interest
sigPromoCpgLst <- list()
sigPromoCpgVec <- c()
for(i in 1:length(allSigPromoDfLst)){
  tmpG <- allSigPromoDfLst[[i]]
  tmpName <- names(allSigPromoDfLst)[i]
  tmpVec <- c()
  for(j in 1:length(tmpG)){
    tmpDf <- tmpG[[j]]
    tmpVec <- append(tmpVec, rownames(tmpDf))
  }
  sigPromoCpgVec <- append(sigPromoCpgVec,tmpVec)
  sigPromoCpgLst[[i]] <- tmpVec
  names(sigPromoCpgLst)[i] <- tmpName
}
makeCsvSave(sigPromoCpgLst, "sigPromoCpgVecs")
write.csv(sigPromoCpgVec,paste(outPath, "allSigPromoCpgVec.csv", sep=""))

sigPromoCpgFiles <- list.files(path = outPath, 
                               pattern = "sigPromoCpgVecs",
                               recursive = TRUE,
                               full.names = TRUE)
sigPromoCpgLst <- makeDfLst(sigPromoCpgFiles)
for(i in 1:length(sigPromoCpgLst)){
  sigPromoCpgLst[[i]] <- sigPromoCpgLst[[i]][[1]]
}
sigPromoCpgVec <- read.csv(paste(outPath, "allSigPromoCpgVec.csv", sep=""), row.names = 1)
sigPromoCpgVec <- sigPromoCpgVec[[1]]

################################################################################
# Create a list of dataframes containing the CpGs for the significant promoter regions
################################################################################
sigPromoBetaDfLst <- list()
for(i in 1:length(allSigPromoCpGSumDfLst)){
  tmpH <- names(allSigPromoCpGSumDfLst)[i]
  tmpCpg <- allSigPromoCpGSumDfLst[[i]]
  hBetaDf <- makeCpgDf(tmpCpg, 
                       trainBeta, 
                       trainPheno, 
                       tmpH)
  sigPromoBetaDfLst[[i]] <- hBetaDf
  names(sigPromoBetaDfLst)[i] <- tmpH
}

# Get cpg locations for best-probes
cpgLocsLst <- list()
for(i in 1:length(sigPromoCpgLst)){
  tmpCpgs <- sigPromoCpgLst[[i]]
  tmpH <- names(sigPromoCpgLst)[i]
  tmpLocs <- hg38CpgLocsEPIC[tmpCpgs,]
  cpgLocsLst[[i]] <- tmpLocs
  names(cpgLocsLst)[i] <- tmpH
}
makeCsvSave(cpgLocsLst, "sigPromoCpgDf")

################################################################################
# Get dataframes for each histotype reflective of the promoters passing threshold criterions
################################################################################
keepPromoCpgPos <- list() 
for(i in 1:length(allSigPromo)){
  tmpPromoGenes <- allSigPromo[[i]]
  tmpName <- names(allSigPromo)[i]
  tmpPromo <- hg38GenePromoterCpgs[which(hg38GenePromoterCpgs$Gene %in% tmpPromoGenes$ensembl_gene_id),]
  tmpInf <- geneInf[which(geneInf$ensembl_gene_id %in%  tmpPromoGenes$ensembl_gene_id),]
  keepPromoCpgPos[[i]] <- tmpPromo
  names(keepPromoCpgPos)[i] <- tmpName
}

makeCsvSave(keepPromoCpgPos, "sigPromoCpGLocDf")

################################################################################
# Create list of gene-inf dataframes
################################################################################
keepPromoInfLst <- list()
for(i in 1:length(allSigPromo)){
  tmpSep <- allSigPromo[[i]]
  tmpH <- names(allSigPromo)[i]
  tmpH38 <- geneInf[geneInf$ensembl_gene_id %in% tmpSep$ensembl_gene_id, ]
  keepPromoInfLst[[i]] <- tmpH38
  names(keepPromoInfLst)[i] <- tmpH
}
makeCsvSave(keepPromoInfLst, "sigPromoInfDf")

################################################################################
# Plot Trimean dotplots and boxplots for for significant promoter regions
################################################################################

for(i in 1:length(allSigPromo)){
  tmpS <- allSigPromo[[i]]
  tmpG <- geneInf[geneInf$ensembl_gene_id %in% tmpS$ensembl_gene_id,]
  tmpH <- names(allSigPromo)[i]
  geneDotPlots <- makeGeneLstCpgDotPlot(inpGenes = tmpG, 
                                        inpBeta = trainBeta, 
                                        inpPheno = trainPheno, 
                                        cpgInp = hg38CpgLocsEPIC, 
                                        geneInfInp = geneInf,
                                        fileExt=paste(tmpH,"SigPromo", sep="_"), 
                                        dirExt=paste(tmpH,"allSigPromoHg38_triMean", sep="_"),
                                        promoBool = TRUE)
  geneBoxPlots <- makeGeneLstCpgBoxPlot(inpGenes = tmpG, 
                                        inpBeta = trainBeta, 
                                        inpPheno = trainPheno, 
                                        cpgInp = hg38CpgLocsEPIC, 
                                        geneInfInp = geneInf, 
                                        fileExt=paste(tmpH,"SigPromo", sep="_"),  
                                        dirExt=paste(tmpH,"allSigPromoHg38_triMean", sep="_"), promoBool = TRUE)
}

################################################################################
################################################################################
# Get methylation type in percent of significant regions
message("Summarizing HSPs methylation data")
################################################################################
################################################################################

################################################################################
# Get methylation type stats for all CpGs in the promoter regions
################################################################################
sigPromoAllMethTypeLst <- list()
for(i in 1:length(sigPromoCpgLst)){
  tmpB <- trainBeta[sigPromoCpgLst[[i]], ]
  tmpMethType <- makeBetaDfMethType_MULT(tmpB,
                                         trainPheno, 
                                         nCores = 6)
  sigPromoAllMethTypeLst[[i]] <- tmpMethType
  names(sigPromoAllMethTypeLst)[i] <- names(sigPromoCpgLst)[i]
}
makeCsvSave(sigPromoAllMethTypeLst, "allCpgSigPromoStats")

################################################################################
# Get the percentage of samples of hypo, hemi, hypermethylation type in individual promoter regions
################################################################################
hspMethTypes <- list()
for(i in 1:length(allSigPromo)){
  tmpH <- names(allSigPromo)[i] 
  hspMethTypes[[i]] <- makeCpgMethTypePerc(inpGenes = allSigPromo[[i]], 
                                           inpCpgs = allSigPromoCpGSumDfLst[[i]], 
                                           inpPheno = trainPheno, 
                                           inpBeta = trainBeta,
                                           inpFoc = tmpH)
  names(hspMethTypes)[i] <- tmpH
}
makeCsvSave(hspMethTypes, "sigPromoMethTypesPerc")

################################################################################
# Get mean value methylation type for phenotypes in all promoter regions
################################################################################
promoBetaMeanTypeStatLst <- list()
for(i in 1:length(allSigPromoDfLst)){
  tmpStats <- makePromoBetaMeanTypeStats(inpPromoPos = allSigPromoDfLst[[i]], 
                                     inpBeta = trainBeta, 
                                     inpPheno = trainPheno)
  promoBetaMeanTypeStatLst[[i]] <- tmpStats
  names(promoBetaMeanTypeStatLst)[i] <- names(allSigPromoDfLst)[i]
}
makeCsvSave(promoBetaMeanTypeStatLst, "SigPromoBetaTypeStats")

################################################################################
# Calculate distribution types of CpG's in the promoter region
################################################################################
sigPromoDistTypeLst <- list()
for(i in 1:length(sigPromoBetaDfLst)){
  tmpDf <- sigPromoBetaDfLst[[i]]
  tmpH <- names(sigPromoBetaDfLst)[i]
  tmpDist <- makeCpGDistTypeDf_MULT(tmpDf,
                                        trainPheno,
                                        nCores = 6)
  sigPromoDistTypeLst[[i]] <- tmpDist
  names(sigPromoDistTypeLst)[i] <- tmpH
  makeDistBar(tmpDist, 
              tmpDf, 
              fileExt = paste(tmpH, "_SigPromo_DistType",sep=""))
}
makeCsvSave(sigPromoDistTypeLst, "sigPromoDistTypes")
#sigPromoDistStatFiles <- list.files(path = outPath, 
#                               pattern = "sigPromoDistTypes.csv",
#                               recursive = TRUE,
#                               full.names = TRUE)
#sigPromoDistTypeLst <- makeDfLst(sigPromoDistStatFiles)

# Deprecated
# train_cpg_diff <- makeBetaTMDiff(inpCpg = allSigPromoCpGSumDfLst, 
#                                  inpPheno = trainPheno, 
#                                  inpBeta = trainBeta, 
#                                  inpGeneInf = hg19GeneInf)
# gse226823_cpg_diff <- makeBetaTMDiff(inpCpg = allSigPromoCpGSumDfLst, inpPheno = gse226823_pheno, inpBeta = gse226823_beta, inpGeneInf = hg19GeneInf)
# gse51820_cpg_diff <- makeBetaTMDiff(inpCpg = allSigPromoCpGSumDfLst, inpPheno = gse51820_pheno, inpBeta = gse51820_beta, inpGeneInf = hg19GeneInf)
# makeCsvSave(list("TRAIN"= train_cpg_diff, "GSE226823"= gse226823_cpg_diff, "GSE51820"= gse51820_cpg_diff), "Diff_TM_Mean_HSP")

################################################################################
################################################################################
# Calculate TM differences across cohort, also include students t-test
message("Comparing Trimean vs. Mean for HSPs")
################################################################################
################################################################################

################################################################################
#  Get variance of mean vs. trimean in HSPs
################################################################################

tm_mean_var_stat <- makeBetaTM_Mean_Var(inpCpg=allSigPromoCpGSumDfLst, 
                                        inpPheno=trainPheno, 
                                        inpBeta = trainBeta, 
                                        inpGeneInf = geneInf)
makeCsvSave(tm_mean_var_stat, "sig_PromoMeanVar")

# Calculate for histotypes serparately
sigPromoVarDfLst <- list()
for(i in 1:length(sigPromoBetaDfLst)){
  tmpH <- names(sigPromoBetaDfLst)[i]
  tmpDf <- sigPromoBetaDfLst[[i]]
  tmpVar <- makeTMVarDiff(inpBeta = tmpDf, 
                          inpPheno =  trainPheno, 
                          fileExt = paste(tmpH,"_SigPromoVarDf", sep=""),
                          pltBool = TRUE)
  sigPromoVarDfLst[[i]] <- tmpVar
  names(sigPromoVarDfLst)[i] <- tmpH
}

################################################################################
# Compare Trimean, mean and median to one another for the same CpG sites
# In group comparisons
################################################################################
allSigBetaDf <- trainBeta[which(rownames(trainBeta) %in% sigPromoCpgVec), ]
allSigBetaDf <- na.omit(allSigBetaDf)
allHSPVarDf <- makeTMVarDiff(inpBeta = allSigBetaDf, 
                             inpPheno = trainPheno, 
                             fileExt = paste("All_SigPromoTMVarDiff", sep=""),
                             pltBool = TRUE)
write.csv(allHSPVarDf, paste(outPath, "All_SigPromoTMVarDiff.csv", sep="/"))

################################################################################
################################################################################
# Compare and plot mean and trimean to one another
# For significant promoter region CpGs
################################################################################
################################################################################
# allSigCpgPromo <- list()
# for(i in 1:length(allSigPromo)){
#   tmpSigs <- allSigPromo[[i]]
#   tmpH <- names(allSigPromo)[i]
#   sigPromoVec <- c()
#   for(j in 1:nrow(tmpSigs)){
#     tmpEns <- tmpSigs$ensembl_gene_id[j]
#     tmpReg <- hg38GenePromoterCpgs[hg38GenePromoterCpgs$Gene %in% tmpEns,]
#     tmpDf <-  makeCpgDf(data.frame(tmpReg),
#                         trainBeta, 
#                         trainPheno, 
#                         tmpH)
#     tmpCpg <- makeRegSigCpg(inpReg = tmpDf,
#                             inpBeta = trainBeta,
#                             inpPheno = trainPheno,
#                             inpH = tmpH,
#                             minCpg = 2)
#     if(length(tmpCpg) > 0){
#       sigPromoVec <- append(sigPromoVec, tmpCpg)
#     }
#   }
#   allSigCpgPromo[[i]] <- sigPromoVec
#   names(allSigCpgPromo)[i] <- tmpH
# }

allCpgPromo <- list()
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
    sigPromoVec <- append(sigPromoVec, rownames(tmpDf))
  }
  allCpgPromo[[i]] <- sigPromoVec
  names(allCpgPromo)[i] <- tmpH
}

################################################################################
# Calculate differences in mean vs. tukeys trimean
################################################################################
diffMeanLst <- list()
for(i in 1:length(allCpgPromo)){
  tmpPromo <- allCpgPromo[[i]]
  tmpH <- names(allCpgPromo)[i]
  tmpB  <- trainBeta[which(rownames(trainM) %in% tmpPromo), ]
  sigPromoDiffMean <- makeMeanTMDiffBetaDf(inpMat = tmpB, 
                                             inpPheno = trainPheno, 
                                             inpH = tmpH,
                                             limRows = NULL,
                                             pltBool = TRUE,
                                             fileExt = "sigPromo_Diff_Mean_TM_Beta")
  diffMeanLst[[i]] <- sigPromoDiffMean[[1]]
  names(diffMeanLst)[i] <- tmpH
}

################################################################################
# Calculate differences in T-test vs. Brunner-Munnz
################################################################################
diffPLst_B <- list()
for(i in 1:length(allCpgPromo)){
  tmpPromo <- allCpgPromo[[i]]
  tmpH <- names(allCpgPromo)[i]
  tmpB  <- trainBeta[which(rownames(trainM) %in% tmpPromo), ]
  sigPromoDiffTtBm <- makePairPvalDiffBetaDf(inpMat = tmpB, 
                                             inpPheno = trainPheno, 
                                             inpH = tmpH,
                                             limRows = NULL,
                                             pltBool = TRUE,
                                             fileExt = "sigPromo_Diff_TT_BM_Beta")
  diffPLst_B[[i]] <- sigPromoDiffTtBm[[1]]
  names(diffPLst_B)[i] <- tmpH
}

diffMultLst_B <- list()
for(i in 1:length(allCpgPromo)){
  tmpPromo <- allCpgPromo[[i]]
  tmpH <- names(allCpgPromo)[i]
  tmpB  <- trainBeta[which(rownames(trainM) %in% tmpPromo), ]
  sigPromoDiffTtBm <- makeMultPvalDiffBetaDf(inpMat = tmpB, 
                                             inpPheno = trainPheno, 
                                             inpH = tmpH,
                                             limRows = NULL,
                                             pltBool = TRUE,
                                             fileExt = "sigPromo_MultDiff_TT_BM_Beta")
  diffMultLst_B[[i]] <- sigPromoDiffTtBm[[1]]
  names(diffMultLst_B)[i] <- tmpH
}

################################################################################
# Repeat for M-values
################################################################################
diffPLst_M <- list()
for(i in 1:length(allCpgPromo)){
  tmpPromo <- allCpgPromo[[i]]
  tmpH <- names(allCpgPromo)[i]
  tmpM  <- trainM[which(rownames(trainM) %in% tmpPromo), ]
  sigPromoDiffTtBm <- makePairPvalDiffBetaDf(inpMat = tmpM, 
                                             inpPheno = trainPheno, 
                                             inpH = tmpH,
                                             limRows = NULL,
                                             pltBool = TRUE,
                                             fileExt = "sigPromo_PairDiff_TT_BM_M")
  diffPLst_M[[i]] <- sigPromoDiffTtBm[[1]]
  names(diffPLst_M)[i] <- tmpH
}

diffMLst_M <- list()
for(i in 1:length(allCpgPromo)){
  tmpPromo <- allCpgPromo[[i]]
  tmpH <- names(allCpgPromo)[i]
  tmpM  <- trainM[which(rownames(trainM) %in% tmpPromo), ]
  sigPromoDiffTtBm <- makeMultPvalDiffBetaDf(inpMat = tmpM, 
                                             inpPheno = trainPheno, 
                                             inpH = tmpH,
                                             limRows = NULL,
                                             pltBool = TRUE,
                                             fileExt = "sigPromo_MultDiff_TT_BM_M")
  diffMLst_M[[i]] <- sigPromoDiffTtBm[[1]]
  names(diffMLst_M)[i] <- tmpH
}

################################################################################
# Calculate differences in F-statistic vs. Dunn
################################################################################
sigPromoDiffFsDunn <- makepMultPvalDiffBetaDf(inpBeta = allSigBetaDf, 
                                         inpPheno = trainPheno, 
                                         inpM = trainM, 
                                         sigBool = NULL,
                                         limRows = NULL,
                                         fileExt = "sigPromo_Diff_TT_BM_M")

################################################################################
# Repeat for the full cohort
################################################################################
diffMeanTM_ALL <- makeMeanTMDiffBetaDf(inpBeta = allSigBetaDf, 
                                            inpPheno = trainPheno, 
                                            inpM = NULL, 
                                            sigBool = NULL,
                                            limRows = NULL,
                                            fileExt = "Mean_TM_Diff_Beta_sigPromo_CpGs")

diffTTDunnDf_ALL <- makepValDiffBetaDf_MULT(inpBeta = trainBeta, 
                                       inpPheno = trainPheno, 
                                       inpM = trainM, 
                                       sigBool = TRUE,
                                       limRows = NULL,
                                       fileExt = "allCpG_Beta",
                                       pltBool = TRUE)

################################################################################
# Plot differences between Dunns and t-test
################################################################################

#makeTTestDunnCompPlt(allSigBetaDf, 
#                trainPheno, 
#                fileAdd = "AllSigPromo_Beta", 
#                mBool = NULL)

makeTTestDunnCompPlt(allSigBetaDf, 
                trainPheno, 
                fileAdd = "AllSigPromo_M", 
                mBool = TRUE)

################################################################################
################################################################################
# Get the oncogenic potential for HSPs
message("Running Oncoscore for HSP promoters")
################################################################################
################################################################################

################################################################################
# CCC SigPromo
################################################################################
# Get HGNC symbols for the significant promoter regions
cccSigSyms <- geneInf$external_gene_name[geneInf$ensembl_gene_id %in% allSigPromo$CCC$ensembl_gene_id]
cccSigSyms <- cccSigSyms[!cccSigSyms %in% c("", "BAD", "GPT")]
cccSigSyms <- append(cccSigSyms, "ALT1")
cccOsQueryAll <-  OncoScore::perform.query(cccSigSyms)
# Calculate OncoScore, keep only significant hits 
cccOsRes <- OncoScore::compute.oncoscore(cccOsQueryAll)
cccOsDf <- data.frame(cccOsRes[order(cccOsRes[,"OncoScore"], decreasing = TRUE),])
sigCccOsDf <- cccOsDf[cccOsDf$OncoScore >= 21.09, ]
write.csv(sigCccOsDf,paste(outPath,"CCC/CCC_SigPromo_OncoScoreResults.csv", sep=""), row.names = TRUE)

################################################################################
# HGSC SigPromo
################################################################################
hgscSigSyms <- geneInf$external_gene_name[geneInf$ensembl_gene_id %in% allSigPromo$HGSC$ensembl_gene_id]
hgscSigSyms <- hgscSigSyms[!hgscSigSyms %in% ""]
hgscOsQueryAll <-  OncoScore::perform.query(hgscSigSyms)
hgscOsRes <- OncoScore::compute.oncoscore(hgscOsQueryAll)
hgscOsDf <- data.frame(hgscOsRes[order(hgscOsRes[,"OncoScore"], decreasing = TRUE),])
sigHgscOsDf <- hgscOsDf[cccOsDf$OncoScore >= 21.09, ]
write.csv(sigHgscOsDf,paste(outPath,"/HGSC/HGSC_SigPromo_OncoScoreResults.csv", sep=""), row.names = TRUE)

################################################################################
# MC SigPromo
################################################################################
mcSigSyms <- geneInf$external_gene_name[geneInf$ensembl_gene_id %in% allSigPromo$MC$ensembl_gene_id]
# Remove genes which can not be used in package, add synonyms
mcSigSyms <- mcSigSyms[!mcSigSyms %in% c("", "CD40")]
mcSigSyms <- append(mcSigSyms, "TNFRSF5")
# Perform OS query for all significant promoter regions
mcOsQueryAll <-  OncoScore::perform.query(mcSigSyms)
mcOsRes <- OncoScore::compute.oncoscore(mcOsQueryAll)
mcOsDf <- data.frame(mcOsRes[order(mcOsRes[,"OncoScore"], decreasing = TRUE),])
sigMcOsDf <- mcOsDf[mcOsDf$OncoScore >= 21.09, ]
write.csv(sigMcOsDf,paste(outPath,"/MC/MC_SigPromo_OncoScoreResults.csv", sep=""), row.names = TRUE)

################################################################################
################################################################################
# Load test-datasets
message("Loading test cohort datasets")
################################################################################
################################################################################

################################################################################
# Read in phenotype data
################################################################################
gse51820_pheno <- read.csv(paste(dataPath, "GSE51820/gse51820_pheno.csv",sep=""))
gse226823_pheno <- read.csv(paste(dataPath, "GSE226823/gse226823_pheno.csv", sep=""))

################################################################################
# Read in beta-matrices
################################################################################
gse51820_beta <- read.csv(paste(dataPath, "GSE51820/gse51820_Beta.csv",sep=""), row.names = 1)
gse226823_beta <- read.csv(paste(dataPath, "GSE226823/gse226823_Beta.csv",sep=""), row.names = 1)

################################################################################
# Intersect samples, keep only those in both pheno and beta-matrix
################################################################################
gse51820_filt <- makePhenoBetaFilter(gse51820_beta, gse51820_pheno)
gse51820_beta <- gse51820_filt[[1]]
gse51820_pheno <- gse51820_filt[[2]]

gse226823_filt <- makePhenoBetaFilter(gse226823_beta, gse226823_pheno)
gse226823_beta <- gse226823_filt[[1]]
gse226823_pheno <- gse226823_filt[[2]]

extPhenoLst <- list("GSE51820" = gse51820_pheno,
                    "GSE226823" = gse226823_pheno)
extBetaLst <- list("GSE51820" = gse51820_beta,
                   "GSE226823" = gse226823_beta)

################################################################################
# Plot PCA
################################################################################
for(i in 1:length(extPhenoLst)){
  tmpB <- extBetaLst[[i]]
  tmpName <- names(extPhenoLst )[i]
  tmpPheno <- extPhenoLst[[tmpName]]
  makeMethPcaPlot(tmpB, tmpPheno, nameBool = tmpName, noProbes = 1000)
}

################################################################################
# Plot Beta densities of test-cohorts
################################################################################
for(i in 1:length(extBetaLst)){
  tmpB <- extBetaLst[[i]]
  tmpName <- names(extBetaLst)[i]
  tmpPheno <- extPhenoLst[[tmpName]]
  makeMethBetaFreqPlt(tmpB, tmpPheno, inpSid = tmpName)
}

################################################################################
# Plot Heatmaps of test-cohorts
################################################################################
for(i in 1:length(extBetaLst)){
  tmpB <- extBetaLst[[i]]
  tmpSid <- names(extBetaLst)[i]
  tmpP <- extPhenoLst[[tmpSid]]
  makeTopHMap(tmpB, tmpP, noGenes = 1000, fileExt = paste(tmpSid), pltBool = TRUE, 
              titleVal = "Top 1000 variable probes", rowAnnInpDf = NULL, nameVal = "Beta", heatBool = FALSE)
}

################################################################################
# Calculate CpG distribution distribution of CpG sites in test cohorts
# With respect to each histotype
################################################################################
cpgTypeDfExtLst <- list()
cpgTypeDfExtLstM <- list()
for(i in 1:length(extBetaLst)){
  tmpB <- extBetaLst[[i]]
  tmpM <- log2(tmpB/(1-tmpB))
  tmpP <- extPhenoLst[[names(extBetaLst)[i]]]
  tmpDistType <- makeCpGDistTypeDf_MULT(tmpB, 
                                     tmpP, 
                                     nCores = 6)
  tmpDistTypeM <- makeCpGDistTypeDf_MULT(tmpM, 
                                      tmpP, 
                                      nCores = 6)
  cpgTypeDfExtLst[[i]] <- tmpDistType
  names(cpgTypeDfExtLst)[i] <- names(extBetaLst)[i]
  cpgTypeDfExtLstM[[i]] <- tmpDistTypeM
  names(cpgTypeDfExtLstM)[i] <- names(extBetaLst)[i]
}

makeCsvSave(cpgTypeDfExtLst, "cpgDistType_Beta_ExtCohort")
makeCsvSave(cpgTypeDfExtLst, "cpgDistType_M_ExtCohort")

################################################################################
# Plot the cpg-site distribution type of cohorts next to one another
# With respect to each histotype
################################################################################
compBetaLst <- extBetaLst
compBetaLst[[3]] <- trainBeta
names(compBetaLst)[3] <- "Train"

trainCDist <- read.csv(paste(outPath, "trainCpgDistribution.csv", sep=""), sep=",", row.names = 1)
trainMDist <- read.csv(paste(outPath, "trainCpgDistribution_M_Vals.csv", sep=""), sep=",", row.names = 1)
distTypeFiles_Beta <- list.files(path = outPath, 
                               pattern = "cpgDistType_Beta",
                               recursive = TRUE,
                               full.names = TRUE)
cpgTypeDfExtLst_Beta <- makeDfLst(distTypeFiles_Beta)
cpgTypeDfExtLst_Beta[[3]] <- trainCDist
names(cpgTypeDfExtLst_Beta)[3] <- "Train"

distTypeFiles_M <- list.files(path = outPath, 
                                 pattern = "cpgDistType_M",
                                 recursive = TRUE,
                                 full.names = TRUE)
cpgTypeDfExtLst_M <- makeDfLst(distTypeFiles_M)
cpgTypeDfExtLst_M[[3]] <- trainMDist
names(cpgTypeDfExtLst_M)[3] <- "Train"

makeDistBarLst(inpDistLst = cpgTypeDfExtLst_Beta, 
               inpMLst = cpgTypeDfExtLst_M, 
               inpBetaLst = compBetaLst)

extBTypeLst <- list()
for(i in 1:length(extBetaLst)){
  tmpB <- extBetaLst[[i]]
  tmpP <- extPhenoLst[[names(extBetaLst)[i]]]
  bTypeTmp <- makeBetaDfMethType_MULT(tmpB, 
                                      tmpP,
                                      nCores = 6)
  # bTypeTmp <- makeBetaDfMethType(tmpB, tmpP)
  extBTypeLst[[i]] <- bTypeTmp
  names(extBTypeLst)[i] <- names(extBetaLst)[i]
}
makeCsvSave(extBTypeLst, "BetaType_Ext")

cpgTypePercDfLst_B <- list()
cpgTypePercDfLst_M <- list()
for(i in 1:length(cpgTypeDfExtLst_Beta)){
  tmpSid <- names(cpgTypeDfExtLst_Beta)[i]
  tmpB <- compBetaLst[[tmpSid]]
  tmpD <- cpgTypeDfExtLst_Beta[[i]]
  tmpM <- cpgTypeDfExtLst_M[[i]]
  percDf_B <- round(tmpD/nrow(tmpB),3)
  percDf_M <- round(tmpM/nrow(tmpB),3)
  cpgTypePercDfLst_B[[i]] <- percDf_B
  names(cpgTypePercDfLst_B)[i] <- tmpSid
  cpgTypePercDfLst_M[[i]] <- percDf_M
  names(cpgTypePercDfLst_M)[i] <- tmpSid
}

makeCsvSave(cpgTypePercDfLst_B, "cpgTypePercDf_B")
makeCsvSave(cpgTypePercDfLst_M, "cpgTypePercDf_M")

################################################################################
# Promoter methylation of test-cohorts
message("Analysing signifcant promoter regions in test-cohorts")
################################################################################
promoBetaMeanTypeStatLst_51820 <- list()
promoBetaMeanTypeStatLst_226823 <- list()
for(i in 1:length(allSigPromoDfLst)){
  tmpH <- names(allSigPromoDfLst)[i]
  tmpStats_51820 <- makePromoBetaMeanTypeStats(inpPromoPos = allSigPromoDfLst[[i]], 
                                         inpBeta = extBetaLst$GSE51820, 
                                         inpPheno = extPhenoLst$GSE51820)
  promoBetaMeanTypeStatLst_51820[[i]] <- tmpStats_51820
  names(promoBetaMeanTypeStatLst_51820)[i] <- tmpH
  if(!tmpH %in% "MC"){
    tmpStats_226823 <- makePromoBetaMeanTypeStats(inpPromoPos = allSigPromoDfLst[[i]], 
                                               inpBeta = extBetaLst$GSE226823, 
                                               inpPheno = extPhenoLst$GSE226823)
    promoBetaMeanTypeStatLst_226823[[i]] <- tmpStats_226823
    names(promoBetaMeanTypeStatLst_226823)[i] <- tmpH
  }
}
makeCsvSave(promoBetaMeanTypeStatLst_51820, "SigPromoMeanStatsExt")
makeCsvSave(promoBetaMeanTypeStatLst_226823, "SigPromoMeanStatsExt")

################################################################################
# Plot trimean dotplots for genes of interest in external cohorts
################################################################################
extDotPltLst <- list()
for(i in 1:length(allSigPromo)){
  hPltLst <- list()
  tmpS <- allSigPromo[[i]]
  tmpG <- geneInf[match(tmpS$ensembl_gene_id, 
                        geneInf$ensembl_gene_id),]
  tmpH <- names(allSigPromo)[i]
  tmpLocs <- hg38CpgLocsEPIC
  # Plot the gene plots for each entry in the test cohort
  for(j in 1:length(extBetaLst)){
    tmpB <- extBetaLst[[j]]
    tmpSid <- names(extBetaLst)[j]
    tmpPheno <- extPhenoLst[[tmpSid]]
    message(paste("Creating dot plots for cohort: ", tmpSid, ", Histotype: ", tmpH, sep=""))
    pltLst <- makeGeneLstCpgDotPlot(inpGenes = tmpG, 
                                    inpBeta = tmpB, 
                                    inpPheno = tmpPheno, 
                                    cpgInp = tmpLocs, 
                                    geneInfInp = geneInf,
                                    fileExt=paste(tmpH, tmpSid, sep="_"),
                                    dirExt=paste(tmpSid, "_TestCohort", tmpH, sep=""),
                                    promoBool = TRUE)
    hPltLst[[j]] <- pltLst
    names(hPltLst)[j] <- tmpSid
    #message(paste("Creating box plots for cohort: ", tmpSid, ", Histotype: ", tmpH, sep=""))
    #makeGeneLstCpgBoxPlot(tmpG, tmpB, tmpPheno, hg38CpgLocs, geneInf, fileExt=paste(tmpH, tmpSid, sep="_"), dirExt=paste(tmpSid, "_", tmpH, sep=""))
  }
  extDotPltLst[[i]] <- hPltLst
  names(extDotPltLst)[i] <- tmpH
}

################################################################################
# Create "merged" ggarrange plot of genes in different cohorts
# i.e. comparative plot of the same gene
################################################################################
for(i in 1:length(extDotPltLst)){
  tmpHLst <- extDotPltLst[[i]]
  tmpH <- names(extDotPltLst)[i]
  tmpGenes <- allSigPromo[[tmpH]]
  for(j in 1:nrow(tmpGenes)){
    extLst <- list()
    tmpG <- tmpGenes[j,]
    for(k in 1:length(tmpHLst)){
      tmpExt <- tmpHLst[[k]]
      if(!tmpG$external_gene_name %in% names(tmpExt)){
        extLst[[k]] <- ""
        names(extLst)[[k]] <- ""
      }else{
        tmpPlt <- tmpExt[[tmpG$external_gene_name]]
        extLst[[k]] <- tmpPlt
        names(extLst)[k] <- names(tmpHLst)[k]
      }
    }
    extLst <- extLst[!names(extLst) %in% ""]
    if(length(extLst) < 2){
      message("Few or no matches for gene: ", tmpG$external_gene_name, " Gene will be skipped for comparative plot")
      next()
    }else if(length(extLst) == 2){
      pltCols <- 2
      pltRows <- 1
    }else{
      pltCols <- ceiling(length(extLst)/2)
      pltRows <- pltCols
    }
    gridPlt <- ggarrange(plotlist=extLst, widths = c(1,1), ncol=pltCols, nrow=pltRows, labels = c(names(extLst)))
    outDir <- paste(plotPath, "cpgDotPlot/SigPromo/", tmpH, sep="")
    ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir)), FALSE)
    outFile <- paste(outDir, "/", tmpG$external_gene_name, "_SigPromo_extCompPlot.pdf",sep="")
    ggsave(outFile, plot=gridPlt, width=60, height=40, units = "cm")
  }
}

save.image(paste(wd, "/main.RData", sep=""))
#load.image(paste(wd, "/main.RData", sep=""))
