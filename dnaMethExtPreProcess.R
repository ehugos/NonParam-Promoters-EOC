################################################################################
################################################################################
################################################################################
# Script for preprocessing datasets used in study
################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
# Load packages
################################################################################
################################################################################
# Check bioconductor packages, install if missing
bioc_packs <- c("ChAMP", "GEOquery", "minfi", "minfiData")
miss_bc_packs <- bioc_packs[!(bioc_packs %in% installed.packages()[,"Package"])]
if(length(miss_bc_packs)){
  BiocManager::install(miss_bc_packs)
} 

invisible(lapply(r_packs, library, character.only = TRUE))
invisible(lapply(bioc_packs, library, character.only = TRUE))
# Install maxprobes (not on biocondcuctor)
# remotes::install_github("markgene/maxprobes")
library(maxprobes)

################################################################################
################################################################################
# Define global parameters
################################################################################
################################################################################
histotypes <- c("CCC", "EC","HGSC", "MC")
catColVec <- viridis(length(histotypes))
names(catColVec) <- histotypes

# Check r-packages, install if missing
r_packs <- c("dplyr", "remotes", "tidyverse", "viridis")
miss_r_packs <- r_packs[!(r_packs %in% installed.packages()[,"Package"])]
if(length(miss_r_packs)){
  install.packages(miss_r_packs)
} 
invisible(lapply(r_packs, library, character.only = TRUE))

# Download options for GEO (win-based)
options(download.file.method.GEOquery = 'curl')
options(timeout = max(999999, getOption("timeout")))

# Datasets used in study
#GSE51820
#GSE133556
#GSE155760 
#GSE185008 
#GSE211686
#GSE226823 
#GSE263434 
#GSE267068
#GSE326000

################################################################################
################################################################################
# Define paths
################################################################################
################################################################################

wd <- getwd()
dataPath <- paste(wd, "Data/", sep="/") 
source(paste(wd, "Source/methPreProcessScripts.R", sep="/"))

GSE51820_path <- paste(dataPath, "GSE51820/GSE51820_series_matrix.txt.gz", sep="")
GSE133556_path <- paste(dataPath, "GSE133556/GSE133556_series_matrix.txt.gz", sep="")
GSE155760_path <- paste(dataPath, "GSE155760/GSE155760_series_matrix.txt.gz", sep="")
GSE185008_1_path <- paste(dataPath, "GSE185008/GSE185008-GPL13534_series_matrix.txt.gz", sep="")
GSE185008_2_path <- paste(dataPath, "GSE185008/GSE185008-GPL21145_series_matrix.txt.gz", sep="")
GSE211686_1_path <- paste(dataPath, "GSE211686/GSE211686-GPL13534_series_matrix.txt.gz", sep="")
GSE211686_2_path <- paste(dataPath, "GSE211686/GSE211686-GPL21145_series_matrix.txt.gz", sep="")
GSE226823_path <- paste(dataPath, "GSE226823/GSE226823_series_matrix.txt.gz", sep="")
GSE226872_1_path <- paste(dataPath, "GSE226872/GSE226872-GPL13534_series_matrix.txt.gz", sep="")
GSE226872_2_path <- paste(dataPath, "GSE226872/GSE226872-GPL16791_series_matrix.txt.gz", sep="")
GSE263434_path <- paste(dataPath, "GSE263434/GSE263434_series_matrix.txt.gz", sep="")
GSE267068_path <- paste(dataPath, "GSE267068/GSE267068_series_matrix.txt.gz", sep="")
trainPath <- paste(dataPath, "Train/GSE326000_series_matrix.txt.gz", sep="")

################################################################################
################################################################################
# Cohort preparation/Preprocessing
################################################################################
################################################################################

################################################################################
# GSE51820
# 13 CCC, 11 EC, 8 MC, 53 Ser
################################################################################
gse51820 <- getGEO(filename=GSE51820_path, getGPL = FALSE)
gse51820_supp <- getGEOSuppFiles("GSE51820", makeDirectory = TRUE, baseDir = paste(dataPath, "GSE51820/", sep=""))
gse51820_pheno <- gse51820@phenoData@data

# 51820 in .bpm format instead of .idat, use non-normalized values instead

# Read in genomicmethylset
gse51820_beta_raw <- readGEORawFile(paste(dataPath, "GSE51820/GSE51820_non_normalized.txt.gz", sep=""), 
                                    sep = "\t", 
                                    Uname = "Unmethylated Signal",
                                    Mname = "Methylated signal", 
                                    row.names = 1, 
                                    #pData =gse51820_pheno,
                                    array = "IlluminaHumanMethylation450k",
                                    # annotation = .default.450k.annotation,
                                    mergeManifest = FALSE,
                                    showProgress = TRUE)
gse51820_beta <- getBeta(gse51820_beta_raw, type="Illumina")
gse51820_beta <- gse51820_beta[,colnames(gse51820_beta) %in% gse51820_pheno$title]

# gse51820_pheno
gse51820_pheno <- tibble::rownames_to_column(gse51820_pheno, "Sample_ID")
gse51820_pheno$Histotype <- gse51820_pheno$characteristics_ch1.1
# Rename histotypes to work with scripts 
gse51820_pheno$Histotype <- plyr::mapvalues(
  gse51820_pheno$Histotype,
  c("histologic type: clear cell", "histologic type: endometrioid", "histologic type: serous", "histologic type: mucinous"),
  c('CCC', 'EC', 'HGSC', 'MC'))
colnames(gse51820_beta) <- gse51820_pheno$Sample_ID[match(colnames(gse51820_beta), gse51820_pheno$title)]

gse51820_pheno <- gse51820_pheno[gse51820_pheno$Histotype %in% histotypes, ]
gse51820_pheno$barcode <- gse51820_pheno$Sample_ID
gse51820_beta <- gse51820_beta[,colnames(gse51820_beta) %in% gse51820_pheno$barcode]
gse51820_pheno$Sample_Name <- gse51820_pheno$Sample_ID

# Filter using dmrcate
gse51820_beta <- DMRcate::rmSNPandCH(gse51820_beta)

# Retrieve cross-reactive probes 
xloci <- maxprobes::xreactive_probes(array_type = "450K")
gse51820_beta <- gse51820_beta[!rownames(gse51820_beta) %in% xloci, ]
gse51820_beta <- gse51820_beta[ ,colnames(gse51820_beta) %in% gse51820_pheno$barcode]
gse51820_filter <- champ.filter(beta = gse51820_beta,
                                pd=gse51820_pheno,
                                detP=NULL,
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
                                arraytype = "450K")
message("Probe filtering completed!")

# Perform BMIQ normalisation
gse51820_beta <- champ.norm(beta=gse51820_Filter$beta,
                            resultsDir=paste(dataPath, "51820/CHAMP_Normalization/", sep=""),
                            method="BMIQ",
                            plotBMIQ=FALSE,
                            arraytype="450K",
                            cores=6)

# Plot beta intensities for the different histotypes in the dataset
makeMethBetaFreqPlt(inpBeta = as.data.frame(gse51820_beta), inpPheno = gse51820_pheno, inpSid = "GSE51820", inpDir = paste(dataPath, "GSE51820/",sep=""))

write.csv(gse51820_pheno, paste(dataPath, "GSE51820/GSE51820_pheno.csv", sep=""), row.names = TRUE)
write.csv(gse51820_beta, paste(dataPath, "GSE51820/GSE51820_Beta.csv", sep=""))

################################################################################
# GSE133556 - EPIC
# 99 HGSC, 13 "Normal"
################################################################################
gse133556 = getGEO(filename=GSE133556_path, getGPL = FALSE)
GSE133556_pheno <- gse133556@phenoData@data
GSE133556_pheno$barcode <- gsub("_[a-zA-Z0-9_]{0,3}.idat.gz", "", basename(GSE133556_pheno$supplementary_file))
# Change basename to work with local files
GSE133556_pheno$Basename <- GSE133556_pheno$supplementary_file
GSE133556_pheno$Basename <- basename(GSE133556_pheno$Basename)
GSE133556_pheno$Basename <- paste(dataPath, "GSE133556/", GSE133556_pheno$Basename, sep="")
GSE133556_pheno$Basename <- gsub("_Grn.idat.gz", "", GSE133556_pheno$Basename)
GSE133556_pheno <- tibble::rownames_to_column(GSE133556_pheno, "Sample_ID")
GSE133556_pheno$Histotype <- GSE133556_pheno$'tissue:ch1'
# Rename histotypes to work with scripts 
GSE133556_pheno$Histotype <- plyr::mapvalues(
  GSE133556_pheno$Histotype,
  c("High grade serous ovarian cancer","Normal Fallopian tube control"),
  c("HGSC", "CONTROL"))
gse133556_pheno <- GSE133556_pheno[GSE133556_pheno$Histotype %in% c(histotypes, "CONTROL"),]
# Save phenotypic manifest to local drive for future use
write.csv(gse133556_pheno , paste(dataPath, "GSE133556/GSE133556_pheno.csv", sep=""))

################################################################################
# GSE155760 - EPIC
# 23 HGSC, 37 non HGSC, 36 "Normal"
################################################################################
gse155760 <- getGEO(filename=GSE155760_path, getGPL = FALSE)
gse155760_pheno  <- gse155760@phenoData@data
gse155760_pheno$barcode <- gsub("_[a-zA-Z0-9_]{0,3}.idat.gz", "", basename(gse155760_pheno$supplementary_file))
# Change basename to work with local files
gse155760_pheno$Basename <- gse155760_pheno$supplementary_file
gse155760_pheno$Basename <- basename(gse155760_pheno$Basename)
gse155760_pheno$Basename <- paste(dataPath, "GSE155760/", gse155760_pheno$Basename, sep="")
gse155760_pheno$Basename <- gsub("_Grn.idat.gz", "", gse155760_pheno$Basename)
gse155760_pheno <- tibble::rownames_to_column(gse155760_pheno, "Sample_ID")
gse155760_pheno$Histotype <- gse155760_pheno$'histology:ch1'
# Rename histotypes to work with scripts 
gse155760_pheno$Histotype <- plyr::mapvalues(
  gse155760_pheno$Histotype,
  c('cervical mucosa from cancer-free normal control', 'fallopian tube mucosa from cancer-free normal control', 'endometrial mucosa from cancer-free normal control', 
    'high-grade serous ovarian carcinoma'),
  c("CERV-CONTROL", "CONTROL", "END-CONTROL", "HGSC"))
# Filter to only keep desired files
gse155760_pheno <- gse155760_pheno[gse155760_pheno$Histotype %in% c(histotypes, "CONTROL"), ]
write.csv(gse155760_pheno, paste(dataPath, "GSE155760/gse155760_pheno.csv",sep=""))

################################################################################
# GSE185008
# All (270) CCC
################################################################################
# gse185008_1 - 450K
gse185008_1 = getGEO(filename=GSE185008_1_path)
GSE185008_1_pheno <- gse185008_1@phenoData@data
GSE185008_1_pheno$barcode <- gsub("_[a-zA-Z0-9_]{0,3}.idat.gz", "", basename(GSE185008_1_pheno$supplementary_file))
# Change basename to work with local files
GSE185008_1_pheno$Basename <- GSE185008_1_pheno$supplementary_file
GSE185008_1_pheno$Basename <- basename(GSE185008_1_pheno$Basename)
GSE185008_1_pheno$Basename <- paste(dataPath, "GSE185008/", GSE185008_1_pheno$Basename, sep="")
GSE185008_1_pheno$Basename<- gsub("_Grn.idat.gz", "", GSE185008_1_pheno$Basename)
GSE185008_1_pheno <- tibble::rownames_to_column(GSE185008_1_pheno, "Sample_ID")
GSE185008_1_pheno$Histotype <- GSE185008_1_pheno$`histology:ch1`
# Rename histotypes to work with scripts 
GSE185008_1_pheno$Histotype <- plyr::mapvalues(
  GSE185008_1_pheno$Histotype,
  c('Ovarian Clear Cell Carcinoma (OCCC)'),
  c("CCC"))
# Save phenotypic manifest to local drive for future use
gse185008_1_pheno <- GSE185008_1_pheno
write.csv(gse185008_1_pheno, paste(dataPath, "GSE185008/GSE185008_1_pheno.csv", sep=""))

################################################################################
# gse185008_2 - EPIC
################################################################################
gse185008_2 = getGEO(filename=GSE185008_2_path)
GSE185008_2_pheno <- gse185008_2@phenoData@data
GSE185008_2_pheno$barcode <- gsub("_[a-zA-Z0-9_]{0,3}.idat.gz", "", basename(GSE185008_2_pheno$supplementary_file))
# Change basename to work with local files
GSE185008_2_pheno$Basename <- GSE185008_2_pheno$supplementary_file
GSE185008_2_pheno$Basename <- basename(GSE185008_2_pheno$Basename)
GSE185008_2_pheno$Basename <- paste(dataPath, "GSE185008/", GSE185008_2_pheno$Basename, sep="")
GSE185008_2_pheno$Basename <- gsub("_Grn.idat.gz", "", GSE185008_2_pheno$Basename)
GSE185008_2_pheno <- tibble::rownames_to_column(GSE185008_2_pheno, "Sample_ID")
GSE185008_2_pheno$Histotype <- GSE185008_2_pheno$`histology:ch1`
# Rename histotypes to work with scripts 
GSE185008_2_pheno$Histotype <- plyr::mapvalues(
  GSE185008_2_pheno$Histotype,
  c('Ovarian Clear Cell Carcinoma (OCCC)'),
  c("CCC"))
gse185008_2_pheno <- GSE185008_2_pheno
# Save phenotypic manifest to local drive for future use
write.csv(gse185008_2_pheno, paste(dataPath, "GSE185008/GSE185008_2_pheno.csv", sep=""))

################################################################################
# GSE211686_1 - 450K
# All HGSC
################################################################################
# From dataset GSE65820, rename when doing comparative analysis
gse211686_1= getGEO(filename=GSE211686_1_path, getGPL = FALSE)
GSE211686_1_pheno <- gse211686_1@phenoData@data
GSE211686_1_pheno$barcode <- gsub("_[a-zA-Z0-9_]{0,3}.idat.gz", "", basename(GSE211686_1_pheno$supplementary_file))
# Change basename to work with local files
GSE211686_1_pheno$Basename <- GSE211686_1_pheno$supplementary_file
GSE211686_1_pheno$Basename <- basename(GSE211686_1_pheno$Basename)
GSE211686_1_pheno$Basename <- paste(dataPath, "GSE211686/", GSE211686_1_pheno$Basename, sep="")  
GSE211686_1_pheno$Basename <- gsub("_Grn.idat.gz", "", GSE211686_1_pheno$Basename)
GSE211686_1_pheno <- tibble::rownames_to_column(GSE211686_1_pheno, "Sample_ID")
GSE211686_1_pheno$Histotype <- GSE211686_1_pheno$'characteristics_ch1.3'
# Rename histotypes to work with scripts 
GSE211686_1_pheno$Histotype <- plyr::mapvalues(
  GSE211686_1_pheno$Histotype,
  c('cell type: high grade serous ovarian cancer'),
  c("HGSC"))
# Save phenotypic manifest to local drive for future use
write.csv(GSE211686_1_pheno, paste(dataPath, "GSE211686/GSE211686_1_pheno.csv", sep=""))

################################################################################
# GSE211686_2 - EPIC
# All HGSC
################################################################################
gse211686_2 = getGEO(filename=GSE211686_2_path, getGPL = FALSE)
GSE211686_2_pheno <- gse211686_2@phenoData@data
GSE211686_2_pheno$barcode <- gsub("_[a-zA-Z0-9_]{0,3}.idat.gz", "", basename(GSE211686_2_pheno$supplementary_file))
# Change basename to work with local files
GSE211686_2_pheno$Basename <- GSE211686_2_pheno$supplementary_file
GSE211686_2_pheno$Basename <- basename(GSE211686_2_pheno$Basename)
GSE211686_2_pheno$Basename <- paste(dataPath, "GSE211686/", GSE211686_2_pheno$Basename, sep="")
GSE211686_2_pheno$Basename <- gsub("_Grn.idat.gz", "", GSE211686_2_pheno$Basename)
GSE211686_2_pheno <- tibble::rownames_to_column(GSE211686_2_pheno, "Sample_ID")
GSE211686_2_pheno$Histotype <- GSE211686_2_pheno$'characteristics_ch1.3'
# Rename histotypes to work with scripts 
GSE211686_2_pheno$Histotype <- plyr::mapvalues(
  GSE211686_2_pheno$Histotype,
  c('cell type: high grade serous ovarian cancer'),
  c("HGSC"))
# Save phenotypic manifest to local drive for future use
write.csv(GSE211686_2_pheno, paste(dataPath, "GSE211686/GSE211686_2_pheno.csv", sep=""))

################################################################################
# GSE226823
# 19 CCC, 48 EC, 60 HGSC 
################################################################################

gse226823 = getGEO(filename=GSE226823_path, getGPL = FALSE)
gse226823_pheno <- gse226823@phenoData@data
gse226823_pheno$Basename <- gse226823_pheno$supplementary_file
gse226823_pheno$Basename <- basename(gse226823_pheno$Basename)
gse226823_pheno$Basename <- paste(dataPath, "GSE226823/", gse226823_pheno$Basename, sep="")
gse226823_pheno <- tibble::rownames_to_column(gse226823_pheno, "Sample_ID")
gse226823_pheno$Histotype <- gse226823_pheno$characteristics_ch1.1
# Rename histotypes to work with scripts 
gse226823_pheno$Histotype <- plyr::mapvalues(
  gse226823_pheno$Histotype,
  c("histology: CCOC","histology: ENOC", "histology: HGSC"),
  c('CCC', 'EC', 'HGSC'))
# Filter to only keep desired files
gse226823_pheno <- gse226823_pheno[gse226823_pheno$Histotype %in% histotypes, ]
gse226823_pheno$barcode <- gsub("_[a-zA-Z0-9_]{0,3}.idat.gz", "", basename(gse226823_pheno$supplementary_file))
# Save phenotypic manifest to local drive for future use
write.csv(gse226823_pheno, paste(dataPath, "GSE226823/GSE226823_pheno.csv", sep=""))

################################################################################
# GSE263434 - EPIC
################################################################################

# 285 HGSC 
gse263434=getGEO(filename=GSE263434_path, getGPL = FALSE)
gse263434_pheno  <- gse263434@phenoData@data
gse263434_pheno$barcode <- gsub("_[a-zA-Z0-9_]{0,3}.idat.gz", "", basename(gse263434_pheno$supplementary_file))
# Change basename to work with local files
gse263434_pheno$Basename <- gse263434_pheno$supplementary_file
gse263434_pheno$Basename <- basename(gse263434_pheno$Basename)
gse263434_pheno$Basename <- paste(dataPath, "GSE263434/", gse263434_pheno$Basename, sep="")
gse263434_pheno$Basename <- gsub("_Grn.idat.gz", "", gse263434_pheno$Basename)
gse263434_pheno <- tibble::rownames_to_column(gse263434_pheno, "Sample_ID")
gse263434_pheno$Histotype <- gse263434_pheno$'histology_iti:ch1'
# Rename histotypes to work with scripts 
gse263434_pheno$Histotype <- plyr::mapvalues(
  gse263434_pheno$Histotype,
  c('Clear cell carcinoma', "Endometrioid carcinoma", 'Serous carcinoma'),
  c("RM", "EC", "HGSC"))
# Filter to only keep desired files
gse263434_pheno <- gse263434_pheno[gse263434_pheno$Histotype %in% histotypes, ]

write.csv(gse263434_pheno, paste(dataPath, "GSE263434/gse263434_pheno.csv",sep=""))

################################################################################
# GSE217673 - EPIC
################################################################################

# # 29 HGSC
# gse217673_path<- "C:/Users/xswehu/OneDrive - University of Gothenburg/xswehu_project/Data/Pek_2/GSE217673/Manifest/GSE217673_series_matrix.txt.gz"
# gse217673=getGEO(filename=gse217673_path, getGPL = FALSE)
# # gse263434 <- getGEO(GEO="GSE263434", GSEMatrix = TRUE)
# # , GSEMatrix = TRUE)
# gse217673_pheno  <- gse217673@phenoData@data
# # gse263434_pheno <- read.table("C:/Users/xswehu/OneDrive - University of Gothenburg/xswehu_project/Data/Pek_2/GSE263434/Manifest/GSE263434_series_matrix.txt", sep= "\t")
# gse217673_pheno$barcode <- gsub("_[a-zA-Z0-9_]{0,3}.idat.gz", "", basename(gse217673_pheno$supplementary_file))
# # Change basename to work with local files
# gse217673_pheno$Basename <- gse217673_pheno$supplementary_file
# gse217673_pheno$Basename <- basename(gse217673_pheno$Basename)
# gse217673_pheno$Basename <- paste("C:/Users/xswehu/OneDrive - University of Gothenburg/xswehu_project/Data/Pek_2/GSE217673/Data", gse217673_pheno$Basename, sep="/")
# gse217673_pheno$Basename <- gsub("_Grn.idat.gz", "", gse217673_pheno$Basename)
# gse217673_pheno <- tibble::rownames_to_column(gse217673_pheno, "Sample_ID")
# gse217673_pheno$Histotype <- gse217673_pheno$'characteristics_ch1.1'
# # Rename histotypes to work with scripts 
# gse217673_pheno$Histotype <- plyr::mapvalues(
#   gse217673_pheno$Histotype,
#   c("tissue: Autopsy Tumour"),
#   c("HGSC"))
# # Filter to only keep desired files
# gse217673_pheno <- gse217673_pheno[gse217673_pheno$Histotype %in% histotypes, ]
# # Save phenotypic manifest to local drive for future use
# write.csv(gse217673_pheno, "C:/Users/xswehu/OneDrive - University of Gothenburg/xswehu_project/Data/Pek_2/GSE217673/Manifest/gse217673_pheno.csv")

################################################################################
# GSE267068 - EPIC
################################################################################

# 92 HGSC
gse267068=getGEO(filename=GSE267068_path, getGPL = FALSE)
gse267068_pheno  <- gse267068@phenoData@data
gse267068_pheno$barcode <- gsub("_[a-zA-Z0-9_]{0,3}.idat.gz", "", basename(gse267068_pheno$supplementary_file))
# Change basename to work with local files
gse267068_pheno$Basename <- gse267068_pheno$supplementary_file
gse267068_pheno$Basename <- basename(gse267068_pheno$Basename)
gse267068_pheno$Basename <- paste(dataPath, "GSE267068/", gse267068_pheno$Basename, sep="")
gse267068_pheno$Basename <- gsub("_Grn.idat.gz", "", gse267068_pheno$Basename)
gse267068_pheno <- tibble::rownames_to_column(gse267068_pheno, "Sample_ID")
gse267068_pheno$Histotype <- gse267068_pheno$'description.1'
# Rename histotypes to work with scripts 
gse267068_pheno$Histotype <- plyr::mapvalues(
  gse267068_pheno$Histotype,
  c('high-grade serous ovarian carcinoma'),
  c("HGSC"))
# Filter to only keep desired files
gse267068_pheno <- gse267068_pheno[gse267068_pheno$Histotype %in% histotypes, ]
gse267068_pheno <- gse267068_pheno[gse267068_pheno$source_name_ch1 %in% "Frozen",]
write.csv(gse267068_pheno, paste(dataPath, "GSE267068/gse267068_pheno.csv", sep=""))

################################################################################
# GSE326000 (Train) - EPIC
# 13 CCC, 21 EC, 45 HGSC, 7 MC 
################################################################################

# Training cohort
trainPheno <- getGEO(filename=trainPath , getGPL = FALSE)
trainPheno  <- trainPheno@phenoData@data
# Change basename to work with local files
trainPheno$Basename <- trainPheno$supplementary_file
trainPheno$Basename <- basename(trainPheno$Basename)
trainPheno$Basename <- paste(dataPath, "Train/", trainPheno$Basename, sep="")
trainPheno$Basename <- gsub("_Grn.idat.gz", "", trainPheno$Basename)
trainPheno <- tibble::rownames_to_column(trainPheno, "Sample_ID")
trainPheno$Histotype <- trainPheno$'histological subtype:ch1'
write.csv(trainPheno, paste(dataPath, "Train/trainPheno.csv",sep=""))

################################################################################
# Generate beta-value matrices from .idat files after quality control
################################################################################

extPhenoLst <- list(GSE133556_pheno, gse155760_pheno, 
                    GSE185008_1_pheno, gse185008_2_pheno,
                    GSE211686_1_pheno, GSE211686_2_pheno,
                    gse226823_pheno, gse263434_pheno, 
                    gse267068_pheno, trainPheno)

names(extPhenoLst) <- c("GSE133556", "GSE155760",
                        "GSE185008_1", "GSE185008_2",
                        "GSE211686_1", "GSE211686_2",
                        "GSE226823", "GSE263434",
                        "GSE267068", "Train")

# Process idat files into filtered, normalized beta matrices
for(i in 1:length(extPhenoLst)){
  tmpSid <- names(extPhenoLst)[i]
  tmpPheno <- extPhenoLst[[i]]
  if(tmpSid %in% c("GSE133556","GSE155760","GSE185008_2", "GSE211686_2", "GSE263434", "GSE267068", "Train")){
    arrType <- "EPIC"
  }else{
    arrType <- "450K"
  }
  message("Creating beta-matrix (nbmiq) for cohort;", tmpSid)
  tmpBeta <- makeBetaMatrix(inpPheno = tmpPheno ,
                            inpArrType = arrType, 
                            inpName =  tmpSid, 
                            inpDir = dataPath,  
                            writeBool = TRUE)
}

