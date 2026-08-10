################################################################################
################################################################################
################################################################################
# Scripts associated with predictive classification analysis
################################################################################
################################################################################
################################################################################

################################################################################
# Misc scripts
################################################################################

# Script for retrieving the multiclass model and saving it as a dataframe
getHspModel <- function(inpGeneDf, inpHsp){
  geneVec <- inpGeneDf$GModel
  geneVec <- strsplit(geneVec, "\\+")[[1]]
  modelHsp <- data.frame(matrix(nrow=0, ncol=5))
  colnames(modelHsp) <- c("Histotype", "ensembl_gene_id", "external_gene_name", "noCpg", "sumDist")
  for(i in 1:length(inpHsp)){
    tmpHsp <- inpHsp[[i]]
    tmpH <- names(inpHsp)[i]
    if(nrow(tmpHsp) == 0){
      next()
    }
    tmpMod <- tmpHsp[which(tmpHsp$ensembl_gene_id %in% geneVec),  ]
    for(j in 1:nrow(tmpMod)){
      modelHsp[nrow(modelHsp)+1, ] <- c(tmpH, 
                                        tmpMod$ensembl_gene_id[j],
                                        tmpMod$external_gene_name[j], 
                                        tmpMod$noCpg[j], 
                                        tmpMod$sumDist[j])
    }
  }
  return(modelHsp)
}

getAnnotatedCpgModel <- function(inpCpg, inpPromoLocs, spltBool = NULL){
  # Automated version,issues with duplicates
  tmpCpg <- inpCpg
  # If spltbool is not null, we split the string based on the provided character
  if(!is.null(spltBool)){
    tmpCpg <-  strsplit(tmpCpg[[1]], "\\+")[[1]]
  }
  # Remove any numbers from the str
  names(tmpCpg) <- gsub('[0-9]+', '', names(tmpCpg))
  ensModel <- c()
  inpPromoLocs <- inpPromoLocs[names(inpPromoLocs) %in% names(table(names(tmpCpg)))]
  # Annotate the cpg-model
  skipVec <- c()
  
  for(i in 1:length(inpPromoLocs)){
    tmpLoc <- inpPromoLocs[[i]]
    tmpH <- names(inpPromoLocs)[i]
    matchCpg <- tmpCpg[names(tmpCpg) %in% tmpH]
    for(j in 1:length(matchCpg)){
      tmpC <- matchCpg[j]
      matchGene <- subset(
        tmpLoc,
        rowSums(sapply(tmpLoc, grepl, pattern = tmpC)) > 0
      )
      if(nrow(matchGene) == 0){
        skipVec <- append(skipVec, tmpC)
      }else if(nrow(matchGene > 1)){
        matchGene <- matchGene[1,]
      }
      ensModel <- append(ensModel, matchGene$Gene)
      names(ensModel)[length(ensModel)] <- tmpH
    }
  }
  if(length(skipVec) > 1){
    tmpCpg <- tmpCpg[!tmpCpg %in% skipVec]
  }
  sumVec <- c()
  for(k in 1:length(names(table(gsub('[0-9]+', '', names(tmpCpg)))))){
    tmpH <- names(table(gsub('[0-9]+', '', names(tmpCpg))))[k]
    hInds <- grep(tmpH, names(tmpCpg))
    hNames <- paste(tmpH, 1:length(hInds), sep="")
    names(tmpCpg)[hInds] <- hNames
  }
  names(ensModel) <- names(tmpCpg)
  return(list("CpG" = tmpCpg,"Ens" = ensModel))
}

################################################################################
# Grid selection / Parameter tuning scripts
################################################################################

makeXGBoostCaretGridSMulti <- function(gridDat){
  # Script for xgboost parameter tuning using Caret (superml has incredibly long runtime)
  # Ratio of negative to positive classes
  # class_imba <- c(length(which(gridLabel %in% "Other"))/length(gridLabel), length(which(gridLabel %in% gridCats[!gridCats %in% "Other"]))/length(gridLabel))
  
  # Parameter tuning advicetaken from https://machinelearningmastery.com/configure-gradient-boosting-algorithm/
  # minStart <- round(1/sqrt(table(gridDat$Histotype)[2]/table(gridDat$Histotype)[1]))[[1]]
  # learning_rate and num_boost_round are fixed at 0.1 and 1000 respectively at round 1,2
  # GROUP 1: max_depth , min_child_weight - fixed no.trees and learn-rate
  # GROUP 2: subsample, colsample_bytree - fixed no.trees, learn-rate, depth and child-weight
  # GROUP 3: learning_rate, num_boost_round
  
  # Round 1
  # gridDat$Histotype <- ifelse(gridDat$Histotype == 'Other', 0, 1)
  
  gridTrain <- gridDat[,-grep("Histotype", colnames(gridDat))]
  gridLabel <- gridDat[,grep("Histotype", colnames(gridDat))]
  gridCats <- names(table(gridLabel))
  
  #gridLabel <- as.character(gridLabel)
  #for(i in 1:length(table(gridLabel))){
  #  gridLabel <- replace(gridLabel, gridLabel==names(table(gridDat$Histotype))[[i]], i-1)
  #}
  #gridLabel <- as.numeric(trainLabel)
  #gridLabel <- factor(gridLabel,levels=names(table(gridLabel)))
  
  # Round 1 parameter tuning
  xgb_grid_1 = expand.grid(
    nrounds = c(1000),
    max_depth = c(3, 5, 7, 9),
    eta = c(0.1),
    min_child_weight = c(2, 5, 10, 15, 50),
    gamma = c(0),
    colsample_bytree = c(1),
    subsample=c(1)
  )
  
  # pack the training control parameters
  xgb_trcontrol_1 = caret::trainControl(
    method = "cv",
    number = 3,
    verboseIter = FALSE,
    returnData = FALSE,
    returnResamp = "all",                                                        
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    allowParallel = TRUE
  )
  
  # train the model for each parameter combination in the grid, 
  #   using CV to evaluate
  xgb_train_1 = caret::train(
    x = as.matrix(gridTrain),
    y = gridLabel,
    trControl = xgb_trcontrol_1,
    tuneGrid = xgb_grid_1,
    method = "xgbTree",
    num.threads = 10
  )
  
  # We optimize based on F1 score, as we are dealign with imbalanced datasets
  best_param_1 <- xgb_train_1$results[which.max(xgb_train_1$results$Mean_F1),]
  if(length(which(xgb_train_1$results$Mean_F1 %in% best_param_1$Mean_F1)) > 1){
    roc_rows <- xgb_train_1$results[which(xgb_train_1$results$Mean_F1 %in% best_param_1$Mean_F1),]
    best_param_1 <- roc_rows[which.max(roc_rows$prAUC),]
  }
  best_d <- best_param_1$max_depth
  best_cw <- best_param_1$min_child_weight
  # Round 2
  xgb_grid_2 = expand.grid(
    nrounds = c(1000),
    max_depth = c(best_d),
    eta = c(0.1),
    min_child_weight = c(best_cw),
    gamma = c(0),
    subsample = c(0.5, 0.75, 1.0),
    colsample_bytree = c(0.5, 0.6, 0.8, 1.0)
  )
  xgb_train_2 = caret::train(
    x = as.matrix(gridTrain),
    y = gridLabel,
    trControl = xgb_trcontrol_1,
    tuneGrid = xgb_grid_2,
    method = "xgbTree",
    num.threads = 10
  )
  best_param_2 <- xgb_train_2$results[which.max(xgb_train_2$results$Mean_F1),]
  if(length(which(xgb_train_2$results$Mean_F1 %in% best_param_2$Mean_F1)) > 1){
    roc_rows_2 <- xgb_train_2$results[which(xgb_train_2$results$Mean_F1 %in% best_param_2$Mean_F1),]
    best_param_2 <- roc_rows_2[which.max(roc_rows_2$prAUC),]
  }
  best_ss <- best_param_2$subsample
  best_cst <- best_param_2$colsample_bytree
  # Round 3
  xgb_grid_3 = expand.grid(
    nrounds = c(100, 150, 300, 500, 1000),
    max_depth = c(best_d),
    eta = c(0.01, 0.025, 0.1, 0.15, 0.2, 0.3),
    min_child_weight = c(best_cw),
    gamma = c(0),
    subsample = c(best_ss),
    colsample_bytree = c(best_cst)
  )
  
  # train the model for each parameter combination in the grid, 
  #   using CV to evaluate
  xgb_train_3 = caret::train(
    x = as.matrix(gridTrain),
    y = gridLabel,
    trControl = xgb_trcontrol_1,
    tuneGrid = xgb_grid_3,
    method = "xgbTree",
    num.threads = 10
  )
  best_param_3 <- xgb_train_3$results[which.max(xgb_train_3$results$Mean_F1),]
  if(length(which(xgb_train_3$results$Mean_F1 %in% best_param_3$Mean_F1)) > 1){
    roc_rows_3 <- xgb_train_3$results[which(xgb_train_3$results$Mean_F1 %in% best_param_3$Mean_F1),]
    best_param_3 <- roc_rows_3[which.max(roc_rows_3$prAUC),]
  }
  best_eta <- best_param_3$eta
  best_nr <- best_param_3$nrounds
  return(c("eta"= best_eta, "n_estimators" = best_nr, 
           "subsample" = best_ss, "colsample_bytree"= best_cst,
           "max_depth" = best_d, "min_child_weight"=best_cw))
}

makeXGBoostCaretGridS <- function(gridDat){
  # Script for xgboost parameter tuning using Caret (superml has incredibly long runtime)
  
  # Ratio of negative to positive classes
  
  gridTrain <- gridDat[,-which(colnames(gridDat) %like% "Histotype")]
  gridLabel <- gridDat[,which(colnames(gridDat) %like% "Histotype")]
  
  gridCats <- names(table(gridLabel))
  class_imba <- c(length(which(gridLabel %in% "Other"))/length(gridLabel), length(which(gridLabel %in% gridCats[!gridCats %in% "Other"]))/length(gridLabel))
  
  # Parameter tuning advicetaken from https://machinelearningmastery.com/configure-gradient-boosting-algorithm/
  minStart <- round(1/sqrt(table(gridDat$Histotype)[2]/table(gridDat$Histotype)[1]))[[1]]
  # learning_rate and num_boost_round are fixed at 0.1 and 1000 respectively at round 1,2
  # GROUP 1: max_depth , min_child_weight - fixed no.trees and learn-rate
  # GROUP 2: subsample, colsample_bytree - fixed no.trees, learn-rate, depth and child-weight
  # GROUP 3: learning_rate, num_boost_round
  
  # Create weights for random forest classification (balance)
  xgWeights <-  1/table(gridDat$Histotype)
  xgWeights <- xgWeights/sum(xgWeights)
  xgWVec <-  ifelse(gridDat$Histotype == 'Other', xgWeights[1], xgWeights[2])
  
  # Round 1
  # gridDat$Histotype <- ifelse(gridDat$Histotype == 'Other', 0, 1)
  
  gridTrain <- gridDat[,-which(colnames(gridDat) %like% "Histotype")]
  gridLabel <- gridDat[,which(colnames(gridDat) %like% "Histotype")]
  
  # Round 1 parameter tuning
  xgb_grid_1 = expand.grid(
    nrounds = c(1000),
    max_depth = c(3, 5, 7, 9),
    eta = c(0.1),
    min_child_weight = c(minStart, 5, 10, 15, 50),
    gamma = c(0),
    colsample_bytree = c(1),
    subsample=c(1)
  )
  
  # pack the training control parameters
  xgb_trcontrol_1 = caret::trainControl(
    method = "cv",
    number = 3,
    verboseIter = FALSE,
    returnData = FALSE,
    returnResamp = "all",                                                        
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    allowParallel = TRUE
  )
  
  # train the model for each parameter combination in the grid, 
  #   using CV to evaluate
  xgb_train_1 = caret::train(
    x = as.matrix(gridTrain),
    y = gridLabel,
    trControl = xgb_trcontrol_1,
    tuneGrid = xgb_grid_1,
    method = "xgbTree",
    num.threads = 10
  )
  best_param_1 <- xgb_train_1$results[which.max(xgb_train_1$results$ROC),]
  if(length(which(xgb_train_1$results$ROC %in% best_param_1$ROC)) > 1){
    roc_rows <- xgb_train_1$results[which(xgb_train_1$results$ROC %in% best_param_1$ROC),]
    best_param_1 <- roc_rows[which.max(roc_rows$Sens),]
  }
  best_d <- best_param_1$max_depth
  best_cw <- best_param_1$min_child_weight
  # Round 2
  xgb_grid_2 = expand.grid(
    nrounds = c(1000),
    max_depth = c(best_d),
    eta = c(0.1),
    min_child_weight = c(best_cw),
    gamma = c(0),
    subsample = c(0.5, 0.75, 1.0),
    colsample_bytree = c(0.5, 0.6, 0.8, 1.0)
  )
  xgb_train_2 = caret::train(
    x = as.matrix(gridTrain),
    y = gridLabel,
    trControl = xgb_trcontrol_1,
    tuneGrid = xgb_grid_2,
    method = "xgbTree",
    num.threads = 10
  )
  best_param_2 <- xgb_train_2$results[which.max(xgb_train_2$results$ROC),]
  if(length(which(xgb_train_2$results$ROC %in% best_param_2$ROC)) > 1){
    roc_rows_2 <- xgb_train_2$results[which(xgb_train_2$results$ROC %in% best_param_2$ROC),]
    best_param_2 <- roc_rows_2[which.max(roc_rows_2$Sens),]
  }
  best_ss <- best_param_2$subsample
  best_cst <- best_param_2$colsample_bytree
  # Round 3
  xgb_grid_3 = expand.grid(
    nrounds = c(100, 150, 300, 500, 1000),
    max_depth = c(best_d),
    eta = c(0.01, 0.025, 0.1, 0.15, 0.2, 0.3),
    min_child_weight = c(best_cw),
    gamma = c(0),
    subsample = c(best_ss),
    colsample_bytree = c(best_cst)
  )
  
  # train the model for each parameter combination in the grid, 
  #   using CV to evaluate
  xgb_train_3 = caret::train(
    x = as.matrix(gridTrain),
    y = gridLabel,
    trControl = xgb_trcontrol_1,
    tuneGrid = xgb_grid_3,
    method = "xgbTree",
    num.threads = 10
  )
  best_param_3 <- xgb_train_3$results[which.max(xgb_train_3$results$ROC),]
  if(length(which(xgb_train_3$results$ROC %in% best_param_3$ROC)) > 1){
    roc_rows_3 <- xgb_train_3$results[which(xgb_train_3$results$ROC %in% best_param_3$ROC),]
    best_param_3 <- roc_rows_3[which.max(roc_rows_3$Sens),]
  }
  best_eta <- best_param_3$eta
  best_nr <- best_param_3$nrounds
  return(c("eta"= best_eta, "n_estimators" = best_nr, 
           "subsample" = best_ss, "colsample_bytree"= best_cst,
           "max_depth" = best_d, "min_child_weight"=best_cw))
}

makeXGBoostGridS <- function(gridDat){
  # Function for performing sequential grid search for xgboost parameters
  
  # Grid final three step sequential tuning;
  # https://blog.cambridgespark.com/hyperparameter-tuning-in-xgboost-4ff9100a3b2f
  # learning_rate and num_boost_round are fixed at 0.1 and 1000 respectively at round 1,2
  # GROUP 1: max_depth , min_child_weight - fixed no.trees and learn-rate
  # GROUP 2: subsample, colsample_bytree - fixed no.trees, learn-rate, depth and child-weight
  # GROUP 3: learning_rate, num_boost_round
  
  # colsample_bytree. Lower ratios avoid over-fitting.
  # max_depth. Lower values avoid over-fitting. 
  # gamma. Larger values avoid over-fitting
  # eta. Lower values avoid over-fitting.
  # min_child_weight. Larger values avoid over-fitting.
  # Ratio of negative to positive classes
  # weight. Vector of class-weights to be used for classification
  
  # Create weights for random forest classification (balance)
  xgWeights <-  1/table(gridDat$Histotype)
  xgWeights <- xgWeights/sum(xgWeights)
  xgWVec <-  ifelse(gridDat$Histotype == 'Other', xgWeights[1], xgWeights[2])
  # Parameter tuning advicetaken from https://machinelearningmastery.com/configure-gradient-boosting-algorithm/
  minStart <- round(1/sqrt(table(gridDat$Histotype)[2]/table(gridDat$Histotype)[1]))[[1]]
  # Round 1
  gridDat$Histotype <- ifelse(gridDat$Histotype == 'Other', 0, 1)
  #class_imba <- length(which(gridDat$Histotype %in% 0))/length(which(gridDat$Histotype %in% 1))
  xgb <-XGBTrainer$new()
  gst_1 <- superml::GridSearchCV$new(trainer = xgb,
                                     parameters = list(n_estimators = 100,
                                                       learning_rate = 0.1,
                                                       max_depth =  c(3,5,7,9),
                                                       min_child_weight = c(minStart, 5, 10, 15, 50),
                                                       verbose = 0,
                                                       early_stopping = 10,
                                                       nthread= 10
                                                       # weight =  xgWVec
                                     ),
                                     n_folds = 3,
                                     scoring = c('auc'))
  gridCv_1 <- gst_1$fit(gridDat, "Histotype")
  best_param_1 <-  gridCv_1[which.max(gridCv_1$auc_avg),]
  best_d <- best_param_1$max_depth
  best_cw <- best_param_1$min_child_weight
  # Round 2
  gst_2 <- superml::GridSearchCV$new(trainer = xgb,
                                     parameters = list(n_estimators = 100,
                                                       learning_rate = 0.1,
                                                       # scale_pos_weight = class_imba,
                                                       max_depth =  best_d,
                                                       min_child_weight = best_cw,
                                                       subsample = c(0.4, 0.6, 0.8, 1.0),
                                                       colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
                                                       verbose = 0,
                                                       early_stopping = 10,
                                                       nthread= 10
                                                       # weight =  xgWVec
                                     ),
                                     n_folds = 3,
                                     scoring = c('auc'))
  gridCv_2 <- gst_2$fit(gridDat, "Histotype")
  best_param_2 <-  gridCv_2[which.max(gridCv_2$auc_avg),]
  best_ss <- best_param_2$subsample
  best_cst <- best_param_2$colsample_bytree 
  # Round 3
  gst_3 <- superml::GridSearchCV$new(trainer = xgb,
                                     parameters = list(n_estimators = c(100, 150, 300, 500, 1000),
                                                       learning_rate = c(0.01, 0.025, 0.1, 0.15, 0.2, 0.3),
                                                       # scale_pos_weight = class_imba,
                                                       max_depth =  best_d,
                                                       min_child_weight = best_cw,
                                                       subsample = best_ss,
                                                       colsample_bytree = best_cst,
                                                       verbose = 0,
                                                       early_stopping = 10,
                                                       nthread= 10
                                                       # weight =  xgWVec
                                     ),
                                     n_folds = 3,
                                     scoring = c('auc'))
  gridCv_3 <- gst_3$fit(gridDat, "Histotype")
  best_param_3 <-  gridCv_3[which.max(gridCv_3$auc_avg),]
  best_eta <- best_param_3$learning_rate
  best_n_est <- best_param_3$n_estimators
  return(c("eta"= best_eta, "n_estimators" = best_n_est, 
           "subsample" = best_ss, "colsample_bytree"= best_cst,
           "max_depth" = best_d, "min_child_weight"=best_cw))
}

################################################################################
# Model selection scripts
################################################################################

makeStepwiseMultModel <- function(inpCpgDfLst, refBeta, refPheno, inpCpgPos, 
                                  inpCpgLocs=NULL, noEpochs = NULL, nFold = NULL, 
                                  paramBool = NULL, rmHist = NULL, minCpg = NULL, 
                                  itMax = NULL, fullPromoBool=NULL, nCatBool = NULL){
  # Script for creating a model for multiple classification, given multiple defined genomic regions
  if(is.null(nFold)){
    nFold <- 5
  }
  if(is.null(noEpochs)){
    noEpochs <- 1
  }
  if(is.null(minCpg)){
    minCpg <- 2
  }
  
  if(is.null(itMax)){
    itMax <- 5
  }
  
  nFold<- nFold
  if(!is.null(rmHist)){
    refPheno <- refPheno[!refPheno$Histotype %in% rmHist, ]
    refBeta <- refBeta[,refPheno$barcode]
  }
  refPheno <- refPheno[which(refPheno$barcode %in% colnames(refBeta)),]
  refBeta <- refBeta[,which(colnames(refBeta) %in% refPheno$barcode)]
  refM <- log2(refBeta/(1-refBeta))
  
  genePredLst <- list()
  for(j in 1:length(inpCpgDfLst)){
    # Second version, perform pred-class on individiual genes
    genePredDf <- data.frame(matrix(nrow=0, ncol=10))
    colnames(genePredDf) <- c("ensembl_gene_id", 
                              "Sensitivity", "Specificity","Accuracy", "F1",
                              "CI_Lower", "CI_Upper",
                              "AUC","AUPRC", 
                              "CpGs")
    inpCpgDf <- inpCpgDfLst[[j]]
    inpHist <- names(inpCpgDfLst)[j]
    # Specify the number of phenotype categories the site should be significant against
    if(is.null(nCatBool)){
      nCat <- length(table(refPheno$Histotype))-1
    }else{
      nCat <- nCatBool
    }
    if(is.null(nrow(inpCpgDf))){
      next()
    }else if(nrow(inpCpgDf) == 0){
      next()
    }
    for(k in 1:nrow(inpCpgDf)){
      if(abs(k)%%round(nrow(inpCpgDf)/10) == 0){
        message(paste("Processing row: [", k, "/", nrow(inpCpgDf),"]"))
      }
      tmpCpgDf <- data.frame(inpCpgDf[k,])
      tmpEns <- tmpCpgDf$Gene
      # Prepare reference dataframe
      promoDf <- makeCpgDf(tmpCpgDf, refBeta, refPheno, inpHist)
      if(is.null(nrow(promoDf))){
        message(paste("Gene: ", tmpEns, " Not found in array-data and will be skipped", sep=""))
        next()
      }else if(nrow(promoDf) <= minCpg){
        message(paste("Gene: ", tmpEns, " Has less then or equal to: ",  minCpg, " total CpG sites and will be skipped", sep=""))
        next()
      }else{
        promoDf <- promoDf[rownames(promoDf) %in% rownames(refBeta),]
        # Either get the significant CpGs, or model the full promoter region
        if(!is.null(fullPromoBool)){
          minDists <- nrow(promoDf)
        }else{
          fakeDf <- data.frame(nrow=1, ncol=2)
          fakeDf[1,] <- c(tmpEns, inpHist)
          colnames(fakeDf) <- c("ensembl_gene_id", "Histotype")
          fakeLst <- list(fakeDf)
          names(fakeLst)[length(fakeLst)] <- inpHist
          fakeCpgLst <- list(promoDf)
          names(fakeCpgLst) <- tmpEns
          minDists <- makeRegSigCp(inpBeta = refBeta,
                                      inpM = refM ,
                                      inpPheno = refPheno,
                                      inpGenes = fakeLst,
                                      inpPromoCpgs = fakeCpgLst,
                                      noCats = nCat)
          minDists <- minDists[[1]]
        }
        if(length(minDists) < 2){
          next()
        }else{
          # Prepare training-data  based on the less stringent filtering
          refPredDf <- refBeta[minDists, ]
          refPredDf <- as.data.frame(t(refPredDf))
          refPredDf  <- refPredDf[rowSums(is.na(refPredDf)) == 0,]
          refPredDf$Histotype <- refPheno$Histotype[match(rownames(refPredDf), refPheno$barcode)]
          refPredDf$Histotype <- ifelse(refPredDf$Histotype %in% inpHist, inpHist, "Other")
          refPredDf$Histotype <- factor(refPredDf$Histotype, levels = c("Other", inpHist))
          # Create models to be used with data
          minModel <- paste(minDists, collapse = " + ")
          minModel <- reformulate(minModel, "Histotype")
          # Shuffle rows to remove any bias when sampling
          refPredDf <- refPredDf[sample(nrow(refPredDf)),]
          kfDfLst <- list()
          if(!is.null(paramBool)){
            # Perform grid-search to identify best parameters for model
            best_params <- makeXGBoostCaretGridS(refPredDf)
          }
          # Perform xgboost on training cohort
          sumCvDf <- data.frame(matrix(NA, nrow=noEpochs, ncol=8))
          colnames(sumCvDf) <- c("Sensitivity", "Specificity", "Accuracy","F1", 
                                 "CI_Lower", "CI_Upper", 
                                 "AUC", "AUPRC")
          inpRef <- refPredDf
          kfDfLst <- list()
          for(l in 1:noEpochs){
            folds <- caret::createFolds(factor(inpRef$Histotype),
                                        k = nFold, 
                                        list = FALSE) 
            # predStatLst <- list()
            kfDf <- data.frame(matrix(nrow=0, ncol=8))
            colnames(kfDf) <- c("Sensitivity", "Specificity", "Accuracy", "F1",
                                "CI_Lower","CI_Upper", 
                                "AUC", "AUPRC")
            for(m in 1:nFold){
              # Segment train-data by folds, with respect to the fold index no. 
              testIndexes <- which(folds==m, arr.ind=TRUE)
              trainData <- inpRef[-testIndexes, ]
              testData <- inpRef[testIndexes, ]
              # best_params <- makeXGBoostCaretGridS(trainData)
              # Create weights for random forest classification (balance)
              trainVals <- trainData[,-which(grepl("Histotype", colnames(trainData)))]
              testVals <- testData[,-which(grepl("Histotype", colnames(testData)))]
              # Set up labels for test and training data
              trainLabel <- trainData[,which(grepl("Histotype",colnames(trainData)))]
              trainLabel <- ifelse(trainLabel == inpHist, 1, 0)
              trainLabel <- factor(trainLabel, levels = c(0,1))
              testLabel <- testData[,which(grepl("Histotype",colnames(testData)))]
              testLabel <- ifelse(testLabel == inpHist, 1, 0)
              testLabel <- factor(testLabel, levels = c(0,1))
              # Ratio of negative to positive classes
              # class_imba <- length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
              class_imba <- length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
              # best_params <- makeXGBoostGridS(trainData)
              # Perform grid-search on training data to retrieve optimal model parameters
              # best_params <- makeXGBoostGridS(trainData)
              # Build gradient boosting forest model
              if(is.null(paramBool)){
                tmpBoost <- xgboost::xgboost(x = as.matrix(trainVals),
                                             y = trainLabel,
                                             booster = "gbtree",
                                             objective = "binary:logistic",
                                             scale_pos_weight = class_imba,
                                             nrounds = 100,
                                             nthread = 10)
              }else{
                tmpBoost <- xgboost(data = as.matrix(trainVals),
                                    label = trainLabel,
                                    booster = "gbtree",
                                    objective = "binary:logistic",
                                    eta =  best_params["eta"] ,
                                    max_depth =  best_params["max_depth"],
                                    scale_pos_weight = class_imba,
                                    min_child_weight = best_params["min_child_weight"],
                                    subsample = best_params["subsample"],
                                    colsample_bytree = best_params["colsample_bytree"],
                                    nrounds = best_params["n_estimators"],
                                    verbose = 0,
                                    nthread = 10)
              }
              # Apply RF model on test-data (external dataset)
              tmpExtPred <- predict(object=tmpBoost, 
                                    newdata = as.matrix(testVals))
              # Categorize predictions based on their corresponding classification seen to model
              conMat <- confusionMatrix(data = factor(as.numeric(tmpExtPred > 0.5), levels = c(1,0)), 
                                        reference = factor(testLabel, levels=c(1,0)))
              # As we are dealing with heavily imbalanced data, optimize on f1 score and AUPRC
              sens <- conMat$byClass[[1]]
              spec <- conMat$byClass[[2]]
              f1Score <- conMat$byClass[[7]]
              acc <- conMat$overall[[1]]
              # Get confidence interval of statistics
              confIntLower <- conMat$overall[[3]]
              confIntUpper <- conMat$overall[[4]]
              # Get AUC for model
              roc <- pROC::roc(testLabel, tmpExtPred, smooth=FALSE,  quiet = TRUE)
              # Get AUPRC
              # Indice locations of pos and neg-classes
              indice_POS=which(testLabel %in% 1)
              indice_NEG=which(testLabel %in% 0)
              # Get scores of pos and neg classes
              clas_score_POS=tmpExtPred[indice_POS]
              clas_score_NEG=tmpExtPred[indice_NEG]
              # Calculate area under precision recall curve
              # testLabelFreqs <- table(testLabel)/length(testLabel)
              # testLabelW <-ifelse(testLabel %in% names(testLabelFreqs)[1], testLableFreqs[1], testLableFreqs[2])
              auprc <-PRROC::pr.curve(scores.class0 = clas_score_POS,
                                      scores.class1 = clas_score_NEG)
              # Save to list/DF
              kfDf[nrow(kfDf)+1,] <- c(sens, spec, acc, f1Score, 
                                       confIntLower, confIntUpper, 
                                       roc$auc[1], auprc$auc.integral[1])
              rownames(kfDf)[nrow(kfDf)] <- paste("k=",m, sep="")
              
            }
            kfDf[sapply(kfDf, is.nan)] <- 0
            kfDf[sapply(kfDf, is.na)] <- 0
            kfDf[nrow(kfDf)+1, ] <- colMeans(kfDf)
            sumCvDf[l,] <- colMeans(kfDf)
          }
          meanPred <- colMeans(sumCvDf)
          genePredDf[nrow(genePredDf)+1,] <- c(tmpEns, 
                                               meanPred[1], meanPred[2], meanPred[3], meanPred[4],
                                               meanPred[5], meanPred[6],
                                               meanPred[7], meanPred[8],
                                               paste(minDists, collapse = ","))
        }
      }
    }
    # Select only models with a good combined AUC
    # genePredDf <- genePredDf[which(genePredDf$F1 > 0.7), ]
    genePredLst[[length(genePredLst)+1]] <- genePredDf
    names(genePredLst)[length(genePredLst)] <- inpHist
  }
  message("Individual gene-pred completed")
  
  # Remove empty entries
  genePredLst <- genePredLst[!names(genePredLst) %in% ""]
  genePredLst <- genePredLst[sapply(genePredLst, nrow) > 0]
  modelLst <- list()
  for(q in 1:length(genePredLst)){
    bestRow <- NULL
    genePredLst2 <- list()
    sumDf <- data.frame(matrix(nrow=0, ncol=11))
    colnames(sumDf) <- c("Histotype","Model", 
                         "Sensitivity", "Specificity", "Accuracy", "F1",
                         "CI_Lower","CI_Upper",
                         "AUC", "AUPRC",
                         "CpGs")
    tmpH <- names(genePredLst)[q]
    gDf <- genePredLst[[tmpH]]
    nSamp <- nrow(gDf)
    if(is.null(nSamp)){
      next()
    }else if(nSamp == 1){
      nIt <- 1 
    }else if(nSamp < itMax){
      nIt <- nSamp-1
    }else{
      nIt <- itMax-1
    }
    for(m in 1:nIt){
      genePredDf2 <- data.frame(matrix(nrow=0, ncol=11))
      colnames(genePredDf2) <- c("Histotype", "Model", 
                                 "Sensitivity", "Specificity", "Accuracy", "F1",
                                 "CI_Lower","CI_Upper",
                                 "AUC", "AUPRC",
                                 "CpGs")
      # Select best model best on F1 score, if tied base it on AUPRC for model
      if(is.null(bestRow)){
        bestRow <- gDf[which.max(gDf$F1),]
        bestInd <- which.max(gDf$F1)
        if(nrow(bestRow) > 1){
          bestRow <- bestRow[which.max(bestRow$AUPRC),]
          bestInd <- bestInd[-which.min(bestRow$AUPRC)]
        }
        ens1 <- bestRow$ensembl_gene_id
        if(nIt == 1){
          tmpDf <- gDf
        }else{
          tmpDf <- gDf[-bestInd, ]
        }
        # Create new model
        pred1 <- bestRow$CpGs
        pred1 <- strsplit(pred1, ",")[[1]]
        pred1str <- paste(pred1, collapse = " + ")
      }else{
        bestRow <- sumDf[nrow(sumDf),]
        ens1 <- bestRow$Model
        ens1 <- strsplit(ens1, "\\+")[[1]]
        tmpDf <- gDf[-which(gDf$ensembl_gene_id %in% ens1), ]
        # Create new model
        bPreds <- gDf$CpGs[which(gDf$ensembl_gene_id %in% ens1)]
        bPreds <- paste(bPreds, collapse = ",")
        pred1 <- strsplit(bPreds, ",")[[1]]
        pred1str <- paste(pred1, collapse = " + ")
      }
      
      for(n in 1:nrow(tmpDf)){
        nRow <- tmpDf[n,]
        ens2 <- nRow$ensembl_gene_id
        pred2 <- nRow$CpGs
        pred2 <- strsplit(pred2, ",")[[1]]
        pred2str <- paste(pred2, collapse = " + ")
        nMod <- paste(pred1str, pred2str, sep = " + ")
        nMod <- reformulate(nMod, "Histotype")  
        modPreds <- c(pred1, pred2)
        
        # Try new model
        modDf <- refBeta[modPreds, ]
        modDf <- as.data.frame(t(modDf))
        modDf  <- modDf[rowSums(is.na(modDf)) == 0,]
        modDf$Histotype <- refPheno$Histotype[match(rownames(modDf), refPheno$barcode)]
        modDf$Histotype <- ifelse(modDf$Histotype %in% tmpH, 1, 0)
        modDf$Histotype <- factor(modDf$Histotype, levels = c(1, 0))
        # Shuffle rows to remove any bias when sampling
        modDf <- modDf[sample(nrow(modDf)),]
        # Perform xgboost on training cohort
        sumCvDf2 <- data.frame(matrix(NA, nrow=0, ncol=8))
        colnames(sumCvDf2) <- c("Sensitivity", "Specificity", "Accuracy", "F1",
                                "CI_Lower", "CI_Upper",
                                "AUC", "AUPRC")
        folds <- caret::createFolds(modDf$Histotype, k = nFold, list = FALSE) 
        # predStatLst <- list()
        kfDf2 <- data.frame(matrix(nrow=0, ncol=8))
        colnames(kfDf2) <- colnames(kfDf2) <- c("Sensitivity", "Specificity", "Accuracy", "F1",
                                                "CI_Lower", "CI_Upper",
                                                "AUC", "AUPRC")
        for(p in 1:nFold){
          # Segment train-data by folds, with respect to the fold index no. 
          testIndexes <- which(folds==p, arr.ind=TRUE)
          testData <- modDf[testIndexes, ]
          trainData <- modDf[-testIndexes, ]
          # best_params <- makeXGBoostCaretGridS(trainData)
          # Create weights for random forest classification (balance)
          trainVals <- trainData[,-which(grepl("Histotype",colnames(trainData)))]
          trainLabel <- unfactor(trainData[,which(grepl("Histotype",colnames(trainData)))])
          #trainLabel <- factor(trainLabel, levels=c(1, 0))
          testVals <- testData[,-which(grepl("Histotype",colnames(testData)))]
          testLabel <- unfactor(testData[,which(grepl("Histotype", colnames(testData)))])
          #testLabel <- factor(testLabel, levels=c(1, 0))
          # Ratio of negative to positive classes
          class_imba <- length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
          # length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
          # best_params <- makeXGBoostGridS(trainData)
          # Perform grid-search on training data to retrieve optimal model parameters
          # best_params <- makeXGBoostGridS(trainData)
          # Build gradient boosting forest model
          if(is.null(paramBool)){
            tmpBoost <- xgboost(data = as.matrix(trainVals),
                                label = trainLabel,
                                booster = "gbtree",
                                objective = "binary:logistic",
                                scale_pos_weight = class_imba,
                                nrounds = 100,
                                verbose = 0,
                                nthread = 10)
          }else{
            tmpBoost <- xgboost(data = as.matrix(trainVals),
                                label = trainLabel,
                                booster = "gbtree",
                                objective = "binary:logistic",
                                eta =  best_params["eta"] ,
                                max_depth =  best_params["max_depth"],
                                scale_pos_weight = class_imba,
                                min_child_weight = best_params["min_child_weight"],
                                subsample = best_params["subsample"],
                                colsample_bytree = best_params["colsample_bytree"],
                                nrounds = best_params["n_estimators"],
                                verbose = 0,
                                nthread = 10)
          }
          
          # Apply RF model on test-data (external dataset)
          tmpExtPred <- predict(object=tmpBoost, 
                                newdata = as.matrix(testVals))
          conMat <- confusionMatrix(data = factor(as.numeric(tmpExtPred >0.5), levels = c(1,0)), 
                                    reference = factor(testLabel, levels = c(1,0)))
          
          # Get AUC for model
          roc <- pROC::roc(testLabel, tmpExtPred, smooth=FALSE,  quiet = TRUE)
          # As we are dealing with heavily imbalanced data, optimize on f1 score and AUPRC
          sens <- conMat$byClass[[1]]
          spec <- conMat$byClass[[2]]
          f1Score <- conMat$byClass[[7]]
          acc <- conMat$overall[[1]]
          # Get confidence interval of statistics
          confIntLower <- conMat$overall[[3]]
          confIntUpper <- conMat$overall[[4]]
          # Get AUC for model
          roc <- pROC::roc(testLabel, tmpExtPred, smooth=FALSE,  quiet = TRUE)
          
          # Get AUPRC
          # Indice locations of pos and neg-classes
          indice_POS=which(testLabel %in% 1)
          indice_NEG=which(testLabel %in% 0)
          # Get scores of pos and neg classes
          clas_score_POS=tmpExtPred[indice_POS]
          clas_score_NEG=tmpExtPred[indice_NEG]
          # Calculate area under precision recall curve
          # Scores class 0 is equal to the positive class (i.e. histotype)
          auprc <-PRROC::pr.curve(scores.class1 = clas_score_NEG, 
                                  scores.class0 = clas_score_POS)
          
          # Save to dataframe
          kfDf2[nrow(kfDf2)+1,] <- c(sens, spec, acc, f1Score,
                                     confIntLower, confIntUpper,
                                     roc$auc[1], auprc$auc.integral[1])
          kfDf2[sapply(kfDf2, is.nan)] <- 0
          #kfDf2[nrow(kfDf2)+1, ] <- colMeans(kfDf2)
          #rownames(kfDf2)[nrow(kfDf2)] <- "MeanVal"
          sumCvDf2[nrow(sumCvDf2)+1,] <- colMeans(kfDf2)
          #rownames(kfDf2)[nrow(kfDf2)] <- paste("k=",p, sep="")
        }
        meanPred2 <- colMeans(sumCvDf2)
        genePredDf2[nrow(genePredDf2)+1,] <- c(tmpH, paste(paste(ens1, collapse="+"), ens2, sep="+"), 
                                               meanPred2[1], meanPred2[2], meanPred2[3], meanPred2[4],
                                               meanPred2[5],meanPred2[6],
                                               meanPred2[7],meanPred2[8],
                                               paste(modPreds, collapse=","))
      }
      genePredDf2 <- genePredDf2[order(genePredDf2$F1, decreasing=TRUE), ]
      if(length(which.max(genePredDf2$F1)) > 2){
        genePredDf2AUPRC <- genePredDf2[which.max(genePredDf2$F1),]
        bestModel <- genePredDf2AUPRC[which.max(genePredDf2AUPRC$AUPRC),]
      }else{
        if(length(which(is.na(genePredDf2$F1))) > 1){
          bestModel <- genePredDf2[which.max(genePredDf2$AUPRC),]
        }else{
          bestModel <- genePredDf2[which.max(genePredDf2$F1),]
        }
      }
      sumDf[nrow(sumDf)+1,] <- bestModel 
      bestRow <- sumDf[nrow(sumDf),]
    }
    modelLst[[q]] <- sumDf
    names(modelLst)[q] <- tmpH
  }
  message("Model combination analysis completed")
  
  # Combine all different CpG-sets, 
  # Use the most beneficial combination of promoters as a gene-panel for histotype stratification
  fullModel <- c()
  geneModel <- c()
  bModelDf <- data.frame(matrix(nrow=0, ncol=ncol(modelLst[[1]])))
  colnames(bModelDf) <- colnames(modelLst[[1]])
  for(r in 1:length(modelLst)){
    tmpModels <- modelLst[[r]]
    tmpH <- names(modelLst)[r]
    # tmpModels$SensSpec <- as.numeric(tmpModels$Sensitivity) + as.numeric(tmpModels$Specificity)
    bModel <- tmpModels[which.max(tmpModels$F1), ]
    if(nrow(bModel) >1){
      bModel <- bModel[which.max(tmpModels$AUPRC), ]
    }
    tmpGM <- bModel$Model
    tmpGM <- strsplit(tmpGM , "\\+")[[1]]
    names(tmpGM) <- paste(tmpH, 1:length(tmpGM), sep="")
    tmpFM <- bModel$CpGs
    tmpFM <- strsplit(tmpFM , ",")[[1]]
    # Need to chang the names here (ens-id's over CpGs, Histotypes over Ens-genes), alternatively call script for getting the model form the CpG's
    names(tmpFM) <-  paste(tmpH, 1:length(tmpFM), sep="")
    fullModel <- append(fullModel, tmpFM)
    geneModel <- append(geneModel, tmpGM)
    bModelDf[nrow(bModelDf)+1,] <- bModel 
  }
  bModelDf$Histotype <- names(modelLst)
  message("Combined best model created")
  return(list("CpGModel" = fullModel, "ensModel" = geneModel, "statDf" = bModelDf, "modelLst" = modelLst))
}

makeGenePredLstComb_MULT <- function(genePredLst, refBeta, refPheno, 
                                     noEpochs = NULL, nFold = NULL, 
                                     paramBool = NULL, itMax = NULL){
  # Function for finding the optimal combination of models from genomic regions of interest, for predictive classification purposes
  if(is.null(nFold)){
    nFold <- 5
  }
  if(is.null(noEpochs)){
    noEpochs <- 1
  }
  genePredLst <- genePredLst[!names(genePredLst) %in% ""]
  genePredLst <- genePredLst[sapply(genePredLst, nrow) > 0]
  modelLst <- list()
  for(i in 1:length(genePredLst)){
    bestRow <- NULL
    genePredLst2 <- list()
    sumDf <- data.frame(matrix(nrow=0, ncol=11))
    colnames(sumDf) <- c("Histotype","Model", 
                         "Sensitivity", "Specificity", "Accuracy", "F1",
                         "CI_Lower","CI_Upper",
                         "AUC", "AUPRC",
                         "CpGs")
    tmpH <- names(genePredLst)[i]
    gDf <- genePredLst[[tmpH]]
    nSamp <- nrow(gDf)
    if(is.null(itMax)){
      itMax <- nSamp - 1
    }
    if(is.null(nSamp)){
      next()
    }else if(nSamp == 1){
      nIt <- 1 
    }else if(nSamp < itMax){
      nIt <- nSamp-1
    }else{
      nIt <- itMax-1
    }
    for(j in 1:nIt){
      message(paste("Running iteration: [", j,"/", nIt,"] of the model combination", sep=""))
      genePredDf <- data.frame(matrix(nrow=0, ncol=11))
      colnames(genePredDf) <- c("Histotype", "Model", 
                                 "Sensitivity", "Specificity", "Accuracy", "F1",
                                 "CI_Lower","CI_Upper",
                                 "AUC", "AUPRC",
                                 "CpGs")
      # Select best model best on F1 score, if tied base it on AUPRC for model
      if(is.null(bestRow)){
        bestRow <- gDf[which.max(gDf$F1),]
        bestInd <- which.max(gDf$F1)
        if(nrow(bestRow) > 1){
          bestRow <- bestRow[which.max(bestRow$AUPRC),]
          bestInd <- bestInd[-which.min(bestRow$AUPRC)]
        }
        ens1 <- bestRow$ensembl_gene_id
        if(nIt == 1){
          tmpDf <- gDf
        }else{
          tmpDf <- gDf[-bestInd, ]
        }
        # Create new model
        pred1 <- bestRow$CpGs
        pred1 <- strsplit(pred1, ",")[[1]]
        pred1str <- paste(pred1, collapse = " + ")
      }else{
        bestRow <- sumDf[nrow(sumDf),]
        ens1 <- bestRow$Model
        ens1 <- strsplit(ens1, "\\+")[[1]]
        tmpDf <- gDf[-which(gDf$ensembl_gene_id %in% ens1), ]
        # Create new model
        bPreds <- gDf$CpGs[which(gDf$ensembl_gene_id %in% ens1)]
        bPreds <- paste(bPreds, collapse = ",")
        pred1 <- strsplit(bPreds, ",")[[1]]
        pred1str <- paste(pred1, collapse = " + ")
      }
      
      for(k in 1:nrow(tmpDf)){
        message(paste("Processing gene: ", k ," out of: ", nrow(tmpDf), sep=""))
        nRow <- tmpDf[k,]
        ens2 <- nRow$ensembl_gene_id
        pred2 <- nRow$CpGs
        pred2 <- strsplit(pred2, ",")[[1]]
        pred2str <- paste(pred2, collapse = " + ")
        nMod <- paste(pred1str, pred2str, sep = " + ")
        nMod <- reformulate(nMod, "Histotype")  
        modPreds <- c(pred1, pred2)
        # Try new model
        modDf <- refBeta[modPreds, ]
        modDf <- as.data.frame(t(modDf))
        modDf  <- modDf[rowSums(is.na(modDf)) == 0,]
        modDf$Histotype <- refPheno$Histotype[match(rownames(modDf), refPheno$barcode)]
        modDf$Histotype <- ifelse(modDf$Histotype %in% tmpH, 1, 0)
        modDf$Histotype <- factor(modDf$Histotype, levels = c(1, 0))
        # Shuffle rows to remove any bias when sampling
        modDf <- modDf[sample(nrow(modDf)),]
        # Perform xgboost on training cohort
        epochLst <- list()
        # Perform multicore classification
        epochLst <- future_lapply(1:noEpochs, future.seed=TRUE, function(l){
          folds <- caret::createFolds(modDf$Histotype, 
                                      k = nFold, 
                                      list = FALSE) 
          # predStatLst <- list()
          kfDf <- data.frame(matrix(nrow=0, ncol=8))
          colnames(kfDf) <- c("Sensitivity", "Specificity", "Accuracy", "F1",
                              "CI_Lower","CI_Upper", 
                              "AUC", "AUPRC")
          for(m in 1:nFold){
            # Segment train-data by folds, with respect to the fold index no. 
            testIndexes <- which(folds==m, arr.ind=TRUE)
            trainData <- modDf[-testIndexes, ]
            testData <- modDf[testIndexes, ]
            # Create weights for random forest classification (balance)
            trainVals <- trainData[,-which(grepl("Histotype",colnames(trainData)))]
            trainLabel <- trainData[,which(grepl("Histotype",colnames(trainData)))]
            trainLabel <- factor(trainLabel, levels = c(0,1))
            testVals <- testData[,-which(grepl("Histotype",colnames(testData)))]
            testLabel <- testData[,which(grepl("Histotype",colnames(testData)))]
            testLabel <- factor(testLabel, levels = c(0,1))
            # Get the ratio of negative (0) to positive (1) classes
            class_imba <- length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
            if(is.null(paramBool)){
              tmpBoost <- xgboost(x = as.matrix(trainVals),
                                  y = trainLabel,
                                  booster = "gbtree",
                                  objective = "binary:logistic",
                                  scale_pos_weight = class_imba,
                                  nrounds = 100,
                                  nthread = 10)
            }else{
              tmpBoost <- xgboost(data = as.matrix(trainVals),
                                  label = trainLabel,
                                  booster = "gbtree",
                                  objective = "binary:logistic",
                                  eta =  best_params["eta"] ,
                                  max_depth =  best_params["max_depth"],
                                  scale_pos_weight = class_imba,
                                  min_child_weight = best_params["min_child_weight"],
                                  subsample = best_params["subsample"],
                                  colsample_bytree = best_params["colsample_bytree"],
                                  nrounds = best_params["n_estimators"],
                                  verbose = 0,
                                  nthread = 10)
            }
            tmpExtPred <- predict(object=tmpBoost, 
                                  newdata = as.matrix(testVals),
                                  type="response")
            # Categorize predictions based on their corresponding classification seen to model
            conMat <- confusionMatrix(data = factor(as.numeric(tmpExtPred > 0.5), levels = c(1,0)), 
                                      reference = factor(testLabel, levels = c(1,0)))
            # Get AUC for model
            roc <- pROC::roc(testLabel, tmpExtPred, smooth=FALSE,  quiet = TRUE)
            # As we are dealing with heavily imbalanced data, optimize on f1 score and AUPRC
            sens <- conMat$byClass[[1]]
            spec <- conMat$byClass[[2]]
            f1Score <- conMat$byClass[[7]]
            acc <- conMat$overall[[1]]
            # Get confidence interval of statistics
            confIntLower <- conMat$overall[[3]]
            confIntUpper <- conMat$overall[[4]]
            # Get AUC for model
            roc <- pROC::roc(testLabel, tmpExtPred, smooth=FALSE,  quiet = TRUE)
            # Get AUPRC
            # Indice locations of pos and neg-classes
            indice_POS=which(testLabel %in% 1)
            indice_NEG=which(testLabel %in% 0)
            # Get scores of pos and neg classes
            clas_score_POS=tmpExtPred[indice_POS]
            clas_score_NEG=tmpExtPred[indice_NEG]
            # Calculate area under precision recall curve
            # Scores class 0 is equal to the positive class (i.e. histotype)
            auprc <-PRROC::pr.curve(scores.class1 = clas_score_NEG, 
                                    scores.class0 = clas_score_POS)
            
            # Save to dataframe
            kfDf[nrow(kfDf)+1,] <- c(sens, spec, acc, f1Score,
                                       confIntLower, confIntUpper,
                                       roc$auc[1], auprc$auc.integral[1])
            kfDf[sapply(kfDf, is.nan)] <- 0
            #kfDf[nrow(kfDf)+1, ] <- colMeans(kfDf2)
            #rownames(kfD)[nrow(kfDf)] <- "MeanVal"
            #rownames(kfDf)[nrow(kfDf)] <- paste("k=",p, sep="")
          }
          epochLst[[l]] <- colMeans(kfDf, na.rm = TRUE)
        })
        sumCvDf <- do.call(rbind.data.frame, epochLst)
        colnames(sumCvDf) <- c("Sensitivity", "Specificity", "Accuracy","F1", 
                                "CI_Lower", "CI_Upper", 
                                "AUC", "AUPRC")
        meanPred <- colMeans(sumCvDf, na.rm = TRUE)
        genePredDf[nrow(genePredDf)+1,] <- c(tmpH, paste(paste(ens1, collapse="+"), ens2, sep="+"), 
                                               meanPred[1], meanPred[2], meanPred[3], meanPred[4],
                                               meanPred[5],meanPred[6],
                                               meanPred[7],meanPred[8],
                                               paste(modPreds, collapse=","))
      }
      genePredDf <- genePredDf[order(genePredDf$F1, decreasing=TRUE), ]
      if(length(which.max(genePredDf$F1)) > 2){
        genePredDfAUPRC <- genePredDf[which.max(genePredDf$F1),]
        bestModel <- genePredDfAUPRC[which.max(genePredDfAUPRC$AUPRC),]
      }else{
        if(length(which(is.na(genePredDf$F1))) > 1){
          bestModel <- genePredDf[which.max(genePredDf$AUPRC),]
        }else{
          bestModel <- genePredDf[which.max(genePredDf$F1),]
        }
      }
      sumDf[nrow(sumDf)+1,] <- bestModel 
      bestRow <- sumDf[nrow(sumDf),]
    }
    modelLst[[i]] <- sumDf
    names(modelLst)[i] <- tmpH
  }
  message("Model combination analysis completed")
  
  # Combine all different CpG-sets, 
  # Use the most beneficial combination of promoters as a gene-panel for histotype stratification
  fullModel <- c()
  geneModel <- c()
  bModelDf <- data.frame(matrix(nrow=0, ncol=ncol(modelLst[[1]])))
  colnames(bModelDf) <- colnames(modelLst[[1]])
  for(r in 1:length(modelLst)){
    tmpModels <- modelLst[[r]]
    tmpH <- names(modelLst)[r]
    # tmpModels$SensSpec <- as.numeric(tmpModels$Sensitivity) + as.numeric(tmpModels$Specificity)
    bModel <- tmpModels[which.max(tmpModels$F1), ]
    if(nrow(bModel) >1){
      bModel <- bModel[which.max(tmpModels$AUPRC), ]
    }
    tmpGM <- bModel$Model
    tmpGM <- strsplit(tmpGM , "\\+")[[1]]
    names(tmpGM) <- paste(tmpH, 1:length(tmpGM), sep="")
    tmpFM <- bModel$CpGs
    tmpFM <- strsplit(tmpFM , ",")[[1]]
    #TODO: Change names to ens-id's over CpGs, and Histotypes over Ens-genes, alternatively call script for getting the model using the CpG's
    names(tmpFM) <-  paste(tmpH, 1:length(tmpFM), sep="")
    fullModel <- append(fullModel, tmpFM)
    geneModel <- append(geneModel, tmpGM)
    bModelDf[nrow(bModelDf)+1,] <- bModel 
  }
  bModelDf$Histotype <- names(modelLst)
  message("Combined best model created")
  outLst <- list("CpGModel" = fullModel, "ensModel" = geneModel, "statDf" = bModelDf, "modelLst" = modelLst)
  return(outLst)
}

makeStepwiseMultModel_MULT <- function(inpCpgDfLst, refBeta, refPheno, inpCpgPos, 
                                       inpCpgLocs=NULL, noEpochs = NULL, nFold = NULL, 
                                       paramBool = NULL, rmHist = NULL, minCpg = NULL, 
                                       itMax = NULL, fullPromoBool=NULL, nCatBool = NULL){
  # Driver function for selecting the best multiclass model using single genes
  if(is.null(nFold)){
    nFold <- 5
  }
  if(is.null(noEpochs)){
    noEpochs <- 10
  }
  if(is.null(minCpg)){
    minCpg <- 2
  }
  
  if(is.null(itMax)){
    itMax <- 5
  }
  nFold <- nFold
  if(!is.null(rmHist)){
    refPheno <- refPheno[!refPheno$Histotype %in% rmHist, ]
    refBeta <- refBeta[,refPheno$barcode]
  }
  genePredLst <- makeSingleGenePred_MULT(inpCpgDfLst = keepPromoCpgPos, 
                                         refBeta = refBeta, 
                                         refPheno = refPheno,
                                         inpCpgPos = inpCpgPos,
                                         paramBool = NULL,
                                         noEpochs = noEpochs,
                                         minCpg = minCpg,
                                         nFold = nFold)
  bestModel <- makeGenePredLstComb_MULT(genePredLst = genePredLst, 
                                        refBeta = refBeta,
                                        refPheno = refPheno, 
                                        noEpochs = 50, 
                                        paramBool = NULL, 
                                        itMax = itMax,
                                        nFold = nFold)
  return(bestModel)
}

makeGridBestModel <- function(inpModelLst, inpBeta, inpPheno, sigCpgLocs){
  # allModelDf <- do.call(rbind.data.frame, bestModel$modelLst)
  allModelDf <- do.call(rbind.data.frame, 
                        inpModelLst)
  allModelGrid <- expand.grid(x = allModelDf$CpGs[which(allModelDf$Histotype %in% "CCC")],
                              y = allModelDf$CpGs[which(allModelDf$Histotype %in% "HGSC")],
                              z = allModelDf$CpGs[which(allModelDf$Histotype %in% "MC")])
  # Create dataframe for storing data
  sumCvDfMult <- data.frame(matrix(nrow=0, ncol=7+(4*(length(table(trainP$Histotype))))))
  colnames(sumCvDfMult) <- c(paste(names(table(trainP$Histotype)), "Sensitivity",sep="_"), 
                             paste(names(table(trainP$Histotype)), "Specificity",sep="_"),
                             paste(names(table(trainP$Histotype)), "Accuracy",sep="_"),
                             paste(names(table(trainP$Histotype)), "F1Score",sep="_"),
                             "Acc_Overall","Acc_Lower", "Acc_Upper",
                             "AUC", "AUPRC", "Model", "Fold")
  
  modelLst <- list()
  for(i in 1:nrow(allModelGrid)){
    message(paste("Running model: ", i, " out of total models: ", nrow(allModelGrid), sep=""))
    cccMod <- strsplit(as.character(allModelGrid$x[i]), ",")[[1]]
    names(cccMod) <- paste("CCC", 1:length(cccMod), sep="") 
    hgscMod <- strsplit(as.character(allModelGrid$y[i]), ",")[[1]]
    names(hgscMod) <- paste("HGSC", 1:length(hgscMod), sep="") 
    mcMod <- strsplit(as.character(allModelGrid$z[i]), ",")[[1]]
    names(mcMod) <- paste("MC", 1:length(mcMod), sep="") 
    tmpMod <- c(cccMod, hgscMod, mcMod)
    # Set up starting data for model
    cpgModel <- tmpMod
    startData <- inpBeta[cpgModel,]
    startData <- data.frame(t(startData))
    trainP <- refPheno
    annMod <- getAnnotatedCpgModel(inpCpg=cpgModel, 
                                     inpPromoLocs = sigCpgLocs)
    modelLst[[i]] <- annMod$CpG
    folds <- caret::createFolds(trainP$Histotype, 
                                k = nFold, 
                                list = FALSE) 
    for(k in 1:nFold){
      testIndexes <- which(folds == k, 
                           arr.ind=TRUE)
      trainData <- startData[-testIndexes,]
      testData <- startData[testIndexes,]
      trainData$Histotype <- trainP$Histotype[match(rownames(trainData), trainP$barcode)]
      testData$Histotype <- trainP$Histotype[match(rownames(testData), trainP$barcode)]
      trainVals <- trainData[,-which(colnames(trainData) %in% "Histotype")]
      testVals <- testData[,-which(colnames(testData) %in% "Histotype")]
      
      # Create phenotypic annotation
      trainLabel <- trainData$Histotype
      testLabel <- testData$Histotype
      
      repVec <- c()
      for(v in 1:length(table(trainLabel))){
        repVec <- append(repVec, v-1)
        names(repVec)[length(repVec)] <- names(table(trainLabel))[v]
        trainLabel <- replace(trainLabel, trainLabel==names(table(trainLabel))[[v]], v-1)
      }
      # Replace values in testlabel using repvec
      for(w in 1:length(repVec)){
        testLabel <- ifelse(testLabel %in% names(repVec)[w], repVec[w], testLabel)
      }
      # Ratio of classes in multiclass
      class_imba <- table(trainLabel)/length(trainLabel)
      cW <- min(class_imba)/class_imba
      cWArr <- rep(NA, length(trainLabel))
      for(x in 1:length(trainLabel)){
        cWArr[x] <- cW[which(names(cW) %in% trainLabel[x])]
      }
      tmpBoost <- xgboost(x = as.matrix(trainVals),
                                y = trainLabel,
                                booster = "gbtree",
                                objective = "multi:softprob",
                                weights = cWArr,
                                nrounds =100,
                                nthread = 10)

      tmpExtPred <- predict(object=tmpBoost, 
                                newdata = as.matrix(testVals))
      colnames(tmpExtPred) <- names(table(trainLabel))
          
      # May not be needed
      predicted_labels <- factor(colnames(tmpExtPred)[max.col(tmpExtPred)], levels=repVec)
      testLabelFac <- factor(testLabel, levels=repVec)
      conMat <- caret::confusionMatrix(data = predicted_labels, 
                                           reference = testLabelFac)
      # Correction for instances where we only have 2 classes 
      if(length(levels(testLabelFac))<=2){
        sens <- conMat$byClass[1]
        sens <- append(sens, sens)
        names(sens) <- names(table(testLabel))
        spec <- conMat$byClass[2]
        spec <- append(spec, spec)
        names(spec) <- names(table(testLabel))
        f1Score <- conMat$byClass[7]
        f1Score  <- append(f1Score, f1Score)
        names(f1Score) <- names(table(testLabel))
        acc <- conMat$byClass[11]
        acc  <- append(acc, acc)
        names(acc) <- names(table(testLabel))
      }else{
        sens <- conMat$byClass[,1]
        spec <- conMat$byClass[,2]
        f1Score <- conMat$byClass[,7]
        acc <- conMat$byClass[,11]
      }
      # Get overall statistics
      accAll <- conMat$overall[[1]]
      # Get confidence interval of accuracy (lower, upper)
      accLower <- conMat$overall[[3]]
      accUpper <- conMat$overall[[4]]
          
      # Get AUC for model
      rocPred <- tmpExtPred
      colnames(rocPred) <- repVec
      # Get AUC for model
      # If we dont have more then 1 class, we cant get these values and instead assign them to 0 
      if(length(table(testLabel)) > 1){
        # Try testlabelfac instead for the instances with no or few classes 
        roc <- pROC::multiclass.roc(predictor = rocPred,
                                        response= testLabel, 
                                        smooth=FALSE,  
                                        quiet = TRUE)
        aucVal <- roc$auc[1]
            
        # Get PRAUC for model
        prAucDf <- as.data.frame(tmpExtPred)
        prLevFac <- droplevels(testLabelFac)
        prAucDf <- prAucDf[,colnames(prAucDf) %in% levels(prLevFac)]
        prAucDf$ref <-  prLevFac
        prauc <- prAucDf %>% yardstick::pr_auc(estimator = "macro",
                                                   ref,
                                                   min(levels(prLevFac)):max(levels(prLevFac)))
        prAucVal <- prauc$.estimate[1]
      }else{
        aucVal <- 0
        prAucVal <- 0
      }
      for(u in 1:length(names(table(trainP$Histotype)))){
        tmpH <- names(table(trainP$Histotype))[u]
        if(!tmpH %in% names(repVec)){
          sens <- append(sens, 0, after=u-1)
          names(sens)[u] <- tmpH
          spec <- append(spec, 0, after=u-1)
          names(spec)[u] <- tmpH
          f1Score <- append(f1Score, 0, after=u-1)
          names(f1Score)[u] <- tmpH
          acc <- append(acc, 0, after=u-1)
          names(acc)[u] <- tmpH
        }else{
          names(sens)[u] <- tmpH
          names(spec)[u] <- tmpH
          names(acc)[u] <- tmpH
          names(f1Score)[u] <- tmpH
        }
      }
      sens <- sens[names(sens) %in% names(table(trainP$Histotype))]
      spec <- spec[names(spec) %in% names(table(trainP$Histotype))]
      acc <- sens[names(acc) %in% names(table(trainP$Histotype))]
      f1Score <- f1Score[names(f1Score) %in% names(table(trainP$Histotype))]
          
      sumCvDfMult[nrow(sumCvDfMult)+1,] <- c(sens, spec, acc, f1Score,
                                                 accAll, accLower, accUpper,
                                                 aucVal, prAucVal,
                                                 i, k)
          
    }
  }
  sumCvDfMult[sapply(sumCvDfMult, is.nan)] <- 0
  sumCvDfMult[sapply(sumCvDfMult, is.na)] <- 0
  meanItDf <- sumCvDfMult %>%
  group_by(Model) %>%
  summarise_all("mean")
  meanItDf <- data.frame(meanItDf)
  bestModel <- modelLst[[which.max(meanItDf$AUPRC)]]
  return(bestModel)
}
    
makeModelLstCombs <- function(modelLst, refBeta, refPheno, noEpochs = NULL, nFold = NULL, selCat = NULL){
  # Combine list of models in all possible configurations to see which best classifies the dataset
  if(is.null(nFold)){
    nFold <- 5
  }
  if(is.null(noEpochs)){
    noEpochs <- 25
  }
  if(is.null(selCat)){
    selCat <- "PRAUC" 
  }
  
  if(!"Histotype" %in% colnames(modelLst[[1]])){
    for(i in 1:length(modelLst)){
      modelLst[[i]]$Histotype <- names(modelLst)[i]
    }
  }
  
  # Merge all model-lists
  allMods <- Reduce(function(...) merge(..., all=T), modelLst)
  
  modLst <- list()
  for(x in 1:length(table(allMods$Histotype))){
    tmpH <- names(table(allMods$Histotype))[x]
    hMods <- allMods[allMods$Histotype %in% tmpH,]
    modLst[[length(modLst)+1]] <- hMods$CpGs
    names(modLst)[length(modLst)] <- tmpH
  }
  modCombs <- expand.grid(modLst)
  modMeanDf <- data.frame(matrix(nrow=0, ncol=7+(4*(length(table(refPheno$Histotype))))))
  colnames(modMeanDf) <- c(paste(names(table(refPheno$Histotype)), "Sensitivity",sep="_"), 
                           paste(names(table(refPheno$Histotype)), "Specificity",sep="_"),
                           paste(names(table(refPheno$Histotype)), "Accuracy",sep="_"),
                           paste(names(table(refPheno$Histotype)), "F1Score",sep="_"),
                           "Acc_Overall",
                           "CI_Lower",
                           "CI_Upper",
                           "AUC",
                           "PRAUC",
                           "Model_Row",
                           "CpGs")
  for(z in 1:nrow(modCombs)){
    tmpMod <- modCombs[z,]
    # Create model based on model-combinations
    modVec <- c()
    for(t in 1:ncol(tmpMod)){
      modVec <- append(modVec, sapply(as.character(tmpMod[,t]), function(x) strsplit(x, ","))[[1]])
    }
    # Prepare data for training the model
    # Choose between training the model on our cohort and applying it directly on ext,
    # Create dataframe for storing data
    modSumDf <- data.frame(matrix(nrow=0, ncol=7+(4*(length(table(refPheno$Histotype))))))
    colnames(modSumDf) <- c(paste(names(table(refPheno$Histotype)), "Sensitivity",sep="_"), 
                            paste(names(table(refPheno$Histotype)), "Specificity",sep="_"),
                            paste(names(table(refPheno$Histotype)), "Accuracy",sep="_"),
                            paste(names(table(refPheno$Histotype)), "F1Score",sep="_"),
                            "Acc_Overall",
                            "CI_Lower",
                            "CI_Upper",
                            "AUC",
                            "PRAUC",
                            "Epoch", 
                            "Fold")
    for(e in 1:noEpochs){
      folds <- caret::createFolds(factor(refPheno$Histotype), k = nFold, list = FALSE)
      for(f in 1:nFold){
        multModBeta <- refBeta[modVec,]
        # Segment train-data by folds, with respect to the fold index no. 
        testIndexes <- which(folds==f, arr.ind=TRUE)
        multModBeta <- as.data.frame(t(multModBeta))
        trainData <- multModBeta[-testIndexes,]
        testData <- multModBeta[testIndexes,]
        # Add Histotype label to DF
        trainData$Histotype <- refPheno$Histotype[match(rownames(trainData), refPheno$barcode)]
        testData$Histotype <- refPheno$Histotype[match(rownames(testData), refPheno$barcode)]
        trainVals <- trainData[,-which(colnames(trainData) %in% "Histotype")]
        testVals <- testData[,-which(colnames(testData) %in% "Histotype")]
        # Create phenotypic annotation
        trainLabel <- trainData$Histotype
        testLabel <- testData$Histotype
        repVec <- c()
        for(v in 1:length(table(trainLabel))){
          repVec <- append(repVec, v-1)
          names(repVec)[length(repVec)] <- names(table(trainLabel))[v]
          trainLabel <- replace(trainLabel, trainLabel==names(table(trainLabel))[[v]], v-1)
        }
        # Replace values in testlabel using repvec
        for(w in 1:length(repVec)){
          testLabel <- ifelse(testLabel %in% names(repVec)[w], repVec[w], testLabel)
        }
        # Ratio of classes in multiclass
        class_imba <- table(trainLabel)/length(trainLabel)
        cW <- min(class_imba)/class_imba
        cWArr <- rep(NA, length(trainLabel))
        for(x in 1:length(trainLabel)){
          cWArr[x] <- cW[which(names(cW) %in% trainLabel[x])]
        }
        
        # Perform predictive classification
        tmpBoost <- xgboost(data = as.matrix(trainVals),
                            label = trainLabel,
                            booster = "gbtree",
                            objective = "multi:softprob",
                            weight = cWArr,
                            verbose = 0,
                            nthread = 10,
                            nrounds = 100,
                            num_class = length(table(trainLabel)))
        
        # Apply RF model on test-data (external dataset)
        tmpExtPred <- parsnip::xgb_predict(object=tmpBoost, 
                                           new_data = as.matrix(testVals))
        colnames(tmpExtPred) <- names(table(trainLabel))
        # May not be needed
        predicted_labels <- factor(colnames(tmpExtPred)[max.col(tmpExtPred)],levels=repVec)
        testLabelFac <- factor(testLabel, levels=repVec)
        conMat <- caret::confusionMatrix(data = predicted_labels, 
                                         reference = testLabelFac)
        # If we accidentally get binary dueto lack of classes for multiclass
        if(length(levels(testLabelFac))<=2){
          sens <- conMat$byClass[1]
          sens <- append(sens, sens)
          names(sens) <- names(table(testLabel))
          spec <- conMat$byClass[2]
          spec <- append(spec, spec)
          names(spec) <- names(table(testLabel))
          f1Score <- conMat$byClass[7]
          f1Score  <- append(f1Score, f1Score)
          names(f1Score) <- names(table(testLabel))
          acc <- conMat$byClass[11]
          acc  <- append(acc, acc)
          names(acc) <- names(table(testLabel))
        }else{
          sens <- conMat$byClass[,1]
          spec <- conMat$byClass[,2]
          f1Score <- conMat$byClass[,7]
          acc <- conMat$byClass[,11]
        }
        # Get overall statistics
        accAll <- conMat$overall[[1]]
        # Get confidence interval of accuracy (lower, upper)
        accLower <- conMat$overall[[3]]
        accUpper <- conMat$overall[[4]]
        
        # Get AUC for model
        rocPred <- tmpExtPred
        colnames(rocPred) <- names(table(testLabel))
        # Get AUC for model
        roc <- pROC::multiclass.roc(predictor = rocPred,
                                    response= testLabel, 
                                    smooth=FALSE,  
                                    quiet = TRUE)
        # Get PRAUC for model
        prAucDf <- as.data.frame(tmpExtPred)
        prAucDf$ref <- testLabelFac
        prauc <- prAucDf %>% yardstick::pr_auc(estimator = "macro",
                                               ref,
                                               min(testLabel):max(testLabel))
        for(u in 1:length(names(table(refPheno$Histotype)))){
          tmpH <- names(table(refPheno$Histotype))[u]
          if(!tmpH %in% names(repVec)){
            sens <- append(sens, 0, after=u-1)
            names(sens)[u] <- tmpH
            spec <- append(spec, 0, after=u-1)
            names(spec)[u] <- tmpH
            f1Score <- append(f1Score, 0, after=u-1)
            names(f1Score)[u] <- tmpH
            acc <- append(acc, 0, after=u-1)
            names(acc)[u] <- tmpH
          }else{
            names(sens)[u] <- tmpH
            names(spec)[u] <- tmpH
            names(acc)[u] <- tmpH
            names(f1Score)[u] <- tmpH
          }
        }
        modSumDf[nrow(modSumDf)+1,] <- c(sens, spec, acc, f1Score,
                                         accAll, accLower, accUpper,
                                         roc$auc[1], prauc$.estimate[1],
                                         e,f)
        
      }
    }
    modSumDf[sapply(modSumDf, is.nan)] <- 0
    modSumDf[sapply(modSumDf, is.na)] <- 0
    meanItDf <- modSumDf %>%
      group_by(Epoch) %>%
      summarise_all("mean")
    meanItDf <- data.frame(meanItDf)
    meanItDf <- meanItDf[!colnames(meanItDf) %in% c("Epoch", "Fold")]
    modelRow <- colMeans(meanItDf)
    modelRow[length(modelRow)+1] <- z
    names(modelRow)[length(modelRow)] <- "Model_Row"
    modelRow[length(modelRow)+1] <- paste(modVec, collapse=",")
    names(modelRow)[length(modelRow)] <- "CpGs"
    modMeanDf[nrow(modMeanDf)+1,] <- modelRow
  }
  # Select best model based on F1 score for categories we are investigating, and PR AUC
  if(selCat == "PRAUC"){
    bestMod <- modMeanDf[which.max(modMeanDf$PRAUC),]
    if(nrow(bestMod) > 0){
      sigCats <- names(modelLst)
      f1Cols <- paste(sigCats, "F1Score", sep="_")
      allF1 <- bestMod[, f1Cols]
      allF1 <- rowSums(mutate_all(allF1, function(x) as.numeric(as.character(x))))
      bestMod  <- bestMod[which.max(allF1 ),]
    }
  }else if(selCat == "F1"){
    sigCats <- names(modelLst)
    f1Cols <- paste(sigCats, "F1Score", sep="_")
    allF1 <- modMeanDf[, f1Cols]
    allF1 <- rowSums(mutate_all(allF1, function(x) as.numeric(as.character(x))))
    bestMod  <- modMeanDf[which.max(allF1 ),]
  }else{
    "Invalid category, please choose from F1 or PRAUC"
  }
  return(list("ModelCombs" = modCombs, "StatDf" = modMeanDf, "Model"=bestMod$CpGs))
}

################################################################################
# Predictive classification scripts
################################################################################

makeExtPredMult <- function(inpCpgDfLst, inpBeta, refBeta, inpPheno, refPheno, inpSid, inpGeneInf, noEpochs = NULL, nFold = NULL, paramBool = NULL, extBool = NULL){
  # Create/Select model in source dataset, apply on external dataset
  if(is.null(nFold)){
    nFold <- 5
  }
  if(is.null(noEpochs)){
    noEpochs <- 25
  }
  # Make sure that samples in the beta dataframe are present in the phenotypic data
  refPheno <- refPheno[which(refPheno$barcode %in% colnames(refBeta)),]
  refBeta <- refBeta[,which(colnames(refBeta) %in% refPheno$barcode)]
  inpBeta <- inpBeta[,which(colnames(inpBeta) %in% inpPheno$barcode)]
  inpPheno <- inpPheno[which(inpPheno$barcode %in% colnames(inpBeta)),]
  # Get CpG sites which are shared between the two input arrays
  sharedCpgs <- intersect(rownames(inpBeta), rownames(refBeta))
  # Make sure to only include CpG sites present in both dataframes (i.e. drop EPIC only probes to make dataframes compatible)
  refBeta <- refBeta[sharedCpgs,]
  inpBeta <- inpBeta[sharedCpgs,]
  # Create dataframes for the one vs. all approach (for grid-search and feature selection)
  refPredLst <- list()
  extPredLst <- list()
  manFiltLst <- list()
  for(i in 1:length(inpCpgDfLst)){
    inpCpgDf <- inpCpgDfLst[[i]]
    if(nrow(inpCpgDf) == 0){
      next()
    }
    inpHist <- names(inpCpgDfLst)[i]
    # Prepare reference dataframe
    refPredDf <- makeCpgDf(inpCpgDf, refBeta, refPheno, inpHist)
    refPredDf <- as.data.frame(t(refPredDf))
    refPredDf[is.na(refPredDf)] <- 0
    # refMat <- round(refMat, 5)
    refPredDf$Histotype <- refPheno$Histotype[match(rownames(refPredDf), refPheno$barcode)]
    refPredDf$Histotype <- ifelse(refPredDf$Histotype %in% inpHist, inpHist, "Other")
    # if(refPreds$Histotype)
    phenoBetas <- refPredDf[which(refPredDf$Histotype %in% inpHist),]
    phenoBetas <- phenoBetas[,-which(colnames(phenoBetas) %in% "Histotype")]
    phenoBetas <- t(phenoBetas)
    colnames(phenoBetas) <- refPheno$Sample_ID[match(colnames(phenoBetas), refPheno$barcode)]
    nPhenoBetas <- refPredDf[-which(refPredDf$Histotype %in% inpHist),]
    nPhenoBetas <- nPhenoBetas[,-which(colnames(nPhenoBetas) %in% "Histotype")]
    nPhenoBetas <- t(nPhenoBetas)
    colnames(nPhenoBetas) <- refPheno$Sample_ID[match(colnames(nPhenoBetas), refPheno$barcode)]
    # Calculate the mean* distance for predictors 
    # * TriMean or Geometric mean
    tmpHistotypes <- names(table(refPheno$Histotype))[!names(table(refPheno$Histotype)) %in% inpHist]
    betaDists <- makeBetaGrpDists(phenoBetas, nPhenoBetas, refPheno,  tmpHistotypes, inpHist, histBool=TRUE, distType="TRIMEAN")
    betaDists <- betaDists[order(rowSums(betaDists), decreasing = TRUE), ]
    # Filter CpG sites passing delta beta >= 0.2 for 2,3 histotypes for predictors
    minDists <- betaDists[which(rowSums(betaDists >0.2) >= 2),]
    filtDists <- minDists[which(rowSums(minDists >0.2) >= 3),]
    manFiltLst[[i]] <- rownames(filtDists)
    names(manFiltLst)[i] <- inpHist
    # Prepare training-data  based on the less stringent filtering
    refPredDf <- refBeta[rownames(betaDists), ]
    refPredDf <- as.data.frame(t(refPredDf))
    refPredDf  <- refPredDf[rowSums(is.na(refPredDf)) == 0,]
    refPredDf$Histotype <- refPheno$Histotype[match(rownames(refPredDf), refPheno$barcode)]
    refPredDf$Histotype <- ifelse(refPredDf$Histotype %in% inpHist, inpHist, "Other")
    refPredDf$Histotype <- factor(refPredDf$Histotype, levels = c("Other", inpHist))
    # Do the same for external cohort
    extPredDf <- inpBeta[rownames(minDists ), ]
    extPredDf <- as.data.frame(t(extPredDf))
    extPredDf  <- extPredDf[rowSums(is.na(extPredDf)) == 0,]
    extPredDf$Histotype <- inpPheno$Histotype[match(rownames(extPredDf), inpPheno$barcode)]
    extPredDf$Histotype <- ifelse(extPredDf$Histotype %in% inpHist, inpHist, "Other")
    extPredDf$Histotype <- factor(extPredDf$Histotype, levels = c("Other", inpHist))
    # Save each processed DF to list
    refPredLst[[length(refPredLst)+1]] <- refPredDf
    names(refPredLst)[length(refPredLst)] <- inpHist
    extPredLst[[length(extPredLst)+1]] <- extPredDf
    names(extPredLst)[length(extPredLst)] <- inpHist
  }
  allManPreds <- unlist(manFiltLst)
  allManPreds  <- allManPreds[!duplicated(allManPreds )]
  
  # Perform grid-search and feature selection using the one vs all approach for each histotype
  impPredLst <- list()
  borPredLst <- list()
  for(j in 1:length(refPredLst)){
    inpRef <- refPredLst[[j]]
    tmpHist <- names(refPredLst)[j]
    #inpExt <- extPredLst[[tmpHist]]
    # Shuffle rows to remove any bias when sampling
    inpRef <- inpRef[sample(nrow(inpRef)),]
    #inpExt <- inpExt[sample(nrow(inpExt)),]
    # Perform grid-search to identify best parameters for model
    fsTrain <- inpRef[,-which(grepl("Histotype",colnames(inpRef)))]
    fsLabel <- inpRef[,which(grepl("Histotype",colnames(inpRef)))]
    fsLabel <- ifelse(fsLabel == 'Other', 0, 1)
    # Ratio of negative to positive classes
    fs_imba <- length(which(fsLabel %in% 0))/length(which(fsLabel %in% 1))
    if(is.null(paramBool)){
      best_params <- NULL
      boostFS <- xgboost(data = as.matrix(fsTrain),
                         label = fsLabel,
                         scale_pos_weight = fs_imba,
                         verbose = 0,
                         nrounds = 100,
                         nthread = 10)
    }else{
      best_params <- makeXGBoostCaretGridS(inpRef)
      boostFS <- xgboost(data = as.matrix(fsTrain),
                         label = fsLabel,
                         booster = "gbtree",
                         objective = "binary:logistic",
                         eta =  best_params["eta"] ,
                         max_depth =  best_params["max_depth"],
                         scale_pos_weight = fs_imba,
                         min_child_weight = best_params["min_child_weight"],
                         subsample = best_params["subsample"],
                         colsample_bytree = best_params["colsample_bytree"],
                         nrounds = best_params["n_estimators"],
                         verbose = 0,
                         nthread = 10)
    }
    
    impPreds <- NULL
    borPreds <- NULL
    importance_matrix = xgb.importance(colnames(as.matrix(fsTrain)), model = boostFS)
    impPreds <- importance_matrix$Feature[which(importance_matrix$Gain > 0.01)]
    if(length(impPreds) < 5){
      if(length(importance_matrix$Feature) < 5){
        ifLim <- length(importance_matrix$Feature) 
      }else{
        ifLim <- 5
      }
      impPreds <- importance_matrix$Feature[1:ifLim]
    }
    if(is.null(paramBool)){
      borutaRes = Boruta::Boruta(as.matrix(fsTrain),
                                 y= fsLabel,
                                 maxRuns=100, 
                                 doTrace=2,
                                 nrounds=100,
                                 holdHistory=TRUE,
                                 getImp=getImpXgboost,
                                 scale_pos_weight= fs_imba, 
                                 eval_metric="auc", 
                                 eval_metric="rmse", 
                                 eval_metric="logloss",
                                 objective="binary:logistic",
                                 tree_method="hist",
                                 nthread=10)
    }else{
      borutaRes = Boruta::Boruta(as.matrix(fsTrain),
                                 y= fsLabel,
                                 maxRuns=100, 
                                 doTrace=2,
                                 holdHistory=TRUE,
                                 getImp=getImpXgboost,
                                 max.depth=best_params["max_depth"], 
                                 eta= best_params["eta"], 
                                 min_child_weight=best_params["min_child_weight"],
                                 scale_pos_weight= fs_imba, 
                                 subsample = best_params["subsample"],
                                 colsample_bytree = best_params["colsample_bytree"],
                                 nrounds= best_params["n_estimators"], 
                                 eval_metric="auc", 
                                 eval_metric="rmse", 
                                 eval_metric="logloss",
                                 objective="binary:logistic",
                                 tree_method="hist",
                                 nthread=10)
    }
    #extract Boruta's decision
    boruta_dec = attStats(borutaRes)
    borPreds <- rownames(boruta_dec[!boruta_dec$decision %in% c("Rejected"),])
    if(length(borPreds) < 5){
      if(nrow(boruta_dec) < 5){
        bcLim <- nrow(boruta_dec)
      }else{
        bcLim <- 5
      }
      boruta_dec <- boruta_dec[order(boruta_dec$meanImp, decreasing = TRUE),] 
      borPreds <- rownames(boruta_dec[1:bcLim,])
    }
    impPredLst[[j]] <- impPreds
    names(impPredLst)[j] <- tmpHist
    borPredLst[[j]] <- borPreds
    names(borPredLst)[j] <- tmpHist
  }
  allBorPreds <- unlist(borPredLst)
  allImpPreds <- unlist(impPredLst)
  predLst <- list("BORUTA"= allBorPreds, "IMPORTANCE"=allImpPreds, "MANUAL" = allManPreds)
  if(is.null(best_params)){
    best_params <- c("eta"= 0.15, "n_estimators" = 500, 
                     "subsample" = 0.7, "colsample_bytree"= 0.6,
                     "max_depth" = 6, "min_child_weight"=1)
  }
  
  predTypeLst <- list()
  for(k in 1:length(predLst)){
    kfDfLst <- list()
    tmpPreds <- predLst[[k]]
    # Prepare reference dataframe
    refDf <- data.frame(t(refBeta[tmpPreds ,]))
    refDf$Histotype <- refPheno$Histotype[match(rownames(refDf), refPheno$barcode)]
    refDf$Histotype <- factor(refDf$Histotype, levels = names(table(refPheno$Histotype)))
    if(!is.null(paramBool)){
      best_params <- makeXGBoostCaretGridSMulti(refDf)
    }
    # Repeat for external dataframe
    extDf <- data.frame(t(inpBeta[tmpPreds,]))
    extDf$Histotype <- inpPheno$Histotype[match(rownames(extDf), inpPheno$barcode)]
    extDf$Histotype <- factor(extDf$Histotype, levels = names(table(inpPheno$Histotype)))     
    
    sumCvDf <- data.frame(matrix(nrow=0, ncol=3+(2*(length(table(refDf$Histotype))))))
    colnames(sumCvDf) <- c(paste(names(table(refDf$Histotype)), "Sensitivity",sep="_"), 
                           paste(names(table(refDf$Histotype)), "Specificity",sep="_"),
                           "Overall_Accuracy","Accuracy_p", "Overall_AUC")
    
    for(l in 1:noEpochs){
      folds <- caret::createFolds(factor(extDf$Histotype), k = nFold, list = FALSE) 
      # predStatLst <- list()
      kfDf <- data.frame(matrix(nrow=0, ncol=3+(2*(length(table(refDf$Histotype))))))
      colnames(kfDf) <- c(paste(names(table(refDf$Histotype)), "Sensitivity",sep="_"), 
                          paste(names(table(refDf$Histotype)), "Specificity",sep="_"),
                          "Overall_Accuracy","Accuracy_p", "Overall_AUC")
      for(m in 1:nFold){
        # Segment train-data by folds, with respect to the fold index no. 
        testIndexes <- which(folds==m,arr.ind=TRUE)
        testData <- extDf[testIndexes, ]
        if(!is.null(extBool)){
          trainData <- extDf[-testIndexes, ]
        }else{
          trainData <- refDf
          trainData <- trainData[trainData$Histotype %in% extDf$Histotype, ]
          trainData$Histotype <- droplevels(trainData$Histotype)
        }
        
        # Create weights for random forest classification (balance)
        trainTrain <- trainData[,-which(grepl("Histotype",colnames(trainData)))]
        trainLabel <- trainData[,which(grepl("Histotype",colnames(trainData)))]
        trainLabel <- as.character(trainLabel)
        for(n in 1:length(table(trainLabel))){
          trainLabel <- replace(trainLabel, trainLabel==names(table(trainLabel))[[n]], n-1)
        }
        trainLabel <- as.numeric(trainLabel)
        # trainLabel <- as.factor(trainLabel)
        
        # Ratio of negative to positive classes
        # class_imba <- table(trainLabel)/length(trainLabel)
        # cW <- min(class_imba)/class_imba
        # cWArr <- rep(NA, length(trainLabel))
        # for(p in 1:length(trainLabel)){
        #   cWArr[p] <- cW[which(names(cW) %in% trainLabel[p])]
        # }
        
        # best_params <- makeXGBoostGridS(trainData)
        # Perform grid-search on training data to retrieve optimal model parameters
        # best_params <- makeXGBoostGridS(trainData)
        # Build gradient boosting forest model
        if(is.null(paramBool)){
          tmpBoost <- xgboost(data = as.matrix(trainTrain),
                              label = trainLabel,
                              booster = "gbtree",
                              objective = "multi:softprob",
                              #weight = cWArr,
                              verbose = 0,
                              nthread = 10,
                              nrounds = 100,
                              num_class = length(table(extDf$Histotype)))
        }else{
          tmpBoost <- xgboost(data = as.matrix(trainTrain),
                              label = trainLabel,
                              booster = "gbtree",
                              objective = "multi:softprob",
                              eta =  best_params["eta"] ,
                              max_depth =  best_params["max_depth"],
                              #weight = cWArr,
                              min_child_weight = best_params["min_child_weight"],
                              subsample = best_params["subsample"],
                              colsample_bytree = best_params["colsample_bytree"],
                              nrounds = best_params["n_estimators"],
                              verbose = 0,
                              nthread = 10,
                              num_class = length(table(extDf$Histotype)))
        }
        
        # Apply RF model on test-data (external dataset)
        tmpExtPred <- parsnip::xgb_predict(object=tmpBoost, 
                                           new_data = as.matrix(testData[,-which(grepl("Histotype",colnames(testData)))]))
        testLabel <- as.character(testData$Histotype)
        for(o in 1:length(table(refDf$Histotype))){
          testLabel <- replace(testLabel, testLabel==names(table(refDf$Histotype))[[o]], o-1)
        }
        testLabel <-as.numeric(testLabel)
        rocPred <- tmpExtPred
        colnames(rocPred) <- names(table(trainLabel))
        # Get AUC for model
        roc <- pROC::multiclass.roc(response=testLabel, predictor = rocPred, smooth=FALSE,  quiet = TRUE)
        
        colnames(tmpExtPred) <- names(table(trainLabel))
        predicted_labels=factor(colnames(tmpExtPred)[max.col(tmpExtPred)],levels=0:(length(table(trainLabel))-1))
        testLabel <- factor(testLabel, levels = levels(predicted_labels))
        conMat <- caret::confusionMatrix(testLabel,predicted_labels)
        conAll <- conMat$overall
        conClass <- data.frame(conMat$byClass)
        rownames(conClass) <- names(table(trainData$Histotype))
        sensVec <- conClass$Sensitivity
        specVec <- conClass$Specificity
        for(n in 1:length(names(table(refDf$Histotype)))){
          tmpH <- names(table(refDf$Histotype))[n]
          if(!tmpH %in% rownames(conClass)){
            sensVec <- append(sensVec, 0, after=n-1)
            specVec <- append(specVec, 0, after=n-1)
          }
        }
        kfDf[nrow(kfDf)+1,] <- c(sensVec, specVec, conAll[[1]], conAll[[6]], roc$auc[1])
        rownames(kfDf)[nrow(kfDf)] <- paste("k=",m, sep="")
      }
      kfDf[sapply(kfDf, is.nan)] <- 0
      kfDf[sapply(kfDf, is.na)] <- 0
      kfDf[nrow(kfDf)+1, ] <- colMeans(kfDf)
      rownames(kfDf)[nrow(kfDf)] <- "MeanVal"
      kfDfLst[[l]] <- kfDf
      names(kfDfLst)[l] <- l
      sumCvDf[l,] <- kfDf["MeanVal",]
      rownames(sumCvDf)[l] <- l
    }
    sumCvDf[nrow(sumCvDf)+1, ] <- colMeans(sumCvDf)
    rownames(sumCvDf)[nrow(sumCvDf)] <- "MeanVal"
    predTypeLst[[k]] <- sumCvDf
    names(predTypeLst)[k] <- names(predLst)[k] 
  }
  return(list("BORUTA"= predTypeLst$BORUTA,
              "IMPORTANCE" = predTypeLst$IMPORTANCE,
              "MANUAL" = predTypeLst$MANUAL,
              "BORUTA_PREDS" = allBorPreds,
              "IMPORTANCE_PREDS" = allImpPreds,
              "MANUAL_PREDS" = allManPreds))
}


makeModelExtPredMult <- function(inpExtBetaLst, inpExtPhenoLst, inpCModel, refPheno, sigCpgLocs ,refBeta = NULL, extBool = NULL, noEpochs = NULL, nFold = NULL, paramBool = NULL, remHistBool=NULL){
  # Function for carrying out predictive classification on external dataset after initial model selection on source dataset 
  # Given an initial model
  if(is.null(refBeta)){
    extBool <- NULL
  }
  if(is.null(noEpochs)){
    noEpochs <- 100
  }
  if(is.null(nFold)){
    nFold <- 5
  }
  # cpgModel <-  sapply(inpCModel, function(x) strsplit(x, "\\+")[[1]][1], USE.NAMES=FALSE)
  # ensModel <- sapply(inpEModel, function(x) strsplit(x, "\\+")[[1]][1], USE.NAMES=FALSE)
  # Run multiclass prediction on test-datasets
  extPredLst <- list()
  cvMultLst <- list()
  cvMultMeanLst <- list()
  cpgModVec <- c()
  annCpgModVec <- c()
  ensModVec <- c()
  missPredVec <- c()
  hCatVec <- c()
  for(s in 1:length(inpExtBetaLst)){
    cpgModel <- inpCModel
    tmpSid <- names(inpExtBetaLst)[s]
    message(paste("Performing multiclass classification on cohort: ", tmpSid, sep=""))
    tmpExtB <- inpExtBetaLst[[tmpSid]]
    
    # Control what CpGs are in the dataframe compared to the model
    if(length(which(!cpgModel %in% rownames(tmpExtB))) > 1){
      message(paste("Warning, CpGs: ", paste(cpgModel[which(!cpgModel %in% rownames(tmpExtB))], collapse=", "), 
                    " not found in beta-matrix, and will be excluded from analysis", sep=""))
      remCpg <- cpgModel[which(!cpgModel %in% rownames(tmpExtB))]
      missPredVec <- append(missPredVec, paste(remCpg, collapse="+"))
      cpgModel <- cpgModel[!cpgModel %in% remCpg]
    }else{
      missPredVec <- append(missPredVec, 0)
    }
    names(missPredVec)[length(missPredVec)] <- tmpSid
    tmpExtB <- tmpExtB[cpgModel,]
    
    # Save model for later use
    cpgModVec <- append(cpgModVec, paste(cpgModel, collapse="+"))
    names(cpgModVec)[length(cpgModVec)] <- tmpSid
    # Annotate model
    annMod <- getAnnotatedCpgModel(inpCpg= cpgModel, 
                                   inpPromoLocs = sigCpgLocs)
    # Save annotated model
    annCpgModVec <- append(annCpgModVec, paste(annMod$CpG, collapse="+"))
    names(annCpgModVec)[length(annCpgModVec)] <- tmpSid
    
    ensModVec <- append(ensModVec, paste(annMod$Ens, collapse="+"))
    names(ensModVec)[length(ensModVec)] <- tmpSid
    
    hCatVec <-  append(hCatVec, paste(names(annMod$CpG), collapse="+"))
    names(hCatVec)[length(hCatVec)] <- tmpSid
    
    tmpExtP <- inpExtPhenoLst[[tmpSid]]
    # Make sure all samples are present in both pheno and beta
    tmpExtB  <- tmpExtB[,which(colnames(tmpExtB) %in% tmpExtP$barcode)]
    tmpExtP <- tmpExtP[which(tmpExtP$barcode %in% colnames(tmpExtB)),]
    
    # Create dataframe for storing data
    sumCvDfMult <- data.frame(matrix(nrow=0, ncol=7+(4*(length(table(refPheno$Histotype))))))
    colnames(sumCvDfMult) <- c(paste(names(table(refPheno$Histotype)), "Sensitivity",sep="_"), 
                               paste(names(table(refPheno$Histotype)), "Specificity",sep="_"),
                               paste(names(table(refPheno$Histotype)), "Accuracy",sep="_"),
                               paste(names(table(refPheno$Histotype)), "F1Score",sep="_"),
                               "Acc_Overall","Acc_Lower", "Acc_Upper",
                               "AUC", "AUPRC", "Epoch", "Fold")
    
    # Should we use internal or external data for train (i.e. FULL holdout)?
    if(!is.null(extBool)){
      startData <- refBeta[cpgModel,]
      startData <- t(trainData)
      trainP <- refPheno
      testP <- tmpExtP
    }else{
      startData <- t(tmpExtB)
      trainP <- tmpExtP
      testP <- tmpExtP
    }
    
    # If not null remove histotype and associated samples from model 
    if(!is.null(remHistBool)){
      startData <- startData[!rownames(startData) %in% trainP$barcode[which(trainP$Histotype %in% remHistBool)], ]
      trainP <- trainP[!trainP$Histotype %in% remHistBool,]
    }
    
    # Should we optimize the parameters of the data when performing multiclass prediction
    # Increases runtime substantially
    if(!is.null(paramBool)){
      paramDf <- startData
      paramDf <- data.frame(paramDf)
      paramDf$Histotype <- trainP$Histotype[match(rownames(paramDf), trainP$barcode)]
      paramDf$Histotype <- factor(paramDf$Histotype, levels = names(table(trainP$Histotype)))
      best_params <- makeXGBoostCaretGridSMulti(gridDat = paramDf)
    }
    
    # Prepare data for training the model
    # Choose between training the model on our cohort and applying it directly on ext,
    for(t in 1:noEpochs){
      folds <- caret::createFolds(factor(trainP$Histotype), k = nFold, list = FALSE)
      for(k in 1:nFold){
        # Segment train-data by folds, with respect to the fold index no. 
        testIndexes <- which(folds==k, arr.ind=TRUE)
        trainData <- data.frame(startData[-testIndexes, ])
        if(!is.null(extBool)){
          testData <- data.frame(tmpExtB)
        }else{
          testData <- data.frame(startData[testIndexes, ])
        }
        # Add Histotype label to DF
        trainData$Histotype <- trainP$Histotype[match(rownames(trainData), trainP$barcode)]
        testData$Histotype <- testP$Histotype[match(rownames(testData), testP$barcode)]
        trainVals <- trainData[,-which(colnames(trainData) %in% "Histotype")]
        testVals <- testData[,-which(colnames(testData) %in% "Histotype")]
        # Create phenotypic annotation
        trainLabel <- trainData$Histotype
        testLabel <- testData$Histotype
        repVec <- c()
        for(v in 1:length(table(trainLabel))){
          repVec <- append(repVec, v-1)
          names(repVec)[length(repVec)] <- names(table(trainLabel))[v]
          trainLabel <- replace(trainLabel, trainLabel==names(table(trainLabel))[[v]], v-1)
        }
        # Replace values in testlabel using repvec
        for(w in 1:length(repVec)){
          testLabel <- ifelse(testLabel %in% names(repVec)[w], repVec[w], testLabel)
        }
        # Ratio of classes in multiclass
        class_imba <- table(trainLabel)/length(trainLabel)
        cW <- min(class_imba)/class_imba
        cWArr <- rep(NA, length(trainLabel))
        for(x in 1:length(trainLabel)){
          cWArr[x] <- cW[which(names(cW) %in% trainLabel[x])]
        }
        #inpBoost <- as.matrix(trainVals)
        #inpBoost <- sapply(inpBoost, as.numeric)
        if(is.null(paramBool)){
          tmpBoost <- xgboost(data = as.matrix(trainVals),
                              label = trainLabel,
                              booster = "gbtree",
                              objective = "multi:softprob",
                              # weight = cWArr,
                              verbose = 0,
                              nthread = 10,
                              nrounds = 100,
                              num_class = length(table(trainLabel)))
        }else{
          tmpBoost <- xgboost(data = as.matrix(trainVals),
                              label = trainLabel,
                              booster = "gbtree",
                              objective = "multi:softprob",
                              eta =  best_params["eta"] ,
                              max_depth =  best_params["max_depth"],
                              weight = cWArr,
                              min_child_weight = best_params["min_child_weight"],
                              subsample = best_params["subsample"],
                              colsample_bytree = best_params["colsample_bytree"],
                              nrounds = best_params["n_estimators"],
                              verbose = 0,
                              nthread = 10,
                              num_class = length(table(trainLabel)))
        }
        # Apply RF model on test-data (external dataset)
        tmpExtPred <- parsnip::xgb_predict(object=tmpBoost, 
                                           new_data = as.matrix(testVals))
        colnames(tmpExtPred) <- names(table(trainLabel))
        # May not be needed
        predicted_labels <- factor(colnames(tmpExtPred)[max.col(tmpExtPred)],levels=repVec)
        testLabelFac <- factor(testLabel, levels=repVec)
        conMat <- caret::confusionMatrix(data = predicted_labels, 
                                         reference = testLabelFac)
        if(length(levels(testLabelFac))<=2){
          sens <- conMat$byClass[1]
          sens <- append(sens, sens)
          names(sens) <- names(table(testLabel))
          spec <- conMat$byClass[2]
          spec <- append(spec, spec)
          names(spec) <- names(table(testLabel))
          f1Score <- conMat$byClass[7]
          f1Score  <- append(f1Score, f1Score)
          names(f1Score) <- names(table(testLabel))
          acc <- conMat$byClass[11]
          acc  <- append(acc, acc)
          names(acc) <- names(table(testLabel))
        }else{
          sens <- conMat$byClass[,1]
          spec <- conMat$byClass[,2]
          f1Score <- conMat$byClass[,7]
          acc <- conMat$byClass[,11]
          
        }
        # Get overall statistics
        accAll <- conMat$overall[[1]]
        # Get confidence interval of accuracy (lower, upper)
        accLower <- conMat$overall[[3]]
        accUpper <- conMat$overall[[4]]
        
        # Get AUC for model
        rocPred <- tmpExtPred
        colnames(rocPred) <- names(table(testLabel))
        # Get AUC for model
        roc <- pROC::multiclass.roc(predictor = rocPred,
                                    response= testLabel, smooth=FALSE,  quiet = TRUE)
        # Get PRAUC for model
        prAucDf <- as.data.frame(tmpExtPred)
        prAucDf$ref <- testLabelFac
        prauc <- prAucDf %>% yardstick::pr_auc(estimator = "macro",
                                               ref,
                                               min(testLabel):max(testLabel))
        for(u in 1:length(names(table(refPheno$Histotype)))){
          tmpH <- names(table(refPheno$Histotype))[u]
          if(!tmpH %in% names(repVec)){
            sens <- append(sens, 0, after=u-1)
            names(sens)[u] <- tmpH
            spec <- append(spec, 0, after=u-1)
            names(spec)[u] <- tmpH
            f1Score <- append(f1Score, 0, after=u-1)
            names(f1Score)[u] <- tmpH
            acc <- append(acc, 0, after=u-1)
            names(acc)[u] <- tmpH
          }else{
            names(sens)[u] <- tmpH
            names(spec)[u] <- tmpH
            names(acc)[u] <- tmpH
            names(f1Score)[u] <- tmpH
          }
        }
        sumCvDfMult[nrow(sumCvDfMult)+1,] <- c(sens, spec, acc, f1Score,
                                               accAll, accLower, accUpper,
                                               roc$auc[1], prauc$.estimate[1],
                                               t,k)
        
      }
    }
    sumCvDfMult[sapply(sumCvDfMult, is.nan)] <- 0
    sumCvDfMult[sapply(sumCvDfMult, is.na)] <- 0
    meanItDf <- sumCvDfMult %>%
      group_by(Epoch) %>%
      summarise_all("mean")
    
    sumCvDfMult$SID <- tmpSid
    cvMultLst[[s]] <- sumCvDfMult
    names(cvMultLst)[s] <- tmpSid
    
    meanItDf <- data.frame(meanItDf)
    meanItDf[nrow(meanItDf)+1, ] <- colMeans(meanItDf)
    rownames(meanItDf)[nrow(meanItDf)] <- "MeanVal"
    meanItDf$SID <- tmpSid
    cvMultMeanLst[[s]] <- meanItDf
    names(cvMultMeanLst)[s] <- tmpSid
  }
  
  sumDf <- data.frame(matrix(nrow=0, ncol=ncol(cvMultMeanLst[[1]])))
  colnames(sumDf) <- colnames(cvMultMeanLst[[1]])
  for(v in 1:length(cvMultMeanLst)){
    tmpDf <- cvMultMeanLst[[v]]
    meanRow <- tmpDf["MeanVal",]
    sumDf[nrow(sumDf)+1,] <- unlist(meanRow)
  }
  sumDf <- sumDf[,!colnames(sumDf) %in% c("Epoch", "Fold")]
  # Add cpg-model to dataframe
  sumDf$EModel <- ensModVec
  sumDf$CpGModel <- annCpgModVec
  sumDf$ModelCat <- hCatVec
  return(list("Merge"=sumDf, 
              "MeanCvLst"= cvMultMeanLst, 
              "cvAllLst" = cvMultLst, 
              "CpGModel" = annCpgModVec, 
              "ensModel" = ensModVec,
              "modelCat" = hCatVec))
}

makeExtPredMultStepwise <- function(inpCpgDfLst, inpExtBLst, refBeta, inpExtPLst, refPheno, inpCpgLocs=NULL, noEpochs = NULL, nFold = NULL, paramBool = NULL, rmHist = NULL, minCpg = NULL, itMax = NULL, manModLst=NULL, fullModLst=NULL){
  if(is.null(nFold)){
    nFold <- 5
  }
  if(is.null(noEpochs)){
    noEpochs <- 1
  }
  if(is.null(minCpg)){
    minCpg <- 2
  }
  
  if(is.null(itMax)){
    itMax <- 5
  }
  nFold <- nFold
  refPheno <- refPheno[which(refPheno$barcode %in% colnames(refBeta)),]
  refBeta <- refBeta[,which(colnames(refBeta) %in% refPheno$barcode)]
  
  if(!is.null(rmHist)){
    refPheno <- refPheno[!refPheno$Histotype %in% rmHist, ]
    refBeta <- refBeta[,refPheno$barcode]
    for(a in 1:length(inpExtPLst)){
      rmExtP <- inpExtPLst[[a]]
      rmSid <- names(inpExtPLst)[a]
      rmExtB <- inpExtBLst[[rmSid]]
      rmExtP <- rmExtP[!rmExtP$Histotype %in% rmHist, ]
      rmExtB <- rmExtB[,rmExtP$barcode]
      inpExtPLst[[a]] <- rmExtP
      inpExtBLst[[rmSid]] <- rmExtB
    }
  }
  
  if(!is.null(fullModLst)){
    fullModel <- fullModLst[[1]]
    geneModel <- fullModLst[[2]]
  }else if(!is.null(manModLst)){
    fullModel <- c()
    geneModel <- c()
    for(b in 1:length(inpCpgLocs)){
      tmpDf <- inpCpgLocs[[b]]
      if(nrow(tmpDf) == 0){
        next()
      }
      tmpH <- names(inpCpgLocs)[b]
      ensMatches <- tmpDf[which(manModLst %in% tmpDf$Gene), ]
      if(nrow(ensMatches) == 0){
        next()
      }
      ensCpgs <- ensMatches[, grep("cpg", colnames(ensMatches))]
      for(c in 1:nrow(ensCpgs)){
        tmpC <- ensCpgs[c,]
        tmpC <- tmpC[!is.na(tmpC)]
        tmpC <- tmpC[tmpC %in% rownames(refBeta)]
        hDf <- refBeta[tmpC, refPheno$barcode[which(refPheno$Histotype %in% tmpH)]]
        colnames(hDf) <- refPheno$Sample_ID[match(colnames(hDf), refPheno$barcode)]
        nhDf <- refBeta[tmpC, refPheno$barcode[which(!refPheno$Histotype %in% tmpH)]]
        colnames(nhDf) <- refPheno$Sample_ID[match(colnames(nhDf), refPheno$barcode)]
        tmpHistotypes <- names(table(refPheno$Histotype))[!names(table(refPheno$Histotype)) %in% tmpH]
        betaDists <- makeBetaGrpDists(hDf, nhDf, refPheno,  tmpHistotypes, tmpH, histBool=TRUE, distType="TRIMEAN")
        tmpC <- rownames(betaDists)[which(rowSums(betaDists > 0.2) >= 3)]
        names(tmpC) <- paste(tmpH, 1:length(tmpC), sep="")
        fullModel <- append(fullModel, tmpC)
        geneModel <- append(geneModel, ensMatches$Gene[c])
        names(geneModel)[length(geneModel)] <- paste(tmpH, c, sep="")
      }
    }
  }else{
    # Perform prediction for individual genes
    genePredLst <- list()
    for(j in 1:length(inpCpgDfLst)){
      # Second version, perform pred-class on individiual genes
      genePredDf <- data.frame(matrix(nrow=0, ncol=5))
      colnames(genePredDf) <- c("ensembl_gene_id", "Sensitivity", "Specificity", "AUC", "preds")
      inpCpgDf <- inpCpgDfLst[[j]]
      inpHist <- names(inpCpgDfLst)[j]
      if(nrow(inpCpgDf) == 0){
        next()
      }
      for(k in 1:nrow(inpCpgDf)){
        if(abs(k)%%round(nrow(inpCpgDf)/10) == 0){
          message(paste("Processing row: [", k, "/", nrow(inpCpgDf),"]"))
        }
        tmpCpgDf <- data.frame(inpCpgDf[k,])
        tmpEns <- tmpCpgDf$Gene
        # Prepare reference dataframe
        promoDf <- makeCpgDf(tmpCpgDf, refBeta, refPheno, inpHist)
        if(nrow(promoDf) == 0){
          next()
        }
        promoDf <- promoDf[rownames(promoDf) %in% rownames(refBeta)]
        if(nrow(promoDf) <= minCpg){
          message(paste("Gene: ", tmpEns, " Has less then or equal to: ",  minCpg, " total CpG sites and will be skipped", sep=""))
          next()
        }else{
          tmpHistotypes <- names(table(refPheno$Histotype))[!names(table(refPheno$Histotype)) %in% inpHist]
          promoPredDf <- as.data.frame(t(promoDf))
          promoPredDf[is.na(promoPredDf)] <- 0
          # refMat <- round(refMat, 5)
          promoPredDf$Histotype <- refPheno$Histotype[match(rownames(promoPredDf), refPheno$barcode)]
          promoPredDf$Histotype <- ifelse(promoPredDf$Histotype %in% inpHist, inpHist, "Other")
          # if(refPreds$Histotype)
          promoPhenoBetas <- promoPredDf[which(promoPredDf$Histotype %in% inpHist),]
          promoPhenoBetas <- promoPhenoBetas[,-which(colnames(promoPhenoBetas) %in% "Histotype")]
          promoPhenoBetas <- t(promoPhenoBetas)
          colnames(promoPhenoBetas) <- refPheno$Sample_ID[match(colnames(promoPhenoBetas), refPheno$barcode)]
          promoNPhenoBetas <- promoPredDf[-which(promoPredDf$Histotype %in% inpHist),]
          promoNPhenoBetas <- promoNPhenoBetas[,-which(colnames(promoNPhenoBetas) %in% "Histotype")]
          promoNPhenoBetas <- t(promoNPhenoBetas)
          colnames(promoNPhenoBetas) <- refPheno$Sample_ID[match(colnames(promoNPhenoBetas), refPheno$barcode)]
          # Calculate the mean* distance for predictors 
          # * TriMean or Geometric mean
          tmpHistotypes <- names(table(refPheno$Histotype))[!names(table(refPheno$Histotype)) %in% inpHist]
          betaDists <- makeBetaGrpDists(promoPhenoBetas, promoNPhenoBetas, refPheno,  tmpHistotypes, inpHist, histBool=TRUE, distType="TRIMEAN")
          # Filter CpG sites passing delta beta >= 0.2 for 2,3 histotypes for predictors
          minDists <- betaDists[which(rowSums(betaDists >= 0.2) >= 3),]
          if(nrow(minDists) < 2){
            next()
          }else{
            # Prepare training-data  based on the less stringent filtering
            refPredDf <- refBeta[rownames(minDists), ]
            refPredDf <- as.data.frame(t(refPredDf))
            refPredDf  <- refPredDf[rowSums(is.na(refPredDf)) == 0,]
            refPredDf$Histotype <- refPheno$Histotype[match(rownames(refPredDf), refPheno$barcode)]
            refPredDf$Histotype <- ifelse(refPredDf$Histotype %in% inpHist, inpHist, "Other")
            refPredDf$Histotype <- factor(refPredDf$Histotype, levels = c("Other", inpHist))
            # Create models to be used with data
            minModel <- paste(rownames(minDists), collapse = " + ")
            minModel <- reformulate(minModel, "Histotype")
            # Shuffle rows to remove any bias when sampling
            refPredDf <- refPredDf[sample(nrow(refPredDf)),]
            kfDfLst <- list()
            if(!is.null(paramBool)){
              # Perform grid-search to identify best parameters for model
              best_params <- makeXGBoostCaretGridS(refPredDf)
            }
            # Perform xgboost on training cohort
            sumCvDf <- data.frame(matrix(NA, nrow=noEpochs, ncol=3))
            colnames(sumCvDf) <- c("Sensitivity", "Specificity", "AUC")
            inpRef <- refPredDf
            kfDfLst <- list()
            for(l in 1:noEpochs){
              folds <- caret::createFolds(factor(inpRef$Histotype), k = nFold, list = FALSE) 
              # predStatLst <- list()
              kfDf <- data.frame(matrix(nrow=0, ncol=3))
              colnames(kfDf) <- c("Sensitivity", "Specificity", "AUC")
              for(m in 1:nFold){
                # Segment train-data by folds, with respect to the fold index no. 
                testIndexes <- which(folds==m, arr.ind=TRUE)
                testData <- inpRef[testIndexes, ]
                trainData <- inpRef[-testIndexes, ]
                # best_params <- makeXGBoostCaretGridS(trainData)
                # Create weights for random forest classification (balance)
                trainTrain <- trainData[,-which(grepl("Histotype",colnames(trainData)))]
                trainLabel <- trainData[,which(grepl("Histotype",colnames(trainData)))]
                trainLabel <- ifelse(trainLabel == 'Other', 0, 1)
                # Ratio of negative to positive classes
                # class_imba <- length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
                class_imba <- length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
                # best_params <- makeXGBoostGridS(trainData)
                # Perform grid-search on training data to retrieve optimal model parameters
                # best_params <- makeXGBoostGridS(trainData)
                # Build gradient boosting forest model
                inpBoost <- as.matrix(trainTrain)
                class(inpBoost) <- "numeric"
                
                if(is.null(paramBool)){
                  tmpBoost <- xgboost(data = inpBoost,
                                      label = trainLabel,
                                      booster = "gbtree",
                                      objective = "binary:logistic",
                                      scale_pos_weight = class_imba,
                                      nrounds = 100,
                                      verbose = 0,
                                      nthread = 10)
                }else{
                  tmpBoost <- xgboost(data = inpBoost,
                                      label = trainLabel,
                                      booster = "gbtree",
                                      objective = "binary:logistic",
                                      eta =  best_params["eta"] ,
                                      max_depth =  best_params["max_depth"],
                                      scale_pos_weight = class_imba,
                                      min_child_weight = best_params["min_child_weight"],
                                      subsample = best_params["subsample"],
                                      colsample_bytree = best_params["colsample_bytree"],
                                      nrounds = best_params["n_estimators"],
                                      verbose = 0,
                                      nthread = 10)
                }
                # Apply RF model on test-data (external dataset)
                tmpExtPred <- predict(object=tmpBoost, 
                                      newdata = as.matrix(testData[,-which(grepl("Histotype",colnames(testData)))]))
                extStats <- testData
                extStats$PredStats <- tmpExtPred 
                # Categorize predictions based on their corresponding classification seen to model
                conMat <- confusionMatrix(data = factor(as.numeric(extStats$PredStats >0.5), levels = c(1,0)), 
                                          reference = factor(as.numeric(ifelse(extStats$Histotype %in% "Other", 0,1)), levels=c(1,0)))
                sens <- conMat[[4]][[2]]
                spec <- conMat[[4]][[1]]
                # Get AUC for model
                roc <- pROC::roc(testData$Histotype, tmpExtPred, smooth=FALSE,  quiet = TRUE)
                kfDf[nrow(kfDf)+1,] <- c(sens, spec, roc$auc[1])
                rownames(kfDf)[nrow(kfDf)] <- paste("k=",m, sep="")
              }
              kfDf[sapply(kfDf, is.nan)] <- 0
              kfDf[nrow(kfDf)+1, ] <- colMeans(kfDf)
              sumCvDf[l,] <- colMeans(kfDf)
            }
            meanPred <- colMeans(sumCvDf)
            genePredDf[nrow(genePredDf)+1,] <- c(tmpEns, meanPred[1], meanPred[2], meanPred[3], paste(rownames(minDists), collapse = ","))
          }
        }
      }
      # Select only models with a good combined AUC
      genePredDf <- genePredDf[which(genePredDf$AUC > 0.8), ]
      genePredLst[[length(genePredLst)+1]] <- genePredDf
      names(genePredLst)[length(genePredLst)] <- inpHist
    }
    message("Individual gene-pred completed")
    
    # Remove empty entries
    genePredLst <- genePredLst[!names(genePredLst) %in% ""]
    genePredLst <- genePredLst[sapply(genePredLst, nrow) > 0]
    modelLst <- list()
    for(q in 1:length(genePredLst)){
      bestRow <- NULL
      genePredLst2 <- list()
      sumDf <- data.frame(matrix(nrow=0, ncol=5))
      colnames(sumDf) <- c("Model", "Sensitivity", "Specificity", "AUC", "CpG")
      tmpH <- names(genePredLst)[q]
      gDf <- genePredLst[[tmpH]]
      nSamp <- nrow(gDf)
      if(nSamp == 1){
        nIt <- 1 
      }else if(nSamp < itMax){
        nIt <- nSamp-1
      }else{
        nIt <- itMax-1
      }
      for(m in 1:nIt){
        genePredDf2 <- data.frame(matrix(nrow=0, ncol=5))
        colnames(genePredDf2) <- c("Model", "Sensitivity", "Specificity", "AUC", "CpG")
        if(is.null(bestRow)){
          bestRow <- gDf[which.max(gDf$AUC),]
          ens1 <- bestRow$ensembl_gene_id
          if(nIt == 1){
            tmpDf <- gDf
          }else{
            tmpDf <- gDf[-which.max(gDf$AUC), ]
          }
          # Create new model
          pred1 <- bestRow$preds
          pred1 <- strsplit(pred1, ",")[[1]]
          pred1str <- paste(pred1, collapse = " + ")
        }else{
          bestRow <- sumDf[nrow(sumDf),]
          ens1 <- bestRow$Model
          ens1 <- strsplit(ens1, "\\+")[[1]]
          tmpDf <- gDf[-which(gDf$ensembl_gene_id %in% ens1), ]
          # Create new model
          bPreds <- gDf$preds[which(gDf$ensembl_gene_id %in% ens1)]
          bPreds <- paste(bPreds, collapse = ",")
          pred1 <- strsplit(bPreds, ",")[[1]]
          pred1str <- paste(pred1, collapse = " + ")
        }
        for(n in 1:nrow(tmpDf)){
          nRow <- tmpDf[n,]
          ens2 <- nRow$ensembl_gene_id
          pred2 <- nRow$preds
          pred2 <- strsplit(pred2, ",")[[1]]
          pred2str <- paste(pred2, collapse = " + ")
          nMod <- paste(pred1str, pred2str, sep = " + ")
          nMod <- reformulate(nMod, "Histotype")  
          modPreds <- c(pred1, pred2)
          # Try new model
          modDf <- refBeta[modPreds, ]
          modDf <- as.data.frame(t(modDf))
          modDf  <- modDf[rowSums(is.na(modDf)) == 0,]
          modDf$Histotype <- refPheno$Histotype[match(rownames(modDf), refPheno$barcode)]
          modDf$Histotype <- ifelse(modDf$Histotype %in% tmpH, tmpH, "Other")
          modDf$Histotype <- factor(modDf$Histotype, levels = c("Other", tmpH))
          # Shuffle rows to remove any bias when sampling
          modDf <- modDf[sample(nrow(modDf)),]
          # Perform xgboost on training cohort
          sumCvDf2 <- data.frame(matrix(NA, nrow=0, ncol=3))
          colnames(sumCvDf2) <- c("Sensitivity", "Specificity", "AUC")
          folds <- caret::createFolds(factor(modDf$Histotype), k = nFold, list = FALSE) 
          # predStatLst <- list()
          kfDf2 <- data.frame(matrix(nrow=0, ncol=3))
          colnames(kfDf2) <- c("Sensitivity", "Specificity", "AUC")
          for(p in 1:nFold){
            # Segment train-data by folds, with respect to the fold index no. 
            testIndexes <- which(folds==p, arr.ind=TRUE)
            testData <- modDf[testIndexes, ]
            trainData <- modDf
            trainData <- modDf[-testIndexes, ]
            # best_params <- makeXGBoostCaretGridS(trainData)
            # Create weights for random forest classification (balance)
            trainTrain <- trainData[,-which(grepl("Histotype",colnames(trainData)))]
            trainLabel <- trainData[,which(grepl("Histotype",colnames(trainData)))]
            trainLabel <- ifelse(trainLabel == 'Other', 0, 1)
            # Ratio of negative to positive classes
            class_imba <- length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
            # length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
            # best_params <- makeXGBoostGridS(trainData)
            # Perform grid-search on training data to retrieve optimal model parameters
            # best_params <- makeXGBoostGridS(trainData)
            # Build gradient boosting forest model
            inpBoost <- as.matrix(trainTrain)
            class(inpBoost) <- "numeric"
            if(is.null(paramBool)){
              tmpBoost <- xgboost(data = inpBoost,
                                  label = trainLabel,
                                  booster = "gbtree",
                                  objective = "binary:logistic",
                                  scale_pos_weight = class_imba,
                                  nrounds = 100,
                                  verbose = 0,
                                  nthread = 10)
            }else{
              tmpBoost <- xgboost(data = inpBoost,
                                  label = trainLabel,
                                  booster = "gbtree",
                                  objective = "binary:logistic",
                                  eta =  best_params["eta"] ,
                                  max_depth =  best_params["max_depth"],
                                  scale_pos_weight = class_imba,
                                  min_child_weight = best_params["min_child_weight"],
                                  subsample = best_params["subsample"],
                                  colsample_bytree = best_params["colsample_bytree"],
                                  nrounds = best_params["n_estimators"],
                                  verbose = 0,
                                  nthread = 10)
            }
            # Apply RF model on test-data (external dataset)
            tmpExtPred <- predict(object=tmpBoost, 
                                  newdata = as.matrix(testData[,-which(grepl("Histotype",colnames(testData)))]))
            extStats <- testData
            extStats$PredStats <- tmpExtPred 
            conMat <- confusionMatrix(data = factor(as.numeric(extStats$PredStats >0.5), levels = c(1,0)), 
                                      reference = factor(as.numeric(ifelse(extStats$Histotype %in% "Other", 0,1)), levels=c(1,0)))
            sens <- conMat[[4]][[2]]
            spec <- conMat[[4]][[1]]
            # Get AUC for model
            roc <- pROC::roc(testData$Histotype, tmpExtPred, smooth=FALSE,  quiet = TRUE)
            kfDf2[nrow(kfDf2)+1,] <- c(sens, spec, roc$auc[1])
            kfDf2[sapply(kfDf2, is.nan)] <- 0
            #kfDf2[nrow(kfDf2)+1, ] <- colMeans(kfDf2)
            #rownames(kfDf2)[nrow(kfDf2)] <- "MeanVal"
            sumCvDf2[nrow(sumCvDf2)+1,] <- colMeans(kfDf2)
            #rownames(kfDf2)[nrow(kfDf2)] <- paste("k=",p, sep="")
          }
          meanPred2 <- colMeans(sumCvDf2)
          genePredDf2[nrow(genePredDf2)+1,] <- c(paste(paste(ens1, collapse="+"), ens2, sep="+"), meanPred2[1], meanPred2[2], meanPred2[3], paste(modPreds, collapse=","))
        }
        genePredDfSel <- genePredDf2
        genePredDf2$SensSpec <- as.numeric(genePredDf2$Sensitivity) + as.numeric(genePredDf2$Specificity)
        genePredDf2 <- genePredDf2[order(genePredDf2$SensSpec, decreasing=TRUE), ]
        bestModel <-  genePredDf2[1,]
        bestModel <- bestModel[!names(bestModel) %in% "SensSpec"]
        # nBestModel <- genePredDf2[2,]
        sumDf[nrow(sumDf)+1,] <- bestModel 
        bestRow <- sumDf[nrow(sumDf),]
      }
      modelLst[[q]] <- sumDf
      names(modelLst)[q] <- tmpH
    }
    message("Model combination analysis completed")
    
    # Combine all different CpG-sets, 
    # Use the most beneficial combination of promoters as a gene-panel for histotype stratification
    fullModel <- c()
    geneModel <- c()
    for(r in 1:length(modelLst)){
      tmpModels <- modelLst[[r]]
      tmpH <- names(modelLst)[r]
      tmpModels$SensSpec <- as.numeric(tmpModels$Sensitivity) + as.numeric(tmpModels$Specificity)
      bModel <- tmpModels[which.max(tmpModels$AUC), ]
      tmpGM <- bModel$Model
      tmpGM <- strsplit(tmpGM , "\\+")[[1]]
      names(tmpGM) <- paste(tmpH, 1:length(tmpGM), sep="")
      tmpFM <- bModel$CpG
      tmpFM <- strsplit(tmpFM , ",")[[1]]
      names(tmpFM) <-  paste(tmpH, 1:length(tmpFM), sep="")
      fullModel <- append(fullModel, tmpFM)
      geneModel <- append(geneModel, tmpGM)
    }
  }
  message("Combined best model created")
  
  # Run multiclass prediction on test-datasets
  extPredLst <- list()
  cvMultLst <- list()
  for(s in 1:length(inpExtBLst)){
    tmpExtB <- inpExtBLst[[s]]
    tmpSid <- names(inpExtBLst)[s]
    tmpExtP <- inpExtPLst[[tmpSid]]
    tmpExtB  <- tmpExtB[,which(colnames(tmpExtB) %in% tmpExtP$barcode)]
    tmpExtP <- tmpExtP[which(tmpExtP$barcode %in% colnames(tmpExtB)),]
    # Get CpG sites which are shared between the two input arrays
    sharedCpgs <- intersect(rownames(tmpExtB), rownames(refBeta))
    # Make sure to only include CpG sites present in both dataframes (i.e. drop EPIC only probes to make dataframes compatible)
    tmpExtB <- tmpExtB[sharedCpgs,]
    tmpEFM <- fullModel[which(fullModel %in% rownames(tmpExtB))]
    # Repeat for external dataframe
    tmpExtB <- data.frame(t(tmpExtB[tmpEFM,]))
    tmpExtB$Histotype <- tmpExtP$Histotype[match(rownames(tmpExtB), tmpExtP$barcode)]
    tmpExtB$Histotype <- factor(tmpExtB$Histotype, levels = names(table(tmpExtP$Histotype))) 
    
    if(!is.null(paramBool)){
      # Perform grid-search to identify best parameters for model
      best_params2 <- makeXGBoostCaretGridS(tmpExtB)
    }
    
    sumCvDfMult <- data.frame(matrix(nrow=0, ncol=3+(2*(length(table(trainPheno$Histotype))))))
    colnames(sumCvDfMult) <- c(paste(names(table(trainPheno$Histotype)), "Sensitivity",sep="_"), 
                               paste(names(table(trainPheno$Histotype)), "Specificity",sep="_"),
                               "Overall_Accuracy","Accuracy_p", "Overall_AUC")
    
    
    for(t in 1:100){
      folds <- caret::createFolds(factor(tmpExtB$Histotype), k = nFold, list = FALSE) 
      # predStatLst <- list()
      kfDf3 <- data.frame(matrix(nrow=0, ncol=3+(2*(length(table(refPheno$Histotype))))))
      colnames(kfDf3) <- c(paste(names(table(refPheno$Histotype)), "Sensitivity",sep="_"), 
                           paste(names(table(refPheno$Histotype)), "Specificity",sep="_"),
                           "Overall_Accuracy","Accuracy_p", "Overall_AUC")
      for(u in 1:nFold){
        # Segment train-data by folds, with respect to the fold index no. 
        testIndexes <- which(folds==u,arr.ind=TRUE)
        tmpExtBTest <- tmpExtB[testIndexes, ]
        tmpExtBTrain <- tmpExtB[-testIndexes, ]
        # Create weights for random forest classification (balance)
        trainTrain <- tmpExtBTrain[,-which(grepl("Histotype",colnames(tmpExtBTrain)))]
        trainLabel <- tmpExtBTrain[,which(grepl("Histotype",colnames(tmpExtBTrain)))]
        trainLabel <- as.character(trainLabel)
        for(v in 1:length(table(trainLabel))){
          trainLabel <- replace(trainLabel, trainLabel==names(table(trainLabel))[[v]], v-1)
        }
        trainLabel <- as.numeric(trainLabel)
        
        # Ratio of negative to positive classes
        class_imba <- table(trainLabel)/length(trainLabel)
        cW <- min(class_imba)/class_imba
        cWArr <- rep(NA, length(trainLabel))
        for(x in 1:length(trainLabel)){
          cWArr[x] <- cW[which(names(cW) %in% trainLabel[x])]
        }
        
        inpBoost <- as.matrix(trainTrain)
        class(inpBoost) <- "numeric"
        if(is.null(paramBool)){
          tmpBoost <- xgboost(data = inpBoost,
                              label = trainLabel,
                              booster = "gbtree",
                              objective = "multi:softprob",
                              weight = cWArr,
                              verbose = 0,
                              nthread = 10,
                              nrounds = 100,
                              num_class = length(table(tmpExtB$Histotype)))
        }else{
          tmpBoost <- xgboost(data = inpBoost,
                              label = trainLabel,
                              booster = "gbtree",
                              objective = "multi:softprob",
                              eta =  best_params2["eta"] ,
                              max_depth =  best_params2["max_depth"],
                              weight = cWArr,
                              min_child_weight = best_params2["min_child_weight"],
                              subsample = best_params2["subsample"],
                              colsample_bytree = best_params2["colsample_bytree"],
                              nrounds = best_params2["n_estimators"],
                              verbose = 0,
                              nthread = 10,
                              num_class = length(table(tmpExtB$Histotype)))
        }
        
        # Apply RF model on test-data (external dataset)
        tmpExtPred <- parsnip::xgb_predict(object=tmpBoost, 
                                           new_data = as.matrix(tmpExtBTest[,-which(grepl("Histotype",colnames(tmpExtBTest)))]))
        testLabel <- as.character(tmpExtBTest$Histotype)
        for(z in 1:length(table(tmpExtB$Histotype))){
          testLabel <- replace(testLabel, testLabel==names(table(tmpExtB$Histotype))[[z]], z-1)
        }
        testLabel <-as.numeric(testLabel)
        rocPred <- tmpExtPred
        colnames(rocPred) <- names(table(trainLabel))
        # Get AUC for model
        roc <- pROC::multiclass.roc(response=testLabel, predictor = rocPred, smooth=FALSE,  quiet = TRUE)
        colnames(tmpExtPred) <- names(table(trainLabel))
        predicted_labels=factor(colnames(tmpExtPred)[max.col(tmpExtPred)],levels=0:(length(table(trainLabel))-1))
        testLabel <- factor(testLabel, levels = levels(predicted_labels))
        conMat <- caret::confusionMatrix(testLabel,predicted_labels)
        conAll <- conMat$overall
        conClass <- data.frame(conMat$byClass)
        rownames(conClass) <- names(table(tmpExtBTrain$Histotype))
        sensVec <- conClass$Sensitivity
        specVec <- conClass$Specificity
        for(a in 1:length(names(table(trainPheno$Histotype)))){
          tmpH <- names(table(trainPheno$Histotype))[a]
          if(!tmpH %in% rownames(conClass)){
            sensVec <- append(sensVec, 0, after=a-1)
            specVec <- append(specVec, 0, after=a-1)
          }
        }
        kfDf3[nrow(kfDf3)+1,] <- c(sensVec, specVec, conAll[[1]], conAll[[6]], roc$auc[1])
      }
      kfDf3[sapply(kfDf3, is.nan)] <- 0
      kfDf3[sapply(kfDf3, is.na)] <- 0
      kfDf3[nrow(kfDf3)+1, ] <- colMeans(kfDf3)
      rownames(kfDf3)[nrow(kfDf3)] <- "MeanVal"
      sumCvDfMult[t,] <- kfDf3["MeanVal",]
      rownames(sumCvDfMult)[t] <- t
    }
    cvMultLst[[s]] <- sumCvDfMult
    names(cvMultLst)[s] <- tmpSid
    sumCvDfMult[nrow(sumCvDfMult)+1, ] <- colMeans(sumCvDfMult)
    rownames(sumCvDfMult)[nrow(sumCvDfMult)] <- "MeanVal"
    sumCvDfMult$SID <- tmpSid
    extPredLst[[s]] <- sumCvDfMult[nrow(sumCvDfMult),]
    names(extPredLst)[s] <- tmpSid 
  }
  mergeExtDf <- Reduce(function(...) merge(..., all=T), extPredLst)
  mergeExtDf$Model <- paste(fullModel, collapse="+")
  mergeExtDf$GModel <-paste(geneModel, collapse="+")
  mergeExtDf$HistotypeModel <- paste(names(geneModel), collapse="+")
  return(list("Merge"=mergeExtDf, "PredLst"= cvMultLst))
}

# Script for binary classification
makeBinPred <- function(refBeta, refPheno, inpExtBLst, inpExtPLst, fullModel, inpLocs, nFold=NULL, noEpochs=NULL, extBool =NULL){
  # Function for binay classification using a pre-determined model of CpG sites
  if(is.null(nFold)){
    nFold <- 5
  }else{
    nFold <- nFold
  }
  if(is.null(noEpochs)){
    noEpochs <- 100
  }
  
  ensModel <- c()
  # Annotate the cpg-model
  for(a in 1:length(inpLocs)){
    tmpLoc <- inpLocs[[a]]
    tmpH <- names(inpLocs)[a]
    for(b in 1:nrow(tmpLoc)){
      if(length(which(tmpLoc[b,] %in% fullModel)) > 0){
        tmpCpgs <- tmpLoc[b,which(tmpLoc[b,] %in% fullModel)]
        names(fullModel)[which(fullModel %in% tmpCpgs)] <- paste(tmpH, 1:length(tmpCpgs), sep="")
        for(c in 1:length(tmpCpgs)){
          ensModel <- append(ensModel, tmpLoc$Gene[b])
        }
      }
    }
  }
  
  # Run binary prediction on test-datasets
  extPredDf <- data.frame(matrix(nrow=length(inpExtBLst), ncol=3*(length(table(refPheno$Histotype)))))
  colnames(extPredDf) <- c(paste(names(table(refPheno$Histotype)), "Sensitivity",sep="_"), 
                           paste(names(table(refPheno$Histotype)), "Specificity",sep="_"),
                           paste(names(table(refPheno$Histotype)), "AUC",sep="_")) 
  cvMultLst <- list()
  for(i in 1:length(inpExtBLst)){
    sumCvDfBin <- data.frame(matrix(nrow=noEpochs, ncol=3*(length(table(refPheno$Histotype)))))
    colnames(sumCvDfBin) <- c(paste(names(table(refPheno$Histotype)), "Sensitivity",sep="_"), 
                              paste(names(table(refPheno$Histotype)), "Specificity",sep="_"),
                              paste(names(table(refPheno$Histotype)), "AUC",sep="_"))
    
    extB <- inpExtBLst[[i]]
    tmpSid <- names(inpExtBLst)[i]
    tmpExtP <- inpExtPLst[[tmpSid]]
    extB  <- extB[,which(colnames(extB) %in% tmpExtP$barcode)]
    tmpExtP <- tmpExtP[which(tmpExtP$barcode %in% colnames(extB)),]
    # Get CpG sites which are shared between the two input arrays
    sharedCpgs <- intersect(rownames(extB), rownames(refBeta))
    # Make sure to only include CpG sites present in both dataframes (i.e. drop EPIC only probes to make dataframes compatible)
    extB <- extB[sharedCpgs,]
    tmpEFM <- fullModel[which(fullModel %in% rownames(extB))]
    
    extB <- data.frame(t(extB[tmpEFM,]))
    extB$Histotype <- tmpExtP$Histotype[match(rownames(extB), tmpExtP$barcode)]
    extB$Histotype <- factor(extB$Histotype, levels = names(table(tmpExtP$Histotype))) 
    
    # Repeat for external dataframe
    tB <- data.frame(t(refBeta[tmpEFM,]))
    tB$Histotype <- refPheno$Histotype[match(rownames(tB), refPheno$barcode)]
    tB$Histotype <- factor(tB$Histotype, levels = names(table(refPheno$Histotype)))
    tB <- tB[which(tB$Histotype %in% names(table(extB$Histotype))),]
    
    extB$Histotype <- droplevels(extB$Histotype)
    tB$Histotype <- droplevels(tB$Histotype)
    
    for(j in 1:length(names(table(extB$Histotype)))){
      tmpExtB <- extB
      tmpTB <- tB
      tmpH <- names(table(extB$Histotype))[j]
      if(length(grep(tmpH, names(fullModel))) == 0){
        next()
      }else{
        tmpExtB <- tmpExtB[, which(colnames(tmpExtB) %in% c(fullModel[grep(tmpH, names(fullModel))], "Histotype"))]
        tmpTB <- tmpTB[, which(colnames(tmpTB) %in% c(fullModel[grep(tmpH, names(fullModel))], "Histotype"))]
      }
      tmpExtB$Histotype <- ifelse(tmpExtB$Histotype %in% tmpH, tmpH, "Other")
      tmpTB$Histotype <- ifelse(tmpTB$Histotype %in% tmpH, tmpH, "Other")
      for(t in 1:noEpochs){
        folds <- caret::createFolds(factor(tmpExtB$Histotype), k = nFold, list = FALSE) 
        # predStatLst <- list()
        kfDf3 <- data.frame(matrix(nrow=0, ncol=3))
        colnames(kfDf3) <- c("Sensitivity", 
                             "Specificity",
                             "AUC")
        for(u in 1:nFold){
          # Segment train-data by folds, with respect to the fold index no. 
          testIndexes <- which(folds==u, arr.ind=TRUE)
          tmpExtBTest <- tmpExtB[testIndexes, ]
          tmpExtBTrain <- tmpExtB[-testIndexes, ]
          
          # Create weights for random forest classification (balance)
          trainTrain <- tmpExtBTrain[,-which(grepl("Histotype",colnames(tmpExtBTrain)))]
          trainLabel <- tmpExtBTrain[,which(grepl("Histotype",colnames(tmpExtBTrain)))]
          trainLabel <- ifelse(trainLabel %in% tmpH, 1,0)
          trainLabel <- as.numeric(trainLabel)
          
          # Create weights for random forest classification (balance)
          #trainTrain <- tmpTB[,-which(grepl("Histotype",colnames(tmpTB)))]
          #trainLabel <- tmpTB[,which(grepl("Histotype",colnames(tmpTB)))]
          #trainLabel <- ifelse(trainLabel %in% tmpH, 1,0)
          #trainLabel <- as.numeric(trainLabel)
          
          # Ratio of negative to positive classes
          class_imba <- table(trainLabel)/length(trainLabel)
          #cW <- min(class_imba)/class_imba
          #cWArr <- rep(NA, length(trainLabel))
          #for(x in 1:length(trainLabel)){
          #  cWArr[x] <- cW[which(names(cW) %in% trainLabel[x])]
          #}
          tBoost <- as.matrix(trainTrain)
          class(tBoost) <- "numeric"
          tmpBoost <- xgboost(data = tBoost,
                              label = trainLabel,
                              booster = "gbtree",
                              objective = "binary:logistic",
                              eval_metric = "logloss",
                              scale_pos_weight = class_imba, 
                              verbose = 0,
                              nthread = 10,
                              nrounds = 100)
          
          # Apply RF model on test-data (external dataset)
          tmpExtPred <- predict(object=tmpBoost, 
                                newdata = as.matrix(tmpExtBTest[,-which(grepl("Histotype",colnames(tmpExtBTest)))]))
          extStats <- tmpExtBTest
          extStats$PredStats <- tmpExtPred 
          conMat <- confusionMatrix(data = factor(as.numeric(extStats$PredStats >0.5), levels = c(1,0)), 
                                    reference = factor(as.numeric(ifelse(extStats$Histotype %in% tmpH, 1, 0)), levels=c(1,0)))
          conAll <- conMat$overall
          sens <- conMat[[4]][[1]]
          spec <- conMat[[4]][[2]]
          # Get AUC for model
          roc <- pROC::roc(tmpExtBTest$Histotype, tmpExtPred, smooth=FALSE,  quiet = TRUE)
          kfDf3[nrow(kfDf3)+1,] <- c(sens, spec, roc$auc[1])
        }
        kfDf3[sapply(kfDf3, is.nan)] <- 0
        kfDf3[sapply(kfDf3, is.na)] <- 0
        kfDf3[nrow(kfDf3)+1, ] <- colMeans(kfDf3)
        rownames(kfDf3)[nrow(kfDf3)] <- "MeanVal"
        colnames(kfDf3) <- paste(tmpH, colnames(kfDf3), sep="_")
        colInds <- which(colnames(sumCvDfBin) %in% colnames(kfDf3))
        sumCvDfBin[t,colInds] <- kfDf3["MeanVal",]
      }
    }
    sumCvDfBin[is.na(sumCvDfBin)] <- 0
    sumCvDfBin[nrow(sumCvDfBin)+1, ] <- colMeans(sumCvDfBin)
    rownames(sumCvDfBin)[nrow(sumCvDfBin)] <- "MeanVal"
    cvMultLst[[i]] <- sumCvDfBin
    names(cvMultLst)[i] <- tmpSid
    extPredDf[i,] <- sumCvDfBin[nrow(sumCvDfBin), ]
    rownames(extPredDf)[i] <- tmpSid
  }
  extPredDf <- rownames_to_column(extPredDf, "SID")
  return(list("Merge"=extPredDf, "PredLst"= cvMultLst, "CpG_Model" = fullModel, "EnsModel" = ensModel))
}

makeDirExtMultPred <- function(inpExtBetaLst, inpExtPhenoLst, inpCModel, refPheno, sigCpgLocs, refBeta, noEpochs = NULL, nFold = NULL, paramBool = NULL, remHistBool=NULL){
  # Script for applying model directly on external datasets after training on source dataset
  if(is.null(noEpochs)){
    noEpochs <- 100
  }
  #if(is.null(nFold)){
  #  nFold <- 5
  #}
  # cpgModel <-  sapply(inpCModel, function(x) strsplit(x, "\\+")[[1]][1], USE.NAMES=FALSE)
  # ensModel <- sapply(inpEModel, function(x) strsplit(x, "\\+")[[1]][1], USE.NAMES=FALSE)
  # Run multiclass prediction on test-datasets
  extPredLst <- list()
  cvMultLst <- list()
  cvMultMeanLst <- list()
  cpgModVec <- c()
  annCpgModVec <- c()
  ensModVec <- c()
  missPredVec <- c()
  hCatVec <- c()
  
  # If not null remove histotype and associated samples from model 
  if(!is.null(remHistBool)){
    refBeta<- refBeta[, !colnames(refBeta) %in% trainPheno$barcode[which(trainPheno$Histotype %in% remHistBool)]]
    refPheno <- refPheno[!refPheno$Histotype %in% remHistBool,]
    for(x in 1:length(extBetaLst)){
      tmpB <- inpExtBetaLst[[x]]
      tmpP <- inpExtPhenoLst[[x]]
      tmpB <- tmpB[, !colnames(tmpB) %in% tmpP$barcode[which(tmpP$Histotype %in% remHistBool)]]
      tmpP <- tmpP[!tmpP$Histotype %in% remHistBool,]
      inpExtBetaLst[[x]] <- tmpB
      inpExtPhenoLst[[x]] <- tmpP
    }
  }
  
  # Set up starting data for model
  cpgModel <- inpCModel
  startData <- refBeta[cpgModel,]
  startData <- t(startData)
  trainP <- refPheno
  if(is.null(paramBool)){
    # Tune parameters for XGB model before we apply it
    tuneData <- data.frame(startData)
    tuneData$Histotype <- trainP$Histotype[match(rownames(tuneData), trainP$barcode)]
    best_params <- makeXGBoostCaretGridSMulti(gridDat = tuneData)
  }else if(!isFALSE(paramBool)){
    best_params <- paramBool
  }
  
  for(s in 1:length(inpExtBetaLst)){
    cpgModel <- inpCModel
    tmpSid <- names(inpExtBetaLst)[s]
    message(paste("Performing multiclass classification on cohort: ", tmpSid, sep=""))
    tmpExtB <- inpExtBetaLst[[tmpSid]]
    # Control what CpGs are in the dataframe compared to the model
    if(length(which(!cpgModel %in% rownames(tmpExtB))) > 0 | length(which(!cpgModel %in% rownames(refBeta))) > 0){
      # Get missing CpG's
      refMiss <- paste(cpgModel[which(!cpgModel %in% rownames(refBeta))])
      extMiss <- paste(cpgModel[which(!cpgModel %in% rownames(tmpExtB))])
      allMiss <- c(refMiss, extMiss)
      allMiss <- allMiss[!allMiss %in% ""]
      message(paste("Warning, CpGs: ", paste(allMiss,collapse=", "), 
                    " not found in external beta-matrix, and will be excluded from analysis", sep=""))
      
      missPredVec <- append(missPredVec, paste(allMiss, collapse="+"))
      cpgModel <- cpgModel[!cpgModel %in% allMiss]
    }else{
      missPredVec <- append(missPredVec, 0)
    }
    names(missPredVec)[length(missPredVec)] <- tmpSid
    # Annotate model with ensembl-id's
    annMod <- getAnnotatedCpgModel(inpCpg=cpgModel, 
                                   inpPromoLocs = sigCpgLocs)
    
    # Save model for later use
    cpgModVec <- append(cpgModVec, paste(annMod$CpG, collapse="+"))
    names(cpgModVec)[length(cpgModVec)] <- tmpSid
    
    # Save annotated model
    annCpgModVec <- append(annCpgModVec, paste(annMod$CpG, collapse="+"))
    names(annCpgModVec)[length(annCpgModVec)] <- tmpSid
    
    ensModVec <- append(ensModVec, paste(annMod$Ens, collapse="+"))
    names(ensModVec)[length(ensModVec)] <- tmpSid
    
    hCatVec <-  append(hCatVec, paste(names(annMod$CpG), collapse="+"))
    names(hCatVec)[length(hCatVec)] <- tmpSid
    
    tmpExtP <- inpExtPhenoLst[[tmpSid]]
    # Make sure all samples are present in both pheno and beta
    tmpExtB  <- tmpExtB[,which(colnames(tmpExtB) %in% tmpExtP$barcode)]
    tmpExtP <- tmpExtP[which(tmpExtP$barcode %in% colnames(tmpExtB)),]
    # Create dataframe for storing data
    sumCvDfMult <- data.frame(matrix(nrow=0, ncol=7+(4*(length(table(trainP$Histotype))))))
    colnames(sumCvDfMult) <- c(paste(names(table(trainP$Histotype)), "Sensitivity",sep="_"), 
                               paste(names(table(trainP$Histotype)), "Specificity",sep="_"),
                               paste(names(table(trainP$Histotype)), "Accuracy",sep="_"),
                               paste(names(table(trainP$Histotype)), "F1Score",sep="_"),
                               "Acc_Overall","Acc_Lower", "Acc_Upper",
                               "AUC", "AUPRC", "Epoch", "Fold")
    # Prepare data for training the model
    # Choose between training the model on our cohort and applying it directly on ext,
    tBData <- tmpExtB[which(rownames(tmpExtB) %in% annMod$CpG),] 
    tBData <- t(tBData)
    startData <- startData[,which(colnames(startData) %in% colnames(tBData))]
    for(t in 1:noEpochs){
      if(!is.null(nFold)){
        folds <- caret::createFolds(factor(tmpExtP$Histotype), k = nFold, list = FALSE)
      }else{
        nFold <- 1
      }
      for(k in 1:nFold){
        if(nFold > 1){
          testIndexes <- which(folds == k, 
                               arr.ind=TRUE)
        }else{
          # Alternative approach to k-fold 
          # We utilize stratified sampling here to ensure randomness
          # 70/30 split, with the 30 group used for testing
          # i.e. we remove 30% of the samples randomly, and test on the remainder
          folds <- tmpExtP %>%
            group_by(Histotype) %>%
            slice_sample(prop = 0.7)
          folds <- data.frame(folds)
          testIndexes <- which(tmpExtP$barcode %in% folds$barcode, 
                               arr.ind=TRUE)
        }
        trainData <- data.frame(startData)
        testData <- tBData[testIndexes, ]
        testData <- data.frame(testData)
        testP <- tmpExtP[testIndexes, ]
        # Add Histotype label to DF
        trainData$Histotype <- trainP$Histotype[match(rownames(trainData), trainP$barcode)]
        testData$Histotype <- testP$Histotype[match(rownames(testData), testP$barcode)]
        trainVals <- trainData[,-which(colnames(trainData) %in% "Histotype")]
        testVals <- testData[,-which(colnames(testData) %in% "Histotype")]
        # Create phenotypic annotation
        trainLabel <- trainData$Histotype
        testLabel <- testData$Histotype
        # If other group is missing a category, we add this one as a new category at the end of the model
        modCats <- names(table(trainLabel))
        if(length(which(!names(table(testLabel)) %in% names(table(trainLabel)))) > 0){
          missCats <- names(table(testLabel))[which(!names(table(testLabel)) %in% names(table(trainLabel)))]
          for(o in 1:length(missCats)){
            modCats <- append(modCats, missCats[o])
          }
        }
        repVec <- c()
        for(v in 1:length(modCats)){
          repVec <- append(repVec, v-1)
          names(repVec)[length(repVec)] <- modCats[v]
          trainLabel <- replace(trainLabel, trainLabel==modCats[v], v-1)
        }
        # Replace values in testlabel using repvec
        for(w in 1:length(repVec)){
          testLabel <- ifelse(testLabel %in% names(repVec)[w], repVec[w], testLabel)
        }
        
        # Ratio of classes in multiclass
        class_imba <- table(trainLabel)/length(trainLabel)
        cW <- min(class_imba)/class_imba
        cWArr <- rep(NA, length(trainLabel))
         for(x in 1:length(trainLabel)){
           cWArr[x] <- cW[which(names(cW) %in% trainLabel[x])]
        }
        #inpBoost <- as.matrix(trainVals)
        #inpBoost <- sapply(inpBoost, as.numeric)
        if(isFALSE(paramBool)){
          tmpBoost <- xgboost(x = as.matrix(trainVals),
                              y = trainLabel,
                              booster = "gbtree",
                              objective = "multi:softprob",
                              weights = cWArr,
                              nrounds =100,
                              nthread = 10)
        }else{
          tmpBoost <- xgboost(data = as.matrix(trainVals),
                              label = trainLabel,
                              booster = "gbtree",
                              objective = "multi:softprob",
                              eta =  best_params["eta"] ,
                              max_depth =  best_params["max_depth"],
                              min_child_weight = best_params["min_child_weight"],
                              subsample = best_params["subsample"],
                              colsample_bytree = best_params["colsample_bytree"],
                              nrounds = best_params["n_estimators"],
                              verbose = 0,
                              nthread = 10,
                              num_class = length(modCats))
        }
        tmpExtPred <- predict(object=tmpBoost, 
                                           newdata = as.matrix(testVals))
        colnames(tmpExtPred) <- repVec
        
        # May not be needed
        predicted_labels <- factor(colnames(tmpExtPred)[max.col(tmpExtPred)], levels=repVec)
        testLabelFac <- factor(testLabel, levels=repVec)
        conMat <- caret::confusionMatrix(data = predicted_labels, 
                                         reference = testLabelFac)
        # Correction for instances where we only have 2 classes 
        if(length(levels(testLabelFac))<=2){
          sens <- conMat$byClass[1]
          sens <- append(sens, sens)
          names(sens) <- names(table(testLabel))
          spec <- conMat$byClass[2]
          spec <- append(spec, spec)
          names(spec) <- names(table(testLabel))
          f1Score <- conMat$byClass[7]
          f1Score  <- append(f1Score, f1Score)
          names(f1Score) <- names(table(testLabel))
          acc <- conMat$byClass[11]
          acc  <- append(acc, acc)
          names(acc) <- names(table(testLabel))
        }else{
          sens <- conMat$byClass[,1]
          spec <- conMat$byClass[,2]
          f1Score <- conMat$byClass[,7]
          acc <- conMat$byClass[,11]
          
        }
        # Get overall statistics
        accAll <- conMat$overall[[1]]
        # Get confidence interval of accuracy (lower, upper)
        accLower <- conMat$overall[[3]]
        accUpper <- conMat$overall[[4]]
        
        # Get AUC for model
        rocPred <- tmpExtPred
        colnames(rocPred) <- repVec
        # Get AUC for model
        # If we dont have more then 1 class, we cant get these values and instead assign them to 0 
        if(length(table(testLabel)) > 1){
          # Try testlabelfac instead for the instances with no or few classes 
          roc <- pROC::multiclass.roc(predictor = rocPred,
                                      response= testLabel, 
                                      smooth=FALSE,  
                                      quiet = TRUE)
          aucVal <- roc$auc[1]
          
          # Get PRAUC for model
          prAucDf <- as.data.frame(tmpExtPred)
          prLevFac <- droplevels(testLabelFac)
          prAucDf <- prAucDf[,colnames(prAucDf) %in% levels(prLevFac)]
          prAucDf$ref <-  prLevFac
          prauc <- prAucDf %>% yardstick::pr_auc(estimator = "macro",
                                                 ref,
                                                 min(levels(prLevFac)):max(levels(prLevFac)))
          prAucVal <- prauc$.estimate[1]
          
        }else{
          aucVal <- 0
          prAucVal <- 0
        }
        
        for(u in 1:length(names(table(trainP$Histotype)))){
          tmpH <- names(table(trainP$Histotype))[u]
          if(!tmpH %in% names(repVec)){
            sens <- append(sens, 0, after=u-1)
            names(sens)[u] <- tmpH
            spec <- append(spec, 0, after=u-1)
            names(spec)[u] <- tmpH
            f1Score <- append(f1Score, 0, after=u-1)
            names(f1Score)[u] <- tmpH
            acc <- append(acc, 0, after=u-1)
            names(acc)[u] <- tmpH
          }else{
            names(sens)[u] <- tmpH
            names(spec)[u] <- tmpH
            names(acc)[u] <- tmpH
            names(f1Score)[u] <- tmpH
          }
        }
        
        sens <- sens[names(sens) %in% names(table(trainP$Histotype))]
        spec <- spec[names(spec) %in% names(table(trainP$Histotype))]
        acc <- sens[names(acc) %in% names(table(trainP$Histotype))]
        f1Score <- f1Score[names(f1Score) %in% names(table(trainP$Histotype))]
        
        sumCvDfMult[nrow(sumCvDfMult)+1,] <- c(sens, spec, acc, f1Score,
                                               accAll, accLower, accUpper,
                                               aucVal, prAucVal,
                                               t,k)
        
      }
    }
    sumCvDfMult[sapply(sumCvDfMult, is.nan)] <- 0
    sumCvDfMult[sapply(sumCvDfMult, is.na)] <- 0
    meanItDf <- sumCvDfMult %>%
      group_by(Epoch) %>%
      summarise_all("mean")
    
    sumCvDfMult$SID <- tmpSid
    cvMultLst[[s]] <- sumCvDfMult
    names(cvMultLst)[s] <- tmpSid
    
    meanItDf <- data.frame(meanItDf)
    meanItDf[nrow(meanItDf)+1, ] <- colMeans(meanItDf)
    rownames(meanItDf)[nrow(meanItDf)] <- "MeanVal"
    meanItDf$SID <- tmpSid
    cvMultMeanLst[[s]] <- meanItDf
    names(cvMultMeanLst)[s] <- tmpSid
  }
  
  sumDf <- data.frame(matrix(nrow=0, ncol=ncol(cvMultMeanLst[[1]])))
  colnames(sumDf) <- colnames(cvMultMeanLst[[1]])
  for(v in 1:length(cvMultMeanLst)){
    tmpDf <- cvMultMeanLst[[v]]
    meanRow <- tmpDf["MeanVal",]
    sumDf[nrow(sumDf)+1,] <- unlist(meanRow)
  }
  sumDf <- sumDf[,!colnames(sumDf) %in% c("Epoch", "Fold")]
  # Add cpg-model to dataframe
  sumDf$EModel <- ensModVec
  sumDf$CpGModel <- annCpgModVec
  sumDf$ModelCat <- hCatVec
  resDfLst <- list("Merge"=sumDf, 
                "MeanCvLst"= cvMultMeanLst, 
                "cvAllLst" = cvMultLst, 
                "CpGModel" = annCpgModVec, 
                "ensModel" = ensModVec,
                "modelCat" = hCatVec)
  return(resDfLst)
}

makeDirExtBinPred <- function(inpExtBetaLst, inpExtPhenoLst, inpCModel, refPheno, sigCpgLocs, refBeta, 
                              noEpochs = NULL, nFold = NULL, paramBool = NULL, remHistBool=NULL){
  # Directly apply model trained on our data on external data, binary classification
  if(is.null(noEpochs)){
    noEpochs <- 100
  }
  # cpgModel <-  sapply(inpCModel, function(x) strsplit(x, "\\+")[[1]][1], USE.NAMES=FALSE)
  # ensModel <- sapply(inpEModel, function(x) strsplit(x, "\\+")[[1]][1], USE.NAMES=FALSE)
  # Run multiclass prediction on test-datasets
  extPredLst <- list()
  
  cvAllLst <- list()
  meanCvLst <- list()
  
  cpgModVec <- c()
  annCpgModVec <- c()
  ensModVec <- c()
  missPredVec <- c()
  hCatVec <- c()
  
  # If not null remove histotype and associated samples from model 
  if(!is.null(remHistBool)){
    refBeta <- refBeta[, !colnames(refBeta) %in% trainPheno$barcode[which(trainPheno$Histotype %in% remHistBool)]]
    refPheno <- refPheno[!refPheno$Histotype %in% remHistBool,]
    for(x in 1:length(extBetaLst)){
      tmpB <- inpExtBetaLst[[x]]
      tmpP <- inpExtPhenoLst[[x]]
      tmpB <- tmpB[, !colnames(tmpB) %in% tmpP$barcode[which(tmpP$Histotype %in% remHistBool)]]
      tmpP <- tmpP[!tmpP$Histotype %in% remHistBool,]
      inpExtBetaLst[[x]] <- tmpB
      inpExtPhenoLst[[x]] <- tmpP
    }
  }
  
  cpgModel <- inpCModel
  # Set up starting data for model
  startData <- refBeta[cpgModel,]
  startData <- data.frame(t(startData))
  trainP <- refPheno
  
  if(is.null(paramBool)){
    # Tune parameters for XGB model before we apply it
    tuneData <- data.frame(startData)
    tuneData$Histotype <- trainP$Histotype[match(rownames(tuneData), 
                                                 trainP$barcode)]
    best_params <- makeXGBoostCaretGridSMulti(gridDat = tuneData)
  }else if(!isFALSE(paramBool)){
    best_params <- paramBool
  }
  
  for(s in 1:length(inpExtBetaLst)){
    tmpSid <- names(inpExtBetaLst)[s]
    message(paste("Performing binary classification on cohort: ", tmpSid, sep=""))
    tmpExtB <- inpExtBetaLst[[tmpSid]]
    
    # Control what CpGs are in the dataframe compared to the model
    if(length(which(!cpgModel %in% rownames(tmpExtB))) > 1){
      message(paste("Warning, CpGs: ", paste(cpgModel[which(!cpgModel %in% rownames(tmpExtB))], collapse=", "), 
                    " not found in beta-matrix, and will be excluded from analysis", sep=""))
      remCpg <- cpgModel[which(!cpgModel %in% rownames(tmpExtB))]
      missPredVec <- append(missPredVec, paste(remCpg, collapse="+"))
      cpgModel <- cpgModel[!cpgModel %in% remCpg]
    }else{
      missPredVec <- append(missPredVec, 0)
    }
    names(missPredVec)[length(missPredVec)] <- tmpSid
    tmpExtB <- tmpExtB[cpgModel,]
    # Annotate model
    annMod <- getAnnotatedCpgModel(inpCpg=cpgModel, 
                                   inpPromoLocs=sigCpgLocs)
    # Save model for later use
    cpgModVec <- append(cpgModVec, paste(annMod$CpG, collapse="+"))
    names(cpgModVec)[length(cpgModVec)] <- tmpSid
    # Save annotated model
    annCpgModVec <- append(annCpgModVec, paste(annMod$CpG, collapse="+"))
    names(annCpgModVec)[length(annCpgModVec)] <- tmpSid
    ensModVec <- append(ensModVec, paste(annMod$Ens, collapse="+"))
    names(ensModVec)[length(ensModVec)] <- tmpSid
    hCatVec <-  append(hCatVec, paste(names(annMod$CpG), collapse="+"))
    names(hCatVec)[length(hCatVec)] <- tmpSid
    
    tmpExtP <- inpExtPhenoLst[[tmpSid]]
    # Make sure all samples are present in both pheno and beta
    tmpExtB  <- tmpExtB[,which(colnames(tmpExtB) %in% tmpExtP$barcode)]
    tmpExtP <- tmpExtP[which(tmpExtP$barcode %in% colnames(tmpExtB)),]
    # Create dataframe for storing data
    sumCvDfMult <- data.frame(matrix(nrow=0, ncol=2+(8*(length(table(trainP$Histotype))))))
    colnames(sumCvDfMult) <- c(paste(names(table(trainP$Histotype)), "Sensitivity",sep="_"), 
                               paste(names(table(trainP$Histotype)), "Specificity",sep="_"),
                               paste(names(table(trainP$Histotype)), "Accuracy",sep="_"),
                               paste(names(table(trainP$Histotype)), "F1Score",sep="_"),
                               paste(names(table(trainP$Histotype)), "CI_Lower",sep="_"),
                               paste(names(table(trainP$Histotype)), "CI_Upper",sep="_"),
                               paste(names(table(trainP$Histotype)), "AUC",sep="_"),
                               paste(names(table(trainP$Histotype)), "PRAUC",sep="_"),
                               "Epoch", "Fold")
    # Create dataframe for test-dataset
    tmpExtB <- tmpExtB[annMod$CpG,]
    tmpExtB <- na.omit(tmpExtB)
    tmpExtB <- data.frame(t(tmpExtB))
    # Repeat for train-dataset
    startData <- startData[,which(colnames(startData) %in% colnames(tmpExtB))]
    # Add phenotype column to data
    startData$Histotype <- trainP$Histotype[match(rownames(startData), trainP$barcode)]
    tmpExtB$Histotype <- tmpExtP$Histotype[match(rownames(tmpExtB), tmpExtP$barcode)]
    histKfLst <- list()
    for(i in 1:length(table(trainP$Histotype))){
      trainData <- startData
      testData <- tmpExtB
      tmpH <- names(table(trainP$Histotype))[i]
      tmpCModel <- annMod$CpG[grep(tmpH, names(annMod$CpG))]
      kfDf <- data.frame(matrix(nrow=0, ncol=11))
      colnames(kfDf) <- c("Sensitivity", 
                          "Specificity",
                          "Accuracy",
                          "F1Score",
                          "CI_Lower",
                          "CI_Upper",
                          "AUC",
                          "PRAUC",
                          "Epoch",
                          "Fold",
                          "Histotype")
      if(length(tmpCModel) == 0){
        next()
      }
      trainLabel <- trainData$Histotype
      testLabel <- testData$Histotype
      if(!tmpH %in% testLabel){
        next()
      }
      # Convert to numerical categories
      trainLabel <- ifelse(trainLabel == tmpH, 1, 0)
      trainLabel <- factor(trainLabel, levels = c(0,1))
      testLabel <- ifelse(testLabel == tmpH, 1, 0)
      testLabel <- factor(testLabel, levels = c(0,1))
      # Get values
      trainVals <- trainData[,-which(colnames(trainData) %in% "Histotype")]
      testVals <- testData[,-which(colnames(testData) %in% "Histotype")]
      # Include only values in model for histotype of interest
      trainVals <- trainVals[,tmpCModel]
      testVals <- testVals[,tmpCModel]
      # Loop through t epochs over k folds
      for(t in 1:noEpochs){
        if(!is.null(nFold)){
          folds <- caret::createFolds(testLabel, 
                                     k = nFold, 
                                     list = FALSE)
        }else{
          nFold = 1
        }
        for(k in 1:nFold){
          if(nFold > 1){
            # Segment train-data by folds, with respect to the fold index no. 
            testIndexes <- which(!folds == k, 
                                   arr.ind=TRUE)
          }else{
            # We utilize stratified sampling here to ensure randomness
            # 70/30 split, with the 70 group used for testing
            # i.e. we remove 30% of the samples randomly, and test on the remainder
            folds <- tmpExtP %>%
              group_by(Histotype) %>%
              slice_sample(prop = 0.7)
            folds <- data.frame(folds)
            testIndexes <- which(rownames(testData) %in% folds$barcode, 
                                 arr.ind=TRUE)
          }
          testFold <- testVals[testIndexes, ]
          labelFold <- testLabel[testIndexes]
          dfVec <- c()
          if(isFALSE(paramBool)){
            class_imba <- length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
            tmpBoost <- xgboost(x = as.matrix(trainVals),
                                y = trainLabel,
                                booster = "gbtree",
                                objective = "binary:logistic",
                                scale_pos_weight = class_imba,
                                nrounds = 100,
                                nthread = 10)
          }else{
            tmpBoost <- xgboost(data = as.matrix(trainVals),
                                label = trainLabel,
                                booster = "gbtree",
                                objective = "binary:logistic",
                                eta =  best_params["eta"] ,
                                max_depth =  best_params["max_depth"],
                                min_child_weight = best_params["min_child_weight"],
                                subsample = best_params["subsample"],
                                colsample_bytree = best_params["colsample_bytree"],
                                nrounds = best_params["n_estimators"],
                                verbose = 0,
                                nthread = 10)
          }
          # Apply model on test-data (external dataset)
          tmpExtPred <- predict(object=tmpBoost, 
                                newdata = as.matrix(testFold))
          # Get confusion matrix
          predicted_labels <- factor(ifelse(tmpExtPred > 0.5, 1,0), 
                                     levels=c(1,0))
          conMat <- caret::confusionMatrix(data = predicted_labels, 
                                           reference = factor(labelFold, levels=c(1,0)))
          # Get initial statistics
          sens <- conMat$byClass[1]
          spec <- conMat$byClass[2]
          f1Score <- conMat$byClass[7]
          acc <- conMat$byClass[11]
          # Get overall statistics
          accAll <- conMat$overall[[1]]
          # Get confidence interval of accuracy (lower, upper)
          accLower <- conMat$overall[[3]]
          accUpper <- conMat$overall[[4]]
          # Get AUC for model
          roc <- pROC::roc(labelFold, 
                           tmpExtPred, smooth=FALSE,  quiet = TRUE)
          # Get AUPRC
          # Indice locations of pos and neg-classes
          indice_POS=which(labelFold %in% 1)
          indice_NEG=which(labelFold %in% 0)
          # Get scores of pos and neg classes
          clas_score_POS=tmpExtPred[indice_POS]
          clas_score_NEG=tmpExtPred[indice_NEG]
          # Calculate area under precision recall curve
          # testLabelFreqs <- table(testLabel)/length(testLabel)
          # testLabelW <-ifelse(testLabel %in% names(testLabelFreqs)[1], testLableFreqs[1], testLableFreqs[2])
          # Class 0 is equal to the positive class
          # Class 1 is equal to the negative class
          auprc <- PRROC::pr.curve(scores.class0 = clas_score_POS, 
                                  scores.class1 = clas_score_NEG)
          statVec <- c(sens, spec, acc, f1Score, accLower, accUpper, roc$auc[1], auprc$auc.integral[1], t,k, tmpH)
          names(statVec) <- c("Sensitivity", "Specificity", "Accuracy", "F1Score", 
                              "CI_Lower", "CI_Upper", "AUC", "PRAUC", 
                              "Epoch", "Fold", "Histotype")
          kfDf[nrow(kfDf)+1, ] <- statVec
        }
      }
      histKfLst[[i]] <- kfDf
      names(histKfLst)[i] <- tmpH
    }
    # Merge list of stats
    histKfLst <- histKfLst[!lengths(histKfLst) == 0]
    cvSumKDf <- do.call(cbind.data.frame, histKfLst)
    # Replace Na/NAN with 0
    cvSumKDf[sapply(cvSumKDf, is.nan)] <- 0
    cvSumKDf[sapply(cvSumKDf, is.na)] <- 0
    # Rename columns, match against result DF
    colnames(cvSumKDf) <- gsub("\\.", "_", colnames(cvSumKDf))
    epochDf <- cvSumKDf[,colnames(cvSumKDf) %in% colnames(sumCvDfMult)]
    # Ugly fix due to cbind forcibly changing the value from numeric to character
    for(r in 1:ncol(epochDf)){
      epochDf[,r] <- as.numeric(epochDf[,r])
    }
    # Summarize statistics based on fold
    epochDf$Epoch <- as.numeric(cvSumKDf[,grep("Epoch", colnames(cvSumKDf))[1]])
    # Get dataframe with ALL stats
    cvSumOutDf <- epochDf
    cvSumOutDf$Fold <- as.numeric(cvSumKDf[,grep("Fold", colnames(cvSumKDf))[1]])
    cvSumOutDf$SID <- tmpSid
    
    # Get dataframe with mean for each Epoch
    foldMeans <- epochDf %>%
      group_by(Epoch) %>%
      summarise_all("mean")  
    foldMeans <- data.frame(foldMeans)
    foldMeans[nrow(foldMeans)+1, ] <- colMeans(foldMeans)
    rownames(foldMeans)[nrow(foldMeans)] <- "MeanVal"
    foldMeans$SID <- tmpSid
    # Save to lists
    cvAllLst[[s]] <- epochDf
    names(cvAllLst)[s] <- tmpSid
    meanCvLst[[s]] <- foldMeans
    names(meanCvLst)[s] <- tmpSid
  }
  sumDf <- data.frame(matrix(nrow=length(meanCvLst), ncol=ncol(meanCvLst[[1]])))
  colnames(sumDf) <- colnames(meanCvLst[[1]])
  for(v in 1:length(meanCvLst)){
    tmpDf <- meanCvLst[[v]]
    colIdx <- match(colnames(tmpDf), colnames(sumDf))
    meanRow <- tmpDf["MeanVal",]
    sumDf[v, colIdx] <- unlist(meanRow)
  }
  sumDf <- sumDf[,!colnames(sumDf) %in% c("Epoch", "Fold")]
  # Add cpg-model to dataframe
  sumDf$EModel <- ensModVec
  sumDf$CpGModel <- annCpgModVec
  sumDf$ModelCat <- hCatVec
  return(list("Merge"= sumDf, 
              "MeanCvLst"= meanCvLst, 
              "cvAllLst" = cvAllLst, 
              "CpGModel" = annCpgModVec, 
              "ensModel" = ensModVec,
              "modelCat" = hCatVec))
}

makeSingleGenePred_MULT <- function(inpCpgDfLst, refBeta, refPheno, inpCpgPos, 
                                    inpCpgLocs=NULL, noEpochs = NULL, nFold = NULL, 
                                    paramBool = NULL, rmHist = NULL, minCpg = NULL, 
                                    fullPromoBool=NULL, nCatBool = NULL){
  # Function for predictive classification using a defined genomic region of interest 
  if(is.null(nFold)){
    nFold <- 5
  }
  if(is.null(noEpochs)){
    noEpochs <- 1
  }
  if(is.null(minCpg)){
    minCpg <- 2
  }
  if(!is.null(rmHist)){
    refPheno <- refPheno[!refPheno$Histotype %in% rmHist, ]
    refBeta <- refBeta[,refPheno$barcode]
  }
  refPheno <- refPheno[which(refPheno$barcode %in% colnames(refBeta)),]
  refBeta <- refBeta[,which(colnames(refBeta) %in% refPheno$barcode)]
  refM <- log2(refBeta/(1-refBeta))
  genePredLst <- list()
  for(j in 1:length(inpCpgDfLst)){
    # Second version, perform pred-class on individiual genes
    genePredDf <- data.frame(matrix(nrow=0, ncol=10))
    colnames(genePredDf) <- c("ensembl_gene_id", 
                              "Sensitivity", "Specificity","Accuracy", "F1",
                              "CI_Lower", "CI_Upper",
                              "AUC","AUPRC", 
                              "CpGs")
    inpCpgDf <- inpCpgDfLst[[j]]
    inpHist <- names(inpCpgDfLst)[j]
    # Specify the number of phenotype categories the site should be significant against
    if(is.null(nCatBool)){
      nCat <- length(table(refPheno$Histotype))-1
    }else{
      nCat <- nCatBool
    }
    if(is.null(nrow(inpCpgDf))){
      next()
    }else if(nrow(inpCpgDf) == 0){
      next()
    }
    for(k in 1:nrow(inpCpgDf)){
      if(abs(k)%%round(nrow(inpCpgDf)/10) == 0){
        message(paste("Processing row: [", k, "/", nrow(inpCpgDf),"]"))
      }
      tmpCpgDf <- data.frame(inpCpgDf[k,])
      tmpEns <- tmpCpgDf$Gene
      # Prepare reference dataframe
      promoDf <- makeCpgDf(tmpCpgDf, 
                           refBeta, 
                           refPheno, 
                           inpHist)
      if(is.null(nrow(promoDf))){
        message(paste("Gene: ", tmpEns, " Not found in array-data and will be skipped", sep=""))
        next()
      }else if(nrow(promoDf) <= minCpg){
        message(paste("Gene: ", tmpEns, " Has less then or equal to: ",  minCpg, " total CpG sites and will be skipped", sep=""))
        next()
      }else{
        promoDf <- promoDf[rownames(promoDf) %in% rownames(refBeta),]
        # Either get the significant CpGs, or model the full promoter region
        if(!is.null(fullPromoBool)){
          minDists <- nrow(promoDf)
        }else{
          minDists <- makeRegSigCpg(inpReg = promoDf,
                   inpBeta = refBeta,
                   inpPheno = refPheno,
                   inpH = inpHist,
                   noCats = nCat,
                   minCpg = minCpg)
        }
        if(length(minDists) < 2){
          next()
        }else{
          # Prepare training-data  based on the less stringent filtering
          refPredDf <- refBeta[minDists, ]
          refPredDf <- as.data.frame(t(refPredDf))
          refPredDf  <- refPredDf[rowSums(is.na(refPredDf)) == 0,]
          refPredDf$Histotype <- refPheno$Histotype[match(rownames(refPredDf), refPheno$barcode)]
          refPredDf$Histotype <- ifelse(refPredDf$Histotype %in% inpHist, inpHist, "Other")
          refPredDf$Histotype <- factor(refPredDf$Histotype, levels = c("Other", inpHist))
          # Create models to be used with data
          minModel <- paste(minDists, collapse = " + ")
          minModel <- reformulate(minModel, "Histotype")
          # Shuffle rows to remove any bias when sampling
          refPredDf <- refPredDf[sample(nrow(refPredDf)),]
          kfDfLst <- list()
          if(!is.null(paramBool)){
            # Perform grid-search to identify best parameters for model
            # Long runtime, should only be performed when necessary
            best_params <- makeXGBoostCaretGridS(refPredDf)
          }
          inpRef <- refPredDf
          epochLst <- list()
          # Perform multicore classification
          epochLst <- future_lapply(1:noEpochs, future.seed=TRUE, function(l){
            folds <- caret::createFolds(factor(inpRef$Histotype), 
                                        k = nFold, 
                                        list = FALSE) 
            # predStatLst <- list()
            kfDf <- data.frame(matrix(nrow=0, ncol=8))
            colnames(kfDf) <- c("Sensitivity", "Specificity", "Accuracy", "F1",
                                "CI_Lower","CI_Upper", 
                                "AUC", "AUPRC")
            
            for(m in 1:nFold){
              # Segment train-data by n folds, with respect to the fold index no. 
              testIndexes <- which(folds==m, arr.ind=TRUE)
              trainData <- inpRef[-testIndexes, ]
              testData <- inpRef[testIndexes, ]
              # best_params <- makeXGBoostCaretGridS(trainData)
              # Create weights for random forest classification (balance)
              trainVals <- trainData[,-which(grepl("Histotype", colnames(trainData)))]
              testVals <- testData[,-which(grepl("Histotype", colnames(testData)))]
              # Set up labels for test and training data
              # XGBoost now requires binary labels to be in factor instead of numeric
              # Set positive class (1) and negative (0), i.e. Yes=1 and No=0 as "Other" is the largest class 
              trainLabel <- trainData[,which(grepl("Histotype",colnames(trainData)))]
              trainLabel <- ifelse(trainLabel == inpHist, 1, 0)
              trainLabel <- factor(trainLabel, levels = c(0,1))
              testLabel <- testData[,which(grepl("Histotype",colnames(testData)))]
              testLabel <- ifelse(testLabel == inpHist, 1, 0)
              testLabel <- factor(testLabel, levels = c(0,1))
              # Get the ratio of negative (0) to positive (1) classes
              class_imba <- length(which(trainLabel %in% 0))/length(which(trainLabel %in% 1))
              
              # Build gradient boosting forest model
              if(is.null(paramBool)){
                tmpBoost <- xgboost::xgboost(x = as.matrix(trainVals),
                                             y = trainLabel,
                                             booster = "gbtree",
                                             objective = "binary:logistic",
                                             scale_pos_weight = class_imba,
                                             nrounds = 100,
                                             nthread = 10)
              }else{
                tmpBoost <- xgboost(data = as.matrix(trainVals),
                                    label = trainLabel,
                                    booster = "gbtree",
                                    objective = "binary:logistic",
                                    eta =  best_params["eta"] ,
                                    max_depth =  best_params["max_depth"],
                                    scale_pos_weight = class_imba,
                                    min_child_weight = best_params["min_child_weight"],
                                    subsample = best_params["subsample"],
                                    colsample_bytree = best_params["colsample_bytree"],
                                    nrounds = best_params["n_estimators"],
                                    verbose = 0,
                                    nthread = 10)
              }
              # Apply RF model on test-data (external dataset)
              # XGBoosts prediction class returns the probabilities of a sample belonging to the last class in binary classification (i.e. class 1)
              # This is confusing when viewing the confusionMatrix, as the classes will seem to be "swapped"
              tmpExtPred <- predict(object=tmpBoost, 
                                    newdata = as.matrix(testVals),
                                    type="response")
              # Categorize predictions based on their corresponding classification seen to model
              conMat <- confusionMatrix(data = factor(as.numeric(tmpExtPred > 0.5), levels = c(1,0)), 
                                        reference = factor(testLabel, levels = c(1,0)))
              # As we are dealing with heavily imbalanced data, optimize on f1 score and AUPRC
              sens <- conMat$byClass[[1]]
              spec <- conMat$byClass[[2]]
              f1Score <- conMat$byClass[[7]]
              acc <- conMat$overall[[1]]
              # Get confidence interval of statistics
              confIntLower <- conMat$overall[[3]]
              confIntUpper <- conMat$overall[[4]]
              # Get AUC for model
              roc <- pROC::roc(testLabel, 
                               tmpExtPred, 
                               smooth=FALSE,  
                               quiet = TRUE)
              # Get AUPRC
              # Indice locations of pos and neg-classes
              indice_POS=which(testLabel %in% 1)
              indice_NEG=which(testLabel %in% 0)
              # Get scores of pos and neg classes
              clas_score_POS=tmpExtPred[indice_POS]
              clas_score_NEG=tmpExtPred[indice_NEG]
              # Calculate area under precision recall curve
              # testLabelFreqs <- table(testLabel)/length(testLabel)
              # testLabelW <-ifelse(testLabel %in% names(testLabelFreqs)[1], testLableFreqs[1], testLableFreqs[2])
              auprc <-PRROC::pr.curve(scores.class0 = clas_score_POS,
                                      scores.class1 = clas_score_NEG)
              # Save to list/DF
              kfDf[nrow(kfDf)+1,] <- c(sens, spec, acc, f1Score, 
                                       confIntLower, confIntUpper, 
                                       roc$auc[1], auprc$auc.integral[1])
              rownames(kfDf)[nrow(kfDf)] <- paste("k=",m, sep="")
              
            }
            kfDf[sapply(kfDf, is.nan)] <- 0
            kfDf[sapply(kfDf, is.na)] <- 0
            kfDf[nrow(kfDf)+1, ] <- colMeans(kfDf)
            epochLst[[l]]<- colMeans(kfDf)
          })
          sumCvDf2 <- do.call(rbind.data.frame, epochLst)
          colnames(sumCvDf2) <- c("Sensitivity", "Specificity", "Accuracy","F1", 
                                  "CI_Lower", "CI_Upper", 
                                  "AUC", "AUPRC")
          meanPred <- colMeans(sumCvDf2)
          genePredDf[nrow(genePredDf)+1,] <- c(tmpEns, 
                                               meanPred[1], meanPred[2], meanPred[3], meanPred[4],
                                               meanPred[5], meanPred[6],
                                               meanPred[7], meanPred[8],
                                               paste(minDists, collapse = ","))
        }
      }
    }
    # Select only models with a good combined AUC
    # genePredDf <- genePredDf[which(genePredDf$F1 > 0.7), ]
    genePredLst[[length(genePredLst)+1]] <- genePredDf
    names(genePredLst)[length(genePredLst)] <- inpHist
  }
  message("Individual gene-pred completed")
  return(genePredLst)
}
  
  ##############################################################################
  # Deprecated
  ##############################################################################
  
# Directly apply model trained on our data on external data
makeDirExtMultPred_OLD <- function(inpExtBetaLst, inpExtPhenoLst, inpCModel, refPheno, sigCpgLocs, refBeta, noEpochs = NULL, nFold = NULL, paramBool = NULL, remHistBool=NULL){
    if(is.null(noEpochs)){
      noEpochs <- 100
    }
    if(is.null(nFold)){
      nFold <- 5
    }
    # cpgModel <-  sapply(inpCModel, function(x) strsplit(x, "\\+")[[1]][1], USE.NAMES=FALSE)
    # ensModel <- sapply(inpEModel, function(x) strsplit(x, "\\+")[[1]][1], USE.NAMES=FALSE)
    # Run multiclass prediction on test-datasets
    extPredLst <- list()
    cvMultLst <- list()
    cvMultMeanLst <- list()
    cpgModVec <- c()
    annCpgModVec <- c()
    ensModVec <- c()
    missPredVec <- c()
    hCatVec <- c()
    
    # If not null remove histotype and associated samples from model 
    if(!is.null(remHistBool)){
      refBeta<- refBeta[, !colnames(refBeta) %in% trainPheno$barcode[which(trainPheno$Histotype %in% remHistBool)]]
      refPheno <- refPheno[!refPheno$Histotype %in% remHistBool,]
      for(x in 1:length(extBetaLst)){
        tmpB <- inpExtBetaLst[[x]]
        tmpP <- inpExtPhenoLst[[x]]
        tmpB <- tmpB[, !colnames(tmpB) %in% tmpP$barcode[which(tmpP$Histotype %in% remHistBool)]]
        tmpP <- tmpP[!tmpP$Histotype %in% remHistBool,]
        inpExtBetaLst[[x]] <- tmpB
        inpExtPhenoLst[[x]] <- tmpP
      }
    }
    
    # Set up starting data for model
    startData <- refBeta[inpCModel,]
    startData <- t(startData)
    trainP <- refPheno
    
    if(is.null(paramBool)){
      # Tune parameters for XGB model before we apply it
      tuneData <- data.frame(startData)
      tuneData$Histotype <- trainP$Histotype[match(rownames(tuneData), trainP$barcode)]
      best_params <- makeXGBoostCaretGridSMulti(gridDat = tuneData)
    }else if(!isFALSE(paramBool)){
      best_params <- paramBool
    }
    
    for(s in 1:length(inpExtBetaLst)){
      cpgModel <- inpCModel
      tmpSid <- names(inpExtBetaLst)[s]
      message(paste("Performing multiclass classification on cohort: ", tmpSid, sep=""))
      tmpExtB <- inpExtBetaLst[[tmpSid]]
      # Control what CpGs are in the dataframe compared to the model
      if(length(which(!cpgModel %in% rownames(tmpExtB))) > 0 | length(which(!cpgModel %in% rownames(refBeta))) > 0){
        # Get missing CpG's
        refMiss <- paste(cpgModel[which(!cpgModel %in% rownames(refBeta))], collapse=", ")
        extMiss <- paste(cpgModel[which(!cpgModel %in% rownames(tmpExtB))], collapse=", ")
        allMiss <- c(refMiss, extMiss)
        allMiss <- allMiss[!allMiss %in% ""]
        message(paste("Warning, CpGs: ", paste(allMiss,collapse=", "), 
                      " not found in external beta-matrix, and will be excluded from analysis", sep=""))
        missPredVec <- append(missPredVec, paste(allMiss, collapse="+"))
        cpgModel <- cpgModel[!cpgModel %in% allMiss]
        startData <- startData[,colnames(startData) %in% cpgModel]
      }else{
        missPredVec <- append(missPredVec, 0)
      }
      names(missPredVec)[length(missPredVec)] <- tmpSid
      
      tmpExtB <- tmpExtB[cpgModel,]
      # Save model for later use
      cpgModVec <- append(cpgModVec, paste(cpgModel, collapse="+"))
      names(cpgModVec)[length(cpgModVec)] <- tmpSid
      # Annotate model
      annMod <- getAnnotatedCpgModel(inpCpg = cpgModel, 
                                     inpPromoLocs = sigCpgLocs)
      # Save annotated model
      annCpgModVec <- append(annCpgModVec, paste(annMod$CpG, collapse="+"))
      names(annCpgModVec)[length(annCpgModVec)] <- tmpSid
      
      ensModVec <- append(ensModVec, paste(annMod$Ens, collapse="+"))
      names(ensModVec)[length(ensModVec)] <- tmpSid
      
      hCatVec <-  append(hCatVec, paste(names(annMod$CpG), collapse="+"))
      names(hCatVec)[length(hCatVec)] <- tmpSid
      
      tmpExtP <- inpExtPhenoLst[[tmpSid]]
      # Make sure all samples are present in both pheno and beta
      tmpExtB  <- tmpExtB[,which(colnames(tmpExtB) %in% tmpExtP$barcode)]
      tmpExtP <- tmpExtP[which(tmpExtP$barcode %in% colnames(tmpExtB)),]
      # Create dataframe for storing data
      sumCvDfMult <- data.frame(matrix(nrow=0, ncol=7+(4*(length(table(trainP$Histotype))))))
      colnames(sumCvDfMult) <- c(paste(names(table(trainP$Histotype)), "Sensitivity",sep="_"), 
                                 paste(names(table(trainP$Histotype)), "Specificity",sep="_"),
                                 paste(names(table(trainP$Histotype)), "Accuracy",sep="_"),
                                 paste(names(table(trainP$Histotype)), "F1Score",sep="_"),
                                 "Acc_Overall","Acc_Lower", "Acc_Upper",
                                 "AUC", "AUPRC", "Epoch", "Fold")
      # Prepare data for training the model
      # Choose between training the model on our cohort and applying it directly on ext,
      for(t in 1:noEpochs){
        folds <- caret::createFolds(factor(tmpExtP$Histotype), k = nFold, list = FALSE)
        for(k in 1:nFold){
          testIndexes <- which(folds==k, arr.ind=TRUE)
          trainData <- data.frame(startData)
          # Training data remains constant (i.e. the full source dataset)
          # Due to low sample size, we use the entire test dataset
          # To ensure variation but maximize size, we instead but remove the holdout set samples
          # i.e. The non-fold samples are used as test set, instead of the smaller number of fold-samples
          testData <- tmpExtB[, -testIndexes]
          testData <- data.frame(t(testData))
          testP <- tmpExtP[-testIndexes, ]
          # Add Histotype label to DF
          trainData$Histotype <- trainP$Histotype[match(rownames(trainData), trainP$barcode)]
          testData$Histotype <- testP$Histotype[match(rownames(testData), testP$barcode)]
          trainVals <- trainData[,-which(colnames(trainData) %in% "Histotype")]
          testVals <- testData[,-which(colnames(testData) %in% "Histotype")]
          # Create phenotypic annotation
          trainLabel <- trainData$Histotype
          testLabel <- testData$Histotype
          
          # If other group is missing a category, we add this one as a new category at the end of the model
          modCats <- names(table(trainLabel))
          if(length(which(!names(table(testLabel)) %in% names(table(trainLabel)))) > 0){
            missCats <- names(table(testLabel))[which(!names(table(testLabel)) %in% names(table(trainLabel)))]
            for(o in 1:length(missCats)){
              modCats <- append(modCats, missCats[o])
            }
          }
          
          repVec <- c()
          for(v in 1:length(modCats)){
            repVec <- append(repVec, v-1)
            names(repVec)[length(repVec)] <- modCats[v]
            trainLabel <- replace(trainLabel, trainLabel==modCats[v], v-1)
          }
          # Replace values in testlabel using repvec
          for(w in 1:length(repVec)){
            testLabel <- ifelse(testLabel %in% names(repVec)[w], repVec[w], testLabel)
          }
          # # Ratio of classes in multiclass
          # class_imba <- table(trainLabel)/length(trainLabel)
          # cW <- min(class_imba)/class_imba
          # cWArr <- rep(NA, length(trainLabel))
          # for(x in 1:length(trainLabel)){
          #   cWArr[x] <- cW[which(names(cW) %in% trainLabel[x])]
          # }
          #inpBoost <- as.matrix(trainVals)
          #inpBoost <- sapply(inpBoost, as.numeric)
          if(isFALSE(paramBool)){
            tmpBoost <- xgboost(data = as.matrix(trainVals),
                                label = trainLabel,
                                booster = "gbtree",
                                objective = "multi:softprob",
                                nrounds =100,
                                verbose = 0,
                                nthread = 10,
                                num_class = length(modCats))
          }else{
            tmpBoost <- xgboost(data = as.matrix(trainVals),
                                label = trainLabel,
                                booster = "gbtree",
                                objective = "multi:softprob",
                                eta =  best_params["eta"] ,
                                max_depth =  best_params["max_depth"],
                                min_child_weight = best_params["min_child_weight"],
                                subsample = best_params["subsample"],
                                colsample_bytree = best_params["colsample_bytree"],
                                nrounds = best_params["n_estimators"],
                                verbose = 0,
                                nthread = 10,
                                num_class = length(modCats))
          }
          tmpExtPred <- predict(object=tmpBoost, 
                                newdata = as.matrix(testVals))
          colnames(tmpExtPred) <- repVec
          # May not be needed
          predicted_labels <- factor(colnames(tmpExtPred)[max.col(tmpExtPred)], levels=repVec)
          testLabelFac <- factor(testLabel, levels=repVec)
          conMat <- caret::confusionMatrix(data = predicted_labels, 
                                           reference = testLabelFac)
          # Correction for instances where we only have 2 classes 
          if(length(levels(testLabelFac))<=2){
            sens <- conMat$byClass[1]
            sens <- append(sens, sens)
            names(sens) <- names(table(testLabel))
            spec <- conMat$byClass[2]
            spec <- append(spec, spec)
            names(spec) <- names(table(testLabel))
            f1Score <- conMat$byClass[7]
            f1Score  <- append(f1Score, f1Score)
            names(f1Score) <- names(table(testLabel))
            acc <- conMat$byClass[11]
            acc  <- append(acc, acc)
            names(acc) <- names(table(testLabel))
          }else{
            sens <- conMat$byClass[,1]
            spec <- conMat$byClass[,2]
            f1Score <- conMat$byClass[,7]
            acc <- conMat$byClass[,11]
            
          }
          # Get overall statistics
          accAll <- conMat$overall[[1]]
          # Get confidence interval of accuracy (lower, upper)
          accLower <- conMat$overall[[3]]
          accUpper <- conMat$overall[[4]]
          
          # Get AUC for model
          rocPred <- tmpExtPred
          colnames(rocPred) <- repVec
          # Get AUC for model
          # If we dont have more then 1 class, we cant get these values and instead assign them to 0 
          if(length(table(testLabel)) > 1){
            # Try testlabelfac instead for the instances with no or few classes 
            roc <- pROC::multiclass.roc(predictor = rocPred,
                                        response= testLabel, 
                                        smooth=FALSE,  
                                        quiet = TRUE)
            aucVal <- roc$auc[1]
            
            # Get PRAUC for model
            prAucDf <- as.data.frame(tmpExtPred)
            prLevFac <- droplevels(testLabelFac)
            prAucDf <- prAucDf[,colnames(prAucDf) %in% levels(prLevFac)]
            prAucDf$ref <-  prLevFac
            prauc <- prAucDf %>% yardstick::pr_auc(estimator = "macro",
                                                   ref,
                                                   min(levels(prLevFac)):max(levels(prLevFac)))
            prAucVal <- prauc$.estimate[1]
            
          }else{
            aucVal <- 0
            prAucVal <- 0
          }
          
          for(u in 1:length(names(table(trainP$Histotype)))){
            tmpH <- names(table(trainP$Histotype))[u]
            if(!tmpH %in% names(repVec)){
              sens <- append(sens, 0, after=u-1)
              names(sens)[u] <- tmpH
              spec <- append(spec, 0, after=u-1)
              names(spec)[u] <- tmpH
              f1Score <- append(f1Score, 0, after=u-1)
              names(f1Score)[u] <- tmpH
              acc <- append(acc, 0, after=u-1)
              names(acc)[u] <- tmpH
            }else{
              names(sens)[u] <- tmpH
              names(spec)[u] <- tmpH
              names(acc)[u] <- tmpH
              names(f1Score)[u] <- tmpH
            }
          }
          
          sens <- sens[names(sens) %in% names(table(trainP$Histotype))]
          spec <- spec[names(spec) %in% names(table(trainP$Histotype))]
          acc <- sens[names(acc) %in% names(table(trainP$Histotype))]
          f1Score <- f1Score[names(f1Score) %in% names(table(trainP$Histotype))]
          
          sumCvDfMult[nrow(sumCvDfMult)+1,] <- c(sens, spec, acc, f1Score,
                                                 accAll, accLower, accUpper,
                                                 aucVal, prAucVal,
                                                 t,k)
          
        }
      }
      sumCvDfMult[sapply(sumCvDfMult, is.nan)] <- 0
      sumCvDfMult[sapply(sumCvDfMult, is.na)] <- 0
      meanItDf <- sumCvDfMult %>%
        group_by(Epoch) %>%
        summarise_all("mean")
      
      sumCvDfMult$SID <- tmpSid
      cvMultLst[[s]] <- sumCvDfMult
      names(cvMultLst)[s] <- tmpSid
      
      meanItDf <- data.frame(meanItDf)
      meanItDf[nrow(meanItDf)+1, ] <- colMeans(meanItDf)
      rownames(meanItDf)[nrow(meanItDf)] <- "MeanVal"
      meanItDf$SID <- tmpSid
      cvMultMeanLst[[s]] <- meanItDf
      names(cvMultMeanLst)[s] <- tmpSid
    }
    
    sumDf <- data.frame(matrix(nrow=0, ncol=ncol(cvMultMeanLst[[1]])))
    colnames(sumDf) <- colnames(cvMultMeanLst[[1]])
    for(v in 1:length(cvMultMeanLst)){
      tmpDf <- cvMultMeanLst[[v]]
      meanRow <- tmpDf["MeanVal",]
      sumDf[nrow(sumDf)+1,] <- unlist(meanRow)
    }
    sumDf <- sumDf[,!colnames(sumDf) %in% c("Epoch", "Fold")]
    # Add cpg-model to dataframe
    sumDf$EModel <- ensModVec
    sumDf$CpGModel <- annCpgModVec
    sumDf$ModelCat <- hCatVec
    return(list("Merge"=sumDf, 
                "MeanCvLst"= cvMultMeanLst, 
                "cvAllLst" = cvMultLst, 
                "CpGModel" = annCpgModVec, 
                "ensModel" = ensModVec,
                "modelCat" = hCatVec))
}
