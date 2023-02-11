

#		-- predictorIndices: indicates the rows of the matrix that should be considered as predictors for a regression problem.
#		-- targetIndices: indicates which genes should be considered as targets for a regression problem.
#	Parameters (optional):
# 		
#		-- rankThreshold[DEFAULT=round((length(predictorIndices)/20))] : the amount of top ranked features that should be awarded 1 instead of 0 during a feature ranking in an iteration step
#		-- ensembleSize[DEFAULT=2000] : the amount of models in the ensemble
#		-- minSampleSize[DEFAULT=round((dim(expressionMatrix)[1])/5)] : mininimum sample size of the experiments
#		-- maxSampleSize[DEFAULT=4*round((dim(expressionMatrix)[1])/5)] : maximum sample size of the experiments
#		-- predictorSampleSizeMin[DEFAULT=round(length(predictorIndices)/5)] : mininimum sample size of predictors
#		-- predictorSampleSizeMax[DEFAULT=4*round(length(predictorIndices)/5)] : maximum sample size of the predictors
#		-- alpha[DEFAULT=0.3]: the mixing parameter of the elastic net. Range [0=L2norm,1=L1norm].
#		-- trace[DEFAULT=TRUE]: index of target currently computed is reported to stdout
#       -- traceDetail[DEFAULT=FALSE] : index of subproblen is reported to stdout
#		-- ... : the other parameters are passed to the glmnet functions
#
#
#	Returns:								
#		
#		A matrix where each cell(row i, column j) specifies the importance of variable i to target j.

#' ElVariableEnsembleSolve
#' Solves the feature selection problem using L0 regularization.
#'
#' @param expressionMatrix Expression matrix as returned by "constructExpressionMatrixFromFile".
#'                         Format: samples x genes ; colnames are unique gene labels ; rownames could be ignored
#' @param predictorIndices 
#' @param targetIndices 
#' @param rankThreshold 
#' @param ensembleSize 
#' @param minSampleSize 
#' @param maxSampleSize 
#' @param predictorSampleSizeMin 
#' @param predictorSampleSizeMax 
#' @param alpha 
#' @param trace 
#' @param traceDetail 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
ElVariableEnsembleSolve <- function(expressionMatrix,
                                    predictorIndices,
                                    targetIndices,
                                    rankThreshold = round((length(predictorIndices) / 20)),
                                    ensembleSize = 2000,
                                    minSampleSize = round((dim(expressionMatrix)[1]) / 5),
                                    maxSampleSize = 4 * round((dim(expressionMatrix)[1]) / 5),
                                    predictorSampleSizeMin = round(length(predictorIndices) / 5),
                                    predictorSampleSizeMax = 4 * round(length(predictorIndices) / 5),
                                    alpha = 0.3,
                                    trace = TRUE,
                                    traceDetail = FALSE,
                                    ...) {
  
  # Check input
  if ((minSampleSize <= 0 || maxSampleSize > dim(expressionMatrix)[1] || maxSampleSize < minSampleSize)) {
    stop("Please specify a valid sample sample size minimum and maximum.")
  }
  
  # Preallocate the return matrix
  resultMatrix <- matrix(0.0, length(predictorIndices), length(targetIndices))
  rownames(resultMatrix) <- colnames(expressionMatrix)[predictorIndices]
  colnames(resultMatrix) <- colnames(expressionMatrix)[targetIndices]
  
  # Take a sample this many times
  for (i in 1:ensembleSize) {
    
    # sampling
    if (minSampleSize == maxSampleSize) {
      sampleSize <- minSampleSize
    } else {
      sampleSize <- sample(minSampleSize:maxSampleSize, 1)
    }
    
    # Take sample of specified size 
    ind <- sample(dim(expressionMatrix)[1], sampleSize)
    sampleMatrix <- expressionMatrix[ind,]
    
    # report on progress
    if (trace) {
      print(paste("Sample: ", i, " of size: ", sampleSize, " out of ", ensembleSize, " iterations and a rank threshold of ", rankThreshold))
    }
    
    # Call EL procedure
    resultMatrix <- resultMatrix + elasticNetRankedSolve(sampleMatrix, 
                                                         predictorIndices = predictorIndices, 
                                                         targetIndices = targetIndices, 
                                                         alpha = alpha, 
                                                         traceDetail = traceDetail, 
                                                         rankThreshold = rankThreshold, 
                                                         predictorSampleSizeMin = predictorSampleSizeMin, 
                                                         predictorSampleSizeMax = predictorSampleSizeMax, 
                                                         ...)
    
  }
  # return
  return(resultMatrix)
  
}

# 	elasticNetRankedSolve
#
#	
# 	Solves the feature selection problem using the Elastic Net/LASSO/Ridge Regression using a rank threshold.
#
#
#   This method should not be called directly but by using ElVariableEnsembleSolve.
#
#	Parameters	(required):
#		-- expressionMatrix: expression matrix as returned by constructExpressionMatrixFromFile.
#							 Format: samples x genes ; colnames are unique gene labels ; rownames are ignored
#		-- predictorIndices: indicates the rows of the matrix that should be considered as predictors for a regression problem.
#		-- targetIndices   : indicates which genes should be considered as targets for a regression problem.
#		-- rankThreshold : the amount of top ranked features that should be awarded 1 instead of 0 during a feature ranking in an iteration step
#		-- predictorSampleSizeMin : mininimum sample size of predictors
#		-- predictorSampleSizeMax : maximum sample size of the predictors
#		-- alpha: the mixing parameter of the elastic net. Range [0=L2norm,1=L1norm].
#		-- traceDetail: index of target currently computed is reported to stdout				
#	
#	Parameters (optional):
# 	

#		-- ... : the other parameters are passed to the glmnet functions
#
#
#Returns:
#	A matrix where each cell(row i, column j) specifies the importance of variable i to target j.
#

elasticNetRankedSolve <- function(expressionMatrix,
                                  predictorIndices,
                                  targetIndices,
                                  alpha,
                                  rankThreshold,
                                  traceDetail,
                                  predictorSampleSizeMin,
                                  predictorSampleSizeMax,
                                  penalty = "L0",
                                  ...) {
  
  
  # Check arguments
  if ((predictorSampleSizeMin <= 0 || predictorSampleSizeMax >= length(predictorIndices) || predictorSampleSizeMin > predictorSampleSizeMax)) {
    stop("Please specify a valid predictor sample size minimum and maximum (between 1 inclusive and the length of the predictorIndices exclusive).")
  }
  
  # Only the predictors should be evaluated
  predictors <- expressionMatrix[, predictorIndices]
  # Pre allocation
  resultMatrix <- matrix(0.0, length(colnames(predictors)), length(targetIndices))
  rownames(resultMatrix) <- colnames(expressionMatrix)[predictorIndices]
  colnames(resultMatrix) <- colnames(expressionMatrix)[targetIndices]
  
  # Loop over all possible targets
  i <- 0
  for (targetIndex in targetIndices) {
    # Report on progress 
    if (traceDetail) {
      i <- i + 1
      cat(paste("Computing target ", targetIndex, " iteration ", i, "\\", length(targetIndices), "\n"))
      flush.console()
    }
    # extract the target
    target <- expressionMatrix[, targetIndex, drop = FALSE]
    # Check if the target has already been removed from the predictors, if not, remove it
    if (targetIndex %in% predictorIndices) {
      nameCol <- colnames(expressionMatrix)[targetIndex]
      indexInPredictorsToRemove <- which(nameCol == colnames(predictors))
      predictorsWithoutTarget <- predictors[, - indexInPredictorsToRemove, drop = FALSE]
    }
    else {
      predictorsWithoutTarget <- predictors
    }
    
    # Sample on possible predictors
    if (predictorSampleSizeMin == predictorSampleSizeMax) {
      predictorSampleSize <- predictorSampleSizeMin
      pIndices <- sample(1:(dim(predictorsWithoutTarget)[2]), predictorSampleSize)
      predictorsWithoutTarget <- predictorsWithoutTarget[, pIndices, drop = FALSE]
    }
    else {
      predictorSampleSize <- sample(predictorSampleSizeMin:predictorSampleSizeMax, 1)
      pIndices <- sample(1:(dim(predictorsWithoutTarget)[2]), predictorSampleSize)
      predictorsWithoutTarget <- predictorsWithoutTarget[, pIndices, drop = FALSE]
    }
    
    if (F) {
      nfolds <- round(sqrt(dim(predictorsWithoutTarget)[1]))
      # Cross validate
      cross <- glmnet::cv.glmnet(predictorsWithoutTarget,
                                 target,
                                 nfolds = nfolds,
                                 alpha = alpha,
                                 family = "gaussian",
                                 standardize = FALSE,
                                 ...)
      # Retrain optimal
      bestModel <- glmnet::glmnet(predictorsWithoutTarget, target, family = "gaussian", alpha = alpha, lambda = cross$lambda.min, standardize = FALSE, ...)
      # get the weights
      wghts <- abs(as.vector(bestModel$beta))
    } else {
      # message("----- Run ", penalty, " model for ", names(target), "! -----")
      L0_Model <- L0Learn::L0Learn.fit(predictorsWithoutTarget,
                                       target,
                                       penalty = penalty,
                                       maxSuppSize = ncol(predictorsWithoutTarget))
      L0_Model_Information <- as.data.frame(print(L0_Model))
      L0_Model_Information <- L0_Model_Information[order(L0_Model_Information$suppSize,
                                                         decreasing = TRUE), ]
      lambda_L0 <- L0_Model_Information$lambda[1]
      gamma_L0 <- L0_Model_Information$gamma[1]
      temp <- coef(L0_Model,
                   lambda = lambda_L0,
                   gamma = gamma_L0)
      
      # suppressPackageStartupMessages(
      #   temp <- LO_fit(predictorsWithoutTarget,
      #                  target,
      #                  penalty = penalty,
      #                  nFolds = 10,
      #                  seed = 1,
      #                  maxSuppSize = ncol(predictorsWithoutTarget),
      #                  nGamma = 5,
      #                  gammaMin = 0.0001,
      #                  gammaMax = 10
      #   )
      # )
      
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      wghts <- wghts / sum(wghts)
    }
    
    # Now sort the wghts
    indices <- sort.list(wghts, decreasing = TRUE)
    # Check for zero entries
    zeros <- which(wghts == 0)
    # Now replace by ones that are in the top and are non-zero
    wghts[1:length(wghts)] <- 0
    wghts[indices[1:rankThreshold]] <- 1
    # Set the ones that were zero to zero anyway
    wghts[zeros] <- 0
    
    # write result to matrix
    resultMatrix[colnames(predictorsWithoutTarget), colnames(expressionMatrix)[targetIndex]] <- wghts
    
  }
  return(resultMatrix)
}
