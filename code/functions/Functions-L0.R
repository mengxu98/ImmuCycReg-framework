#' L0VariableEnsembleSolve
#'  Solves the feature selection problem using L0 regularization.
#'
#' @param expressionMatrix Expression matrix as returned by "constructExpressionMatrixFromFile"
#'                         Format: samples x genes ; colnames are unique gene labels ; rownames could be ignored
#' @param predictorIndices Indicates the rows of the matrix that should be considered as predictors for a regression problem
#' @param targetIndices Indicates which genes should be considered as targets for a regression problem
#' @param rankThreshold [DEFAULT=round((length(predictorIndices)/20))]
#'    The amount of top ranked features that should be awarded 1 instead of 0 during a feature ranking in an iteration step
#' @param ensembleSize [DEFAULT=2000] The amount of models in the ensemble
#' @param minSampleSize [DEFAULT=round((dim(expressionMatrix)[1])/5)] Mininimum sample size of the experiments
#' @param maxSampleSize [DEFAULT=4*round((dim(expressionMatrix)[1])/5)] Maximum sample size of the experiments
#' @param predictorSampleSizeMin [DEFAULT=round(length(predictorIndices)/5)] Mininimum sample size of predictors
#' @param predictorSampleSizeMax [DEFAULT=4*round(length(predictorIndices)/5)] Maximum sample size of the predictors
#' @param trace [DEFAULT=TRUE] Index of target currently computed is reported to stdout
#' @param traceDetail [DEFAULT=FALSE] Index of subproblen is reported to stdout
#' @param ...
#'
#' @return A matrix where each cell (row i, column j) specifies the importance of variable i to target j
#' @export
#'
L0VariableEnsembleSolve <- function(expressionMatrix,
                                    predictorIndices,
                                    targetIndices,
                                    rankThreshold = round((length(predictorIndices) / 20)),
                                    ensembleSize = 2000,
                                    minSampleSize = round((dim(expressionMatrix)[1]) / 5),
                                    maxSampleSize = 4 * round((dim(expressionMatrix)[1]) / 5),
                                    predictorSampleSizeMin = round(length(predictorIndices) / 5),
                                    predictorSampleSizeMax = 4 * round(length(predictorIndices) / 5),
                                    trace = TRUE,
                                    traceDetail = FALSE,
                                    ...) {
  
  # Check input
  if ((minSampleSize <= 0 || maxSampleSize > dim(expressionMatrix)[1] || maxSampleSize < minSampleSize)) {
    stop("Please specify a valid sample sample size minimum and maximum......")
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
    if(trace){
      message(paste(
        "Sampling:", i, "of", ensembleSize,
        "times, and select", sampleSize, "of", nrow(expressionMatrix),"samples......"
      ))
    }
    
    # Call procedure
    resultMatrix <- resultMatrix + L0RankedSolve(sampleMatrix, 
                                                 predictorIndices = predictorIndices, 
                                                 targetIndices = targetIndices,
                                                 traceDetail = traceDetail, 
                                                 rankThreshold = rankThreshold, 
                                                 predictorSampleSizeMin = predictorSampleSizeMin, 
                                                 predictorSampleSizeMax = predictorSampleSizeMax)
    
  }
  return(resultMatrix)
}

#' L0RankedSolve
#'
#' @param expressionMatrix Expression matrix as returned by constructExpressionMatrixFromFile
#'  Format: samples x genes ; colnames are unique gene labels ; rownames are ignored
#' @param predictorIndices Indicates the rows of the matrix that should be considered as predictors for a regression problem
#' @param targetIndices Indicates which genes should be considered as targets for a regression problem
#' @param rankThreshold The amount of top ranked features that should be awarded 1 instead of 0 during a feature ranking in an iteration step
#' @param traceDetail 
#' @param predictorSampleSizeMin Mininimum sample size of predictors
#' @param predictorSampleSizeMax Maximum sample size of the predictors
#' @param penalty 
#' @param ...
#'
#' @return A matrix where each cell(row i, column j) specifies the importance of variable i to target j
#' @export
#'
L0RankedSolve <- function(expressionMatrix,
                          predictorIndices,
                          targetIndices,
                          rankThreshold,
                          traceDetail,
                          predictorSampleSizeMin,
                          predictorSampleSizeMax,
                          penalty = "L0") {
  # Check arguments
  if ((predictorSampleSizeMin <= 0 || predictorSampleSizeMax >= length(predictorIndices) || predictorSampleSizeMin > predictorSampleSizeMax)) {
    stop("Please specify a valid predictor sample size minimum and maximum......")
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
    
    # Extract the target
    target <- expressionMatrix[, targetIndex, drop = FALSE]
    # Check if the target has already been removed from the predictors, if not, remove it
    if (targetIndex %in% predictorIndices) {
      nameCol <- colnames(expressionMatrix)[targetIndex]
      indexInPredictorsToRemove <- which(nameCol == colnames(predictors))
      predictorsWithoutTarget <- predictors[, - indexInPredictorsToRemove, drop = FALSE]
    } else {
      predictorsWithoutTarget <- predictors
    }
    
    # Sample on possible predictors
    if (predictorSampleSizeMin == predictorSampleSizeMax) {
      predictorSampleSize <- predictorSampleSizeMin
      pIndices <- sample(1:dim(predictorsWithoutTarget)[2], predictorSampleSize)
      predictorsWithoutTarget <- predictorsWithoutTarget[, pIndices, drop = FALSE]
    } else {
      predictorSampleSize <- sample(predictorSampleSizeMin:predictorSampleSizeMax, 1)
      pIndices <- sample(1:dim(predictorsWithoutTarget)[2], predictorSampleSize)
      predictorsWithoutTarget <- predictorsWithoutTarget[, pIndices, drop = FALSE]
    }
    
    L0Model <- inferCSN::inferCSN.fit(predictorsWithoutTarget,
                                      target,
                                      penalty = penalty,
                                      maxSuppSize = ncol(predictorsWithoutTarget))
    L0ModelInfor <- as.data.frame(print(L0Model))
    L0ModelInfor <- L0ModelInfor[order(L0ModelInfor$suppSize,
                                       decreasing = TRUE), ]
    
    temp <- coef(L0Model,
                 lambda = L0ModelInfor$lambda[1],
                 gamma = L0ModelInfor$gamma[1]) %>% as.vector()
    
    wghts <- temp[-1] %>% abs()
    wghts <- wghts / sum(wghts)
    
    # Now sort the wghts
    indices <- sort.list(wghts, decreasing = TRUE)
    # Check for zero entries
    zeros <- which(wghts == 0)
    # Now replace by ones that are in the top and are non-zero
    wghts[1:length(wghts)] <- 0
    wghts[indices[1:rankThreshold]] <- 1
    # Set the ones that were zero to zero anyway
    wghts[zeros] <- 0
    
    # Write result to matrix
    resultMatrix[colnames(predictorsWithoutTarget),
                 colnames(expressionMatrix)[targetIndex]] <- wghts
  }
  return(resultMatrix)
}
