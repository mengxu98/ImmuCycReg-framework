

#' L0DWGRN
#'
#' @param matrix The rows are samples and the columns are genes of the matrix
#' @param penalty
#' @param regulators
#' @param targets
#' @param maxSuppSize
#' @param cores
#' @param crossInteraction
#' @param ensembleSize 
#'
#' @return
#' @export
#'
#' @examples
L0DWGRN <- function(matrix,
                    penalty = NULL,
                    regulators = NULL,
                    targets = NULL,
                    maxSuppSize = NULL,
                    cores = 1,
                    crossInteraction = FALSE,
                    ensembleSize = 100) {
  matrix <- as.data.frame(t(matrix))
  geneNum <- dim(matrix)[2]

  if (is.null(penalty)) penalty <- "L0"
  if (is.null(maxSuppSize)) maxSuppSize <- geneNum - 1
  if (is.null(targets)) targets <- colnames(matrix)
  if (is.null(regulators)) {
    regulators <- colnames(matrix)
  } else {
    matrix <- matrix[, regulators]
  }

  if (crossInteraction == FALSE) {
    weightList <- L0DWGRN_core(matrix,
      penalty = penalty,
      regulators = regulators,
      targets = targets,
      maxSuppSize = maxSuppSize,
      cores = cores
    )
  } else {
    weightList <- L0DWGRN_cross(matrix,
      penalty = penalty,
      cores = cores,
      ensembleSize = ensembleSize
    )
  }
  return(weightList)
}

#' L0DWGRN_core
#'
#' @param matrix The rows are samples and the columns are genes of the matrix
#' @param penalty
#' @param regulators
#' @param targets
#' @param maxSuppSize
#' @param cores 
#'
#' @return
#' @export
#'
#' @examples
L0DWGRN_core <- function(matrix,
                         penalty = penalty,
                         regulators = regulators,
                         targets = targets,
                         maxSuppSize = maxSuppSize,
                         cores = 1) {
  if (cores == 1) {
    weightList <- c()
    for (i in 1:length(regulators)) {
      X <- as.matrix(matrix[, -which(colnames(matrix) == regulators[i])])
      Y <- matrix[, regulators[i]]
      temp <- LO_fit(X, Y,
                     penalty = penalty,
                     nFolds = 10,
                     seed = 1,
                     maxSuppSize = maxSuppSize,
                     nGamma = 5,
                     gammaMin = 0.0001,
                     gammaMax = 10
      )
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      wghts <- wghts / sum(wghts)
      if (F) {
        wghts <- wghts / max(wghts)
        indices <- sort.list(wghts, decreasing = TRUE)
        zeros <- which(wghts <= 0.8)
        # wghts[1:length(wghts)] <- 1
        wghts[zeros] <- 0
      }
      if (length(wghts) != dim(X)[2]) {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = 0)
      } else {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = wghts)
      }
      weightList <- rbind.data.frame(weightList, weightd)
      if (i == length(regulators)) {
        weightList <- weightList[order(weightList$weight, decreasing = TRUE), ]
      }
    }
  } else {
    cores <- min(parallel::detectCores(logical = F), cores)
    # cl <- parallel::makeCluster(cores)
    # doParallel::registerDoParallel(cl)
    doParallel::registerDoParallel(cores = cores)
    message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
    "%dopar%" <- foreach::"%dopar%"
    suppressPackageStartupMessages(
      weightList <- doRNG::"%dorng%"(foreach::foreach(regulator = regulators, .combine = "rbind", .export = "LO_fit"), {
        X <- as.matrix(matrix[, -which(colnames(matrix) == regulator)])
        Y <- matrix[, regulator]
        temp <- LO_fit(X, Y,
                       penalty = penalty,
                       nFolds = 10,
                       seed = 1,
                       maxSuppSize = maxSuppSize,
                       nGamma = 5,
                       gammaMin = 0.0001,
                       gammaMax = 10
        )
        temp <- as.vector(temp)
        wghts <- temp[-1]
        wghts <- abs(wghts)
        wghts <- wghts / sum(wghts)
        
        if (F) {
          # indices <- sort.list(wghts, decreasing = TRUE)
          # zeros <- which(wghts <= 0.8)
          # # wghts[1:length(wghts)] <- 1
          # wghts[zeros] <- 0
          
          # Now sort the wghts
          indices <- sort.list(wghts, decreasing = TRUE)
          # Check for zero entries
          zeros <- which(wghts == 0)
          # Now replace by ones that are in the top and are non-zero
          wghts[1:length(wghts)] <- 0
          rankThreshold <- 5
          wghts[indices[1:rankThreshold]] <- 1
          # Set the ones that were zero to zero anyway
          wghts[zeros] <- 0
        }
        
        if (length(wghts) != dim(X)[2]) {
          weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulator, weight = 0)
        } else {
          weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulator, weight = wghts)
        }
      })
    )
    # parallel::stopCluster(cl)
  }
  weightList <- weightList[order(weightList$weight, decreasing = TRUE), ]
  return(weightList)
}


#' Title
#'
#' @param matrix
#' @param regulatoryGene
#' @param targetGene
#' @param trace
#' @param predSampleMin
#' @param predSampleMax
#' @param sampleMin
#' @param sampleMax
#' @param rankThreshold
#' @param ensembleSize
#' @param candidateSplitCount
#' @param penalty 
#' @param cores 
#'
#' @return
#' @export
#'
#' @examples
L0DWGRN_cross <- function(matrix,
                          regulatoryGene = NULL,
                          targetGene = NULL,
                          trace = TRUE,
                          predSampleMin = 20,
                          predSampleMax = 80,
                          sampleMin = 20,
                          sampleMax = 80,
                          rankThreshold = 5,
                          ensembleSize = ensembleSize,
                          candidateSplitCount = "ALL",
                          penalty = penalty,
                          cores = cores) {
  
  # Check arguments
  # The following arguments should be specified as percentages
  if (!(isPercentage(predSampleMin) &&
    isPercentage(predSampleMax) &&
    isPercentage(sampleMin) &&
    isPercentage(sampleMax) &&
    isPercentage(rankThreshold))) {
    stop("Sample sizes and rank thresholds should be specified as percentages.")
  }
  
  # Ensemble size should be positive and non zero
  if (!(isWholeNumber(ensembleSize) && ensembleSize > 0)) {
    stop("Ensemble sizes should be specified as non-zero integers.")
  }

  # Set the regulatoryGene to all if they are not set
  if (is.null(regulatoryGene)) {
    regulatoryGene <- getIndicesOfGenesInMatrix(matrix)
  }
  
  # Set targetGene to all if they are are not set
  if (is.null(targetGene)) {
    targetGene <- getIndicesOfGenesInMatrix(matrix)
  }

  # Convert the percentages to the numbers
  lenPredictors <- dim(matrix)[2]
  lenSamples <- dim(matrix)[1]
  predSampleMin <- round((predSampleMin * lenPredictors) / 100)
  predSampleMax <- round((predSampleMax * lenPredictors) / 100)
  sampleMin <- round((sampleMin * lenSamples) / 100)
  sampleMax <- round((sampleMax * lenSamples) / 100)
  rankThreshold <- round((rankThreshold * lenPredictors) / 100)

  # Preallocate the return matrix
  resultMatrix <- matrix(0.0, length(regulatoryGene), length(targetGene))
  rownames(resultMatrix) <- colnames(matrix)[regulatoryGene]
  colnames(resultMatrix) <- colnames(matrix)[targetGene]
  resultMatrix <- sortMatrix(resultMatrix)

  weightList <- L0_fit_cross(matrix,
    regulatoryGene,
    targetGene,
    trace = trace,
    ensembleSize = ensembleSize,
    rankThreshold = rankThreshold,
    predSampleMinSize = predSampleMin,
    predSampleMaxSize = predSampleMax,
    minSampleSize = sampleMin,
    maxSampleSize = sampleMax,
    penalty = penalty,
    cores = cores
  )

  weightList <- matrix2List(inputFile = weightList)
  weightList[, 3] <- 1:(dim(weightList)[1])
  weightList <- list2Matrix(input = weightList)
  weightList <- sortMatrix(weightList)
  resultMatrix <- resultMatrix + weightList
  resultMatrix <- matrix2List(inputFile = resultMatrix, desc = FALSE)
}


#' LO_fit
#'
#' @param X The rows are samples and the columns are genes of the matrix
#' @param Y
#' @param penalty
#' @param nFolds
#' @param seed
#' @param maxSuppSize
#' @param nGamma
#' @param gammaMin
#' @param gammaMax
#'
#' @return
#' @export
#'
#' @examples
LO_fit <- function(X, Y,
                   penalty = penalty,
                   nFolds = 10,
                   seed = 1,
                   maxSuppSize = maxSuppSize,
                   nGamma = 5,
                   gammaMin = 0.0001,
                   gammaMax = 10) {
  tryCatch(
    {
      fit <- L0Learn::L0Learn.cvfit(X, Y,
        penalty = penalty,
        maxSuppSize = maxSuppSize,
        nFolds = 10,
        seed = 1,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
      )
      fit_inf <- print(fit)
      optimalGammaIndex <- which(unlist(lapply(fit$cvMeans, min)) == min(unlist(lapply(fit$cvMeans, min))))
      gamma <- fit$fit$gamma[optimalGammaIndex]
      lambda_list <- fit_inf[which(fit_inf$gamma == gamma), ]
      if (is.null(maxSuppSize)) {
        lambda <- min(lambda_list$lambda)
      } else {
        if (maxSuppSize %in% lambda_list$maxSuppSize) {
          lambda <- lambda_list$maxSuppSize[which(lambda_list$maxSuppSize == maxSuppSize)]
        } else {
          lambda <- min(lambda_list$lambda)
        }
      }
      temp <- coef(fit,
        lambda = lambda,
        gamma = gamma
      )
    },
    error = function(e) {
      fit <- L0Learn::L0Learn.fit(X, Y,
        penalty = penalty,
        maxSuppSize = maxSuppSize,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
      )
      fit_inf <- print(fit)
      fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
      lambda <- fit_inf$lambda[1]
      gamma <- fit_inf$gamma[1]
      temp <- coef(fit,
        lambda = lambda,
        gamma = gamma
      )
    }
  )
}

# LO_fit <- function(X, Y,
#                    penalty = penalty,
#                    nFolds = 10,
#                    seed = 1,
#                    maxSuppSize = maxSuppSize,
#                    nGamma = 5,
#                    gammaMin = 0.0001,
#                    gammaMax = 10) {
#
#       fit <- L0Learn::L0Learn.fit(X, Y,
#                                   penalty = penalty,
#                                   maxSuppSize = maxSuppSize,
#                                   nGamma = 5,
#                                   gammaMin = 0.0001,
#                                   gammaMax = 10
#       )
#       fit_inf <- print(fit)
#       fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
#       lambda <- fit_inf$lambda[1]
#       gamma <- fit_inf$gamma[1]
#       temp <- coef(fit,
#                    lambda = lambda,
#                    gamma = gamma
#
#   )
# }


#' L0_fit_cross
#' Solves the feature selection problem using L0 regularization.
#'
#' @param matrix Expression matrix as returned by "constructmatrixFromFile".
#'               Format: samples X genes; colnames are unique gene labels; rownames could be ignored
#' @param regulatoryGene Indicates the rows of the matrix that should be considered as predictors for a regression problem.
#' @param targetGene Indicates which genes should be considered as targets for a regression problem.
#' @param rankThreshold [DEFAULT=round((length(regulatoryGene)/20))] The amount of top ranked features that should be awarded 1 instead of 0 during a feature ranking in an iteration step
#' @param ensembleSize [DEFAULT=2000] The amount of models in the ensemble
#' @param minSampleSize [DEFAULT=round((dim(matrix)[1])/5)] Mininimum sample size of the experiments
#' @param maxSampleSize [DEFAULT=4*round((dim(matrix)[1])/5)] Maximum sample size of the experiments
#' @param predSampleMinSize [DEFAULT=round(length(regulatoryGene)/5)] Mininimum sample size of predictors
#' @param predSampleMaxSize [DEFAULT=4*round(length(regulatoryGene)/5)] Maximum sample size of the predictors
#' @param trace [DEFAULT=TRUE] Index of target currently computed is reported to stdout
#' @param penalty 
#' @param cores 
#' @param traceDetail [DEFAULT=FALSE] Index of subproblen is reported to stdout
#'
#' @return A matrix where each cell(row i, column j) specifies the importance of variable i to target j.
#' @export
#'
#' @examples
L0_fit_cross <- function(matrix,
                         regulatoryGene,
                         targetGene,
                         rankThreshold = round((length(regulatoryGene) / 20)),
                         ensembleSize = ensembleSize,
                         minSampleSize = round((dim(matrix)[1]) / 5),
                         maxSampleSize = 4 * round((dim(matrix)[1]) / 5),
                         predSampleMinSize = round(length(regulatoryGene) / 5),
                         predSampleMaxSize = 4 * round(length(regulatoryGene) / 5),
                         trace = trace,
                         traceDetail = traceDetail,
                         penalty = penalty,
                         cores = cores) {
  # Check input
  if ((minSampleSize <= 0 || 
       maxSampleSize > dim(matrix)[1] || 
       maxSampleSize < minSampleSize)) {
    stop("Please specify a valid sample sample size minimum and maximum.")
  }

  # Preallocate the return matrix
  resultMatrix <- matrix(0.0, length(regulatoryGene), length(targetGene))
  rownames(resultMatrix) <- colnames(matrix)[regulatoryGene]
  colnames(resultMatrix) <- colnames(matrix)[targetGene]

  # Multiple sampling
  for (i in 1:ensembleSize) {
    # sampling
    if (minSampleSize == maxSampleSize) {
      sampleSize <- minSampleSize
    } else {
      sampleSize <- sample(minSampleSize:maxSampleSize, 1)
    }

    # Take sample of specified size
    sampleMatrix <- matrix[sample(dim(matrix)[1], sampleSize), ]

    # report on progress
    if (trace) {
      message(paste0(
        "Sample: ",
        i,
        " of size: ",
        sampleSize,
        " out of ",
        ensembleSize,
        " iterations and a rank threshold of ",
        rankThreshold
      ))
    }

    resultMatrix <- resultMatrix + L0_fit_cross_core(sampleMatrix,
      regulatoryGene = regulatoryGene,
      targetGene = targetGene,
      rankThreshold = rankThreshold,
      traceDetail = traceDetail,
      predSampleMinSize = predSampleMinSize,
      predSampleMaxSize = predSampleMaxSize,
      penalty = penalty
    )
  }
  return(resultMatrix)
}


L0_fit_cross_core <- function(matrix,
                              regulatoryGene,
                              targetGene,
                              rankThreshold,
                              traceDetail,
                              predSampleMinSize,
                              predSampleMaxSize,
                              penalty) {
  # Check arguments
  if ((predSampleMinSize <= 0 || 
       predSampleMaxSize >= length(regulatoryGene) || 
       predSampleMinSize > predSampleMaxSize)) {
    stop("Please specify a valid predictor sample size minimum and maximum 
         (between 1 inclusive and the length of the regulatoryGene exclusive).")
  }
  
  # Only the predictors should be evaluated
  predictors <- matrix[, regulatoryGene]
  # Pre allocation
  resultMatrix <- matrix(0.0, length(colnames(predictors)), length(targetGene))
  rownames(resultMatrix) <- colnames(matrix)[regulatoryGene]
  colnames(resultMatrix) <- colnames(matrix)[targetGene]
  
  # Loop over all possible targets
  i <- 0
  for (targetIndex in targetGene) {
    # Report on progress
    if (traceDetail) {
      i <- i + 1
      cat(paste("Computing target ", targetIndex, " iteration ", i, "\\", length(targetGene), "\n"))
      flush.console()
    }
    # extract the target
    matrix <- as.matrix(matrix)
    Y <- matrix[, targetIndex, drop = FALSE]
    
    # Check if the target has already been removed from the predictors, if not, remove it
    if (targetIndex %in% regulatoryGene) {
      nameCol <- colnames(matrix)[targetIndex]
      X <- predictors[, -which(nameCol == colnames(predictors)), drop = FALSE]
    } else {
      X <- predictors
    }
    
    # Sample on possible predictors
    if (predSampleMinSize == predSampleMaxSize) {
      predictorSampleSize <- predSampleMinSize
      pIndices <- sample(1:(dim(X)[2]), predictorSampleSize)
      X <- X[, pIndices, drop = FALSE]
    } else {
      predictorSampleSize <- sample(predSampleMinSize:predSampleMaxSize, 1)
      pIndices <- sample(1:(dim(X)[2]), predictorSampleSize)
      X <- X[, pIndices, drop = FALSE]
    }
    
    L0_Model <- L0Learn::L0Learn.fit(as.matrix(X),
                                     Y,
                                     penalty = penalty,
                                     maxSuppSize = ncol(X)
    )
    L0_Model_Information <- as.data.frame(print(L0_Model))
    L0_Model_Information <- L0_Model_Information[order(L0_Model_Information$suppSize,
                                                       decreasing = TRUE
    ), ]
    lambda_L0 <- L0_Model_Information$lambda[1]
    gamma_L0 <- L0_Model_Information$gamma[1]
    temp <- coef(L0_Model,
                 lambda = lambda_L0,
                 gamma = gamma_L0
    )
    
    # suppressPackageStartupMessages(
    #   temp <- LO_fit(X,
    #                  target,
    #                  penalty = penalty,
    #                  nFolds = 10,
    #                  seed = 1,
    #                  maxSuppSize = ncol(X),
    #                  nGamma = 5,
    #                  gammaMin = 0.0001,
    #                  gammaMax = 10
    #   )
    # )
    
    
    temp <- as.vector(temp)
    wghts <- temp[-1]
    wghts <- abs(wghts)
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
    
    # write result to matrix
    resultMatrix[colnames(X), colnames(matrix)[targetIndex]] <- wghts
  }
  weightList <- matrix2List(inputFile = resultMatrix)
  # weightList[, 3] <- 1:(dim(weightList)[1])
  
  return(weightList)
}




#' L0_fit_cross_core
#'  Solves the feature selection problem using the Elastic Net/LASSO/Ridge Regression using a rank threshold.
#'
#' @param matrix Expression matrix as returned by constructmatrixFromFile.
#'               Format: samples x genes ; colnames are unique gene labels ; rownames are ignored
#' @param regulatoryGene Indicates the rows of the matrix that should be considered as predictors for a regression problem.
#' @param targetGene Indicates which genes should be considered as targets for a regression problem.
#' @param rankThreshold The amount of top ranked features that should be awarded 1 instead of 0 during a feature ranking in an iteration step
#' @param traceDetail index of target currently computed is reported to stdout
#' @param predSampleMinSize mininimum sample size of predictors
#' @param predSampleMaxSize maximum sample size of the predictors
#' @param penalty
#' @param ...
#'
#' @return A matrix where each cell(row i, column j) specifies the importance of variable i to target j.
#' @export
#'
#' @examples
L0_fit_cross_core1 <- function(matrix,
                              regulatoryGene,
                              targetGene,
                              rankThreshold,
                              traceDetail,
                              predSampleMinSize,
                              predSampleMaxSize,
                              penalty = penalty,
                              cores = cores,
                              ...) {
  # Check arguments
  if ((predSampleMinSize <= 0 || predSampleMaxSize >= length(regulatoryGene) || predSampleMinSize > predSampleMaxSize)) {
    stop("Please specify a valid predictor sample size minimum and maximum (between 1 inclusive and the length of the regulatoryGene exclusive).")
  }

  # Only the predictors should be evaluated
  predictors <- matrix[, regulatoryGene]
  # Pre allocation
  resultMatrix <- matrix(0.0, length(colnames(predictors)), length(targetGene))
  rownames(resultMatrix) <- colnames(matrix)[regulatoryGene]
  colnames(resultMatrix) <- colnames(matrix)[targetGene]
  
  if (cores == 1) {
    # Loop over all possible targets
    i <- 0
    for (targetIndex in targetGene) {
      # Report on progress
      if (traceDetail) {
        i <- i + 1
        cat(paste("Computing target ", targetIndex, " iteration ", i, "\\", length(targetGene), "\n"))
        flush.console()
      }
      # extract the target
      matrix <- as.matrix(matrix)
      Y <- matrix[, targetIndex, drop = FALSE]

      # Check if the target has already been removed from the predictors, if not, remove it
      if (targetIndex %in% regulatoryGene) {
        nameCol <- colnames(matrix)[targetIndex]
        X <- predictors[, -which(nameCol == colnames(predictors)), drop = FALSE]
      } else {
        X <- predictors
      }

      # Sample on possible predictors
      if (predSampleMinSize == predSampleMaxSize) {
        predictorSampleSize <- predSampleMinSize
        pIndices <- sample(1:(dim(X)[2]), predictorSampleSize)
        X <- X[, pIndices, drop = FALSE]
      } else {
        predictorSampleSize <- sample(predSampleMinSize:predSampleMaxSize, 1)
        pIndices <- sample(1:(dim(X)[2]), predictorSampleSize)
        X <- X[, pIndices, drop = FALSE]
      }

      L0_Model <- L0Learn::L0Learn.fit(as.matrix(X),
        Y,
        penalty = penalty,
        maxSuppSize = ncol(X)
      )
      L0_Model_Information <- as.data.frame(print(L0_Model))
      L0_Model_Information <- L0_Model_Information[order(L0_Model_Information$suppSize, decreasing = TRUE), ]
      lambda_L0 <- L0_Model_Information$lambda[1]
      gamma_L0 <- L0_Model_Information$gamma[1]
      temp <- coef(L0_Model,
        lambda = lambda_L0,
        gamma = gamma_L0
      )
      
      # suppressPackageStartupMessages(
      #   temp <- LO_fit(X,
      #                  target,
      #                  penalty = penalty,
      #                  nFolds = 10,
      #                  seed = 1,
      #                  maxSuppSize = ncol(X),
      #                  nGamma = 5,
      #                  gammaMin = 0.0001,
      #                  gammaMax = 10
      #   )
      # )


      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
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

      # write result to matrix
      resultMatrix[colnames(X), colnames(matrix)[targetIndex]] <- wghts
    }
  } else {
    cores <- min(parallel::detectCores(logical = F), cores)
    # cl <- parallel::makeCluster(cores)
    # doParallel::registerDoParallel(cl)
    doParallel::registerDoParallel(cores = cores)
    message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
    "%dopar%" <- foreach::"%dopar%"

    # for (targetIndex in targetGene) {
    #   # Report on progress
    #   if (traceDetail) {
    #     i <- i + 1
    #     cat(paste("Computing target ", targetIndex, " iteration ", i, "\\", length(targetGene), "\n"))
    #     flush.console()
    #   }

    suppressPackageStartupMessages(
      resultMatrix <- doRNG::"%dorng%"(foreach::foreach(
        targetIndex = targetGene #,
        # .export = "LO_fit",
        # .combine = "cbind"
      ), {
        # extract the target
        matrix <- as.matrix(matrix)
        Y <- matrix[, targetIndex, drop = FALSE]

        # Check if the target has already been removed from the predictors, if not, remove it
        if (targetIndex %in% regulatoryGene) {
          nameCol <- colnames(matrix)[targetIndex]
          X <- predictors[, -which(nameCol == colnames(predictors)), drop = FALSE]
        } else {
          X <- predictors
        }

        # Sample on possible predictors
        if (predSampleMinSize == predSampleMaxSize) {
          predictorSampleSize <- predSampleMinSize
          pIndices <- sample(1:(dim(X)[2]), predictorSampleSize)
          X <- X[, pIndices, drop = FALSE]
        } else {
          predictorSampleSize <- sample(predSampleMinSize:predSampleMaxSize, 1)
          pIndices <- sample(1:(dim(X)[2]), predictorSampleSize)
          X <- X[, pIndices, drop = FALSE]
        }

        L0_Model <- L0Learn::L0Learn.fit(as.matrix(X),
          Y,
          penalty = penalty,
          maxSuppSize = ncol(X)
        )
        L0_Model_Information <- as.data.frame(print(L0_Model))
        L0_Model_Information <- L0_Model_Information[order(L0_Model_Information$suppSize,
          decreasing = TRUE
        ), ]
        lambda_L0 <- L0_Model_Information$lambda[1]
        gamma_L0 <- L0_Model_Information$gamma[1]
        temp <- coef(L0_Model,
          lambda = lambda_L0,
          gamma = gamma_L0
        )

        # suppressPackageStartupMessages(
        #   temp <- LO_fit(X,
        #                  target,
        #                  penalty = penalty,
        #                  nFolds = 10,
        #                  seed = 1,
        #                  maxSuppSize = ncol(X),
        #                  nGamma = 5,
        #                  gammaMin = 0.0001,
        #                  gammaMax = 10
        #   )
        # )
        
        
        temp <- as.vector(temp)
        wghts <- temp[-1]
        wghts <- abs(wghts)
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

        # write result to matrix
        # resultMatrix[colnames(X), colnames(matrix)[targetIndex]] <- wghts
        setNames(list(setNames(wghts, colnames(X))), colnames(matrix)[targetIndex])
      })
    )
    attr(resultMatrix, "rng") <- NULL
    attr(resultMatrix, "doRNG_version") <- NULL
    resultMatrix <- unlist(resultMatrix, recursive=FALSE)

    weightMatrix <- dplyr::bind_rows(resultMatrix)
    weightMatrix <- as.matrix(weightMatrix)
    rownames(weightMatrix) <- names(resultMatrix)
    weightMatrix <- t(weightMatrix)
    weightMatrix[which(is.na(weightMatrix), arr.ind = TRUE)] <- 0 # optional?  
    resultMatrix <- weightMatrix
    
  }
  return(resultMatrix)
}


# matrix2List
#
# 	Converts an adjacency matrix (genes x genes)  to a TSV format where each line
#   is formatted as [GENE_X][tab][GENE_Y][tab][RANK].
# 		-- inputFilename: name or full path to the file containing the matrix with the adjancency scores.
# 		-- ouputFilename: name or full path to the outputfile
#
# 	Parameters (optional):
#
# 	Returns: NULL
matrix2List <- function(inputFilename = NULL, inputFile = NULL, outputFilename = NULL, desc = TRUE) {
  if (is.null(inputFile)) {
    # Read the table from file
    tbl <- read.table(inputFilename)
  } else {
    tbl <- inputFile
  }
  # First column -> repeat the predictors
  firstCol <- rep(rownames(tbl), length(colnames(tbl)))
  # Second column -> repeat the targets
  secondCol <- c()
  for (i in 1:length(colnames(tbl))) {
    secondCol <- append(secondCol, rep(colnames(tbl)[i], length(rownames(tbl))))
  }
  thirdCol <- as.matrix(tbl)
  thirdCol <- as.vector(thirdCol)
  # Gets the indices from a desc sort on the adjancy measures
  if (desc) {
    or <- sort.list(-thirdCol)
  } else {
    or <- sort.list(thirdCol)
  }
  # Glue everything together
  result <- cbind(firstCol, secondCol, thirdCol)
  # Convert to dataframe
  result <- as.data.frame(result)
  # Sort it using the indices obtained before
  result <- result[or, , drop = FALSE]
  # Write to file if filename is given
  if (!is.null(outputFilename)) {
    write.table(result, file = outputFilename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  # else write to function output
  else {
    return(result)
  }
}





# 	list2Matrix
#
# 	Reads a sorted TSV format from file. Each line is formatted as [GENE_X][tab][GENE_Y][tab][Adjacency_Measure]. The lines are sorted by decreasing adjancency scores. This
# 	format is frequently used as standard for DREAM challenge evaluation scripts (http://wiki.c2b2.columbia.edu/dream/index.php/The_DREAM_Project).
# 	Returns a file formatted as an adjacency matrix (genes x genes), the colums and rows are sorted(ascending) (Format: spaces as seperators, colnames and rownames present to file.
# 	Parameters	(required):
# 		-- inputFilename: name or full path to the TSV file
# 		-- ouputFilename: name or full path to the outputfile
#
# 	Parameters (optional):
#
# 	Returns: NULL
#
#   Warning: Temporary version, either I'm missing something with all the data type restrictions or this could be a lot shorter.
list2Matrix <- function(inputFilename = NULL, input = NULL) {
  # Read the table from file
  if (!is.null(inputFilename)) {
    tbl <- read.table(inputFilename)
  } else {
    tbl <- input
  }
  # Sort by second column
  tbl <- tbl[sort.list(tbl[, 2]), ]
  # Get the number of unique elements in the predictors
  pred <- unique(tbl[, 1])
  # Get the number of unique elements in the targets
  targ <- unique(tbl[, 2])
  # Pre allocate return matrix
  m <- matrix(0.0, length(pred), length(targ))
  # Set the col/row-names
  colnames(m) <- targ
  rownames(m) <- pred
  # Get the duplicates
  dups <- duplicated(tbl[, 2])
  # Get the starIndices of another column
  startIndices <- which(FALSE == dups)
  for (i in 1:(length(startIndices) - 1)) {
    targetToAdd <- tbl[startIndices[i], 2]
    predToAdd <- tbl[startIndices[i]:(startIndices[i + 1] - 1), 1]
    valuesToAdd <- tbl[startIndices[i]:(startIndices[i + 1] - 1), 3]
    colIndex <- which(colnames(m) %in% targetToAdd)
    tmp <- c()
    for (i in predToAdd) {
      tmp <- c(tmp, which(rownames(m) == i))
    }
    rowIndexes <- tmp
    m[rowIndexes, colIndex] <- valuesToAdd
  }
  targetToAdd <- tbl[startIndices[length(startIndices)], 2]
  predToAdd <- tbl[startIndices[length(startIndices)]:length(tbl[, 1]), 1]
  valuesToAdd <- tbl[startIndices[length(startIndices)]:length(tbl[, 1]), 3]
  tmp <- c()
  for (i in predToAdd) {
    tmp <- c(tmp, which(rownames(m) == i))
  }
  rowIndexes <- tmp
  colIndex <- which(colnames(m) %in% targetToAdd)
  m[rowIndexes, colIndex] <- valuesToAdd
  return(m)
}


#' isPercentage
#'
#' @param number
#'
#' @return
#' @export
#'
#' @examples
isPercentage <- function(number) {
  if (!is.numeric(number) || number <= 0 || number > 100 || length(number) != 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}


#' getIndicesOfGenesInMatrix
#'  Returns and checks the column indices to be used as targets or predictors.
#'
#' @param expressionMatrix Expression matrix as returned by constructExpressionMatrixFromFile.
#'                         Format: samples x genes ; colnames are unique gene labels ; rownames are ignored
#' @param genes [DEFAULT = NULL] Vector indicating the genes whose column indices need to be returned.
#'                               Parameters genes and inputFile can not be both non-null.
#'                               Format: 1- A vector of indices, 2- A vector of gene names (strings)
#' @param inputFile [DEFAULT = NULL] Path to a file containing a tab or new line seperated file with the genes specified as in inputfileFormat.
#' @param seperator [DEFAULT = NULL]
#' @param inputFileFormat  [DEFAULT = character()] Indicates whether the genes are specified as strings or indices in the file.
#'
#' @return Sorted indices list of columns
#' @export
#'
#' @examples
getIndicesOfGenesInMatrix <- function(expressionMatrix,
                                      genes = NULL,
                                      inputFile = NULL,
                                      seperator = "\t",
                                      inputFileFormat = character()) {
  # Check for conflicting arguments
  if ((!is.null(genes)) && (!is.null(inputFile))) {
    stop("Aborting. Parameters \"genes\" and \"inputFile\" cannot be both non-null")
  }
  # Read the genes from file if requested
  if (!is.null(inputFile)) {
    genes <- (scan(inputFile, character(), sep = seperator, what = inputFileFormat))
  }
  # Check if genes are provided, if not all genes are considered
  if (is.null(genes)) {
    genes <- 1:length(colnames(expressionMatrix))
  } else {
    # Check if indices or gene names are given
    if (!isPositiveNonZeroIntegerVector(genes)) { # Do conversions if genes are specified by gene name
      # Check first if all geneTargets are in the expressionMatrix
      tmp <- setdiff(genes, colnames(expressionMatrix))
      if (length(tmp) > 0) {
        for (i in tmp) {
          warning(paste0("Gene ", i, " not in expression matrix\n"))
        }
        warning("Certain genes not present in expression matrix")
      }
      # If correct, convert into indices array
      tmp <- mat.or.vec(length(genes), 1)
      for (i in 1:length(genes)) {
        tmp[i] <- match(genes[i], colnames(expressionMatrix))
      }

      tmp <- tmp[!is.na(tmp)]
      tmp <- tmp[sort.list(as.numeric(tmp))]
    } else { # already an indices array, just sort
      genes <- sort(geneTargets)
      # Check for invalid indexes
      if (genes[length(genes)] > length(colnames(expressionMatrix))) {
        warning("Target index exceeds amount of genes/columns.")
      }
    }
    # Remove doubles from the list
    genes <- unique(genes)
  }
  return(genes)
}

# 	isPositiveNonZeroIntegerVector - checks if a vector is entirely composed of positive non zero integers
isPositiveNonZeroIntegerVector <- function(integerVector) {
  return(isIntegerVector(integerVector) && !any(integerVector <= 0))
}

# Sorts a matrix according to col and rownames
sortMatrix <- function(mat) {
  mat <- mat[sort.list(rownames(mat)), sort.list(colnames(mat))]
  return(mat)
}
