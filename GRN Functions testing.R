
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
      if (T) {
        doParallel::registerDoParallel(cores=6)
        message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
        "%dopar%"<- foreach::"%dopar%"
        suppressPackageStartupMessages(fit <- doRNG::"%dorng%"(foreach::foreach(X=X,Y=Y,
                                                                                penalty = penalty,
                                                                                maxSuppSize = maxSuppSize,
                                                                                nFolds = 10,
                                                                                seed = 1,
                                                                                nGamma = 5,
                                                                                gammaMin = 0.0001,
                                                                                gammaMax = 10),
                                                               {
                                                                 L0Learn::L0Learn.cvfit(X, Y,
                                                                                        penalty = penalty,
                                                                                        maxSuppSize = maxSuppSize,
                                                                                        nFolds = 10,
                                                                                        seed = 1,
                                                                                        nGamma = 5,
                                                                                        gammaMin = 0.0001,
                                                                                        gammaMax = 10
                                                                 )
                                                                 
                                                               }))
        attr(fit, "rng") <- NULL
        attr(fit,"doRNG_version") <- NULL
        fit <- unlist(fit, recursive=FALSE)
      }else{
        fit <- L0Learn.cvfit(X, Y,
                             penalty = penalty,
                             maxSuppSize = maxSuppSize,
                             nFolds = 10,
                             seed = 1,
                             nGamma = 5,
                             gammaMin = 0.0001,
                             gammaMax = 10
        )
      }
      
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
      temp <- coef(fit, lambda = lambda, gamma = gamma)
    },
    error = function(e) {
      if (T) {
        doParallel::registerDoParallel(cores=6)
        message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
        "%dopar%"<- foreach::"%dopar%"
        suppressPackageStartupMessages(fit <- doRNG::"%dorng%"(foreach::foreach(Y),
                                                               {
                                                                 L0Learn::L0Learn.fit(X, Y,
                                                                                      penalty = penalty,
                                                                                      maxSuppSize = maxSuppSize,
                                                                                      nGamma = 5,
                                                                                      gammaMin = 0.0001,
                                                                                      gammaMax = 10
                                                                 )
                                                                 
                                                               }))
        attr(fit, "rng") <- NULL
        attr(fit, "doRNG_version") <- NULL
        fit <- unlist(fit, recursive=FALSE)
      }else{
        fit <- L0Learn.fit(X, Y,
                           penalty = penalty,
                           maxSuppSize = maxSuppSize,
                           nGamma = 5,
                           gammaMin = 0.0001,
                           gammaMax = 10
        )
      }
      
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
      fit <- L0Learn.cvfit(X, Y,
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
      temp <- coef(fit, lambda = lambda, gamma = gamma)
    },
    error = function(e) {
      fit <- L0Learn.fit(X, Y,
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

#' Title
#'
#' @param matrix The rows are samples and the columns are genes of the matrix
#' @param penalty
#' @param regulators
#' @param targets
#' @param maxSuppSize
#'
#' @return
#' @export
#'
#' @examples
L0DWGRN <- function(matrix,
                    penalty = NULL,
                    regulators = NULL,
                    targets = NULL,
                    maxSuppSize = NULL) {
  library(L0Learn)
  matrix <- as.data.frame(t(matrix))
  weightdf <- c()
  if (is.null(penalty)) {
    penalty <- "L0"
  }
  if (is.null(maxSuppSize)) {
    maxSuppSize <- dim(matrix)[2]
  }
  if (is.null(targets)) {
    targets <- colnames(matrix)
  }
  if (!is.null(regulators)) {
    matrix_reg <- matrix[, regulators]
    for (i in 1:length(targets)) {
      if (targets[i] %in% regulators) {
        X <- as.matrix(matrix_reg[, -which(colnames(matrix_reg) == targets[i])])
      } else {
        X <- as.matrix(matrix_reg)
      }
      Y <- matrix[, targets[i]]
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
      # wghts <- wghts / max(wghts)
      if (F) {
        wghts <- wghts / max(wghts)
        indices <- sort.list(wghts, decreasing = TRUE)
        zeros <- which(wghts <= 0.8)
        # wghts[1:length(wghts)] <- 1
        wghts[zeros] <- 0
      }
      if (length(wghts) != dim(X)[2]) {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = targets[i], weight = 0)
        # weightd <- data.frame(regulatoryGene = targets[i], targetGene = colnames(X), weight = 0)
      } else {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = targets[i], weight = wghts)
        # weightd <- data.frame(regulatoryGene = targets[i], targetGene = colnames(X), weight = wghts)
      }
      # weightd$weight <- weightd$weight / max(weightd$weight)
      weightdf <- rbind.data.frame(weightdf, weightd)
      if (i == length(regulators)) {
        weightdf <- weightdf[order(weightdf$weight, decreasing = TRUE), ]
      }
    }
  } else {
    
    regulators <- colnames(matrix)
    
    
    doParallel::registerDoParallel(cores=6)
    message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
    "%dopar%"<- foreach::"%dopar%"
    suppressPackageStartupMessages(weightdf <- doRNG::"%dorng%"(foreach::foreach(regulator = regulators),
                                                           {
                                                             X <- as.matrix(matrix[, -which(colnames(matrix) == regulator)])
                                                             Y <- matrix[, regulator]
                                                             
                                                             ####
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
                                                                   temp <- coef(fit, lambda = lambda, gamma = gamma)
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
                                                             ####
                                                             
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
                                                             # wghts <- wghts / max(wghts)
                                                             if (F) {
                                                               wghts <- wghts / max(wghts)
                                                               wghts <- wghts / sum(wghts)
                                                               indices <- sort.list(wghts, decreasing = TRUE)
                                                               zeros <- which(wghts <= 0.8)
                                                               # wghts[1:length(wghts)] <- 1
                                                               wghts[zeros] <- 0
                                                             }
                                                             if (length(wghts) != dim(X)[2]) {
                                                               weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulator, weight = 0)
                                                               # weightd <- data.frame(regulatoryGene = regulators[i], targetGene = colnames(X), weight = 0)
                                                             } else {
                                                               weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulator, weight = wghts)
                                                               # weightd <- data.frame(regulatoryGene = regulators[i], targetGene = colnames(X), weight = wghts)
                                                             }
                                                             # weightd$weight <- weightd$weight / max(weightd$weight)
                                                             weightdf <- rbind.data.frame(weightdf, weightd)
                                                             
                                                             # 后续需要一个排序函数处理regulators
                                                             # if (i == length(regulators)) {
                                                             #   weightdf <- weightdf[order(weightdf$weight, decreasing = TRUE), ]
                                                             # }
                                                             # setNames(list(setNames(weightdf, colnames(X))), regulator)
                                                           }))
    attr(weightdf, "rng") <- NULL
    attr(weightdf, "doRNG_version") <- NULL
    # weightlist <- unlist(weightdf, recursive=FALSE)
    # weightlist <- weightdf
    # weightdf <- c()
    # for (i in 1:length(weightlist)) {
    #   weightdf <- rbind.data.frame(weightdf, weightlist[[i]])
    # }
    
    # doParallel::stopImplicitCluster()
    
    
    # for (i in 1:length(regulators)) {
    #   X <- as.matrix(matrix[, -which(colnames(matrix) == regulators[i])])
    #   Y <- matrix[, regulators[i]]
    #   temp <- LO_fit(X, Y,
    #                  penalty = penalty,
    #                  nFolds = 10,
    #                  seed = 1,
    #                  maxSuppSize = maxSuppSize,
    #                  nGamma = 5,
    #                  gammaMin = 0.0001,
    #                  gammaMax = 10
    #   )
    #   temp <- as.vector(temp)
    #   wghts <- temp[-1]
    #   wghts <- abs(wghts)
    #   # wghts <- wghts / max(wghts)
    #   if (F) {
    #     wghts <- wghts / max(wghts)
    #     indices <- sort.list(wghts, decreasing = TRUE)
    #     zeros <- which(wghts <= 0.8)
    #     # wghts[1:length(wghts)] <- 1
    #     wghts[zeros] <- 0
    #   }
    #   if (length(wghts) != dim(X)[2]) {
    #     weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = 0)
    #     # weightd <- data.frame(regulatoryGene = regulators[i], targetGene = colnames(X), weight = 0)
    #   } else {
    #     weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = wghts)
    #     # weightd <- data.frame(regulatoryGene = regulators[i], targetGene = colnames(X), weight = wghts)
    #   }
    #   # weightd$weight <- weightd$weight / max(weightd$weight)
    #   weightdf <- rbind.data.frame(weightdf, weightd)
    #   if (i == length(regulators)) {
    #     weightdf <- weightdf[order(weightdf$weight, decreasing = TRUE), ]
    #   }
    # }
  }
  return(weightdf)
}

