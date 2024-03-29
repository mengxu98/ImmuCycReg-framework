library(magrittr)
source("functions/Functions-basic_evaluation.R")
source("functions/Functions-L0.R")
source("functions/Functions-svm.R")
source("functions/Functions-io.R")
source("functions/Functions-util.R")

#' L0REG
#' 	Constructs regression problems from an expression matrix using the columns specified as targets
#' 	and solves the regression problems using ensemble feature selection techniques
#'
#' @param expressionMatrix Expression matrix as returned by constructExpressionMatrixFromFile.
#'                         Format: samples x genes; colnames are unique gene labels; rownames are ignored.
#' @param predictorIndices [DEFAULT=NULL] Indicates the rows of the matrix that should be considered as predictors
#'                                        for a regression problem. It is highly recommened to use the output of
#'                                        getPredictorIndices. If null, getIndicesOfGenesInMatrixs will be called
#'                                        on the matrix and all rows will be considered predictors if needed for the algorithm.
#' @param targetIndices [DEFAULT=NULL] Indicates which genes should be considered as targets for a regression problem.
#'                                     It is highly recommened to use the output of getIndicesOfGenesInMatrixs.
#'                                     If null, getIndidicesForRegressionProblems will be called on the matrix and
#'                                     All columns will be considered targets if needed for the algorithm.
#' @param SVM [DEFAULT=TRUE] Indicates if SVM should be run and be considered in the ensemble
#' @param SVMPredSampleMin [DEFAULT=20] Lower bound of the size of the predictor sample, expressed in %
#' @param SVMPredSampleMax [DEFAULT=80 ] Upper bound of the size of the predictor sample, expressed in %
#' @param SVMExpSampleMin [DEFAULT=20] Lower bound of the size of the experiment sample, expressed in %
#' @param SVMExpSampleMax [DEFAULT=80] Upper bound of the size of the experiment sample, expressed in %
#' @param SVMRankThreshold [DEFAULT=5]  Amount of top predictors at the top of ranking which are awarded a score, expressed in %
#' @param SVMEnsembleSize [DEFAULT=2000] Amount of iterations
#' @param SVMTrace [DEFAULT=TRUE] Indicates if progress should be reported during execution of SVM
#' @param L0 [DEFAULT=TRUE] Indicates if L0 should be run and be considered in the ensemble
#' @param L0Trace [DEFAULT=TRUE] Indicates if progress should be reported during execution of L0
#' @param L0PredSampleMin [DEFAULT=20] : Lower bound of the size of the predictor sample, expressed in %
#' @param L0PredSampleMax [DEFAULT=80] : Upper bound of the size of the predictor sample, expressed in %
#' @param L0ExpSampleMin [DEFAULT=20] : Lower bound of the size of the experiment sample, expressed in %
#' @param L0ExpSampleMax [DEFAULT=80] : Upper bound of the size of the experiment sample, expressed in %
#' @param L0RankThreshold [DEFAULT=5] Amount of top predictors at the top of ranking which are awarded a score, expressed in %
#' @param L0EnsembleSize [DEFAULT=2000] Amount of iterations
#' @param outputFileName [DEFAULT="output.txt"] The result will be written to this file.
#' @param ... The other parameters are passed to the specific method to be used.
#'
#' @return A matrix where each cell (row i, column j) specifies the importance of variable i to target j.
#' @export
#'
L0REG <- function(expressionMatrix,
                  predictorIndices = NULL,
                  targetIndices = NULL,
                  SVM = TRUE,
                  SVMPredSampleMin = 20,
                  SVMPredSampleMax = 80,
                  SVMExpSampleMin = 20,
                  SVMExpSampleMax = 80,
                  SVMRankThreshold = 5,
                  SVMEnsembleSize = 600,
                  SVMTrace = TRUE,
                  L0 = TRUE,
                  L0Trace = TRUE,
                  L0PredSampleMin = 20,
                  L0PredSampleMax = 80,
                  L0ExpSampleMin = 20,
                  L0ExpSampleMax = 80,
                  L0RankThreshold = 5,
                  L0EnsembleSize = 600,
                  outputFileName = "output.txt",
                  ...) {
  # Check arguments
  # The following arguments should be specified as percentages
  if (!(isPercentage(SVMPredSampleMin) && isPercentage(SVMPredSampleMax) &&
        isPercentage(SVMExpSampleMax) && isPercentage(SVMExpSampleMin) && isPercentage(SVMRankThreshold) &&
        isPercentage(L0PredSampleMin) && isPercentage(L0PredSampleMax) &&
        isPercentage(L0ExpSampleMin) && isPercentage(L0ExpSampleMax) && isPercentage(L0RankThreshold))) {
    stop("Sample sizes and rank thresholds should be specified as percentages.")
  }
  # Ensemble size should be positive and non zero
  if (!(isWholeNumber(L0EnsembleSize) && isWholeNumber(SVMEnsembleSize) && L0EnsembleSize > 0 && SVMEnsembleSize > 0)) {
    stop("Ensemble sizes should be specified as non-zero integers.")
  }
  
  # Set the predictorIndices to all if they are not set
  if (is.null(predictorIndices)) {
    predictorIndices <- getIndicesOfGenesInMatrix(expressionMatrix)
  }
  # Set targetIndices to all if they are are not set
  if (is.null(targetIndices)) {
    targetIndices <- getIndicesOfGenesInMatrix(expressionMatrix)
  }
  
  # Convert the percentages to the numbers
  lenPredictors <- dim(expressionMatrix)[2]
  lenSamples <- dim(expressionMatrix)[1]
  
  SVMPredSampleMin <- round((SVMPredSampleMin * lenPredictors) / 100)
  SVMPredSampleMax <- round((SVMPredSampleMax * lenPredictors) / 100)
  SVMExpSampleMin <- round((SVMExpSampleMin * lenSamples) / 100)
  SVMExpSampleMax <- round((SVMExpSampleMax * lenSamples) / 100)
  SVMRankThreshold <- round((SVMRankThreshold * lenPredictors) / 100)
  
  L0PredSampleMin <- round((L0PredSampleMin * lenPredictors) / 100)
  L0PredSampleMax <- round((L0PredSampleMax * lenPredictors) / 100)
  L0ExpSampleMin <- round((L0ExpSampleMin * lenSamples) / 100)
  L0ExpSampleMax <- round((L0ExpSampleMax * lenSamples) / 100)
  L0RankThreshold <- round((L0RankThreshold * lenPredictors) / 100)
  
  # Call the algorithms
  L0Result <- NULL
  svmResult <- NULL
  
  # Preallocate the return matrix
  resultMatrix <- matrix(0.0, length(predictorIndices), length(targetIndices))
  rownames(resultMatrix) <- colnames(expressionMatrix)[predictorIndices]
  colnames(resultMatrix) <- colnames(expressionMatrix)[targetIndices]
  resultMatrix <- sortMatrix(resultMatrix)
  
  if (SVM) {
    message("Starting E-SVM network inference......")
    svmResult <- SVMVariableEnsembleSolve(expressionMatrix,
                                          predictorIndices,
                                          targetIndices,
                                          trace = SVMTrace,
                                          ensembleSize = SVMEnsembleSize,
                                          rankThreshold = SVMRankThreshold,
                                          predictorSampleSizeMin = SVMPredSampleMin,
                                          predictorSampleSizeMax = SVMPredSampleMax,
                                          minSampleSize = SVMExpSampleMin,
                                          maxSampleSize = SVMExpSampleMax,
                                          ...)
    svmResult <- convertAdjMatrixToSortedRankTSV(inputFile = svmResult)
    svmResult[, 3] <- 1:(dim(svmResult)[1])
    svmResult <- convertSortedRankTSVToAdjMatrix(input = svmResult)
    svmResult <- sortMatrix(svmResult)
    resultMatrix <- resultMatrix + svmResult
  }
  
  if (L0) {
    message("Starting L0 inference......")
    L0Result <- L0VariableEnsembleSolve(expressionMatrix,
                                        predictorIndices,
                                        targetIndices,
                                        trace = L0Trace,
                                        ensembleSize = L0EnsembleSize,
                                        rankThreshold = L0RankThreshold,
                                        predictorSampleSizeMin = L0PredSampleMin,
                                        predictorSampleSizeMax = L0PredSampleMax,
                                        minSampleSize = L0ExpSampleMin,
                                        maxSampleSize = L0ExpSampleMax, ...)
    L0Result <- convertAdjMatrixToSortedRankTSV(inputFile = L0Result)
    L0Result[, 3] <- 1:(dim(L0Result)[1])
    L0Result <- convertSortedRankTSVToAdjMatrix(input = L0Result)
    L0Result <- sortMatrix(L0Result)
    resultMatrix <- resultMatrix + L0Result
  }
  
  resultMatrix <- convertAdjMatrixToSortedRankTSV(inputFile = resultMatrix, desc = FALSE)
  
  if (is.null(outputFileName)) {
    return(resultMatrix)
  } else {
    write.table(resultMatrix,
                outputFileName,
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE)
  }
}
