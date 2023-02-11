

source("Functions-basic_evaluation.R")
source("Functions-elasticNet.R")
source("Functions-genie3.R")
source("Functions-svm.R")
source("Functions-io.R")
source("Functions-util.R")

#' NIMEFI
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
#' @param GENIE [DEFAULT=TRUE] Indicates if GENIE3 should be run and be considered in the ensemble
#' @param treesInEnsembleCount [DEFAULT=1000] GENIE parameter (see. genie3.R)
#' @param candidateSplitCount [DEFAULT="ALL"] GENIE parameter (see. genie3.R)
#' @param GENIETrace [DEFAULT=TRUE] Indicates if progress should be reported during execution of GENIE
#' @param SVM [DEFAULT=TRUE] Indicates if SVM should be run and be considered in the ensemble
#' @param SVMPredSampleMin [DEFAULT=20] Lower bound of the size of the predictor sample, expressed in %
#' @param SVMPredSampleMax [DEFAULT=80 ] Upper bound of the size of the predictor sample, expressed in %
#' @param SVMExpSampleMin [DEFAULT=20] Lower bound of the size of the experiment sample, expressed in %
#' @param SVMExpSampleMax [DEFAULT=80] Upper bound of the size of the experiment sample, expressed in %
#' @param SVMRankThreshold [DEFAULT=5]  Amount of top predictors at the top of ranking which are awarded a score, expressed in %
#' @param SVMEnsembleSize [DEFAULT=2000 ] Amount of iterations
#' @param SVMTrace [DEFAULT=TRUE] Indicates if progress should be reported during execution of SVM
#' @param EL [DEFAULT=TRUE] Indicates if EL should be run and be considered in the ensemble
#' @param ELTrace [DEFAULT=TRUE] Indicates if progress should be reported during execution of EL
#' @param ELPredSampleMin [DEFAULT=20] : Lower bound of the size of the predictor sample, expressed in %
#' @param ELPredSampleMax [DEFAULT=80] : Upper bound of the size of the predictor sample, expressed in %
#' @param ELExpSampleMin [DEFAULT=20] : Lower bound of the size of the experiment sample, expressed in %
#' @param ELExpSampleMax [DEFAULT=80] : Upper bound of the size of the experiment sample, expressed in %
#' @param ELRankThreshold [DEFAULT=5] Amount of top predictors at the top of ranking which are awarded a score, expressed in %
#' @param ELEnsembleSize [DEFAULT=2000] Amount of iterations
#' @param outputFileName [DEFAULT="output.txt"] The result will be written to this file.
#' @param ... The other parameters are passed to the specific method to be used.
#'
#' @return A matrix where each cell (row i, column j) specifies the importance of variable i to target j.
#' @export
#'
#' @examples
NIMEFI <- function(expressionMatrix,
                   predictorIndices = NULL,
                   targetIndices = NULL,
                   GENIE = TRUE,
                   treesInEnsembleCount = 1000,
                   candidateSplitCount = "ALL",
                   GENIETrace = TRUE,
                   SVM = TRUE,
                   SVMPredSampleMin = 20,
                   SVMPredSampleMax = 80,
                   SVMExpSampleMin = 20,
                   SVMExpSampleMax = 80,
                   SVMRankThreshold = 5,
                   SVMEnsembleSize = 600,
                   SVMTrace = TRUE,
                   EL = TRUE,
                   ELTrace = TRUE,
                   ELPredSampleMin = 20,
                   ELPredSampleMax = 80,
                   ELExpSampleMin = 20,
                   ELExpSampleMax = 80,
                   ELRankThreshold = 5,
                   ELEnsembleSize = 600,
                   outputFileName = "output.txt",
                   ...) {
  # Check arguments
  # The following arguments should be specified as percentages
  if (!(isPercentage(SVMPredSampleMin) && isPercentage(SVMPredSampleMax) &&
    isPercentage(SVMExpSampleMax) && isPercentage(SVMExpSampleMin) && isPercentage(SVMRankThreshold) &&
    isPercentage(ELPredSampleMin) && isPercentage(ELPredSampleMax) &&
    isPercentage(ELExpSampleMin) && isPercentage(ELExpSampleMax) && isPercentage(ELRankThreshold))) {
    stop("Sample sizes and rank thresholds should be specified as percentages.")
  }
  # Ensemble size should be positive and non zero
  if (!(isWholeNumber(ELEnsembleSize) && isWholeNumber(SVMEnsembleSize) && ELEnsembleSize > 0 && SVMEnsembleSize > 0)) {
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
  ELPredSampleMin <- round((ELPredSampleMin * lenPredictors) / 100)
  ELPredSampleMax <- round((ELPredSampleMax * lenPredictors) / 100)
  ELExpSampleMin <- round((ELExpSampleMin * lenSamples) / 100)
  ELExpSampleMax <- round((ELExpSampleMax * lenSamples) / 100)
  ELRankThreshold <- round((ELRankThreshold * lenPredictors) / 100)

  # Call the algorithms
  genieResult <- NULL
  elResult <- NULL
  svmResult <- NULL

  # Preallocate the return matrix
  resultMatrix <- matrix(0.0, length(predictorIndices), length(targetIndices))
  rownames(resultMatrix) <- colnames(expressionMatrix)[predictorIndices]
  colnames(resultMatrix) <- colnames(expressionMatrix)[targetIndices]
  resultMatrix <- sortMatrix(resultMatrix)

  if (SVM) {
    message("Starting E-SVM network inference...")
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
      ...
    )
    svmResult <- convertAdjMatrixToSortedRankTSV(inputFile = svmResult)
    svmResult[, 3] <- 1:(dim(svmResult)[1])
    svmResult <- convertSortedRankTSVToAdjMatrix(input = svmResult)
    svmResult <- sortMatrix(svmResult)
    resultMatrix <- resultMatrix + svmResult
  }

  if (EL) {
    message("Starting E-EL network inference...")
    elResult <- ElVariableEnsembleSolve(expressionMatrix, predictorIndices, targetIndices,
      trace = ELTrace, ensembleSize = ELEnsembleSize, rankThreshold = ELRankThreshold,
      predictorSampleSizeMin = ELPredSampleMin, predictorSampleSizeMax = ELPredSampleMax,
      minSampleSize = ELExpSampleMin, maxSampleSize = ELExpSampleMax, ...
    )
    elResult <- convertAdjMatrixToSortedRankTSV(inputFile = elResult)
    elResult[, 3] <- 1:(dim(elResult)[1])
    elResult <- convertSortedRankTSVToAdjMatrix(input = elResult)
    elResult <- sortMatrix(elResult)
    resultMatrix <- resultMatrix + elResult
  }

  if (GENIE) {
    message("Starting GENIE3 network inference...")
    genieResult <- genieSolve(expressionMatrix,
      predictorIndices,
      targetIndices,
      trace = GENIETrace,
      candidateSplitCount = candidateSplitCount,
      treesInEnsembleCount = treesInEnsembleCount,
      ...
    )
    genieResult <- convertAdjMatrixToSortedRankTSV(inputFile = genieResult)
    genieResult[, 3] <- 1:(dim(genieResult)[1])
    genieResult <- convertSortedRankTSVToAdjMatrix(input = genieResult)
    genieResult <- sortMatrix(genieResult)
    resultMatrix <- resultMatrix + genieResult
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
      quote = FALSE
    )
  }
}
