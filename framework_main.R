

source("elasticNet.R")
source("genie3.R")
source("svm.R")
source("io.R")
source("util.R")
source("basic_evaluation.R")

# 	NIMEFI
# 	Constructs regression problems from an expression matrix using the columns specified as targets and solves the regression problems
# 	using ensemble feature selection techniques

# 	Parameters	(required):
# 		-- expressionMatrix: expression matrix as returned by constructExpressionMatrixFromFile.
# 							 Format: samples x genes ; colnames are unique gene labels ; rownames are ignored

#
# 	Parameters (optional):
# 		-- predictorIndices: [DEFAULT=NULL]: indicates the rows of the matrix that should be considered as predictors for a regression problem.
# 										 It is highly recommened to use the output of getPredictorIndices
# 										 If null, getIndicesOfGenesInMatrixs will be called on the matrix and all rows will be considered predictors if needed for the algorithm.
# 		-- targetIndices [DEFAULT=NULL]: indicates which genes should be considered as targets for a regression problem.
# 										 It is highly recommened to use the output of getIndicesOfGenesInMatrixs
# 										 If null, getIndidicesForRegressionProblems will be called on the matrix and all columns will be considered targets if needed for the algorithm.
# 		-- GENIE [DEFAULT=TRUE]        : indicates if GENIE3 should be run and be considered in the ensemble
# 		-- SVM [DEFAULT=TRUE]                    : indicates if SVM should be run and be considered in the ensemble
# 		-- EL [DEFAULT=TRUE]                     : indicates if EL  should be run and be considered in the ensemble
#       -- GENIETrace [DEFAULT=TRUE]             : indicates if progress should be reported during execution of GENIE
#       -- SVMTrace [DEFAULT=TRUE]               : indicates if progress should be reported during execution of SVM
#       -- ELTrace [DEFAULT=TRUE]                : indicates if progress should be reported during execution of EL
#       -- candidateSplitCount [DEFAULT="ALL]    : GENIE parameter (see. genie3.R)
#       -- treesInEnsembleCount [ DEFAULT: 1000] :  GENIE parameter (see. genie3.R)
#       -- SVMPredSampleMin [DEFAULT=20] : lower bound of the size of the predictor sample, expressed in %
#       -- SVMPredSampleMax [DEFAULT=80 ] : upper bound of the size of the predictor sample, expressed in %
#       -- SVMExpSampleMin [DEFAULT=20] : lower bound of the size of the experiment sample, expressed in %
#       -- SVMExpSampleMax [DEFAULT=80] : upper bound of the size of the experiment sample, expressed in %
#       -- SVMRankThreshold [DEFAULT=5] :  amount of top predictors at the top of ranking which are awarded a score, expressed in %
#       -- SVMEnsembleSize [DEFAULT=2000 ] : amount of iterations
#       -- ELPredSampleMin [DEFAULT=20] : lower bound of the size of the predictor sample, expressed in %
#       -- ELPredSampleMax [DEFAULT=80] : upper bound of the size of the predictor sample, expressed in %
#       -- ELExpSampleMin  [DEFAULT=20] : lower bound of the size of the experiment sample, expressed in %
#       -- ELExpSampleMax [DEFAULT=80] : upper bound of the size of the experiment sample, expressed in %
#       -- ELRankThreshold [DEFAULT=5] : amount of top predictors at the top of ranking which are awarded a score, expressed in %
#       -- ELEnsembleSize [DEFAULT=2000] :amount of iterations
# 		-- outputFileName [DEFAULT: "output.txt"] : The result will be written to this file.
#
# 		-- ... : the other parameters are passed to the specific method to be used.
#
# 	Returns:
#
# 		A matrix where each cell(row i, column j) specifies the importance of variable i to target j.
#

NIMEFI <- function(expressionMatrix,
                   predictorIndices = NULL,
                   targetIndices = NULL,
                   GENIE = TRUE,
                   treesInEnsembleCount = 1000,
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
                   candidateSplitCount = "ALL",
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
    print("Starting E-SVM network inference...")
    svmResult <- SVMVariableEnsembleSolve(expressionMatrix, predictorIndices, targetIndices,
      trace = SVMTrace, ensembleSize = SVMEnsembleSize, rankThreshold = SVMRankThreshold,
      predictorSampleSizeMin = SVMPredSampleMin, predictorSampleSizeMax = SVMPredSampleMax,
      minSampleSize = SVMExpSampleMin, maxSampleSize = SVMExpSampleMax, ...
    )
    svmResult <- convertAdjMatrixToSortedRankTSV(inputFile = svmResult)
    svmResult[, 3] <- 1:(dim(svmResult)[1])
    svmResult <- convertSortedRankTSVToAdjMatrix(input = svmResult)
    svmResult <- sortMatrix(svmResult)
    resultMatrix <- resultMatrix + svmResult
  }

  if (EL) {
    print("Starting E-EL network inference...")
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
    print("Starting GENIE3 network inference...")
    genieResult <- genieSolve(expressionMatrix, 
                              predictorIndices, 
                              targetIndices, 
                              trace = GENIETrace, 
                              candidateSplitCount = candidateSplitCount, 
                              treesInEnsembleCount = treesInEnsembleCount,
                              ...)
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
    write.table(resultMatrix, outputFileName, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}
